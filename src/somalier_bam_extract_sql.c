/* Parallel sparse BAM/CRAM panel-site counting over one prepared panel. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <ctype.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/faidx.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include "include/bam_site_counts.h"
#include "duckhts_registration.h"

#define BAM_EXTRACT_MAX_PATH_BYTES 16384u
#define BAM_EXTRACT_MAX_IDENTITY_BYTES 1024u
#define BAM_EXTRACT_MAX_SITES UINT64_C(100000000)
#define BAM_EXTRACT_MAX_WORKERS 64u
#define BAM_EXTRACT_ERROR_BYTES 1024u
#define BAM_EXTRACT_PANEL_TABLE \
    "temp.main.__duckhts_somalier_bam_panel"

enum bam_extract_output_field {
    BAM_OUT_SAMPLE_ID = 0,
    BAM_OUT_SOURCE_PATH,
    BAM_OUT_ASSEMBLY,
    BAM_OUT_SITE_INDEX,
    BAM_OUT_REGION,
    BAM_OUT_POSITION,
    BAM_OUT_ALLELE_A,
    BAM_OUT_ALLELE_B,
    BAM_OUT_A,
    BAM_OUT_B,
    BAM_OUT_OTHER,
    BAM_OUT_SOURCE_METHOD,
    BAM_OUT_COUNT_SCOPE,
    BAM_OUT_OVERLAP_POLICY,
    BAM_OUT_STATUS,
    BAM_OUT_FIELD_COUNT
};

static const char *const bam_extract_output_names[BAM_OUT_FIELD_COUNT] = {
    "sample_id", "source_path", "assembly", "site_index", "region",
    "position", "allele_a", "allele_b", "a", "b", "other",
    "source_method", "count_scope", "overlap_policy", "status"
};

typedef enum bam_extract_site_status {
    BAM_SITE_MEASURED = 0,
    BAM_SITE_REFERENCE_CONTIG_UNAVAILABLE,
    BAM_SITE_REFERENCE_POSITION_UNAVAILABLE,
    BAM_SITE_REFERENCE_MISMATCH,
    BAM_SITE_ALIGNMENT_CONTIG_UNAVAILABLE,
    BAM_SITE_ALIGNMENT_POSITION_UNAVAILABLE
} bam_extract_site_status_t;

typedef struct bam_extract_registry {
    pthread_mutex_t query_mutex;
    duckdb_connection query_connection;
} bam_extract_registry_t;

static _Thread_local unsigned bam_extract_preparation_depth;

typedef struct bam_extract_bind {
    bam_extract_registry_t *registry;
    char *source_path;
    char *panel_table;
    char *panel_parquet;
    char *sample_id;
    char *reference_path;
    char *index_path;
    char *reference_index_path;
    char *cram_reference;
    uint32_t min_mapq;
    uint32_t min_baseq;
    uint32_t require_flags;
    uint32_t exclude_flags;
    uint32_t decompression_threads;
    uint32_t worker_count;
    uint32_t max_depth;
    uint32_t max_overlap_qnames;
    uint64_t max_sites;
    uint64_t max_region_bytes;
    uint64_t remote_block_bytes;
    uint64_t remote_cache_bytes;
    uint64_t reference_cache_bytes;
    duckhts_bam_site_overlap_policy_t overlap_policy;
} bam_extract_bind_t;

typedef struct bam_extract_site {
    uint64_t site_index;
    char *region;
    uint32_t region_length;
    uint64_t position1;
    char allele_a;
    char allele_b;
    int32_t tid;
    bam_extract_site_status_t status;
    duckhts_bam_site_counts_t counts;
} bam_extract_site_t;

typedef struct bam_extract_order {
    size_t site_offset;
    int32_t tid;
    hts_pos_t pos0;
} bam_extract_order_t;

typedef struct bam_extract_job {
    size_t begin;
    size_t end;
} bam_extract_job_t;

typedef struct bam_extract_global {
    bam_extract_site_t *sites;
    char *region_storage;
    char *assembly;
    uint32_t assembly_length;
    size_t site_count;
    bam_extract_job_t jobs[BAM_EXTRACT_MAX_WORKERS];
    uint32_t job_count;
    uint32_t next_job;
} bam_extract_global_t;

typedef struct bam_extract_local {
    size_t cursor;
    size_t end;
    uint8_t have_job;
    uint8_t done;
} bam_extract_local_t;

typedef struct bam_extract_reader {
    samFile *fp;
    hts_itr_t *iterator;
    uint32_t min_mapq;
    uint16_t require_flags;
    uint16_t exclude_flags;
    int io_error;
} bam_extract_reader_t;

static bool vector_row_valid(duckdb_vector vector, idx_t row) {
    uint64_t *validity = duckdb_vector_get_validity(vector);
    return validity == NULL || duckdb_validity_row_is_valid(validity, row);
}

static bool checked_size_product(size_t count, size_t width,
                                 size_t *product) {
    if (width != 0u && count > SIZE_MAX / width) return false;
    *product = count * width;
    return true;
}

static void copy_error(char *error, size_t error_size, const char *message) {
    if (error_size == 0u) return;
    snprintf(error, error_size, "%s", message ? message :
             "duckhts_somalier_bam_counts: unknown error");
}

static const char *site_status_name(bam_extract_site_status_t status) {
    switch (status) {
        case BAM_SITE_MEASURED:
            return "measured";
        case BAM_SITE_REFERENCE_CONTIG_UNAVAILABLE:
            return "unavailable_reference_contig";
        case BAM_SITE_REFERENCE_POSITION_UNAVAILABLE:
            return "unavailable_reference_position";
        case BAM_SITE_REFERENCE_MISMATCH:
            return "unavailable_reference_mismatch";
        case BAM_SITE_ALIGNMENT_CONTIG_UNAVAILABLE:
            return "unavailable_alignment_contig";
        case BAM_SITE_ALIGNMENT_POSITION_UNAVAILABLE:
            return "unavailable_alignment_position";
    }
    return "unknown";
}

static int compare_order(const void *left_ptr, const void *right_ptr) {
    const bam_extract_order_t *left = left_ptr;
    const bam_extract_order_t *right = right_ptr;
    if (left->tid != right->tid) return left->tid < right->tid ? -1 : 1;
    if (left->pos0 != right->pos0) return left->pos0 < right->pos0 ? -1 : 1;
    return 0;
}

static int bam_extract_fetch(void *context, bam1_t *record) {
    bam_extract_reader_t *reader = context;
    for (;;) {
        int status = sam_itr_next(reader->fp, reader->iterator, record);
        uint16_t flags;
        if (status < 0) {
            if (status < -1) reader->io_error = 1;
            return status;
        }
        flags = record->core.flag;
        if ((flags & reader->exclude_flags) != 0u) continue;
        if ((flags & reader->require_flags) != reader->require_flags) continue;
        if (record->core.qual < reader->min_mapq) continue;
        return status;
    }
}

static bool read_reference_base(const faidx_t *fai,
                                const bam_extract_site_t *site,
                                char *base,
                                bam_extract_site_status_t *status) {
    char buffer[4];
    hts_pos_t length = -1;
    hts_pos_t sequence_length;
    size_t required = 0u;

    if (!faidx_has_seq(fai, site->region)) {
        *status = BAM_SITE_REFERENCE_CONTIG_UNAVAILABLE;
        return true;
    }
    sequence_length = faidx_seq_len64(fai, site->region);
    if (sequence_length < 0 || (uint64_t)sequence_length < site->position1) {
        *status = BAM_SITE_REFERENCE_POSITION_UNAVAILABLE;
        return true;
    }
    if (faidx_fetch_seq64_into(fai, site->region,
                               (hts_pos_t)(site->position1 - 1u),
                               (hts_pos_t)(site->position1 - 1u), buffer,
                               sizeof(buffer), &length, &required) != 0) {
        if (length == -1) return false;
        *status = length == -2 ? BAM_SITE_REFERENCE_CONTIG_UNAVAILABLE
                               : BAM_SITE_REFERENCE_POSITION_UNAVAILABLE;
        return true;
    }
    if (length != 1) {
        *status = BAM_SITE_REFERENCE_POSITION_UNAVAILABLE;
        return true;
    }
    *base = (char)toupper((unsigned char)buffer[0]);
    return true;
}

static bool build_region_list(const bam_extract_order_t *order,
                              size_t count, uint64_t max_region_bytes,
                              hts_reglist_t **regions_out,
                              unsigned int *region_count_out,
                              const char **error) {
    hts_reglist_t *regions = NULL;
    size_t region_count = 0u;
    size_t required;

    if (count == 0u) {
        *regions_out = NULL;
        *region_count_out = 0u;
        return true;
    }
    if (count > INT_MAX || count > SIZE_MAX / sizeof(hts_pair_pos_t)) {
        *error = "duckhts_somalier_bam_counts: panel has too many indexed regions";
        return false;
    }
    for (size_t i = 0u; i < count; i++) {
        if (i == 0u || order[i].tid != order[i - 1u].tid) region_count++;
    }
    if (!checked_size_product(region_count, sizeof(*regions), &required) ||
        count > (SIZE_MAX - required) / sizeof(hts_pair_pos_t) ||
        required + count * sizeof(hts_pair_pos_t) > max_region_bytes) {
        *error = "duckhts_somalier_bam_counts: indexed region workspace exceeds max_region_bytes";
        return false;
    }
    regions = calloc(region_count, sizeof(*regions));
    if (!regions) {
        *error = "duckhts_somalier_bam_counts: out of memory allocating multi-region workspace";
        return false;
    }
    for (size_t begin = 0u, group = 0u; begin < count; group++) {
        size_t end = begin + 1u;
        size_t interval_count;

        while (end < count && order[end].tid == order[begin].tid) end++;
        interval_count = end - begin;
        regions[group].intervals = calloc(interval_count,
                                          sizeof(*regions[group].intervals));
        if (!regions[group].intervals) {
            *error = "duckhts_somalier_bam_counts: out of memory allocating region intervals";
            hts_reglist_free(regions, (int)region_count);
            return false;
        }
        regions[group].tid = order[begin].tid;
        regions[group].count = (uint32_t)interval_count;
        regions[group].min_beg = order[begin].pos0;
        regions[group].max_end = order[end - 1u].pos0 + 1;
        for (size_t i = begin; i < end; i++) {
            regions[group].intervals[i - begin].beg = order[i].pos0;
            regions[group].intervals[i - begin].end = order[i].pos0 + 1;
        }
        begin = end;
    }
    *regions_out = regions;
    *region_count_out = (unsigned int)region_count;
    return true;
}

static bool scan_sites(samFile *fp, sam_hdr_t *header, hts_idx_t *index,
                       bam_extract_site_t *sites,
                       const bam_extract_order_t *order, size_t available_count,
                       const bam_extract_bind_t *bind, const char **error) {
    bam_extract_reader_t reader = {0};
    bam_plp_t pileup = NULL;
    duckhts_bam_site_overlap_slot_t *overlap_slots = NULL;
    duckhts_bam_site_overlap_scratch_t overlap_scratch = {0};
    duckhts_bam_site_count_config_t count_config;
    hts_reglist_t *regions = NULL;
    unsigned int region_count = 0u;
    size_t cursor = 0u;
    bool ok = false;

    if (available_count == 0u) return true;
    if (!build_region_list(order, available_count, bind->max_region_bytes,
                           &regions, &region_count, error)) {
        return false;
    }
    reader.fp = fp;
    reader.min_mapq = bind->min_mapq;
    reader.require_flags = (uint16_t)bind->require_flags;
    reader.exclude_flags = (uint16_t)bind->exclude_flags;
    reader.iterator = sam_itr_regions(index, header, regions, region_count);
    if (!reader.iterator) {
        *error = "duckhts_somalier_bam_counts: could not construct one multi-region iterator";
        goto cleanup;
    }
    regions = NULL;
    if (bind->overlap_policy == DUCKHTS_BAM_SITE_OVERLAP_HILEUP_V0_1_0 &&
        bind->max_overlap_qnames != 0u) {
        size_t overlap_bytes;
        size_t slot_capacity;
        if (!duckhts_bam_site_overlap_slot_capacity(
                (size_t)bind->max_overlap_qnames, &slot_capacity) ||
            !checked_size_product(slot_capacity, sizeof(*overlap_slots),
                                  &overlap_bytes)) {
            *error = "duckhts_somalier_bam_counts: overlap scratch size overflows size_t";
            goto cleanup;
        }
        overlap_slots = duckdb_malloc(overlap_bytes);
        if (!overlap_slots) {
            *error = "duckhts_somalier_bam_counts: out of memory allocating overlap scratch";
            goto cleanup;
        }
        memset(overlap_slots, 0, overlap_bytes);
        overlap_scratch.slot_capacity = slot_capacity;
    }
    overlap_scratch.slots = overlap_slots;
    overlap_scratch.max_entries = bind->max_overlap_qnames;
    count_config.min_baseq = (uint8_t)bind->min_baseq;
    count_config.overlap_policy = bind->overlap_policy;

    pileup = bam_plp_init(bam_extract_fetch, &reader);
    if (!pileup) {
        *error = "duckhts_somalier_bam_counts: could not initialize pileup";
        goto cleanup;
    }
    /* htslib's active-node count includes its tail sentinel and drops a
     * same-start record only when that count is already greater than maxcnt.
     * Admit one observation beyond the public limit so exhaustion is visible. */
    bam_plp_set_maxcnt(pileup, (int)bind->max_depth + 1);
    for (;;) {
        int tid = -1;
        int depth = 0;
        hts_pos_t pos0 = -1;
        const bam_pileup1_t *observations =
            bam_plp64_auto(pileup, &tid, &pos0, &depth);
        duckhts_bam_site_status_t count_status;
        duckhts_bam_site_t kernel_site;

        if (!observations) {
            if (depth < 0 || reader.io_error) {
                *error = "duckhts_somalier_bam_counts: alignment or pileup read failed";
                goto cleanup;
            }
            break;
        }
        if ((uint32_t)depth > bind->max_depth) {
            *error = "duckhts_somalier_bam_counts: observed pileup exceeds max_depth";
            goto cleanup;
        }
        while (cursor < available_count &&
               (order[cursor].tid < tid ||
                (order[cursor].tid == tid && order[cursor].pos0 < pos0))) {
            cursor++;
        }
        if (cursor == available_count) break;
        if (order[cursor].tid != tid || order[cursor].pos0 != pos0) continue;

        kernel_site.pos0 = pos0;
        kernel_site.allele_a = sites[order[cursor].site_offset].allele_a;
        kernel_site.allele_b = sites[order[cursor].site_offset].allele_b;
        count_status = duckhts_bam_site_count_pileup(
            &kernel_site, observations, (size_t)depth, &count_config,
            bind->overlap_policy == DUCKHTS_BAM_SITE_OVERLAP_NONE
                ? NULL : &overlap_scratch,
            &sites[order[cursor].site_offset].counts);
        if (count_status != DUCKHTS_BAM_SITE_OK) {
            *error = duckhts_bam_site_status_string(count_status);
            goto cleanup;
        }
        cursor++;
    }
    ok = true;

cleanup:
    if (pileup) bam_plp_destroy(pileup);
    if (reader.iterator) hts_itr_destroy(reader.iterator);
    if (overlap_slots) duckdb_free(overlap_slots);
    if (regions) hts_reglist_free(regions, (int)region_count);
    return ok;
}

static bool set_remote_tuning(htsFile *fp, const char *path,
                              uint64_t block_bytes, uint64_t cache_bytes,
                              const char **error) {
    if (!hisremote(path)) return true;
    if ((block_bytes != 0u &&
         hts_set_opt(fp, HTS_OPT_BLOCK_SIZE, (int)block_bytes) < 0) ||
        (cache_bytes != 0u &&
         hts_set_opt(fp, HTS_OPT_CACHE_SIZE, (int)cache_bytes) < 0)) {
        *error = "duckhts_somalier_bam_counts: could not apply remote I/O settings";
        return false;
    }
    return true;
}

static char *cram_reference_locator(const char *reference_path,
                                    const char *reference_index_path,
                                    const char **error) {
    size_t reference_length;
    size_t index_length;
    size_t delimiter_length = sizeof(HTS_IDX_DELIM) - 1u;
    size_t required;
    char *locator;

    if (!reference_index_path) return NULL;
    reference_length = strlen(reference_path);
    index_length = strlen(reference_index_path);
    if (reference_length > SIZE_MAX - delimiter_length ||
        reference_length + delimiter_length > SIZE_MAX - index_length - 1u) {
        *error = "duckhts_somalier_bam_counts: CRAM reference locator is too large";
        return NULL;
    }
    required = reference_length + delimiter_length + index_length + 1u;
    locator = duckdb_malloc(required);
    if (!locator) {
        *error = "duckhts_somalier_bam_counts: out of memory building CRAM reference locator";
        return NULL;
    }
    memcpy(locator, reference_path, reference_length);
    memcpy(locator + reference_length, HTS_IDX_DELIM, delimiter_length);
    memcpy(locator + reference_length + delimiter_length,
           reference_index_path, index_length);
    locator[required - 1u] = '\0';
    return locator;
}

static void bam_extract_registry_destroy(void *pointer) {
    bam_extract_registry_t *registry = pointer;
    if (!registry) return;
    if (registry->query_connection) {
        duckdb_disconnect(&registry->query_connection);
    }
    pthread_mutex_destroy(&registry->query_mutex);
    duckdb_free(registry);
}

static void bam_extract_bind_destroy(void *pointer) {
    bam_extract_bind_t *bind = pointer;
    if (!bind) return;
    duckdb_free(bind->cram_reference);
    duckdb_free(bind->reference_index_path);
    duckdb_free(bind->index_path);
    duckdb_free(bind->reference_path);
    duckdb_free(bind->sample_id);
    duckdb_free(bind->panel_parquet);
    duckdb_free(bind->panel_table);
    duckdb_free(bind->source_path);
    duckdb_free(bind);
}

static void bam_extract_global_destroy(void *pointer) {
    bam_extract_global_t *global = pointer;
    if (!global) return;
    duckdb_free(global->assembly);
    duckdb_free(global->region_storage);
    duckdb_free(global->sites);
    duckdb_free(global);
}

static void bam_extract_local_destroy(void *pointer) {
    duckdb_free(pointer);
}

static char *copy_value_string(duckdb_value value, bool nullable,
                               size_t max_bytes, const char *name,
                               char *error, size_t error_size) {
    char *text;
    size_t length;
    if (!value || duckdb_is_null_value(value)) {
        if (nullable) return NULL;
        snprintf(error, error_size,
                 "duckhts_somalier_bam_counts: %s is required", name);
        return NULL;
    }
    text = duckdb_get_varchar(value);
    if (!text) {
        snprintf(error, error_size,
                 "duckhts_somalier_bam_counts: could not read %s", name);
        return NULL;
    }
    length = strlen(text);
    if (length == 0u || length > max_bytes) {
        snprintf(error, error_size,
                 "duckhts_somalier_bam_counts: %s must be 1..%zu bytes",
                 name, max_bytes);
        duckdb_free(text);
        return NULL;
    }
    return text;
}

static char *bind_positional_string(duckdb_bind_info info, idx_t index,
                                    bool nullable, size_t max_bytes,
                                    const char *name, char *error,
                                    size_t error_size) {
    duckdb_value value = duckdb_bind_get_parameter(info, index);
    char *text = copy_value_string(value, nullable, max_bytes, name,
                                   error, error_size);
    duckdb_destroy_value(&value);
    return text;
}

static char *bind_named_string(duckdb_bind_info info, const char *name,
                               bool nullable, size_t max_bytes,
                               char *error, size_t error_size) {
    duckdb_value value = duckdb_bind_get_named_parameter(info, name);
    char *text;
    if (!value) return NULL;
    text = copy_value_string(value, nullable, max_bytes, name,
                             error, error_size);
    duckdb_destroy_value(&value);
    return text;
}

static bool bind_named_uint64(duckdb_bind_info info, const char *name,
                              uint64_t default_value, uint64_t minimum,
                              uint64_t maximum, uint64_t *result,
                              char *error, size_t error_size) {
    duckdb_value value = duckdb_bind_get_named_parameter(info, name);
    uint64_t parsed;
    if (!value) {
        *result = default_value;
        return true;
    }
    if (duckdb_is_null_value(value)) {
        snprintf(error, error_size,
                 "duckhts_somalier_bam_counts: %s cannot be NULL", name);
        duckdb_destroy_value(&value);
        return false;
    }
    parsed = duckdb_get_uint64(value);
    duckdb_destroy_value(&value);
    if (parsed < minimum || parsed > maximum) {
        snprintf(error, error_size,
                 "duckhts_somalier_bam_counts: %s is outside [%" PRIu64
                 ", %" PRIu64 "]", name, minimum, maximum);
        return false;
    }
    *result = parsed;
    return true;
}

static bool bind_worker_count(duckdb_bind_info info, uint32_t *result,
                              char *error, size_t error_size) {
    duckdb_value value = duckdb_bind_get_named_parameter(info, "worker_count");
    double parsed;
    if (!value) {
        *result = 1u;
        return true;
    }
    if (duckdb_is_null_value(value)) {
        copy_error(error, error_size,
            "duckhts_somalier_bam_counts: worker_count cannot be NULL");
        duckdb_destroy_value(&value);
        return false;
    }
    parsed = duckdb_get_double(value);
    duckdb_destroy_value(&value);
    if (!isfinite(parsed) || parsed < 1.0 ||
        parsed > (double)BAM_EXTRACT_MAX_WORKERS || floor(parsed) != parsed) {
        copy_error(error, error_size,
            "duckhts_somalier_bam_counts: worker_count must be an integer in [1, 64]");
        return false;
    }
    *result = (uint32_t)parsed;
    return true;
}

static void add_result_column(duckdb_bind_info info, const char *name,
                              duckdb_type type) {
    duckdb_logical_type logical = duckdb_create_logical_type(type);
    duckdb_bind_add_result_column(info, name, logical);
    duckdb_destroy_logical_type(&logical);
}

static void bam_extract_bind(duckdb_bind_info info) {
    bam_extract_bind_t *bind = duckdb_malloc(sizeof(*bind));
    char error[BAM_EXTRACT_ERROR_BYTES] = {0};
    const char *locator_error = NULL;
    uint64_t value;

    if (!bind) {
        duckdb_bind_set_error(info,
            "duckhts_somalier_bam_counts: out of memory allocating bind state");
        return;
    }
    memset(bind, 0, sizeof(*bind));
    bind->registry = duckdb_bind_get_extra_info(info);
    bind->source_path = bind_positional_string(info, 0u, false,
        BAM_EXTRACT_MAX_PATH_BYTES, "source_path", error, sizeof(error));
    bind->panel_table = bind_positional_string(info, 1u, true,
        BAM_EXTRACT_MAX_PATH_BYTES, "panel_table", error, sizeof(error));
    bind->sample_id = bind_positional_string(info, 2u, false,
        BAM_EXTRACT_MAX_IDENTITY_BYTES, "sample_id", error, sizeof(error));
    bind->reference_path = bind_positional_string(info, 3u, false,
        BAM_EXTRACT_MAX_PATH_BYTES, "reference_path", error, sizeof(error));
    bind->panel_parquet = bind_named_string(info, "panel_parquet", true,
        BAM_EXTRACT_MAX_PATH_BYTES, error, sizeof(error));
    bind->index_path = bind_named_string(info, "index_path", true,
        BAM_EXTRACT_MAX_PATH_BYTES, error, sizeof(error));
    bind->reference_index_path = bind_named_string(info,
        "reference_index_path", true, BAM_EXTRACT_MAX_PATH_BYTES,
        error, sizeof(error));
    if (error[0] != '\0') goto fail;
    if ((bind->panel_table == NULL) == (bind->panel_parquet == NULL)) {
        copy_error(error, sizeof(error),
            "duckhts_somalier_bam_counts: provide exactly one of panel_table or panel_parquet");
        goto fail;
    }

    if (!bind_named_uint64(info, "min_mapq", 1u, 0u, UINT8_MAX,
                           &value, error, sizeof(error))) goto fail;
    bind->min_mapq = (uint32_t)value;
    if (!bind_named_uint64(info, "min_baseq", 0u, 0u, UINT8_MAX,
                           &value, error, sizeof(error))) goto fail;
    bind->min_baseq = (uint32_t)value;
    if (!bind_named_uint64(info, "require_flags", 0u, 0u, UINT16_MAX,
                           &value, error, sizeof(error))) goto fail;
    bind->require_flags = (uint32_t)value;
    if (!bind_named_uint64(info, "exclude_flags", 1796u, 0u, UINT16_MAX,
                           &value, error, sizeof(error))) goto fail;
    bind->exclude_flags = (uint32_t)value;
    if (!bind_named_uint64(info, "decompression_threads", 0u, 0u, INT_MAX,
                           &value, error, sizeof(error))) goto fail;
    bind->decompression_threads = (uint32_t)value;
    if (!bind_worker_count(info, &bind->worker_count,
                           error, sizeof(error))) goto fail;
    if (!bind_named_uint64(info, "max_depth", 100000u, 1u, INT_MAX - 1u,
                           &value, error, sizeof(error))) goto fail;
    bind->max_depth = (uint32_t)value;
    if (!bind_named_uint64(info, "max_overlap_qnames", 100000u, 0u,
                           UINT32_MAX, &value, error, sizeof(error))) goto fail;
    bind->max_overlap_qnames = (uint32_t)value;
    if (!bind_named_uint64(info, "max_sites", 1000000u, 1u,
                           BAM_EXTRACT_MAX_SITES, &bind->max_sites,
                           error, sizeof(error))) goto fail;
    if (!bind_named_uint64(info, "max_region_bytes", 67108864u, 1u,
                           SIZE_MAX, &bind->max_region_bytes,
                           error, sizeof(error))) goto fail;
    if (!bind_named_uint64(info, "remote_block_bytes", 1048576u, 0u,
                           INT_MAX, &bind->remote_block_bytes,
                           error, sizeof(error))) goto fail;
    if (!bind_named_uint64(info, "remote_cache_bytes", 67108864u, 0u,
                           INT_MAX, &bind->remote_cache_bytes,
                           error, sizeof(error))) goto fail;
    if (!bind_named_uint64(info, "reference_cache_bytes", 67108864u, 0u,
                           INT_MAX, &bind->reference_cache_bytes,
                           error, sizeof(error))) goto fail;
    {
        duckdb_value overlap_value = duckdb_bind_get_named_parameter(
            info, "overlap_policy");
        char *overlap = NULL;
        if (overlap_value && duckdb_is_null_value(overlap_value)) {
            copy_error(error, sizeof(error),
                "duckhts_somalier_bam_counts: overlap_policy cannot be NULL");
        } else if (overlap_value) {
            overlap = copy_value_string(overlap_value, false, 32u,
                "overlap_policy", error, sizeof(error));
        }
        if (overlap_value) duckdb_destroy_value(&overlap_value);
        if (error[0] != '\0') {
            duckdb_free(overlap);
            goto fail;
        }
        if (!overlap || strcmp(overlap, "hileup_v0.1.0") == 0) {
            bind->overlap_policy = DUCKHTS_BAM_SITE_OVERLAP_HILEUP_V0_1_0;
        } else if (strcmp(overlap, "none") == 0) {
            bind->overlap_policy = DUCKHTS_BAM_SITE_OVERLAP_NONE;
        } else {
            copy_error(error, sizeof(error),
                "duckhts_somalier_bam_counts: overlap_policy must be 'none' or 'hileup_v0.1.0'");
        }
        duckdb_free(overlap);
        if (error[0] != '\0') goto fail;
    }
    if (bind->reference_index_path) {
        bind->cram_reference = cram_reference_locator(bind->reference_path,
            bind->reference_index_path, &locator_error);
        if (!bind->cram_reference) {
            copy_error(error, sizeof(error), locator_error);
            goto fail;
        }
    }

    for (unsigned i = 0u; i < BAM_OUT_FIELD_COUNT; i++) {
        duckdb_type type = i == BAM_OUT_SITE_INDEX || i == BAM_OUT_POSITION ||
                           (i >= BAM_OUT_A && i <= BAM_OUT_OTHER)
            ? DUCKDB_TYPE_UBIGINT : DUCKDB_TYPE_VARCHAR;
        add_result_column(info, bam_extract_output_names[i], type);
    }
    duckdb_bind_set_bind_data(info, bind, bam_extract_bind_destroy);
    return;

fail:
    duckdb_bind_set_error(info, error[0] != '\0' ? error :
        "duckhts_somalier_bam_counts: invalid arguments");
    bam_extract_bind_destroy(bind);
}

static bool execute_prepared_source(duckdb_connection connection,
                                    const char *sql, const char *source,
                                    uint64_t row_limit,
                                    char *error, size_t error_size) {
    duckdb_prepared_statement statement = NULL;
    duckdb_result result = {0};
    bool ok = false;
    if (duckdb_prepare(connection, sql, &statement) != DuckDBSuccess) {
        copy_error(error, error_size, duckdb_prepare_error(statement));
        goto cleanup;
    }
    if (duckdb_bind_varchar(statement, 1u, source) != DuckDBSuccess ||
        duckdb_bind_uint64(statement, 2u, row_limit) != DuckDBSuccess ||
        duckdb_execute_prepared(statement, &result) != DuckDBSuccess) {
        copy_error(error, error_size, duckdb_result_error(&result));
        goto cleanup;
    }
    ok = true;
cleanup:
    duckdb_destroy_result(&result);
    if (statement) duckdb_destroy_prepare(&statement);
    return ok;
}

static bool execute_query(duckdb_connection connection, const char *sql,
                          duckdb_result *result, char *error,
                          size_t error_size) {
    if (duckdb_query(connection, sql, result) != DuckDBSuccess) {
        copy_error(error, error_size, duckdb_result_error(result));
        return false;
    }
    return true;
}

static bool load_panel_rows(duckdb_connection connection,
                            const bam_extract_bind_t *bind,
                            bam_extract_global_t *global,
                            char *error, size_t error_size) {
    static const char create_table[] =
        "CREATE OR REPLACE TEMP TABLE " BAM_EXTRACT_PANEL_TABLE
        " AS SELECT assembly, site_index, region, position, allele_a, allele_b "
        "FROM query_table(?) LIMIT ?";
    static const char create_parquet[] =
        "CREATE OR REPLACE TEMP TABLE " BAM_EXTRACT_PANEL_TABLE
        " AS SELECT assembly, site_index, region, position, allele_a, allele_b "
        "FROM read_parquet(?) LIMIT ?";
    static const char validate[] =
        "SELECT duckhts_somalier_panel_sha256('" BAM_EXTRACT_PANEL_TABLE "')";
    static const char stats[] =
        "SELECT count(*)::UBIGINT, "
        "coalesce(sum(octet_length(encode(region))), 0)::UBIGINT FROM "
        BAM_EXTRACT_PANEL_TABLE;
    static const char rows[] =
        "SELECT assembly::VARCHAR, site_index::UBIGINT, region::VARCHAR, "
        "position::UBIGINT, allele_a::VARCHAR, allele_b::VARCHAR FROM "
        BAM_EXTRACT_PANEL_TABLE " ORDER BY site_index";
    static const char drop[] = "DROP TABLE IF EXISTS " BAM_EXTRACT_PANEL_TABLE;
    duckdb_result result = {0};
    duckdb_data_chunk chunk = NULL;
    size_t region_bytes;
    size_t region_capacity;
    size_t region_offset = 0u;
    size_t site_offset = 0u;
    bool created = false;
    bool ok = false;

    if (!execute_prepared_source(connection,
            bind->panel_table ? create_table : create_parquet,
            bind->panel_table ? bind->panel_table : bind->panel_parquet,
            bind->max_sites + 1u, error, error_size)) {
        goto cleanup;
    }
    created = true;
    if (!execute_query(connection, stats, &result, error, error_size)) {
        goto cleanup;
    }
    chunk = duckdb_fetch_chunk(result);
    if (!chunk || duckdb_data_chunk_get_size(chunk) != 1u) {
        copy_error(error, error_size,
            "duckhts_somalier_bam_counts: could not read panel cardinality");
        goto cleanup;
    }
    {
        uint64_t site_count64 = ((uint64_t *)duckdb_vector_get_data(
            duckdb_data_chunk_get_vector(chunk, 0u)))[0];
        uint64_t region_bytes64 = ((uint64_t *)duckdb_vector_get_data(
            duckdb_data_chunk_get_vector(chunk, 1u)))[0];
        if (site_count64 > SIZE_MAX || region_bytes64 > SIZE_MAX) {
            copy_error(error, error_size,
                "duckhts_somalier_bam_counts: panel workspace exceeds this platform");
            goto cleanup;
        }
        global->site_count = (size_t)site_count64;
        region_bytes = (size_t)region_bytes64;
    }
    duckdb_destroy_data_chunk(&chunk);
    duckdb_destroy_result(&result);
    if (global->site_count == 0u || global->site_count > bind->max_sites ||
        global->site_count > SIZE_MAX / sizeof(*global->sites) ||
        region_bytes > SIZE_MAX - global->site_count) {
        copy_error(error, error_size,
            "duckhts_somalier_bam_counts: panel is empty, exceeds max_sites, or overflows workspace");
        goto cleanup;
    }
    if (!execute_query(connection, validate, &result, error, error_size)) {
        goto cleanup;
    }
    duckdb_destroy_result(&result);

    region_capacity = region_bytes + global->site_count;
    global->sites = duckdb_malloc(global->site_count * sizeof(*global->sites));
    global->region_storage = duckdb_malloc(region_capacity);
    if (!global->sites || !global->region_storage) {
        copy_error(error, error_size,
            "duckhts_somalier_bam_counts: out of memory preparing panel");
        goto cleanup;
    }
    memset(global->sites, 0, global->site_count * sizeof(*global->sites));
    if (!execute_query(connection, rows, &result, error, error_size)) {
        goto cleanup;
    }
    while ((chunk = duckdb_fetch_chunk(result)) != NULL) {
        duckdb_vector vectors[6];
        idx_t chunk_rows = duckdb_data_chunk_get_size(chunk);
        for (unsigned i = 0u; i < 6u; i++) {
            vectors[i] = duckdb_data_chunk_get_vector(chunk, i);
        }
        for (idx_t row = 0u; row < chunk_rows; row++) {
            duckdb_string_t *assembly;
            duckdb_string_t *region;
            duckdb_string_t *allele_a;
            duckdb_string_t *allele_b;
            bam_extract_site_t *site;
            uint32_t assembly_length;
            uint32_t region_length;
            size_t region_item;
            if (site_offset >= global->site_count) {
                copy_error(error, error_size,
                    "duckhts_somalier_bam_counts: panel cardinality changed during preparation");
                goto cleanup;
            }
            for (unsigned i = 0u; i < 6u; i++) {
                if (!vector_row_valid(vectors[i], row)) {
                    copy_error(error, error_size,
                        "duckhts_somalier_bam_counts: panel identity fields cannot be NULL");
                    goto cleanup;
                }
            }
            assembly = &((duckdb_string_t *)duckdb_vector_get_data(
                vectors[0]))[row];
            region = &((duckdb_string_t *)duckdb_vector_get_data(
                vectors[2]))[row];
            allele_a = &((duckdb_string_t *)duckdb_vector_get_data(
                vectors[4]))[row];
            allele_b = &((duckdb_string_t *)duckdb_vector_get_data(
                vectors[5]))[row];
            assembly_length = duckdb_string_t_length(*assembly);
            region_length = duckdb_string_t_length(*region);
            if (site_offset == 0u) {
                global->assembly = duckdb_malloc((size_t)assembly_length + 1u);
                if (!global->assembly) {
                    copy_error(error, error_size,
                        "duckhts_somalier_bam_counts: out of memory copying assembly");
                    goto cleanup;
                }
                memcpy(global->assembly, duckdb_string_t_data(assembly),
                       assembly_length);
                global->assembly[assembly_length] = '\0';
                global->assembly_length = assembly_length;
            }
            region_item = (size_t)region_length + 1u;
            if (region_offset > region_capacity ||
                region_item > region_capacity - region_offset) {
                copy_error(error, error_size,
                    "duckhts_somalier_bam_counts: panel region storage overflow");
                goto cleanup;
            }
            site = &global->sites[site_offset];
            site->site_index = ((uint64_t *)duckdb_vector_get_data(
                vectors[1]))[row];
            site->position1 = ((uint64_t *)duckdb_vector_get_data(
                vectors[3]))[row];
            site->region = global->region_storage + region_offset;
            site->region_length = region_length;
            memcpy(site->region, duckdb_string_t_data(region), region_length);
            site->region[region_length] = '\0';
            region_offset += region_item;
            site->allele_a = duckdb_string_t_data(allele_a)[0];
            site->allele_b = duckdb_string_t_data(allele_b)[0];
            site->status = BAM_SITE_MEASURED;
            site_offset++;
        }
        duckdb_destroy_data_chunk(&chunk);
    }
    if (duckdb_result_error(&result) != NULL) {
        copy_error(error, error_size, duckdb_result_error(&result));
        goto cleanup;
    }
    if (site_offset != global->site_count) {
        copy_error(error, error_size,
            "duckhts_somalier_bam_counts: panel cardinality changed during preparation");
        goto cleanup;
    }
    ok = true;

cleanup:
    if (chunk) duckdb_destroy_data_chunk(&chunk);
    duckdb_destroy_result(&result);
    if (created) {
        duckdb_result drop_result = {0};
        if (duckdb_query(connection, drop, &drop_result) != DuckDBSuccess && ok) {
            copy_error(error, error_size, duckdb_result_error(&drop_result));
            ok = false;
        }
        duckdb_destroy_result(&drop_result);
    }
    return ok;
}

static void bam_extract_global_init(duckdb_init_info info) {
    const bam_extract_bind_t *bind = duckdb_init_get_bind_data(info);
    bam_extract_global_t *global = duckdb_malloc(sizeof(*global));
    char error[BAM_EXTRACT_ERROR_BYTES] = {0};
    uint32_t jobs;

    if (!global) {
        duckdb_init_set_error(info,
            "duckhts_somalier_bam_counts: out of memory allocating global state");
        return;
    }
    memset(global, 0, sizeof(*global));
    if (bam_extract_preparation_depth != 0u ||
        pthread_mutex_trylock(&bind->registry->query_mutex) != 0) {
        duckdb_init_set_error(info,
            "duckhts_somalier_bam_counts: panel preparation is busy or recursive");
        bam_extract_global_destroy(global);
        return;
    }
    bam_extract_preparation_depth++;
    if (!load_panel_rows(bind->registry->query_connection, bind, global,
                         error, sizeof(error))) {
        bam_extract_preparation_depth--;
        pthread_mutex_unlock(&bind->registry->query_mutex);
        duckdb_init_set_error(info, error);
        bam_extract_global_destroy(global);
        return;
    }
    bam_extract_preparation_depth--;
    pthread_mutex_unlock(&bind->registry->query_mutex);

    jobs = bind->worker_count;
    if ((size_t)jobs > global->site_count) jobs = (uint32_t)global->site_count;
    for (uint32_t i = 0u; i < jobs; i++) {
        size_t width = global->site_count / jobs;
        size_t extra = global->site_count % jobs;
        global->jobs[i].begin = (size_t)i * width + (i < extra ? i : extra);
        global->jobs[i].end = global->jobs[i].begin + width +
                              (i < extra ? 1u : 0u);
    }
    global->job_count = jobs;
    duckdb_init_set_max_threads(info, jobs == 0u ? 1u : jobs);
    duckdb_init_set_init_data(info, global, bam_extract_global_destroy);
}

static void bam_extract_local_init(duckdb_init_info info) {
    bam_extract_local_t *local = duckdb_malloc(sizeof(*local));
    if (!local) {
        duckdb_init_set_error(info,
            "duckhts_somalier_bam_counts: out of memory allocating worker state");
        return;
    }
    memset(local, 0, sizeof(*local));
    duckdb_init_set_init_data(info, local, bam_extract_local_destroy);
}

static bool process_job(const bam_extract_bind_t *bind,
                        bam_extract_global_t *global,
                        const bam_extract_job_t *job,
                        char *error_buffer, size_t error_size) {
    samFile *fp = NULL;
    sam_hdr_t *header = NULL;
    hts_idx_t *index = NULL;
    faidx_t *fai = NULL;
    bam_extract_order_t *order = NULL;
    size_t available_count = 0u;
    size_t site_count = job->end - job->begin;
    const char *error = NULL;
    bool ok = false;

    fp = sam_open(bind->source_path, "r");
    if (!fp) {
        error = "duckhts_somalier_bam_counts: failed to open BAM/CRAM source";
        goto cleanup;
    }
    if (!set_remote_tuning(fp, bind->source_path,
                           bind->remote_block_bytes,
                           bind->remote_cache_bytes, &error)) goto cleanup;
    if (bind->decompression_threads != 0u &&
        hts_set_threads(fp, (int)bind->decompression_threads) < 0) {
        error = "duckhts_somalier_bam_counts: failed to configure decompression threads";
        goto cleanup;
    }
    header = sam_hdr_read(fp);
    if (!header) {
        error = "duckhts_somalier_bam_counts: failed to read alignment header";
        goto cleanup;
    }
    fai = fai_load3_format(bind->reference_path, bind->reference_index_path,
                           NULL, 0, FAI_FASTA);
    if (!fai) {
        error = "duckhts_somalier_bam_counts: failed to open reference and existing faidx";
        goto cleanup;
    }
    /* Read @SQ and require an existing FAI before attaching a CRAM reference,
     * so FASTA lengths cannot replace the alignment coordinate domain and the
     * CRAM loader cannot create an inferred reference sidecar. */
    if (hts_set_opt(fp, CRAM_OPT_REFERENCE,
                    bind->cram_reference ? bind->cram_reference :
                                           bind->reference_path) < 0) {
        error = "duckhts_somalier_bam_counts: failed to configure CRAM reference";
        goto cleanup;
    }
    index = sam_index_load3(fp, bind->source_path, bind->index_path,
                            HTS_IDX_SILENT_FAIL);
    if (!index) {
        error = "duckhts_somalier_bam_counts: indexed BAM/CRAM source is required";
        goto cleanup;
    }
    if (hisremote(bind->reference_path) && bind->reference_cache_bytes != 0u) {
        fai_set_cache_size(fai, (int)bind->reference_cache_bytes);
    }
    order = duckdb_malloc(site_count * sizeof(*order));
    if (!order) {
        error = "duckhts_somalier_bam_counts: out of memory allocating worker order";
        goto cleanup;
    }
    for (size_t offset = job->begin; offset < job->end; offset++) {
        bam_extract_site_t *site = &global->sites[offset];
        char reference_base = '\0';
        hts_pos_t alignment_length;
        memset(&site->counts, 0, sizeof(site->counts));
        site->status = BAM_SITE_MEASURED;
        site->tid = -1;
        if (!read_reference_base(fai, site, &reference_base, &site->status)) {
            error = "duckhts_somalier_bam_counts: reference lookup failed";
            goto cleanup;
        }
        if (site->status != BAM_SITE_MEASURED) continue;
        if (reference_base != site->allele_a &&
            reference_base != site->allele_b) {
            site->status = BAM_SITE_REFERENCE_MISMATCH;
            continue;
        }
        site->tid = sam_hdr_name2tid(header, site->region);
        if (site->tid < 0) {
            site->status = BAM_SITE_ALIGNMENT_CONTIG_UNAVAILABLE;
            continue;
        }
        alignment_length = sam_hdr_tid2len(header, site->tid);
        if (alignment_length < 0 ||
            (uint64_t)alignment_length < site->position1) {
            site->status = BAM_SITE_ALIGNMENT_POSITION_UNAVAILABLE;
            continue;
        }
        order[available_count].site_offset = offset;
        order[available_count].tid = site->tid;
        order[available_count].pos0 = (hts_pos_t)(site->position1 - 1u);
        available_count++;
    }
    qsort(order, available_count, sizeof(*order), compare_order);
    if (!scan_sites(fp, header, index, global->sites, order,
                    available_count, bind, &error)) goto cleanup;
    ok = true;

cleanup:
    if (!ok) copy_error(error_buffer, error_size, error);
    duckdb_free(order);
    if (fai) fai_destroy(fai);
    if (index) hts_idx_destroy(index);
    if (header) sam_hdr_destroy(header);
    if (fp) sam_close(fp);
    return ok;
}

static void write_output_row(duckdb_vector *vectors, idx_t row,
                             const bam_extract_bind_t *bind,
                             const bam_extract_global_t *global,
                             const bam_extract_site_t *site) {
    uint64_t *validity_a = duckdb_vector_get_validity(vectors[BAM_OUT_A]);
    uint64_t *validity_b = duckdb_vector_get_validity(vectors[BAM_OUT_B]);
    uint64_t *validity_other =
        duckdb_vector_get_validity(vectors[BAM_OUT_OTHER]);

    duckdb_vector_assign_string_element(vectors[BAM_OUT_SAMPLE_ID], row,
                                        bind->sample_id);
    duckdb_vector_assign_string_element(vectors[BAM_OUT_SOURCE_PATH], row,
                                        bind->source_path);
    duckdb_vector_assign_string_element_len(vectors[BAM_OUT_ASSEMBLY], row,
        global->assembly, global->assembly_length);
    ((uint64_t *)duckdb_vector_get_data(
        vectors[BAM_OUT_SITE_INDEX]))[row] = site->site_index;
    duckdb_vector_assign_string_element_len(vectors[BAM_OUT_REGION], row,
        site->region, site->region_length);
    ((uint64_t *)duckdb_vector_get_data(
        vectors[BAM_OUT_POSITION]))[row] = site->position1;
    duckdb_vector_assign_string_element_len(vectors[BAM_OUT_ALLELE_A], row,
                                            &site->allele_a, 1u);
    duckdb_vector_assign_string_element_len(vectors[BAM_OUT_ALLELE_B], row,
                                            &site->allele_b, 1u);
    duckdb_vector_assign_string_element(vectors[BAM_OUT_SOURCE_METHOD], row,
                                        "bam_cram_pileup");
    duckdb_vector_assign_string_element(vectors[BAM_OUT_COUNT_SCOPE], row,
        "observed_query_bases_after_explicit_filters");
    duckdb_vector_assign_string_element(vectors[BAM_OUT_OVERLAP_POLICY], row,
        bind->overlap_policy == DUCKHTS_BAM_SITE_OVERLAP_NONE
            ? "none" : "hileup_v0.1.0");
    duckdb_vector_assign_string_element(vectors[BAM_OUT_STATUS], row,
                                        site_status_name(site->status));
    if (site->status == BAM_SITE_MEASURED) {
        ((uint64_t *)duckdb_vector_get_data(vectors[BAM_OUT_A]))[row] =
            site->counts.allele_a;
        ((uint64_t *)duckdb_vector_get_data(vectors[BAM_OUT_B]))[row] =
            site->counts.allele_b;
        ((uint64_t *)duckdb_vector_get_data(vectors[BAM_OUT_OTHER]))[row] =
            site->counts.other;
        duckdb_validity_set_row_valid(validity_a, row);
        duckdb_validity_set_row_valid(validity_b, row);
        duckdb_validity_set_row_valid(validity_other, row);
    } else {
        duckdb_validity_set_row_invalid(validity_a, row);
        duckdb_validity_set_row_invalid(validity_b, row);
        duckdb_validity_set_row_invalid(validity_other, row);
    }
}

static void bam_extract_scan(duckdb_function_info info,
                             duckdb_data_chunk output) {
    const bam_extract_bind_t *bind = duckdb_function_get_bind_data(info);
    bam_extract_global_t *global = duckdb_function_get_init_data(info);
    bam_extract_local_t *local = duckdb_function_get_local_init_data(info);
    duckdb_vector vectors[BAM_OUT_FIELD_COUNT];
    idx_t row_count = 0u;
    idx_t vector_size = duckdb_vector_size();
    char error[BAM_EXTRACT_ERROR_BYTES] = {0};

    if (!bind || !global || !local || local->done) {
        duckdb_data_chunk_set_size(output, 0u);
        return;
    }
    for (unsigned i = 0u; i < BAM_OUT_FIELD_COUNT; i++) {
        vectors[i] = duckdb_data_chunk_get_vector(output, i);
    }
    duckdb_vector_ensure_validity_writable(vectors[BAM_OUT_A]);
    duckdb_vector_ensure_validity_writable(vectors[BAM_OUT_B]);
    duckdb_vector_ensure_validity_writable(vectors[BAM_OUT_OTHER]);

    while (row_count < vector_size && !local->done) {
        if (!local->have_job) {
            uint32_t job_index = __sync_fetch_and_add(&global->next_job, 1u);
            if (job_index >= global->job_count) {
                local->done = 1u;
                break;
            }
            if (!process_job(bind, global, &global->jobs[job_index],
                             error, sizeof(error))) {
                duckdb_function_set_error(info, error);
                local->done = 1u;
                duckdb_data_chunk_set_size(output, 0u);
                return;
            }
            local->cursor = global->jobs[job_index].begin;
            local->end = global->jobs[job_index].end;
            local->have_job = 1u;
        }
        while (row_count < vector_size && local->cursor < local->end) {
            write_output_row(vectors, row_count, bind, global,
                             &global->sites[local->cursor]);
            local->cursor++;
            row_count++;
        }
        if (local->cursor == local->end) local->have_job = 0u;
    }
    duckdb_data_chunk_set_size(output, row_count);
}

bool register_duckhts_somalier_bam_extract_functions(
    duckhts_registration_t *registration, duckdb_database database) {
    bam_extract_registry_t *registry = duckdb_malloc(sizeof(*registry));
    duckdb_table_function function = NULL;
    duckdb_logical_type varchar_type = NULL;
    duckdb_logical_type ubigint_type = NULL;
    duckdb_logical_type double_type = NULL;
    bool ok = false;

    if (!registry) {
        return duckhts_registration_error(registration,
            "DuckHTS could not allocate Somalier BAM registration state");
    }
    memset(registry, 0, sizeof(*registry));
    if (pthread_mutex_init(&registry->query_mutex, NULL) != 0) {
        duckdb_free(registry);
        return duckhts_registration_error(registration,
            "DuckHTS could not initialize the Somalier BAM registration mutex");
    }
    if (duckdb_connect(database, &registry->query_connection) != DuckDBSuccess ||
        !registry->query_connection) {
        bam_extract_registry_destroy(registry);
        return duckhts_registration_error(registration,
            "DuckHTS could not open the Somalier BAM registration connection");
    }
    function = duckdb_create_table_function();
    varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    double_type = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_table_function_set_name(function, "duckhts_somalier_bam_counts");
    for (unsigned i = 0u; i < 4u; i++) {
        duckdb_table_function_add_parameter(function, varchar_type);
    }
    duckdb_table_function_add_named_parameter(function, "panel_parquet",
                                               varchar_type);
    duckdb_table_function_add_named_parameter(function, "index_path",
                                               varchar_type);
    duckdb_table_function_add_named_parameter(function, "reference_index_path",
                                               varchar_type);
    duckdb_table_function_add_named_parameter(function, "min_mapq",
                                               ubigint_type);
    duckdb_table_function_add_named_parameter(function, "min_baseq",
                                               ubigint_type);
    duckdb_table_function_add_named_parameter(function, "require_flags",
                                               ubigint_type);
    duckdb_table_function_add_named_parameter(function, "exclude_flags",
                                               ubigint_type);
    duckdb_table_function_add_named_parameter(function, "overlap_policy",
                                               varchar_type);
    duckdb_table_function_add_named_parameter(function, "decompression_threads",
                                               ubigint_type);
    duckdb_table_function_add_named_parameter(function, "worker_count",
                                               double_type);
    duckdb_table_function_add_named_parameter(function, "max_depth",
                                               ubigint_type);
    duckdb_table_function_add_named_parameter(function, "max_overlap_qnames",
                                               ubigint_type);
    duckdb_table_function_add_named_parameter(function, "max_sites",
                                               ubigint_type);
    duckdb_table_function_add_named_parameter(function, "max_region_bytes",
                                               ubigint_type);
    duckdb_table_function_add_named_parameter(function, "remote_block_bytes",
                                               ubigint_type);
    duckdb_table_function_add_named_parameter(function, "remote_cache_bytes",
                                               ubigint_type);
    duckdb_table_function_add_named_parameter(function, "reference_cache_bytes",
                                               ubigint_type);
    duckdb_table_function_set_extra_info(function, registry,
                                         bam_extract_registry_destroy);
    duckdb_table_function_set_bind(function, bam_extract_bind);
    duckdb_table_function_set_init(function, bam_extract_global_init);
    duckdb_table_function_set_local_init(function, bam_extract_local_init);
    duckdb_table_function_set_function(function, bam_extract_scan);
    ok = duckdb_register_table_function(registration->connection, function) == DuckDBSuccess;

    duckdb_destroy_logical_type(&double_type);
    duckdb_destroy_logical_type(&ubigint_type);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_table_function(&function);
    if (!ok) {
        return duckhts_registration_error(registration,
            "DuckHTS could not register duckhts_somalier_bam_counts");
    }
    return true;
}
