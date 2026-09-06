/* Record-major typed GT/PS and its shared, original-header sample catalog. */
#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN
#include "include/bcf_genotypes.h"
#include "include/duckdb_alloc.h"
#include "include/duckdb_list.h"
#include "include/region_list.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { GENO_RECORD_INDEX, GENO_CHROM, GENO_POS, GENO_ID, GENO_REF, GENO_ALT, GENO_CALLS, GENO_COLUMNS };

typedef struct {
    char *path;
    char *index_path;
    char **regions;
    unsigned int region_count;
    duckhts_bcf_index_t index;
    duckhts_bcf_samples_t samples;
    duckhts_bcf_decode_policy_t policy;
    int decompression_threads;
    int non_reference_only;
    int catalog;
} geno_bind_t;

typedef struct {
    duckhts_bcf_scan_t scan;
    bcf1_t *record;
    duckhts_bcf_genotypes_t genotypes;
    idx_t columns[GENO_COLUMNS];
    idx_t column_count;
    uint64_t next_index;
    int ordinal_exhausted;
    int needs_calls;
    int unpack_mask;
    int done;
} geno_local_t;

static void geno_bind_destroy(void *data) {
    geno_bind_t *bind = data;
    duckdb_free(bind->path);
    duckdb_free(bind->index_path);
    free(bind->regions);
    duckhts_bcf_index_destroy(&bind->index);
    duckhts_bcf_samples_destroy(&bind->samples);
    duckdb_free(bind);
}

static void geno_local_destroy(void *data) {
    geno_local_t *local = data;
    duckhts_bcf_genotypes_destroy(&local->genotypes);
    if (local->record) bcf_destroy(local->record);
    duckhts_bcf_scan_close(&local->scan);
    duckdb_free(local);
}

/* Named NULL means the default. Allocation failure must not also mean default. */
static int geno_named_string(duckdb_bind_info info, const char *name, char **out) {
    *out = NULL;
    duckdb_value value = duckdb_bind_get_named_parameter(info, name);
    int present = value && !duckdb_is_null_value(value);
    if (present) *out = duckdb_get_varchar(value);
    if (value) duckdb_destroy_value(&value);
    return !present || *out;
}

static duckdb_logical_type geno_calls_type(void) {
    duckdb_logical_type integer = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type boolean = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type children[] = {
        duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER),
        duckdb_create_list_type(integer), duckdb_create_list_type(boolean),
        duckdb_create_logical_type(DUCKDB_TYPE_BIGINT)
    };
    const char *names[] = {"sample_index", "alleles", "phase_before", "phase_set"};
    duckdb_logical_type call = duckdb_create_struct_type(children, names, 4);
    duckdb_logical_type calls = duckdb_create_list_type(call);
    duckdb_destroy_logical_type(&call);
    for (int i = 0; i < 4; i++) duckdb_destroy_logical_type(&children[i]);
    duckdb_destroy_logical_type(&integer);
    duckdb_destroy_logical_type(&boolean);
    return calls;
}

static void geno_bind_schema(duckdb_bind_info info, int catalog) {
    duckdb_logical_type text = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    if (catalog) {
        duckdb_logical_type sample = duckdb_create_logical_type(DUCKDB_TYPE_UINTEGER);
        duckdb_bind_add_result_column(info, "sample_index", sample);
        duckdb_bind_add_result_column(info, "sample_name", text);
        duckdb_destroy_logical_type(&sample);
    } else {
        duckdb_logical_type ordinal = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
        duckdb_logical_type position = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
        duckdb_logical_type alts = duckdb_create_list_type(text);
        duckdb_logical_type calls = geno_calls_type();
        duckdb_bind_add_result_column(info, "record_index", ordinal);
        duckdb_bind_add_result_column(info, "CHROM", text);
        duckdb_bind_add_result_column(info, "POS", position);
        duckdb_bind_add_result_column(info, "ID", text);
        duckdb_bind_add_result_column(info, "REF", text);
        duckdb_bind_add_result_column(info, "ALT", alts);
        duckdb_bind_add_result_column(info, "calls", calls);
        duckdb_destroy_logical_type(&ordinal);
        duckdb_destroy_logical_type(&position);
        duckdb_destroy_logical_type(&alts);
        duckdb_destroy_logical_type(&calls);
    }
    duckdb_destroy_logical_type(&text);
}

static void geno_bind(duckdb_bind_info info, int catalog) {
    geno_bind_t *bind = duckhts_alloc_array(1, sizeof(*bind));
    duckhts_bcf_scan_t metadata = {0};
    char *selector = NULL, *region = NULL, *mode = NULL, *policy = NULL;
    char error[512] = "read_geno: out of memory allocating bind state";
    if (!bind) goto fail;
    bind->catalog = catalog;
    duckdb_value value = duckdb_bind_get_parameter(info, 0);
    int path_null = !value || duckdb_is_null_value(value);
    if (!path_null) bind->path = duckdb_get_varchar(value);
    if (value) duckdb_destroy_value(&value);
    if (path_null || (bind->path && !*bind->path)) {
        snprintf(error, sizeof(error), "%s requires a file path", catalog ? "read_bcf_samples" : "read_geno");
        goto fail;
    }
    if (!bind->path || !geno_named_string(info, "samples", &selector)) goto fail;
    int sequential = 0;
    if (!catalog) {
        if (!geno_named_string(info, "region", &region) ||
            !geno_named_string(info, "index_path", &bind->index_path) ||
            !geno_named_string(info, "scan_mode", &mode) ||
            !geno_named_string(info, "decode_error_policy", &policy)) goto fail;
        if (!duckhts_bcf_parse_scan_mode(mode, &sequential)) {
            snprintf(error, sizeof(error), "read_geno: scan_mode must be 'auto' or 'sequential'");
            goto fail;
        }
        if (!duckhts_bcf_parse_decode_policy(policy, &bind->policy)) {
            snprintf(error, sizeof(error), "read_geno: decode_error_policy must be 'null', 'warn', or 'error'");
            goto fail;
        }
        if (!duckhts_region_list_parse(region, &bind->regions, &bind->region_count, error, sizeof(error))) goto fail;
        if (sequential && bind->region_count) {
            snprintf(error, sizeof(error), "read_geno: scan_mode := 'sequential' is incompatible with region queries");
            goto fail;
        }
        value = duckdb_bind_get_named_parameter(info, "decompression_threads");
        int64_t threads = value && !duckdb_is_null_value(value) ? duckdb_get_int64(value) : 0;
        if (value) duckdb_destroy_value(&value);
        if (threads < 0 || threads > INT_MAX) {
            snprintf(error, sizeof(error), "read_geno: decompression_threads must be between 0 and INT_MAX");
            goto fail;
        }
        bind->decompression_threads = (int)threads;
        value = duckdb_bind_get_named_parameter(info, "non_reference_only");
        bind->non_reference_only = value && !duckdb_is_null_value(value) && duckdb_get_bool(value);
        if (value) duckdb_destroy_value(&value);
    }
    if (!duckhts_bcf_scan_open(&metadata, bind->path, NULL, 0, DUCKHTS_HTS_IO_PROFILE_METADATA,
            catalog ? "read_bcf_samples" : "read_geno", error, sizeof(error)) ||
        !duckhts_bcf_samples_build(&bind->samples, metadata.hdr, selector, error, sizeof(error))) goto fail;
    /* A full scan always follows the input stream, independent of index
     * availability. Only explicit regions need an immutable index snapshot. */
    if (bind->region_count) {
        int loaded = duckhts_bcf_index_load(&bind->index, hts_get_format(metadata.fp)->format,
                                            bind->path, bind->index_path, HTS_IDX_SILENT_FAIL);
        if (loaded <= 0) {
            snprintf(error, sizeof(error), "read_geno: %s", loaded < 0
                ? "out of memory finalizing index" : "region query requires an index file (.tbi or .csi)");
            goto fail;
        }
    }
    geno_bind_schema(info, catalog);
    if (catalog) duckdb_bind_set_cardinality(info, (idx_t)bind->samples.count, true);
    duckdb_bind_set_bind_data(info, bind, geno_bind_destroy);
    bind = NULL;
    goto cleanup;
fail:
    duckdb_bind_set_error(info, error);
    if (bind) geno_bind_destroy(bind);
cleanup:
    duckhts_bcf_scan_close(&metadata);
    duckdb_free(selector);
    duckdb_free(region);
    duckdb_free(mode);
    duckdb_free(policy);
}

static void geno_read_bind(duckdb_bind_info info) { geno_bind(info, 0); }
static void geno_samples_bind(duckdb_bind_info info) { geno_bind(info, 1); }

static void geno_global_init(duckdb_init_info info) {
    /* No atomic arrival-order ordinal and no shared mutable HTSlib handle. */
    duckdb_init_set_max_threads(info, 1);
}

static void geno_local_init(duckdb_init_info info) {
    const geno_bind_t *bind = duckdb_init_get_bind_data(info);
    geno_local_t *local = duckhts_alloc_array(1, sizeof(*local));
    char error[512] = "read_geno: out of memory allocating scan state";
    if (!local) goto fail;
    local->column_count = duckdb_init_get_column_count(info);
    if (local->column_count > GENO_COLUMNS) {
        snprintf(error, sizeof(error), "read_geno: invalid projection count");
        goto fail;
    }
    for (idx_t i = 0; i < local->column_count; i++) {
        idx_t col = duckdb_init_get_column_index(info, i);
        local->columns[i] = col;
        if (!bind->catalog) {
            if (col == GENO_CALLS) local->needs_calls = 1;
            if (col == GENO_ID || col == GENO_REF || col == GENO_ALT) local->unpack_mask |= BCF_UN_STR;
        }
    }
    if (!bind->catalog) {
        if (!duckhts_bcf_scan_open(&local->scan, bind->path, &bind->index, bind->decompression_threads,
                bind->region_count ? DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION : DUCKHTS_HTS_IO_PROFILE_STREAMING,
                "read_geno", error, sizeof(error)) ||
            !duckhts_bcf_samples_apply(&bind->samples, local->scan.hdr, error, sizeof(error))) goto fail;
        local->record = bcf_init();
        if (!local->record) {
            snprintf(error, sizeof(error), "read_geno: out of memory allocating BCF/VCF record");
            goto fail;
        }
        if (bind->region_count && !duckhts_bcf_scan_regions(&local->scan, bind->regions,
                bind->region_count, error, sizeof(error))) goto fail;
    }
    duckdb_init_set_init_data(info, local, geno_local_destroy);
    return;
fail:
    duckdb_init_set_error(info, error);
    if (local) geno_local_destroy(local);
}

static void geno_set_null(duckdb_vector vector, idx_t row) {
    duckdb_vector_ensure_validity_writable(vector);
    duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vector), row);
}

static int geno_write_calls(duckdb_vector vector, idx_t row, const geno_bind_t *bind,
                            const duckhts_bcf_genotypes_t *values) {
    idx_t count = 0, slots = 0;
    for (int i = 0; i < values->samples; i++) {
        const int32_t *gt = values->gt_stride ? values->gt + (size_t)i * values->gt_stride : NULL;
        int ploidy = duckhts_bcf_genotype_ploidy(gt, values->gt_stride);
        if (bind->non_reference_only && !duckhts_bcf_genotype_has_alt(gt, ploidy)) continue;
        if ((idx_t)ploidy > UINT64_MAX - slots) return 0;
        count++;
        slots += (idx_t)ploidy;
    }
    duckdb_list_entry calls, alleles, phases;
    if (!duckhts_list_extend(vector, count, &calls)) return 0;
    duckdb_vector call = duckdb_list_vector_get_child(vector);
    duckdb_vector sample = duckdb_struct_vector_get_child(call, 0);
    duckdb_vector allele_list = duckdb_struct_vector_get_child(call, 1);
    duckdb_vector phase_list = duckdb_struct_vector_get_child(call, 2);
    duckdb_vector ps = duckdb_struct_vector_get_child(call, 3);
    if (!duckhts_list_extend(allele_list, slots, &alleles) ||
        !duckhts_list_extend(phase_list, slots, &phases)) return 0;
    /* No child data access until all three reserves succeed. */
    duckdb_vector allele_child = duckdb_list_vector_get_child(allele_list);
    int32_t *allele_data = duckdb_vector_get_data(allele_child);
    bool *phase_data = duckdb_vector_get_data(duckdb_list_vector_get_child(phase_list));
    uint32_t *sample_data = duckdb_vector_get_data(sample);
    int64_t *ps_data = duckdb_vector_get_data(ps);
    duckdb_list_entry *allele_entries = duckdb_vector_get_data(allele_list);
    duckdb_list_entry *phase_entries = duckdb_vector_get_data(phase_list);
    idx_t output_call = calls.offset, output_slot = 0;
    for (int i = 0; i < values->samples; i++) {
        const int32_t *gt = values->gt_stride ? values->gt + (size_t)i * values->gt_stride : NULL;
        int ploidy = duckhts_bcf_genotype_ploidy(gt, values->gt_stride);
        if (bind->non_reference_only && !duckhts_bcf_genotype_has_alt(gt, ploidy)) continue;
        sample_data[output_call] = bind->samples.indices[i];
        allele_entries[output_call] = (duckdb_list_entry){alleles.offset + output_slot, (idx_t)ploidy};
        phase_entries[output_call] = (duckdb_list_entry){phases.offset + output_slot, (idx_t)ploidy};
        if (!values->gt_stride) {
            geno_set_null(allele_list, output_call);
            geno_set_null(phase_list, output_call);
        }
        if (!values->ps_present || values->ps[i] == bcf_int32_missing || values->ps[i] == bcf_int32_vector_end) {
            geno_set_null(ps, output_call);
        } else ps_data[output_call] = values->ps[i];
        for (int slot = 0; slot < ploidy; slot++, output_slot++) {
            if (bcf_gt_is_missing(gt[slot])) geno_set_null(allele_child, alleles.offset + output_slot);
            else allele_data[alleles.offset + output_slot] = bcf_gt_allele(gt[slot]);
            phase_data[phases.offset + output_slot] = bcf_gt_is_phased(gt[slot]);
        }
        output_call++;
    }
    ((duckdb_list_entry *)duckdb_vector_get_data(vector))[row] = calls;
    return 1;
}

static int geno_write_alts(duckdb_vector vector, idx_t row, const bcf1_t *record) {
    duckdb_list_entry entry;
    idx_t count = record->n_allele > 1 ? record->n_allele - 1 : 0;
    if (!duckhts_list_extend(vector, count, &entry)) return 0;
    duckdb_vector child = duckdb_list_vector_get_child(vector);
    for (idx_t i = 0; i < count; i++) {
        duckdb_vector_assign_string_element(child, entry.offset + i, record->d.allele[i + 1]);
    }
    ((duckdb_list_entry *)duckdb_vector_get_data(vector))[row] = entry;
    return 1;
}

static void geno_read(duckdb_function_info info, duckdb_data_chunk output) {
    const geno_bind_t *bind = duckdb_function_get_bind_data(info);
    geno_local_t *local = duckdb_function_get_local_init_data(info);
    idx_t rows = 0, capacity = duckdb_vector_size();
    char error[512];
    while (!local->done && rows < capacity) {
        if (bind->catalog) {
            if (local->next_index == (uint64_t)bind->samples.count) { local->done = 1; break; }
        } else {
            int ret = duckhts_bcf_scan_next(&local->scan, local->record);
            if (ret == -1) { local->done = 1; break; }
            if (ret < -1) {
                snprintf(error, sizeof(error), "read_geno: failed to read BCF/VCF record from %s (htslib return %d)",
                         bind->path, ret);
                goto fail;
            }
            if (local->ordinal_exhausted) {
                snprintf(error, sizeof(error), "read_geno: record_index exceeds UBIGINT capacity");
                goto fail;
            }
            if (bcf_unpack(local->record, local->unpack_mask) < 0) {
                snprintf(error, sizeof(error), "read_geno: failed to unpack BCF/VCF record");
                goto fail;
            }
            if (local->needs_calls && !duckhts_bcf_genotypes_decode(&local->genotypes,
                    local->scan.hdr, local->record, bind->policy, error, sizeof(error))) goto fail;
        }
        for (idx_t i = 0; i < local->column_count; i++) {
            duckdb_vector vector = duckdb_data_chunk_get_vector(output, i);
            idx_t col = local->columns[i];
            if (bind->catalog) {
                if (col == 0) ((uint32_t *)duckdb_vector_get_data(vector))[rows] = bind->samples.indices[local->next_index];
                else duckdb_vector_assign_string_element(vector, rows, bind->samples.names[local->next_index]);
                continue;
            }
            bcf1_t *record = local->record;
            switch (col) {
            case GENO_RECORD_INDEX:
                ((uint64_t *)duckdb_vector_get_data(vector))[rows] = local->next_index;
                break;
            case GENO_CHROM: {
                const char *name = record->rid >= 0 ? bcf_hdr_id2name(local->scan.hdr, record->rid) : NULL;
                if (name) duckdb_vector_assign_string_element(vector, rows, name);
                else geno_set_null(vector, rows);
                break;
            }
            case GENO_POS:
                if (record->pos == INT64_MAX) {
                    snprintf(error, sizeof(error), "read_geno: one-based position exceeds BIGINT");
                    goto fail;
                }
                ((int64_t *)duckdb_vector_get_data(vector))[rows] = record->pos + 1;
                break;
            case GENO_ID:
                if (!record->d.id || !*record->d.id || strcmp(record->d.id, ".") == 0) geno_set_null(vector, rows);
                else duckdb_vector_assign_string_element(vector, rows, record->d.id);
                break;
            case GENO_REF:
                if (!record->n_allele) geno_set_null(vector, rows);
                else duckdb_vector_assign_string_element(vector, rows, record->d.allele[0]);
                break;
            case GENO_ALT:
                if (!geno_write_alts(vector, rows, record)) goto list_error;
                break;
            case GENO_CALLS:
                if (!geno_write_calls(vector, rows, bind, &local->genotypes)) goto list_error;
                break;
            default:
                snprintf(error, sizeof(error), "read_geno: invalid projected column");
                goto fail;
            }
        }
        if (local->next_index == UINT64_MAX) local->ordinal_exhausted = 1;
        else local->next_index++;
        rows++;
    }
    duckdb_data_chunk_set_size(output, rows);
    return;
list_error:
    snprintf(error, sizeof(error), "read_geno: failed to grow output list at record_index %llu in %s",
             (unsigned long long)local->next_index, bind->path);
fail:
    local->done = 1;
    duckdb_function_set_error(info, error);
    duckdb_data_chunk_set_size(output, 0);
}

void register_read_geno_functions(duckdb_connection connection) {
    duckdb_logical_type text = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type boolean = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type bigint = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    for (int catalog = 0; catalog < 2; catalog++) {
        duckdb_table_function function = duckdb_create_table_function();
        duckdb_table_function_set_name(function, catalog ? "read_bcf_samples" : "read_geno");
        duckdb_table_function_add_parameter(function, text);
        duckdb_table_function_add_named_parameter(function, "samples", text);
        if (!catalog) {
            duckdb_table_function_add_named_parameter(function, "region", text);
            duckdb_table_function_add_named_parameter(function, "index_path", text);
            duckdb_table_function_add_named_parameter(function, "scan_mode", text);
            duckdb_table_function_add_named_parameter(function, "decode_error_policy", text);
            duckdb_table_function_add_named_parameter(function, "decompression_threads", bigint);
            duckdb_table_function_add_named_parameter(function, "non_reference_only", boolean);
        }
        duckdb_table_function_set_bind(function, catalog ? geno_samples_bind : geno_read_bind);
        duckdb_table_function_set_init(function, geno_global_init);
        duckdb_table_function_set_local_init(function, geno_local_init);
        duckdb_table_function_set_function(function, geno_read);
        duckdb_table_function_supports_projection_pushdown(function, true);
        duckdb_register_table_function(connection, function);
        duckdb_destroy_table_function(&function);
    }
    duckdb_destroy_logical_type(&text);
    duckdb_destroy_logical_type(&boolean);
    duckdb_destroy_logical_type(&bigint);
}
