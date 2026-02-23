/**
 * DuckHTS BAM/SAM/CRAM Reader
 *
 * Table function for reading alignment files via htslib.
 * Provides columns matching SAM spec: QNAME, FLAG, RNAME, POS, MAPQ,
 * CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL.
 *
 * Parallelism strategy (indexed BAM/CRAM, no user region):
 *   - Global init advertises max_threads = min(n_contigs, 16)
 *   - Each DuckDB thread opens its own samFile + hts_idx_t
 *   - Threads claim contigs via __sync_fetch_and_add on a shared counter
 *   - Per-thread hts_set_threads for htslib I/O decompression
 *
 * For user-supplied region queries the multi-region iterator
 *   sam_itr_regarray() is used (handles overlap dedup internally).
 *
 * API reference: htslib-1.23 samples/read_bam.c, samples/index_multireg_read.c,
 *                samples/split_thread2.c, samples/read_aux.c
 */

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <inttypes.h>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/kstring.h>

/* ================================================================
 * Helpers
 * ================================================================ */

static inline void set_null(duckdb_vector vec, idx_t row) {
    duckdb_vector_ensure_validity_writable(vec);
    uint64_t *v = duckdb_vector_get_validity(vec);
    v[row / 64] &= ~((uint64_t)1 << (row % 64));
}

/* ================================================================
 * Standard SAM AUX tags (SAMtags) for typed columns
 * ================================================================ */

typedef struct {
    const char *tag;
    char type;     /* A, i, f, Z, H, B */
    char subtype;  /* for B arrays: c,C,s,S,i,I,f */
} bam_std_tag_t;

static const bam_std_tag_t BAM_STD_TAGS[] = {
    {"AM", 'i', 0}, {"AS", 'i', 0}, {"BC", 'Z', 0}, {"BQ", 'Z', 0},
    {"BZ", 'Z', 0}, {"CB", 'Z', 0}, {"CC", 'Z', 0}, {"CG", 'B', 'I'},
    {"CM", 'i', 0}, {"CO", 'Z', 0}, {"CP", 'i', 0}, {"CQ", 'Z', 0},
    {"CR", 'Z', 0}, {"CS", 'Z', 0}, {"CT", 'Z', 0}, {"CY", 'Z', 0},
    {"E2", 'Z', 0}, {"FI", 'i', 0}, {"FS", 'Z', 0}, {"FZ", 'B', 'S'},
    {"H0", 'i', 0}, {"H1", 'i', 0}, {"H2", 'i', 0}, {"HI", 'i', 0},
    {"IH", 'i', 0}, {"LB", 'Z', 0}, {"MC", 'Z', 0}, {"MD", 'Z', 0},
    {"MI", 'Z', 0}, {"ML", 'B', 'C'}, {"MM", 'Z', 0}, {"MN", 'i', 0},
    {"MQ", 'i', 0}, {"NH", 'i', 0}, {"NM", 'i', 0}, {"OA", 'Z', 0},
    {"OC", 'Z', 0}, {"OP", 'i', 0}, {"OQ", 'Z', 0}, {"OX", 'Z', 0},
    {"PG", 'Z', 0}, {"PQ", 'i', 0}, {"PT", 'Z', 0}, {"PU", 'Z', 0},
    {"Q2", 'Z', 0}, {"QT", 'Z', 0}, {"QX", 'Z', 0}, {"R2", 'Z', 0},
    {"RG", 'Z', 0}, {"RX", 'Z', 0}, {"SA", 'Z', 0}, {"SM", 'i', 0},
    {"TC", 'i', 0}, {"TS", 'A', 0}, {"U2", 'Z', 0}, {"UQ", 'i', 0},
    {NULL, 0, 0}
};

static int bam_std_tag_index(const char *tag) {
    for (int i = 0; BAM_STD_TAGS[i].tag; i++) {
        if (tag[0] == BAM_STD_TAGS[i].tag[0] &&
            tag[1] == BAM_STD_TAGS[i].tag[1] &&
            tag[2] == '\0') {
            return i;
        }
    }
    return -1;
}

static duckdb_logical_type bam_std_tag_type(const bam_std_tag_t *t) {
    switch (t->type) {
        case 'A':
        case 'Z':
        case 'H':
            return duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
        case 'i':
            return duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
        case 'f':
            return duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
        case 'B': {
            duckdb_logical_type child =
                (t->subtype == 'f') ? duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE)
                                    : duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
            duckdb_logical_type list = duckdb_create_list_type(child);
            duckdb_destroy_logical_type(&child);
            return list;
        }
        default:
            return duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    }
}

static void bam_assign_list_int(duckdb_vector vec, idx_t row, const uint8_t *aux) {
    uint32_t len = bam_auxB_len(aux);
    duckdb_list_entry entry;
    entry.offset = duckdb_list_vector_get_size(vec);
    entry.length = len;
    duckdb_list_vector_reserve(vec, entry.offset + entry.length);
    duckdb_list_vector_set_size(vec, entry.offset + entry.length);

    duckdb_vector child = duckdb_list_vector_get_child(vec);
    int64_t *data = (int64_t *)duckdb_vector_get_data(child);
    for (uint32_t i = 0; i < len; i++) {
        data[entry.offset + i] = bam_auxB2i(aux, i);
    }
    duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(vec);
    list_data[row] = entry;
}

static void bam_assign_list_double(duckdb_vector vec, idx_t row, const uint8_t *aux) {
    uint32_t len = bam_auxB_len(aux);
    duckdb_list_entry entry;
    entry.offset = duckdb_list_vector_get_size(vec);
    entry.length = len;
    duckdb_list_vector_reserve(vec, entry.offset + entry.length);
    duckdb_list_vector_set_size(vec, entry.offset + entry.length);

    duckdb_vector child = duckdb_list_vector_get_child(vec);
    double *data = (double *)duckdb_vector_get_data(child);
    for (uint32_t i = 0; i < len; i++) {
        data[entry.offset + i] = bam_auxB2f(aux, i);
    }
    duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(vec);
    list_data[row] = entry;
}

static void bam_aux_to_string(const uint8_t *aux, kstring_t *ks) {
    ks->l = 0;
    char type = bam_aux_type(aux);
    switch (type) {
        case 'A': {
            kputc(bam_aux2A(aux), ks);
            break;
        }
        case 'i':
        case 'I':
        case 's':
        case 'S':
        case 'c':
        case 'C':
            ksprintf(ks, "%" PRId64, bam_aux2i(aux));
            break;
        case 'f':
        case 'd':
            ksprintf(ks, "%g", bam_aux2f(aux));
            break;
        case 'Z':
        case 'H': {
            const char *z = bam_aux2Z(aux);
            if (z) kputs(z, ks);
            break;
        }
        case 'B': {
            char subtype = aux[1];
            uint32_t len = bam_auxB_len(aux);
            kputc(subtype, ks);
            for (uint32_t i = 0; i < len; i++) {
                kputc(',', ks);
                if (subtype == 'f' || subtype == 'd') {
                    ksprintf(ks, "%g", bam_auxB2f(aux, i));
                } else {
                    ksprintf(ks, "%" PRId64, bam_auxB2i(aux, i));
                }
            }
            break;
        }
        default:
            break;
    }
}

/* ================================================================
 * Column indices — matches SAM spec field order
 * ================================================================ */

enum {
    BAM_COL_QNAME = 0,
    BAM_COL_FLAG,
    BAM_COL_RNAME,
    BAM_COL_POS,
    BAM_COL_MAPQ,
    BAM_COL_CIGAR,
    BAM_COL_RNEXT,
    BAM_COL_PNEXT,
    BAM_COL_TLEN,
    BAM_COL_SEQ,
    BAM_COL_QUAL,
    BAM_COL_READ_GROUP_ID,
    BAM_COL_SAMPLE_ID,
    BAM_COL_CORE_COUNT
};

/* ================================================================
 * Bind Data — shared across all threads (immutable after bind)
 * ================================================================ */

typedef struct {
    char *file_path;
    char *index_path;
    char *reference;

    /* Parsed from the "region" named parameter.
     * May contain comma-separated multi-region specs. */
    char *region;       /* original string (owned) */
    char **regions;     /* split array for sam_itr_regarray (owned) */
    unsigned int n_regions;

    int has_index;
    int n_contigs;      /* sam_hdr_nref(); used for parallel partitioning */
    int standard_tags;
    int auxiliary_tags;
    int std_col_start;
    int std_col_count;
    int aux_col_idx;
} bam_bind_data_t;

/* ================================================================
 * Global Init Data — shared mutable state for parallel contig claim
 * ================================================================ */

typedef struct {
    int current_contig; /* atomic counter — threads fetch-and-add */
    int n_contigs;
    int has_region;
} bam_global_init_data_t;

/* ================================================================
 * Local Init Data — per-thread scanning state
 * ================================================================ */

typedef struct {
    samFile *fp;
    sam_hdr_t *hdr;
    bam1_t *rec;
    hts_idx_t *idx;
    hts_itr_t *itr;

    int done;
    int is_parallel;
    int assigned_contig;
    int needs_next_contig;

    idx_t column_count;
    idx_t *column_ids;

    /* Reusable buffers for SEQ / QUAL / CIGAR conversion */
    char *seq_buf;
    size_t seq_buf_cap;
    char *qual_buf;
    size_t qual_buf_cap;
    char cigar_buf[8192];

    /* Read group caching */
    char *last_rg_id;
    char *last_sample_id;
    kstring_t rg_tmp;
} bam_local_init_data_t;

/* ================================================================
 * Destructors
 * ================================================================ */

static void destroy_bam_bind(void *data) {
    bam_bind_data_t *b = (bam_bind_data_t *)data;
    if (!b) return;
    if (b->file_path) duckdb_free(b->file_path);
    if (b->index_path) duckdb_free(b->index_path);
    if (b->reference) duckdb_free(b->reference);
    if (b->region) duckdb_free(b->region);
    if (b->regions) {
        for (unsigned int i = 0; i < b->n_regions; i++)
            duckdb_free(b->regions[i]);
        duckdb_free(b->regions);
    }
    duckdb_free(b);
}

static void destroy_bam_global(void *data) {
    if (data) duckdb_free(data);
}

static void destroy_bam_local(void *data) {
    bam_local_init_data_t *l = (bam_local_init_data_t *)data;
    if (!l) return;
    if (l->itr) hts_itr_destroy(l->itr);
    if (l->idx) hts_idx_destroy(l->idx);
    if (l->rec) bam_destroy1(l->rec);
    if (l->hdr) sam_hdr_destroy(l->hdr);
    if (l->fp) sam_close(l->fp);
    if (l->column_ids) duckdb_free(l->column_ids);
    if (l->seq_buf) free(l->seq_buf);
    if (l->qual_buf) free(l->qual_buf);
    if (l->last_rg_id) free(l->last_rg_id);
    if (l->last_sample_id) free(l->last_sample_id);
    free(l->rg_tmp.s);
    duckdb_free(l);
}

/* ================================================================
 * Region parsing helper — split comma-separated regions
 * ================================================================ */

static void parse_regions(const char *region_str, char ***out_regions, unsigned int *out_count) {
    *out_regions = NULL;
    *out_count = 0;
    if (!region_str || strlen(region_str) == 0) return;

    /* Count commas to determine array size */
    unsigned int count = 1;
    for (const char *p = region_str; *p; p++)
        if (*p == ',') count++;

    char **arr = (char **)duckdb_malloc(sizeof(char *) * count);
    char *dup = (char *)duckdb_malloc(strlen(region_str) + 1);
    strcpy(dup, region_str);

    unsigned int idx = 0;
    char *tok = strtok(dup, ",");
    while (tok && idx < count) {
        size_t len = strlen(tok) + 1;
        arr[idx] = (char *)duckdb_malloc(len);
        memcpy(arr[idx], tok, len);
        idx++;
        tok = strtok(NULL, ",");
    }
    duckdb_free(dup);

    *out_regions = arr;
    *out_count = idx;
}

/* ================================================================
 * Ensure reusable buffers are large enough
 * ================================================================ */

static int ensure_seq_buf(bam_local_init_data_t *l, int seq_len) {
    size_t need = (size_t)(seq_len + 1);
    if (need > l->seq_buf_cap) {
        size_t new_cap = need * 2;
        char *new_seq = (char *)realloc(l->seq_buf, new_cap);
        if (!new_seq) return 0;
        char *new_qual = (char *)realloc(l->qual_buf, new_cap);
        if (!new_qual) {
            l->seq_buf = new_seq;
            return 0;
        }
        l->seq_buf = new_seq;
        l->qual_buf = new_qual;
        l->seq_buf_cap = new_cap;
    }
    return 1;
}

/* ================================================================
 * CIGAR → string
 * Uses bam_cigar_op / bam_cigar_oplen / bam_cigar_opchr from sam.h
 * ================================================================ */

static void cigar_to_string(const uint32_t *cigar, int n_cigar,
                            char *buf, size_t buf_size) {
    size_t pos = 0;
    for (int i = 0; i < n_cigar && pos + 12 < buf_size; i++) {
        pos += (size_t)snprintf(buf + pos, buf_size - pos, "%d%c",
                                bam_cigar_oplen(cigar[i]),
                                bam_cigar_opchr(bam_cigar_op(cigar[i])));
    }
    buf[pos] = '\0';
}

/* ================================================================
 * SEQ → string
 * Uses seq_nt16_str[] and bam_seqi() from htslib (sam.h)
 * ================================================================ */

static void seq_to_string(const uint8_t *seq_data, int len, char *buf) {
    for (int i = 0; i < len; i++)
        buf[i] = seq_nt16_str[bam_seqi(seq_data, i)];
    buf[len] = '\0';
}

/* ================================================================
 * QUAL → string (Phred+33, matching SAM text output)
 * ================================================================ */

static void qual_to_string(const uint8_t *qual, int len, char *buf) {
    for (int i = 0; i < len; i++)
        buf[i] = (char)(qual[i] + 33);
    buf[len] = '\0';
}

/* ================================================================
 * Bind
 * ================================================================ */

static void bam_read_bind(duckdb_bind_info info) {
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    char *file_path = duckdb_get_varchar(path_val);
    duckdb_destroy_value(&path_val);

    if (!file_path || strlen(file_path) == 0) {
        duckdb_bind_set_error(info, "read_bam requires a file path");
        if (file_path) duckdb_free(file_path);
        return;
    }

    /* Parse optional region parameter */
    char *region = NULL;
    duckdb_value region_val = duckdb_bind_get_named_parameter(info, "region");
    if (region_val && !duckdb_is_null_value(region_val))
        region = duckdb_get_varchar(region_val);
    if (region_val) duckdb_destroy_value(&region_val);

    /* Parse optional explicit index path */
    char *index_path = NULL;
    duckdb_value index_val = duckdb_bind_get_named_parameter(info, "index_path");
    if (index_val && !duckdb_is_null_value(index_val))
        index_path = duckdb_get_varchar(index_val);
    if (index_val) duckdb_destroy_value(&index_val);

    /* Parse optional reference for CRAM */
    char *reference = NULL;
    duckdb_value ref_val = duckdb_bind_get_named_parameter(info, "reference");
    if (ref_val && !duckdb_is_null_value(ref_val))
        reference = duckdb_get_varchar(ref_val);
    if (ref_val) duckdb_destroy_value(&ref_val);

    /* Probe the file: open, read header, check index */
    samFile *fp = sam_open(file_path, "r");
    if (!fp) {
        char err[512];
        snprintf(err, sizeof(err), "Failed to open SAM/BAM/CRAM file: %s", file_path);
        duckdb_bind_set_error(info, err);
        duckdb_free(file_path);
        if (index_path) duckdb_free(index_path);
        if (region) duckdb_free(region);
        if (reference) duckdb_free(reference);
        return;
    }

    if (reference) {
        hts_set_opt(fp, CRAM_OPT_REFERENCE, reference);
    }
    sam_hdr_t *hdr = sam_hdr_read(fp);
    if (!hdr) {
        sam_close(fp);
        duckdb_bind_set_error(info, "Failed to read SAM/BAM/CRAM header");
        duckdb_free(file_path);
        if (index_path) duckdb_free(index_path);
        if (region) duckdb_free(region);
        if (reference) duckdb_free(reference);
        return;
    }

    bam_bind_data_t *bind = (bam_bind_data_t *)duckdb_malloc(sizeof(bam_bind_data_t));
    memset(bind, 0, sizeof(bam_bind_data_t));
    bind->file_path = file_path;
    bind->index_path = index_path;
    bind->reference = reference;
    bind->region = region;
    bind->n_contigs = sam_hdr_nref(hdr);
    bind->standard_tags = 0;
    bind->auxiliary_tags = 0;
    bind->std_col_start = BAM_COL_CORE_COUNT;
    bind->std_col_count = 0;
    bind->aux_col_idx = -1;

    /* Parse comma-separated regions (if any) */
    parse_regions(region, &bind->regions, &bind->n_regions);

    /* Optional tag controls */
    duckdb_value std_val = duckdb_bind_get_named_parameter(info, "standard_tags");
    if (std_val && !duckdb_is_null_value(std_val)) {
        bind->standard_tags = duckdb_get_bool(std_val) ? 1 : 0;
    }
    if (std_val) duckdb_destroy_value(&std_val);

    duckdb_value aux_val = duckdb_bind_get_named_parameter(info, "auxiliary_tags");
    if (aux_val && !duckdb_is_null_value(aux_val)) {
        bind->auxiliary_tags = duckdb_get_bool(aux_val) ? 1 : 0;
    }
    if (aux_val) duckdb_destroy_value(&aux_val);

    /* Check for index availability */
    hts_idx_t *idx = sam_index_load3(fp, file_path, index_path, HTS_IDX_SILENT_FAIL);
    if (idx) {
        bind->has_index = 1;
        hts_idx_destroy(idx);
    }

    sam_hdr_destroy(hdr);
    sam_close(fp);

    /* ----- Define output schema (SAM spec columns) ----- */
    duckdb_logical_type varchar_type   = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type int32_type     = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type bigint_type    = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type usmallint_type = duckdb_create_logical_type(DUCKDB_TYPE_USMALLINT);

    duckdb_bind_add_result_column(info, "QNAME", varchar_type);
    duckdb_bind_add_result_column(info, "FLAG",  usmallint_type);
    duckdb_bind_add_result_column(info, "RNAME", varchar_type);
    duckdb_bind_add_result_column(info, "POS",   bigint_type);
    duckdb_bind_add_result_column(info, "MAPQ",  int32_type);
    duckdb_bind_add_result_column(info, "CIGAR", varchar_type);
    duckdb_bind_add_result_column(info, "RNEXT", varchar_type);
    duckdb_bind_add_result_column(info, "PNEXT", bigint_type);
    duckdb_bind_add_result_column(info, "TLEN",  bigint_type);
    duckdb_bind_add_result_column(info, "SEQ",   varchar_type);
    duckdb_bind_add_result_column(info, "QUAL",  varchar_type);
    duckdb_bind_add_result_column(info, "READ_GROUP_ID", varchar_type);
    duckdb_bind_add_result_column(info, "SAMPLE_ID", varchar_type);

    if (bind->standard_tags) {
        bind->std_col_start = BAM_COL_CORE_COUNT;
        int count = 0;
        for (int i = 0; BAM_STD_TAGS[i].tag; i++) {
            duckdb_logical_type t = bam_std_tag_type(&BAM_STD_TAGS[i]);
            duckdb_bind_add_result_column(info, BAM_STD_TAGS[i].tag, t);
            duckdb_destroy_logical_type(&t);
            count++;
        }
        bind->std_col_count = count;
    }

    if (bind->auxiliary_tags) {
        duckdb_logical_type key_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
        duckdb_logical_type val_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
        duckdb_logical_type map_type = duckdb_create_map_type(key_type, val_type);
        bind->aux_col_idx = BAM_COL_CORE_COUNT + bind->std_col_count;
        duckdb_bind_add_result_column(info, "AUXILIARY_TAGS", map_type);
        duckdb_destroy_logical_type(&key_type);
        duckdb_destroy_logical_type(&val_type);
        duckdb_destroy_logical_type(&map_type);
    }

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&int32_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&usmallint_type);

    duckdb_bind_set_bind_data(info, bind, destroy_bam_bind);
}

/* ================================================================
 * Global Init — set up parallel contig partitioning
 * ================================================================ */

static void bam_read_global_init(duckdb_init_info info) {
    bam_bind_data_t *bind = (bam_bind_data_t *)duckdb_init_get_bind_data(info);

    bam_global_init_data_t *global = (bam_global_init_data_t *)duckdb_malloc(
        sizeof(bam_global_init_data_t));
    memset(global, 0, sizeof(bam_global_init_data_t));

    global->current_contig = 0;
    global->has_region = (bind->n_regions > 0);

    /*
     * Parallel scan: each thread claims contigs via sam_itr_queryi().
     * Only when indexed and no user-specified region.
     */
    if (bind->has_index && bind->n_contigs > 1 && !global->has_region) {
        global->n_contigs = bind->n_contigs;
        idx_t max_threads = (idx_t)bind->n_contigs;
        if (max_threads > 16) max_threads = 16;
        duckdb_init_set_max_threads(info, max_threads);
    } else {
        global->n_contigs = 0;
        duckdb_init_set_max_threads(info, 1);
    }

    duckdb_init_set_init_data(info, global, destroy_bam_global);
}

/* ================================================================
 * Local Init — per-thread: own file handle, index, iterator
 * ================================================================ */

static void bam_read_local_init(duckdb_init_info info) {
    bam_bind_data_t *bind = (bam_bind_data_t *)duckdb_init_get_bind_data(info);

    bam_local_init_data_t *local = (bam_local_init_data_t *)duckdb_malloc(
        sizeof(bam_local_init_data_t));
    memset(local, 0, sizeof(bam_local_init_data_t));

    int is_parallel = (bind->has_index && bind->n_contigs > 1 && bind->n_regions == 0);
    local->is_parallel = is_parallel;
    local->assigned_contig = -1;
    local->needs_next_contig = is_parallel;

    /* Each thread opens its own file handle (required for parallel seeks) */
    local->fp = sam_open(bind->file_path, "r");
    if (!local->fp) {
        duckdb_init_set_error(info, "Failed to open SAM/BAM/CRAM file");
        duckdb_free(local);
        return;
    }

    if (bind->reference) {
        if (hts_set_opt(local->fp, CRAM_OPT_REFERENCE, bind->reference) < 0) {
            duckdb_init_set_error(info, "Failed to set CRAM reference");
            sam_close(local->fp);
            duckdb_free(local);
            return;
        }
    }

    /* Enable htslib I/O threads for BAM/CRAM decompression.
     * hts_set_threads creates non-shared threads for this file handle. */
    hts_set_threads(local->fp, 2);

    /* Read header — each thread needs its own copy */
    local->hdr = sam_hdr_read(local->fp);
    if (!local->hdr) {
        sam_close(local->fp); local->fp = NULL;
        duckdb_init_set_error(info, "Failed to read SAM/BAM/CRAM header");
        duckdb_free(local);
        return;
    }

    /* Allocate a reusable bam1_t record */
    local->rec = bam_init1();
    local->rg_tmp.l = 0;
    local->rg_tmp.m = 0;
    local->rg_tmp.s = NULL;

    /* Load index if needed for parallel scanning or region queries */
    if (is_parallel || bind->n_regions > 0) {
        local->idx = sam_index_load3(local->fp, bind->file_path, bind->index_path, HTS_IDX_SILENT_FAIL);
        if (!local->idx) {
            if (bind->n_regions > 0) {
                duckdb_init_set_error(info,
                    "Region query requires an index (.bai/.csi/.crai)");
                destroy_bam_local(local);
                return;
            }
            /* No index but parallel was requested — fall back to sequential */
            local->is_parallel = 0;
            local->needs_next_contig = 0;
        }
    }

    /* User-supplied region(s): use sam_itr_regarray for multi-region support.
     * htslib handles overlap deduplication internally. */
    if (!is_parallel && bind->n_regions > 0 && local->idx) {
        local->itr = sam_itr_regarray(local->idx, local->hdr,
                                       bind->regions, bind->n_regions);
        if (!local->itr) {
            char err[512];
            snprintf(err, sizeof(err), "No reads found for region(s): %s",
                     bind->region);
            duckdb_init_set_error(info, err);
            destroy_bam_local(local);
            return;
        }
    }

    local->done = 0;

    /* Projection pushdown: record which columns DuckDB actually needs */
    local->column_count = duckdb_init_get_column_count(info);
    local->column_ids = (idx_t *)duckdb_malloc(sizeof(idx_t) * local->column_count);
    for (idx_t i = 0; i < local->column_count; i++)
        local->column_ids[i] = duckdb_init_get_column_index(info, i);

    duckdb_init_set_init_data(info, local, destroy_bam_local);
}

/* ================================================================
 * Claim next contig for parallel scanning
 * Returns 1 if a contig was claimed, 0 if done
 * ================================================================ */

static int claim_next_contig(bam_local_init_data_t *local,
                             bam_global_init_data_t *global) {
    if (!local->is_parallel || !global || global->n_contigs == 0)
        return 0;

    int next = __sync_fetch_and_add(&global->current_contig, 1);
    if (next >= global->n_contigs)
        return 0;  /* all contigs claimed */

    /* Destroy previous iterator */
    if (local->itr) {
        hts_itr_destroy(local->itr);
        local->itr = NULL;
    }

    local->assigned_contig = next;

    /* sam_itr_queryi: iterate reads overlapping the entire contig.
     * tid=next, beg=0, end=HTS_POS_MAX covers the full reference. */
    local->itr = sam_itr_queryi(local->idx, next, 0, HTS_POS_MAX);
    if (!local->itr) {
        /* Contig may have zero aligned reads — skip to next */
        return claim_next_contig(local, global);
    }

    local->needs_next_contig = 0;
    return 1;
}

/* ================================================================
 * Scan — fill one DuckDB data chunk
 * ================================================================ */

static void bam_read_function(duckdb_function_info info, duckdb_data_chunk output) {
    bam_bind_data_t *bind =
        (bam_bind_data_t *)duckdb_function_get_bind_data(info);
    bam_global_init_data_t *global =
        (bam_global_init_data_t *)duckdb_function_get_init_data(info);
    bam_local_init_data_t *local =
        (bam_local_init_data_t *)duckdb_function_get_local_init_data(info);

    if (!local || local->done) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    /* For parallel scans, claim first/next contig if needed */
    if (local->needs_next_contig) {
        if (!claim_next_contig(local, global)) {
            local->done = 1;
            duckdb_data_chunk_set_size(output, 0);
            return;
        }
    }

    idx_t vector_size = duckdb_vector_size();
    idx_t row_count = 0;

    while (row_count < vector_size) {
        int ret;
        if (local->itr)
            ret = sam_itr_next(local->fp, local->itr, local->rec);
        else
            ret = sam_read1(local->fp, local->hdr, local->rec);

        if (ret < 0) {
            /* ret == -1: EOF/end-of-region.  ret < -1: error. */
            if (local->is_parallel && ret == -1) {
                /* Try next contig */
                local->needs_next_contig = 1;
                if (!claim_next_contig(local, global)) {
                    local->done = 1;
                }
                break;
            }
            local->done = 1;
            break;
        }

        bam1_t *b = local->rec;
        int seq_len = b->core.l_qseq;

        /* Grow SEQ/QUAL conversion buffers if needed */
        if (!ensure_seq_buf(local, seq_len)) {
            duckdb_function_set_error(info, "read_bam: out of memory allocating sequence buffers");
            local->done = 1;
            duckdb_data_chunk_set_size(output, 0);
            return;
        }

        for (idx_t i = 0; i < local->column_count; i++) {
            idx_t col_id = local->column_ids[i];
            duckdb_vector vec = duckdb_data_chunk_get_vector(output, i);

            switch (col_id) {

            case BAM_COL_QNAME: {
                const char *qname = bam_get_qname(b);
                duckdb_vector_assign_string_element(vec, row_count,
                                                     qname ? qname : "*");
                break;
            }

            case BAM_COL_FLAG: {
                uint16_t *data = (uint16_t *)duckdb_vector_get_data(vec);
                data[row_count] = b->core.flag;
                break;
            }

            case BAM_COL_RNAME: {
                const char *rname = (b->core.tid >= 0)
                    ? sam_hdr_tid2name(local->hdr, b->core.tid)
                    : NULL;
                duckdb_vector_assign_string_element(vec, row_count,
                                                     rname ? rname : "*");
                break;
            }

            case BAM_COL_POS: {
                int64_t *data = (int64_t *)duckdb_vector_get_data(vec);
                data[row_count] = b->core.pos + 1;  /* 1-based */
                break;
            }

            case BAM_COL_MAPQ: {
                int32_t *data = (int32_t *)duckdb_vector_get_data(vec);
                data[row_count] = (int32_t)b->core.qual;
                break;
            }

            case BAM_COL_CIGAR: {
                if (b->core.n_cigar > 0) {
                    cigar_to_string(bam_get_cigar(b), (int)b->core.n_cigar,
                                    local->cigar_buf, sizeof(local->cigar_buf));
                    duckdb_vector_assign_string_element(vec, row_count,
                                                         local->cigar_buf);
                } else {
                    duckdb_vector_assign_string_element(vec, row_count, "*");
                }
                break;
            }

            case BAM_COL_RNEXT: {
                const char *rnext = (b->core.mtid >= 0)
                    ? sam_hdr_tid2name(local->hdr, b->core.mtid)
                    : NULL;
                duckdb_vector_assign_string_element(vec, row_count,
                                                     rnext ? rnext : "*");
                break;
            }

            case BAM_COL_PNEXT: {
                int64_t *data = (int64_t *)duckdb_vector_get_data(vec);
                data[row_count] = b->core.mpos + 1;
                break;
            }

            case BAM_COL_TLEN: {
                int64_t *data = (int64_t *)duckdb_vector_get_data(vec);
                data[row_count] = b->core.isize;
                break;
            }

            case BAM_COL_SEQ: {
                if (seq_len > 0) {
                    seq_to_string(bam_get_seq(b), seq_len, local->seq_buf);
                    duckdb_vector_assign_string_element(vec, row_count,
                                                         local->seq_buf);
                } else {
                    duckdb_vector_assign_string_element(vec, row_count, "*");
                }
                break;
            }

            case BAM_COL_QUAL: {
                if (seq_len > 0 && bam_get_qual(b)[0] != 255) {
                    qual_to_string(bam_get_qual(b), seq_len, local->qual_buf);
                    duckdb_vector_assign_string_element(vec, row_count,
                                                         local->qual_buf);
                } else {
                    duckdb_vector_assign_string_element(vec, row_count, "*");
                }
                break;
            }

            case BAM_COL_READ_GROUP_ID: {
                uint8_t *aux = bam_aux_get(b, "RG");
                if (aux) {
                    const char *rg = bam_aux2Z(aux);
                    if (rg) {
                        duckdb_vector_assign_string_element(vec, row_count, rg);
                    } else {
                        set_null(vec, row_count);
                    }
                } else {
                    set_null(vec, row_count);
                }
                break;
            }

            case BAM_COL_SAMPLE_ID: {
                uint8_t *aux = bam_aux_get(b, "RG");
                const char *rg = aux ? bam_aux2Z(aux) : NULL;
                if (!rg) {
                    set_null(vec, row_count);
                    break;
                }
                if (!local->last_rg_id || strcmp(local->last_rg_id, rg) != 0) {
                    free(local->last_rg_id);
                    free(local->last_sample_id);
                    local->last_rg_id = strdup(rg);
                    local->last_sample_id = NULL;
                    local->rg_tmp.l = 0;
                    if (sam_hdr_find_tag_id(local->hdr, "RG", "ID", rg, "SM", &local->rg_tmp) == 0 &&
                        local->rg_tmp.s && local->rg_tmp.l > 0) {
                        local->last_sample_id = strdup(local->rg_tmp.s);
                    }
                }
                if (local->last_sample_id) {
                    duckdb_vector_assign_string_element(vec, row_count, local->last_sample_id);
                } else {
                    set_null(vec, row_count);
                }
                break;
            }

            default: {
                if (col_id >= BAM_COL_CORE_COUNT) {
                    if (bind->standard_tags &&
                        col_id < (idx_t)(BAM_COL_CORE_COUNT + bind->std_col_count)) {
                        int std_idx = (int)(col_id - BAM_COL_CORE_COUNT);
                        const bam_std_tag_t *st = &BAM_STD_TAGS[std_idx];
                        uint8_t *aux = bam_aux_get(b, st->tag);
                        if (!aux) {
                            set_null(vec, row_count);
                            break;
                        }
                        switch (st->type) {
                            case 'A': {
                                char tmp[2] = { bam_aux2A(aux), '\0' };
                                duckdb_vector_assign_string_element(vec, row_count, tmp);
                                break;
                            }
                            case 'Z':
                            case 'H': {
                                const char *z = bam_aux2Z(aux);
                                if (z) duckdb_vector_assign_string_element(vec, row_count, z);
                                else set_null(vec, row_count);
                                break;
                            }
                            case 'i': {
                                int64_t *data = (int64_t *)duckdb_vector_get_data(vec);
                                data[row_count] = bam_aux2i(aux);
                                break;
                            }
                            case 'f': {
                                double *data = (double *)duckdb_vector_get_data(vec);
                                data[row_count] = bam_aux2f(aux);
                                break;
                            }
                            case 'B': {
                                char subtype = aux[1];
                                if (subtype == 'f' || subtype == 'd') {
                                    bam_assign_list_double(vec, row_count, aux);
                                } else {
                                    bam_assign_list_int(vec, row_count, aux);
                                }
                                break;
                            }
                            default:
                                set_null(vec, row_count);
                                break;
                        }
                    } else if (bind->auxiliary_tags &&
                               (int)col_id == bind->aux_col_idx) {
                        uint8_t *aux = bam_aux_first(b);
                        int count = 0;
                        while (aux) {
                            char tagbuf[3] = {0, 0, 0};
                            const char *tag = bam_aux_tag(aux);
                            tagbuf[0] = tag[0];
                            tagbuf[1] = tag[1];
                            if (!(bind->standard_tags && bam_std_tag_index(tagbuf) >= 0)) {
                                count++;
                            }
                            aux = bam_aux_next(b, aux);
                        }

                        if (count == 0) {
                            duckdb_vector_ensure_validity_writable(vec);
                            uint64_t *validity = duckdb_vector_get_validity(vec);
                            duckdb_validity_set_row_invalid(validity, row_count);
                            duckdb_list_entry entry = {duckdb_list_vector_get_size(vec), 0};
                            duckdb_list_entry *list_data =
                                (duckdb_list_entry *)duckdb_vector_get_data(vec);
                            list_data[row_count] = entry;
                            break;
                        }

                        duckdb_list_entry entry;
                        entry.offset = duckdb_list_vector_get_size(vec);
                        entry.length = count;
                        duckdb_list_vector_reserve(vec, entry.offset + entry.length);
                        duckdb_list_vector_set_size(vec, entry.offset + entry.length);

                        duckdb_vector child = duckdb_list_vector_get_child(vec);
                        duckdb_vector key_vec = duckdb_struct_vector_get_child(child, 0);
                        duckdb_vector val_vec = duckdb_struct_vector_get_child(child, 1);

                        kstring_t ks = {0, 0, NULL};
                        aux = bam_aux_first(b);
                        int write_idx = 0;
                        while (aux && write_idx < count) {
                            char tagbuf[3] = {0, 0, 0};
                            const char *tag = bam_aux_tag(aux);
                            tagbuf[0] = tag[0];
                            tagbuf[1] = tag[1];
                            if (bind->standard_tags && bam_std_tag_index(tagbuf) >= 0) {
                                aux = bam_aux_next(b, aux);
                                continue;
                            }
                            duckdb_vector_assign_string_element(
                                key_vec, entry.offset + write_idx, tagbuf);
                            bam_aux_to_string(aux, &ks);
                            duckdb_vector_assign_string_element(
                                val_vec, entry.offset + write_idx, ks.s ? ks.s : "");
                            write_idx++;
                            aux = bam_aux_next(b, aux);
                        }
                        free(ks.s);

                        duckdb_list_entry *list_data =
                            (duckdb_list_entry *)duckdb_vector_get_data(vec);
                        list_data[row_count] = entry;
                    }
                }
                break;
            }
            } /* switch */
        } /* for columns */

        row_count++;
    } /* while rows */

    duckdb_data_chunk_set_size(output, row_count);
}

/* ================================================================
 * Registration
 * ================================================================ */

void register_read_bam_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "read_bam");

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "region", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "reference", varchar_type);
    duckdb_destroy_logical_type(&varchar_type);

    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_table_function_add_named_parameter(tf, "standard_tags", bool_type);
    duckdb_table_function_add_named_parameter(tf, "auxiliary_tags", bool_type);
    duckdb_destroy_logical_type(&bool_type);

    duckdb_table_function_set_bind(tf, bam_read_bind);
    duckdb_table_function_set_init(tf, bam_read_global_init);
    duckdb_table_function_set_local_init(tf, bam_read_local_init);
    duckdb_table_function_set_function(tf, bam_read_function);
    duckdb_table_function_supports_projection_pushdown(tf, true);

    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}
