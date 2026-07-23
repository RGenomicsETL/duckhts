/**
 * DuckDB BCF/VCF Reader Extension
 *
 * A properly-typed VCF/BCF reader for DuckDB that matches the type system
 * of the nanoarrow vcf_arrow_stream implementation.
 *
 * Features:
 *   - VCF spec-compliant type validation with warnings
 *   - Proper DuckDB types: INT32, INT64, FLOAT, DOUBLE, VARCHAR, LIST, STRUCT
 *   - Boolean support for FLAG fields
 *   - Nullable fields with validity tracking
 *   - Parallel scan support for indexed files (CSI/TBI)
 *   - Region filtering
 *   - Projection pushdown
 *
 * Usage:
 *   LOAD 'bcf_reader.duckdb_extension';
 *   SELECT * FROM bcf_read('path/to/file.vcf.gz');
 *   SELECT * FROM bcf_read('path/to/file.bcf', region := 'chr1:1000-2000');
 *
 * Build:
 *   make (uses package htslib from RBCFTools)
 *
 * Copyright (c) 2026 RBCFTools Authors
 * Licensed under MIT License
 */

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include "include/vcf_types.h"
#include "include/vep_parser.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <inttypes.h>
#include <limits.h>
#include <strings.h>
#include <pthread.h>

// htslib headers
#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <htslib/hts_log.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>
#include <htslib/bgzf.h>

#include "include/hts_io_tuning.h"

// =============================================================================
// Constants
// =============================================================================

#define BCF_READER_DEFAULT_BATCH_SIZE 2048
#define VEP_TRANSCRIPT_ALL 0
#define VEP_TRANSCRIPT_FIRST 1

// Debug/progress tracking
#define BCF_READER_PROGRESS_INTERVAL 100000  // Print progress every N records
#define BCF_READER_ENABLE_PROGRESS 0
#define BCF_APPENDER_MAX_REGION_THREADS 16

// Decode policy for dirty/corrupt BCF header-vs-payload mismatches.
typedef enum {
    BCF_DECODE_ERROR_NULL = 0,  // Materialize NULL for decode failures.
    BCF_DECODE_ERROR_WARN = 1,  // Emit a DuckHTS warning, then materialize NULL.
    BCF_DECODE_ERROR_ERROR = 2  // Surface decode failures as DuckDB errors.
} bcf_decode_error_policy_t;

// Column indices for core VCF fields
enum {
    COL_CHROM = 0,
    COL_POS,
    COL_ID,
    COL_REF,
    COL_ALT,
    COL_QUAL,
    COL_FILTER,
    COL_CORE_COUNT  // Number of core columns (7)
};

// =============================================================================
// Field Metadata Structure
// =============================================================================

typedef struct {
    char* name;              // Field name (owned)
    int header_id;           // Header ID for bcf_get_* functions
    int header_type;         // BCF_HT_* from header (used for reading data)
    int schema_type;         // BCF_HT_* for schema (may be corrected)
    int vl_type;             // BCF_VL_* (corrected per VCF spec)
    int fixed_count;         // Exact fixed cardinality for Number=N fields, else <= 1
    int is_list;             // Whether this is a list type
    int duckdb_col_idx;      // Column index in DuckDB result
} field_meta_t;

typedef enum {
    BCF_OUT_INVALID = 0,
    BCF_OUT_CHROM,
    BCF_OUT_POS,
    BCF_OUT_ID,
    BCF_OUT_REF,
    BCF_OUT_ALT,
    BCF_OUT_QUAL,
    BCF_OUT_FILTER,
    BCF_OUT_VEP,
    BCF_OUT_INFO,
    BCF_OUT_SAMPLE_ID,
    BCF_OUT_FORMAT_GT,
    BCF_OUT_FORMAT_INT,
    BCF_OUT_FORMAT_FLOAT,
    BCF_OUT_FORMAT_STRING
} bcf_out_kind_t;

typedef struct {
    idx_t out_idx;                  // Projected output vector index in DuckDB chunk
    idx_t col_id;                   // Bound table-function column id
    bcf_out_kind_t kind;
    int field_idx;                  // INFO/FORMAT/VEP output field index, if applicable
    int schema_field_idx;           // VEP schema field index, if applicable
    int sample_idx;                 // Wide FORMAT sample index; -1 means current tidy sample
    int is_list;
    int header_type;
    const field_meta_t *field;
    const vep_field_t *vep_field;
} bcf_projected_col_t;

typedef struct {
    int field_idx;
    bcf_out_kind_t kind;
    idx_t column_count;
    idx_t *projected_indices;       // Indexes into bcf_init_data_t.projected_cols
} bcf_format_group_t;

// =============================================================================
// Bind Data - stores parameters and schema info
// =============================================================================

typedef struct {
    char* file_path;
    char* index_path;          // Optional explicit index path
    char* region;              // Optional region filter
    char* additional_csq_column_types;  // Optional bcftools-style override rules
    char* sample_filter_spec;  // Optional htslib sample include/exclude spec for v2
    int sample_filter_is_file; // Whether sample_filter_spec points at a file/list file
    int sample_filter_active;  // Whether FORMAT records need iterator subsetting
    int force_samples;         // Ignore unknown sample names in v2 selectors
    int decompression_threads; // htslib decompression worker threads per file handle
    int scan_sequential;       // Force full-file sequential streaming (no index count/parallel scan)
    int reader_version;        // 1 = legacy read_bcf, 2 = experimental read_bcf_v2 optimizations
    bcf_decode_error_policy_t decode_error_policy; // null|warn|error for BCF decode mismatches
    char** regions;            // Parsed comma-separated regions
    unsigned int n_regions;
    int include_info;          // Include INFO fields
    int include_format;        // Include FORMAT/sample fields
    int include_vep;           // Include VEP/CSQ/ANN/BCSQ fields
    int n_samples;             // Number of samples
    char** sample_names;       // Sample names (owned)

    // Tidy format options
    int tidy_format;           // If true, emit one row per variant-sample with SAMPLE_ID column
    int sample_id_col_idx;     // Column index for SAMPLE_ID (when tidy_format=true)

    // Field metadata
    int n_info_fields;
    field_meta_t* info_fields;

    int n_format_fields;
    field_meta_t* format_fields;

    // VEP/CSQ/BCSQ/ANN schema
    int n_vep_fields;
    int vep_col_start;       // Starting column index for VEP fields
    int *vep_field_indices;  // Output VEP field index -> schema field index
    vep_schema_t* vep_schema;
    int vep_transcript_mode; // VEP_TRANSCRIPT_FIRST (scalar) for now
    int info_col_start;
    int format_col_start;

    // Total column count
    int total_columns;

    // Parallel scan info (populated if index exists)
    int has_index;             // Whether an index was found
    int n_contigs;             // Number of contigs for parallel scan
    char** contig_names;       // Contig names (owned)
    uint64_t index_row_count;
    int index_row_count_valid;
} bcf_bind_data_t;

// =============================================================================
// Global Init Data - shared across all threads
// =============================================================================

typedef struct {
    volatile int current_contig;  // Next contig to assign (use atomic ops!)
    int n_contigs;                // Total number of contigs
    char** contig_names;          // Contig names (reference to bind data)
    int has_region;               // User specified a region
} bcf_global_init_data_t;

// =============================================================================
// Init Data - per-thread scanning state (now used as local init)
// =============================================================================

typedef struct {
    htsFile* fp;
    bcf_hdr_t* hdr;
    bcf1_t* rec;

    // Index support
    hts_idx_t* idx;           // BCF index (CSI)
    tbx_t* tbx;               // VCF tabix index (TBI)
    hts_itr_t* itr;           // Iterator
    kstring_t kstr;           // String buffer for VCF text parsing
    kstring_t gt_kstr;        // Thread-local reusable GT formatter buffer

    int64_t current_row;
    int done;

    // Projection pushdown
    idx_t column_count;
    idx_t* column_ids;
    bcf_projected_col_t* projected_cols;
    idx_t site_column_count;
    idx_t* site_column_indices;
    int n_format_groups;
    bcf_format_group_t* format_groups;
    int n_projected_vep_fields;
    int* projected_vep_field_indices;
    int unpack_mask;
    int need_vep;
    int count_only;
    uint64_t count_remaining;

    // Parallel scan state
    int is_parallel;
    int assigned_contig;       // Which contig this thread is scanning (-1 = all)
    const char* contig_name;   // Name of assigned contig (reference, don't free)
    int needs_next_contig;     // Flag to request next contig assignment
    // Tidy format state: tracks which sample we're emitting for current record
    int tidy_current_sample;   // Current sample index in tidy mode (-1 = need to read next record)
    int tidy_record_valid;     // Whether we have a valid record buffered for tidy mode

    // Debug/progress tracking
    int64_t total_records_processed;  // Total records processed by this thread
    struct timespec batch_start_time;  // Start time for performance measurement
    struct timespec last_progress_time;  // Last time progress was logged
    int timing_initialized;           // Flag to indicate timing is set up

    // Persistent decode cache for read_bcf_v2. The legacy read_bcf path keeps
    // per-output-chunk caches so v1 remains behaviorally isolated for benchmarks.
    int cache_persistent;
    int cache_n_format_fields;
    int cache_n_info_fields;
    int *fmt_loaded;
    int *fmt_ret;
    int *fmt_n_values;
    int32_t **fmt_i32;
    float **fmt_f32;
    char ***fmt_str;
    int gt_loaded;
    int gt_ret;
    int gt_n_values;
    int32_t *gt_arr;
    int *info_loaded;
    int *info_ret;
    int *info_n_values;
    int *info_flag;
    int32_t **info_i32;
    float **info_f32;
    char **info_str;
    vep_record_t *vep_rec;
    int vep_loaded;
} bcf_init_data_t;

// =============================================================================
// Warning Callback for DuckDB
// =============================================================================

static void duckdb_vcf_warning(const char* msg, void* ctx) {
    (void)ctx;
    // In DuckDB extensions, we can't easily emit warnings
    // For now, print to stderr
    fprintf(stderr, "[bcf_reader] %s\n", msg);
}

static int bcf_parse_decode_error_policy(const char *policy, bcf_decode_error_policy_t *out) {
    if (!out) return 0;
    *out = BCF_DECODE_ERROR_NULL;
    if (!policy || policy[0] == '\0' || strcasecmp(policy, "null") == 0) {
        return 1;
    }
    if (strcasecmp(policy, "warn") == 0) {
        *out = BCF_DECODE_ERROR_WARN;
        return 1;
    }
    if (strcasecmp(policy, "error") == 0) {
        *out = BCF_DECODE_ERROR_ERROR;
        return 1;
    }
    return 0;
}

static int bcf_parse_scan_mode(const char* mode, int* scan_sequential) {
    if (!scan_sequential) return 0;
    *scan_sequential = 0;
    if (!mode || mode[0] == '\0' || strcasecmp(mode, "auto") == 0) {
        return 1;
    }
    if (strcasecmp(mode, "sequential") == 0 ||
        strcasecmp(mode, "streaming") == 0 ||
        strcasecmp(mode, "stream") == 0 ||
        strcasecmp(mode, "seq") == 0) {
        *scan_sequential = 1;
        return 1;
    }
    return 0;
}

// =============================================================================
// Memory Management
// =============================================================================

static void free_bcf_format_string_array(char **values);

static void destroy_bind_data(void* data) {
    bcf_bind_data_t* bind = (bcf_bind_data_t*)data;
    if (!bind) return;

    if (bind->file_path) duckdb_free(bind->file_path);
    if (bind->index_path) duckdb_free(bind->index_path);
    if (bind->region) duckdb_free(bind->region);
    if (bind->additional_csq_column_types) duckdb_free(bind->additional_csq_column_types);
    if (bind->sample_filter_spec) duckdb_free(bind->sample_filter_spec);
    if (bind->regions) {
        for (unsigned int i = 0; i < bind->n_regions; i++) {
            if (bind->regions[i]) duckdb_free(bind->regions[i]);
        }
        duckdb_free(bind->regions);
    }

    if (bind->sample_names) {
        for (int i = 0; i < bind->n_samples; i++) {
            if (bind->sample_names[i]) duckdb_free(bind->sample_names[i]);
        }
        duckdb_free(bind->sample_names);
    }

    if (bind->info_fields) {
        for (int i = 0; i < bind->n_info_fields; i++) {
            if (bind->info_fields[i].name) duckdb_free(bind->info_fields[i].name);
        }
        duckdb_free(bind->info_fields);
    }

    if (bind->format_fields) {
        for (int i = 0; i < bind->n_format_fields; i++) {
            if (bind->format_fields[i].name) duckdb_free(bind->format_fields[i].name);
        }
        duckdb_free(bind->format_fields);
    }

    if (bind->contig_names) {
        for (int i = 0; i < bind->n_contigs; i++) {
            if (bind->contig_names[i]) duckdb_free(bind->contig_names[i]);
        }
        duckdb_free(bind->contig_names);
    }

    if (bind->vep_field_indices) duckdb_free(bind->vep_field_indices);
    if (bind->vep_schema) {
        vep_schema_destroy(bind->vep_schema);
    }

    duckdb_free(bind);
}

static void destroy_global_init_data(void* data) {
    bcf_global_init_data_t* global = (bcf_global_init_data_t*)data;
    if (global) {
        // contig_names is a reference, don't free
        duckdb_free(global);
    }
}

static void bcf_decode_cache_free(bcf_init_data_t *init) {
    if (!init) return;
    for (int i = 0; i < init->cache_n_format_fields; i++) {
        if (init->fmt_i32 && init->fmt_i32[i]) free(init->fmt_i32[i]);
        if (init->fmt_f32 && init->fmt_f32[i]) free(init->fmt_f32[i]);
        if (init->fmt_str && init->fmt_str[i]) free_bcf_format_string_array(init->fmt_str[i]);
    }
    if (init->gt_arr) free(init->gt_arr);
    if (init->fmt_loaded) duckdb_free(init->fmt_loaded);
    if (init->fmt_ret) duckdb_free(init->fmt_ret);
    if (init->fmt_n_values) duckdb_free(init->fmt_n_values);
    if (init->fmt_i32) duckdb_free(init->fmt_i32);
    if (init->fmt_f32) duckdb_free(init->fmt_f32);
    if (init->fmt_str) duckdb_free(init->fmt_str);

    for (int i = 0; i < init->cache_n_info_fields; i++) {
        if (init->info_i32 && init->info_i32[i]) free(init->info_i32[i]);
        if (init->info_f32 && init->info_f32[i]) free(init->info_f32[i]);
        if (init->info_str && init->info_str[i]) free(init->info_str[i]);
    }
    if (init->info_loaded) duckdb_free(init->info_loaded);
    if (init->info_ret) duckdb_free(init->info_ret);
    if (init->info_n_values) duckdb_free(init->info_n_values);
    if (init->info_flag) duckdb_free(init->info_flag);
    if (init->info_i32) duckdb_free(init->info_i32);
    if (init->info_f32) duckdb_free(init->info_f32);
    if (init->info_str) duckdb_free(init->info_str);
    if (init->vep_rec) vep_record_destroy(init->vep_rec);

    init->fmt_loaded = NULL;
    init->fmt_ret = NULL;
    init->fmt_n_values = NULL;
    init->fmt_i32 = NULL;
    init->fmt_f32 = NULL;
    init->fmt_str = NULL;
    init->gt_arr = NULL;
    init->gt_loaded = 0;
    init->gt_ret = 0;
    init->gt_n_values = 0;
    init->info_loaded = NULL;
    init->info_ret = NULL;
    init->info_n_values = NULL;
    init->info_flag = NULL;
    init->info_i32 = NULL;
    init->info_f32 = NULL;
    init->info_str = NULL;
    init->vep_rec = NULL;
    init->vep_loaded = 0;
    init->cache_n_format_fields = 0;
    init->cache_n_info_fields = 0;
    init->cache_persistent = 0;
}

static int bcf_decode_cache_init(bcf_init_data_t *init, const bcf_bind_data_t *bind) {
    if (!init || !bind) return 0;
    init->cache_persistent = 1;
    init->cache_n_format_fields = bind->n_format_fields;
    init->cache_n_info_fields = bind->n_info_fields;

    if (bind->n_format_fields > 0) {
        size_t nfmt = (size_t)bind->n_format_fields;
        init->fmt_loaded = (int *)duckdb_malloc(nfmt * sizeof(int));
        init->fmt_ret = (int *)duckdb_malloc(nfmt * sizeof(int));
        init->fmt_n_values = (int *)duckdb_malloc(nfmt * sizeof(int));
        init->fmt_i32 = (int32_t **)duckdb_malloc(nfmt * sizeof(int32_t *));
        init->fmt_f32 = (float **)duckdb_malloc(nfmt * sizeof(float *));
        init->fmt_str = (char ***)duckdb_malloc(nfmt * sizeof(char **));
        if (!init->fmt_loaded || !init->fmt_ret || !init->fmt_n_values ||
            !init->fmt_i32 || !init->fmt_f32 || !init->fmt_str) {
            bcf_decode_cache_free(init);
            return 0;
        }
        memset(init->fmt_loaded, 0, nfmt * sizeof(int));
        memset(init->fmt_ret, 0, nfmt * sizeof(int));
        memset(init->fmt_n_values, 0, nfmt * sizeof(int));
        memset(init->fmt_i32, 0, nfmt * sizeof(int32_t *));
        memset(init->fmt_f32, 0, nfmt * sizeof(float *));
        memset(init->fmt_str, 0, nfmt * sizeof(char **));
    }

    if (bind->n_info_fields > 0) {
        size_t ninfo = (size_t)bind->n_info_fields;
        init->info_loaded = (int *)duckdb_malloc(ninfo * sizeof(int));
        init->info_ret = (int *)duckdb_malloc(ninfo * sizeof(int));
        init->info_n_values = (int *)duckdb_malloc(ninfo * sizeof(int));
        init->info_flag = (int *)duckdb_malloc(ninfo * sizeof(int));
        init->info_i32 = (int32_t **)duckdb_malloc(ninfo * sizeof(int32_t *));
        init->info_f32 = (float **)duckdb_malloc(ninfo * sizeof(float *));
        init->info_str = (char **)duckdb_malloc(ninfo * sizeof(char *));
        if (!init->info_loaded || !init->info_ret || !init->info_n_values ||
            !init->info_flag || !init->info_i32 || !init->info_f32 || !init->info_str) {
            bcf_decode_cache_free(init);
            return 0;
        }
        memset(init->info_loaded, 0, ninfo * sizeof(int));
        memset(init->info_ret, 0, ninfo * sizeof(int));
        memset(init->info_n_values, 0, ninfo * sizeof(int));
        memset(init->info_flag, 0, ninfo * sizeof(int));
        memset(init->info_i32, 0, ninfo * sizeof(int32_t *));
        memset(init->info_f32, 0, ninfo * sizeof(float *));
        memset(init->info_str, 0, ninfo * sizeof(char *));
    }

    return 1;
}

static void destroy_init_data(void* data) {
    bcf_init_data_t* init = (bcf_init_data_t*)data;
    if (!init) return;

    bcf_decode_cache_free(init);
    if (init->itr) hts_itr_destroy(init->itr);
    if (init->tbx) tbx_destroy(init->tbx);
    if (init->idx) hts_idx_destroy(init->idx);
    if (init->rec) bcf_destroy(init->rec);
    if (init->hdr) bcf_hdr_destroy(init->hdr);
    if (init->fp) hts_close(init->fp);
    if (init->column_ids) duckdb_free(init->column_ids);
    if (init->projected_cols) duckdb_free(init->projected_cols);
    if (init->site_column_indices) duckdb_free(init->site_column_indices);
    if (init->projected_vep_field_indices) duckdb_free(init->projected_vep_field_indices);
    if (init->format_groups) {
        for (int i = 0; i < init->n_format_groups; i++) {
            if (init->format_groups[i].projected_indices) {
                duckdb_free(init->format_groups[i].projected_indices);
            }
        }
        duckdb_free(init->format_groups);
    }
    ks_free(&init->kstr);
    ks_free(&init->gt_kstr);

    duckdb_free(init);
}

// =============================================================================
// String Utilities
// =============================================================================

static char* strdup_duckdb(const char* s) {
    if (!s) return NULL;
    size_t len = strlen(s) + 1;
    char* copy = (char*)duckdb_malloc(len);
    if (copy) memcpy(copy, s, len);
    return copy;
}

static void parse_regions_duckdb(const char *region_str, char ***out_regions, unsigned int *out_count) {
    *out_regions = NULL;
    *out_count = 0;
    if (!region_str || region_str[0] == '\0') return;

    unsigned int count = 1;
    for (const char *p = region_str; *p; p++) {
        if (*p == ',') count++;
    }

    char **arr = (char **)duckdb_malloc(sizeof(char *) * count);
    char *dup = (char *)duckdb_malloc(strlen(region_str) + 1);
    if (!arr || !dup) {
        if (arr) duckdb_free(arr);
        if (dup) duckdb_free(dup);
        return;
    }
    strcpy(dup, region_str);

    unsigned int idx = 0;
    char *tok = strtok(dup, ",");
    while (tok && idx < count) {
        while (*tok == ' ' || *tok == '\t') tok++;
        int len = (int)strlen(tok);
        while (len > 0 && (tok[len - 1] == ' ' || tok[len - 1] == '\t')) {
            tok[--len] = '\0';
        }
        if (len > 0) {
            arr[idx] = strdup_duckdb(tok);
            if (!arr[idx]) {
                for (unsigned int i = 0; i < idx; i++) duckdb_free(arr[i]);
                duckdb_free(arr);
                duckdb_free(dup);
                return;
            }
            idx++;
        }
        tok = strtok(NULL, ",");
    }

    duckdb_free(dup);
    *out_regions = arr;
    *out_count = idx;
}

/* Open one htslib iterator for the complete requested region set.  HTSlib
 * 1.24's multi-region iterators merge overlapping chunks and return a record
 * once even when several requested regions overlap it. */
static hts_itr_t *bcf_open_region_iterator(hts_idx_t *idx, bcf_hdr_t *hdr,
                                           tbx_t *tbx, char **regions,
                                           unsigned int n_regions) {
    if (!regions || n_regions == 0) return NULL;

    if (idx) {
        return n_regions == 1
            ? bcf_itr_querys(idx, hdr, regions[0])
            : bcf_itr_regarray(idx, hdr, regions, n_regions);
    }
    if (tbx) {
        return n_regions == 1
            ? tbx_itr_querys(tbx, regions[0])
            : tbx_itr_regarray(tbx, regions, n_regions);
    }
    return NULL;
}

static void free_string_list_duckdb(char **items, int n_items) {
    if (!items) return;
    for (int i = 0; i < n_items; i++) {
        if (items[i]) duckdb_free(items[i]);
    }
    duckdb_free(items);
}

static void free_bcf_format_string_array(char **values) {
    if (!values) return;
    if (values[0]) free(values[0]);
    free(values);
}

static int parse_field_filter_list(const char *spec, const char *param_name,
                                   char ***out_items, int *out_n_items,
                                   char *err, size_t err_size) {
    *out_items = NULL;
    *out_n_items = 0;
    if (!spec) return 1;

    unsigned int max_items = 1;
    for (const char *p = spec; *p; p++) {
        if (*p == ',') max_items++;
    }

    char **items = (char **)duckdb_malloc(sizeof(char *) * max_items);
    if (!items) {
        snprintf(err, err_size, "read_bcf_v2: out of memory parsing %s", param_name);
        return 0;
    }
    memset(items, 0, sizeof(char *) * max_items);

    int n_items = 0;
    const char *start = spec;
    for (;;) {
        const char *end = start;
        while (*end && *end != ',') end++;

        const char *tok_start = start;
        const char *tok_end = end;
        while (tok_start < tok_end && (*tok_start == ' ' || *tok_start == '\t' || *tok_start == '\n' || *tok_start == '\r')) {
            tok_start++;
        }
        while (tok_end > tok_start && (tok_end[-1] == ' ' || tok_end[-1] == '\t' || tok_end[-1] == '\n' || tok_end[-1] == '\r')) {
            tok_end--;
        }

        size_t len = (size_t)(tok_end - tok_start);
        if (len == 0) {
            snprintf(err, err_size, "read_bcf_v2: %s must be a non-empty comma-separated field list", param_name);
            free_string_list_duckdb(items, n_items);
            return 0;
        }

        char *copy = (char *)duckdb_malloc(len + 1);
        if (!copy) {
            snprintf(err, err_size, "read_bcf_v2: out of memory parsing %s", param_name);
            free_string_list_duckdb(items, n_items);
            return 0;
        }
        memcpy(copy, tok_start, len);
        copy[len] = '\0';
        items[n_items++] = copy;

        if (*end == '\0') break;
        start = end + 1;
    }

    *out_items = items;
    *out_n_items = n_items;
    return 1;
}

static int field_filter_item_matches(const char *item, const char *name, const char *prefix) {
    if (!item || !name) return 0;
    if (strcmp(item, name) == 0) return 1;
    if (prefix) {
        size_t prefix_len = strlen(prefix);
        if (strncmp(item, prefix, prefix_len) == 0 && strcmp(item + prefix_len, name) == 0) {
            return 1;
        }
    }
    return 0;
}

static int field_filter_matches(char **items, int n_items, const char *name,
                                const char *prefix, int *matched) {
    if (n_items == 0) return 1;

    int found = 0;
    for (int i = 0; i < n_items; i++) {
        if (field_filter_item_matches(items[i], name, prefix)) {
            if (matched) matched[i] = 1;
            found = 1;
        }
    }
    return found;
}

static int field_filter_check_all_matched(char **items, int n_items, int *matched,
                                          const char *param_name, char *err, size_t err_size) {
    if (!items || !matched) return 1;
    for (int i = 0; i < n_items; i++) {
        if (!matched[i]) {
            snprintf(err, err_size, "read_bcf_v2: %s contains unknown field: %s", param_name, items[i]);
            return 0;
        }
    }
    return 1;
}

static int bcf_projection_unpack_mask(const bcf_bind_data_t* bind, const idx_t* column_ids, idx_t column_count) {
    int mask = 0;
    if (!bind || !column_ids) {
        return mask;
    }

    for (idx_t i = 0; i < column_count; i++) {
        idx_t col_id = column_ids[i];

        if (col_id == COL_ID || col_id == COL_REF || col_id == COL_ALT) {
            mask |= BCF_UN_STR;
            continue;
        }
        if (col_id == COL_FILTER) {
            mask |= BCF_UN_FLT;
            continue;
        }
        if (bind->vep_schema &&
            col_id >= (idx_t)bind->vep_col_start &&
            col_id < (idx_t)(bind->vep_col_start + bind->n_vep_fields)) {
            mask |= BCF_UN_INFO;
            continue;
        }
        if (col_id >= (idx_t)bind->info_col_start &&
            col_id < (idx_t)(bind->info_col_start + bind->n_info_fields)) {
            mask |= BCF_UN_INFO;
            continue;
        }
        if (col_id >= (idx_t)bind->format_col_start &&
            col_id < (idx_t)bind->total_columns) {
            mask |= BCF_UN_FMT;
        }
    }

    return mask;
}

static int bcf_projection_needs_vep(const bcf_bind_data_t* bind, const idx_t* column_ids, idx_t column_count) {
    if (!bind || !bind->vep_schema || !column_ids) {
        return 0;
    }

    for (idx_t i = 0; i < column_count; i++) {
        idx_t col_id = column_ids[i];
        if (col_id >= (idx_t)bind->vep_col_start &&
            col_id < (idx_t)(bind->vep_col_start + bind->n_vep_fields)) {
            return 1;
        }
    }

    return 0;
}

static bcf_out_kind_t bcf_projected_format_kind(const field_meta_t *field) {
    if (!field) {
        return BCF_OUT_INVALID;
    }
    if (field->header_type == BCF_HT_INT) {
        return BCF_OUT_FORMAT_INT;
    }
    if (field->header_type == BCF_HT_REAL) {
        return BCF_OUT_FORMAT_FLOAT;
    }
    if (field->name && strcmp(field->name, "GT") == 0) {
        return BCF_OUT_FORMAT_GT;
    }
    return BCF_OUT_FORMAT_STRING;
}

static int bcf_out_kind_is_format(bcf_out_kind_t kind) {
    return kind == BCF_OUT_FORMAT_GT ||
           kind == BCF_OUT_FORMAT_INT ||
           kind == BCF_OUT_FORMAT_FLOAT ||
           kind == BCF_OUT_FORMAT_STRING;
}

static int bcf_record_has_info_field(bcf1_t *rec, const field_meta_t *field) {
    return rec && field && bcf_get_info_id(rec, field->header_id) != NULL;
}

static int bcf_record_has_format_field(bcf1_t *rec, const field_meta_t *field) {
    return rec && field && bcf_get_fmt_id(rec, field->header_id) != NULL;
}

static const char *bcf_reader_name(const bcf_bind_data_t *bind) {
    return bind && bind->reader_version >= 2 ? "read_bcf_v2" : "read_bcf";
}

static int bcf_reader_input_is_bcf(htsFile *fp) {
    const htsFormat *fmt = fp ? hts_get_format(fp) : NULL;
    return fmt && fmt->format == bcf;
}

static const char *bcf_encoded_type_name(int type) {
    switch (type) {
    case BCF_BT_NULL: return "NULL";
    case BCF_BT_INT8: return "INT8";
    case BCF_BT_INT16: return "INT16";
    case BCF_BT_INT32: return "INT32";
    case BCF_BT_INT64: return "INT64";
    case BCF_BT_FLOAT: return "FLOAT";
    case BCF_BT_CHAR: return "CHAR";
    default: return "unknown";
    }
}

static const char *bcf_header_type_name(int type) {
    switch (type) {
    case BCF_HT_FLAG: return "Flag";
    case BCF_HT_INT: return "Integer";
    case BCF_HT_REAL: return "Float";
    case BCF_HT_STR: return "String";
    default: return "unknown";
    }
}

static void bcf_record_location(bcf_hdr_t *hdr, bcf1_t *rec, char *buf, size_t buf_size) {
    if (!buf || buf_size == 0) return;
    const char *chrom = (hdr && rec && rec->rid >= 0) ? bcf_hdr_id2name(hdr, rec->rid) : NULL;
    if (!chrom) chrom = "?";
    long long pos = rec ? (long long)rec->pos + 1 : 0;
    snprintf(buf, buf_size, "%s:%lld", chrom, pos);
}

static int bcf_encoded_type_matches_header(int header_type, int encoded_type) {
    switch (header_type) {
    case BCF_HT_INT:
        return encoded_type == BCF_BT_INT8 ||
               encoded_type == BCF_BT_INT16 ||
               encoded_type == BCF_BT_INT32;
    case BCF_HT_REAL:
        return encoded_type == BCF_BT_FLOAT;
    case BCF_HT_STR:
        return encoded_type == BCF_BT_CHAR;
    default:
        return 1;
    }
}

static int bcf_encoded_type_matches_format_kind(bcf_out_kind_t kind, int encoded_type) {
    switch (kind) {
    case BCF_OUT_FORMAT_GT:
    case BCF_OUT_FORMAT_INT:
        return encoded_type == BCF_BT_INT8 ||
               encoded_type == BCF_BT_INT16 ||
               encoded_type == BCF_BT_INT32;
    case BCF_OUT_FORMAT_FLOAT:
        return encoded_type == BCF_BT_FLOAT;
    case BCF_OUT_FORMAT_STRING:
        return encoded_type == BCF_BT_CHAR;
    default:
        return 1;
    }
}

static int bcf_check_decode_ret(const char *reader_name,
                                const char *field_class,
                                const char *tag,
                                bcf_hdr_t *hdr,
                                bcf1_t *rec,
                                int ret,
                                bcf_decode_error_policy_t policy,
                                char *err,
                                size_t err_size) {
    if (ret >= 0 || ret == -3) {
        return 1;
    }

    char loc[128];
    bcf_record_location(hdr, rec, loc, sizeof(loc));
    const char *name = reader_name ? reader_name : "read_bcf";
    const char *klass = field_class ? field_class : "field";
    const char *field = tag ? tag : "?";
    char msg[512];

    if (ret == -1) {
        snprintf(msg, sizeof(msg), "%s: %s/%s is not defined in the BCF/VCF header at %s", name, klass, field, loc);
    } else if (ret == -2) {
        snprintf(msg, sizeof(msg), "%s: %s/%s encoded BCF type does not match the header at %s", name, klass, field, loc);
        if (policy == BCF_DECODE_ERROR_NULL) {
            return 1;
        }
        if (policy == BCF_DECODE_ERROR_WARN) {
            vcf_emit_warning(msg);
            return 1;
        }
    } else if (ret == -4) {
        snprintf(msg, sizeof(msg), "%s: out of memory decoding %s/%s at %s", name, klass, field, loc);
    } else {
        snprintf(msg, sizeof(msg), "%s: failed to decode %s/%s at %s (htslib return %d)", name, klass, field, loc, ret);
    }
    if (err && err_size > 0) {
        snprintf(err, err_size, "%s", msg);
    }
    return 0;
}

static int bcf_handle_decode_diagnostic(bcf_decode_error_policy_t policy, const char *msg) {
    if (policy == BCF_DECODE_ERROR_ERROR) {
        return 0;
    }
    if (policy == BCF_DECODE_ERROR_WARN && msg && msg[0]) {
        vcf_emit_warning(msg);
    }
    return 1;
}

static int bcf_check_format_width(const char *reader_name,
                                  const char *tag,
                                  bcf_hdr_t *hdr,
                                  bcf1_t *rec,
                                  int ret,
                                  int n_samples,
                                  char *err,
                                  size_t err_size) {
    if (ret <= 0) {
        return 1;
    }
    if (n_samples <= 0 || ret % n_samples == 0) {
        return 1;
    }

    char loc[128];
    bcf_record_location(hdr, rec, loc, sizeof(loc));
    snprintf(err, err_size,
             "%s: FORMAT/%s decoded value count %d is not divisible by sample count %d at %s",
             reader_name ? reader_name : "read_bcf",
             tag ? tag : "?",
             ret, n_samples, loc);
    return 0;
}

static int bcf_handle_format_width(const char *reader_name,
                                   const char *tag,
                                   bcf_hdr_t *hdr,
                                   bcf1_t *rec,
                                   int *ret,
                                   int n_samples,
                                   bcf_decode_error_policy_t policy,
                                   char *err,
                                   size_t err_size) {
    if (!ret || bcf_check_format_width(reader_name, tag, hdr, rec, *ret, n_samples, err, err_size)) {
        return 1;
    }
    if (!bcf_handle_decode_diagnostic(policy, err)) {
        return 0;
    }
    *ret = -3;
    return 1;
}

static int bcf_preflight_info_encoded_type(bcf_hdr_t *hdr,
                                           bcf1_t *rec,
                                           const field_meta_t *field,
                                           const char *reader_name,
                                           char *err,
                                           size_t err_size) {
    if (!hdr || !rec || !field || !bcf_record_has_info_field(rec, field)) {
        return 1;
    }
    if (field->header_type == BCF_HT_FLAG) {
        return 1;
    }

    bcf_unpack(rec, BCF_UN_INFO);
    bcf_info_t *info = bcf_get_info_id(rec, field->header_id);
    if (!info || !info->vptr) {
        return 1;
    }
    if (bcf_encoded_type_matches_header(field->header_type, info->type)) {
        return 1;
    }

    char loc[128];
    bcf_record_location(hdr, rec, loc, sizeof(loc));
    snprintf(err, err_size,
             "%s: INFO/%s encoded BCF type %s does not match header Type=%s at %s",
             reader_name ? reader_name : "read_bcf",
             field->name ? field->name : "?",
             bcf_encoded_type_name(info->type),
             bcf_header_type_name(field->header_type),
             loc);
    return 0;
}

static int bcf_preflight_format_encoded_type(bcf_hdr_t *hdr,
                                             bcf1_t *rec,
                                             int header_id,
                                             const char *tag,
                                             bcf_out_kind_t kind,
                                             int header_type,
                                             const char *reader_name,
                                             char *err,
                                             size_t err_size) {
    if (!hdr || !rec || header_id < 0) {
        return 1;
    }

    bcf_unpack(rec, BCF_UN_FMT);
    bcf_fmt_t *fmt = bcf_get_fmt_id(rec, header_id);
    if (!fmt || !fmt->p) {
        return 1;
    }
    if (bcf_encoded_type_matches_format_kind(kind, fmt->type)) {
        return 1;
    }

    char loc[128];
    bcf_record_location(hdr, rec, loc, sizeof(loc));
    snprintf(err, err_size,
             "%s: FORMAT/%s encoded BCF type %s does not match header Type=%s at %s",
             reader_name ? reader_name : "read_bcf",
             tag ? tag : "?",
             bcf_encoded_type_name(fmt->type),
             kind == BCF_OUT_FORMAT_GT ? "String/GT-encoded-integer" : bcf_header_type_name(header_type),
             loc);
    return 0;
}

static int bcf_projected_columns_init(bcf_init_data_t *local, const bcf_bind_data_t *bind,
                                      char *err, size_t err_size) {
    if (!local || !bind || local->column_count == 0) {
        return 1;
    }

    const char *reader_name = bind->reader_version >= 2 ? "read_bcf_v2" : "read_bcf";

    local->projected_cols = (bcf_projected_col_t *)duckdb_malloc(
        (size_t)local->column_count * sizeof(bcf_projected_col_t));
    if (!local->projected_cols) {
        snprintf(err, err_size, "%s: out of memory allocating projected column descriptors", reader_name);
        return 0;
    }
    memset(local->projected_cols, 0, (size_t)local->column_count * sizeof(bcf_projected_col_t));

    int tidy_mode = bind->tidy_format && bind->include_format && bind->n_samples > 0;

    for (idx_t i = 0; i < local->column_count; i++) {
        idx_t col_id = local->column_ids[i];
        bcf_projected_col_t *col = &local->projected_cols[i];
        col->out_idx = i;
        col->col_id = col_id;
        col->field_idx = -1;
        col->schema_field_idx = -1;
        col->sample_idx = -1;
        col->kind = BCF_OUT_INVALID;

        switch (col_id) {
            case COL_CHROM:  col->kind = BCF_OUT_CHROM; break;
            case COL_POS:    col->kind = BCF_OUT_POS; break;
            case COL_ID:     col->kind = BCF_OUT_ID; break;
            case COL_REF:    col->kind = BCF_OUT_REF; break;
            case COL_ALT:    col->kind = BCF_OUT_ALT; break;
            case COL_QUAL:   col->kind = BCF_OUT_QUAL; break;
            case COL_FILTER: col->kind = BCF_OUT_FILTER; break;
            default: break;
        }
        if (col->kind != BCF_OUT_INVALID) {
            continue;
        }

        if (bind->vep_schema &&
            col_id >= (idx_t)bind->vep_col_start &&
            col_id < (idx_t)(bind->vep_col_start + bind->n_vep_fields)) {
            int field_idx = (int)(col_id - (idx_t)bind->vep_col_start);
            int schema_field_idx = bind->vep_field_indices ? bind->vep_field_indices[field_idx] : field_idx;
            const vep_field_t *field = vep_schema_get_field(bind->vep_schema, schema_field_idx);
            if (!field) {
                snprintf(err, err_size, "%s: invalid projected VEP field index", reader_name);
                return 0;
            }
            col->kind = BCF_OUT_VEP;
            col->field_idx = field_idx;
            col->schema_field_idx = schema_field_idx;
            col->is_list = 1;
            col->header_type = field->type;
            col->vep_field = field;
            continue;
        }

        if (col_id >= (idx_t)bind->info_col_start &&
            col_id < (idx_t)(bind->info_col_start + bind->n_info_fields)) {
            int field_idx = (int)(col_id - (idx_t)bind->info_col_start);
            const field_meta_t *field = &bind->info_fields[field_idx];
            col->kind = BCF_OUT_INFO;
            col->field_idx = field_idx;
            col->is_list = field->is_list;
            col->header_type = field->header_type;
            col->field = field;
            continue;
        }

        if (tidy_mode && col_id == (idx_t)bind->sample_id_col_idx) {
            col->kind = BCF_OUT_SAMPLE_ID;
            continue;
        }

        if (col_id >= (idx_t)bind->format_col_start && col_id < (idx_t)bind->total_columns) {
            int format_col_idx = (int)(col_id - (idx_t)bind->format_col_start);
            int field_idx;
            int sample_idx;

            if (bind->n_format_fields <= 0) {
                snprintf(err, err_size, "%s: projected FORMAT column without FORMAT fields", reader_name);
                return 0;
            }

            if (tidy_mode) {
                field_idx = format_col_idx;
                sample_idx = -1;
            } else {
                sample_idx = format_col_idx / bind->n_format_fields;
                field_idx = format_col_idx % bind->n_format_fields;
            }

            if (field_idx < 0 || field_idx >= bind->n_format_fields ||
                (!tidy_mode && (sample_idx < 0 || sample_idx >= bind->n_samples))) {
                snprintf(err, err_size, "%s: invalid projected FORMAT column index", reader_name);
                return 0;
            }

            const field_meta_t *field = &bind->format_fields[field_idx];
            col->kind = bcf_projected_format_kind(field);
            col->field_idx = field_idx;
            col->sample_idx = sample_idx;
            col->is_list = field->is_list;
            col->header_type = field->header_type;
            col->field = field;
            continue;
        }

        snprintf(err, err_size, "%s: invalid projected column index: %llu",
                 reader_name, (unsigned long long)col_id);
        return 0;
    }

    local->site_column_count = 0;
    local->n_format_groups = 0;
    local->n_projected_vep_fields = 0;

    int *format_counts = NULL;
    if (bind->n_format_fields > 0) {
        format_counts = (int *)duckdb_malloc((size_t)bind->n_format_fields * sizeof(int));
        if (!format_counts) {
            snprintf(err, err_size, "%s: out of memory allocating FORMAT projection counts", reader_name);
            return 0;
        }
        memset(format_counts, 0, (size_t)bind->n_format_fields * sizeof(int));
    }

    for (idx_t i = 0; i < local->column_count; i++) {
        bcf_projected_col_t *col = &local->projected_cols[i];
        if (bcf_out_kind_is_format(col->kind)) {
            if (col->field_idx >= 0 && col->field_idx < bind->n_format_fields) {
                if (format_counts[col->field_idx] == 0) {
                    local->n_format_groups++;
                }
                format_counts[col->field_idx]++;
            }
        } else {
            local->site_column_count++;
            if (col->kind == BCF_OUT_VEP) {
                int seen = 0;
                for (idx_t j = 0; j < i; j++) {
                    bcf_projected_col_t *prev = &local->projected_cols[j];
                    if (prev->kind == BCF_OUT_VEP && prev->schema_field_idx == col->schema_field_idx) {
                        seen = 1;
                        break;
                    }
                }
                if (!seen) {
                    local->n_projected_vep_fields++;
                }
            }
        }
    }

    if (local->n_projected_vep_fields > 0) {
        local->projected_vep_field_indices = (int *)duckdb_malloc(
            (size_t)local->n_projected_vep_fields * sizeof(int));
        if (!local->projected_vep_field_indices) {
            snprintf(err, err_size, "%s: out of memory allocating VEP projection descriptors", reader_name);
            if (format_counts) duckdb_free(format_counts);
            return 0;
        }
        int vep_write = 0;
        for (idx_t i = 0; i < local->column_count; i++) {
            bcf_projected_col_t *col = &local->projected_cols[i];
            if (col->kind != BCF_OUT_VEP) continue;
            int seen = 0;
            for (int j = 0; j < vep_write; j++) {
                if (local->projected_vep_field_indices[j] == col->schema_field_idx) {
                    seen = 1;
                    break;
                }
            }
            if (!seen) {
                local->projected_vep_field_indices[vep_write++] = col->schema_field_idx;
            }
        }
    }

    if (local->site_column_count > 0) {
        local->site_column_indices = (idx_t *)duckdb_malloc(
            (size_t)local->site_column_count * sizeof(idx_t));
        if (!local->site_column_indices) {
            snprintf(err, err_size, "%s: out of memory allocating site column descriptors", reader_name);
            if (format_counts) duckdb_free(format_counts);
            return 0;
        }
    }

    int *format_group_for_field = NULL;
    int *format_group_offsets = NULL;
    if (local->n_format_groups > 0) {
        local->format_groups = (bcf_format_group_t *)duckdb_malloc(
            (size_t)local->n_format_groups * sizeof(bcf_format_group_t));
        format_group_for_field = (int *)duckdb_malloc((size_t)bind->n_format_fields * sizeof(int));
        format_group_offsets = (int *)duckdb_malloc((size_t)local->n_format_groups * sizeof(int));
        if (!local->format_groups || !format_group_for_field || !format_group_offsets) {
            snprintf(err, err_size, "%s: out of memory allocating FORMAT projection groups", reader_name);
            if (format_counts) duckdb_free(format_counts);
            if (format_group_for_field) duckdb_free(format_group_for_field);
            if (format_group_offsets) duckdb_free(format_group_offsets);
            return 0;
        }
        memset(local->format_groups, 0, (size_t)local->n_format_groups * sizeof(bcf_format_group_t));
        for (int i = 0; i < bind->n_format_fields; i++) {
            format_group_for_field[i] = -1;
        }
        memset(format_group_offsets, 0, (size_t)local->n_format_groups * sizeof(int));

        int group_idx = 0;
        for (int field_idx = 0; field_idx < bind->n_format_fields; field_idx++) {
            if (format_counts[field_idx] == 0) {
                continue;
            }
            bcf_format_group_t *group = &local->format_groups[group_idx];
            group->field_idx = field_idx;
            group->kind = bcf_projected_format_kind(&bind->format_fields[field_idx]);
            group->column_count = (idx_t)format_counts[field_idx];
            group->projected_indices = (idx_t *)duckdb_malloc(
                (size_t)group->column_count * sizeof(idx_t));
            if (!group->projected_indices) {
                snprintf(err, err_size, "%s: out of memory allocating FORMAT projection group", reader_name);
                if (format_counts) duckdb_free(format_counts);
                if (format_group_for_field) duckdb_free(format_group_for_field);
                if (format_group_offsets) duckdb_free(format_group_offsets);
                return 0;
            }
            format_group_for_field[field_idx] = group_idx;
            group_idx++;
        }
    }

    idx_t site_write = 0;
    for (idx_t i = 0; i < local->column_count; i++) {
        bcf_projected_col_t *col = &local->projected_cols[i];
        if (bcf_out_kind_is_format(col->kind)) {
            int group_idx = format_group_for_field ? format_group_for_field[col->field_idx] : -1;
            if (group_idx >= 0) {
                int offset = format_group_offsets[group_idx]++;
                local->format_groups[group_idx].projected_indices[offset] = i;
            }
        } else if (local->site_column_indices) {
            local->site_column_indices[site_write++] = i;
        }
    }

    if (format_counts) duckdb_free(format_counts);
    if (format_group_for_field) duckdb_free(format_group_for_field);
    if (format_group_offsets) duckdb_free(format_group_offsets);
    return 1;
}

static int bcf_try_get_index_row_count(hts_idx_t *idx, uint64_t row_multiplier, uint64_t *out_total) {
    if (!idx || !out_total || row_multiplier == 0) {
        return 0;
    }

    int nseq = hts_idx_nseq(idx);
    uint64_t total = hts_idx_get_n_no_coor(idx);
    for (int tid = 0; tid < nseq; tid++) {
        uint64_t mapped = 0, unmapped = 0;
        if (hts_idx_get_stat(idx, tid, &mapped, &unmapped) != 0) {
            return 0;
        }
        total += mapped + unmapped;
    }

    *out_total = total * row_multiplier;
    return 1;
}

static inline void set_validity_bit(uint64_t* validity, idx_t row, int is_valid);
static void process_comma_separated_list(duckdb_vector vec, idx_t row, const char* value);

// =============================================================================
// DuckDB Type Creation Helpers
// =============================================================================

/**
 * Create a DuckDB logical type for a BCF field.
 */
static duckdb_logical_type create_bcf_field_type(int bcf_type, int is_list) {
    duckdb_logical_type element_type;

    switch (bcf_type) {
        case BCF_HT_FLAG:
            element_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
            break;
        case BCF_HT_INT:
            element_type = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
            break;
        case BCF_HT_REAL:
            element_type = duckdb_create_logical_type(DUCKDB_TYPE_FLOAT);
            break;
        case BCF_HT_STR:
        default:
            element_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
            break;
    }

    if (is_list) {
        duckdb_logical_type list_type = duckdb_create_list_type(element_type);
        duckdb_destroy_logical_type(&element_type);
        return list_type;
    }

    return element_type;
}

static duckdb_logical_type create_vep_field_type(vep_field_type_t vep_type, int is_list) {
    duckdb_logical_type element_type;

    switch (vep_type) {
        case VEP_TYPE_INTEGER:
            element_type = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
            break;
        case VEP_TYPE_FLOAT:
            element_type = duckdb_create_logical_type(DUCKDB_TYPE_FLOAT);
            break;
        case VEP_TYPE_FLAG:
            element_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
            break;
        case VEP_TYPE_STRING:
        default:
            element_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
            break;
    }

    if (is_list) {
        duckdb_logical_type list_type = duckdb_create_list_type(element_type);
        duckdb_destroy_logical_type(&element_type);
        return list_type;
    }

    return element_type;
}

static void set_field_null(duckdb_vector vec, idx_t row, int is_list) {
    duckdb_vector_ensure_validity_writable(vec);
    uint64_t* validity = duckdb_vector_get_validity(vec);
    set_validity_bit(validity, row, 0);

    if (is_list) {
        duckdb_list_entry entry = {duckdb_list_vector_get_size(vec), 0};
        duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
        list_data[row] = entry;
    }
}

static void assign_int32_field(duckdb_vector vec, idx_t row, const int32_t* values, int n_values, int is_list) {
    if (!values || n_values <= 0) {
        set_field_null(vec, row, is_list);
        return;
    }

    if (!is_list) {
        if (values[0] != bcf_int32_missing && values[0] != bcf_int32_vector_end) {
            int32_t* data = (int32_t*)duckdb_vector_get_data(vec);
            data[row] = values[0];
        } else {
            set_field_null(vec, row, 0);
        }
        return;
    }

    duckdb_list_entry entry;
    entry.offset = duckdb_list_vector_get_size(vec);
    entry.length = 0;

    for (int i = 0; i < n_values; i++) {
        if (values[i] != bcf_int32_missing && values[i] != bcf_int32_vector_end) {
            entry.length++;
        }
    }

    if (entry.length > 0) {
        duckdb_list_vector_reserve(vec, entry.offset + entry.length);
        duckdb_list_vector_set_size(vec, entry.offset + entry.length);

        duckdb_vector child_vec = duckdb_list_vector_get_child(vec);
        int32_t* child_data = (int32_t*)duckdb_vector_get_data(child_vec);
        int write_idx = 0;
        for (int i = 0; i < n_values; i++) {
            if (values[i] != bcf_int32_missing && values[i] != bcf_int32_vector_end) {
                child_data[entry.offset + write_idx] = values[i];
                write_idx++;
            }
        }
    }

    duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
    list_data[row] = entry;
}

static void assign_float_field(duckdb_vector vec, idx_t row, const float* values, int n_values, int is_list) {
    if (!values || n_values <= 0) {
        set_field_null(vec, row, is_list);
        return;
    }

    if (!is_list) {
        if (!bcf_float_is_missing(values[0]) && !bcf_float_is_vector_end(values[0])) {
            float* data = (float*)duckdb_vector_get_data(vec);
            data[row] = values[0];
        } else {
            set_field_null(vec, row, 0);
        }
        return;
    }

    duckdb_list_entry entry;
    entry.offset = duckdb_list_vector_get_size(vec);
    entry.length = 0;

    for (int i = 0; i < n_values; i++) {
        if (!bcf_float_is_missing(values[i]) && !bcf_float_is_vector_end(values[i])) {
            entry.length++;
        }
    }

    if (entry.length > 0) {
        duckdb_list_vector_reserve(vec, entry.offset + entry.length);
        duckdb_list_vector_set_size(vec, entry.offset + entry.length);

        duckdb_vector child_vec = duckdb_list_vector_get_child(vec);
        float* child_data = (float*)duckdb_vector_get_data(child_vec);
        int write_idx = 0;
        for (int i = 0; i < n_values; i++) {
            if (!bcf_float_is_missing(values[i]) && !bcf_float_is_vector_end(values[i])) {
                child_data[entry.offset + write_idx] = values[i];
                write_idx++;
            }
        }
    }

    duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
    list_data[row] = entry;
}

static void assign_string_field(duckdb_vector vec, idx_t row, const char* value, int is_list) {
    if (!value || strcmp(value, ".") == 0) {
        set_field_null(vec, row, is_list);
        return;
    }

    if (is_list) {
        process_comma_separated_list(vec, row, value);
    } else {
        duckdb_vector_assign_string_element(vec, row, value);
    }
}

static int format_gt_fast(const int32_t *gt, int max_ploidy,
                          char small[4], const char **out, size_t *out_len) {
    if (!gt || max_ploidy <= 0 || !small || !out || !out_len) {
        return 0;
    }

    int n = 0;
    while (n < max_ploidy && gt[n] != bcf_int32_vector_end) {
        n++;
    }
    if (n == 0 || n > 2) {
        return 0;
    }

    for (int i = 0; i < n; i++) {
        if (!bcf_gt_is_missing(gt[i])) {
            int allele = bcf_gt_allele(gt[i]);
            if (allele < 0 || allele > 9) {
                return 0;
            }
        }
    }

    small[0] = bcf_gt_is_missing(gt[0]) ? '.' : (char)('0' + bcf_gt_allele(gt[0]));
    if (n == 1) {
        *out = small;
        *out_len = 1;
        return 1;
    }

    small[1] = bcf_gt_is_phased(gt[1]) ? '|' : '/';
    small[2] = bcf_gt_is_missing(gt[1]) ? '.' : (char)('0' + bcf_gt_allele(gt[1]));
    *out = small;
    *out_len = 3;
    return 1;
}

// =============================================================================
// Schema Building - Bind Function
// =============================================================================

static void bcf_read_bind_common(duckdb_bind_info info, int reader_version) {
    const char *reader_name = reader_version >= 2 ? "read_bcf_v2" : "read_bcf";

    // Set up warning callback
    vcf_set_warning_callback(duckdb_vcf_warning, NULL);

    // Get the file path parameter
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    char* file_path = duckdb_get_varchar(path_val);
    duckdb_destroy_value(&path_val);

    if (!file_path || strlen(file_path) == 0) {
        char err[96];
        snprintf(err, sizeof(err), "%s requires a file path", reader_name);
        duckdb_bind_set_error(info, err);
        if (file_path) duckdb_free(file_path);
        return;
    }

    // Get optional region named parameter
    char* region = NULL;
    duckdb_value region_val = duckdb_bind_get_named_parameter(info, "region");
    if (region_val && !duckdb_is_null_value(region_val)) {
        region = duckdb_get_varchar(region_val);
    }
    if (region_val) duckdb_destroy_value(&region_val);

    // Optional explicit index path
    char* index_path = NULL;
    duckdb_value idx_val = duckdb_bind_get_named_parameter(info, "index_path");
    if (idx_val && !duckdb_is_null_value(idx_val)) {
        index_path = duckdb_get_varchar(idx_val);
    }
    if (idx_val) duckdb_destroy_value(&idx_val);

    // Optional scan mode: auto (current index-aware behavior) or sequential streaming.
    int scan_sequential = 0;
    duckdb_value scan_mode_val = duckdb_bind_get_named_parameter(info, "scan_mode");
    if (scan_mode_val && !duckdb_is_null_value(scan_mode_val)) {
        char* scan_mode = duckdb_get_varchar(scan_mode_val);
        if (!bcf_parse_scan_mode(scan_mode, &scan_sequential)) {
            char err[128];
            snprintf(err, sizeof(err), "%s: scan_mode must be 'auto' or 'sequential'", reader_name);
            duckdb_bind_set_error(info, err);
            if (scan_mode) duckdb_free(scan_mode);
            duckdb_destroy_value(&scan_mode_val);
            duckdb_free(file_path);
            if (index_path) duckdb_free(index_path);
            if (region) duckdb_free(region);
            return;
        }
        if (scan_mode) duckdb_free(scan_mode);
    }
    if (scan_mode_val) duckdb_destroy_value(&scan_mode_val);
    if (scan_sequential && region && region[0] != '\0') {
        char err[160];
        snprintf(err, sizeof(err), "%s: scan_mode := 'sequential' is incompatible with region queries", reader_name);
        duckdb_bind_set_error(info, err);
        duckdb_free(file_path);
        if (index_path) duckdb_free(index_path);
        if (region) duckdb_free(region);
        return;
    }

    bcf_decode_error_policy_t decode_error_policy = BCF_DECODE_ERROR_NULL;
    duckdb_value decode_policy_val = duckdb_bind_get_named_parameter(info, "decode_error_policy");
    if (decode_policy_val && !duckdb_is_null_value(decode_policy_val)) {
        char *decode_policy = duckdb_get_varchar(decode_policy_val);
        if (!bcf_parse_decode_error_policy(decode_policy, &decode_error_policy)) {
            char err[192];
            snprintf(err, sizeof(err), "%s: decode_error_policy must be 'null', 'warn', or 'error'", reader_name);
            duckdb_bind_set_error(info, err);
            if (decode_policy) duckdb_free(decode_policy);
            duckdb_destroy_value(&decode_policy_val);
            duckdb_free(file_path);
            if (index_path) duckdb_free(index_path);
            if (region) duckdb_free(region);
            return;
        }
        if (decode_policy) duckdb_free(decode_policy);
    }
    if (decode_policy_val) duckdb_destroy_value(&decode_policy_val);

    // Get optional tidy_format named parameter (default: false)
    int tidy_format = 0;
    duckdb_value tidy_val = duckdb_bind_get_named_parameter(info, "tidy_format");
    if (tidy_val && !duckdb_is_null_value(tidy_val)) {
        tidy_format = duckdb_get_bool(tidy_val);
    }
    if (tidy_val) duckdb_destroy_value(&tidy_val);

    // Optional bcftools-style CSQ/ANN/BCSQ type overrides
    char* additional_csq_column_types = NULL;
    int decompression_threads = 0;
    duckdb_value csq_types_val = duckdb_bind_get_named_parameter(info, "additional_csq_column_types");
    if (csq_types_val && !duckdb_is_null_value(csq_types_val)) {
        additional_csq_column_types = duckdb_get_varchar(csq_types_val);
    }
    if (csq_types_val) duckdb_destroy_value(&csq_types_val);
    if (additional_csq_column_types) {
        char err[256];
        if (!vep_validate_column_type_rules(additional_csq_column_types, err, sizeof(err))) {
            duckdb_bind_set_error(info, err);
            duckdb_free(file_path);
            if (index_path) duckdb_free(index_path);
            if (region) duckdb_free(region);
            duckdb_free(additional_csq_column_types);
            return;
        }
    }

    duckdb_value dthreads_val = duckdb_bind_get_named_parameter(info, "decompression_threads");
    if (dthreads_val && !duckdb_is_null_value(dthreads_val)) {
        int64_t requested_threads = duckdb_get_int64(dthreads_val);
        if (requested_threads < 0 || requested_threads > INT_MAX) {
            char err[160];
            snprintf(err, sizeof(err), "%s: decompression_threads must be between 0 and INT_MAX", reader_name);
            duckdb_bind_set_error(info, err);
            duckdb_free(file_path);
            if (index_path) duckdb_free(index_path);
            if (region) duckdb_free(region);
            if (additional_csq_column_types) duckdb_free(additional_csq_column_types);
            duckdb_destroy_value(&dthreads_val);
            return;
        }
        decompression_threads = (int)requested_threads;
    }
    if (dthreads_val) duckdb_destroy_value(&dthreads_val);

    // Optional read_bcf_v2 sample pushdown. These are intentionally v2-only so
    // legacy read_bcf remains a stable benchmark/control surface.
    char *sample_filter_spec = NULL;
    int sample_filter_is_file = 0;
    int sample_filter_active = 0;
    int force_samples = 0;
    if (reader_version >= 2) {
        char *samples = NULL;
        char *samples_file = NULL;
        char *exclude_samples = NULL;
        char *exclude_samples_file = NULL;

        duckdb_value samples_val = duckdb_bind_get_named_parameter(info, "samples");
        if (samples_val && !duckdb_is_null_value(samples_val)) samples = duckdb_get_varchar(samples_val);
        if (samples_val) duckdb_destroy_value(&samples_val);

        duckdb_value samples_file_val = duckdb_bind_get_named_parameter(info, "samples_file");
        if (samples_file_val && !duckdb_is_null_value(samples_file_val)) samples_file = duckdb_get_varchar(samples_file_val);
        if (samples_file_val) duckdb_destroy_value(&samples_file_val);

        duckdb_value exclude_samples_val = duckdb_bind_get_named_parameter(info, "exclude_samples");
        if (exclude_samples_val && !duckdb_is_null_value(exclude_samples_val)) exclude_samples = duckdb_get_varchar(exclude_samples_val);
        if (exclude_samples_val) duckdb_destroy_value(&exclude_samples_val);

        duckdb_value exclude_samples_file_val = duckdb_bind_get_named_parameter(info, "exclude_samples_file");
        if (exclude_samples_file_val && !duckdb_is_null_value(exclude_samples_file_val)) exclude_samples_file = duckdb_get_varchar(exclude_samples_file_val);
        if (exclude_samples_file_val) duckdb_destroy_value(&exclude_samples_file_val);

        duckdb_value force_samples_val = duckdb_bind_get_named_parameter(info, "force_samples");
        if (force_samples_val && !duckdb_is_null_value(force_samples_val)) force_samples = duckdb_get_bool(force_samples_val);
        if (force_samples_val) duckdb_destroy_value(&force_samples_val);

        int selector_count = 0;
        if (samples) selector_count++;
        if (samples_file) selector_count++;
        if (exclude_samples) selector_count++;
        if (exclude_samples_file) selector_count++;
        if (selector_count > 1) {
            duckdb_bind_set_error(info, "read_bcf_v2: specify only one of samples, samples_file, exclude_samples, exclude_samples_file");
            duckdb_free(file_path);
            if (index_path) duckdb_free(index_path);
            if (region) duckdb_free(region);
            if (additional_csq_column_types) duckdb_free(additional_csq_column_types);
            if (samples) duckdb_free(samples);
            if (samples_file) duckdb_free(samples_file);
            if (exclude_samples) duckdb_free(exclude_samples);
            if (exclude_samples_file) duckdb_free(exclude_samples_file);
            return;
        }

        const char *raw_sample_spec = NULL;
        int is_exclusion = 0;
        if (samples) {
            raw_sample_spec = samples;
            sample_filter_is_file = 0;
        } else if (samples_file) {
            raw_sample_spec = samples_file;
            sample_filter_is_file = 1;
        } else if (exclude_samples) {
            raw_sample_spec = exclude_samples;
            sample_filter_is_file = 0;
            is_exclusion = 1;
        } else if (exclude_samples_file) {
            raw_sample_spec = exclude_samples_file;
            sample_filter_is_file = 1;
            is_exclusion = 1;
        }

        if (raw_sample_spec) {
            if (strlen(raw_sample_spec) == 0) {
                duckdb_bind_set_error(info, "read_bcf_v2: sample selectors must not be empty");
                duckdb_free(file_path);
                if (index_path) duckdb_free(index_path);
                if (region) duckdb_free(region);
                if (additional_csq_column_types) duckdb_free(additional_csq_column_types);
                if (samples) duckdb_free(samples);
                if (samples_file) duckdb_free(samples_file);
                if (exclude_samples) duckdb_free(exclude_samples);
                if (exclude_samples_file) duckdb_free(exclude_samples_file);
                return;
            }
            size_t raw_len = strlen(raw_sample_spec);
            sample_filter_spec = (char *)duckdb_malloc(raw_len + (is_exclusion ? 2 : 1));
            if (!sample_filter_spec) {
                duckdb_bind_set_error(info, "read_bcf_v2: out of memory allocating sample selector");
                duckdb_free(file_path);
                if (index_path) duckdb_free(index_path);
                if (region) duckdb_free(region);
                if (additional_csq_column_types) duckdb_free(additional_csq_column_types);
                if (samples) duckdb_free(samples);
                if (samples_file) duckdb_free(samples_file);
                if (exclude_samples) duckdb_free(exclude_samples);
                if (exclude_samples_file) duckdb_free(exclude_samples_file);
                return;
            }
            if (is_exclusion) {
                sample_filter_spec[0] = '^';
                memcpy(sample_filter_spec + 1, raw_sample_spec, raw_len + 1);
            } else {
                memcpy(sample_filter_spec, raw_sample_spec, raw_len + 1);
            }
            sample_filter_active = strcmp(sample_filter_spec, "-") != 0;
        }

        if (samples) duckdb_free(samples);
        if (samples_file) duckdb_free(samples_file);
        if (exclude_samples) duckdb_free(exclude_samples);
        if (exclude_samples_file) duckdb_free(exclude_samples_file);
    }

    int include_info = 1;
    int include_format = 1;
    int include_vep = 1;
    char *info_fields_spec = NULL;
    char *format_fields_spec = NULL;
    char *vep_fields_spec = NULL;
    char **info_field_filter = NULL;
    char **format_field_filter = NULL;
    char **vep_field_filter = NULL;
    int n_info_field_filter = 0;
    int n_format_field_filter = 0;
    int n_vep_field_filter = 0;

    if (reader_version >= 2) {
        duckdb_value include_info_val = duckdb_bind_get_named_parameter(info, "include_info");
        if (include_info_val && !duckdb_is_null_value(include_info_val)) include_info = duckdb_get_bool(include_info_val);
        if (include_info_val) duckdb_destroy_value(&include_info_val);

        duckdb_value include_format_val = duckdb_bind_get_named_parameter(info, "include_format");
        if (include_format_val && !duckdb_is_null_value(include_format_val)) include_format = duckdb_get_bool(include_format_val);
        if (include_format_val) duckdb_destroy_value(&include_format_val);

        duckdb_value include_vep_val = duckdb_bind_get_named_parameter(info, "include_vep");
        if (include_vep_val && !duckdb_is_null_value(include_vep_val)) include_vep = duckdb_get_bool(include_vep_val);
        if (include_vep_val) duckdb_destroy_value(&include_vep_val);

        duckdb_value info_fields_val = duckdb_bind_get_named_parameter(info, "info_fields");
        if (info_fields_val && !duckdb_is_null_value(info_fields_val)) info_fields_spec = duckdb_get_varchar(info_fields_val);
        if (info_fields_val) duckdb_destroy_value(&info_fields_val);

        duckdb_value format_fields_val = duckdb_bind_get_named_parameter(info, "format_fields");
        if (format_fields_val && !duckdb_is_null_value(format_fields_val)) format_fields_spec = duckdb_get_varchar(format_fields_val);
        if (format_fields_val) duckdb_destroy_value(&format_fields_val);

        duckdb_value vep_fields_val = duckdb_bind_get_named_parameter(info, "vep_fields");
        if (vep_fields_val && !duckdb_is_null_value(vep_fields_val)) vep_fields_spec = duckdb_get_varchar(vep_fields_val);
        if (vep_fields_val) duckdb_destroy_value(&vep_fields_val);

        if (!include_info && info_fields_spec) {
            duckdb_bind_set_error(info, "read_bcf_v2: info_fields requires include_info := true");
            duckdb_free(file_path);
            if (index_path) duckdb_free(index_path);
            if (region) duckdb_free(region);
            if (additional_csq_column_types) duckdb_free(additional_csq_column_types);
            if (sample_filter_spec) duckdb_free(sample_filter_spec);
            if (info_fields_spec) duckdb_free(info_fields_spec);
            if (format_fields_spec) duckdb_free(format_fields_spec);
            if (vep_fields_spec) duckdb_free(vep_fields_spec);
            return;
        }
        if (!include_format && format_fields_spec) {
            duckdb_bind_set_error(info, "read_bcf_v2: format_fields requires include_format := true");
            duckdb_free(file_path);
            if (index_path) duckdb_free(index_path);
            if (region) duckdb_free(region);
            if (additional_csq_column_types) duckdb_free(additional_csq_column_types);
            if (sample_filter_spec) duckdb_free(sample_filter_spec);
            if (info_fields_spec) duckdb_free(info_fields_spec);
            if (format_fields_spec) duckdb_free(format_fields_spec);
            if (vep_fields_spec) duckdb_free(vep_fields_spec);
            return;
        }
        if (!include_vep && vep_fields_spec) {
            duckdb_bind_set_error(info, "read_bcf_v2: vep_fields requires include_vep := true");
            duckdb_free(file_path);
            if (index_path) duckdb_free(index_path);
            if (region) duckdb_free(region);
            if (additional_csq_column_types) duckdb_free(additional_csq_column_types);
            if (sample_filter_spec) duckdb_free(sample_filter_spec);
            if (info_fields_spec) duckdb_free(info_fields_spec);
            if (format_fields_spec) duckdb_free(format_fields_spec);
            if (vep_fields_spec) duckdb_free(vep_fields_spec);
            return;
        }

        char filter_err[256];
        if ((info_fields_spec && !parse_field_filter_list(info_fields_spec, "info_fields", &info_field_filter, &n_info_field_filter, filter_err, sizeof(filter_err))) ||
            (format_fields_spec && !parse_field_filter_list(format_fields_spec, "format_fields", &format_field_filter, &n_format_field_filter, filter_err, sizeof(filter_err))) ||
            (vep_fields_spec && !parse_field_filter_list(vep_fields_spec, "vep_fields", &vep_field_filter, &n_vep_field_filter, filter_err, sizeof(filter_err)))) {
            duckdb_bind_set_error(info, filter_err);
            duckdb_free(file_path);
            if (index_path) duckdb_free(index_path);
            if (region) duckdb_free(region);
            if (additional_csq_column_types) duckdb_free(additional_csq_column_types);
            if (sample_filter_spec) duckdb_free(sample_filter_spec);
            if (info_fields_spec) duckdb_free(info_fields_spec);
            if (format_fields_spec) duckdb_free(format_fields_spec);
            if (vep_fields_spec) duckdb_free(vep_fields_spec);
            free_string_list_duckdb(info_field_filter, n_info_field_filter);
            free_string_list_duckdb(format_field_filter, n_format_field_filter);
            free_string_list_duckdb(vep_field_filter, n_vep_field_filter);
            return;
        }
    }

    if (info_fields_spec) duckdb_free(info_fields_spec);
    if (format_fields_spec) duckdb_free(format_fields_spec);
    if (vep_fields_spec) duckdb_free(vep_fields_spec);

    // Open the file to read header
    htsFile* fp = hts_open(file_path, "r");
    if (!fp) {
        char err[512];
        snprintf(err, sizeof(err), "%s: failed to open BCF/VCF file: %s", reader_name, file_path);
        duckdb_bind_set_error(info, err);
        duckdb_free(file_path);
        if (index_path) duckdb_free(index_path);
        if (region) duckdb_free(region);
        if (additional_csq_column_types) duckdb_free(additional_csq_column_types);
        if (sample_filter_spec) duckdb_free(sample_filter_spec);
        free_string_list_duckdb(info_field_filter, n_info_field_filter);
        free_string_list_duckdb(format_field_filter, n_format_field_filter);
        free_string_list_duckdb(vep_field_filter, n_vep_field_filter);
        return;
    }

    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        char err[128];
        hts_close(fp);
        snprintf(err, sizeof(err), "%s: failed to read BCF/VCF header", reader_name);
        duckdb_bind_set_error(info, err);
        duckdb_free(file_path);
        if (index_path) duckdb_free(index_path);
        if (region) duckdb_free(region);
        if (additional_csq_column_types) duckdb_free(additional_csq_column_types);
        if (sample_filter_spec) duckdb_free(sample_filter_spec);
        free_string_list_duckdb(info_field_filter, n_info_field_filter);
        free_string_list_duckdb(format_field_filter, n_format_field_filter);
        free_string_list_duckdb(vep_field_filter, n_vep_field_filter);
        return;
    }

    if (reader_version >= 2 && sample_filter_spec) {
        int sample_ret = bcf_hdr_set_samples(hdr, sample_filter_spec, sample_filter_is_file);
        if (sample_ret < 0 || (sample_ret > 0 && !force_samples)) {
            char err[512];
            if (sample_ret > 0) {
                snprintf(err, sizeof(err),
                         "read_bcf_v2: sample selector contains an unknown sample at position %d; set force_samples := true to ignore unknown names",
                         sample_ret);
            } else {
                snprintf(err, sizeof(err), "read_bcf_v2: failed to apply sample selector");
            }
            bcf_hdr_destroy(hdr);
            hts_close(fp);
            duckdb_bind_set_error(info, err);
            duckdb_free(file_path);
            if (index_path) duckdb_free(index_path);
            if (region) duckdb_free(region);
            if (additional_csq_column_types) duckdb_free(additional_csq_column_types);
            if (sample_filter_spec) duckdb_free(sample_filter_spec);
            free_string_list_duckdb(info_field_filter, n_info_field_filter);
            free_string_list_duckdb(format_field_filter, n_format_field_filter);
            free_string_list_duckdb(vep_field_filter, n_vep_field_filter);
            return;
        }
    }

    if (reader_version >= 2 && include_format && n_format_field_filter > 0 && bcf_hdr_nsamples(hdr) == 0) {
        duckdb_bind_set_error(info, "read_bcf_v2: format_fields specified but no samples are selected");
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        duckdb_free(file_path);
        if (index_path) duckdb_free(index_path);
        if (region) duckdb_free(region);
        if (additional_csq_column_types) duckdb_free(additional_csq_column_types);
        if (sample_filter_spec) duckdb_free(sample_filter_spec);
        free_string_list_duckdb(info_field_filter, n_info_field_filter);
        free_string_list_duckdb(format_field_filter, n_format_field_filter);
        free_string_list_duckdb(vep_field_filter, n_vep_field_filter);
        return;
    }

    // Create bind data
    bcf_bind_data_t* bind = (bcf_bind_data_t*)duckdb_malloc(sizeof(bcf_bind_data_t));
    memset(bind, 0, sizeof(bcf_bind_data_t));
    bind->file_path = file_path;
    bind->index_path = index_path;
    bind->region = region;
    bind->additional_csq_column_types = additional_csq_column_types;
    bind->sample_filter_spec = sample_filter_spec;
    bind->sample_filter_is_file = sample_filter_is_file;
    bind->sample_filter_active = sample_filter_active;
    bind->force_samples = force_samples;
    bind->decompression_threads = decompression_threads;
    bind->scan_sequential = scan_sequential;
    bind->reader_version = reader_version;
    bind->decode_error_policy = decode_error_policy;
    parse_regions_duckdb(region, &bind->regions, &bind->n_regions);
    bind->include_info = include_info;
    bind->include_format = include_format;
    bind->include_vep = include_vep;
    bind->n_samples = bcf_hdr_nsamples(hdr);
    bind->tidy_format = tidy_format;
    bind->sample_id_col_idx = -1;  // Will be set if tidy_format=true
    bind->n_vep_fields = 0;
    bind->vep_col_start = COL_CORE_COUNT;
    bind->info_col_start = COL_CORE_COUNT;
    bind->format_col_start = COL_CORE_COUNT;
    bind->vep_field_indices = NULL;
    bind->vep_schema = NULL;
    bind->vep_transcript_mode = VEP_TRANSCRIPT_FIRST;

    // Copy sample names
    if (bind->n_samples > 0) {
        bind->sample_names = (char**)duckdb_malloc(bind->n_samples * sizeof(char*));
        for (int i = 0; i < bind->n_samples; i++) {
            bind->sample_names[i] = strdup_duckdb(hdr->samples[i]);
        }
    }

    // Create logical types for schema
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type double_type = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type varchar_list_type = duckdb_create_list_type(varchar_type);

    int col_idx = 0;

    // -------------------------------------------------------------------------
    // Core VCF columns (matching nanoarrow schema)
    // -------------------------------------------------------------------------

    // CHROM - VARCHAR (not null)
    duckdb_bind_add_result_column(info, "CHROM", varchar_type);
    col_idx++;

    // POS - BIGINT (1-based position)
    duckdb_bind_add_result_column(info, "POS", bigint_type);
    col_idx++;

    // ID - VARCHAR (nullable)
    duckdb_bind_add_result_column(info, "ID", varchar_type);
    col_idx++;

    // REF - VARCHAR
    duckdb_bind_add_result_column(info, "REF", varchar_type);
    col_idx++;

    // ALT - LIST(VARCHAR) - list of alternate alleles
    duckdb_bind_add_result_column(info, "ALT", varchar_list_type);
    col_idx++;

    // QUAL - DOUBLE (nullable, matching nanoarrow FLOAT64)
    duckdb_bind_add_result_column(info, "QUAL", double_type);
    col_idx++;

    // FILTER - LIST(VARCHAR) - list of filter names
    duckdb_bind_add_result_column(info, "FILTER", varchar_list_type);
    col_idx++;

    // -------------------------------------------------------------------------
    // VEP/CSQ/BCSQ/ANN fields (auto-detected)
    // -------------------------------------------------------------------------
    if (bind->include_vep) {
        bind->vep_schema = vep_schema_parse(hdr, NULL, bind->additional_csq_column_types);
    }
    if (bind->vep_schema) {
        int *vep_filter_matched = NULL;
        if (n_vep_field_filter > 0) {
            vep_filter_matched = (int *)duckdb_malloc((size_t)n_vep_field_filter * sizeof(int));
            if (!vep_filter_matched) {
                duckdb_bind_set_error(info, "read_bcf_v2: out of memory validating vep_fields");
                free_string_list_duckdb(info_field_filter, n_info_field_filter);
                free_string_list_duckdb(format_field_filter, n_format_field_filter);
                free_string_list_duckdb(vep_field_filter, n_vep_field_filter);
                duckdb_destroy_logical_type(&varchar_type);
                duckdb_destroy_logical_type(&bigint_type);
                duckdb_destroy_logical_type(&double_type);
                duckdb_destroy_logical_type(&varchar_list_type);
                bcf_hdr_destroy(hdr);
                hts_close(fp);
                destroy_bind_data(bind);
                return;
            }
            memset(vep_filter_matched, 0, (size_t)n_vep_field_filter * sizeof(int));
        }

        int selected_vep_fields = 0;
        for (int v = 0; v < bind->vep_schema->n_fields; v++) {
            const vep_field_t* field = vep_schema_get_field(bind->vep_schema, v);
            if (field && field_filter_matches(vep_field_filter, n_vep_field_filter, field->name, "VEP_", vep_filter_matched)) {
                selected_vep_fields++;
            }
        }

        if (n_vep_field_filter > 0) {
            char filter_err[256];
            if (!field_filter_check_all_matched(vep_field_filter, n_vep_field_filter, vep_filter_matched,
                                                "vep_fields", filter_err, sizeof(filter_err))) {
                duckdb_bind_set_error(info, filter_err);
                duckdb_free(vep_filter_matched);
                free_string_list_duckdb(info_field_filter, n_info_field_filter);
                free_string_list_duckdb(format_field_filter, n_format_field_filter);
                free_string_list_duckdb(vep_field_filter, n_vep_field_filter);
                duckdb_destroy_logical_type(&varchar_type);
                duckdb_destroy_logical_type(&bigint_type);
                duckdb_destroy_logical_type(&double_type);
                duckdb_destroy_logical_type(&varchar_list_type);
                bcf_hdr_destroy(hdr);
                hts_close(fp);
                destroy_bind_data(bind);
                return;
            }
        }

        bind->n_vep_fields = selected_vep_fields;
        bind->vep_col_start = col_idx;
        if (selected_vep_fields > 0) {
            bind->vep_field_indices = (int *)duckdb_malloc((size_t)selected_vep_fields * sizeof(int));
            if (!bind->vep_field_indices) {
                duckdb_bind_set_error(info, "read_bcf_v2: out of memory allocating VEP field mapping");
                if (vep_filter_matched) duckdb_free(vep_filter_matched);
                free_string_list_duckdb(info_field_filter, n_info_field_filter);
                free_string_list_duckdb(format_field_filter, n_format_field_filter);
                free_string_list_duckdb(vep_field_filter, n_vep_field_filter);
                duckdb_destroy_logical_type(&varchar_type);
                duckdb_destroy_logical_type(&bigint_type);
                duckdb_destroy_logical_type(&double_type);
                duckdb_destroy_logical_type(&varchar_list_type);
                bcf_hdr_destroy(hdr);
                hts_close(fp);
                destroy_bind_data(bind);
                return;
            }
        }

        int out_v = 0;
        for (int v = 0; v < bind->vep_schema->n_fields; v++) {
            const vep_field_t* field = vep_schema_get_field(bind->vep_schema, v);
            if (!field || !field_filter_matches(vep_field_filter, n_vep_field_filter, field->name, "VEP_", NULL)) {
                continue;
            }
            bind->vep_field_indices[out_v] = v;

            char col_name[256];
            snprintf(col_name, sizeof(col_name), "VEP_%s", field->name);

            // Expose all transcripts as list columns for full preservation
            duckdb_logical_type field_type = create_vep_field_type(field->type, 1);
            duckdb_bind_add_result_column(info, col_name, field_type);
            duckdb_destroy_logical_type(&field_type);

            col_idx++;
            out_v++;
        }
        if (vep_filter_matched) duckdb_free(vep_filter_matched);

        bind->vep_transcript_mode = VEP_TRANSCRIPT_ALL;
    } else if (n_vep_field_filter > 0) {
        duckdb_bind_set_error(info, "read_bcf_v2: vep_fields specified but no VEP/CSQ/ANN/BCSQ schema was found");
        free_string_list_duckdb(info_field_filter, n_info_field_filter);
        free_string_list_duckdb(format_field_filter, n_format_field_filter);
        free_string_list_duckdb(vep_field_filter, n_vep_field_filter);
        duckdb_destroy_logical_type(&varchar_type);
        duckdb_destroy_logical_type(&bigint_type);
        duckdb_destroy_logical_type(&double_type);
        duckdb_destroy_logical_type(&varchar_list_type);
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        destroy_bind_data(bind);
        return;
    }

    // -------------------------------------------------------------------------
    // INFO fields (with type validation)
    // -------------------------------------------------------------------------
    bind->info_col_start = col_idx;

    // Count INFO fields
    bind->n_info_fields = 0;
    int *info_filter_matched = NULL;
    if (n_info_field_filter > 0) {
        info_filter_matched = (int *)duckdb_malloc((size_t)n_info_field_filter * sizeof(int));
        if (!info_filter_matched) {
            duckdb_bind_set_error(info, "read_bcf_v2: out of memory validating info_fields");
            free_string_list_duckdb(info_field_filter, n_info_field_filter);
            free_string_list_duckdb(format_field_filter, n_format_field_filter);
            free_string_list_duckdb(vep_field_filter, n_vep_field_filter);
            duckdb_destroy_logical_type(&varchar_type);
            duckdb_destroy_logical_type(&bigint_type);
            duckdb_destroy_logical_type(&double_type);
            duckdb_destroy_logical_type(&varchar_list_type);
            bcf_hdr_destroy(hdr);
            hts_close(fp);
            destroy_bind_data(bind);
            return;
        }
        memset(info_filter_matched, 0, (size_t)n_info_field_filter * sizeof(int));
    }
    if (bind->include_info) {
        for (int i = 0; i < hdr->n[BCF_DT_ID]; i++) {
            if (hdr->id[BCF_DT_ID][i].val &&
                hdr->id[BCF_DT_ID][i].val->hrec[BCF_HL_INFO]) {
                const char* field_name = hdr->id[BCF_DT_ID][i].key;
                if (field_filter_matches(info_field_filter, n_info_field_filter, field_name, "INFO_", info_filter_matched)) {
                    bind->n_info_fields++;
                }
            }
        }
    }
    if (n_info_field_filter > 0) {
        char filter_err[256];
        if (!field_filter_check_all_matched(info_field_filter, n_info_field_filter, info_filter_matched,
                                            "info_fields", filter_err, sizeof(filter_err))) {
            duckdb_bind_set_error(info, filter_err);
            duckdb_free(info_filter_matched);
            free_string_list_duckdb(info_field_filter, n_info_field_filter);
            free_string_list_duckdb(format_field_filter, n_format_field_filter);
            free_string_list_duckdb(vep_field_filter, n_vep_field_filter);
            duckdb_destroy_logical_type(&varchar_type);
            duckdb_destroy_logical_type(&bigint_type);
            duckdb_destroy_logical_type(&double_type);
            duckdb_destroy_logical_type(&varchar_list_type);
            bcf_hdr_destroy(hdr);
            hts_close(fp);
            destroy_bind_data(bind);
            return;
        }
        duckdb_free(info_filter_matched);
        info_filter_matched = NULL;
    }

    if (bind->n_info_fields > 0) {
        bind->info_fields = (field_meta_t*)duckdb_malloc(bind->n_info_fields * sizeof(field_meta_t));
        memset(bind->info_fields, 0, bind->n_info_fields * sizeof(field_meta_t));

        int info_idx = 0;
        for (int i = 0; i < hdr->n[BCF_DT_ID] && info_idx < bind->n_info_fields; i++) {
            if (hdr->id[BCF_DT_ID][i].val &&
                hdr->id[BCF_DT_ID][i].val->hrec[BCF_HL_INFO]) {
                const char* field_name = hdr->id[BCF_DT_ID][i].key;
                if (!field_filter_matches(info_field_filter, n_info_field_filter, field_name, "INFO_", NULL)) {
                    continue;
                }
                int header_type = bcf_hdr_id2type(hdr, BCF_HL_INFO, i);
                int header_vl_type = bcf_hdr_id2length(hdr, BCF_HL_INFO, i);
                int header_count = bcf_hdr_id2number(hdr, BCF_HL_INFO, i);

                // Validate against VCF spec (emits warnings)
                int corrected_type;
                int corrected_count;
                int corrected_vl_type = vcf_validate_info_field(field_name, header_vl_type,
                                                                 header_count, header_type,
                                                                 &corrected_type, &corrected_count);

                field_meta_t* field = &bind->info_fields[info_idx];
                field->name = strdup_duckdb(field_name);
                field->header_id = i;
                field->header_type = header_type;
                field->schema_type = header_type;  // Use header type for data
                field->vl_type = corrected_vl_type;
                field->fixed_count = corrected_vl_type == BCF_VL_FIXED ? corrected_count : 1;
                field->is_list = vcf_is_list_type(corrected_vl_type, field->fixed_count);
                field->duckdb_col_idx = col_idx;

                // Create column name: INFO_<fieldname>
                char col_name[256];
                snprintf(col_name, sizeof(col_name), "INFO_%s", field_name);

                // Create DuckDB type
                duckdb_logical_type field_type = create_bcf_field_type(header_type, field->is_list);
                duckdb_bind_add_result_column(info, col_name, field_type);
                duckdb_destroy_logical_type(&field_type);

                col_idx++;
                info_idx++;
            }
        }
    }

    // -------------------------------------------------------------------------
    // FORMAT fields per sample (with type validation)
    // -------------------------------------------------------------------------

    bind->format_col_start = col_idx;

    if (bind->include_format && bind->n_samples > 0) {
        // Count FORMAT fields
        bind->n_format_fields = 0;
        int *format_filter_matched = NULL;
        if (n_format_field_filter > 0) {
            format_filter_matched = (int *)duckdb_malloc((size_t)n_format_field_filter * sizeof(int));
            if (!format_filter_matched) {
                duckdb_bind_set_error(info, "read_bcf_v2: out of memory validating format_fields");
                free_string_list_duckdb(info_field_filter, n_info_field_filter);
                free_string_list_duckdb(format_field_filter, n_format_field_filter);
                free_string_list_duckdb(vep_field_filter, n_vep_field_filter);
                duckdb_destroy_logical_type(&varchar_type);
                duckdb_destroy_logical_type(&bigint_type);
                duckdb_destroy_logical_type(&double_type);
                duckdb_destroy_logical_type(&varchar_list_type);
                bcf_hdr_destroy(hdr);
                hts_close(fp);
                destroy_bind_data(bind);
                return;
            }
            memset(format_filter_matched, 0, (size_t)n_format_field_filter * sizeof(int));
        }
        int header_format_count = 0;
        for (int i = 0; i < hdr->n[BCF_DT_ID]; i++) {
            if (hdr->id[BCF_DT_ID][i].val &&
                hdr->id[BCF_DT_ID][i].val->hrec[BCF_HL_FMT]) {
                header_format_count++;
                const char* field_name = hdr->id[BCF_DT_ID][i].key;
                if (field_filter_matches(format_field_filter, n_format_field_filter, field_name, "FORMAT_", format_filter_matched)) {
                    bind->n_format_fields++;
                }
            }
        }
        if (bind->n_format_fields == 0 && header_format_count == 0 &&
            field_filter_matches(format_field_filter, n_format_field_filter, "GT", "FORMAT_", format_filter_matched)) {
            bind->n_format_fields = 1;
        }
        if (n_format_field_filter > 0) {
            char filter_err[256];
            if (!field_filter_check_all_matched(format_field_filter, n_format_field_filter, format_filter_matched,
                                                "format_fields", filter_err, sizeof(filter_err))) {
                duckdb_bind_set_error(info, filter_err);
                duckdb_free(format_filter_matched);
                free_string_list_duckdb(info_field_filter, n_info_field_filter);
                free_string_list_duckdb(format_field_filter, n_format_field_filter);
                free_string_list_duckdb(vep_field_filter, n_vep_field_filter);
                duckdb_destroy_logical_type(&varchar_type);
                duckdb_destroy_logical_type(&bigint_type);
                duckdb_destroy_logical_type(&double_type);
                duckdb_destroy_logical_type(&varchar_list_type);
                bcf_hdr_destroy(hdr);
                hts_close(fp);
                destroy_bind_data(bind);
                return;
            }
            duckdb_free(format_filter_matched);
            format_filter_matched = NULL;
        }

        if (bind->n_format_fields == 1 && header_format_count == 0 &&
            field_filter_matches(format_field_filter, n_format_field_filter, "GT", "FORMAT_", NULL)) {
            // Add GT as default for old header-light fixtures.
            bind->format_fields = (field_meta_t*)duckdb_malloc(sizeof(field_meta_t));
            memset(bind->format_fields, 0, sizeof(field_meta_t));
            bind->format_fields[0].name = strdup_duckdb("GT");
            bind->format_fields[0].header_type = BCF_HT_STR;
            bind->format_fields[0].schema_type = BCF_HT_STR;
            bind->format_fields[0].vl_type = BCF_VL_FIXED;
            bind->format_fields[0].fixed_count = 1;
            bind->format_fields[0].is_list = 0;
        }
        if (!bind->format_fields && bind->n_format_fields > 0) {
            bind->format_fields = (field_meta_t*)duckdb_malloc(bind->n_format_fields * sizeof(field_meta_t));
            memset(bind->format_fields, 0, bind->n_format_fields * sizeof(field_meta_t));

            int fmt_idx = 0;
            for (int i = 0; i < hdr->n[BCF_DT_ID] && fmt_idx < bind->n_format_fields; i++) {
                if (hdr->id[BCF_DT_ID][i].val &&
                    hdr->id[BCF_DT_ID][i].val->hrec[BCF_HL_FMT]) {
                    const char* field_name = hdr->id[BCF_DT_ID][i].key;
                    if (!field_filter_matches(format_field_filter, n_format_field_filter, field_name, "FORMAT_", NULL)) {
                        continue;
                    }
                    int header_type = bcf_hdr_id2type(hdr, BCF_HL_FMT, i);
                    int header_vl_type = bcf_hdr_id2length(hdr, BCF_HL_FMT, i);
                    int header_count = bcf_hdr_id2number(hdr, BCF_HL_FMT, i);

                    // Validate against VCF spec (emits warnings, only once)
                    int corrected_type;
                    int corrected_count;
                    int corrected_vl_type = vcf_validate_format_field(field_name, header_vl_type,
                                                                       header_count, header_type,
                                                                       &corrected_type, &corrected_count);

                    field_meta_t* field = &bind->format_fields[fmt_idx];
                    field->name = strdup_duckdb(field_name);
                    field->header_id = i;
                    field->header_type = header_type;
                    field->schema_type = header_type;
                    field->vl_type = corrected_vl_type;
                    field->fixed_count = corrected_vl_type == BCF_VL_FIXED ? corrected_count : 1;
                    field->is_list = vcf_is_list_type(corrected_vl_type, field->fixed_count);

                    fmt_idx++;
                }
            }
        }

        // Add FORMAT columns for each sample (or single set for tidy format)
        if (bind->tidy_format) {
            // Tidy format: Add SAMPLE_ID column, then FORMAT_<field> (no sample suffix)
            bind->sample_id_col_idx = col_idx;
            duckdb_bind_add_result_column(info, "SAMPLE_ID", varchar_type);
            col_idx++;

            // Update format_col_start to be after SAMPLE_ID
            bind->format_col_start = col_idx;

            // Add FORMAT columns once (no sample suffix)
            for (int f = 0; f < bind->n_format_fields; f++) {
                field_meta_t* field = &bind->format_fields[f];

                char col_name[256];
                snprintf(col_name, sizeof(col_name), "FORMAT_%s", field->name);

                duckdb_logical_type field_type = create_bcf_field_type(field->header_type, field->is_list);
                duckdb_bind_add_result_column(info, col_name, field_type);
                duckdb_destroy_logical_type(&field_type);

                col_idx++;
            }
        } else {
            // Wide format: Add FORMAT columns for each sample
            for (int s = 0; s < bind->n_samples; s++) {
                for (int f = 0; f < bind->n_format_fields; f++) {
                    field_meta_t* field = &bind->format_fields[f];

                    // Column name: FORMAT_<fieldname>_<samplename>
                    char col_name[512];
                    snprintf(col_name, sizeof(col_name), "FORMAT_%s_%s",
                             field->name, bind->sample_names[s]);

                    duckdb_logical_type field_type = create_bcf_field_type(field->header_type, field->is_list);
                    duckdb_bind_add_result_column(info, col_name, field_type);
                    duckdb_destroy_logical_type(&field_type);

                    col_idx++;
                }
            }
        }
    }

    bind->total_columns = col_idx;

    // -------------------------------------------------------------------------
    // Check for index and extract contig names for parallel scanning
    // -------------------------------------------------------------------------

    bind->has_index = 0;
    bind->n_contigs = 0;
    bind->contig_names = NULL;

    // Only set up parallel scan if no user-specified region
    if (bind->n_regions == 0 && !bind->scan_sequential) {
        // Try to load index using *_load3 with minimal flags to avoid network timeouts
        // Only use HTS_IDX_SAVE_REMOTE for actual remote protocols
        int is_remote = (strncmp(file_path, "http://", 7) == 0 ||
                         strncmp(file_path, "https://", 8) == 0 ||
                         strncmp(file_path, "ftp://", 6) == 0 ||
                         strncmp(file_path, "s3://", 5) == 0 ||
                         strncmp(file_path, "gs://", 5) == 0);

        hts_idx_t* idx = NULL;
        tbx_t* tbx = NULL;
        enum htsExactFormat fmt = hts_get_format(fp)->format;
        int flags = HTS_IDX_SILENT_FAIL;
        if (is_remote) {
            flags |= HTS_IDX_SAVE_REMOTE;
        }

        if (fmt == bcf) {
            idx = bcf_index_load3(file_path, index_path, flags);
        } else {
            tbx = tbx_index_load3(file_path, index_path, flags);
            if (!tbx) {
                idx = bcf_index_load3(file_path, index_path, flags);
            }
        }

        if (idx || tbx) {
            bind->has_index = 1;
            int tidy_rows = bind->tidy_format && bind->include_format && bind->n_samples > 0;
            uint64_t row_multiplier = tidy_rows ? (uint64_t)bind->n_samples : 1;
            bind->index_row_count_valid = bcf_try_get_index_row_count(idx ? idx : tbx->idx,
                                                                      row_multiplier,
                                                                      &bind->index_row_count);

            // Get contig names from header for legacy automatic parallel contig scans.
            // read_bcf_v2 streams full-file scans by default, so it only needs index
            // metadata for count-only shortcuts and avoids creating seek-heavy iterators.
            if (bind->reader_version < 2) {
                int n_seqs = hdr->n[BCF_DT_CTG];
                if (n_seqs > 0) {
                    bind->n_contigs = n_seqs;
                    bind->contig_names = (char**)duckdb_malloc(n_seqs * sizeof(char*));

                    for (int i = 0; i < n_seqs; i++) {
                        bind->contig_names[i] = strdup_duckdb(hdr->id[BCF_DT_CTG][i].key);
                    }
                }
            }

            if (idx) hts_idx_destroy(idx);
            if (tbx) tbx_destroy(tbx);
        }
    }

    // Cleanup
    free_string_list_duckdb(info_field_filter, n_info_field_filter);
    free_string_list_duckdb(format_field_filter, n_format_field_filter);
    free_string_list_duckdb(vep_field_filter, n_vep_field_filter);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&double_type);
    duckdb_destroy_logical_type(&varchar_list_type);

    bcf_hdr_destroy(hdr);
    hts_close(fp);

    duckdb_bind_set_bind_data(info, bind, destroy_bind_data);
}

static void bcf_read_bind(duckdb_bind_info info) {
    bcf_read_bind_common(info, 1);
}

static void bcf_read_v2_bind(duckdb_bind_info info) {
    bcf_read_bind_common(info, 2);
}

static int bcf_v2_vcf_max_unpack(int unpack_mask, idx_t column_count) {
    if (column_count == 0) {
        return BCF_UN_STR;
    }
    if (unpack_mask & BCF_UN_FMT) {
        return 0;
    }
    if (unpack_mask & BCF_UN_INFO) {
        return BCF_UN_INFO;
    }
    if (unpack_mask & BCF_UN_FLT) {
        return BCF_UN_FLT;
    }
    return BCF_UN_STR;
}

// =============================================================================
// Global Init Function - Set up parallel scanning
// =============================================================================

static void bcf_read_global_init(duckdb_init_info info) {
    bcf_bind_data_t* bind = (bcf_bind_data_t*)duckdb_init_get_bind_data(info);
    idx_t column_count = duckdb_init_get_column_count(info);

    bcf_global_init_data_t* global = (bcf_global_init_data_t*)duckdb_malloc(sizeof(bcf_global_init_data_t));
    memset(global, 0, sizeof(bcf_global_init_data_t));

    global->current_contig = 0;
    global->has_region = (bind->n_regions > 0);

    // Enable parallel scan if:
    // 1. Index exists
    // 2. Multiple contigs available
    // 3. No user-specified region (region queries are already filtered)
    if (column_count == 0 && !bind->scan_sequential && bind->index_row_count_valid && !global->has_region) {
        global->n_contigs = 0;
        global->contig_names = NULL;
        duckdb_init_set_max_threads(info, 1);
    } else if (!bind->scan_sequential && bind->reader_version < 2 && bind->has_index && bind->n_contigs > 1 && !global->has_region) {
        global->n_contigs = bind->n_contigs;
        global->contig_names = bind->contig_names;  // Reference only

        // Cap threads at number of contigs or reasonable max
        idx_t max_threads = bind->n_contigs;
        if (max_threads > 16) max_threads = 16;
        duckdb_init_set_max_threads(info, max_threads);
    } else {
        // Single-threaded scan (single contig or no index)
        global->n_contigs = 0;
        global->contig_names = NULL;
        duckdb_init_set_max_threads(info, 1);
    }

    duckdb_init_set_init_data(info, global, destroy_global_init_data);
}

// =============================================================================
// Local Init Function - Per-thread scanning state
// =============================================================================

static void bcf_read_local_init(duckdb_init_info info) {
    bcf_bind_data_t* bind = (bcf_bind_data_t*)duckdb_init_get_bind_data(info);
    const char *reader_name = bind && bind->reader_version >= 2 ? "read_bcf_v2" : "read_bcf";

    bcf_init_data_t* local = (bcf_init_data_t*)duckdb_malloc(sizeof(bcf_init_data_t));
    memset(local, 0, sizeof(bcf_init_data_t));

    local->column_count = duckdb_init_get_column_count(info);
    if (local->column_count > 0) {
        local->column_ids = (idx_t*)duckdb_malloc(sizeof(idx_t) * local->column_count);
        for (idx_t i = 0; i < local->column_count; i++) {
            local->column_ids[i] = duckdb_init_get_column_index(info, i);
        }
        local->unpack_mask = bcf_projection_unpack_mask(bind, local->column_ids, local->column_count);
        local->need_vep = bcf_projection_needs_vep(bind, local->column_ids, local->column_count);
        char projection_err[256];
        if (!bcf_projected_columns_init(local, bind, projection_err, sizeof(projection_err))) {
            duckdb_init_set_error(info, projection_err);
            destroy_init_data(local);
            return;
        }
    }

    if (local->column_count == 0 && !bind->scan_sequential && bind->n_regions == 0 && bind->index_row_count_valid) {
        local->count_only = 1;
        local->count_remaining = bind->index_row_count;
        local->done = (local->count_remaining == 0);
        duckdb_init_set_init_data(info, local, destroy_init_data);
        return;
    }

    if (bind->reader_version >= 2 && !bcf_decode_cache_init(local, bind)) {
        duckdb_init_set_error(info, "read_bcf_v2: out of memory allocating decode cache");
        destroy_init_data(local);
        return;
    }

    // Check if we're in parallel mode based on bind data
    int is_parallel = (!bind->scan_sequential && bind->reader_version < 2 && bind->has_index && bind->n_contigs > 1 && bind->n_regions == 0);

    // Initialize parallel scan state
    local->is_parallel = is_parallel;
    local->assigned_contig = -1;
    local->contig_name = NULL;
    local->needs_next_contig = is_parallel;  // Start by requesting first contig

    // Open file (each thread gets its own file handle)
    local->fp = hts_open(bind->file_path, "r");
    if (!local->fp) {
        char err[512];
        snprintf(err, sizeof(err), "%s: failed to open BCF/VCF file: %s", reader_name, bind->file_path);
        duckdb_init_set_error(info, err);
        destroy_init_data(local);
        return;
    }
    duckhts_apply_remote_hts_tuning(
        local->fp,
        bind->file_path,
        (is_parallel || bind->n_regions > 0)
            ? DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION
            : DUCKHTS_HTS_IO_PROFILE_STREAMING
    );

    if (bind->decompression_threads > 0 &&
        hts_get_format(local->fp)->compression != no_compression &&
        hts_set_threads(local->fp, bind->decompression_threads) < 0) {
        char err[160];
        snprintf(err, sizeof(err), "%s: failed to configure BCF/VCF decompression threads", reader_name);
        duckdb_init_set_error(info, err);
        destroy_init_data(local);
        return;
    }

    // Read header
    local->hdr = bcf_hdr_read(local->fp);
    if (!local->hdr) {
        char err[128];
        snprintf(err, sizeof(err), "%s: failed to read BCF/VCF header", reader_name);
        duckdb_init_set_error(info, err);
        destroy_init_data(local);
        return;
    }
    if (bind->reader_version >= 2 && bind->sample_filter_spec) {
        int sample_ret = bcf_hdr_set_samples(local->hdr, bind->sample_filter_spec, bind->sample_filter_is_file);
        if (sample_ret < 0 || (sample_ret > 0 && !bind->force_samples)) {
            duckdb_init_set_error(info, "read_bcf_v2: failed to apply sample selector during scan initialization");
            destroy_init_data(local);
            return;
        }
    }

    // Allocate record
    local->rec = bcf_init();
    if (!local->rec) {
        char err[128];
        snprintf(err, sizeof(err), "%s: failed to allocate BCF/VCF record", reader_name);
        duckdb_init_set_error(info, err);
        destroy_init_data(local);
        return;
    }
    if (bind->reader_version >= 2) {
        local->rec->max_unpack = bcf_v2_vcf_max_unpack(local->unpack_mask, local->column_count);
    }

    // Load index for parallel scanning or region queries
    // Use *_load3 with HTS_IDX_SAVE_REMOTE for remote file support
    if (is_parallel || bind->n_regions > 0) {
        enum htsExactFormat fmt = hts_get_format(local->fp)->format;

        if (fmt == bcf) {
            local->idx = bcf_index_load3(bind->file_path, bind->index_path, HTS_IDX_SAVE_REMOTE | HTS_IDX_SILENT_FAIL);
        } else {
            local->tbx = tbx_index_load3(bind->file_path, bind->index_path, HTS_IDX_SAVE_REMOTE | HTS_IDX_SILENT_FAIL);
            if (!local->tbx) {
                local->idx = bcf_index_load3(bind->file_path, bind->index_path, HTS_IDX_SAVE_REMOTE | HTS_IDX_SILENT_FAIL);
            }
        }
    }

    // Set up region query if user specified a region (non-parallel case)
    if (!is_parallel && bind->n_regions > 0) {
        // First check if we have an index
        if (!local->idx && !local->tbx) {
            char err[512];
            snprintf(err, sizeof(err),
                     "%s: region query requires an index file (.tbi or .csi). Region: %s",
                     reader_name, bind->region);
            duckdb_init_set_error(info, err);
            destroy_init_data(local);
            return;
        }

        local->itr = bcf_open_region_iterator(local->idx, local->hdr, local->tbx,
                                               bind->regions, bind->n_regions);

        if (!local->itr) {
            local->done = 1;
            duckdb_init_set_init_data(info, local, destroy_init_data);
            return;
        }
    }

    local->current_row = 0;
    local->done = 0;

// Initialize debug/progress tracking
    local->total_records_processed = 0;
    memset(&local->batch_start_time, 0, sizeof(local->batch_start_time));
    memset(&local->last_progress_time, 0, sizeof(local->last_progress_time));
    local->timing_initialized = 0;

    // Store as local init data
    duckdb_init_set_init_data(info, local, destroy_init_data);
}

// =============================================================================
// Helper: Set validity bit
// =============================================================================

static inline void set_validity_bit(uint64_t* validity, idx_t row, int is_valid) {
    if (!validity) return;
    idx_t entry_idx = row / 64;
    idx_t bit_idx = row % 64;
    if (is_valid) {
        validity[entry_idx] |= ((uint64_t)1 << bit_idx);
    } else {
        validity[entry_idx] &= ~((uint64_t)1 << bit_idx);
    }
}

// Single-pass comma-separated string list processing
static void process_comma_separated_list(duckdb_vector vec, idx_t row, const char* value) {
    if (!value || strcmp(value, ".") == 0) {
        // NULL value - empty list
        duckdb_vector_ensure_validity_writable(vec);
        uint64_t* validity = duckdb_vector_get_validity(vec);
        set_validity_bit(validity, row, 0);
        duckdb_list_entry entry = {duckdb_list_vector_get_size(vec), 0};
        duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
        list_data[row] = entry;
        return;
    }

    duckdb_list_entry entry;
    entry.offset = duckdb_list_vector_get_size(vec);
    entry.length = 0;

    // Single-pass: count tokens and assign in one go
    const char* p = value;
    const char* token_start = p;
    int token_count = 0;

    // First pass: count tokens
    while (*p) {
        if (*p == ',') {
            token_count++;
            token_start = p + 1;
        }
        p++;
    }
    if (p > token_start) token_count++;  // Last token

    entry.length = token_count;

    // Reserve and fill
    if (entry.length > 0) {
        duckdb_list_vector_reserve(vec, entry.offset + entry.length);
        duckdb_list_vector_set_size(vec, entry.offset + entry.length);
        duckdb_vector child_vec = duckdb_list_vector_get_child(vec);

        // Second pass: assign tokens
        p = value;
        token_start = p;
        int write_idx = 0;

        while (*p) {
            if (*p == ',') {
                // Assign current token
                duckdb_vector_assign_string_element_len(child_vec, entry.offset + write_idx,
                                                     token_start, p - token_start);
                write_idx++;
                token_start = p + 1;
            }
            p++;
        }

        // Last token
        if (p > token_start) {
            duckdb_vector_assign_string_element_len(child_vec, entry.offset + write_idx,
                                                 token_start, p - token_start);
        }
    }

    duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
    list_data[row] = entry;
}

#if BCF_READER_ENABLE_PROGRESS
// =============================================================================
// Helper: Print debug progress
// =============================================================================

static void print_progress(bcf_init_data_t* init, const char* context) {
    if (init->total_records_processed % BCF_READER_PROGRESS_INTERVAL == 0) {
        // Calculate records per second if we have timing
        double records_per_sec = 0.0;
        if (init->timing_initialized) {
            struct timespec now;
            clock_gettime(CLOCK_MONOTONIC, &now);

            // Calculate elapsed time since last progress report
            double elapsed = (now.tv_sec - init->last_progress_time.tv_sec) +
                          (now.tv_nsec - init->last_progress_time.tv_nsec) / 1e9;

            // Calculate rate - for very fast processing, use total elapsed time
            if (elapsed > 0.001) {  // Only calculate after 1ms to avoid division by zero
                records_per_sec = BCF_READER_PROGRESS_INTERVAL / elapsed;
            } else {
                // If processing is very fast (<1ms per 100 records), use total time
                double total_elapsed = (now.tv_sec - init->batch_start_time.tv_sec) +
                                     (now.tv_nsec - init->batch_start_time.tv_nsec) / 1e9;
                if (total_elapsed > 0.001) {
                    records_per_sec = init->total_records_processed / total_elapsed;
                }
            }

            // Update last progress time
            init->last_progress_time = now;

            // Debug: print elapsed time for first few intervals
            if (init->total_records_processed <= 300) {
                fprintf(stderr, "[bcf_reader] DEBUG: elapsed=%.6f sec, timing_initialized=%d\n",
                        elapsed, init->timing_initialized);
            }
        }

        fprintf(stderr, "[bcf_reader] %s: Processed %" PRId64 " records (%.0f rec/s)\n",
                context, init->total_records_processed, records_per_sec);
    }
}
#endif

// =============================================================================
// Helper: Claim next contig for parallel scanning
// Returns 1 if a new contig was claimed, 0 if no more contigs
// =============================================================================

static int claim_next_contig(bcf_init_data_t* init, bcf_global_init_data_t* global) {
    if (!init->is_parallel || !global || global->n_contigs == 0) {
        return 0;
    }

    // Destroy old iterator if exists
    if (init->itr) {
        hts_itr_destroy(init->itr);
        init->itr = NULL;
    }

    for (;;) {
        // Atomically claim next contig using fetch-and-add.
        int next = __sync_fetch_and_add(&global->current_contig, 1);
        if (next >= global->n_contigs) {
            return 0;  // No more contigs
        }

        // Set up iterator for this contig.
        const char* contig = global->contig_names[next];
        init->assigned_contig = next;
        init->contig_name = contig;

        if (init->idx) {
            init->itr = bcf_itr_querys(init->idx, init->hdr, contig);
        } else if (init->tbx) {
            init->itr = tbx_itr_querys(init->tbx, contig);
        }

        if (init->itr) {
            init->needs_next_contig = 0;
            return 1;
        }

        // Empty contig: keep claiming until we find a live iterator.
    }
}

static int bcf_decode_format_projected_field(bcf_init_data_t *init,
                                             const bcf_bind_data_t *bind,
                                             const bcf_projected_col_t *proj_col,
                                             int *fmt_loaded, int *fmt_ret,
                                             int *fmt_n_values,
                                             int32_t **fmt_i32, float **fmt_f32,
                                             char ***fmt_str,
                                             int *gt_loaded, int *gt_ret,
                                             int *gt_n_values, int32_t **gt_arr,
                                             char *err, size_t err_size) {
    if (!init || !bind || !proj_col || proj_col->field_idx < 0 || !proj_col->field) {
        return 1;
    }

    int field_idx = proj_col->field_idx;
    const char *tag = proj_col->field->name;
    const char *reader_name = bcf_reader_name(bind);
    bcf_decode_error_policy_t policy = bind->decode_error_policy;

    switch (proj_col->kind) {
    case BCF_OUT_FORMAT_INT:
        if (fmt_loaded && !fmt_loaded[field_idx]) {
            if (!bcf_record_has_format_field(init->rec, proj_col->field)) {
                fmt_ret[field_idx] = -3;
            } else if (bcf_reader_input_is_bcf(init->fp) &&
                       !bcf_preflight_format_encoded_type(init->hdr, init->rec,
                                                          proj_col->field->header_id,
                                                          tag, proj_col->kind,
                                                          proj_col->field->header_type,
                                                          reader_name, err, err_size)) {
                if (!bcf_handle_decode_diagnostic(policy, err)) {
                    return 0;
                }
                fmt_ret[field_idx] = -3;
            } else {
                fmt_ret[field_idx] = bcf_get_format_int32(init->hdr, init->rec, tag,
                                                           &fmt_i32[field_idx],
                                                           &fmt_n_values[field_idx]);
            }
            fmt_loaded[field_idx] = 1;
        }
        if (!bcf_check_decode_ret(reader_name, "FORMAT", tag, init->hdr, init->rec,
                                  fmt_ret[field_idx], policy, err, err_size)) {
            return 0;
        }
        return bcf_handle_format_width(reader_name, tag, init->hdr, init->rec,
                                       &fmt_ret[field_idx], bind->n_samples,
                                       policy, err, err_size);
    case BCF_OUT_FORMAT_FLOAT:
        if (fmt_loaded && !fmt_loaded[field_idx]) {
            if (!bcf_record_has_format_field(init->rec, proj_col->field)) {
                fmt_ret[field_idx] = -3;
            } else if (bcf_reader_input_is_bcf(init->fp) &&
                       !bcf_preflight_format_encoded_type(init->hdr, init->rec,
                                                          proj_col->field->header_id,
                                                          tag, proj_col->kind,
                                                          proj_col->field->header_type,
                                                          reader_name, err, err_size)) {
                if (!bcf_handle_decode_diagnostic(policy, err)) {
                    return 0;
                }
                fmt_ret[field_idx] = -3;
            } else {
                fmt_ret[field_idx] = bcf_get_format_float(init->hdr, init->rec, tag,
                                                           &fmt_f32[field_idx],
                                                           &fmt_n_values[field_idx]);
            }
            fmt_loaded[field_idx] = 1;
        }
        if (!bcf_check_decode_ret(reader_name, "FORMAT", tag, init->hdr, init->rec,
                                  fmt_ret[field_idx], policy, err, err_size)) {
            return 0;
        }
        return bcf_handle_format_width(reader_name, tag, init->hdr, init->rec,
                                       &fmt_ret[field_idx], bind->n_samples,
                                       policy, err, err_size);
    case BCF_OUT_FORMAT_GT:
        if (gt_loaded && !*gt_loaded) {
            if (!bcf_record_has_format_field(init->rec, proj_col->field)) {
                *gt_ret = -3;
            } else if (bcf_reader_input_is_bcf(init->fp) &&
                       !bcf_preflight_format_encoded_type(init->hdr, init->rec,
                                                          proj_col->field->header_id,
                                                          tag, proj_col->kind,
                                                          proj_col->field->header_type,
                                                          reader_name, err, err_size)) {
                if (!bcf_handle_decode_diagnostic(policy, err)) {
                    return 0;
                }
                *gt_ret = -3;
            } else {
                *gt_ret = bcf_get_genotypes(init->hdr, init->rec, gt_arr, gt_n_values);
            }
            *gt_loaded = 1;
        }
        if (!bcf_check_decode_ret(reader_name, "FORMAT", tag, init->hdr, init->rec,
                                  *gt_ret, policy, err, err_size)) {
            return 0;
        }
        return bcf_handle_format_width(reader_name, tag, init->hdr, init->rec,
                                       gt_ret, bind->n_samples,
                                       policy, err, err_size);
    case BCF_OUT_FORMAT_STRING:
        if (fmt_loaded && !fmt_loaded[field_idx]) {
            if (!bcf_record_has_format_field(init->rec, proj_col->field)) {
                fmt_ret[field_idx] = -3;
            } else if (bcf_reader_input_is_bcf(init->fp) &&
                       !bcf_preflight_format_encoded_type(init->hdr, init->rec,
                                                          proj_col->field->header_id,
                                                          tag, proj_col->kind,
                                                          proj_col->field->header_type,
                                                          reader_name, err, err_size)) {
                if (!bcf_handle_decode_diagnostic(policy, err)) {
                    return 0;
                }
                fmt_ret[field_idx] = -3;
            } else {
                fmt_ret[field_idx] = bcf_get_format_string(init->hdr, init->rec, tag,
                                                            &fmt_str[field_idx],
                                                            &fmt_n_values[field_idx]);
            }
            fmt_loaded[field_idx] = 1;
        }
        return bcf_check_decode_ret(reader_name, "FORMAT", tag, init->hdr, init->rec,
                                    fmt_ret[field_idx], policy, err, err_size);
    default:
        break;
    }
    return 1;
}

static void bcf_fill_gt_string(duckdb_vector vec, idx_t row_count,
                               bcf_init_data_t *init, const int32_t *sample_gt,
                               int ploidy) {
    char gt_small[4];
    const char *gt_text = NULL;
    size_t gt_len = 0;
    if (format_gt_fast(sample_gt, ploidy, gt_small, &gt_text, &gt_len)) {
        duckdb_vector_assign_string_element_len(vec, row_count, gt_text, gt_len);
        return;
    }

    init->gt_kstr.l = 0;
    if (init->gt_kstr.s) init->gt_kstr.s[0] = '\0';

    for (int p = 0; p < ploidy; p++) {
        if (sample_gt[p] == bcf_int32_vector_end) {
            break;
        }
        if (p > 0) {
            if (kputc(bcf_gt_is_phased(sample_gt[p]) ? '|' : '/', &init->gt_kstr) < 0) {
                break;
            }
        }

        if (bcf_gt_is_missing(sample_gt[p])) {
            if (kputc('.', &init->gt_kstr) < 0) {
                break;
            }
        } else {
            int allele = bcf_gt_allele(sample_gt[p]);
            if (kputw(allele, &init->gt_kstr) < 0) {
                break;
            }
        }
    }

    if (init->gt_kstr.l > 0 && init->gt_kstr.s) {
        duckdb_vector_assign_string_element_len(vec, row_count,
                                                init->gt_kstr.s,
                                                init->gt_kstr.l);
    } else {
        duckdb_vector_ensure_validity_writable(vec);
        uint64_t* validity = duckdb_vector_get_validity(vec);
        set_validity_bit(validity, row_count, 0);
    }
}

static void bcf_fill_format_projected_col(bcf_init_data_t *init,
                                          const bcf_bind_data_t *bind,
                                          const bcf_projected_col_t *proj_col,
                                          duckdb_vector vec,
                                          idx_t row_count,
                                          int current_sample,
                                          int *fmt_ret,
                                          int32_t **fmt_i32,
                                          float **fmt_f32,
                                          char ***fmt_str,
                                          int gt_ret,
                                          int32_t *gt_arr) {
    int field_idx = proj_col->field_idx;
    int sample_idx = (bind->tidy_format && bind->include_format && bind->n_samples > 0)
        ? current_sample
        : proj_col->sample_idx;

    if (sample_idx < 0 || sample_idx >= bind->n_samples || field_idx < 0 || field_idx >= bind->n_format_fields) {
        return;
    }

    switch (proj_col->kind) {
    case BCF_OUT_FORMAT_INT: {
        int ret_fmt = fmt_ret[field_idx];
        int32_t *values = fmt_i32[field_idx];
        if (ret_fmt > 0 && values) {
            int vals_per_sample = ret_fmt / bind->n_samples;
            int32_t *sample_vals = values + sample_idx * vals_per_sample;
            assign_int32_field(vec, row_count, sample_vals, vals_per_sample, proj_col->is_list);
        } else {
            set_field_null(vec, row_count, proj_col->is_list);
        }
        break;
    }
    case BCF_OUT_FORMAT_FLOAT: {
        int ret_fmt = fmt_ret[field_idx];
        float *values = fmt_f32[field_idx];
        if (ret_fmt > 0 && values) {
            int vals_per_sample = ret_fmt / bind->n_samples;
            float *sample_vals = values + sample_idx * vals_per_sample;
            assign_float_field(vec, row_count, sample_vals, vals_per_sample, proj_col->is_list);
        } else {
            set_field_null(vec, row_count, proj_col->is_list);
        }
        break;
    }
    case BCF_OUT_FORMAT_GT:
        if (gt_ret > 0 && gt_arr) {
            int ploidy = gt_ret / bind->n_samples;
            int32_t *sample_gt = gt_arr + sample_idx * ploidy;
            bcf_fill_gt_string(vec, row_count, init, sample_gt, ploidy);
        } else {
            duckdb_vector_ensure_validity_writable(vec);
            uint64_t* validity = duckdb_vector_get_validity(vec);
            set_validity_bit(validity, row_count, 0);
        }
        break;
    case BCF_OUT_FORMAT_STRING: {
        int ret_fmt = fmt_ret[field_idx];
        char **values = fmt_str[field_idx];
        assign_string_field(
            vec,
            row_count,
            (ret_fmt > 0 && values) ? values[sample_idx] : NULL,
            proj_col->is_list
        );
        break;
    }
    default:
        break;
    }
}

// =============================================================================
// Main Scan Function
// =============================================================================

static void bcf_read_function(duckdb_function_info info, duckdb_data_chunk output) {
    bcf_bind_data_t* bind = (bcf_bind_data_t*)duckdb_function_get_bind_data(info);
    bcf_global_init_data_t* global = (bcf_global_init_data_t*)duckdb_function_get_init_data(info);

    // Try to get local init data first (for parallel scans)
    bcf_init_data_t* init = (bcf_init_data_t*)duckdb_function_get_local_init_data(info);
    if (!init) {
        // Fall back to regular init data (shouldn't happen with our setup)
        init = (bcf_init_data_t*)duckdb_function_get_init_data(info);
    }

    if (!init || init->done) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    if (init->count_only) {
        idx_t row_count = (init->count_remaining > (uint64_t)duckdb_vector_size())
            ? duckdb_vector_size()
            : (idx_t)init->count_remaining;
        init->count_remaining -= row_count;
        if (init->count_remaining == 0) {
            init->done = 1;
        }
        duckdb_data_chunk_set_size(output, row_count);
        return;
    }

    // For parallel scans, claim first/next contig if needed
    if (init->needs_next_contig) {
        if (!claim_next_contig(init, global)) {
            // No more contigs to process
            init->done = 1;
            duckdb_data_chunk_set_size(output, 0);
            return;
        }
    }

    idx_t vector_size = duckdb_vector_size();
    idx_t row_count = 0;

    // Cache vector pointers to reduce repeated calls
    duckdb_vector* vectors = NULL;
    if (init->column_count > 0) {
        vectors = (duckdb_vector*)duckdb_malloc(init->column_count * sizeof(duckdb_vector));
        for (idx_t i = 0; i < init->column_count; i++) {
            vectors[i] = duckdb_data_chunk_get_vector(output, i);
        }
    }

    // Tidy format variables
    int tidy_mode = bind->tidy_format && bind->include_format && bind->n_samples > 0;
    int current_sample = 0;  // Which sample we're emitting (only used in tidy mode)

    // Per-record FORMAT cache to avoid repeated bcf_get_format_* calls across
    // projected columns/samples for the same record. read_bcf_v2 keeps these
    // buffers in local init data so tidy records split across DuckDB chunks do
    // not re-decode the same BCF/VCF record on the next function call.
    int use_persistent_cache = bind->reader_version >= 2 && init->cache_persistent;
    int *fmt_loaded = NULL;
    int *fmt_ret = NULL;
    int *fmt_n_values = NULL;
    int32_t **fmt_i32 = NULL;
    float **fmt_f32 = NULL;
    char ***fmt_str = NULL;
    int gt_loaded = 0;
    int gt_ret = 0;
    int gt_n_values = 0;
    int32_t *gt_arr = NULL;
    if (use_persistent_cache) {
        fmt_loaded = init->fmt_loaded;
        fmt_ret = init->fmt_ret;
        fmt_n_values = init->fmt_n_values;
        fmt_i32 = init->fmt_i32;
        fmt_f32 = init->fmt_f32;
        fmt_str = init->fmt_str;
        gt_loaded = init->gt_loaded;
        gt_ret = init->gt_ret;
        gt_n_values = init->gt_n_values;
        gt_arr = init->gt_arr;
    } else if (bind->n_format_fields > 0) {
        size_t nfmt = (size_t)bind->n_format_fields;
        fmt_loaded = (int *)duckdb_malloc(nfmt * sizeof(int));
        fmt_ret = (int *)duckdb_malloc(nfmt * sizeof(int));
        fmt_n_values = (int *)duckdb_malloc(nfmt * sizeof(int));
        fmt_i32 = (int32_t **)duckdb_malloc(nfmt * sizeof(int32_t *));
        fmt_f32 = (float **)duckdb_malloc(nfmt * sizeof(float *));
        fmt_str = (char ***)duckdb_malloc(nfmt * sizeof(char **));
        if (!fmt_loaded || !fmt_ret || !fmt_n_values || !fmt_i32 || !fmt_f32 || !fmt_str) {
            duckdb_function_set_error(info, "read_bcf: out of memory allocating format cache");
            if (fmt_loaded) duckdb_free(fmt_loaded);
            if (fmt_ret) duckdb_free(fmt_ret);
            if (fmt_n_values) duckdb_free(fmt_n_values);
            if (fmt_i32) duckdb_free(fmt_i32);
            if (fmt_f32) duckdb_free(fmt_f32);
            if (fmt_str) duckdb_free(fmt_str);
            duckdb_free(vectors);
            duckdb_data_chunk_set_size(output, 0);
            return;
        }
        memset(fmt_loaded, 0, nfmt * sizeof(int));
        memset(fmt_ret, 0, nfmt * sizeof(int));
        memset(fmt_n_values, 0, nfmt * sizeof(int));
        memset(fmt_i32, 0, nfmt * sizeof(int32_t *));
        memset(fmt_f32, 0, nfmt * sizeof(float *));
        memset(fmt_str, 0, nfmt * sizeof(char **));
    }

    // Per-record INFO cache (same rationale as FORMAT cache above).
    int *info_loaded = NULL;
    int *info_ret = NULL;
    int *info_n_values = NULL;
    int *info_flag = NULL;
    int32_t **info_i32 = NULL;
    float **info_f32 = NULL;
    char **info_str = NULL;
    if (use_persistent_cache) {
        info_loaded = init->info_loaded;
        info_ret = init->info_ret;
        info_n_values = init->info_n_values;
        info_flag = init->info_flag;
        info_i32 = init->info_i32;
        info_f32 = init->info_f32;
        info_str = init->info_str;
    } else if (bind->n_info_fields > 0) {
        size_t ninfo = (size_t)bind->n_info_fields;
        info_loaded = (int *)duckdb_malloc(ninfo * sizeof(int));
        info_ret = (int *)duckdb_malloc(ninfo * sizeof(int));
        info_n_values = (int *)duckdb_malloc(ninfo * sizeof(int));
        info_flag = (int *)duckdb_malloc(ninfo * sizeof(int));
        info_i32 = (int32_t **)duckdb_malloc(ninfo * sizeof(int32_t *));
        info_f32 = (float **)duckdb_malloc(ninfo * sizeof(float *));
        info_str = (char **)duckdb_malloc(ninfo * sizeof(char *));
        if (!info_loaded || !info_ret || !info_n_values || !info_flag || !info_i32 || !info_f32 || !info_str) {
            duckdb_function_set_error(info, "read_bcf: out of memory allocating info cache");
            if (info_loaded) duckdb_free(info_loaded);
            if (info_ret) duckdb_free(info_ret);
            if (info_n_values) duckdb_free(info_n_values);
            if (info_flag) duckdb_free(info_flag);
            if (info_i32) duckdb_free(info_i32);
            if (info_f32) duckdb_free(info_f32);
            if (info_str) duckdb_free(info_str);
            if (fmt_loaded) duckdb_free(fmt_loaded);
            if (fmt_ret) duckdb_free(fmt_ret);
            if (fmt_n_values) duckdb_free(fmt_n_values);
            if (fmt_i32) duckdb_free(fmt_i32);
            if (fmt_f32) duckdb_free(fmt_f32);
            if (fmt_str) duckdb_free(fmt_str);
            duckdb_free(vectors);
            duckdb_data_chunk_set_size(output, 0);
            return;
        }
        memset(info_loaded, 0, ninfo * sizeof(int));
        memset(info_ret, 0, ninfo * sizeof(int));
        memset(info_n_values, 0, ninfo * sizeof(int));
        memset(info_flag, 0, ninfo * sizeof(int));
        memset(info_i32, 0, ninfo * sizeof(int32_t *));
        memset(info_f32, 0, ninfo * sizeof(float *));
        memset(info_str, 0, ninfo * sizeof(char *));
    }

    int scan_error = 0;
    char scan_err[512] = {0};

    // Read records
    while (row_count < vector_size) {
        // In tidy mode, only read a new record when we've emitted all samples
        int need_read = 1;
        if (tidy_mode && init->tidy_record_valid) {
            // We have a buffered record - check if we still have samples to emit
            if (init->tidy_current_sample < bind->n_samples) {
                need_read = 0;
                current_sample = init->tidy_current_sample;
            }
        }

        if (need_read) {
            int ret;

            if (init->itr) {
                if (init->tbx) {
                    // VCF with tabix: read text line then parse. read_bcf_v2
                    // count-only scans can count returned lines directly.
                    ret = tbx_itr_next(init->fp, init->tbx, init->itr, &init->kstr);
                    if (ret >= 0) {
                        if (bind->reader_version >= 2 && init->column_count == 0) {
                            init->kstr.l = 0;
                            ret = 0;
                        } else {
                            ret = vcf_parse1(&init->kstr, init->hdr, init->rec);
                            init->kstr.l = 0;
                        }
                    }
                } else {
                    // BCF with index
                    ret = bcf_itr_next(init->fp, init->itr, init->rec);
                }
            } else {
                if (bind->reader_version >= 2 && init->column_count == 0 &&
                    hts_get_format(init->fp)->format == vcf) {
                    ret = hts_getline(init->fp, '\n', &init->kstr);
                    if (ret >= 0) {
                        init->kstr.l = 0;
                        ret = 0;
                    }
                } else {
                    ret = bcf_read(init->fp, init->hdr, init->rec);
                }
            }

            if (ret < -1) {
                char err[256];
                const char *reader_name = bind->reader_version >= 2 ? "read_bcf_v2" : "read_bcf";
                snprintf(err, sizeof(err), "%s: failed to read or parse BCF/VCF record", reader_name);
                duckdb_function_set_error(info, err);
                init->done = 1;
                break;
            }

            if (ret < 0) {
                // End of current contig/file
                if (init->is_parallel) {
                    // Try to claim next contig
                    if (claim_next_contig(init, global)) {
                        continue;  // Continue reading from new contig
                    }
                }
                init->done = 1;
                break;
            }

            if (bind->sample_filter_active && init->itr && init->idx &&
                (init->unpack_mask & BCF_UN_FMT) &&
                bcf_subset_format(init->hdr, init->rec) < 0) {
                duckdb_function_set_error(info, "read_bcf_v2: failed to subset FORMAT fields for indexed BCF record");
                init->done = 1;
                break;
            }

            if (init->unpack_mask) {
                bcf_unpack(init->rec, init->unpack_mask);
            }

            // For tidy mode, reset sample counter
            if (tidy_mode) {
                init->tidy_current_sample = 0;
                init->tidy_record_valid = 1;
                current_sample = 0;
            }

            // Update debug/progress counters (only when reading a new record)
#if BCF_READER_ENABLE_PROGRESS
            if (!init->timing_initialized) {
                // First record - start timing
                clock_gettime(CLOCK_MONOTONIC, &init->batch_start_time);
                init->last_progress_time = init->batch_start_time;
                init->timing_initialized = 1;
            }
#endif
            if (bind->n_format_fields > 0) {
                memset(fmt_loaded, 0, (size_t)bind->n_format_fields * sizeof(int));
            }
            if (bind->n_info_fields > 0) {
                memset(info_loaded, 0, (size_t)bind->n_info_fields * sizeof(int));
            }
            gt_loaded = 0;
            if (use_persistent_cache) {
                if (init->vep_rec) {
                    vep_record_destroy(init->vep_rec);
                    init->vep_rec = NULL;
                }
                init->vep_loaded = 0;
            }
        }

        // Parse VEP annotation once per record if needed. read_bcf_v2 caches the
        // parsed record across tidy sample rows and DuckDB chunk boundaries.
        vep_record_t* vep_rec = NULL;
        if (init->need_vep) {
            if (use_persistent_cache) {
                if (!init->vep_loaded) {
                    if (bind->reader_version >= 2) {
                        init->vep_rec = vep_record_parse_selected_bcf(
                            bind->vep_schema, init->hdr, init->rec,
                            init->projected_vep_field_indices,
                            init->n_projected_vep_fields
                        );
                    } else {
                        init->vep_rec = vep_record_parse_bcf(bind->vep_schema, init->hdr, init->rec);
                    }
                    init->vep_loaded = 1;
                }
                vep_rec = init->vep_rec;
            } else if (!tidy_mode || current_sample == 0) {
                if (bind->reader_version >= 2) {
                    vep_rec = vep_record_parse_selected_bcf(
                        bind->vep_schema, init->hdr, init->rec,
                        init->projected_vep_field_indices,
                        init->n_projected_vep_fields
                    );
                } else {
                    vep_rec = vep_record_parse_bcf(bind->vep_schema, init->hdr, init->rec);
                }
            }
        }

        idx_t emit_run_len = 1;
        if (bind->reader_version >= 2 && tidy_mode) {
            int samples_remaining = bind->n_samples - current_sample;
            idx_t rows_remaining = vector_size - row_count;
            if (samples_remaining > 0 && rows_remaining > 0) {
                emit_run_len = (idx_t)samples_remaining < rows_remaining
                    ? (idx_t)samples_remaining
                    : rows_remaining;
            }
        }

        // Decode each projected FORMAT field once for this record before
        // filling all sample rows in the current tidy run.
        for (int group_idx = 0; group_idx < init->n_format_groups; group_idx++) {
            bcf_format_group_t *group = &init->format_groups[group_idx];
            if (group->column_count == 0) {
                continue;
            }
            const bcf_projected_col_t *first_col = &init->projected_cols[group->projected_indices[0]];
            if (!bcf_decode_format_projected_field(init, bind, first_col,
                                                   fmt_loaded, fmt_ret, fmt_n_values,
                                                   fmt_i32, fmt_f32, fmt_str,
                                                   &gt_loaded, &gt_ret,
                                                   &gt_n_values, &gt_arr,
                                                   scan_err, sizeof(scan_err))) {
                duckdb_function_set_error(info, scan_err[0] ? scan_err : "read_bcf: failed to decode FORMAT field");
                init->done = 1;
                scan_error = 1;
                break;
            }
        }
        if (scan_error) {
            break;
        }

        if (bind->reader_version >= 2 && emit_run_len > 1) {
            for (idx_t site_i = 0; site_i < init->site_column_count; site_i++) {
                const bcf_projected_col_t *proj_col = &init->projected_cols[init->site_column_indices[site_i]];
                idx_t child_values_per_row = 0;
                if (proj_col->kind == BCF_OUT_ALT) {
                    child_values_per_row = init->rec->n_allele > 1 ? (idx_t)(init->rec->n_allele - 1) : 0;
                } else if (proj_col->kind == BCF_OUT_FILTER) {
                    child_values_per_row = init->rec->d.n_flt == 0 ? 1 : (idx_t)init->rec->d.n_flt;
                } else if (proj_col->kind == BCF_OUT_VEP && vep_rec) {
                    child_values_per_row = (idx_t)vep_rec->n_transcripts;
                }
                if (child_values_per_row > 0) {
                    duckdb_vector vec = vectors[proj_col->out_idx];
                    idx_t current_child_size = duckdb_list_vector_get_size(vec);
                    duckdb_list_vector_reserve(vec, current_child_size + child_values_per_row * emit_run_len);
                }
            }
        }

        for (idx_t emit_i = 0; emit_i < emit_run_len; emit_i++) {
            idx_t row_idx = row_count + emit_i;
            int sample_for_row = tidy_mode ? current_sample + (int)emit_i : current_sample;

            // Process site-level columns using descriptors prepared in local init.
            for (idx_t site_i = 0; site_i < init->site_column_count; site_i++) {
                idx_t i = init->site_column_indices[site_i];
                const bcf_projected_col_t *proj_col = &init->projected_cols[i];
                duckdb_vector vec = vectors[proj_col->out_idx];

                switch (proj_col->kind) {
                case BCF_OUT_CHROM: {
                    const char* chrom = bcf_hdr_id2name(init->hdr, init->rec->rid);
                    duckdb_vector_assign_string_element(vec, row_idx, chrom ? chrom : ".");
                    break;
                }
                case BCF_OUT_POS: {
                    int64_t* data = (int64_t*)duckdb_vector_get_data(vec);
                    data[row_idx] = init->rec->pos + 1;  // 1-based
                    break;
                }
                case BCF_OUT_ID: {
                    const char* id = init->rec->d.id;
                    if (id && strcmp(id, ".") != 0) {
                        duckdb_vector_assign_string_element(vec, row_idx, id);
                    } else {
                        duckdb_vector_ensure_validity_writable(vec);
                        uint64_t* validity = duckdb_vector_get_validity(vec);
                        set_validity_bit(validity, row_idx, 0);
                    }
                    break;
                }
                case BCF_OUT_REF: {
                    const char* ref = init->rec->d.allele[0];
                    duckdb_vector_assign_string_element(vec, row_idx, ref ? ref : ".");
                    break;
                }
                case BCF_OUT_ALT: {
                    // ALT is a LIST(VARCHAR)
                    duckdb_list_entry entry;
                    entry.offset = duckdb_list_vector_get_size(vec);
                    entry.length = init->rec->n_allele > 1 ? init->rec->n_allele - 1 : 0;

                    duckdb_vector child_vec = duckdb_list_vector_get_child(vec);

                    // Reserve and set space for all ALT alleles
                    if (entry.length > 0) {
                        duckdb_list_vector_reserve(vec, entry.offset + entry.length);
                        duckdb_list_vector_set_size(vec, entry.offset + entry.length);

                        for (int a = 1; a < init->rec->n_allele; a++) {
                            duckdb_vector_assign_string_element(child_vec, entry.offset + a - 1,
                                                                init->rec->d.allele[a]);
                        }
                    }

                    duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
                    list_data[row_idx] = entry;
                    break;
                }
                case BCF_OUT_QUAL: {
                    double* data = (double*)duckdb_vector_get_data(vec);
                    if (bcf_float_is_missing(init->rec->qual)) {
                        duckdb_vector_ensure_validity_writable(vec);
                        uint64_t* validity = duckdb_vector_get_validity(vec);
                        set_validity_bit(validity, row_idx, 0);
                        data[row_idx] = 0.0;
                    } else {
                        data[row_idx] = init->rec->qual;
                    }
                    break;
                }
                case BCF_OUT_FILTER: {
                    // FILTER is a LIST(VARCHAR)
                    duckdb_list_entry entry;
                    entry.offset = duckdb_list_vector_get_size(vec);

                    duckdb_vector child_vec = duckdb_list_vector_get_child(vec);

                    if (init->rec->d.n_flt == 0) {
                        // No filters means PASS
                        entry.length = 1;
                        duckdb_list_vector_reserve(vec, entry.offset + entry.length);
                        duckdb_list_vector_set_size(vec, entry.offset + entry.length);
                        duckdb_vector_assign_string_element(child_vec, entry.offset, "PASS");
                    } else {
                        entry.length = init->rec->d.n_flt;
                        // Reserve space for all filters at once
                        duckdb_list_vector_reserve(vec, entry.offset + entry.length);
                        duckdb_list_vector_set_size(vec, entry.offset + entry.length);
                        for (int f = 0; f < init->rec->d.n_flt; f++) {
                            const char* flt_name = bcf_hdr_int2id(init->hdr, BCF_DT_ID,
                                                                  init->rec->d.flt[f]);
                            duckdb_vector_assign_string_element(child_vec, entry.offset + f,
                                                                flt_name ? flt_name : ".");
                        }
                    }

                    duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
                    list_data[row_idx] = entry;
                    break;
                }
                case BCF_OUT_VEP: {
                    int schema_field_idx = proj_col->schema_field_idx;
                    const vep_field_t* field = proj_col->vep_field;

                    duckdb_list_entry entry;
                    entry.offset = duckdb_list_vector_get_size(vec);
                    entry.length = (vep_rec && field) ? vep_rec->n_transcripts : 0;

                    duckdb_vector child_vec = duckdb_list_vector_get_child(vec);

                    if (entry.length > 0) {
                        duckdb_vector_ensure_validity_writable(vec);
                        uint64_t* parent_validity = duckdb_vector_get_validity(vec);
                        set_validity_bit(parent_validity, row_idx, 1);

                        duckdb_list_vector_reserve(vec, entry.offset + entry.length);
                        duckdb_list_vector_set_size(vec, entry.offset + entry.length);

                        // Ensure child validity writable for nulls
                        duckdb_vector_ensure_validity_writable(child_vec);
                        uint64_t* child_validity = duckdb_vector_get_validity(child_vec);

                        if (field->type == VEP_TYPE_STRING) {
                            for (idx_t t = 0; t < entry.length; t++) {
                                const vep_value_t* val = vep_record_get_value(vep_rec, t, schema_field_idx);
                                if (val && !val->is_missing && val->str_value) {
                                    duckdb_vector_assign_string_element(child_vec, entry.offset + t, val->str_value);
                                    set_validity_bit(child_validity, entry.offset + t, 1);
                                } else {
                                    set_validity_bit(child_validity, entry.offset + t, 0);
                                }
                            }
                        } else if (field->type == VEP_TYPE_INTEGER) {
                            int32_t* data = (int32_t*)duckdb_vector_get_data(child_vec);
                            for (idx_t t = 0; t < entry.length; t++) {
                                const vep_value_t* val = vep_record_get_value(vep_rec, t, schema_field_idx);
                                if (val && !val->is_missing) {
                                    data[entry.offset + t] = val->int_value;
                                    set_validity_bit(child_validity, entry.offset + t, 1);
                                } else {
                                    set_validity_bit(child_validity, entry.offset + t, 0);
                                }
                            }
                        } else if (field->type == VEP_TYPE_FLOAT) {
                            float* data = (float*)duckdb_vector_get_data(child_vec);
                            for (idx_t t = 0; t < entry.length; t++) {
                                const vep_value_t* val = vep_record_get_value(vep_rec, t, schema_field_idx);
                                if (val && !val->is_missing) {
                                    data[entry.offset + t] = val->float_value;
                                    set_validity_bit(child_validity, entry.offset + t, 1);
                                } else {
                                    set_validity_bit(child_validity, entry.offset + t, 0);
                                }
                            }
                        } else if (field->type == VEP_TYPE_FLAG) {
                            bool* data = (bool*)duckdb_vector_get_data(child_vec);
                            for (idx_t t = 0; t < entry.length; t++) {
                                const vep_value_t* val = vep_record_get_value(vep_rec, t, schema_field_idx);
                                if (val && !val->is_missing) {
                                    data[entry.offset + t] = 1;
                                    set_validity_bit(child_validity, entry.offset + t, 1);
                                } else {
                                    set_validity_bit(child_validity, entry.offset + t, 0);
                                }
                            }
                        }
                    } else {
                        // No VEP data for this record
                        duckdb_vector_ensure_validity_writable(vec);
                        uint64_t* validity = duckdb_vector_get_validity(vec);
                        set_validity_bit(validity, row_idx, 0);
                        duckdb_list_vector_set_size(vec, entry.offset);
                    }

                    duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
                    list_data[row_idx] = entry;
                    break;
                }
                case BCF_OUT_INFO: {
                    // INFO field
                    int field_idx = proj_col->field_idx;
                    const field_meta_t* field = proj_col->field;
                    const char* tag = field->name;

                    if (proj_col->header_type == BCF_HT_FLAG) {
                        // Boolean field
                        bool* data = (bool*)duckdb_vector_get_data(vec);
                        if (!info_loaded[field_idx]) {
                            if (!bcf_record_has_info_field(init->rec, field)) {
                                info_ret[field_idx] = 0;
                            } else {
                                int* dummy = NULL;
                                int ndummy = 0;
                                info_ret[field_idx] = bcf_get_info_flag(init->hdr, init->rec, tag, &dummy, &ndummy);
                                if (dummy) free(dummy);  // Only free if allocated
                            }
                            info_flag[field_idx] = (info_ret[field_idx] == 1);
                            info_loaded[field_idx] = 1;
                        }
                        data[row_idx] = (info_flag[field_idx] != 0);
                    }
                    else if (proj_col->header_type == BCF_HT_INT) {
                        int32_t* values = NULL;
                        int ret_info = 0;
                        if (!info_loaded[field_idx]) {
                            if (!bcf_record_has_info_field(init->rec, field)) {
                                info_ret[field_idx] = -3;
                            } else if (bcf_reader_input_is_bcf(init->fp) &&
                                       !bcf_preflight_info_encoded_type(init->hdr, init->rec, field,
                                                                       bcf_reader_name(bind),
                                                                       scan_err, sizeof(scan_err))) {
                                if (!bcf_handle_decode_diagnostic(bind->decode_error_policy, scan_err)) {
                                    duckdb_function_set_error(info, scan_err[0] ? scan_err : "read_bcf: failed to decode INFO field");
                                    init->done = 1;
                                    scan_error = 1;
                                    break;
                                }
                                info_ret[field_idx] = -3;
                            } else {
                                info_ret[field_idx] = bcf_get_info_int32(init->hdr, init->rec, tag,
                                                                          &info_i32[field_idx],
                                                                          &info_n_values[field_idx]);
                            }
                            info_loaded[field_idx] = 1;
                        }
                        values = info_i32[field_idx];
                        ret_info = info_ret[field_idx];
                        if (!bcf_check_decode_ret(bcf_reader_name(bind), "INFO", tag, init->hdr, init->rec,
                                                  ret_info, bind->decode_error_policy, scan_err, sizeof(scan_err))) {
                            duckdb_function_set_error(info, scan_err[0] ? scan_err : "read_bcf: failed to decode INFO field");
                            init->done = 1;
                            scan_error = 1;
                            break;
                        }

                        assign_int32_field(vec, row_idx, values, ret_info, proj_col->is_list);
                    }
                    else if (proj_col->header_type == BCF_HT_REAL) {
                        float* values = NULL;
                        int ret_info = 0;
                        if (!info_loaded[field_idx]) {
                            if (!bcf_record_has_info_field(init->rec, field)) {
                                info_ret[field_idx] = -3;
                            } else if (bcf_reader_input_is_bcf(init->fp) &&
                                       !bcf_preflight_info_encoded_type(init->hdr, init->rec, field,
                                                                       bcf_reader_name(bind),
                                                                       scan_err, sizeof(scan_err))) {
                                if (!bcf_handle_decode_diagnostic(bind->decode_error_policy, scan_err)) {
                                    duckdb_function_set_error(info, scan_err[0] ? scan_err : "read_bcf: failed to decode INFO field");
                                    init->done = 1;
                                    scan_error = 1;
                                    break;
                                }
                                info_ret[field_idx] = -3;
                            } else {
                                info_ret[field_idx] = bcf_get_info_float(init->hdr, init->rec, tag,
                                                                          &info_f32[field_idx],
                                                                          &info_n_values[field_idx]);
                            }
                            info_loaded[field_idx] = 1;
                        }
                        values = info_f32[field_idx];
                        ret_info = info_ret[field_idx];
                        if (!bcf_check_decode_ret(bcf_reader_name(bind), "INFO", tag, init->hdr, init->rec,
                                                  ret_info, bind->decode_error_policy, scan_err, sizeof(scan_err))) {
                            duckdb_function_set_error(info, scan_err[0] ? scan_err : "read_bcf: failed to decode INFO field");
                            init->done = 1;
                            scan_error = 1;
                            break;
                        }

                        assign_float_field(vec, row_idx, values, ret_info, proj_col->is_list);
                    }
                    else {
                        // String type
                        char* value = NULL;
                        int ret_info = 0;
                        if (!info_loaded[field_idx]) {
                            if (!bcf_record_has_info_field(init->rec, field)) {
                                info_ret[field_idx] = -3;
                            } else if (bcf_reader_input_is_bcf(init->fp) &&
                                       !bcf_preflight_info_encoded_type(init->hdr, init->rec, field,
                                                                       bcf_reader_name(bind),
                                                                       scan_err, sizeof(scan_err))) {
                                if (!bcf_handle_decode_diagnostic(bind->decode_error_policy, scan_err)) {
                                    duckdb_function_set_error(info, scan_err[0] ? scan_err : "read_bcf: failed to decode INFO field");
                                    init->done = 1;
                                    scan_error = 1;
                                    break;
                                }
                                info_ret[field_idx] = -3;
                            } else {
                                info_ret[field_idx] = bcf_get_info_string(init->hdr, init->rec, tag,
                                                                           &info_str[field_idx],
                                                                           &info_n_values[field_idx]);
                            }
                            info_loaded[field_idx] = 1;
                        }
                        value = info_str[field_idx];
                        ret_info = info_ret[field_idx];
                        if (!bcf_check_decode_ret(bcf_reader_name(bind), "INFO", tag, init->hdr, init->rec,
                                                  ret_info, bind->decode_error_policy, scan_err, sizeof(scan_err))) {
                            duckdb_function_set_error(info, scan_err[0] ? scan_err : "read_bcf: failed to decode INFO field");
                            init->done = 1;
                            scan_error = 1;
                            break;
                        }

                        assign_string_field(vec, row_idx, ret_info > 0 ? value : NULL, proj_col->is_list);
                    }
                    break;
                }
                case BCF_OUT_SAMPLE_ID:
                    // SAMPLE_ID column in tidy mode
                    duckdb_vector_assign_string_element(vec, row_idx, bind->sample_names[sample_for_row]);
                    break;

                case BCF_OUT_INVALID:
                default:
                    break;
                }
                if (scan_error) {
                    break;
                }
            }

            if (scan_error) {
                break;
            }

            // FORMAT columns are grouped by FORMAT field; fill every projected
            // sample column from the decoded group buffers.
            for (int group_idx = 0; group_idx < init->n_format_groups; group_idx++) {
                bcf_format_group_t *group = &init->format_groups[group_idx];
                for (idx_t group_col_idx = 0; group_col_idx < group->column_count; group_col_idx++) {
                    const bcf_projected_col_t *proj_col = &init->projected_cols[group->projected_indices[group_col_idx]];
                    duckdb_vector vec = vectors[proj_col->out_idx];
                    bcf_fill_format_projected_col(init, bind, proj_col, vec,
                                                  row_idx, sample_for_row,
                                                  fmt_ret, fmt_i32, fmt_f32,
                                                  fmt_str, gt_ret, gt_arr);
                }
            }

        }

        if (scan_error) {
            if (vep_rec && !use_persistent_cache) {
                vep_record_destroy(vep_rec);
            }
            break;
        }

        if (vep_rec && !use_persistent_cache) {
            vep_record_destroy(vep_rec);
        }

        row_count += emit_run_len;
        init->current_row += emit_run_len;

        // In tidy mode, advance by the emitted sample run (or mark record as consumed)
        if (tidy_mode) {
            init->tidy_current_sample += (int)emit_run_len;
            if (init->tidy_current_sample >= bind->n_samples) {
                // All samples emitted for this record - next iteration will read new record
                init->tidy_record_valid = 0;
                init->total_records_processed++;
            }
        } else {
            init->total_records_processed++;
        }

        // Print progress every N records (only count actual VCF records, not per-sample rows)
#if BCF_READER_ENABLE_PROGRESS
        if (!tidy_mode || !init->tidy_record_valid) {
            if (init->is_parallel && init->contig_name) {
                char context[256];
                snprintf(context, sizeof(context), "scan (contig: %s)", init->contig_name);
                print_progress(init, context);
            } else {
                print_progress(init, "scan");
            }
        }
#endif
    }

    // Cleanup cached vectors
    duckdb_free(vectors);
    if (use_persistent_cache) {
        init->gt_loaded = gt_loaded;
        init->gt_ret = gt_ret;
        init->gt_n_values = gt_n_values;
        init->gt_arr = gt_arr;
    } else {
        if (bind->n_format_fields > 0) {
            for (int i = 0; i < bind->n_format_fields; i++) {
                if (fmt_i32 && fmt_i32[i]) free(fmt_i32[i]);
                if (fmt_f32 && fmt_f32[i]) free(fmt_f32[i]);
                if (fmt_str && fmt_str[i]) free_bcf_format_string_array(fmt_str[i]);
            }
            if (gt_arr) free(gt_arr);
            if (fmt_loaded) duckdb_free(fmt_loaded);
            if (fmt_ret) duckdb_free(fmt_ret);
            if (fmt_n_values) duckdb_free(fmt_n_values);
            if (fmt_i32) duckdb_free(fmt_i32);
            if (fmt_f32) duckdb_free(fmt_f32);
            if (fmt_str) duckdb_free(fmt_str);
        }
        if (bind->n_info_fields > 0) {
            for (int i = 0; i < bind->n_info_fields; i++) {
                if (info_i32 && info_i32[i]) free(info_i32[i]);
                if (info_f32 && info_f32[i]) free(info_f32[i]);
                if (info_str && info_str[i]) free(info_str[i]);
            }
            if (info_loaded) duckdb_free(info_loaded);
            if (info_ret) duckdb_free(info_ret);
            if (info_n_values) duckdb_free(info_n_values);
            if (info_flag) duckdb_free(info_flag);
            if (info_i32) duckdb_free(info_i32);
            if (info_f32) duckdb_free(info_f32);
            if (info_str) duckdb_free(info_str);
        }
    }

    duckdb_data_chunk_set_size(output, row_count);
}

// =============================================================================
// Narrow BCF/VCF appender benchmark table function
// =============================================================================

typedef struct {
    duckdb_connection connection;
} bcf_appender_extra_t;

typedef struct {
    duckdb_connection connection;
    char *file_path;
    char *target_table;
    char *region;
    char **regions;
    unsigned int n_regions;
    int tidy_format;
    int overwrite;
    int include_file_offset;
    int region_threads;
    int decompression_threads;
    bcf_decode_error_policy_t decode_error_policy;
} bcf_appender_bind_t;

typedef struct {
    int done;
} bcf_appender_init_t;

typedef struct {
    duckdb_data_chunk chunk;
    duckdb_logical_type types[7];
    duckdb_vector chrom;
    duckdb_vector pos;
    duckdb_vector ref;
    duckdb_vector alt;
    duckdb_vector sample_id;
    duckdb_vector file_offset;
    duckdb_vector format_gt;
    idx_t row_count;
    idx_t column_count;
    int include_file_offset;
} bcf_appender_chunk_t;

typedef struct {
    int tid;
    char **regions;
    unsigned int n_regions;
    unsigned int capacity;
} bcf_appender_contig_task_t;

typedef struct {
    const bcf_appender_bind_t *bind;
    bcf_appender_contig_task_t *tasks;
    int n_tasks;
    volatile int next_task_idx;
    volatile int failed;
    char error[1024];
} bcf_appender_region_shared_t;

typedef struct {
    bcf_appender_region_shared_t *shared;
    duckdb_appender appender;
    uint64_t rows_written;
} bcf_appender_region_worker_t;

static void destroy_bcf_appender_extra(void *data) {
    bcf_appender_extra_t *extra = (bcf_appender_extra_t *)data;
    if (!extra) return;
    if (extra->connection) duckdb_disconnect(&extra->connection);
    duckdb_free(extra);
}

static void destroy_bcf_appender_bind(void *data) {
    bcf_appender_bind_t *bind = (bcf_appender_bind_t *)data;
    if (!bind) return;
    if (bind->file_path) duckdb_free(bind->file_path);
    if (bind->target_table) duckdb_free(bind->target_table);
    if (bind->region) duckdb_free(bind->region);
    if (bind->regions) {
        for (unsigned int i = 0; i < bind->n_regions; i++) {
            if (bind->regions[i]) duckdb_free(bind->regions[i]);
        }
        duckdb_free(bind->regions);
    }
    duckdb_free(bind);
}

static void destroy_bcf_appender_init(void *data) {
    if (data) duckdb_free(data);
}

static int bcf_appender_valid_table_name(const char *s) {
    if (!s || !s[0]) return 0;
    if (!((s[0] >= 'A' && s[0] <= 'Z') || (s[0] >= 'a' && s[0] <= 'z') || s[0] == '_')) {
        return 0;
    }
    for (const char *p = s + 1; *p; p++) {
        if ((*p >= 'A' && *p <= 'Z') || (*p >= 'a' && *p <= 'z') ||
            (*p >= '0' && *p <= '9') || *p == '_') {
            continue;
        }
        return 0;
    }
    return 1;
}

static int bcf_appender_exec_sql(duckdb_connection con, const char *sql,
                                 char *err, size_t err_size) {
    duckdb_result res;
    duckdb_state state = duckdb_query(con, sql, &res);
    if (state == DuckDBError) {
        const char *msg = duckdb_result_error(&res);
        snprintf(err, err_size, "%s", msg ? msg : "DuckDB query failed");
        duckdb_destroy_result(&res);
        return 0;
    }
    duckdb_destroy_result(&res);
    return 1;
}

static int bcf_appender_create_target(duckdb_connection con, const bcf_appender_bind_t *bind,
                                      char *err, size_t err_size) {
    char sql[2048];

    if (bind->overwrite) {
        snprintf(sql, sizeof(sql), "DROP TABLE IF EXISTS \"%s\"", bind->target_table);
        if (!bcf_appender_exec_sql(con, sql, err, err_size)) return 0;
    }

    if (bind->include_file_offset) {
        snprintf(sql, sizeof(sql),
                 "CREATE TABLE IF NOT EXISTS \"%s\" ("
                 "CHROM VARCHAR, POS BIGINT, REF VARCHAR, ALT VARCHAR, "
                 "SAMPLE_ID VARCHAR, FILE_OFFSET UBIGINT, FORMAT_GT VARCHAR)",
                 bind->target_table);
    } else {
        snprintf(sql, sizeof(sql),
                 "CREATE TABLE IF NOT EXISTS \"%s\" ("
                 "CHROM VARCHAR, POS BIGINT, REF VARCHAR, ALT VARCHAR, "
                 "SAMPLE_ID VARCHAR, FORMAT_GT VARCHAR)",
                 bind->target_table);
    }
    return bcf_appender_exec_sql(con, sql, err, err_size);
}

static void bcf_appender_chunk_refresh_vectors(bcf_appender_chunk_t *chunk) {
    idx_t col = 0;
    chunk->chrom = duckdb_data_chunk_get_vector(chunk->chunk, col++);
    chunk->pos = duckdb_data_chunk_get_vector(chunk->chunk, col++);
    chunk->ref = duckdb_data_chunk_get_vector(chunk->chunk, col++);
    chunk->alt = duckdb_data_chunk_get_vector(chunk->chunk, col++);
    chunk->sample_id = duckdb_data_chunk_get_vector(chunk->chunk, col++);
    if (chunk->include_file_offset) {
        chunk->file_offset = duckdb_data_chunk_get_vector(chunk->chunk, col++);
    } else {
        chunk->file_offset = NULL;
    }
    chunk->format_gt = duckdb_data_chunk_get_vector(chunk->chunk, col++);
}

static int bcf_appender_chunk_init(bcf_appender_chunk_t *chunk,
                                   int include_file_offset) {
    idx_t col = 0;
    memset(chunk, 0, sizeof(*chunk));
    chunk->include_file_offset = include_file_offset ? 1 : 0;
    chunk->column_count = 6 + (chunk->include_file_offset ? 1 : 0);
    chunk->types[col++] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    chunk->types[col++] = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    chunk->types[col++] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    chunk->types[col++] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    chunk->types[col++] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    if (chunk->include_file_offset) {
        chunk->types[col++] = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    }
    chunk->types[col++] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    chunk->chunk = duckdb_create_data_chunk(chunk->types, chunk->column_count);
    if (!chunk->chunk) return 0;
    bcf_appender_chunk_refresh_vectors(chunk);
    return 1;
}

static void bcf_appender_chunk_destroy(bcf_appender_chunk_t *chunk) {
    if (!chunk) return;
    if (chunk->chunk) duckdb_destroy_data_chunk(&chunk->chunk);
    for (idx_t i = 0; i < chunk->column_count; i++) {
        if (chunk->types[i]) duckdb_destroy_logical_type(&chunk->types[i]);
    }
    memset(chunk, 0, sizeof(*chunk));
}

static int bcf_appender_chunk_flush(duckdb_appender appender, bcf_appender_chunk_t *chunk,
                                    char *err, size_t err_size) {
    if (chunk->row_count == 0) return 1;
    duckdb_data_chunk_set_size(chunk->chunk, chunk->row_count);
    if (duckdb_append_data_chunk(appender, chunk->chunk) == DuckDBError) {
        const char *msg = duckdb_appender_error(appender);
        snprintf(err, err_size, "%s", msg ? msg : "duckdb_append_data_chunk failed");
        return 0;
    }
    duckdb_data_chunk_reset(chunk->chunk);
    bcf_appender_chunk_refresh_vectors(chunk);
    chunk->row_count = 0;
    return 1;
}

static void bcf_appender_set_null(duckdb_vector vec, idx_t row) {
    duckdb_vector_ensure_validity_writable(vec);
    uint64_t *validity = duckdb_vector_get_validity(vec);
    set_validity_bit(validity, row, 0);
}

static int bcf_appender_format_gt(const int32_t *sample_gt, int ploidy, kstring_t *gt) {
    gt->l = 0;
    if (gt->s) gt->s[0] = '\0';
    for (int p = 0; p < ploidy; p++) {
        if (sample_gt[p] == bcf_int32_vector_end) break;
        if (p > 0 && kputc(bcf_gt_is_phased(sample_gt[p]) ? '|' : '/', gt) < 0) return 0;
        if (bcf_gt_is_missing(sample_gt[p])) {
            if (kputc('.', gt) < 0) return 0;
        } else {
            if (kputw(bcf_gt_allele(sample_gt[p]), gt) < 0) return 0;
        }
    }
    return 1;
}

static int bcf_appender_append_row(duckdb_appender appender, bcf_appender_chunk_t *chunk,
                                   const char *chrom, int64_t pos, const char *ref,
                                   const char *alt, const char *sample_id,
                                   int has_file_offset, uint64_t file_offset,
                                   const char *gt, size_t gt_len,
                                   uint64_t *rows_written,
                                   char *err, size_t err_size) {
    if (chunk->row_count >= duckdb_vector_size()) {
        if (!bcf_appender_chunk_flush(appender, chunk, err, err_size)) return 0;
    }

    idx_t row = chunk->row_count;
    int64_t *pos_data = (int64_t *)duckdb_vector_get_data(chunk->pos);
    uint64_t *offset_data = chunk->include_file_offset
        ? (uint64_t *)duckdb_vector_get_data(chunk->file_offset)
        : NULL;

    duckdb_vector_assign_string_element(chunk->chrom, row, chrom ? chrom : "");
    pos_data[row] = pos;
    duckdb_vector_assign_string_element(chunk->ref, row, ref ? ref : "");
    duckdb_vector_assign_string_element(chunk->alt, row, alt ? alt : "");
    if (sample_id) {
        duckdb_vector_assign_string_element(chunk->sample_id, row, sample_id);
    } else {
        bcf_appender_set_null(chunk->sample_id, row);
    }
    if (chunk->include_file_offset) {
        if (has_file_offset) {
            offset_data[row] = file_offset;
        } else {
            bcf_appender_set_null(chunk->file_offset, row);
        }
    }
    if (gt && gt_len > 0) {
        duckdb_vector_assign_string_element_len(chunk->format_gt, row, gt, gt_len);
    } else {
        bcf_appender_set_null(chunk->format_gt, row);
    }

    chunk->row_count++;
    (*rows_written)++;
    return 1;
}

static int bcf_appender_next_record(htsFile *fp, bcf_hdr_t *hdr, bcf1_t *rec,
                                    tbx_t *tbx, hts_itr_t *itr, kstring_t *line,
                                    int want_file_offset,
                                    uint64_t *file_offset, int *has_file_offset,
                                    char *err, size_t err_size) {
    int ret;
    if (itr) {
        if (tbx) {
            ret = tbx_itr_next(fp, tbx, itr, line);
            if (ret >= 0) {
                ret = vcf_parse1(line, hdr, rec);
                line->l = 0;
            }
        } else {
            ret = bcf_itr_next(fp, itr, rec);
        }
    } else {
        ret = bcf_read(fp, hdr, rec);
    }

    if (ret >= 0) {
        if (want_file_offset) {
            BGZF *bgzf = hts_get_bgzfp(fp);
            if (bgzf) {
                *file_offset = (uint64_t)bgzf_tell(bgzf);
                *has_file_offset = 1;
            } else {
                *file_offset = 0;
                *has_file_offset = 0;
            }
        } else {
            *file_offset = 0;
            *has_file_offset = 0;
        }
        return 1;
    }

    if (ret < -1) {
        snprintf(err, err_size, "read_bcf_appender: failed to read or parse BCF/VCF record");
        return -1;
    }
    return 0;
}

static int bcf_appender_append_record(duckdb_appender appender,
                                      bcf_appender_chunk_t *chunk,
                                      bcf_hdr_t *hdr,
                                      bcf1_t *rec,
                                      int tidy_format,
                                      int has_file_offset,
                                      uint64_t file_offset,
                                      int input_is_bcf,
                                      bcf_decode_error_policy_t decode_error_policy,
                                      kstring_t *alt,
                                      kstring_t *gt,
                                      int32_t **gt_arr,
                                      int *gt_arr_n,
                                      uint64_t *rows_written,
                                      char *err,
                                      size_t err_size) {
    bcf_unpack(rec, BCF_UN_STR);
    const char *chrom = bcf_hdr_id2name(hdr, rec->rid);
    const char *ref = (rec->n_allele > 0 && rec->d.allele) ? rec->d.allele[0] : "";

    alt->l = 0;
    if (alt->s) alt->s[0] = '\0';
    if (rec->n_allele <= 1 || !rec->d.allele) {
        if (kputc('.', alt) < 0) {
            snprintf(err, err_size, "read_bcf_appender: failed to format ALT");
            return 0;
        }
    } else {
        for (int a = 1; a < rec->n_allele; a++) {
            if (a > 1 && kputc(',', alt) < 0) {
                snprintf(err, err_size, "read_bcf_appender: failed to format ALT");
                return 0;
            }
            if (kputs(rec->d.allele[a], alt) < 0) {
                snprintf(err, err_size, "read_bcf_appender: failed to format ALT");
                return 0;
            }
        }
    }

    int n_samples = bcf_hdr_nsamples(hdr);
    int gt_ret = -3;
    int ploidy = 0;
    if (n_samples > 0) {
        int gt_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "GT");
        if (gt_id >= 0 && bcf_hdr_idinfo_exists(hdr, BCF_HL_FMT, gt_id)) {
            if (input_is_bcf && bcf_get_fmt_id(rec, gt_id) &&
                !bcf_preflight_format_encoded_type(hdr, rec, gt_id, "GT",
                                                   BCF_OUT_FORMAT_GT, BCF_HT_STR,
                                                   "read_bcf_appender",
                                                   err, err_size)) {
                if (!bcf_handle_decode_diagnostic(decode_error_policy, err)) {
                    return 0;
                }
                gt_ret = -3;
            } else {
                gt_ret = bcf_get_genotypes(hdr, rec, gt_arr, gt_arr_n);
            }
            if (!bcf_check_decode_ret("read_bcf_appender", "FORMAT", "GT",
                                      hdr, rec, gt_ret, decode_error_policy,
                                      err, err_size)) {
                return 0;
            }
            if (!bcf_handle_format_width("read_bcf_appender", "GT",
                                         hdr, rec, &gt_ret, n_samples,
                                         decode_error_policy,
                                         err, err_size)) {
                return 0;
            }
            if (gt_ret > 0) ploidy = gt_ret / n_samples;
        }
    }

    if (tidy_format && n_samples > 0) {
        for (int s = 0; s < n_samples; s++) {
            const char *gt_s = NULL;
            size_t gt_len = 0;
            if (gt_ret > 0 && *gt_arr && ploidy > 0) {
                if (!bcf_appender_format_gt(*gt_arr + s * ploidy, ploidy, gt)) {
                    snprintf(err, err_size, "read_bcf_appender: failed to format GT");
                    return 0;
                }
                if (gt->l > 0) {
                    gt_s = gt->s;
                    gt_len = gt->l;
                }
            }
            if (!bcf_appender_append_row(appender, chunk,
                                         chrom, (int64_t)rec->pos + 1,
                                         ref, alt->s ? alt->s : "", hdr->samples[s],
                                         has_file_offset, file_offset,
                                         gt_s, gt_len, rows_written,
                                         err, err_size)) {
                return 0;
            }
        }
    } else {
        if (!bcf_appender_append_row(appender, chunk,
                                     chrom, (int64_t)rec->pos + 1,
                                     ref, alt->s ? alt->s : "", NULL,
                                     has_file_offset, file_offset,
                                     NULL, 0, rows_written,
                                     err, err_size)) {
            return 0;
        }
    }

    return 1;
}

static int bcf_appender_bcf_name2id(void *data, const char *name) {
    return bcf_hdr_name2id((bcf_hdr_t *)data, name);
}

static int bcf_appender_tbx_name2id(void *data, const char *name) {
    return tbx_name2id((tbx_t *)data, name);
}

static void bcf_appender_destroy_contig_tasks(bcf_appender_contig_task_t *tasks,
                                               int n_tasks) {
    if (!tasks) return;
    for (int i = 0; i < n_tasks; i++) {
        if (!tasks[i].regions) continue;
        for (unsigned int j = 0; j < tasks[i].n_regions; j++) {
            free(tasks[i].regions[j]);
        }
        free(tasks[i].regions);
    }
    free(tasks);
}

static int bcf_appender_contig_task_add(bcf_appender_contig_task_t *task,
                                        const char *region) {
    if (task->n_regions == task->capacity) {
        unsigned int new_capacity = task->capacity ? task->capacity * 2 : 4;
        size_t new_bytes;
        if (new_capacity < task->capacity) return 0;
        new_bytes = (size_t)new_capacity * sizeof(*task->regions);
        if (new_capacity != 0 &&
            new_bytes / (size_t)new_capacity != sizeof(*task->regions)) return 0;
        char **new_regions = (char **)realloc(
            task->regions, new_bytes
        );
        if (!new_regions) return 0;
        task->regions = new_regions;
        task->capacity = new_capacity;
    }
    task->regions[task->n_regions] = strdup(region);
    if (!task->regions[task->n_regions]) return 0;
    task->n_regions++;
    return 1;
}

/* Build one job per primary contig.  Every job retains one native htslib
 * multi-region iterator, so overlapping or repeated intervals on that contig
 * are compacted and a spanning record is emitted once. */
static int bcf_appender_build_contig_tasks(const bcf_appender_bind_t *bind,
                                           bcf_appender_contig_task_t **out_tasks,
                                           int *out_n_tasks,
                                           char *err,
                                           size_t err_size) {
    htsFile *fp = NULL;
    bcf_hdr_t *hdr = NULL;
    hts_idx_t *idx = NULL;
    tbx_t *tbx = NULL;
    bcf_appender_contig_task_t *tasks = NULL;
    int n_tasks = 0;
    int ok = 0;

    *out_tasks = NULL;
    *out_n_tasks = 0;

    fp = hts_open(bind->file_path, "r");
    if (!fp) {
        snprintf(err, err_size, "read_bcf_appender: failed to open BCF/VCF file");
        goto cleanup;
    }
    duckhts_apply_remote_hts_tuning(fp, bind->file_path,
                                    DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION);
    hdr = bcf_hdr_read(fp);
    if (!hdr) {
        snprintf(err, err_size, "read_bcf_appender: failed to read BCF/VCF header");
        goto cleanup;
    }

    if (hts_get_format(fp)->format == bcf) {
        idx = bcf_index_load3(bind->file_path, NULL,
                              HTS_IDX_SAVE_REMOTE | HTS_IDX_SILENT_FAIL);
    } else {
        tbx = tbx_index_load3(bind->file_path, NULL,
                              HTS_IDX_SAVE_REMOTE | HTS_IDX_SILENT_FAIL);
        if (!tbx) {
            idx = bcf_index_load3(bind->file_path, NULL,
                                  HTS_IDX_SAVE_REMOTE | HTS_IDX_SILENT_FAIL);
        }
    }
    if (!idx && !tbx) {
        snprintf(err, err_size, "read_bcf_appender: parallel region query requires an index");
        goto cleanup;
    }

    tasks = (bcf_appender_contig_task_t *)calloc(bind->n_regions, sizeof(*tasks));
    if (!tasks) {
        snprintf(err, err_size, "read_bcf_appender: out of memory creating contig jobs");
        goto cleanup;
    }

    for (unsigned int i = 0; i < bind->n_regions; i++) {
        int tid = -1;
        hts_pos_t beg;
        hts_pos_t end;
        const char *parsed;

        if (strcmp(bind->regions[i], ".") == 0) {
            tid = HTS_IDX_START;
            parsed = bind->regions[i] + 1;
        } else if (strcmp(bind->regions[i], "*") == 0) {
            tid = HTS_IDX_NOCOOR;
            parsed = bind->regions[i] + 1;
        } else {
            parsed = hts_parse_region(
                bind->regions[i], &tid, &beg, &end,
                tbx ? bcf_appender_tbx_name2id : bcf_appender_bcf_name2id,
                tbx ? (void *)tbx : (void *)hdr,
                HTS_PARSE_THOUSANDS_SEP
            );
        }
        if (!parsed) {
            if (tid < -1) {
                snprintf(err, err_size, "read_bcf_appender: failed to parse region %s",
                         bind->regions[i]);
                goto cleanup;
            }
            hts_log_warning("Region '%s' specifies an unknown reference name; ignoring it",
                            bind->regions[i]);
            continue;
        }

        if (tid == HTS_IDX_START) {
            bcf_appender_destroy_contig_tasks(tasks, n_tasks);
            tasks = (bcf_appender_contig_task_t *)calloc(1, sizeof(*tasks));
            n_tasks = tasks ? 1 : 0;
            if (!tasks) {
                snprintf(err, err_size, "read_bcf_appender: out of memory creating full-file job");
                goto cleanup;
            }
            tasks[0].tid = HTS_IDX_START;
            if (!bcf_appender_contig_task_add(&tasks[0], ".")) {
                snprintf(err, err_size, "read_bcf_appender: out of memory creating full-file job");
                goto cleanup;
            }
            ok = 1;
            goto cleanup;
        }

        int task_idx = -1;
        for (int j = 0; j < n_tasks; j++) {
            if (tasks[j].tid == tid) {
                task_idx = j;
                break;
            }
        }
        if (task_idx < 0) {
            task_idx = n_tasks++;
            tasks[task_idx].tid = tid;
        }
        if (!bcf_appender_contig_task_add(&tasks[task_idx], bind->regions[i])) {
            snprintf(err, err_size, "read_bcf_appender: out of memory copying region");
            goto cleanup;
        }
    }
    ok = 1;

cleanup:
    if (idx) hts_idx_destroy(idx);
    if (tbx) tbx_destroy(tbx);
    if (hdr) bcf_hdr_destroy(hdr);
    if (fp) hts_close(fp);
    if (!ok) {
        bcf_appender_destroy_contig_tasks(tasks, n_tasks);
        return 0;
    }
    if (n_tasks == 0) {
        free(tasks);
        tasks = NULL;
    }
    *out_tasks = tasks;
    *out_n_tasks = n_tasks;
    return 1;
}

static void bcf_appender_region_set_error(bcf_appender_region_shared_t *shared, const char *msg) {
    if (!shared || !msg) return;
    if (__sync_bool_compare_and_swap(&shared->failed, 0, 1)) {
        snprintf(shared->error, sizeof(shared->error), "%s", msg);
    }
}

static int bcf_appender_region_failed(const bcf_appender_region_shared_t *shared) {
    return shared && __atomic_load_n(&shared->failed, __ATOMIC_ACQUIRE) != 0;
}

static void *bcf_appender_region_worker_main(void *arg) {
    bcf_appender_region_worker_t *worker = (bcf_appender_region_worker_t *)arg;
    bcf_appender_region_shared_t *shared = worker ? worker->shared : NULL;
    const bcf_appender_bind_t *bind = shared ? shared->bind : NULL;
    htsFile *fp = NULL;
    bcf_hdr_t *hdr = NULL;
    bcf1_t *rec = NULL;
    hts_idx_t *idx = NULL;
    tbx_t *tbx = NULL;
    hts_itr_t *itr = NULL;
    kstring_t line = {0, 0, NULL};
    kstring_t alt = {0, 0, NULL};
    kstring_t gt = {0, 0, NULL};
    int32_t *gt_arr = NULL;
    int gt_arr_n = 0;
    bcf_appender_chunk_t chunk;
    int chunk_initialized = 0;
    char err[1024];

    if (!shared || !bind) return NULL;
    memset(&chunk, 0, sizeof(chunk));

    if (!bcf_appender_chunk_init(&chunk, bind->include_file_offset)) {
        bcf_appender_region_set_error(shared, "read_bcf_appender: failed to create worker data chunk");
        return NULL;
    }
    chunk_initialized = 1;

    fp = hts_open(bind->file_path, "r");
    if (!fp) {
        bcf_appender_region_set_error(shared, "read_bcf_appender: failed to open BCF/VCF file");
        goto cleanup;
    }
    duckhts_apply_remote_hts_tuning(fp, bind->file_path,
                                    DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION);
    if (bind->decompression_threads > 0 &&
        hts_get_format(fp)->compression != no_compression &&
        hts_set_threads(fp, bind->decompression_threads) < 0) {
        bcf_appender_region_set_error(shared, "read_bcf_appender: failed to configure decompression threads");
        goto cleanup;
    }

    hdr = bcf_hdr_read(fp);
    if (!hdr) {
        bcf_appender_region_set_error(shared, "read_bcf_appender: failed to read BCF/VCF header");
        goto cleanup;
    }
    rec = bcf_init();
    if (!rec) {
        bcf_appender_region_set_error(shared, "read_bcf_appender: failed to allocate BCF record");
        goto cleanup;
    }

    enum htsExactFormat fmt = hts_get_format(fp)->format;
    if (fmt == bcf) {
        idx = bcf_index_load3(bind->file_path, NULL, HTS_IDX_SAVE_REMOTE | HTS_IDX_SILENT_FAIL);
    } else {
        tbx = tbx_index_load3(bind->file_path, NULL, HTS_IDX_SAVE_REMOTE | HTS_IDX_SILENT_FAIL);
        if (!tbx) idx = bcf_index_load3(bind->file_path, NULL, HTS_IDX_SAVE_REMOTE | HTS_IDX_SILENT_FAIL);
    }
    if (!idx && !tbx) {
        bcf_appender_region_set_error(shared, "read_bcf_appender: parallel region query requires an index");
        goto cleanup;
    }

    /* Claim jobs in a loop.  Empty indexed contigs consume a job and advance;
     * they never recurse and never acquire another file/index handle. */
    for (;;) {
        int claim = __sync_fetch_and_add(&shared->next_task_idx, 1);
        if (bcf_appender_region_failed(shared)) break;
        if (claim >= shared->n_tasks) break;

        bcf_appender_contig_task_t *task = &shared->tasks[claim];
        itr = bcf_open_region_iterator(idx, hdr, tbx,
                                       task->regions, task->n_regions);
        if (!itr) continue;

        for (;;) {
            uint64_t file_offset = 0;
            int has_file_offset = 0;
            int next_record = bcf_appender_next_record(
                fp, hdr, rec, tbx, itr, &line,
                bind->include_file_offset,
                &file_offset, &has_file_offset,
                err, sizeof(err)
            );
            if (next_record == 0) break;
            if (next_record < 0) {
                bcf_appender_region_set_error(shared, err);
                break;
            }
            if (!bcf_appender_append_record(worker->appender, &chunk,
                                            hdr, rec,
                                            bind->tidy_format,
                                            has_file_offset, file_offset,
                                            bcf_reader_input_is_bcf(fp),
                                            bind->decode_error_policy,
                                            &alt, &gt, &gt_arr, &gt_arr_n,
                                            &worker->rows_written,
                                            err, sizeof(err))) {
                bcf_appender_region_set_error(shared, err);
                break;
            }
            if (bcf_appender_region_failed(shared)) break;
        }

        hts_itr_destroy(itr);
        itr = NULL;
        if (bcf_appender_region_failed(shared)) break;
    }

    if (!bcf_appender_region_failed(shared) &&
        !bcf_appender_chunk_flush(worker->appender, &chunk,
                                  err, sizeof(err))) {
        bcf_appender_region_set_error(shared, err);
    }
    if (!bcf_appender_region_failed(shared) &&
        duckdb_appender_close(worker->appender) == DuckDBError) {
        const char *msg = duckdb_appender_error(worker->appender);
        bcf_appender_region_set_error(shared, msg ? msg : "duckdb_appender_close failed");
    }

cleanup:
    if (itr) hts_itr_destroy(itr);
    if (gt_arr) free(gt_arr);
    ks_free(&line);
    ks_free(&alt);
    ks_free(&gt);
    if (idx) hts_idx_destroy(idx);
    if (tbx) tbx_destroy(tbx);
    if (rec) bcf_destroy(rec);
    if (hdr) bcf_hdr_destroy(hdr);
    if (fp) hts_close(fp);
    if (chunk_initialized) bcf_appender_chunk_destroy(&chunk);
    return NULL;
}

static int bcf_appender_run_parallel_regions(const bcf_appender_bind_t *bind,
                                             uint64_t *rows_written,
                                             char *err,
                                             size_t err_size) {
    bcf_appender_region_shared_t shared;
    pthread_t *threads = NULL;
    bcf_appender_region_worker_t *workers = NULL;
    int n_threads;
    int started = 0;

    *rows_written = 0;
    memset(&shared, 0, sizeof(shared));
    shared.bind = bind;
    if (!bcf_appender_build_contig_tasks(bind, &shared.tasks, &shared.n_tasks,
                                         err, err_size)) {
        return 0;
    }
    if (shared.n_tasks == 0) return 1;

    /* Worker and open-handle count depends on requested contig jobs, not on
     * the number of contig declarations in the file header. */
    n_threads = bind->region_threads;
    if (n_threads > shared.n_tasks) n_threads = shared.n_tasks;
    if (n_threads > BCF_APPENDER_MAX_REGION_THREADS) n_threads = BCF_APPENDER_MAX_REGION_THREADS;
    if (n_threads < 1) n_threads = 1;

    threads = (pthread_t *)malloc(sizeof(pthread_t) * (size_t)n_threads);
    workers = (bcf_appender_region_worker_t *)calloc((size_t)n_threads, sizeof(*workers));
    if (!threads || !workers) {
        snprintf(err, err_size, "read_bcf_appender: out of memory creating region workers");
        free(threads);
        free(workers);
        bcf_appender_destroy_contig_tasks(shared.tasks, shared.n_tasks);
        return 0;
    }

    for (int i = 0; i < n_threads; i++) {
        workers[i].shared = &shared;
        if (duckdb_appender_create(bind->connection, NULL, bind->target_table,
                                   &workers[i].appender) == DuckDBError) {
            bcf_appender_region_set_error(&shared, "read_bcf_appender: failed to create worker appender");
            break;
        }
        if (pthread_create(&threads[i], NULL, bcf_appender_region_worker_main, &workers[i]) != 0) {
            bcf_appender_region_set_error(&shared, "read_bcf_appender: failed to start region worker");
            break;
        }
        started++;
    }
    for (int i = 0; i < started; i++) {
        pthread_join(threads[i], NULL);
        *rows_written += workers[i].rows_written;
    }

    for (int i = 0; i < n_threads; i++) {
        if (workers[i].appender) duckdb_appender_destroy(&workers[i].appender);
    }
    free(threads);
    free(workers);
    bcf_appender_destroy_contig_tasks(shared.tasks, shared.n_tasks);

    if (bcf_appender_region_failed(&shared)) {
        snprintf(err, err_size, "%s", shared.error[0] ? shared.error : "read_bcf_appender: region worker failure");
        return 0;
    }
    return 1;
}

static int bcf_appender_run(const bcf_appender_bind_t *bind, uint64_t *rows_written,
                            double *seconds, char *err, size_t err_size) {
    duckdb_connection con = bind->connection;
    duckdb_appender appender = NULL;
    htsFile *fp = NULL;
    bcf_hdr_t *hdr = NULL;
    bcf1_t *rec = NULL;
    hts_idx_t *idx = NULL;
    tbx_t *tbx = NULL;
    hts_itr_t *itr = NULL;
    kstring_t line = {0, 0, NULL};
    kstring_t alt = {0, 0, NULL};
    kstring_t gt = {0, 0, NULL};
    int32_t *gt_arr = NULL;
    int gt_arr_n = 0;
    bcf_appender_chunk_t chunk;
    int chunk_initialized = 0;
    int transaction_started = 0;
    struct timespec t0, t1;
    int ok = 0;

    memset(&chunk, 0, sizeof(chunk));
    *rows_written = 0;
    *seconds = 0.0;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    if (!con) {
        snprintf(err, err_size, "read_bcf_appender: missing internal DuckDB connection");
        goto cleanup;
    }
    if (!bcf_appender_exec_sql(con, "BEGIN TRANSACTION", err, err_size)) goto cleanup;
    transaction_started = 1;
    if (!bcf_appender_create_target(con, bind, err, err_size)) goto cleanup;

    if (bind->n_regions > 1 && bind->region_threads > 1) {
        if (!bcf_appender_run_parallel_regions(bind, rows_written, err, err_size)) goto cleanup;
        if (transaction_started) {
            if (!bcf_appender_exec_sql(con, "COMMIT", err, err_size)) goto cleanup;
            transaction_started = 0;
        }
        ok = 1;
        goto cleanup;
    }

    if (duckdb_appender_create(con, NULL, bind->target_table, &appender) == DuckDBError) {
        snprintf(err, err_size, "read_bcf_appender: failed to create appender for table %s", bind->target_table);
        goto cleanup;
    }
    if (!bcf_appender_chunk_init(&chunk, bind->include_file_offset)) {
        snprintf(err, err_size, "read_bcf_appender: failed to create data chunk");
        goto cleanup;
    }
    chunk_initialized = 1;

    fp = hts_open(bind->file_path, "r");
    if (!fp) {
        snprintf(err, err_size, "read_bcf_appender: failed to open %s", bind->file_path);
        goto cleanup;
    }
    duckhts_apply_remote_hts_tuning(fp, bind->file_path,
                                    bind->n_regions > 0 ? DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION
                                                        : DUCKHTS_HTS_IO_PROFILE_STREAMING);
    if (bind->decompression_threads > 0 &&
        hts_get_format(fp)->compression != no_compression &&
        hts_set_threads(fp, bind->decompression_threads) < 0) {
        snprintf(err, err_size, "read_bcf_appender: failed to configure decompression threads");
        goto cleanup;
    }

    hdr = bcf_hdr_read(fp);
    if (!hdr) {
        snprintf(err, err_size, "read_bcf_appender: failed to read BCF/VCF header");
        goto cleanup;
    }
    rec = bcf_init();
    if (!rec) {
        snprintf(err, err_size, "read_bcf_appender: failed to allocate BCF record");
        goto cleanup;
    }

    if (bind->n_regions > 0) {
        enum htsExactFormat fmt = hts_get_format(fp)->format;
        if (fmt == bcf) {
            idx = bcf_index_load3(bind->file_path, NULL, HTS_IDX_SAVE_REMOTE | HTS_IDX_SILENT_FAIL);
        } else {
            tbx = tbx_index_load3(bind->file_path, NULL, HTS_IDX_SAVE_REMOTE | HTS_IDX_SILENT_FAIL);
            if (!tbx) idx = bcf_index_load3(bind->file_path, NULL, HTS_IDX_SAVE_REMOTE | HTS_IDX_SILENT_FAIL);
        }
        if (!idx && !tbx) {
            snprintf(err, err_size, "read_bcf_appender: region query requires an index");
            goto cleanup;
        }
    }

    int region_scan_empty = 0;
    if (bind->n_regions > 0) {
        itr = bcf_open_region_iterator(idx, hdr, tbx,
                                       bind->regions, bind->n_regions);
        // No requested region produced an iterator (absent contig or
        // out-of-range coordinates). Skip the read loop, but still finalize and
        // commit so the empty target table is created/replaced rather than
        // rolled back by cleanup with the transaction left open.
        if (!itr) region_scan_empty = 1;
    }

    if (!region_scan_empty) {
        for (;;) {
            uint64_t file_offset = 0;
            int has_file_offset = 0;
            int next_record = bcf_appender_next_record(fp, hdr, rec, tbx, itr, &line,
                                                       bind->include_file_offset,
                                                       &file_offset, &has_file_offset,
                                                       err, err_size);
            if (next_record < 0) goto cleanup;
            if (next_record == 0) {
                break;
            }
            if (!bcf_appender_append_record(appender, &chunk, hdr, rec,
                                            bind->tidy_format,
                                            has_file_offset, file_offset,
                                            bcf_reader_input_is_bcf(fp),
                                            bind->decode_error_policy,
                                            &alt, &gt, &gt_arr, &gt_arr_n,
                                            rows_written, err, err_size)) {
                goto cleanup;
            }
        }
    }

    if (!bcf_appender_chunk_flush(appender, &chunk, err, err_size)) goto cleanup;
    if (duckdb_appender_close(appender) == DuckDBError) {
        const char *msg = duckdb_appender_error(appender);
        snprintf(err, err_size, "%s", msg ? msg : "duckdb_appender_close failed");
        goto cleanup;
    }
    duckdb_appender_destroy(&appender);
    if (transaction_started) {
        if (!bcf_appender_exec_sql(con, "COMMIT", err, err_size)) goto cleanup;
        transaction_started = 0;
    }
    ok = 1;

cleanup:
    clock_gettime(CLOCK_MONOTONIC, &t1);
    *seconds = (double)(t1.tv_sec - t0.tv_sec) + (double)(t1.tv_nsec - t0.tv_nsec) / 1e9;
    if (gt_arr) free(gt_arr);
    ks_free(&line);
    ks_free(&alt);
    ks_free(&gt);
    if (itr) hts_itr_destroy(itr);
    if (idx) hts_idx_destroy(idx);
    if (tbx) tbx_destroy(tbx);
    if (rec) bcf_destroy(rec);
    if (hdr) bcf_hdr_destroy(hdr);
    if (fp) hts_close(fp);
    if (chunk_initialized) bcf_appender_chunk_destroy(&chunk);
    if (appender) duckdb_appender_destroy(&appender);
    if (transaction_started) {
        char rollback_err[256];
        rollback_err[0] = '\0';
        (void)bcf_appender_exec_sql(con, "ROLLBACK", rollback_err, sizeof(rollback_err));
    }
    return ok;
}

static void bcf_appender_bind(duckdb_bind_info info) {
    bcf_appender_extra_t *extra = (bcf_appender_extra_t *)duckdb_bind_get_extra_info(info);
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    duckdb_value table_val = duckdb_bind_get_parameter(info, 1);
    char *file_path = duckdb_get_varchar(path_val);
    char *target_table = duckdb_get_varchar(table_val);
    duckdb_destroy_value(&path_val);
    duckdb_destroy_value(&table_val);

    if (!file_path || !file_path[0] || !target_table || !target_table[0]) {
        duckdb_bind_set_error(info, "read_bcf_appender requires path and target_table");
        if (file_path) duckdb_free(file_path);
        if (target_table) duckdb_free(target_table);
        return;
    }
    if (!bcf_appender_valid_table_name(target_table)) {
        duckdb_bind_set_error(info, "read_bcf_appender target_table must be an unqualified SQL identifier");
        duckdb_free(file_path);
        duckdb_free(target_table);
        return;
    }

    bcf_appender_bind_t *bind = (bcf_appender_bind_t *)duckdb_malloc(sizeof(*bind));
    memset(bind, 0, sizeof(*bind));
    bind->connection = extra ? extra->connection : NULL;
    bind->file_path = file_path;
    bind->target_table = target_table;
    bind->tidy_format = 1;
    bind->overwrite = 1;
    bind->include_file_offset = 0;
    bind->region_threads = 1;
    bind->decompression_threads = 0;
    bind->decode_error_policy = BCF_DECODE_ERROR_NULL;

    duckdb_value region_val = duckdb_bind_get_named_parameter(info, "region");
    if (region_val && !duckdb_is_null_value(region_val)) {
        bind->region = duckdb_get_varchar(region_val);
        parse_regions_duckdb(bind->region, &bind->regions, &bind->n_regions);
    }
    if (region_val) duckdb_destroy_value(&region_val);

    duckdb_value tidy_val = duckdb_bind_get_named_parameter(info, "tidy_format");
    if (tidy_val && !duckdb_is_null_value(tidy_val)) bind->tidy_format = duckdb_get_bool(tidy_val) ? 1 : 0;
    if (tidy_val) duckdb_destroy_value(&tidy_val);

    duckdb_value overwrite_val = duckdb_bind_get_named_parameter(info, "overwrite");
    if (overwrite_val && !duckdb_is_null_value(overwrite_val)) bind->overwrite = duckdb_get_bool(overwrite_val) ? 1 : 0;
    if (overwrite_val) duckdb_destroy_value(&overwrite_val);

    duckdb_value offset_val = duckdb_bind_get_named_parameter(info, "include_file_offset");
    if (offset_val && !duckdb_is_null_value(offset_val)) bind->include_file_offset = duckdb_get_bool(offset_val) ? 1 : 0;
    if (offset_val) duckdb_destroy_value(&offset_val);

    duckdb_value region_threads_val = duckdb_bind_get_named_parameter(info, "region_threads");
    if (region_threads_val && !duckdb_is_null_value(region_threads_val)) {
        int64_t v = duckdb_get_int64(region_threads_val);
        if (v < 1 || v > BCF_APPENDER_MAX_REGION_THREADS) {
            duckdb_bind_set_error(info, "region_threads must be an integer between 1 and 16");
            duckdb_destroy_value(&region_threads_val);
            destroy_bcf_appender_bind(bind);
            return;
        }
        bind->region_threads = (int)v;
    }
    if (region_threads_val) duckdb_destroy_value(&region_threads_val);
    duckdb_value threads_val = duckdb_bind_get_named_parameter(info, "decompression_threads");
    if (threads_val && !duckdb_is_null_value(threads_val)) {
        int64_t v = duckdb_get_int64(threads_val);
        if (v < 0 || v > INT_MAX) {
            duckdb_bind_set_error(info, "decompression_threads must be a non-negative integer");
            duckdb_destroy_value(&threads_val);
            destroy_bcf_appender_bind(bind);
            return;
        }
        bind->decompression_threads = (int)v;
    }
    if (threads_val) duckdb_destroy_value(&threads_val);

    duckdb_value decode_policy_val = duckdb_bind_get_named_parameter(info, "decode_error_policy");
    if (decode_policy_val && !duckdb_is_null_value(decode_policy_val)) {
        char *decode_policy = duckdb_get_varchar(decode_policy_val);
        if (!bcf_parse_decode_error_policy(decode_policy, &bind->decode_error_policy)) {
            duckdb_bind_set_error(info, "read_bcf_appender: decode_error_policy must be 'null', 'warn', or 'error'");
            if (decode_policy) duckdb_free(decode_policy);
            duckdb_destroy_value(&decode_policy_val);
            destroy_bcf_appender_bind(bind);
            return;
        }
        if (decode_policy) duckdb_free(decode_policy);
    }
    if (decode_policy_val) duckdb_destroy_value(&decode_policy_val);

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type ubigint_type = duckdb_create_logical_type(DUCKDB_TYPE_UBIGINT);
    duckdb_logical_type double_type = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_bind_add_result_column(info, "target_table", varchar_type);
    duckdb_bind_add_result_column(info, "rows_written", ubigint_type);
    duckdb_bind_add_result_column(info, "seconds", double_type);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&ubigint_type);
    duckdb_destroy_logical_type(&double_type);

    duckdb_bind_set_bind_data(info, bind, destroy_bcf_appender_bind);
}

static void bcf_appender_init(duckdb_init_info info) {
    bcf_appender_init_t *init = (bcf_appender_init_t *)duckdb_malloc(sizeof(*init));
    memset(init, 0, sizeof(*init));
    duckdb_init_set_max_threads(info, 1);
    duckdb_init_set_init_data(info, init, destroy_bcf_appender_init);
}

static void bcf_appender_function(duckdb_function_info info, duckdb_data_chunk output) {
    bcf_appender_bind_t *bind = (bcf_appender_bind_t *)duckdb_function_get_bind_data(info);
    bcf_appender_init_t *init = (bcf_appender_init_t *)duckdb_function_get_init_data(info);
    if (!bind || !init || init->done) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    uint64_t rows_written = 0;
    double seconds = 0.0;
    char err[1024] = {0};
    if (!bcf_appender_run(bind, &rows_written, &seconds, err, sizeof(err))) {
        duckdb_function_set_error(info, err[0] ? err : "read_bcf_appender failed");
        init->done = 1;
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    duckdb_vector target_vec = duckdb_data_chunk_get_vector(output, 0);
    duckdb_vector rows_vec = duckdb_data_chunk_get_vector(output, 1);
    duckdb_vector seconds_vec = duckdb_data_chunk_get_vector(output, 2);
    uint64_t *rows_data = (uint64_t *)duckdb_vector_get_data(rows_vec);
    double *seconds_data = (double *)duckdb_vector_get_data(seconds_vec);
    duckdb_vector_assign_string_element(target_vec, 0, bind->target_table);
    rows_data[0] = rows_written;
    seconds_data[0] = seconds;
    duckdb_data_chunk_set_size(output, 1);
    init->done = 1;
}

void register_read_bcf_appender_function(duckdb_connection connection, duckdb_database database) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    bcf_appender_extra_t *extra = (bcf_appender_extra_t *)duckdb_malloc(sizeof(*extra));
    if (!extra) return;
    if (duckdb_connect(database, &extra->connection) == DuckDBError || !extra->connection) {
        duckdb_free(extra);
        duckdb_destroy_table_function(&tf);
        duckdb_destroy_logical_type(&varchar_type);
        duckdb_destroy_logical_type(&bool_type);
        duckdb_destroy_logical_type(&bigint_type);
        return;
    }

    duckdb_table_function_set_name(tf, "read_bcf_appender");
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "region", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "tidy_format", bool_type);
    duckdb_table_function_add_named_parameter(tf, "overwrite", bool_type);
    duckdb_table_function_add_named_parameter(tf, "include_file_offset", bool_type);
    duckdb_table_function_add_named_parameter(tf, "region_threads", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "decompression_threads", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "decode_error_policy", varchar_type);
    duckdb_table_function_set_extra_info(tf, extra, destroy_bcf_appender_extra);
    duckdb_table_function_set_bind(tf, bcf_appender_bind);
    duckdb_table_function_set_init(tf, bcf_appender_init);
    duckdb_table_function_set_function(tf, bcf_appender_function);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_logical_type(&bigint_type);
}

// =============================================================================
// Register the bcf_read Table Functions
// =============================================================================

static void register_read_bcf_named_function(duckdb_connection connection, const char *name, int reader_version) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, name);

    // Parameters
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_table_function_add_parameter(tf, varchar_type);  // file_path
    duckdb_table_function_add_named_parameter(tf, "region", varchar_type);  // optional region
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);  // optional explicit index path
    duckdb_table_function_add_named_parameter(tf, "tidy_format", bool_type);  // optional tidy format
    duckdb_table_function_add_named_parameter(tf, "additional_csq_column_types", varchar_type);  // optional bcftools-style CSQ/ANN/BCSQ type overrides
    duckdb_table_function_add_named_parameter(tf, "decompression_threads", bigint_type);  // optional htslib worker threads
    duckdb_table_function_add_named_parameter(tf, "scan_mode", varchar_type);  // 'auto' or 'sequential'
    duckdb_table_function_add_named_parameter(tf, "decode_error_policy", varchar_type);  // 'null', 'warn', or 'error'
    if (reader_version >= 2) {
        duckdb_table_function_add_named_parameter(tf, "samples", varchar_type);  // optional comma-separated sample inclusion list
        duckdb_table_function_add_named_parameter(tf, "samples_file", varchar_type);  // optional sample inclusion file
        duckdb_table_function_add_named_parameter(tf, "exclude_samples", varchar_type);  // optional comma-separated sample exclusion list
        duckdb_table_function_add_named_parameter(tf, "exclude_samples_file", varchar_type);  // optional sample exclusion file
        duckdb_table_function_add_named_parameter(tf, "force_samples", bool_type);  // ignore unknown sample names
        duckdb_table_function_add_named_parameter(tf, "include_info", bool_type);  // bind INFO columns
        duckdb_table_function_add_named_parameter(tf, "include_format", bool_type);  // bind FORMAT/sample columns
        duckdb_table_function_add_named_parameter(tf, "include_vep", bool_type);  // bind VEP/CSQ/ANN/BCSQ columns
        duckdb_table_function_add_named_parameter(tf, "info_fields", varchar_type);  // comma-separated INFO field IDs
        duckdb_table_function_add_named_parameter(tf, "format_fields", varchar_type);  // comma-separated FORMAT field IDs
        duckdb_table_function_add_named_parameter(tf, "vep_fields", varchar_type);  // comma-separated VEP field names
    }
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_logical_type(&bigint_type);

    // Callbacks - use global init + local init for parallel scan support
    duckdb_table_function_set_bind(tf, reader_version >= 2 ? bcf_read_v2_bind : bcf_read_bind);
    duckdb_table_function_set_init(tf, bcf_read_global_init);       // Global init
    duckdb_table_function_set_local_init(tf, bcf_read_local_init);  // Per-thread init
    duckdb_table_function_set_function(tf, bcf_read_function);

    // Enable projection pushdown
    duckdb_table_function_supports_projection_pushdown(tf, true);

    // Register
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}

void register_read_bcf_function(duckdb_connection connection) {
    register_read_bcf_named_function(connection, "read_bcf", 1);
}

void register_read_bcf_v2_function(duckdb_connection connection) {
    register_read_bcf_named_function(connection, "read_bcf_v2", 2);
}
