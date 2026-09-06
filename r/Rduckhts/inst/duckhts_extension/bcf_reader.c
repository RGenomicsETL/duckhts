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

// htslib headers
#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <htslib/hts_log.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>

#include "include/hts_io_tuning.h"
#include "include/duckdb_alloc.h"
#include "include/region_list.h"

// =============================================================================
// Constants
// =============================================================================

#define BCF_READER_DEFAULT_BATCH_SIZE 2048
#define VEP_TRANSCRIPT_ALL 0
#define VEP_TRANSCRIPT_FIRST 1

// Debug/progress tracking
#define BCF_READER_PROGRESS_INTERVAL 100000  // Print progress every N records
#define BCF_READER_ENABLE_PROGRESS 0

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
    int decompression_threads; // htslib decompression worker threads per file handle
    int scan_sequential;       // Force full-file sequential streaming (no index count/parallel scan)
    bcf_decode_error_policy_t decode_error_policy; // null|warn|error for BCF decode mismatches
    char** regions;            // Parsed comma-separated regions
    unsigned int n_regions;
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

    // Worker-owned decode buffers. Loaded values and CSQ descriptors belong to
    // the current record, including when tidy samples span output chunks.
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
    free(bind->regions); /* one region-list allocation, including strings */

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
}

static int bcf_decode_cache_init(bcf_init_data_t *init, const bcf_bind_data_t *bind) {
    if (!init || !bind) return 0;
    init->cache_n_format_fields = bind->n_format_fields;
    init->cache_n_info_fields = bind->n_info_fields;

    if (bind->n_format_fields > 0) {
        size_t nfmt = (size_t)bind->n_format_fields;
        init->fmt_loaded = duckhts_alloc_array(nfmt, sizeof(*init->fmt_loaded));
        init->fmt_ret = duckhts_alloc_array(nfmt, sizeof(*init->fmt_ret));
        init->fmt_n_values = duckhts_alloc_array(nfmt, sizeof(*init->fmt_n_values));
        init->fmt_i32 = duckhts_alloc_array(nfmt, sizeof(*init->fmt_i32));
        init->fmt_f32 = duckhts_alloc_array(nfmt, sizeof(*init->fmt_f32));
        init->fmt_str = duckhts_alloc_array(nfmt, sizeof(*init->fmt_str));
        if (!init->fmt_loaded || !init->fmt_ret || !init->fmt_n_values ||
            !init->fmt_i32 || !init->fmt_f32 || !init->fmt_str) {
            bcf_decode_cache_free(init);
            return 0;
        }
    }

    if (bind->n_info_fields > 0) {
        size_t ninfo = (size_t)bind->n_info_fields;
        init->info_loaded = duckhts_alloc_array(ninfo, sizeof(*init->info_loaded));
        init->info_ret = duckhts_alloc_array(ninfo, sizeof(*init->info_ret));
        init->info_n_values = duckhts_alloc_array(ninfo, sizeof(*init->info_n_values));
        init->info_flag = duckhts_alloc_array(ninfo, sizeof(*init->info_flag));
        init->info_i32 = duckhts_alloc_array(ninfo, sizeof(*init->info_i32));
        init->info_f32 = duckhts_alloc_array(ninfo, sizeof(*init->info_f32));
        init->info_str = duckhts_alloc_array(ninfo, sizeof(*init->info_str));
        if (!init->info_loaded || !init->info_ret || !init->info_n_values ||
            !init->info_flag || !init->info_i32 || !init->info_f32 || !init->info_str) {
            bcf_decode_cache_free(init);
            return 0;
        }
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

static int bcf_region_name2id(void *hdr, const char *name) {
    return bcf_hdr_name2id((bcf_hdr_t *)hdr, name);
}

static int bcf_tabix_region_name2id(void *tbx, const char *name) {
    return tbx_name2id((tbx_t *)tbx, name);
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

static void free_bcf_format_string_array(char **values) {
    if (!values) return;
    if (values[0]) free(values[0]);
    free(values);
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

    const char *reader_name = "read_bcf";

    local->projected_cols = duckhts_alloc_array(local->column_count, sizeof(*local->projected_cols));
    if (!local->projected_cols) {
        snprintf(err, err_size, "%s: out of memory allocating projected column descriptors", reader_name);
        return 0;
    }

    int tidy_mode = bind->tidy_format && bind->n_samples > 0;

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
            int schema_field_idx = field_idx;
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

    int *format_counts = NULL;
    if (bind->n_format_fields > 0) {
        format_counts = duckhts_alloc_array(bind->n_format_fields, sizeof(*format_counts));
        if (!format_counts) {
            snprintf(err, err_size, "%s: out of memory allocating FORMAT projection counts", reader_name);
            return 0;
        }
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
        }
    }

    if (local->site_column_count > 0) {
        local->site_column_indices = duckhts_alloc_array(local->site_column_count, sizeof(*local->site_column_indices));
        if (!local->site_column_indices) {
            snprintf(err, err_size, "%s: out of memory allocating site column descriptors", reader_name);
            if (format_counts) duckdb_free(format_counts);
            return 0;
        }
    }

    int *format_group_for_field = NULL;
    int *format_group_offsets = NULL;
    if (local->n_format_groups > 0) {
        local->format_groups = duckhts_alloc_array(local->n_format_groups, sizeof(*local->format_groups));
        format_group_for_field = duckhts_alloc_array(bind->n_format_fields, sizeof(*format_group_for_field));
        format_group_offsets = duckhts_alloc_array(local->n_format_groups, sizeof(*format_group_offsets));
        if (!local->format_groups || !format_group_for_field || !format_group_offsets) {
            snprintf(err, err_size, "%s: out of memory allocating FORMAT projection groups", reader_name);
            if (format_counts) duckdb_free(format_counts);
            if (format_group_for_field) duckdb_free(format_group_for_field);
            if (format_group_offsets) duckdb_free(format_group_offsets);
            return 0;
        }
        for (int i = 0; i < bind->n_format_fields; i++) {
            format_group_for_field[i] = -1;
        }

        int group_idx = 0;
        for (int field_idx = 0; field_idx < bind->n_format_fields; field_idx++) {
            if (format_counts[field_idx] == 0) {
                continue;
            }
            bcf_format_group_t *group = &local->format_groups[group_idx];
            group->field_idx = field_idx;
            group->kind = bcf_projected_format_kind(&bind->format_fields[field_idx]);
            group->column_count = (idx_t)format_counts[field_idx];
            group->projected_indices = duckhts_alloc_array(group->column_count, sizeof(*group->projected_indices));
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

static void bcf_read_bind(duckdb_bind_info info) {
    const char *reader_name = "read_bcf";

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
        if (!region) {
            duckdb_destroy_value(&region_val);
            duckdb_free(file_path);
            duckdb_bind_set_error(info, "region list: out of memory reading region");
            return;
        }
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
        return;
    }


    duckdb_logical_type varchar_type = NULL, bigint_type = NULL;
    duckdb_logical_type double_type = NULL, varchar_list_type = NULL;
    bcf_bind_data_t* bind = duckhts_alloc_array(1, sizeof(*bind));
    if (!bind) goto bind_oom;
    bind->file_path = file_path;
    bind->index_path = index_path;
    bind->region = region;
    bind->additional_csq_column_types = additional_csq_column_types;
    bind->decompression_threads = decompression_threads;
    bind->scan_sequential = scan_sequential;
    bind->decode_error_policy = decode_error_policy;
    char region_error[256];
    if (!duckhts_region_list_parse(region, &bind->regions, &bind->n_regions,
                                   region_error, sizeof(region_error))) {
        duckdb_bind_set_error(info, region_error);
        goto bind_error;
    }
    bind->n_samples = bcf_hdr_nsamples(hdr);
    bind->tidy_format = tidy_format;
    bind->sample_id_col_idx = -1;  // Will be set if tidy_format=true
    bind->n_vep_fields = 0;
    bind->vep_col_start = COL_CORE_COUNT;
    bind->info_col_start = COL_CORE_COUNT;
    bind->format_col_start = COL_CORE_COUNT;
    bind->vep_schema = NULL;
    bind->vep_transcript_mode = VEP_TRANSCRIPT_FIRST;

    // Copy sample names
    if (bind->n_samples > 0) {
        bind->sample_names = duckhts_alloc_array(bind->n_samples, sizeof(*bind->sample_names));
        if (!bind->sample_names) goto bind_oom;
        for (int i = 0; i < bind->n_samples; i++) {
            bind->sample_names[i] = duckhts_copy_string(hdr->samples[i]);
            if (!bind->sample_names[i]) goto bind_oom;
        }
    }

    // Create logical types for schema
    varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    double_type = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    varchar_list_type = duckdb_create_list_type(varchar_type);

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
    bind->vep_schema = vep_schema_parse(hdr, NULL, bind->additional_csq_column_types);
    if (bind->vep_schema) {
        bind->n_vep_fields = bind->vep_schema->n_fields;
        if (bind->n_vep_fields > INT_MAX - col_idx) goto bind_too_many_columns;
        bind->vep_col_start = col_idx;
        for (int v = 0; v < bind->vep_schema->n_fields; v++) {
            const vep_field_t* field = vep_schema_get_field(bind->vep_schema, v);
            if (!field) {
                continue;
            }

            char col_name[256];
            snprintf(col_name, sizeof(col_name), "VEP_%s", field->name);

            // Expose all transcripts as list columns for full preservation
            duckdb_logical_type field_type = create_vep_field_type(field->type, 1);
            duckdb_bind_add_result_column(info, col_name, field_type);
            duckdb_destroy_logical_type(&field_type);

            col_idx++;
        }

        bind->vep_transcript_mode = VEP_TRANSCRIPT_ALL;
    }
    // -------------------------------------------------------------------------
    // INFO fields (with type validation)
    // -------------------------------------------------------------------------
    bind->info_col_start = col_idx;

    // Count INFO fields
    bind->n_info_fields = 0;
    for (int i = 0; i < hdr->n[BCF_DT_ID]; i++) {
        if (hdr->id[BCF_DT_ID][i].val &&
            hdr->id[BCF_DT_ID][i].val->hrec[BCF_HL_INFO]) {
            bind->n_info_fields++;
        }
    }

    if (bind->n_info_fields > INT_MAX - col_idx) goto bind_too_many_columns;
    if (bind->n_info_fields > 0) {
        bind->info_fields = duckhts_alloc_array(bind->n_info_fields, sizeof(*bind->info_fields));
        if (!bind->info_fields) goto bind_oom;

        int info_idx = 0;
        for (int i = 0; i < hdr->n[BCF_DT_ID] && info_idx < bind->n_info_fields; i++) {
            if (hdr->id[BCF_DT_ID][i].val &&
                hdr->id[BCF_DT_ID][i].val->hrec[BCF_HL_INFO]) {
                const char* field_name = hdr->id[BCF_DT_ID][i].key;
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
                field->name = duckhts_copy_string(field_name);
                if (!field->name) goto bind_oom;
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

    if (bind->n_samples > 0) {
        // Count FORMAT fields
        bind->n_format_fields = 0;
        int header_format_count = 0;
        for (int i = 0; i < hdr->n[BCF_DT_ID]; i++) {
            if (hdr->id[BCF_DT_ID][i].val &&
                hdr->id[BCF_DT_ID][i].val->hrec[BCF_HL_FMT]) {
                header_format_count++;
                bind->n_format_fields++;
            }
        }
        if (bind->n_format_fields == 0 && header_format_count == 0) {
            bind->n_format_fields = 1;
        }
        uint64_t format_columns = bind->tidy_format
            ? (uint64_t)bind->n_format_fields + 1
            : (uint64_t)bind->n_format_fields * (uint64_t)bind->n_samples;
        if (format_columns > (uint64_t)(INT_MAX - col_idx)) goto bind_too_many_columns;

        if (bind->n_format_fields == 1 && header_format_count == 0) {
            // Add GT as default for old header-light fixtures.
            bind->format_fields = duckhts_alloc_array(1, sizeof(*bind->format_fields));
            if (!bind->format_fields) goto bind_oom;
            bind->format_fields[0].name = duckhts_copy_string("GT");
            if (!bind->format_fields[0].name) goto bind_oom;
            bind->format_fields[0].header_type = BCF_HT_STR;
            bind->format_fields[0].schema_type = BCF_HT_STR;
            bind->format_fields[0].vl_type = BCF_VL_FIXED;
            bind->format_fields[0].fixed_count = 1;
            bind->format_fields[0].is_list = 0;
        }
        if (!bind->format_fields && bind->n_format_fields > 0) {
            bind->format_fields = duckhts_alloc_array(bind->n_format_fields, sizeof(*bind->format_fields));
            if (!bind->format_fields) goto bind_oom;

            int fmt_idx = 0;
            for (int i = 0; i < hdr->n[BCF_DT_ID] && fmt_idx < bind->n_format_fields; i++) {
                if (hdr->id[BCF_DT_ID][i].val &&
                    hdr->id[BCF_DT_ID][i].val->hrec[BCF_HL_FMT]) {
                    const char* field_name = hdr->id[BCF_DT_ID][i].key;
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
                    field->name = duckhts_copy_string(field_name);
                    if (!field->name) goto bind_oom;
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
            int tidy_rows = bind->tidy_format && bind->n_samples > 0;
            uint64_t row_multiplier = tidy_rows ? (uint64_t)bind->n_samples : 1;
            bind->index_row_count_valid = bcf_try_get_index_row_count(idx ? idx : tbx->idx,
                                                                      row_multiplier,
                                                                      &bind->index_row_count);
            if (idx) hts_idx_destroy(idx);
            if (tbx) tbx_destroy(tbx);

            // Get contig names for automatic parallel contig scans.
            int n_seqs = hdr->n[BCF_DT_CTG];
            if (n_seqs > 0) {
                bind->n_contigs = n_seqs;
                bind->contig_names = duckhts_alloc_array(n_seqs, sizeof(*bind->contig_names));
                if (!bind->contig_names) goto bind_oom;

                for (int i = 0; i < n_seqs; i++) {
                    bind->contig_names[i] = duckhts_copy_string(hdr->id[BCF_DT_CTG][i].key);
                    if (!bind->contig_names[i]) goto bind_oom;
                }
            }
        }
    }

    goto bind_cleanup;

bind_too_many_columns:
    duckdb_bind_set_error(info, "read_bcf: header exceeds the supported column count");
    goto bind_error;
bind_oom:
    duckdb_bind_set_error(info, "read_bcf: out of memory allocating bind schema");
bind_error:
    if (bind) {
        destroy_bind_data(bind);
        bind = NULL;
    } else {
        duckdb_free(file_path);
        if (index_path) duckdb_free(index_path);
        if (region) duckdb_free(region);
        if (additional_csq_column_types) duckdb_free(additional_csq_column_types);
    }
bind_cleanup:
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&double_type);
    duckdb_destroy_logical_type(&varchar_list_type);

    bcf_hdr_destroy(hdr);
    hts_close(fp);

    if (bind) duckdb_bind_set_bind_data(info, bind, destroy_bind_data);
}

// =============================================================================
// Global Init Function - Set up parallel scanning
// =============================================================================

static void bcf_read_global_init(duckdb_init_info info) {
    bcf_bind_data_t* bind = (bcf_bind_data_t*)duckdb_init_get_bind_data(info);
    idx_t column_count = duckdb_init_get_column_count(info);

    bcf_global_init_data_t* global = duckhts_alloc_array(1, sizeof(*global));
    if (!global) {
        duckdb_init_set_error(info, "read_bcf: out of memory allocating global state");
        return;
    }

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
    } else if (!bind->scan_sequential && bind->has_index && bind->n_contigs > 1 && !global->has_region) {
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
    const char *reader_name = "read_bcf";

    bcf_init_data_t* local = duckhts_alloc_array(1, sizeof(*local));
    if (!local) {
        duckdb_init_set_error(info, "read_bcf: out of memory allocating worker state");
        return;
    }

    local->column_count = duckdb_init_get_column_count(info);
    if (local->column_count > 0) {
        local->column_ids = duckhts_alloc_array(local->column_count, sizeof(*local->column_ids));
        if (!local->column_ids) {
            duckdb_init_set_error(info, "read_bcf: out of memory allocating projected columns");
            destroy_init_data(local);
            return;
        }
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

    if (!bcf_decode_cache_init(local, bind)) {
        char err[128];
        snprintf(err, sizeof(err), "%s: out of memory allocating decode cache", reader_name);
        duckdb_init_set_error(info, err);
        destroy_init_data(local);
        return;
    }

    // Check if we're in parallel mode based on bind data
    int is_parallel = (!bind->scan_sequential && bind->has_index && bind->n_contigs > 1 && bind->n_regions == 0);

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

    // Allocate record
    local->rec = bcf_init();
    if (!local->rec) {
        char err[128];
        snprintf(err, sizeof(err), "%s: failed to allocate BCF/VCF record", reader_name);
        duckdb_init_set_error(info, err);
        destroy_init_data(local);
        return;
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

    /* Bind already chose one-owner-per-contig scanning. Without the worker's
     * index, every contig would look empty; streaming fallback would duplicate
     * the file across workers. Fail this plan instead. */
    if (is_parallel && !local->idx && !local->tbx) {
        duckdb_init_set_error(info, "read_bcf: failed to reload index for parallel scan");
        destroy_init_data(local);
        return;
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

        /* VCF indexes can name contigs absent from the text header, especially
         * with an explicit non-colocated index. Use the iterator's dictionary. */
        char region_error[256];
        if (!duckhts_region_list_validate(bind->regions, bind->n_regions,
                                          local->idx ? bcf_region_name2id : bcf_tabix_region_name2id,
                                          local->idx ? (void *)local->hdr : (void *)local->tbx,
                                          region_error, sizeof(region_error))) {
            duckdb_init_set_error(info, region_error);
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
    const char *reader_name = "read_bcf";
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
    int sample_idx = (bind->tidy_format && bind->n_samples > 0)
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
        if (init->column_count <= SIZE_MAX / sizeof(*vectors)) {
            vectors = duckdb_malloc((size_t)init->column_count * sizeof(*vectors));
        }
        if (!vectors) {
            duckdb_function_set_error(info, "read_bcf: out of memory allocating output vector pointers");
            duckdb_data_chunk_set_size(output, 0);
            return;
        }
        for (idx_t i = 0; i < init->column_count; i++) {
            vectors[i] = duckdb_data_chunk_get_vector(output, i);
        }
    }

    // Tidy format variables
    int tidy_mode = bind->tidy_format && bind->n_samples > 0;
    int current_sample = 0;  // Which sample we're emitting (only used in tidy mode)

    // Borrow the worker cache; only the record transition resets loaded state.
    int *fmt_loaded = init->fmt_loaded;
    int *fmt_ret = init->fmt_ret;
    int *fmt_n_values = init->fmt_n_values;
    int32_t **fmt_i32 = init->fmt_i32;
    float **fmt_f32 = init->fmt_f32;
    char ***fmt_str = init->fmt_str;
    int gt_loaded = init->gt_loaded;
    int gt_ret = init->gt_ret;
    int gt_n_values = init->gt_n_values;
    int32_t *gt_arr = init->gt_arr;
    int *info_loaded = init->info_loaded;
    int *info_ret = init->info_ret;
    int *info_n_values = init->info_n_values;
    int *info_flag = init->info_flag;
    int32_t **info_i32 = init->info_i32;
    float **info_f32 = init->info_f32;
    char **info_str = init->info_str;

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
                    // VCF with tabix: read text line then parse.
                    ret = tbx_itr_next(init->fp, init->tbx, init->itr, &init->kstr);
                    if (ret >= 0) {
                        ret = vcf_parse1(&init->kstr, init->hdr, init->rec);
                        init->kstr.l = 0;
                    }
                } else {
                    // BCF with index
                    ret = bcf_itr_next(init->fp, init->itr, init->rec);
                }
            } else {
                ret = bcf_read(init->fp, init->hdr, init->rec);
            }

            if (ret < -1) {
                char err[256];
                const char *reader_name = "read_bcf";
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
            if (init->vep_rec) {
                vep_record_destroy(init->vep_rec);
                init->vep_rec = NULL;
            }
            init->vep_loaded = 0;
        }

        // Retain parsed annotation until every sample of this record is emitted.
        vep_record_t* vep_rec = NULL;
        if (init->need_vep) {
            if (!init->vep_loaded) {
                init->vep_rec = vep_record_parse_bcf(bind->vep_schema, init->hdr, init->rec);
                init->vep_loaded = 1;
            }
            vep_rec = init->vep_rec;
        }

        // Decode each projected FORMAT field once per input record; later
        // tidy samples borrow the loaded values.
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

        idx_t row_idx = row_count;
        int sample_for_row = current_sample;

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
                                                                   "read_bcf",
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
                    if (!bcf_check_decode_ret("read_bcf", "INFO", tag, init->hdr, init->rec,
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
                                                                   "read_bcf",
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
                    if (!bcf_check_decode_ret("read_bcf", "INFO", tag, init->hdr, init->rec,
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
                                                                   "read_bcf",
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
                    if (!bcf_check_decode_ret("read_bcf", "INFO", tag, init->hdr, init->rec,
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

        if (scan_error) {
            break;
        }

        row_count++;
        init->current_row++;

        // Advance the tidy sample, or mark the record as consumed.
        if (tidy_mode) {
            init->tidy_current_sample++;
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

    duckdb_free(vectors);
    init->gt_loaded = gt_loaded;
    init->gt_ret = gt_ret;
    init->gt_n_values = gt_n_values;
    init->gt_arr = gt_arr;

    duckdb_data_chunk_set_size(output, row_count);
}


// =============================================================================
// Register the bcf_read Table Functions
// =============================================================================

void register_read_bcf_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "read_bcf");

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
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_logical_type(&bigint_type);

    // Callbacks - use global init + local init for parallel scan support
    duckdb_table_function_set_bind(tf, bcf_read_bind);
    duckdb_table_function_set_init(tf, bcf_read_global_init);       // Global init
    duckdb_table_function_set_local_init(tf, bcf_read_local_init);  // Per-thread init
    duckdb_table_function_set_function(tf, bcf_read_function);

    // Enable projection pushdown
    duckdb_table_function_supports_projection_pushdown(tf, true);

    // Register
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}
