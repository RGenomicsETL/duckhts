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

// htslib headers
#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <htslib/hts_log.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>

// =============================================================================
// Constants
// =============================================================================

#define BCF_READER_DEFAULT_BATCH_SIZE 2048
#define VEP_TRANSCRIPT_ALL 0
#define VEP_TRANSCRIPT_FIRST 1

// Debug/progress tracking
#define BCF_READER_PROGRESS_INTERVAL 100000  // Print progress every N records
#define BCF_READER_ENABLE_PROGRESS 0

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
    int is_list;             // Whether this is a list type
    int duckdb_col_idx;      // Column index in DuckDB result
} field_meta_t;

// =============================================================================
// Bind Data - stores parameters and schema info
// =============================================================================

typedef struct {
    char* file_path;
    char* index_path;          // Optional explicit index path
    char* region;              // Optional region filter
    char** regions;            // Parsed comma-separated regions
    unsigned int n_regions;
    int include_info;          // Include INFO fields
    int include_format;        // Include FORMAT/sample fields
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
    
    int64_t current_row;
    int done;
    
    // Projection pushdown
    idx_t column_count;
    idx_t* column_ids;
    
    // Parallel scan state
    int is_parallel;
    int assigned_contig;       // Which contig this thread is scanning (-1 = all)
    const char* contig_name;   // Name of assigned contig (reference, don't free)
    int needs_next_contig;     // Flag to request next contig assignment
    unsigned int next_region_idx; // Region cursor for chained region scans
    
    // Tidy format state: tracks which sample we're emitting for current record
    int tidy_current_sample;   // Current sample index in tidy mode (-1 = need to read next record)
    int tidy_record_valid;     // Whether we have a valid record buffered for tidy mode
    
    // Debug/progress tracking
    int64_t total_records_processed;  // Total records processed by this thread
    struct timespec batch_start_time;  // Start time for performance measurement
    struct timespec last_progress_time;  // Last time progress was logged
    int timing_initialized;           // Flag to indicate timing is set up
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

// =============================================================================
// Memory Management
// =============================================================================

static void destroy_bind_data(void* data) {
    bcf_bind_data_t* bind = (bcf_bind_data_t*)data;
    if (!bind) return;
    
    if (bind->file_path) duckdb_free(bind->file_path);
    if (bind->index_path) duckdb_free(bind->index_path);
    if (bind->region) duckdb_free(bind->region);
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

static void destroy_init_data(void* data) {
    bcf_init_data_t* init = (bcf_init_data_t*)data;
    if (!init) return;
    
    if (init->itr) hts_itr_destroy(init->itr);
    if (init->tbx) tbx_destroy(init->tbx);
    if (init->idx) hts_idx_destroy(init->idx);
    if (init->rec) bcf_destroy(init->rec);
    if (init->hdr) bcf_hdr_destroy(init->hdr);
    if (init->fp) hts_close(init->fp);
    if (init->column_ids) duckdb_free(init->column_ids);
    ks_free(&init->kstr);
    
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

// =============================================================================
// Schema Building - Bind Function
// =============================================================================

static void bcf_read_bind(duckdb_bind_info info) {
    // Set up warning callback
    vcf_set_warning_callback(duckdb_vcf_warning, NULL);
    
    // Get the file path parameter
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    char* file_path = duckdb_get_varchar(path_val);
    duckdb_destroy_value(&path_val);
    
    if (!file_path || strlen(file_path) == 0) {
        duckdb_bind_set_error(info, "read_bcf requires a file path");
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
    
    // Get optional tidy_format named parameter (default: false)
    int tidy_format = 0;
    duckdb_value tidy_val = duckdb_bind_get_named_parameter(info, "tidy_format");
    if (tidy_val && !duckdb_is_null_value(tidy_val)) {
        tidy_format = duckdb_get_bool(tidy_val);
    }
    if (tidy_val) duckdb_destroy_value(&tidy_val);
    
    // Open the file to read header
    htsFile* fp = hts_open(file_path, "r");
    if (!fp) {
        char err[512];
        snprintf(err, sizeof(err), "Failed to open BCF/VCF file: %s", file_path);
        duckdb_bind_set_error(info, err);
        duckdb_free(file_path);
        if (index_path) duckdb_free(index_path);
        if (region) duckdb_free(region);
        return;
    }
    
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        duckdb_bind_set_error(info, "Failed to read BCF/VCF header");
        duckdb_free(file_path);
        if (index_path) duckdb_free(index_path);
        if (region) duckdb_free(region);
        return;
    }
    
    // Create bind data
    bcf_bind_data_t* bind = (bcf_bind_data_t*)duckdb_malloc(sizeof(bcf_bind_data_t));
    memset(bind, 0, sizeof(bcf_bind_data_t));
    bind->file_path = file_path;
    bind->index_path = index_path;
    bind->region = region;
    parse_regions_duckdb(region, &bind->regions, &bind->n_regions);
    bind->include_info = 1;
    bind->include_format = 1;
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
    bind->vep_schema = vep_schema_parse(hdr, NULL);
    if (bind->vep_schema) {
        bind->n_vep_fields = bind->vep_schema->n_fields;
        bind->vep_col_start = col_idx;
        
        for (int v = 0; v < bind->n_vep_fields; v++) {
            const vep_field_t* field = vep_schema_get_field(bind->vep_schema, v);
            if (!field) continue;
            
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
    
    if (bind->n_info_fields > 0) {
        bind->info_fields = (field_meta_t*)duckdb_malloc(bind->n_info_fields * sizeof(field_meta_t));
        memset(bind->info_fields, 0, bind->n_info_fields * sizeof(field_meta_t));
        
        int info_idx = 0;
        for (int i = 0; i < hdr->n[BCF_DT_ID] && info_idx < bind->n_info_fields; i++) {
            if (hdr->id[BCF_DT_ID][i].val && 
                hdr->id[BCF_DT_ID][i].val->hrec[BCF_HL_INFO]) {
                const char* field_name = hdr->id[BCF_DT_ID][i].key;
                int header_type = bcf_hdr_id2type(hdr, BCF_HL_INFO, i);
                int header_vl_type = bcf_hdr_id2length(hdr, BCF_HL_INFO, i);
                
                // Validate against VCF spec (emits warnings)
                int corrected_type;
                int corrected_vl_type = vcf_validate_info_field(field_name, header_vl_type, 
                                                                 header_type, &corrected_type);
                
                field_meta_t* field = &bind->info_fields[info_idx];
                field->name = strdup_duckdb(field_name);
                field->header_id = i;
                field->header_type = header_type;
                field->schema_type = header_type;  // Use header type for data
                field->vl_type = corrected_vl_type;
                field->is_list = vcf_is_list_type(corrected_vl_type);
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
        for (int i = 0; i < hdr->n[BCF_DT_ID]; i++) {
            if (hdr->id[BCF_DT_ID][i].val && 
                hdr->id[BCF_DT_ID][i].val->hrec[BCF_HL_FMT]) {
                bind->n_format_fields++;
            }
        }
        
        if (bind->n_format_fields == 0) {
            // Add GT as default
            bind->n_format_fields = 1;
            bind->format_fields = (field_meta_t*)duckdb_malloc(sizeof(field_meta_t));
            memset(bind->format_fields, 0, sizeof(field_meta_t));
            bind->format_fields[0].name = strdup_duckdb("GT");
            bind->format_fields[0].header_type = BCF_HT_STR;
            bind->format_fields[0].schema_type = BCF_HT_STR;
            bind->format_fields[0].vl_type = BCF_VL_FIXED;
            bind->format_fields[0].is_list = 0;
        } else {
            bind->format_fields = (field_meta_t*)duckdb_malloc(bind->n_format_fields * sizeof(field_meta_t));
            memset(bind->format_fields, 0, bind->n_format_fields * sizeof(field_meta_t));
            
            int fmt_idx = 0;
            for (int i = 0; i < hdr->n[BCF_DT_ID] && fmt_idx < bind->n_format_fields; i++) {
                if (hdr->id[BCF_DT_ID][i].val && 
                    hdr->id[BCF_DT_ID][i].val->hrec[BCF_HL_FMT]) {
                    const char* field_name = hdr->id[BCF_DT_ID][i].key;
                    int header_type = bcf_hdr_id2type(hdr, BCF_HL_FMT, i);
                    int header_vl_type = bcf_hdr_id2length(hdr, BCF_HL_FMT, i);
                    
                    // Validate against VCF spec (emits warnings, only once)
                    int corrected_type;
                    int corrected_vl_type = vcf_validate_format_field(field_name, header_vl_type,
                                                                       header_type, &corrected_type);
                    
                    field_meta_t* field = &bind->format_fields[fmt_idx];
                    field->name = strdup_duckdb(field_name);
                    field->header_id = i;
                    field->header_type = header_type;
                    field->schema_type = header_type;
                    field->vl_type = corrected_vl_type;
                    field->is_list = vcf_is_list_type(corrected_vl_type);
                    
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
    if (bind->n_regions == 0) {
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
            
            // Get contig names from header for parallel scan
            int n_seqs = hdr->n[BCF_DT_CTG];
            if (n_seqs > 0) {
                bind->n_contigs = n_seqs;
                bind->contig_names = (char**)duckdb_malloc(n_seqs * sizeof(char*));
                
                for (int i = 0; i < n_seqs; i++) {
                    bind->contig_names[i] = strdup_duckdb(hdr->id[BCF_DT_CTG][i].key);
                }
            }
            
            if (idx) hts_idx_destroy(idx);
            if (tbx) tbx_destroy(tbx);
        }
    }
    
    // Cleanup
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&double_type);
    duckdb_destroy_logical_type(&varchar_list_type);
    
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    
    duckdb_bind_set_bind_data(info, bind, destroy_bind_data);
}

// =============================================================================
// Global Init Function - Set up parallel scanning
// =============================================================================

static void bcf_read_global_init(duckdb_init_info info) {
    bcf_bind_data_t* bind = (bcf_bind_data_t*)duckdb_init_get_bind_data(info);
    
    bcf_global_init_data_t* global = (bcf_global_init_data_t*)duckdb_malloc(sizeof(bcf_global_init_data_t));
    memset(global, 0, sizeof(bcf_global_init_data_t));
    
    global->current_contig = 0;
    global->has_region = (bind->n_regions > 0);
    
    // Enable parallel scan if:
    // 1. Index exists
    // 2. Multiple contigs available
    // 3. No user-specified region (region queries are already filtered)
    if (bind->has_index && bind->n_contigs > 1 && !global->has_region) {
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
    
    bcf_init_data_t* local = (bcf_init_data_t*)duckdb_malloc(sizeof(bcf_init_data_t));
    memset(local, 0, sizeof(bcf_init_data_t));
    
    // Check if we're in parallel mode based on bind data
    int is_parallel = (bind->has_index && bind->n_contigs > 1 && bind->n_regions == 0);
    
    // Initialize parallel scan state
    local->is_parallel = is_parallel;
    local->assigned_contig = -1;
    local->contig_name = NULL;
    local->needs_next_contig = is_parallel;  // Start by requesting first contig
    
    // Open file (each thread gets its own file handle)
    local->fp = hts_open(bind->file_path, "r");
    if (!local->fp) {
        duckdb_init_set_error(info, "Failed to open BCF/VCF file");
        duckdb_free(local);
        return;
    }
    
    // Read header
    local->hdr = bcf_hdr_read(local->fp);
    if (!local->hdr) {
        hts_close(local->fp);
        duckdb_init_set_error(info, "Failed to read BCF/VCF header");
        duckdb_free(local);
        return;
    }
    
    // Allocate record
    local->rec = bcf_init();
    
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
                     "Region query requires an index file (.tbi or .csi). Region: %s", bind->region);
            duckdb_init_set_error(info, err);
            destroy_init_data(local);
            return;
        }

        if (bind->n_regions > 1) {
            vcf_emit_warning("Multi-region BCF/VCF queries are executed as a chained union of single-region iterators; overlapping regions may return duplicate rows");
        }
        
        local->next_region_idx = 0;
        while (local->next_region_idx < bind->n_regions) {
            const char *query_region = bind->regions[local->next_region_idx++];
            if (local->idx) {
                local->itr = bcf_itr_querys(local->idx, local->hdr, query_region);
            } else if (local->tbx) {
                local->itr = tbx_itr_querys(local->tbx, query_region);
            }
            if (local->itr) break;

            // Avoid hard failure for non-conforming headers.
            char msg[512];
            snprintf(msg, sizeof(msg),
                     "Region query returned no iterator; skipping region: %s",
                     query_region);
            vcf_emit_warning(msg);
        }
        
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
    
    // Get projection pushdown info
    local->column_count = duckdb_init_get_column_count(info);
    local->column_ids = (idx_t*)duckdb_malloc(sizeof(idx_t) * local->column_count);
    for (idx_t i = 0; i < local->column_count; i++) {
        local->column_ids[i] = duckdb_init_get_column_index(info, i);
    }
    
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
    
    duckdb_vector child_vec = duckdb_list_vector_get_child(vec);
    
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
    
    // Atomically claim next contig using fetch-and-add
    // This prevents race conditions where two threads grab the same contig
    int next = __sync_fetch_and_add(&global->current_contig, 1);
    if (next >= global->n_contigs) {
        return 0;  // No more contigs
    }
    
    // Destroy old iterator if exists
    if (init->itr) {
        hts_itr_destroy(init->itr);
        init->itr = NULL;
    }
    
    // Set up iterator for this contig
    const char* contig = global->contig_names[next];
    init->assigned_contig = next;
    init->contig_name = contig;
    
    if (init->idx) {
        init->itr = bcf_itr_querys(init->idx, init->hdr, contig);
    } else if (init->tbx) {
        init->itr = tbx_itr_querys(init->tbx, contig);
    }
    
    if (!init->itr) {
        // This contig might not have any records - try next
        return claim_next_contig(init, global);
    }
    
    init->needs_next_contig = 0;
    return 1;
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

    // Determine if any VEP columns are requested for this scan
    int need_vep = (bind->vep_schema != NULL);
    if (need_vep) {
        need_vep = 0;
        for (idx_t i = 0; i < init->column_count; i++) {
            idx_t col_id = init->column_ids[i];
            if (col_id >= (idx_t)bind->vep_col_start &&
                col_id < (idx_t)(bind->vep_col_start + bind->n_vep_fields)) {
                need_vep = 1;
                break;
            }
        }
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
    duckdb_vector* vectors = (duckdb_vector*)duckdb_malloc(init->column_count * sizeof(duckdb_vector));
    
    for (idx_t i = 0; i < init->column_count; i++) {
        vectors[i] = duckdb_data_chunk_get_vector(output, i);
    }
    
    // Tidy format variables
    int tidy_mode = bind->tidy_format && bind->n_samples > 0;
    int current_sample = 0;  // Which sample we're emitting (only used in tidy mode)
    
    // Per-record FORMAT cache to avoid repeated bcf_get_format_* calls across
    // projected columns/samples for the same record.
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
    if (bind->n_format_fields > 0) {
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
    if (bind->n_info_fields > 0) {
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
                    // VCF with tabix: read text line then parse
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
            
            if (ret < 0) {
                // End of current contig/file
                if (init->is_parallel) {
                    // Try to claim next contig
                    if (claim_next_contig(init, global)) {
                        continue;  // Continue reading from new contig
                    }
                } else if (bind->n_regions > 0) {
                    if (init->itr) {
                        hts_itr_destroy(init->itr);
                        init->itr = NULL;
                    }
                    while (init->next_region_idx < bind->n_regions) {
                        const char *query_region = bind->regions[init->next_region_idx++];
                        if (init->idx) {
                            init->itr = bcf_itr_querys(init->idx, init->hdr, query_region);
                        } else if (init->tbx) {
                            init->itr = tbx_itr_querys(init->tbx, query_region);
                        }
                        if (init->itr) break;
                    }
                    if (init->itr) {
                        continue;
                    }
                }
                init->done = 1;
                break;
            }
            
            // Unpack record
            bcf_unpack(init->rec, BCF_UN_ALL);
            
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
        }

        // Parse VEP annotation once per record if needed (only on first sample in tidy mode)
        vep_record_t* vep_rec = NULL;
        if (need_vep && (!tidy_mode || current_sample == 0)) {
            vep_rec = vep_record_parse_bcf(bind->vep_schema, init->hdr, init->rec);
        }
        
        // Process each requested column (using cached vectors)
        for (idx_t i = 0; i < init->column_count; i++) {
            idx_t col_id = init->column_ids[i];
            duckdb_vector vec = vectors[i];
            
            // Core VCF columns
            if (col_id == COL_CHROM) {
                const char* chrom = bcf_hdr_id2name(init->hdr, init->rec->rid);
                duckdb_vector_assign_string_element(vec, row_count, chrom ? chrom : ".");
            }
            else if (col_id == COL_POS) {
                int64_t* data = (int64_t*)duckdb_vector_get_data(vec);
                data[row_count] = init->rec->pos + 1;  // 1-based
            }
            else if (col_id == COL_ID) {
                const char* id = init->rec->d.id;
                if (id && strcmp(id, ".") != 0) {
                    duckdb_vector_assign_string_element(vec, row_count, id);
                } else {
                    duckdb_vector_ensure_validity_writable(vec);
                    uint64_t* validity = duckdb_vector_get_validity(vec);
                    set_validity_bit(validity, row_count, 0);
                }
            }
            else if (col_id == COL_REF) {
                const char* ref = init->rec->d.allele[0];
                duckdb_vector_assign_string_element(vec, row_count, ref ? ref : ".");
            }
            else if (col_id == COL_ALT) {
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
                list_data[row_count] = entry;
            }
            else if (col_id == COL_QUAL) {
                double* data = (double*)duckdb_vector_get_data(vec);
                if (bcf_float_is_missing(init->rec->qual)) {
                    duckdb_vector_ensure_validity_writable(vec);
                    uint64_t* validity = duckdb_vector_get_validity(vec);
                    set_validity_bit(validity, row_count, 0);
                    data[row_count] = 0.0;
                } else {
                    data[row_count] = init->rec->qual;
                }
            }
            else if (col_id == COL_FILTER) {
                // FILTER is a LIST(VARCHAR)
                duckdb_list_entry entry;
                entry.offset = duckdb_list_vector_get_size(vec);
                
                duckdb_vector child_vec = duckdb_list_vector_get_child(vec);
                
                if (init->rec->d.n_flt == 0) {
                    // No filters means PASS
                    entry.length = 1;
                    duckdb_list_vector_set_size(vec, entry.offset + 1);
                    duckdb_vector_assign_string_element(child_vec, entry.offset, "PASS");
                } else {
                    entry.length = init->rec->d.n_flt;
                    // Reserve space for all filters at once
                    duckdb_list_vector_set_size(vec, entry.offset + entry.length);
                    for (int f = 0; f < init->rec->d.n_flt; f++) {
                        const char* flt_name = bcf_hdr_int2id(init->hdr, BCF_DT_ID, 
                                                              init->rec->d.flt[f]);
                        duckdb_vector_assign_string_element(child_vec, entry.offset + f,
                                                            flt_name ? flt_name : ".");
                    }
                }
                
                duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
                list_data[row_count] = entry;
            }
            else if (bind->vep_schema &&
                     col_id >= (idx_t)bind->vep_col_start &&
                     col_id < (idx_t)(bind->vep_col_start + bind->n_vep_fields)) {
                int field_idx = col_id - bind->vep_col_start;
                const vep_field_t* field = vep_schema_get_field(bind->vep_schema, field_idx);

                duckdb_list_entry entry;
                entry.offset = duckdb_list_vector_get_size(vec);
                entry.length = (vep_rec && field) ? vep_rec->n_transcripts : 0;

                duckdb_vector child_vec = duckdb_list_vector_get_child(vec);

                if (entry.length > 0) {
                    duckdb_vector_ensure_validity_writable(vec);
                    uint64_t* parent_validity = duckdb_vector_get_validity(vec);
                    set_validity_bit(parent_validity, row_count, 1);

                    duckdb_list_vector_reserve(vec, entry.offset + entry.length);
                    duckdb_list_vector_set_size(vec, entry.offset + entry.length);

                    // Ensure child validity writable for nulls
                    duckdb_vector_ensure_validity_writable(child_vec);
                    uint64_t* child_validity = duckdb_vector_get_validity(child_vec);

                    if (field->type == VEP_TYPE_STRING) {
                        for (idx_t t = 0; t < entry.length; t++) {
                            const vep_value_t* val = vep_record_get_value(vep_rec, t, field_idx);
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
                            const vep_value_t* val = vep_record_get_value(vep_rec, t, field_idx);
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
                            const vep_value_t* val = vep_record_get_value(vep_rec, t, field_idx);
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
                            const vep_value_t* val = vep_record_get_value(vep_rec, t, field_idx);
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
                    set_validity_bit(validity, row_count, 0);
                    duckdb_list_vector_set_size(vec, entry.offset);
                }

                duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
                list_data[row_count] = entry;
            }
            else if (col_id >= (idx_t)bind->info_col_start && 
                     col_id < (idx_t)(bind->info_col_start + bind->n_info_fields)) {
                // INFO field
                int field_idx = col_id - bind->info_col_start;
                field_meta_t* field = &bind->info_fields[field_idx];
                const char* tag = field->name;
                
                if (field->header_type == BCF_HT_FLAG) {
                    // Boolean field
                    bool* data = (bool*)duckdb_vector_get_data(vec);
                    if (!info_loaded[field_idx]) {
                        int* dummy = NULL;
                        int ndummy = 0;
                        info_ret[field_idx] = bcf_get_info_flag(init->hdr, init->rec, tag, &dummy, &ndummy);
                        if (dummy) free(dummy);  // Only free if allocated
                        info_flag[field_idx] = (info_ret[field_idx] == 1);
                        info_loaded[field_idx] = 1;
                    }
                    data[row_count] = (info_flag[field_idx] != 0);
                }
                else if (field->header_type == BCF_HT_INT) {
                    int32_t* values = NULL;
                    int ret_info = 0;
                    if (!info_loaded[field_idx]) {
                        info_ret[field_idx] = bcf_get_info_int32(init->hdr, init->rec, tag,
                                                                  &info_i32[field_idx],
                                                                  &info_n_values[field_idx]);
                        info_loaded[field_idx] = 1;
                    }
                    values = info_i32[field_idx];
                    ret_info = info_ret[field_idx];
                    
                    if (ret_info > 0 && values) {
                        if (field->is_list) {
                            // List of integers
                            duckdb_list_entry entry;
                            entry.offset = duckdb_list_vector_get_size(vec);
                            entry.length = 0;
                            
                            duckdb_vector child_vec = duckdb_list_vector_get_child(vec);
                            
                            // First count valid values
                            for (int v = 0; v < ret_info; v++) {
                                if (values[v] != bcf_int32_missing && values[v] != bcf_int32_vector_end) {
                                    entry.length++;
                                }
                            }
                            
                            // Reserve and set size
                            if (entry.length > 0) {
                                duckdb_list_vector_reserve(vec, entry.offset + entry.length);
                                duckdb_list_vector_set_size(vec, entry.offset + entry.length);
                                
                                // Now fill in the values
                                int32_t* child_data = (int32_t*)duckdb_vector_get_data(child_vec);
                                int write_idx = 0;
                                for (int v = 0; v < ret_info; v++) {
                                    if (values[v] != bcf_int32_missing && values[v] != bcf_int32_vector_end) {
                                        child_data[entry.offset + write_idx] = values[v];
                                        write_idx++;
                                    }
                                }
                            }
                            
                            duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
                            list_data[row_count] = entry;
                        } else {
                            // Scalar integer
                            int32_t* data = (int32_t*)duckdb_vector_get_data(vec);
                            if (values[0] != bcf_int32_missing) {
                                data[row_count] = values[0];
                            } else {
                                duckdb_vector_ensure_validity_writable(vec);
                                uint64_t* validity = duckdb_vector_get_validity(vec);
                                set_validity_bit(validity, row_count, 0);
                            }
                        }
                    } else {
                        // No data - NULL
                        duckdb_vector_ensure_validity_writable(vec);
                        uint64_t* validity = duckdb_vector_get_validity(vec);
                        set_validity_bit(validity, row_count, 0);
                        
                        if (field->is_list) {
                            duckdb_list_entry entry = {duckdb_list_vector_get_size(vec), 0};
                            duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
                            list_data[row_count] = entry;
                        }
                    }
                }
                else if (field->header_type == BCF_HT_REAL) {
                    float* values = NULL;
                    int ret_info = 0;
                    if (!info_loaded[field_idx]) {
                        info_ret[field_idx] = bcf_get_info_float(init->hdr, init->rec, tag,
                                                                  &info_f32[field_idx],
                                                                  &info_n_values[field_idx]);
                        info_loaded[field_idx] = 1;
                    }
                    values = info_f32[field_idx];
                    ret_info = info_ret[field_idx];
                    
                    if (ret_info > 0 && values) {
                        if (field->is_list) {
                            // List of floats
                            duckdb_list_entry entry;
                            entry.offset = duckdb_list_vector_get_size(vec);
                            entry.length = 0;
                            
                            duckdb_vector child_vec = duckdb_list_vector_get_child(vec);
                            
                            // First count valid values
                            for (int v = 0; v < ret_info; v++) {
                                if (!bcf_float_is_missing(values[v]) && !bcf_float_is_vector_end(values[v])) {
                                    entry.length++;
                                }
                            }
                            
                            // Reserve and set size
                            if (entry.length > 0) {
                                duckdb_list_vector_reserve(vec, entry.offset + entry.length);
                                duckdb_list_vector_set_size(vec, entry.offset + entry.length);
                                
                                // Now fill in the values
                                float* child_data = (float*)duckdb_vector_get_data(child_vec);
                                int write_idx = 0;
                                for (int v = 0; v < ret_info; v++) {
                                    if (!bcf_float_is_missing(values[v]) && !bcf_float_is_vector_end(values[v])) {
                                        child_data[entry.offset + write_idx] = values[v];
                                        write_idx++;
                                    }
                                }
                            }
                            
                            duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
                            list_data[row_count] = entry;
                        } else {
                            // Scalar float
                            float* data = (float*)duckdb_vector_get_data(vec);
                            if (!bcf_float_is_missing(values[0])) {
                                data[row_count] = values[0];
                            } else {
                                duckdb_vector_ensure_validity_writable(vec);
                                uint64_t* validity = duckdb_vector_get_validity(vec);
                                set_validity_bit(validity, row_count, 0);
                            }
                        }
                    } else {
                        duckdb_vector_ensure_validity_writable(vec);
                        uint64_t* validity = duckdb_vector_get_validity(vec);
                        set_validity_bit(validity, row_count, 0);
                        
                        if (field->is_list) {
                            duckdb_list_entry entry = {duckdb_list_vector_get_size(vec), 0};
                            duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
                            list_data[row_count] = entry;
                        }
                    }
                }
                else {
                    // String type
                    char* value = NULL;
                    int ret_info = 0;
                    if (!info_loaded[field_idx]) {
                        info_ret[field_idx] = bcf_get_info_string(init->hdr, init->rec, tag,
                                                                   &info_str[field_idx],
                                                                   &info_n_values[field_idx]);
                        info_loaded[field_idx] = 1;
                    }
                    value = info_str[field_idx];
                    ret_info = info_ret[field_idx];
                    
                    if (ret_info > 0 && value && strcmp(value, ".") != 0) {
                        if (field->is_list) {
                            // Use optimized single-pass comma-separated list processing
                            process_comma_separated_list(vec, row_count, value);
                        } else {
                            duckdb_vector_assign_string_element(vec, row_count, value);
                        }
                    } else {
                        duckdb_vector_ensure_validity_writable(vec);
                        uint64_t* validity = duckdb_vector_get_validity(vec);
                        set_validity_bit(validity, row_count, 0);
                        
                        if (field->is_list) {
                            duckdb_list_entry entry = {duckdb_list_vector_get_size(vec), 0};
                            duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
                            list_data[row_count] = entry;
                        }
                    }
                }
            }
            else if (tidy_mode && col_id == (idx_t)bind->sample_id_col_idx) {
                // SAMPLE_ID column in tidy mode
                duckdb_vector_assign_string_element(vec, row_count, bind->sample_names[current_sample]);
            }
            else if (col_id >= (idx_t)bind->format_col_start) {
                // FORMAT field for a sample
                int format_col_idx = col_id - bind->format_col_start;
                int sample_idx, field_idx;
                
                if (tidy_mode) {
                    // In tidy mode: column index directly maps to field (no sample suffix)
                    // SAMPLE_ID column is at sample_id_col_idx, FORMAT columns start after that
                    field_idx = format_col_idx;
                    sample_idx = current_sample;
                } else {
                    // Wide mode: column index encodes both sample and field
                    sample_idx = format_col_idx / bind->n_format_fields;
                    field_idx = format_col_idx % bind->n_format_fields;
                }
                
                if (sample_idx < bind->n_samples && field_idx < bind->n_format_fields) {
                    field_meta_t* field = &bind->format_fields[field_idx];
                    const char* tag = field->name;
                    
                    if (field->header_type == BCF_HT_INT) {
                        int32_t* values = NULL;
                        int ret_fmt = 0;
                        if (!fmt_loaded[field_idx]) {
                            fmt_ret[field_idx] = bcf_get_format_int32(init->hdr, init->rec, tag,
                                                                       &fmt_i32[field_idx],
                                                                       &fmt_n_values[field_idx]);
                            fmt_loaded[field_idx] = 1;
                        }
                        values = fmt_i32[field_idx];
                        ret_fmt = fmt_ret[field_idx];
                        
                        if (ret_fmt > 0 && values) {
                            int vals_per_sample = ret_fmt / bind->n_samples;
                            int32_t* sample_vals = values + sample_idx * vals_per_sample;
                            
                            if (field->is_list) {
                                duckdb_list_entry entry;
                                entry.offset = duckdb_list_vector_get_size(vec);
                                entry.length = 0;
                                
                                duckdb_vector child_vec = duckdb_list_vector_get_child(vec);
                                
                                // First count valid values
                                for (int v = 0; v < vals_per_sample; v++) {
                                    if (sample_vals[v] != bcf_int32_missing && 
                                        sample_vals[v] != bcf_int32_vector_end) {
                                        entry.length++;
                                    }
                                }
                                
                                // Reserve and set size
                                if (entry.length > 0) {
                                    duckdb_list_vector_reserve(vec, entry.offset + entry.length);
                                    duckdb_list_vector_set_size(vec, entry.offset + entry.length);
                                    
                                    // Now fill in the values
                                    int32_t* child_data = (int32_t*)duckdb_vector_get_data(child_vec);
                                    int write_idx = 0;
                                    for (int v = 0; v < vals_per_sample; v++) {
                                        if (sample_vals[v] != bcf_int32_missing && 
                                            sample_vals[v] != bcf_int32_vector_end) {
                                            child_data[entry.offset + write_idx] = sample_vals[v];
                                            write_idx++;
                                        }
                                    }
                                }
                                
                                duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
                                list_data[row_count] = entry;
                            } else {
                                int32_t* data = (int32_t*)duckdb_vector_get_data(vec);
                                if (sample_vals[0] != bcf_int32_missing) {
                                    data[row_count] = sample_vals[0];
                                } else {
                                    duckdb_vector_ensure_validity_writable(vec);
                                    uint64_t* validity = duckdb_vector_get_validity(vec);
                                    set_validity_bit(validity, row_count, 0);
                                }
                            }
                        } else {
                            duckdb_vector_ensure_validity_writable(vec);
                            uint64_t* validity = duckdb_vector_get_validity(vec);
                            set_validity_bit(validity, row_count, 0);
                            
                            if (field->is_list) {
                                duckdb_list_entry entry = {duckdb_list_vector_get_size(vec), 0};
                                duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
                                list_data[row_count] = entry;
                            }
                        }
                    }
                    else if (field->header_type == BCF_HT_REAL) {
                        float* values = NULL;
                        int ret_fmt = 0;
                        if (!fmt_loaded[field_idx]) {
                            fmt_ret[field_idx] = bcf_get_format_float(init->hdr, init->rec, tag,
                                                                       &fmt_f32[field_idx],
                                                                       &fmt_n_values[field_idx]);
                            fmt_loaded[field_idx] = 1;
                        }
                        values = fmt_f32[field_idx];
                        ret_fmt = fmt_ret[field_idx];
                        
                        if (ret_fmt > 0 && values) {
                            int vals_per_sample = ret_fmt / bind->n_samples;
                            float* sample_vals = values + sample_idx * vals_per_sample;
                            
                            if (field->is_list) {
                                duckdb_list_entry entry;
                                entry.offset = duckdb_list_vector_get_size(vec);
                                entry.length = 0;
                                
                                duckdb_vector child_vec = duckdb_list_vector_get_child(vec);
                                
                                // First count valid values
                                for (int v = 0; v < vals_per_sample; v++) {
                                    if (!bcf_float_is_missing(sample_vals[v]) && 
                                        !bcf_float_is_vector_end(sample_vals[v])) {
                                        entry.length++;
                                    }
                                }
                                
                                // Reserve and set size
                                if (entry.length > 0) {
                                    duckdb_list_vector_reserve(vec, entry.offset + entry.length);
                                    duckdb_list_vector_set_size(vec, entry.offset + entry.length);
                                    
                                    // Now fill in the values
                                    float* child_data = (float*)duckdb_vector_get_data(child_vec);
                                    int write_idx = 0;
                                    for (int v = 0; v < vals_per_sample; v++) {
                                        if (!bcf_float_is_missing(sample_vals[v]) && 
                                            !bcf_float_is_vector_end(sample_vals[v])) {
                                            child_data[entry.offset + write_idx] = sample_vals[v];
                                            write_idx++;
                                        }
                                    }
                                }
                                
                                duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
                                list_data[row_count] = entry;
                            } else {
                                float* data = (float*)duckdb_vector_get_data(vec);
                                if (!bcf_float_is_missing(sample_vals[0])) {
                                    data[row_count] = sample_vals[0];
                                } else {
                                    duckdb_vector_ensure_validity_writable(vec);
                                    uint64_t* validity = duckdb_vector_get_validity(vec);
                                    set_validity_bit(validity, row_count, 0);
                                }
                            }
                        } else {
                            duckdb_vector_ensure_validity_writable(vec);
                            uint64_t* validity = duckdb_vector_get_validity(vec);
                            set_validity_bit(validity, row_count, 0);
                            
                            if (field->is_list) {
                                duckdb_list_entry entry = {duckdb_list_vector_get_size(vec), 0};
                                duckdb_list_entry* list_data = (duckdb_list_entry*)duckdb_vector_get_data(vec);
                                list_data[row_count] = entry;
                            }
                        }
                    }
                    else {
                        // String type - GT needs special handling
                        if (strcmp(tag, "GT") == 0) {
                            // GT is stored as encoded integers, use bcf_get_genotypes()
                            if (!gt_loaded) {
                                gt_ret = bcf_get_genotypes(init->hdr, init->rec, &gt_arr, &gt_n_values);
                                gt_loaded = 1;
                            }
                            int ret_gt = gt_ret;
                            
                            if (ret_gt > 0 && gt_arr) {
                                int ploidy = ret_gt / bind->n_samples;
                                int32_t* sample_gt = gt_arr + sample_idx * ploidy;
                                
                                // Build GT string (e.g., "0/1", "1|1", "./.")
                                char gt_str[64];
                                int pos = 0;
                                
                                for (int p = 0; p < ploidy && pos < 60; p++) {
                                    if (p > 0) {
                                        // Add separator: '|' for phased, '/' for unphased
                                        gt_str[pos++] = bcf_gt_is_phased(sample_gt[p]) ? '|' : '/';
                                    }
                                    
                                    if (sample_gt[p] == bcf_int32_vector_end) {
                                        break;  // End of genotype
                                    } else if (bcf_gt_is_missing(sample_gt[p])) {
                                        gt_str[pos++] = '.';
                                    } else {
                                        int allele = bcf_gt_allele(sample_gt[p]);
                                        pos += snprintf(gt_str + pos, sizeof(gt_str) - pos, "%d", allele);
                                    }
                                }
                                gt_str[pos] = '\0';
                                
                                if (pos > 0) {
                                    duckdb_vector_assign_string_element(vec, row_count, gt_str);
                                } else {
                                    duckdb_vector_ensure_validity_writable(vec);
                                    uint64_t* validity = duckdb_vector_get_validity(vec);
                                    set_validity_bit(validity, row_count, 0);
                                }
                            } else {
                                duckdb_vector_ensure_validity_writable(vec);
                                uint64_t* validity = duckdb_vector_get_validity(vec);
                                set_validity_bit(validity, row_count, 0);
                            }
                        } else {
                            // Other string FORMAT fields
                            char** values = NULL;
                            int ret_fmt = 0;
                            if (!fmt_loaded[field_idx]) {
                                fmt_ret[field_idx] = bcf_get_format_string(init->hdr, init->rec, tag,
                                                                            &fmt_str[field_idx],
                                                                            &fmt_n_values[field_idx]);
                                fmt_loaded[field_idx] = 1;
                            }
                            values = fmt_str[field_idx];
                            ret_fmt = fmt_ret[field_idx];
                            
                            if (ret_fmt > 0 && values && values[sample_idx]) {
                                duckdb_vector_assign_string_element(vec, row_count, values[sample_idx]);
                            } else {
                                duckdb_vector_ensure_validity_writable(vec);
                                uint64_t* validity = duckdb_vector_get_validity(vec);
                                set_validity_bit(validity, row_count, 0);
                            }
                            
                        }
                    }
                }
            }
        }
        if (vep_rec) {
            vep_record_destroy(vep_rec);
        }

        row_count++;
        init->current_row++;
        
        // In tidy mode, advance to next sample (or mark record as consumed)
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
    
    // Cleanup cached vectors
    duckdb_free(vectors);
    if (bind->n_format_fields > 0) {
        for (int i = 0; i < bind->n_format_fields; i++) {
            if (fmt_i32 && fmt_i32[i]) free(fmt_i32[i]);
            if (fmt_f32 && fmt_f32[i]) free(fmt_f32[i]);
            if (fmt_str && fmt_str[i]) free(fmt_str[i]);
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
    
    duckdb_data_chunk_set_size(output, row_count);
}

// =============================================================================
// Register the bcf_read Table Function
// =============================================================================

void register_read_bcf_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "read_bcf");
    
    // Parameters
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_table_function_add_parameter(tf, varchar_type);  // file_path
    duckdb_table_function_add_named_parameter(tf, "region", varchar_type);  // optional region
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);  // optional explicit index path
    duckdb_table_function_add_named_parameter(tf, "tidy_format", bool_type);  // optional tidy format
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bool_type);
    
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
