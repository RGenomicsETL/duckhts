/**
 * DuckHTS interval readers and FASTA interval composition helpers.
 *
 * read_bed(path, region := NULL, index_path := NULL)
 *   -> BED3-BED12 reader with canonical typed columns and trailing extras.
 *
 * fasta_nuc(fasta_path, bed_path := NULL, bin_width := NULL, region := NULL,
 *           index_path := NULL, bed_index_path := NULL, include_seq := FALSE)
 *   -> bedtools nuc-style interval composition metrics over either supplied
 *      BED intervals or generated fixed-width bins.
 */

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <ctype.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/tbx.h>

#define INTERVAL_BATCH_SIZE 2048

enum {
    BED_COL_CHROM = 0,
    BED_COL_START,
    BED_COL_END,
    BED_COL_NAME,
    BED_COL_SCORE,
    BED_COL_STRAND,
    BED_COL_THICK_START,
    BED_COL_THICK_END,
    BED_COL_ITEM_RGB,
    BED_COL_BLOCK_COUNT,
    BED_COL_BLOCK_SIZES,
    BED_COL_BLOCK_STARTS,
    BED_COL_EXTRA,
    BED_COL_COUNT
};

enum {
    NUC_COL_CHROM = 0,
    NUC_COL_START,
    NUC_COL_END,
    NUC_COL_PCT_AT,
    NUC_COL_PCT_GC,
    NUC_COL_NUM_A,
    NUC_COL_NUM_C,
    NUC_COL_NUM_G,
    NUC_COL_NUM_T,
    NUC_COL_NUM_N,
    NUC_COL_NUM_OTHER,
    NUC_COL_SEQ_LEN,
    NUC_COL_SEQ,
    NUC_COL_COUNT
};

typedef struct {
    char *file_path;
    char *index_path;
    char *region;
} bed_bind_data_t;

typedef struct {
    htsFile *fp;
    tbx_t *tbx;
    hts_itr_t *itr;
    kstring_t line;
    bool finished;
    idx_t *column_ids;
    idx_t n_projected_cols;
} bed_init_data_t;

typedef enum {
    FASTA_NUC_MODE_BED = 0,
    FASTA_NUC_MODE_BINS = 1
} fasta_nuc_mode_t;

typedef struct {
    char *fasta_path;
    char *index_path;
    char *bed_path;
    char *bed_index_path;
    char *region;
    int64_t bin_width;
    bool include_seq;
    fasta_nuc_mode_t mode;
} fasta_nuc_bind_data_t;

typedef struct {
    faidx_t *fai;
    fasta_nuc_bind_data_t *bind;

    htsFile *bed_fp;
    tbx_t *bed_tbx;
    hts_itr_t *bed_itr;
    kstring_t bed_line;
    bool bed_finished;

    int region_tid;
    hts_pos_t region_beg;
    hts_pos_t region_end;
    const char *region_seq;
    bool has_region;

    int current_tid;
    hts_pos_t next_bin_start;
    hts_pos_t current_seq_end_exclusive;
    bool done;

    idx_t *column_ids;
    idx_t n_projected_cols;
} fasta_nuc_init_data_t;

static inline void set_null(duckdb_vector vec, idx_t row) {
    duckdb_vector_ensure_validity_writable(vec);
    uint64_t *v = duckdb_vector_get_validity(vec);
    v[row / 64] &= ~((uint64_t)1 << (row % 64));
}

static char *duckdb_strdup_local(const char *s) {
    if (!s) return NULL;
    size_t n = strlen(s) + 1;
    char *out = (char *)duckdb_malloc(n);
    if (!out) return NULL;
    memcpy(out, s, n);
    return out;
}

static int parse_int64_span_local(const char *s, int len, int64_t *out) {
    if (!s || len <= 0) return 0;
    char *tmp = (char *)malloc((size_t)len + 1);
    if (!tmp) return 0;
    memcpy(tmp, s, (size_t)len);
    tmp[len] = '\0';
    char *end = NULL;
    long long v = strtoll(tmp, &end, 10);
    int ok = (end && *end == '\0');
    if (ok) *out = (int64_t)v;
    free(tmp);
    return ok;
}

static bool is_meta_bed_line(const char *s) {
    if (!s || !*s) return true;
    if (s[0] == '#') return true;
    if (strncmp(s, "track", 5) == 0) return true;
    if (strncmp(s, "browser", 7) == 0) return true;
    return false;
}

static int count_tab_fields(const char *s) {
    if (!s || !*s) return 0;
    int n = 1;
    while (*s) {
        if (*s == '\t') n++;
        s++;
    }
    return n;
}

static const char *get_field_span(const char *s, int idx, int *len) {
    int cur = 0;
    const char *start = s;
    const char *p = s;
    while (*p) {
        if (*p == '\t') {
            if (cur == idx) {
                *len = (int)(p - start);
                return start;
            }
            cur++;
            start = p + 1;
        }
        p++;
    }
    if (cur == idx) {
        *len = (int)(p - start);
        return start;
    }
    *len = 0;
    return NULL;
}

static const char *get_extra_span(const char *s, int start_field, int *len) {
    int cur = 0;
    const char *p = s;
    while (*p) {
        if (cur == start_field) {
            *len = (int)strlen(p);
            return p;
        }
        if (*p == '\t') cur++;
        p++;
    }
    *len = 0;
    return NULL;
}

static void destroy_bed_bind(void *data) {
    bed_bind_data_t *bind = (bed_bind_data_t *)data;
    if (!bind) return;
    if (bind->file_path) duckdb_free(bind->file_path);
    if (bind->index_path) duckdb_free(bind->index_path);
    if (bind->region) duckdb_free(bind->region);
    duckdb_free(bind);
}

static void destroy_bed_init(void *data) {
    bed_init_data_t *init = (bed_init_data_t *)data;
    if (!init) return;
    if (init->itr) hts_itr_destroy(init->itr);
    if (init->tbx) tbx_destroy(init->tbx);
    if (init->fp) hts_close(init->fp);
    free(init->line.s);
    if (init->column_ids) duckdb_free(init->column_ids);
    duckdb_free(init);
}

static void add_bed_result_columns(duckdb_bind_info info) {
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_bind_add_result_column(info, "chrom", varchar_type);
    duckdb_bind_add_result_column(info, "start", bigint_type);
    duckdb_bind_add_result_column(info, "end", bigint_type);
    duckdb_bind_add_result_column(info, "name", varchar_type);
    duckdb_bind_add_result_column(info, "score", varchar_type);
    duckdb_bind_add_result_column(info, "strand", varchar_type);
    duckdb_bind_add_result_column(info, "thick_start", bigint_type);
    duckdb_bind_add_result_column(info, "thick_end", bigint_type);
    duckdb_bind_add_result_column(info, "item_rgb", varchar_type);
    duckdb_bind_add_result_column(info, "block_count", bigint_type);
    duckdb_bind_add_result_column(info, "block_sizes", varchar_type);
    duckdb_bind_add_result_column(info, "block_starts", varchar_type);
    duckdb_bind_add_result_column(info, "extra", varchar_type);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
}

static void read_bed_bind(duckdb_bind_info info) {
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    char *file_path = duckdb_get_varchar(path_val);
    duckdb_destroy_value(&path_val);

    if (!file_path || file_path[0] == '\0') {
        duckdb_bind_set_error(info, "read_bed requires a file path");
        if (file_path) duckdb_free(file_path);
        return;
    }

    char *region = NULL;
    duckdb_value region_val = duckdb_bind_get_named_parameter(info, "region");
    if (region_val && !duckdb_is_null_value(region_val)) {
        region = duckdb_get_varchar(region_val);
    }
    if (region_val) duckdb_destroy_value(&region_val);

    char *index_path = NULL;
    duckdb_value index_val = duckdb_bind_get_named_parameter(info, "index_path");
    if (index_val && !duckdb_is_null_value(index_val)) {
        index_path = duckdb_get_varchar(index_val);
    }
    if (index_val) duckdb_destroy_value(&index_val);

    htsFile *fp = hts_open(file_path, "r");
    if (!fp) {
        char err[512];
        snprintf(err, sizeof(err), "read_bed: failed to open file: %s", file_path);
        duckdb_bind_set_error(info, err);
        duckdb_free(file_path);
        if (region) duckdb_free(region);
        if (index_path) duckdb_free(index_path);
        return;
    }
    hts_close(fp);

    if (region) {
        tbx_t *tbx = tbx_index_load3(file_path, index_path, HTS_IDX_SILENT_FAIL);
        if (!tbx) {
            duckdb_bind_set_error(info, "read_bed: region queries require a tabix index");
            duckdb_free(file_path);
            duckdb_free(region);
            if (index_path) duckdb_free(index_path);
            return;
        }
        tbx_destroy(tbx);
    }

    add_bed_result_columns(info);

    bed_bind_data_t *bind = (bed_bind_data_t *)duckdb_malloc(sizeof(bed_bind_data_t));
    bind->file_path = file_path;
    bind->index_path = index_path;
    bind->region = region;
    duckdb_bind_set_bind_data(info, bind, destroy_bed_bind);
}

static void read_bed_init(duckdb_init_info info) {
    bed_bind_data_t *bind = (bed_bind_data_t *)duckdb_init_get_bind_data(info);
    bed_init_data_t *init = (bed_init_data_t *)duckdb_malloc(sizeof(bed_init_data_t));
    memset(init, 0, sizeof(*init));

    init->fp = hts_open(bind->file_path, "r");
    if (!init->fp) {
        duckdb_init_set_error(info, "read_bed: failed to open file during init");
        destroy_bed_init(init);
        return;
    }

    if (bind->region) {
        init->tbx = tbx_index_load3(bind->file_path, bind->index_path, HTS_IDX_SILENT_FAIL);
        if (!init->tbx) {
            duckdb_init_set_error(info, "read_bed: failed to load tabix index");
            destroy_bed_init(init);
            return;
        }
        init->itr = tbx_itr_querys(init->tbx, bind->region);
        if (!init->itr) {
            duckdb_init_set_error(info, "read_bed: failed to create region iterator");
            destroy_bed_init(init);
            return;
        }
    }

    init->n_projected_cols = duckdb_init_get_column_count(info);
    init->column_ids = (idx_t *)duckdb_malloc(sizeof(idx_t) * init->n_projected_cols);
    for (idx_t i = 0; i < init->n_projected_cols; i++) {
        init->column_ids[i] = duckdb_init_get_column_index(info, i);
    }
    duckdb_init_set_init_data(info, init, destroy_bed_init);
}

static int next_bed_line(bed_init_data_t *init) {
    while (!init->finished) {
        int ret = init->itr ? tbx_itr_next(init->fp, init->tbx, init->itr, &init->line)
                            : hts_getline(init->fp, '\n', &init->line);
        if (ret < 0) {
            init->finished = true;
            return 0;
        }
        if (init->line.l == 0 || is_meta_bed_line(init->line.s)) continue;
        return 1;
    }
    return 0;
}

static void read_bed_scan(duckdb_function_info info, duckdb_data_chunk output) {
    bed_init_data_t *init = (bed_init_data_t *)duckdb_function_get_init_data(info);
    if (!init || init->finished) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    idx_t row_count = 0;
    idx_t col_count = duckdb_data_chunk_get_column_count(output);
    duckdb_vector vectors[BED_COL_COUNT];
    for (idx_t c = 0; c < col_count; c++) {
        vectors[c] = duckdb_data_chunk_get_vector(output, c);
    }

    while (row_count < INTERVAL_BATCH_SIZE && next_bed_line(init)) {
        int n_fields = count_tab_fields(init->line.s);
        if (n_fields < 3) {
            duckdb_function_set_error(info, "read_bed: BED line has fewer than 3 tab-delimited fields");
            init->finished = true;
            duckdb_data_chunk_set_size(output, 0);
            return;
        }

        for (idx_t c = 0; c < col_count; c++) {
            int logical_col = (int)init->column_ids[c];
            int len = 0;
            const char *field = NULL;
            int64_t ival = 0;
            switch (logical_col) {
                case BED_COL_CHROM:
                case BED_COL_NAME:
                case BED_COL_SCORE:
                case BED_COL_STRAND:
                case BED_COL_ITEM_RGB:
                case BED_COL_BLOCK_SIZES:
                case BED_COL_BLOCK_STARTS:
                    field = get_field_span(init->line.s, logical_col, &len);
                    if (!field || len == 0) {
                        set_null(vectors[c], row_count);
                    } else {
                        duckdb_vector_assign_string_element_len(vectors[c], row_count, field, len);
                    }
                    break;
                case BED_COL_START:
                case BED_COL_END:
                case BED_COL_THICK_START:
                case BED_COL_THICK_END:
                case BED_COL_BLOCK_COUNT:
                    field = get_field_span(
                        init->line.s,
                        logical_col == BED_COL_START ? 1 :
                        logical_col == BED_COL_END ? 2 :
                        logical_col == BED_COL_THICK_START ? 6 :
                        logical_col == BED_COL_THICK_END ? 7 : 9,
                        &len
                    );
                    if (!field || len == 0) {
                        set_null(vectors[c], row_count);
                    } else if (!parse_int64_span_local(field, len, &ival)) {
                        set_null(vectors[c], row_count);
                    } else {
                        int64_t *data = (int64_t *)duckdb_vector_get_data(vectors[c]);
                        data[row_count] = ival;
                    }
                    break;
                case BED_COL_EXTRA:
                    field = get_extra_span(init->line.s, 12, &len);
                    if (!field || len == 0) {
                        set_null(vectors[c], row_count);
                    } else {
                        duckdb_vector_assign_string_element_len(vectors[c], row_count, field, len);
                    }
                    break;
                default:
                    set_null(vectors[c], row_count);
                    break;
            }
        }
        row_count++;
    }

    duckdb_data_chunk_set_size(output, row_count);
}

static void destroy_fasta_nuc_bind(void *data) {
    fasta_nuc_bind_data_t *bind = (fasta_nuc_bind_data_t *)data;
    if (!bind) return;
    if (bind->fasta_path) duckdb_free(bind->fasta_path);
    if (bind->index_path) duckdb_free(bind->index_path);
    if (bind->bed_path) duckdb_free(bind->bed_path);
    if (bind->bed_index_path) duckdb_free(bind->bed_index_path);
    if (bind->region) duckdb_free(bind->region);
    duckdb_free(bind);
}

static void destroy_fasta_nuc_init(void *data) {
    fasta_nuc_init_data_t *init = (fasta_nuc_init_data_t *)data;
    if (!init) return;
    if (init->fai) fai_destroy(init->fai);
    if (init->bed_itr) hts_itr_destroy(init->bed_itr);
    if (init->bed_tbx) tbx_destroy(init->bed_tbx);
    if (init->bed_fp) hts_close(init->bed_fp);
    free(init->bed_line.s);
    if (init->column_ids) duckdb_free(init->column_ids);
    duckdb_free(init);
}

static void add_fasta_nuc_columns(duckdb_bind_info info, bool include_seq) {
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type double_type = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_bind_add_result_column(info, "chrom", varchar_type);
    duckdb_bind_add_result_column(info, "start", bigint_type);
    duckdb_bind_add_result_column(info, "end", bigint_type);
    duckdb_bind_add_result_column(info, "pct_at", double_type);
    duckdb_bind_add_result_column(info, "pct_gc", double_type);
    duckdb_bind_add_result_column(info, "num_a", bigint_type);
    duckdb_bind_add_result_column(info, "num_c", bigint_type);
    duckdb_bind_add_result_column(info, "num_g", bigint_type);
    duckdb_bind_add_result_column(info, "num_t", bigint_type);
    duckdb_bind_add_result_column(info, "num_n", bigint_type);
    duckdb_bind_add_result_column(info, "num_other", bigint_type);
    duckdb_bind_add_result_column(info, "seq_len", bigint_type);
    if (include_seq) {
        duckdb_bind_add_result_column(info, "seq", varchar_type);
    }
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&double_type);
}

static void fasta_nuc_bind(duckdb_bind_info info) {
    duckdb_value fasta_val = duckdb_bind_get_parameter(info, 0);
    char *fasta_path = duckdb_get_varchar(fasta_val);
    duckdb_destroy_value(&fasta_val);
    if (!fasta_path || fasta_path[0] == '\0') {
        duckdb_bind_set_error(info, "fasta_nuc requires a FASTA path");
        if (fasta_path) duckdb_free(fasta_path);
        return;
    }

    char *bed_path = NULL;
    duckdb_value bed_val = duckdb_bind_get_named_parameter(info, "bed_path");
    if (bed_val && !duckdb_is_null_value(bed_val)) bed_path = duckdb_get_varchar(bed_val);
    if (bed_val) duckdb_destroy_value(&bed_val);

    int64_t bin_width = 0;
    bool has_bin_width = false;
    duckdb_value bin_val = duckdb_bind_get_named_parameter(info, "bin_width");
    if (bin_val && !duckdb_is_null_value(bin_val)) {
        bin_width = duckdb_get_int64(bin_val);
        has_bin_width = true;
    }
    if (bin_val) duckdb_destroy_value(&bin_val);

    if ((bed_path && has_bin_width) || (!bed_path && !has_bin_width)) {
        duckdb_bind_set_error(info, "fasta_nuc requires exactly one of bed_path or bin_width");
        duckdb_free(fasta_path);
        if (bed_path) duckdb_free(bed_path);
        return;
    }
    if (has_bin_width && bin_width <= 0) {
        duckdb_bind_set_error(info, "fasta_nuc bin_width must be > 0");
        duckdb_free(fasta_path);
        if (bed_path) duckdb_free(bed_path);
        return;
    }

    char *region = NULL;
    duckdb_value region_val = duckdb_bind_get_named_parameter(info, "region");
    if (region_val && !duckdb_is_null_value(region_val)) region = duckdb_get_varchar(region_val);
    if (region_val) duckdb_destroy_value(&region_val);

    char *index_path = NULL;
    duckdb_value index_val = duckdb_bind_get_named_parameter(info, "index_path");
    if (index_val && !duckdb_is_null_value(index_val)) index_path = duckdb_get_varchar(index_val);
    if (index_val) duckdb_destroy_value(&index_val);

    char *bed_index_path = NULL;
    duckdb_value bed_index_val = duckdb_bind_get_named_parameter(info, "bed_index_path");
    if (bed_index_val && !duckdb_is_null_value(bed_index_val)) bed_index_path = duckdb_get_varchar(bed_index_val);
    if (bed_index_val) duckdb_destroy_value(&bed_index_val);

    bool include_seq = false;
    duckdb_value include_seq_val = duckdb_bind_get_named_parameter(info, "include_seq");
    if (include_seq_val && !duckdb_is_null_value(include_seq_val)) include_seq = duckdb_get_bool(include_seq_val);
    if (include_seq_val) duckdb_destroy_value(&include_seq_val);

    faidx_t *fai = fai_load3_format(fasta_path, index_path, NULL, 0, FAI_FASTA);
    if (!fai) {
        duckdb_bind_set_error(info, "fasta_nuc: failed to open FASTA index");
        duckdb_free(fasta_path);
        if (bed_path) duckdb_free(bed_path);
        if (region) duckdb_free(region);
        if (index_path) duckdb_free(index_path);
        if (bed_index_path) duckdb_free(bed_index_path);
        return;
    }
    fai_destroy(fai);

    add_fasta_nuc_columns(info, include_seq);

    fasta_nuc_bind_data_t *bind = (fasta_nuc_bind_data_t *)duckdb_malloc(sizeof(fasta_nuc_bind_data_t));
    bind->fasta_path = fasta_path;
    bind->index_path = index_path;
    bind->bed_path = bed_path;
    bind->bed_index_path = bed_index_path;
    bind->region = region;
    bind->bin_width = bin_width;
    bind->include_seq = include_seq;
    bind->mode = bed_path ? FASTA_NUC_MODE_BED : FASTA_NUC_MODE_BINS;
    duckdb_bind_set_bind_data(info, bind, destroy_fasta_nuc_bind);
}

static bool init_fasta_region(fasta_nuc_init_data_t *init, const char *region) {
    if (!region || !*region) return true;
    int tid = -1;
    hts_pos_t beg = 0, end = 0;
    const char *remaining = fai_parse_region(init->fai, region, &tid, &beg, &end, 0);
    if (!remaining || *remaining != '\0') return false;
    if (fai_adjust_region(init->fai, tid, &beg, &end) != 0) return false;
    init->region_tid = tid;
    init->region_beg = beg;
    init->region_end = end;
    init->region_seq = faidx_iseq(init->fai, tid);
    init->has_region = true;
    return true;
}

static void fasta_nuc_init(duckdb_init_info info) {
    fasta_nuc_bind_data_t *bind = (fasta_nuc_bind_data_t *)duckdb_init_get_bind_data(info);
    fasta_nuc_init_data_t *init = (fasta_nuc_init_data_t *)duckdb_malloc(sizeof(fasta_nuc_init_data_t));
    memset(init, 0, sizeof(*init));
    init->bind = bind;
    init->fai = fai_load3_format(bind->fasta_path, bind->index_path, NULL, 0, FAI_FASTA);
    if (!init->fai) {
        duckdb_init_set_error(info, "fasta_nuc: failed to load FASTA index");
        destroy_fasta_nuc_init(init);
        return;
    }
    if (!init_fasta_region(init, bind->region)) {
        duckdb_init_set_error(info, "fasta_nuc: invalid FASTA region");
        destroy_fasta_nuc_init(init);
        return;
    }

    if (bind->mode == FASTA_NUC_MODE_BED) {
        init->bed_fp = hts_open(bind->bed_path, "r");
        if (!init->bed_fp) {
            duckdb_init_set_error(info, "fasta_nuc: failed to open BED file");
            destroy_fasta_nuc_init(init);
            return;
        }
        if (bind->region) {
            init->bed_tbx = tbx_index_load3(bind->bed_path, bind->bed_index_path, HTS_IDX_SILENT_FAIL);
            if (init->bed_tbx) {
                init->bed_itr = tbx_itr_querys(init->bed_tbx, bind->region);
                if (!init->bed_itr) {
                    duckdb_init_set_error(info, "fasta_nuc: failed to create BED region iterator");
                    destroy_fasta_nuc_init(init);
                    return;
                }
            }
        }
    } else {
        if (init->has_region) {
            init->current_tid = init->region_tid;
            init->next_bin_start = init->region_beg;
            init->current_seq_end_exclusive = init->region_end;
        } else {
            init->current_tid = 0;
            init->next_bin_start = 0;
            const char *seqname = faidx_iseq(init->fai, 0);
            init->current_seq_end_exclusive = faidx_seq_len64(init->fai, seqname);
        }
    }

    init->n_projected_cols = duckdb_init_get_column_count(info);
    init->column_ids = (idx_t *)duckdb_malloc(sizeof(idx_t) * init->n_projected_cols);
    for (idx_t i = 0; i < init->n_projected_cols; i++) {
        init->column_ids[i] = duckdb_init_get_column_index(info, i);
    }
    duckdb_init_set_init_data(info, init, destroy_fasta_nuc_init);
}

static void count_nucleotides(const char *seq, hts_pos_t len,
                              int64_t *a, int64_t *c, int64_t *g, int64_t *t,
                              int64_t *n, int64_t *other) {
    *a = *c = *g = *t = *n = *other = 0;
    for (hts_pos_t i = 0; i < len; i++) {
        switch (toupper((unsigned char)seq[i])) {
            case 'A': (*a)++; break;
            case 'C': (*c)++; break;
            case 'G': (*g)++; break;
            case 'T': (*t)++; break;
            case 'N': (*n)++; break;
            default: (*other)++; break;
        }
    }
}

static bool bed_overlap_region(const fasta_nuc_init_data_t *init, const char *chrom, int64_t start, int64_t end) {
    if (!init->has_region) return true;
    if (!chrom || strcmp(chrom, init->region_seq) != 0) return false;
    return end > (int64_t)init->region_beg && start < (int64_t)init->region_end;
}

static int next_fasta_nuc_bed_interval(fasta_nuc_init_data_t *init, const char **chrom, int64_t *start, int64_t *end) {
    while (!init->bed_finished) {
        int ret = init->bed_itr ? tbx_itr_next(init->bed_fp, init->bed_tbx, init->bed_itr, &init->bed_line)
                                : hts_getline(init->bed_fp, '\n', &init->bed_line);
        if (ret < 0) {
            init->bed_finished = true;
            return 0;
        }
        if (init->bed_line.l == 0 || is_meta_bed_line(init->bed_line.s)) continue;
        int len0 = 0, len1 = 0, len2 = 0;
        const char *f0 = get_field_span(init->bed_line.s, 0, &len0);
        const char *f1 = get_field_span(init->bed_line.s, 1, &len1);
        const char *f2 = get_field_span(init->bed_line.s, 2, &len2);
        int64_t s = 0, e = 0;
        if (!f0 || !f1 || !f2 || !parse_int64_span_local(f1, len1, &s) || !parse_int64_span_local(f2, len2, &e)) {
            continue;
        }
        char *chrom_buf = (char *)malloc((size_t)len0 + 1);
        if (!chrom_buf) return 0;
        memcpy(chrom_buf, f0, (size_t)len0);
        chrom_buf[len0] = '\0';
        if (!bed_overlap_region(init, chrom_buf, s, e)) {
            free(chrom_buf);
            continue;
        }
        *chrom = chrom_buf;
        *start = s;
        *end = e;
        return 1;
    }
    return 0;
}

static int next_fasta_nuc_bin_interval(fasta_nuc_init_data_t *init, const char **chrom, int64_t *start, int64_t *end) {
    int nseq = faidx_nseq(init->fai);
    while (!init->done) {
        if (init->current_tid >= nseq) {
            init->done = true;
            return 0;
        }
        const char *seqname = faidx_iseq(init->fai, init->current_tid);
        if (!seqname) {
            init->done = true;
            return 0;
        }
        hts_pos_t seq_len = init->has_region && init->current_tid == init->region_tid
            ? init->current_seq_end_exclusive
            : faidx_seq_len64(init->fai, seqname);
        if (seq_len < 0) {
            init->done = true;
            return 0;
        }
        if (init->next_bin_start >= seq_len) {
            if (init->has_region) {
                init->done = true;
                return 0;
            }
            init->current_tid++;
            if (init->current_tid >= nseq) {
                init->done = true;
                return 0;
            }
            init->next_bin_start = 0;
            init->current_seq_end_exclusive = faidx_seq_len64(init->fai, faidx_iseq(init->fai, init->current_tid));
            continue;
        }
        hts_pos_t bin_end = init->next_bin_start + init->bind->bin_width;
        if (bin_end > seq_len) bin_end = seq_len;
        *chrom = seqname;
        *start = init->next_bin_start;
        *end = bin_end;
        init->next_bin_start = bin_end;
        return 1;
    }
    return 0;
}

static void fasta_nuc_scan(duckdb_function_info info, duckdb_data_chunk output) {
    fasta_nuc_init_data_t *init = (fasta_nuc_init_data_t *)duckdb_function_get_init_data(info);
    if (!init || init->done) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    idx_t row_count = 0;
    idx_t col_count = duckdb_data_chunk_get_column_count(output);
    duckdb_vector vectors[NUC_COL_COUNT];
    for (idx_t c = 0; c < col_count; c++) {
        vectors[c] = duckdb_data_chunk_get_vector(output, c);
    }

    while (row_count < INTERVAL_BATCH_SIZE) {
        const char *chrom = NULL;
        int64_t start = 0, end = 0;
        int have_interval = init->bind->mode == FASTA_NUC_MODE_BED
            ? next_fasta_nuc_bed_interval(init, &chrom, &start, &end)
            : next_fasta_nuc_bin_interval(init, &chrom, &start, &end);
        if (!have_interval) break;

        hts_pos_t seq_len = (hts_pos_t)(end - start);
        char *seq = NULL;
        hts_pos_t fetch_len = 0;
        int64_t num_a = 0, num_c = 0, num_g = 0, num_t = 0, num_n = 0, num_other = 0;
        double pct_at = 0.0, pct_gc = 0.0;

        if (seq_len > 0) {
            seq = faidx_fetch_seq64(init->fai, chrom, (hts_pos_t)start, (hts_pos_t)end - 1, &fetch_len);
            if (!seq || fetch_len < 0) {
                if (init->bind->mode == FASTA_NUC_MODE_BED) free((void *)chrom);
                free(seq);
                continue;
            }
            seq_len = fetch_len;
            count_nucleotides(seq, seq_len, &num_a, &num_c, &num_g, &num_t, &num_n, &num_other);
            if (seq_len > 0) {
                pct_at = (double)(num_a + num_t) / (double)seq_len;
                pct_gc = (double)(num_c + num_g) / (double)seq_len;
            }
        }

        for (idx_t c = 0; c < col_count; c++) {
            int logical_col = (int)init->column_ids[c];
            switch (logical_col) {
                case NUC_COL_CHROM:
                    duckdb_vector_assign_string_element(vectors[c], row_count, chrom);
                    break;
                case NUC_COL_START: {
                    int64_t *data = (int64_t *)duckdb_vector_get_data(vectors[c]);
                    data[row_count] = start;
                    break;
                }
                case NUC_COL_END: {
                    int64_t *data = (int64_t *)duckdb_vector_get_data(vectors[c]);
                    data[row_count] = end;
                    break;
                }
                case NUC_COL_PCT_AT: {
                    double *data = (double *)duckdb_vector_get_data(vectors[c]);
                    data[row_count] = pct_at;
                    break;
                }
                case NUC_COL_PCT_GC: {
                    double *data = (double *)duckdb_vector_get_data(vectors[c]);
                    data[row_count] = pct_gc;
                    break;
                }
                case NUC_COL_NUM_A:
                case NUC_COL_NUM_C:
                case NUC_COL_NUM_G:
                case NUC_COL_NUM_T:
                case NUC_COL_NUM_N:
                case NUC_COL_NUM_OTHER:
                case NUC_COL_SEQ_LEN: {
                    int64_t *data = (int64_t *)duckdb_vector_get_data(vectors[c]);
                    data[row_count] =
                        logical_col == NUC_COL_NUM_A ? num_a :
                        logical_col == NUC_COL_NUM_C ? num_c :
                        logical_col == NUC_COL_NUM_G ? num_g :
                        logical_col == NUC_COL_NUM_T ? num_t :
                        logical_col == NUC_COL_NUM_N ? num_n :
                        logical_col == NUC_COL_NUM_OTHER ? num_other :
                        (int64_t)seq_len;
                    break;
                }
                case NUC_COL_SEQ:
                    if (init->bind->include_seq && seq) {
                        duckdb_vector_assign_string_element_len(vectors[c], row_count, seq, (idx_t)seq_len);
                    } else {
                        set_null(vectors[c], row_count);
                    }
                    break;
                default:
                    set_null(vectors[c], row_count);
                    break;
            }
        }

        if (init->bind->mode == FASTA_NUC_MODE_BED) free((void *)chrom);
        free(seq);
        row_count++;
    }

    if (row_count == 0 && init->bind->mode == FASTA_NUC_MODE_BINS) init->done = true;
    if (row_count == 0 && init->bind->mode == FASTA_NUC_MODE_BED && init->bed_finished) init->done = true;
    duckdb_data_chunk_set_size(output, row_count);
}

void register_read_bed_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "read_bed");
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "region", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_table_function_set_bind(tf, read_bed_bind);
    duckdb_table_function_set_init(tf, read_bed_init);
    duckdb_table_function_set_function(tf, read_bed_scan);
    duckdb_table_function_supports_projection_pushdown(tf, true);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}

void register_fasta_nuc_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "fasta_nuc");
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "bed_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "bin_width", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "region", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "bed_index_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "include_seq", bool_type);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_table_function_set_bind(tf, fasta_nuc_bind);
    duckdb_table_function_set_init(tf, fasta_nuc_init);
    duckdb_table_function_set_function(tf, fasta_nuc_scan);
    duckdb_table_function_supports_projection_pushdown(tf, true);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}
