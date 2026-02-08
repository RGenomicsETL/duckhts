/**
 * DuckHTS Tabix/GTF/GFF Reader
 *
 * Table functions for reading tabix-indexed tab-delimited files,
 * with specialised parsers for GTF and GFF3 formats.
 *
 * read_tabix(path, [region])  → generic: splits lines by \t into columns
 * read_gtf(path, [region])    → GTF-aware: typed SEQNAME..ATTRIBUTES + parsed attrs
 * read_gff(path, [region])    → GFF3-aware: same structure as GTF
 *
 * GTF/GFF columns (1-indexed in spec):
 *   1. seqname   VARCHAR
 *   2. source    VARCHAR
 *   3. feature   VARCHAR
 *   4. start     INTEGER
 *   5. end       INTEGER
 *   6. score     DOUBLE
 *   7. strand    VARCHAR
 *   8. frame     VARCHAR
 *   9. attributes VARCHAR  (raw attribute string)
 *
 * All three functions support an optional 'region' named parameter for
 * tabix-indexed queries (e.g. region := 'chr1:1000-2000').
 * Without an index, the file is scanned sequentially.
 *
 * API reference: htslib tbx.h, hts_getline()
 */

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>

/* ================================================================
 * Constants
 * ================================================================ */

#define TABIX_BATCH_SIZE 2048
#define TABIX_MAX_GENERIC_COLS 256

/* GTF/GFF column indices */
enum {
    GXF_COL_SEQNAME = 0,
    GXF_COL_SOURCE,
    GXF_COL_FEATURE,
    GXF_COL_START,
    GXF_COL_END,
    GXF_COL_SCORE,
    GXF_COL_STRAND,
    GXF_COL_FRAME,
    GXF_COL_ATTRIBUTES,
    GXF_COL_COUNT
};

/* Reader mode */
typedef enum {
    TABIX_MODE_GENERIC = 0,
    TABIX_MODE_GTF,
    TABIX_MODE_GFF
} tabix_mode_t;

/* ================================================================
 * Bind Data
 * ================================================================ */

typedef struct {
    char *file_path;
    char *region;          /* NULL = full scan */
    tabix_mode_t mode;
    int  n_cols;           /* detected or fixed (9 for GTF/GFF) */
} tabix_bind_data_t;

static void tabix_bind_data_destroy(void *data) {
    tabix_bind_data_t *bd = (tabix_bind_data_t *)data;
    if (bd) {
        free(bd->file_path);
        free(bd->region);
        free(bd);
    }
}

/* ================================================================
 * Init Data  (single-threaded scan for tabix)
 * ================================================================ */

typedef struct {
    htsFile *fp;
    tbx_t   *tbx;          /* NULL when file is not indexed */
    hts_itr_t *itr;        /* NULL when doing sequential scan */
    kstring_t  line;
    bool       finished;
} tabix_init_data_t;

static void tabix_init_data_destroy(void *data) {
    tabix_init_data_t *id = (tabix_init_data_t *)data;
    if (id) {
        if (id->itr) hts_itr_destroy(id->itr);
        if (id->tbx) tbx_destroy(id->tbx);
        if (id->fp)  hts_close(id->fp);
        free(id->line.s);
        free(id);
    }
}

/* ================================================================
 * Helpers
 * ================================================================ */

/* Count tab-separated fields in a line */
static int count_fields(const char *s) {
    int n = 1;
    while (*s) {
        if (*s == '\t') n++;
        s++;
    }
    return n;
}

/* Get n-th tab-separated field (0-based). Returns pointer into s, sets *len.
 * Returns NULL if field index is out of range. */
static const char *get_field(const char *s, int idx, int *len) {
    int cur = 0;
    const char *start = s;
    while (*s) {
        if (*s == '\t') {
            if (cur == idx) {
                *len = (int)(s - start);
                return start;
            }
            cur++;
            start = s + 1;
        }
        s++;
    }
    if (cur == idx) {
        *len = (int)(s - start);
        return start;
    }
    *len = 0;
    return NULL;
}

/* ================================================================
 * Bind  (shared logic for all three modes)
 * ================================================================ */

static void tabix_bind(duckdb_bind_info info, tabix_mode_t mode) {
    tabix_bind_data_t *bd = calloc(1, sizeof(tabix_bind_data_t));
    if (!bd) {
        duckdb_bind_set_error(info, "Out of memory");
        return;
    }
    bd->mode = mode;

    /* positional param: file path */
    duckdb_value val = duckdb_bind_get_parameter(info, 0);
    if (duckdb_get_type_id(duckdb_get_value_type(val)) == DUCKDB_TYPE_VARCHAR) {
        const char *path = duckdb_get_varchar(val);
        bd->file_path = strdup(path);
        duckdb_free((void *)path);
    }
    duckdb_destroy_value(&val);

    if (!bd->file_path || bd->file_path[0] == '\0') {
        const char *names[] = {"read_tabix", "read_gtf", "read_gff"};
        char msg[128];
        snprintf(msg, sizeof(msg), "%s requires a file path", names[mode]);
        duckdb_bind_set_error(info, msg);
        tabix_bind_data_destroy(bd);
        return;
    }

    /* named param: region */
    val = duckdb_bind_get_named_parameter(info, "region");
    if (val) {
        if (duckdb_get_type_id(duckdb_get_value_type(val)) == DUCKDB_TYPE_VARCHAR) {
            const char *r = duckdb_get_varchar(val);
            if (r && r[0]) bd->region = strdup(r);
            duckdb_free((void *)r);
        }
        duckdb_destroy_value(&val);
    }

    /* Schema: GTF/GFF have fixed 9 columns; generic auto-detects */
    if (mode == TABIX_MODE_GTF || mode == TABIX_MODE_GFF) {
        bd->n_cols = GXF_COL_COUNT;
        duckdb_bind_add_result_column(info, "seqname",    duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR));
        duckdb_bind_add_result_column(info, "source",     duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR));
        duckdb_bind_add_result_column(info, "feature",    duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR));
        duckdb_bind_add_result_column(info, "start",      duckdb_create_logical_type(DUCKDB_TYPE_BIGINT));
        duckdb_bind_add_result_column(info, "end",        duckdb_create_logical_type(DUCKDB_TYPE_BIGINT));
        duckdb_bind_add_result_column(info, "score",      duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE));
        duckdb_bind_add_result_column(info, "strand",     duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR));
        duckdb_bind_add_result_column(info, "frame",      duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR));
        duckdb_bind_add_result_column(info, "attributes", duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR));
    } else {
        /* Generic tabix: peek at first line to determine column count */
        htsFile *fp = hts_open(bd->file_path, "r");
        if (!fp) {
            duckdb_bind_set_error(info, "Cannot open file");
            tabix_bind_data_destroy(bd);
            return;
        }
        kstring_t line = {0, 0, NULL};
        int n_cols = 0;

        /* Skip comment/header lines (starting with # or meta_char) */
        while (hts_getline(fp, '\n', &line) >= 0) {
            if (line.l > 0 && line.s[0] != '#') {
                n_cols = count_fields(line.s);
                break;
            }
        }
        free(line.s);
        hts_close(fp);

        if (n_cols == 0) n_cols = 1;
        if (n_cols > TABIX_MAX_GENERIC_COLS) n_cols = TABIX_MAX_GENERIC_COLS;
        bd->n_cols = n_cols;

        duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
        for (int i = 0; i < n_cols; i++) {
            char col_name[32];
            snprintf(col_name, sizeof(col_name), "column%d", i);
            duckdb_bind_add_result_column(info, col_name, varchar_type);
        }
        duckdb_destroy_logical_type(&varchar_type);
    }

    duckdb_bind_set_bind_data(info, bd, tabix_bind_data_destroy);
}

static void read_tabix_bind(duckdb_bind_info info) { tabix_bind(info, TABIX_MODE_GENERIC); }
static void read_gtf_bind(duckdb_bind_info info)   { tabix_bind(info, TABIX_MODE_GTF); }
static void read_gff_bind(duckdb_bind_info info)    { tabix_bind(info, TABIX_MODE_GFF); }

/* ================================================================
 * Init
 * ================================================================ */

static void tabix_init(duckdb_init_info info) {
    tabix_bind_data_t *bd = (tabix_bind_data_t *)duckdb_init_get_bind_data(info);
    tabix_init_data_t *id = calloc(1, sizeof(tabix_init_data_t));
    if (!id) {
        duckdb_init_set_error(info, "Out of memory");
        return;
    }

    id->fp = hts_open(bd->file_path, "r");
    if (!id->fp) {
        char msg[512];
        snprintf(msg, sizeof(msg), "Cannot open file: %s", bd->file_path);
        duckdb_init_set_error(info, msg);
        free(id);
        return;
    }

    /* Try to load tabix index */
    id->tbx = tbx_index_load(bd->file_path);

    if (bd->region && bd->region[0]) {
        if (!id->tbx) {
            char msg[512];
            snprintf(msg, sizeof(msg),
                     "Region query requested but no tabix index found for: %s",
                     bd->file_path);
            duckdb_init_set_error(info, msg);
            hts_close(id->fp);
            free(id);
            return;
        }
        id->itr = tbx_itr_querys(id->tbx, bd->region);
        if (!id->itr) {
            /* Region doesn't match any sequences – return empty result */
            id->finished = true;
        }
    }
    /* else: sequential scan (no iterator) */

    id->line.l = 0;
    id->line.m = 0;
    id->line.s = NULL;

    duckdb_init_set_init_data(info, id, tabix_init_data_destroy);
}

/* ================================================================
 * Scan function
 * ================================================================ */

static void tabix_scan(duckdb_function_info info, duckdb_data_chunk output) {
    tabix_bind_data_t *bd = (tabix_bind_data_t *)duckdb_function_get_bind_data(info);
    tabix_init_data_t *id = (tabix_init_data_t *)duckdb_function_get_init_data(info);

    if (id->finished) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    idx_t row_count = 0;
    int n_cols = bd->n_cols;

    /* Pre-fetch column vectors */
    duckdb_vector *vectors = (duckdb_vector *)alloca(sizeof(duckdb_vector) * n_cols);
    for (int c = 0; c < n_cols; c++) {
        vectors[c] = duckdb_data_chunk_get_vector(output, c);
    }

    while (row_count < TABIX_BATCH_SIZE) {
        int ret;
        if (id->itr) {
            ret = tbx_itr_next(id->fp, id->tbx, id->itr, &id->line);
        } else {
            ret = hts_getline(id->fp, '\n', &id->line);
        }

        if (ret < 0) {
            id->finished = true;
            break;
        }

        /* Skip comment/header lines */
        if (id->line.l == 0 || id->line.s[0] == '#') continue;

        if (bd->mode == TABIX_MODE_GTF || bd->mode == TABIX_MODE_GFF) {
            /* Parse GTF/GFF 9 columns with proper types */
            for (int c = 0; c < GXF_COL_COUNT; c++) {
                int flen = 0;
                const char *fld = get_field(id->line.s, c, &flen);

                if (!fld || flen == 0 || (flen == 1 && fld[0] == '.')) {
                    /* Missing value */
                    switch (c) {
                        case GXF_COL_START:
                        case GXF_COL_END:
                            duckdb_vector_assign_string_element_len(vectors[c], row_count, "0", 1);
                            break;
                        case GXF_COL_SCORE: {
                            /* Set NULL for missing score */
                            uint64_t *validity = duckdb_vector_get_validity(vectors[c]);
                            if (validity) {
                                duckdb_validity_set_row_invalid(validity, row_count);
                            }
                            break;
                        }
                        default:
                            duckdb_vector_assign_string_element_len(vectors[c], row_count, ".", 1);
                            break;
                    }
                    continue;
                }

                switch (c) {
                    case GXF_COL_START:
                    case GXF_COL_END: {
                        /* Parse as int64 */
                        char buf[32];
                        int copy_len = flen < 31 ? flen : 31;
                        memcpy(buf, fld, copy_len);
                        buf[copy_len] = '\0';
                        int64_t ival = strtoll(buf, NULL, 10);
                        int64_t *data = (int64_t *)duckdb_vector_get_data(vectors[c]);
                        data[row_count] = ival;
                        break;
                    }
                    case GXF_COL_SCORE: {
                        /* Parse as double, '.' = NULL */
                        char buf[64];
                        int copy_len = flen < 63 ? flen : 63;
                        memcpy(buf, fld, copy_len);
                        buf[copy_len] = '\0';
                        double dval = strtod(buf, NULL);
                        double *data = (double *)duckdb_vector_get_data(vectors[c]);
                        data[row_count] = dval;
                        break;
                    }
                    default:
                        /* VARCHAR columns */
                        duckdb_vector_assign_string_element_len(vectors[c], row_count, fld, flen);
                        break;
                }
            }
        } else {
            /* Generic mode: all VARCHAR */
            for (int c = 0; c < n_cols; c++) {
                int flen = 0;
                const char *fld = get_field(id->line.s, c, &flen);
                if (fld) {
                    duckdb_vector_assign_string_element_len(vectors[c], row_count, fld, flen);
                } else {
                    duckdb_vector_assign_string_element_len(vectors[c], row_count, "", 0);
                }
            }
        }

        row_count++;
    }

    duckdb_data_chunk_set_size(output, row_count);
}

/* ================================================================
 * Registration helpers
 * ================================================================ */

static duckdb_table_function create_tabix_tf(const char *name,
                                              void (*bind_fn)(duckdb_bind_info)) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, name);

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "region", varchar_type);
    duckdb_destroy_logical_type(&varchar_type);

    duckdb_table_function_set_bind(tf, bind_fn);
    duckdb_table_function_set_init(tf, tabix_init);
    duckdb_table_function_set_function(tf, tabix_scan);
    duckdb_table_function_supports_projection_pushdown(tf, true);

    return tf;
}

void register_read_tabix_function(duckdb_connection connection) {
    duckdb_table_function tf = create_tabix_tf("read_tabix", read_tabix_bind);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}

void register_read_gtf_function(duckdb_connection connection) {
    duckdb_table_function tf = create_tabix_tf("read_gtf", read_gtf_bind);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}

void register_read_gff_function(duckdb_connection connection) {
    duckdb_table_function tf = create_tabix_tf("read_gff", read_gff_bind);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}
