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
    GXF_COL_ATTRIBUTES_MAP,
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
    int  include_attr_map;
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
    idx_t     *column_ids;      /* logical column indices (for projection pushdown) */
    idx_t      n_projected_cols;
} tabix_init_data_t;

static void tabix_init_data_destroy(void *data) {
    tabix_init_data_t *id = (tabix_init_data_t *)data;
    if (id) {
        if (id->itr) hts_itr_destroy(id->itr);
        if (id->tbx) tbx_destroy(id->tbx);
        if (id->fp)  hts_close(id->fp);
        free(id->line.s);
        free(id->column_ids);
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

static void trim_span(const char **start, int *len) {
    const char *s = *start;
    int l = *len;
    while (l > 0 && (*s == ' ' || *s == '\t')) {
        s++;
        l--;
    }
    while (l > 0 && (s[l - 1] == ' ' || s[l - 1] == '\t')) {
        l--;
    }
    *start = s;
    *len = l;
}

static int count_gff_pairs(const char *s) {
    int count = 0;
    const char *p = s;
    while (*p) {
        while (*p == ';' || *p == ' ' || *p == '\t') p++;
        if (!*p) break;
        const char *key = p;
        while (*p && *p != '=' && *p != ';') p++;
        if (*p != '=') {
            while (*p && *p != ';') p++;
            continue;
        }
        int key_len = (int)(p - key);
        p++; /* skip '=' */
        const char *val = p;
        while (*p && *p != ';') p++;
        int val_len = (int)(p - val);
        trim_span(&key, &key_len);
        trim_span(&val, &val_len);
        if (key_len > 0) count++;
        if (*p == ';') p++;
    }
    return count;
}

static int count_gtf_pairs(const char *s) {
    int count = 0;
    const char *p = s;
    while (*p) {
        while (*p == ';' || *p == ' ' || *p == '\t') p++;
        if (!*p) break;
        const char *key = p;
        while (*p && *p != ' ' && *p != '\t' && *p != ';') p++;
        int key_len = (int)(p - key);
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '"') {
            p++;
            while (*p && *p != '"') p++;
            if (*p == '"') p++;
        } else {
            while (*p && *p != ';') p++;
        }
        trim_span(&key, &key_len);
        if (key_len > 0) count++;
        while (*p && *p != ';') p++;
        if (*p == ';') p++;
    }
    return count;
}

static void fill_attr_map(duckdb_vector vec, idx_t row, const char *s, bool is_gff) {
    if (!s || s[0] == '\0' || (s[0] == '.' && s[1] == '\0')) {
        duckdb_vector_ensure_validity_writable(vec);
        uint64_t *validity = duckdb_vector_get_validity(vec);
        duckdb_validity_set_row_invalid(validity, row);
        duckdb_list_entry entry = {duckdb_list_vector_get_size(vec), 0};
        duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(vec);
        list_data[row] = entry;
        return;
    }

    int pair_count = is_gff ? count_gff_pairs(s) : count_gtf_pairs(s);
    duckdb_list_entry entry;
    entry.offset = duckdb_list_vector_get_size(vec);
    entry.length = pair_count;

    if (pair_count == 0) {
        duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(vec);
        list_data[row] = entry;
        return;
    }

    duckdb_list_vector_reserve(vec, entry.offset + entry.length);
    duckdb_list_vector_set_size(vec, entry.offset + entry.length);

    duckdb_vector child = duckdb_list_vector_get_child(vec);
    duckdb_vector key_vec = duckdb_struct_vector_get_child(child, 0);
    duckdb_vector val_vec = duckdb_struct_vector_get_child(child, 1);

    const char *p = s;
    int write_idx = 0;
    while (*p && write_idx < pair_count) {
        while (*p == ';' || *p == ' ' || *p == '\t') p++;
        if (!*p) break;

        const char *key = p;
        int key_len = 0;
        const char *val = NULL;
        int val_len = 0;

        if (is_gff) {
            while (*p && *p != '=' && *p != ';') p++;
            if (*p != '=') {
                while (*p && *p != ';') p++;
                continue;
            }
            key_len = (int)(p - key);
            p++; /* '=' */
            val = p;
            while (*p && *p != ';') p++;
            val_len = (int)(p - val);
        } else {
            while (*p && *p != ' ' && *p != '\t' && *p != ';') p++;
            key_len = (int)(p - key);
            while (*p == ' ' || *p == '\t') p++;
            if (*p == '"') {
                p++;
                val = p;
                while (*p && *p != '"') p++;
                val_len = (int)(p - val);
                if (*p == '"') p++;
            } else {
                val = p;
                while (*p && *p != ';') p++;
                val_len = (int)(p - val);
            }
        }

        trim_span(&key, &key_len);
        trim_span(&val, &val_len);
        if (key_len > 0) {
            duckdb_vector_assign_string_element_len(key_vec, entry.offset + write_idx, key, key_len);
            duckdb_vector_assign_string_element_len(val_vec, entry.offset + write_idx, val, val_len);
            write_idx++;
        }
        while (*p && *p != ';') p++;
        if (*p == ';') p++;
    }

    entry.length = write_idx;
    duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(vec);
    list_data[row] = entry;
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
        /* File has 9 fixed columns; attributes_map is a derived column. */
        bd->n_cols = GXF_COL_COUNT - 1;
        duckdb_value attr_map_val = duckdb_bind_get_named_parameter(info, "attributes_map");
        if (attr_map_val && !duckdb_is_null_value(attr_map_val)) {
            bd->include_attr_map = duckdb_get_bool(attr_map_val) ? 1 : 0;
        }
        if (attr_map_val) duckdb_destroy_value(&attr_map_val);

        duckdb_bind_add_result_column(info, "seqname",    duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR));
        duckdb_bind_add_result_column(info, "source",     duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR));
        duckdb_bind_add_result_column(info, "feature",    duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR));
        duckdb_bind_add_result_column(info, "start",      duckdb_create_logical_type(DUCKDB_TYPE_BIGINT));
        duckdb_bind_add_result_column(info, "end",        duckdb_create_logical_type(DUCKDB_TYPE_BIGINT));
        duckdb_bind_add_result_column(info, "score",      duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE));
        duckdb_bind_add_result_column(info, "strand",     duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR));
        duckdb_bind_add_result_column(info, "frame",      duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR));
        duckdb_bind_add_result_column(info, "attributes", duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR));
        if (bd->include_attr_map) {
            duckdb_logical_type key_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
            duckdb_logical_type val_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
            duckdb_logical_type map_type = duckdb_create_map_type(key_type, val_type);
            duckdb_bind_add_result_column(info, "attributes_map", map_type);
            duckdb_destroy_logical_type(&key_type);
            duckdb_destroy_logical_type(&val_type);
            duckdb_destroy_logical_type(&map_type);
        }
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

    /* Store projection pushdown column mapping */
    id->n_projected_cols = duckdb_init_get_column_count(info);
    if (id->n_projected_cols > 0) {
        id->column_ids = (idx_t *)malloc(sizeof(idx_t) * id->n_projected_cols);
        for (idx_t i = 0; i < id->n_projected_cols; i++) {
            id->column_ids[i] = duckdb_init_get_column_index(info, i);
        }
    } else {
        id->column_ids = NULL;
    }

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
    idx_t chunk_col_count = duckdb_data_chunk_get_column_count(output);

    /* Pre-fetch column vectors */
    duckdb_vector *vectors = NULL;
    if (chunk_col_count > 0) {
        vectors = (duckdb_vector *)alloca(sizeof(duckdb_vector) * chunk_col_count);
        for (idx_t c = 0; c < chunk_col_count; c++) {
            vectors[c] = duckdb_data_chunk_get_vector(output, c);
        }
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

        if (chunk_col_count > 0) {
        if (bd->mode == TABIX_MODE_GTF || bd->mode == TABIX_MODE_GFF) {
            /* Parse GTF/GFF columns with projection pushdown:
             * Vector index c maps to logical column id->column_ids[c],
             * which is the field index in the TSV line. */
            for (idx_t c = 0; c < chunk_col_count; c++) {
                int logical_col = (int)id->column_ids[c];
                if (logical_col == GXF_COL_ATTRIBUTES_MAP && bd->include_attr_map) {
                    const char *fld = NULL;
                    int flen = 0;
                    fld = get_field(id->line.s, GXF_COL_ATTRIBUTES, &flen);
                    char *tmp = NULL;
                    if (fld && flen > 0) {
                        tmp = (char *)alloca((size_t)flen + 1);
                        memcpy(tmp, fld, (size_t)flen);
                        tmp[flen] = '\0';
                    }
                    fill_attr_map(vectors[c], row_count, tmp ? tmp : ".", bd->mode == TABIX_MODE_GFF);
                    continue;
                }
                if (logical_col >= GXF_COL_COUNT - 1) continue;

                int flen = 0;
                const char *fld = get_field(id->line.s, logical_col, &flen);

                if (!fld || flen == 0 || (flen == 1 && fld[0] == '.')) {
                    /* Missing value */
                    switch (logical_col) {
                        case GXF_COL_START:
                        case GXF_COL_END: {
                            int64_t *data = (int64_t *)duckdb_vector_get_data(vectors[c]);
                            data[row_count] = 0;
                            break;
                        }
                        case GXF_COL_SCORE: {
                            /* Set NULL for missing score */
                            duckdb_vector_ensure_validity_writable(vectors[c]);
                            uint64_t *validity = duckdb_vector_get_validity(vectors[c]);
                            duckdb_validity_set_row_invalid(validity, row_count);
                            break;
                        }
                        default:
                            duckdb_vector_assign_string_element_len(vectors[c], row_count, ".", 1);
                            break;
                    }
                    continue;
                }

                switch (logical_col) {
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
            /* Generic mode: all VARCHAR, with projection pushdown */
            for (idx_t c = 0; c < chunk_col_count; c++) {
                int logical_col = (int)id->column_ids[c];
                if (logical_col >= n_cols) continue;
                int flen = 0;
                const char *fld = get_field(id->line.s, logical_col, &flen);
                if (fld) {
                    duckdb_vector_assign_string_element_len(vectors[c], row_count, fld, flen);
                } else {
                    duckdb_vector_assign_string_element_len(vectors[c], row_count, "", 0);
                }
            }
        }
        } /* chunk_col_count > 0 */

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

    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_table_function_add_named_parameter(tf, "attributes_map", bool_type);
    duckdb_destroy_logical_type(&bool_type);

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
