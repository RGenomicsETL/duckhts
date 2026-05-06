/**
 * DuckHTS Tabix/GTF/GFF Reader
 *
 * Table functions for reading tabix-indexed tab-delimited files,
 * with specialised parsers for GTF and GFF3 formats.
 *
 * read_tabix(path, [region, header, header_names])  → generic: splits lines by \t into columns
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
#include "include/vcf_types.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <strings.h>
#include <ctype.h>
#include <errno.h>

#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <htslib/kstring.h>

/* ================================================================
 * Constants
 * ================================================================ */

#define TABIX_BATCH_SIZE 2048
#define TABIX_MAX_GENERIC_COLS 256
#define TABIX_NUM_PARSE_BUF 128

static inline void set_null(duckdb_vector vec, idx_t row) {
    duckdb_vector_ensure_validity_writable(vec);
    uint64_t *v = duckdb_vector_get_validity(vec);
    v[row / 64] &= ~((uint64_t)1 << (row % 64));
}

static int parse_int64_span(const char *s, int len, int64_t *out) {
    if (!s || len <= 0 || len >= TABIX_NUM_PARSE_BUF) return 0;
    char buf[TABIX_NUM_PARSE_BUF];
    memcpy(buf, s, (size_t)len);
    buf[len] = '\0';
    char *end = NULL;
    errno = 0;
    int64_t v = strtoll(buf, &end, 10);
    if (errno == ERANGE || !end || *end != '\0') return 0;
    *out = v;
    return 1;
}

static int parse_double_span(const char *s, int len, double *out) {
    if (!s || len <= 0 || len >= TABIX_NUM_PARSE_BUF) return 0;
    char buf[TABIX_NUM_PARSE_BUF];
    memcpy(buf, s, (size_t)len);
    buf[len] = '\0';
    char *end = NULL;
    double v = strtod(buf, &end);
    if (!end || *end != '\0') return 0;
    *out = v;
    return 1;
}

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
    GXF_BASE_COL_COUNT
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
    char *index_path;
    char *region;          /* NULL = full scan */
    char **regions;        /* parsed comma-separated regions */
    unsigned int n_regions;
    tabix_mode_t mode;
    int  n_cols;           /* detected or fixed (9 for GTF/GFF) */
    int  include_attr_map;
    int  include_attr_list;
    int  include_attr_pairs;
    int  attr_map_col;
    int  attr_list_col;
    int  attr_pairs_col;
    int  strict;
    int  header;
    int  skip_header_line;
    int  header_names_count;
    char **header_names;
    char meta_char;
    int  line_skip;
    int  auto_detect;
    int  col_types_provided;
    int *col_types;
    uint64_t index_row_count;
    int index_row_count_valid;
} tabix_bind_data_t;

static void tabix_bind_data_destroy(void *data) {
    tabix_bind_data_t *bd = (tabix_bind_data_t *)data;
    if (bd) {
        free(bd->file_path);
        free(bd->index_path);
        free(bd->region);
        if (bd->regions) {
            for (unsigned int i = 0; i < bd->n_regions; i++) {
                free(bd->regions[i]);
            }
            free(bd->regions);
        }
        if (bd->header_names) {
            for (int i = 0; i < bd->header_names_count; i++) {
                free(bd->header_names[i]);
            }
            free(bd->header_names);
        }
        free(bd->col_types);
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
    int        count_only;
    uint64_t   count_remaining;
    unsigned int next_region_idx;
    int        skip_remaining;
    int        skipped_header;
    uint64_t   line_number;
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

static int is_integer_field(const char *s, int len) {
    if (len == 0) return 0;
    int i = 0;
    if (s[i] == '-' || s[i] == '+') i++;
    if (i >= len) return 0;
    for (; i < len; i++) {
        if (s[i] < '0' || s[i] > '9') return 0;
    }
    return 1;
}

static int is_float_field(const char *s, int len) {
    if (len == 0) return 0;
    char *tmp = (char *)malloc((size_t)len + 1);
    if (!tmp) return 0;
    memcpy(tmp, s, (size_t)len);
    tmp[len] = '\0';
    char *end = NULL;
    strtod(tmp, &end);
    int ok = (end && *end == '\0');
    free(tmp);
    return ok;
}

static int parse_type_name(const char *s) {
    if (!s) return DUCKDB_TYPE_VARCHAR;
    if (strcasecmp(s, "INT") == 0 || strcasecmp(s, "INTEGER") == 0)
        return DUCKDB_TYPE_INTEGER;
    if (strcasecmp(s, "BIGINT") == 0 || strcasecmp(s, "LONG") == 0)
        return DUCKDB_TYPE_BIGINT;
    if (strcasecmp(s, "DOUBLE") == 0 || strcasecmp(s, "FLOAT") == 0 ||
        strcasecmp(s, "REAL") == 0)
        return DUCKDB_TYPE_DOUBLE;
    if (strcasecmp(s, "VARCHAR") == 0 || strcasecmp(s, "STRING") == 0)
        return DUCKDB_TYPE_VARCHAR;
    return DUCKDB_TYPE_VARCHAR;
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

static int span_equals_lit(const char *s, int len, const char *lit) {
    size_t lit_len = strlen(lit);
    return len == (int)lit_len && memcmp(s, lit, lit_len) == 0;
}

static int span_has_whitespace(const char *s, int len) {
    for (int i = 0; i < len; i++) {
        if (isspace((unsigned char)s[i])) return 1;
    }
    return 0;
}

static int gff3_attr_segments_valid(const char *s, int len) {
    if (!s || len <= 0) return 0;
    const char *p = s;
    const char *end = s + len;
    int saw_pair = 0;
    while (p < end) {
        while (p < end && (*p == ';' || *p == ' ' || *p == '\t')) p++;
        if (p >= end) break;

        const char *seg = p;
        while (p < end && *p != ';') p++;
        int seg_len = (int)(p - seg);
        trim_span(&seg, &seg_len);
        if (seg_len <= 0) {
            if (p < end && *p == ';') p++;
            continue;
        }

        const char *eq = memchr(seg, '=', (size_t)seg_len);
        if (!eq) return 0;
        const char *key = seg;
        int key_len = (int)(eq - seg);
        trim_span(&key, &key_len);
        if (key_len <= 0) return 0;
        saw_pair = 1;

        if (p < end && *p == ';') p++;
    }
    return saw_pair;
}

static int hex_value(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}

static char *dup_trimmed_span(const char *s, int len, int *out_len) {
    if (!s || len < 0) return NULL;
    trim_span(&s, &len);
    char *out = (char *)malloc((size_t)len + 1);
    if (!out) return NULL;
    if (len > 0) memcpy(out, s, (size_t)len);
    out[len] = '\0';
    if (out_len) *out_len = len;
    return out;
}

static char *percent_decode_span(const char *s, int len, int *out_len) {
    if (!s || len < 0) return NULL;
    trim_span(&s, &len);
    char *out = (char *)malloc((size_t)len + 1);
    if (!out) return NULL;
    int w = 0;
    for (int i = 0; i < len; i++) {
        if (s[i] == '%' && i + 2 < len) {
            int hi = hex_value(s[i + 1]);
            int lo = hex_value(s[i + 2]);
            if (hi >= 0 && lo >= 0) {
                out[w++] = (char)((hi << 4) | lo);
                i += 2;
                continue;
            }
        }
        out[w++] = s[i];
    }
    out[w] = '\0';
    if (out_len) *out_len = w;
    return out;
}

typedef struct {
    char *key;
    int key_len;
    char *value;
    int value_len;
    int idx;
} gxf_attr_pair_t;

static void free_attr_pairs(gxf_attr_pair_t *pairs, int n_pairs) {
    if (!pairs) return;
    for (int i = 0; i < n_pairs; i++) {
        free(pairs[i].key);
        free(pairs[i].value);
    }
    free(pairs);
}

static int attr_pair_next_idx(gxf_attr_pair_t *pairs, int n_pairs, const char *key, int key_len) {
    int idx = 0;
    for (int i = 0; i < n_pairs; i++) {
        if (pairs[i].key_len == key_len && memcmp(pairs[i].key, key, (size_t)key_len) == 0) idx++;
    }
    return idx;
}

static int append_attr_pair(gxf_attr_pair_t **pairs, int *n_pairs, int *cap,
                            const char *key, int key_len,
                            char *value, int value_len) {
    if (key_len <= 0) {
        free(value);
        return 1;
    }
    if (*n_pairs >= *cap) {
        int new_cap = *cap ? (*cap * 2) : 8;
        gxf_attr_pair_t *tmp = (gxf_attr_pair_t *)realloc(*pairs, sizeof(gxf_attr_pair_t) * (size_t)new_cap);
        if (!tmp) {
            free(value);
            return 0;
        }
        *pairs = tmp;
        *cap = new_cap;
    }
    int stored_key_len = 0;
    char *stored_key = dup_trimmed_span(key, key_len, &stored_key_len);
    if (!stored_key) {
        free(value);
        return 0;
    }
    int idx = attr_pair_next_idx(*pairs, *n_pairs, stored_key, stored_key_len);
    (*pairs)[*n_pairs].key = stored_key;
    (*pairs)[*n_pairs].key_len = stored_key_len;
    (*pairs)[*n_pairs].value = value;
    (*pairs)[*n_pairs].value_len = value_len;
    (*pairs)[*n_pairs].idx = idx;
    (*n_pairs)++;
    return 1;
}

static int parse_gff3_attr_pairs(const char *s, gxf_attr_pair_t **out_pairs, int *out_count) {
    *out_pairs = NULL;
    *out_count = 0;
    if (!s || s[0] == '\0' || (s[0] == '.' && s[1] == '\0')) return 1;

    gxf_attr_pair_t *pairs = NULL;
    int n_pairs = 0, cap = 0;
    const char *p = s;
    while (*p) {
        while (*p == ';' || *p == ' ' || *p == '\t') p++;
        if (!*p) break;
        const char *key = p;
        while (*p && *p != '=' && *p != ';') p++;
        if (*p != '=') {
            while (*p && *p != ';') p++;
            if (*p == ';') p++;
            continue;
        }
        int key_len = (int)(p - key);
        const char *key_trim = key;
        trim_span(&key_trim, &key_len);
        p++;
        const char *val = p;
        while (*p && *p != ';') p++;
        int val_len = (int)(p - val);
        const char *seg = val;
        const char *seg_end = val + val_len;
        while (seg <= seg_end) {
            const char *comma = seg;
            while (comma < seg_end && *comma != ',') comma++;
            int part_len = (int)(comma - seg);
            int decoded_len = 0;
            char *decoded = percent_decode_span(seg, part_len, &decoded_len);
            if (!decoded || !append_attr_pair(&pairs, &n_pairs, &cap, key_trim, key_len, decoded, decoded_len)) {
                free_attr_pairs(pairs, n_pairs);
                return 0;
            }
            if (comma >= seg_end) break;
            seg = comma + 1;
        }
        if (*p == ';') p++;
    }
    *out_pairs = pairs;
    *out_count = n_pairs;
    return 1;
}

static int parse_gtf_attr_pairs(const char *s, gxf_attr_pair_t **out_pairs, int *out_count) {
    *out_pairs = NULL;
    *out_count = 0;
    if (!s || s[0] == '\0' || (s[0] == '.' && s[1] == '\0')) return 1;

    gxf_attr_pair_t *pairs = NULL;
    int n_pairs = 0, cap = 0;
    const char *p = s;
    while (*p) {
        while (*p == ';' || *p == ' ' || *p == '\t') p++;
        if (!*p) break;
        const char *key = p;
        while (*p && *p != ' ' && *p != '\t' && *p != ';') p++;
        int key_len = (int)(p - key);
        const char *key_trim = key;
        trim_span(&key_trim, &key_len);
        while (*p == ' ' || *p == '\t') p++;
        const char *val = p;
        int val_len = 0;
        if (*p == '"') {
            p++;
            val = p;
            while (*p && *p != '"') p++;
            val_len = (int)(p - val);
            if (*p == '"') p++;
        } else {
            while (*p && *p != ';') p++;
            val_len = (int)(p - val);
        }
        int stored_len = 0;
        char *stored = dup_trimmed_span(val, val_len, &stored_len);
        if (!stored || !append_attr_pair(&pairs, &n_pairs, &cap, key_trim, key_len, stored, stored_len)) {
            free_attr_pairs(pairs, n_pairs);
            return 0;
        }
        while (*p && *p != ';') p++;
        if (*p == ';') p++;
    }
    *out_pairs = pairs;
    *out_count = n_pairs;
    return 1;
}

static int parse_attr_pairs(const char *s, bool is_gff, gxf_attr_pair_t **out_pairs, int *out_count) {
    return is_gff ? parse_gff3_attr_pairs(s, out_pairs, out_count)
                  : parse_gtf_attr_pairs(s, out_pairs, out_count);
}

static int validate_gff3_line_strict(const char *line, char *err, size_t err_len) {
    int n_fields = count_fields(line);
    if (n_fields != 9) {
        if (n_fields < 9) {
            snprintf(err, err_len, "TooFewFields: expected exactly 9 tab-separated fields, found %d", n_fields);
        } else {
            snprintf(err, err_len, "TooManyFields: expected exactly 9 tab-separated fields, found %d", n_fields);
        }
        return 0;
    }

    int len = 0;
    const char *seqid = get_field(line, GXF_COL_SEQNAME, &len);
    if (!seqid || len == 0) {
        snprintf(err, err_len, "EmptySeqid: seqid (column 1) is empty");
        return 0;
    }

    const char *feature = get_field(line, GXF_COL_FEATURE, &len);
    if (!feature || len == 0) {
        snprintf(err, err_len, "EmptyFeaturetype: feature type (column 3) is empty");
        return 0;
    }
    if (span_has_whitespace(feature, len)) {
        snprintf(err, err_len, "InvalidFeaturetype: feature type contains whitespace");
        return 0;
    }

    int start_len = 0, end_len = 0;
    const char *start_s = get_field(line, GXF_COL_START, &start_len);
    const char *end_s = get_field(line, GXF_COL_END, &end_len);
    int64_t start = 0, endv = 0;
    int have_start = start_s && start_len > 0 && !(start_len == 1 && start_s[0] == '.');
    int have_end = end_s && end_len > 0 && !(end_len == 1 && end_s[0] == '.');
    if (have_start) {
        if (!parse_int64_span(start_s, start_len, &start)) {
            snprintf(err, err_len, "InvalidCoordinate: start coordinate is not an integer");
            return 0;
        }
        if (start < 1) {
            snprintf(err, err_len, "InvalidCoordinate: start coordinate must be >= 1");
            return 0;
        }
    }
    if (have_end) {
        if (!parse_int64_span(end_s, end_len, &endv)) {
            snprintf(err, err_len, "InvalidCoordinate: end coordinate is not an integer");
            return 0;
        }
        if (endv < 1) {
            snprintf(err, err_len, "InvalidCoordinate: end coordinate must be >= 1");
            return 0;
        }
    }
    if (have_start && have_end && endv < start) {
        snprintf(err, err_len, "InvalidCoordinate: end coordinate is less than start coordinate");
        return 0;
    }

    const char *score = get_field(line, GXF_COL_SCORE, &len);
    if (score && len > 0 && !(len == 1 && score[0] == '.')) {
        double d = 0.0;
        if (!parse_double_span(score, len, &d)) {
            snprintf(err, err_len, "InvalidScore: score must be a float or '.'");
            return 0;
        }
    }

    const char *strand = get_field(line, GXF_COL_STRAND, &len);
    if (!strand || len != 1 || !(strand[0] == '+' || strand[0] == '-' || strand[0] == '?' || strand[0] == '.')) {
        snprintf(err, err_len, "InvalidStrand: strand must be one of '+', '-', '?', '.'");
        return 0;
    }

    int frame_len = 0;
    const char *frame = get_field(line, GXF_COL_FRAME, &frame_len);
    if (!frame || frame_len != 1 || !(frame[0] == '.' || frame[0] == '0' || frame[0] == '1' || frame[0] == '2')) {
        snprintf(err, err_len, "InvalidPhase: phase must be 0, 1, 2, or '.'");
        return 0;
    }
    int feature_len = 0;
    const char *feature_check = get_field(line, GXF_COL_FEATURE, &feature_len);
    if (feature_check && span_equals_lit(feature_check, feature_len, "CDS") && frame[0] == '.') {
        snprintf(err, err_len, "InvalidPhase: CDS row missing required phase");
        return 0;
    }

    int attr_len = 0;
    const char *attrs = get_field(line, GXF_COL_ATTRIBUTES, &attr_len);
    if (attrs) trim_span(&attrs, &attr_len);
    if (attrs && attr_len > 0 && !(attr_len == 1 && attrs[0] == '.')) {
        if (!gff3_attr_segments_valid(attrs, attr_len)) {
            snprintf(err, err_len, "InvalidAttribute: every non-empty attribute segment must be key=value");
            return 0;
        }
    }

    return 1;
}

static char *dup_field_name(const char *start, int len) {
    trim_span(&start, &len);
    char *name = (char *)malloc((size_t)len + 1);
    if (!name) return NULL;
    memcpy(name, start, (size_t)len);
    name[len] = '\0';
    return name;
}

static int parse_header_names(const char *line, char ***out_names) {
    int n = count_fields(line);
    char **names = (char **)malloc(sizeof(char *) * (size_t)n);
    if (!names) return 0;
    for (int i = 0; i < n; i++) {
        int len = 0;
        const char *start = get_field(line, i, &len);
        if (!start) {
            names[i] = dup_field_name("", 0);
        } else {
            names[i] = dup_field_name(start, len);
        }
        if (!names[i]) {
            for (int j = 0; j < i; j++) free(names[j]);
            free(names);
            return 0;
        }
    }
    *out_names = names;
    return n;
}

static void parse_regions(const char *region_str, char ***out_regions, unsigned int *out_count) {
    *out_regions = NULL;
    *out_count = 0;
    if (!region_str || region_str[0] == '\0') return;

    unsigned int count = 1;
    for (const char *p = region_str; *p; p++) {
        if (*p == ',') count++;
    }

    char **arr = (char **)malloc(sizeof(char *) * count);
    if (!arr) return;

    char *dup = strdup(region_str);
    if (!dup) {
        free(arr);
        return;
    }

    unsigned int idx = 0;
    char *tok = strtok(dup, ",");
    while (tok && idx < count) {
        while (*tok == ' ' || *tok == '\t') tok++;
        int len = (int)strlen(tok);
        while (len > 0 && (tok[len - 1] == ' ' || tok[len - 1] == '\t')) {
            tok[--len] = '\0';
        }
        if (len > 0) {
            arr[idx] = strdup(tok);
            if (!arr[idx]) {
                for (unsigned int i = 0; i < idx; i++) free(arr[i]);
                free(arr);
                free(dup);
                return;
            }
            idx++;
        }
        tok = strtok(NULL, ",");
    }

    free(dup);
    *out_regions = arr;
    *out_count = idx;
}

static int tabix_advance_region_iterator(tabix_init_data_t *id, tabix_bind_data_t *bd) {
    if (!id || !bd || !id->tbx || !bd->regions) return 0;

    if (id->itr) {
        hts_itr_destroy(id->itr);
        id->itr = NULL;
    }

    while (id->next_region_idx < bd->n_regions) {
        const char *region = bd->regions[id->next_region_idx++];
        id->itr = tbx_itr_querys(id->tbx, region);
        if (id->itr) return 1;
    }
    return 0;
}

static int tabix_try_get_index_row_count(tbx_t *tbx, uint64_t *out_total) {
    if (!tbx || !tbx->idx || !out_total) return 0;

    int n = 0;
    const char **names = tbx_seqnames(tbx, &n);
    uint64_t total = 0;
    for (int tid = 0; tid < n; tid++) {
        uint64_t mapped = 0, unmapped = 0;
        if (hts_idx_get_stat(tbx->idx, tid, &mapped, &unmapped) != 0) {
            free(names);
            return 0;
        }
        total += mapped + unmapped;
    }
    free(names);
    *out_total = total;
    return 1;
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

static int attr_group_index(char **keys, int *key_lens, int n_keys, const char *key, int key_len);

static int fill_attr_map(duckdb_vector vec, idx_t row, const char *s, bool is_gff) {
    if (!s || s[0] == '\0' || (s[0] == '.' && s[1] == '\0')) {
        duckdb_vector_ensure_validity_writable(vec);
        uint64_t *validity = duckdb_vector_get_validity(vec);
        duckdb_validity_set_row_invalid(validity, row);
        duckdb_list_entry entry = {duckdb_list_vector_get_size(vec), 0};
        duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(vec);
        list_data[row] = entry;
        return 1;
    }

    int pair_count = is_gff ? count_gff_pairs(s) : count_gtf_pairs(s);
    duckdb_list_entry entry;
    entry.offset = duckdb_list_vector_get_size(vec);
    entry.length = pair_count;

    if (pair_count == 0) {
        duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(vec);
        list_data[row] = entry;
        return 1;
    }

    if (duckdb_list_vector_reserve(vec, entry.offset + entry.length) == DuckDBError ||
        duckdb_list_vector_set_size(vec, entry.offset + entry.length) == DuckDBError) {
        return 0;
    }

    duckdb_vector child = duckdb_list_vector_get_child(vec);
    duckdb_vector key_vec = duckdb_struct_vector_get_child(child, 0);
    duckdb_vector val_vec = duckdb_struct_vector_get_child(child, 1);
    char **seen_keys = (char **)calloc((size_t)pair_count, sizeof(char *));
    int *seen_lens = (int *)calloc((size_t)pair_count, sizeof(int));
    if (!seen_keys || !seen_lens) {
        free(seen_keys);
        free(seen_lens);
        return 0;
    }

    const char *p = s;
    int write_idx = 0;
    int seen_count = 0;
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
        if (key_len > 0 && attr_group_index(seen_keys, seen_lens, seen_count, key, key_len) < 0) {
            int stored_key_len = 0;
            char *stored_key = dup_trimmed_span(key, key_len, &stored_key_len);
            if (!stored_key) {
                for (int i = 0; i < seen_count; i++) free(seen_keys[i]);
                free(seen_keys);
                free(seen_lens);
                return 0;
            }
            seen_keys[seen_count] = stored_key;
            seen_lens[seen_count] = stored_key_len;
            seen_count++;
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
    for (int i = 0; i < seen_count; i++) free(seen_keys[i]);
    free(seen_keys);
    free(seen_lens);
    return 1;
}

static void set_null_list_like(duckdb_vector vec, idx_t row) {
    duckdb_vector_ensure_validity_writable(vec);
    uint64_t *validity = duckdb_vector_get_validity(vec);
    duckdb_validity_set_row_invalid(validity, row);
    duckdb_list_entry entry = {duckdb_list_vector_get_size(vec), 0};
    duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(vec);
    list_data[row] = entry;
}

static int attr_group_index(char **keys, int *key_lens, int n_keys, const char *key, int key_len) {
    for (int i = 0; i < n_keys; i++) {
        if (key_lens[i] == key_len && memcmp(keys[i], key, (size_t)key_len) == 0) return i;
    }
    return -1;
}

static int fill_attr_list_map(duckdb_vector vec, idx_t row, const char *s, bool is_gff) {
    gxf_attr_pair_t *pairs = NULL;
    int n_pairs = 0;
    if (!parse_attr_pairs(s, is_gff, &pairs, &n_pairs)) return 0;
    if (n_pairs == 0) {
        set_null_list_like(vec, row);
        free_attr_pairs(pairs, n_pairs);
        return 1;
    }

    char **keys = (char **)calloc((size_t)n_pairs, sizeof(char *));
    int *key_lens = (int *)calloc((size_t)n_pairs, sizeof(int));
    int *counts = (int *)calloc((size_t)n_pairs, sizeof(int));
    if (!keys || !key_lens || !counts) {
        free(keys); free(key_lens); free(counts);
        free_attr_pairs(pairs, n_pairs);
        return 0;
    }

    int n_keys = 0;
    for (int i = 0; i < n_pairs; i++) {
        int g = attr_group_index(keys, key_lens, n_keys, pairs[i].key, pairs[i].key_len);
        if (g < 0) {
            keys[n_keys] = pairs[i].key;
            key_lens[n_keys] = pairs[i].key_len;
            counts[n_keys] = 1;
            n_keys++;
        } else {
            counts[g]++;
        }
    }

    duckdb_list_entry entry;
    entry.offset = duckdb_list_vector_get_size(vec);
    entry.length = (idx_t)n_keys;
    if (duckdb_list_vector_reserve(vec, entry.offset + entry.length) == DuckDBError ||
        duckdb_list_vector_set_size(vec, entry.offset + entry.length) == DuckDBError) {
        free(keys); free(key_lens); free(counts);
        free_attr_pairs(pairs, n_pairs);
        return 0;
    }

    duckdb_vector child = duckdb_list_vector_get_child(vec);
    duckdb_vector key_vec = duckdb_struct_vector_get_child(child, 0);
    duckdb_vector val_list_vec = duckdb_struct_vector_get_child(child, 1);

    for (int g = 0; g < n_keys; g++) {
        idx_t child_row = entry.offset + (idx_t)g;
        duckdb_vector_assign_string_element_len(key_vec, child_row, keys[g], key_lens[g]);

        duckdb_list_entry value_entry;
        value_entry.offset = duckdb_list_vector_get_size(val_list_vec);
        value_entry.length = (idx_t)counts[g];
        if (duckdb_list_vector_reserve(val_list_vec, value_entry.offset + value_entry.length) == DuckDBError ||
            duckdb_list_vector_set_size(val_list_vec, value_entry.offset + value_entry.length) == DuckDBError) {
            free(keys); free(key_lens); free(counts);
            free_attr_pairs(pairs, n_pairs);
            return 0;
        }
        duckdb_list_entry *value_entries = (duckdb_list_entry *)duckdb_vector_get_data(val_list_vec);
        duckdb_vector value_child = duckdb_list_vector_get_child(val_list_vec);
        idx_t out_idx = 0;
        for (int i = 0; i < n_pairs; i++) {
            if (pairs[i].key_len == key_lens[g] && memcmp(pairs[i].key, keys[g], (size_t)key_lens[g]) == 0) {
                duckdb_vector_assign_string_element_len(value_child, value_entry.offset + out_idx,
                                                        pairs[i].value, pairs[i].value_len);
                out_idx++;
            }
        }
        value_entries[child_row] = value_entry;
    }

    duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(vec);
    list_data[row] = entry;
    free(keys); free(key_lens); free(counts);
    free_attr_pairs(pairs, n_pairs);
    return 1;
}

static int fill_attr_pairs_list(duckdb_vector vec, idx_t row, const char *s, bool is_gff) {
    gxf_attr_pair_t *pairs = NULL;
    int n_pairs = 0;
    if (!parse_attr_pairs(s, is_gff, &pairs, &n_pairs)) return 0;
    if (n_pairs == 0) {
        set_null_list_like(vec, row);
        free_attr_pairs(pairs, n_pairs);
        return 1;
    }

    duckdb_list_entry entry;
    entry.offset = duckdb_list_vector_get_size(vec);
    entry.length = (idx_t)n_pairs;
    if (duckdb_list_vector_reserve(vec, entry.offset + entry.length) == DuckDBError ||
        duckdb_list_vector_set_size(vec, entry.offset + entry.length) == DuckDBError) {
        free_attr_pairs(pairs, n_pairs);
        return 0;
    }

    duckdb_vector child = duckdb_list_vector_get_child(vec);
    duckdb_vector key_vec = duckdb_struct_vector_get_child(child, 0);
    duckdb_vector val_vec = duckdb_struct_vector_get_child(child, 1);
    duckdb_vector idx_vec = duckdb_struct_vector_get_child(child, 2);
    int32_t *idx_data = (int32_t *)duckdb_vector_get_data(idx_vec);

    for (int i = 0; i < n_pairs; i++) {
        idx_t child_row = entry.offset + (idx_t)i;
        duckdb_vector_assign_string_element_len(key_vec, child_row, pairs[i].key, pairs[i].key_len);
        duckdb_vector_assign_string_element_len(val_vec, child_row, pairs[i].value, pairs[i].value_len);
        idx_data[child_row] = (int32_t)pairs[i].idx;
    }

    duckdb_list_entry *list_data = (duckdb_list_entry *)duckdb_vector_get_data(vec);
    list_data[row] = entry;
    free_attr_pairs(pairs, n_pairs);
    return 1;
}

static duckdb_logical_type create_attr_pairs_type(void) {
    duckdb_logical_type member_types[3];
    const char *member_names[3] = {"key", "value", "idx"};
    member_types[0] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    member_types[1] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    member_types[2] = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type struct_type = duckdb_create_struct_type(member_types, member_names, 3);
    duckdb_logical_type list_type = duckdb_create_list_type(struct_type);
    for (int i = 0; i < 3; i++) duckdb_destroy_logical_type(&member_types[i]);
    duckdb_destroy_logical_type(&struct_type);
    return list_type;
}

static duckdb_logical_type create_attr_list_map_type(void) {
    duckdb_logical_type key_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type val_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type val_list_type = duckdb_create_list_type(val_type);
    duckdb_logical_type map_type = duckdb_create_map_type(key_type, val_list_type);
    duckdb_destroy_logical_type(&key_type);
    duckdb_destroy_logical_type(&val_type);
    duckdb_destroy_logical_type(&val_list_type);
    return map_type;
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
    bd->meta_char = '#';
    bd->line_skip = 0;
    bd->col_types_provided = 0;
    bd->attr_map_col = -1;
    bd->attr_list_col = -1;
    bd->attr_pairs_col = -1;

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
            if (r && r[0]) {
                bd->region = strdup(r);
                parse_regions(bd->region, &bd->regions, &bd->n_regions);
            }
            duckdb_free((void *)r);
        }
        duckdb_destroy_value(&val);
    }

    /* named param: explicit index path */
    val = duckdb_bind_get_named_parameter(info, "index_path");
    if (val) {
        if (duckdb_get_type_id(duckdb_get_value_type(val)) == DUCKDB_TYPE_VARCHAR) {
            const char *idx = duckdb_get_varchar(val);
            if (idx && idx[0]) bd->index_path = strdup(idx);
            duckdb_free((void *)idx);
        }
        duckdb_destroy_value(&val);
    }

    /* Schema: GTF/GFF have fixed 9 columns; generic auto-detects */
    if (mode == TABIX_MODE_GTF || mode == TABIX_MODE_GFF) {
        /* File has 9 fixed columns; attributes_* columns are derived. */
        bd->n_cols = GXF_BASE_COL_COUNT;
        duckdb_value attr_map_val = duckdb_bind_get_named_parameter(info, "attributes_map");
        if (attr_map_val && !duckdb_is_null_value(attr_map_val)) {
            bd->include_attr_map = duckdb_get_bool(attr_map_val) ? 1 : 0;
        }
        if (attr_map_val) duckdb_destroy_value(&attr_map_val);

        duckdb_value attr_list_val = duckdb_bind_get_named_parameter(info, "attributes_list");
        if (attr_list_val && !duckdb_is_null_value(attr_list_val)) {
            bd->include_attr_list = duckdb_get_bool(attr_list_val) ? 1 : 0;
        }
        if (attr_list_val) duckdb_destroy_value(&attr_list_val);

        duckdb_value attr_pairs_val = duckdb_bind_get_named_parameter(info, "attributes_pairs");
        if (attr_pairs_val && !duckdb_is_null_value(attr_pairs_val)) {
            bd->include_attr_pairs = duckdb_get_bool(attr_pairs_val) ? 1 : 0;
        }
        if (attr_pairs_val) duckdb_destroy_value(&attr_pairs_val);

        if (mode == TABIX_MODE_GFF) {
            duckdb_value strict_val = duckdb_bind_get_named_parameter(info, "strict");
            if (strict_val && !duckdb_is_null_value(strict_val)) {
                bd->strict = duckdb_get_bool(strict_val) ? 1 : 0;
            }
            if (strict_val) duckdb_destroy_value(&strict_val);
        }

        duckdb_logical_type gxf_varchar = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
        duckdb_logical_type gxf_bigint = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
        duckdb_logical_type gxf_double = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
        duckdb_bind_add_result_column(info, "seqname", gxf_varchar);
        duckdb_bind_add_result_column(info, "source", gxf_varchar);
        duckdb_bind_add_result_column(info, "feature", gxf_varchar);
        duckdb_bind_add_result_column(info, "start", gxf_bigint);
        duckdb_bind_add_result_column(info, "end", gxf_bigint);
        duckdb_bind_add_result_column(info, "score", gxf_double);
        duckdb_bind_add_result_column(info, "strand", gxf_varchar);
        duckdb_bind_add_result_column(info, "frame", gxf_varchar);
        duckdb_bind_add_result_column(info, "attributes", gxf_varchar);
        duckdb_destroy_logical_type(&gxf_varchar);
        duckdb_destroy_logical_type(&gxf_bigint);
        duckdb_destroy_logical_type(&gxf_double);
        int next_attr_col = GXF_BASE_COL_COUNT;
        if (bd->include_attr_map) {
            duckdb_logical_type key_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
            duckdb_logical_type val_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
            duckdb_logical_type map_type = duckdb_create_map_type(key_type, val_type);
            bd->attr_map_col = next_attr_col++;
            duckdb_bind_add_result_column(info, "attributes_map", map_type);
            duckdb_destroy_logical_type(&key_type);
            duckdb_destroy_logical_type(&val_type);
            duckdb_destroy_logical_type(&map_type);
        }
        if (bd->include_attr_list) {
            duckdb_logical_type attr_list_type = create_attr_list_map_type();
            bd->attr_list_col = next_attr_col++;
            duckdb_bind_add_result_column(info, "attributes_list", attr_list_type);
            duckdb_destroy_logical_type(&attr_list_type);
        }
        if (bd->include_attr_pairs) {
            duckdb_logical_type attr_pairs_type = create_attr_pairs_type();
            bd->attr_pairs_col = next_attr_col++;
            duckdb_bind_add_result_column(info, "attributes_pairs", attr_pairs_type);
            duckdb_destroy_logical_type(&attr_pairs_type);
        }
    } else {
        /* Optional header handling for generic tabix */
        val = duckdb_bind_get_named_parameter(info, "header");
        if (val && !duckdb_is_null_value(val)) {
            bd->header = duckdb_get_bool(val) ? 1 : 0;
        }
        if (val) duckdb_destroy_value(&val);

        val = duckdb_bind_get_named_parameter(info, "auto_detect");
        if (val && !duckdb_is_null_value(val)) {
            bd->auto_detect = duckdb_get_bool(val) ? 1 : 0;
        }
        if (val) duckdb_destroy_value(&val);

        val = duckdb_bind_get_named_parameter(info, "header_names");
        if (val && !duckdb_is_null_value(val)) {
            idx_t n = duckdb_get_list_size(val);
            if (n > 0) {
                bd->header_names = (char **)malloc(sizeof(char *) * (size_t)n);
                bd->header_names_count = (int)n;
                for (idx_t i = 0; i < n; i++) {
                    duckdb_value elem = duckdb_get_list_child(val, i);
                    bd->header_names[i] = duckdb_get_varchar(elem);
                    duckdb_destroy_value(&elem);
                }
            }
        }
        if (val) duckdb_destroy_value(&val);

        val = duckdb_bind_get_named_parameter(info, "column_types");
        if (val && !duckdb_is_null_value(val)) {
            idx_t n = duckdb_get_list_size(val);
            if (n > 0) {
                bd->col_types = (int *)malloc(sizeof(int) * (size_t)n);
                bd->col_types_provided = 1;
                for (idx_t i = 0; i < n; i++) {
                    duckdb_value elem = duckdb_get_list_child(val, i);
                    const char *tname = duckdb_get_varchar(elem);
                    bd->col_types[i] = parse_type_name(tname);
                    duckdb_free((void *)tname);
                    duckdb_destroy_value(&elem);
                }
                bd->n_cols = (int)n;
            }
        }
        if (val) duckdb_destroy_value(&val);

        /* Generic tabix: if indexed, use header/meta settings; otherwise peek */
        htsFile *fp = hts_open(bd->file_path, "r");
        if (!fp) {
            duckdb_bind_set_error(info, "Cannot open file");
            tabix_bind_data_destroy(bd);
            return;
        }
        kstring_t line = {0, 0, NULL};
        int n_cols = 0;
        char meta_char = '#';
        int skip_lines = 0;
        int header_from_skip = 0;
        char *header_candidate = NULL;

        tbx_t *tbx = tbx_index_load2(bd->file_path, bd->index_path);
        if (tbx) {
            tbx_conf_t conf = tbx->conf;
            meta_char = conf.meta_char ? conf.meta_char : '#';
            skip_lines = conf.line_skip;
        }
        bd->meta_char = meta_char;
        bd->line_skip = skip_lines;

        /* Skip comment/header lines (meta char or line_skip prefix) */
        while (hts_getline(fp, '\n', &line) >= 0) {
            if (line.l == 0) continue;
            if (skip_lines > 0) {
                if (bd->header && !bd->header_names && line.l > 0) {
                    free(header_candidate);
                    header_candidate = strdup(line.s);
                    header_from_skip = 1;
                }
                skip_lines--;
                continue;
            }
            if (meta_char && line.s[0] == meta_char) {
                continue;
            }
            if (bd->header && !bd->header_names && !header_candidate) {
                header_candidate = strdup(line.s);
            } else {
                n_cols = count_fields(line.s);
                break;
            }
        }
        free(line.s);
        hts_close(fp);
        if (tbx) tbx_destroy(tbx);

        if (bd->header_names && bd->header_names_count > 0) {
            n_cols = bd->header_names_count;
            bd->skip_header_line = bd->header ? 1 : 0;
        } else if (bd->header && header_candidate) {
            bd->header_names_count = parse_header_names(header_candidate, &bd->header_names);
            n_cols = bd->header_names_count;
            bd->skip_header_line = header_from_skip ? 0 : 1;
        }
        free(header_candidate);

        if (n_cols == 0) n_cols = 1;
        if (n_cols > TABIX_MAX_GENERIC_COLS) n_cols = TABIX_MAX_GENERIC_COLS;
        if (bd->n_cols == 0) bd->n_cols = n_cols;

        if (bd->col_types && bd->n_cols != n_cols) {
            duckdb_bind_set_error(info, "column_types length does not match detected column count");
            tabix_bind_data_destroy(bd);
            return;
        }

        if (!bd->col_types) {
            bd->col_types = (int *)malloc(sizeof(int) * (size_t)bd->n_cols);
            for (int i = 0; i < bd->n_cols; i++) bd->col_types[i] = DUCKDB_TYPE_VARCHAR;
        }

        if (bd->auto_detect && !bd->col_types_provided) {
            int *type_state = (int *)malloc(sizeof(int) * (size_t)bd->n_cols);
            for (int i = 0; i < bd->n_cols; i++) type_state[i] = DUCKDB_TYPE_INTEGER;

            /* Reopen to scan for type inference */
            htsFile *fp2 = hts_open(bd->file_path, "r");
            if (fp2) {
                kstring_t l2 = {0, 0, NULL};
                int skip2 = bd->line_skip;
                int skip_header = bd->skip_header_line;
                int seen = 0;
                while (seen < 100 && hts_getline(fp2, '\n', &l2) >= 0) {
                    if (l2.l == 0) continue;
                    if (skip2 > 0) { skip2--; continue; }
                    if (bd->meta_char && l2.s[0] == bd->meta_char) continue;
                    if (skip_header) { skip_header = 0; continue; }

                    for (int i = 0; i < bd->n_cols; i++) {
                        int flen = 0;
                        const char *fld = get_field(l2.s, i, &flen);
                        if (!fld || flen == 0 || (flen == 1 && fld[0] == '.')) continue;
                        if (is_integer_field(fld, flen)) {
                            continue;
                        } else if (is_float_field(fld, flen)) {
                            if (type_state[i] != DUCKDB_TYPE_VARCHAR)
                                type_state[i] = DUCKDB_TYPE_DOUBLE;
                        } else {
                            type_state[i] = DUCKDB_TYPE_VARCHAR;
                        }
                    }
                    seen++;
                }
                free(l2.s);
                hts_close(fp2);
            }

            for (int i = 0; i < bd->n_cols; i++) {
                if (type_state[i] == DUCKDB_TYPE_INTEGER) {
                    bd->col_types[i] = DUCKDB_TYPE_BIGINT;
                } else if (type_state[i] == DUCKDB_TYPE_DOUBLE) {
                    bd->col_types[i] = DUCKDB_TYPE_DOUBLE;
                } else {
                    bd->col_types[i] = DUCKDB_TYPE_VARCHAR;
                }
            }
            free(type_state);
        }

        duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
        for (int i = 0; i < bd->n_cols; i++) {
            char col_name[32];
            duckdb_logical_type col_type = duckdb_create_logical_type((duckdb_type)bd->col_types[i]);
            const char *name = NULL;
            if (bd->header_names && i < bd->header_names_count &&
                bd->header_names[i] && bd->header_names[i][0] != '\0') {
                name = bd->header_names[i];
            } else {
                snprintf(col_name, sizeof(col_name), "column%d", i);
                name = col_name;
            }
            duckdb_bind_add_result_column(info, name, col_type);
            duckdb_destroy_logical_type(&col_type);
        }
        duckdb_destroy_logical_type(&varchar_type);
    }

    if (bd->n_regions == 0) {
        tbx_t *tbx_stats = tbx_index_load2(bd->file_path, bd->index_path);
        if (tbx_stats) {
            bd->index_row_count_valid =
                tabix_try_get_index_row_count(tbx_stats, &bd->index_row_count);
            if (bd->index_row_count_valid && bd->skip_header_line && bd->index_row_count > 0) {
                bd->index_row_count--;
            }
            tbx_destroy(tbx_stats);
        }
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

    id->n_projected_cols = duckdb_init_get_column_count(info);
    if (id->n_projected_cols > 0) {
        id->column_ids = (idx_t *)malloc(sizeof(idx_t) * id->n_projected_cols);
        for (idx_t i = 0; i < id->n_projected_cols; i++) {
            id->column_ids[i] = duckdb_init_get_column_index(info, i);
        }
    } else {
        id->column_ids = NULL;
    }

    if (id->n_projected_cols == 0 && bd->n_regions == 0 && bd->index_row_count_valid && !bd->strict) {
        id->count_only = 1;
        id->count_remaining = bd->index_row_count;
        id->finished = (id->count_remaining == 0);
        duckdb_init_set_init_data(info, id, tabix_init_data_destroy);
        return;
    }

    id->fp = hts_open(bd->file_path, "r");
    if (!id->fp) {
        char msg[512];
        snprintf(msg, sizeof(msg), "Cannot open file: %s", bd->file_path);
        duckdb_init_set_error(info, msg);
        tabix_init_data_destroy(id);
        return;
    }

    /* Try to load tabix index */
    id->tbx = tbx_index_load2(bd->file_path, bd->index_path);

    if (bd->n_regions > 0) {
        if (!id->tbx) {
            char msg[512];
            snprintf(msg, sizeof(msg),
                     "Region query requested but no tabix index found for: %s",
                     bd->file_path);
            duckdb_init_set_error(info, msg);
            tabix_init_data_destroy(id);
            return;
        }
        if (bd->n_regions > 1) {
            vcf_emit_warning("Multi-region tabix queries are executed as a chained union of single-region iterators; overlapping regions may return duplicate rows");
        }
        id->next_region_idx = 0;
        if (!tabix_advance_region_iterator(id, bd)) {
            /* Region doesn't match any sequences – return empty result */
            id->finished = true;
        }
    }
    /* else: sequential scan (no iterator) */

    id->line.l = 0;
    id->line.m = 0;
    id->line.s = NULL;
    id->skip_remaining = bd->line_skip;
    id->skipped_header = 0;

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

    if (id->count_only) {
        idx_t row_count = (id->count_remaining > (uint64_t)duckdb_vector_size())
            ? duckdb_vector_size()
            : (idx_t)id->count_remaining;
        id->count_remaining -= row_count;
        if (id->count_remaining == 0) {
            id->finished = true;
        }
        duckdb_data_chunk_set_size(output, row_count);
        return;
    }

    idx_t row_count = 0;
    int n_cols = bd->n_cols;
    idx_t chunk_col_count = duckdb_data_chunk_get_column_count(output);

    /* Pre-fetch column vectors */
    duckdb_vector *vectors = NULL;
    if (chunk_col_count > 0) {
        vectors = (duckdb_vector *)malloc(sizeof(duckdb_vector) * chunk_col_count);
        if (!vectors) {
            duckdb_function_set_error(info, "tabix_scan: out of memory allocating vectors");
            id->finished = true;
            duckdb_data_chunk_set_size(output, 0);
            return;
        }
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
            if (id->itr && tabix_advance_region_iterator(id, bd)) {
                continue;
            }
            id->finished = true;
            break;
        }
        id->line_number++;

        /* Skip comment/header lines */
        if (id->line.l == 0) continue;
        if (!id->itr && id->skip_remaining > 0) {
            id->skip_remaining--;
            continue;
        }
        if (bd->mode == TABIX_MODE_GFF && strncmp(id->line.s, "##FASTA", 7) == 0) {
            id->finished = true;
            break;
        }
        if (bd->meta_char && id->line.s[0] == bd->meta_char) continue;
        if (!id->itr && bd->skip_header_line && !id->skipped_header) {
            id->skipped_header = 1;
            continue;
        }

        if (bd->strict && bd->mode == TABIX_MODE_GFF) {
            char err[256];
            if (!validate_gff3_line_strict(id->line.s, err, sizeof(err))) {
                char msg[384];
                if (id->itr) {
                    snprintf(msg, sizeof(msg), "read_gff strict validation failed during indexed region scan: %s", err);
                } else if (id->line_number > 0) {
                    snprintf(msg, sizeof(msg), "read_gff strict validation failed at line %llu: %s",
                             (unsigned long long)id->line_number, err);
                } else {
                    snprintf(msg, sizeof(msg), "read_gff strict validation failed: %s", err);
                }
                duckdb_function_set_error(info, msg);
                id->finished = true;
                if (vectors) free(vectors);
                duckdb_data_chunk_set_size(output, 0);
                return;
            }
        }

        if (chunk_col_count > 0) {
        if (bd->mode == TABIX_MODE_GTF || bd->mode == TABIX_MODE_GFF) {
            /* Parse GTF/GFF columns with projection pushdown:
             * Vector index c maps to logical column id->column_ids[c],
             * which is the field index in the TSV line. */
            for (idx_t c = 0; c < chunk_col_count; c++) {
                int logical_col = (int)id->column_ids[c];
                if (logical_col == bd->attr_map_col && bd->include_attr_map) {
                    const char *fld = NULL;
                    int flen = 0;
                    fld = get_field(id->line.s, GXF_COL_ATTRIBUTES, &flen);
                        char *tmp = NULL;
                        if (fld && flen > 0) {
                            tmp = (char *)malloc((size_t)flen + 1);
                            if (!tmp) {
                                duckdb_function_set_error(info, "tabix_scan: out of memory parsing attributes");
                                id->finished = true;
                                if (vectors) free(vectors);
                                duckdb_data_chunk_set_size(output, 0);
                                return;
                            }
                            memcpy(tmp, fld, (size_t)flen);
                            tmp[flen] = '\0';
                        }
                        if (!fill_attr_map(vectors[c], row_count, tmp ? tmp : ".", bd->mode == TABIX_MODE_GFF)) {
                            if (tmp) free(tmp);
                            duckdb_function_set_error(info, "tabix_scan: out of memory parsing attributes_map");
                            id->finished = true;
                            if (vectors) free(vectors);
                            duckdb_data_chunk_set_size(output, 0);
                            return;
                        }
                        if (tmp) free(tmp);
                        continue;
                    }
                if (logical_col == bd->attr_list_col && bd->include_attr_list) {
                    int flen = 0;
                    const char *fld = get_field(id->line.s, GXF_COL_ATTRIBUTES, &flen);
                    char *tmp = NULL;
                    if (fld && flen > 0) {
                        tmp = (char *)malloc((size_t)flen + 1);
                        if (!tmp) {
                            duckdb_function_set_error(info, "tabix_scan: out of memory parsing attributes_list");
                            id->finished = true;
                            if (vectors) free(vectors);
                            duckdb_data_chunk_set_size(output, 0);
                            return;
                        }
                        memcpy(tmp, fld, (size_t)flen);
                        tmp[flen] = '\0';
                    }
                    if (!fill_attr_list_map(vectors[c], row_count, tmp ? tmp : ".", bd->mode == TABIX_MODE_GFF)) {
                        if (tmp) free(tmp);
                        duckdb_function_set_error(info, "tabix_scan: out of memory parsing attributes_list");
                        id->finished = true;
                        if (vectors) free(vectors);
                        duckdb_data_chunk_set_size(output, 0);
                        return;
                    }
                    if (tmp) free(tmp);
                    continue;
                }
                if (logical_col == bd->attr_pairs_col && bd->include_attr_pairs) {
                    int flen = 0;
                    const char *fld = get_field(id->line.s, GXF_COL_ATTRIBUTES, &flen);
                    char *tmp = NULL;
                    if (fld && flen > 0) {
                        tmp = (char *)malloc((size_t)flen + 1);
                        if (!tmp) {
                            duckdb_function_set_error(info, "tabix_scan: out of memory parsing attributes_pairs");
                            id->finished = true;
                            if (vectors) free(vectors);
                            duckdb_data_chunk_set_size(output, 0);
                            return;
                        }
                        memcpy(tmp, fld, (size_t)flen);
                        tmp[flen] = '\0';
                    }
                    if (!fill_attr_pairs_list(vectors[c], row_count, tmp ? tmp : ".", bd->mode == TABIX_MODE_GFF)) {
                        if (tmp) free(tmp);
                        duckdb_function_set_error(info, "tabix_scan: out of memory parsing attributes_pairs");
                        id->finished = true;
                        if (vectors) free(vectors);
                        duckdb_data_chunk_set_size(output, 0);
                        return;
                    }
                    if (tmp) free(tmp);
                    continue;
                }
                if (logical_col >= GXF_BASE_COL_COUNT) {
                    duckdb_function_set_error(info, "tabix_scan: internal derived-column projection mismatch");
                    id->finished = true;
                    if (vectors) free(vectors);
                    duckdb_data_chunk_set_size(output, 0);
                    return;
                }

                int flen = 0;
                const char *fld = get_field(id->line.s, logical_col, &flen);

                if (!fld || flen == 0 || (flen == 1 && fld[0] == '.')) {
                    /* Missing value */
                    switch (logical_col) {
                        case GXF_COL_START:
                        case GXF_COL_END: {
                            if (bd->strict && bd->mode == TABIX_MODE_GFF) {
                                set_null(vectors[c], row_count);
                            } else {
                                int64_t *data = (int64_t *)duckdb_vector_get_data(vectors[c]);
                                data[row_count] = 0;
                            }
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
                        int64_t ival = 0;
                        if (!parse_int64_span(fld, flen, &ival)) {
                            set_null(vectors[c], row_count);
                            break;
                        }
                        int64_t *data = (int64_t *)duckdb_vector_get_data(vectors[c]);
                        data[row_count] = ival;
                        break;
                    }
                    case GXF_COL_SCORE: {
                        double dval = 0.0;
                        if (!parse_double_span(fld, flen, &dval)) {
                            set_null(vectors[c], row_count);
                            break;
                        }
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
            /* Generic mode with optional typing */
            for (idx_t c = 0; c < chunk_col_count; c++) {
                int logical_col = (int)id->column_ids[c];
                if (logical_col >= n_cols) continue;
                int flen = 0;
                const char *fld = get_field(id->line.s, logical_col, &flen);
                if (!fld || flen == 0 || (flen == 1 && fld[0] == '.')) {
                    set_null(vectors[c], row_count);
                    continue;
                }
                int t = bd->col_types ? bd->col_types[logical_col] : DUCKDB_TYPE_VARCHAR;
                if (t == DUCKDB_TYPE_INTEGER || t == DUCKDB_TYPE_BIGINT) {
                    int64_t v = 0;
                    if (parse_int64_span(fld, flen, &v)) {
                        if (t == DUCKDB_TYPE_INTEGER) {
                            int32_t *data = (int32_t *)duckdb_vector_get_data(vectors[c]);
                            data[row_count] = (int32_t)v;
                        } else {
                            int64_t *data = (int64_t *)duckdb_vector_get_data(vectors[c]);
                            data[row_count] = v;
                        }
                    } else {
                        set_null(vectors[c], row_count);
                    }
                } else if (t == DUCKDB_TYPE_DOUBLE) {
                    double v = 0.0;
                    if (parse_double_span(fld, flen, &v)) {
                        double *data = (double *)duckdb_vector_get_data(vectors[c]);
                        data[row_count] = v;
                    } else {
                        set_null(vectors[c], row_count);
                    }
                } else {
                    duckdb_vector_assign_string_element_len(vectors[c], row_count, fld, flen);
                }
            }
        }
        } /* chunk_col_count > 0 */

        row_count++;
    }

    if (vectors) free(vectors);
    duckdb_data_chunk_set_size(output, row_count);
}

/* ================================================================
 * Registration helpers
 * ================================================================ */

static duckdb_table_function create_tabix_tf(const char *name,
                                              void (*bind_fn)(duckdb_bind_info),
                                              int include_attributes_map,
                                              int include_strict) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, name);

    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "region", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);
    duckdb_destroy_logical_type(&varchar_type);

    if (include_attributes_map) {
        duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
        duckdb_table_function_add_named_parameter(tf, "attributes_map", bool_type);
        duckdb_table_function_add_named_parameter(tf, "attributes_list", bool_type);
        duckdb_table_function_add_named_parameter(tf, "attributes_pairs", bool_type);
        duckdb_destroy_logical_type(&bool_type);
    }

    if (include_strict) {
        duckdb_logical_type strict_bool = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
        duckdb_table_function_add_named_parameter(tf, "strict", strict_bool);
        duckdb_destroy_logical_type(&strict_bool);
    }

    duckdb_logical_type header_bool = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_table_function_add_named_parameter(tf, "header", header_bool);
    duckdb_destroy_logical_type(&header_bool);

    duckdb_logical_type list_child = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type list_type = duckdb_create_list_type(list_child);
    duckdb_table_function_add_named_parameter(tf, "header_names", list_type);
    duckdb_destroy_logical_type(&list_child);
    duckdb_destroy_logical_type(&list_type);

    duckdb_logical_type auto_bool = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_table_function_add_named_parameter(tf, "auto_detect", auto_bool);
    duckdb_destroy_logical_type(&auto_bool);

    duckdb_logical_type types_child = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type types_list = duckdb_create_list_type(types_child);
    duckdb_table_function_add_named_parameter(tf, "column_types", types_list);
    duckdb_destroy_logical_type(&types_child);
    duckdb_destroy_logical_type(&types_list);

    duckdb_table_function_set_bind(tf, bind_fn);
    duckdb_table_function_set_init(tf, tabix_init);
    duckdb_table_function_set_function(tf, tabix_scan);
    duckdb_table_function_supports_projection_pushdown(tf, true);

    return tf;
}

void register_read_tabix_function(duckdb_connection connection) {
    duckdb_table_function tf = create_tabix_tf("read_tabix", read_tabix_bind, 0, 0);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}

void register_read_gtf_function(duckdb_connection connection) {
    duckdb_table_function tf = create_tabix_tf("read_gtf", read_gtf_bind, 1, 0);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}

void register_read_gff_function(duckdb_connection connection) {
    duckdb_table_function tf = create_tabix_tf("read_gff", read_gff_bind, 1, 1);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
}
