/**
 * DuckHTS GenBank reader: DuckDB glue over genbank_core.
 *
 * read_genbank(path, attributes_map := FALSE, attributes := []) emits one row per feature
 * segment in read_gff's column shape. Bind declares the schema, init opens
 * the hFILE, and each scan pulls the next record from the core as the chunk
 * drains, so memory is bounded by the largest record. Only projected columns
 * are materialized; the GFF3 attribute string is built once per feature and
 * only when attributes or attributes_map is requested. attributes := [...]
 * adds one VARCHAR column per qualifier key, equal to attributes_map[key];
 * a value is computed only for a projected key.
 *
 * genbank_to_fasta(path, output_path := NULL, line_width := 70, overwrite := FALSE)
 * writes each record's ORIGIN under the name read_genbank reports as seqname.
 * Bind validates parameters only. Init opens the input, claims the destination
 * (O_CREAT|O_EXCL unless overwrite) and opens a temp file beside it; the first
 * scan streams residues to the temp file and renames it into place after a
 * clean end of input. The init destructor owns every cleanup.
 *
 * The GenBank to GFF3 mapping is documented in functions.yaml and genbank_core.h.
 */

/* mkstemp, fdopen and fchmod are POSIX, not ISO C; declare them under a strict -std=c11. */
#if !defined(_WIN32) && !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE 1
#endif

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include "duckdb_alloc.h"
#include "duckdb_list.h"
#include "genbank_core.h"
#include "include/named_attribute_columns.h"

#include <errno.h>
#include <fcntl.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef _WIN32
#include <windows.h>
#endif

#include <htslib/hts.h>
#include <htslib/kseq.h>

enum {
    GB_COL_SEQNAME = 0,
    GB_COL_SOURCE,
    GB_COL_FEATURE,
    GB_COL_START,
    GB_COL_END,
    GB_COL_SCORE,
    GB_COL_STRAND,
    GB_COL_FRAME,
    GB_COL_ATTRIBUTES,
    GB_COL_ATTRIBUTES_MAP
};

/* ------------------------------------------------------------------------ */
/* Shared helpers                                                            */
/* ------------------------------------------------------------------------ */

/* "<function>: <path>: line N: <message>", omitting the line when unknown. */
static void gb_format_error(char *buf, size_t n, const char *function, const char *path, const gb_error_t *e) {
    if (e->line > 0)
        snprintf(buf, n, "%s: %s: line %ld: %s", function, path, e->line, e->msg);
    else
        snprintf(buf, n, "%s: %s: %s", function, path, e->msg);
}

/* Non-NULL, non-empty VARCHAR parameter; NULL on anything else. Caller frees. */
static char *gb_bind_string(duckdb_value val) {
    if (!val || duckdb_is_null_value(val)) return NULL;
    if (duckdb_get_type_id(duckdb_get_value_type(val)) != DUCKDB_TYPE_VARCHAR) return NULL;
    char *s = duckdb_get_varchar(val);
    if (s && s[0] == '\0') {
        duckdb_free(s);
        return NULL;
    }
    return s;
}

static gb_span_t gb_contig(const gb_parser_t *p) {
    return p->version.len ? p->version : p->accession.len ? p->accession : p->locus;
}

/* ------------------------------------------------------------------------ */
/* read_genbank                                                              */
/* ------------------------------------------------------------------------ */

typedef struct {
    char *path;
    bool include_attr_map;
    duckhts_attribute_columns named_attributes;
} gb_bind_t;

static void gb_bind_destroy(void *data) {
    gb_bind_t *bd = (gb_bind_t *)data;
    if (!bd) return;
    duckhts_attribute_columns_destroy(&bd->named_attributes);
    duckdb_free(bd->path);
    duckdb_free(bd);
}

static void read_genbank_bind(duckdb_bind_info info) {
    duckdb_value val = duckdb_bind_get_parameter(info, 0);
    char *path = gb_bind_string(val);
    if (val) duckdb_destroy_value(&val);
    if (!path) {
        duckdb_bind_set_error(info, "read_genbank requires a file path");
        return;
    }

    gb_bind_t *bd = (gb_bind_t *)duckdb_malloc(sizeof(*bd));
    if (!bd) {
        duckdb_free(path);
        duckdb_bind_set_error(info, "read_genbank: out of memory");
        return;
    }
    memset(bd, 0, sizeof(*bd));
    bd->path = path;

    val = duckdb_bind_get_named_parameter(info, "attributes_map");
    if (val && !duckdb_is_null_value(val)) bd->include_attr_map = duckdb_get_bool(val);
    if (val) duckdb_destroy_value(&val);

    const char *reserved[] = {"seqname", "source", "feature", "start", "end", "score",
                              "strand", "frame", "attributes", "attributes_map"};
    size_t reserved_count = GB_COL_ATTRIBUTES_MAP + (bd->include_attr_map ? 1 : 0);
    if (!duckhts_attribute_columns_bind(info, "attributes", reserved, reserved_count,
                                        &bd->named_attributes)) {
        gb_bind_destroy(bd);
        return;
    }

    duckdb_logical_type vc = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bi = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type db = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_bind_add_result_column(info, "seqname", vc);
    duckdb_bind_add_result_column(info, "source", vc);
    duckdb_bind_add_result_column(info, "feature", vc);
    duckdb_bind_add_result_column(info, "start", bi);
    duckdb_bind_add_result_column(info, "end", bi);
    duckdb_bind_add_result_column(info, "score", db);
    duckdb_bind_add_result_column(info, "strand", vc);
    duckdb_bind_add_result_column(info, "frame", vc);
    duckdb_bind_add_result_column(info, "attributes", vc);
    duckdb_destroy_logical_type(&vc);
    duckdb_destroy_logical_type(&bi);
    duckdb_destroy_logical_type(&db);
    if (bd->include_attr_map) {
        duckdb_logical_type kt = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
        duckdb_logical_type vt = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
        duckdb_logical_type mt = duckdb_create_map_type(kt, vt);
        duckdb_bind_add_result_column(info, "attributes_map", mt);
        duckdb_destroy_logical_type(&kt);
        duckdb_destroy_logical_type(&vt);
        duckdb_destroy_logical_type(&mt);
    }
    duckhts_attribute_columns_declare(info, &bd->named_attributes, reserved_count);
    duckdb_bind_set_bind_data(info, bd, gb_bind_destroy);
}

typedef enum {
    GB_SCAN_PULL = 0, /* no record in the parser: read until the next // */
    GB_SCAN_EMIT,     /* the parser holds a resolved record; rows remain */
    GB_SCAN_DONE      /* input consumed to a clean end */
} gb_scan_state_t;

typedef struct {
    htsFile *fp;
    kstring_t line;
    gb_parser_t parser;
    gb_scan_state_t state;
    size_t row_cursor; /* next parser row to emit */
    gb_span_t contig;
    kstring_t attr;     /* GFF3 attribute string cache */
    size_t attr_feature; /* feature the cache was built for; SIZE_MAX when none */
    kstring_t key, value;
    idx_t *column_ids;
    idx_t n_projected;
    duckhts_projected_attribute *projected_attributes;
    idx_t n_projected_attributes;
    kstring_t *attribute_values; /* one buffer per projected attribute */
    duckdb_vector *vectors;      /* output vectors of the current chunk */
} gb_init_t;

static void gb_init_destroy(void *data) {
    gb_init_t *st = (gb_init_t *)data;
    if (!st) return;
    if (st->fp) hts_close(st->fp);
    free(st->line.s);
    free(st->attr.s);
    free(st->key.s);
    free(st->value.s);
    gb_parser_destroy(&st->parser);
    free(st->column_ids);
    free(st->projected_attributes);
    for (idx_t i = 0; st->attribute_values && i < st->n_projected_attributes; i++) {
        free(st->attribute_values[i].s);
    }
    free(st->attribute_values);
    free(st->vectors);
    free(st);
}

static void read_genbank_init(duckdb_init_info info) {
    gb_bind_t *bd = (gb_bind_t *)duckdb_init_get_bind_data(info);
    gb_init_t *st = (gb_init_t *)calloc(1, sizeof(*st));
    if (!st) {
        duckdb_init_set_error(info, "read_genbank: out of memory");
        return;
    }
    gb_parser_init(&st->parser, GB_MODE_FEATURES);
    st->parser.flags = GB_FLAG_DROP_TRANSLATION;
    st->attr_feature = SIZE_MAX;

    st->fp = hts_open(bd->path, "r");
    if (!st->fp) {
        char msg[512];
        snprintf(msg, sizeof(msg), "read_genbank: cannot open file: %s", bd->path);
        duckdb_init_set_error(info, msg);
        gb_init_destroy(st);
        return;
    }

    st->n_projected = duckdb_init_get_column_count(info);
    if (st->n_projected) {
        st->column_ids = (idx_t *)malloc(sizeof(idx_t) * st->n_projected);
        if (!st->column_ids) {
            duckdb_init_set_error(info, "read_genbank: out of memory");
            gb_init_destroy(st);
            return;
        }
        for (idx_t i = 0; i < st->n_projected; i++) st->column_ids[i] = duckdb_init_get_column_index(info, i);
    }
    if (bd->named_attributes.count && st->n_projected) {
        st->projected_attributes = duckhts_attribute_columns_project(
            &bd->named_attributes, st->column_ids, st->n_projected, &st->n_projected_attributes);
        if (st->n_projected_attributes) {
            st->attribute_values = (kstring_t *)calloc(st->n_projected_attributes,
                                                       sizeof(*st->attribute_values));
            st->vectors = (duckdb_vector *)calloc(st->n_projected, sizeof(*st->vectors));
        }
        if (st->n_projected_attributes &&
            (!st->projected_attributes || !st->attribute_values || !st->vectors)) {
            duckdb_init_set_error(info, "read_genbank: out of memory");
            gb_init_destroy(st);
            return;
        }
    }
    duckdb_init_set_max_threads(info, 1);
    duckdb_init_set_init_data(info, st, gb_init_destroy);
}

/* Pull the next record into the parser. 1 record, 0 clean end of input, -1 error (msg set). */
static int gb_next_record(gb_init_t *st, const char *path, char *msg, size_t msg_len) {
    st->state = GB_SCAN_PULL;
    st->attr_feature = SIZE_MAX;
    for (;;) {
        int n = hts_getline(st->fp, KS_SEP_LINE, &st->line);
        if (n < -1) {
            snprintf(msg, msg_len, "read_genbank: %s: read failed after line %ld%s%s", path, st->parser.line_no,
                     errno ? ": " : "", errno ? strerror(errno) : "");
            return -1;
        }
        if (n == -1) {
            if (gb_parser_finish(&st->parser) != GB_OK) {
                gb_format_error(msg, msg_len, "read_genbank", path, &st->parser.err);
                return -1;
            }
            st->state = GB_SCAN_DONE;
            return 0;
        }
        gb_feed_t r = gb_parser_feed(&st->parser, st->line.s, st->line.l);
        if (r == GB_FEED_ERROR) {
            gb_format_error(msg, msg_len, "read_genbank", path, &st->parser.err);
            return -1;
        }
        if (r == GB_FEED_RECORD) {
            st->state = GB_SCAN_EMIT;
            st->row_cursor = 0;
            st->contig = gb_contig(&st->parser);
            return 1;
        }
    }
}

static bool gb_ensure_attributes(gb_init_t *st, size_t feature) {
    if (st->attr_feature == feature) return true;
    if (gb_feature_attributes(&st->parser, feature, &st->attr) != GB_OK) return false;
    st->attr_feature = feature;
    return true;
}

static bool gb_fill_attr_map(gb_init_t *st, duckdb_vector vec, idx_t row, size_t feature) {
    size_t n = gb_feature_attr_count(&st->parser, feature);
    duckdb_list_entry entry;
    if (!duckhts_list_extend(vec, (idx_t)n, &entry)) return false;
    duckdb_vector child = duckdb_list_vector_get_child(vec);
    duckdb_vector key_vec = duckdb_struct_vector_get_child(child, 0);
    duckdb_vector val_vec = duckdb_struct_vector_get_child(child, 1);
    for (size_t i = 0; i < n; i++) {
        if (gb_feature_attr_at(&st->parser, feature, i, &st->key, &st->value) != GB_OK) return false;
        duckdb_vector_assign_string_element_len(key_vec, entry.offset + i, st->key.s, st->key.l);
        duckdb_vector_assign_string_element_len(val_vec, entry.offset + i, st->value.s, st->value.l);
    }
    ((duckdb_list_entry *)duckdb_vector_get_data(vec))[row] = entry;
    return true;
}

/* Point each projected attribute at the value attributes_map[key] would hold, or
 * NULL when the key is absent. The first pair with a key wins, as in a MAP lookup.
 * Values are encoded only for matching keys. False only on allocation failure. */
static bool gb_find_projected_attributes(gb_init_t *st, size_t feature) {
    duckhts_projected_attribute *projected = st->projected_attributes;
    idx_t remaining = st->n_projected_attributes;
    for (idx_t i = 0; i < remaining; i++) projected[i].value = NULL;

    size_t n = gb_feature_attr_count(&st->parser, feature);
    for (size_t index = 0; index < n && remaining; index++) {
        if (gb_feature_attr_key(&st->parser, feature, index, &st->key) != GB_OK) return false;
        bool value_ready = false;
        for (idx_t i = 0; i < st->n_projected_attributes; i++) {
            if (projected[i].value || projected[i].key->length != st->key.l ||
                memcmp(projected[i].key->name, st->key.s, st->key.l) != 0) {
                continue;
            }
            if (!value_ready) {
                if (gb_feature_attr_value(&st->parser, feature, index, &st->value) != GB_OK) return false;
                value_ready = true;
            }
            kstring_t *owned = &st->attribute_values[i];
            owned->l = 0;
            if (kputsn(st->value.s ? st->value.s : "", st->value.l, owned) == EOF) return false;
            projected[i].value = owned->s;
            projected[i].value_length = owned->l;
            remaining--;
        }
    }
    return true;
}

/* Fill one projected column of one output row. False only on allocation failure. */
static bool gb_fill_column(gb_init_t *st, duckdb_vector vec, idx_t col, idx_t row, const gb_row_t *r) {
    const gb_parser_t *p = &st->parser;
    const gb_feature_t *f = &p->feats[r->feature];
    switch (col) {
    case GB_COL_SEQNAME: {
        gb_span_t name = r->remote.len ? r->remote : st->contig;
        duckdb_vector_assign_string_element_len(vec, row, gb_str(p, name), name.len);
        return true;
    }
    case GB_COL_SOURCE:
        duckdb_vector_assign_string_element(vec, row, "GenBank");
        return true;
    case GB_COL_FEATURE:
        duckdb_vector_assign_string_element_len(vec, row, gb_str(p, f->key), f->key.len);
        return true;
    case GB_COL_START:
        ((int64_t *)duckdb_vector_get_data(vec))[row] = r->start;
        return true;
    case GB_COL_END:
        ((int64_t *)duckdb_vector_get_data(vec))[row] = r->end;
        return true;
    case GB_COL_SCORE:
        duckdb_vector_ensure_validity_writable(vec);
        duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec), row);
        return true;
    case GB_COL_STRAND:
        duckdb_vector_assign_string_element_len(vec, row, &r->strand, 1);
        return true;
    case GB_COL_FRAME: {
        char frame = r->phase < 0 ? '.' : (char)('0' + r->phase);
        duckdb_vector_assign_string_element_len(vec, row, &frame, 1);
        return true;
    }
    case GB_COL_ATTRIBUTES:
        if (!gb_ensure_attributes(st, r->feature)) return false;
        duckdb_vector_assign_string_element_len(vec, row, st->attr.s ? st->attr.s : "", st->attr.l);
        return true;
    case GB_COL_ATTRIBUTES_MAP:
        return gb_fill_attr_map(st, vec, row, r->feature);
    default:
        return true;
    }
}

static void read_genbank_scan(duckdb_function_info info, duckdb_data_chunk output) {
    gb_bind_t *bd = (gb_bind_t *)duckdb_function_get_bind_data(info);
    gb_init_t *st = (gb_init_t *)duckdb_function_get_init_data(info);
    const gb_parser_t *p = &st->parser;
    idx_t capacity = duckdb_vector_size();
    idx_t chunk_cols = duckdb_data_chunk_get_column_count(output);
    char msg[768];

    idx_t row = 0;
    while (row < capacity && st->state != GB_SCAN_DONE) {
        if (st->state == GB_SCAN_PULL || st->row_cursor >= p->n_rows) {
            int next = gb_next_record(st, bd->path, msg, sizeof(msg));
            if (next < 0) {
                duckdb_function_set_error(info, msg);
                duckdb_data_chunk_set_size(output, 0);
                return;
            }
            continue;
        }

        const gb_row_t *r = &p->rows[st->row_cursor++];
        const gb_feature_t *f = &p->feats[r->feature];
        if (f->key.len == 6 && memcmp(gb_str(p, f->key), "source", 6) == 0) continue;

        if (st->n_projected_attributes) {
            for (idx_t c = 0; c < chunk_cols; c++) {
                st->vectors[c] = duckdb_data_chunk_get_vector(output, c);
            }
            if (!gb_find_projected_attributes(st, r->feature)) {
                duckdb_function_set_error(info, "read_genbank: out of memory building attributes");
                duckdb_data_chunk_set_size(output, 0);
                return;
            }
            duckhts_attribute_columns_write(st->projected_attributes, st->n_projected_attributes,
                                            st->vectors, row);
        }
        for (idx_t c = 0; c < chunk_cols; c++) {
            idx_t column = st->column_ids[c];
            if (column >= bd->named_attributes.first_column) continue; /* named column, written above */
            if (!gb_fill_column(st, duckdb_data_chunk_get_vector(output, c), column, row, r)) {
                duckdb_function_set_error(info, "read_genbank: out of memory building attributes");
                duckdb_data_chunk_set_size(output, 0);
                return;
            }
        }
        row++;
    }
    duckdb_data_chunk_set_size(output, row);
}

void register_read_genbank_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "read_genbank");
    duckdb_logical_type vc = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bl = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_table_function_add_parameter(tf, vc);
    duckdb_table_function_add_named_parameter(tf, "attributes_map", bl);
    duckdb_logical_type keys = duckdb_create_list_type(vc);
    duckdb_table_function_add_named_parameter(tf, "attributes", keys);
    duckdb_destroy_logical_type(&keys);
    duckdb_table_function_set_bind(tf, read_genbank_bind);
    duckdb_table_function_set_init(tf, read_genbank_init);
    duckdb_table_function_set_function(tf, read_genbank_scan);
    duckdb_table_function_supports_projection_pushdown(tf, true);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
    duckdb_destroy_logical_type(&vc);
    duckdb_destroy_logical_type(&bl);
}

/* ------------------------------------------------------------------------ */
/* genbank_to_fasta                                                          */
/* ------------------------------------------------------------------------ */

typedef struct {
    char *input;
    char *output;
    int64_t line_width;
    bool overwrite;
} gb2fa_bind_t;

static void gb2fa_bind_destroy(void *data) {
    gb2fa_bind_t *bd = (gb2fa_bind_t *)data;
    if (!bd) return;
    duckdb_free(bd->input);
    duckdb_free(bd->output);
    duckdb_free(bd);
}

static void genbank_to_fasta_bind(duckdb_bind_info info) {
    duckdb_value val = duckdb_bind_get_parameter(info, 0);
    char *input = gb_bind_string(val);
    if (val) duckdb_destroy_value(&val);
    if (!input) {
        duckdb_bind_set_error(info, "genbank_to_fasta requires a file path");
        return;
    }

    gb2fa_bind_t *bd = (gb2fa_bind_t *)duckdb_malloc(sizeof(*bd));
    if (!bd) {
        duckdb_free(input);
        duckdb_bind_set_error(info, "genbank_to_fasta: out of memory");
        return;
    }
    memset(bd, 0, sizeof(*bd));
    bd->input = input;
    bd->line_width = 70;

    val = duckdb_bind_get_named_parameter(info, "output_path");
    if (val && !duckdb_is_null_value(val)) {
        bd->output = gb_bind_string(val);
        if (!bd->output) {
            duckdb_destroy_value(&val);
            duckdb_bind_set_error(info, "genbank_to_fasta: output_path must not be empty");
            gb2fa_bind_destroy(bd);
            return;
        }
    }
    if (val) duckdb_destroy_value(&val);
    if (!bd->output) {
        size_t n = strlen(input);
        bd->output = (char *)duckdb_malloc(n + 4);
        if (!bd->output) {
            duckdb_bind_set_error(info, "genbank_to_fasta: out of memory");
            gb2fa_bind_destroy(bd);
            return;
        }
        memcpy(bd->output, input, n);
        memcpy(bd->output + n, ".fa", 4);
    }

    val = duckdb_bind_get_named_parameter(info, "line_width");
    if (val && !duckdb_is_null_value(val)) bd->line_width = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);
    if (bd->line_width < 1) {
        duckdb_bind_set_error(info, "genbank_to_fasta: line_width must be at least 1");
        gb2fa_bind_destroy(bd);
        return;
    }

    val = duckdb_bind_get_named_parameter(info, "overwrite");
    if (val && !duckdb_is_null_value(val)) bd->overwrite = duckdb_get_bool(val);
    if (val) duckdb_destroy_value(&val);

    duckdb_logical_type bl = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type vc = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bi = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_bind_add_result_column(info, "success", bl);
    duckdb_bind_add_result_column(info, "output_path", vc);
    duckdb_bind_add_result_column(info, "records_written", bi);
    duckdb_destroy_logical_type(&bl);
    duckdb_destroy_logical_type(&vc);
    duckdb_destroy_logical_type(&bi);
    duckdb_bind_set_bind_data(info, bd, gb2fa_bind_destroy);
}

typedef enum {
    GB2FA_STAGED = 0, /* temp file open beside the destination; the scan has not run */
    GB2FA_PUBLISHED,  /* temp file renamed over the destination */
    GB2FA_FAILED      /* error raised; the destructor cleans up */
} gb2fa_stage_t;

typedef enum {
    GB2FA_BETWEEN_RECORDS = 0,
    GB2FA_IN_SEQUENCE /* a defline is written; residues of this record follow */
} gb2fa_record_state_t;

typedef struct {
    htsFile *fp;
    FILE *out;         /* temp file beside the destination */
    char *tmp_path;    /* owned; NULL once published or removed */
    char *dest_path;   /* owned copy of the destination, for the destructor */
    bool created_dest; /* the claim created the destination, so a failure removes it */
    gb2fa_stage_t stage;
    gb2fa_record_state_t record_state;
    kstring_t line;
    gb_parser_t parser;
    int64_t records;
    int64_t col;       /* letters on the current output line */
} gb2fa_init_t;

/* The destructor owns every cleanup, so a query that never reaches the scan
 * (LIMIT 0, an interrupt) still leaves neither a temp file nor an empty claim. */
static void gb2fa_init_destroy(void *data) {
    gb2fa_init_t *st = (gb2fa_init_t *)data;
    if (!st) return;
    if (st->out) fclose(st->out);
    if (st->tmp_path) {
        unlink(st->tmp_path);
        free(st->tmp_path);
    }
    if (st->created_dest && st->stage != GB2FA_PUBLISHED && st->dest_path) unlink(st->dest_path);
    duckdb_free(st->dest_path);
    if (st->fp) hts_close(st->fp);
    free(st->line.s);
    gb_parser_destroy(&st->parser);
    free(st);
}

/* Refuse to write onto the input. Remote inputs do not stat and cannot collide. */
static bool gb2fa_same_file(const char *input, const char *output) {
    struct stat si, so;
    if (stat(input, &si) != 0 || stat(output, &so) != 0) return false;
    return si.st_ino != 0 && si.st_dev == so.st_dev && si.st_ino == so.st_ino;
}

/* Claim the destination: create it exclusively, or open the existing file when
 * overwriting (without truncating it). Returns the descriptor, or -1 with msg. */
static int gb2fa_claim(const gb2fa_bind_t *bd, bool *created, char *msg, size_t msg_len) {
    int flags = O_WRONLY | O_CREAT | O_EXCL;
#ifdef O_BINARY
    flags |= O_BINARY;
#endif
    int fd = open(bd->output, flags, 0666);
    if (fd >= 0) {
        *created = true;
        return fd;
    }
    if (errno != EEXIST) {
        snprintf(msg, msg_len, "genbank_to_fasta: cannot create output %s: %s", bd->output, strerror(errno));
        return -1;
    }
    if (!bd->overwrite) {
        snprintf(msg, msg_len, "genbank_to_fasta: output '%s' already exists (use overwrite := TRUE)", bd->output);
        return -1;
    }
    fd = open(bd->output, flags & ~(O_CREAT | O_EXCL));
    if (fd < 0) {
        snprintf(msg, msg_len, "genbank_to_fasta: cannot open output %s: %s", bd->output, strerror(errno));
        return -1;
    }
    *created = false;
    return fd;
}

static void genbank_to_fasta_init(duckdb_init_info info) {
    gb2fa_bind_t *bd = (gb2fa_bind_t *)duckdb_init_get_bind_data(info);
    char msg[768];
    gb2fa_init_t *st = (gb2fa_init_t *)calloc(1, sizeof(*st));
    if (!st) {
        duckdb_init_set_error(info, "genbank_to_fasta: out of memory");
        return;
    }
    gb_parser_init(&st->parser, GB_MODE_SEQUENCE);
    st->dest_path = duckhts_copy_string(bd->output);
    if (!st->dest_path) {
        duckdb_init_set_error(info, "genbank_to_fasta: out of memory");
        gb2fa_init_destroy(st);
        return;
    }

    st->fp = hts_open(bd->input, "r");
    if (!st->fp) {
        snprintf(msg, sizeof(msg), "genbank_to_fasta: cannot open input %s", bd->input);
        duckdb_init_set_error(info, msg);
        gb2fa_init_destroy(st);
        return;
    }
    if (gb2fa_same_file(bd->input, bd->output)) {
        snprintf(msg, sizeof(msg), "genbank_to_fasta: output '%s' is the input file", bd->output);
        duckdb_init_set_error(info, msg);
        gb2fa_init_destroy(st);
        return;
    }

    int claim = gb2fa_claim(bd, &st->created_dest, msg, sizeof(msg));
    if (claim < 0) {
        duckdb_init_set_error(info, msg);
        gb2fa_init_destroy(st);
        return;
    }
    struct stat claimed;
    bool have_mode = fstat(claim, &claimed) == 0;
    close(claim);

    size_t n = strlen(bd->output);
    st->tmp_path = (char *)malloc(n + 8);
    if (!st->tmp_path) {
        duckdb_init_set_error(info, "genbank_to_fasta: out of memory");
        gb2fa_init_destroy(st);
        return;
    }
    memcpy(st->tmp_path, bd->output, n);
    memcpy(st->tmp_path + n, ".XXXXXX", 8);
    int fd = mkstemp(st->tmp_path);
    if (fd < 0) {
        snprintf(msg, sizeof(msg), "genbank_to_fasta: cannot create temporary output beside %s: %s", bd->output,
                 strerror(errno));
        duckdb_init_set_error(info, msg);
        free(st->tmp_path);
        st->tmp_path = NULL;
        gb2fa_init_destroy(st);
        return;
    }
#ifndef _WIN32
    /* mkstemp creates 0600; the published file should carry the destination's mode. */
    if (have_mode) (void)fchmod(fd, claimed.st_mode & 07777);
#else
    (void)have_mode;
#endif
    st->out = fdopen(fd, "wb");
    if (!st->out) {
        snprintf(msg, sizeof(msg), "genbank_to_fasta: cannot open temporary output: %s", strerror(errno));
        duckdb_init_set_error(info, msg);
        close(fd);
        gb2fa_init_destroy(st);
        return;
    }
    duckdb_init_set_max_threads(info, 1);
    duckdb_init_set_init_data(info, st, gb2fa_init_destroy);
}

static bool gb2fa_put(gb2fa_init_t *st, const char *s, size_t n) {
    return fwrite(s, 1, n, st->out) == n;
}

/* Write residues wrapped at line_width, continuing the current line. */
static bool gb2fa_write_residues(gb2fa_init_t *st, int64_t line_width, const char *s, size_t n) {
    while (n) {
        size_t room = (size_t)(line_width - st->col);
        size_t w = n < room ? n : room;
        if (!gb2fa_put(st, s, w)) return false;
        s += w;
        n -= w;
        st->col += (int64_t)w;
        if (st->col == line_width) {
            if (!gb2fa_put(st, "\n", 1)) return false;
            st->col = 0;
        }
    }
    return true;
}

static bool gb2fa_write_defline(gb2fa_init_t *st) {
    const gb_parser_t *p = &st->parser;
    gb_span_t contig = gb_contig(p);
    if (!gb2fa_put(st, ">", 1)) return false;
    if (contig.len ? !gb2fa_put(st, gb_str(p, contig), contig.len) : !gb2fa_put(st, "unknown", 7)) return false;
    /* NCBI's FASTA export drops the period that ends every DEFINITION. */
    size_t def_len = p->definition.len;
    if (def_len && gb_str(p, p->definition)[def_len - 1] == '.') def_len--;
    if (def_len) {
        if (!gb2fa_put(st, " ", 1) || !gb2fa_put(st, gb_str(p, p->definition), def_len)) return false;
    }
    return gb2fa_put(st, "\n", 1);
}

/* Stream the whole input into the temp file. Returns false with msg on any failure. */
static bool gb2fa_convert(gb2fa_init_t *st, const gb2fa_bind_t *bd, char *msg, size_t msg_len) {
    for (;;) {
        int n = hts_getline(st->fp, KS_SEP_LINE, &st->line);
        if (n < -1) {
            snprintf(msg, msg_len, "genbank_to_fasta: %s: read failed after line %ld%s%s", bd->input,
                     st->parser.line_no, errno ? ": " : "", errno ? strerror(errno) : "");
            return false;
        }
        if (n == -1) {
            if (gb_parser_finish(&st->parser) != GB_OK) {
                gb_format_error(msg, msg_len, "genbank_to_fasta", bd->input, &st->parser.err);
                return false;
            }
            return true;
        }
        bool ok = true;
        switch (gb_parser_feed(&st->parser, st->line.s, st->line.l)) {
        case GB_FEED_ORIGIN:
            ok = gb2fa_write_defline(st);
            st->record_state = GB2FA_IN_SEQUENCE;
            st->col = 0;
            break;
        case GB_FEED_RESIDUES:
            ok = gb2fa_write_residues(st, bd->line_width, st->parser.residues.s, st->parser.residues.l);
            break;
        case GB_FEED_RECORD:
            if (st->record_state == GB2FA_IN_SEQUENCE) {
                if (st->col > 0) ok = gb2fa_put(st, "\n", 1);
                st->records++;
                st->record_state = GB2FA_BETWEEN_RECORDS;
                st->col = 0;
            }
            break;
        case GB_FEED_ERROR:
            gb_format_error(msg, msg_len, "genbank_to_fasta", bd->input, &st->parser.err);
            return false;
        default:
            break;
        }
        if (!ok) {
            snprintf(msg, msg_len, "genbank_to_fasta: write failed on %s: %s", bd->output, strerror(errno));
            return false;
        }
    }
}

/* Move the finished temp file over the destination. */
static bool gb2fa_publish(gb2fa_init_t *st, const gb2fa_bind_t *bd, char *msg, size_t msg_len) {
    /* fclose flushes: a disk that filled during the last buffered write fails here. */
    FILE *out = st->out;
    st->out = NULL;
    if (fclose(out) != 0) {
        snprintf(msg, msg_len, "genbank_to_fasta: write failed on %s: %s", bd->output, strerror(errno));
        return false;
    }
#ifdef _WIN32
    bool moved = MoveFileExA(st->tmp_path, bd->output, MOVEFILE_REPLACE_EXISTING) != 0;
#else
    bool moved = rename(st->tmp_path, bd->output) == 0;
#endif
    if (!moved) {
        snprintf(msg, msg_len, "genbank_to_fasta: cannot publish %s: %s", bd->output, strerror(errno));
        return false;
    }
    free(st->tmp_path);
    st->tmp_path = NULL;
    st->stage = GB2FA_PUBLISHED;
    return true;
}

static void genbank_to_fasta_scan(duckdb_function_info info, duckdb_data_chunk output) {
    gb2fa_bind_t *bd = (gb2fa_bind_t *)duckdb_function_get_bind_data(info);
    gb2fa_init_t *st = (gb2fa_init_t *)duckdb_function_get_init_data(info);
    if (st->stage != GB2FA_STAGED) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    char msg[768];
    bool ok = gb2fa_convert(st, bd, msg, sizeof(msg));
    if (ok && st->records == 0) {
        snprintf(msg, sizeof(msg), "genbank_to_fasta: no sequence records found in %s", bd->input);
        ok = false;
    }
    if (ok) ok = gb2fa_publish(st, bd, msg, sizeof(msg));
    if (!ok) {
        /* The destructor removes the temp file and any destination this query created. */
        st->stage = GB2FA_FAILED;
        duckdb_function_set_error(info, msg);
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    ((bool *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 0)))[0] = true;
    duckdb_vector_assign_string_element(duckdb_data_chunk_get_vector(output, 1), 0, bd->output);
    ((int64_t *)duckdb_vector_get_data(duckdb_data_chunk_get_vector(output, 2)))[0] = st->records;
    duckdb_data_chunk_set_size(output, 1);
}

void register_genbank_to_fasta_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "genbank_to_fasta");
    duckdb_logical_type vc = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type it = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type bl = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_table_function_add_parameter(tf, vc);
    duckdb_table_function_add_named_parameter(tf, "output_path", vc);
    duckdb_table_function_add_named_parameter(tf, "line_width", it);
    duckdb_table_function_add_named_parameter(tf, "overwrite", bl);
    duckdb_table_function_set_bind(tf, genbank_to_fasta_bind);
    duckdb_table_function_set_init(tf, genbank_to_fasta_init);
    duckdb_table_function_set_function(tf, genbank_to_fasta_scan);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
    duckdb_destroy_logical_type(&vc);
    duckdb_destroy_logical_type(&it);
    duckdb_destroy_logical_type(&bl);
}
