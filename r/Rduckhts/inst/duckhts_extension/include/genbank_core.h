#ifndef DUCKHTS_GENBANK_CORE_H
#define DUCKHTS_GENBANK_CORE_H

/* GenBank flat-file parsing core. No DuckDB, no file I/O: the caller pushes
 * one line at a time and reads the record DOM after the terminator.
 *
 * Lifecycle: gb_parser_init, then gb_parser_feed per line, then
 * gb_parser_finish at end of input, then gb_parser_destroy. After a feed
 * returns GB_FEED_RECORD the DOM arrays and header spans describe that
 * record until the next feed call, which resets them. Every span is a
 * borrowed, NUL-terminated byte range inside `arena`, so gb_str() returns a
 * C string. Arena pointers move as the record grows; hold spans, not
 * pointers, until GB_FEED_RECORD.
 *
 * Memory is proportional to one record's kept FEATURES text plus the DOM.
 * ORIGIN letters are never stored. In GB_MODE_SEQUENCE the feature table is
 * not stored either; in GB_MODE_FEATURES sequence lines are skipped.
 *
 * Table layout follows BioPython's GenBank scanner: a feature key may start
 * anywhere before column 22, everything from column 22 is location text,
 * qualifiers or continuation, and a FEATURES table must reach a sequence
 * section (ORIGIN, CONTIG, BASE COUNT, WGS, TSA or TLS) before the record
 * terminator; any other keyword in that position is a syntax error. */

#include <stddef.h>
#include <stdint.h>

#include <htslib/kstring.h>

typedef struct {
    size_t off, len;
} gb_span_t;

typedef enum {
    GB_OK = 0,
    GB_ERR_NOMEM,
    GB_ERR_SYNTAX,      /* malformed line, qualifier or location */
    GB_ERR_UNSUPPORTED, /* legal INSDC form this reader declines: one-of, bond, gap, a.b */
    GB_ERR_TRUNCATED,   /* end of input inside a record */
    GB_ERR_EMPTY        /* end of input with no LOCUS record */
} gb_status_t;

typedef struct {
    gb_status_t code;
    long line; /* 1-based input line, 0 when not line-specific */
    char msg[256];
} gb_error_t;

typedef enum {
    GB_MODE_FEATURES = 0,
    GB_MODE_SEQUENCE = 1
} gb_mode_t;

/* Set in gb_parser_t.flags before the first feed. */
enum {
    GB_FLAG_DROP_TRANSLATION = 1u /* never store /translation values */
};

typedef enum {
    GB_FEED_CONTINUE = 0,
    GB_FEED_ORIGIN,   /* ORIGIN keyword: header spans are final, sequence follows */
    GB_FEED_RESIDUES, /* GB_MODE_SEQUENCE: `residues` holds this line's letters */
    GB_FEED_RECORD,   /* // seen: the record DOM is readable */
    GB_FEED_ERROR     /* `err` is set; the parser accepts no further input */
} gb_feed_t;

typedef enum {
    GB_OP_NONE = 0,
    GB_OP_JOIN,
    GB_OP_ORDER
} gb_op_t;

typedef enum {
    GB_TOPOLOGY_UNKNOWN = 0,
    GB_TOPOLOGY_LINEAR,
    GB_TOPOLOGY_CIRCULAR
} gb_topology_t;

/* A span whose end precedes its start on a circular record wraps the origin
 * and yields two segments, start..length then 1..end (reversed under an
 * element-level complement), as BioPython resolves it. On a linear record
 * it is an error. A between site n^m requires m == n + 1, or m == 1 with n
 * the sequence length, and is reported as the zero-length site start == end
 * == n, which is how GFF3 places an insertion site. */
typedef struct {
    int64_t start, end; /* 1-based inclusive; n^m gives start == end == n */
    gb_span_t remote;   /* accession before ':'; len 0 when local */
    uint8_t complement; /* element-level complement( ) */
    uint8_t partial5;   /* '<' on start */
    uint8_t partial3;   /* '>' on end */
    uint8_t between;    /* n^m */
} gb_seg_t;

typedef struct {
    gb_span_t key;
    gb_span_t value;   /* quotes stripped, "" unescaped, lines joined by one space */
    uint8_t has_value; /* 0 for a valueless qualifier such as /pseudo */
} gb_qual_t;

typedef struct {
    gb_span_t key;
    size_t first_value, n_values; /* into attr_values: qualifier indexes */
} gb_attr_t;

typedef struct {
    size_t feature;
    int64_t start, end;
    char strand;  /* '+' or '-' */
    int8_t phase; /* 0..2 on CDS rows, -1 otherwise */
    gb_span_t remote;
} gb_row_t;

typedef struct {
    gb_span_t key;      /* feature key, e.g. CDS */
    gb_span_t loc_text; /* location as written, continuation lines concatenated */
    long line;          /* input line of the feature start */
    long ordinal;       /* 0-based count across the whole input, for GFF3 IDs */
    gb_op_t op;
    uint8_t outer_complement;
    uint8_t is_gene;
    size_t first_seg, n_segs;   /* into segs, file order */
    size_t first_qual, n_quals; /* into quals, file order */
    size_t first_attr, n_attrs; /* into attrs, first-occurrence order */
    size_t first_row, n_rows;   /* into rows, biological order */
    gb_span_t locus;            /* /locus_tag, else /gene; len 0 when neither */
    gb_span_t name;             /* GFF3 Name: gene, product, label, locus_tag, else key */
    gb_span_t parent_locus;     /* locus of the linked gene; len 0 when none */
} gb_feature_t;

/* Parser state axes. Private to the core; exposed only so the struct is flat. */
typedef enum {
    GB_SECTION_BETWEEN = 0, /* before LOCUS or after // */
    GB_SECTION_HEADER,
    GB_SECTION_FEATURES,
    GB_SECTION_OTHER,       /* CONTIG, BASE COUNT and similar after FEATURES */
    GB_SECTION_SEQUENCE,    /* after ORIGIN */
    GB_SECTION_FAILED
} gb_section_t;

typedef enum { GB_HEADER_NONE = 0, GB_HEADER_DEFINITION } gb_header_field_t;
typedef enum { GB_PART_NONE = 0, GB_PART_LOCATION, GB_PART_QUALIFIERS } gb_feature_part_t;
typedef enum { GB_VALUE_NONE = 0, GB_VALUE_UNQUOTED, GB_VALUE_QUOTE_OPEN, GB_VALUE_QUOTE_CLOSED } gb_value_form_t;
typedef enum { GB_SINK_STORE = 0, GB_SINK_DISCARD } gb_sink_t;

typedef struct {
    kstring_t arena;

    gb_feature_t *feats;
    size_t n_feats, cap_feats;
    gb_qual_t *quals;
    size_t n_quals, cap_quals;
    gb_seg_t *segs;
    size_t n_segs, cap_segs;
    gb_attr_t *attrs;
    size_t n_attrs, cap_attrs;
    size_t *attr_values;
    size_t n_attr_values, cap_attr_values;
    gb_row_t *rows;
    size_t n_rows, cap_rows;

    gb_span_t locus, accession, version, definition;
    int64_t seq_length;     /* LOCUS length in bp, 0 when absent */
    gb_topology_t topology; /* LOCUS topology */
    kstring_t residues; /* GB_MODE_SEQUENCE: uppercase letters of the last sequence line */

    gb_mode_t mode;
    unsigned flags;
    gb_section_t section;
    gb_header_field_t header_field;
    gb_feature_part_t feature_part;
    gb_value_form_t value_form;
    gb_sink_t sink;
    long line_no;
    long next_ordinal;
    size_t n_records;
    gb_error_t err;
} gb_parser_t;

void gb_parser_init(gb_parser_t *p, gb_mode_t mode);
void gb_parser_destroy(gb_parser_t *p);

/* `line` excludes the newline; a trailing CR is ignored. */
gb_feed_t gb_parser_feed(gb_parser_t *p, const char *line, size_t len);

/* End of input. GB_ERR_TRUNCATED inside a record, GB_ERR_EMPTY with no record. */
gb_status_t gb_parser_finish(gb_parser_t *p);

static inline const char *gb_str(const gb_parser_t *p, gb_span_t s) {
    return s.len ? p->arena.s + s.off : "";
}

/* GFF3 column-9 output. Emitted attributes are ID, Name, Parent when linked,
 * then the grouped qualifiers with comma-joined, percent-encoded values.
 * A valueless qualifier contributes "true". */
size_t gb_feature_attr_count(const gb_parser_t *p, size_t feature);
gb_status_t gb_feature_attr_at(const gb_parser_t *p, size_t feature, size_t index,
                               kstring_t *key, kstring_t *value);
/* Key and value of one attribute, emitted separately so a caller that needs only
 * some pairs can skip encoding the values it discards. attr_at is both. */
gb_status_t gb_feature_attr_key(const gb_parser_t *p, size_t feature, size_t index, kstring_t *key);
gb_status_t gb_feature_attr_value(const gb_parser_t *p, size_t feature, size_t index, kstring_t *value);
gb_status_t gb_feature_attributes(const gb_parser_t *p, size_t feature, kstring_t *out);

#endif
