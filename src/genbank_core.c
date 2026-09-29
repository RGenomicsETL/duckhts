/* GenBank flat-file parsing core. See genbank_core.h for the contract.
 *
 * The parser is a line classifier driven by five independent enums:
 *   section      which top-level block of the record we are in
 *   header_field which header keyword owns indented continuation lines
 *   feature_part whether a feature is open and which part continuation extends
 *   value_form   how the current qualifier value was written
 *   sink         whether the current qualifier value is stored or discarded
 * Every stored token is appended to one arena and referenced by span. A token
 * still being extended by continuation lines is always the arena tail. */

#include "genbank_core.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/hts_defs.h>
#include <htslib/khash_str2int.h>

#define GB_TEXT_COLUMN 21u /* INSDC column 22: location, qualifiers, continuations */

/* ------------------------------------------------------------------------ */
/* Small helpers                                                             */
/* ------------------------------------------------------------------------ */

static int gb_reserve(void **data, size_t *cap, size_t need, size_t width) {
    if (need <= *cap) return 1;
    size_t new_cap = *cap ? *cap : 16;
    while (new_cap < need) {
        if (new_cap > SIZE_MAX / 2) return 0;
        new_cap *= 2;
    }
    if (new_cap > SIZE_MAX / width) return 0;
    void *grown = realloc(*data, new_cap * width);
    if (!grown) return 0;
    *data = grown;
    *cap = new_cap;
    return 1;
}

#define GB_PUSH(p, arr, n, cap, item)                                                       \
    (gb_reserve((void **)&(p)->arr, &(p)->cap, (p)->n + 1, sizeof(*(p)->arr))              \
         ? ((p)->arr[(p)->n++] = (item), 1)                                                 \
         : 0)

static size_t gb_rtrim(const char *s, size_t n) {
    while (n && (s[n - 1] == ' ' || s[n - 1] == '\t')) n--;
    return n;
}

static size_t gb_skip_spaces(const char *s, size_t n, size_t i) {
    while (i < n && (s[i] == ' ' || s[i] == '\t')) i++;
    return i;
}

static size_t gb_token_end(const char *s, size_t n, size_t i) {
    while (i < n && s[i] != ' ' && s[i] != '\t') i++;
    return i;
}

static int gb_span_eq(const gb_parser_t *p, gb_span_t a, const char *lit) {
    size_t n = strlen(lit);
    return a.len == n && memcmp(gb_str(p, a), lit, n) == 0;
}

static int gb_span_eq_span(const gb_parser_t *p, gb_span_t a, gb_span_t b) {
    return a.len == b.len && memcmp(gb_str(p, a), gb_str(p, b), a.len) == 0;
}

static gb_feed_t gb_fail(gb_parser_t *p, gb_status_t code, const char *fmt, ...) HTS_FORMAT(HTS_PRINTF_FMT, 3, 4);

static gb_feed_t gb_fail(gb_parser_t *p, gb_status_t code, const char *fmt, ...) {
    va_list ap;
    p->section = GB_SECTION_FAILED;
    p->err.code = code;
    p->err.line = p->line_no;
    va_start(ap, fmt);
    vsnprintf(p->err.msg, sizeof(p->err.msg), fmt, ap);
    va_end(ap);
    return GB_FEED_ERROR;
}

static gb_feed_t gb_oom(gb_parser_t *p) {
    return gb_fail(p, GB_ERR_NOMEM, "out of memory");
}

/* Append a new NUL-terminated token to the arena. */
static int gb_arena_put(gb_parser_t *p, const char *s, size_t n, gb_span_t *out) {
    size_t off = p->arena.l;
    if (kputsn(s, n, &p->arena) == EOF || kputc('\0', &p->arena) == EOF) return 0;
    out->off = off;
    out->len = n;
    return 1;
}

/* Extend the token at the arena tail. */
static int gb_arena_extend(gb_parser_t *p, gb_span_t *span, const char *s, size_t n) {
    if (p->arena.l != span->off + span->len + 1) return 0; /* not the tail: invariant broken */
    p->arena.l -= 1;
    if (kputsn(s, n, &p->arena) == EOF || kputc('\0', &p->arena) == EOF) return 0;
    span->len += n;
    return 1;
}

/* ------------------------------------------------------------------------ */
/* Lifecycle                                                                 */
/* ------------------------------------------------------------------------ */

void gb_parser_init(gb_parser_t *p, gb_mode_t mode) {
    memset(p, 0, sizeof(*p));
    p->mode = mode;
    p->section = GB_SECTION_BETWEEN;
}

void gb_parser_destroy(gb_parser_t *p) {
    free(p->arena.s);
    free(p->residues.s);
    free(p->feats);
    free(p->quals);
    free(p->segs);
    free(p->attrs);
    free(p->attr_values);
    free(p->rows);
    memset(p, 0, sizeof(*p));
    p->section = GB_SECTION_FAILED;
}

static void gb_reset_record(gb_parser_t *p) {
    p->arena.l = 0;
    p->residues.l = 0;
    p->n_feats = p->n_quals = p->n_segs = p->n_attrs = p->n_attr_values = p->n_rows = 0;
    p->locus.off = p->locus.len = 0;
    p->accession.off = p->accession.len = 0;
    p->version.off = p->version.len = 0;
    p->definition.off = p->definition.len = 0;
    p->seq_length = 0;
    p->topology = GB_TOPOLOGY_UNKNOWN;
    p->header_field = GB_HEADER_NONE;
    p->feature_part = GB_PART_NONE;
    p->value_form = GB_VALUE_NONE;
    p->sink = GB_SINK_STORE;
}

/* ------------------------------------------------------------------------ */
/* Feature table                                                             */
/* ------------------------------------------------------------------------ */

static gb_feature_t *gb_cur_feature(gb_parser_t *p) {
    return p->n_feats ? &p->feats[p->n_feats - 1] : NULL;
}

static gb_qual_t *gb_cur_qual(gb_parser_t *p) {
    return p->n_quals ? &p->quals[p->n_quals - 1] : NULL;
}

/* ------------------------------------------------------------------------ */
/* Location grammar (INSDC 3.4.2)                                            */
/*                                                                           */
/*   location := 'complement(' inner ')' | inner                             */
/*   inner    := ('join(' | 'order(') element (',' element)* ')' | element   */
/*   element  := 'complement(' span ')' | span                               */
/*   span     := [accession ':'] point ('..' point | '^' point)?             */
/*   point    := ['<' | '>'] digits                                          */
/*                                                                           */
/* join/order never nest, and complement never nests, so the grammar itself  */
/* bounds recursion. The text lives in the arena, so the cursor holds an     */
/* offset rather than a pointer: appending segments may move the arena.      */
/* ------------------------------------------------------------------------ */

typedef struct {
    size_t off, n, i; /* text at arena + off, length n, cursor i */
    gb_feature_t *f;
} gb_loc_cursor_t;

static char gb_lc_at(const gb_parser_t *p, const gb_loc_cursor_t *c, size_t k) {
    return c->i + k < c->n ? p->arena.s[c->off + c->i + k] : '\0';
}

static int gb_lc_starts(const gb_parser_t *p, const gb_loc_cursor_t *c, const char *lit) {
    size_t n = strlen(lit);
    return c->i + n <= c->n && memcmp(p->arena.s + c->off + c->i, lit, n) == 0;
}

/* Location errors surface when the feature closes; report the feature's own line. */
static gb_feed_t gb_loc_fail(gb_parser_t *p, const gb_loc_cursor_t *c, gb_status_t code, const char *what) {
    /* A location can run to thousands of characters; the message keeps a prefix. */
    gb_fail(p, code, "feature %s at line %ld: location '%.80s%s': %s", gb_str(p, c->f->key), c->f->line,
            gb_str(p, c->f->loc_text), c->f->loc_text.len > 80 ? "..." : "", what);
    p->err.line = c->f->line;
    return GB_FEED_ERROR;
}

static gb_feed_t gb_parse_point(gb_parser_t *p, gb_loc_cursor_t *c, int64_t *out, uint8_t *partial) {
    char ch = gb_lc_at(p, c, 0);
    if (ch == '<' || ch == '>') {
        *partial = 1;
        c->i++;
        ch = gb_lc_at(p, c, 0);
    }
    if (ch < '0' || ch > '9') return gb_loc_fail(p, c, GB_ERR_SYNTAX, "expected a coordinate");
    int64_t v = 0;
    while (ch >= '0' && ch <= '9') {
        int64_t d = ch - '0';
        if (v > (INT64_MAX - d) / 10) return gb_loc_fail(p, c, GB_ERR_SYNTAX, "coordinate overflows 64 bits");
        v = v * 10 + d;
        c->i++;
        ch = gb_lc_at(p, c, 0);
    }
    if (v == 0) return gb_loc_fail(p, c, GB_ERR_SYNTAX, "coordinates are 1-based");
    *out = v;
    return GB_FEED_CONTINUE;
}

/* BioPython's reference pattern: a letter, then letters, digits, '_', '.' or '|'. */
static int gb_is_accession_char(char ch) {
    return (ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') || (ch >= '0' && ch <= '9') || ch == '_' ||
           ch == '.' || ch == '|';
}

static int gb_is_letter(char ch) {
    return (ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z');
}

static gb_feed_t gb_parse_span(gb_parser_t *p, gb_loc_cursor_t *c, uint8_t complement) {
    gb_seg_t seg;
    memset(&seg, 0, sizeof(seg));
    seg.complement = complement;

    size_t j = 0;
    if (gb_is_letter(gb_lc_at(p, c, 0)))
        while (gb_is_accession_char(gb_lc_at(p, c, j))) j++;
    if (gb_lc_at(p, c, j) == ':') {
        if (j == 0) return gb_loc_fail(p, c, GB_ERR_SYNTAX, "empty remote accession");
        /* The accession is copied from the arena into the arena: reserve first so the
         * copy cannot reallocate its own source. */
        if (ks_resize(&p->arena, p->arena.l + j + 2) < 0) return gb_oom(p);
        if (!gb_arena_put(p, p->arena.s + c->off + c->i, j, &seg.remote)) return gb_oom(p);
        c->i += j + 1;
    }

    gb_feed_t r = gb_parse_point(p, c, &seg.start, &seg.partial5);
    if (r == GB_FEED_ERROR) return r;
    seg.end = seg.start;
    if (gb_lc_at(p, c, 0) == '.' && gb_lc_at(p, c, 1) == '.') {
        c->i += 2;
    } else if (gb_lc_at(p, c, 0) == '^') {
        seg.between = 1;
        c->i += 1;
    } else if (gb_lc_at(p, c, 0) == '.') {
        return gb_loc_fail(p, c, GB_ERR_UNSUPPORTED, "single base within a range (a.b) is not supported");
    } else {
        if (!GB_PUSH(p, segs, n_segs, cap_segs, seg)) return gb_oom(p);
        return GB_FEED_CONTINUE;
    }
    r = gb_parse_point(p, c, &seg.end, &seg.partial3);
    if (r == GB_FEED_ERROR) return r;
    if (seg.between) {
        int junction = seg.end == 1 && p->seq_length > 0 && seg.start == p->seq_length;
        if (seg.end != seg.start + 1 && !junction)
            return gb_loc_fail(p, c, GB_ERR_SYNTAX, "between positions must be adjacent");
        /* GFF3: a zero-length site has start == end at the base to its left. */
        seg.end = seg.start;
        if (!GB_PUSH(p, segs, n_segs, cap_segs, seg)) return gb_oom(p);
        return GB_FEED_CONTINUE;
    }
    if (seg.end >= seg.start) {
        if (!GB_PUSH(p, segs, n_segs, cap_segs, seg)) return gb_oom(p);
        return GB_FEED_CONTINUE;
    }
    /* end before start: an origin-spanning span on a circular record. */
    if (p->topology != GB_TOPOLOGY_CIRCULAR || p->seq_length <= 0 || seg.start > p->seq_length)
        return gb_loc_fail(p, c, GB_ERR_SYNTAX, "end precedes start on a non-circular record");
    gb_seg_t head = seg, tail = seg;
    head.end = p->seq_length;
    head.partial3 = 0;
    tail.start = 1;
    tail.partial5 = 0;
    gb_seg_t first = seg.complement ? tail : head, second = seg.complement ? head : tail;
    if (!GB_PUSH(p, segs, n_segs, cap_segs, first) || !GB_PUSH(p, segs, n_segs, cap_segs, second)) return gb_oom(p);
    return GB_FEED_CONTINUE;
}

static gb_feed_t gb_parse_element(gb_parser_t *p, gb_loc_cursor_t *c, int allow_complement) {
    if (gb_lc_starts(p, c, "complement(")) {
        if (!allow_complement) return gb_loc_fail(p, c, GB_ERR_UNSUPPORTED, "nested complement");
        c->i += 11;
        gb_feed_t r = gb_parse_span(p, c, 1);
        if (r == GB_FEED_ERROR) return r;
        if (gb_lc_at(p, c, 0) != ')') return gb_loc_fail(p, c, GB_ERR_SYNTAX, "expected ')' after complement");
        c->i++;
        return GB_FEED_CONTINUE;
    }
    if (gb_lc_starts(p, c, "join(") || gb_lc_starts(p, c, "order("))
        return gb_loc_fail(p, c, GB_ERR_UNSUPPORTED, "nested join/order");
    if (gb_lc_starts(p, c, "one-of(") || gb_lc_starts(p, c, "gap(") || gb_lc_starts(p, c, "bond("))
        return gb_loc_fail(p, c, GB_ERR_UNSUPPORTED, "one-of, gap and bond are not supported");
    return gb_parse_span(p, c, 0);
}

static gb_feed_t gb_parse_inner(gb_parser_t *p, gb_loc_cursor_t *c, int allow_complement) {
    if (gb_lc_starts(p, c, "join(")) {
        c->f->op = GB_OP_JOIN;
        c->i += 5;
    } else if (gb_lc_starts(p, c, "order(")) {
        c->f->op = GB_OP_ORDER;
        c->i += 6;
    } else {
        return gb_parse_element(p, c, allow_complement);
    }
    for (;;) {
        gb_feed_t r = gb_parse_element(p, c, allow_complement);
        if (r == GB_FEED_ERROR) return r;
        char ch = gb_lc_at(p, c, 0);
        c->i++;
        if (ch == ',') continue;
        if (ch == ')') return GB_FEED_CONTINUE;
        return gb_loc_fail(p, c, GB_ERR_SYNTAX, "expected ',' or ')' in list");
    }
}

static gb_feed_t gb_parse_location(gb_parser_t *p, gb_feature_t *f) {
    gb_loc_cursor_t c = {f->loc_text.off, f->loc_text.len, 0, f};
    f->first_seg = p->n_segs;
    gb_feed_t r;
    if (gb_lc_starts(p, &c, "complement(")) {
        f->outer_complement = 1;
        c.i += 11;
        r = gb_parse_inner(p, &c, 0);
        if (r == GB_FEED_ERROR) return r;
        if (gb_lc_at(p, &c, 0) != ')') return gb_loc_fail(p, &c, GB_ERR_SYNTAX, "expected ')' after complement");
        c.i++;
    } else {
        r = gb_parse_inner(p, &c, 1);
        if (r == GB_FEED_ERROR) return r;
    }
    if (c.i != c.n) return gb_loc_fail(p, &c, GB_ERR_SYNTAX, "unexpected text after location");
    f->n_segs = p->n_segs - f->first_seg;
    return GB_FEED_CONTINUE;
}

static gb_feed_t gb_close_feature(gb_parser_t *p) {
    if (p->feature_part == GB_PART_NONE) return GB_FEED_CONTINUE;
    gb_feature_t *f = gb_cur_feature(p);
    if (p->value_form == GB_VALUE_QUOTE_OPEN) {
        gb_fail(p, GB_ERR_SYNTAX, "feature %s at line %ld: unterminated quoted qualifier value", gb_str(p, f->key),
                f->line);
        p->err.line = f->line;
        return GB_FEED_ERROR;
    }
    p->feature_part = GB_PART_NONE;
    return gb_parse_location(p, f);
}

/* Append value bytes to the current qualifier, or discard them. */
static int gb_value_append(gb_parser_t *p, const char *s, size_t n) {
    if (p->sink == GB_SINK_DISCARD) return 1;
    return gb_arena_extend(p, &gb_cur_qual(p)->value, s, n);
}

/* Scan quoted text: "" is a literal quote, a lone quote closes and must end the line. */
static gb_feed_t gb_scan_quoted(gb_parser_t *p, const char *s, size_t n) {
    for (size_t i = 0; i < n; i++) {
        if (s[i] != '"') {
            if (!gb_value_append(p, s + i, 1)) return gb_oom(p);
            continue;
        }
        if (i + 1 < n && s[i + 1] == '"') {
            if (!gb_value_append(p, "\"", 1)) return gb_oom(p);
            i++;
            continue;
        }
        if (i + 1 != n) return gb_fail(p, GB_ERR_SYNTAX, "text after the closing quote of a qualifier value");
        p->value_form = GB_VALUE_QUOTE_CLOSED;
        return GB_FEED_CONTINUE;
    }
    p->value_form = GB_VALUE_QUOTE_OPEN;
    return GB_FEED_CONTINUE;
}

static gb_feed_t gb_feature_line(gb_parser_t *p, const char *s, size_t n, size_t key_start) {
    gb_feed_t closed = gb_close_feature(p);
    if (closed == GB_FEED_ERROR) return closed;

    size_t key_end = gb_token_end(s, n, key_start);
    size_t loc_start = gb_skip_spaces(s, n, key_end);
    size_t loc_end = gb_rtrim(s, n);
    if (loc_start >= loc_end)
        return gb_fail(p, GB_ERR_SYNTAX, "feature %.*s has no location", (int)(key_end - key_start), s + key_start);

    gb_feature_t f;
    memset(&f, 0, sizeof(f));
    f.line = p->line_no;
    f.ordinal = p->next_ordinal++;
    f.first_qual = p->n_quals;
    f.first_seg = p->n_segs;
    if (!gb_arena_put(p, s + key_start, key_end - key_start, &f.key)) return gb_oom(p);
    if (!gb_arena_put(p, s + loc_start, loc_end - loc_start, &f.loc_text)) return gb_oom(p);
    if (!GB_PUSH(p, feats, n_feats, cap_feats, f)) return gb_oom(p);
    p->feature_part = GB_PART_LOCATION;
    p->value_form = GB_VALUE_NONE;
    p->sink = GB_SINK_STORE;
    return GB_FEED_CONTINUE;
}

/* `s` starts just after the '/'. */
static gb_feed_t gb_qualifier_line(gb_parser_t *p, const char *s, size_t n) {
    if (p->feature_part == GB_PART_NONE) return gb_fail(p, GB_ERR_SYNTAX, "qualifier before any feature");
    p->feature_part = GB_PART_QUALIFIERS;

    size_t eq = 0;
    while (eq < n && s[eq] != '=') eq++;
    size_t key_len = gb_rtrim(s, eq);
    if (key_len == 0) return gb_fail(p, GB_ERR_SYNTAX, "qualifier without a name");

    int drop = (p->flags & GB_FLAG_DROP_TRANSLATION) && key_len == 11 && memcmp(s, "translation", 11) == 0;
    p->sink = drop ? GB_SINK_DISCARD : GB_SINK_STORE;
    if (!drop) {
        gb_qual_t q;
        memset(&q, 0, sizeof(q));
        if (!gb_arena_put(p, s, key_len, &q.key)) return gb_oom(p);
        if (!gb_arena_put(p, "", 0, &q.value)) return gb_oom(p);
        q.has_value = (uint8_t)(eq < n);
        if (!GB_PUSH(p, quals, n_quals, cap_quals, q)) return gb_oom(p);
        gb_cur_feature(p)->n_quals++;
    }

    if (eq >= n) {
        p->value_form = GB_VALUE_NONE;
        return GB_FEED_CONTINUE;
    }
    size_t vstart = gb_skip_spaces(s, n, eq + 1); /* "/key= value" is written by some tools */
    const char *v = s + vstart;
    size_t vn = gb_rtrim(v, n - vstart);
    if (vn && v[0] == '"') return gb_scan_quoted(p, v + 1, vn - 1);
    p->value_form = GB_VALUE_UNQUOTED;
    if (!gb_value_append(p, v, vn)) return gb_oom(p);
    return GB_FEED_CONTINUE;
}

/* `s` is the trimmed text of a continuation line at the text column. */
static gb_feed_t gb_continuation_line(gb_parser_t *p, const char *s, size_t n) {
    switch (p->feature_part) {
    case GB_PART_LOCATION:
        if (!gb_arena_extend(p, &gb_cur_feature(p)->loc_text, s, n)) return gb_oom(p);
        return GB_FEED_CONTINUE;
    case GB_PART_QUALIFIERS:
        switch (p->value_form) {
        case GB_VALUE_QUOTE_OPEN:
            if (!gb_value_append(p, " ", 1)) return gb_oom(p);
            return gb_scan_quoted(p, s, n);
        case GB_VALUE_UNQUOTED:
            if (!gb_value_append(p, " ", 1) || !gb_value_append(p, s, n)) return gb_oom(p);
            return GB_FEED_CONTINUE;
        default:
            return gb_fail(p, GB_ERR_SYNTAX, "continuation line after a complete qualifier");
        }
    default:
        return gb_fail(p, GB_ERR_SYNTAX, "continuation line before any feature");
    }
}

static gb_feed_t gb_features_indented(gb_parser_t *p, const char *s, size_t n) {
    if (p->mode == GB_MODE_SEQUENCE) return GB_FEED_CONTINUE;

    size_t lead = gb_skip_spaces(s, n, 0);
    if (lead >= n) return GB_FEED_CONTINUE; /* blank */
    if (lead >= GB_TEXT_COLUMN) {
        size_t end = gb_rtrim(s, n);
        if (s[lead] == '/' && p->value_form != GB_VALUE_QUOTE_OPEN)
            return gb_qualifier_line(p, s + lead + 1, end - lead - 1);
        return gb_continuation_line(p, s + lead, end - lead);
    }
    return gb_feature_line(p, s, n, lead); /* a key anywhere before the text column */
}

/* ------------------------------------------------------------------------ */
/* Record resolution (at the terminator)                                     */
/* ------------------------------------------------------------------------ */

/* First qualifier with this key that carries a value. */
static int gb_qual_value(const gb_parser_t *p, const gb_feature_t *f, const char *key, gb_span_t *out) {
    for (size_t i = 0; i < f->n_quals; i++) {
        const gb_qual_t *q = &p->quals[f->first_qual + i];
        if (q->has_value && gb_span_eq(p, q->key, key)) {
            *out = q->value;
            return 1;
        }
    }
    return 0;
}

static gb_feed_t gb_feature_fail(gb_parser_t *p, const gb_feature_t *f, gb_status_t code, const char *what,
                                 const char *detail) {
    gb_fail(p, code, "feature %s at line %ld: %s '%.80s'", gb_str(p, f->key), f->line, what, detail);
    p->err.line = f->line;
    return GB_FEED_ERROR;
}

static void gb_resolve_identity(gb_parser_t *p, gb_feature_t *f) {
    gb_span_t gene = {0, 0}, product = {0, 0}, label = {0, 0};
    int has_gene = gb_qual_value(p, f, "gene", &gene);
    int has_product = gb_qual_value(p, f, "product", &product);
    int has_label = gb_qual_value(p, f, "label", &label);
    f->is_gene = (uint8_t)gb_span_eq(p, f->key, "gene");
    if (!gb_qual_value(p, f, "locus_tag", &f->locus) && has_gene) f->locus = gene;
    f->name = has_gene ? gene : has_product ? product : has_label ? label : f->locus.len ? f->locus : f->key;
}

/* Link every non-gene feature to the first gene sharing its locus, in any file order. */
static gb_feed_t gb_resolve_parents(gb_parser_t *p) {
    void *genes = khash_str2int_init();
    if (!genes) return gb_oom(p);
    for (size_t i = 0; i < p->n_feats; i++) {
        const gb_feature_t *f = &p->feats[i];
        if (!f->is_gene || !f->locus.len) continue;
        const char *key = gb_str(p, f->locus);
        if (khash_str2int_has_key(genes, key)) continue;
        if (khash_str2int_set(genes, key, (int)i) < 0) {
            khash_str2int_destroy(genes);
            return gb_oom(p);
        }
    }
    for (size_t i = 0; i < p->n_feats; i++) {
        gb_feature_t *f = &p->feats[i];
        if (!f->is_gene && f->locus.len && khash_str2int_has_key(genes, gb_str(p, f->locus))) f->parent_locus = f->locus;
    }
    khash_str2int_destroy(genes);
    return GB_FEED_CONTINUE;
}

static int gb_is_reserved_key(const gb_parser_t *p, gb_span_t key) {
    return gb_span_eq(p, key, "ID") || gb_span_eq(p, key, "Name") || gb_span_eq(p, key, "Parent");
}

/* Group a feature's qualifiers by key, first-occurrence order, values contiguous. */
static gb_feed_t gb_resolve_attrs(gb_parser_t *p, gb_feature_t *f) {
    f->first_attr = p->n_attrs;
    for (size_t i = 0; i < f->n_quals; i++) {
        gb_span_t key = p->quals[f->first_qual + i].key;
        if (gb_is_reserved_key(p, key))
            return gb_feature_fail(p, f, GB_ERR_UNSUPPORTED, "qualifier collides with a synthesized GFF3 key",
                                   gb_str(p, key));
        size_t g = f->first_attr;
        while (g < p->n_attrs && !gb_span_eq_span(p, p->attrs[g].key, key)) g++;
        if (g < p->n_attrs) continue;
        gb_attr_t a = {key, 0, 0};
        if (!GB_PUSH(p, attrs, n_attrs, cap_attrs, a)) return gb_oom(p);
    }
    f->n_attrs = p->n_attrs - f->first_attr;
    for (size_t g = f->first_attr; g < p->n_attrs; g++) {
        size_t first = p->n_attr_values;
        for (size_t i = 0; i < f->n_quals; i++) {
            size_t qi = f->first_qual + i;
            if (!gb_span_eq_span(p, p->quals[qi].key, p->attrs[g].key)) continue;
            if (!GB_PUSH(p, attr_values, n_attr_values, cap_attr_values, qi)) return gb_oom(p);
        }
        p->attrs[g].first_value = first;
        p->attrs[g].n_values = p->n_attr_values - first;
    }
    return GB_FEED_CONTINUE;
}

/* Emit segments in biological order and carry the CDS phase across them. */
static gb_feed_t gb_resolve_rows(gb_parser_t *p, size_t index) {
    gb_feature_t *f = &p->feats[index];
    int is_cds = gb_span_eq(p, f->key, "CDS");
    int phase = -1;
    if (is_cds) {
        gb_span_t cs;
        phase = 0;
        if (gb_qual_value(p, f, "codon_start", &cs)) {
            const char *v = gb_str(p, cs);
            if (cs.len != 1 || v[0] < '1' || v[0] > '3')
                return gb_feature_fail(p, f, GB_ERR_SYNTAX, "invalid /codon_start", v);
            phase = v[0] - '1';
        }
    }
    f->first_row = p->n_rows;
    for (size_t k = 0; k < f->n_segs; k++) {
        size_t si = f->outer_complement ? f->first_seg + f->n_segs - 1 - k : f->first_seg + k;
        const gb_seg_t *s = &p->segs[si];
        gb_row_t r;
        r.feature = index;
        r.start = s->start;
        r.end = s->end;
        r.remote = s->remote;
        r.strand = (f->outer_complement ^ s->complement) ? '-' : '+';
        r.phase = (int8_t)phase;
        if (!GB_PUSH(p, rows, n_rows, cap_rows, r)) return gb_oom(p);
        if (is_cds) {
            int64_t carry = ((s->end - s->start + 1 - phase) % 3 + 3) % 3;
            phase = (int)((3 - carry) % 3);
        }
    }
    f->n_rows = p->n_rows - f->first_row;
    return GB_FEED_CONTINUE;
}

static gb_feed_t gb_resolve_record(gb_parser_t *p) {
    for (size_t i = 0; i < p->n_feats; i++) gb_resolve_identity(p, &p->feats[i]);
    gb_feed_t r = gb_resolve_parents(p);
    if (r == GB_FEED_ERROR) return r;
    for (size_t i = 0; i < p->n_feats; i++) {
        r = gb_resolve_attrs(p, &p->feats[i]);
        if (r == GB_FEED_ERROR) return r;
        r = gb_resolve_rows(p, i);
        if (r == GB_FEED_ERROR) return r;
    }
    return GB_FEED_CONTINUE;
}

/* ------------------------------------------------------------------------ */
/* Header and sequence                                                       */
/* ------------------------------------------------------------------------ */

static gb_feed_t gb_first_token_field(gb_parser_t *p, const char *s, size_t n, size_t after, gb_span_t *out) {
    size_t start = gb_skip_spaces(s, n, after);
    size_t end = gb_token_end(s, n, start);
    if (!gb_arena_put(p, s + start, end - start, out)) return gb_oom(p);
    return GB_FEED_CONTINUE;
}

static gb_feed_t gb_finish_record(gb_parser_t *p) {
    if (p->section == GB_SECTION_FEATURES)
        return gb_fail(p, GB_ERR_SYNTAX, "record %s: FEATURES table is not followed by a sequence section",
                       gb_str(p, p->locus));
    gb_feed_t closed = gb_close_feature(p);
    if (closed == GB_FEED_ERROR) return closed;
    if (p->mode == GB_MODE_FEATURES) {
        gb_feed_t resolved = gb_resolve_record(p);
        if (resolved == GB_FEED_ERROR) return resolved;
    }
    p->n_records++;
    p->section = GB_SECTION_BETWEEN;
    return GB_FEED_RECORD;
}

/* LOCUS name, then the length token that precedes "bp" or "aa", then the
 * topology token; other LOCUS fields are not needed. */
static gb_feed_t gb_begin_record(gb_parser_t *p, const char *s, size_t n) {
    gb_reset_record(p);
    p->section = GB_SECTION_HEADER;
    gb_feed_t r = gb_first_token_field(p, s, n, 5, &p->locus);
    if (r == GB_FEED_ERROR) return r;
    size_t i = gb_token_end(s, n, gb_skip_spaces(s, n, 5)); /* past the name */
    size_t prev_start = 0, prev_end = 0;
    for (;;) {
        size_t start = gb_skip_spaces(s, n, i);
        if (start >= n) break;
        size_t end = gb_token_end(s, n, start);
        size_t len = end - start;
        if ((len == 2 && (memcmp(s + start, "bp", 2) == 0 || memcmp(s + start, "aa", 2) == 0)) && prev_end > prev_start) {
            int64_t v = 0;
            size_t k = prev_start;
            while (k < prev_end && s[k] >= '0' && s[k] <= '9' && v < INT64_MAX / 10) v = v * 10 + (s[k++] - '0');
            if (k == prev_end) p->seq_length = v;
        } else if (len == 8 && memcmp(s + start, "circular", 8) == 0) {
            p->topology = GB_TOPOLOGY_CIRCULAR;
        } else if (len == 6 && memcmp(s + start, "linear", 6) == 0) {
            p->topology = GB_TOPOLOGY_LINEAR;
        }
        prev_start = start;
        prev_end = end;
        i = end;
    }
    return GB_FEED_CONTINUE;
}

/* Residues are written through a local pointer with ASCII range checks: a
 * per-character kputc would store a NUL terminator after every base, and the
 * locale-aware ctype calls cost a PLT call each, which together dominated the
 * FASTA converter. The line bounds the output, so one reservation suffices. */
static gb_feed_t gb_sequence_line(gb_parser_t *p, const char *s, size_t n) {
    if (p->mode != GB_MODE_SEQUENCE) return GB_FEED_CONTINUE;
    if (ks_resize(&p->residues, n + 1) < 0) return gb_oom(p);
    char *out = p->residues.s;
    size_t m = 0;
    for (size_t i = 0; i < n; i++) {
        char c = s[i];
        if (c == ' ' || c == '\t' || (c >= '0' && c <= '9')) continue;
        if (c >= 'a' && c <= 'z') c = (char)(c - ('a' - 'A'));
        if ((c >= 'A' && c <= 'Z') || c == '-' || c == '*') {
            out[m++] = c;
            continue;
        }
        p->residues.l = 0;
        return gb_fail(p, GB_ERR_SYNTAX, "unexpected character '%c' in sequence data", c);
    }
    out[m] = '\0';
    p->residues.l = m;
    return m ? GB_FEED_RESIDUES : GB_FEED_CONTINUE;
}

static gb_feed_t gb_keyword_line(gb_parser_t *p, const char *s, size_t n) {
    size_t kw = gb_token_end(s, n, 0);
    int is_locus = kw == 5 && memcmp(s, "LOCUS", 5) == 0;
    int is_end = n >= 2 && s[0] == '/' && s[1] == '/';

    if (p->section == GB_SECTION_BETWEEN) {
        if (is_locus) return gb_begin_record(p, s, n);
        return GB_FEED_CONTINUE; /* release-file preamble or other text between records */
    }
    if (is_locus) return gb_fail(p, GB_ERR_SYNTAX, "LOCUS inside record %s", gb_str(p, p->locus));
    if (is_end) return gb_finish_record(p);

    if (p->section == GB_SECTION_SEQUENCE)
        return gb_fail(p, GB_ERR_SYNTAX, "keyword %.*s after ORIGIN", (int)kw, s);

    if (p->section == GB_SECTION_FEATURES) {
        gb_feed_t closed = gb_close_feature(p);
        if (closed == GB_FEED_ERROR) return closed;
    }
    p->header_field = GB_HEADER_NONE;

    if (kw == 6 && memcmp(s, "ORIGIN", 6) == 0) {
        p->section = GB_SECTION_SEQUENCE;
        return GB_FEED_ORIGIN;
    }
    if (kw == 8 && memcmp(s, "FEATURES", 8) == 0) {
        if (p->section != GB_SECTION_HEADER) return gb_fail(p, GB_ERR_SYNTAX, "second FEATURES table");
        p->section = GB_SECTION_FEATURES;
        return GB_FEED_CONTINUE;
    }
    if (p->section == GB_SECTION_FEATURES) {
        int sequence_header =
            (kw == 6 && memcmp(s, "CONTIG", 6) == 0) ||
            (kw == 3 && (memcmp(s, "WGS", 3) == 0 || memcmp(s, "TSA", 3) == 0 || memcmp(s, "TLS", 3) == 0)) ||
            (kw == 4 && n >= 10 && memcmp(s, "BASE COUNT", 10) == 0 &&
             (n == 10 || s[10] == ' ' || s[10] == '\t'));
        if (!sequence_header)
            return gb_fail(p, GB_ERR_SYNTAX, "record %s: FEATURES table is not followed by a sequence section",
                           gb_str(p, p->locus));
        p->section = GB_SECTION_OTHER;
    }
    if (p->section != GB_SECTION_HEADER) return GB_FEED_CONTINUE;

    if (kw == 10 && memcmp(s, "DEFINITION", 10) == 0) {
        size_t start = gb_skip_spaces(s, n, kw);
        size_t end = gb_rtrim(s, n);
        if (!gb_arena_put(p, s + start, end > start ? end - start : 0, &p->definition)) return gb_oom(p);
        p->header_field = GB_HEADER_DEFINITION;
        return GB_FEED_CONTINUE;
    }
    if (kw == 9 && memcmp(s, "ACCESSION", 9) == 0) return gb_first_token_field(p, s, n, kw, &p->accession);
    if (kw == 7 && memcmp(s, "VERSION", 7) == 0) return gb_first_token_field(p, s, n, kw, &p->version);
    return GB_FEED_CONTINUE;
}

static gb_feed_t gb_indented_line(gb_parser_t *p, const char *s, size_t n) {
    switch (p->section) {
    case GB_SECTION_HEADER:
        if (p->header_field == GB_HEADER_DEFINITION) {
            size_t start = gb_skip_spaces(s, n, 0);
            size_t end = gb_rtrim(s, n);
            if (start >= end) return GB_FEED_CONTINUE;
            if (!gb_arena_extend(p, &p->definition, " ", 1) ||
                !gb_arena_extend(p, &p->definition, s + start, end - start))
                return gb_oom(p);
        }
        return GB_FEED_CONTINUE;
    case GB_SECTION_FEATURES:
        return gb_features_indented(p, s, n);
    case GB_SECTION_SEQUENCE:
        return gb_sequence_line(p, s, n);
    default:
        return GB_FEED_CONTINUE;
    }
}

gb_feed_t gb_parser_feed(gb_parser_t *p, const char *line, size_t len) {
    if (p->section == GB_SECTION_FAILED) return GB_FEED_ERROR;
    p->line_no++;
    if (len && line[len - 1] == '\r') len--;
    if (len == 0) return GB_FEED_CONTINUE;
    if (line[0] != ' ' && line[0] != '\t') return gb_keyword_line(p, line, len);
    return gb_indented_line(p, line, len);
}

gb_status_t gb_parser_finish(gb_parser_t *p) {
    if (p->section == GB_SECTION_FAILED) return p->err.code;
    if (p->section != GB_SECTION_BETWEEN) {
        p->err.code = GB_ERR_TRUNCATED;
        p->err.line = p->line_no;
        snprintf(p->err.msg, sizeof(p->err.msg), "record %s is not terminated by //", gb_str(p, p->locus));
        p->section = GB_SECTION_FAILED;
        return GB_ERR_TRUNCATED;
    }
    if (p->n_records == 0) {
        p->err.code = GB_ERR_EMPTY;
        p->err.line = 0;
        snprintf(p->err.msg, sizeof(p->err.msg), "no LOCUS record found");
        p->section = GB_SECTION_FAILED;
        return GB_ERR_EMPTY;
    }
    return GB_OK;
}

/* ------------------------------------------------------------------------ */
/* GFF3 attributes                                                           */
/* ------------------------------------------------------------------------ */

/* GFF3 reserves ; = & , and %; control characters are encoded as well. */
static int gb_put_encoded(kstring_t *k, const char *s, size_t n) {
    for (size_t i = 0; i < n; i++) {
        unsigned char ch = (unsigned char)s[i];
        if (ch == ';' || ch == '=' || ch == '&' || ch == ',' || ch == '%' || ch < 0x20 || ch == 0x7f) {
            char hex[4];
            snprintf(hex, sizeof(hex), "%%%02X", ch);
            if (kputsn(hex, 3, k) == EOF) return 0;
        } else if (kputc(ch, k) == EOF) {
            return 0;
        }
    }
    return 1;
}

static size_t gb_synth_count(const gb_feature_t *f) {
    return 2 + (f->parent_locus.len ? 1 : 0);
}

size_t gb_feature_attr_count(const gb_parser_t *p, size_t feature) {
    const gb_feature_t *f = &p->feats[feature];
    return gb_synth_count(f) + f->n_attrs;
}

static gb_status_t gb_emit_key(const gb_parser_t *p, const gb_feature_t *f, size_t index, kstring_t *k) {
    size_t synth = gb_synth_count(f);
    const char *lit = index == 0 ? "ID" : index == 1 ? "Name" : index < synth ? "Parent" : NULL;
    if (lit) return kputs(lit, k) == EOF ? GB_ERR_NOMEM : GB_OK;
    if (index - synth >= f->n_attrs) return GB_ERR_UNSUPPORTED;
    gb_span_t key = p->attrs[f->first_attr + index - synth].key;
    return gb_put_encoded(k, gb_str(p, key), key.len) ? GB_OK : GB_ERR_NOMEM;
}

static gb_status_t gb_emit_value(const gb_parser_t *p, const gb_feature_t *f, size_t index, kstring_t *k) {
    size_t synth = gb_synth_count(f);
    if (index == 0) {
        if (f->is_gene && f->locus.len) {
            if (kputs("gene-", k) == EOF || !gb_put_encoded(k, gb_str(p, f->locus), f->locus.len)) return GB_ERR_NOMEM;
            return GB_OK;
        }
        char ordinal[32];
        snprintf(ordinal, sizeof(ordinal), "-%ld", f->ordinal);
        if (!gb_put_encoded(k, gb_str(p, f->key), f->key.len) || kputs(ordinal, k) == EOF) return GB_ERR_NOMEM;
        return GB_OK;
    }
    if (index == 1) return gb_put_encoded(k, gb_str(p, f->name), f->name.len) ? GB_OK : GB_ERR_NOMEM;
    if (index < synth) {
        if (kputs("gene-", k) == EOF || !gb_put_encoded(k, gb_str(p, f->parent_locus), f->parent_locus.len))
            return GB_ERR_NOMEM;
        return GB_OK;
    }
    if (index - synth >= f->n_attrs) return GB_ERR_UNSUPPORTED;
    const gb_attr_t *a = &p->attrs[f->first_attr + index - synth];
    for (size_t v = 0; v < a->n_values; v++) {
        const gb_qual_t *q = &p->quals[p->attr_values[a->first_value + v]];
        if (v && kputc(',', k) == EOF) return GB_ERR_NOMEM;
        if (!q->has_value) {
            if (kputs("true", k) == EOF) return GB_ERR_NOMEM;
        } else if (!gb_put_encoded(k, gb_str(p, q->value), q->value.len)) {
            return GB_ERR_NOMEM;
        }
    }
    return GB_OK;
}

gb_status_t gb_feature_attr_key(const gb_parser_t *p, size_t feature, size_t index, kstring_t *key) {
    key->l = 0;
    if (key->s) key->s[0] = '\0';
    return gb_emit_key(p, &p->feats[feature], index, key);
}

gb_status_t gb_feature_attr_value(const gb_parser_t *p, size_t feature, size_t index, kstring_t *value) {
    value->l = 0;
    if (value->s) value->s[0] = '\0';
    return gb_emit_value(p, &p->feats[feature], index, value);
}

gb_status_t gb_feature_attr_at(const gb_parser_t *p, size_t feature, size_t index, kstring_t *key,
                               kstring_t *value) {
    gb_status_t r = gb_feature_attr_key(p, feature, index, key);
    if (r != GB_OK) return r;
    return gb_feature_attr_value(p, feature, index, value);
}

gb_status_t gb_feature_attributes(const gb_parser_t *p, size_t feature, kstring_t *out) {
    const gb_feature_t *f = &p->feats[feature];
    size_t n = gb_feature_attr_count(p, feature);
    out->l = 0;
    for (size_t i = 0; i < n; i++) {
        if (i && kputc(';', out) == EOF) return GB_ERR_NOMEM;
        gb_status_t r = gb_emit_key(p, f, i, out);
        if (r != GB_OK) return r;
        if (kputc('=', out) == EOF) return GB_ERR_NOMEM;
        r = gb_emit_value(p, f, i, out);
        if (r != GB_OK) return r;
    }
    return GB_OK;
}
