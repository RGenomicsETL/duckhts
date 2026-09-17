/**
 * DuckHTS GenBank Reader
 *
 * Two table functions over GenBank flat files (single- or multi-record,
 * optionally bgzipped — I/O goes through htslib's hFILE layer):
 *
 *   read_genbank(path, attributes_map := FALSE)
 *       One row PER FEATURE, in read_gff's exact column shape, so it
 *       substitutes for read_gff anywhere without a schema change:
 *         seqname, source, feature, start BIGINT, end BIGINT, score DOUBLE,
 *         strand, frame, attributes [, attributes_map MAP<VARCHAR,VARCHAR>]
 *       The FEATURES table is parsed; the ORIGIN sequence block is skipped.
 *
 *   genbank_to_fasta(path, output_path := NULL, line_width := 70, overwrite := FALSE)
 *       Write-materializer (shaped like bgzip): emits a canonical multi-record
 *       FASTA to output_path and returns (success, output_path, records_written).
 *       The FEATURES table is skipped; only headers + ORIGIN are read.
 *
 * GenBank -> GFF mapping (read_genbank):
 *   seqname  = VERSION (accession.version) || ACCESSION || LOCUS name
 *              (identical to the genbank_to_fasta defline name, so feature
 *               coordinates land on the FASTA contig of the same name)
 *   source   = "GenBank"
 *   feature  = feature key (gene, CDS, mRNA, tRNA, ...); the whole-record
 *              "source" feature is dropped (it is record metadata, not an interval)
 *   start/end= 1-based inclusive (same as GFF, no offset); join()/order() flatten
 *              to one row per segment; complement(...) -> strand '-'
 *   frame    = /codon_start (1/2/3) mapped to GFF phase (0/1/2) for CDS, else "."
 *   attributes / attributes_map carry synthesized GFF3 keys so the interval
 *   reader's canonical projection works unchanged:
 *     ID     = "gene-<locus_tag>" for genes, "<feature>-<n>" otherwise
 *     Name   = /gene || /product || /label || /locus_tag || feature
 *     Parent = "gene-<locus_tag>" for non-gene features whose /locus_tag matches
 *              a gene seen in the same record (gene = root); absent otherwise
 *   plus the original qualifiers (except /translation, which is large and
 *   redundant with the sequence).
 */

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <ctype.h>
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>
#include <unistd.h>

#include <htslib/hts.h>
#include <htslib/kstring.h>

#define GB_BATCH_SIZE 2048

/* GTF/GFF-compatible column indices (must match read_gff). */
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
    GB_COL_ATTRIBUTES_MAP,
    GB_COL_COUNT
};

/* ================================================================
 * Small vector helpers (validity + GFF3 attribute MAP)
 * ================================================================ */

static void gb_set_null(duckdb_vector vec, idx_t row) {
    duckdb_vector_ensure_validity_writable(vec);
    uint64_t *v = duckdb_vector_get_validity(vec);
    duckdb_validity_set_row_invalid(v, row);
}

static void gb_trim_span(const char **start, int *len) {
    const char *s = *start;
    int l = *len;
    while (l > 0 && (*s == ' ' || *s == '\t')) { s++; l--; }
    while (l > 0 && (s[l - 1] == ' ' || s[l - 1] == '\t')) l--;
    *start = s;
    *len = l;
}

/* Count GFF3 key=value pairs (semicolon-separated). Mirrors tabix_reader.c. */
static int gb_count_pairs(const char *s) {
    int count = 0;
    const char *p = s;
    while (*p) {
        while (*p == ';' || *p == ' ' || *p == '\t') p++;
        if (!*p) break;
        const char *key = p;
        while (*p && *p != '=' && *p != ';') p++;
        if (*p != '=') { while (*p && *p != ';') p++; continue; }
        int key_len = (int)(p - key);
        p++;
        while (*p && *p != ';') p++;
        gb_trim_span(&key, &key_len);
        if (key_len > 0) count++;
        if (*p == ';') p++;
    }
    return count;
}

/* Build a MAP<VARCHAR,VARCHAR> entry for `row` from a GFF3 attribute string. */
static void gb_fill_attr_map(duckdb_vector vec, idx_t row, const char *s) {
    if (!s || s[0] == '\0') {
        duckdb_vector_ensure_validity_writable(vec);
        duckdb_validity_set_row_invalid(duckdb_vector_get_validity(vec), row);
        duckdb_list_entry empty = {duckdb_list_vector_get_size(vec), 0};
        ((duckdb_list_entry *)duckdb_vector_get_data(vec))[row] = empty;
        return;
    }

    int pair_count = gb_count_pairs(s);
    duckdb_list_entry entry;
    entry.offset = duckdb_list_vector_get_size(vec);
    entry.length = (idx_t)pair_count;

    if (pair_count == 0) {
        ((duckdb_list_entry *)duckdb_vector_get_data(vec))[row] = entry;
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
        while (*p && *p != '=' && *p != ';') p++;
        if (*p != '=') { while (*p && *p != ';') p++; continue; }
        int key_len = (int)(p - key);
        p++;
        const char *val = p;
        while (*p && *p != ';') p++;
        int val_len = (int)(p - val);
        gb_trim_span(&key, &key_len);
        gb_trim_span(&val, &val_len);
        if (key_len > 0) {
            duckdb_vector_assign_string_element_len(key_vec, entry.offset + write_idx, key, key_len);
            duckdb_vector_assign_string_element_len(val_vec, entry.offset + write_idx, val, val_len);
            write_idx++;
        }
        if (*p == ';') p++;
    }
    entry.length = (idx_t)write_idx;
    ((duckdb_list_entry *)duckdb_vector_get_data(vec))[row] = entry;
}

/* Append `key=<percent-encoded val>;` to a GFF3 attribute kstring.
 * Reserved chars (;=,&%, control) are percent-encoded per GFF3. */
static void gb_attr_append(kstring_t *k, const char *key, const char *val) {
    if (!key || key[0] == '\0' || !val) return;
    kputs(key, k);
    kputc('=', k);
    for (const char *p = val; *p; p++) {
        unsigned char c = (unsigned char)*p;
        if (c == ';' || c == '=' || c == '&' || c == ',' || c == '%' ||
            c < 0x20 || c == 0x7f) {
            char buf[4];
            snprintf(buf, sizeof(buf), "%%%02X", c);
            kputsn(buf, 3, k);
        } else {
            kputc((char)c, k);
        }
    }
    kputc(';', k);
}

/* ================================================================
 * Parsed-feature emit rows (read_genbank)
 * ================================================================ */

typedef struct {
    char   *seqname;     /* owned */
    char   *feature;     /* owned */
    int64_t start;
    int64_t end;
    char    strand;      /* '+', '-', '.' */
    char    frame;       /* '0'..'2' or '.' */
    char   *attr;        /* owned GFF3 attribute string (never NULL) */
} gb_row_t;

typedef struct {
    gb_row_t *rows;
    size_t    n;
    size_t    cap;
} gb_rows_t;

static int gb_rows_push(gb_rows_t *v, gb_row_t r) {
    if (v->n == v->cap) {
        size_t ncap = v->cap ? v->cap * 2 : 256;
        gb_row_t *nr = (gb_row_t *)realloc(v->rows, ncap * sizeof(gb_row_t));
        if (!nr) return 0;
        v->rows = nr;
        v->cap = ncap;
    }
    v->rows[v->n++] = r;
    return 1;
}

static void gb_rows_free(gb_rows_t *v) {
    if (!v->rows) return;
    for (size_t i = 0; i < v->n; i++) {
        free(v->rows[i].seqname);
        free(v->rows[i].feature);
        free(v->rows[i].attr);
    }
    free(v->rows);
    v->rows = NULL;
    v->n = v->cap = 0;
}

/* A single qualifier (/key=value). */
typedef struct { char *key; char *val; } gb_qual_t;

/* The feature currently being accumulated across lines. */
typedef struct {
    char      key[64];        /* feature key, e.g. "CDS" */
    kstring_t loc;            /* raw location string (may be multi-line) */
    gb_qual_t quals[256];
    int       n_quals;
    kstring_t qbuf;           /* current qualifier value being accumulated */
    char      qkey[64];
    int       in_quoted;      /* inside an unterminated quoted value */
    int       have_qual;      /* a qualifier is being accumulated */
} gb_feature_t;

/* Per-record set of gene /locus_tag values (for Parent linkage). */
typedef struct { char **tags; size_t n, cap; } gb_strset_t;

static int gb_strset_has(gb_strset_t *s, const char *t) {
    if (!t) return 0;
    for (size_t i = 0; i < s->n; i++)
        if (strcmp(s->tags[i], t) == 0) return 1;
    return 0;
}
static void gb_strset_add(gb_strset_t *s, const char *t) {
    if (!t || gb_strset_has(s, t)) return;
    if (s->n == s->cap) {
        size_t ncap = s->cap ? s->cap * 2 : 64;
        char **nt = (char **)realloc(s->tags, ncap * sizeof(char *));
        if (!nt) return;
        s->tags = nt; s->cap = ncap;
    }
    s->tags[s->n++] = strdup(t);
}
static void gb_strset_clear(gb_strset_t *s) {
    for (size_t i = 0; i < s->n; i++) free(s->tags[i]);
    free(s->tags);
    s->tags = NULL; s->n = s->cap = 0;
}

static const char *gb_qual_get(gb_feature_t *f, const char *key) {
    for (int i = 0; i < f->n_quals; i++)
        if (strcmp(f->quals[i].key, key) == 0) return f->quals[i].val;
    return NULL;
}

/* ----- location parsing ----------------------------------------------- */

typedef struct { char *seqname; int64_t start, end; char strand; } gb_seg_t;
typedef struct { gb_seg_t *segs; int n, cap; } gb_seglist_t;

static void gb_seglist_push(gb_seglist_t *l, gb_seg_t s) {
    if (l->n == l->cap) {
        int ncap = l->cap ? l->cap * 2 : 8;
        gb_seg_t *ns = (gb_seg_t *)realloc(l->segs, (size_t)ncap * sizeof(gb_seg_t));
        if (!ns) return;
        l->segs = ns; l->cap = ncap;
    }
    l->segs[l->n++] = s;
}

static char gb_flip(char s) { return s == '-' ? '+' : '-'; }

/* Mutually recursive with gb_parse_location (join/order + complement nesting). */
static void gb_parse_location_n(const char *s, int len, char strand, gb_seglist_t *out);
static void gb_parse_location(const char *s, char strand, gb_seglist_t *out);

/* Parse a single span like "<123..>456", "123", "12^13", or "ACC.1:1..9". */
static void gb_parse_span(const char *s, int len, char strand, gb_seglist_t *out) {
    /* copy + trim into a small buffer */
    char buf[256];
    if (len <= 0) return;
    if (len >= (int)sizeof(buf)) len = (int)sizeof(buf) - 1;
    memcpy(buf, s, (size_t)len);
    buf[len] = '\0';

    char *loc = buf;
    char *seqname_override = NULL;
    /* remote reference: ACCESSION:location */
    char *colon = strchr(loc, ':');
    if (colon) { *colon = '\0'; seqname_override = loc; loc = colon + 1; }

    /* strip partial markers */
    char clean[256];
    int ci = 0;
    for (char *p = loc; *p && ci < (int)sizeof(clean) - 1; p++) {
        if (*p == '<' || *p == '>') continue;
        clean[ci++] = *p;
    }
    clean[ci] = '\0';

    int64_t a = 0, b = 0;
    char *dots = strstr(clean, "..");
    if (dots) {
        *dots = '\0';
        a = strtoll(clean, NULL, 10);
        b = strtoll(dots + 2, NULL, 10);
    } else {
        char *caret = strchr(clean, '^');
        if (caret) { *caret = '\0'; a = strtoll(clean, NULL, 10); b = strtoll(caret + 1, NULL, 10); }
        else { a = b = strtoll(clean, NULL, 10); }
    }
    if (a == 0 && b == 0) return; /* unparseable */
    if (b < a) { int64_t t = a; a = b; b = t; }

    gb_seg_t seg;
    seg.seqname = seqname_override ? strdup(seqname_override) : NULL;
    seg.start = a;
    seg.end = b;
    seg.strand = strand;
    gb_seglist_push(out, seg);
}

/* Recursively parse a location, splitting join()/order() on top-level commas
 * and flipping strand through complement(). */
static void gb_parse_location(const char *s, char strand, gb_seglist_t *out) {
    const char *start = s;
    int len = (int)strlen(s);
    gb_trim_span(&start, &len);
    if (len <= 0) return;

    if (len > 11 && strncmp(start, "complement(", 11) == 0 && start[len - 1] == ')') {
        gb_parse_location_n(start + 11, len - 12, gb_flip(strand), out);
        return;
    }
    if ((len > 5 && strncmp(start, "join(", 5) == 0 && start[len - 1] == ')') ||
        (len > 6 && strncmp(start, "order(", 6) == 0 && start[len - 1] == ')')) {
        int off = (start[0] == 'j') ? 5 : 6;
        const char *inner = start + off;
        int inner_len = len - off - 1;
        /* split inner on top-level commas */
        int depth = 0, seg_start = 0;
        for (int i = 0; i <= inner_len; i++) {
            char c = (i < inner_len) ? inner[i] : ',';
            if (c == '(') depth++;
            else if (c == ')') depth--;
            else if (c == ',' && depth == 0) {
                gb_parse_location_n(inner + seg_start, i - seg_start, strand, out);
                seg_start = i + 1;
            }
        }
        return;
    }
    gb_parse_span(start, len, strand, out);
}

/* length-bounded entry point (operates on a [s, s+len) slice). */
static void gb_parse_location_n(const char *s, int len, char strand, gb_seglist_t *out) {
    char tmp[1024];
    if (len < 0) len = 0;
    if (len >= (int)sizeof(tmp)) len = (int)sizeof(tmp) - 1;
    memcpy(tmp, s, (size_t)len);
    tmp[len] = '\0';
    gb_parse_location(tmp, strand, out);
}

/* ----- feature finalize ----------------------------------------------- */

/* Flush the current qualifier value (qbuf/qkey) into the feature's qual list. */
static void gb_flush_qual(gb_feature_t *f) {
    if (!f->have_qual) return;
    if (f->n_quals < (int)(sizeof(f->quals) / sizeof(f->quals[0]))) {
        f->quals[f->n_quals].key = strdup(f->qkey);
        f->quals[f->n_quals].val = strdup(f->qbuf.l ? f->qbuf.s : "");
        f->n_quals++;
    }
    f->have_qual = 0;
    f->in_quoted = 0;
    f->qbuf.l = 0;
    if (f->qbuf.s) f->qbuf.s[0] = '\0';
    f->qkey[0] = '\0';
}

static void gb_feature_reset(gb_feature_t *f) {
    for (int i = 0; i < f->n_quals; i++) { free(f->quals[i].key); free(f->quals[i].val); }
    f->n_quals = 0;
    f->key[0] = '\0';
    f->loc.l = 0;
    if (f->loc.s) f->loc.s[0] = '\0';
    f->have_qual = 0;
    f->in_quoted = 0;
    f->qbuf.l = 0;
    if (f->qbuf.s) f->qbuf.s[0] = '\0';
    f->qkey[0] = '\0';
}

/* Emit one row per location segment of the finalized feature. */
static void gb_emit_feature(gb_feature_t *f, const char *contig, gb_rows_t *rows,
                            gb_strset_t *gene_tags, long *feat_counter) {
    if (f->key[0] == '\0') return;
    if (f->have_qual) gb_flush_qual(f);
    if (strcmp(f->key, "source") == 0) return; /* whole-record metadata */

    gb_seglist_t segs = {0};
    if (f->loc.l) gb_parse_location(f->loc.s, '+', &segs);
    if (segs.n == 0) { free(segs.segs); return; }

    int is_gene = (strcmp(f->key, "gene") == 0);
    const char *locus = gb_qual_get(f, "locus_tag");
    if (!locus) locus = gb_qual_get(f, "gene");
    const char *gene_name = gb_qual_get(f, "gene");
    const char *product = gb_qual_get(f, "product");
    const char *label = gb_qual_get(f, "label");
    const char *codon_start = gb_qual_get(f, "codon_start");

    if (is_gene && locus) gb_strset_add(gene_tags, locus);

    long this_id = (*feat_counter)++;

    /* Build the shared GFF3 attribute string once for all segments. */
    kstring_t attr = {0, 0, NULL};
    char idbuf[256];
    if (is_gene && locus) snprintf(idbuf, sizeof(idbuf), "gene-%s", locus);
    else snprintf(idbuf, sizeof(idbuf), "%s-%ld", f->key, this_id);
    gb_attr_append(&attr, "ID", idbuf);

    const char *name = gene_name ? gene_name : (product ? product : (label ? label : (locus ? locus : f->key)));
    gb_attr_append(&attr, "Name", name);

    if (!is_gene && locus && gb_strset_has(gene_tags, locus)) {
        char parent[256];
        snprintf(parent, sizeof(parent), "gene-%s", locus);
        gb_attr_append(&attr, "Parent", parent);
    }

    for (int i = 0; i < f->n_quals; i++) {
        if (strcmp(f->quals[i].key, "translation") == 0) continue; /* large + redundant */
        gb_attr_append(&attr, f->quals[i].key, f->quals[i].val);
    }

    char frame = '.';
    if (strcmp(f->key, "CDS") == 0 && codon_start && codon_start[0] >= '1' && codon_start[0] <= '3')
        frame = (char)('0' + (codon_start[0] - '1'));

    for (int i = 0; i < segs.n; i++) {
        gb_row_t r;
        r.seqname = strdup(segs.segs[i].seqname ? segs.segs[i].seqname : contig);
        r.feature = strdup(f->key);
        r.start = segs.segs[i].start;
        r.end = segs.segs[i].end;
        r.strand = segs.segs[i].strand;
        r.frame = frame;
        r.attr = strdup(attr.l ? attr.s : "");
        if (!gb_rows_push(rows, r)) { free(r.seqname); free(r.feature); free(r.attr); }
    }

    for (int i = 0; i < segs.n; i++) free(segs.segs[i].seqname);
    free(segs.segs);
    free(attr.s);
}

/* ----- record header field helpers ------------------------------------ */

/* First whitespace-delimited token after a fixed keyword column. */
static void gb_first_token(const char *line, char *out, size_t outsz) {
    out[0] = '\0';
    const char *p = line;
    while (*p && *p != ' ' && *p != '\t') p++;      /* skip keyword */
    while (*p == ' ' || *p == '\t') p++;            /* skip ws */
    size_t i = 0;
    while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r' && i + 1 < outsz)
        out[i++] = *p++;
    out[i] = '\0';
}

/* ================================================================
 * read_genbank: parse the whole file's FEATURES into emit rows
 * ================================================================ */

static int gb_parse_features(htsFile *fp, gb_rows_t *rows, char *errbuf, size_t errsz) {
    kstring_t line = {0, 0, NULL};
    gb_feature_t feat = {0};
    gb_strset_t gene_tags = {0};
    char locus[256] = "", accession[256] = "", version[256] = "", contig[256] = "";
    long feat_counter = 0;
    int in_features = 0, in_feature = 0;

    while (hts_getline(fp, '\n', &line) >= 0) {
        if (line.l == 0) continue;
        const char *s = line.s;

        /* Record boundary */
        if (s[0] == '/' && s[1] == '/') {
            if (in_feature) { gb_emit_feature(&feat, contig, rows, &gene_tags, &feat_counter); gb_feature_reset(&feat); in_feature = 0; }
            in_features = 0;
            gb_strset_clear(&gene_tags);
            locus[0] = accession[0] = version[0] = contig[0] = '\0';
            continue;
        }

        /* Top-level keyword (column 0 non-space) ends the FEATURES table. */
        if (s[0] != ' ' && s[0] != '\t') {
            if (in_feature) { gb_emit_feature(&feat, contig, rows, &gene_tags, &feat_counter); gb_feature_reset(&feat); in_feature = 0; }
            if (strncmp(s, "LOCUS", 5) == 0) { gb_first_token(s, locus, sizeof(locus)); }
            else if (strncmp(s, "ACCESSION", 9) == 0) { gb_first_token(s, accession, sizeof(accession)); }
            else if (strncmp(s, "VERSION", 7) == 0) { gb_first_token(s, version, sizeof(version)); }
            else if (strncmp(s, "FEATURES", 8) == 0) {
                snprintf(contig, sizeof(contig), "%s",
                         version[0] ? version : (accession[0] ? accession : locus));
                in_features = 1;
            } else {
                in_features = 0; /* ORIGIN, BASE COUNT, CONTIG, etc. */
            }
            continue;
        }

        if (!in_features) continue; /* indented line outside FEATURES (e.g. DEFINITION cont.) */

        /* Within FEATURES: classify by indentation. */
        int i = 0;
        while (s[i] == ' ' || s[i] == '\t') i++;
        if (s[i] == '\0') continue;

        if (s[i] == '/' && i >= 16) {
            /* qualifier line: /key or /key=value */
            if (in_feature) {
                if (feat.in_quoted) {
                    /* continuation of an open quoted value (shouldn't start with '/', but guard) */
                    kputc(' ', &feat.qbuf);
                    kputs(s + i, &feat.qbuf);
                } else {
                    gb_flush_qual(&feat);
                    const char *q = s + i + 1; /* after '/' */
                    const char *eq = strchr(q, '=');
                    feat.have_qual = 1;
                    feat.qbuf.l = 0; if (feat.qbuf.s) feat.qbuf.s[0] = '\0';
                    if (!eq) {
                        snprintf(feat.qkey, sizeof(feat.qkey), "%s", q);
                        /* strip trailing CR/whitespace */
                        for (char *k = feat.qkey; *k; k++) if (*k=='\r'||*k=='\n'||*k==' ') { *k='\0'; break; }
                    } else {
                        int klen = (int)(eq - q);
                        if (klen >= (int)sizeof(feat.qkey)) klen = (int)sizeof(feat.qkey) - 1;
                        memcpy(feat.qkey, q, (size_t)klen); feat.qkey[klen] = '\0';
                        const char *v = eq + 1;
                        if (*v == '"') {
                            v++;
                            const char *endq = strrchr(v, '"');
                            if (endq && endq > v) { kputsn(v, (int)(endq - v), &feat.qbuf); }
                            else if (endq == v) { /* empty "" */ }
                            else { kputs(v, &feat.qbuf); feat.in_quoted = 1; }
                        } else {
                            /* unquoted single-line value (number, token); strip CR */
                            kstring_t t = {0,0,NULL};
                            kputs(v, &t);
                            while (t.l && (t.s[t.l-1]=='\r'||t.s[t.l-1]=='\n'||t.s[t.l-1]==' ')) t.s[--t.l]='\0';
                            kputsn(t.s ? t.s : "", (int)t.l, &feat.qbuf);
                            free(t.s);
                        }
                    }
                }
            }
            continue;
        }

        if (i <= 8) {
            /* new feature line: key + (start of) location */
            if (in_feature) { gb_emit_feature(&feat, contig, rows, &gene_tags, &feat_counter); gb_feature_reset(&feat); }
            in_feature = 1;
            int k = 0;
            while (s[i] && s[i] != ' ' && s[i] != '\t' && k + 1 < (int)sizeof(feat.key)) feat.key[k++] = s[i++];
            feat.key[k] = '\0';
            while (s[i] == ' ' || s[i] == '\t') i++;
            /* rest of line is the (possibly partial) location */
            const char *loc = s + i;
            int llen = (int)strlen(loc);
            while (llen && (loc[llen-1]=='\r'||loc[llen-1]=='\n'||loc[llen-1]==' ')) llen--;
            kputsn(loc, llen, &feat.loc);
            continue;
        }

        /* continuation line (deeper indent, not '/') */
        if (in_feature) {
            const char *cont = s + i;
            int clen = (int)strlen(cont);
            while (clen && (cont[clen-1]=='\r'||cont[clen-1]=='\n'||cont[clen-1]==' ')) clen--;
            if (feat.have_qual && feat.in_quoted) {
                const char *endq = NULL;
                for (int j = clen - 1; j >= 0; j--) if (cont[j]=='"') { endq = cont + j; break; }
                kputc(' ', &feat.qbuf);
                if (endq) { kputsn(cont, (int)(endq - cont), &feat.qbuf); feat.in_quoted = 0; }
                else kputsn(cont, clen, &feat.qbuf);
            } else if (!feat.have_qual) {
                /* still accumulating the location */
                kputsn(cont, clen, &feat.loc);
            }
        }
    }

    if (in_feature) { gb_emit_feature(&feat, contig, rows, &gene_tags, &feat_counter); }
    gb_feature_reset(&feat);
    gb_strset_clear(&gene_tags);
    free(feat.loc.s);
    free(feat.qbuf.s);
    free(line.s);
    (void)errbuf; (void)errsz;
    return 1;
}

/* ----- bind / init / scan --------------------------------------------- */

typedef struct {
    char *file_path;
    int   include_attr_map;
} gb_bind_t;

static void gb_bind_destroy(void *data) {
    gb_bind_t *bd = (gb_bind_t *)data;
    if (!bd) return;
    free(bd->file_path);
    free(bd);
}

static void read_genbank_bind(duckdb_bind_info info) {
    gb_bind_t *bd = (gb_bind_t *)calloc(1, sizeof(gb_bind_t));
    if (!bd) { duckdb_bind_set_error(info, "Out of memory"); return; }

    duckdb_value val = duckdb_bind_get_parameter(info, 0);
    if (val && duckdb_get_type_id(duckdb_get_value_type(val)) == DUCKDB_TYPE_VARCHAR) {
        const char *path = duckdb_get_varchar(val);
        bd->file_path = strdup(path);
        duckdb_free((void *)path);
    }
    if (val) duckdb_destroy_value(&val);
    if (!bd->file_path || bd->file_path[0] == '\0') {
        duckdb_bind_set_error(info, "read_genbank requires a file path");
        gb_bind_destroy(bd);
        return;
    }

    val = duckdb_bind_get_named_parameter(info, "attributes_map");
    if (val && !duckdb_is_null_value(val)) bd->include_attr_map = duckdb_get_bool(val) ? 1 : 0;
    if (val) duckdb_destroy_value(&val);

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

    duckdb_bind_set_bind_data(info, bd, gb_bind_destroy);
}

typedef struct {
    gb_rows_t rows;
    size_t    cursor;
    idx_t    *column_ids;
    idx_t     n_projected;
    int       include_attr_map;
} gb_init_t;

static void gb_init_destroy(void *data) {
    gb_init_t *id = (gb_init_t *)data;
    if (!id) return;
    gb_rows_free(&id->rows);
    free(id->column_ids);
    free(id);
}

static void read_genbank_init(duckdb_init_info info) {
    gb_bind_t *bd = (gb_bind_t *)duckdb_init_get_bind_data(info);
    gb_init_t *id = (gb_init_t *)calloc(1, sizeof(gb_init_t));
    if (!id) { duckdb_init_set_error(info, "Out of memory"); return; }
    id->include_attr_map = bd->include_attr_map;

    htsFile *fp = hts_open(bd->file_path, "r");
    if (!fp) {
        char msg[512];
        snprintf(msg, sizeof(msg), "read_genbank: cannot open file: %s", bd->file_path);
        duckdb_init_set_error(info, msg);
        free(id);
        return;
    }
    char err[256] = "";
    gb_parse_features(fp, &id->rows, err, sizeof(err));
    hts_close(fp);

    id->cursor = 0;
    id->n_projected = duckdb_init_get_column_count(info);
    if (id->n_projected > 0) {
        id->column_ids = (idx_t *)malloc(sizeof(idx_t) * id->n_projected);
        for (idx_t i = 0; i < id->n_projected; i++)
            id->column_ids[i] = duckdb_init_get_column_index(info, i);
    }
    duckdb_init_set_init_data(info, id, gb_init_destroy);
}

static void read_genbank_scan(duckdb_function_info info, duckdb_data_chunk output) {
    gb_init_t *id = (gb_init_t *)duckdb_function_get_init_data(info);
    idx_t chunk_cols = duckdb_data_chunk_get_column_count(output);

    duckdb_vector *vecs = NULL;
    if (chunk_cols > 0) {
        vecs = (duckdb_vector *)malloc(sizeof(duckdb_vector) * chunk_cols);
        if (!vecs) { duckdb_function_set_error(info, "read_genbank: out of memory"); duckdb_data_chunk_set_size(output, 0); return; }
        for (idx_t c = 0; c < chunk_cols; c++) vecs[c] = duckdb_data_chunk_get_vector(output, c);
    }

    idx_t row = 0;
    while (row < GB_BATCH_SIZE && id->cursor < id->rows.n) {
        gb_row_t *r = &id->rows.rows[id->cursor];
        for (idx_t c = 0; c < chunk_cols; c++) {
            int col = (int)id->column_ids[c];
            switch (col) {
                case GB_COL_SEQNAME:
                    duckdb_vector_assign_string_element(vecs[c], row, r->seqname); break;
                case GB_COL_SOURCE:
                    duckdb_vector_assign_string_element(vecs[c], row, "GenBank"); break;
                case GB_COL_FEATURE:
                    duckdb_vector_assign_string_element(vecs[c], row, r->feature); break;
                case GB_COL_START:
                    ((int64_t *)duckdb_vector_get_data(vecs[c]))[row] = r->start; break;
                case GB_COL_END:
                    ((int64_t *)duckdb_vector_get_data(vecs[c]))[row] = r->end; break;
                case GB_COL_SCORE:
                    gb_set_null(vecs[c], row); break;
                case GB_COL_STRAND: {
                    char b[2] = { r->strand, '\0' };
                    duckdb_vector_assign_string_element_len(vecs[c], row, b, 1); break;
                }
                case GB_COL_FRAME: {
                    char b[2] = { r->frame, '\0' };
                    duckdb_vector_assign_string_element_len(vecs[c], row, b, 1); break;
                }
                case GB_COL_ATTRIBUTES:
                    duckdb_vector_assign_string_element(vecs[c], row, r->attr); break;
                case GB_COL_ATTRIBUTES_MAP:
                    if (id->include_attr_map) gb_fill_attr_map(vecs[c], row, r->attr);
                    break;
                default: break;
            }
        }
        id->cursor++;
        row++;
    }

    if (vecs) free(vecs);
    duckdb_data_chunk_set_size(output, row);
}

void register_read_genbank_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "read_genbank");
    duckdb_logical_type vc = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bl = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_table_function_add_parameter(tf, vc);
    duckdb_table_function_add_named_parameter(tf, "attributes_map", bl);
    duckdb_table_function_set_bind(tf, read_genbank_bind);
    duckdb_table_function_set_init(tf, read_genbank_init);
    duckdb_table_function_set_function(tf, read_genbank_scan);
    duckdb_table_function_supports_projection_pushdown(tf, true);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
    duckdb_destroy_logical_type(&vc);
    duckdb_destroy_logical_type(&bl);
}

/* ================================================================
 * genbank_to_fasta: materialize a canonical FASTA from ORIGIN blocks
 * ================================================================ */

typedef struct {
    char   *output_path;   /* duckdb-owned */
    int64_t records_written;
    int     emitted;
} gb2fa_bind_t;

static void gb2fa_destroy(void *data) {
    gb2fa_bind_t *b = (gb2fa_bind_t *)data;
    if (!b) return;
    if (b->output_path) duckdb_free(b->output_path);
    duckdb_free(b);
}

static char *gb2fa_default_output(const char *in) {
    size_t len = strlen(in);
    char *out = (char *)duckdb_malloc(len + 4);
    if (!out) return NULL;
    snprintf(out, len + 4, "%s.fa", in);
    return out;
}

/* Write one record: ">name description" then sequence wrapped at line_width. */
static void gb2fa_write_record(FILE *out, const char *contig, const char *def,
                               kstring_t *seq, int line_width) {
    fputc('>', out);
    fputs(contig[0] ? contig : "unknown", out);
    if (def && def[0]) { fputc(' ', out); fputs(def, out); }
    fputc('\n', out);
    size_t n = seq->l;
    for (size_t i = 0; i < n; i += (size_t)line_width) {
        size_t w = (n - i < (size_t)line_width) ? (n - i) : (size_t)line_width;
        fwrite(seq->s + i, 1, w, out);
        fputc('\n', out);
    }
}

static void genbank_to_fasta_bind(duckdb_bind_info info) {
    duckdb_value pv = duckdb_bind_get_parameter(info, 0);
    char *input_path = (pv && duckdb_get_type_id(duckdb_get_value_type(pv)) == DUCKDB_TYPE_VARCHAR)
                         ? duckdb_get_varchar(pv) : NULL;
    if (pv) duckdb_destroy_value(&pv);
    if (!input_path || input_path[0] == '\0') {
        duckdb_bind_set_error(info, "genbank_to_fasta requires a file path");
        if (input_path) duckdb_free(input_path);
        return;
    }

    char *output_path = NULL;
    int line_width = 70;
    int overwrite = 0;
    duckdb_value v;
    v = duckdb_bind_get_named_parameter(info, "output_path");
    if (v && !duckdb_is_null_value(v)) output_path = duckdb_get_varchar(v);
    if (v) duckdb_destroy_value(&v);
    v = duckdb_bind_get_named_parameter(info, "line_width");
    if (v && !duckdb_is_null_value(v)) { line_width = (int)duckdb_get_int64(v); if (line_width < 1) line_width = 70; }
    if (v) duckdb_destroy_value(&v);
    v = duckdb_bind_get_named_parameter(info, "overwrite");
    if (v && !duckdb_is_null_value(v)) overwrite = duckdb_get_bool(v) ? 1 : 0;
    if (v) duckdb_destroy_value(&v);

    if (!output_path) {
        output_path = gb2fa_default_output(input_path);
        if (!output_path) { duckdb_bind_set_error(info, "Out of memory"); duckdb_free(input_path); return; }
    }
    if (!overwrite) {
        struct stat st;
        if (stat(output_path, &st) == 0) {
            char err[512];
            snprintf(err, sizeof(err), "genbank_to_fasta: output '%s' already exists (use overwrite := TRUE)", output_path);
            duckdb_bind_set_error(info, err);
            duckdb_free(input_path); duckdb_free(output_path);
            return;
        }
    }

    htsFile *fp = hts_open(input_path, "r");
    if (!fp) {
        char err[512];
        snprintf(err, sizeof(err), "genbank_to_fasta: cannot open input %s", input_path);
        duckdb_bind_set_error(info, err);
        duckdb_free(input_path); duckdb_free(output_path);
        return;
    }
    FILE *out = fopen(output_path, "wb");
    if (!out) {
        char err[512];
        snprintf(err, sizeof(err), "genbank_to_fasta: cannot open output %s: %s", output_path, strerror(errno));
        duckdb_bind_set_error(info, err);
        hts_close(fp);
        duckdb_free(input_path); duckdb_free(output_path);
        return;
    }

    kstring_t line = {0, 0, NULL};
    kstring_t seq = {0, 0, NULL};
    kstring_t def = {0, 0, NULL};
    char locus[256] = "", accession[256] = "", version[256] = "", contig[256] = "";
    int in_origin = 0, in_def = 0, have_record = 0;
    int64_t records = 0;

    while (hts_getline(fp, '\n', &line) >= 0) {
        if (line.l == 0) continue;
        const char *s = line.s;

        if (s[0] == '/' && s[1] == '/') {
            if (have_record && seq.l > 0) {
                snprintf(contig, sizeof(contig), "%s", version[0] ? version : (accession[0] ? accession : locus));
                gb2fa_write_record(out, contig, def.s ? def.s : "", &seq, line_width);
                records++;
            }
            seq.l = 0; if (seq.s) seq.s[0] = '\0';
            def.l = 0; if (def.s) def.s[0] = '\0';
            locus[0] = accession[0] = version[0] = '\0';
            in_origin = in_def = have_record = 0;
            continue;
        }

        if (s[0] != ' ' && s[0] != '\t') {
            in_def = 0;
            if (strncmp(s, "LOCUS", 5) == 0) { gb_first_token(s, locus, sizeof(locus)); have_record = 1; }
            else if (strncmp(s, "ACCESSION", 9) == 0) gb_first_token(s, accession, sizeof(accession));
            else if (strncmp(s, "VERSION", 7) == 0) gb_first_token(s, version, sizeof(version));
            else if (strncmp(s, "DEFINITION", 10) == 0) {
                const char *p = s + 10; while (*p == ' ' || *p == '\t') p++;
                int l = (int)strlen(p); while (l && (p[l-1]=='\r'||p[l-1]=='\n'||p[l-1]==' ')) l--;
                kputsn(p, l, &def);
                in_def = 1;
            }
            else if (strncmp(s, "ORIGIN", 6) == 0) in_origin = 1;
            else in_origin = 0;
            continue;
        }

        /* indented line */
        if (in_origin) {
            for (const char *p = s; *p; p++) {
                char c = *p;
                if ((c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || c == '*' || c == '-')
                    kputc((char)toupper((unsigned char)c), &seq);
            }
        } else if (in_def) {
            const char *p = s; while (*p == ' ' || *p == '\t') p++;
            int l = (int)strlen(p); while (l && (p[l-1]=='\r'||p[l-1]=='\n'||p[l-1]==' ')) l--;
            if (l) { kputc(' ', &def); kputsn(p, l, &def); }
        }
    }
    if (have_record && seq.l > 0) {
        snprintf(contig, sizeof(contig), "%s", version[0] ? version : (accession[0] ? accession : locus));
        gb2fa_write_record(out, contig, def.s ? def.s : "", &seq, line_width);
        records++;
    }

    free(line.s); free(seq.s); free(def.s);
    fclose(out);
    hts_close(fp);

    if (records == 0) {
        duckdb_bind_set_error(info, "genbank_to_fasta: no sequence records found in input");
        unlink(output_path);
        duckdb_free(input_path); duckdb_free(output_path);
        return;
    }

    duckdb_logical_type bl = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type vc = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bi = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_bind_add_result_column(info, "success", bl);
    duckdb_bind_add_result_column(info, "output_path", vc);
    duckdb_bind_add_result_column(info, "records_written", bi);
    duckdb_destroy_logical_type(&bl);
    duckdb_destroy_logical_type(&vc);
    duckdb_destroy_logical_type(&bi);

    gb2fa_bind_t *bind = (gb2fa_bind_t *)duckdb_malloc(sizeof(gb2fa_bind_t));
    bind->output_path = output_path;
    bind->records_written = records;
    bind->emitted = 0;
    duckdb_bind_set_bind_data(info, bind, gb2fa_destroy);
    duckdb_free(input_path);
}

static void genbank_to_fasta_init(duckdb_init_info info) {
    gb2fa_bind_t *b = (gb2fa_bind_t *)duckdb_init_get_bind_data(info);
    b->emitted = 0;
}

static void genbank_to_fasta_scan(duckdb_function_info info, duckdb_data_chunk output) {
    gb2fa_bind_t *b = (gb2fa_bind_t *)duckdb_function_get_bind_data(info);
    if (b->emitted) { duckdb_data_chunk_set_size(output, 0); return; }
    duckdb_vector sv = duckdb_data_chunk_get_vector(output, 0);
    duckdb_vector pv = duckdb_data_chunk_get_vector(output, 1);
    duckdb_vector rv = duckdb_data_chunk_get_vector(output, 2);
    ((bool *)duckdb_vector_get_data(sv))[0] = true;
    duckdb_vector_assign_string_element(pv, 0, b->output_path ? b->output_path : "");
    ((int64_t *)duckdb_vector_get_data(rv))[0] = b->records_written;
    b->emitted = 1;
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
