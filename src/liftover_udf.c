#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <ctype.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/khash.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <htslib/regidx.h>

#define LIFTOVER_WARNING_BUF 256
#define LIFTOVER_CACHE_MAX_ENTRIES 8
typedef struct {
    int t_start;
    int q_start;
    int size;
    int chain_ind;
    int t_start_gap;
    int t_end_gap;
    int q_start_gap;
    int q_end_gap;
} lo_block_t;

typedef struct {
    uint64_t score;
    char *t_name;
    int t_size;
    int t_start;
    int t_end;
    char *q_name;
    int q_size;
    int q_strand;
    int q_start;
    int q_end;
    int id;
    int block_ind;
    int n_blocks;
} lo_chain_t;

typedef struct {
    char *chain_path;
    char *dst_fasta_ref;
    char *src_fasta_ref;
    int max_snp_gap;
    int max_indel_inc;
    lo_chain_t *chains;
    int n_chains;
    lo_block_t *blocks;
    regidx_t *idx;
    faidx_t *src_fai;
    faidx_t *dst_fai;
} liftover_bind_t;

typedef struct liftover_cache_entry {
    liftover_bind_t *bind;
    struct liftover_cache_entry *next;
} liftover_cache_entry_t;

static void liftover_bind_destroy(void *ptr);

/*
 * faidx_t owns mutable internal I/O state and cannot be shared safely across
 * concurrent DuckDB worker threads. Keep one liftover context cache per thread.
 */
static pthread_key_t g_liftover_cache_key;
static pthread_once_t g_liftover_cache_key_once = PTHREAD_ONCE_INIT;

static void destroy_liftover_cache(void *ptr) {
    liftover_cache_entry_t *entry = (liftover_cache_entry_t *)ptr;
    while (entry) {
        liftover_cache_entry_t *next = entry->next;
        liftover_bind_destroy(entry->bind);
        free(entry);
        entry = next;
    }
}

static void init_liftover_cache_key(void) {
    pthread_key_create(&g_liftover_cache_key, destroy_liftover_cache);
}

static liftover_cache_entry_t *get_liftover_cache_head(void) {
    pthread_once(&g_liftover_cache_key_once, init_liftover_cache_key);
    return (liftover_cache_entry_t *)pthread_getspecific(g_liftover_cache_key);
}

static void set_liftover_cache_head(liftover_cache_entry_t *head) {
    pthread_once(&g_liftover_cache_key_once, init_liftover_cache_key);
    pthread_setspecific(g_liftover_cache_key, head);
}

static void trim_liftover_cache_head(liftover_cache_entry_t *head) {
    size_t count = 0;
    liftover_cache_entry_t *entry = head;
    liftover_cache_entry_t *prev = NULL;
    while (entry) {
        count++;
        if (count > LIFTOVER_CACHE_MAX_ENTRIES) {
            if (prev) prev->next = NULL;
            destroy_liftover_cache(entry);
            return;
        }
        prev = entry;
        entry = entry->next;
    }
}

enum {
    OUT_SRC_CHROM = 0,
    OUT_SRC_POS,
    OUT_SRC_REF,
    OUT_SRC_ALT,
    OUT_DEST_CHROM,
    OUT_DEST_POS,
    OUT_DEST_END,
    OUT_DEST_REF,
    OUT_DEST_ALT,
    OUT_MAPPED,
    OUT_REVERSE_COMPLEMENTED,
    OUT_SWAP,
    OUT_REJECT_REASON,
    OUT_NOTE,
    OUT_COUNT
};

static const char *LIFTOVER_FIELD_NAMES[OUT_COUNT] = {
    "src_chrom", "src_pos", "src_ref", "src_alt",
    "dest_chrom", "dest_pos", "dest_end", "dest_ref", "dest_alt",
    "mapped", "reverse_complemented", "swap", "reject_reason", "note"
};

static inline int row_is_valid(duckdb_vector vector, idx_t row) {
    uint64_t *validity = duckdb_vector_get_validity(vector);
    return !validity || ((validity[row / 64] >> (row % 64)) & 1ULL);
}

static inline void set_null_at(duckdb_vector vector, idx_t row) {
    duckdb_vector_ensure_validity_writable(vector);
    uint64_t *validity = duckdb_vector_get_validity(vector);
    validity[row / 64] &= ~((uint64_t)1 << (row % 64));
}

static inline const char *get_string_at(duckdb_vector vector, idx_t row, idx_t *len) {
    duckdb_string_t *data = (duckdb_string_t *)duckdb_vector_get_data(vector);
    duckdb_string_t *val = &data[row];
    *len = duckdb_string_t_length(*val);
    return duckdb_string_t_data(val);
}

static inline int64_t get_int64_at(duckdb_vector vector, idx_t row) {
    return ((int64_t *)duckdb_vector_get_data(vector))[row];
}

static inline int32_t get_int32_at(duckdb_vector vector, idx_t row) {
    return ((int32_t *)duckdb_vector_get_data(vector))[row];
}

static char *dup_cstr(const char *s) {
    if (!s) return NULL;
    size_t n = strlen(s);
    char *out = (char *)malloc(n + 1);
    if (!out) return NULL;
    memcpy(out, s, n + 1);
    return out;
}

static char *dup_span(const char *s, size_t n) {
    char *out = (char *)malloc(n + 1);
    if (!out) return NULL;
    memcpy(out, s, n);
    out[n] = '\0';
    return out;
}

static void free_allele_array(char **alleles, int n_allele) {
    if (!alleles) return;
    for (int i = 0; i < n_allele; i++) free(alleles[i]);
    free(alleles);
}

static int split_alt_csv(const char *alt_csv, char ***out_alts, int *out_n_alt) {
    if (strcmp(alt_csv, ".") == 0) {
        *out_alts = NULL;
        *out_n_alt = 0;
        return 1;
    }
    char **alts = NULL;
    int n_alt = 0;
    int m_alt = 0;
    const char *p = alt_csv;
    const char *tok = alt_csv;

    while (1) {
        if (*p == ',' || *p == '\0') {
            size_t tok_len = (size_t)(p - tok);
            char *a = dup_span(tok, tok_len);
            if (!a) {
                free_allele_array(alts, n_alt);
                return 0;
            }
            if (n_alt == m_alt) {
                int new_m = m_alt ? m_alt * 2 : 4;
                char **tmp = (char **)realloc(alts, sizeof(char *) * (size_t)new_m);
                if (!tmp) {
                    free(a);
                    free_allele_array(alts, n_alt);
                    return 0;
                }
                alts = tmp;
                m_alt = new_m;
            }
            alts[n_alt++] = a;
            if (*p == '\0') break;
            tok = p + 1;
        }
        p++;
    }

    *out_alts = alts;
    *out_n_alt = n_alt;
    return 1;
}

static int make_allele_array(const char *ref, const char *alt_csv, char ***out_alleles, int *out_n_allele) {
    char **alleles = NULL;
    int n_allele = 0;
    int n_alt = 0;
    char **alts = NULL;

    if (ref && *ref) n_allele++;
    if (alt_csv && *alt_csv) {
        if (!split_alt_csv(alt_csv, &alts, &n_alt)) return 0;
        n_allele += n_alt;
    }

    if (n_allele == 0) {
        *out_alleles = NULL;
        *out_n_allele = 0;
        return 1;
    }

    alleles = (char **)calloc((size_t)n_allele, sizeof(char *));
    if (!alleles) {
        free_allele_array(alts, n_alt);
        return 0;
    }

    int idx = 0;
    if (ref && *ref) {
        alleles[idx] = dup_cstr(ref);
        if (!alleles[idx]) {
            free(alleles);
            free_allele_array(alts, n_alt);
            return 0;
        }
        idx++;
    }
    for (int i = 0; i < n_alt; i++, idx++) {
        alleles[idx] = alts[i];
        alts[i] = NULL;
    }
    free_allele_array(alts, n_alt);

    *out_alleles = alleles;
    *out_n_allele = n_allele;
    return 1;
}

static int scalar_is_symbolic_alleles(char **alleles, int n_allele) {
    for (int i = 0; i < n_allele; i++) {
        if (!alleles[i]) continue;
        if (alleles[i][0] == '*' || alleles[i][0] == '<') return 1;
    }
    return 0;
}

static int scalar_is_snp_alleles(char **alleles, int n_allele) {
    if (!alleles || n_allele < 2) return 0;
    if (scalar_is_symbolic_alleles(alleles, n_allele)) return 0;
    for (int i = 0; i < n_allele; i++) {
        if (!alleles[i]) return 0;
        if (strlen(alleles[i]) != 1) return 0;
    }
    return 1;
}

static void scalar_normalize_alleles(char ***alleles_ptr, int *n_allele_ptr) {
    char **alleles = *alleles_ptr;
    int n_allele = *n_allele_ptr;
    if (!alleles || n_allele < 2) return;

    char **norm = (char **)calloc((size_t)n_allele, sizeof(char *));
    if (!norm) return;
    norm[0] = dup_cstr(alleles[0]);
    if (!norm[0]) {
        free(norm);
        return;
    }
    int n_norm = 1;

    for (int i = 1; i < n_allele; i++) {
        if (!alleles[i]) continue;
        if (strcmp(alleles[i], alleles[0]) == 0) continue;
        int seen = 0;
        for (int j = 1; j < n_norm; j++) {
            if (strcmp(alleles[i], norm[j]) == 0) {
                seen = 1;
                break;
            }
        }
        if (seen) continue;
        norm[n_norm] = dup_cstr(alleles[i]);
        if (!norm[n_norm]) {
            free_allele_array(norm, n_norm);
            return;
        }
        n_norm++;
    }

    free_allele_array(alleles, n_allele);
    *alleles_ptr = norm;
    *n_allele_ptr = n_norm;
}

static char *join_alt_csv(char **alleles, int n_allele) {
    if (!alleles || n_allele <= 1) return NULL;

    size_t total = 0;
    for (int i = 1; i < n_allele; i++) {
        if (!alleles[i]) continue;
        total += strlen(alleles[i]);
        if (i + 1 < n_allele) total++;
    }

    if (total == 0) return NULL;
    char *out = (char *)malloc(total + 1);
    if (!out) return NULL;

    size_t off = 0;
    for (int i = 1; i < n_allele; i++) {
        if (!alleles[i]) continue;
        size_t len = strlen(alleles[i]);
        memcpy(out + off, alleles[i], len);
        off += len;
        if (i + 1 < n_allele) out[off++] = ',';
    }
    out[off] = '\0';
    return out;
}

static char *fmt_message(const char *prefix, const char *detail) {
    size_t prefix_len;
    size_t detail_len;
    char *out;
    if (!prefix) return dup_cstr(detail);
    if (!detail || !*detail) return dup_cstr(prefix);
    prefix_len = strlen(prefix);
    detail_len = strlen(detail);
    out = (char *)malloc(prefix_len + 2 + detail_len + 1);
    if (!out) return NULL;
    memcpy(out, prefix, prefix_len);
    out[prefix_len] = ':';
    out[prefix_len + 1] = ' ';
    memcpy(out + prefix_len + 2, detail, detail_len + 1);
    return out;
}

static char *fmt_path_message(const char *prefix, const char *path) {
    size_t prefix_len;
    size_t path_len;
    char *out;
    if (!path || !*path) return dup_cstr(prefix);
    prefix_len = strlen(prefix);
    path_len = strlen(path);
    out = (char *)malloc(prefix_len + 2 + path_len + 1);
    if (!out) return NULL;
    memcpy(out, prefix, prefix_len);
    out[prefix_len] = ':';
    out[prefix_len + 1] = ' ';
    memcpy(out + prefix_len + 2, path, path_len + 1);
    return out;
}

static char *canonical_contig_name(const char *chr) {
    if (!chr || !*chr) return NULL;
    const char *s = chr;
    if (strncasecmp(s, "chr", 3) == 0 && s[3] != '\0') s += 3;
    /* Upstream bcf_hdr_name2id_flexible(): "26" and "MT" and "chrM" all map to MT */
    if (strcasecmp(s, "M") == 0 || strcasecmp(s, "MT") == 0 || strcmp(s, "26") == 0) return dup_cstr("MT");
    /* Upstream: "23", "25", "XY", "XX", "PAR1", "PAR2" all map to X */
    if (strcmp(s, "23") == 0 || strcmp(s, "25") == 0 || strcasecmp(s, "X") == 0 ||
        strcasecmp(s, "XY") == 0 || strcasecmp(s, "XX") == 0 ||
        strcasecmp(s, "PAR1") == 0 || strcasecmp(s, "PAR2") == 0) {
        return dup_cstr("X");
    }
    if (strcmp(s, "24") == 0 || strcasecmp(s, "Y") == 0) return dup_cstr("Y");
    return dup_cstr(s);
}

static void seq_to_upper_ascii(char *seq) {
    if (!seq) return;
    while (*seq) {
        *seq = (char)toupper((unsigned char)*seq);
        seq++;
    }
}

/* Exact port of upstream replace_iupac_codes() from bcftools/vcfnorm.c via liftover.c.
   Replaces any non-ACGTN ambiguity codes with N after uppercasing.
   Returns the count of replaced characters. */
static inline int replace_iupac_codes(char *seq, int nseq) {
    int i, n = 0;
    for (i = 0; i < nseq; i++) {
        char c = toupper(seq[i]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
            seq[i] = 'N';
            n++;
        }
    }
    return n;
}

/* Full IUPAC complement table — exact port of upstream rev_nt() from liftover.c.
   Handles all IUPAC ambiguity codes (R↔Y, S↔S, W↔W, K↔M, B↔V, D↔H, N↔N),
   bracket chars (]↔[), and special chars (-, /).
   Non-IUPAC positions are identity-mapped (no silent fallback to N). */
static inline char dna_complement(char iupac) {
    static const char iupac_complement[128] = {
        0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0A, 0x0B, 0x0C, 0x0D, 0x0E, 0x0F,
        0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1A, 0x1B, 0x1C, 0x1D, 0x1E, 0x1F,
        0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27, 0x28, 0x29, 0x2A, 0x2B, 0x2C, '-',  0x2E, '/',
        0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3A, 0x3B, 0x3C, 0x3D, 0x3E, 0x3F,
        0x40, 'T',  'V',  'G',  'H',  0x45, 0x46, 'C',  'D',  0x49, 0x4A, 'M',  0x4C, 'K',  'N',  0x4F,
        0x50, 0x51, 'Y',  'S',  'A',  0x55, 'B',  'W',  0x58, 'R',  0x5A, ']',  0x5C, '[',  0x5E, 0x5F,
        0x60, 't',  'v',  'g',  'h',  0x65, 0x66, 'c',  'd',  0x69, 0x6A, 'm',  0x6C, 'k',  'n',  0x6F,
        0x70, 0x71, 'y',  's',  'a',  0x75, 'b',  'w',  0x78, 'r',  0x7A, 0x7B, 0x7C, 0x7D, 0x7E, 0x7F,
    };
    return iupac_complement[(int)(iupac & 0x7F)];
}

static char *reverse_complement_copy(const char *s) {
    size_t len;
    char *out;
    if (!s) return NULL;
    len = strlen(s);
    out = (char *)malloc(len + 1);
    if (!out) return NULL;
    for (size_t i = 0; i < len; i++) {
        out[i] = dna_complement(s[len - i - 1]);
    }
    out[len] = '\0';
    return out;
}

static void tag_append(char *buf, size_t buf_size, const char *tag) {
    size_t cur;
    if (!buf || !tag || !*tag) return;
    cur = strlen(buf);
    if (cur > 0 && cur + 1 < buf_size) {
        buf[cur++] = ';';
        buf[cur] = '\0';
    }
    if (cur + strlen(tag) + 1 >= buf_size) return;
    memcpy(buf + cur, tag, strlen(tag) + 1);
}

/* Check if a canonical contig name is mitochondrial.
 * canonical_contig_name() already normalizes MT/chrM/M/26 → "MT". */
static int is_mt_contig(const char *canonical) {
    return canonical && strcmp(canonical, "MT") == 0;
}

/* Find the destination MT contig name in the FASTA index.
 * Returns a malloc'd string on success, NULL if not found. */
static char *find_dst_mt_name(faidx_t *fai) {
    static const char *mt_names[] = {"chrM", "MT", "M", NULL};
    for (int i = 0; mt_names[i]; i++) {
        if (faidx_has_seq(fai, mt_names[i]) > 0) return dup_cstr(mt_names[i]);
    }
    return NULL;
}

/* Get the source MT contig size from chain data.
 * Returns 0 if no MT chain found. */
static int get_chain_mt_size(const liftover_bind_t *bind) {
    for (int i = 0; i < bind->n_chains; i++) {
        char *canon = canonical_contig_name(bind->chains[i].t_name);
        if (canon && strcmp(canon, "MT") == 0) {
            int sz = bind->chains[i].t_size;
            free(canon);
            return sz;
        }
        free(canon);
    }
    return 0;
}

static int has_source_contig(const liftover_bind_t *bind, const char *chrom) {
    char *canon = canonical_contig_name(chrom);
    int found = 0;
    if (!bind || !canon) return 0;
    for (int i = 0; i < bind->n_chains; i++) {
        if (strcmp(bind->chains[i].t_name, canon) == 0) {
            found = 1;
            break;
        }
    }
    free(canon);
    return found;
}

static inline const lo_block_t *prev_block(const lo_block_t *block, const lo_block_t *blocks, const lo_chain_t *chains) {
    int ind = (int)((block - blocks) - chains[block->chain_ind].block_ind);
    return ind == 0 ? NULL : block - 1;
}

static inline const lo_block_t *next_block(const lo_block_t *block, const lo_block_t *blocks, const lo_chain_t *chains) {
    int ind = (int)((block - blocks) - chains[block->chain_ind].block_ind);
    return ind == chains[block->chain_ind].n_blocks - 1 ? NULL : block + 1;
}

static char *fetch_sequence_flexible(faidx_t *fai, const char *chrom, hts_pos_t start, hts_pos_t end) {
    static const char *mt_aliases[] = {"MT", "chrM", "M", NULL};
    char with_chr[512];
    const char *aliases[8];
    hts_pos_t len = 0;
    char *ref = NULL;
    int idx = 0;
    if (!fai || !chrom || start < 1 || end < start) return NULL;

    aliases[idx++] = chrom;
    if (strncasecmp(chrom, "chr", 3) != 0) {
        snprintf(with_chr, sizeof(with_chr), "chr%s", chrom);
        aliases[idx++] = with_chr;
    } else if (chrom[3] != '\0') {
        aliases[idx++] = chrom + 3;
    }
    if (strcasecmp(chrom, "MT") == 0 || strcasecmp(chrom, "M") == 0 || strcasecmp(chrom, "chrM") == 0) {
        for (int i = 0; mt_aliases[i]; i++) aliases[idx++] = mt_aliases[i];
    }
    aliases[idx] = NULL;

    for (int i = 0; aliases[i]; i++) {
        if (faidx_has_seq(fai, aliases[i]) <= 0) {
            continue;
        }
        ref = faidx_fetch_seq64(fai, aliases[i], start - 1, end - 1, &len);
        if (ref && len == end - start + 1) {
            seq_to_upper_ascii(ref);
            replace_iupac_codes(ref, (int)len);
            return ref;
        }
        free(ref);
        ref = NULL;
    }
    return NULL;
}

static int read_chains(htsFile *fp, int max_snp_gap, lo_chain_t **chains, lo_block_t **blocks, char **errbuf) {
    int n_chains = 0, n_blocks = 0, m_chains = 0, m_blocks = 0;
    int moff = 0, *off = NULL;
    char *tmp = NULL;
    kstring_t str = {0, 0, NULL};

    while (hts_getline(fp, KS_SEP_LINE, &str) >= 0) {
        lo_chain_t *chain;
        int t_start = 0, q_start = 0, dt = -1, dq = -1;
        int merge_in_progress = 0;
        if (str.l == 0) continue;

        if (n_chains >= m_chains) {
            lo_chain_t *new_chains;
            m_chains = m_chains ? m_chains * 2 : 64;
            new_chains = (lo_chain_t *)realloc(*chains, (size_t)m_chains * sizeof(lo_chain_t));
            if (!new_chains) {
                *errbuf = dup_cstr("bcftools_liftover: out of memory while growing chain list");
                free(off);
                free(str.s);
                return -1;
            }
            *chains = new_chains;
        }
        chain = &(*chains)[n_chains];
        memset(chain, 0, sizeof(*chain));
        int ncols = ksplit_core(str.s, 0, &moff, &off);
        if (ncols != 13 || strcmp(&str.s[off[0]], "chain") != 0) {
            *errbuf = dup_cstr("bcftools_liftover: malformed chain header");
            free(off);
            free(str.s);
            return -1;
        }

        chain->score = (uint64_t)strtoull(&str.s[off[1]], &tmp, 10);
        if (*tmp) {
            *errbuf = dup_cstr("bcftools_liftover: failed to parse chain score");
            free(off);
            free(str.s);
            return -1;
        }
        chain->t_name = canonical_contig_name(&str.s[off[2]]);
        chain->t_size = (int)strtol(&str.s[off[3]], &tmp, 10);
        if (*tmp || !chain->t_name) {
            *errbuf = dup_cstr("bcftools_liftover: failed to parse source contig");
            free(off);
            free(str.s);
            return -1;
        }
        if (str.s[off[4]] != '+') {
            *errbuf = dup_cstr("bcftools_liftover: chain source strand must be +");
            free(off);
            free(str.s);
            return -1;
        }
        chain->t_start = (int)strtol(&str.s[off[5]], &tmp, 10);
        if (*tmp) {
            *errbuf = dup_cstr("bcftools_liftover: failed to parse source start");
            free(off);
            free(str.s);
            return -1;
        }
        chain->t_end = (int)strtol(&str.s[off[6]], &tmp, 10);
        if (*tmp) {
            *errbuf = dup_cstr("bcftools_liftover: failed to parse source end");
            free(off);
            free(str.s);
            return -1;
        }
        chain->q_name = dup_cstr(&str.s[off[7]]);
        chain->q_size = (int)strtol(&str.s[off[8]], &tmp, 10);
        if (*tmp || !chain->q_name) {
            *errbuf = dup_cstr("bcftools_liftover: failed to parse destination contig");
            free(off);
            free(str.s);
            return -1;
        }
        if (str.s[off[9]] != '+' && str.s[off[9]] != '-') {
            *errbuf = dup_cstr("bcftools_liftover: chain destination strand must be +/-");
            free(off);
            free(str.s);
            return -1;
        }
        chain->q_strand = (str.s[off[9]] == '-');
        chain->q_start = (int)strtol(&str.s[off[10]], &tmp, 10);
        if (*tmp) {
            *errbuf = dup_cstr("bcftools_liftover: failed to parse destination start");
            free(off);
            free(str.s);
            return -1;
        }
        chain->q_end = (int)strtol(&str.s[off[11]], &tmp, 10);
        if (*tmp) {
            *errbuf = dup_cstr("bcftools_liftover: failed to parse destination end");
            free(off);
            free(str.s);
            return -1;
        }
        chain->id = (int)strtol(&str.s[off[12]], &tmp, 10);
        if (*tmp) {
            *errbuf = dup_cstr("bcftools_liftover: failed to parse chain id");
            free(off);
            free(str.s);
            return -1;
        }
        chain->block_ind = n_blocks;
        chain->n_blocks = 0;

        while (hts_getline(fp, KS_SEP_LINE, &str) > 0) {
            lo_block_t *block;
            int bncols = ksplit_core(str.s, 0, &moff, &off);
            int size;
            if (bncols != 1 && bncols != 3) {
                *errbuf = dup_cstr("bcftools_liftover: malformed chain block");
                free(off);
                free(str.s);
                return -1;
            }
            size = (int)strtol(&str.s[off[0]], &tmp, 10);
            if (*tmp) {
                *errbuf = dup_cstr("bcftools_liftover: failed to parse block size");
                free(off);
                free(str.s);
                return -1;
            }
            if (merge_in_progress) {
                block = &(*blocks)[n_blocks - 1];
                block->size += size;
            } else {
                if (n_blocks >= m_blocks) {
                    lo_block_t *new_blocks;
                    m_blocks = m_blocks ? m_blocks * 2 : 256;
                    new_blocks = (lo_block_t *)realloc(*blocks, (size_t)m_blocks * sizeof(lo_block_t));
                    if (!new_blocks) {
                        *errbuf = dup_cstr("bcftools_liftover: out of memory while growing block list");
                        free(off);
                        free(str.s);
                        return -1;
                    }
                    *blocks = new_blocks;
                }
                block = &(*blocks)[n_blocks++];
                block->t_start = t_start;
                block->q_start = q_start;
                block->size = size;
                block->chain_ind = n_chains;
                block->t_start_gap = dt;
                block->t_end_gap = -1;
                block->q_start_gap = dq;
                block->q_end_gap = -1;
                chain->n_blocks++;
            }

            if (bncols == 1) {
                t_start += size;
                q_start += size;
                /* Upstream chain coverage validation (liftover.c:307-313):
                   verify that accumulated block sizes fully cover the chain interval */
                if (chain->t_start + t_start != chain->t_end) {
                    *errbuf = dup_cstr("bcftools_liftover: chain malformed, target interval not fully covered");
                    free(off);
                    free(str.s);
                    return -1;
                }
                if (chain->q_start + q_start != chain->q_end) {
                    *errbuf = dup_cstr("bcftools_liftover: chain malformed, query interval not fully covered");
                    free(off);
                    free(str.s);
                    return -1;
                }
                break;
            }

            dt = (int)strtol(&str.s[off[1]], &tmp, 10);
            if (*tmp) {
                *errbuf = dup_cstr("bcftools_liftover: failed to parse source gap");
                free(off);
                free(str.s);
                return -1;
            }
            dq = (int)strtol(&str.s[off[2]], &tmp, 10);
            if (*tmp) {
                *errbuf = dup_cstr("bcftools_liftover: failed to parse destination gap");
                free(off);
                free(str.s);
                return -1;
            }
            t_start += size + dt;
            q_start += size + dq;
            merge_in_progress = (dt == dq && dt <= max_snp_gap);
            if (merge_in_progress) {
                block->size += dt;
            } else {
                block->t_end_gap = dt;
                block->q_end_gap = dq;
            }
        }
        n_chains++;
    }

    free(off);
    free(str.s);
    return n_chains;
}

/* Exact port of upstream regidx_init_chains() from liftover.c:360-388.
   - Uses uint64_t payload (not int) matching upstream sizeof(uint64_t).
   - Deduplicates chains by id using khash (Ensembl chain files contain duplicates).
   - Chains with id==0 are never deduplicated (upstream convention). */
KHASH_MAP_INIT_INT(lo_chaindedup, char)
static regidx_t *regidx_init_chains(const lo_chain_t *chains, int n_chains, const lo_block_t *blocks) {
    khash_t(lo_chaindedup) *h = kh_init(lo_chaindedup);
    regidx_t *idx = regidx_init(NULL, NULL, NULL, sizeof(uint64_t), NULL);
    if (!idx) { kh_destroy(lo_chaindedup, h); return NULL; }
    for (int i = 0; i < n_chains; i++) {
        const lo_chain_t *chain = &chains[i];

        /* Check whether the chain has already been added (Ensembl duplicates chains).
           Chains with id==0 are not deduplicated — exact upstream behavior. */
        if (chain->id) {
            khiter_t k = kh_get(lo_chaindedup, h, chain->id);
            int is_missing = (k == kh_end(h));
            if (!is_missing) continue;
            int ret;
            kh_put(lo_chaindedup, h, chain->id, &ret);
        }

        int len = (int)strlen(chain->t_name);
        for (int j = 0; j < chain->n_blocks; j++) {
            uint64_t block_ind = (uint64_t)(chain->block_ind + j);
            const lo_block_t *block = &blocks[block_ind];
            regidx_push(idx, (char *)chain->t_name, (char *)chain->t_name + len,
                        chain->t_start + block->t_start + 1,
                        chain->t_start + block->t_start + block->size, (void *)&block_ind);
        }
    }
    kh_destroy(lo_chaindedup, h);
    return idx;
}

static void liftover_bind_destroy(void *ptr) {
    liftover_bind_t *bind = (liftover_bind_t *)ptr;
    if (!bind) return;
    free(bind->chain_path);
    free(bind->dst_fasta_ref);
    free(bind->src_fasta_ref);
    if (bind->chains) {
        for (int i = 0; i < bind->n_chains; i++) {
            free(bind->chains[i].t_name);
            free(bind->chains[i].q_name);
        }
    }
    free(bind->chains);
    free(bind->blocks);
    if (bind->idx) regidx_destroy(bind->idx);
    if (bind->src_fai) fai_destroy(bind->src_fai);
    if (bind->dst_fai) fai_destroy(bind->dst_fai);
    free(bind);
}

static int liftover_bp(liftover_bind_t *bind, const char *t_chr, hts_pos_t t_pos,
                       const char **q_chr, hts_pos_t *q_pos, int *q_strand, int *multi_match) {
    int block_ind = -1;
    char *canon = canonical_contig_name(t_chr);
    regitr_t *itr = NULL;
    if (!canon) return -1;
    itr = regitr_init(bind->idx);
    if (!itr) {
        free(canon);
        return -1;
    }
    *multi_match = 0;
    if (regidx_overlap(bind->idx, canon, (uint32_t)t_pos, (uint32_t)t_pos, itr)) {
        for (int i = 0; regitr_overlap(itr); i++) {
            const lo_block_t *block;
            const lo_chain_t *chain;
            int block_pos;
            if (i > 0) *multi_match = 1;
            block_ind = (int)regitr_payload(itr, uint64_t);
            block = &bind->blocks[block_ind];
            chain = &bind->chains[block->chain_ind];
            *q_chr = chain->q_name;
            *q_strand = chain->q_strand;
            block_pos = (int)(t_pos - itr->beg);
            if (*q_strand) {
                *q_pos = (hts_pos_t)(chain->q_size - chain->q_start - block->q_start - block_pos);
            } else {
                *q_pos = (hts_pos_t)(chain->q_start + block->q_start + block_pos + 1);
            }
        }
    }
    regitr_destroy(itr);
    free(canon);
    return block_ind;
}

static int liftover_indel(liftover_bind_t *bind, const char *src_chr, hts_pos_t src_pos5, hts_pos_t src_pos3,
                          int *dst_rid_unused, const char **dst_chr, hts_pos_t *dst_pos5, hts_pos_t *dst_pos3,
                          int *strand, int *npad) {
    int rid5 = -1, rid3 = -1, strand5 = 0, strand3 = 0;
    hts_pos_t pos5 = -1, pos3 = -1;
    const char *chr5 = NULL, *chr3 = NULL;
    int multi5 = 0, multi3 = 0;
    int block_ind5 = liftover_bp(bind, src_chr, src_pos5, &chr5, &pos5, &strand5, &multi5);
    const lo_block_t *block5 = block_ind5 < 0 ? NULL : &bind->blocks[block_ind5];
    int block_ind3 = liftover_bp(bind, src_chr, src_pos3, &chr3, &pos3, &strand3, &multi3);
    const lo_block_t *block3 = block_ind3 < 0 ? NULL : &bind->blocks[block_ind3];
    (void)rid5;
    (void)rid3;
    (void)dst_rid_unused;
    if (!block5 && !block3) return -1;

    if (!block5) {
        const lo_chain_t *aux_chain = &bind->chains[block3->chain_ind];
        int dst_chr_size = aux_chain->q_size;
        const lo_block_t *aux_block;
        *dst_chr = chr3;
        *strand = strand3;
        aux_block = prev_block(block3, bind->blocks, bind->chains);
        *npad = aux_block ? aux_chain->t_start + aux_block->t_start + aux_block->size - (int)src_pos5 : -bind->max_indel_inc;
        if (*npad < -bind->max_indel_inc || *npad >= 0) {
            *npad = -bind->max_indel_inc;
            if (*npad > src_pos5 - src_pos3) return -2;
        }
        if (src_pos5 + *npad < 1) *npad = (int)(1 - src_pos5);
        block_ind5 = liftover_bp(bind, src_chr, src_pos5 + *npad, &chr5, &pos5, &strand5, &multi5);
        (void)block_ind5;
        int dst_npad = *strand ? (int)(pos5 - pos3) : (int)(pos3 - pos5);
        if (!chr5 || strcmp(chr5, chr3) != 0 || strand5 != strand3 || dst_npad < 0 || dst_npad > bind->max_indel_inc) {
            if (*strand) {
                if (pos3 - *npad > dst_chr_size) *npad = (int)(pos3 - dst_chr_size);
                pos5 = pos3 - *npad;
            } else {
                if (pos3 - (src_pos3 - src_pos5) + *npad < 1) *npad = (int)(1 - pos3 + (src_pos3 - src_pos5));
                pos5 = pos3 - (src_pos3 - src_pos5) + *npad;
            }
        }
    } else if (!block3) {
        const lo_chain_t *aux_chain = &bind->chains[block5->chain_ind];
        int src_chr_size = aux_chain->t_size;
        int dst_chr_size = aux_chain->q_size;
        const lo_block_t *aux_block;
        *dst_chr = chr5;
        *strand = strand5;
        aux_block = next_block(block5, bind->blocks, bind->chains);
        *npad = aux_block ? aux_chain->t_start + aux_block->t_start + 1 - (int)src_pos3 : bind->max_indel_inc;
        if (*npad > bind->max_indel_inc || *npad <= 0) {
            *npad = bind->max_indel_inc;
            if (*npad < src_pos3 - src_pos5) return -3;
        }
        if (src_pos3 + *npad > src_chr_size) *npad = (int)(src_chr_size - src_pos3);
        block_ind3 = liftover_bp(bind, src_chr, src_pos3 + *npad, &chr3, &pos3, &strand3, &multi3);
        (void)block_ind3;
        int dst_npad = *strand ? (int)(pos5 - pos3) : (int)(pos3 - pos5);
        if (!chr3 || strcmp(chr5, chr3) != 0 || strand5 != strand3 || dst_npad < 0 || dst_npad > bind->max_indel_inc) {
            if (*strand) {
                if (pos5 - *npad < 1) *npad = (int)(pos5 - 1);
                pos3 = pos5 - *npad;
            } else {
                if (pos5 + (src_pos3 - src_pos5) + *npad > dst_chr_size) {
                    *npad = (int)(dst_chr_size - pos5 - (src_pos3 - src_pos5));
                }
                pos3 = pos5 + (src_pos3 - src_pos5) + *npad;
            }
        }
    } else {
        if (block5->chain_ind != block3->chain_ind) return -4;
        if (strcmp(chr5, chr3) != 0) return -4;
        *dst_chr = chr5;
        if (strand5 != strand3) return -4;
        *strand = strand5;
        if (abs((int)(pos3 - pos5)) > src_pos3 - src_pos5 + bind->max_indel_inc) return -5;
    }

    *dst_pos5 = *strand == 0 ? pos5 : pos3;
    *dst_pos3 = *strand == 0 ? pos3 : pos5;
    return 0;
}

/* ========================================================================
 * TIER 3: Allele Extension, Needleman-Wunsch, find_reference, left-alignment
 * Exact port of upstream bcftools liftover.c indel pipeline.
 *
 * Adapted for scalar (biallelic) UDF: n_allele is always 2 (ref + alt).
 * Instead of mutating bcf1_t records, we work with mutable string copies.
 * ======================================================================== */

/* Scalar version of upstream bcf1_realign_t.
 * Holds a dynamic allele array (ref=als[0], alt(s)=als[1..]) plus the cached
 * reference sequence window, FASTA index, chromosome name, and variant pos.
 * pos is 0-based (same as bcf1_t.pos in upstream). */
typedef struct {
    const char *chrom;      /* chromosome name for FASTA lookups */
    faidx_t *fai;           /* FASTA index handle */
    kstring_t *als;         /* als[0]=ref, als[1..]=alt(s) */
    int n_allele;           /* >= 1 */
    hts_pos_t pos;          /* 0-based position */
    int aln_win;            /* window size for reference fetch (default 100) */
    char *ref;              /* cached reference sequence */
    hts_pos_t beg;          /* 0-based start of cached ref window */
    hts_pos_t end;          /* 0-based end of cached ref window */
} scalar_realign_t;

#define SCALAR_ALN_WIN 100

/* Adapted from upstream safe_fetch_sequence().
 * Uses chrom string directly instead of bcf_seqname(hdr, rec).
 * Clamps beg/end to valid range and fetches via fetch_sequence_flexible(). */
static void scalar_safe_fetch_sequence(scalar_realign_t *ra) {
    int seq_len = faidx_seq_len64(ra->fai, ra->chrom);
    if (seq_len < 0) {
        /* Try alias probing — fetch_sequence_flexible handles this */
        seq_len = 0;
    }
    if (ra->beg + 1 < 1) ra->beg = 0;
    if (ra->end + 1 > seq_len) ra->end = seq_len - 1;
    free(ra->ref);
    /* fetch_sequence_flexible uses 1-based inclusive coords */
    ra->ref = fetch_sequence_flexible(ra->fai, ra->chrom, ra->beg + 1, ra->end + 1);
}

/* Adapted from upstream initialize_alleles().
 * Takes explicit ref_str/alt_str and stores in ra->als[0..1].
 * Sets up the reference window and validates ref matches reference. */
static int scalar_initialize_alleles(scalar_realign_t *ra, char **alleles, int n_allele) {
    if (!alleles || n_allele < 1) return 0;
    ra->als = (kstring_t *)calloc((size_t)n_allele, sizeof(kstring_t));
    if (!ra->als) return 0;
    ra->n_allele = n_allele;
    for (int i = 0; i < n_allele; i++) {
        if (!alleles[i]) continue;
        kputsn_(alleles[i], strlen(alleles[i]), &ra->als[i]);
    }

    ra->beg = ra->pos - 1;
    ra->end = ra->pos + (hts_pos_t)ra->als[0].l;
    free(ra->ref);
    ra->ref = NULL;
    scalar_safe_fetch_sequence(ra);

    /* Upstream validates ref matches reference — we skip fatal error in scalar UDF,
       but the mismatch will be caught downstream by find_reference. */
    return 1;
}

/* Exact port of upstream pad_left() for scalar_realign_t */
static void scalar_pad_left(scalar_realign_t *ra, int npad) {
    if (ra->beg == 0 || ra->pos - npad < ra->beg) {
        ra->beg = ra->pos - npad;
        free(ra->ref);
        ra->ref = NULL;
        scalar_safe_fetch_sequence(ra);
    }
    if (!ra->ref) return;
    const char *ptr = &ra->ref[ra->pos - ra->beg - npad];
    for (int i = 0; i < ra->n_allele; i++) {
        if (ra->als[i].s[0] == '*') continue;
        ks_resize(&ra->als[i], ra->als[i].l + (size_t)npad);
        memmove(ra->als[i].s + npad, ra->als[i].s, ra->als[i].l);
        memcpy(ra->als[i].s, ptr, (size_t)npad);
        ra->als[i].l += (size_t)npad;
    }
    ra->pos -= npad;
}

/* Exact port of upstream pad_right() for scalar_realign_t */
static void scalar_pad_right(scalar_realign_t *ra, int npad) {
    if (ra->end < ra->pos + (hts_pos_t)ra->als[0].l + npad - 1) {
        ra->end = ra->pos + (hts_pos_t)ra->als[0].l + npad - 1;
        free(ra->ref);
        ra->ref = NULL;
        scalar_safe_fetch_sequence(ra);
    }
    if (!ra->ref) return;
    /* Bounds check: ensure we have enough reference to pad from */
    hts_pos_t ref_offset = ra->pos - ra->beg + (hts_pos_t)ra->als[0].l;
    hts_pos_t ref_available = ra->end - ra->beg + 1;
    if (ref_offset + npad > ref_available) return; /* not enough reference */
    const char *ptr = &ra->ref[ref_offset];
    for (int i = 0; i < ra->n_allele; i++) {
        if (ra->als[i].s[0] == '*') continue;
        kputsn(ptr, (size_t)npad, &ra->als[i]);
    }
}

/* Exact port of upstream pad_from_right() for scalar_realign_t */
static void scalar_pad_from_right(scalar_realign_t *ra, const char *s_ptr, int d) {
    int npad = 0;
    hts_pos_t ref_offset = ra->pos - ra->beg + (hts_pos_t)ra->als[0].l;
    hts_pos_t ref_available = ra->end - ra->beg + 1;
    if (ref_offset >= ref_available) return; /* no reference available to pad from */
    const char *l_ptr = &ra->ref[ref_offset];
    while (1) {
        if (ra->pos + (hts_pos_t)ra->als[0].l + npad > ra->end) {
            ra->end += ra->aln_win;
            free(ra->ref);
            ra->ref = NULL;
            scalar_safe_fetch_sequence(ra);
            if (!ra->ref) break;
            ref_available = ra->end - ra->beg + 1;
            l_ptr = &ra->ref[ra->pos - ra->beg + (hts_pos_t)ra->als[0].l + npad];
            if (npad >= d) s_ptr = &ra->ref[ra->pos - ra->beg + (hts_pos_t)ra->als[0].l + npad - d];
        }
        /* Bounds check */
        if (ra->pos - ra->beg + (hts_pos_t)ra->als[0].l + npad >= ref_available) break;
        if (npad == d) s_ptr = &ra->ref[ra->pos - ra->beg + (hts_pos_t)ra->als[0].l];
        npad++;
        if (*l_ptr != *s_ptr) break;
        l_ptr++;
        s_ptr++;
    }
    if (!ra->ref || npad == 0) return;
    ref_offset = ra->pos - ra->beg + (hts_pos_t)ra->als[0].l;
    if (ref_offset + npad > ref_available) npad = (int)(ref_available - ref_offset);
    if (npad <= 0) return;
    const char *ptr = &ra->ref[ref_offset];
    for (int i = 0; i < ra->n_allele; i++) {
        if (ra->als[i].s[0] == '*') continue;
        kputsn(ptr, (size_t)npad, &ra->als[i]);
    }
}

/* Exact port of upstream pad_from_left() for scalar_realign_t */
static void scalar_pad_from_left(scalar_realign_t *ra, const char *s_ptr, int d) {
    int npad = 0;
    const char *l_ptr = &ra->ref[ra->pos - ra->beg - 1];
    while (1) {
        if (ra->pos - npad <= ra->beg) {
            ra->beg -= ra->aln_win;
            free(ra->ref);
            ra->ref = NULL;
            scalar_safe_fetch_sequence(ra);
            if (!ra->ref) break;
            l_ptr = &ra->ref[ra->pos - ra->beg - npad - 1];
            if (npad >= d) s_ptr = &ra->ref[ra->pos - ra->beg - npad + d - 1];
        }
        if (npad == d) s_ptr = &ra->ref[ra->pos - ra->beg - 1];
        npad++;
        if (*l_ptr != *s_ptr) break;
        l_ptr--;
        s_ptr--;
    }
    if (!ra->ref) return;
    for (int i = 0; i < ra->n_allele; i++) {
        if (ra->als[i].s[0] == '*') continue;
        ks_resize(&ra->als[i], ra->als[i].l + (size_t)npad);
        memmove(ra->als[i].s + npad, ra->als[i].s, ra->als[i].l);
        memcpy(ra->als[i].s, l_ptr, (size_t)npad);
        ra->als[i].l += (size_t)npad;
    }
    ra->pos -= npad;
}

/* Scalar version of upstream is_extension_needed().
 * For biallelic (ref_str, alt_str), returns 1 if alleles need extension. */
static int scalar_is_extension_needed(char **alleles, int n_allele) {
    if (!alleles || n_allele < 1) return 0;
    if (n_allele == 1 && alleles[0] && strlen(alleles[0]) == 1) return 1;

    for (int i = 0; i < n_allele; i++) {
        if (!alleles[i] || alleles[i][0] == '*') continue;
        int len_i = (int)strlen(alleles[i]);
        for (int j = i + 1; j < n_allele; j++) {
            if (!alleles[j] || alleles[j][0] == '*') continue;
            int len_j = (int)strlen(alleles[j]);
            if (alleles[i][0] != alleles[j][0]) return 1;
            if (alleles[i][len_i - 1] != alleles[j][len_j - 1]) return 1;
            if (len_i == len_j) continue;
            if (len_i > len_j) {
                if (strncasecmp(alleles[i], alleles[j], (size_t)len_j) == 0) return 1;
                if (strncasecmp(alleles[i] + len_i - len_j, alleles[j], (size_t)len_j) == 0) return 1;
            } else {
                if (strncasecmp(alleles[i], alleles[j], (size_t)len_i) == 0) return 1;
                if (strncasecmp(alleles[i], alleles[j] + len_j - len_i, (size_t)len_i) == 0) return 1;
            }
        }
    }
    return 0;
}

/* Exact port of upstream extend_alleles() for scalar_realign_t.
 * Inspired by Petr Danecek's code from realign() in bcftools/vcfnorm.c. */
static void scalar_extend_alleles(scalar_realign_t *ra, char **alleles, int n_allele) {
    if (!scalar_initialize_alleles(ra, alleles, n_allele) || !ra->ref) return;

    for (int i = 0; i < ra->n_allele; i++) {
        if (ra->als[i].s[0] == '*') continue;
        for (int j = i + 1; j < ra->n_allele; j++) {
            if (ra->als[j].s[0] == '*') continue;
            if (ra->als[i].s[0] != ra->als[j].s[0]) scalar_pad_left(ra, 1);
            if (ra->als[i].s[ra->als[i].l - 1] != ra->als[j].s[ra->als[j].l - 1]) scalar_pad_right(ra, 1);
            if (ra->als[i].l == ra->als[j].l) continue;
            kstring_t *l_str = ra->als[i].l > ra->als[j].l ? &ra->als[i] : &ra->als[j];
            kstring_t *s_str = ra->als[i].l > ra->als[j].l ? &ra->als[j] : &ra->als[i];
            if (strncasecmp(l_str->s, s_str->s, s_str->l) == 0)
                scalar_pad_from_right(ra, l_str->s + s_str->l, (int)(l_str->l - s_str->l));
            if (strncasecmp(l_str->s + l_str->l - s_str->l, s_str->s, s_str->l) == 0)
                scalar_pad_from_left(ra, l_str->s + l_str->l - s_str->l - 1, (int)(l_str->l - s_str->l));
        }
    }
}

/* Scalar update_alleles: copies result from ra->als[0..1] into output strings.
 * Returns new (pos, ref, alt) via out params. Caller must free output alleles.
 * Returns the 0-based position (ra->pos). */
static hts_pos_t scalar_update_alleles(scalar_realign_t *ra, char ***out_alleles, int *out_n_allele) {
    char **alleles = (char **)calloc((size_t)ra->n_allele, sizeof(char *));
    if (!alleles) return -1;
    for (int i = 0; i < ra->n_allele; i++) {
        alleles[i] = dup_span(ra->als[i].s, ra->als[i].l);
        if (!alleles[i]) {
            free_allele_array(alleles, ra->n_allele);
            return -1;
        }
    }
    *out_alleles = alleles;
    *out_n_allele = ra->n_allele;
    return ra->pos;
}

static void scalar_realign_cleanup(scalar_realign_t *ra) {
    if (ra->als) {
        for (int i = 0; i < ra->n_allele; i++) free(ra->als[i].s);
        free(ra->als);
    }
    free(ra->ref);
    ra->als = NULL;
    ra->n_allele = 0;
    ra->ref = NULL;
}

/* ========================================================================
 * L-C2: Needleman-Wunsch Re-alignment
 * Exact port of upstream nw(), get_shift(), clip_pad() from liftover.c.
 * ======================================================================== */

#define _M 0
#define _D 1
#define _I 2
/* This version left-aligns deletions and insertions (exact upstream macro) */
#define MAX_IND(v) ((v[_D]) < (v[_I]) ? (v[_I]) < (v[_M]) ? (_M) : (_I) : (v[_D]) < (v[_M]) ? (_M) : (_D))

/* Exact port of upstream nw() — Needleman-Wunsch with affine gap penalty.
 * Scoring values drawn from bwa mem: match=1, mismatch=-4, gap_open=-6, gap_ext=-1.
 * Returns the alignment score. Path is stored in path->s. */
static int scalar_nw(const char *s, size_t s_l, const char *t, size_t t_l, kstring_t *path) {
    if (!path) return -1;
    path->l = 0;
    if (path->s) path->s[0] = '\0';
    if (s_l == 0 || t_l == 0) return -1;
    if (s_l > (size_t)(INT_MAX - 1) || t_l > (size_t)(INT_MAX - 1)) return -1;
    int a = 1;   /* score for a sequence match */
    int b = -4;  /* penalty for a mismatch */
    int o = -6;  /* gap open penalty */
    int e = -1;  /* gap extension penalty */
    int n = (int)s_l + 1;
    int m = (int)t_l + 1;
    int i, j, cell, ind, score[3];
    int *A[3];
    char *B[3];
    if (ks_resize(path, (size_t)(n + m - 1)) < 0) return -1;
    for (i = 0; i < 3; i++) {
        A[i] = (int *)malloc(sizeof(int) * (size_t)(n * m));
        B[i] = (char *)malloc(sizeof(char) * (size_t)(n * m));
        if (!A[i] || !B[i]) {
            for (int k = 0; k <= i; k++) {
                free(A[k]);
                free(B[k]);
            }
            return -1;
        }
    }

    /* Build three scoring matrices corresponding to the three possible states */
    A[_M][0] = 0;
    A[_D][0] = INT32_MIN / 2;
    A[_I][0] = INT32_MIN / 2;
    for (j = 1; j < m; j++) {
        A[_M][j] = INT32_MIN / 2;
        A[_D][j] = INT32_MIN / 2;
        A[_I][j] = o + e * j;
        B[_I][j] = _I;
    }
    for (i = 1; i < n; i++) {
        A[_M][i * m] = INT32_MIN / 2;
        A[_D][i * m] = o + e * i;
        B[_D][i * m] = _D;
        A[_I][i * m] = INT32_MIN / 2;
        for (j = 1; j < m; j++) {
            cell = i * m + j;
            score[_M] = A[_M][cell - m - 1];
            score[_D] = A[_D][cell - m - 1];
            score[_I] = A[_I][cell - m - 1];
            ind = MAX_IND(score);
            A[_M][cell] = (s[i - 1] == t[j - 1] || s[i - 1] == 'N' || t[j - 1] == 'N' ? a : b) + score[ind];
            B[_M][cell] = (char)ind;

            score[_M] = o + A[_M][cell - m];
            score[_D] = A[_D][cell - m];
            score[_I] = o + A[_I][cell - m];
            ind = MAX_IND(score);
            A[_D][cell] = e + score[ind];
            B[_D][cell] = (char)ind;

            score[_M] = o + A[_M][cell - 1];
            score[_D] = o + A[_D][cell - 1];
            score[_I] = A[_I][cell - 1];
            ind = MAX_IND(score);
            A[_I][cell] = e + score[ind];
            B[_I][cell] = (char)ind;
        }
    }

    /* Adjust the values at the end to account for 3' clipping */
    cell = m * n - 1;
    score[_M] = A[_M][cell];
    score[_D] = A[_D][cell];
    score[_I] = A[_I][cell];
    ind = MAX_IND(score);
    int ret = score[ind];

    /* Reconstruct the path by backtracking */
    i = n + m - 2;
    path->s[i] = '\0';
    while (cell > 0) {
        path->s[--i] = 1 + ind;
        int shift = (ind == _M || ind == _D ? m : 0) + (ind == _M || ind == _I ? 1 : 0);
        ind = B[ind][cell];
        cell -= shift;
    }

    /* i here represents the number of matches */
    if (i) memmove(path->s, &path->s[i], (size_t)(n + m - 1 - i));
    path->l = (size_t)(n + m - i);

    for (i = 0; i < 3; i++) {
        free(A[i]);
        free(B[i]);
    }

    return ret;
}

/* Exact port of upstream get_shift() */
static int scalar_get_shift(char *path, int npad) {
    if (!path || npad == 0) return 0;
    size_t path_len = strlen(path);
    if (path_len == 0) return 0;

    int shift = 0;
    if (npad < 0) {
        char *ptr = path;
        while (npad < 0 && *ptr) {
            int ind = (*ptr++) - 1;
            if (ind == _M || ind == _D) shift++;
            if (ind == _M || ind == _I) npad++;
        }
    } else {
        char *ptr = &path[path_len - 1];
        while (npad > 0) {
            if (ptr < path) break;
            int ind = (*ptr--) - 1;
            if (ind == _M || ind == _D) shift--;
            if (ind == _M || ind == _I) npad--;
        }
    }
    return shift;
}

/* Scalar version of upstream clip_pad().
 * For biallelic, iterates over ref and alt alleles.
 * ref_str is the destination reference sequence spanning [pos5..pos3].
 * pad is the source pad sequence, npad is its signed length.
 * Adjusts pos5/pos3 in-place. Uses allele strings from extend_alleles output. */
static void scalar_clip_pad(char **alleles, int n_allele,
                            const char *dst_ref, char *pad, int npad,
                            hts_pos_t *pos5, hts_pos_t *pos3) {
    int best_score = INT32_MIN;
    int best_shift = 0;
    kstring_t path = {0, 0, NULL};
    kstring_t src_seq = {0, 0, NULL};
    for (int i = 0; i < n_allele; i++) {
        if (!alleles[i] || alleles[i][0] == '*') continue;

        /* Generate source reference sequence: pad + allele or allele + pad */
        src_seq.l = 0;
        if (npad < 0) {
            kputsn_(&pad[0], (size_t)(-npad), &src_seq);
            kputs(alleles[i], &src_seq);
        } else {
            kputsn_(alleles[i], strlen(alleles[i]), &src_seq);
            kputsn(pad, (size_t)npad, &src_seq);
        }

        /* Pairwise align source and destination sequences excluding the anchors */
        int new_score = scalar_nw(dst_ref + 1, (size_t)(*pos3 - *pos5 - 1),
                                  src_seq.s + 1, src_seq.l - 2, &path);
        if (new_score < 0 || !path.s || path.s[0] == '\0') continue;
        if (new_score > best_score) {
            best_score = new_score;
            best_shift = scalar_get_shift(path.s, npad);
        }
    }

    if (npad < 0)
        *pos5 += best_shift;
    else
        *pos3 += best_shift;

    free(path.s);
    free(src_seq.s);
}

/* ========================================================================
 * L-C4: Sophisticated find_reference
 * Scalar version of upstream find_reference() from liftover.c:1436-1526.
 * Takes mutable allele strings and the destination reference.
 * Returns: swap value (0=no swap, >0=allele index swapped, -1=ref added).
 * Modifies allele strings in-place (anchor update, allele reordering).
 * Also updates *pos_ptr for right-match SNP rescue.
 * ======================================================================== */
static int scalar_find_reference(char ***alleles_ptr, int *n_allele_ptr,
                                 const char *dst_ref, int is_snp,
                                 hts_pos_t *pos_ptr) {
    char **alleles = *alleles_ptr;
    int n_allele = *n_allele_ptr;
    int swap = -1;
    int ref_len = (int)strlen(dst_ref);

    /* Phase 1: Direct matching — interior comparison for each allele */
    for (int i = 0; i < n_allele; i++) {
        if (!alleles[i] || alleles[i][0] == '*') continue;
        int len = (int)strlen(alleles[i]);
        if (len != ref_len) continue;
        int match;
        if (ref_len == 1) {
            match = (dst_ref[0] == alleles[i][0]);
        } else {
            match = (strncasecmp(dst_ref + 1, alleles[i] + 1, (size_t)(len - 2)) == 0);
        }
        if (match) {
            if (swap >= 0) {
                /* Two matches — ambiguous, can only happen if alleles not maximally extended */
                break;
            }
            swap = i;
        }
    }

    /* Phase 2: SNP rescue hack — when direct match fails for SNPs with ref_len > 2
     * and the allele lengths differ from ref_len, tries left-match and right-match
     * of shorter alleles (saves ~30% of SNPs from becoming complex variants) */
    if (swap < 0 && is_snp && ref_len > 2 && ref_len != (int)strlen(alleles[0])) {
        int n_matches = 0;
        for (int i = 0; i < n_allele; i++) {
            if (!alleles[i] || alleles[i][0] == '*') continue;
            int len = (int)strlen(alleles[i]);
            if (len < ref_len) {
                int left_match = (strncasecmp(dst_ref, alleles[i], (size_t)(len - 1)) == 0);
                int right_match = (strncasecmp(dst_ref + 1 + ref_len - len, alleles[i] + 1, (size_t)len) == 0);
                n_matches += left_match + right_match;
                if (n_matches > 1) break; /* too many matches */
                if (left_match || right_match) swap = i;
            }
        }
        if (n_matches == 1) {
            int len = (int)strlen(alleles[swap]);
            int left_match = (strncasecmp(dst_ref, alleles[swap], (size_t)(len - 1)) == 0);
            if (left_match) {
                ref_len = len;
            } else { /* right match */
                *pos_ptr += ref_len - len;
                dst_ref += ref_len - len;
                ref_len = len;
            }
        } else {
            swap = -1;
        }
    }

    /* Phase 3: Update left and right anchors.
     * Sometimes anchors can include SNPs between assemblies that should
     * not be incorporated in the liftover. */
    if (ref_len > 1) { /* don't update if the variant is a SNP */
        char ref5 = dst_ref[0];
        char ref3 = dst_ref[ref_len - 1];
        /* Update anchors in both alleles */
        for (int i = 0; i < n_allele; i++) {
            if (!alleles[i] || alleles[i][0] == '*') continue;
            int len = (int)strlen(alleles[i]);
            alleles[i][len - 1] = ref3;
            alleles[i][0] = ref5;
        }
    }

    /* Update alleles if necessary:
     * swap == 0: ref already matches destination ref — no swap needed
     * swap > 0: alt matches destination ref — swap ref and alt
     * swap < 0: no match found — prepend destination ref as new reference */
    if (swap != 0) {
        if (swap < 0) {
            char **new_alleles = (char **)calloc((size_t)(n_allele + 1), sizeof(char *));
            if (!new_alleles) return swap;
            new_alleles[0] = dup_span(dst_ref, (size_t)ref_len);
            if (!new_alleles[0]) {
                free(new_alleles);
                return swap;
            }
            for (int i = 0; i < n_allele; i++) {
                new_alleles[i + 1] = dup_cstr(alleles[i]);
                if (!new_alleles[i + 1]) {
                    free_allele_array(new_alleles, n_allele + 1);
                    return swap;
                }
            }
            free_allele_array(alleles, n_allele);
            *alleles_ptr = new_alleles;
            *n_allele_ptr = n_allele + 1;
        } else {
            char *tmp = alleles[0];
            alleles[0] = alleles[swap];
            alleles[swap] = tmp;
        }
        /* swap == 0 is handled implicitly — ref already matches, nothing to do */
    }

    return swap;
}

/* ========================================================================
 * L-C3: Left-alignment (Tan et al. 2015)
 * Scalar versions of is_left_aligned, trim_right, trim_left, shift_left.
 * ======================================================================== */

/* Scalar version of upstream is_left_aligned().
 * Returns 1 if the variant is already left-aligned, 0 if it needs alignment. */
static int scalar_is_left_aligned(char **alleles, int n_allele) {
    if (!alleles || n_allele < 1) return 1;
    if (n_allele == 1) return strlen(alleles[0]) <= 1;

    int *len = (int *)calloc((size_t)n_allele, sizeof(int));
    if (!len) return 1;
    int min_len_is_one = 0;
    for (int i = 0; i < n_allele; i++) {
        if (!alleles[i] || alleles[i][0] == '*') continue;
        len[i] = (int)strlen(alleles[i]);
        if (len[i] == 1) min_len_is_one = 1;
    }

    char last_ref = alleles[0][len[0] - 1];
    int i;
    for (i = 0; i < n_allele; i++) {
        if (!alleles[i] || alleles[i][0] == '*') continue;
        if (alleles[i][len[i] - 1] != last_ref) break;
    }
    if (i == n_allele) {
        free(len);
        return 0;
    }

    if (!min_len_is_one) {
        char first_ref = alleles[0][0];
        for (i = 0; i < n_allele; i++) {
            if (!alleles[i] || alleles[i][0] == '*') continue;
            if (alleles[i][0] != first_ref) break;
        }
        if (i == n_allele) {
            free(len);
            return 0;
        }
    }

    free(len);
    return 1;
}

/* Exact port of upstream trim_right() for scalar_realign_t */
static void scalar_trim_right(scalar_realign_t *ra) {
    if (ra->als[0].l <= 1) return;
    int i_trim;
    for (i_trim = 0; i_trim < (int)ra->als[0].l - 1; i_trim++) {
        char last_ref = ra->als[0].s[ra->als[0].l - 1 - i_trim];
        int all_match = 1;
        for (int j = 1; j < ra->n_allele; j++) {
            if (ra->als[j].s[0] == '*') continue;
            if (i_trim >= (int)ra->als[j].l - 1 || ra->als[j].s[ra->als[j].l - 1 - i_trim] != last_ref) {
                all_match = 0;
                break;
            }
        }
        if (!all_match) break;
    }
    for (int j = 0; j < ra->n_allele; j++) {
        if (ra->als[j].s[0] == '*') continue;
        ra->als[j].l -= (size_t)i_trim;
    }
}

/* Exact port of upstream trim_left() for scalar_realign_t */
static void scalar_trim_left(scalar_realign_t *ra) {
    if (ra->als[0].l <= 1) return;
    int i_trim;
    for (i_trim = 0; i_trim < (int)ra->als[0].l - 1; i_trim++) {
        char first_ref = ra->als[0].s[i_trim];
        int all_match = 1;
        for (int j = 1; j < ra->n_allele; j++) {
            if (ra->als[j].s[0] == '*') continue;
            if (i_trim >= (int)ra->als[j].l - 1 || ra->als[j].s[i_trim] != first_ref) {
                all_match = 0;
                break;
            }
        }
        if (!all_match) break;
    }
    for (int j = 0; j < ra->n_allele; j++) {
        if (ra->als[j].s[0] == '*') continue;
        ra->als[j].l -= (size_t)i_trim;
        memmove(ra->als[j].s, ra->als[j].s + i_trim, ra->als[j].l);
    }
    ra->pos += i_trim;
}

/* Exact port of upstream shift_left() for scalar_realign_t */
static void scalar_shift_left(scalar_realign_t *ra) {
    int n_allele_non_star = 0;
    for (int i = 1; i < ra->n_allele; i++)
        if (ra->als[i].s[0] != '*') n_allele_non_star++;
    if (n_allele_non_star == 0) return;

    while (ra->pos > 0) {
        char last_ref = ra->als[0].s[ra->als[0].l - 1];
        int all_match = 1;
        for (int i = 1; i < ra->n_allele; i++) {
            if (ra->als[i].s[0] == '*') continue;
            if (ra->als[i].s[ra->als[i].l - 1] != last_ref) { all_match = 0; break; }
        }
        if (!all_match) return;
        if (ra->beg == 0 || ra->pos - 1 < ra->beg) {
            ra->beg = ra->pos - ra->aln_win;
            free(ra->ref);
            ra->ref = NULL;
            scalar_safe_fetch_sequence(ra);
            if (!ra->ref) return;
        }
        char first_ref = ra->ref[ra->pos - ra->beg - 1];
        for (int i = 0; i < ra->n_allele; i++) {
            if (ra->als[i].s[0] == '*') continue;
            memmove(ra->als[i].s + 1, ra->als[i].s, ra->als[i].l);
            ra->als[i].s[0] = first_ref;
        }
        ra->pos--;
    }
}

/* ======================================================================== */

static liftover_bind_t *load_liftover_context(const char *chain_path, const char *dst_fasta_ref,
                                              const char *src_fasta_ref, int max_snp_gap, int max_indel_inc,
                                              char **errbuf_out) {
    liftover_bind_t *bind = NULL;
    htsFile *chain_fp = NULL;
    char *errbuf = NULL;

    bind = (liftover_bind_t *)calloc(1, sizeof(*bind));
    if (!bind) return NULL;
    bind->chain_path = dup_cstr(chain_path);
    bind->dst_fasta_ref = dup_cstr(dst_fasta_ref);
    bind->src_fasta_ref = dup_cstr(src_fasta_ref);
    bind->max_snp_gap = max_snp_gap;
    bind->max_indel_inc = max_indel_inc;

    chain_fp = hts_open(bind->chain_path, "r");
    if (!chain_fp) {
        if (errbuf_out) *errbuf_out = fmt_path_message("bcftools_liftover: failed to open chain file", bind->chain_path);
        liftover_bind_destroy(bind);
        return NULL;
    }
    bind->n_chains = read_chains(chain_fp, bind->max_snp_gap, &bind->chains, &bind->blocks, &errbuf);
    hts_close(chain_fp);
    if (bind->n_chains < 0) {
        if (errbuf_out) *errbuf_out = errbuf ? errbuf : dup_cstr("bcftools_liftover: failed to parse chain file");
        else free(errbuf);
        liftover_bind_destroy(bind);
        return NULL;
    }
    bind->idx = regidx_init_chains(bind->chains, bind->n_chains, bind->blocks);
    if (!bind->idx) {
        if (errbuf_out) *errbuf_out = dup_cstr("bcftools_liftover: failed to build chain interval index");
        liftover_bind_destroy(bind);
        return NULL;
    }
    bind->dst_fai = fai_load(bind->dst_fasta_ref);
    if (!bind->dst_fai) {
        if (errbuf_out) *errbuf_out = fmt_path_message("bcftools_liftover: failed to load destination FASTA index", bind->dst_fasta_ref);
        liftover_bind_destroy(bind);
        return NULL;
    }
    if (bind->src_fasta_ref) {
        bind->src_fai = fai_load(bind->src_fasta_ref);
        if (!bind->src_fai) {
            if (errbuf_out) *errbuf_out = fmt_path_message("bcftools_liftover: failed to load source FASTA index", bind->src_fasta_ref);
            liftover_bind_destroy(bind);
            return NULL;
        }
    }
    return bind;
}

static liftover_bind_t *get_liftover_context(const char *chain_path, const char *dst_fasta_ref,
                                             const char *src_fasta_ref, int max_snp_gap, int max_indel_inc,
                                             char **errbuf_out) {
    liftover_cache_entry_t *head = get_liftover_cache_head();
    liftover_cache_entry_t *entry = NULL;
    liftover_cache_entry_t *prev = NULL;
    liftover_bind_t *loaded = NULL;

    for (entry = head; entry; prev = entry, entry = entry->next) {
        liftover_bind_t *bind = entry->bind;
        if (strcmp(bind->chain_path, chain_path) == 0 &&
            strcmp(bind->dst_fasta_ref, dst_fasta_ref) == 0 &&
            ((bind->src_fasta_ref == NULL && src_fasta_ref == NULL) ||
             (bind->src_fasta_ref && src_fasta_ref && strcmp(bind->src_fasta_ref, src_fasta_ref) == 0)) &&
            bind->max_snp_gap == max_snp_gap &&
            bind->max_indel_inc == max_indel_inc) {
            if (prev) {
                prev->next = entry->next;
                entry->next = head;
                head = entry;
                set_liftover_cache_head(head);
            }
            return bind;
        }
    }

    loaded = load_liftover_context(chain_path, dst_fasta_ref, src_fasta_ref, max_snp_gap, max_indel_inc, errbuf_out);
    if (!loaded) return NULL;

    entry = (liftover_cache_entry_t *)calloc(1, sizeof(*entry));
    if (!entry) {
        liftover_bind_destroy(loaded);
        return NULL;
    }
    entry->bind = loaded;
    entry->next = head;
    head = entry;
    trim_liftover_cache_head(head);
    set_liftover_cache_head(entry);
    return loaded;
}

static void set_liftover_error(duckdb_function_info info, const char *msg) {
    duckdb_scalar_function_set_error(info, msg ? msg : "bcftools_liftover: unknown error");
}

static void bcftools_liftover_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector chrom_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector pos_vec = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector ref_vec = duckdb_data_chunk_get_vector(input, 2);
    duckdb_vector alt_vec = duckdb_data_chunk_get_vector(input, 3);
    duckdb_vector chain_vec = duckdb_data_chunk_get_vector(input, 4);
    duckdb_vector dst_fasta_vec = duckdb_data_chunk_get_vector(input, 5);
    duckdb_vector src_fasta_vec = duckdb_data_chunk_get_vector(input, 6);
    duckdb_vector max_snp_gap_vec = duckdb_data_chunk_get_vector(input, 7);
    duckdb_vector max_indel_inc_vec = duckdb_data_chunk_get_vector(input, 8);
    duckdb_vector lift_mt_vec = duckdb_data_chunk_get_vector(input, 9);
    duckdb_vector end_pos_vec = duckdb_data_chunk_get_vector(input, 10);
    idx_t n_input_cols = duckdb_data_chunk_get_column_count(input);
    duckdb_vector no_left_align_vec = (n_input_cols > 11) ? duckdb_data_chunk_get_vector(input, 11) : NULL;
    duckdb_vector child_vecs[OUT_COUNT];
    int64_t *src_pos_data;
    int64_t *dest_pos_data;
    int64_t *dest_end_data;
    bool *mapped_data;
    bool *revcomp_data;
    int32_t *swap_data;
    idx_t row_count = duckdb_data_chunk_get_size(input);

    for (int i = 0; i < OUT_COUNT; i++) child_vecs[i] = duckdb_struct_vector_get_child(output, (idx_t)i);
    src_pos_data = (int64_t *)duckdb_vector_get_data(child_vecs[OUT_SRC_POS]);
    dest_pos_data = (int64_t *)duckdb_vector_get_data(child_vecs[OUT_DEST_POS]);
    dest_end_data = (int64_t *)duckdb_vector_get_data(child_vecs[OUT_DEST_END]);
    mapped_data = (bool *)duckdb_vector_get_data(child_vecs[OUT_MAPPED]);
    revcomp_data = (bool *)duckdb_vector_get_data(child_vecs[OUT_REVERSE_COMPLEMENTED]);
    swap_data = (int32_t *)duckdb_vector_get_data(child_vecs[OUT_SWAP]);

    for (idx_t row = 0; row < row_count; row++) {
        idx_t chrom_len = 0, ref_len = 0, alt_len = 0;
        idx_t chain_len = 0, dst_fasta_len = 0, src_fasta_len = 0;
        const char *chrom;
        const char *ref = NULL;
        const char *alt = NULL;
        const char *chain_path = NULL;
        const char *dst_fasta_ref = NULL;
        const char *src_fasta_ref = NULL;
        char *chain_path_copy = NULL;
        char *dst_fasta_ref_copy = NULL;
        char *src_fasta_ref_copy = NULL;
        char *load_err = NULL;
        char reject_reason[LIFTOVER_WARNING_BUF] = {0};
        char note[LIFTOVER_WARNING_BUF] = {0};
        char *src_ref_copy = NULL, *src_alt_copy = NULL;
        char *lift_ref = NULL, *lift_alt = NULL, *dst_ref = NULL;
        const char *dst_chr = NULL;
        char *dst_chr_alloc = NULL; /* for MT passthrough: allocated dest contig name */
        hts_pos_t dst_pos5 = -1, dst_pos3 = -1;
        int swap = 0;
        int is_reverse = 0;
        int mapped = 0;
        int32_t max_snp_gap = 1;
        int32_t max_indel_inc = 250;
        int lift_mt = 0;
        int no_left_align = 0;
        int has_end_pos = 0;
        int64_t end_pos = 0;
        hts_pos_t dest_end = -1;
        liftover_bind_t *bind = NULL;

        if (!row_is_valid(chrom_vec, row)) {
            set_liftover_error(info, "bcftools_liftover: chrom must be non-null");
            return;
        }
        if (!row_is_valid(pos_vec, row)) {
            set_liftover_error(info, "bcftools_liftover: pos must be non-null");
            return;
        }
        if (!row_is_valid(chain_vec, row)) {
            set_liftover_error(info, "bcftools_liftover: chain_path must be non-null");
            return;
        }
        if (!row_is_valid(dst_fasta_vec, row)) {
            set_liftover_error(info, "bcftools_liftover: dst_fasta_ref must be non-null");
            return;
        }

        chrom = get_string_at(chrom_vec, row, &chrom_len);
        chain_path = get_string_at(chain_vec, row, &chain_len);
        dst_fasta_ref = get_string_at(dst_fasta_vec, row, &dst_fasta_len);
        if (row_is_valid(src_fasta_vec, row)) src_fasta_ref = get_string_at(src_fasta_vec, row, &src_fasta_len);
        if (row_is_valid(max_snp_gap_vec, row)) max_snp_gap = get_int32_at(max_snp_gap_vec, row);
        if (row_is_valid(max_indel_inc_vec, row)) max_indel_inc = get_int32_at(max_indel_inc_vec, row);
        if (row_is_valid(lift_mt_vec, row)) lift_mt = ((bool *)duckdb_vector_get_data(lift_mt_vec))[row] ? 1 : 0;
        if (no_left_align_vec && row_is_valid(no_left_align_vec, row)) no_left_align = ((bool *)duckdb_vector_get_data(no_left_align_vec))[row] ? 1 : 0;
        if (row_is_valid(end_pos_vec, row)) {
            end_pos = ((int64_t *)duckdb_vector_get_data(end_pos_vec))[row];
            has_end_pos = 1;
        }
        if (!chrom || chrom_len == 0) {
            set_liftover_error(info, "bcftools_liftover: chrom must be non-empty");
            return;
        }
        if (!chain_path || chain_len == 0) {
            set_liftover_error(info, "bcftools_liftover: chain_path must be non-empty");
            return;
        }
        if (!dst_fasta_ref || dst_fasta_len == 0) {
            set_liftover_error(info, "bcftools_liftover: dst_fasta_ref must be non-empty");
            return;
        }
        if (max_snp_gap < 0) {
            set_liftover_error(info, "bcftools_liftover: max_snp_gap must be >= 0");
            return;
        }
        if (max_indel_inc < 0) {
            set_liftover_error(info, "bcftools_liftover: max_indel_inc must be >= 0");
            return;
        }
        chain_path_copy = dup_span(chain_path, chain_len);
        dst_fasta_ref_copy = dup_span(dst_fasta_ref, dst_fasta_len);
        if (src_fasta_ref && src_fasta_len > 0) src_fasta_ref_copy = dup_span(src_fasta_ref, src_fasta_len);
        bind = get_liftover_context(chain_path_copy, dst_fasta_ref_copy, src_fasta_ref_copy, max_snp_gap, max_indel_inc,
                                    &load_err);
        free(chain_path_copy);
        free(dst_fasta_ref_copy);
        free(src_fasta_ref_copy);
        if (!bind) {
            char *full_err = fmt_message("bcftools_liftover: failed to load chain or FASTA context", load_err);
            set_liftover_error(info, full_err ? full_err : "bcftools_liftover: failed to load chain or FASTA context");
            free(load_err);
            free(full_err);
            return;
        }

        src_pos_data[row] = get_int64_at(pos_vec, row);
        if (src_pos_data[row] < 1) {
            set_liftover_error(info, "bcftools_liftover: pos must be >= 1");
            return;
        }
        duckdb_vector_assign_string_element_len(child_vecs[OUT_SRC_CHROM], row, chrom, chrom_len);
        if (row_is_valid(ref_vec, row)) {
            ref = get_string_at(ref_vec, row, &ref_len);
            if (ref_len > 0) {
                src_ref_copy = dup_span(ref, ref_len);
                seq_to_upper_ascii(src_ref_copy);
                duckdb_vector_assign_string_element_len(child_vecs[OUT_SRC_REF], row, src_ref_copy, ref_len);
            } else {
                set_null_at(child_vecs[OUT_SRC_REF], row);
            }
        } else {
            set_null_at(child_vecs[OUT_SRC_REF], row);
        }
        if (row_is_valid(alt_vec, row)) {
            alt = get_string_at(alt_vec, row, &alt_len);
            if (alt_len > 0) {
                src_alt_copy = dup_span(alt, alt_len);
                seq_to_upper_ascii(src_alt_copy);
                duckdb_vector_assign_string_element_len(child_vecs[OUT_SRC_ALT], row, src_alt_copy, alt_len);
            } else {
                set_null_at(child_vecs[OUT_SRC_ALT], row);
            }
        } else {
            set_null_at(child_vecs[OUT_SRC_ALT], row);
        }

        if (!src_ref_copy && src_alt_copy) tag_append(note, sizeof(note), "MissingSourceRef");

        /* Mitochondria handling (L-M1): upstream liftover.c:915-922,2398-2401.
         * When lift_mt is false (default), MT variants with matching source/dest
         * contig lengths are passed through with only contig rename — they share
         * the same MT assembly and do not need chain-based liftover. */
        {
            char *chrom_canon = canonical_contig_name(chrom);
            if (!lift_mt && is_mt_contig(chrom_canon)) {
                char *dst_mt_name = find_dst_mt_name(bind->dst_fai);
                if (dst_mt_name) {
                    int src_mt_sz = get_chain_mt_size(bind);
                    hts_pos_t dst_mt_len = faidx_seq_len64(bind->dst_fai, dst_mt_name);
                    if (src_mt_sz == 0 || src_mt_sz == (int)dst_mt_len) {
                        /* Pass through: rename contig, keep position */
                        mapped = 1;
                        tag_append(note, sizeof(note), "MitochondriaPassthrough");
                        dst_chr_alloc = dst_mt_name;
                        dst_chr = dst_chr_alloc;
                        dst_pos5 = src_pos_data[row];
                        lift_ref = src_ref_copy ? dup_cstr(src_ref_copy) : NULL;
                        lift_alt = src_alt_copy ? dup_cstr(src_alt_copy) : NULL;
                        if (has_end_pos) dest_end = end_pos;
                        free(chrom_canon);
                        goto emit_row;
                    }
                    free(dst_mt_name);
                    /* Sizes differ: fall through to normal chain liftover */
                }
            }
            free(chrom_canon);
        }

        if (!has_source_contig(bind, chrom)) {
            tag_append(reject_reason, sizeof(reject_reason), "MissingContig");
        } else {
            /* ================================================================
             * Full upstream pipeline (process() from liftover.c:2341-2514):
             *
             * Step 1: Try liftover_bp() for the 5' position of ALL variants.
             *         If SNP fails bp, set is_difficult_snp → fall to indel path.
             * Step 2: For non-SNPs or difficult SNPs: extend_alleles() if src_fai.
             * Step 3: For non-SNPs or difficult SNPs: liftover_indel().
             * Step 4: If npad && src_fai: fetch pad from source ref.
             * Step 5: Flip alleles if strand is reversed.
             * Step 6: If npad && dst_pos3 > dst_pos5: clip_pad() via NW.
             * Step 7: find_reference() against destination ref.
             * Step 8: Left-align if needed.
             * ================================================================ */
            int npad = 0;
            int ret = 0;
            int is_difficult_snp = 0;
            int multi_match = 0;
            char *pad_seq = NULL;

            /* Mutable copies for the allele extension / left-alignment pipeline.
             * These start as copies of src_ref/alt and get modified in-place. */
            char **work_alleles = NULL;
            int work_n_allele = 0;
            /* Working 0-based position (upstream rec->pos convention) */
            hts_pos_t work_pos = src_pos_data[row] - 1; /* convert 1-based → 0-based */

            int is_symbolic = 0;
            int is_snp = 0;

            if (src_ref_copy && src_alt_copy) {
                if (!make_allele_array(src_ref_copy, src_alt_copy, &work_alleles, &work_n_allele)) {
                    tag_append(reject_reason, sizeof(reject_reason), "OutOfMemory");
                }
            }
            if (work_alleles && work_n_allele > 0) {
                is_symbolic = scalar_is_symbolic_alleles(work_alleles, work_n_allele);
                is_snp = scalar_is_snp_alleles(work_alleles, work_n_allele);
                if (is_symbolic) tag_append(note, sizeof(note), "SymbolicAlleles");
            }

            /* Step 1: Try liftover_bp() for 5' anchor */
            {
                int block_ind = liftover_bp(bind, chrom, src_pos_data[row],
                                            &dst_chr, &dst_pos5, &is_reverse, &multi_match);
                dst_pos3 = dst_pos5;
                if (block_ind < 0) {
                    is_difficult_snp = 1;
                } else {
                    if (multi_match) tag_append(note, sizeof(note), "MultiBlock");
                }
            }

            /* Step 2: For non-SNPs or difficult SNPs, extend alleles if needed */
            if ((!is_snp || is_difficult_snp) && !is_symbolic && work_alleles && work_n_allele > 1) {
                if (bind->src_fai) {
                    if (scalar_is_extension_needed(work_alleles, work_n_allele)) {
                        scalar_realign_t ra = {0};
                        ra.chrom = chrom;
                        ra.fai = bind->src_fai;
                        ra.pos = work_pos;
                        ra.aln_win = SCALAR_ALN_WIN;
                        scalar_extend_alleles(&ra, work_alleles, work_n_allele);
                        if (ra.ref) {
                            /* Update working alleles and position from extension result */
                            char **new_alleles = NULL;
                            int new_n_allele = 0;
                            work_pos = scalar_update_alleles(&ra, &new_alleles, &new_n_allele);
                            if (work_pos >= 0) {
                                free_allele_array(work_alleles, work_n_allele);
                                work_alleles = new_alleles;
                                work_n_allele = new_n_allele;
                            }
                        }
                        scalar_realign_cleanup(&ra);
                    }
                }
            }

            /* Step 3: For non-SNPs or difficult SNPs, call liftover_indel() */
            if ((!is_snp || is_difficult_snp) && !is_symbolic && work_alleles && work_n_allele > 0) {
                hts_pos_t src_pos5 = work_pos + 1; /* 1-based */
                hts_pos_t src_pos3 = work_pos + (hts_pos_t)strlen(work_alleles[0]); /* 1-based */
                ret = liftover_indel(bind, chrom, src_pos5, src_pos3,
                                     NULL, &dst_chr, &dst_pos5, &dst_pos3, &is_reverse, &npad);
                if (ret < 0) {
                    tag_append(reject_reason, sizeof(reject_reason),
                               ret == -1 ? "UnmappedAnchors" :
                               ret == -2 ? "UnmappedAnchor5" :
                               ret == -3 ? "UnmappedAnchor3" :
                               ret == -4 ? "MismatchAnchors" : "ApartAnchors");
                } else if (npad != 0 && !bind->src_fai) {
                    tag_append(reject_reason, sizeof(reject_reason), "MissingFasta");
                    ret = -1; /* mark as rejected */
                } else {
                    mapped = 1;
                    if (npad != 0) tag_append(note, sizeof(note), "Padded");
                }
            } else if (!is_difficult_snp && !reject_reason[0]) {
                /* SNP that succeeded bp in Step 1 */
                if (dst_pos5 > 0) {
                    mapped = 1;
                } else {
                    tag_append(reject_reason, sizeof(reject_reason), "UnmappedAnchors");
                }
            }

            /* Step 4: Collect pad sequence from source reference if needed */
            if (mapped && npad != 0 && bind->src_fai && work_alleles && work_n_allele > 0) {
                hts_pos_t src_pos5_1 = work_pos + 1; /* 1-based */
                hts_pos_t src_pos3_1 = work_pos + (hts_pos_t)strlen(work_alleles[0]); /* 1-based */
                hts_pos_t npad_pos5 = npad < 0 ? src_pos5_1 + npad : src_pos3_1 + 1;
                hts_pos_t npad_pos3 = npad < 0 ? src_pos5_1 - 1 : src_pos3_1 + npad;
                pad_seq = fetch_sequence_flexible(bind->src_fai, chrom, npad_pos5, npad_pos3);
                if (!pad_seq) {
                    tag_append(note, sizeof(note), "MissingPadSequence");
                }
            }

            /* Step 5: Flip alleles if strand is reversed */
            if (mapped && is_reverse) {
                for (int i = 0; i < work_n_allele; i++) {
                    if (!work_alleles[i]) continue;
                    if (work_alleles[i][0] != '<' && work_alleles[i][0] != '*' && work_alleles[i][0] != ',') {
                        char *rc = reverse_complement_copy(work_alleles[i]);
                        free(work_alleles[i]);
                        work_alleles[i] = rc;
                    }
                }
                if (pad_seq) {
                    /* Reverse complement pad and negate npad (exact upstream behavior) */
                    char *rc = reverse_complement_copy(pad_seq);
                    free(pad_seq);
                    pad_seq = rc;
                    npad = -npad;
                }
            }

            /* Step 6: If padded and dst_pos3 > dst_pos5, clip pad via NW re-alignment */
            if (mapped && npad != 0 && pad_seq && dst_pos3 > dst_pos5 && bind->dst_fai) {
                char *dst_ref_nw = fetch_sequence_flexible(bind->dst_fai, dst_chr, dst_pos5, dst_pos3);
                if (dst_ref_nw) {
                    scalar_clip_pad(work_alleles, work_n_allele, dst_ref_nw, pad_seq, npad, &dst_pos5, &dst_pos3);
                    free(dst_ref_nw);
                }
            }

            /* Step 7: find_reference() against destination reference */
            if (mapped && !is_symbolic && bind->dst_fai && work_alleles && work_n_allele > 0) {
                /* Update working position to destination 0-based */
                hts_pos_t dst_work_pos = dst_pos5 - 1; /* 1-based → 0-based */
                dst_ref = fetch_sequence_flexible(bind->dst_fai, dst_chr, dst_pos5, dst_pos3);
                if (!dst_ref) {
                    tag_append(note, sizeof(note), "MissingDestinationRef");
                } else {
                    swap = scalar_find_reference(&work_alleles, &work_n_allele, dst_ref, is_snp, &dst_work_pos);
                    /* Update dst_pos5 if find_reference adjusted position (right-match SNP rescue) */
                    dst_pos5 = dst_work_pos + 1; /* 0-based → 1-based */
                    free(dst_ref);
                    dst_ref = NULL;
                }
            }

            /* Step 8: Left-align indels if needed */
            if (mapped && (!is_snp || is_difficult_snp) && !is_symbolic
                && work_alleles && work_n_allele > 0
                && !no_left_align
                && !scalar_is_left_aligned(work_alleles, work_n_allele)
                && bind->dst_fai) {
                scalar_realign_t ra = {0};
                ra.chrom = dst_chr;
                ra.fai = bind->dst_fai;
                ra.pos = dst_pos5 - 1; /* 1-based → 0-based */
                ra.aln_win = SCALAR_ALN_WIN;
                if (scalar_initialize_alleles(&ra, work_alleles, work_n_allele) && ra.ref) {
                    char **new_alleles = NULL;
                    int new_n_allele = 0;
                    scalar_trim_right(&ra);
                    scalar_trim_left(&ra);
                    scalar_shift_left(&ra);
                    hts_pos_t new_pos = scalar_update_alleles(&ra, &new_alleles, &new_n_allele);
                    if (new_pos >= 0) {
                        dst_pos5 = new_pos + 1; /* 0-based → 1-based */
                        free_allele_array(work_alleles, work_n_allele);
                        work_alleles = new_alleles;
                        work_n_allele = new_n_allele;
                    }
                }
                scalar_realign_cleanup(&ra);
            }

            /* Normalize alleles to upstream VCF semantics before emitting */
            if (work_alleles && work_n_allele > 0) {
                scalar_normalize_alleles(&work_alleles, &work_n_allele);
            }

            /* Assign outputs */
            if (work_alleles && work_n_allele > 0) {
                lift_ref = dup_cstr(work_alleles[0]);
                lift_alt = join_alt_csv(work_alleles, work_n_allele);
            }
            free_allele_array(work_alleles, work_n_allele);
            free(pad_seq);

            /* Lift INFO/END if provided (L-M2): upstream liftover.c:2516-2534.
             * Lift end_pos through liftover_bp(); validate same contig, same strand,
             * and ordering (forward: end >= pos; reverse: end <= pos). */
            if (has_end_pos && mapped) {
                const char *end_dst_chr = NULL;
                hts_pos_t end_dst_pos = -1;
                int end_is_reverse = 0;
                int end_multi = 0;
                int end_ret = liftover_bp(bind, chrom, end_pos,
                                          &end_dst_chr, &end_dst_pos, &end_is_reverse, &end_multi);
                if (end_ret >= 0 && end_dst_chr && dst_chr &&
                    strcmp(end_dst_chr, dst_chr) == 0 && end_is_reverse == is_reverse) {
                    /* Validate ordering: forward strand → end >= pos; reverse → end <= pos */
                    if ((!is_reverse && end_dst_pos >= dst_pos5) ||
                        (is_reverse && end_dst_pos <= dst_pos5)) {
                        dest_end = end_dst_pos;
                    } else {
                        dest_end = -1; /* ordering violation */
                        tag_append(note, sizeof(note), "EndLiftOrderViolation");
                    }
                } else {
                    dest_end = -1; /* lift failed or different contig/strand */
                    tag_append(note, sizeof(note), "EndLiftFailed");
                }
            }
        }

    emit_row:
        mapped_data[row] = (bool)mapped;
        revcomp_data[row] = (bool)is_reverse;
        swap_data[row] = swap;

        if (!mapped) {
            set_null_at(child_vecs[OUT_DEST_CHROM], row);
            set_null_at(child_vecs[OUT_DEST_POS], row);
            set_null_at(child_vecs[OUT_DEST_END], row);
            set_null_at(child_vecs[OUT_DEST_REF], row);
            set_null_at(child_vecs[OUT_DEST_ALT], row);
            if (reject_reason[0]) duckdb_vector_assign_string_element(child_vecs[OUT_REJECT_REASON], row, reject_reason);
            else set_null_at(child_vecs[OUT_REJECT_REASON], row);
            if (note[0]) duckdb_vector_assign_string_element(child_vecs[OUT_NOTE], row, note);
            else set_null_at(child_vecs[OUT_NOTE], row);
            free(src_ref_copy);
            free(src_alt_copy);
            free(lift_ref);
            free(lift_alt);
            continue;
        }

        duckdb_vector_assign_string_element(child_vecs[OUT_DEST_CHROM], row, dst_chr);
        dest_pos_data[row] = dst_pos5;

        if (has_end_pos && dest_end >= 0) dest_end_data[row] = dest_end;
        else set_null_at(child_vecs[OUT_DEST_END], row);

        if (lift_ref) duckdb_vector_assign_string_element(child_vecs[OUT_DEST_REF], row, lift_ref);
        else set_null_at(child_vecs[OUT_DEST_REF], row);
        if (lift_alt) duckdb_vector_assign_string_element(child_vecs[OUT_DEST_ALT], row, lift_alt);
        else set_null_at(child_vecs[OUT_DEST_ALT], row);
        set_null_at(child_vecs[OUT_REJECT_REASON], row);
        if (note[0]) duckdb_vector_assign_string_element(child_vecs[OUT_NOTE], row, note);
        else set_null_at(child_vecs[OUT_NOTE], row);

        free(src_ref_copy);
        free(src_alt_copy);
        free(lift_ref);
        free(lift_alt);
        free(dst_ref);
        free(dst_chr_alloc);
    }
}

/* Build a scalar function object for one bcftools_liftover overload.
   include_no_left_align=0 → 11-param legacy signature (no_left_align absent,
   defaults to false at runtime); =1 → 12-param signature with explicit no_left_align. */
static duckdb_scalar_function build_liftover_overload(duckdb_logical_type varchar_type,
                                                       duckdb_logical_type bigint_type,
                                                       duckdb_logical_type int_type,
                                                       duckdb_logical_type bool_type,
                                                       duckdb_logical_type struct_type,
                                                       int include_no_left_align) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_scalar_function_add_parameter(fn, varchar_type);   /* chrom */
    duckdb_scalar_function_add_parameter(fn, bigint_type);    /* pos */
    duckdb_scalar_function_add_parameter(fn, varchar_type);   /* ref */
    duckdb_scalar_function_add_parameter(fn, varchar_type);   /* alt */
    duckdb_scalar_function_add_parameter(fn, varchar_type);   /* chain_path */
    duckdb_scalar_function_add_parameter(fn, varchar_type);   /* dst_fasta_ref */
    duckdb_scalar_function_add_parameter(fn, varchar_type);   /* src_fasta_ref */
    duckdb_scalar_function_add_parameter(fn, int_type);       /* max_snp_gap */
    duckdb_scalar_function_add_parameter(fn, int_type);       /* max_indel_inc */
    duckdb_scalar_function_add_parameter(fn, bool_type);      /* lift_mt */
    duckdb_scalar_function_add_parameter(fn, bigint_type);    /* end_pos */
    if (include_no_left_align)
        duckdb_scalar_function_add_parameter(fn, bool_type);  /* no_left_align */
    duckdb_scalar_function_set_return_type(fn, struct_type);
    duckdb_scalar_function_set_special_handling(fn);
    duckdb_scalar_function_set_function(fn, bcftools_liftover_scalar);
    return fn;
}

void register_liftover_functions(duckdb_connection connection) {
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type int_type = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type fields[OUT_COUNT];
    duckdb_logical_type struct_type;
    duckdb_scalar_function fn12;

    fields[OUT_SRC_CHROM] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[OUT_SRC_POS] = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    fields[OUT_SRC_REF] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[OUT_SRC_ALT] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[OUT_DEST_CHROM] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[OUT_DEST_POS] = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    fields[OUT_DEST_END] = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    fields[OUT_DEST_REF] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[OUT_DEST_ALT] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[OUT_MAPPED] = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    fields[OUT_REVERSE_COMPLEMENTED] = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    fields[OUT_SWAP] = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    fields[OUT_REJECT_REASON] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[OUT_NOTE] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    struct_type = duckdb_create_struct_type(fields, LIFTOVER_FIELD_NAMES, OUT_COUNT);

    /* Register only the 12-param overload; no_left_align defaults to false in
       the duckdb_liftover SQL macro wrapper for callers that don't need it. */
    fn12 = build_liftover_overload(varchar_type, bigint_type, int_type, bool_type, struct_type, 1);
    duckdb_scalar_function_set_name(fn12, "bcftools_liftover");
    duckdb_register_scalar_function(connection, fn12);

    duckdb_destroy_scalar_function(&fn12);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&int_type);
    duckdb_destroy_logical_type(&bool_type);
    for (int i = 0; i < OUT_COUNT; i++) duckdb_destroy_logical_type(&fields[i]);
    duckdb_destroy_logical_type(&struct_type);
}
