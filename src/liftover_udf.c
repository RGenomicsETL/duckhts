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
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <htslib/regidx.h>

#define LIFTOVER_WARNING_BUF 256
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

static liftover_cache_entry_t *g_liftover_cache_head = NULL;
static pthread_mutex_t g_liftover_cache_mutex = PTHREAD_MUTEX_INITIALIZER;

enum {
    OUT_SRC_CHROM = 0,
    OUT_SRC_POS,
    OUT_SRC_REF,
    OUT_SRC_ALT,
    OUT_DEST_CHROM,
    OUT_DEST_POS,
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
    "dest_chrom", "dest_pos", "dest_ref", "dest_alt",
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
    if (strcasecmp(s, "M") == 0 || strcasecmp(s, "MT") == 0) return dup_cstr("MT");
    if (strcmp(s, "23") == 0 || strcasecmp(s, "X") == 0 || strcasecmp(s, "XY") == 0 ||
        strcasecmp(s, "XX") == 0 || strcasecmp(s, "PAR1") == 0 || strcasecmp(s, "PAR2") == 0) {
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

static char dna_complement(char base) {
    switch (toupper((unsigned char)base)) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'N': return 'N';
        default: return 'N';
    }
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

static regidx_t *regidx_init_chains(const lo_chain_t *chains, int n_chains, const lo_block_t *blocks) {
    regidx_t *idx = regidx_init(NULL, NULL, NULL, sizeof(int), NULL);
    if (!idx) return NULL;
    for (int i = 0; i < n_chains; i++) {
        const lo_chain_t *chain = &chains[i];
        for (int j = 0; j < chain->n_blocks; j++) {
            int block_ind = chain->block_ind + j;
            const lo_block_t *block = &blocks[block_ind];
            regidx_push(idx, (char *)chain->t_name, (char *)chain->t_name + strlen(chain->t_name),
                        chain->t_start + block->t_start + 1,
                        chain->t_start + block->t_start + block->size, (void *)&block_ind);
        }
    }
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

static liftover_bind_t *liftover_bind_copy_data(const liftover_bind_t *src) {
    liftover_bind_t *bind;
    htsFile *chain_fp = NULL;
    char *errbuf = NULL;
    if (!src) return NULL;
    bind = (liftover_bind_t *)calloc(1, sizeof(*bind));
    if (!bind) return NULL;
    bind->chain_path = dup_cstr(src->chain_path);
    bind->dst_fasta_ref = dup_cstr(src->dst_fasta_ref);
    bind->src_fasta_ref = dup_cstr(src->src_fasta_ref);
    bind->max_snp_gap = src->max_snp_gap;
    bind->max_indel_inc = src->max_indel_inc;

    chain_fp = hts_open(bind->chain_path, "r");
    if (!chain_fp) {
        liftover_bind_destroy(bind);
        return NULL;
    }
    bind->n_chains = read_chains(chain_fp, bind->max_snp_gap, &bind->chains, &bind->blocks, &errbuf);
    hts_close(chain_fp);
    if (bind->n_chains < 0) {
        free(errbuf);
        liftover_bind_destroy(bind);
        return NULL;
    }
    bind->idx = regidx_init_chains(bind->chains, bind->n_chains, bind->blocks);
    bind->dst_fai = fai_load(bind->dst_fasta_ref);
    if (!bind->dst_fai) {
        liftover_bind_destroy(bind);
        return NULL;
    }
    if (bind->src_fasta_ref) {
        bind->src_fai = fai_load(bind->src_fasta_ref);
        if (!bind->src_fai) {
            liftover_bind_destroy(bind);
            return NULL;
        }
    }
    return bind;
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
            block_ind = regitr_payload(itr, int);
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
    liftover_cache_entry_t *entry = NULL;
    liftover_bind_t *loaded = NULL;

    pthread_mutex_lock(&g_liftover_cache_mutex);
    for (entry = g_liftover_cache_head; entry; entry = entry->next) {
        liftover_bind_t *bind = entry->bind;
        if (strcmp(bind->chain_path, chain_path) == 0 &&
            strcmp(bind->dst_fasta_ref, dst_fasta_ref) == 0 &&
            ((bind->src_fasta_ref == NULL && src_fasta_ref == NULL) ||
             (bind->src_fasta_ref && src_fasta_ref && strcmp(bind->src_fasta_ref, src_fasta_ref) == 0)) &&
            bind->max_snp_gap == max_snp_gap &&
            bind->max_indel_inc == max_indel_inc) {
            pthread_mutex_unlock(&g_liftover_cache_mutex);
            return bind;
        }
    }
    pthread_mutex_unlock(&g_liftover_cache_mutex);

    loaded = load_liftover_context(chain_path, dst_fasta_ref, src_fasta_ref, max_snp_gap, max_indel_inc, errbuf_out);
    if (!loaded) return NULL;

    pthread_mutex_lock(&g_liftover_cache_mutex);
    for (entry = g_liftover_cache_head; entry; entry = entry->next) {
        liftover_bind_t *bind = entry->bind;
        if (strcmp(bind->chain_path, chain_path) == 0 &&
            strcmp(bind->dst_fasta_ref, dst_fasta_ref) == 0 &&
            ((bind->src_fasta_ref == NULL && src_fasta_ref == NULL) ||
             (bind->src_fasta_ref && src_fasta_ref && strcmp(bind->src_fasta_ref, src_fasta_ref) == 0)) &&
            bind->max_snp_gap == max_snp_gap &&
            bind->max_indel_inc == max_indel_inc) {
            pthread_mutex_unlock(&g_liftover_cache_mutex);
            liftover_bind_destroy(loaded);
            return bind;
        }
    }

    entry = (liftover_cache_entry_t *)calloc(1, sizeof(*entry));
    if (!entry) {
        pthread_mutex_unlock(&g_liftover_cache_mutex);
        liftover_bind_destroy(loaded);
        return NULL;
    }
    entry->bind = loaded;
    entry->next = g_liftover_cache_head;
    g_liftover_cache_head = entry;
    pthread_mutex_unlock(&g_liftover_cache_mutex);
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
    duckdb_vector child_vecs[OUT_COUNT];
    int64_t *src_pos_data;
    int64_t *dest_pos_data;
    bool *mapped_data;
    bool *revcomp_data;
    int32_t *swap_data;
    idx_t row_count = duckdb_data_chunk_get_size(input);

    for (int i = 0; i < OUT_COUNT; i++) child_vecs[i] = duckdb_struct_vector_get_child(output, (idx_t)i);
    src_pos_data = (int64_t *)duckdb_vector_get_data(child_vecs[OUT_SRC_POS]);
    dest_pos_data = (int64_t *)duckdb_vector_get_data(child_vecs[OUT_DEST_POS]);
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
        hts_pos_t dst_pos5 = -1, dst_pos3 = -1;
        int swap = 0;
        int is_reverse = 0;
        int mapped = 0;
        int32_t max_snp_gap = 1;
        int32_t max_indel_inc = 250;
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
        if (row_is_valid(max_snp_gap_vec, row)) max_snp_gap = (int32_t)get_int64_at(max_snp_gap_vec, row);
        if (row_is_valid(max_indel_inc_vec, row)) max_indel_inc = (int32_t)get_int64_at(max_indel_inc_vec, row);
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
        if (src_alt_copy && (strchr(src_alt_copy, '<') || strchr(src_alt_copy, '*') || strchr(src_alt_copy, ','))) {
            tag_append(note, sizeof(note), "SymbolicAlleles");
        }

        if (!has_source_contig(bind, chrom)) {
            tag_append(reject_reason, sizeof(reject_reason), "MissingContig");
        } else if (src_ref_copy && strlen(src_ref_copy) > 1) {
            int npad = 0;
            int ret = liftover_indel(bind, chrom, src_pos_data[row],
                                     src_pos_data[row] + (hts_pos_t)strlen(src_ref_copy) - 1,
                                     NULL, &dst_chr, &dst_pos5, &dst_pos3, &is_reverse, &npad);
            if (ret < 0) {
                tag_append(reject_reason, sizeof(reject_reason),
                           ret == -1 ? "UnmappedAnchors" :
                           ret == -2 ? "UnmappedAnchor5" :
                           ret == -3 ? "UnmappedAnchor3" :
                           ret == -4 ? "MismatchAnchors" : "ApartAnchors");
            } else {
                if (npad != 0 && !bind->src_fai) {
                    tag_append(reject_reason, sizeof(reject_reason), "MissingFasta");
                } else {
                    mapped = 1;
                    if (npad != 0) tag_append(note, sizeof(note), "Padded");
                }
            }
        } else {
            int multi_match = 0;
            int block = liftover_bp(bind, chrom, src_pos_data[row], &dst_chr, &dst_pos5, &is_reverse, &multi_match);
            dst_pos3 = dst_pos5;
            if (block >= 0) {
                mapped = 1;
                if (multi_match) tag_append(note, sizeof(note), "MultiBlock");
            } else {
                tag_append(reject_reason, sizeof(reject_reason), "UnmappedAnchors");
            }
        }

        mapped_data[row] = (bool)mapped;
        revcomp_data[row] = (bool)is_reverse;
        swap_data[row] = swap;

        if (!mapped) {
            set_null_at(child_vecs[OUT_DEST_CHROM], row);
            set_null_at(child_vecs[OUT_DEST_POS], row);
            set_null_at(child_vecs[OUT_DEST_REF], row);
            set_null_at(child_vecs[OUT_DEST_ALT], row);
            if (reject_reason[0]) duckdb_vector_assign_string_element(child_vecs[OUT_REJECT_REASON], row, reject_reason);
            else set_null_at(child_vecs[OUT_REJECT_REASON], row);
            if (note[0]) duckdb_vector_assign_string_element(child_vecs[OUT_NOTE], row, note);
            else set_null_at(child_vecs[OUT_NOTE], row);
            free(src_ref_copy);
            free(src_alt_copy);
            continue;
        }

        duckdb_vector_assign_string_element(child_vecs[OUT_DEST_CHROM], row, dst_chr);
        dest_pos_data[row] = dst_pos5;

        if (src_ref_copy) {
            lift_ref = is_reverse ? reverse_complement_copy(src_ref_copy) : dup_cstr(src_ref_copy);
        }
        if (src_alt_copy) {
            lift_alt = is_reverse ? reverse_complement_copy(src_alt_copy) : dup_cstr(src_alt_copy);
        }

        if (dst_chr && bind->dst_fai) {
            dst_ref = fetch_sequence_flexible(bind->dst_fai, dst_chr, dst_pos5, dst_pos3);
        }
        if (!dst_ref) {
            tag_append(note, sizeof(note), "MissingDestinationRef");
        } else {
            if (lift_ref && strcmp(dst_ref, lift_ref) == 0) {
                free(lift_ref);
                lift_ref = dup_cstr(dst_ref);
            } else if (lift_alt && strcmp(dst_ref, lift_alt) == 0) {
                char *tmp = lift_ref;
                lift_ref = dup_cstr(dst_ref);
                lift_alt = tmp;
                swap = 1;
            } else if (lift_ref) {
                tag_append(note, sizeof(note), "RefMismatch");
                free(lift_ref);
                lift_ref = dup_cstr(dst_ref);
            } else {
                lift_ref = dup_cstr(dst_ref);
            }
        }
        swap_data[row] = swap;

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
    }
}

void register_liftover_functions(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type int_type = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type fields[OUT_COUNT];
    duckdb_logical_type struct_type;

    duckdb_scalar_function_set_name(fn, "bcftools_liftover");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, int_type);
    duckdb_scalar_function_add_parameter(fn, int_type);

    fields[OUT_SRC_CHROM] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[OUT_SRC_POS] = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    fields[OUT_SRC_REF] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[OUT_SRC_ALT] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[OUT_DEST_CHROM] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[OUT_DEST_POS] = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    fields[OUT_DEST_REF] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[OUT_DEST_ALT] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[OUT_MAPPED] = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    fields[OUT_REVERSE_COMPLEMENTED] = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    fields[OUT_SWAP] = duckdb_create_logical_type(DUCKDB_TYPE_INTEGER);
    fields[OUT_REJECT_REASON] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[OUT_NOTE] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    struct_type = duckdb_create_struct_type(fields, LIFTOVER_FIELD_NAMES, OUT_COUNT);

    duckdb_scalar_function_set_return_type(fn, struct_type);
    duckdb_scalar_function_set_special_handling(fn);
    duckdb_scalar_function_set_function(fn, bcftools_liftover_scalar);
    duckdb_register_scalar_function(connection, fn);

    duckdb_destroy_scalar_function(&fn);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&int_type);
    duckdb_destroy_logical_type(&bool_type);
    for (int i = 0; i < OUT_COUNT; i++) duckdb_destroy_logical_type(&fields[i]);
    duckdb_destroy_logical_type(&struct_type);
}
