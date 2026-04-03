#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <ctype.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/faidx.h>

#define MUNGE_TAG_BUF 128

typedef struct munge_cache_entry {
    char *fasta_path;
    faidx_t *fai;
    struct munge_cache_entry *next;
} munge_cache_entry_t;

/*
 * faidx_t owns a mutable BGZF cursor, so a single handle cannot be shared
 * across concurrent DuckDB worker threads. Keep one cache per thread instead.
 */
static pthread_key_t g_munge_cache_key;
static pthread_once_t g_munge_cache_key_once = PTHREAD_ONCE_INIT;
static pthread_mutex_t g_munge_fai_fetch_mutex = PTHREAD_MUTEX_INITIALIZER;

static void destroy_munge_cache(void *ptr) {
    munge_cache_entry_t *entry = (munge_cache_entry_t *)ptr;
    while (entry) {
        munge_cache_entry_t *next = entry->next;
        free(entry->fasta_path);
        if (entry->fai) fai_destroy(entry->fai);
        free(entry);
        entry = next;
    }
}

static void init_munge_cache_key(void) {
    pthread_key_create(&g_munge_cache_key, destroy_munge_cache);
}

static munge_cache_entry_t *get_munge_cache_head(void) {
    pthread_once(&g_munge_cache_key_once, init_munge_cache_key);
    return (munge_cache_entry_t *)pthread_getspecific(g_munge_cache_key);
}

static void set_munge_cache_head(munge_cache_entry_t *head) {
    pthread_once(&g_munge_cache_key_once, init_munge_cache_key);
    pthread_setspecific(g_munge_cache_key, head);
}

enum {
    MUNGE_OUT_CHROM = 0,
    MUNGE_OUT_POS,
    MUNGE_OUT_ID,
    MUNGE_OUT_REF,
    MUNGE_OUT_ALT,
    MUNGE_OUT_ALLELES_SWAPPED,
    MUNGE_OUT_FILTER,
    MUNGE_OUT_NS,
    MUNGE_OUT_EZ,
    MUNGE_OUT_NC,
    MUNGE_OUT_ES,
    MUNGE_OUT_SE,
    MUNGE_OUT_LP,
    MUNGE_OUT_AF,
    MUNGE_OUT_AC,
    MUNGE_OUT_NE,
    MUNGE_OUT_SI,
    MUNGE_OUT_I2,
    MUNGE_OUT_CQ,
    MUNGE_OUT_ED,
    MUNGE_OUT_COUNT
};

static const char *MUNGE_FIELD_NAMES[MUNGE_OUT_COUNT] = {
    "chrom", "pos", "id", "ref", "alt", "alleles_swapped", "filter",
    "ns", "ez", "nc", "es", "se", "lp", "af", "ac", "ne",
    "si", "i2", "cq", "ed"
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

static inline double get_double_at(duckdb_vector vector, idx_t row) {
    return ((double *)duckdb_vector_get_data(vector))[row];
}

static char *dup_cstr(const char *s) {
    size_t n;
    char *out;
    if (!s) return NULL;
    n = strlen(s);
    out = (char *)malloc(n + 1);
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

static void set_munge_error(duckdb_function_info info, const char *msg) {
    duckdb_scalar_function_set_error(info, msg ? msg : "bcftools_munge_row: unknown error");
}

static void seq_to_upper_ascii(char *seq) {
    if (!seq) return;
    while (*seq) {
        *seq = (char)toupper((unsigned char)*seq);
        seq++;
    }
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
    if (strcasecmp(chrom, "M") == 0 || strcasecmp(chrom, "MT") == 0 || strcasecmp(chrom, "chrM") == 0) {
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

static char *normalize_allele(const char *s, size_t n) {
    int symbolic = 0;
    char *out;
    if (!s || n == 0) return NULL;
    if (s[0] == '<' || s[0] == '*') return dup_span(s, n);
    if (s[0] == 'D' || s[0] == 'I' || s[0] == 'd' || s[0] == 'i') symbolic = 1;
    for (size_t i = 0; i < n; i++) {
        if (s[i] == '+') symbolic = 1;
    }
    if (symbolic) {
        out = (char *)malloc(n + 3);
        if (!out) return NULL;
        out[0] = '<';
        memcpy(out + 1, s, n);
        out[n + 1] = '>';
        out[n + 2] = '\0';
        return out;
    }
    if (n == 1 && ((s[0] - '1') & 0xFC) == 0) {
        out = (char *)malloc(2);
        if (!out) return NULL;
        out[0] = "ACGT"[s[0] - '1'];
        out[1] = '\0';
        return out;
    }
    out = dup_span(s, n);
    seq_to_upper_ascii(out);
    return out;
}

static faidx_t *get_munge_fai(const char *fasta_path, char **err) {
    munge_cache_entry_t *head = get_munge_cache_head();
    munge_cache_entry_t *entry;
    faidx_t *loaded;
    char msg[1024];
    for (entry = head; entry; entry = entry->next) {
        if (strcmp(entry->fasta_path, fasta_path) == 0) {
            return entry->fai;
        }
    }

    loaded = fai_load3(fasta_path, NULL, NULL, FAI_CREATE);
    if (!loaded) {
        snprintf(msg, sizeof(msg), "bcftools_munge_row: failed to load FASTA index: %s", fasta_path);
        *err = dup_cstr(msg);
        return NULL;
    }

    entry = (munge_cache_entry_t *)calloc(1, sizeof(*entry));
    if (!entry) {
        fai_destroy(loaded);
        *err = dup_cstr("bcftools_munge_row: out of memory");
        return NULL;
    }
    entry->fasta_path = dup_cstr(fasta_path);
    if (!entry->fasta_path) {
        fai_destroy(loaded);
        free(entry);
        *err = dup_cstr("bcftools_munge_row: out of memory");
        return NULL;
    }
    entry->fai = loaded;

    entry->next = head;
    set_munge_cache_head(entry);
    return loaded;
}

static double safe_minus_log10(double x) {
    if (isnan(x) || x <= 0.0) return NAN;
    return -log10(x);
}

static double safe_log(double x) {
    if (isnan(x) || x <= 0.0) return NAN;
    return log(x);
}

static void assign_double_or_null(duckdb_vector vector, idx_t row, double value) {
    if (isnan(value)) set_null_at(vector, row);
    else ((double *)duckdb_vector_get_data(vector))[row] = value;
}

static void bcftools_munge_row_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector chrom_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector pos_vec = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector a1_vec = duckdb_data_chunk_get_vector(input, 2);
    duckdb_vector a2_vec = duckdb_data_chunk_get_vector(input, 3);
    duckdb_vector id_vec = duckdb_data_chunk_get_vector(input, 4);
    duckdb_vector p_vec = duckdb_data_chunk_get_vector(input, 5);
    duckdb_vector z_vec = duckdb_data_chunk_get_vector(input, 6);
    duckdb_vector or_vec = duckdb_data_chunk_get_vector(input, 7);
    duckdb_vector beta_vec = duckdb_data_chunk_get_vector(input, 8);
    duckdb_vector n_vec = duckdb_data_chunk_get_vector(input, 9);
    duckdb_vector n_cas_vec = duckdb_data_chunk_get_vector(input, 10);
    duckdb_vector n_con_vec = duckdb_data_chunk_get_vector(input, 11);
    duckdb_vector info_vec = duckdb_data_chunk_get_vector(input, 12);
    duckdb_vector frq_vec = duckdb_data_chunk_get_vector(input, 13);
    duckdb_vector se_vec = duckdb_data_chunk_get_vector(input, 14);
    duckdb_vector lp_vec = duckdb_data_chunk_get_vector(input, 15);
    duckdb_vector ac_vec = duckdb_data_chunk_get_vector(input, 16);
    duckdb_vector neff_vec = duckdb_data_chunk_get_vector(input, 17);
    duckdb_vector neffdiv2_vec = duckdb_data_chunk_get_vector(input, 18);
    duckdb_vector het_i2_vec = duckdb_data_chunk_get_vector(input, 19);
    duckdb_vector het_p_vec = duckdb_data_chunk_get_vector(input, 20);
    duckdb_vector het_lp_vec = duckdb_data_chunk_get_vector(input, 21);
    duckdb_vector dire_vec = duckdb_data_chunk_get_vector(input, 22);
    duckdb_vector fasta_ref_vec = duckdb_data_chunk_get_vector(input, 23);
    duckdb_vector iffy_tag_vec = duckdb_data_chunk_get_vector(input, 24);
    duckdb_vector mismatch_tag_vec = duckdb_data_chunk_get_vector(input, 25);
    duckdb_vector ns_override_vec = duckdb_data_chunk_get_vector(input, 26);
    duckdb_vector nc_override_vec = duckdb_data_chunk_get_vector(input, 27);
    duckdb_vector ne_override_vec = duckdb_data_chunk_get_vector(input, 28);
    duckdb_vector child_vecs[MUNGE_OUT_COUNT];
    idx_t row_count = duckdb_data_chunk_get_size(input);
    char *cached_fasta_path = NULL;
    idx_t cached_fasta_len = 0;
    faidx_t *cached_fai = NULL;

    for (int i = 0; i < MUNGE_OUT_COUNT; i++) child_vecs[i] = duckdb_struct_vector_get_child(output, (idx_t)i);

    for (idx_t row = 0; row < row_count; row++) {
        idx_t chrom_len = 0, a1_len = 0, a2_len = 0, id_len = 0, fasta_len = 0, tag_len = 4, mismatch_len = 12;
        const char *chrom;
        const char *a1;
        const char *a2;
        const char *id = NULL;
        const char *fasta_ref;
        const char *iffy_tag = "IFFY";
        const char *mismatch_tag = "REF_MISMATCH";
        char *fasta_path = NULL;
        char *err = NULL;
        char *norm_alt = NULL;
        char *norm_ref = NULL;
        char *ref_fetch = NULL;
        char *chrom_cstr = NULL;
        faidx_t *fai = NULL;
        double ns = NAN, ez = NAN, nc = NAN, es = NAN, se = NAN, lp = NAN, af = NAN, ac = NAN, ne = NAN;
        double si = NAN, i2 = NAN, cq = NAN;
        const char *ed = NULL;
        idx_t ed_len = 0;
        char *ed_flipped = NULL;
        bool swapped = false;
        int32_t pos;

        if (!row_is_valid(chrom_vec, row)) {
            set_munge_error(info, "bcftools_munge_row: chrom must be non-null");
            free(cached_fasta_path);
            return;
        }
        if (!row_is_valid(pos_vec, row)) {
            set_munge_error(info, "bcftools_munge_row: pos must be non-null");
            free(cached_fasta_path);
            return;
        }
        if (!row_is_valid(a1_vec, row) || !row_is_valid(a2_vec, row)) {
            set_munge_error(info, "bcftools_munge_row: a1 and a2 must be non-null");
            free(cached_fasta_path);
            return;
        }

        chrom = get_string_at(chrom_vec, row, &chrom_len);
        a1 = get_string_at(a1_vec, row, &a1_len);
        a2 = get_string_at(a2_vec, row, &a2_len);
        if (row_is_valid(fasta_ref_vec, row)) fasta_ref = get_string_at(fasta_ref_vec, row, &fasta_len);
        else { fasta_ref = NULL; fasta_len = 0; }
        if (row_is_valid(id_vec, row)) id = get_string_at(id_vec, row, &id_len);
        if (row_is_valid(iffy_tag_vec, row)) iffy_tag = get_string_at(iffy_tag_vec, row, &tag_len);
        if (row_is_valid(mismatch_tag_vec, row)) mismatch_tag = get_string_at(mismatch_tag_vec, row, &mismatch_len);

        if (!chrom || chrom_len == 0) {
            set_munge_error(info, "bcftools_munge_row: chrom must be non-empty");
            free(cached_fasta_path);
            return;
        }
        if (!a1 || a1_len == 0 || !a2 || a2_len == 0) {
            set_munge_error(info, "bcftools_munge_row: a1 and a2 must be non-empty");
            free(cached_fasta_path);
            return;
        }

        pos = (int32_t)get_int64_at(pos_vec, row);
        if (pos < 1) {
            set_munge_error(info, "bcftools_munge_row: pos must be >= 1");
            free(cached_fasta_path);
            return;
        }

        norm_alt = normalize_allele(a1, a1_len);
        norm_ref = normalize_allele(a2, a2_len);
        if (!norm_alt || !norm_ref) {
            free(norm_alt);
            free(norm_ref);
            set_munge_error(info, "bcftools_munge_row: failed to normalize alleles");
            free(cached_fasta_path);
            return;
        }

        if (row_is_valid(z_vec, row)) ez = get_double_at(z_vec, row);
        if (row_is_valid(n_cas_vec, row)) nc = get_double_at(n_cas_vec, row);
        if (row_is_valid(se_vec, row)) se = get_double_at(se_vec, row);
        if (row_is_valid(frq_vec, row)) af = get_double_at(frq_vec, row);
        if (row_is_valid(ac_vec, row)) ac = get_double_at(ac_vec, row);
        if (row_is_valid(n_vec, row)) ns = get_double_at(n_vec, row);
        if (row_is_valid(neff_vec, row)) ne = get_double_at(neff_vec, row);
        else if (row_is_valid(neffdiv2_vec, row)) ne = get_double_at(neffdiv2_vec, row) * 2.0;
        if (row_is_valid(beta_vec, row)) es = get_double_at(beta_vec, row);
        else if (row_is_valid(or_vec, row)) es = safe_log(get_double_at(or_vec, row));
        if (row_is_valid(lp_vec, row)) lp = get_double_at(lp_vec, row);
        else if (row_is_valid(p_vec, row)) lp = safe_minus_log10(get_double_at(p_vec, row));
        if (row_is_valid(info_vec, row)) si = get_double_at(info_vec, row);
        if (row_is_valid(het_i2_vec, row)) i2 = get_double_at(het_i2_vec, row);
        if (row_is_valid(het_lp_vec, row)) cq = get_double_at(het_lp_vec, row);
        else if (row_is_valid(het_p_vec, row)) cq = safe_minus_log10(get_double_at(het_p_vec, row));
        if (row_is_valid(dire_vec, row)) ed = get_string_at(dire_vec, row, &ed_len);
        if (row_is_valid(nc_override_vec, row)) nc = get_double_at(nc_override_vec, row);
        if (row_is_valid(ns_override_vec, row)) ns = get_double_at(ns_override_vec, row);
        else if (!isnan(nc) && row_is_valid(n_con_vec, row)) ns = nc + get_double_at(n_con_vec, row);
        if (row_is_valid(ne_override_vec, row)) ne = get_double_at(ne_override_vec, row);
        if (isnan(ac) && !isnan(af) && !isnan(ns)) ac = 2.0 * ns * af;

        /* M-M9: fai-only mode — when fasta_ref is NULL/empty, skip ref-matching
         * entirely (upstream munge.c:590 gates swap logic on ref_fname). Output
         * alleles as-is: A2→REF, A1→ALT, no swap, no FILTER tag. */
        if (!fasta_ref || fasta_len == 0) {
            duckdb_vector_assign_string_element_len(child_vecs[MUNGE_OUT_CHROM], row, chrom, chrom_len);
            ((int64_t *)duckdb_vector_get_data(child_vecs[MUNGE_OUT_POS]))[row] = pos;
            if (id && id_len > 0) duckdb_vector_assign_string_element_len(child_vecs[MUNGE_OUT_ID], row, id, id_len);
            else set_null_at(child_vecs[MUNGE_OUT_ID], row);
            duckdb_vector_assign_string_element(child_vecs[MUNGE_OUT_REF], row, norm_ref);
            duckdb_vector_assign_string_element(child_vecs[MUNGE_OUT_ALT], row, norm_alt);
            ((bool *)duckdb_vector_get_data(child_vecs[MUNGE_OUT_ALLELES_SWAPPED]))[row] = false;
            set_null_at(child_vecs[MUNGE_OUT_FILTER], row);
            assign_double_or_null(child_vecs[MUNGE_OUT_NS], row, ns);
            assign_double_or_null(child_vecs[MUNGE_OUT_EZ], row, ez);
            assign_double_or_null(child_vecs[MUNGE_OUT_NC], row, nc);
            assign_double_or_null(child_vecs[MUNGE_OUT_ES], row, es);
            assign_double_or_null(child_vecs[MUNGE_OUT_SE], row, se);
            assign_double_or_null(child_vecs[MUNGE_OUT_LP], row, lp);
            assign_double_or_null(child_vecs[MUNGE_OUT_AF], row, af);
            assign_double_or_null(child_vecs[MUNGE_OUT_AC], row, ac);
            assign_double_or_null(child_vecs[MUNGE_OUT_NE], row, ne);
            assign_double_or_null(child_vecs[MUNGE_OUT_SI], row, si);
            assign_double_or_null(child_vecs[MUNGE_OUT_I2], row, i2);
            assign_double_or_null(child_vecs[MUNGE_OUT_CQ], row, cq);
            if (ed && ed_len > 0) duckdb_vector_assign_string_element_len(child_vecs[MUNGE_OUT_ED], row, ed, ed_len);
            else set_null_at(child_vecs[MUNGE_OUT_ED], row);
            free(norm_alt);
            free(norm_ref);
            continue;
        }

        if (cached_fai && cached_fasta_path && cached_fasta_len == fasta_len &&
            memcmp(cached_fasta_path, fasta_ref, (size_t)fasta_len) == 0) {
            fai = cached_fai;
        } else {
            fasta_path = dup_span(fasta_ref, fasta_len);
            if (!fasta_path) {
                set_munge_error(info, "bcftools_munge_row: out of memory");
                free(norm_alt);
                free(norm_ref);
                free(cached_fasta_path);
                return;
            }
            fai = get_munge_fai(fasta_path, &err);
            if (fai) {
                free(cached_fasta_path);
                cached_fasta_path = fasta_path;
                cached_fasta_len = fasta_len;
                cached_fai = fai;
                fasta_path = NULL;
            }
            free(fasta_path);
        }
        if (!fai) {
            set_munge_error(info, err ? err : "bcftools_munge_row: failed to load FASTA");
            free(err);
            free(norm_alt);
            free(norm_ref);
            free(cached_fasta_path);
            return;
        }

        chrom_cstr = dup_span(chrom, chrom_len);
        if (!chrom_cstr) {
            set_munge_error(info, "bcftools_munge_row: out of memory");
            free(norm_alt);
            free(norm_ref);
            free(cached_fasta_path);
            return;
        }
        pthread_mutex_lock(&g_munge_fai_fetch_mutex);
        ref_fetch = fetch_sequence_flexible(
            fai,
            chrom_cstr,
            pos,
            pos + (hts_pos_t)((strlen(norm_ref) > strlen(norm_alt) ? strlen(norm_ref) : strlen(norm_alt))) - 1
        );
        pthread_mutex_unlock(&g_munge_fai_fetch_mutex);
        free(chrom_cstr);
        if (!ref_fetch) {
            /* Upstream hard-errors here, but in SQL context we emit the row
               with a FILTER tag so the query doesn't abort on scaffolds or
               contigs missing from the FASTA. */
            duckdb_vector_assign_string_element_len(child_vecs[MUNGE_OUT_CHROM], row, chrom, chrom_len);
            ((int64_t *)duckdb_vector_get_data(child_vecs[MUNGE_OUT_POS]))[row] = pos;
            if (id && id_len > 0) duckdb_vector_assign_string_element_len(child_vecs[MUNGE_OUT_ID], row, id, id_len);
            else set_null_at(child_vecs[MUNGE_OUT_ID], row);
            duckdb_vector_assign_string_element(child_vecs[MUNGE_OUT_REF], row, norm_ref);
            duckdb_vector_assign_string_element(child_vecs[MUNGE_OUT_ALT], row, norm_alt);
            ((bool *)duckdb_vector_get_data(child_vecs[MUNGE_OUT_ALLELES_SWAPPED]))[row] = false;
            duckdb_vector_assign_string_element(child_vecs[MUNGE_OUT_FILTER], row, "MissingContig");
            assign_double_or_null(child_vecs[MUNGE_OUT_NS], row, ns);
            assign_double_or_null(child_vecs[MUNGE_OUT_EZ], row, ez);
            assign_double_or_null(child_vecs[MUNGE_OUT_NC], row, nc);
            assign_double_or_null(child_vecs[MUNGE_OUT_ES], row, es);
            assign_double_or_null(child_vecs[MUNGE_OUT_SE], row, se);
            assign_double_or_null(child_vecs[MUNGE_OUT_LP], row, lp);
            assign_double_or_null(child_vecs[MUNGE_OUT_AF], row, af);
            assign_double_or_null(child_vecs[MUNGE_OUT_AC], row, ac);
            assign_double_or_null(child_vecs[MUNGE_OUT_NE], row, ne);
            assign_double_or_null(child_vecs[MUNGE_OUT_SI], row, si);
            assign_double_or_null(child_vecs[MUNGE_OUT_I2], row, i2);
            assign_double_or_null(child_vecs[MUNGE_OUT_CQ], row, cq);
            if (ed && ed_len > 0) duckdb_vector_assign_string_element_len(child_vecs[MUNGE_OUT_ED], row, ed, ed_len);
            else set_null_at(child_vecs[MUNGE_OUT_ED], row);
            free(norm_alt);
            free(norm_ref);
            continue;
        }

        {
            int ref_match = strncasecmp(ref_fetch, norm_ref, strlen(norm_ref)) == 0;
            int alt_match = strncasecmp(ref_fetch, norm_alt, strlen(norm_alt)) == 0;
            if (!ref_match && alt_match) {
                char *tmp = norm_ref;
                norm_ref = norm_alt;
                norm_alt = tmp;
                swapped = true;
                set_null_at(child_vecs[MUNGE_OUT_FILTER], row);
                if (!isnan(ez)) ez = -ez;
                if (!isnan(es)) es = -es;
                if (!isnan(af)) af = 1.0 - af;
                if (!isnan(ac) && !isnan(ns)) ac = 2.0 * ns - ac;
                else if (!isnan(ac)) ac = NAN;
                /* upstream munge.c:617-624 — flip ED direction on swap */
                if (ed && ed_len > 0) {
                    ed_flipped = dup_span(ed, ed_len);
                    if (ed_flipped) {
                        for (idx_t k = 0; k < ed_len; k++) {
                            if (ed_flipped[k] == '+') ed_flipped[k] = '-';
                            else if (ed_flipped[k] == '-') ed_flipped[k] = '+';
                        }
                    }
                }
            } else if (ref_match && alt_match) {
                duckdb_vector_assign_string_element_len(child_vecs[MUNGE_OUT_FILTER], row, iffy_tag, tag_len);
            } else if (!ref_match && !alt_match) {
                duckdb_vector_assign_string_element_len(child_vecs[MUNGE_OUT_FILTER], row, mismatch_tag, mismatch_len);
            } else {
                set_null_at(child_vecs[MUNGE_OUT_FILTER], row);
            }
        }

        duckdb_vector_assign_string_element_len(child_vecs[MUNGE_OUT_CHROM], row, chrom, chrom_len);
        ((int64_t *)duckdb_vector_get_data(child_vecs[MUNGE_OUT_POS]))[row] = pos;
        if (id && id_len > 0) duckdb_vector_assign_string_element_len(child_vecs[MUNGE_OUT_ID], row, id, id_len);
        else set_null_at(child_vecs[MUNGE_OUT_ID], row);
        duckdb_vector_assign_string_element(child_vecs[MUNGE_OUT_REF], row, norm_ref);
        duckdb_vector_assign_string_element(child_vecs[MUNGE_OUT_ALT], row, norm_alt);
        ((bool *)duckdb_vector_get_data(child_vecs[MUNGE_OUT_ALLELES_SWAPPED]))[row] = swapped;
        assign_double_or_null(child_vecs[MUNGE_OUT_NS], row, ns);
        assign_double_or_null(child_vecs[MUNGE_OUT_EZ], row, ez);
        assign_double_or_null(child_vecs[MUNGE_OUT_NC], row, nc);
        assign_double_or_null(child_vecs[MUNGE_OUT_ES], row, es);
        assign_double_or_null(child_vecs[MUNGE_OUT_SE], row, se);
        assign_double_or_null(child_vecs[MUNGE_OUT_LP], row, lp);
        assign_double_or_null(child_vecs[MUNGE_OUT_AF], row, af);
        assign_double_or_null(child_vecs[MUNGE_OUT_AC], row, ac);
        assign_double_or_null(child_vecs[MUNGE_OUT_NE], row, ne);
        assign_double_or_null(child_vecs[MUNGE_OUT_SI], row, si);
        assign_double_or_null(child_vecs[MUNGE_OUT_I2], row, i2);
        assign_double_or_null(child_vecs[MUNGE_OUT_CQ], row, cq);
        if (ed_flipped) {
            duckdb_vector_assign_string_element_len(child_vecs[MUNGE_OUT_ED], row, ed_flipped, ed_len);
        } else if (ed && ed_len > 0) {
            duckdb_vector_assign_string_element_len(child_vecs[MUNGE_OUT_ED], row, ed, ed_len);
        } else {
            set_null_at(child_vecs[MUNGE_OUT_ED], row);
        }

        free(ed_flipped);
        free(norm_alt);
        free(norm_ref);
        free(ref_fetch);
    }

    free(cached_fasta_path);
}

void register_munge_functions(duckdb_connection connection) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type double_type = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type fields[MUNGE_OUT_COUNT];
    duckdb_logical_type struct_type;

    duckdb_scalar_function_set_name(fn, "bcftools_munge_row");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    for (int i = 0; i < 17; i++) duckdb_scalar_function_add_parameter(fn, double_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, double_type);
    duckdb_scalar_function_add_parameter(fn, double_type);
    duckdb_scalar_function_add_parameter(fn, double_type);

    fields[MUNGE_OUT_CHROM] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[MUNGE_OUT_POS] = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    fields[MUNGE_OUT_ID] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[MUNGE_OUT_REF] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[MUNGE_OUT_ALT] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[MUNGE_OUT_ALLELES_SWAPPED] = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    fields[MUNGE_OUT_FILTER] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    for (int i = MUNGE_OUT_NS; i <= MUNGE_OUT_CQ; i++) fields[i] = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);
    fields[MUNGE_OUT_ED] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);

    struct_type = duckdb_create_struct_type(fields, MUNGE_FIELD_NAMES, MUNGE_OUT_COUNT);
    duckdb_scalar_function_set_return_type(fn, struct_type);
    duckdb_scalar_function_set_special_handling(fn);
    duckdb_scalar_function_set_function(fn, bcftools_munge_row_scalar);
    duckdb_register_scalar_function(connection, fn);
    duckdb_destroy_scalar_function(&fn);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&double_type);
    duckdb_destroy_logical_type(&bool_type);
    for (int i = 0; i < MUNGE_OUT_COUNT; i++) duckdb_destroy_logical_type(&fields[i]);
    duckdb_destroy_logical_type(&struct_type);
}
