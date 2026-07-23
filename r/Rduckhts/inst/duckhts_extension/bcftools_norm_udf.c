#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <ctype.h>
#include <limits.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>

#include "hts_io_tuning.h"

#define NORM_ALN_WIN 100
#define NORM_CACHE_MAX_ENTRIES 8
#define NORM_REF_FETCH_QUANTUM 65536
#define NORM_REF_FETCH_PAD 4096

typedef struct norm_cache_entry {
    char *fasta_path;
    char *fai_path;
    char *gzi_path;
    faidx_t *fai;
    char *window_contig;
    hts_pos_t window_beg;
    hts_pos_t window_end;
    char *window_seq;
    struct norm_cache_entry *next;
} norm_cache_entry_t;

/*
 * faidx_t and the cached reference window own mutable htslib/file state.  Never
 * share a norm cache entry across DuckDB worker threads; issue #17/PR #18 fixed
 * the same class of FASTA-cache race in munge/liftover.  The cache below is
 * deliberately pthread-local, and reference fetches are serialized defensively.
 */
static pthread_key_t g_norm_cache_key;
static pthread_once_t g_norm_cache_key_once = PTHREAD_ONCE_INIT;
static pthread_mutex_t g_norm_fai_fetch_mutex = PTHREAD_MUTEX_INITIALIZER;

enum {
    NORM_OUT_POS = 0,
    NORM_OUT_END_POS,
    NORM_OUT_REF,
    NORM_OUT_ALT,
    NORM_OUT_NORMED,
    NORM_OUT_STATUS,
    NORM_OUT_COUNT
};

static const char *NORM_FIELD_NAMES[NORM_OUT_COUNT] = {
    "pos_normed",
    "end_pos_normed",
    "ref_normed",
    "alt_normed",
    "normed",
    "norm_status"
};

typedef struct {
    int applicable;
    int changed;
    const char *status;
    int64_t pos1;
    int64_t end_pos1;
    char *ref;
    char **alts;
    int n_alt;
} norm_result_t;

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

static char *dup_span(const char *s, size_t n) {
    char *out;
    out = (char *)malloc(n + 1);
    if (!out) return NULL;
    if (n) memcpy(out, s, n);
    out[n] = '\0';
    return out;
}

static char *dup_cstr(const char *s) {
    if (!s) return NULL;
    return dup_span(s, strlen(s));
}

static void free_string_array(char **items, int n_items) {
    if (!items) return;
    for (int i = 0; i < n_items; i++) free(items[i]);
    free(items);
}

static void free_norm_result(norm_result_t *result) {
    if (!result) return;
    free(result->ref);
    free_string_array(result->alts, result->n_alt);
    memset(result, 0, sizeof(*result));
}

static void destroy_norm_cache(void *ptr) {
    norm_cache_entry_t *entry = (norm_cache_entry_t *)ptr;
    while (entry) {
        norm_cache_entry_t *next = entry->next;
        free(entry->fasta_path);
        free(entry->fai_path);
        free(entry->gzi_path);
        if (entry->fai) fai_destroy(entry->fai);
        free(entry->window_contig);
        free(entry->window_seq);
        free(entry);
        entry = next;
    }
}

static void init_norm_cache_key(void) {
    pthread_key_create(&g_norm_cache_key, destroy_norm_cache);
}

static norm_cache_entry_t *get_norm_cache_head(void) {
    pthread_once(&g_norm_cache_key_once, init_norm_cache_key);
    return (norm_cache_entry_t *)pthread_getspecific(g_norm_cache_key);
}

static void set_norm_cache_head(norm_cache_entry_t *head) {
    pthread_once(&g_norm_cache_key_once, init_norm_cache_key);
    pthread_setspecific(g_norm_cache_key, head);
}

static void trim_norm_cache_head(norm_cache_entry_t *head) {
    size_t count = 0;
    norm_cache_entry_t *entry = head;
    norm_cache_entry_t *prev = NULL;
    while (entry) {
        count++;
        if (count > NORM_CACHE_MAX_ENTRIES) {
            if (prev) prev->next = NULL;
            destroy_norm_cache(entry);
            return;
        }
        prev = entry;
        entry = entry->next;
    }
}

static int same_nullable_string(const char *a, const char *b) {
    if (!a && !b) return 1;
    if (!a || !b) return 0;
    return strcmp(a, b) == 0;
}

static int same_nullable_span(const char *a, const char *b, size_t b_len) {
    if (!a && !b) return 1;
    if (!a || !b) return 0;
    return strlen(a) == b_len && memcmp(a, b, b_len) == 0;
}

static void seq_to_upper_ascii(char *seq) {
    if (!seq) return;
    while (*seq) {
        *seq = (char)toupper((unsigned char)*seq);
        seq++;
    }
}

static inline int replace_iupac_codes(char *seq, int nseq) {
    int i, n = 0;
    for (i = 0; i < nseq; i++) {
        char c = (char)toupper((unsigned char)seq[i]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') {
            seq[i] = 'N';
            n++;
        } else {
            seq[i] = c;
        }
    }
    return n;
}

static inline int has_non_acgtn(const char *seq, int nseq) {
    for (int i = 0; i < nseq; i++) {
        char c = (char)toupper((unsigned char)seq[i]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T' && c != 'N') return 1;
    }
    return 0;
}

static int is_breakend_alt(const char *alt) {
    return alt && (strchr(alt, '[') || strchr(alt, ']'));
}

static int is_spanning_deletion_alt(const char *alt) {
    return alt && alt[0] == '*' && alt[1] == '\0';
}

static int is_symbolic_alt(const char *alt) {
    size_t len;
    if (!alt) return 0;
    len = strlen(alt);
    return len >= 3 && alt[0] == '<' && alt[len - 1] == '>';
}

static int is_gvcf_symbolic_alt(const char *alt) {
    return alt && (strcmp(alt, "<NON_REF>") == 0 || strcmp(alt, "<*>") == 0);
}

static int is_gvcf_ignored_alt(const char *alt) {
    return is_gvcf_symbolic_alt(alt) || is_spanning_deletion_alt(alt);
}

static char *resolve_fasta_contig_name(faidx_t *fai, const char *chrom) {
    static const char *mt_aliases[] = {"MT", "chrM", "M", NULL};
    char with_chr[512];
    const char *aliases[8];
    int idx = 0;

    if (!fai || !chrom || !*chrom) return NULL;

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
        if (faidx_has_seq(fai, aliases[i]) > 0) return dup_cstr(aliases[i]);
    }
    return NULL;
}

static norm_cache_entry_t *get_norm_cache_entry(const char *fasta_path,
                                                const char *fai_path,
                                                const char *gzi_path,
                                                char **err) {
    norm_cache_entry_t *head = get_norm_cache_head();
    norm_cache_entry_t *entry = NULL;
    norm_cache_entry_t *prev = NULL;
    faidx_t *loaded = NULL;

    for (entry = head; entry; prev = entry, entry = entry->next) {
        if (same_nullable_string(entry->fasta_path, fasta_path)
            && same_nullable_string(entry->fai_path, fai_path)
            && same_nullable_string(entry->gzi_path, gzi_path)) {
            if (prev) {
                prev->next = entry->next;
                entry->next = head;
                head = entry;
                set_norm_cache_head(head);
            }
            return entry;
        }
    }

    loaded = fai_load3_format(fasta_path,
                              (fai_path && *fai_path) ? fai_path : NULL,
                              (gzi_path && *gzi_path) ? gzi_path : NULL,
                              FAI_CREATE,
                              FAI_FASTA);
    if (!loaded) {
        *err = dup_cstr("bcftools_norm_row: failed to load FASTA index");
        return NULL;
    }
    duckhts_apply_remote_faidx_tuning(loaded, fasta_path, DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION);

    entry = (norm_cache_entry_t *)calloc(1, sizeof(*entry));
    if (!entry) {
        fai_destroy(loaded);
        *err = dup_cstr("bcftools_norm_row: out of memory");
        return NULL;
    }
    entry->fasta_path = dup_cstr(fasta_path);
    entry->fai_path = dup_cstr(fai_path);
    entry->gzi_path = dup_cstr(gzi_path);
    if (!entry->fasta_path || ((fai_path && *fai_path) && !entry->fai_path)
        || ((gzi_path && *gzi_path) && !entry->gzi_path)) {
        fai_destroy(loaded);
        free(entry->fasta_path);
        free(entry->fai_path);
        free(entry->gzi_path);
        free(entry);
        *err = dup_cstr("bcftools_norm_row: out of memory");
        return NULL;
    }
    entry->fai = loaded;
    entry->next = head;
    trim_norm_cache_head(entry);
    set_norm_cache_head(entry);
    return entry;
}

static char *fetch_windowed_sequence(norm_cache_entry_t *entry,
                                     const char *chrom,
                                     hts_pos_t start0,
                                     hts_pos_t end0) {
    hts_pos_t seq_len;
    hts_pos_t fetch_beg;
    hts_pos_t fetch_end;
    hts_pos_t contig_len;
    char *resolved = NULL;
    char *fetched = NULL;
    char *out = NULL;

    if (!entry || !entry->fai || !chrom || start0 < 0 || end0 < start0) return NULL;

    if (entry->window_seq && entry->window_contig && strcmp(entry->window_contig, chrom) == 0
        && start0 >= entry->window_beg && end0 <= entry->window_end) {
        return dup_span(entry->window_seq + (size_t)(start0 - entry->window_beg), (size_t)(end0 - start0 + 1));
    }

    resolved = resolve_fasta_contig_name(entry->fai, chrom);
    if (!resolved) return NULL;

    if (entry->window_seq && entry->window_contig && strcmp(entry->window_contig, resolved) == 0
        && start0 >= entry->window_beg && end0 <= entry->window_end) {
        out = dup_span(entry->window_seq + (size_t)(start0 - entry->window_beg), (size_t)(end0 - start0 + 1));
        free(resolved);
        return out;
    }

    contig_len = faidx_seq_len64(entry->fai, resolved);
    if (contig_len <= 0) {
        free(resolved);
        return NULL;
    }
    if (start0 >= contig_len) {
        free(resolved);
        return NULL;
    }
    if (end0 >= contig_len) end0 = contig_len - 1;

    fetch_beg = start0 > NORM_REF_FETCH_PAD ? start0 - NORM_REF_FETCH_PAD : 0;
    fetch_end = fetch_beg + NORM_REF_FETCH_QUANTUM - 1;
    if (fetch_end < end0) fetch_end = end0;
    if (fetch_end >= contig_len) fetch_end = contig_len - 1;

    pthread_mutex_lock(&g_norm_fai_fetch_mutex);
    fetched = faidx_fetch_seq64(entry->fai, resolved, fetch_beg, fetch_end, &seq_len);
    pthread_mutex_unlock(&g_norm_fai_fetch_mutex);
    if (!fetched || seq_len != fetch_end - fetch_beg + 1) {
        free(fetched);
        free(resolved);
        return NULL;
    }
    replace_iupac_codes(fetched, (int)seq_len);

    free(entry->window_contig);
    free(entry->window_seq);
    entry->window_contig = resolved;
    entry->window_seq = fetched;
    entry->window_beg = fetch_beg;
    entry->window_end = fetch_end;

    return dup_span(entry->window_seq + (size_t)(start0 - entry->window_beg), (size_t)(end0 - start0 + 1));
}

static int split_alt_text(const char *alt_text, size_t alt_len, char ***out_alts, int *out_n_alt) {
    char **alts = NULL;
    int n_alt = 0;
    int cap = 0;
    size_t start = 0;

    if (!alt_text) {
        *out_alts = NULL;
        *out_n_alt = 0;
        return 1;
    }
    if (alt_len == 1 && alt_text[0] == '.') {
        *out_alts = NULL;
        *out_n_alt = 0;
        return 1;
    }

    while (start <= alt_len) {
        size_t end = start;
        char *alt;
        while (end < alt_len && alt_text[end] != ',') end++;
        alt = dup_span(alt_text + start, end - start);
        if (!alt) {
            free_string_array(alts, n_alt);
            return 0;
        }
        if (n_alt == cap) {
            int new_cap = cap ? cap * 2 : 4;
            char **tmp = (char **)realloc(alts, sizeof(char *) * (size_t)new_cap);
            if (!tmp) {
                free(alt);
                free_string_array(alts, n_alt);
                return 0;
            }
            alts = tmp;
            cap = new_cap;
        }
        alts[n_alt++] = alt;
        if (end == alt_len) break;
        start = end + 1;
    }

    *out_alts = alts;
    *out_n_alt = n_alt;
    return 1;
}

static int extract_alt_list(duckdb_vector alt_vec,
                            idx_t row,
                            char ***out_alts,
                            int *out_n_alt,
                            char **err) {
    duckdb_logical_type alt_type;
    duckdb_type type_id;

    *out_alts = NULL;
    *out_n_alt = 0;

    if (!row_is_valid(alt_vec, row)) return 1;

    alt_type = duckdb_vector_get_column_type(alt_vec);
    type_id = duckdb_get_type_id(alt_type);

    if (type_id == DUCKDB_TYPE_VARCHAR) {
        idx_t alt_len = 0;
        const char *alt_text = get_string_at(alt_vec, row, &alt_len);
        duckdb_destroy_logical_type(&alt_type);
        if (!split_alt_text(alt_text, (size_t)alt_len, out_alts, out_n_alt)) {
            *err = dup_cstr("bcftools_norm_row: out of memory");
            return 0;
        }
        return 1;
    }

    if (type_id == DUCKDB_TYPE_LIST) {
        duckdb_list_entry *entries = (duckdb_list_entry *)duckdb_vector_get_data(alt_vec);
        duckdb_list_entry entry = entries[row];
        duckdb_vector child = duckdb_list_vector_get_child(alt_vec);
        char **alts = NULL;
        if (entry.length > 0) {
            alts = (char **)calloc((size_t)entry.length, sizeof(char *));
            if (!alts) {
                duckdb_destroy_logical_type(&alt_type);
                *err = dup_cstr("bcftools_norm_row: out of memory");
                return 0;
            }
        }
        for (idx_t i = 0; i < entry.length; i++) {
            idx_t child_row = entry.offset + i;
            idx_t alt_len = 0;
            const char *alt_text;
            if (!row_is_valid(child, child_row)) {
                continue;
            }
            alt_text = get_string_at(child, child_row, &alt_len);
            alts[i] = dup_span(alt_text, (size_t)alt_len);
            if (!alts[i]) {
                duckdb_destroy_logical_type(&alt_type);
                free_string_array(alts, (int)entry.length);
                *err = dup_cstr("bcftools_norm_row: out of memory");
                return 0;
            }
        }
        *out_alts = alts;
        *out_n_alt = (int)entry.length;
        duckdb_destroy_logical_type(&alt_type);
        return 1;
    }

    duckdb_destroy_logical_type(&alt_type);
    *err = dup_cstr("bcftools_norm_row: ALT must be VARCHAR or VARCHAR[]");
    return 0;
}

static int assign_alt_list(duckdb_vector vec, idx_t row, char **alts, int n_alt) {
    duckdb_list_entry *list_data;
    duckdb_list_entry entry;
    duckdb_vector child;

    entry.offset = duckdb_list_vector_get_size(vec);
    entry.length = (idx_t)n_alt;
    if (duckdb_list_vector_reserve(vec, entry.offset + entry.length) != DuckDBSuccess) return 0;
    if (duckdb_list_vector_set_size(vec, entry.offset + entry.length) != DuckDBSuccess) return 0;
    child = duckdb_list_vector_get_child(vec);
    list_data = (duckdb_list_entry *)duckdb_vector_get_data(vec);
    list_data[row] = entry;
    for (int i = 0; i < n_alt; i++) {
        idx_t child_row = entry.offset + (idx_t)i;
        if (!alts || !alts[i]) {
            set_null_at(child, child_row);
            continue;
        }
        duckdb_vector_assign_string_element(child, child_row, alts[i]);
    }
    return 1;
}

static void alt_to_list_scalar(duckdb_function_info info, duckdb_data_chunk input, duckdb_vector output) {
    duckdb_vector alt_vec = duckdb_data_chunk_get_vector(input, 0);
    idx_t row_count = duckdb_data_chunk_get_size(input);

    duckdb_list_vector_set_size(output, 0);
    for (idx_t row = 0; row < row_count; row++) {
        char **alt_items = NULL;
        int n_alt_items = 0;
        char *err = NULL;
        if (!extract_alt_list(alt_vec, row, &alt_items, &n_alt_items, &err)) {
            duckdb_scalar_function_set_error(info, err ? err : "duckhts_alt_to_list: failed to parse ALT");
            free(err);
            return;
        }
        if (!row_is_valid(alt_vec, row)) {
            set_null_at(output, row);
        } else if (!assign_alt_list(output, row, alt_items, n_alt_items)) {
            free_string_array(alt_items, n_alt_items);
            duckdb_scalar_function_set_error(info, "duckhts_alt_to_list: failed to assign ALT output");
            return;
        }
        free_string_array(alt_items, n_alt_items);
    }
}

static int copy_alt_array(char **src, int n_src, char ***out, int *out_n) {
    char **dst = NULL;
    if (n_src > 0) {
        dst = (char **)calloc((size_t)n_src, sizeof(char *));
        if (!dst) return 0;
    }
    for (int i = 0; i < n_src; i++) {
        if (!src || !src[i]) continue;
        dst[i] = dup_cstr(src[i]);
        if (!dst[i]) {
            free_string_array(dst, n_src);
            return 0;
        }
    }
    *out = dst;
    *out_n = n_src;
    return 1;
}

static int compare_alt_arrays(char **a, int n_a, char **b, int n_b) {
    if (n_a != n_b) return 0;
    for (int i = 0; i < n_a; i++) {
        const char *sa = a ? a[i] : NULL;
        const char *sb = b ? b[i] : NULL;
        if (!same_nullable_string(sa, sb)) return 0;
    }
    return 1;
}

static int prepare_sequence_allele(const char *src, kstring_t *dst) {
    dst->l = 0;
    if (!src) return 1;
    if (kputs(src, dst) < 0) return 0;
    seq_to_upper_ascii(dst->s);
    return 1;
}

static void free_kstring_array(kstring_t *arr, int n) {
    if (!arr) return;
    for (int i = 0; i < n; i++) free(arr[i].s);
    free(arr);
}

static hts_pos_t realign_left_cached(norm_cache_entry_t *cache,
                                     const char *chrom,
                                     kstring_t *als,
                                     int n_allele,
                                     hts_pos_t pos0,
                                     int aln_win) {
    char *ref = NULL;
    hts_pos_t new_pos = pos0;

    while (1) {
        int min_len = (int)als[0].l;
        int i;
        for (i = 1; i < n_allele; i++) {
            if (toupper((unsigned char)als[0].s[als[0].l - 1]) != toupper((unsigned char)als[i].s[als[i].l - 1])) break;
            if ((int)als[i].l < min_len) min_len = (int)als[i].l;
        }
        if (i != n_allele) break;
        if (min_len <= 1 && new_pos == 0) break;

        int pad_from_left = 0;
        for (i = 0; i < n_allele; i++) {
            als[i].l--;
            if (!als[i].l) pad_from_left = 1;
        }
        if (pad_from_left) {
            int npad = new_pos >= aln_win ? aln_win : (int)new_pos;
            free(ref);
            ref = fetch_windowed_sequence(cache, chrom, new_pos - npad, new_pos - 1);
            if (!ref) break;
            for (i = 0; i < n_allele; i++) {
                if (ks_resize(&als[i], als[i].l + (size_t)npad) < 0) {
                    free(ref);
                    return new_pos;
                }
                if (als[i].l) memmove(als[i].s + npad, als[i].s, als[i].l);
                memcpy(als[i].s, ref, (size_t)npad);
                als[i].l += (size_t)npad;
            }
            new_pos -= npad;
        }
    }
    free(ref);

    {
        int ntrim_left = 0;
        while (1) {
            int min_len = (int)als[0].l - ntrim_left;
            int i;
            for (i = 1; i < n_allele; i++) {
                if (toupper((unsigned char)als[0].s[ntrim_left]) != toupper((unsigned char)als[i].s[ntrim_left])) break;
                if (min_len > (int)als[i].l - ntrim_left) min_len = (int)als[i].l - ntrim_left;
            }
            if (i != n_allele || min_len <= 1) break;
            ntrim_left++;
        }
        if (ntrim_left) {
            for (int i = 0; i < n_allele; i++) {
                memmove(als[i].s, als[i].s + ntrim_left, als[i].l - (size_t)ntrim_left);
                als[i].l -= (size_t)ntrim_left;
                als[i].s[als[i].l] = '\0';
            }
            new_pos += ntrim_left;
        }
    }
    return new_pos;
}

static void set_result_passthrough(norm_result_t *result,
                                   int64_t pos1,
                                   int64_t end_pos1,
                                   const char *ref,
                                   char **alts,
                                   int n_alt,
                                   const char *status) {
    result->applicable = 0;
    result->changed = 0;
    result->status = status;
    result->pos1 = pos1;
    result->end_pos1 = end_pos1;
    result->ref = dup_cstr(ref);
    copy_alt_array(alts, n_alt, &result->alts, &result->n_alt);
}

static void set_result_applicable(norm_result_t *result,
                                  int changed,
                                  const char *status,
                                  int64_t pos1,
                                  int64_t end_pos1,
                                  char *ref,
                                  char **alts,
                                  int n_alt) {
    result->applicable = 1;
    result->changed = changed;
    result->status = status;
    result->pos1 = pos1;
    result->end_pos1 = end_pos1;
    result->ref = ref;
    result->alts = alts;
    result->n_alt = n_alt;
}

static int normalize_variant_row(norm_cache_entry_t *cache,
                                 const char *chrom,
                                 int64_t pos1,
                                 const char *ref_in,
                                 char **alt_in,
                                 int n_alt_in,
                                 int has_end_pos,
                                 int64_t end_pos_in,
                                 int has_svlen,
                                 int64_t svlen_in,
                                 norm_result_t *result,
                                 char **err) {
    char *orig_ref = NULL;
    char **orig_alt = NULL;
    int n_orig_alt = 0;
    char *ref_fetch = NULL;
    kstring_t *als = NULL;
    kstring_t *sym = NULL;
    int n_allele;
    int symbolic_alts = 1;
    int any_breakend = 0;
    int any_spanning = 0;
    int64_t original_end;
    hts_pos_t new_pos0;
    int full_ref_len;
    char *out_ref = NULL;
    char **out_alts = NULL;
    int changed;

    memset(result, 0, sizeof(*result));

    if (!chrom || !*chrom || !ref_in || !*ref_in || pos1 < 1) {
        original_end = (pos1 >= 1 && ref_in) ? pos1 + (int64_t)strlen(ref_in) - 1 : 0;
        set_result_passthrough(result, pos1, original_end, ref_in, alt_in, n_alt_in, "NullInput");
        return 1;
    }

    orig_ref = dup_cstr(ref_in);
    if (!orig_ref) goto oom;
    seq_to_upper_ascii(orig_ref);
    if (!copy_alt_array(alt_in, n_alt_in, &orig_alt, &n_orig_alt)) goto oom;
    for (int i = 0; i < n_orig_alt; i++) {
        if (orig_alt[i]) seq_to_upper_ascii(orig_alt[i]);
    }

    original_end = has_end_pos ? end_pos_in : pos1 + (int64_t)strlen(orig_ref) - 1;

    if (n_orig_alt == 0) {
        set_result_applicable(result, 0, "RefOnly", pos1, original_end, orig_ref, orig_alt, n_orig_alt);
        return 1;
    }

    for (int i = 0; i < n_orig_alt; i++) {
        if (!orig_alt[i] || !*orig_alt[i]) {
            set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt, "NullInput");
            free(orig_ref);
            free_string_array(orig_alt, n_orig_alt);
            return 1;
        }
        any_breakend |= is_breakend_alt(orig_alt[i]);
        any_spanning |= is_spanning_deletion_alt(orig_alt[i]);
    }

    int n_gvcf_ignored = 0;
    int has_gvcf_symbolic = 0;
    for (int i = 0; i < n_orig_alt; i++) {
        if (!is_gvcf_ignored_alt(orig_alt[i])) continue;
        n_gvcf_ignored++;
        if (is_gvcf_symbolic_alt(orig_alt[i])) has_gvcf_symbolic = 1;
    }
    if (n_gvcf_ignored > 0) {
        if (n_gvcf_ignored == n_orig_alt) {
            if (has_gvcf_symbolic) {
                set_result_applicable(result, 0, "GVCFReferenceBlock", pos1, original_end, orig_ref, orig_alt, n_orig_alt);
            } else {
                set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt, "SpanningDeletion");
                free(orig_ref);
                free_string_array(orig_alt, n_orig_alt);
            }
            return 1;
        }

        int n_real_alt = n_orig_alt - n_gvcf_ignored;
        char **real_alt = (char **)calloc((size_t)n_real_alt, sizeof(char *));
        int *real_to_orig = (int *)calloc((size_t)n_real_alt, sizeof(int));
        norm_result_t real_result;
        memset(&real_result, 0, sizeof(real_result));
        if (!real_alt || !real_to_orig) {
            free(real_alt);
            free(real_to_orig);
            goto oom;
        }
        int real_i = 0;
        for (int i = 0; i < n_orig_alt; i++) {
            if (is_gvcf_ignored_alt(orig_alt[i])) continue;
            real_alt[real_i] = dup_cstr(orig_alt[i]);
            if (!real_alt[real_i]) {
                free_string_array(real_alt, n_real_alt);
                free(real_to_orig);
                goto oom;
            }
            real_to_orig[real_i] = i;
            real_i++;
        }

        if (!normalize_variant_row(cache, chrom, pos1, orig_ref, real_alt, n_real_alt,
                                   has_end_pos, end_pos_in, has_svlen, svlen_in,
                                   &real_result, err)) {
            free_string_array(real_alt, n_real_alt);
            free(real_to_orig);
            free(orig_ref);
            free_string_array(orig_alt, n_orig_alt);
            return 0;
        }
        free_string_array(real_alt, n_real_alt);
        free(real_to_orig);

        if (!real_result.applicable) {
            set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt,
                                   real_result.status ? real_result.status : "GVCFMixedPassthrough");
            free_norm_result(&real_result);
            free(orig_ref);
            free_string_array(orig_alt, n_orig_alt);
            return 1;
        }

        char **merged_alt = (char **)calloc((size_t)n_orig_alt, sizeof(char *));
        char *merged_ref = dup_cstr(real_result.ref);
        if (!merged_alt || !merged_ref) {
            free(merged_ref);
            free_string_array(merged_alt, n_orig_alt);
            free_norm_result(&real_result);
            goto oom;
        }
        real_i = 0;
        for (int i = 0; i < n_orig_alt; i++) {
            const char *src_alt = NULL;
            if (is_gvcf_ignored_alt(orig_alt[i])) {
                src_alt = orig_alt[i];
            } else {
                src_alt = real_i < real_result.n_alt ? real_result.alts[real_i] : NULL;
                real_i++;
            }
            merged_alt[i] = dup_cstr(src_alt);
            if (!merged_alt[i]) {
                free(merged_ref);
                free_string_array(merged_alt, n_orig_alt);
                free_norm_result(&real_result);
                goto oom;
            }
        }
        int64_t merged_end = (has_end_pos && has_gvcf_symbolic) ? original_end : real_result.end_pos1;
        changed = real_result.pos1 != pos1 || merged_end != original_end ||
                  strcmp(merged_ref, orig_ref) != 0 ||
                  !compare_alt_arrays(merged_alt, n_orig_alt, orig_alt, n_orig_alt);
        set_result_applicable(result, changed,
                              changed ? "Normalized" : "Unchanged",
                              real_result.pos1, merged_end,
                              merged_ref, merged_alt, n_orig_alt);
        free_norm_result(&real_result);
        free(orig_ref);
        free_string_array(orig_alt, n_orig_alt);
        return 1;
    }

    if (any_breakend) {
        set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt, "Breakend");
        free(orig_ref);
        free_string_array(orig_alt, n_orig_alt);
        return 1;
    }
    if (any_spanning) {
        set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt, "SpanningDeletion");
        free(orig_ref);
        free_string_array(orig_alt, n_orig_alt);
        return 1;
    }
    if (has_non_acgtn(orig_ref, (int)strlen(orig_ref))) {
        set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt, "InvalidRefAllele");
        free(orig_ref);
        free_string_array(orig_alt, n_orig_alt);
        return 1;
    }

    ref_fetch = fetch_windowed_sequence(cache, chrom, pos1 - 1, pos1 + (int64_t)strlen(orig_ref) - 2);
    if (!ref_fetch) {
        set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt, "MissingContig");
        free(orig_ref);
        free_string_array(orig_alt, n_orig_alt);
        return 1;
    }
    if (strcasecmp(ref_fetch, orig_ref) != 0) {
        free(ref_fetch);
        set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt, "RefMismatch");
        free(orig_ref);
        free_string_array(orig_alt, n_orig_alt);
        return 1;
    }
    free(ref_fetch);
    ref_fetch = NULL;

    {
        size_t ref_len = strlen(orig_ref);
        size_t min_len = ref_len;
        int tail = toupper((unsigned char)orig_ref[ref_len - 1]);
        int head = toupper((unsigned char)orig_ref[0]);
        int common_tail = 1;
        int common_head = 1;
        int can_skip_realign = 1;

        for (int i = 0; i < n_orig_alt; i++) {
            char *alt = orig_alt[i];
            size_t alt_len;

            if (is_symbolic_alt(alt)) {
                can_skip_realign = 0;
                break;
            }
            alt_len = strlen(alt);
            if (has_non_acgtn(alt, (int)alt_len)) {
                set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt, "InvalidAltAllele");
                free(orig_ref);
                free_string_array(orig_alt, n_orig_alt);
                return 1;
            }
            if (alt_len == ref_len && strcmp(alt, orig_ref) == 0) {
                set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt, "DuplicateAllele");
                free(orig_ref);
                free_string_array(orig_alt, n_orig_alt);
                return 1;
            }
            if (alt_len < min_len) min_len = alt_len;
            if (toupper((unsigned char)alt[alt_len - 1]) != tail) common_tail = 0;
            if (toupper((unsigned char)alt[0]) != head) common_head = 0;
        }

        if (can_skip_realign && !common_tail && (!common_head || min_len <= 1)) {
            set_result_applicable(result, 0, "Unchanged", pos1, pos1 + (int64_t)ref_len - 1,
                                  orig_ref, orig_alt, n_orig_alt);
            return 1;
        }
    }

    n_allele = n_orig_alt + 1;
    als = (kstring_t *)calloc((size_t)n_allele, sizeof(kstring_t));
    sym = (kstring_t *)calloc((size_t)n_allele, sizeof(kstring_t));
    if (!als || !sym) goto oom;
    if (!prepare_sequence_allele(orig_ref, &als[0])) goto oom;

    for (int i = 1; i < n_allele; i++) {
        char *alt = orig_alt[i - 1];
        if (is_symbolic_alt(alt)) {
            if (!strncmp(alt, "<DEL", 4)) {
                int64_t rlen = has_end_pos && end_pos_in >= pos1 ? end_pos_in - pos1 + 1 : (int64_t)strlen(orig_ref);
                char *expanded_ref;
                if (rlen <= 0) {
                    set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt,
                                           "SymbolicInsufficientMetadata");
                    free_kstring_array(als, n_allele);
                    free_kstring_array(sym, n_allele);
                    free(orig_ref);
                    free_string_array(orig_alt, n_orig_alt);
                    return 1;
                }
                expanded_ref = fetch_windowed_sequence(cache, chrom, pos1 - 1, pos1 + rlen - 2);
                if (!expanded_ref) {
                    set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt, "MissingContig");
                    free_kstring_array(als, n_allele);
                    free_kstring_array(sym, n_allele);
                    free(orig_ref);
                    free_string_array(orig_alt, n_orig_alt);
                    return 1;
                }
                free(als[0].s);
                als[0].s = expanded_ref;
                als[0].l = strlen(expanded_ref);
                als[0].m = als[0].l + 1;
                if (kputsn(als[0].s, 1, &als[i]) < 0) goto oom;
            } else if (!strncmp(alt, "<DUP", 4)) {
                char *dup_ref;
                if (!has_svlen || svlen_in <= 0) {
                    set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt,
                                           "SymbolicInsufficientMetadata");
                    free_kstring_array(als, n_allele);
                    free_kstring_array(sym, n_allele);
                    free(orig_ref);
                    free_string_array(orig_alt, n_orig_alt);
                    return 1;
                }
                dup_ref = fetch_windowed_sequence(cache, chrom, pos1 - 1, pos1 + svlen_in - 1);
                if (!dup_ref) {
                    set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt, "MissingContig");
                    free_kstring_array(als, n_allele);
                    free_kstring_array(sym, n_allele);
                    free(orig_ref);
                    free_string_array(orig_alt, n_orig_alt);
                    return 1;
                }
                if (kputs(dup_ref, &als[i]) < 0) {
                    free(dup_ref);
                    goto oom;
                }
                free(dup_ref);
            } else {
                set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt,
                                       "UnsupportedSymbolicAllele");
                free_kstring_array(als, n_allele);
                free_kstring_array(sym, n_allele);
                free(orig_ref);
                free_string_array(orig_alt, n_orig_alt);
                return 1;
            }
            if (kputs(alt, &sym[i]) < 0) goto oom;
            continue;
        }
        symbolic_alts = 0;
        if (has_non_acgtn(alt, (int)strlen(alt))) {
            set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt, "InvalidAltAllele");
            free_kstring_array(als, n_allele);
            free_kstring_array(sym, n_allele);
            free(orig_ref);
            free_string_array(orig_alt, n_orig_alt);
            return 1;
        }
        if (!prepare_sequence_allele(alt, &als[i])) goto oom;
        if (als[i].l == als[0].l && strcasecmp(als[i].s, als[0].s) == 0) {
            set_result_passthrough(result, pos1, original_end, orig_ref, orig_alt, n_orig_alt, "DuplicateAllele");
            free_kstring_array(als, n_allele);
            free_kstring_array(sym, n_allele);
            free(orig_ref);
            free_string_array(orig_alt, n_orig_alt);
            return 1;
        }
    }

    new_pos0 = realign_left_cached(cache, chrom, als, n_allele, pos1 - 1, NORM_ALN_WIN);
    full_ref_len = (int)als[0].l;
    if (symbolic_alts && als[0].l > 1) {
        als[0].l = 1;
        als[0].s[1] = '\0';
    }

    out_ref = dup_span(als[0].s, als[0].l);
    if (!out_ref) goto oom;
    if (n_orig_alt > 0) {
        out_alts = (char **)calloc((size_t)n_orig_alt, sizeof(char *));
        if (!out_alts) goto oom;
    }
    for (int i = 0; i < n_orig_alt; i++) {
        if (sym[i + 1].l) {
            out_alts[i] = dup_span(sym[i + 1].s, sym[i + 1].l);
        } else {
            out_alts[i] = dup_span(als[i + 1].s, als[i + 1].l);
        }
        if (!out_alts[i]) goto oom;
    }

    changed = new_pos0 + 1 != pos1 || strcmp(out_ref, orig_ref) != 0 || !compare_alt_arrays(out_alts, n_orig_alt, orig_alt, n_orig_alt);
    set_result_applicable(result,
                          changed,
                          changed ? "Normalized" : "Unchanged",
                          new_pos0 + 1,
                          new_pos0 + full_ref_len,
                          out_ref,
                          out_alts,
                          n_orig_alt);
    out_ref = NULL;
    out_alts = NULL;
    free_kstring_array(als, n_allele);
    free_kstring_array(sym, n_allele);
    free(orig_ref);
    free_string_array(orig_alt, n_orig_alt);
    return 1;

oom:
    free(out_ref);
    free_string_array(out_alts, n_orig_alt);
    free_kstring_array(als, n_allele);
    free_kstring_array(sym, n_allele);
    free(ref_fetch);
    free(orig_ref);
    free_string_array(orig_alt, n_orig_alt);
    *err = dup_cstr("bcftools_norm_row: out of memory");
    return 0;
}

static void bcftools_norm_row_scalar(duckdb_function_info info,
                                     duckdb_data_chunk input,
                                     duckdb_vector output) {
    duckdb_vector chrom_vec = duckdb_data_chunk_get_vector(input, 0);
    duckdb_vector pos_vec = duckdb_data_chunk_get_vector(input, 1);
    duckdb_vector ref_vec = duckdb_data_chunk_get_vector(input, 2);
    duckdb_vector alt_vec = duckdb_data_chunk_get_vector(input, 3);
    duckdb_vector fasta_vec = duckdb_data_chunk_get_vector(input, 4);
    duckdb_vector end_pos_vec = duckdb_data_chunk_get_vector(input, 5);
    duckdb_vector svlen_vec = duckdb_data_chunk_get_vector(input, 6);
    duckdb_vector fai_vec = duckdb_data_chunk_get_vector(input, 7);
    duckdb_vector gzi_vec = duckdb_data_chunk_get_vector(input, 8);
    duckdb_vector child_vecs[NORM_OUT_COUNT];
    idx_t row_count = duckdb_data_chunk_get_size(input);
    char *cached_fasta_path = NULL;
    char *cached_fai_path = NULL;
    char *cached_gzi_path = NULL;
    norm_cache_entry_t *cached_entry = NULL;

    for (int i = 0; i < NORM_OUT_COUNT; i++) child_vecs[i] = duckdb_struct_vector_get_child(output, (idx_t)i);
    duckdb_list_vector_set_size(child_vecs[NORM_OUT_ALT], 0);

    for (idx_t row = 0; row < row_count; row++) {
        idx_t chrom_len = 0, ref_len = 0, fasta_len = 0, fai_len = 0, gzi_len = 0;
        const char *chrom = NULL;
        const char *ref = NULL;
        const char *fasta_ref = NULL;
        const char *fai_path = NULL;
        const char *gzi_path = NULL;
        char *chrom_key = NULL;
        char *ref_key = NULL;
        char *fasta_key = NULL;
        char *fai_key = NULL;
        char *gzi_key = NULL;
        char **alt_items = NULL;
        int n_alt_items = 0;
        char *err = NULL;
        norm_result_t result;
        int64_t pos1 = row_is_valid(pos_vec, row) ? get_int64_at(pos_vec, row) : 0;
        int64_t end_pos = row_is_valid(end_pos_vec, row) ? get_int64_at(end_pos_vec, row) : 0;
        int64_t svlen = row_is_valid(svlen_vec, row) ? get_int64_at(svlen_vec, row) : 0;
        int has_end_pos = row_is_valid(end_pos_vec, row);
        int has_svlen = row_is_valid(svlen_vec, row);

        memset(&result, 0, sizeof(result));

        if (row_is_valid(chrom_vec, row)) chrom = get_string_at(chrom_vec, row, &chrom_len);
        if (row_is_valid(ref_vec, row)) ref = get_string_at(ref_vec, row, &ref_len);
        if (row_is_valid(fasta_vec, row)) fasta_ref = get_string_at(fasta_vec, row, &fasta_len);
        if (row_is_valid(fai_vec, row)) fai_path = get_string_at(fai_vec, row, &fai_len);
        if (row_is_valid(gzi_vec, row)) gzi_path = get_string_at(gzi_vec, row, &gzi_len);

        if (!fasta_ref || fasta_len == 0) {
            duckdb_scalar_function_set_error(info, "bcftools_norm_row: fasta_ref must be non-empty");
            free(cached_fasta_path);
            free(cached_fai_path);
            free(cached_gzi_path);
            return;
        }

        if (chrom && chrom_len > 0) {
            chrom_key = dup_span(chrom, chrom_len);
            if (!chrom_key) {
                duckdb_scalar_function_set_error(info, "bcftools_norm_row: out of memory");
                free(cached_fasta_path);
                free(cached_fai_path);
                free(cached_gzi_path);
                return;
            }
        }
        if (ref && ref_len > 0) {
            ref_key = dup_span(ref, ref_len);
            if (!ref_key) {
                free(chrom_key);
                duckdb_scalar_function_set_error(info, "bcftools_norm_row: out of memory");
                free(cached_fasta_path);
                free(cached_fai_path);
                free(cached_gzi_path);
                return;
            }
        }

        if (!cached_entry || !same_nullable_span(cached_fasta_path, fasta_ref, (size_t)fasta_len)
            || !same_nullable_span(cached_fai_path, (fai_path && fai_len > 0) ? fai_path : NULL, (size_t)fai_len)
            || !same_nullable_span(cached_gzi_path, (gzi_path && gzi_len > 0) ? gzi_path : NULL, (size_t)gzi_len)) {
            fasta_key = dup_span(fasta_ref, fasta_len);
            if (!fasta_key) {
                free(chrom_key);
                free(ref_key);
                duckdb_scalar_function_set_error(info, "bcftools_norm_row: out of memory");
                free(cached_fasta_path);
                free(cached_fai_path);
                free(cached_gzi_path);
                return;
            }
            if (fai_path && fai_len > 0) {
                fai_key = dup_span(fai_path, fai_len);
                if (!fai_key) {
                    free(chrom_key);
                    free(ref_key);
                    free(fasta_key);
                    duckdb_scalar_function_set_error(info, "bcftools_norm_row: out of memory");
                    free(cached_fasta_path);
                    free(cached_fai_path);
                    free(cached_gzi_path);
                    return;
                }
            }
            if (gzi_path && gzi_len > 0) {
                gzi_key = dup_span(gzi_path, gzi_len);
                if (!gzi_key) {
                    free(chrom_key);
                    free(ref_key);
                    free(fasta_key);
                    free(fai_key);
                    duckdb_scalar_function_set_error(info, "bcftools_norm_row: out of memory");
                    free(cached_fasta_path);
                    free(cached_fai_path);
                    free(cached_gzi_path);
                    return;
                }
            }

            cached_entry = get_norm_cache_entry(fasta_key, fai_key, gzi_key, &err);
            if (!cached_entry) {
                duckdb_scalar_function_set_error(info, err ? err : "bcftools_norm_row: failed to load FASTA index");
                free(err);
                free(chrom_key);
                free(ref_key);
                free(fasta_key);
                free(fai_key);
                free(gzi_key);
                free(cached_fasta_path);
                free(cached_fai_path);
                free(cached_gzi_path);
                return;
            }
            free(cached_fasta_path);
            free(cached_fai_path);
            free(cached_gzi_path);
            cached_fasta_path = fasta_key;
            cached_fai_path = fai_key;
            cached_gzi_path = gzi_key;
            fasta_key = NULL;
            fai_key = NULL;
            gzi_key = NULL;
        }
        free(fasta_key);
        free(fai_key);
        free(gzi_key);

        if (!extract_alt_list(alt_vec, row, &alt_items, &n_alt_items, &err)) {
            duckdb_scalar_function_set_error(info, err ? err : "bcftools_norm_row: failed to read ALT");
            free(err);
            free(chrom_key);
            free(ref_key);
            free(cached_fasta_path);
            free(cached_fai_path);
            free(cached_gzi_path);
            return;
        }
        if (!normalize_variant_row(cached_entry,
                                   chrom_key,
                                   pos1,
                                   ref_key,
                                   alt_items,
                                   n_alt_items,
                                   has_end_pos,
                                   end_pos,
                                   has_svlen,
                                   svlen,
                                   &result,
                                   &err)) {
            duckdb_scalar_function_set_error(info, err ? err : "bcftools_norm_row: normalization failed");
            free(err);
            free(chrom_key);
            free(ref_key);
            free_string_array(alt_items, n_alt_items);
            free(cached_fasta_path);
            free(cached_fai_path);
            free(cached_gzi_path);
            return;
        }

        if (result.ref) duckdb_vector_assign_string_element(child_vecs[NORM_OUT_REF], row, result.ref);
        else set_null_at(child_vecs[NORM_OUT_REF], row);
        if (!assign_alt_list(child_vecs[NORM_OUT_ALT], row, result.alts, result.n_alt)) {
            duckdb_scalar_function_set_error(info, "bcftools_norm_row: failed to assign ALT output");
            free_norm_result(&result);
            free(chrom_key);
            free(ref_key);
            free_string_array(alt_items, n_alt_items);
            free(cached_fasta_path);
            free(cached_fai_path);
            free(cached_gzi_path);
            return;
        }
        if (result.pos1 > 0) ((int64_t *)duckdb_vector_get_data(child_vecs[NORM_OUT_POS]))[row] = result.pos1;
        else set_null_at(child_vecs[NORM_OUT_POS], row);
        if (result.end_pos1 > 0) ((int64_t *)duckdb_vector_get_data(child_vecs[NORM_OUT_END_POS]))[row] = result.end_pos1;
        else set_null_at(child_vecs[NORM_OUT_END_POS], row);
        if (result.applicable) ((bool *)duckdb_vector_get_data(child_vecs[NORM_OUT_NORMED]))[row] = result.changed ? true : false;
        else set_null_at(child_vecs[NORM_OUT_NORMED], row);
        if (result.status) duckdb_vector_assign_string_element(child_vecs[NORM_OUT_STATUS], row, result.status);
        else set_null_at(child_vecs[NORM_OUT_STATUS], row);

        free_norm_result(&result);
        free(chrom_key);
        free(ref_key);
        free_string_array(alt_items, n_alt_items);
    }

    free(cached_fasta_path);
    free(cached_fai_path);
    free(cached_gzi_path);
}

static void add_bcftools_norm_scalar(
    duckdb_scalar_function_set set,
    duckdb_logical_type        alt_type) {
    duckdb_scalar_function fn = duckdb_create_scalar_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type fields[NORM_OUT_COUNT];
    duckdb_logical_type struct_type;

    duckdb_scalar_function_set_name(fn, "bcftools_norm_row");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, alt_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, bigint_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_add_parameter(fn, varchar_type);

    fields[NORM_OUT_POS] = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    fields[NORM_OUT_END_POS] = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    fields[NORM_OUT_REF] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    fields[NORM_OUT_ALT] = duckdb_create_list_type(varchar_type);
    fields[NORM_OUT_NORMED] = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    fields[NORM_OUT_STATUS] = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    struct_type = duckdb_create_struct_type(fields, NORM_FIELD_NAMES, NORM_OUT_COUNT);

    duckdb_scalar_function_set_return_type(fn, struct_type);
    duckdb_scalar_function_set_special_handling(fn);
    duckdb_scalar_function_set_function(fn, bcftools_norm_row_scalar);
    duckdb_add_scalar_function_to_set(set, fn);
    duckdb_destroy_scalar_function(&fn);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    for (int i = 0; i < NORM_OUT_COUNT; i++) duckdb_destroy_logical_type(&fields[i]);
    duckdb_destroy_logical_type(&struct_type);
}

void register_bcftools_norm_functions(duckdb_connection connection) {
    duckdb_scalar_function fn;
    duckdb_scalar_function_set norm_set;
    duckdb_scalar_function_set alt_set;
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type varchar_list_type = duckdb_create_list_type(varchar_type);

    norm_set = duckdb_create_scalar_function_set("bcftools_norm_row");
    add_bcftools_norm_scalar(norm_set, varchar_type);
    add_bcftools_norm_scalar(norm_set, varchar_list_type);
    duckdb_register_scalar_function_set(connection, norm_set);
    duckdb_destroy_scalar_function_set(&norm_set);

    alt_set = duckdb_create_scalar_function_set("duckhts_alt_to_list");
    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "duckhts_alt_to_list");
    duckdb_scalar_function_add_parameter(fn, varchar_type);
    duckdb_scalar_function_set_return_type(fn, varchar_list_type);
    duckdb_scalar_function_set_special_handling(fn);
    duckdb_scalar_function_set_function(fn, alt_to_list_scalar);
    duckdb_add_scalar_function_to_set(alt_set, fn);
    duckdb_destroy_scalar_function(&fn);

    fn = duckdb_create_scalar_function();
    duckdb_scalar_function_set_name(fn, "duckhts_alt_to_list");
    duckdb_scalar_function_add_parameter(fn, varchar_list_type);
    duckdb_scalar_function_set_return_type(fn, varchar_list_type);
    duckdb_scalar_function_set_special_handling(fn);
    duckdb_scalar_function_set_function(fn, alt_to_list_scalar);
    duckdb_add_scalar_function_to_set(alt_set, fn);
    duckdb_destroy_scalar_function(&fn);

    duckdb_register_scalar_function_set(connection, alt_set);
    duckdb_destroy_scalar_function_set(&alt_set);

    duckdb_destroy_logical_type(&varchar_list_type);
    duckdb_destroy_logical_type(&varchar_type);
}
