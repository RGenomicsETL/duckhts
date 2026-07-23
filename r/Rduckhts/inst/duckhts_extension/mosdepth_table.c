/**
 * DuckHTS native mosdepth-compatible rewrite subset.
 *
 * Current scope:
 *   - indexed BAM/CRAM input
 *   - fast mode, fragment mode, and default mode
 *   - summary/global distribution/per-base output
 *   - optional by := <window> or by := <bed>
 *   - optional quantize and thresholds outputs
 *   - optional read-group and fragment-length filters
 *
 * Deferred options still error loudly instead of silently diverging from
 * mosdepth.
 */

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <pthread.h>

#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/khash.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>

#include "include/hts_io_tuning.h"

#define MOSDEPTH_DEFAULT_PRECISION 2
#define MOSDEPTH_MAX_COVERAGE 400000
#define MOSDEPTH_CSI_MIN_SHIFT 14
/* Upstream mosdepth uses 512 buckets for region/global distribution arrays. */
#define MOSDEPTH_REGION_DIST_BUCKETS 512
/* Upstream CountStat for --use-median is initialized with 65536 coverage bins. */
#define MOSDEPTH_MEDIAN_COUNTS 65536

KHASH_MAP_INIT_STR(mdmate, bam1_t *)
KHASH_SET_INIT_STR(mdstr)

typedef struct {
    int64_t start;
    int64_t stop;
    char *name;
} mosdepth_region_t;

typedef struct {
    mosdepth_region_t *data;
    size_t count;
    size_t cap;
} mosdepth_region_list_t;

typedef struct {
    int64_t *data;
    size_t len;
} mosdepth_dist_t;

typedef struct {
    int64_t *data;
    size_t count;
} mosdepth_threshold_list_t;

typedef struct {
    int64_t *cuts;
    char **labels;
    size_t count;
} mosdepth_quantize_t;

typedef struct {
    hts_pos_t start;
    hts_pos_t stop;
} mosdepth_span_t;

typedef struct {
    mosdepth_span_t *data;
    size_t count;
    size_t cap;
} mosdepth_span_list_t;

typedef struct {
    uint64_t cum_depth;
    int64_t cum_length;
    uint32_t min_depth;
    uint32_t max_depth;
} mosdepth_stat_t;

typedef struct {
    uint64_t *counts;
    size_t len;
    uint64_t n;
} mosdepth_count_stat_t;

typedef struct {
    char *prefix;
    char *path;
    char *chrom;
    char *by;
    char *fasta;
    char *read_groups;
    char *index_path;
    char *summary_path;
    char *global_dist_path;
    char *per_base_path;
    char *regions_path;
    char *region_dist_path;
    char *quantized_path;
    char *thresholds_path;
    int64_t threads;
    int64_t processing_threads;
    int64_t flag;
    int64_t include_flag;
    int64_t mapq;
    int64_t min_frag_len;
    int64_t max_frag_len;
    int64_t precision;
    mosdepth_threshold_list_t thresholds;
    mosdepth_quantize_t quantize;
    int no_per_base;
    int fast_mode;
    int fragment_mode;
    int use_median;
    int overwrite;
    int emitted;
} mosdepth_bind_t;

static char *append_suffix(const char *prefix, const char *suffix) {
    size_t a = strlen(prefix);
    size_t b = strlen(suffix);
    char *out = (char *)duckdb_malloc(a + b + 1);
    if (!out) return NULL;
    memcpy(out, prefix, a);
    memcpy(out + a, suffix, b + 1);
    return out;
}

static char *append_csi_suffix(const char *path) {
    return append_suffix(path, ".csi");
}

static int file_exists(const char *path) {
    struct stat st;
    return path && stat(path, &st) == 0;
}

static void set_null(duckdb_vector vec, idx_t row) {
    duckdb_vector_ensure_validity_writable(vec);
    uint64_t *validity = duckdb_vector_get_validity(vec);
    validity[row / 64] &= ~((uint64_t)1 << (row % 64));
}

static int is_digits_only(const char *s) {
    if (!s || !*s) return 0;
    for (const unsigned char *p = (const unsigned char *)s; *p; p++) {
        if (!isdigit(*p)) return 0;
    }
    return 1;
}

static int parse_int64_span_local(const char *s, int len, int64_t *out) {
    if (!s || len <= 0 || !out) return 0;
    char *tmp = (char *)malloc((size_t)len + 1);
    if (!tmp) return 0;
    memcpy(tmp, s, (size_t)len);
    tmp[len] = '\0';
    char *end = NULL;
    long long v = strtoll(tmp, &end, 10);
    int ok = (end && *end == '\0');
    if (ok) *out = (int64_t)v;
    free(tmp);
    return ok;
}

static int is_meta_bed_line(const char *s) {
    if (!s || !*s) return 1;
    if (s[0] == '#') return 1;
    if (strncmp(s, "track", 5) == 0) return 1;
    if (strncmp(s, "browser", 7) == 0) return 1;
    return 0;
}

static const char *get_field_span(const char *s, int idx, int *len) {
    int cur = 0;
    const char *start = s;
    const char *p = s;
    while (*p) {
        if (*p == '\t') {
            if (cur == idx) {
                *len = (int)(p - start);
                return start;
            }
            cur++;
            start = p + 1;
        }
        p++;
    }
    if (cur == idx) {
        *len = (int)(p - start);
        return start;
    }
    *len = 0;
    return NULL;
}

static void stat_clear(mosdepth_stat_t *stat) {
    stat->cum_depth = 0;
    stat->cum_length = 0;
    stat->min_depth = UINT32_MAX;
    stat->max_depth = 0;
}

static void stat_add_depth(mosdepth_stat_t *stat, uint32_t depth) {
    stat->cum_depth += (uint64_t)depth;
    stat->cum_length++;
    if (depth < stat->min_depth) stat->min_depth = depth;
    if (depth > stat->max_depth) stat->max_depth = depth;
}

static void stat_merge(mosdepth_stat_t *dst, const mosdepth_stat_t *src) {
    if (!dst || !src) return;
    if (src->cum_length == 0) return;
    dst->cum_depth += src->cum_depth;
    dst->cum_length += src->cum_length;
    if (src->min_depth < dst->min_depth) dst->min_depth = src->min_depth;
    if (src->max_depth > dst->max_depth) dst->max_depth = src->max_depth;
}

static int depth_bucket_value(uint32_t depth) {
    if ((int)depth > MOSDEPTH_MAX_COVERAGE) return MOSDEPTH_MAX_COVERAGE - 10;
    return (int)depth;
}

/*
 * Upstream mosdepth buckets window-region distributions from the mean depth
 * using the integer bucket that matches its current implementation behavior.
 * Keep this helper trivial and inline so the compiler can fold it directly
 * into the hot region loop.
 */
static inline int mosdepth_window_region_bucket(double region_value) {
    int bucket = (int)(region_value + 0.5);
    int max_bucket = MOSDEPTH_REGION_DIST_BUCKETS - 1;
    return bucket > max_bucket ? max_bucket : bucket;
}

static int count_stat_init(mosdepth_count_stat_t *stat, size_t len) {
    if (!stat) return -1;
    memset(stat, 0, sizeof(*stat));
    if (len == 0) return 0;
    stat->counts = (uint64_t *)calloc(len, sizeof(uint64_t));
    if (!stat->counts) return -1;
    stat->len = len;
    return 0;
}

static void count_stat_destroy(mosdepth_count_stat_t *stat) {
    if (!stat) return;
    free(stat->counts);
    stat->counts = NULL;
    stat->len = 0;
    stat->n = 0;
}

static void count_stat_clear(mosdepth_count_stat_t *stat) {
    if (!stat || !stat->counts || stat->n == 0) {
        if (stat) stat->n = 0;
        return;
    }
    memset(stat->counts, 0, stat->len * sizeof(uint64_t));
    stat->n = 0;
}

static inline void count_stat_add(mosdepth_count_stat_t *stat, uint32_t depth) {
    size_t idx;
    if (!stat || !stat->counts || stat->len == 0) return;
    idx = (size_t)depth;
    if (idx >= stat->len) idx = stat->len - 1;
    stat->counts[idx] += 1;
    stat->n += 1;
}

static uint32_t count_stat_median(const mosdepth_count_stat_t *stat) {
    uint64_t stop_n;
    uint64_t cum = 0;

    if (!stat || !stat->counts || stat->len == 0 || stat->n == 0) return 0;
    stop_n = (uint64_t)(0.5 + (double)stat->n * 0.5);
    for (size_t i = 0; i < stat->len; i++) {
        cum += stat->counts[i];
        if (cum >= stop_n) return (uint32_t)i;
    }
    return 0;
}

static int dist_ensure(mosdepth_dist_t *dist, size_t idx) {
    if (!dist) return -1;
    if (idx < dist->len) return 0;
    size_t new_len = idx + 10;
    int64_t *new_data = (int64_t *)realloc(dist->data, new_len * sizeof(int64_t));
    if (!new_data) return -1;
    for (size_t i = dist->len; i < new_len; i++) new_data[i] = 0;
    dist->data = new_data;
    dist->len = new_len;
    return 0;
}

static void dist_destroy(mosdepth_dist_t *dist) {
    if (!dist) return;
    free(dist->data);
    dist->data = NULL;
    dist->len = 0;
}

static int threshold_cmp(const void *a, const void *b) {
    const int64_t aa = *(const int64_t *)a;
    const int64_t bb = *(const int64_t *)b;
    if (aa < bb) return -1;
    if (aa > bb) return 1;
    return 0;
}

static void threshold_list_destroy(mosdepth_threshold_list_t *thresholds) {
    if (!thresholds) return;
    free(thresholds->data);
    thresholds->data = NULL;
    thresholds->count = 0;
}

static void quantize_destroy(mosdepth_quantize_t *quantize) {
    if (!quantize) return;
    free(quantize->cuts);
    if (quantize->labels) {
        for (size_t i = 0; i + 1 < quantize->count; i++) free(quantize->labels[i]);
    }
    free(quantize->labels);
    quantize->cuts = NULL;
    quantize->labels = NULL;
    quantize->count = 0;
}

static int span_list_reserve(mosdepth_span_list_t *list, size_t need) {
    if (!list) return -1;
    if (need <= list->cap) return 0;
    size_t new_cap = list->cap ? list->cap : 8;
    while (new_cap < need) new_cap *= 2;
    mosdepth_span_t *new_data =
        (mosdepth_span_t *)realloc(list->data, new_cap * sizeof(mosdepth_span_t));
    if (!new_data) return -1;
    list->data = new_data;
    list->cap = new_cap;
    return 0;
}

static void span_list_destroy(mosdepth_span_list_t *list) {
    if (!list) return;
    free(list->data);
    list->data = NULL;
    list->count = 0;
    list->cap = 0;
}

static int span_list_append(mosdepth_span_list_t *list, hts_pos_t start, hts_pos_t stop) {
    if (!list || stop <= start) return 0;
    if (span_list_reserve(list, list->count + 1) != 0) return -1;
    list->data[list->count].start = start;
    list->data[list->count].stop = stop;
    list->count++;
    return 0;
}

static void read_group_set_destroy(khash_t(mdstr) *set) {
    if (!set) return;
    for (khiter_t it = kh_begin(set); it != kh_end(set); ++it) {
        if (kh_exist(set, it)) free((char *)kh_key(set, it));
    }
    kh_destroy(mdstr, set);
}

static int parse_read_group_set(const char *spec, khash_t(mdstr) **out,
                                char *err, size_t errlen) {
    khash_t(mdstr) *set = NULL;
    char *copy = NULL;
    char *saveptr = NULL;
    char *token = NULL;

    if (!spec || !*spec) return 0;

    set = kh_init(mdstr);
    copy = strdup(spec);
    if (!set || !copy) {
        snprintf(err, errlen, "duckhts_mosdepth: out of memory");
        read_group_set_destroy(set);
        free(copy);
        return -1;
    }

    token = strtok_r(copy, ",", &saveptr);
    while (token) {
        int ret = 0;
        char *key = NULL;
        if (*token == '\0') {
            snprintf(err, errlen, "duckhts_mosdepth: invalid read_groups string");
            read_group_set_destroy(set);
            free(copy);
            return -1;
        }
        key = strdup(token);
        if (!key) {
            snprintf(err, errlen, "duckhts_mosdepth: out of memory");
            read_group_set_destroy(set);
            free(copy);
            return -1;
        }
        (void)kh_put(mdstr, set, key, &ret);
        if (ret < 0) {
            free(key);
            snprintf(err, errlen, "duckhts_mosdepth: out of memory");
            read_group_set_destroy(set);
            free(copy);
            return -1;
        }
        if (ret == 0) free(key);
        token = strtok_r(NULL, ",", &saveptr);
    }

    free(copy);
    *out = set;
    return 0;
}

static int record_in_read_groups(const bam1_t *rec, khash_t(mdstr) *read_groups) {
    uint8_t *aux;
    const char *rg;

    if (!read_groups) return 1;
    aux = bam_aux_get(rec, "RG");
    if (!aux) return 0;
    rg = bam_aux2Z(aux);
    if (!rg) return 0;
    return kh_get(mdstr, read_groups, rg) != kh_end(read_groups);
}

static inline void add_coverage_delta(int32_t *coverage, int64_t chrom_len,
                                      hts_pos_t start, hts_pos_t stop, int32_t delta) {
    if (!coverage || delta == 0 || stop <= start) return;
    if (start < 0) start = 0;
    if (stop > chrom_len) stop = chrom_len;
    if (start >= chrom_len || stop <= start) return;
    coverage[start] += delta;
    coverage[stop] -= delta;
}

static int collect_query_ref_spans(const bam1_t *rec, int64_t chrom_len,
                                   mosdepth_span_list_t *spans) {
    hts_pos_t pos;
    hts_pos_t span_start = -1;
    hts_pos_t last_stop = -1;
    uint32_t *cigar;

    if (!rec || !spans) return -1;
    spans->count = 0;
    if (rec->core.n_cigar == 0) return 0;
    if (span_list_reserve(spans, rec->core.n_cigar) != 0) return -1;

    pos = rec->core.pos;
    cigar = bam_get_cigar(rec);
    for (uint32_t i = 0; i < rec->core.n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        hts_pos_t oplen = (hts_pos_t)bam_cigar_oplen(cigar[i]);
        int type = bam_cigar_type(op);

        if (!(type & 2)) continue;
        if (type & 1) {
            if (span_start < 0 || pos != last_stop) {
                if (span_start >= 0 &&
                    span_list_append(spans, span_start, last_stop) != 0) {
                    return -1;
                }
                span_start = pos;
            }
            last_stop = pos + oplen;
        }
        pos += oplen;
    }

    if (span_start >= 0 &&
        span_list_append(spans, span_start, last_stop) != 0) {
        return -1;
    }

    for (size_t i = 0; i < spans->count; i++) {
        if (spans->data[i].start < 0) spans->data[i].start = 0;
        if (spans->data[i].stop > chrom_len) spans->data[i].stop = chrom_len;
    }
    return 0;
}

static void add_record_cigar_coverage(const bam1_t *rec, int32_t *coverage, int64_t chrom_len) {
    hts_pos_t pos;
    hts_pos_t span_start = -1;
    hts_pos_t last_stop = -1;
    uint32_t *cigar;

    if (!rec || !coverage || rec->core.n_cigar == 0) return;

    cigar = bam_get_cigar(rec);
    if (rec->core.n_cigar == 1) {
        int op = bam_cigar_op(cigar[0]);
        if ((bam_cigar_type(op) & 3) == 3) {
            add_coverage_delta(coverage, chrom_len, rec->core.pos,
                               rec->core.pos + (hts_pos_t)bam_cigar_oplen(cigar[0]), 1);
            return;
        }
    }

    pos = rec->core.pos;
    for (uint32_t i = 0; i < rec->core.n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        hts_pos_t oplen = (hts_pos_t)bam_cigar_oplen(cigar[i]);
        int type = bam_cigar_type(op);

        if (!(type & 2)) continue;
        if (type & 1) {
            if (span_start < 0 || pos != last_stop) {
                if (span_start >= 0) add_coverage_delta(coverage, chrom_len, span_start, last_stop, 1);
                span_start = pos;
            }
            last_stop = pos + oplen;
        }
        pos += oplen;
    }
    if (span_start >= 0) add_coverage_delta(coverage, chrom_len, span_start, last_stop, 1);
}

static int subtract_pair_overlap(const bam1_t *rec, const bam1_t *mate,
                                 int32_t *coverage, int64_t chrom_len) {
    mosdepth_span_list_t rec_spans = {0};
    mosdepth_span_list_t mate_spans = {0};
    int rc = -1;
    size_t i = 0, j = 0;

    if (!rec || !mate || !coverage) return 0;
    if (rec->core.n_cigar == 1 && mate->core.n_cigar == 1) {
        add_coverage_delta(coverage, chrom_len, rec->core.pos, bam_endpos(mate), -1);
        return 0;
    }

    if (collect_query_ref_spans(rec, chrom_len, &rec_spans) != 0) goto cleanup;
    if (collect_query_ref_spans(mate, chrom_len, &mate_spans) != 0) goto cleanup;

    while (i < rec_spans.count && j < mate_spans.count) {
        hts_pos_t rec_stop = rec_spans.data[i].stop;
        hts_pos_t mate_stop = mate_spans.data[j].stop;
        hts_pos_t start = rec_spans.data[i].start > mate_spans.data[j].start
                              ? rec_spans.data[i].start
                              : mate_spans.data[j].start;
        hts_pos_t stop = rec_stop < mate_stop ? rec_stop : mate_stop;
        if (stop > start) add_coverage_delta(coverage, chrom_len, start, stop, -1);
        if (rec_stop <= mate_stop) i++;
        if (mate_stop <= rec_stop) j++;
    }

    rc = 0;

cleanup:
    span_list_destroy(&rec_spans);
    span_list_destroy(&mate_spans);
    return rc;
}

static void seen_mates_destroy(khash_t(mdmate) *seen) {
    if (!seen) return;
    for (khiter_t it = kh_begin(seen); it != kh_end(seen); ++it) {
        if (kh_exist(seen, it)) {
            free((char *)kh_key(seen, it));
            bam_destroy1(kh_val(seen, it));
        }
    }
    kh_destroy(mdmate, seen);
}

static int handle_default_overlap(const bam1_t *rec, khash_t(mdmate) *seen,
                                  int32_t *coverage, int64_t chrom_len) {
    const char *qname;
    hts_pos_t start;
    hts_pos_t stop;
    hts_pos_t matepos;

    if (!rec || !seen) return 0;
    if (!(rec->core.flag & BAM_FPROPER_PAIR) || (rec->core.flag & BAM_FSUPPLEMENTARY)) return 0;

    qname = bam_get_qname(rec);
    start = rec->core.pos;
    stop = bam_endpos(rec);
    matepos = rec->core.mpos;

    if (rec->core.tid == rec->core.mtid && stop > matepos &&
        (start < matepos ||
         (start == matepos && kh_get(mdmate, seen, qname) == kh_end(seen)))) {
        bam1_t *copy = bam_dup1(rec);
        char *key = strdup(qname);
        khiter_t it;
        int ret = 0;
        if (!copy || !key) {
            free(key);
            bam_destroy1(copy);
            return -1;
        }
        it = kh_put(mdmate, seen, key, &ret);
        if (ret < 0) {
            free(key);
            bam_destroy1(copy);
            return -1;
        }
        if (ret == 0) {
            free(key);
            bam_destroy1(kh_val(seen, it));
        }
        kh_val(seen, it) = copy;
        return 0;
    }

    {
        khiter_t it = kh_get(mdmate, seen, qname);
        if (it != kh_end(seen)) {
            bam1_t *mate = kh_val(seen, it);
            int rc = subtract_pair_overlap(rec, mate, coverage, chrom_len);
            free((char *)kh_key(seen, it));
            bam_destroy1(mate);
            kh_del(mdmate, seen, it);
            return rc;
        }
    }

    return 0;
}

static int parse_thresholds_string(const char *spec, mosdepth_threshold_list_t *out,
                                   char *err, size_t errlen) {
    char *copy = NULL;
    char *saveptr = NULL;
    char *token = NULL;
    size_t count = 0;
    size_t cap = 1;
    int64_t *values = NULL;

    if (!spec || !*spec) return 0;

    for (const char *p = spec; *p; p++) {
        if (*p == ',') cap++;
    }

    copy = strdup(spec);
    values = (int64_t *)malloc(cap * sizeof(int64_t));
    if (!copy || !values) {
        snprintf(err, errlen, "duckhts_mosdepth: out of memory");
        free(copy);
        free(values);
        return -1;
    }

    token = strtok_r(copy, ",", &saveptr);
    while (token) {
        char *endptr = NULL;
        long long value;
        if (*token == '\0') {
            snprintf(err, errlen, "duckhts_mosdepth: invalid thresholds string");
            free(copy);
            free(values);
            return -1;
        }
        errno = 0;
        value = strtoll(token, &endptr, 10);
        if (errno != 0 || !endptr || *endptr != '\0') {
            snprintf(err, errlen, "duckhts_mosdepth: invalid thresholds string");
            free(copy);
            free(values);
            return -1;
        }
        values[count++] = (int64_t)value;
        token = strtok_r(NULL, ",", &saveptr);
    }

    qsort(values, count, sizeof(int64_t), threshold_cmp);
    out->data = values;
    out->count = count;
    free(copy);
    return 0;
}

static int int64_cmp(const void *a, const void *b) {
    const int64_t aa = *(const int64_t *)a;
    const int64_t bb = *(const int64_t *)b;
    if (aa < bb) return -1;
    if (aa > bb) return 1;
    return 0;
}

static int append_literal(char **buf, size_t *len, size_t *cap, const char *src, size_t src_len) {
    if (!buf || !len || !cap || (!src && src_len > 0)) return -1;
    if (*len + src_len + 1 > *cap) {
        size_t new_cap = *cap ? *cap : 32;
        while (*len + src_len + 1 > new_cap) new_cap *= 2;
        char *new_buf = (char *)realloc(*buf, new_cap);
        if (!new_buf) return -1;
        *buf = new_buf;
        *cap = new_cap;
    }
    if (src_len > 0) memcpy(*buf + *len, src, src_len);
    *len += src_len;
    (*buf)[*len] = '\0';
    return 0;
}

static int append_char(char **buf, size_t *len, size_t *cap, char ch) {
    return append_literal(buf, len, cap, &ch, 1);
}

static int append_int64(char **buf, size_t *len, size_t *cap, int64_t value) {
    char tmp[64];
    int written = snprintf(tmp, sizeof(tmp), "%" PRId64, value);
    if (written < 0) return -1;
    return append_literal(buf, len, cap, tmp, (size_t)written);
}

static int parse_quantize_string(const char *spec, mosdepth_quantize_t *out,
                                 char *err, size_t errlen) {
    char *normalized = NULL;
    char *copy = NULL;
    char *saveptr = NULL;
    char *token = NULL;
    size_t normalized_len = 0;
    size_t normalized_cap = 0;
    size_t cut_count = 0;
    size_t label_count = 0;
    int64_t *cuts = NULL;
    char **labels = NULL;

    if (!spec || !*spec) return 0;

    if (strchr(spec, ':') == NULL) {
        if (append_char(&normalized, &normalized_len, &normalized_cap, ':') != 0 ||
            append_literal(&normalized, &normalized_len, &normalized_cap, spec, strlen(spec)) != 0 ||
            append_char(&normalized, &normalized_len, &normalized_cap, ':') != 0) {
            snprintf(err, errlen, "duckhts_mosdepth: out of memory");
            goto fail;
        }
    } else if (append_literal(&normalized, &normalized_len, &normalized_cap, spec, strlen(spec)) != 0) {
        snprintf(err, errlen, "duckhts_mosdepth: out of memory");
        goto fail;
    }

    if (normalized[0] == ':') {
        char *tmp = NULL;
        size_t tmp_len = 0;
        size_t tmp_cap = 0;
        if (append_char(&tmp, &tmp_len, &tmp_cap, '0') != 0 ||
            append_literal(&tmp, &tmp_len, &tmp_cap, normalized, normalized_len) != 0) {
            free(tmp);
            snprintf(err, errlen, "duckhts_mosdepth: out of memory");
            goto fail;
        }
        free(normalized);
        normalized = tmp;
        normalized_len = tmp_len;
        normalized_cap = tmp_cap;
    }
    if (normalized[normalized_len - 1] == ':') {
        if (append_int64(&normalized, &normalized_len, &normalized_cap, INT64_MAX) != 0) {
            snprintf(err, errlen, "duckhts_mosdepth: out of memory");
            goto fail;
        }
    }
    if (strstr(normalized, "::") != NULL) {
        snprintf(err, errlen, "duckhts_mosdepth: invalid quantize string");
        goto fail;
    }

    for (size_t i = 0; i < normalized_len; i++) {
        if (normalized[i] == ':') cut_count++;
    }
    cut_count += 1;
    cuts = (int64_t *)malloc(cut_count * sizeof(int64_t));
    if (!cuts) {
        snprintf(err, errlen, "duckhts_mosdepth: out of memory");
        goto fail;
    }

    copy = strdup(normalized);
    if (!copy) {
        snprintf(err, errlen, "duckhts_mosdepth: out of memory");
        goto fail;
    }
    token = strtok_r(copy, ":", &saveptr);
    while (token) {
        char *endptr = NULL;
        long long value;
        errno = 0;
        value = strtoll(token, &endptr, 10);
        if (errno != 0 || !endptr || *endptr != '\0') {
            snprintf(err, errlen, "duckhts_mosdepth: invalid quantize string");
            goto fail;
        }
        cuts[label_count++] = (int64_t)value;
        token = strtok_r(NULL, ":", &saveptr);
    }
    if (label_count < 2) {
        snprintf(err, errlen, "duckhts_mosdepth: invalid quantize string");
        goto fail;
    }
    qsort(cuts, label_count, sizeof(int64_t), int64_cmp);

    labels = (char **)calloc(label_count - 1, sizeof(char *));
    if (!labels) {
        snprintf(err, errlen, "duckhts_mosdepth: out of memory");
        goto fail;
    }
    for (size_t i = 0; i + 1 < label_count; i++) {
        char label[128];
        int written;
        if (cuts[i + 1] == INT64_MAX) {
            written = snprintf(label, sizeof(label), "%" PRId64 ":inf", cuts[i]);
        } else {
            written = snprintf(label, sizeof(label), "%" PRId64 ":%" PRId64, cuts[i], cuts[i + 1]);
        }
        if (written < 0) {
            snprintf(err, errlen, "duckhts_mosdepth: invalid quantize string");
            goto fail;
        }
        labels[i] = strdup(label);
        if (!labels[i]) {
            snprintf(err, errlen, "duckhts_mosdepth: out of memory");
            goto fail;
        }
    }

    out->cuts = cuts;
    out->labels = labels;
    out->count = label_count;
    free(copy);
    free(normalized);
    return 0;

fail:
    free(copy);
    free(normalized);
    free(cuts);
    if (labels) {
        for (size_t i = 0; i + 1 < label_count; i++) free(labels[i]);
    }
    free(labels);
    return -1;
}

static inline int quantize_bucket(const mosdepth_quantize_t *quantize, int32_t depth) {
    if (!quantize || !quantize->cuts || quantize->count < 2) return -1;
    if ((int64_t)depth < quantize->cuts[0] || (int64_t)depth > quantize->cuts[quantize->count - 1]) {
        return -1;
    }
    for (size_t i = 0; i < quantize->count; i++) {
        if (quantize->cuts[i] > (int64_t)depth) return (int)i - 1;
        if (quantize->cuts[i] == (int64_t)depth) return (int)i;
    }
    return (int)quantize->count - 1;
}

static int dist_add_count(mosdepth_dist_t *dist, int depth, int64_t count) {
    if (!dist || depth < 0 || count <= 0) return 0;
    if (depth > MOSDEPTH_MAX_COVERAGE) depth = MOSDEPTH_MAX_COVERAGE - 10;
    if (dist_ensure(dist, (size_t)depth) != 0) return -1;
    dist->data[depth] += count;
    return 0;
}

static int dist_sum_into(const mosdepth_dist_t *from, mosdepth_dist_t *to) {
    if (!from || !to || !from->data || from->len == 0) return 0;
    if (dist_ensure(to, from->len - 1) != 0) return -1;
    for (size_t i = 0; i < from->len; i++) {
        to->data[i] += from->data[i];
    }
    return 0;
}

static int write_distribution(FILE *fh, const char *chrom, const mosdepth_dist_t *dist,
                              int precision) {
    if (!fh || !chrom || !dist || !dist->data || dist->len == 0) return 0;
    int64_t sum = 0;
    for (size_t i = 0; i < dist->len; i++) sum += dist->data[i];
    if (sum < 1) return 0;

    double cum = 0.0;
    for (size_t i = dist->len; i > 0; i--) {
        size_t depth = i - 1;
        int64_t count = dist->data[depth];
        if (depth > 300 && count == 0) continue;
        cum += (double)count / (double)sum;
        if (cum < 8e-5) continue;
        if (fprintf(fh, "%s\t%zu\t%.*f\n", chrom, depth, precision, cum) < 0) {
            return -1;
        }
    }
    return 0;
}

static int write_summary(FILE *fh, int *header_written, const char *label,
                         const mosdepth_stat_t *stat, int precision) {
    if (!fh || !header_written || !label || !stat) return -1;
    if (!*header_written) {
        if (fprintf(fh, "chrom\tlength\tbases\tmean\tmin\tmax\n") < 0) return -1;
        *header_written = 1;
    }
    double mean = stat->cum_length > 0 ? (double)stat->cum_depth / (double)stat->cum_length : 0.0;
    uint32_t min_depth = stat->cum_length > 0 ? stat->min_depth : 0;
    if (fprintf(fh, "%s\t%" PRId64 "\t%" PRIu64 "\t%.*f\t%u\t%u\n",
                label, stat->cum_length, stat->cum_depth, precision, mean,
                min_depth, stat->max_depth) < 0) {
        return -1;
    }
    return 0;
}

static int bgzf_write_line(BGZF *fp, kstring_t *line) {
    if (!fp || !line) return -1;
    if (line->l > 0 && bgzf_write(fp, line->s, (size_t)line->l) < 0) return -1;
    if (bgzf_write(fp, "\n", 1) < 0) return -1;
    return 0;
}

static void region_list_destroy(mosdepth_region_list_t *list) {
    if (!list) return;
    for (size_t i = 0; i < list->count; i++) {
        free(list->data[i].name);
    }
    free(list->data);
    list->data = NULL;
    list->count = 0;
    list->cap = 0;
}

static int region_cmp(const void *a, const void *b) {
    const mosdepth_region_t *ra = (const mosdepth_region_t *)a;
    const mosdepth_region_t *rb = (const mosdepth_region_t *)b;
    if (ra->start < rb->start) return -1;
    if (ra->start > rb->start) return 1;
    if (ra->stop < rb->stop) return -1;
    if (ra->stop > rb->stop) return 1;
    return 0;
}

static int region_list_append(mosdepth_region_list_t *list,
                              int64_t start, int64_t stop, const char *name) {
    if (!list) return -1;
    if (list->count == list->cap) {
        size_t new_cap = list->cap ? list->cap * 2 : 16;
        mosdepth_region_t *new_data =
            (mosdepth_region_t *)realloc(list->data, new_cap * sizeof(mosdepth_region_t));
        if (!new_data) return -1;
        list->data = new_data;
        list->cap = new_cap;
    }
    list->data[list->count].start = start;
    list->data[list->count].stop = stop;
    list->data[list->count].name = name ? strdup(name) : NULL;
    if (name && !list->data[list->count].name) return -1;
    list->count++;
    return 0;
}

static int parse_bed_regions(const char *path, sam_hdr_t *hdr,
                             mosdepth_region_list_t *per_tid,
                             char *err, size_t errlen) {
    htsFile *fp = NULL;
    kstring_t line = {0, 0, NULL};
    int rc = -1;

    fp = hts_open(path, "r");
    if (!fp) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to open BED file '%s'", path);
        goto cleanup;
    }

    while (hts_getline(fp, '\n', &line) >= 0) {
        if (line.l == 0 || is_meta_bed_line(line.s)) continue;

        int len0 = 0, len1 = 0, len2 = 0, len3 = 0;
        const char *f0 = get_field_span(line.s, 0, &len0);
        const char *f1 = get_field_span(line.s, 1, &len1);
        const char *f2 = get_field_span(line.s, 2, &len2);
        const char *f3 = get_field_span(line.s, 3, &len3);
        int64_t start = 0, stop = 0;

        if (!f0 || !f1 || !f2) {
            snprintf(err, errlen, "duckhts_mosdepth: BED line has fewer than 3 tab-delimited fields");
            goto cleanup;
        }
        if (!parse_int64_span_local(f1, len1, &start) ||
            !parse_int64_span_local(f2, len2, &stop) ||
            start < 0 || stop < start) {
            snprintf(err, errlen, "duckhts_mosdepth: invalid BED coordinates");
            goto cleanup;
        }

        char *chrom = (char *)malloc((size_t)len0 + 1);
        if (!chrom) {
            snprintf(err, errlen, "duckhts_mosdepth: out of memory");
            goto cleanup;
        }
        memcpy(chrom, f0, (size_t)len0);
        chrom[len0] = '\0';
        int tid = sam_hdr_name2tid(hdr, chrom);
        free(chrom);
        if (tid < 0) continue;

        char *name = NULL;
        if (f3 && len3 > 0) {
            name = (char *)malloc((size_t)len3 + 1);
            if (!name) {
                snprintf(err, errlen, "duckhts_mosdepth: out of memory");
                goto cleanup;
            }
            memcpy(name, f3, (size_t)len3);
            name[len3] = '\0';
        }
        if (region_list_append(&per_tid[tid], start, stop, name) != 0) {
            free(name);
            snprintf(err, errlen, "duckhts_mosdepth: out of memory");
            goto cleanup;
        }
        free(name);
    }

    for (int tid = 0; tid < sam_hdr_nref(hdr); tid++) {
        if (per_tid[tid].count > 1) {
            qsort(per_tid[tid].data, per_tid[tid].count, sizeof(mosdepth_region_t), region_cmp);
        }
    }

    rc = 0;

cleanup:
    free(line.s);
    if (fp) hts_close(fp);
    return rc;
}

static int remove_output_if_needed(const char *path) {
    if (!path) return 0;
    if (file_exists(path) && unlink(path) != 0) return -1;
    char *csi = append_csi_suffix(path);
    if (!csi) return -1;
    if (file_exists(csi) && unlink(csi) != 0) {
        duckdb_free(csi);
        return -1;
    }
    duckdb_free(csi);
    return 0;
}

static int ensure_output_available(const char *path, int overwrite, char *err, size_t errlen) {
    if (!path) return 0;
    char *csi = append_csi_suffix(path);
    if (!csi) {
        snprintf(err, errlen, "duckhts_mosdepth: out of memory");
        return -1;
    }
    if (!overwrite && (file_exists(path) || file_exists(csi))) {
        snprintf(err, errlen,
                 "duckhts_mosdepth: output '%s' already exists (use overwrite := TRUE to replace)",
                 path);
        duckdb_free(csi);
        return -1;
    }
    duckdb_free(csi);
    if (overwrite && remove_output_if_needed(path) != 0) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to replace existing output '%s'", path);
        return -1;
    }
    return 0;
}

static int build_bed_csi(const char *path, char *err, size_t errlen) {
    if (!path) return 0;
    if (tbx_index_build3(path, NULL, MOSDEPTH_CSI_MIN_SHIFT, 1, &tbx_conf_bed) != 0) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to build CSI index for '%s'", path);
        return -1;
    }
    return 0;
}

static int open_text_or_error(FILE **fh, const char *path, char *err, size_t errlen) {
    *fh = fopen(path, "wb");
    if (!*fh) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to open '%s': %s", path, strerror(errno));
        return -1;
    }
    return 0;
}

static int open_bgzf_or_error(BGZF **fp, const char *path, char *err, size_t errlen) {
    *fp = bgzf_open(path, "w1");
    if (!*fp) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to open '%s' for BGZF output", path);
        return -1;
    }
    return 0;
}

/* Fast BED-line helpers using kstring.h direct formatters instead of ksprintf.
 * kputs/kputc/kputll/kputw avoid printf format parsing overhead, which is
 * significant when emitting hundreds of millions of BED lines. */
static inline int bed_put_chrom_start_end(kstring_t *s, const char *chrom,
                                          int64_t start, int64_t end) {
    return (kputs(chrom, s) < 0 || kputc('\t', s) < 0 ||
            kputll(start, s) < 0 || kputc('\t', s) < 0 ||
            kputll(end, s) < 0) ? -1 : 0;
}

static inline int bed_put_precision_float(kstring_t *s, int precision, double val) {
    char buf[64];
    int n = snprintf(buf, sizeof(buf), "%.*f", precision, val);
    if (n < 0) return -1;
    return kputsn(buf, (size_t)n, s);
}

static int write_zero_per_base(BGZF *fp, kstring_t *line, const char *chrom, int64_t len) {
    line->l = 0;
    if (kputs(chrom, line) < 0 || kputsn("\t0\t", 3, line) < 0 ||
        kputll(len, line) < 0 || kputsn("\t0", 2, line) < 0) return -1;
    return bgzf_write_line(fp, line);
}

static int write_thresholds_header(BGZF *fp, const mosdepth_threshold_list_t *thresholds) {
    static const char header[] = "#chrom\tstart\tend\tregion";
    if (!fp || !thresholds || thresholds->count == 0) return 0;
    if (bgzf_write(fp, header, strlen(header)) < 0) return -1;
    for (size_t i = 0; i < thresholds->count; i++) {
        char buf[64];
        int len = snprintf(buf, sizeof(buf), "\t%" PRId64 "X", thresholds->data[i]);
        if (len < 0 || bgzf_write(fp, buf, (size_t)len) < 0) return -1;
    }
    if (bgzf_write(fp, "\n", 1) < 0) return -1;
    return 0;
}

static int write_thresholds_row(BGZF *fp, kstring_t *line, const char *chrom,
                                int64_t start, int64_t stop, const char *name,
                                const int64_t *counts,
                                const mosdepth_threshold_list_t *thresholds) {
    if (!fp || !line || !chrom || !thresholds || thresholds->count == 0) return 0;
    line->l = 0;
    if (bed_put_chrom_start_end(line, chrom, start, stop) < 0 ||
        kputc('\t', line) < 0 ||
        kputs((name && *name) ? name : "unknown", line) < 0) {
        return -1;
    }
    for (size_t i = 0; i < thresholds->count; i++) {
        if (kputc('\t', line) < 0 || kputll(counts ? counts[i] : 0, line) < 0) return -1;
    }
    return bgzf_write_line(fp, line);
}

static int write_per_base_rle(BGZF *fp, kstring_t *line, const char *chrom,
                              const int32_t *coverage, int64_t len) {
    if (len <= 0) return 0;
    int32_t last_depth = coverage[0];
    int64_t last_start = 0;
    for (int64_t i = 1; i < len; i++) {
        int32_t depth = coverage[i];
        if (depth == last_depth) continue;
        line->l = 0;
        if (bed_put_chrom_start_end(line, chrom, last_start, i) < 0 ||
            kputc('\t', line) < 0 || kputw((int)last_depth, line) < 0) {
            return -1;
        }
        if (bgzf_write_line(fp, line) != 0) return -1;
        last_start = i;
        last_depth = depth;
    }
    line->l = 0;
    if (bed_put_chrom_start_end(line, chrom, last_start, len) < 0 ||
        kputc('\t', line) < 0 || kputw((int)last_depth, line) < 0) {
        return -1;
    }
    return bgzf_write_line(fp, line);
}

static int write_quantized_rle(BGZF *fp, kstring_t *line, const char *chrom,
                               const int32_t *coverage, int64_t len,
                               const mosdepth_quantize_t *quantize) {
    int last_bucket;
    int64_t last_start;

    if (!fp || !line || !chrom || !coverage || len <= 0 || !quantize || quantize->count < 2) {
        return 0;
    }

    last_bucket = quantize_bucket(quantize, coverage[0]);
    last_start = 0;
    for (int64_t i = 1; i < len; i++) {
        int bucket = quantize_bucket(quantize, coverage[i]);
        if (bucket == last_bucket) continue;
        if (last_bucket >= 0 && (size_t)last_bucket + 1 < quantize->count) {
            line->l = 0;
            if (bed_put_chrom_start_end(line, chrom, last_start, i) < 0 ||
                kputc('\t', line) < 0 ||
                kputs(quantize->labels[last_bucket], line) < 0) {
                return -1;
            }
            if (bgzf_write_line(fp, line) != 0) return -1;
        }
        last_bucket = bucket;
        last_start = i;
    }

    if (last_bucket >= 0 && (size_t)last_bucket + 1 < quantize->count && last_start < len) {
        line->l = 0;
        if (bed_put_chrom_start_end(line, chrom, last_start, len) < 0 ||
            kputc('\t', line) < 0 ||
            kputs(quantize->labels[last_bucket], line) < 0) {
            return -1;
        }
        if (bgzf_write_line(fp, line) != 0) return -1;
    }
    return 0;
}

static int write_window_regions(BGZF *fp, kstring_t *line, const char *chrom,
                                const int32_t *coverage, int64_t len, int64_t window,
                                mosdepth_stat_t *region_stat,
                                mosdepth_dist_t *region_dist, int precision,
                                int has_reads, int use_median,
                                BGZF *thresholds_fp,
                                const mosdepth_threshold_list_t *thresholds) {
    int64_t *threshold_counts = NULL;
    mosdepth_count_stat_t median_stat = {0};
    int rc = -1;

    if (use_median && count_stat_init(&median_stat, MOSDEPTH_MEDIAN_COUNTS) != 0) return -1;
    if (thresholds_fp && thresholds && thresholds->count > 0) {
        threshold_counts = (int64_t *)calloc(thresholds->count, sizeof(int64_t));
        if (!threshold_counts) goto cleanup;
    }
    for (int64_t start = 0; start < len; start += window) {
        int64_t stop = start + window;
        if (stop > len) stop = len;
        double region_value = 0.0;
        if (threshold_counts) memset(threshold_counts, 0, thresholds->count * sizeof(int64_t));
        if (use_median) count_stat_clear(&median_stat);
        if (has_reads) {
            for (int64_t i = start; i < stop; i++) {
                uint32_t depth = (uint32_t)coverage[i];
                if (use_median) count_stat_add(&median_stat, depth);
                else region_value += (double)depth / (double)(stop - start);
                stat_add_depth(region_stat, depth);
                if (threshold_counts) {
                    for (size_t j = 0; j < thresholds->count; j++) {
                        if ((int64_t)depth < thresholds->data[j]) break;
                        threshold_counts[j] += 1;
                    }
                }
            }
            if (use_median) region_value = (double)count_stat_median(&median_stat);
            if (dist_add_count(region_dist, mosdepth_window_region_bucket(region_value), 1) != 0) goto cleanup;
        }
        line->l = 0;
        if (bed_put_chrom_start_end(line, chrom, start, stop) < 0 ||
            kputc('\t', line) < 0 ||
            bed_put_precision_float(line, precision, region_value) < 0) goto cleanup;
        if (bgzf_write_line(fp, line) != 0) goto cleanup;
        if (threshold_counts &&
            write_thresholds_row(thresholds_fp, line, chrom, start, stop, NULL,
                                 threshold_counts, thresholds) != 0) goto cleanup;
    }
    rc = 0;

cleanup:
    free(threshold_counts);
    count_stat_destroy(&median_stat);
    return rc;
}

static int write_bed_regions(BGZF *fp, kstring_t *line, const char *chrom,
                             const int32_t *coverage, int64_t len,
                             const mosdepth_region_list_t *regions,
                             mosdepth_stat_t *region_stat,
                             mosdepth_dist_t *region_dist, int precision,
                             int has_reads, int use_median,
                             BGZF *thresholds_fp,
                             const mosdepth_threshold_list_t *thresholds) {
    int64_t *threshold_counts = NULL;
    mosdepth_count_stat_t median_stat = {0};
    int rc = -1;

    if (!regions) return 0;
    if (use_median && count_stat_init(&median_stat, MOSDEPTH_MEDIAN_COUNTS) != 0) return -1;
    if (thresholds_fp && thresholds && thresholds->count > 0) {
        threshold_counts = (int64_t *)calloc(thresholds->count, sizeof(int64_t));
        if (!threshold_counts) goto cleanup;
    }
    for (size_t i = 0; i < regions->count; i++) {
        const mosdepth_region_t *r = &regions->data[i];
        int64_t in_start = r->start < 0 ? 0 : r->start;
        int64_t in_stop = r->stop > len ? len : r->stop;
        double region_value = 0.0;
        if (threshold_counts) memset(threshold_counts, 0, thresholds->count * sizeof(int64_t));
        if (use_median) count_stat_clear(&median_stat);
        if (has_reads && in_start < in_stop) {
            for (int64_t pos = in_start; pos < in_stop; pos++) {
                uint32_t depth = (uint32_t)coverage[pos];
                if (use_median) count_stat_add(&median_stat, depth);
                else region_value += (double)depth / (double)(r->stop - r->start);
                stat_add_depth(region_stat, depth);
                if (dist_add_count(region_dist, depth_bucket_value(depth), 1) != 0) {
                    goto cleanup;
                }
                if (threshold_counts) {
                    for (size_t j = 0; j < thresholds->count; j++) {
                        if ((int64_t)depth < thresholds->data[j]) break;
                        threshold_counts[j] += 1;
                    }
                }
            }
            if (use_median) region_value = (double)count_stat_median(&median_stat);
        }
        line->l = 0;
        if (r->name && *r->name) {
            if (bed_put_chrom_start_end(line, chrom, r->start, r->stop) < 0 ||
                kputc('\t', line) < 0 || kputs(r->name, line) < 0 ||
                kputc('\t', line) < 0 ||
                bed_put_precision_float(line, precision, region_value) < 0) goto cleanup;
        } else {
            if (bed_put_chrom_start_end(line, chrom, r->start, r->stop) < 0 ||
                kputc('\t', line) < 0 ||
                bed_put_precision_float(line, precision, region_value) < 0) goto cleanup;
        }
        if (bgzf_write_line(fp, line) != 0) goto cleanup;
        if (threshold_counts &&
            write_thresholds_row(thresholds_fp, line, chrom, r->start, r->stop, r->name,
                                 threshold_counts, thresholds) != 0) goto cleanup;
    }
    rc = 0;

cleanup:
    free(threshold_counts);
    count_stat_destroy(&median_stat);
    return rc;
}

/* ─── Parallel contig processing ─────────────────────────────────────── */

typedef struct {
    /* Shared read-only configuration */
    const mosdepth_bind_t *bind;
    sam_hdr_t *hdr;
    mosdepth_region_list_t *region_lists;
    khash_t(mdstr) *read_group_set;
    int by_is_window;
    int64_t window;
    int n_targets;
    int precision;

    /* Work-stealing counter (atomic via __sync builtins) */
    volatile int next_claim;

    /* Ordered output assembly */
    pthread_mutex_t write_mutex;
    pthread_cond_t write_cond;
    int next_write_tid;

    /* Error state */
    volatile int failed;
    char error[512];

    /* Output handles — written under write_mutex */
    FILE *fh_summary;
    FILE *fh_global;
    FILE *fh_region;
    BGZF *bgzf_per_base;
    BGZF *bgzf_regions;
    BGZF *bgzf_quantized;
    BGZF *bgzf_thresholds;
    int summary_header_written;

    /* Accumulated totals — merged under write_mutex */
    mosdepth_dist_t total_global_dist;
    mosdepth_dist_t total_region_dist;
    mosdepth_stat_t total_stat;
    mosdepth_stat_t total_region_stat;
} mosdepth_parallel_ctx_t;

static void mosdepth_parallel_set_error(mosdepth_parallel_ctx_t *ctx, const char *msg) {
    if (__sync_bool_compare_and_swap(&ctx->failed, 0, 1)) {
        snprintf(ctx->error, sizeof(ctx->error), "%s", msg);
    }
}

static void *mosdepth_parallel_worker(void *arg) {
    mosdepth_parallel_ctx_t *ctx = (mosdepth_parallel_ctx_t *)arg;
    const mosdepth_bind_t *bind = ctx->bind;
    samFile *fp = NULL;
    sam_hdr_t *hdr = NULL;
    hts_idx_t *idx = NULL;
    bam1_t *rec = NULL;
    int32_t *coverage = NULL;
    size_t coverage_cap = 0;
    kstring_t line = {0, 0, NULL};

    /* Each worker opens its own BAM file handle for independent I/O */
    fp = sam_open(bind->path, "r");
    if (!fp) {
        mosdepth_parallel_set_error(ctx, "duckhts_mosdepth: worker failed to open alignment file");
        goto worker_done;
    }
    duckhts_apply_remote_hts_tuning(fp, bind->path,
                                    DUCKHTS_HTS_IO_PROFILE_STREAMING);
    if (bind->fasta && hts_set_opt(fp, CRAM_OPT_REFERENCE, bind->fasta) < 0) {
        mosdepth_parallel_set_error(ctx, "duckhts_mosdepth: worker failed to set CRAM reference");
        goto worker_done;
    }
    if (bind->threads > 0 && hts_set_threads(fp, (int)bind->threads) < 0) {
        mosdepth_parallel_set_error(ctx, "duckhts_mosdepth: worker failed to set decompression threads");
        goto worker_done;
    }
    hdr = sam_hdr_read(fp);
    if (!hdr) {
        mosdepth_parallel_set_error(ctx, "duckhts_mosdepth: worker failed to read header");
        goto worker_done;
    }
    idx = sam_index_load3(fp, bind->path, bind->index_path,
                          HTS_IDX_SILENT_FAIL | HTS_IDX_SAVE_REMOTE);
    if (!idx) {
        mosdepth_parallel_set_error(ctx, "duckhts_mosdepth: worker failed to load index");
        goto worker_done;
    }
    rec = bam_init1();
    if (!rec) {
        mosdepth_parallel_set_error(ctx, "duckhts_mosdepth: worker out of memory");
        goto worker_done;
    }

    for (;;) {
        int tid = __sync_fetch_and_add(&ctx->next_claim, 1);
        const char *chrom;
        int64_t chrom_len;
        int skipped = 0;
        int saw_record = 0;
        int write_rc = 0;
        hts_itr_t *itr = NULL;
        khash_t(mdmate) *seen = NULL;
        mosdepth_dist_t chrom_global_dist = {0};
        mosdepth_dist_t chrom_region_dist = {0};
        mosdepth_stat_t chrom_stat;
        mosdepth_stat_t chrom_region_stat;
        char local_err[512];

        if (tid >= ctx->n_targets) break;
        if (ctx->failed) break;

        chrom = sam_hdr_tid2name(ctx->hdr, tid);
        chrom_len = (int64_t)sam_hdr_tid2len(ctx->hdr, tid);
        stat_clear(&chrom_stat);
        stat_clear(&chrom_region_stat);

        /* Skip checks — must match sequential path exactly */
        if (bind->chrom && bind->chrom[0] != '\0' &&
            strcmp(bind->chrom, chrom) != 0) {
            skipped = 1;
        } else if (bind->no_per_base && bind->regions_path &&
                   !ctx->by_is_window &&
                   (!ctx->region_lists ||
                    ctx->region_lists[tid].count == 0)) {
            skipped = 1;
        } else if (chrom_len <= 0) {
            skipped = 1;
        }

        if (!skipped) {
            /* Allocate/reuse coverage array */
            if (coverage_cap < (size_t)chrom_len + 1) {
                int32_t *new_cov = (int32_t *)realloc(
                    coverage, ((size_t)chrom_len + 1) * sizeof(int32_t));
                if (!new_cov) {
                    snprintf(local_err, sizeof(local_err),
                             "duckhts_mosdepth: worker OOM for '%s'", chrom);
                    mosdepth_parallel_set_error(ctx, local_err);
                    goto contig_error;
                }
                coverage = new_cov;
                coverage_cap = (size_t)chrom_len + 1;
            }
            memset(coverage, 0,
                   ((size_t)chrom_len + 1) * sizeof(int32_t));

            /* Mate tracking for default mode */
            if (!bind->fast_mode && !bind->fragment_mode) {
                seen = kh_init(mdmate);
                if (!seen) {
                    mosdepth_parallel_set_error(
                        ctx, "duckhts_mosdepth: worker OOM");
                    goto contig_error;
                }
            }

            /* Read BAM records for this contig */
            itr = sam_itr_queryi(idx, tid, 0, chrom_len);
            if (itr) {
                int iter_rc;
                while ((iter_rc = sam_itr_next(fp, itr, rec)) >= 0) {
                    hts_pos_t start = rec->core.pos;
                    hts_pos_t stop = bam_endpos(rec);
                    int64_t abs_isize = rec->core.isize < 0
                                            ? -(int64_t)rec->core.isize
                                            : (int64_t)rec->core.isize;

                    saw_record = 1;
                    if ((int64_t)rec->core.qual < bind->mapq) continue;
                    if (bind->min_frag_len >= 0 &&
                        abs_isize < bind->min_frag_len)
                        continue;
                    if (bind->max_frag_len >= 0 &&
                        abs_isize > bind->max_frag_len)
                        continue;
                    if (((uint16_t)rec->core.flag &
                         (uint16_t)bind->flag) != 0)
                        continue;
                    if (bind->include_flag != 0 &&
                        (((uint16_t)rec->core.flag &
                          (uint16_t)bind->include_flag) == 0))
                        continue;
                    if (!record_in_read_groups(rec, ctx->read_group_set))
                        continue;

                    if (bind->fragment_mode) {
                        hts_pos_t fragment_start, fragment_stop;
                        if ((rec->core.flag & BAM_FREAD2) ||
                            !(rec->core.flag & BAM_FPROPER_PAIR) ||
                            (rec->core.flag & BAM_FSUPPLEMENTARY))
                            continue;
                        fragment_start = start < rec->core.mpos
                                             ? start
                                             : rec->core.mpos;
                        fragment_stop = fragment_start + abs_isize;
                        add_coverage_delta(coverage, chrom_len,
                                           fragment_start, fragment_stop,
                                           1);
                        continue;
                    }

                    if (!bind->fast_mode) {
                        if (handle_default_overlap(rec, seen, coverage,
                                                   chrom_len) != 0) {
                            mosdepth_parallel_set_error(
                                ctx,
                                "duckhts_mosdepth: worker OOM in overlap");
                            hts_itr_destroy(itr);
                            if (seen) seen_mates_destroy(seen);
                            goto contig_error;
                        }
                        add_record_cigar_coverage(rec, coverage, chrom_len);
                        continue;
                    }

                    add_coverage_delta(coverage, chrom_len, start, stop, 1);
                }
                if (iter_rc < -1) {
                    if (fp->format.format == cram && !bind->fasta) {
                        snprintf(local_err, sizeof(local_err),
                                 "duckhts_mosdepth: CRAM decode failed "
                                 "scanning '%s'; provide fasta",
                                 chrom);
                    } else {
                        snprintf(local_err, sizeof(local_err),
                                 "duckhts_mosdepth: worker failed scanning "
                                 "'%s'",
                                 chrom);
                    }
                    mosdepth_parallel_set_error(ctx, local_err);
                    hts_itr_destroy(itr);
                    if (seen) seen_mates_destroy(seen);
                    goto contig_error;
                }
                hts_itr_destroy(itr);
                itr = NULL;
            }
            if (seen) {
                seen_mates_destroy(seen);
                seen = NULL;
            }

            /* Prefix-sum and stats (only when records were seen) */
            if (saw_record) {
                for (int64_t i = 1; i <= chrom_len; i++)
                    coverage[i] += coverage[i - 1];
                for (int64_t i = 0; i < chrom_len; i++) {
                    uint32_t depth = (uint32_t)coverage[i];
                    stat_add_depth(&chrom_stat, depth);
                    if (dist_add_count(&chrom_global_dist,
                                       depth_bucket_value(depth), 1) != 0) {
                        mosdepth_parallel_set_error(
                            ctx, "duckhts_mosdepth: worker OOM in dist");
                        goto contig_error;
                    }
                }
            }
        }

        /* ── Ordered write phase ────────────────────────────────── */
        pthread_mutex_lock(&ctx->write_mutex);
        while (ctx->next_write_tid < tid && !ctx->failed)
            pthread_cond_wait(&ctx->write_cond, &ctx->write_mutex);

        if (ctx->failed) {
            pthread_mutex_unlock(&ctx->write_mutex);
            dist_destroy(&chrom_global_dist);
            dist_destroy(&chrom_region_dist);
            break;
        }

        if (!skipped && !saw_record) {
            /* Zero-coverage output for contigs with no reads */
            if (ctx->bgzf_per_base &&
                write_zero_per_base(ctx->bgzf_per_base, &line, chrom,
                                    chrom_len) != 0) {
                write_rc = -1;
            }
            if (!write_rc && ctx->bgzf_quantized &&
                bind->quantize.count > 1 &&
                bind->quantize.cuts[0] == 0) {
                line.l = 0;
                if (bed_put_chrom_start_end(&line, chrom, 0, chrom_len) <
                        0 ||
                    kputc('\t', &line) < 0 ||
                    kputs(bind->quantize.labels[0], &line) < 0 ||
                    bgzf_write_line(ctx->bgzf_quantized, &line) != 0) {
                    write_rc = -1;
                }
            }
            if (!write_rc && ctx->bgzf_regions) {
                if (ctx->by_is_window) {
                    if (write_window_regions(
                            ctx->bgzf_regions, &line, chrom, NULL,
                            chrom_len, ctx->window,
                            &chrom_region_stat, &chrom_region_dist,
                            ctx->precision, 0, bind->use_median,
                            ctx->bgzf_thresholds,
                            &bind->thresholds) != 0)
                        write_rc = -1;
                } else if (ctx->region_lists) {
                    if (write_bed_regions(
                            ctx->bgzf_regions, &line, chrom, NULL,
                            chrom_len, &ctx->region_lists[tid],
                            &chrom_region_stat, &chrom_region_dist,
                            ctx->precision, 0, bind->use_median,
                            ctx->bgzf_thresholds,
                            &bind->thresholds) != 0)
                        write_rc = -1;
                }
            }
        } else if (!skipped && saw_record) {
            /* Normal output for contigs with reads */
            if (write_summary(ctx->fh_summary,
                              &ctx->summary_header_written, chrom,
                              &chrom_stat, ctx->precision) != 0)
                write_rc = -1;
            if (!write_rc &&
                write_distribution(ctx->fh_global, chrom,
                                   &chrom_global_dist,
                                   ctx->precision) != 0)
                write_rc = -1;
            if (!write_rc &&
                dist_sum_into(&chrom_global_dist,
                              &ctx->total_global_dist) != 0)
                write_rc = -1;
            if (!write_rc) stat_merge(&ctx->total_stat, &chrom_stat);

            if (!write_rc && ctx->bgzf_per_base &&
                write_per_base_rle(ctx->bgzf_per_base, &line, chrom,
                                   coverage, chrom_len) != 0)
                write_rc = -1;
            if (!write_rc && ctx->bgzf_quantized &&
                write_quantized_rle(ctx->bgzf_quantized, &line, chrom,
                                    coverage, chrom_len,
                                    &bind->quantize) != 0)
                write_rc = -1;

            if (!write_rc && ctx->bgzf_regions) {
                if (ctx->by_is_window) {
                    if (write_window_regions(
                            ctx->bgzf_regions, &line, chrom, coverage,
                            chrom_len, ctx->window,
                            &chrom_region_stat, &chrom_region_dist,
                            ctx->precision, 1, bind->use_median,
                            ctx->bgzf_thresholds,
                            &bind->thresholds) != 0)
                        write_rc = -1;
                } else if (ctx->region_lists) {
                    if (write_bed_regions(
                            ctx->bgzf_regions, &line, chrom, coverage,
                            chrom_len, &ctx->region_lists[tid],
                            &chrom_region_stat, &chrom_region_dist,
                            ctx->precision, 1, bind->use_median,
                            ctx->bgzf_thresholds,
                            &bind->thresholds) != 0)
                        write_rc = -1;
                }
                if (!write_rc && chrom_region_stat.cum_length > 0) {
                    char *label = append_suffix(chrom, "_region");
                    if (!label) {
                        write_rc = -1;
                    } else {
                        if (write_summary(
                                ctx->fh_summary,
                                &ctx->summary_header_written, label,
                                &chrom_region_stat,
                                ctx->precision) != 0)
                            write_rc = -1;
                        if (!write_rc) duckdb_free(label);
                        else {
                            duckdb_free(label);
                            /* write_rc already -1 */
                        }
                    }
                    if (!write_rc &&
                        write_distribution(ctx->fh_region, chrom,
                                           &chrom_region_dist,
                                           ctx->precision) != 0)
                        write_rc = -1;
                    if (!write_rc &&
                        dist_sum_into(&chrom_region_dist,
                                      &ctx->total_region_dist) != 0)
                        write_rc = -1;
                    if (!write_rc)
                        stat_merge(&ctx->total_region_stat,
                                   &chrom_region_stat);
                }
            }
        }

        if (write_rc != 0) {
            snprintf(local_err, sizeof(local_err),
                     "duckhts_mosdepth: worker write failure for '%s'",
                     chrom);
            mosdepth_parallel_set_error(ctx, local_err);
        }

        ctx->next_write_tid = tid + 1;
        pthread_cond_broadcast(&ctx->write_cond);
        pthread_mutex_unlock(&ctx->write_mutex);

        dist_destroy(&chrom_global_dist);
        dist_destroy(&chrom_region_dist);
        if (ctx->failed) break;
        continue;

    contig_error:
        /* Error already set via mosdepth_parallel_set_error.
         * Broadcast to unblock any workers waiting for their write turn. */
        pthread_mutex_lock(&ctx->write_mutex);
        pthread_cond_broadcast(&ctx->write_cond);
        pthread_mutex_unlock(&ctx->write_mutex);
        dist_destroy(&chrom_global_dist);
        dist_destroy(&chrom_region_dist);
        break;
    }

worker_done:
    free(coverage);
    free(line.s);
    if (rec) bam_destroy1(rec);
    if (idx) hts_idx_destroy(idx);
    if (hdr) sam_hdr_destroy(hdr);
    if (fp) sam_close(fp);
    return NULL;
}

static int run_duckhts_mosdepth(mosdepth_bind_t *bind, char *err, size_t errlen) {
    samFile *fp = NULL;
    sam_hdr_t *hdr = NULL;
    hts_idx_t *idx = NULL;
    bam1_t *rec = NULL;
    FILE *fh_summary = NULL;
    FILE *fh_global = NULL;
    FILE *fh_region = NULL;
    BGZF *bgzf_per_base = NULL;
    BGZF *bgzf_regions = NULL;
    BGZF *bgzf_quantized = NULL;
    BGZF *bgzf_thresholds = NULL;
    kstring_t line = {0, 0, NULL};
    mosdepth_region_list_t *region_lists = NULL;
    khash_t(mdstr) *read_group_set = NULL;
    mosdepth_dist_t total_global_dist = {0};
    mosdepth_dist_t total_region_dist = {0};
    mosdepth_stat_t total_stat;
    mosdepth_stat_t total_region_stat;
    int32_t *coverage = NULL;
    size_t coverage_cap = 0;
    int summary_header_written = 0;
    int precision = (int)bind->precision;
    int rc = -1;
    int n_targets = 0;
    int by_is_window = 0;
    int64_t window = 0;

    stat_clear(&total_stat);
    stat_clear(&total_region_stat);

    fp = sam_open(bind->path, "r");
    if (!fp) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to open alignment file '%s'", bind->path);
        goto cleanup;
    }
    duckhts_apply_remote_hts_tuning(fp, bind->path,
                                    DUCKHTS_HTS_IO_PROFILE_STREAMING);
    if (bind->fasta && hts_set_opt(fp, CRAM_OPT_REFERENCE, bind->fasta) < 0) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to set CRAM reference '%s'", bind->fasta);
        goto cleanup;
    }
    if (bind->threads > 0 && hts_set_threads(fp, (int)bind->threads) < 0) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to configure decompression threads");
        goto cleanup;
    }

    hdr = sam_hdr_read(fp);
    if (!hdr) {
        if (fp->format.format == cram && !bind->fasta) {
            snprintf(err, errlen,
                     "duckhts_mosdepth: failed to read CRAM header; provide fasta := '<reference.fa>' when required");
        } else {
            snprintf(err, errlen, "duckhts_mosdepth: failed to read alignment header");
        }
        goto cleanup;
    }
    idx = sam_index_load3(fp, bind->path, bind->index_path, HTS_IDX_SILENT_FAIL | HTS_IDX_SAVE_REMOTE);
    if (!idx) {
        snprintf(err, errlen,
                 "duckhts_mosdepth: indexed BAM/CRAM input is required (failed to load .bai/.csi/.crai)");
        goto cleanup;
    }
    rec = bam_init1();
    if (!rec) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to allocate BAM record");
        goto cleanup;
    }

    n_targets = sam_hdr_nref(hdr);
    if (bind->chrom && bind->chrom[0] != '\0' && sam_hdr_name2tid(hdr, bind->chrom) < 0) {
        snprintf(err, errlen, "duckhts_mosdepth: chromosome '%s' not found in BAM header", bind->chrom);
        goto cleanup;
    }

    if (bind->by && bind->by[0] != '\0') {
        if (is_digits_only(bind->by)) {
            by_is_window = 1;
            window = strtoll(bind->by, NULL, 10);
            if (window <= 0) {
                snprintf(err, errlen, "duckhts_mosdepth: by must be a positive window size or a BED path");
                goto cleanup;
            }
        } else {
            region_lists = (mosdepth_region_list_t *)calloc((size_t)n_targets, sizeof(mosdepth_region_list_t));
            if (!region_lists) {
                snprintf(err, errlen, "duckhts_mosdepth: out of memory");
                goto cleanup;
            }
            if (parse_bed_regions(bind->by, hdr, region_lists, err, errlen) != 0) {
                goto cleanup;
            }
        }
    }

    if (parse_read_group_set(bind->read_groups, &read_group_set, err, errlen) != 0) {
        goto cleanup;
    }

    if (open_text_or_error(&fh_summary, bind->summary_path, err, errlen) != 0) goto cleanup;
    if (open_text_or_error(&fh_global, bind->global_dist_path, err, errlen) != 0) goto cleanup;
    if (bind->region_dist_path && open_text_or_error(&fh_region, bind->region_dist_path, err, errlen) != 0) {
        goto cleanup;
    }
    if (bind->per_base_path && open_bgzf_or_error(&bgzf_per_base, bind->per_base_path, err, errlen) != 0) {
        goto cleanup;
    }
    if (bind->regions_path && open_bgzf_or_error(&bgzf_regions, bind->regions_path, err, errlen) != 0) {
        goto cleanup;
    }
    if (bind->quantized_path &&
        open_bgzf_or_error(&bgzf_quantized, bind->quantized_path, err, errlen) != 0) {
        goto cleanup;
    }
    if (bind->thresholds_path && open_bgzf_or_error(&bgzf_thresholds, bind->thresholds_path, err, errlen) != 0) {
        goto cleanup;
    }
    if (bgzf_thresholds && write_thresholds_header(bgzf_thresholds, &bind->thresholds) != 0) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to write thresholds header");
        goto cleanup;
    }

    if (bind->processing_threads > 1) {
        /* ── Parallel contig processing ─────────────────────────── */
        mosdepth_parallel_ctx_t pctx;
        pthread_t *threads_arr = NULL;
        int n_workers = (int)bind->processing_threads;
        int i;

        memset(&pctx, 0, sizeof(pctx));
        pctx.bind = bind;
        pctx.hdr = hdr;
        pctx.region_lists = region_lists;
        pctx.read_group_set = read_group_set;
        pctx.by_is_window = by_is_window;
        pctx.window = window;
        pctx.n_targets = n_targets;
        pctx.precision = precision;
        pctx.fh_summary = fh_summary;
        pctx.fh_global = fh_global;
        pctx.fh_region = fh_region;
        pctx.bgzf_per_base = bgzf_per_base;
        pctx.bgzf_regions = bgzf_regions;
        pctx.bgzf_quantized = bgzf_quantized;
        pctx.bgzf_thresholds = bgzf_thresholds;
        stat_clear(&pctx.total_stat);
        stat_clear(&pctx.total_region_stat);
        pthread_mutex_init(&pctx.write_mutex, NULL);
        pthread_cond_init(&pctx.write_cond, NULL);

        threads_arr = (pthread_t *)malloc(sizeof(pthread_t) * (size_t)n_workers);
        if (!threads_arr) {
            snprintf(err, errlen, "duckhts_mosdepth: failed to allocate worker threads");
            pthread_mutex_destroy(&pctx.write_mutex);
            pthread_cond_destroy(&pctx.write_cond);
            goto cleanup;
        }

        for (i = 0; i < n_workers; i++) {
            if (pthread_create(&threads_arr[i], NULL,
                               mosdepth_parallel_worker, &pctx) != 0) {
                pctx.failed = 1;
                snprintf(pctx.error, sizeof(pctx.error),
                         "duckhts_mosdepth: failed to start worker thread");
                n_workers = i;
                break;
            }
        }
        for (i = 0; i < n_workers; i++) {
            pthread_join(threads_arr[i], NULL);
        }
        free(threads_arr);

        /* Copy parallel results back to local variables */
        total_global_dist = pctx.total_global_dist;
        total_region_dist = pctx.total_region_dist;
        total_stat = pctx.total_stat;
        total_region_stat = pctx.total_region_stat;
        summary_header_written = pctx.summary_header_written;

        /* Zero out pctx dists so they aren't double-freed */
        memset(&pctx.total_global_dist, 0, sizeof(mosdepth_dist_t));
        memset(&pctx.total_region_dist, 0, sizeof(mosdepth_dist_t));

        pthread_mutex_destroy(&pctx.write_mutex);
        pthread_cond_destroy(&pctx.write_cond);

        if (pctx.failed) {
            snprintf(err, errlen, "%s",
                     pctx.error[0] ? pctx.error
                                   : "duckhts_mosdepth: parallel worker failure");
            goto cleanup;
        }
    } else {
    /* ── Sequential contig processing (original path) ───────── */

    for (int tid = 0; tid < n_targets; tid++) {
        const char *chrom = sam_hdr_tid2name(hdr, tid);
        int64_t chrom_len = (int64_t)sam_hdr_tid2len(hdr, tid);
        hts_itr_t *itr = NULL;
        khash_t(mdmate) *seen = NULL;
        int saw_record = 0;
        int iter_rc = -1;
        mosdepth_dist_t chrom_global_dist = {0};
        mosdepth_dist_t chrom_region_dist = {0};
        mosdepth_stat_t chrom_stat;
        mosdepth_stat_t chrom_region_stat;

        if (bind->chrom && bind->chrom[0] != '\0' && strcmp(bind->chrom, chrom) != 0) continue;
        if (bind->no_per_base && bind->regions_path && !by_is_window &&
            (!region_lists || region_lists[tid].count == 0)) {
            continue;
        }
        if (chrom_len <= 0) continue;

        stat_clear(&chrom_stat);
        stat_clear(&chrom_region_stat);

        if (coverage_cap < (size_t)chrom_len + 1) {
            int32_t *new_coverage =
                (int32_t *)realloc(coverage, ((size_t)chrom_len + 1) * sizeof(int32_t));
            if (!new_coverage) {
                snprintf(err, errlen, "duckhts_mosdepth: failed to allocate coverage array for '%s'", chrom);
                goto contig_cleanup;
            }
            coverage = new_coverage;
            coverage_cap = (size_t)chrom_len + 1;
        }
        memset(coverage, 0, ((size_t)chrom_len + 1) * sizeof(int32_t));

        if (!bind->fast_mode && !bind->fragment_mode) {
            seen = kh_init(mdmate);
            if (!seen) {
                snprintf(err, errlen, "duckhts_mosdepth: out of memory");
                goto contig_cleanup;
            }
        }

        itr = sam_itr_queryi(idx, tid, 0, chrom_len);
        if (itr) {
            while ((iter_rc = sam_itr_next(fp, itr, rec)) >= 0) {
                hts_pos_t start = rec->core.pos;
                hts_pos_t stop = bam_endpos(rec);
                int64_t abs_isize = rec->core.isize < 0 ? -(int64_t)rec->core.isize : (int64_t)rec->core.isize;

                saw_record = 1;
                if ((int64_t)rec->core.qual < bind->mapq) continue;
                if (bind->min_frag_len >= 0 && abs_isize < bind->min_frag_len) continue;
                if (bind->max_frag_len >= 0 && abs_isize > bind->max_frag_len) continue;
                if (((uint16_t)rec->core.flag & (uint16_t)bind->flag) != 0) continue;
                if (bind->include_flag != 0 &&
                    (((uint16_t)rec->core.flag & (uint16_t)bind->include_flag) == 0)) {
                    continue;
                }
                if (!record_in_read_groups(rec, read_group_set)) continue;

                if (bind->fragment_mode) {
                    hts_pos_t fragment_start;
                    hts_pos_t fragment_stop;

                    if ((rec->core.flag & BAM_FREAD2) ||
                        !(rec->core.flag & BAM_FPROPER_PAIR) ||
                        (rec->core.flag & BAM_FSUPPLEMENTARY)) {
                        continue;
                    }
                    fragment_start = start < rec->core.mpos ? start : rec->core.mpos;
                    fragment_stop = fragment_start + abs_isize;
                    add_coverage_delta(coverage, chrom_len, fragment_start, fragment_stop, 1);
                    continue;
                }

                if (!bind->fast_mode) {
                    if (handle_default_overlap(rec, seen, coverage, chrom_len) != 0) {
                        snprintf(err, errlen, "duckhts_mosdepth: out of memory");
                        goto contig_cleanup;
                    }
                    add_record_cigar_coverage(rec, coverage, chrom_len);
                    continue;
                }

                add_coverage_delta(coverage, chrom_len, start, stop, 1);
            }
            if (iter_rc < -1) {
                if (fp->format.format == cram && !bind->fasta) {
                    snprintf(err, errlen,
                             "duckhts_mosdepth: CRAM decode failed while scanning '%s'; provide fasta := '<reference.fa>'",
                             chrom);
                } else {
                    snprintf(err, errlen, "duckhts_mosdepth: failed while scanning '%s'", chrom);
                }
                goto contig_cleanup;
            }
        }

        if (!saw_record) {
            if (bgzf_per_base && write_zero_per_base(bgzf_per_base, &line, chrom, chrom_len) != 0) {
                snprintf(err, errlen, "duckhts_mosdepth: failed to write per-base output");
                goto contig_cleanup;
            }
            if (bgzf_quantized && bind->quantize.count > 1 && bind->quantize.cuts[0] == 0) {
                line.l = 0;
                if (bed_put_chrom_start_end(&line, chrom, 0, chrom_len) < 0 ||
                    kputc('\t', &line) < 0 ||
                    kputs(bind->quantize.labels[0], &line) < 0 ||
                    bgzf_write_line(bgzf_quantized, &line) != 0) {
                    snprintf(err, errlen, "duckhts_mosdepth: failed to write quantized output");
                    goto contig_cleanup;
                }
            }
            if (bgzf_regions) {
                if (by_is_window) {
                    if (write_window_regions(bgzf_regions, &line, chrom, NULL, chrom_len, window,
                                             &chrom_region_stat, &chrom_region_dist,
                                             precision, 0, bind->use_median, bgzf_thresholds,
                                             &bind->thresholds) != 0) {
                        snprintf(err, errlen, "duckhts_mosdepth: failed to write region output");
                        goto contig_cleanup;
                    }
                } else if (region_lists) {
                    if (write_bed_regions(bgzf_regions, &line, chrom, NULL, chrom_len, &region_lists[tid],
                                          &chrom_region_stat, &chrom_region_dist,
                                          precision, 0, bind->use_median, bgzf_thresholds,
                                          &bind->thresholds) != 0) {
                        snprintf(err, errlen, "duckhts_mosdepth: failed to write region output");
                        goto contig_cleanup;
                    }
                }
            }
            rc = 0;
            goto contig_cleanup;
        }

        for (int64_t i = 1; i <= chrom_len; i++) coverage[i] += coverage[i - 1];

        for (int64_t i = 0; i < chrom_len; i++) {
            uint32_t depth = (uint32_t)coverage[i];
            stat_add_depth(&chrom_stat, depth);
            if (dist_add_count(&chrom_global_dist, depth_bucket_value(depth), 1) != 0) {
                snprintf(err, errlen, "duckhts_mosdepth: out of memory");
                goto contig_cleanup;
            }
        }
        if (write_summary(fh_summary, &summary_header_written, chrom, &chrom_stat, precision) != 0) {
            snprintf(err, errlen, "duckhts_mosdepth: failed to write summary");
            goto contig_cleanup;
        }
        if (write_distribution(fh_global, chrom, &chrom_global_dist, precision) != 0) {
            snprintf(err, errlen, "duckhts_mosdepth: failed to write global distribution");
            goto contig_cleanup;
        }
        if (dist_sum_into(&chrom_global_dist, &total_global_dist) != 0) {
            snprintf(err, errlen, "duckhts_mosdepth: out of memory");
            goto contig_cleanup;
        }
        stat_merge(&total_stat, &chrom_stat);

        if (bgzf_per_base && write_per_base_rle(bgzf_per_base, &line, chrom, coverage, chrom_len) != 0) {
            snprintf(err, errlen, "duckhts_mosdepth: failed to write per-base output");
            goto contig_cleanup;
        }
        if (bgzf_quantized &&
            write_quantized_rle(bgzf_quantized, &line, chrom, coverage, chrom_len, &bind->quantize) != 0) {
            snprintf(err, errlen, "duckhts_mosdepth: failed to write quantized output");
            goto contig_cleanup;
        }

        if (bgzf_regions) {
            if (by_is_window) {
                if (write_window_regions(bgzf_regions, &line, chrom, coverage, chrom_len, window,
                                         &chrom_region_stat, &chrom_region_dist,
                                         precision, 1, bind->use_median, bgzf_thresholds,
                                         &bind->thresholds) != 0) {
                    snprintf(err, errlen, "duckhts_mosdepth: failed to write region output");
                    goto contig_cleanup;
                }
            } else if (region_lists) {
                if (write_bed_regions(bgzf_regions, &line, chrom, coverage, chrom_len, &region_lists[tid],
                                      &chrom_region_stat, &chrom_region_dist,
                                      precision, 1, bind->use_median, bgzf_thresholds,
                                      &bind->thresholds) != 0) {
                    snprintf(err, errlen, "duckhts_mosdepth: failed to write region output");
                    goto contig_cleanup;
                }
            }
            if (chrom_region_stat.cum_length > 0) {
                char *label = append_suffix(chrom, "_region");
                if (!label) {
                    snprintf(err, errlen, "duckhts_mosdepth: out of memory");
                    goto contig_cleanup;
                }
                if (write_summary(fh_summary, &summary_header_written, label, &chrom_region_stat, precision) != 0) {
                    duckdb_free(label);
                    snprintf(err, errlen, "duckhts_mosdepth: failed to write region summary");
                    goto contig_cleanup;
                }
                duckdb_free(label);
                if (write_distribution(fh_region, chrom, &chrom_region_dist, precision) != 0) {
                    snprintf(err, errlen, "duckhts_mosdepth: failed to write region distribution");
                    goto contig_cleanup;
                }
                if (dist_sum_into(&chrom_region_dist, &total_region_dist) != 0) {
                    snprintf(err, errlen, "duckhts_mosdepth: out of memory");
                    goto contig_cleanup;
                }
                stat_merge(&total_region_stat, &chrom_region_stat);
            }
        }

        rc = 0;

contig_cleanup:
        if (itr) hts_itr_destroy(itr);
        if (seen) seen_mates_destroy(seen);
        dist_destroy(&chrom_global_dist);
        dist_destroy(&chrom_region_dist);
        if (rc != 0) goto cleanup;
    }

    } /* end sequential else block */

    if (write_summary(fh_summary, &summary_header_written, "total", &total_stat, precision) != 0) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to write total summary");
        goto cleanup;
    }
    if (write_distribution(fh_global, "total", &total_global_dist, precision) != 0) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to write total global distribution");
        goto cleanup;
    }
    if (fh_region) {
        if (write_summary(fh_summary, &summary_header_written, "total_region", &total_region_stat, precision) != 0) {
            snprintf(err, errlen, "duckhts_mosdepth: failed to write total region summary");
            goto cleanup;
        }
        if (write_distribution(fh_region, "total", &total_region_dist, precision) != 0) {
            snprintf(err, errlen, "duckhts_mosdepth: failed to write total region distribution");
            goto cleanup;
        }
    }

    rc = 0;

cleanup:
    free(line.s);
    if (bgzf_per_base) {
        if (bgzf_close(bgzf_per_base) != 0 && rc == 0) {
            snprintf(err, errlen, "duckhts_mosdepth: failed to close per-base output");
            rc = -1;
        }
    }
    if (bgzf_regions) {
        if (bgzf_close(bgzf_regions) != 0 && rc == 0) {
            snprintf(err, errlen, "duckhts_mosdepth: failed to close region output");
            rc = -1;
        }
    }
    if (bgzf_quantized) {
        if (bgzf_close(bgzf_quantized) != 0 && rc == 0) {
            snprintf(err, errlen, "duckhts_mosdepth: failed to close quantized output");
            rc = -1;
        }
    }
    if (bgzf_thresholds) {
        if (bgzf_close(bgzf_thresholds) != 0 && rc == 0) {
            snprintf(err, errlen, "duckhts_mosdepth: failed to close thresholds output");
            rc = -1;
        }
    }
    if (fh_summary && fclose(fh_summary) != 0 && rc == 0) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to close summary output");
        rc = -1;
    }
    if (fh_global && fclose(fh_global) != 0 && rc == 0) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to close global distribution output");
        rc = -1;
    }
    if (fh_region && fclose(fh_region) != 0 && rc == 0) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to close region distribution output");
        rc = -1;
    }
    if (rec) bam_destroy1(rec);
    if (idx) hts_idx_destroy(idx);
    if (hdr) sam_hdr_destroy(hdr);
    if (fp) sam_close(fp);
    free(coverage);
    if (region_lists) {
        for (int i = 0; i < n_targets; i++) region_list_destroy(&region_lists[i]);
        free(region_lists);
    }
    read_group_set_destroy(read_group_set);
    dist_destroy(&total_global_dist);
    dist_destroy(&total_region_dist);

    if (rc == 0) {
        if (bind->per_base_path && build_bed_csi(bind->per_base_path, err, errlen) != 0) rc = -1;
        if (rc == 0 && bind->regions_path && build_bed_csi(bind->regions_path, err, errlen) != 0) rc = -1;
        if (rc == 0 && bind->quantized_path && build_bed_csi(bind->quantized_path, err, errlen) != 0) rc = -1;
        if (rc == 0 && bind->thresholds_path && build_bed_csi(bind->thresholds_path, err, errlen) != 0) rc = -1;
    }
    return rc;
}

static void destroy_mosdepth_bind(void *data) {
    mosdepth_bind_t *bind = (mosdepth_bind_t *)data;
    if (!bind) return;
    if (bind->prefix) duckdb_free(bind->prefix);
    if (bind->path) duckdb_free(bind->path);
    if (bind->chrom) duckdb_free(bind->chrom);
    if (bind->by) duckdb_free(bind->by);
    if (bind->fasta) duckdb_free(bind->fasta);
    if (bind->read_groups) duckdb_free(bind->read_groups);
    if (bind->index_path) duckdb_free(bind->index_path);
    if (bind->summary_path) duckdb_free(bind->summary_path);
    if (bind->global_dist_path) duckdb_free(bind->global_dist_path);
    if (bind->per_base_path) duckdb_free(bind->per_base_path);
    if (bind->regions_path) duckdb_free(bind->regions_path);
    if (bind->region_dist_path) duckdb_free(bind->region_dist_path);
    if (bind->quantized_path) duckdb_free(bind->quantized_path);
    if (bind->thresholds_path) duckdb_free(bind->thresholds_path);
    threshold_list_destroy(&bind->thresholds);
    quantize_destroy(&bind->quantize);
    duckdb_free(bind);
}

static void add_result_columns(duckdb_bind_info info) {
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_bind_add_result_column(info, "success", bool_type);
    duckdb_bind_add_result_column(info, "prefix", varchar_type);
    duckdb_bind_add_result_column(info, "summary_path", varchar_type);
    duckdb_bind_add_result_column(info, "global_dist_path", varchar_type);
    duckdb_bind_add_result_column(info, "per_base_path", varchar_type);
    duckdb_bind_add_result_column(info, "regions_path", varchar_type);
    duckdb_bind_add_result_column(info, "region_dist_path", varchar_type);
    duckdb_bind_add_result_column(info, "quantized_path", varchar_type);
    duckdb_bind_add_result_column(info, "thresholds_path", varchar_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_logical_type(&varchar_type);
}

static void duckhts_mosdepth_bind(duckdb_bind_info info) {
    duckdb_value prefix_val = duckdb_bind_get_parameter(info, 0);
    duckdb_value path_val = duckdb_bind_get_parameter(info, 1);
    mosdepth_bind_t *bind = NULL;
    duckdb_value val;
    char err[512];

    memset(err, 0, sizeof(err));

    char *prefix = duckdb_get_varchar(prefix_val);
    char *path = duckdb_get_varchar(path_val);
    duckdb_destroy_value(&prefix_val);
    duckdb_destroy_value(&path_val);

    if (!prefix || prefix[0] == '\0') {
        duckdb_bind_set_error(info, "duckhts_mosdepth requires a non-empty prefix");
        if (prefix) duckdb_free(prefix);
        if (path) duckdb_free(path);
        return;
    }
    if (!path || path[0] == '\0') {
        duckdb_bind_set_error(info, "duckhts_mosdepth requires a BAM file path");
        if (prefix) duckdb_free(prefix);
        if (path) duckdb_free(path);
        return;
    }

    bind = (mosdepth_bind_t *)duckdb_malloc(sizeof(mosdepth_bind_t));
    if (!bind) {
        duckdb_bind_set_error(info, "duckhts_mosdepth: out of memory");
        duckdb_free(prefix);
        duckdb_free(path);
        return;
    }
    memset(bind, 0, sizeof(mosdepth_bind_t));
    bind->prefix = prefix;
    bind->path = path;
    bind->threads = 2;
    bind->processing_threads = 2;
    bind->flag = 1796;
    bind->include_flag = 0;
    bind->mapq = 0;
    bind->min_frag_len = -1;
    bind->max_frag_len = -1;
    bind->precision = MOSDEPTH_DEFAULT_PRECISION;
    bind->fast_mode = 0;
    bind->fragment_mode = 0;
    bind->use_median = 0;
    bind->no_per_base = 0;
    bind->overwrite = 0;

    val = duckdb_bind_get_named_parameter(info, "chrom");
    if (val && !duckdb_is_null_value(val)) bind->chrom = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "by");
    if (val && !duckdb_is_null_value(val)) bind->by = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "fasta");
    if (val && !duckdb_is_null_value(val)) bind->fasta = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "read_groups");
    if (val && !duckdb_is_null_value(val)) bind->read_groups = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "index_path");
    if (val && !duckdb_is_null_value(val)) bind->index_path = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "threads");
    if (val && !duckdb_is_null_value(val)) bind->threads = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "processing_threads");
    if (val && !duckdb_is_null_value(val)) bind->processing_threads = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "flag");
    if (val && !duckdb_is_null_value(val)) bind->flag = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "include_flag");
    if (val && !duckdb_is_null_value(val)) bind->include_flag = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "mapq");
    if (val && !duckdb_is_null_value(val)) bind->mapq = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "min_frag_len");
    if (val && !duckdb_is_null_value(val)) bind->min_frag_len = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "max_frag_len");
    if (val && !duckdb_is_null_value(val)) bind->max_frag_len = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "precision_digits");
    if (val && !duckdb_is_null_value(val)) bind->precision = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "thresholds");
    if (val && !duckdb_is_null_value(val)) {
        char *thresholds = duckdb_get_varchar(val);
        if (thresholds) {
            if (parse_thresholds_string(thresholds, &bind->thresholds, err, sizeof(err)) != 0) {
                duckdb_free(thresholds);
                duckdb_destroy_value(&val);
                destroy_mosdepth_bind(bind);
                duckdb_bind_set_error(info, err);
                return;
            }
            duckdb_free(thresholds);
        }
    }
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "quantize");
    if (val && !duckdb_is_null_value(val)) {
        char *quantize = duckdb_get_varchar(val);
        if (quantize) {
            if (parse_quantize_string(quantize, &bind->quantize, err, sizeof(err)) != 0) {
                duckdb_free(quantize);
                duckdb_destroy_value(&val);
                destroy_mosdepth_bind(bind);
                duckdb_bind_set_error(info, err);
                return;
            }
            duckdb_free(quantize);
        }
    }
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "no_per_base");
    if (val && !duckdb_is_null_value(val)) bind->no_per_base = duckdb_get_bool(val) ? 1 : 0;
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "fast_mode");
    if (val && !duckdb_is_null_value(val)) bind->fast_mode = duckdb_get_bool(val) ? 1 : 0;
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "fragment_mode");
    if (val && !duckdb_is_null_value(val)) bind->fragment_mode = duckdb_get_bool(val) ? 1 : 0;
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "use_median");
    if (val && !duckdb_is_null_value(val)) bind->use_median = duckdb_get_bool(val) ? 1 : 0;
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "overwrite");
    if (val && !duckdb_is_null_value(val)) bind->overwrite = duckdb_get_bool(val) ? 1 : 0;
    if (val) duckdb_destroy_value(&val);

    if (bind->threads < 0) {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info, "duckhts_mosdepth: threads must be >= 0");
        return;
    }
    if (bind->processing_threads < 0) {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info, "duckhts_mosdepth: processing_threads must be >= 0");
        return;
    }
    if (bind->fasta && bind->fasta[0] == '\0') {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info, "duckhts_mosdepth: fasta must be non-empty when provided");
        return;
    }
    if (bind->read_groups && bind->read_groups[0] == '\0') {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info, "duckhts_mosdepth: read_groups must be non-empty when provided");
        return;
    }
    if (bind->flag < 0 || bind->flag > 65535 || bind->include_flag < 0 || bind->include_flag > 65535) {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info, "duckhts_mosdepth: flag and include_flag must be between 0 and 65535");
        return;
    }
    if (bind->mapq < 0 || bind->mapq > INT32_MAX) {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info, "duckhts_mosdepth: mapq must be >= 0");
        return;
    }
    if (bind->precision < 0 || bind->precision > 18) {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info, "duckhts_mosdepth: precision_digits must be between 0 and 18");
        return;
    }
    if (bind->fast_mode && bind->fragment_mode) {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info,
                              "duckhts_mosdepth: only one of fast_mode and fragment_mode can be TRUE");
        return;
    }
    if (bind->min_frag_len < -1 || bind->max_frag_len < -1) {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info, "duckhts_mosdepth: min_frag_len and max_frag_len must be >= -1");
        return;
    }
    if (bind->max_frag_len >= 0 && bind->min_frag_len >= 0 &&
        bind->max_frag_len < bind->min_frag_len) {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info,
                              "duckhts_mosdepth: max_frag_len must be >= min_frag_len");
        return;
    }
    if (bind->thresholds.count > 0 && (!bind->by || bind->by[0] == '\0')) {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info,
                              "duckhts_mosdepth: thresholds requires by := '<window|bed>'");
        return;
    }

    bind->summary_path = append_suffix(bind->prefix, ".mosdepth.summary.txt");
    bind->global_dist_path = append_suffix(bind->prefix, ".mosdepth.global.dist.txt");
    bind->per_base_path = bind->no_per_base ? NULL : append_suffix(bind->prefix, ".per-base.bed.gz");
    bind->regions_path = (bind->by && bind->by[0] != '\0') ? append_suffix(bind->prefix, ".regions.bed.gz") : NULL;
    bind->region_dist_path = (bind->by && bind->by[0] != '\0') ? append_suffix(bind->prefix, ".mosdepth.region.dist.txt") : NULL;
    bind->quantized_path = (bind->quantize.count > 0) ? append_suffix(bind->prefix, ".quantized.bed.gz") : NULL;
    bind->thresholds_path = (bind->thresholds.count > 0) ? append_suffix(bind->prefix, ".thresholds.bed.gz") : NULL;
    if (!bind->summary_path || !bind->global_dist_path ||
        (!bind->no_per_base && !bind->per_base_path) ||
        ((bind->by && bind->by[0] != '\0') && (!bind->regions_path || !bind->region_dist_path)) ||
        (bind->quantize.count > 0 && !bind->quantized_path) ||
        (bind->thresholds.count > 0 && !bind->thresholds_path)) {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info, "duckhts_mosdepth: out of memory");
        return;
    }

    if (ensure_output_available(bind->summary_path, bind->overwrite, err, sizeof(err)) != 0 ||
        ensure_output_available(bind->global_dist_path, bind->overwrite, err, sizeof(err)) != 0 ||
        ensure_output_available(bind->per_base_path, bind->overwrite, err, sizeof(err)) != 0 ||
        ensure_output_available(bind->regions_path, bind->overwrite, err, sizeof(err)) != 0 ||
        ensure_output_available(bind->region_dist_path, bind->overwrite, err, sizeof(err)) != 0 ||
        ensure_output_available(bind->quantized_path, bind->overwrite, err, sizeof(err)) != 0 ||
        ensure_output_available(bind->thresholds_path, bind->overwrite, err, sizeof(err)) != 0) {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info, err);
        return;
    }

    if (run_duckhts_mosdepth(bind, err, sizeof(err)) != 0) {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info, err);
        return;
    }

    add_result_columns(info);
    bind->emitted = 0;
    duckdb_bind_set_bind_data(info, bind, destroy_mosdepth_bind);
}

static void duckhts_mosdepth_init(duckdb_init_info info) {
    mosdepth_bind_t *bind = (mosdepth_bind_t *)duckdb_init_get_bind_data(info);
    bind->emitted = 0;
    duckdb_init_set_max_threads(info, 1);
}

static void duckhts_mosdepth_scan(duckdb_function_info info, duckdb_data_chunk output) {
    mosdepth_bind_t *bind = (mosdepth_bind_t *)duckdb_function_get_bind_data(info);
    if (!bind || bind->emitted) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    duckdb_vector success_vec = duckdb_data_chunk_get_vector(output, 0);
    duckdb_vector prefix_vec = duckdb_data_chunk_get_vector(output, 1);
    duckdb_vector summary_vec = duckdb_data_chunk_get_vector(output, 2);
    duckdb_vector global_vec = duckdb_data_chunk_get_vector(output, 3);
    duckdb_vector per_base_vec = duckdb_data_chunk_get_vector(output, 4);
    duckdb_vector regions_vec = duckdb_data_chunk_get_vector(output, 5);
    duckdb_vector region_dist_vec = duckdb_data_chunk_get_vector(output, 6);
    duckdb_vector quantized_vec = duckdb_data_chunk_get_vector(output, 7);
    duckdb_vector thresholds_vec = duckdb_data_chunk_get_vector(output, 8);

    bool *success = (bool *)duckdb_vector_get_data(success_vec);
    success[0] = true;
    duckdb_vector_assign_string_element(prefix_vec, 0, bind->prefix);
    duckdb_vector_assign_string_element(summary_vec, 0, bind->summary_path);
    duckdb_vector_assign_string_element(global_vec, 0, bind->global_dist_path);
    if (bind->per_base_path) duckdb_vector_assign_string_element(per_base_vec, 0, bind->per_base_path);
    else set_null(per_base_vec, 0);
    if (bind->regions_path) duckdb_vector_assign_string_element(regions_vec, 0, bind->regions_path);
    else set_null(regions_vec, 0);
    if (bind->region_dist_path) duckdb_vector_assign_string_element(region_dist_vec, 0, bind->region_dist_path);
    else set_null(region_dist_vec, 0);
    if (bind->quantized_path) duckdb_vector_assign_string_element(quantized_vec, 0, bind->quantized_path);
    else set_null(quantized_vec, 0);
    if (bind->thresholds_path) duckdb_vector_assign_string_element(thresholds_vec, 0, bind->thresholds_path);
    else set_null(thresholds_vec, 0);

    bind->emitted = 1;
    duckdb_data_chunk_set_size(output, 1);
}

void register_duckhts_mosdepth_function(duckdb_connection connection) {
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);

    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_table_function_set_name(tf, "duckhts_mosdepth");
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "chrom", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "by", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "fasta", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "read_groups", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "no_per_base", bool_type);
    duckdb_table_function_add_named_parameter(tf, "threads", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "processing_threads", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "flag", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "include_flag", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "fast_mode", bool_type);
    duckdb_table_function_add_named_parameter(tf, "fragment_mode", bool_type);
    duckdb_table_function_add_named_parameter(tf, "use_median", bool_type);
    duckdb_table_function_add_named_parameter(tf, "mapq", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "min_frag_len", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "max_frag_len", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "precision_digits", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "quantize", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "thresholds", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "overwrite", bool_type);
    duckdb_table_function_set_bind(tf, duckhts_mosdepth_bind);
    duckdb_table_function_set_init(tf, duckhts_mosdepth_init);
    duckdb_table_function_set_function(tf, duckhts_mosdepth_scan);
    duckdb_register_table_function(connection, tf);

    duckdb_destroy_table_function(&tf);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_logical_type(&bigint_type);
}
