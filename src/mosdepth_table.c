/**
 * DuckHTS native mosdepth-compatible fast-mode subset.
 *
 * Current scope:
 *   - indexed BAM input
 *   - fast_mode := TRUE only
 *   - summary/global distribution/per-base output
 *   - optional by := <window> or by := <bed>
 *
 * Deferred options error loudly instead of silently diverging from mosdepth.
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

#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/tbx.h>

#define MOSDEPTH_DEFAULT_PRECISION 2
#define MOSDEPTH_MAX_COVERAGE 400000
#define MOSDEPTH_CSI_MIN_SHIFT 14

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
    uint64_t cum_depth;
    int64_t cum_length;
    uint32_t min_depth;
    uint32_t max_depth;
} mosdepth_stat_t;

typedef struct {
    char *prefix;
    char *path;
    char *chrom;
    char *by;
    char *index_path;
    char *summary_path;
    char *global_dist_path;
    char *per_base_path;
    char *regions_path;
    char *region_dist_path;
    int64_t threads;
    int64_t flag;
    int64_t include_flag;
    int64_t mapq;
    int no_per_base;
    int fast_mode;
    int overwrite;
    int emitted;
} mosdepth_bind_t;

static char *dup_string(const char *s) {
    if (!s) return NULL;
    size_t len = strlen(s) + 1;
    char *out = (char *)duckdb_malloc(len);
    if (!out) return NULL;
    memcpy(out, s, len);
    return out;
}

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

static int get_mosdepth_precision(void) {
    const char *env = getenv("MOSDEPTH_PRECISION");
    if (!env || !*env) return MOSDEPTH_DEFAULT_PRECISION;
    char *end = NULL;
    long v = strtol(env, &end, 10);
    if (!end || *end != '\0' || v < 0 || v > 18) return MOSDEPTH_DEFAULT_PRECISION;
    return (int)v;
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

static void dist_reset(mosdepth_dist_t *dist) {
    if (!dist || !dist->data) return;
    memset(dist->data, 0, dist->len * sizeof(int64_t));
}

static void dist_destroy(mosdepth_dist_t *dist) {
    if (!dist) return;
    free(dist->data);
    dist->data = NULL;
    dist->len = 0;
}

static int dist_add_count(mosdepth_dist_t *dist, int depth, int64_t count) {
    if (!dist || depth < 0 || count <= 0) return 0;
    if (depth > MOSDEPTH_MAX_COVERAGE) depth = MOSDEPTH_MAX_COVERAGE - 10;
    if (dist_ensure(dist, (size_t)depth) != 0) return -1;
    dist->data[depth] += count;
    return 0;
}

static int dist_add_span_coverage(mosdepth_dist_t *dist, const int32_t *coverage,
                                  int64_t start, int64_t stop, int64_t len) {
    if (!dist || !coverage || start >= stop) return 0;
    if (start < 0) start = 0;
    if (stop > len) stop = len;
    for (int64_t i = start; i < stop; i++) {
        if (dist_add_count(dist, depth_bucket_value((uint32_t)coverage[i]), 1) != 0) {
            return -1;
        }
    }
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

static int write_zero_per_base(BGZF *fp, kstring_t *line, const char *chrom, int64_t len) {
    line->l = 0;
    if (ksprintf(line, "%s\t0\t%" PRId64 "\t0", chrom, len) < 0) return -1;
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
        if (ksprintf(line, "%s\t%" PRId64 "\t%" PRId64 "\t%d",
                     chrom, last_start, i, (int)last_depth) < 0) {
            return -1;
        }
        if (bgzf_write_line(fp, line) != 0) return -1;
        last_start = i;
        last_depth = depth;
    }
    line->l = 0;
    if (ksprintf(line, "%s\t%" PRId64 "\t%" PRId64 "\t%d",
                 chrom, last_start, len, (int)last_depth) < 0) {
        return -1;
    }
    return bgzf_write_line(fp, line);
}

static int write_window_regions(BGZF *fp, kstring_t *line, const char *chrom,
                                const int32_t *coverage, int64_t len, int64_t window,
                                mosdepth_stat_t *region_stat,
                                mosdepth_dist_t *region_dist, int precision,
                                int has_reads) {
    for (int64_t start = 0; start < len; start += window) {
        int64_t stop = start + window;
        if (stop > len) stop = len;
        double mean = 0.0;
        if (has_reads) {
            uint64_t sum = 0;
            for (int64_t i = start; i < stop; i++) {
                uint32_t depth = (uint32_t)coverage[i];
                sum += depth;
                stat_add_depth(region_stat, depth);
            }
            mean = (double)sum / (double)(stop - start);
            if (dist_add_count(region_dist, (int)(mean + 0.5), 1) != 0) return -1;
        }
        line->l = 0;
        if (ksprintf(line, "%s\t%" PRId64 "\t%" PRId64 "\t%.*f",
                     chrom, start, stop, precision, mean) < 0) {
            return -1;
        }
        if (bgzf_write_line(fp, line) != 0) return -1;
    }
    return 0;
}

static int write_bed_regions(BGZF *fp, kstring_t *line, const char *chrom,
                             const int32_t *coverage, int64_t len,
                             const mosdepth_region_list_t *regions,
                             mosdepth_stat_t *region_stat,
                             mosdepth_dist_t *region_dist, int precision,
                             int has_reads) {
    if (!regions) return 0;
    for (size_t i = 0; i < regions->count; i++) {
        const mosdepth_region_t *r = &regions->data[i];
        int64_t in_start = r->start < 0 ? 0 : r->start;
        int64_t in_stop = r->stop > len ? len : r->stop;
        double mean = 0.0;
        if (has_reads && in_start < in_stop) {
            uint64_t sum = 0;
            for (int64_t pos = in_start; pos < in_stop; pos++) {
                uint32_t depth = (uint32_t)coverage[pos];
                sum += depth;
                stat_add_depth(region_stat, depth);
                if (dist_add_count(region_dist, depth_bucket_value(depth), 1) != 0) return -1;
            }
            mean = (double)sum / (double)(r->stop - r->start);
        }
        line->l = 0;
        if (r->name && *r->name) {
            if (ksprintf(line, "%s\t%" PRId64 "\t%" PRId64 "\t%s\t%.*f",
                         chrom, r->start, r->stop, r->name, precision, mean) < 0) {
                return -1;
            }
        } else {
            if (ksprintf(line, "%s\t%" PRId64 "\t%" PRId64 "\t%.*f",
                         chrom, r->start, r->stop, precision, mean) < 0) {
                return -1;
            }
        }
        if (bgzf_write_line(fp, line) != 0) return -1;
    }
    return 0;
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
    kstring_t line = {0, 0, NULL};
    mosdepth_region_list_t *region_lists = NULL;
    mosdepth_dist_t total_global_dist = {0};
    mosdepth_dist_t total_region_dist = {0};
    mosdepth_stat_t total_stat;
    mosdepth_stat_t total_region_stat;
    int summary_header_written = 0;
    int precision = get_mosdepth_precision();
    int rc = -1;
    int n_targets;
    int by_is_window = 0;
    int64_t window = 0;

    stat_clear(&total_stat);
    stat_clear(&total_region_stat);

    fp = sam_open(bind->path, "r");
    if (!fp) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to open alignment file '%s'", bind->path);
        goto cleanup;
    }
    if (bind->threads > 0 && hts_set_threads(fp, (int)bind->threads) < 0) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to configure decompression threads");
        goto cleanup;
    }
    if (fp->format.format == cram) {
        snprintf(err, errlen,
                 "duckhts_mosdepth: CRAM input is not supported yet in the native rewrite; use BAM for now");
        goto cleanup;
    }

    hdr = sam_hdr_read(fp);
    if (!hdr) {
        snprintf(err, errlen, "duckhts_mosdepth: failed to read BAM header");
        goto cleanup;
    }
    idx = sam_index_load3(fp, bind->path, bind->index_path, HTS_IDX_SILENT_FAIL | HTS_IDX_SAVE_REMOTE);
    if (!idx) {
        snprintf(err, errlen,
                 "duckhts_mosdepth: indexed BAM input is required (failed to load .bai/.csi)");
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

    for (int tid = 0; tid < n_targets; tid++) {
        const char *chrom = sam_hdr_tid2name(hdr, tid);
        int64_t chrom_len = (int64_t)sam_hdr_tid2len(hdr, tid);
        hts_itr_t *itr = NULL;
        int32_t *coverage = NULL;
        int has_reads = 0;
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

        coverage = (int32_t *)calloc((size_t)chrom_len + 1, sizeof(int32_t));
        if (!coverage) {
            snprintf(err, errlen, "duckhts_mosdepth: failed to allocate coverage array for '%s'", chrom);
            goto contig_cleanup;
        }

        itr = sam_itr_queryi(idx, tid, 0, chrom_len);
        if (itr) {
            while (sam_itr_next(fp, itr, rec) >= 0) {
                hts_pos_t start = rec->core.pos;
                hts_pos_t stop = bam_endpos(rec);

                if ((int64_t)rec->core.qual < bind->mapq) continue;
                if (((uint16_t)rec->core.flag & (uint16_t)bind->flag) != 0) continue;
                if (bind->include_flag != 0 &&
                    (((uint16_t)rec->core.flag & (uint16_t)bind->include_flag) == 0)) {
                    continue;
                }
                if (start < 0) start = 0;
                if (stop < start) stop = start;
                if (start >= chrom_len) continue;
                if (stop > chrom_len) stop = chrom_len;
                coverage[start] += 1;
                coverage[stop] -= 1;
                has_reads = 1;
            }
        }

        if (!has_reads) {
            if (bgzf_per_base && write_zero_per_base(bgzf_per_base, &line, chrom, chrom_len) != 0) {
                snprintf(err, errlen, "duckhts_mosdepth: failed to write per-base output");
                goto contig_cleanup;
            }
            if (bgzf_regions) {
                if (by_is_window) {
                    if (write_window_regions(bgzf_regions, &line, chrom, NULL, chrom_len, window,
                                             &chrom_region_stat, &chrom_region_dist,
                                             precision, 0) != 0) {
                        snprintf(err, errlen, "duckhts_mosdepth: failed to write region output");
                        goto contig_cleanup;
                    }
                } else if (region_lists) {
                    if (write_bed_regions(bgzf_regions, &line, chrom, NULL, chrom_len, &region_lists[tid],
                                          &chrom_region_stat, &chrom_region_dist,
                                          precision, 0) != 0) {
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

        if (bgzf_regions) {
            if (by_is_window) {
                if (write_window_regions(bgzf_regions, &line, chrom, coverage, chrom_len, window,
                                         &chrom_region_stat, &chrom_region_dist,
                                         precision, 1) != 0) {
                    snprintf(err, errlen, "duckhts_mosdepth: failed to write region output");
                    goto contig_cleanup;
                }
            } else if (region_lists) {
                if (write_bed_regions(bgzf_regions, &line, chrom, coverage, chrom_len, &region_lists[tid],
                                      &chrom_region_stat, &chrom_region_dist,
                                      precision, 1) != 0) {
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
        free(coverage);
        dist_destroy(&chrom_global_dist);
        dist_destroy(&chrom_region_dist);
        if (rc != 0) goto cleanup;
    }

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
    if (region_lists) {
        for (int i = 0; i < n_targets; i++) region_list_destroy(&region_lists[i]);
        free(region_lists);
    }
    dist_destroy(&total_global_dist);
    dist_destroy(&total_region_dist);

    if (rc == 0) {
        if (bind->per_base_path && build_bed_csi(bind->per_base_path, err, errlen) != 0) rc = -1;
        if (rc == 0 && bind->regions_path && build_bed_csi(bind->regions_path, err, errlen) != 0) rc = -1;
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
    if (bind->index_path) duckdb_free(bind->index_path);
    if (bind->summary_path) duckdb_free(bind->summary_path);
    if (bind->global_dist_path) duckdb_free(bind->global_dist_path);
    if (bind->per_base_path) duckdb_free(bind->per_base_path);
    if (bind->regions_path) duckdb_free(bind->regions_path);
    if (bind->region_dist_path) duckdb_free(bind->region_dist_path);
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
    bind->threads = 0;
    bind->flag = 1796;
    bind->include_flag = 0;
    bind->mapq = 0;
    bind->fast_mode = 1;
    bind->no_per_base = 0;
    bind->overwrite = 0;

    val = duckdb_bind_get_named_parameter(info, "chrom");
    if (val && !duckdb_is_null_value(val)) bind->chrom = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "by");
    if (val && !duckdb_is_null_value(val)) bind->by = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "index_path");
    if (val && !duckdb_is_null_value(val)) bind->index_path = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "threads");
    if (val && !duckdb_is_null_value(val)) bind->threads = duckdb_get_int64(val);
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

    val = duckdb_bind_get_named_parameter(info, "no_per_base");
    if (val && !duckdb_is_null_value(val)) bind->no_per_base = duckdb_get_bool(val) ? 1 : 0;
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "fast_mode");
    if (val && !duckdb_is_null_value(val)) bind->fast_mode = duckdb_get_bool(val) ? 1 : 0;
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "overwrite");
    if (val && !duckdb_is_null_value(val)) bind->overwrite = duckdb_get_bool(val) ? 1 : 0;
    if (val) duckdb_destroy_value(&val);

    if (bind->threads < 0) {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info, "duckhts_mosdepth: threads must be >= 0");
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
    if (!bind->fast_mode) {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info,
                              "duckhts_mosdepth currently implements only fast_mode := TRUE");
        return;
    }

    bind->summary_path = append_suffix(bind->prefix, ".mosdepth.summary.txt");
    bind->global_dist_path = append_suffix(bind->prefix, ".mosdepth.global.dist.txt");
    bind->per_base_path = bind->no_per_base ? NULL : append_suffix(bind->prefix, ".per-base.bed.gz");
    bind->regions_path = (bind->by && bind->by[0] != '\0') ? append_suffix(bind->prefix, ".regions.bed.gz") : NULL;
    bind->region_dist_path = (bind->by && bind->by[0] != '\0') ? append_suffix(bind->prefix, ".mosdepth.region.dist.txt") : NULL;
    if (!bind->summary_path || !bind->global_dist_path ||
        (!bind->no_per_base && !bind->per_base_path) ||
        ((bind->by && bind->by[0] != '\0') && (!bind->regions_path || !bind->region_dist_path))) {
        destroy_mosdepth_bind(bind);
        duckdb_bind_set_error(info, "duckhts_mosdepth: out of memory");
        return;
    }

    if (ensure_output_available(bind->summary_path, bind->overwrite, err, sizeof(err)) != 0 ||
        ensure_output_available(bind->global_dist_path, bind->overwrite, err, sizeof(err)) != 0 ||
        ensure_output_available(bind->per_base_path, bind->overwrite, err, sizeof(err)) != 0 ||
        ensure_output_available(bind->regions_path, bind->overwrite, err, sizeof(err)) != 0 ||
        ensure_output_available(bind->region_dist_path, bind->overwrite, err, sizeof(err)) != 0) {
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
    duckdb_table_function_add_named_parameter(tf, "no_per_base", bool_type);
    duckdb_table_function_add_named_parameter(tf, "threads", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "flag", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "include_flag", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "fast_mode", bool_type);
    duckdb_table_function_add_named_parameter(tf, "mapq", bigint_type);
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
