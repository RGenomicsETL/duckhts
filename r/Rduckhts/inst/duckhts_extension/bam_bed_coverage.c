/**
 * DuckHTS native BAM/CRAM BED regional coverage summary.
 *
 * Initial target: samtools coverage-like per-region summaries with DuckHTS
 * pre/post columns and read-mode strand-specific post summaries.
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

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>

#define DUCKHTS_BEDCOV_DEFAULT_EXCLUDE_FLAGS 1796
#define DUCKHTS_BEDCOV_DEFAULT_MAX_DEPTH 1000000
#define DUCKHTS_BEDCOV_TILE_SIZE 1000000

enum {
    BEDCOV_COL_CHROM = 0,
    BEDCOV_COL_START,
    BEDCOV_COL_END,
    BEDCOV_COL_NAME,
    BEDCOV_COL_NUMREADS_PRE,
    BEDCOV_COL_COVBASES_PRE,
    BEDCOV_COL_COVERAGE_PRE,
    BEDCOV_COL_MEANDEPTH_PRE,
    BEDCOV_COL_MEANBASEQ_PRE,
    BEDCOV_COL_MEANMAPQ_PRE,
    BEDCOV_COL_NUMREADS_POST,
    BEDCOV_COL_COVBASES_POST,
    BEDCOV_COL_COVERAGE_POST,
    BEDCOV_COL_MEANDEPTH_POST,
    BEDCOV_COL_MEANBASEQ_POST,
    BEDCOV_COL_MEANMAPQ_POST,
    BEDCOV_COL_NUMREADS_FWD_POST,
    BEDCOV_COL_NUMREADS_REV_POST,
    BEDCOV_COL_COVBASES_FWD_POST,
    BEDCOV_COL_COVBASES_REV_POST,
    BEDCOV_COL_COVERAGE_FWD_POST,
    BEDCOV_COL_COVERAGE_REV_POST,
    BEDCOV_COL_MEANDEPTH_FWD_POST,
    BEDCOV_COL_MEANDEPTH_REV_POST,
    BEDCOV_COL_COUNT
};

typedef struct {
    char *chrom;
    int64_t start;
    int64_t end;
    char *name;
    int tid;

    uint64_t numreads_pre;
    uint64_t mapq_sum_pre;
    uint64_t numreads_post;
    uint64_t mapq_sum_post;
    uint64_t numreads_fwd_post;
    uint64_t numreads_rev_post;

    uint64_t baseq_sum_pre;
    uint64_t baseq_count_pre;
    uint64_t baseq_sum_post;
    uint64_t baseq_count_post;

    uint64_t covbases_pre;
    uint64_t covbases_post;
    uint64_t covbases_fwd_post;
    uint64_t covbases_rev_post;
    uint64_t summed_cov_pre;
    uint64_t summed_cov_post;
    uint64_t summed_cov_fwd_post;
    uint64_t summed_cov_rev_post;

    uint32_t *depth_pre;
    uint32_t *depth_post;
    uint32_t *depth_fwd_post;
    uint32_t *depth_rev_post;
    int64_t len;
} bedcov_region_t;

typedef struct {
    char *path;
    char *bed_path;
    char *reference;
    char *index_path;
    char *bed_index_path;
    int64_t mapq;
    int64_t min_baseq;
    int64_t min_read_len;
    int64_t require_flags;
    int64_t exclude_flags;
    int64_t min_depth;
    int64_t max_depth;
    int64_t decompression_threads;
    int fragment_mode;
    int strand_outputs;
    int64_t processing_threads;

    bedcov_region_t *regions;
    size_t n_regions;
    size_t emitted;
} bam_bed_cov_bind_t;

typedef struct bedcov_tile_t {
    int64_t start;
    int64_t end;
    int64_t len;
    uint32_t *depth_pre;
    uint32_t *depth_post;
    uint32_t *depth_fwd_post;
    uint32_t *depth_rev_post;
    struct bedcov_tile_t *next;
} bedcov_tile_t;

static inline void set_null(duckdb_vector vec, idx_t row) {
    duckdb_vector_ensure_validity_writable(vec);
    uint64_t *validity = duckdb_vector_get_validity(vec);
    validity[row / 64] &= ~((uint64_t)1 << (row % 64));
}

static char *dup_cstr_bedcov(const char *s) {
    size_t n;
    char *out;
    if (!s) return NULL;
    n = strlen(s);
    out = (char *)duckdb_malloc(n + 1);
    if (!out) return NULL;
    memcpy(out, s, n + 1);
    return out;
}

static int is_meta_bed_line_bedcov(const char *s) {
    if (!s || !*s) return 1;
    if (s[0] == '#') return 1;
    if (strncmp(s, "track", 5) == 0) return 1;
    if (strncmp(s, "browser", 7) == 0) return 1;
    return 0;
}

static int64_t parse_i64_field(const char *s, int len, int *ok) {
    char buf[64];
    char *end = NULL;
    long long v;
    if (ok) *ok = 0;
    if (!s || len <= 0 || (size_t)len >= sizeof(buf)) return 0;
    memcpy(buf, s, (size_t)len);
    buf[len] = '\0';
    v = strtoll(buf, &end, 10);
    if (!end || *end != '\0') return 0;
    if (ok) *ok = 1;
    return (int64_t)v;
}

static const char *field_span(const char *s, int idx, int *len) {
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

static int ensure_region_depth_arrays(bedcov_region_t *region) {
    if (!region || region->len < 0) return -1;
    region->depth_pre = (uint32_t *)calloc((size_t)region->len, sizeof(uint32_t));
    region->depth_post = (uint32_t *)calloc((size_t)region->len, sizeof(uint32_t));
    region->depth_fwd_post = (uint32_t *)calloc((size_t)region->len, sizeof(uint32_t));
    region->depth_rev_post = (uint32_t *)calloc((size_t)region->len, sizeof(uint32_t));
    if (!region->depth_pre || !region->depth_post || !region->depth_fwd_post || !region->depth_rev_post) {
        free(region->depth_pre);
        free(region->depth_post);
        free(region->depth_fwd_post);
        free(region->depth_rev_post);
        region->depth_pre = NULL;
        region->depth_post = NULL;
        region->depth_fwd_post = NULL;
        region->depth_rev_post = NULL;
        return -1;
    }
    return 0;
}

static void destroy_region_depth_arrays(bedcov_region_t *region) {
    if (!region) return;
    free(region->depth_pre);
    free(region->depth_post);
    free(region->depth_fwd_post);
    free(region->depth_rev_post);
    region->depth_pre = NULL;
    region->depth_post = NULL;
    region->depth_fwd_post = NULL;
    region->depth_rev_post = NULL;
}

static int ensure_tile_depth_arrays(bedcov_tile_t *tile, int strand_outputs) {
    if (!tile || tile->len < 0) return -1;
    tile->depth_pre = (uint32_t *)calloc((size_t)tile->len, sizeof(uint32_t));
    tile->depth_post = (uint32_t *)calloc((size_t)tile->len, sizeof(uint32_t));
    if (strand_outputs) {
        tile->depth_fwd_post = (uint32_t *)calloc((size_t)tile->len, sizeof(uint32_t));
        tile->depth_rev_post = (uint32_t *)calloc((size_t)tile->len, sizeof(uint32_t));
    }
    if (!tile->depth_pre || !tile->depth_post || (strand_outputs && (!tile->depth_fwd_post || !tile->depth_rev_post))) {
        free(tile->depth_pre);
        free(tile->depth_post);
        free(tile->depth_fwd_post);
        free(tile->depth_rev_post);
        tile->depth_pre = NULL;
        tile->depth_post = NULL;
        tile->depth_fwd_post = NULL;
        tile->depth_rev_post = NULL;
        return -1;
    }
    return 0;
}

static void destroy_tile_depth_arrays(bedcov_tile_t *tile) {
    if (!tile) return;
    free(tile->depth_pre);
    free(tile->depth_post);
    free(tile->depth_fwd_post);
    free(tile->depth_rev_post);
    tile->depth_pre = NULL;
    tile->depth_post = NULL;
    tile->depth_fwd_post = NULL;
    tile->depth_rev_post = NULL;
}

static void destroy_tile_list(bedcov_tile_t *tile) {
    while (tile) {
        bedcov_tile_t *next = tile->next;
        destroy_tile_depth_arrays(tile);
        free(tile);
        tile = next;
    }
}

static void flush_completed_tiles(bedcov_tile_t **head, bedcov_region_t *region,
                                  const bam_bed_cov_bind_t *bind, int64_t threshold_pos) {
    while (*head && (*head)->end <= threshold_pos) {
        bedcov_tile_t *tile = *head;
        int64_t j;
        for (j = 0; j < tile->len; j++) {
            uint32_t dpre = tile->depth_pre[j];
            uint32_t dpost = tile->depth_post[j];
            uint32_t dfwd = tile->depth_fwd_post ? tile->depth_fwd_post[j] : 0;
            uint32_t drev = tile->depth_rev_post ? tile->depth_rev_post[j] : 0;
            if (dpre >= (uint32_t)bind->min_depth) {
                region->covbases_pre++;
                region->summed_cov_pre += dpre;
            }
            if (dpost >= (uint32_t)bind->min_depth) {
                region->covbases_post++;
                region->summed_cov_post += dpost;
            }
            if (dfwd >= (uint32_t)bind->min_depth) {
                region->covbases_fwd_post++;
                region->summed_cov_fwd_post += dfwd;
            }
            if (drev >= (uint32_t)bind->min_depth) {
                region->covbases_rev_post++;
                region->summed_cov_rev_post += drev;
            }
        }
        *head = tile->next;
        destroy_tile_depth_arrays(tile);
        free(tile);
    }
}

static bedcov_tile_t *get_or_create_tile(bedcov_tile_t **head, int64_t tile_start,
                                         int64_t tile_end, int strand_outputs,
                                         char *err, size_t errlen) {
    bedcov_tile_t *cur = *head;
    bedcov_tile_t *prev = NULL;
    /* Active tiles are flushed in coordinate order as soon as later records can
     * no longer touch them, so this linear scan stays bounded by a tiny live
     * set in the normal streaming path. */
    while (cur) {
        if (cur->start == tile_start) return cur;
        prev = cur;
        cur = cur->next;
    }
    cur = (bedcov_tile_t *)calloc(1, sizeof(*cur));
    if (!cur) {
        snprintf(err, errlen, "duckhts_bam_bed_coverage: out of memory");
        return NULL;
    }
    cur->start = tile_start;
    cur->end = tile_end;
    cur->len = tile_end - tile_start;
    if (ensure_tile_depth_arrays(cur, strand_outputs) != 0) {
        free(cur);
        snprintf(err, errlen, "duckhts_bam_bed_coverage: out of memory");
        return NULL;
    }
    if (prev) prev->next = cur;
    else *head = cur;
    return cur;
}

static void destroy_bam_bed_cov_bind(void *data) {
    bam_bed_cov_bind_t *bind = (bam_bed_cov_bind_t *)data;
    size_t i;
    if (!bind) return;
    if (bind->path) duckdb_free(bind->path);
    if (bind->bed_path) duckdb_free(bind->bed_path);
    if (bind->reference) duckdb_free(bind->reference);
    if (bind->index_path) duckdb_free(bind->index_path);
    if (bind->bed_index_path) duckdb_free(bind->bed_index_path);
    if (bind->regions) {
        for (i = 0; i < bind->n_regions; i++) {
            if (bind->regions[i].chrom) duckdb_free(bind->regions[i].chrom);
            if (bind->regions[i].name) duckdb_free(bind->regions[i].name);
            destroy_region_depth_arrays(&bind->regions[i]);
        }
        free(bind->regions);
    }
    duckdb_free(bind);
}

static int load_bed_regions(const char *path, bedcov_region_t **out_regions, size_t *out_count,
                            char *err, size_t errlen) {
    htsFile *fp = NULL;
    kstring_t line = {0, 0, NULL};
    int linelen;
    size_t cap = 0;
    size_t n = 0;
    bedcov_region_t *regions = NULL;

    if (!path || !out_regions || !out_count) return -1;
    fp = hts_open(path, "r");
    if (!fp) {
        snprintf(err, errlen, "duckhts_bam_bed_coverage: failed to open BED file '%s'", path);
        return -1;
    }

    while ((linelen = hts_getline(fp, '\n', &line)) >= 0) {
        const char *chrom;
        const char *name;
        int chrom_len = 0, start_len = 0, end_len = 0, name_len = 0;
        int ok = 0;
        int64_t start = 0, end = 0;
        bedcov_region_t *r;

        while (linelen > 0 && (line.s[linelen - 1] == '\n' || line.s[linelen - 1] == '\r')) {
            line.s[--linelen] = '\0';
        }
        if (linelen == 0 || is_meta_bed_line_bedcov(line.s)) continue;

        chrom = field_span(line.s, 0, &chrom_len);
        (void)field_span(line.s, 1, &start_len);
        (void)field_span(line.s, 2, &end_len);
        if (!chrom || chrom_len <= 0 || start_len <= 0 || end_len <= 0) {
            snprintf(err, errlen, "duckhts_bam_bed_coverage: BED line has fewer than 3 tab-delimited fields");
            hts_close(fp);
            free(line.s);
            return -1;
        }

        start = parse_i64_field(field_span(line.s, 1, &start_len), start_len, &ok);
        if (!ok) {
            snprintf(err, errlen, "duckhts_bam_bed_coverage: invalid BED start coordinate");
            hts_close(fp);
            free(line.s);
            return -1;
        }
        end = parse_i64_field(field_span(line.s, 2, &end_len), end_len, &ok);
        if (!ok || end < start) {
            snprintf(err, errlen, "duckhts_bam_bed_coverage: invalid BED end coordinate");
            hts_close(fp);
            free(line.s);
            return -1;
        }

        if (n == cap) {
            size_t new_cap = cap ? cap * 2 : 32;
            bedcov_region_t *tmp = (bedcov_region_t *)realloc(regions, new_cap * sizeof(*regions));
            if (!tmp) {
                snprintf(err, errlen, "duckhts_bam_bed_coverage: out of memory");
                hts_close(fp);
                free(line.s);
                return -1;
            }
            regions = tmp;
            memset(regions + cap, 0, (new_cap - cap) * sizeof(*regions));
            cap = new_cap;
        }

        r = &regions[n++];
        memset(r, 0, sizeof(*r));
        r->chrom = (char *)duckdb_malloc((size_t)chrom_len + 1);
        if (!r->chrom) {
            snprintf(err, errlen, "duckhts_bam_bed_coverage: out of memory");
            hts_close(fp);
            free(line.s);
            return -1;
        }
        memcpy(r->chrom, chrom, (size_t)chrom_len);
        r->chrom[chrom_len] = '\0';
        r->start = start;
        r->end = end;
        r->len = end - start;
        r->tid = -1;

        name = field_span(line.s, 3, &name_len);
        if (name && name_len > 0) {
            r->name = (char *)duckdb_malloc((size_t)name_len + 1);
            if (!r->name) {
                snprintf(err, errlen, "duckhts_bam_bed_coverage: out of memory");
                hts_close(fp);
                free(line.s);
                return -1;
            }
            memcpy(r->name, name, (size_t)name_len);
            r->name[name_len] = '\0';
        }

    }

    hts_close(fp);
    free(line.s);
    *out_regions = regions;
    *out_count = n;
    return 0;
}

static int region_set_tids(sam_hdr_t *hdr, bedcov_region_t *regions, size_t n_regions,
                           char *err, size_t errlen) {
    size_t i;
    for (i = 0; i < n_regions; i++) {
        int tid;
        if (regions[i].len <= 0) continue;
        tid = sam_hdr_name2tid(hdr, regions[i].chrom);
        if (tid < 0) {
            snprintf(err, errlen,
                     "duckhts_bam_bed_coverage: chromosome '%s' not found in alignment header",
                     regions[i].chrom);
            return -1;
        }
        regions[i].tid = tid;
        if (regions[i].end > sam_hdr_tid2len(hdr, tid)) {
            regions[i].end = sam_hdr_tid2len(hdr, tid);
            regions[i].len = regions[i].end - regions[i].start;
        }
    }
    return 0;
}

static inline int read_passes_post_filters(const bam1_t *rec, const bam_bed_cov_bind_t *bind) {
    int32_t qlen;
    if (!rec || !bind) return 0;
    if ((rec->core.flag & (uint16_t)bind->exclude_flags) != 0) return 0;
    if ((rec->core.flag & (uint16_t)bind->require_flags) != (uint16_t)bind->require_flags) return 0;
    if (rec->core.qual < bind->mapq) return 0;
    qlen = bam_cigar2qlen(rec->core.n_cigar, bam_get_cigar(rec));
    if (qlen < bind->min_read_len) return 0;
    return 1;
}

static void increment_depth_cap(uint32_t *depth, int64_t idx, int64_t max_depth) {
    if (!depth || idx < 0) return;
    if (max_depth == 0 || depth[idx] < (uint32_t)max_depth) depth[idx]++;
}

static void accumulate_record_summary(const bam1_t *rec, const bam_bed_cov_bind_t *bind,
                                      bedcov_region_t *region) {
    uint32_t *cigar;
    int i;
    int post_ok;
    int is_rev;
    int64_t ref_pos;
    int64_t query_pos;
    uint8_t *qual;
    int has_qual;

    if (!rec || !bind || !region || region->len <= 0) return;

    post_ok = read_passes_post_filters(rec, bind);
    is_rev = (rec->core.flag & BAM_FREVERSE) != 0;
    qual = bam_get_qual(rec);
    has_qual = (rec->core.l_qseq > 0 && qual && qual[0] != 255);

    region->numreads_pre++;
    region->mapq_sum_pre += (uint64_t)rec->core.qual;
    if (post_ok) {
        region->numreads_post++;
        region->mapq_sum_post += (uint64_t)rec->core.qual;
        if (is_rev) region->numreads_rev_post++;
        else region->numreads_fwd_post++;
    }

    cigar = bam_get_cigar(rec);
    ref_pos = rec->core.pos;
    query_pos = 0;
    for (i = 0; i < rec->core.n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int oplen = bam_cigar_oplen(cigar[i]);
        int consumes_ref = bam_cigar_type(op) & 2;
        int consumes_query = bam_cigar_type(op) & 1;

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            int64_t op_start = ref_pos;
            int64_t op_end = ref_pos + oplen;
            int64_t clip_start = op_start < region->start ? region->start : op_start;
            int64_t clip_end = op_end > region->end ? region->end : op_end;
            if (clip_end > clip_start) {
                int64_t query_off = query_pos + (clip_start - op_start);
                int64_t k;
                for (k = clip_start; k < clip_end; k++, query_off++) {
                    if (has_qual && query_off >= 0 && query_off < rec->core.l_qseq) {
                        region->baseq_sum_pre += qual[query_off];
                        region->baseq_count_pre++;
                    }
                    if (post_ok) {
                        int baseq_ok = (!has_qual) || (query_off >= 0 && query_off < rec->core.l_qseq && qual[query_off] >= bind->min_baseq);
                        if (baseq_ok) {
                            if (has_qual && query_off >= 0 && query_off < rec->core.l_qseq) {
                                region->baseq_sum_post += qual[query_off];
                                region->baseq_count_post++;
                            }
                        }
                    }
                }
            }
        }
        if (consumes_ref) ref_pos += oplen;
        if (consumes_query) query_pos += oplen;
    }
}

static void accumulate_record_on_tile(const bam1_t *rec, const bam_bed_cov_bind_t *bind,
                                      int64_t tile_start, int64_t tile_end,
                                      uint32_t *depth_pre, uint32_t *depth_post,
                                      uint32_t *depth_fwd_post, uint32_t *depth_rev_post) {
    uint32_t *cigar;
    int i;
    int post_ok;
    int is_rev;
    int64_t ref_pos;
    int64_t query_pos;
    uint8_t *qual;
    int has_qual;

    if (!rec || !bind || tile_end <= tile_start) return;

    post_ok = read_passes_post_filters(rec, bind);
    is_rev = (rec->core.flag & BAM_FREVERSE) != 0;
    qual = bam_get_qual(rec);
    has_qual = (rec->core.l_qseq > 0 && qual && qual[0] != 255);

    cigar = bam_get_cigar(rec);
    ref_pos = rec->core.pos;
    query_pos = 0;
    for (i = 0; i < rec->core.n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int oplen = bam_cigar_oplen(cigar[i]);
        int consumes_ref = bam_cigar_type(op) & 2;
        int consumes_query = bam_cigar_type(op) & 1;

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            int64_t op_start = ref_pos;
            int64_t op_end = ref_pos + oplen;
            int64_t clip_start = op_start < tile_start ? tile_start : op_start;
            int64_t clip_end = op_end > tile_end ? tile_end : op_end;
            if (clip_end > clip_start) {
                int64_t query_off = query_pos + (clip_start - op_start);
                int64_t k;
                for (k = clip_start; k < clip_end; k++, query_off++) {
                    int64_t local_idx = k - tile_start;
                    increment_depth_cap(depth_pre, local_idx, bind->max_depth);
                    if (post_ok) {
                        int baseq_ok = (!has_qual) || (query_off >= 0 && query_off < rec->core.l_qseq && qual[query_off] >= bind->min_baseq);
                        if (baseq_ok) {
                            increment_depth_cap(depth_post, local_idx, bind->max_depth);
                            if (is_rev) increment_depth_cap(depth_rev_post, local_idx, bind->max_depth);
                            else increment_depth_cap(depth_fwd_post, local_idx, bind->max_depth);
                        }
                    }
                }
            }
        }
        if (consumes_ref) ref_pos += oplen;
        if (consumes_query) query_pos += oplen;
    }
}

static int run_bam_bed_coverage(bam_bed_cov_bind_t *bind, char *err, size_t errlen) {
    samFile *fp = NULL;
    sam_hdr_t *hdr = NULL;
    hts_idx_t *idx = NULL;
    bam1_t *rec = NULL;
    size_t i;

    fp = sam_open(bind->path, "r");
    if (!fp) {
        snprintf(err, errlen, "duckhts_bam_bed_coverage: failed to open alignment file '%s'", bind->path);
        return -1;
    }
    if (bind->decompression_threads > 0 &&
        hts_set_threads(fp, (int)bind->decompression_threads) < 0) {
        snprintf(err, errlen, "duckhts_bam_bed_coverage: failed to configure decompression threads");
        sam_close(fp);
        return -1;
    }
    if (bind->reference && hts_set_fai_filename(fp, bind->reference) != 0) {
        snprintf(err, errlen, "duckhts_bam_bed_coverage: failed to set CRAM reference");
        sam_close(fp);
        return -1;
    }
    hdr = sam_hdr_read(fp);
    if (!hdr) {
        snprintf(err, errlen, "duckhts_bam_bed_coverage: failed to read alignment header");
        sam_close(fp);
        return -1;
    }
    if (region_set_tids(hdr, bind->regions, bind->n_regions, err, errlen) != 0) {
        sam_hdr_destroy(hdr);
        sam_close(fp);
        return -1;
    }

    idx = sam_index_load2(fp, bind->path, bind->index_path);
    if (!idx) {
        snprintf(err, errlen,
                 "duckhts_bam_bed_coverage: indexed BAM/CRAM input is required (failed to load .bai/.csi/.crai)");
        sam_hdr_destroy(hdr);
        sam_close(fp);
        return -1;
    }

    rec = bam_init1();
    if (!rec) {
        snprintf(err, errlen, "duckhts_bam_bed_coverage: failed to allocate BAM record");
        hts_idx_destroy(idx);
        sam_hdr_destroy(hdr);
        sam_close(fp);
        return -1;
    }

    for (i = 0; i < bind->n_regions; i++) {
        hts_itr_t *itr;
        bedcov_region_t *region = &bind->regions[i];
        bedcov_tile_t *tiles = NULL;
        if (region->len <= 0) continue;

        itr = sam_itr_queryi(idx, region->tid, region->start, region->end);
        if (!itr) continue;
        while (sam_itr_next(fp, itr, rec) >= 0) {
            int64_t clip_start;
            int64_t clip_end;
            int64_t tile_start;

            flush_completed_tiles(&tiles, region, bind, rec->core.pos);

            accumulate_record_summary(rec, bind, region);

            clip_start = rec->core.pos < region->start ? region->start : rec->core.pos;
            clip_end = bam_endpos((bam1_t *)rec);
            if (clip_end > region->end) clip_end = region->end;
            if (clip_end <= clip_start) continue;

            for (tile_start = region->start + ((clip_start - region->start) / DUCKHTS_BEDCOV_TILE_SIZE) * DUCKHTS_BEDCOV_TILE_SIZE;
                 tile_start < clip_end;
                 tile_start += DUCKHTS_BEDCOV_TILE_SIZE) {
                int64_t tile_end = tile_start + DUCKHTS_BEDCOV_TILE_SIZE;
                bedcov_tile_t *tile;
                if (tile_end > region->end) tile_end = region->end;
                tile = get_or_create_tile(&tiles, tile_start, tile_end, bind->strand_outputs, err, errlen);
                if (!tile) {
                    hts_itr_destroy(itr);
                    destroy_tile_list(tiles);
                    bam_destroy1(rec);
                    hts_idx_destroy(idx);
                    sam_hdr_destroy(hdr);
                    sam_close(fp);
                    return -1;
                }
                accumulate_record_on_tile(rec, bind, tile_start, tile_end,
                                          tile->depth_pre, tile->depth_post,
                                          tile->depth_fwd_post, tile->depth_rev_post);
            }
        }
        hts_itr_destroy(itr);
        flush_completed_tiles(&tiles, region, bind, INT64_MAX);
        destroy_tile_list(tiles);
    }

    bam_destroy1(rec);
    hts_idx_destroy(idx);
    sam_hdr_destroy(hdr);
    sam_close(fp);
    return 0;
}

static void add_result_columns(duckdb_bind_info info) {
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type double_type = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);

    duckdb_bind_add_result_column(info, "chrom", varchar_type);
    duckdb_bind_add_result_column(info, "start", bigint_type);
    duckdb_bind_add_result_column(info, "end", bigint_type);
    duckdb_bind_add_result_column(info, "name", varchar_type);
    duckdb_bind_add_result_column(info, "numreads_pre", bigint_type);
    duckdb_bind_add_result_column(info, "covbases_pre", bigint_type);
    duckdb_bind_add_result_column(info, "coverage_pre", double_type);
    duckdb_bind_add_result_column(info, "meandepth_pre", double_type);
    duckdb_bind_add_result_column(info, "meanbaseq_pre", double_type);
    duckdb_bind_add_result_column(info, "meanmapq_pre", double_type);
    duckdb_bind_add_result_column(info, "numreads_post", bigint_type);
    duckdb_bind_add_result_column(info, "covbases_post", bigint_type);
    duckdb_bind_add_result_column(info, "coverage_post", double_type);
    duckdb_bind_add_result_column(info, "meandepth_post", double_type);
    duckdb_bind_add_result_column(info, "meanbaseq_post", double_type);
    duckdb_bind_add_result_column(info, "meanmapq_post", double_type);
    duckdb_bind_add_result_column(info, "numreads_fwd_post", bigint_type);
    duckdb_bind_add_result_column(info, "numreads_rev_post", bigint_type);
    duckdb_bind_add_result_column(info, "covbases_fwd_post", bigint_type);
    duckdb_bind_add_result_column(info, "covbases_rev_post", bigint_type);
    duckdb_bind_add_result_column(info, "coverage_fwd_post", double_type);
    duckdb_bind_add_result_column(info, "coverage_rev_post", double_type);
    duckdb_bind_add_result_column(info, "meandepth_fwd_post", double_type);
    duckdb_bind_add_result_column(info, "meandepth_rev_post", double_type);

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&double_type);
}

static void bam_bed_cov_bind(duckdb_bind_info info) {
    bam_bed_cov_bind_t *bind = NULL;
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    duckdb_value bed_val = duckdb_bind_get_parameter(info, 1);
    duckdb_value val;
    char err[512];

    memset(err, 0, sizeof(err));

    bind = (bam_bed_cov_bind_t *)duckdb_malloc(sizeof(*bind));
    if (!bind) {
        duckdb_bind_set_error(info, "duckhts_bam_bed_coverage: out of memory");
        duckdb_destroy_value(&path_val);
        duckdb_destroy_value(&bed_val);
        return;
    }
    memset(bind, 0, sizeof(*bind));
    bind->path = duckdb_get_varchar(path_val);
    bind->bed_path = duckdb_get_varchar(bed_val);
    bind->mapq = 0;
    bind->min_baseq = 0;
    bind->min_read_len = 0;
    bind->require_flags = 0;
    bind->exclude_flags = DUCKHTS_BEDCOV_DEFAULT_EXCLUDE_FLAGS;
    bind->min_depth = 1;
    bind->max_depth = DUCKHTS_BEDCOV_DEFAULT_MAX_DEPTH;
    bind->decompression_threads = 0;
    bind->fragment_mode = 0;
    bind->strand_outputs = 1;
    bind->processing_threads = 0;
    duckdb_destroy_value(&path_val);
    duckdb_destroy_value(&bed_val);

    if (!bind->path || !bind->path[0]) {
        destroy_bam_bed_cov_bind(bind);
        duckdb_bind_set_error(info, "duckhts_bam_bed_coverage requires a non-empty alignment path");
        return;
    }
    if (!bind->bed_path || !bind->bed_path[0]) {
        destroy_bam_bed_cov_bind(bind);
        duckdb_bind_set_error(info, "duckhts_bam_bed_coverage requires a non-empty bed_path");
        return;
    }

    val = duckdb_bind_get_named_parameter(info, "reference");
    if (val && !duckdb_is_null_value(val)) bind->reference = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "index_path");
    if (val && !duckdb_is_null_value(val)) bind->index_path = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "bed_index_path");
    if (val && !duckdb_is_null_value(val)) bind->bed_index_path = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "mapq");
    if (val && !duckdb_is_null_value(val)) bind->mapq = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "min_baseq");
    if (val && !duckdb_is_null_value(val)) bind->min_baseq = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "min_read_len");
    if (val && !duckdb_is_null_value(val)) bind->min_read_len = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "require_flags");
    if (val && !duckdb_is_null_value(val)) bind->require_flags = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "exclude_flags");
    if (val && !duckdb_is_null_value(val)) bind->exclude_flags = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "min_depth");
    if (val && !duckdb_is_null_value(val)) bind->min_depth = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "max_depth");
    if (val && !duckdb_is_null_value(val)) bind->max_depth = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "decompression_threads");
    if (val && !duckdb_is_null_value(val)) bind->decompression_threads = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "fragment_mode");
    if (val && !duckdb_is_null_value(val)) bind->fragment_mode = duckdb_get_bool(val) ? 1 : 0;
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "strand_outputs");
    if (val && !duckdb_is_null_value(val)) bind->strand_outputs = duckdb_get_bool(val) ? 1 : 0;
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "processing_threads");
    if (val && !duckdb_is_null_value(val)) bind->processing_threads = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    if (bind->mapq < 0 || bind->min_baseq < 0 || bind->min_read_len < 0 ||
        bind->require_flags < 0 || bind->exclude_flags < 0 || bind->min_depth < 0 ||
        bind->max_depth < 0 || bind->decompression_threads < 0 || bind->processing_threads < 0) {
        destroy_bam_bed_cov_bind(bind);
        duckdb_bind_set_error(info, "duckhts_bam_bed_coverage: numeric parameters and thread counts must be non-negative");
        return;
    }
    if (bind->min_depth == 0) bind->min_depth = 1;
    if (bind->fragment_mode) {
        destroy_bam_bed_cov_bind(bind);
        duckdb_bind_set_error(info, "duckhts_bam_bed_coverage: fragment_mode is not implemented yet");
        return;
    }
    if (bind->processing_threads != 0) {
        destroy_bam_bed_cov_bind(bind);
        duckdb_bind_set_error(info, "duckhts_bam_bed_coverage: processing_threads is not implemented yet");
        return;
    }

    if (load_bed_regions(bind->bed_path, &bind->regions, &bind->n_regions, err, sizeof(err)) != 0) {
        destroy_bam_bed_cov_bind(bind);
        duckdb_bind_set_error(info, err);
        return;
    }
    if (run_bam_bed_coverage(bind, err, sizeof(err)) != 0) {
        destroy_bam_bed_cov_bind(bind);
        duckdb_bind_set_error(info, err);
        return;
    }

    add_result_columns(info);
    duckdb_bind_set_bind_data(info, bind, destroy_bam_bed_cov_bind);
}

static void bam_bed_cov_init(duckdb_init_info info) {
    bam_bed_cov_bind_t *bind = (bam_bed_cov_bind_t *)duckdb_init_get_bind_data(info);
    bind->emitted = 0;
}

static void emit_double_or_null(duckdb_vector vec, idx_t row, double value, int is_null) {
    double *data = (double *)duckdb_vector_get_data(vec);
    if (is_null) {
        set_null(vec, row);
    } else {
        data[row] = value;
    }
}

static void bam_bed_cov_scan(duckdb_function_info info, duckdb_data_chunk output) {
    bam_bed_cov_bind_t *bind = (bam_bed_cov_bind_t *)duckdb_function_get_bind_data(info);
    idx_t out_idx = 0;
    const idx_t max_rows = duckdb_vector_size();

    while (bind->emitted < bind->n_regions && out_idx < max_rows) {
        bedcov_region_t *r = &bind->regions[bind->emitted++];
        double interval_len = r->len > 0 ? (double)r->len : 0.0;
        duckdb_vector chrom_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_CHROM);
        duckdb_vector start_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_START);
        duckdb_vector end_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_END);
        duckdb_vector name_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_NAME);
        duckdb_vector numreads_pre_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_NUMREADS_PRE);
        duckdb_vector covbases_pre_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_COVBASES_PRE);
        duckdb_vector coverage_pre_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_COVERAGE_PRE);
        duckdb_vector meandepth_pre_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_MEANDEPTH_PRE);
        duckdb_vector meanbaseq_pre_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_MEANBASEQ_PRE);
        duckdb_vector meanmapq_pre_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_MEANMAPQ_PRE);
        duckdb_vector numreads_post_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_NUMREADS_POST);
        duckdb_vector covbases_post_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_COVBASES_POST);
        duckdb_vector coverage_post_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_COVERAGE_POST);
        duckdb_vector meandepth_post_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_MEANDEPTH_POST);
        duckdb_vector meanbaseq_post_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_MEANBASEQ_POST);
        duckdb_vector meanmapq_post_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_MEANMAPQ_POST);
        duckdb_vector numreads_fwd_post_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_NUMREADS_FWD_POST);
        duckdb_vector numreads_rev_post_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_NUMREADS_REV_POST);
        duckdb_vector covbases_fwd_post_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_COVBASES_FWD_POST);
        duckdb_vector covbases_rev_post_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_COVBASES_REV_POST);
        duckdb_vector coverage_fwd_post_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_COVERAGE_FWD_POST);
        duckdb_vector coverage_rev_post_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_COVERAGE_REV_POST);
        duckdb_vector meandepth_fwd_post_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_MEANDEPTH_FWD_POST);
        duckdb_vector meandepth_rev_post_vec = duckdb_data_chunk_get_vector(output, BEDCOV_COL_MEANDEPTH_REV_POST);
        int64_t *start_data = (int64_t *)duckdb_vector_get_data(start_vec);
        int64_t *end_data = (int64_t *)duckdb_vector_get_data(end_vec);
        int64_t *numreads_pre_data = (int64_t *)duckdb_vector_get_data(numreads_pre_vec);
        int64_t *covbases_pre_data = (int64_t *)duckdb_vector_get_data(covbases_pre_vec);
        int64_t *numreads_post_data = (int64_t *)duckdb_vector_get_data(numreads_post_vec);
        int64_t *covbases_post_data = (int64_t *)duckdb_vector_get_data(covbases_post_vec);
        int64_t *numreads_fwd_post_data = (int64_t *)duckdb_vector_get_data(numreads_fwd_post_vec);
        int64_t *numreads_rev_post_data = (int64_t *)duckdb_vector_get_data(numreads_rev_post_vec);
        int64_t *covbases_fwd_post_data = (int64_t *)duckdb_vector_get_data(covbases_fwd_post_vec);
        int64_t *covbases_rev_post_data = (int64_t *)duckdb_vector_get_data(covbases_rev_post_vec);

        duckdb_vector_assign_string_element(chrom_vec, out_idx, r->chrom ? r->chrom : "");
        start_data[out_idx] = r->start;
        end_data[out_idx] = r->end;
        if (r->name) duckdb_vector_assign_string_element(name_vec, out_idx, r->name);
        else set_null(name_vec, out_idx);

        numreads_pre_data[out_idx] = (int64_t)r->numreads_pre;
        covbases_pre_data[out_idx] = (int64_t)r->covbases_pre;
        emit_double_or_null(coverage_pre_vec, out_idx, interval_len > 0 ? (100.0 * (double)r->covbases_pre / interval_len) : 0.0, interval_len <= 0);
        emit_double_or_null(meandepth_pre_vec, out_idx, interval_len > 0 ? ((double)r->summed_cov_pre / interval_len) : 0.0, interval_len <= 0);
        emit_double_or_null(meanbaseq_pre_vec, out_idx,
                            r->baseq_count_pre ? ((double)r->baseq_sum_pre / (double)r->baseq_count_pre) : 0.0,
                            0);
        emit_double_or_null(meanmapq_pre_vec, out_idx,
                            r->numreads_pre ? ((double)r->mapq_sum_pre / (double)r->numreads_pre) : 0.0,
                            0);

        numreads_post_data[out_idx] = (int64_t)r->numreads_post;
        covbases_post_data[out_idx] = (int64_t)r->covbases_post;
        emit_double_or_null(coverage_post_vec, out_idx, interval_len > 0 ? (100.0 * (double)r->covbases_post / interval_len) : 0.0, interval_len <= 0);
        emit_double_or_null(meandepth_post_vec, out_idx, interval_len > 0 ? ((double)r->summed_cov_post / interval_len) : 0.0, interval_len <= 0);
        emit_double_or_null(meanbaseq_post_vec, out_idx,
                            r->baseq_count_post ? ((double)r->baseq_sum_post / (double)r->baseq_count_post) : 0.0,
                            0);
        emit_double_or_null(meanmapq_post_vec, out_idx,
                            r->numreads_post ? ((double)r->mapq_sum_post / (double)r->numreads_post) : 0.0,
                            0);

        numreads_fwd_post_data[out_idx] = (int64_t)r->numreads_fwd_post;
        numreads_rev_post_data[out_idx] = (int64_t)r->numreads_rev_post;
        covbases_fwd_post_data[out_idx] = (int64_t)r->covbases_fwd_post;
        covbases_rev_post_data[out_idx] = (int64_t)r->covbases_rev_post;
        emit_double_or_null(coverage_fwd_post_vec, out_idx, interval_len > 0 ? (100.0 * (double)r->covbases_fwd_post / interval_len) : 0.0, interval_len <= 0);
        emit_double_or_null(coverage_rev_post_vec, out_idx, interval_len > 0 ? (100.0 * (double)r->covbases_rev_post / interval_len) : 0.0, interval_len <= 0);
        emit_double_or_null(meandepth_fwd_post_vec, out_idx, interval_len > 0 ? ((double)r->summed_cov_fwd_post / interval_len) : 0.0, interval_len <= 0);
        emit_double_or_null(meandepth_rev_post_vec, out_idx, interval_len > 0 ? ((double)r->summed_cov_rev_post / interval_len) : 0.0, interval_len <= 0);

        out_idx++;
    }

    duckdb_data_chunk_set_size(output, out_idx);
}

void register_duckhts_bam_bed_coverage_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);

    duckdb_table_function_set_name(tf, "duckhts_bam_bed_coverage");
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_named_parameter(tf, "reference", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "bed_index_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "mapq", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "min_baseq", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "min_read_len", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "require_flags", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "exclude_flags", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "min_depth", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "max_depth", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "decompression_threads", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "fragment_mode", bool_type);
    duckdb_table_function_add_named_parameter(tf, "strand_outputs", bool_type);
    duckdb_table_function_add_named_parameter(tf, "processing_threads", bigint_type);
    duckdb_table_function_set_bind(tf, bam_bed_cov_bind);
    duckdb_table_function_set_init(tf, bam_bed_cov_init);
    duckdb_table_function_set_function(tf, bam_bed_cov_scan);
    duckdb_register_table_function(connection, tf);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&bool_type);
    duckdb_destroy_table_function(&tf);
}
