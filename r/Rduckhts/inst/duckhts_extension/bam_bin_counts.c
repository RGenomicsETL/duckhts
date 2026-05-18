/**
 * DuckHTS native BAM/CRAM fixed-bin counting
 *
 * Minimal generic interface for fixed-width read-start binning with optional
 * WisecondorX-style streaming duplicate suppression and optional per-bin
 * GC/MAPQ sufficient statistics.
 */

#include "duckdb_extension.h"
DUCKDB_EXTENSION_EXTERN

#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/sam.h>

#include "include/hts_io_tuning.h"

#define DUCKHTS_DEFAULT_BINCOUNT_DECOMPRESSION_THREADS 2
#define DUCKHTS_MAX_BINCOUNT_THREADS 16

typedef enum {
    BAM_BIN_RMDUP_NONE = 0,
    BAM_BIN_RMDUP_FLAG,
    BAM_BIN_RMDUP_STREAMING
} bam_bin_rmdup_mode_t;

typedef struct {
    int has_pos;
    int64_t pos0;
    int has_pnext;
    int64_t pnext0;
} bam_bin_stream_state_t;

typedef struct {
    int tid;
    char *chrom_name;
    int64_t ref_length;
    uint64_t weight;

    int64_t max_bin;
    int64_t cap_bins;

    int64_t *count_total;
    int64_t *count_fwd;
    int64_t *count_rev;
    int64_t *count_pre;

    uint64_t *gc_bases_pre;
    uint64_t *bases_pre;
    uint64_t *gc_bases_post;
    uint64_t *bases_post;
    uint64_t *mapq_sum_pre;
    uint64_t *mapq_sum_post;

    int has_first_record;
    int64_t first_bin;
    int first_is_paired;
    int first_is_reverse;
    int first_local_counted;
    int64_t first_pos0;
    int first_has_pnext;
    int64_t first_pnext0;
    uint64_t first_gc_bases;
    uint64_t first_bases;
    uint8_t first_mapq;

    bam_bin_stream_state_t final_state;
} bam_bin_contig_t;

typedef struct {
    const char *chrom;
    int64_t start;
    int64_t end;
    int64_t bin_id;
    int is_unmapped;
    int64_t count_total;
    int64_t count_fwd;
    int64_t count_rev;
    int64_t count_pre;
    uint64_t gc_bases_pre;
    uint64_t bases_pre;
    uint64_t gc_bases_post;
    uint64_t bases_post;
    uint64_t ref_gc_bases;
    uint64_t ref_bases;
    uint64_t mapq_sum_pre;
    uint64_t mapq_sum_post;
} bam_bin_row_t;

typedef struct {
    char *path;
    char *chrom;
    char *reference;
    char *index_path;

    int64_t bin_width;
    int64_t mapq;
    int64_t require_flags;
    int64_t exclude_flags;
    int64_t decompression_threads;

    bam_bin_rmdup_mode_t rmdup;
    int include_unmapped;
    int want_gc;
    int want_mq;

    int64_t unmapped_count_total;
    int64_t unmapped_count_fwd;
    int64_t unmapped_count_rev;
    int64_t unmapped_count_pre;
    uint64_t unmapped_gc_bases_pre;
    uint64_t unmapped_bases_pre;
    uint64_t unmapped_gc_bases_post;
    uint64_t unmapped_bases_post;
    uint64_t unmapped_mapq_sum_pre;
    uint64_t unmapped_mapq_sum_post;

    bam_bin_contig_t *contigs;
    int n_contigs;

    int *claim_order;
    int n_claims;

    bam_bin_row_t *rows;
    size_t row_count;
    size_t emitted;
} bam_bin_bind_t;

typedef struct {
    bam_bin_bind_t *bind;
    bam_bin_contig_t *contigs;
    int *claim_order;
    int n_claims;
    volatile int next_claim;
    volatile int failed;
    char error[512];
} bam_bin_worker_shared_t;

typedef struct {
    bam_bin_worker_shared_t *shared;
} bam_bin_worker_arg_t;

typedef struct {
    int tid;
    uint64_t weight;
    int selected_idx;
} bam_bin_claim_slot_t;

static inline void set_null(duckdb_vector vec, idx_t row) {
    duckdb_vector_ensure_validity_writable(vec);
    uint64_t *validity = duckdb_vector_get_validity(vec);
    validity[row / 64] &= ~((uint64_t)1 << (row % 64));
}

static char *dup_cstr(const char *s) {
    size_t n;
    char *out;
    if (!s) return NULL;
    n = strlen(s);
    out = (char *)duckdb_malloc(n + 1);
    if (!out) return NULL;
    memcpy(out, s, n + 1);
    return out;
}

static void free_contig_arrays(bam_bin_contig_t *contig) {
    if (!contig) return;
    free(contig->count_total);
    free(contig->count_fwd);
    free(contig->count_rev);
    free(contig->count_pre);
    free(contig->gc_bases_pre);
    free(contig->bases_pre);
    free(contig->gc_bases_post);
    free(contig->bases_post);
    free(contig->mapq_sum_pre);
    free(contig->mapq_sum_post);
    contig->count_total = NULL;
    contig->count_fwd = NULL;
    contig->count_rev = NULL;
    contig->count_pre = NULL;
    contig->gc_bases_pre = NULL;
    contig->bases_pre = NULL;
    contig->gc_bases_post = NULL;
    contig->bases_post = NULL;
    contig->mapq_sum_pre = NULL;
    contig->mapq_sum_post = NULL;
    contig->max_bin = -1;
    contig->cap_bins = 0;
}

static void destroy_bam_bin_bind(void *data) {
    bam_bin_bind_t *bind = (bam_bin_bind_t *)data;
    if (!bind) return;
    if (bind->path) duckdb_free(bind->path);
    if (bind->chrom) duckdb_free(bind->chrom);
    if (bind->reference) duckdb_free(bind->reference);
    if (bind->index_path) duckdb_free(bind->index_path);
    if (bind->claim_order) free(bind->claim_order);
    if (bind->rows) free(bind->rows);
    if (bind->contigs) {
        for (int i = 0; i < bind->n_contigs; i++) {
            if (bind->contigs[i].chrom_name) duckdb_free(bind->contigs[i].chrom_name);
            free_contig_arrays(&bind->contigs[i]);
        }
        free(bind->contigs);
    }
    duckdb_free(bind);
}

static int parse_rmdup_mode(const char *s, bam_bin_rmdup_mode_t *mode) {
    if (!s || !mode) return -1;
    if (strcmp(s, "none") == 0) {
        *mode = BAM_BIN_RMDUP_NONE;
        return 0;
    }
    if (strcmp(s, "flag") == 0) {
        *mode = BAM_BIN_RMDUP_FLAG;
        return 0;
    }
    if (strcmp(s, "streaming") == 0) {
        *mode = BAM_BIN_RMDUP_STREAMING;
        return 0;
    }
    return -1;
}

static int parse_stats_string(const char *stats, int *want_gc, int *want_mq) {
    char *copy = NULL;
    char *tok = NULL;
    char *saveptr = NULL;
    if (want_gc) *want_gc = 0;
    if (want_mq) *want_mq = 0;
    if (!stats || !*stats) return 0;

    copy = strdup(stats);
    if (!copy) return -1;
    for (char *p = copy; *p; p++) {
        *p = (char)tolower((unsigned char)*p);
    }

    for (tok = strtok_r(copy, ",", &saveptr); tok; tok = strtok_r(NULL, ",", &saveptr)) {
        while (*tok && isspace((unsigned char)*tok)) tok++;
        if (strcmp(tok, "gc") == 0) {
            if (want_gc) *want_gc = 1;
            continue;
        }
        if (strcmp(tok, "mq") == 0) {
            if (want_mq) *want_mq = 1;
            continue;
        }
        free(copy);
        return -1;
    }

    free(copy);
    return 0;
}

static int claim_slot_cmp(const void *a, const void *b) {
    const bam_bin_claim_slot_t *lhs = (const bam_bin_claim_slot_t *)a;
    const bam_bin_claim_slot_t *rhs = (const bam_bin_claim_slot_t *)b;
    if (lhs->weight > rhs->weight) return -1;
    if (lhs->weight < rhs->weight) return 1;
    if (lhs->tid < rhs->tid) return -1;
    if (lhs->tid > rhs->tid) return 1;
    return 0;
}

static int build_claim_order(bam_bin_bind_t *bind) {
    bam_bin_claim_slot_t *slots = NULL;
    if (!bind || bind->n_contigs <= 0) return 0;

    slots = (bam_bin_claim_slot_t *)malloc(sizeof(bam_bin_claim_slot_t) * (size_t)bind->n_contigs);
    if (!slots) return -1;

    for (int i = 0; i < bind->n_contigs; i++) {
        slots[i].tid = bind->contigs[i].tid;
        slots[i].weight = bind->contigs[i].weight;
        slots[i].selected_idx = i;
    }

    qsort(slots, (size_t)bind->n_contigs, sizeof(bam_bin_claim_slot_t), claim_slot_cmp);

    bind->claim_order = (int *)malloc(sizeof(int) * (size_t)bind->n_contigs);
    if (!bind->claim_order) {
        free(slots);
        return -1;
    }
    for (int i = 0; i < bind->n_contigs; i++) {
        bind->claim_order[i] = slots[i].selected_idx;
    }
    bind->n_claims = bind->n_contigs;
    free(slots);
    return 0;
}

static int64_t *grow_int64_array(int64_t *ptr, size_t new_cap) {
    return (int64_t *)realloc(ptr, sizeof(int64_t) * new_cap);
}

static uint64_t *grow_uint64_array(uint64_t *ptr, size_t new_cap) {
    return (uint64_t *)realloc(ptr, sizeof(uint64_t) * new_cap);
}

static int ensure_contig_capacity(bam_bin_bind_t *bind, bam_bin_contig_t *contig, int64_t need_bins) {
    int64_t old_cap;
    int64_t new_cap;
    size_t zero_count;
    if (!bind || !contig || need_bins <= contig->cap_bins) return 0;

    old_cap = contig->cap_bins;
    new_cap = old_cap > 0 ? old_cap : 16;
    while (new_cap < need_bins) {
        if (new_cap > (INT64_MAX / 2)) return -1;
        new_cap *= 2;
    }

    contig->count_total = grow_int64_array(contig->count_total, (size_t)new_cap);
    contig->count_fwd = grow_int64_array(contig->count_fwd, (size_t)new_cap);
    contig->count_rev = grow_int64_array(contig->count_rev, (size_t)new_cap);
    if (!contig->count_total || !contig->count_fwd || !contig->count_rev) return -1;

    if (bind->want_gc || bind->want_mq) {
        contig->count_pre = grow_int64_array(contig->count_pre, (size_t)new_cap);
        if (!contig->count_pre) return -1;
    }
    if (bind->want_gc) {
        contig->gc_bases_pre = grow_uint64_array(contig->gc_bases_pre, (size_t)new_cap);
        contig->bases_pre = grow_uint64_array(contig->bases_pre, (size_t)new_cap);
        contig->gc_bases_post = grow_uint64_array(contig->gc_bases_post, (size_t)new_cap);
        contig->bases_post = grow_uint64_array(contig->bases_post, (size_t)new_cap);
        if (!contig->gc_bases_pre || !contig->bases_pre ||
            !contig->gc_bases_post || !contig->bases_post) {
            return -1;
        }
    }
    if (bind->want_mq) {
        contig->mapq_sum_pre = grow_uint64_array(contig->mapq_sum_pre, (size_t)new_cap);
        contig->mapq_sum_post = grow_uint64_array(contig->mapq_sum_post, (size_t)new_cap);
        if (!contig->mapq_sum_pre || !contig->mapq_sum_post) return -1;
    }

    zero_count = (size_t)(new_cap - old_cap);
    memset(contig->count_total + old_cap, 0, sizeof(int64_t) * zero_count);
    memset(contig->count_fwd + old_cap, 0, sizeof(int64_t) * zero_count);
    memset(contig->count_rev + old_cap, 0, sizeof(int64_t) * zero_count);
    if (contig->count_pre) memset(contig->count_pre + old_cap, 0, sizeof(int64_t) * zero_count);
    if (contig->gc_bases_pre) memset(contig->gc_bases_pre + old_cap, 0, sizeof(uint64_t) * zero_count);
    if (contig->bases_pre) memset(contig->bases_pre + old_cap, 0, sizeof(uint64_t) * zero_count);
    if (contig->gc_bases_post) memset(contig->gc_bases_post + old_cap, 0, sizeof(uint64_t) * zero_count);
    if (contig->bases_post) memset(contig->bases_post + old_cap, 0, sizeof(uint64_t) * zero_count);
    if (contig->mapq_sum_pre) memset(contig->mapq_sum_pre + old_cap, 0, sizeof(uint64_t) * zero_count);
    if (contig->mapq_sum_post) memset(contig->mapq_sum_post + old_cap, 0, sizeof(uint64_t) * zero_count);

    contig->cap_bins = new_cap;
    return 0;
}

static inline uint64_t count_read_gc_bases(const bam1_t *rec) {
    uint64_t gc = 0;
    const uint8_t *seq = bam_get_seq(rec);
    int len = rec->core.l_qseq;
    for (int i = 0; i < len; i++) {
        int code = bam_seqi(seq, i);
        if (code == 2 || code == 4) gc++;
    }
    return gc;
}

static int record_is_streaming_dup(const bam_bin_stream_state_t *state,
                                   int64_t pos0,
                                   int is_paired,
                                   int has_pnext,
                                   int64_t pnext0) {
    if (!state || !state->has_pos) return 0;
    if (pos0 != state->pos0) return 0;
    if (!is_paired) return 1;
    if (!has_pnext || !state->has_pnext) return 0;
    return pnext0 == state->pnext0;
}

static void update_streaming_state(bam_bin_stream_state_t *state,
                                   int64_t pos0,
                                   int is_paired,
                                   int has_pnext,
                                   int64_t pnext0) {
    if (!state) return;
    state->has_pos = 1;
    state->pos0 = pos0;
    if (is_paired && has_pnext) {
        state->has_pnext = 1;
        state->pnext0 = pnext0;
    }
}

static void subtract_first_post_metrics(bam_bin_bind_t *bind, bam_bin_contig_t *contig) {
    int64_t bin_id;
    if (!bind || !contig || !contig->has_first_record || !contig->first_local_counted) return;
    bin_id = contig->first_bin;
    if (bin_id < 0 || bin_id > contig->max_bin) return;

    contig->count_total[bin_id]--;
    if (contig->first_is_reverse) contig->count_rev[bin_id]--;
    else contig->count_fwd[bin_id]--;

    if (bind->want_gc) {
        contig->gc_bases_post[bin_id] -= contig->first_gc_bases;
        contig->bases_post[bin_id] -= contig->first_bases;
    }
    if (bind->want_mq) {
        contig->mapq_sum_post[bin_id] -= contig->first_mapq;
    }
}

static int process_record_into_contig(bam_bin_bind_t *bind,
                                      bam_bin_contig_t *contig,
                                      bam1_t *rec,
                                      bam_bin_stream_state_t *state) {
    int64_t pos0;
    int64_t pnext0;
    int has_pnext;
    int64_t bin_id;
    int flag;
    int is_paired;
    int is_reverse;
    int pass_post;
    uint64_t gc_bases = 0;
    uint64_t read_bases = 0;

    if (!bind || !contig || !rec) return -1;
    if (rec->core.tid < 0 || rec->core.pos < 0) return 0;

    flag = rec->core.flag;
    if (bind->require_flags > 0 &&
        (((int64_t)flag & bind->require_flags) != bind->require_flags)) {
        return 0;
    }
    if (bind->exclude_flags > 0 &&
        (((int64_t)flag & bind->exclude_flags) != 0)) {
        return 0;
    }

    is_paired = (flag & BAM_FPAIRED) != 0;
    if (bind->rmdup == BAM_BIN_RMDUP_STREAMING && is_paired && !(flag & BAM_FPROPER_PAIR)) {
        return 0;
    }

    pos0 = (int64_t)rec->core.pos;
    pnext0 = (int64_t)rec->core.mpos;
    has_pnext = pnext0 >= 0;
    bin_id = pos0 / bind->bin_width;
    if (bin_id < 0) return 0;
    if (ensure_contig_capacity(bind, contig, bin_id + 1) != 0) return -1;
    if (bin_id > contig->max_bin) contig->max_bin = bin_id;

    is_reverse = (flag & BAM_FREVERSE) != 0;

    if ((bind->want_gc || bind->want_mq) && contig->count_pre) {
        contig->count_pre[bin_id]++;
    }
    if (bind->want_gc) {
        read_bases = (uint64_t)rec->core.l_qseq;
        gc_bases = count_read_gc_bases(rec);
        contig->gc_bases_pre[bin_id] += gc_bases;
        contig->bases_pre[bin_id] += read_bases;
    }
    if (bind->want_mq) {
        contig->mapq_sum_pre[bin_id] += (uint64_t)rec->core.qual;
    }

    pass_post = 1;
    if (bind->rmdup == BAM_BIN_RMDUP_STREAMING) {
        int is_dup = record_is_streaming_dup(state, pos0, is_paired, has_pnext, pnext0);
        if (!contig->has_first_record) {
            contig->has_first_record = 1;
            contig->first_bin = bin_id;
            contig->first_is_paired = is_paired;
            contig->first_is_reverse = is_reverse;
            contig->first_pos0 = pos0;
            contig->first_has_pnext = has_pnext;
            contig->first_pnext0 = pnext0;
            contig->first_gc_bases = gc_bases;
            contig->first_bases = read_bases;
            contig->first_mapq = rec->core.qual;
        }
        update_streaming_state(state, pos0, is_paired, has_pnext, pnext0);
        if (is_dup) pass_post = 0;
    } else if (bind->rmdup == BAM_BIN_RMDUP_FLAG && (flag & BAM_FDUP)) {
        pass_post = 0;
    }

    if ((int64_t)rec->core.qual < bind->mapq) {
        pass_post = 0;
    }

    if (bind->rmdup == BAM_BIN_RMDUP_STREAMING && contig->has_first_record && contig->first_bin == bin_id &&
        contig->first_pos0 == pos0 && !contig->first_local_counted) {
        contig->first_local_counted = pass_post;
    }

    if (pass_post) {
        contig->count_total[bin_id]++;
        if (is_reverse) contig->count_rev[bin_id]++;
        else contig->count_fwd[bin_id]++;
        if (bind->want_gc) {
            contig->gc_bases_post[bin_id] += gc_bases;
            contig->bases_post[bin_id] += read_bases;
        }
        if (bind->want_mq) {
            contig->mapq_sum_post[bin_id] += (uint64_t)rec->core.qual;
        }
    }

    if (bind->rmdup == BAM_BIN_RMDUP_STREAMING) {
        contig->final_state = *state;
    }

    return 0;
}

static int process_unmapped_record(bam_bin_bind_t *bind, bam1_t *rec) {
    int flag;
    int pass_post;
    int is_reverse;
    uint64_t gc_bases = 0;
    uint64_t read_bases = 0;

    if (!bind || !rec) return -1;

    flag = rec->core.flag;
    if (bind->require_flags > 0 &&
        (((int64_t)flag & bind->require_flags) != bind->require_flags)) {
        return 0;
    }
    if (bind->exclude_flags > 0 &&
        (((int64_t)flag & bind->exclude_flags) != 0)) {
        return 0;
    }

    if (bind->want_gc || bind->want_mq) {
        bind->unmapped_count_pre++;
    }
    if (bind->want_gc) {
        read_bases = (uint64_t)rec->core.l_qseq;
        gc_bases = count_read_gc_bases(rec);
        bind->unmapped_gc_bases_pre += gc_bases;
        bind->unmapped_bases_pre += read_bases;
    }
    if (bind->want_mq) {
        bind->unmapped_mapq_sum_pre += (uint64_t)rec->core.qual;
    }

    pass_post = 1;
    if (bind->rmdup == BAM_BIN_RMDUP_FLAG && (flag & BAM_FDUP)) {
        pass_post = 0;
    }
    if ((int64_t)rec->core.qual < bind->mapq) {
        pass_post = 0;
    }

    if (!pass_post) return 0;

    bind->unmapped_count_total++;
    is_reverse = (flag & BAM_FREVERSE) != 0;
    if (is_reverse) bind->unmapped_count_rev++;
    else bind->unmapped_count_fwd++;

    if (bind->want_gc) {
        bind->unmapped_gc_bases_post += gc_bases;
        bind->unmapped_bases_post += read_bases;
    }
    if (bind->want_mq) {
        bind->unmapped_mapq_sum_post += (uint64_t)rec->core.qual;
    }

    return 0;
}

static void worker_set_error(bam_bin_worker_shared_t *shared, const char *msg) {
    if (!shared || !msg) return;
    if (__sync_bool_compare_and_swap(&shared->failed, 0, 1)) {
        snprintf(shared->error, sizeof(shared->error), "%s", msg);
    }
}

static void *bam_bin_worker_main(void *arg) {
    bam_bin_worker_arg_t *worker = (bam_bin_worker_arg_t *)arg;
    bam_bin_worker_shared_t *shared = worker ? worker->shared : NULL;
    bam_bin_bind_t *bind = shared ? shared->bind : NULL;
    samFile *fp = NULL;
    sam_hdr_t *hdr = NULL;
    bam1_t *rec = NULL;
    hts_idx_t *idx = NULL;
    char err[512];

    if (!shared || !bind) return NULL;

    fp = sam_open(bind->path, "r");
    if (!fp) {
        worker_set_error(shared, "bam_bin_counts: failed to open alignment file");
        return NULL;
    }
    duckhts_apply_remote_hts_tuning(fp, bind->path,
                                    DUCKHTS_HTS_IO_PROFILE_STREAMING);
    if (bind->reference) {
        if (hts_set_opt(fp, CRAM_OPT_REFERENCE, bind->reference) < 0) {
            worker_set_error(shared, "bam_bin_counts: failed to set CRAM reference");
            sam_close(fp);
            return NULL;
        }
    }
    if (bind->decompression_threads > 0 &&
        hts_set_threads(fp, (int)bind->decompression_threads) < 0) {
        worker_set_error(shared, "bam_bin_counts: failed to configure decompression threads");
        sam_close(fp);
        return NULL;
    }
    hdr = sam_hdr_read(fp);
    if (!hdr) {
        worker_set_error(shared, "bam_bin_counts: failed to read alignment header");
        sam_close(fp);
        return NULL;
    }
    idx = sam_index_load3(fp, bind->path, bind->index_path, HTS_IDX_SILENT_FAIL);
    if (!idx) {
        sam_hdr_destroy(hdr);
        sam_close(fp);
        worker_set_error(shared, "bam_bin_counts: indexed worker requires a BAM/CRAM index");
        return NULL;
    }
    rec = bam_init1();
    if (!rec) {
        sam_hdr_destroy(hdr);
        hts_idx_destroy(idx);
        sam_close(fp);
        worker_set_error(shared, "bam_bin_counts: out of memory allocating BAM record");
        return NULL;
    }

    for (;;) {
        int claim = __sync_fetch_and_add(&shared->next_claim, 1);
        bam_bin_contig_t *contig;
        hts_itr_t *itr;
        bam_bin_stream_state_t state;

        if (shared->failed) break;
        if (claim >= shared->n_claims) break;

        contig = &shared->contigs[shared->claim_order[claim]];
        itr = sam_itr_queryi(idx, contig->tid, 0, HTS_POS_MAX);
        if (!itr) continue;

        memset(&state, 0, sizeof(state));
        contig->final_state = state;

        for (;;) {
            int ret = sam_itr_next(fp, itr, rec);
            if (ret == -1) break;
            if (ret < -1) {
                snprintf(err, sizeof(err), "bam_bin_counts: iterator failure on '%s'", contig->chrom_name);
                worker_set_error(shared, err);
                break;
            }
            if (process_record_into_contig(bind, contig, rec, &state) != 0) {
                worker_set_error(shared, "bam_bin_counts: out of memory while accumulating bins");
                break;
            }
        }

        hts_itr_destroy(itr);
        if (shared->failed) break;
    }

    bam_destroy1(rec);
    sam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(fp);
    return NULL;
}

static int build_selected_contigs(bam_bin_bind_t *bind,
                                  sam_hdr_t *hdr,
                                  hts_idx_t *idx,
                                  char *err,
                                  size_t errlen) {
    int nref;
    int selected = 0;
    int target_tid = -1;
    if (!bind || !hdr) return -1;

    nref = sam_hdr_nref(hdr);
    if (bind->chrom) {
        target_tid = sam_hdr_name2tid(hdr, bind->chrom);
        if (target_tid < 0) {
            snprintf(err, errlen, "bam_bin_counts: unknown chromosome '%s'", bind->chrom);
            return -1;
        }
        bind->contigs = (bam_bin_contig_t *)calloc(1, sizeof(bam_bin_contig_t));
        if (!bind->contigs) {
            snprintf(err, errlen, "bam_bin_counts: out of memory");
            return -1;
        }
        bind->n_contigs = 1;
        bind->contigs[0].tid = target_tid;
        bind->contigs[0].chrom_name = dup_cstr(sam_hdr_tid2name(hdr, target_tid));
        bind->contigs[0].ref_length = (int64_t)sam_hdr_tid2len(hdr, target_tid);
        bind->contigs[0].max_bin = -1;
        if (!bind->contigs[0].chrom_name) {
            snprintf(err, errlen, "bam_bin_counts: out of memory");
            return -1;
        }
        if (idx) {
            uint64_t mapped = 0, unmapped = 0;
            if (hts_idx_get_stat(idx, target_tid, &mapped, &unmapped) == 0) {
                bind->contigs[0].weight = mapped + unmapped;
            }
        }
        return 0;
    }

    bind->contigs = (bam_bin_contig_t *)calloc((size_t)nref, sizeof(bam_bin_contig_t));
    if (!bind->contigs) {
        snprintf(err, errlen, "bam_bin_counts: out of memory");
        return -1;
    }

    for (int tid = 0; tid < nref; tid++) {
        char *chrom_name = dup_cstr(sam_hdr_tid2name(hdr, tid));
        if (!chrom_name) {
            snprintf(err, errlen, "bam_bin_counts: out of memory");
            return -1;
        }
        bind->contigs[selected].tid = tid;
        bind->contigs[selected].chrom_name = chrom_name;
        bind->contigs[selected].ref_length = (int64_t)sam_hdr_tid2len(hdr, tid);
        bind->contigs[selected].max_bin = -1;
        if (idx) {
            uint64_t mapped = 0, unmapped = 0;
            if (hts_idx_get_stat(idx, tid, &mapped, &unmapped) == 0) {
                bind->contigs[selected].weight = mapped + unmapped;
            }
        }
        selected++;
    }

    bind->n_contigs = selected;
    return 0;
}

static int sequential_scan_no_index(bam_bin_bind_t *bind, char *err, size_t errlen) {
    samFile *fp = NULL;
    sam_hdr_t *hdr = NULL;
    bam1_t *rec = NULL;
    int *tid_to_selected = NULL;

    fp = sam_open(bind->path, "r");
    if (!fp) {
        snprintf(err, errlen, "bam_bin_counts: failed to open alignment file");
        return -1;
    }
    duckhts_apply_remote_hts_tuning(fp, bind->path,
                                    DUCKHTS_HTS_IO_PROFILE_STREAMING);
    if (bind->reference) {
        if (hts_set_opt(fp, CRAM_OPT_REFERENCE, bind->reference) < 0) {
            sam_close(fp);
            snprintf(err, errlen, "bam_bin_counts: failed to set CRAM reference");
            return -1;
        }
    }
    if (bind->decompression_threads > 0 &&
        hts_set_threads(fp, (int)bind->decompression_threads) < 0) {
        sam_close(fp);
        snprintf(err, errlen, "bam_bin_counts: failed to configure decompression threads");
        return -1;
    }

    hdr = sam_hdr_read(fp);
    if (!hdr) {
        sam_close(fp);
        snprintf(err, errlen, "bam_bin_counts: failed to read alignment header");
        return -1;
    }

    rec = bam_init1();
    if (!rec) {
        sam_hdr_destroy(hdr);
        sam_close(fp);
        snprintf(err, errlen, "bam_bin_counts: out of memory allocating BAM record");
        return -1;
    }

    tid_to_selected = (int *)malloc(sizeof(int) * (size_t)sam_hdr_nref(hdr));
    if (!tid_to_selected) {
        bam_destroy1(rec);
        sam_hdr_destroy(hdr);
        sam_close(fp);
        snprintf(err, errlen, "bam_bin_counts: out of memory");
        return -1;
    }
    for (int i = 0; i < sam_hdr_nref(hdr); i++) tid_to_selected[i] = -1;
    for (int i = 0; i < bind->n_contigs; i++) tid_to_selected[bind->contigs[i].tid] = i;

    while (1) {
        int ret = sam_read1(fp, hdr, rec);
        if (ret == -1) break;
        if (ret < -1) {
            free(tid_to_selected);
            bam_destroy1(rec);
            sam_hdr_destroy(hdr);
            sam_close(fp);
            snprintf(err, errlen, "bam_bin_counts: read failure while scanning input");
            return -1;
        }
        if (rec->core.tid < 0 || rec->core.pos < 0) {
            if (bind->include_unmapped && process_unmapped_record(bind, rec) != 0) {
                free(tid_to_selected);
                bam_destroy1(rec);
                sam_hdr_destroy(hdr);
                sam_close(fp);
                snprintf(err, errlen, "bam_bin_counts: out of memory while accumulating unmapped reads");
                return -1;
            }
            continue;
        }
        if (rec->core.tid >= sam_hdr_nref(hdr)) continue;
        if (tid_to_selected[rec->core.tid] < 0) continue;

        {
            bam_bin_contig_t *contig = &bind->contigs[tid_to_selected[rec->core.tid]];
            if (process_record_into_contig(bind, contig, rec, &contig->final_state) != 0) {
                free(tid_to_selected);
                bam_destroy1(rec);
                sam_hdr_destroy(hdr);
                sam_close(fp);
                snprintf(err, errlen, "bam_bin_counts: out of memory while accumulating bins");
                return -1;
            }
        }
    }

    free(tid_to_selected);
    bam_destroy1(rec);
    sam_hdr_destroy(hdr);
    sam_close(fp);
    return 0;
}

static int indexed_parallel_scan(bam_bin_bind_t *bind, char *err, size_t errlen) {
    bam_bin_worker_shared_t shared;
    pthread_t *threads = NULL;
    bam_bin_worker_arg_t *args = NULL;
    int n_threads;

    memset(&shared, 0, sizeof(shared));
    shared.bind = bind;
    shared.contigs = bind->contigs;
    shared.claim_order = bind->claim_order;
    shared.n_claims = bind->n_claims;

    n_threads = bind->n_claims;
    if (n_threads > DUCKHTS_MAX_BINCOUNT_THREADS) n_threads = DUCKHTS_MAX_BINCOUNT_THREADS;
    if (n_threads < 1) n_threads = 1;

    threads = (pthread_t *)malloc(sizeof(pthread_t) * (size_t)n_threads);
    args = (bam_bin_worker_arg_t *)malloc(sizeof(bam_bin_worker_arg_t) * (size_t)n_threads);
    if (!threads || !args) {
        free(threads);
        free(args);
        snprintf(err, errlen, "bam_bin_counts: out of memory creating worker threads");
        return -1;
    }

    for (int i = 0; i < n_threads; i++) {
        args[i].shared = &shared;
        if (pthread_create(&threads[i], NULL, bam_bin_worker_main, &args[i]) != 0) {
            shared.failed = 1;
            snprintf(shared.error, sizeof(shared.error), "bam_bin_counts: failed to start worker thread");
            n_threads = i;
            break;
        }
    }
    for (int i = 0; i < n_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    free(threads);
    free(args);

    if (shared.failed) {
        snprintf(err, errlen, "%s", shared.error[0] ? shared.error : "bam_bin_counts: worker failure");
        return -1;
    }
    return 0;
}

static int indexed_unmapped_scan(bam_bin_bind_t *bind, char *err, size_t errlen) {
    samFile *fp = NULL;
    sam_hdr_t *hdr = NULL;
    bam1_t *rec = NULL;
    hts_idx_t *idx = NULL;
    hts_itr_t *itr = NULL;

    if (!bind || !bind->include_unmapped) return 0;

    fp = sam_open(bind->path, "r");
    if (!fp) {
        snprintf(err, errlen, "bam_bin_counts: failed to open alignment file");
        return -1;
    }
    duckhts_apply_remote_hts_tuning(fp, bind->path,
                                    DUCKHTS_HTS_IO_PROFILE_INDEXED_REGION);
    if (bind->reference) {
        if (hts_set_opt(fp, CRAM_OPT_REFERENCE, bind->reference) < 0) {
            sam_close(fp);
            snprintf(err, errlen, "bam_bin_counts: failed to set CRAM reference");
            return -1;
        }
    }
    if (bind->decompression_threads > 0 &&
        hts_set_threads(fp, (int)bind->decompression_threads) < 0) {
        sam_close(fp);
        snprintf(err, errlen, "bam_bin_counts: failed to configure decompression threads");
        return -1;
    }

    hdr = sam_hdr_read(fp);
    if (!hdr) {
        sam_close(fp);
        snprintf(err, errlen, "bam_bin_counts: failed to read alignment header");
        return -1;
    }
    idx = sam_index_load3(fp, bind->path, bind->index_path, HTS_IDX_SILENT_FAIL);
    if (!idx) {
        sam_hdr_destroy(hdr);
        sam_close(fp);
        snprintf(err, errlen, "bam_bin_counts: failed to open BAM/CRAM index");
        return -1;
    }
    itr = sam_itr_queryi(idx, HTS_IDX_NOCOOR, 0, 0);
    if (!itr) {
        hts_idx_destroy(idx);
        sam_hdr_destroy(hdr);
        sam_close(fp);
        return 0;
    }
    rec = bam_init1();
    if (!rec) {
        hts_itr_destroy(itr);
        hts_idx_destroy(idx);
        sam_hdr_destroy(hdr);
        sam_close(fp);
        snprintf(err, errlen, "bam_bin_counts: out of memory allocating BAM record");
        return -1;
    }

    for (;;) {
        int ret = sam_itr_next(fp, itr, rec);
        if (ret == -1) break;
        if (ret < -1) {
            bam_destroy1(rec);
            hts_itr_destroy(itr);
            hts_idx_destroy(idx);
            sam_hdr_destroy(hdr);
            sam_close(fp);
            snprintf(err, errlen, "bam_bin_counts: iterator failure on unmapped reads");
            return -1;
        }
        if (process_unmapped_record(bind, rec) != 0) {
            bam_destroy1(rec);
            hts_itr_destroy(itr);
            hts_idx_destroy(idx);
            sam_hdr_destroy(hdr);
            sam_close(fp);
            snprintf(err, errlen, "bam_bin_counts: out of memory while accumulating unmapped reads");
            return -1;
        }
    }

    bam_destroy1(rec);
    hts_itr_destroy(itr);
    hts_idx_destroy(idx);
    sam_hdr_destroy(hdr);
    sam_close(fp);
    return 0;
}

static void apply_streaming_boundary_repairs(bam_bin_bind_t *bind) {
    bam_bin_stream_state_t prev;
    if (!bind || bind->rmdup != BAM_BIN_RMDUP_STREAMING) return;
    memset(&prev, 0, sizeof(prev));

    for (int i = 0; i < bind->n_contigs; i++) {
        bam_bin_contig_t *contig = &bind->contigs[i];
        if (contig->has_first_record &&
            record_is_streaming_dup(&prev,
                                    contig->first_pos0,
                                    contig->first_is_paired,
                                    contig->first_has_pnext,
                                    contig->first_pnext0)) {
            subtract_first_post_metrics(bind, contig);
        }
        if (contig->final_state.has_pos) {
            prev = contig->final_state;
        }
    }
}

static int64_t dense_bin_count_for_contig(const bam_bin_contig_t *contig, int64_t bin_width) {
    if (!contig || bin_width <= 0 || contig->ref_length <= 0) return 0;
    return (contig->ref_length + bin_width - 1) / bin_width;
}

static inline int64_t contig_count_at(const int64_t *values, const bam_bin_contig_t *contig, int64_t bin) {
    if (!values || !contig || bin < 0 || bin > contig->max_bin) return 0;
    return values[bin];
}

static inline uint64_t contig_u64_at(const uint64_t *values, const bam_bin_contig_t *contig, int64_t bin) {
    if (!values || !contig || bin < 0 || bin > contig->max_bin) return 0;
    return values[bin];
}

static int count_rows_for_bind(const bam_bin_bind_t *bind, size_t *rows_out) {
    size_t rows = 0;
    if (!bind || !rows_out) return -1;
    for (int i = 0; i < bind->n_contigs; i++) {
        const bam_bin_contig_t *contig = &bind->contigs[i];
        int64_t n_bins = dense_bin_count_for_contig(contig, bind->bin_width);
        if (n_bins < 0 || (size_t)n_bins > (SIZE_MAX - rows)) return -1;
        rows += (size_t)n_bins;
    }
    if (bind->include_unmapped) {
        if (rows == SIZE_MAX) return -1;
        rows += 1;
    }
    *rows_out = rows;
    return 0;
}

static int append_rows_from_contigs(bam_bin_bind_t *bind, char *err, size_t errlen) {
    size_t rows_needed = 0;
    size_t row_idx = 0;
    faidx_t *fai = NULL;

    if (!bind) return -1;
    if (count_rows_for_bind(bind, &rows_needed) != 0) {
        snprintf(err, errlen, "bam_bin_counts: row count overflow");
        return -1;
    }
    bind->rows = (bam_bin_row_t *)calloc(rows_needed ? rows_needed : 1, sizeof(bam_bin_row_t));
    if (!bind->rows) {
        snprintf(err, errlen, "bam_bin_counts: out of memory materializing rows");
        return -1;
    }

    if (bind->want_gc && bind->reference) {
        fai = fai_load3_format(bind->reference, NULL, NULL, 0, FAI_FASTA);
        if (!fai) {
            snprintf(err, errlen, "bam_bin_counts: failed to open FASTA index for reference GC");
            return -1;
        }
    }

    for (int i = 0; i < bind->n_contigs; i++) {
        bam_bin_contig_t *contig = &bind->contigs[i];
        uint64_t *ref_gc_by_bin = NULL;
        int64_t n_bins = dense_bin_count_for_contig(contig, bind->bin_width);

        if (fai && n_bins > 0) {
            hts_pos_t fetch_len = 0;
            char *seq = faidx_fetch_seq64(fai, contig->chrom_name, 0, contig->ref_length - 1, &fetch_len);
            if (!seq || fetch_len < 0) {
                fai_destroy(fai);
                snprintf(err, errlen, "bam_bin_counts: failed to fetch FASTA sequence for '%s'", contig->chrom_name);
                free(seq);
                return -1;
            }
            ref_gc_by_bin = (uint64_t *)calloc((size_t)n_bins, sizeof(uint64_t));
            if (!ref_gc_by_bin) {
                fai_destroy(fai);
                free(seq);
                snprintf(err, errlen, "bam_bin_counts: out of memory computing reference GC");
                return -1;
            }
            for (hts_pos_t pos = 0; pos < fetch_len; pos++) {
                unsigned char c = (unsigned char)toupper((unsigned char)seq[pos]);
                if (c == 'C' || c == 'G') {
                    int64_t bin = (int64_t)pos / bind->bin_width;
                    if (bin >= 0 && bin < n_bins) ref_gc_by_bin[bin]++;
                }
            }
            free(seq);
        }

        for (int64_t bin = 0; bin < n_bins; bin++) {
            int64_t pre = contig_count_at(contig->count_pre, contig, bin);
            bam_bin_row_t *row;
            row = &bind->rows[row_idx++];
            row->chrom = contig->chrom_name;
            row->bin_id = bin;
            row->start = bin * bind->bin_width;
            row->end = row->start + bind->bin_width;
            if (contig->ref_length > 0 && row->end > contig->ref_length) {
                row->end = contig->ref_length;
            }
            row->count_total = contig_count_at(contig->count_total, contig, bin);
            row->count_fwd = contig_count_at(contig->count_fwd, contig, bin);
            row->count_rev = contig_count_at(contig->count_rev, contig, bin);
            row->count_pre = pre;
            if (bind->want_gc) {
                row->gc_bases_pre = contig_u64_at(contig->gc_bases_pre, contig, bin);
                row->bases_pre = contig_u64_at(contig->bases_pre, contig, bin);
                row->gc_bases_post = contig_u64_at(contig->gc_bases_post, contig, bin);
                row->bases_post = contig_u64_at(contig->bases_post, contig, bin);
                if (ref_gc_by_bin) {
                    int64_t remaining = contig->ref_length - row->start;
                    row->ref_bases = remaining > bind->bin_width ? (uint64_t)bind->bin_width :
                        (remaining > 0 ? (uint64_t)remaining : 0);
                    row->ref_gc_bases = ref_gc_by_bin[(size_t)bin];
                }
            }
            if (bind->want_mq) {
                row->mapq_sum_pre = contig_u64_at(contig->mapq_sum_pre, contig, bin);
                row->mapq_sum_post = contig_u64_at(contig->mapq_sum_post, contig, bin);
            }
        }

        free(ref_gc_by_bin);
        free_contig_arrays(contig);
    }

    if (bind->include_unmapped) {
        bam_bin_row_t *row = &bind->rows[row_idx++];
        row->chrom = "*";
        row->is_unmapped = 1;
        row->count_total = bind->unmapped_count_total;
        row->count_fwd = bind->unmapped_count_fwd;
        row->count_rev = bind->unmapped_count_rev;
        row->count_pre = bind->unmapped_count_pre;
        row->gc_bases_pre = bind->unmapped_gc_bases_pre;
        row->bases_pre = bind->unmapped_bases_pre;
        row->gc_bases_post = bind->unmapped_gc_bases_post;
        row->bases_post = bind->unmapped_bases_post;
        row->mapq_sum_pre = bind->unmapped_mapq_sum_pre;
        row->mapq_sum_post = bind->unmapped_mapq_sum_post;
    }

    if (fai) fai_destroy(fai);
    bind->row_count = row_idx;
    return 0;
}

static int run_bam_bin_counts(bam_bin_bind_t *bind, char *err, size_t errlen) {
    samFile *fp = NULL;
    sam_hdr_t *hdr = NULL;
    hts_idx_t *idx = NULL;
    int rc = -1;

    fp = sam_open(bind->path, "r");
    if (!fp) {
        snprintf(err, errlen, "bam_bin_counts: failed to open alignment file");
        return -1;
    }
    if (bind->reference) {
        if (hts_set_opt(fp, CRAM_OPT_REFERENCE, bind->reference) < 0) {
            sam_close(fp);
            snprintf(err, errlen, "bam_bin_counts: failed to set CRAM reference");
            return -1;
        }
    }
    hdr = sam_hdr_read(fp);
    if (!hdr) {
        sam_close(fp);
        snprintf(err, errlen, "bam_bin_counts: failed to read alignment header");
        return -1;
    }

    idx = sam_index_load3(fp, bind->path, bind->index_path, HTS_IDX_SILENT_FAIL);

    if (build_selected_contigs(bind, hdr, idx, err, errlen) != 0) goto cleanup;

    if (idx) {
        if (build_claim_order(bind) != 0) {
            snprintf(err, errlen, "bam_bin_counts: out of memory building contig claim order");
            goto cleanup;
        }
        if (indexed_parallel_scan(bind, err, errlen) != 0) goto cleanup;
        apply_streaming_boundary_repairs(bind);
        if (indexed_unmapped_scan(bind, err, errlen) != 0) goto cleanup;
    } else {
        if (sequential_scan_no_index(bind, err, errlen) != 0) goto cleanup;
    }

    if (append_rows_from_contigs(bind, err, errlen) != 0) goto cleanup;

    rc = 0;

cleanup:
    if (idx) hts_idx_destroy(idx);
    if (hdr) sam_hdr_destroy(hdr);
    if (fp) sam_close(fp);
    return rc;
}

static void add_result_columns(duckdb_bind_info info, const bam_bin_bind_t *bind) {
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type double_type = duckdb_create_logical_type(DUCKDB_TYPE_DOUBLE);

    duckdb_bind_add_result_column(info, "chrom", varchar_type);
    duckdb_bind_add_result_column(info, "start", bigint_type);
    duckdb_bind_add_result_column(info, "end", bigint_type);
    duckdb_bind_add_result_column(info, "bin_id", bigint_type);
    duckdb_bind_add_result_column(info, "count_total", bigint_type);
    duckdb_bind_add_result_column(info, "count_fwd", bigint_type);
    duckdb_bind_add_result_column(info, "count_rev", bigint_type);

    if (bind->want_gc || bind->want_mq) {
        duckdb_bind_add_result_column(info, "count_pre", bigint_type);
    }
    if (bind->want_gc) {
        duckdb_bind_add_result_column(info, "gc_bases_pre", bigint_type);
        duckdb_bind_add_result_column(info, "bases_pre", bigint_type);
        duckdb_bind_add_result_column(info, "gc_perc_pre", double_type);
        duckdb_bind_add_result_column(info, "gc_bases_post", bigint_type);
        duckdb_bind_add_result_column(info, "bases_post", bigint_type);
        duckdb_bind_add_result_column(info, "gc_perc_post", double_type);
        if (bind->reference) {
            duckdb_bind_add_result_column(info, "ref_gc_bases", bigint_type);
            duckdb_bind_add_result_column(info, "ref_bases", bigint_type);
            duckdb_bind_add_result_column(info, "gc_perc_ref", double_type);
        }
    }
    if (bind->want_mq) {
        duckdb_bind_add_result_column(info, "mapq_sum_pre", bigint_type);
        duckdb_bind_add_result_column(info, "mean_mapq_pre", double_type);
        duckdb_bind_add_result_column(info, "mapq_sum_post", bigint_type);
        duckdb_bind_add_result_column(info, "mean_mapq_post", double_type);
    }

    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&double_type);
}

static void bam_bin_counts_bind(duckdb_bind_info info) {
    duckdb_value path_val = duckdb_bind_get_parameter(info, 0);
    duckdb_value bw_val = duckdb_bind_get_parameter(info, 1);
    bam_bin_bind_t *bind = NULL;
    duckdb_value val;
    char err[512];
    char *rmdup = NULL;
    char *stats = NULL;

    memset(err, 0, sizeof(err));

    bind = (bam_bin_bind_t *)duckdb_malloc(sizeof(bam_bin_bind_t));
    if (!bind) {
        duckdb_bind_set_error(info, "bam_bin_counts: out of memory");
        duckdb_destroy_value(&path_val);
        duckdb_destroy_value(&bw_val);
        return;
    }
    memset(bind, 0, sizeof(*bind));
    bind->path = duckdb_get_varchar(path_val);
    bind->bin_width = duckdb_get_int64(bw_val);
    bind->mapq = 0;
    bind->require_flags = 0;
    bind->exclude_flags = 0;
    bind->decompression_threads = DUCKHTS_DEFAULT_BINCOUNT_DECOMPRESSION_THREADS;
    bind->rmdup = BAM_BIN_RMDUP_NONE;
    duckdb_destroy_value(&path_val);
    duckdb_destroy_value(&bw_val);

    if (!bind->path || bind->path[0] == '\0') {
        destroy_bam_bin_bind(bind);
        duckdb_bind_set_error(info, "bam_bin_counts requires a non-empty BAM/CRAM path");
        return;
    }
    if (bind->bin_width <= 0) {
        destroy_bam_bin_bind(bind);
        duckdb_bind_set_error(info, "bam_bin_counts requires a positive bin_width");
        return;
    }

    val = duckdb_bind_get_named_parameter(info, "chrom");
    if (val && !duckdb_is_null_value(val)) bind->chrom = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "include_unmapped");
    if (val && !duckdb_is_null_value(val)) bind->include_unmapped = duckdb_get_bool(val) ? 1 : 0;
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "reference");
    if (val && !duckdb_is_null_value(val)) bind->reference = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "index_path");
    if (val && !duckdb_is_null_value(val)) bind->index_path = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "mapq");
    if (val && !duckdb_is_null_value(val)) bind->mapq = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "require_flags");
    if (val && !duckdb_is_null_value(val)) bind->require_flags = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "exclude_flags");
    if (val && !duckdb_is_null_value(val)) bind->exclude_flags = duckdb_get_int64(val);
    if (val) duckdb_destroy_value(&val);

    val = duckdb_bind_get_named_parameter(info, "rmdup");
    if (val && !duckdb_is_null_value(val)) rmdup = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);
    if (!rmdup) rmdup = dup_cstr("none");

    val = duckdb_bind_get_named_parameter(info, "stats");
    if (val && !duckdb_is_null_value(val)) stats = duckdb_get_varchar(val);
    if (val) duckdb_destroy_value(&val);

    if (bind->mapq < 0 || bind->require_flags < 0 || bind->exclude_flags < 0 ||
        bind->require_flags > INT_MAX || bind->exclude_flags > INT_MAX) {
        if (rmdup) duckdb_free(rmdup);
        if (stats) duckdb_free(stats);
        destroy_bam_bin_bind(bind);
        duckdb_bind_set_error(info, "bam_bin_counts: numeric parameters must be non-negative and fit in INT");
        return;
    }
    if (parse_rmdup_mode(rmdup, &bind->rmdup) != 0) {
        if (rmdup) duckdb_free(rmdup);
        if (stats) duckdb_free(stats);
        destroy_bam_bin_bind(bind);
        duckdb_bind_set_error(info, "bam_bin_counts: rmdup must be 'none', 'flag', or 'streaming'");
        return;
    }
    if (parse_stats_string(stats, &bind->want_gc, &bind->want_mq) != 0) {
        if (rmdup) duckdb_free(rmdup);
        if (stats) duckdb_free(stats);
        destroy_bam_bin_bind(bind);
        duckdb_bind_set_error(info, "bam_bin_counts: stats must be NULL, 'gc', 'mq', or 'gc,mq'");
        return;
    }

    if (rmdup) duckdb_free(rmdup);
    if (stats) duckdb_free(stats);

    if (run_bam_bin_counts(bind, err, sizeof(err)) != 0) {
        destroy_bam_bin_bind(bind);
        duckdb_bind_set_error(info, err);
        return;
    }

    add_result_columns(info, bind);
    bind->emitted = 0;
    duckdb_bind_set_bind_data(info, bind, destroy_bam_bin_bind);
}

static void bam_bin_counts_init(duckdb_init_info info) {
    duckdb_init_set_max_threads(info, 1);
}

static void bam_bin_counts_scan(duckdb_function_info info, duckdb_data_chunk output) {
    bam_bin_bind_t *bind = (bam_bin_bind_t *)duckdb_function_get_bind_data(info);
    idx_t out_size = 0;
    idx_t capacity;
    idx_t col = 0;

    if (!bind || bind->emitted >= bind->row_count) {
        duckdb_data_chunk_set_size(output, 0);
        return;
    }

    capacity = duckdb_vector_size();
    out_size = (idx_t)((bind->row_count - bind->emitted) < (size_t)capacity ?
        (bind->row_count - bind->emitted) : (size_t)capacity);

    duckdb_vector chrom_vec = duckdb_data_chunk_get_vector(output, col++);
    duckdb_vector start_vec = duckdb_data_chunk_get_vector(output, col++);
    duckdb_vector end_vec = duckdb_data_chunk_get_vector(output, col++);
    duckdb_vector bin_vec = duckdb_data_chunk_get_vector(output, col++);
    duckdb_vector total_vec = duckdb_data_chunk_get_vector(output, col++);
    duckdb_vector fwd_vec = duckdb_data_chunk_get_vector(output, col++);
    duckdb_vector rev_vec = duckdb_data_chunk_get_vector(output, col++);

    int64_t *start_data = (int64_t *)duckdb_vector_get_data(start_vec);
    int64_t *end_data = (int64_t *)duckdb_vector_get_data(end_vec);
    int64_t *bin_data = (int64_t *)duckdb_vector_get_data(bin_vec);
    int64_t *total_data = (int64_t *)duckdb_vector_get_data(total_vec);
    int64_t *fwd_data = (int64_t *)duckdb_vector_get_data(fwd_vec);
    int64_t *rev_data = (int64_t *)duckdb_vector_get_data(rev_vec);

    duckdb_vector count_pre_vec = NULL;
    duckdb_vector gc_bases_pre_vec = NULL;
    duckdb_vector bases_pre_vec = NULL;
    duckdb_vector gc_perc_pre_vec = NULL;
    duckdb_vector gc_bases_post_vec = NULL;
    duckdb_vector bases_post_vec = NULL;
    duckdb_vector gc_perc_post_vec = NULL;
    duckdb_vector ref_gc_bases_vec = NULL;
    duckdb_vector ref_bases_vec = NULL;
    duckdb_vector gc_perc_ref_vec = NULL;
    duckdb_vector mapq_sum_pre_vec = NULL;
    duckdb_vector mean_mapq_pre_vec = NULL;
    duckdb_vector mapq_sum_post_vec = NULL;
    duckdb_vector mean_mapq_post_vec = NULL;

    int64_t *count_pre_data = NULL;
    int64_t *gc_bases_pre_data = NULL;
    int64_t *bases_pre_data = NULL;
    double *gc_perc_pre_data = NULL;
    int64_t *gc_bases_post_data = NULL;
    int64_t *bases_post_data = NULL;
    double *gc_perc_post_data = NULL;
    int64_t *ref_gc_bases_data = NULL;
    int64_t *ref_bases_data = NULL;
    double *gc_perc_ref_data = NULL;
    int64_t *mapq_sum_pre_data = NULL;
    double *mean_mapq_pre_data = NULL;
    int64_t *mapq_sum_post_data = NULL;
    double *mean_mapq_post_data = NULL;

    if (bind->want_gc || bind->want_mq) {
        count_pre_vec = duckdb_data_chunk_get_vector(output, col++);
        count_pre_data = (int64_t *)duckdb_vector_get_data(count_pre_vec);
    }
    if (bind->want_gc) {
        gc_bases_pre_vec = duckdb_data_chunk_get_vector(output, col++);
        bases_pre_vec = duckdb_data_chunk_get_vector(output, col++);
        gc_perc_pre_vec = duckdb_data_chunk_get_vector(output, col++);
        gc_bases_post_vec = duckdb_data_chunk_get_vector(output, col++);
        bases_post_vec = duckdb_data_chunk_get_vector(output, col++);
        gc_perc_post_vec = duckdb_data_chunk_get_vector(output, col++);

        gc_bases_pre_data = (int64_t *)duckdb_vector_get_data(gc_bases_pre_vec);
        bases_pre_data = (int64_t *)duckdb_vector_get_data(bases_pre_vec);
        gc_perc_pre_data = (double *)duckdb_vector_get_data(gc_perc_pre_vec);
        gc_bases_post_data = (int64_t *)duckdb_vector_get_data(gc_bases_post_vec);
        bases_post_data = (int64_t *)duckdb_vector_get_data(bases_post_vec);
        gc_perc_post_data = (double *)duckdb_vector_get_data(gc_perc_post_vec);

        if (bind->reference) {
            ref_gc_bases_vec = duckdb_data_chunk_get_vector(output, col++);
            ref_bases_vec = duckdb_data_chunk_get_vector(output, col++);
            gc_perc_ref_vec = duckdb_data_chunk_get_vector(output, col++);
            ref_gc_bases_data = (int64_t *)duckdb_vector_get_data(ref_gc_bases_vec);
            ref_bases_data = (int64_t *)duckdb_vector_get_data(ref_bases_vec);
            gc_perc_ref_data = (double *)duckdb_vector_get_data(gc_perc_ref_vec);
        }
    }
    if (bind->want_mq) {
        mapq_sum_pre_vec = duckdb_data_chunk_get_vector(output, col++);
        mean_mapq_pre_vec = duckdb_data_chunk_get_vector(output, col++);
        mapq_sum_post_vec = duckdb_data_chunk_get_vector(output, col++);
        mean_mapq_post_vec = duckdb_data_chunk_get_vector(output, col++);

        mapq_sum_pre_data = (int64_t *)duckdb_vector_get_data(mapq_sum_pre_vec);
        mean_mapq_pre_data = (double *)duckdb_vector_get_data(mean_mapq_pre_vec);
        mapq_sum_post_data = (int64_t *)duckdb_vector_get_data(mapq_sum_post_vec);
        mean_mapq_post_data = (double *)duckdb_vector_get_data(mean_mapq_post_vec);
    }

    for (idx_t i = 0; i < out_size; i++) {
        bam_bin_row_t *row = &bind->rows[bind->emitted + i];
        duckdb_vector_assign_string_element(chrom_vec, i, row->chrom);
        start_data[i] = row->start;
        end_data[i] = row->end;
        bin_data[i] = row->bin_id;
        if (row->is_unmapped) {
            set_null(start_vec, i);
            set_null(end_vec, i);
            set_null(bin_vec, i);
        }
        total_data[i] = row->count_total;
        fwd_data[i] = row->count_fwd;
        rev_data[i] = row->count_rev;

        if (count_pre_data) count_pre_data[i] = row->count_pre;

        if (bind->want_gc) {
            gc_bases_pre_data[i] = (int64_t)row->gc_bases_pre;
            bases_pre_data[i] = (int64_t)row->bases_pre;
            if (row->bases_pre > 0) gc_perc_pre_data[i] = (double)row->gc_bases_pre / (double)row->bases_pre;
            else set_null(gc_perc_pre_vec, i);

            gc_bases_post_data[i] = (int64_t)row->gc_bases_post;
            bases_post_data[i] = (int64_t)row->bases_post;
            if (row->bases_post > 0) gc_perc_post_data[i] = (double)row->gc_bases_post / (double)row->bases_post;
            else set_null(gc_perc_post_vec, i);

            if (bind->reference) {
                ref_gc_bases_data[i] = (int64_t)row->ref_gc_bases;
                ref_bases_data[i] = (int64_t)row->ref_bases;
                if (row->is_unmapped) {
                    set_null(ref_gc_bases_vec, i);
                    set_null(ref_bases_vec, i);
                    set_null(gc_perc_ref_vec, i);
                } else if (row->ref_bases > 0) {
                    gc_perc_ref_data[i] = (double)row->ref_gc_bases / (double)row->ref_bases;
                } else {
                    set_null(gc_perc_ref_vec, i);
                }
            }
        }

        if (bind->want_mq) {
            mapq_sum_pre_data[i] = (int64_t)row->mapq_sum_pre;
            if (row->count_pre > 0) mean_mapq_pre_data[i] = (double)row->mapq_sum_pre / (double)row->count_pre;
            else set_null(mean_mapq_pre_vec, i);

            mapq_sum_post_data[i] = (int64_t)row->mapq_sum_post;
            if (row->count_total > 0) mean_mapq_post_data[i] = (double)row->mapq_sum_post / (double)row->count_total;
            else set_null(mean_mapq_post_vec, i);
        }
    }

    bind->emitted += out_size;
    duckdb_data_chunk_set_size(output, out_size);
}

void register_bam_bin_counts_function(duckdb_connection connection) {
    duckdb_table_function tf = duckdb_create_table_function();
    duckdb_logical_type varchar_type = duckdb_create_logical_type(DUCKDB_TYPE_VARCHAR);
    duckdb_logical_type bigint_type = duckdb_create_logical_type(DUCKDB_TYPE_BIGINT);
    duckdb_logical_type bool_type = duckdb_create_logical_type(DUCKDB_TYPE_BOOLEAN);

    duckdb_table_function_set_name(tf, "bam_bin_counts");
    duckdb_table_function_add_parameter(tf, varchar_type);
    duckdb_table_function_add_parameter(tf, bigint_type);
    duckdb_table_function_add_named_parameter(tf, "chrom", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "include_unmapped", bool_type);
    duckdb_table_function_add_named_parameter(tf, "reference", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "index_path", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "rmdup", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "stats", varchar_type);
    duckdb_table_function_add_named_parameter(tf, "mapq", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "require_flags", bigint_type);
    duckdb_table_function_add_named_parameter(tf, "exclude_flags", bigint_type);

    duckdb_table_function_set_bind(tf, bam_bin_counts_bind);
    duckdb_table_function_set_init(tf, bam_bin_counts_init);
    duckdb_table_function_set_function(tf, bam_bin_counts_scan);

    duckdb_register_table_function(connection, tf);
    duckdb_destroy_table_function(&tf);
    duckdb_destroy_logical_type(&varchar_type);
    duckdb_destroy_logical_type(&bigint_type);
    duckdb_destroy_logical_type(&bool_type);
}
