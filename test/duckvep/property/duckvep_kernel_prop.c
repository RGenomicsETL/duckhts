/*
 * duckvep_kernel_prop.c — property + smoke tests for libduckvep_kernel.
 *
 * Pure C against the kernel's borrowed-view ABI: NO duckdb_extension.h, NO
 * htslib. That is the point — the engine is exercised exactly the way a CLI,
 * conformance runner, or WASM wrapper would, proving it is independently
 * testable across the boundary. theft drives the generated cases and shrinking;
 * greatest hosts the fixed cases. Ported with the kernel from
 * /root/duckvep-c commit 9f922c8.
 *
 *   make test-duckvep-kernel
 *   make test-duckvep-kernel-asan
 *   make test-duckvep-kernel-ubsan
 *   make test-duckvep-kernel-statistical
 */
#include "duckvep_kernel.h"
#include "duckvep_sweep.h"
#include "duckvep_classify.h"
#include "duckvep_effect.h"
#include "duckvep_event.h"
#include "duckvep_so.h"
#include "duckvep_projection.h"
#include "duckvep_codon.h"
#include "duckvep_coding.h"
#include "duckvep_haplotype.h"
#include "duckvep_variant_tile.h"
#include "duckvep_workspace_internal.h"

#include "greatest.h"
#include "theft.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define KPROP_DEFAULT_TRIALS 3000u
#define KPROP_DEFAULT_SEED   UINT64_C(0xd0c0ffee12345678)

/* Shared zero tx_flags for inline scenes (model_open now requires non-NULL flags;
 * the kernel does not read them yet, so all-zero distilled facts are correct). */
static const uint64_t k_zero_flags[8] = {0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u};
#define KPROP_MAX_VARIANTS   64u

static uint64_t kprop_env_u64(const char *name, uint64_t dflt) {
    const char *v = getenv(name);
    char *end = NULL;
    unsigned long long parsed;
    if (v == NULL || v[0] == '\0') return dflt;
    parsed = strtoull(v, &end, 0);
    if (end == NULL || *end != '\0') return dflt;
    return (uint64_t)parsed;
}

static uint64_t kprop_bounded(struct theft *t, uint64_t bound) {
    if (bound <= 1u) return 0u;
    return theft_random_bits(t, 64) % bound; /* modulo bias is fine for tests */
}

/* A generated, self-owning random variant batch. */
struct kprop_batch {
    duckvep_variant_batch_t view;
    uint16_t *chrom;
    uint32_t *pos1;
    uint32_t *end1;
    uint32_t *roff;
    uint16_t *rlen;
    uint32_t *aoff;
    uint16_t *alen;
    uint8_t  *bytes;
    uint8_t  *kind;
};

static void kprop_batch_free(void *instance, void *env) {
    struct kprop_batch *b = (struct kprop_batch *)instance;
    (void)env;
    if (b == NULL) return;
    free(b->chrom);
    free(b->pos1);
    free(b->end1);
    free(b->roff);
    free(b->rlen);
    free(b->aoff);
    free(b->alen);
    free(b->bytes);
    free(b->kind);
    free(b);
}

static enum theft_alloc_res kprop_batch_alloc(struct theft *t, void *env, void **instance) {
    size_t n = (size_t)kprop_bounded(t, (uint64_t)KPROP_MAX_VARIANTS + 1u);
    struct kprop_batch *b = (struct kprop_batch *)calloc(1u, sizeof *b);
    size_t i;
    (void)env;
    if (b == NULL) return THEFT_ALLOC_ERROR;

    if (n > 0u) {
        b->chrom = (uint16_t *)calloc(n, sizeof *b->chrom);
        b->pos1  = (uint32_t *)calloc(n, sizeof *b->pos1);
        b->end1  = (uint32_t *)calloc(n, sizeof *b->end1);
        b->roff  = (uint32_t *)calloc(n, sizeof *b->roff);
        b->rlen  = (uint16_t *)calloc(n, sizeof *b->rlen);
        b->aoff  = (uint32_t *)calloc(n, sizeof *b->aoff);
        b->alen  = (uint16_t *)calloc(n, sizeof *b->alen);
        b->kind  = (uint8_t  *)calloc(n, sizeof *b->kind);
        b->bytes = (uint8_t  *)calloc(n * 2u, 1u);
        if (b->chrom == NULL || b->pos1 == NULL || b->end1 == NULL ||
            b->roff == NULL || b->rlen == NULL || b->aoff == NULL ||
            b->alen == NULL || b->kind == NULL || b->bytes == NULL) {
            kprop_batch_free(b, NULL);
            return THEFT_ALLOC_ERROR;
        }
        for (i = 0; i < n; i++) {
            uint32_t p = (uint32_t)theft_random_bits(t, 28) + 1u;
            uint32_t span = (uint32_t)kprop_bounded(t, 16u);
            b->chrom[i] = (uint16_t)theft_random_bits(t, 10);
            b->pos1[i]  = p;
            b->end1[i]  = p + span;
            b->kind[i]  = (uint8_t)kprop_bounded(t, 6u); /* duckvep_variant_kind range */
            b->roff[i]  = 0u;
            b->rlen[i]  = 1u;
            b->aoff[i]  = 1u;
            b->alen[i]  = 1u;
        }
    }

    b->view.chrom_id     = b->chrom;
    b->view.pos1         = b->pos1;
    b->view.end1         = b->end1;
    b->view.ref_offset   = b->roff;
    b->view.ref_length   = b->rlen;
    b->view.alt_offset   = b->aoff;
    b->view.alt_length   = b->alen;
    b->view.allele_bytes = b->bytes;
    b->view.allele_bytes_len = n * 2u;
    b->view.variant_kind = b->kind;
    b->view.count        = n;

    *instance = b;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_batch_info = {
    .alloc = kprop_batch_alloc,
    .free  = kprop_batch_free,
};

/* Property: a NULL model is rejected with INVALID_ARG *before* the batch is
 * read, for any generated batch. ASan/UBSan turn this into a no-OOB-read proof
 * regardless of batch contents — a real guard as predicate code is added. */
static enum theft_trial_res prop_null_model_rejected(struct theft *t, void *arg1) {
    const struct kprop_batch *b = (const struct kprop_batch *)arg1;
    duckvep_error_t err;
    duckvep_status_t s;
    (void)t;
    memset(&err, 0, sizeof err);
    s = duckvep_annotate_tile(NULL, &b->view, NULL, NULL, NULL, &err);
    if (s != DUCKVEP_ERR_INVALID_ARG) return THEFT_TRIAL_FAIL;
    if (err.status != DUCKVEP_ERR_INVALID_ARG) return THEFT_TRIAL_FAIL;
    return THEFT_TRIAL_PASS;
}

TEST kernel_version_is_well_formed(void) {
    char expected[16];
    const char *v = duckvep_kernel_version();
    ASSERT(v != NULL);
    (void)snprintf(expected, sizeof expected, "%d.%d.%d",
                   DUCKVEP_KERNEL_VERSION_MAJOR,
                   DUCKVEP_KERNEL_VERSION_MINOR,
                   DUCKVEP_KERNEL_VERSION_PATCH);
    ASSERT_STR_EQ(expected, v);
    PASS();
}

TEST null_args_are_rejected(void) {
    duckvep_error_t err;
    memset(&err, 0, sizeof err);
    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
              duckvep_annotate_tile(NULL, NULL, NULL, NULL, NULL, &err));
    PASS();
}

TEST model_open_rejects_projection_and_sequence_invariant_mutations(void) {
    uint16_t chrom[1] = {0u};
    uint32_t tx_start[1] = {100u};
    uint32_t tx_end[1] = {111u};
    int8_t strand[1] = {1};
    uint64_t flags[1] = {DUCKVEP_TX_HAS_TRANSLATION |
                         DUCKVEP_TX_BIOTYPE_PROTEIN_CODING};
    uint32_t exon_offset[1] = {0u};
    uint16_t exon_count[1] = {2u};
    uint32_t cds_start[1] = {100u};
    uint32_t cds_end[1] = {111u};
    uint32_t exon_start[2] = {100u, 108u};
    uint32_t exon_end[2] = {104u, 111u};
    uint32_t cdna_start[2] = {1u, 6u};
    uint32_t cdna_end[2] = {5u, 9u};
    int8_t phase[2] = {0, 0};
    int8_t end_phase[2] = {0, 0};
    uint8_t cds_bytes[9] = {'A','T','G','A','A','A','T','A','A'};
    uint64_t sequence_offset[1] = {0u};
    uint32_t sequence_length[1] = {9u};
    uint8_t codon_table[1] = {1u};
    uint8_t invalid_tail[DUCKVEP_POST_CDS_BASE_COUNT] = {'A','X',0u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_model_t *model = NULL;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&err, 0, sizeof err);
    tx.chrom_id = chrom; tx.start1 = tx_start; tx.end1 = tx_end;
    tx.strand = strand; tx.flags = flags; tx.exon_offset = exon_offset;
    tx.exon_count = exon_count; tx.cds_start1 = cds_start;
    tx.cds_end1 = cds_end; tx.transcript_count = 1u;
    exons.start1 = exon_start; exons.end1 = exon_end;
    exons.cdna_start1 = cdna_start; exons.cdna_end1 = cdna_end;
    exons.phase = phase; exons.end_phase = end_phase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = sequence_offset; seq.cds_length = sequence_length;
    seq.codon_table = codon_table; seq.transcript_count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    duckvep_model_close(model); model = NULL;

    strand[0] = 0;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(66u, err.where_code); strand[0] = 1;

    tx_start[0] = 99u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(67u, err.where_code); tx_start[0] = 100u;

    cdna_start[1] = 7u; cdna_end[1] = 10u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(68u, err.where_code); cdna_start[1] = 6u; cdna_end[1] = 9u;

    phase[0] = 3;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(69u, err.where_code); phase[0] = 0;

    cds_start[0] = 105u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(70u, err.where_code); cds_start[0] = 100u;

    sequence_length[0] = 8u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(71u, err.where_code); sequence_length[0] = 9u;

    codon_table[0] = 3u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(71u, err.where_code); codon_table[0] = 1u;

    cds_bytes[4] = 'X';
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(71u, err.where_code); cds_bytes[4] = 'A';

    seq.post_cds_bases = invalid_tail;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(71u, err.where_code); seq.post_cds_bases = NULL;

    exons.phase = NULL; exons.end_phase = NULL;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(71u, err.where_code);
    PASS();
}

TEST annotate_tile_rejects_null_model_for_any_batch(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile rejects NULL model without reading the batch";
    cfg.prop1 = prop_null_model_rejected;
    cfg.type_info[0] = &kprop_batch_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

/* ===================================================================== *
 * Sweep candidate generation vs a brute-force O(N*T) oracle.
 *
 * The strongest available check: for EVERY generated (variants, transcripts)
 * scene, the sorted sweep must emit exactly the same candidate-pair SET as the
 * naive double loop. Catches off-by-one in the window, chrom-boundary leaks,
 * active-set evict bugs, and swap-remove ordering hazards.
 * ===================================================================== */

#define KPROP_MAX_TX     32u
#define KPROP_SWEEP_HALO 100u
#define KPROP_COORD_MAX  2000u
#define KPROP_NCHROM     3u
#define KPROP_MAX_PAIRS  (KPROP_MAX_VARIANTS * KPROP_MAX_TX)

struct kprop_scene {
    duckvep_variant_batch_t    v;
    duckvep_transcript_model_t tx;
    uint32_t                   halo;  /* per-scene, so the window size is fuzzed */
    /* variant SoA */
    uint16_t *vchrom; uint32_t *vpos; uint32_t *vend; uint8_t *vkind;
    uint32_t *vroff;  uint16_t *vrlen; uint32_t *vaoff; uint16_t *valen; uint8_t *vbytes;
    /* transcript SoA */
    uint16_t *tchrom; uint32_t *tstart; uint32_t *tend; int8_t *tstrand; uint64_t *tflags;
    uint32_t *texoff; uint16_t *texcnt; uint32_t *tcds_s; uint32_t *tcds_e;
};

struct vrec { uint16_t chrom; uint32_t pos; uint32_t span; uint8_t kind; };
struct trec { uint16_t chrom; uint32_t start; uint32_t len; int8_t strand; };

static int vrec_cmp(const void *a, const void *b) {
    const struct vrec *x = (const struct vrec *)a, *y = (const struct vrec *)b;
    if (x->chrom != y->chrom) return x->chrom < y->chrom ? -1 : 1;
    if (x->pos != y->pos)     return x->pos < y->pos ? -1 : 1;
    return 0;
}
static int trec_cmp(const void *a, const void *b) {
    const struct trec *x = (const struct trec *)a, *y = (const struct trec *)b;
    if (x->chrom != y->chrom) return x->chrom < y->chrom ? -1 : 1;
    if (x->start != y->start) return x->start < y->start ? -1 : 1;
    return 0;
}
static int u64_cmp(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}

static void kprop_scene_free(void *instance, void *env) {
    struct kprop_scene *s = (struct kprop_scene *)instance;
    (void)env;
    if (s == NULL) return;
    free(s->vchrom); free(s->vpos); free(s->vend); free(s->vkind);
    free(s->vroff); free(s->vrlen); free(s->vaoff); free(s->valen); free(s->vbytes);
    free(s->tchrom); free(s->tstart); free(s->tend); free(s->tstrand); free(s->tflags);
    free(s->texoff); free(s->texcnt); free(s->tcds_s); free(s->tcds_e);
    free(s);
}

static enum theft_alloc_res kprop_scene_alloc(struct theft *t, void *env, void **instance) {
    size_t nv = (size_t)kprop_bounded(t, (uint64_t)KPROP_MAX_VARIANTS + 1u);
    size_t ntx = (size_t)kprop_bounded(t, (uint64_t)KPROP_MAX_TX + 1u);
    struct kprop_scene *s = (struct kprop_scene *)calloc(1u, sizeof *s);
    struct vrec vr[KPROP_MAX_VARIANTS];
    struct trec tr[KPROP_MAX_TX];
    size_t i;
    /* Fuzz the window size (incl. 0 and large) and the coordinate scale. A high
     * but overflow-safe base exercises large coordinates; clustering variants and
     * transcripts within [base, base+span] keeps overlaps dense (non-vacuous),
     * while letting transcript length reach the full span gives long transcripts
     * that stay active across many variants (active-set pressure). */
    static const uint32_t halos[5] = {0u, 1u, 50u, 100u, 5000u};
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFFF0000u) + 1u;
    uint32_t span = (uint32_t)kprop_bounded(t, 5000u) + 1u;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->halo = halos[kprop_bounded(t, 5u)];

    for (i = 0; i < nv; i++) {
        vr[i].chrom = (uint16_t)kprop_bounded(t, KPROP_NCHROM);
        vr[i].pos   = base + (uint32_t)kprop_bounded(t, span);
        vr[i].span  = (uint32_t)kprop_bounded(t, 8u);
        vr[i].kind  = kprop_bounded(t, 4u) == 0u
                    ? (uint8_t)DUCKVEP_KIND_SV
                    : (uint8_t)DUCKVEP_KIND_SNV;
    }
    for (i = 0; i < ntx; i++) {
        tr[i].chrom  = (uint16_t)kprop_bounded(t, KPROP_NCHROM);
        tr[i].start  = base + (uint32_t)kprop_bounded(t, span);
        tr[i].len    = (uint32_t)kprop_bounded(t, span); /* up to a full-span transcript */
        tr[i].strand = (kprop_bounded(t, 2u) == 0u) ? (int8_t)1 : (int8_t)-1;
    }
    qsort(vr, nv, sizeof vr[0], vrec_cmp);
    qsort(tr, ntx, sizeof tr[0], trec_cmp);

    s->vchrom = (uint16_t *)calloc(nv ? nv : 1u, sizeof *s->vchrom);
    s->vpos   = (uint32_t *)calloc(nv ? nv : 1u, sizeof *s->vpos);
    s->vend   = (uint32_t *)calloc(nv ? nv : 1u, sizeof *s->vend);
    s->vkind  = (uint8_t  *)calloc(nv ? nv : 1u, sizeof *s->vkind);
    s->vroff  = (uint32_t *)calloc(nv ? nv : 1u, sizeof *s->vroff);
    s->vrlen  = (uint16_t *)calloc(nv ? nv : 1u, sizeof *s->vrlen);
    s->vaoff  = (uint32_t *)calloc(nv ? nv : 1u, sizeof *s->vaoff);
    s->valen  = (uint16_t *)calloc(nv ? nv : 1u, sizeof *s->valen);
    s->vbytes = (uint8_t  *)calloc(nv ? nv : 1u, 1u);
    s->tchrom = (uint16_t *)calloc(ntx ? ntx : 1u, sizeof *s->tchrom);
    s->tstart = (uint32_t *)calloc(ntx ? ntx : 1u, sizeof *s->tstart);
    s->tend   = (uint32_t *)calloc(ntx ? ntx : 1u, sizeof *s->tend);
    s->tstrand= (int8_t   *)calloc(ntx ? ntx : 1u, sizeof *s->tstrand);
    s->tflags = (uint64_t *)calloc(ntx ? ntx : 1u, sizeof *s->tflags);
    s->texoff = (uint32_t *)calloc(ntx ? ntx : 1u, sizeof *s->texoff);
    s->texcnt = (uint16_t *)calloc(ntx ? ntx : 1u, sizeof *s->texcnt);
    s->tcds_s = (uint32_t *)calloc(ntx ? ntx : 1u, sizeof *s->tcds_s);
    s->tcds_e = (uint32_t *)calloc(ntx ? ntx : 1u, sizeof *s->tcds_e);
    if (!s->vchrom || !s->vpos || !s->vend || !s->vkind || !s->vroff || !s->vrlen ||
        !s->vaoff || !s->valen || !s->vbytes || !s->tchrom || !s->tstart || !s->tend ||
        !s->tstrand || !s->tflags || !s->texoff || !s->texcnt || !s->tcds_s || !s->tcds_e) {
        kprop_scene_free(s, NULL);
        return THEFT_ALLOC_ERROR;
    }

    for (i = 0; i < nv; i++) {
        s->vchrom[i] = vr[i].chrom;
        s->vpos[i]   = vr[i].pos;
        s->vkind[i]  = vr[i].kind;
        s->vend[i]   = vr[i].kind == (uint8_t)DUCKVEP_KIND_SV
                     ? vr[i].pos + vr[i].span
                     : vr[i].pos;
        s->vrlen[i]  = 1u; s->valen[i] = 1u; s->vaoff[i] = 1u;
    }
    for (i = 0; i < ntx; i++) {
        s->tchrom[i]  = tr[i].chrom;
        s->tstart[i]  = tr[i].start;
        s->tend[i]    = tr[i].start + tr[i].len;
        s->tstrand[i] = tr[i].strand;
    }

    s->v.chrom_id = s->vchrom; s->v.pos1 = s->vpos; s->v.end1 = s->vend;
    /* Geometry-only scene: omit all allele columns. It generates only SNV/SV rows;
     * allele-aware differing-region behavior is covered by dedicated tests below. */
    s->v.ref_offset = NULL; s->v.ref_length = NULL;
    s->v.alt_offset = NULL; s->v.alt_length = NULL;
    s->v.allele_bytes = NULL; s->v.allele_bytes_len = 0u;
    s->v.variant_kind = s->vkind; s->v.count = nv;

    s->tx.chrom_id = s->tchrom; s->tx.start1 = s->tstart; s->tx.end1 = s->tend;
    s->tx.strand = s->tstrand; s->tx.flags = s->tflags;
    s->tx.exon_offset = s->texoff; s->tx.exon_count = s->texcnt;
    s->tx.cds_start1 = s->tcds_s; s->tx.cds_end1 = s->tcds_e;
    s->tx.transcript_count = ntx;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_scene_info = {
    .alloc = kprop_scene_alloc,
    .free  = kprop_scene_free,
};

struct pair_collector { uint64_t *buf; size_t n; size_t cap; };
static int pair_sink(uint32_t vi, uint32_t ti, void *ctx) {
    struct pair_collector *c = (struct pair_collector *)ctx;
    uint64_t key = ((uint64_t)vi << 32) | (uint64_t)ti;
    if (c->n < c->cap) c->buf[c->n] = key;
    c->n++;
    return 1;
}

/* Naive O(N*T) candidate set using the same window predicate (64-bit math; the
 * generated coords are far from UINT32 limits so there is no saturation to
 * model here). */
static size_t brute_candidates(const struct kprop_scene *s, uint32_t halo,
                               uint64_t *out, size_t cap) {
    size_t n = 0u, vi, ti;
    for (vi = 0u; vi < s->v.count; vi++) {
        uint64_t start = (uint64_t)s->v.pos1[vi];
        uint64_t end = (uint64_t)duckvep_event_effective_end1_at(&s->v, vi);
        for (ti = 0u; ti < s->tx.transcript_count; ti++) {
            if (s->v.chrom_id[vi] != s->tx.chrom_id[ti]) continue;
            if ((uint64_t)s->tx.start1[ti] <= end + halo &&
                (uint64_t)s->tx.end1[ti] + halo >= start) {
                if (n < cap) out[n] = ((uint64_t)vi << 32) | (uint64_t)ti;
                n++;
            }
        }
    }
    return n;
}

static enum theft_trial_res prop_sweep_matches_bruteforce(struct theft *t, void *arg1) {
    const struct kprop_scene *s = (const struct kprop_scene *)arg1;
    uint64_t sweep_buf[KPROP_MAX_PAIRS];
    uint64_t brute_buf[KPROP_MAX_PAIRS];
    uint32_t active[KPROP_MAX_TX];
    uint32_t candidates[KPROP_MAX_TX];
    struct pair_collector col;
    duckvep_status_t st = DUCKVEP_OK;
    size_t n_sweep, n_brute, i;
    (void)t;

    col.buf = sweep_buf; col.n = 0u; col.cap = KPROP_MAX_PAIRS;
    n_sweep = duckvep_sweep_candidates(&s->v, &s->tx, s->halo,
                                       active, KPROP_MAX_TX,
                                       candidates, KPROP_MAX_TX,
                                       pair_sink, &col, &st);
    if (st != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (n_sweep != col.n) return THEFT_TRIAL_FAIL; /* return value vs sink count */

    n_brute = brute_candidates(s, s->halo, brute_buf, KPROP_MAX_PAIRS);
    if (n_sweep != n_brute) return THEFT_TRIAL_FAIL;

    qsort(sweep_buf, n_sweep, sizeof sweep_buf[0], u64_cmp);
    qsort(brute_buf, n_brute, sizeof brute_buf[0], u64_cmp);
    for (i = 0u; i < n_sweep; i++) {
        if (sweep_buf[i] != brute_buf[i]) return THEFT_TRIAL_FAIL;
    }
    return THEFT_TRIAL_PASS;
}

static uint32_t kprop_sat_add_u32(uint32_t a, uint32_t b) {
    return b > UINT32_MAX - a ? UINT32_MAX : a + b;
}

static enum theft_trial_res prop_seeded_sweep_matches_bruteforce(
    struct theft *t, void *arg1) {
    const struct kprop_scene *s = (const struct kprop_scene *)arg1;
    uint64_t sweep_buf[KPROP_MAX_PAIRS];
    uint64_t brute_buf[KPROP_MAX_PAIRS];
    uint32_t active[KPROP_MAX_TX];
    uint32_t candidates[KPROP_MAX_TX];
    uint32_t seed[KPROP_MAX_TX];
    duckvep_sweep_cursor_t cursor;
    size_t seed_count = 0u;
    size_t n_sweep = 0u;
    size_t n_brute;
    size_t i;
    (void)t;

    if (s->v.count == 0u) return THEFT_TRIAL_PASS;
    for (i = 0u; i < s->tx.transcript_count; i++) {
        if (s->tx.chrom_id[i] != s->v.chrom_id[0]) continue;
        if (s->tx.start1[i] <= kprop_sat_add_u32(s->v.pos1[0], s->halo) &&
            kprop_sat_add_u32(s->tx.end1[i], s->halo) >= s->v.pos1[0]) {
            seed[seed_count++] = (uint32_t)i;
        }
    }

    duckvep_sweep_cursor_init(&cursor, &s->v, &s->tx, s->halo,
                              active, KPROP_MAX_TX,
                              candidates, KPROP_MAX_TX);
    if (!duckvep_sweep_cursor_seed(&cursor, seed, seed_count)) {
        return THEFT_TRIAL_FAIL;
    }
    for (;;) {
        uint32_t vi;
        const uint32_t *tx_indices;
        size_t tx_count;
        size_t j;
        int rc = duckvep_sweep_cursor_next(&cursor, &vi, &tx_indices,
                                           &tx_count);
        if (rc < 0) return THEFT_TRIAL_FAIL;
        if (rc == 0) break;
        for (j = 0u; j < tx_count; j++) {
            if (n_sweep >= KPROP_MAX_PAIRS) return THEFT_TRIAL_FAIL;
            sweep_buf[n_sweep++] = ((uint64_t)vi << 32) | tx_indices[j];
        }
    }
    n_brute = brute_candidates(s, s->halo, brute_buf, KPROP_MAX_PAIRS);
    if (n_sweep != n_brute) return THEFT_TRIAL_FAIL;
    qsort(sweep_buf, n_sweep, sizeof sweep_buf[0], u64_cmp);
    qsort(brute_buf, n_brute, sizeof brute_buf[0], u64_cmp);
    for (i = 0u; i < n_sweep; i++) {
        if (sweep_buf[i] != brute_buf[i]) return THEFT_TRIAL_FAIL;
    }
    return THEFT_TRIAL_PASS;
}

TEST seeded_sweep_matches_bruteforce_for_any_scene(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "cgranges-seeded first event + sweep == brute-force candidates";
    cfg.prop1 = prop_seeded_sweep_matches_bruteforce;
    cfg.type_info[0] = &kprop_scene_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

/* Deterministic anchor: a hand-computed scene so the property above cannot pass
 * vacuously (e.g. by always emitting zero pairs). halo = 100.
 *   transcripts: t0 chr0[100,200]  t1 chr0[500,600]  t2 chr1[100,200]
 *   variants:    v0 chr0@150 -> {t0}   v1 chr0@650 -> {t1}   v2 chr1@150 -> {t2}
 * Expected exactly (v0,t0),(v1,t1),(v2,t2). Also exercises chrom boundaries and
 * the evict-then-emit ordering (t0 is admitted then evicted before v1 emits). */
TEST event_load_without_variant_kind_uses_supplied_interval(void) {
    static const uint16_t chrom[1] = {7u};
    static const uint32_t pos[1] = {100u};
    static const uint32_t end[1] = {200u};
    duckvep_variant_batch_t v;
    duckvep_event_t e;
    memset(&v, 0, sizeof v);
    v.chrom_id = chrom;
    v.pos1 = pos;
    v.end1 = end;
    v.count = 1u;

    duckvep_event_load(&v, 0u, &e);
    ASSERT_EQ(7u, e.chrom_id);
    ASSERT_EQ(100u, e.raw_start1);
    ASSERT_EQ(200u, e.raw_end1);
    ASSERT_EQ(100u, e.feature_start1);
    ASSERT_EQ(200u, e.feature_end1);
    ASSERT_EQ(100u, e.start1);
    ASSERT_EQ(200u, e.end1);
    ASSERT_EQ((uint8_t)DUCKVEP_KIND_SV, e.kind);
    ASSERT_EQ((uint8_t)DUCKVEP_SV_UNKNOWN, e.sv_type);
    ASSERT_EQ((uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN, e.copy_change);
    PASS();
}

TEST event_load_trims_small_variant_differing_region(void) {
    static const uint16_t chrom[4] = {0u, 0u, 0u, 0u};
    static const uint32_t pos[4]   = {240u, 240u, 100u, 200u};
    static const uint32_t end[4]   = {241u, 240u, 102u, 201u};
    static const uint8_t  kind[4]  = {
        (uint8_t)DUCKVEP_KIND_DEL,
        (uint8_t)DUCKVEP_KIND_INS,
        (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV
    };
    static const uint8_t  bytes[16] = {
        'T','A', 'T',      /* TA>T: deleted base is POS+1 */
        'T', 'T','G',      /* T>TG: pure insertion after POS */
        'A','C','G', 'A','T','G', /* ACG>ATG: one changed middle base */
        'G','T', 'C','A'   /* GT>CA: unpadded two-base MNV */
    };
    static const uint32_t roff[4] = {0u, 3u, 6u, 12u};
    static const uint32_t aoff[4] = {2u, 4u, 9u, 14u};
    static const uint16_t rlen[4] = {2u, 1u, 3u, 2u};
    static const uint16_t alen[4] = {1u, 2u, 3u, 2u};
    static const uint32_t exp_start[4] = {241u, 240u, 101u, 200u};
    static const uint32_t exp_end[4]   = {241u, 240u, 101u, 201u};
    static const uint32_t exp_feature_start[4] = {241u, 241u, 100u, 200u};
    static const uint32_t exp_feature_end[4]   = {241u, 240u, 102u, 201u};
    static const uint32_t exp_effective_end[4] = {241u, 241u, 102u, 201u};
    static const uint16_t exp_roff[4]  = {1u, 1u, 1u, 0u};
    static const uint16_t exp_aoff[4]  = {1u, 1u, 1u, 0u};
    static const uint16_t exp_rlen[4]  = {1u, 0u, 1u, 2u};
    static const uint16_t exp_alen[4]  = {0u, 1u, 1u, 2u};
    static const uint8_t  exp_inter[4] = {0u, 1u, 0u, 0u};
    static const uint8_t  exp_kind[4] = {
        (uint8_t)DUCKVEP_KIND_DEL,
        (uint8_t)DUCKVEP_KIND_INS,
        (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_MNV
    };
    duckvep_variant_batch_t v;
    duckvep_event_t e;
    size_t i;

    memset(&v, 0, sizeof v);
    v.chrom_id = chrom; v.pos1 = pos; v.end1 = end; v.variant_kind = kind;
    v.allele_bytes = bytes; v.allele_bytes_len = sizeof bytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 4u;

    for (i = 0u; i < 4u; i++) {
        duckvep_event_load(&v, i, &e);
        ASSERT_EQ(pos[i], e.raw_start1);
        ASSERT_EQ(end[i], e.raw_end1);
        ASSERT_EQ(exp_feature_start[i], e.feature_start1);
        ASSERT_EQ(exp_feature_end[i], e.feature_end1);
        ASSERT_EQ(exp_start[i], e.start1);
        ASSERT_EQ(exp_end[i], e.end1);
        ASSERT_EQ(exp_roff[i], e.ref_diff_offset);
        ASSERT_EQ(exp_aoff[i], e.alt_diff_offset);
        ASSERT_EQ(exp_rlen[i], e.ref_diff_length);
        ASSERT_EQ(exp_alen[i], e.alt_diff_length);
        ASSERT_EQ(exp_inter[i], e.interbase);
        ASSERT_EQ(exp_kind[i], e.kind);
        ASSERT_EQ(exp_effective_end[i], duckvep_event_effective_end1_at(&v, i));
    }
    PASS();
}

TEST event_prepare_small_preserves_anchor_side_and_semantic_kind(void) {
    static const uint8_t delete_first_ref[2] = {'A', 'C'};
    static const uint8_t delete_first_alt[1] = {'C'};
    static const uint8_t insert_first_ref[1] = {'C'};
    static const uint8_t insert_first_alt[2] = {'A', 'C'};
    static const uint8_t padded_snv_ref[3] = {'G', 'A', 'C'};
    static const uint8_t padded_snv_alt[3] = {'G', 'A', 'T'};
    duckvep_event_t event;

    ASSERT(duckvep_event_prepare_small(1u, delete_first_ref, 2u,
                                      delete_first_alt, 1u, &event));
    ASSERT_EQ((uint8_t)DUCKVEP_KIND_DEL, event.kind);
    ASSERT_EQ(1u, event.feature_start1);
    ASSERT_EQ(2u, event.feature_end1);
    ASSERT_EQ(1u, event.start1);
    ASSERT_EQ(1u, event.end1);
    ASSERT_EQ(0u, event.ref_diff_offset);
    ASSERT_EQ(1u, event.ref_diff_length);
    ASSERT_EQ(0u, event.alt_diff_length);

    ASSERT(duckvep_event_prepare_small(1u, insert_first_ref, 1u,
                                      insert_first_alt, 2u, &event));
    ASSERT_EQ((uint8_t)DUCKVEP_KIND_INS, event.kind);
    ASSERT_EQ(1u, event.feature_start1);
    ASSERT_EQ(1u, event.feature_end1);
    ASSERT_EQ(0u, event.insertion_boundary0);
    ASSERT_EQ(1u, event.start1);
    ASSERT_EQ((uint8_t)DUCKVEP_EVENT_ANCHOR_RIGHT, event.anchor_side);
    ASSERT_EQ(0u, event.anchor_ref_offset);
    ASSERT_EQ(0u, event.ref_diff_length);
    ASSERT_EQ(1u, event.alt_diff_length);

    ASSERT(duckvep_event_prepare_small(10u, padded_snv_ref, 3u,
                                      padded_snv_alt, 3u, &event));
    ASSERT_EQ((uint8_t)DUCKVEP_KIND_SNV, event.kind);
    ASSERT_EQ(10u, event.feature_start1);
    ASSERT_EQ(12u, event.feature_end1);
    ASSERT_EQ(12u, event.start1);
    ASSERT_EQ(12u, event.end1);
    ASSERT_EQ(2u, event.ref_diff_offset);
    ASSERT_EQ(2u, event.alt_diff_offset);
    ASSERT_EQ(1u, event.ref_diff_length);
    ASSERT_EQ(1u, event.alt_diff_length);
    PASS();
}

TEST sweep_small_variant_differing_region_tail_is_not_persistent(void) {
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vstart[2] = {100u, 101u};
    static const uint32_t vend[2]   = {103u, 101u};
    static const uint8_t  vkind[2]  = {(uint8_t)DUCKVEP_KIND_DEL, (uint8_t)DUCKVEP_KIND_SNV};
    static const uint8_t  bytes[7]  = {'A','C','G','T', 'A', 'C','T'}; /* ACGT>A, C>T */
    static const uint32_t roff[2]   = {0u, 5u};
    static const uint32_t aoff[2]   = {4u, 6u};
    static const uint16_t rlen[2]   = {4u, 1u};
    static const uint16_t alen[2]   = {1u, 1u};
    static const uint16_t tchrom[2] = {0u, 0u};
    static const uint32_t tstart[2] = {100u, 103u};
    static const uint32_t tend[2]   = {100u, 103u};
    duckvep_variant_batch_t v;
    duckvep_transcript_model_t tx;
    duckvep_sweep_cursor_t cursor;
    uint32_t active[2];
    uint32_t candidates[2];
    uint32_t vi;
    const uint32_t *slice;
    size_t n;

    memset(&v, 0, sizeof v); memset(&tx, 0, sizeof tx);
    v.chrom_id = vchrom; v.pos1 = vstart; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = bytes; v.allele_bytes_len = sizeof bytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 2u;
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.transcript_count = 2u;

    duckvep_sweep_cursor_init(&cursor, &v, &tx, 0u, active, 2u, candidates, 2u);
    ASSERT_EQ(1, duckvep_sweep_cursor_next(&cursor, &vi, &slice, &n));
    ASSERT_EQ(0u, vi);
    ASSERT_EQ(2u, n);       /* raw anchor t0 + differing-region tail t1 */
    ASSERT_EQ(1u, cursor.nact); /* t1 was not retained in the point active set */
    ASSERT_EQ(0u, slice[0]);
    ASSERT_EQ(1u, slice[1]);

    ASSERT_EQ(1, duckvep_sweep_cursor_next(&cursor, &vi, &slice, &n));
    ASSERT_EQ(1u, vi);
    ASSERT_EQ(0u, n);       /* t1 did not poison the following point event */
    ASSERT_EQ(0u, cursor.nact);
    ASSERT_EQ(0, duckvep_sweep_cursor_next(&cursor, &vi, &slice, &n));
    PASS();
}

struct kprop_event_norm {
    duckvep_variant_batch_t v;
    uint16_t chrom;
    uint32_t pos;
    uint32_t end;
    uint8_t kind;
    uint8_t bytes[48];
    uint32_t roff;
    uint32_t aoff;
    uint16_t rlen;
    uint16_t alen;
};

static void kprop_event_norm_free(void *instance, void *env) {
    (void)env;
    free(instance);
}

static enum theft_alloc_res kprop_event_norm_alloc(struct theft *t, void *env, void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    struct kprop_event_norm *s = (struct kprop_event_norm *)calloc(1u, sizeof *s);
    uint32_t mode;
    uint32_t prefix;
    uint32_t suffix;
    uint32_t ref_mid;
    uint32_t alt_mid;
    uint32_t i;
    uint32_t o = 0u;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;

    mode = (uint32_t)kprop_bounded(t, 4u);
    prefix = (uint32_t)kprop_bounded(t, 4u);
    suffix = (uint32_t)kprop_bounded(t, 4u);
    if (mode == 0u) {          /* deletion */
        if (prefix == 0u) prefix = 1u;
        ref_mid = (uint32_t)kprop_bounded(t, 4u) + 1u;
        alt_mid = 0u;
        s->kind = (uint8_t)DUCKVEP_KIND_DEL;
    } else if (mode == 1u) {   /* insertion, including suffix-anchored prefix-zero */
        if (prefix == 0u && suffix == 0u) suffix = 1u;
        ref_mid = 0u;
        alt_mid = (uint32_t)kprop_bounded(t, 4u) + 1u;
        s->kind = (uint8_t)DUCKVEP_KIND_INS;
    } else if (mode == 2u) {   /* same-length substitution */
        ref_mid = (uint32_t)kprop_bounded(t, 3u) + 1u;
        alt_mid = ref_mid;
        s->kind = ref_mid == 1u && prefix == 0u && suffix == 0u
                ? (uint8_t)DUCKVEP_KIND_SNV
                : (uint8_t)DUCKVEP_KIND_MNV;
    } else {                   /* unequal replacement / delins-shaped */
        ref_mid = (uint32_t)kprop_bounded(t, 3u) + 1u;
        alt_mid = (uint32_t)kprop_bounded(t, 3u) + 1u;
        if (ref_mid == alt_mid) alt_mid++;
        s->kind = (uint8_t)DUCKVEP_KIND_INDEL;
    }
    if (prefix + ref_mid + suffix == 0u) prefix = 1u;
    if (prefix + alt_mid + suffix == 0u) prefix = 1u;
    if (prefix + ref_mid + suffix > 20u || prefix + alt_mid + suffix > 20u) {
        free(s);
        return THEFT_ALLOC_ERROR;
    }

    s->chrom = 0u;
    s->pos = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1u;
    s->roff = 0u;
    for (i = 0u; i < prefix; i++) s->bytes[o++] = (uint8_t)BASES[i % 4u];
    for (i = 0u; i < ref_mid; i++) s->bytes[o++] = (uint8_t)BASES[(i + 1u) % 4u];
    for (i = 0u; i < suffix; i++) s->bytes[o++] = (uint8_t)BASES[(i + 2u) % 4u];
    s->rlen = (uint16_t)o;
    s->aoff = o;
    for (i = 0u; i < prefix; i++) s->bytes[o++] = (uint8_t)BASES[i % 4u];
    for (i = 0u; i < alt_mid; i++) s->bytes[o++] = (uint8_t)BASES[(i + 3u) % 4u];
    for (i = 0u; i < suffix; i++) s->bytes[o++] = (uint8_t)BASES[(i + 2u) % 4u];
    s->alen = (uint16_t)(o - s->aoff);
    s->end = s->pos + (uint32_t)s->rlen - 1u;

    s->v.chrom_id = &s->chrom;
    s->v.pos1 = &s->pos;
    s->v.end1 = &s->end;
    s->v.variant_kind = &s->kind;
    s->v.allele_bytes = s->bytes;
    s->v.allele_bytes_len = o;
    s->v.ref_offset = &s->roff;
    s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff;
    s->v.alt_length = &s->alen;
    s->v.count = 1u;
    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_event_norm_info = {
    .alloc = kprop_event_norm_alloc,
    .free = kprop_event_norm_free,
};

static struct {
    uint32_t del;
    uint32_t ins;
    uint32_t sub;
    uint32_t indel;
    uint32_t prefix;
    uint32_t suffix;
    uint32_t interbase;
    uint32_t prefix_zero_interbase;
} g_event_norm_cov;

static enum theft_trial_res prop_event_normalization_matches_trim_oracle(
    struct theft *t,
    void         *arg1) {
    const struct kprop_event_norm *s = (const struct kprop_event_norm *)arg1;
    const uint8_t *ref = s->bytes + s->roff;
    const uint8_t *alt = s->bytes + s->aoff;
    uint16_t prefix = 0u;
    uint16_t suffix = 0u;
    uint16_t ref_rem;
    uint16_t alt_rem;
    uint16_t ref_diff_len;
    uint16_t alt_diff_len;
    uint32_t start;
    uint32_t end;
    duckvep_event_t event;
    (void)t;

    while (prefix < s->rlen && prefix < s->alen && ref[prefix] == alt[prefix]) prefix++;
    ref_rem = (uint16_t)(s->rlen - prefix);
    alt_rem = (uint16_t)(s->alen - prefix);
    while (suffix < ref_rem && suffix < alt_rem &&
           ref[(uint16_t)(s->rlen - 1u - suffix)] ==
           alt[(uint16_t)(s->alen - 1u - suffix)]) {
        suffix++;
    }
    ref_diff_len = (uint16_t)(s->rlen - prefix - suffix);
    alt_diff_len = (uint16_t)(s->alen - prefix - suffix);
    if (ref_diff_len == 0u && alt_diff_len > 0u) {
        start = prefix > 0u ? s->pos + (uint32_t)prefix - 1u : s->pos;
        end = start;
    } else {
        start = s->pos + (uint32_t)prefix;
        end = ref_diff_len > 0u ? start + (uint32_t)ref_diff_len - 1u : start;
    }

    duckvep_event_load(&s->v, 0u, &event);
    if (event.raw_start1 != s->pos || event.raw_end1 != s->end) return THEFT_TRIAL_FAIL;
    if (event.feature_start1 !=
            s->pos + (uint32_t)(s->rlen != s->alen && ref[0] == alt[0]) ||
        event.feature_end1 != s->end ||
        event.feature_allele_offset !=
            (uint16_t)(s->rlen != s->alen && ref[0] == alt[0])) {
        return THEFT_TRIAL_FAIL;
    }
    if (event.start1 != start || event.end1 != end) return THEFT_TRIAL_FAIL;
    if (event.ref_diff_offset != prefix || event.alt_diff_offset != prefix) return THEFT_TRIAL_FAIL;
    if (event.ref_diff_length != ref_diff_len ||
        event.alt_diff_length != alt_diff_len) return THEFT_TRIAL_FAIL;
    if (event.interbase != (uint8_t)(ref_diff_len == 0u && alt_diff_len > 0u)) {
        return THEFT_TRIAL_FAIL;
    }

    if (s->kind == (uint8_t)DUCKVEP_KIND_DEL) g_event_norm_cov.del++;
    else if (s->kind == (uint8_t)DUCKVEP_KIND_INS) g_event_norm_cov.ins++;
    else if (s->kind == (uint8_t)DUCKVEP_KIND_MNV ||
             s->kind == (uint8_t)DUCKVEP_KIND_SNV) g_event_norm_cov.sub++;
    else g_event_norm_cov.indel++;
    if (prefix > 0u) g_event_norm_cov.prefix++;
    if (suffix > 0u) g_event_norm_cov.suffix++;
    if (event.interbase) {
        g_event_norm_cov.interbase++;
        if (event.ref_diff_offset == 0u) g_event_norm_cov.prefix_zero_interbase++;
    }
    return THEFT_TRIAL_PASS;
}

TEST event_normalization_matches_trim_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "event differing-region normalization == independent trim oracle";
    cfg.prop1 = prop_event_normalization_matches_trim_oracle;
    cfg.type_info[0] = &kprop_event_norm_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_event_norm_cov, 0, sizeof g_event_norm_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_event_norm_cov.del > 0u);
    ASSERT(g_event_norm_cov.ins > 0u);
    ASSERT(g_event_norm_cov.sub > 0u);
    ASSERT(g_event_norm_cov.indel > 0u);
    ASSERT(g_event_norm_cov.prefix > 0u);
    ASSERT(g_event_norm_cov.suffix > 0u);
    ASSERT(g_event_norm_cov.interbase > 0u);
    ASSERT(g_event_norm_cov.prefix_zero_interbase > 0u);
    fprintf(stderr,
            "[event normalization coverage] del=%u ins=%u sub=%u indel=%u "
            "prefix=%u suffix=%u interbase=%u prefix0_interbase=%u\n",
            g_event_norm_cov.del, g_event_norm_cov.ins, g_event_norm_cov.sub,
            g_event_norm_cov.indel, g_event_norm_cov.prefix,
            g_event_norm_cov.suffix, g_event_norm_cov.interbase,
            g_event_norm_cov.prefix_zero_interbase);
    PASS();
}

struct allele_sweep_vrec {
    uint16_t chrom;
    uint32_t pos;
    uint32_t end;
    uint8_t kind;
    uint8_t ref[24];
    uint8_t alt[24];
    uint16_t rlen;
    uint16_t alen;
};

struct kprop_allele_sweep_scene {
    duckvep_variant_batch_t v;
    duckvep_transcript_model_t tx;
    uint32_t halo;
    uint16_t *vchrom;
    uint32_t *vpos;
    uint32_t *vend;
    uint8_t  *vkind;
    uint32_t *vroff;
    uint16_t *vrlen;
    uint32_t *vaoff;
    uint16_t *valen;
    uint8_t  *vbytes;
    size_t    vbytes_len;
    uint16_t *tchrom;
    uint32_t *tstart;
    uint32_t *tend;
    int8_t   *tstrand;
    uint64_t *tflags;
    uint32_t *texoff;
    uint16_t *texcnt;
    uint32_t *tcds_s;
    uint32_t *tcds_e;
};

static int allele_sweep_vrec_cmp(const void *a, const void *b) {
    const struct allele_sweep_vrec *x = (const struct allele_sweep_vrec *)a;
    const struct allele_sweep_vrec *y = (const struct allele_sweep_vrec *)b;
    if (x->chrom != y->chrom) return x->chrom < y->chrom ? -1 : 1;
    if (x->pos != y->pos) return x->pos < y->pos ? -1 : 1;
    return 0;
}

static uint8_t allele_sweep_base(uint32_t i) {
    static const uint8_t BASES[4] = {'A', 'C', 'G', 'T'};
    return BASES[i % 4u];
}

static void allele_sweep_fill_variant(struct theft *t, struct allele_sweep_vrec *r) {
    uint32_t mode = (uint32_t)kprop_bounded(t, 4u);
    uint32_t prefix = (uint32_t)kprop_bounded(t, 4u);
    uint32_t suffix = (uint32_t)kprop_bounded(t, 3u);
    uint32_t ref_mid;
    uint32_t alt_mid;
    uint32_t i;
    uint32_t nr = 0u;
    uint32_t na = 0u;

    if (mode == 0u) {
        if (prefix == 0u) prefix = 1u;
        ref_mid = (uint32_t)kprop_bounded(t, 4u) + 1u;
        alt_mid = 0u;
        r->kind = (uint8_t)DUCKVEP_KIND_DEL;
    } else if (mode == 1u) {
        if (prefix == 0u) prefix = 1u;
        ref_mid = 0u;
        alt_mid = (uint32_t)kprop_bounded(t, 4u) + 1u;
        r->kind = (uint8_t)DUCKVEP_KIND_INS;
    } else if (mode == 2u) {
        ref_mid = (uint32_t)kprop_bounded(t, 3u) + 1u;
        alt_mid = ref_mid;
        if (prefix + ref_mid + suffix < 2u) suffix = 1u;
        r->kind = (uint8_t)DUCKVEP_KIND_MNV;
    } else {
        ref_mid = (uint32_t)kprop_bounded(t, 3u) + 1u;
        alt_mid = (uint32_t)kprop_bounded(t, 3u) + 1u;
        if (ref_mid == alt_mid) alt_mid++;
        r->kind = (uint8_t)DUCKVEP_KIND_INDEL;
    }

    for (i = 0u; i < prefix; i++) {
        r->ref[nr++] = allele_sweep_base(i);
        r->alt[na++] = allele_sweep_base(i);
    }
    for (i = 0u; i < ref_mid; i++) r->ref[nr++] = allele_sweep_base(i + 1u);
    for (i = 0u; i < alt_mid; i++) r->alt[na++] = allele_sweep_base(i + 3u);
    for (i = 0u; i < suffix; i++) {
        r->ref[nr++] = allele_sweep_base(i + 2u);
        r->alt[na++] = allele_sweep_base(i + 2u);
    }
    r->rlen = (uint16_t)nr;
    r->alen = (uint16_t)na;
    r->end = r->pos + (uint32_t)r->rlen - 1u;
}

static void kprop_allele_sweep_scene_free(void *instance, void *env) {
    struct kprop_allele_sweep_scene *s = (struct kprop_allele_sweep_scene *)instance;
    (void)env;
    if (s == NULL) return;
    free(s->vchrom); free(s->vpos); free(s->vend); free(s->vkind);
    free(s->vroff); free(s->vrlen); free(s->vaoff); free(s->valen); free(s->vbytes);
    free(s->tchrom); free(s->tstart); free(s->tend); free(s->tstrand); free(s->tflags);
    free(s->texoff); free(s->texcnt); free(s->tcds_s); free(s->tcds_e);
    free(s);
}

static enum theft_alloc_res kprop_allele_sweep_scene_alloc(
    struct theft *t,
    void         *env,
    void        **instance) {
    static const uint32_t halos[4] = {0u, 1u, 50u, 100u};
    struct kprop_allele_sweep_scene *s;
    struct allele_sweep_vrec vr[KPROP_MAX_VARIANTS];
    struct trec tr[KPROP_MAX_TX];
    size_t nv = (size_t)kprop_bounded(t, (uint64_t)KPROP_MAX_VARIANTS) + 1u;
    size_t ntx = (size_t)kprop_bounded(t, (uint64_t)KPROP_MAX_TX) + 1u;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFFF0000u) + 1u;
    uint32_t span = (uint32_t)kprop_bounded(t, 5000u) + 64u;
    size_t i;
    size_t off = 0u;
    (void)env;

    s = (struct kprop_allele_sweep_scene *)calloc(1u, sizeof *s);
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->halo = halos[kprop_bounded(t, 4u)];

    for (i = 0u; i < nv; i++) {
        vr[i].chrom = (uint16_t)kprop_bounded(t, KPROP_NCHROM);
        vr[i].pos = base + (uint32_t)kprop_bounded(t, span);
        allele_sweep_fill_variant(t, &vr[i]);
    }
    for (i = 0u; i < ntx; i++) {
        tr[i].chrom = (uint16_t)kprop_bounded(t, KPROP_NCHROM);
        tr[i].start = base + (uint32_t)kprop_bounded(t, span);
        tr[i].len = (uint32_t)kprop_bounded(t, span);
        tr[i].strand = kprop_bounded(t, 2u) == 0u ? (int8_t)1 : (int8_t)-1;
    }
    qsort(vr, nv, sizeof vr[0], allele_sweep_vrec_cmp);
    qsort(tr, ntx, sizeof tr[0], trec_cmp);

    s->vchrom = (uint16_t *)calloc(nv, sizeof *s->vchrom);
    s->vpos   = (uint32_t *)calloc(nv, sizeof *s->vpos);
    s->vend   = (uint32_t *)calloc(nv, sizeof *s->vend);
    s->vkind  = (uint8_t  *)calloc(nv, sizeof *s->vkind);
    s->vroff  = (uint32_t *)calloc(nv, sizeof *s->vroff);
    s->vrlen  = (uint16_t *)calloc(nv, sizeof *s->vrlen);
    s->vaoff  = (uint32_t *)calloc(nv, sizeof *s->vaoff);
    s->valen  = (uint16_t *)calloc(nv, sizeof *s->valen);
    s->vbytes = (uint8_t  *)calloc(nv * 48u, 1u);
    s->tchrom = (uint16_t *)calloc(ntx, sizeof *s->tchrom);
    s->tstart = (uint32_t *)calloc(ntx, sizeof *s->tstart);
    s->tend   = (uint32_t *)calloc(ntx, sizeof *s->tend);
    s->tstrand= (int8_t   *)calloc(ntx, sizeof *s->tstrand);
    s->tflags = (uint64_t *)calloc(ntx, sizeof *s->tflags);
    s->texoff = (uint32_t *)calloc(ntx, sizeof *s->texoff);
    s->texcnt = (uint16_t *)calloc(ntx, sizeof *s->texcnt);
    s->tcds_s = (uint32_t *)calloc(ntx, sizeof *s->tcds_s);
    s->tcds_e = (uint32_t *)calloc(ntx, sizeof *s->tcds_e);
    if (s->vchrom == NULL || s->vpos == NULL || s->vend == NULL ||
        s->vkind == NULL || s->vroff == NULL || s->vrlen == NULL ||
        s->vaoff == NULL || s->valen == NULL || s->vbytes == NULL ||
        s->tchrom == NULL || s->tstart == NULL || s->tend == NULL ||
        s->tstrand == NULL || s->tflags == NULL || s->texoff == NULL ||
        s->texcnt == NULL || s->tcds_s == NULL || s->tcds_e == NULL) {
        kprop_allele_sweep_scene_free(s, NULL);
        return THEFT_ALLOC_ERROR;
    }

    for (i = 0u; i < nv; i++) {
        s->vchrom[i] = vr[i].chrom;
        s->vpos[i] = vr[i].pos;
        s->vend[i] = vr[i].end;
        s->vkind[i] = vr[i].kind;
        s->vroff[i] = (uint32_t)off;
        memcpy(s->vbytes + off, vr[i].ref, vr[i].rlen);
        off += vr[i].rlen;
        s->vrlen[i] = vr[i].rlen;
        s->vaoff[i] = (uint32_t)off;
        memcpy(s->vbytes + off, vr[i].alt, vr[i].alen);
        off += vr[i].alen;
        s->valen[i] = vr[i].alen;
    }
    s->vbytes_len = off;
    for (i = 0u; i < ntx; i++) {
        s->tchrom[i] = tr[i].chrom;
        s->tstart[i] = tr[i].start;
        s->tend[i] = tr[i].start + tr[i].len;
        s->tstrand[i] = tr[i].strand;
    }

    s->v.chrom_id = s->vchrom; s->v.pos1 = s->vpos; s->v.end1 = s->vend;
    s->v.ref_offset = s->vroff; s->v.ref_length = s->vrlen;
    s->v.alt_offset = s->vaoff; s->v.alt_length = s->valen;
    s->v.allele_bytes = s->vbytes; s->v.allele_bytes_len = s->vbytes_len;
    s->v.variant_kind = s->vkind; s->v.count = nv;
    s->tx.chrom_id = s->tchrom; s->tx.start1 = s->tstart; s->tx.end1 = s->tend;
    s->tx.strand = s->tstrand; s->tx.flags = s->tflags;
    s->tx.exon_offset = s->texoff; s->tx.exon_count = s->texcnt;
    s->tx.cds_start1 = s->tcds_s; s->tx.cds_end1 = s->tcds_e;
    s->tx.transcript_count = ntx;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_allele_sweep_scene_info = {
    .alloc = kprop_allele_sweep_scene_alloc,
    .free = kprop_allele_sweep_scene_free,
};

static uint32_t allele_sweep_oracle_end(const duckvep_variant_batch_t *v, size_t i) {
    const uint8_t *ref = v->allele_bytes + v->ref_offset[i];
    const uint8_t *alt = v->allele_bytes + v->alt_offset[i];
    uint16_t ref_len = v->ref_length[i];
    uint16_t alt_len = v->alt_length[i];
    uint32_t feature_start = v->pos1[i] +
        (uint32_t)(ref_len != alt_len && ref[0] == alt[0]);
    uint32_t feature_end = v->pos1[i] + (uint32_t)ref_len - 1u;

    return feature_start > feature_end ? feature_start : feature_end;
}

static size_t brute_allele_sweep_candidates(const struct kprop_allele_sweep_scene *s,
                                            uint64_t *out,
                                            size_t cap) {
    size_t n = 0u;
    size_t vi;
    size_t ti;
    for (vi = 0u; vi < s->v.count; vi++) {
        uint64_t start = (uint64_t)s->v.pos1[vi];
        uint64_t end = (uint64_t)allele_sweep_oracle_end(&s->v, vi);
        for (ti = 0u; ti < s->tx.transcript_count; ti++) {
            if (s->v.chrom_id[vi] != s->tx.chrom_id[ti]) continue;
            if ((uint64_t)s->tx.start1[ti] <= end + s->halo &&
                (uint64_t)s->tx.end1[ti] + s->halo >= start) {
                if (n < cap) out[n] = ((uint64_t)vi << 32) | (uint64_t)ti;
                n++;
            }
        }
    }
    return n;
}

static struct {
    uint32_t del;
    uint32_t ins;
    uint32_t mnv;
    uint32_t indel;
    uint32_t prefix;
    uint32_t suffix;
    uint32_t interbase;
    uint32_t tail;
} g_allele_sweep_cov;

static enum theft_trial_res prop_allele_sweep_matches_trim_oracle(
    struct theft *t,
    void         *arg1) {
    const struct kprop_allele_sweep_scene *s;
    uint64_t sweep_buf[KPROP_MAX_PAIRS];
    uint64_t brute_buf[KPROP_MAX_PAIRS];
    uint32_t active[KPROP_MAX_TX];
    uint32_t candidates[KPROP_MAX_TX];
    struct pair_collector col;
    duckvep_status_t st = DUCKVEP_OK;
    size_t n_sweep;
    size_t n_brute;
    size_t i;
    (void)t;

    s = (const struct kprop_allele_sweep_scene *)arg1;
    col.buf = sweep_buf; col.n = 0u; col.cap = KPROP_MAX_PAIRS;
    n_sweep = duckvep_sweep_candidates(&s->v, &s->tx, s->halo,
                                       active, KPROP_MAX_TX,
                                       candidates, KPROP_MAX_TX,
                                       pair_sink, &col, &st);
    if (st != DUCKVEP_OK || n_sweep != col.n) return THEFT_TRIAL_FAIL;
    n_brute = brute_allele_sweep_candidates(s, brute_buf, KPROP_MAX_PAIRS);
    if (n_sweep != n_brute) return THEFT_TRIAL_FAIL;

    qsort(sweep_buf, n_sweep, sizeof sweep_buf[0], u64_cmp);
    qsort(brute_buf, n_brute, sizeof brute_buf[0], u64_cmp);
    for (i = 0u; i < n_sweep; i++) {
        if (sweep_buf[i] != brute_buf[i]) return THEFT_TRIAL_FAIL;
    }
    for (i = 0u; i < s->v.count; i++) {
        duckvep_event_t event;
        duckvep_event_load(&s->v, i, &event);
        if (duckvep_event_feature_max1(&event) !=
            allele_sweep_oracle_end(&s->v, i)) return THEFT_TRIAL_FAIL;
        if (s->vkind[i] == (uint8_t)DUCKVEP_KIND_DEL) g_allele_sweep_cov.del++;
        else if (s->vkind[i] == (uint8_t)DUCKVEP_KIND_INS) g_allele_sweep_cov.ins++;
        else if (s->vkind[i] == (uint8_t)DUCKVEP_KIND_MNV) g_allele_sweep_cov.mnv++;
        else if (s->vkind[i] == (uint8_t)DUCKVEP_KIND_INDEL) g_allele_sweep_cov.indel++;
        if (event.ref_diff_offset > 0u) g_allele_sweep_cov.prefix++;
        if ((uint16_t)(s->vrlen[i] - event.ref_diff_offset - event.ref_diff_length) > 0u) {
            g_allele_sweep_cov.suffix++;
        }
        if (event.interbase) g_allele_sweep_cov.interbase++;
        if (event.end1 > event.raw_start1) g_allele_sweep_cov.tail++;
    }
    return THEFT_TRIAL_PASS;
}

TEST sweep_vep_feature_span_candidates_match_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "VEP feature-span sweep candidates == independent parser oracle";
    cfg.prop1 = prop_allele_sweep_matches_trim_oracle;
    cfg.type_info[0] = &kprop_allele_sweep_scene_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_allele_sweep_cov, 0, sizeof g_allele_sweep_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_allele_sweep_cov.del > 0u);
    ASSERT(g_allele_sweep_cov.ins > 0u);
    ASSERT(g_allele_sweep_cov.mnv > 0u);
    ASSERT(g_allele_sweep_cov.indel > 0u);
    ASSERT(g_allele_sweep_cov.prefix > 0u);
    ASSERT(g_allele_sweep_cov.suffix > 0u);
    ASSERT(g_allele_sweep_cov.interbase > 0u);
    ASSERT(g_allele_sweep_cov.tail > 0u);
    fprintf(stderr,
            "[allele sweep coverage] del=%u ins=%u mnv=%u indel=%u "
            "prefix=%u suffix=%u interbase=%u tail=%u\n",
            g_allele_sweep_cov.del, g_allele_sweep_cov.ins,
            g_allele_sweep_cov.mnv, g_allele_sweep_cov.indel,
            g_allele_sweep_cov.prefix, g_allele_sweep_cov.suffix,
            g_allele_sweep_cov.interbase, g_allele_sweep_cov.tail);
    PASS();
}

TEST sweep_known_scene_exact_pairs(void) {
    static const uint16_t vchrom[3] = {0u, 0u, 1u};
    static const uint32_t vpos[3]   = {150u, 650u, 150u};
    static const uint32_t vend[3]   = {150u, 650u, 150u};
    static const uint16_t tchrom[3] = {0u, 0u, 1u};
    static const uint32_t tstart[3] = {100u, 500u, 100u};
    static const uint32_t tend[3]   = {200u, 600u, 200u};
    duckvep_variant_batch_t v;
    duckvep_transcript_model_t tx;
    uint64_t buf[8];
    uint32_t active[3];
    uint32_t candidates[3];
    struct pair_collector col;
    duckvep_status_t st = DUCKVEP_OK;
    size_t n;

    memset(&v, 0, sizeof v);
    memset(&tx, 0, sizeof tx);
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.count = 3u;
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.transcript_count = 3u;

    col.buf = buf; col.n = 0u; col.cap = 8u;
    n = duckvep_sweep_candidates(&v, &tx, 100u, active, 3u,
                                 candidates, 3u, pair_sink, &col, &st);
    ASSERT_EQ(DUCKVEP_OK, st);
    ASSERT_EQ(3u, n);
    qsort(buf, n, sizeof buf[0], u64_cmp);
    ASSERT_EQ((((uint64_t)0u << 32) | 0u), buf[0]); /* v0,t0 */
    ASSERT_EQ((((uint64_t)1u << 32) | 1u), buf[1]); /* v1,t1 */
    ASSERT_EQ((((uint64_t)2u << 32) | 2u), buf[2]); /* v2,t2 */
    PASS();
}

/* Event ends are not monotone when events are sorted by start. A wide first
 * event must see the distant transcript, but that span-only tail must not enter
 * the persistent point active set or poison the following short event. */
TEST sweep_span_tail_does_not_poison_point_frontier(void) {
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vstart[2] = {100u, 101u};
    static const uint32_t vend[2] = {1000u, 101u};
    static const uint16_t tchrom[2] = {0u, 0u};
    static const uint32_t tstart[2] = {100u, 900u};
    static const uint32_t tend[2] = {110u, 910u};
    duckvep_variant_batch_t v;
    duckvep_transcript_model_t tx;
    duckvep_sweep_cursor_t cursor;
    uint32_t active[2];
    uint32_t candidates[2];
    uint32_t vi;
    const uint32_t *slice;
    size_t n;

    memset(&v, 0, sizeof v);
    memset(&tx, 0, sizeof tx);
    v.chrom_id = vchrom; v.pos1 = vstart; v.end1 = vend; v.count = 2u;
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.transcript_count = 2u;

    duckvep_sweep_cursor_init(&cursor, &v, &tx, 0u, active, 2u,
                              candidates, 2u);
    ASSERT_EQ(1, duckvep_sweep_cursor_next(&cursor, &vi, &slice, &n));
    ASSERT_EQ(0u, vi);
    ASSERT_EQ(2u, n);
    ASSERT_EQ(1u, cursor.nact); /* distant t1 was a borrowed span tail */
    ASSERT_EQ(0u, slice[0]);
    ASSERT_EQ(1u, slice[1]);

    ASSERT_EQ(1, duckvep_sweep_cursor_next(&cursor, &vi, &slice, &n));
    ASSERT_EQ(1u, vi);
    ASSERT_EQ(1u, n);
    ASSERT_EQ(1u, cursor.nact);
    ASSERT_EQ(0u, slice[0]);
    ASSERT_EQ(0, duckvep_sweep_cursor_next(&cursor, &vi, &slice, &n));
    PASS();
}


/* Raw VCF REF spans are not the small-variant semantic interval. Until the
 * differing-region normalizer lands, only the structural lane may use end1 for
 * overlap. The same supplied [100,1000] geometry is a point for DEL and a span
 * for SV, preventing shared VCF padding from creating boundary over-calls. */
TEST sweep_uses_full_span_only_for_structural_events(void) {
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vstart[2] = {100u, 100u};
    static const uint32_t vend[2] = {1000u, 1000u};
    static const uint8_t vkind[2] = {DUCKVEP_KIND_DEL, DUCKVEP_KIND_SV};
    static const uint16_t tchrom[2] = {0u, 0u};
    static const uint32_t tstart[2] = {100u, 900u};
    static const uint32_t tend[2] = {110u, 910u};
    duckvep_variant_batch_t v;
    duckvep_transcript_model_t tx;
    duckvep_sweep_cursor_t cursor;
    uint32_t active[2];
    uint32_t candidates[2];
    uint32_t vi;
    const uint32_t *slice;
    size_t n;

    memset(&v, 0, sizeof v);
    memset(&tx, 0, sizeof tx);
    v.chrom_id = vchrom; v.pos1 = vstart; v.end1 = vend;
    v.variant_kind = vkind; v.count = 2u;
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.transcript_count = 2u;

    duckvep_sweep_cursor_init(&cursor, &v, &tx, 0u, active, 2u,
                              candidates, 2u);
    ASSERT_EQ(1, duckvep_sweep_cursor_next(&cursor, &vi, &slice, &n));
    ASSERT_EQ(0u, vi);
    ASSERT_EQ(1u, n);
    ASSERT_EQ(0u, slice[0]);

    ASSERT_EQ(1, duckvep_sweep_cursor_next(&cursor, &vi, &slice, &n));
    ASSERT_EQ(1u, vi);
    ASSERT_EQ(2u, n);
    ASSERT_EQ(0u, slice[0]);
    ASSERT_EQ(1u, slice[1]);
    ASSERT_EQ(0, duckvep_sweep_cursor_next(&cursor, &vi, &slice, &n));
    PASS();
}

struct stop_sink_ctx { size_t calls; };
static int stop_immediately_sink(uint32_t vi, uint32_t ti, void *ctx) {
    struct stop_sink_ctx *s = (struct stop_sink_ctx *)ctx;
    (void)vi; (void)ti;
    s->calls++;
    return 0;
}

/* A sink cancellation is control flow, not a hint: once it declines a pair the
 * wrapper must not drain the rest of the candidate relation. */
TEST sweep_sink_stop_is_immediate(void) {
    static const uint16_t vchrom[3] = {0u, 0u, 0u};
    static const uint32_t vpos[3]   = {150u, 160u, 170u};
    static const uint32_t vend[3]   = {150u, 160u, 170u};
    static const uint16_t tchrom[2] = {0u, 0u};
    static const uint32_t tstart[2] = {100u, 100u};
    static const uint32_t tend[2]   = {200u, 200u};
    duckvep_variant_batch_t v;
    duckvep_transcript_model_t tx;
    uint32_t active[2];
    uint32_t candidates[2];
    struct stop_sink_ctx stop;
    duckvep_status_t st = DUCKVEP_OK;
    size_t n;

    memset(&v, 0, sizeof v); memset(&tx, 0, sizeof tx); memset(&stop, 0, sizeof stop);
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.count = 3u;
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.transcript_count = 2u;

    n = duckvep_sweep_candidates(&v, &tx, 0u, active, 2u,
                                 candidates, 2u,
                                 stop_immediately_sink, &stop, &st);
    ASSERT_EQ(DUCKVEP_OK, st);
    ASSERT_EQ(0u, n);          /* the declined pair was not accepted */
    ASSERT_EQ(1u, stop.calls); /* the remaining five candidates were not visited */
    PASS();
}

TEST sweep_matches_bruteforce_for_any_scene(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "sweep candidate set == brute-force candidate set";
    cfg.prop1 = prop_sweep_matches_bruteforce;
    cfg.type_info[0] = &kprop_scene_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

/* Saturation edge: a transcript and variant near UINT32_MAX so that pos+halo and
 * end1+halo both overflow uint32. The saturating add must keep the near-max
 * transcript a candidate and still evict the far one — no wrap-around. */
TEST sweep_saturation_near_uint32_max(void) {
    const uint32_t MAXV = 0xFFFFFFFFu;
    static uint16_t tchrom[2] = {0u, 0u};
    uint32_t tstart[2]; uint32_t tend[2];
    uint16_t vchrom[1] = {0u};
    uint32_t vpos[1];
    uint32_t vend[1];
    duckvep_variant_batch_t v;
    duckvep_transcript_model_t tx;
    uint64_t buf[4];
    uint32_t active[2];
    uint32_t candidates[2];
    struct pair_collector col;
    duckvep_status_t st = DUCKVEP_OK;
    size_t n;

    tstart[0] = 10u;          tend[0] = 20u;            /* far: must be evicted   */
    tstart[1] = MAXV - 50u;   tend[1] = MAXV - 10u;     /* near max: candidate    */
    vpos[0]   = MAXV - 5u;
    vend[0]   = vpos[0];

    memset(&v, 0, sizeof v); memset(&tx, 0, sizeof tx);
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.count = 1u;
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.transcript_count = 2u;

    col.buf = buf; col.n = 0u; col.cap = 4u;
    n = duckvep_sweep_candidates(&v, &tx, 100u, active, 2u,
                                 candidates, 2u, pair_sink, &col, &st);
    ASSERT_EQ(DUCKVEP_OK, st);
    ASSERT_EQ(1u, n);
    ASSERT_EQ((((uint64_t)0u << 32) | 1u), buf[0]); /* v0 -> t1 (near max) */
    PASS();
}

/* active_cap overflow: more simultaneously-active transcripts than scratch can
 * hold must stop with DUCKVEP_ERR_RESULT_FULL, never silently truncate. */
TEST sweep_active_cap_overflow_is_reported(void) {
    static uint16_t tchrom[4] = {0u, 0u, 0u, 0u};
    static uint32_t tstart[4] = {100u, 101u, 102u, 103u};
    static uint32_t tend[4]   = {200u, 200u, 200u, 200u};
    uint16_t vchrom[1] = {0u};
    uint32_t vpos[1]   = {150u};
    uint32_t vend[1]   = {150u};
    duckvep_variant_batch_t v;
    duckvep_transcript_model_t tx;
    uint32_t active[2]; /* capacity 2, but 4 transcripts are simultaneously active */
    uint32_t candidates[2];
    duckvep_status_t st = DUCKVEP_OK;

    memset(&v, 0, sizeof v); memset(&tx, 0, sizeof tx);
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.count = 1u;
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.transcript_count = 4u;

    (void)duckvep_sweep_candidates(&v, &tx, 10u, active, 2u,
                                    candidates, 2u, NULL, NULL, &st);
    ASSERT_EQ(DUCKVEP_ERR_RESULT_FULL, st);
    PASS();
}

/* ===================================================================== *
 * Structural region classification (duckvep_region_mask).
 *
 * Validated two ways: hand-computed deterministic scenes that pin the SPEC, and
 * a property checking genuine structural invariants (not a restatement of the
 * implementation): exactly one primary region bit; outside the span only up/down
 * (and never splice); CDS => coding & in-range; EXON => non-coding; SPLICE iff a
 * boundary is within the splice distance.
 * ===================================================================== */

#define KPROP_MAX_EXONS 6u

static int popcount_u32(uint32_t x) {
    int c = 0;
    while (x != 0u) { x &= (uint32_t)(x - 1u); c++; }
    return c;
}

struct kprop_tx1 {
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t       ex;
    uint16_t chrom; uint32_t tstart; uint32_t tend; int8_t strand; uint64_t flags;
    uint32_t exoff; uint16_t excnt; uint32_t cds_s; uint32_t cds_e;
    uint32_t *es; uint32_t *ee; uint32_t *ecds; uint32_t *ecde; int8_t *eph; int8_t *eeph;
    uint32_t pos; uint32_t splice_exonic; uint32_t splice_intronic;
};

static void kprop_tx1_free(void *instance, void *env) {
    struct kprop_tx1 *s = (struct kprop_tx1 *)instance;
    (void)env;
    if (s == NULL) return;
    free(s->es); free(s->ee); free(s->ecds); free(s->ecde); free(s->eph); free(s->eeph);
    free(s);
}

static enum theft_alloc_res kprop_tx1_alloc(struct theft *t, void *env, void **instance) {
    struct kprop_tx1 *s = (struct kprop_tx1 *)calloc(1u, sizeof *s);
    uint32_t nex = (uint32_t)kprop_bounded(t, KPROP_MAX_EXONS) + 1u;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t span = (uint32_t)kprop_bounded(t, 20000u) + nex * 8u;
    uint32_t seg, i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;

    s->es = (uint32_t *)calloc(nex, sizeof *s->es);
    s->ee = (uint32_t *)calloc(nex, sizeof *s->ee);
    s->ecds = (uint32_t *)calloc(nex, sizeof *s->ecds);
    s->ecde = (uint32_t *)calloc(nex, sizeof *s->ecde);
    s->eph = (int8_t *)calloc(nex, sizeof *s->eph);
    s->eeph = (int8_t *)calloc(nex, sizeof *s->eeph);
    if (!s->es || !s->ee || !s->ecds || !s->ecde || !s->eph || !s->eeph) {
        kprop_tx1_free(s, NULL); return THEFT_ALLOC_ERROR;
    }

    s->chrom = 0u;
    s->tstart = base;
    s->tend = base + span;
    s->strand = (kprop_bounded(t, 2u) == 0u) ? (int8_t)1 : (int8_t)-1;
    seg = span / nex;
    for (i = 0u; i < nex; i++) {
        uint32_t lo = base + i * seg;
        uint32_t a = lo + (uint32_t)kprop_bounded(t, seg / 2u + 1u);
        uint32_t b = a + (uint32_t)kprop_bounded(t, seg / 2u + 1u);
        if (b > s->tend) b = s->tend;
        s->es[i] = a; s->ee[i] = b;
    }
    if (kprop_bounded(t, 2u) == 0u) { /* coding */
        s->cds_s = base + (uint32_t)kprop_bounded(t, span / 2u + 1u) + 1u;
        s->cds_e = s->cds_s + (uint32_t)kprop_bounded(t, span);
        if (s->cds_e > s->tend) s->cds_e = s->tend;
    } else {
        s->cds_s = 0u; s->cds_e = 0u;
    }
    s->splice_exonic = (uint32_t)kprop_bounded(t, 6u);    /* 0..5, incl. 0 = disabled */
    s->splice_intronic = (uint32_t)kprop_bounded(t, 12u); /* 0..11 */
    s->pos = (s->tstart > 100u ? s->tstart - 100u : 0u)
           + (uint32_t)kprop_bounded(t, span + 200u);
    s->exoff = 0u; s->excnt = (uint16_t)nex;

    s->ex.start1 = s->es; s->ex.end1 = s->ee;
    s->ex.cdna_start1 = s->ecds; s->ex.cdna_end1 = s->ecde;
    s->ex.phase = s->eph; s->ex.end_phase = s->eeph; s->ex.exon_count = nex;

    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e;
    s->tx.transcript_count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_tx1_info = {
    .alloc = kprop_tx1_alloc,
    .free  = kprop_tx1_free,
};

static enum theft_trial_res prop_region_mask_invariants(struct theft *t, void *arg1) {
    const struct kprop_tx1 *s = (const struct kprop_tx1 *)arg1;
    uint32_t mask = duckvep_region_mask(&s->tx, &s->ex, 0u, s->pos,
                                        s->splice_exonic, s->splice_intronic);
    uint32_t primary = mask & ~(uint32_t)DUCKVEP_REGION_SPLICE;
    int in_exon = 0, splice_ref = 0;
    size_t e;
    (void)t;

    if (popcount_u32(primary) != 1) return THEFT_TRIAL_FAIL;

    if (s->pos < s->tstart || s->pos > s->tend) {
        uint32_t want = (s->pos < s->tstart)
            ? (s->strand >= 0 ? (uint32_t)DUCKVEP_REGION_UPSTREAM : (uint32_t)DUCKVEP_REGION_DOWNSTREAM)
            : (s->strand >= 0 ? (uint32_t)DUCKVEP_REGION_DOWNSTREAM : (uint32_t)DUCKVEP_REGION_UPSTREAM);
        if (mask != want) return THEFT_TRIAL_FAIL; /* exactly up/down, no splice */
        return THEFT_TRIAL_PASS;
    }

    for (e = 0; e < s->ex.exon_count; e++) {
        uint32_t es = s->es[e], ee = s->ee[e];
        /* Independent per-side splice oracle: exonic reach INSIDE the exon, intronic
         * reach OUTSIDE it (mirrors VEP's asymmetric splice_region definition). */
        if (s->pos >= es && s->pos <= ee) {
            in_exon = 1;
            if (s->pos - es < s->splice_exonic || ee - s->pos < s->splice_exonic) splice_ref = 1;
        } else if (s->pos < es) {
            if (es - s->pos <= s->splice_intronic) splice_ref = 1;
        } else { /* s->pos > ee */
            if (s->pos - ee <= s->splice_intronic) splice_ref = 1;
        }
    }
    if (primary == (uint32_t)DUCKVEP_REGION_INTRON && in_exon) return THEFT_TRIAL_FAIL;
    if (primary != (uint32_t)DUCKVEP_REGION_INTRON && !in_exon) return THEFT_TRIAL_FAIL;
    if (primary == (uint32_t)DUCKVEP_REGION_CDS) {
        if (!(s->cds_s != 0u && s->pos >= s->cds_s && s->pos <= s->cds_e)) return THEFT_TRIAL_FAIL;
    }
    if (primary == (uint32_t)DUCKVEP_REGION_UTR) {
        if (s->cds_s == 0u) return THEFT_TRIAL_FAIL;
        if (s->pos >= s->cds_s && s->pos <= s->cds_e) return THEFT_TRIAL_FAIL;
    }
    if (primary == (uint32_t)DUCKVEP_REGION_EXON && s->cds_s != 0u) return THEFT_TRIAL_FAIL;
    if (((mask & (uint32_t)DUCKVEP_REGION_SPLICE) != 0) != (splice_ref != 0)) return THEFT_TRIAL_FAIL;
    return THEFT_TRIAL_PASS;
}

TEST region_mask_invariants_hold(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "region mask structural invariants";
    cfg.prop1 = prop_region_mask_invariants;
    cfg.type_info[0] = &kprop_tx1_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

/* The production SNV path carries one exon cursor per transcript across sorted
 * tiles. Keep the old all-exon classifiers as an independent oracle: randomized
 * transcript shapes are walked base-by-base, including both strands, short
 * exons/introns, outer transcript padding, CDS/UTR boundaries, and every splice
 * window. */
#define KPROP_POINT_MAX_EXONS 8u

struct kprop_point_scene {
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    uint16_t chrom;
    uint32_t tstart;
    uint32_t tend;
    int8_t strand;
    uint64_t flags;
    uint32_t exoff;
    uint16_t excnt;
    uint32_t cds_s;
    uint32_t cds_e;
    uint32_t es[KPROP_POINT_MAX_EXONS];
    uint32_t ee[KPROP_POINT_MAX_EXONS];
    uint32_t splice_exonic;
    uint32_t splice_intronic;
};

static enum theft_alloc_res kprop_point_scene_alloc(
    struct theft *t, void *env, void **instance) {

    struct kprop_point_scene *s;
    uint32_t asc_s[KPROP_POINT_MAX_EXONS];
    uint32_t asc_e[KPROP_POINT_MAX_EXONS];
    uint32_t pos;
    uint32_t i;
    uint32_t cnt;
    uint32_t left_pad;
    uint32_t right_pad;
    (void)env;

    s = (struct kprop_point_scene *)calloc(1u, sizeof *s);
    if (s == NULL) return THEFT_ALLOC_ERROR;
    cnt = 1u + (uint32_t)kprop_bounded(t, KPROP_POINT_MAX_EXONS);
    pos = 1000u + (uint32_t)kprop_bounded(t, 1000000u);
    for (i = 0u; i < cnt; i++) {
        uint32_t len = 1u + (uint32_t)kprop_bounded(t, 48u);
        uint32_t gap = 1u + (uint32_t)kprop_bounded(t, 48u);
        asc_s[i] = pos;
        asc_e[i] = pos + len - 1u;
        pos = asc_e[i] + gap + 1u;
    }
    s->strand = kprop_bounded(t, 2u) == 0u ? (int8_t)1 : (int8_t)-1;
    for (i = 0u; i < cnt; i++) {
        uint32_t src = s->strand >= 0 ? i : cnt - i - 1u;
        s->es[i] = asc_s[src];
        s->ee[i] = asc_e[src];
    }
    left_pad = (uint32_t)kprop_bounded(t, 20u);
    right_pad = (uint32_t)kprop_bounded(t, 20u);
    s->tstart = asc_s[0] - left_pad;
    s->tend = asc_e[cnt - 1u] + right_pad;
    if (kprop_bounded(t, 3u) != 0u) {
        uint32_t span = s->tend - s->tstart;
        s->cds_s = s->tstart + (uint32_t)kprop_bounded(t, span + 1u);
        s->cds_e = s->cds_s + (uint32_t)kprop_bounded(t,
            s->tend - s->cds_s + 1u);
    }
    s->splice_exonic = (uint32_t)kprop_bounded(t, 40u);
    s->splice_intronic = (uint32_t)kprop_bounded(t, 40u);
    s->excnt = (uint16_t)cnt;

    s->tx.chrom_id = &s->chrom;
    s->tx.start1 = &s->tstart;
    s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand;
    s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff;
    s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s;
    s->tx.cds_end1 = &s->cds_e;
    s->tx.transcript_count = 1u;
    s->ex.start1 = s->es;
    s->ex.end1 = s->ee;
    s->ex.exon_count = cnt;
    *instance = s;
    return THEFT_ALLOC_OK;
}

static void kprop_point_scene_free(void *instance, void *env) {
    (void)env;
    free(instance);
}

static struct theft_type_info kprop_point_scene_info = {
    .alloc = kprop_point_scene_alloc,
    .free = kprop_point_scene_free,
};

static int region_states_equal(const duckvep_region_state_t *a,
                               const duckvep_region_state_t *b) {
    return a->region_mask == b->region_mask &&
           a->within_feature == b->within_feature &&
           a->complete_overlap_feature == b->complete_overlap_feature &&
           a->complete_within_feature == b->complete_within_feature &&
           a->partial_overlap_feature == b->partial_overlap_feature &&
           a->within_cdna == b->within_cdna &&
           a->overlaps_exon == b->overlaps_exon &&
           a->overlaps_intron == b->overlaps_intron &&
           a->overlaps_cds == b->overlaps_cds &&
           a->overlaps_utr5 == b->overlaps_utr5 &&
           a->overlaps_utr3 == b->overlaps_utr3;
}

static int splice_states_equal(const duckvep_splice_state_t *a,
                               const duckvep_splice_state_t *b) {
    return a->splice_donor == b->splice_donor &&
           a->splice_acceptor == b->splice_acceptor &&
           a->splice_donor_5th == b->splice_donor_5th &&
           a->splice_donor_region == b->splice_donor_region &&
           a->splice_polypyrimidine == b->splice_polypyrimidine &&
           a->splice_region == b->splice_region &&
           a->intronic == b->intronic && a->any == b->any;
}

static enum theft_trial_res prop_sorted_point_classifier_matches_exhaustive(
    struct theft *t, void *arg1) {

    const struct kprop_point_scene *s =
        (const struct kprop_point_scene *)arg1;
    uint16_t rank = UINT16_MAX;
    uint32_t first = s->tstart > 24u ? s->tstart - 24u : 1u;
    uint32_t last = s->tend + 24u;
    uint32_t pos;
    (void)t;

    for (pos = first; pos <= last; pos++) {
        duckvep_region_state_t slow_region;
        duckvep_region_state_t fast_region;
        duckvep_splice_state_t slow_splice;
        duckvep_splice_state_t fast_splice;

        slow_region = duckvep_region_classify_span(
            &s->tx, &s->ex, 0u, pos, pos,
            s->splice_exonic, s->splice_intronic);
        slow_splice = duckvep_splice_classify_span(
            &s->tx, &s->ex, 0u, pos, pos, 0u);
        duckvep_classify_point_sorted(
            &s->tx, &s->ex, 0u, pos,
            s->splice_exonic, s->splice_intronic,
            &rank, &fast_region, &fast_splice);
        if (!region_states_equal(&slow_region, &fast_region) ||
            !splice_states_equal(&slow_splice, &fast_splice))
            return THEFT_TRIAL_FAIL;
        if (pos == UINT32_MAX) break;
    }
    return THEFT_TRIAL_PASS;
}

TEST sorted_point_classifier_matches_exhaustive_for_any_transcript(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "sorted point cursor classifier == exhaustive exon/gap scans";
    cfg.prop1 = prop_sorted_point_classifier_matches_exhaustive;
    cfg.type_info[0] = &kprop_point_scene_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

/* Hand-computed: coding transcript [1000,2000], exons [1000,1200] & [1500,2000],
 * CDS [1100,1900], splice_dist 8. tx0 is '+' strand, tx1 is '-' strand. */
TEST region_mask_known_scene(void) {
    static uint16_t chrom[2]  = {0u, 0u};
    static uint32_t tstart[2] = {1000u, 1000u};
    static uint32_t tend[2]   = {2000u, 2000u};
    static int8_t   strand[2] = {(int8_t)1, (int8_t)-1};
    static uint64_t flags[2]  = {0u, 0u};
    static uint32_t exoff[2]  = {0u, 0u};
    static uint16_t excnt[2]  = {2u, 2u};
    static uint32_t cds_s[2]  = {1100u, 1100u};
    static uint32_t cds_e[2]  = {1900u, 1900u};
    static uint32_t es[2]     = {1000u, 1500u};
    static uint32_t ee[2]     = {1200u, 2000u};
    static uint32_t z2[2]     = {0u, 0u};
    static int8_t   zp2[2]    = {0, 0};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;

    memset(&tx, 0, sizeof tx); memset(&ex, 0, sizeof ex);
    tx.chrom_id = chrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = strand;
    tx.flags = flags; tx.exon_offset = exoff; tx.exon_count = excnt;
    tx.cds_start1 = cds_s; tx.cds_end1 = cds_e; tx.transcript_count = 2u;
    ex.start1 = es; ex.end1 = ee; ex.cdna_start1 = z2; ex.cdna_end1 = z2;
    ex.phase = zp2; ex.end_phase = zp2; ex.exon_count = 2u;

    /* '+' strand (tx0) */
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_UPSTREAM,   duckvep_region_mask(&tx, &ex, 0, 900u, 3u, 8u));
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_DOWNSTREAM, duckvep_region_mask(&tx, &ex, 0, 2100u, 3u, 8u));
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_CDS,        duckvep_region_mask(&tx, &ex, 0, 1150u, 3u, 8u));
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_UTR,        duckvep_region_mask(&tx, &ex, 0, 1050u, 3u, 8u));
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_INTRON,     duckvep_region_mask(&tx, &ex, 0, 1300u, 3u, 8u));
    ASSERT_EQ((uint32_t)(DUCKVEP_REGION_INTRON | DUCKVEP_REGION_SPLICE),
                                                   duckvep_region_mask(&tx, &ex, 0, 1205u, 3u, 8u));
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_UTR,        duckvep_region_mask(&tx, &ex, 0, 1950u, 3u, 8u));
    /* Asymmetric splice reach pinned: exon0 ends at 1200 (CDS). 2 bp INSIDE the
     * exon is within the exonic reach (3) -> splice; 5 bp inside is not; the 1205
     * case above is 5 bp into the intron, within the intronic reach (8). Flip the
     * two reaches and the exon-side call gains splice while the intron-side loses
     * it -- proving the two controls are independent and wired. */
    ASSERT_EQ((uint32_t)(DUCKVEP_REGION_CDS | DUCKVEP_REGION_SPLICE),
                                                   duckvep_region_mask(&tx, &ex, 0, 1198u, 3u, 8u));
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_CDS,        duckvep_region_mask(&tx, &ex, 0, 1195u, 3u, 8u));
    /* intron side at 1205 (5 bp in): splice with intronic 8, NOT with intronic 3 */
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_INTRON,     duckvep_region_mask(&tx, &ex, 0, 1205u, 8u, 3u));
    /* '-' strand (tx1): up/down flip */
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_DOWNSTREAM, duckvep_region_mask(&tx, &ex, 1, 900u, 3u, 8u));
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_UPSTREAM,   duckvep_region_mask(&tx, &ex, 1, 2100u, 3u, 8u));
    PASS();
}

/* Regression (found by 300k-trial region-mask fuzz): the exonic splice reach must be
 * clamped to the exon so an exon SHORTER than splice_exonic does not spill its reach past
 * the opposite boundary and flag a flanking intron base as splice. Middle exon is 2 bp
 * [1500,1501] with exonic reach 4; positions 1499 (before) and 1503 (after) are intronic
 * and, with intronic reach 0, must NOT be splice — while in-exon reach and a non-zero
 * intronic reach still fire. */
TEST region_mask_short_exon_splice_reach_clamped(void) {
    static uint16_t chrom[1]  = {0u};
    static uint32_t tstart[1] = {1000u};
    static uint32_t tend[1]   = {2000u};
    static int8_t   strand[1] = {(int8_t)1};
    static uint64_t flags[1]  = {0u};
    static uint32_t exoff[1]  = {0u};
    static uint16_t excnt[1]  = {3u};
    static uint32_t cds_s[1]  = {1100u};
    static uint32_t cds_e[1]  = {1900u};
    static uint32_t es[3]     = {1000u, 1500u, 1800u};
    static uint32_t ee[3]     = {1200u, 1501u, 2000u};
    static uint32_t z3[3]     = {0u, 0u, 0u};
    static int8_t   zp3[3]    = {0, 0, 0};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;

    memset(&tx, 0, sizeof tx); memset(&ex, 0, sizeof ex);
    tx.chrom_id = chrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = strand;
    tx.flags = flags; tx.exon_offset = exoff; tx.exon_count = excnt;
    tx.cds_start1 = cds_s; tx.cds_end1 = cds_e; tx.transcript_count = 1u;
    ex.start1 = es; ex.end1 = ee; ex.cdna_start1 = z3; ex.cdna_end1 = z3;
    ex.phase = zp3; ex.end_phase = zp3; ex.exon_count = 3u;

    /* left spillover (right_lo clamp): 1499 is 1 bp before the 2 bp exon, intronic 0 */
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_INTRON, duckvep_region_mask(&tx, &ex, 0, 1499u, 4u, 0u));
    /* right spillover (left_hi clamp): 1503 is 2 bp after the 2 bp exon, intronic 0 */
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_INTRON, duckvep_region_mask(&tx, &ex, 0, 1503u, 4u, 0u));
    /* legitimate intronic reach still fires when intronic > 0 */
    ASSERT_EQ((uint32_t)(DUCKVEP_REGION_INTRON | DUCKVEP_REGION_SPLICE),
              duckvep_region_mask(&tx, &ex, 0, 1503u, 4u, 8u));
    /* in-exon reach still fires (position inside the short exon, in CDS) */
    ASSERT_EQ((uint32_t)(DUCKVEP_REGION_CDS | DUCKVEP_REGION_SPLICE),
              duckvep_region_mask(&tx, &ex, 0, 1500u, 4u, 0u));
    PASS();
}

TEST region_span_can_cross_cds_intron_and_splice_windows(void) {
    static const uint16_t chrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1] = {300u};
    static const int8_t strand[1] = {1};
    static const uint64_t flags[1] = {0u};
    static const uint32_t exoff[1] = {0u};
    static const uint16_t excnt[1] = {2u};
    static const uint32_t cds_s[1] = {120u};
    static const uint32_t cds_e[1] = {280u};
    static const uint32_t es[2] = {100u, 250u};
    static const uint32_t ee[2] = {150u, 300u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_region_state_t state;
    duckvep_splice_state_t splice;

    memset(&tx, 0, sizeof tx);
    memset(&ex, 0, sizeof ex);
    tx.chrom_id = chrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = strand;
    tx.flags = flags; tx.exon_offset = exoff; tx.exon_count = excnt;
    tx.cds_start1 = cds_s; tx.cds_end1 = cds_e; tx.transcript_count = 1u;
    ex.start1 = es; ex.end1 = ee; ex.exon_count = 2u;

    state = duckvep_region_classify_span(&tx, &ex, 0u, 145u, 255u, 3u, 8u);
    splice = duckvep_splice_classify_span(&tx, &ex, 0u, 145u, 255u, 0u);
    ASSERT(state.within_feature);
    ASSERT(state.complete_within_feature);
    ASSERT(state.within_cdna);
    ASSERT(state.overlaps_exon);
    ASSERT(state.overlaps_intron);
    ASSERT(state.overlaps_cds);
    ASSERT((state.region_mask & (uint32_t)DUCKVEP_REGION_CDS) != 0u);
    ASSERT((state.region_mask & (uint32_t)DUCKVEP_REGION_INTRON) != 0u);
    ASSERT(splice.splice_donor);
    ASSERT(splice.splice_acceptor);
    ASSERT(splice.any);
    PASS();
}

/* ===================================================================== *
 * VEP-source-grounded splice-SITE classification (effect-ctx slice 1).
 *
 * Anchors derived from the VEP 116 _intron_effects / _intron_overlap geometry
 * (see duckvep_classify.h): boundary points (is, is+1, ie-1, ie), the 5th-base
 * and donor-region edges, the polypyrimidine∧splice_region co-fire zone, deep
 * intron, BOTH introns, and the outer transcript ends. These guard regressions;
 * they do NOT prove VEP conformance (no --gff differential yet). Scene: a 3-exon
 * transcript so there are two real introns; tx0 is '+' strand, tx1 is the same
 * model on '-' strand (donor/acceptor swap). Intron0 = [1201,1400] (is=1201,
 * ie=1400), intron1 = [1601,1800]. The outer 5'/3' ends (1000, 2000) MUST NOT be
 * splice sites — the over-call fix.
 * ===================================================================== */
TEST splice_classify_known_scene(void) {
    static uint16_t chrom[2]  = {0u, 0u};
    static uint32_t tstart[2] = {1000u, 1000u};
    static uint32_t tend[2]   = {2000u, 2000u};
    static int8_t   strand[2] = {(int8_t)1, (int8_t)-1};
    static uint64_t flags[2]  = {0u, 0u};
    static uint32_t exoff[2]  = {0u, 0u};
    static uint16_t excnt[2]  = {3u, 3u};
    static uint32_t cds_s[2]  = {1100u, 1100u};
    static uint32_t cds_e[2]  = {1900u, 1900u};
    static uint32_t es[3]     = {1000u, 1401u, 1801u};
    static uint32_t ee[3]     = {1200u, 1600u, 2000u};
    static uint32_t z3[3]     = {0u, 0u, 0u};
    static int8_t   zp3[3]    = {0, 0, 0};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_splice_state_t s;

    memset(&tx, 0, sizeof tx); memset(&ex, 0, sizeof ex);
    tx.chrom_id = chrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = strand;
    tx.flags = flags; tx.exon_offset = exoff; tx.exon_count = excnt;
    tx.cds_start1 = cds_s; tx.cds_end1 = cds_e; tx.transcript_count = 2u;
    ex.start1 = es; ex.end1 = ee; ex.cdna_start1 = z3; ex.cdna_end1 = z3;
    ex.phase = zp3; ex.end_phase = zp3; ex.exon_count = 3u;

    /* --- '+' strand (tx0), intron0 = [1201,1400] --- */
    s = duckvep_splice_classify(&tx, &ex, 0, 1201u); /* essential donor (is)        */
    ASSERT(s.splice_donor && !s.splice_acceptor && !s.intronic && s.any);
    s = duckvep_splice_classify(&tx, &ex, 0, 1202u); /* essential donor (is+1)      */
    ASSERT(s.splice_donor && !s.splice_acceptor && !s.intronic);
    s = duckvep_splice_classify(&tx, &ex, 0, 1399u); /* essential acceptor (ie-1)   */
    ASSERT(s.splice_acceptor && !s.splice_donor && !s.intronic);
    s = duckvep_splice_classify(&tx, &ex, 0, 1400u); /* essential acceptor (ie)     */
    ASSERT(s.splice_acceptor && !s.splice_donor && !s.intronic && s.any);
    s = duckvep_splice_classify(&tx, &ex, 0, 1205u); /* 5th base (is+4) -> 5th+intronic, region/donor_region suppressed */
    ASSERT(s.splice_donor_5th && s.intronic && !s.splice_donor_region && !s.splice_region);
    s = duckvep_splice_classify(&tx, &ex, 0, 1203u); /* donor region (is+2), not 5th */
    ASSERT(s.splice_donor_region && s.intronic && !s.splice_donor_5th && !s.splice_region);
    s = duckvep_splice_classify(&tx, &ex, 0, 1206u); /* donor region edge (is+5)    */
    ASSERT(s.splice_donor_region && !s.splice_donor_5th);
    s = duckvep_splice_classify(&tx, &ex, 0, 1207u); /* past donor region -> splice_region (is+6) */
    ASSERT(!s.splice_donor_region && s.splice_region && s.intronic);
    s = duckvep_splice_classify(&tx, &ex, 0, 1390u); /* polypyrimidine, NOT region [ie-16,ie-8] */
    ASSERT(s.splice_polypyrimidine && s.intronic && !s.splice_region);
    s = duckvep_splice_classify(&tx, &ex, 0, 1395u); /* co-fire: region AND polypyrimidine (VEP does not suppress) */
    ASSERT(s.splice_region && s.splice_polypyrimidine && s.intronic);
    s = duckvep_splice_classify(&tx, &ex, 0, 1199u); /* exon-side splice_region only */
    ASSERT(s.splice_region && !s.intronic && !s.splice_donor && !s.splice_donor_region);
    s = duckvep_splice_classify(&tx, &ex, 0, 1300u); /* deep intron: intronic, no splice */
    ASSERT(s.intronic && !s.any);

    /* VEP compares the feature alleles bytewise and visits only contiguous
     * mismatch islands. Retained bases between two edits must not turn the
     * whole enclosing span into a splice event. */
    {
        uint8_t ref[20];
        uint8_t alt[20];
        memset(ref, 'A', sizeof ref);
        memcpy(alt, ref, sizeof alt);
        alt[0] = 'C';  /* 1190: exon, outside the exon-side splice window */
        alt[19] = 'C'; /* 1209: intronic, past the donor-side region */
        s = duckvep_splice_classify_differing_regions(
            &tx, &ex, 0u, 1190u, ref, 20u, alt, 20u);
        ASSERT(s.intronic && !s.any);

        memcpy(alt, ref, sizeof alt);
        alt[13] = 'C'; /* 1203: donor region */
        alt[15] = 'C'; /* 1205: donor fifth base */
        s = duckvep_splice_classify_differing_regions(
            &tx, &ex, 0u, 1190u, ref, 20u, alt, 20u);
        ASSERT(s.splice_donor_5th && !s.splice_donor_region && s.intronic);
    }

    /* VEP's interval-tree path caches an intron when the REF-shaped feature is
     * within three exonic bases of it. A longer ALT mismatch island may then
     * extend through the essential donor into the intronic interior. Moving
     * that same REF feature one base beyond the cache flank suppresses every
     * intron-derived predicate, even though the ALT island reaches the gap. */
    {
        static const uint8_t ref[3] = {'A', 'A', 'A'};
        static const uint8_t alt[11] = {
            'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'
        };

        s = duckvep_splice_classify_differing_regions(
            &tx, &ex, 0u, 1198u, ref, 3u, alt, 11u);
        ASSERT(s.splice_donor && s.splice_donor_5th && s.intronic);

        s = duckvep_splice_classify_differing_regions(
            &tx, &ex, 0u, 1195u, ref, 3u, alt, 11u);
        ASSERT(!s.any && !s.intronic);
    }

    /* --- second intron1 = [1601,1800] is classified too (not just intron0) --- */
    s = duckvep_splice_classify(&tx, &ex, 0, 1601u); /* donor of intron1            */
    ASSERT(s.splice_donor && !s.splice_acceptor);
    s = duckvep_splice_classify(&tx, &ex, 0, 1800u); /* acceptor of intron1         */
    ASSERT(s.splice_acceptor && !s.splice_donor);

    /* --- over-call fix: outer transcript ends are NOT splice sites --- */
    s = duckvep_splice_classify(&tx, &ex, 0, 1000u); /* outer 5' end                */
    ASSERT(!s.any && !s.intronic);
    s = duckvep_splice_classify(&tx, &ex, 0, 2000u); /* outer 3' end                */
    ASSERT(!s.any && !s.intronic);

    /* --- '-' strand (tx1): donor/acceptor + 5th/region swap to the other end --- */
    s = duckvep_splice_classify(&tx, &ex, 1, 1201u); /* raw start_ss -> acceptor    */
    ASSERT(s.splice_acceptor && !s.splice_donor);
    s = duckvep_splice_classify(&tx, &ex, 1, 1400u); /* raw end_ss -> donor         */
    ASSERT(s.splice_donor && !s.splice_acceptor);
    s = duckvep_splice_classify(&tx, &ex, 1, 1396u); /* ie-4 -> 5th base on '-'      */
    ASSERT(s.splice_donor_5th && !s.splice_donor_region);
    s = duckvep_splice_classify(&tx, &ex, 1, 1397u); /* ie-3 -> donor region on '-' (not 5th) */
    ASSERT(s.splice_donor_region && !s.splice_donor_5th);
    s = duckvep_splice_classify(&tx, &ex, 1, 1210u); /* [is+2,is+16] -> polypyr on '-' */
    ASSERT(s.splice_polypyrimidine);
    s = duckvep_splice_classify(&tx, &ex, 0, 1403u); /* exon-side splice_region, ie+3 arm (acceptor side) */
    ASSERT(s.splice_region && !s.intronic && !s.splice_acceptor);

    /* --- '-' strand stored in TRANSCRIPT order (descending genomic): introns must
     * still be recovered order-independently. With the old is=end[k]+1 /
     * ie=start[k+1]-1 assumption this model would skip every intron (ie<is) and
     * emit NO splice site — this anchor is the regression guard for that bug. The
     * physical introns are the same genomic gaps [1201,1400] and [1601,1800]. --- */
    {
        static uint16_t tchrom[1] = {0u};
        static uint32_t tts[1]    = {1000u};
        static uint32_t tte[1]    = {2000u};
        static int8_t   tstr[1]   = {(int8_t)-1};
        static uint64_t tfl[1]    = {0u};
        static uint32_t teo[1]    = {0u};
        static uint16_t tec[1]    = {3u};
        static uint32_t tcs[1]    = {1100u};
        static uint32_t tce[1]    = {1900u};
        /* exons in transcript (descending genomic) order */
        static uint32_t res[3]    = {1801u, 1401u, 1000u};
        static uint32_t ree[3]    = {2000u, 1600u, 1200u};
        duckvep_transcript_model_t rtx;
        duckvep_exon_model_t rex;
        duckvep_splice_state_t r;

        memset(&rtx, 0, sizeof rtx); memset(&rex, 0, sizeof rex);
        rtx.chrom_id = tchrom; rtx.start1 = tts; rtx.end1 = tte; rtx.strand = tstr;
        rtx.flags = tfl; rtx.exon_offset = teo; rtx.exon_count = tec;
        rtx.cds_start1 = tcs; rtx.cds_end1 = tce; rtx.transcript_count = 1u;
        rex.start1 = res; rex.end1 = ree; rex.cdna_start1 = z3; rex.cdna_end1 = z3;
        rex.phase = zp3; rex.end_phase = zp3; rex.exon_count = 3u;

        r = duckvep_splice_classify(&rtx, &rex, 0, 1400u); /* intron0 ie -> donor on '-' */
        ASSERT(r.splice_donor && !r.splice_acceptor);
        r = duckvep_splice_classify(&rtx, &rex, 0, 1201u); /* intron0 is -> acceptor on '-' */
        ASSERT(r.splice_acceptor && !r.splice_donor);
        r = duckvep_splice_classify(&rtx, &rex, 0, 1396u); /* ie-4 -> 5th base on '-'   */
        ASSERT(r.splice_donor_5th);
        r = duckvep_splice_classify(&rtx, &rex, 0, 1800u); /* intron1 ie -> donor on '-' */
        ASSERT(r.splice_donor && !r.splice_acceptor);
    }

    /* --- short intron: a 6 bp intron [201,206]. Essential donor/acceptor and the
     * intronic interior are pinned; the wider donor-region / PPT windows are
     * deliberately NOT clipped to the intron (VEP's _intron_effects arithmetic is
     * also unclipped), but we only assert the unambiguous sites here. --- */
    {
        static uint16_t schrom[1] = {0u};
        static uint32_t sts[1]    = {100u};
        static uint32_t ste[1]    = {300u};
        static int8_t   sstr[1]   = {(int8_t)1};
        static uint64_t sfl[1]    = {0u};
        static uint32_t seo[1]    = {0u};
        static uint16_t sec[1]    = {2u};
        static uint32_t scs[1]    = {120u};
        static uint32_t sce[1]    = {280u};
        static uint32_t ses[2]    = {100u, 207u};
        static uint32_t see[2]    = {200u, 300u};
        static uint32_t sz2[2]    = {0u, 0u};
        static int8_t   szp2[2]   = {0, 0};
        duckvep_transcript_model_t stx;
        duckvep_exon_model_t sex;
        duckvep_splice_state_t r;

        memset(&stx, 0, sizeof stx); memset(&sex, 0, sizeof sex);
        stx.chrom_id = schrom; stx.start1 = sts; stx.end1 = ste; stx.strand = sstr;
        stx.flags = sfl; stx.exon_offset = seo; stx.exon_count = sec;
        stx.cds_start1 = scs; stx.cds_end1 = sce; stx.transcript_count = 1u;
        sex.start1 = ses; sex.end1 = see; sex.cdna_start1 = sz2; sex.cdna_end1 = sz2;
        sex.phase = szp2; sex.end_phase = szp2; sex.exon_count = 2u;

        r = duckvep_splice_classify(&stx, &sex, 0, 201u); /* is -> donor                */
        ASSERT(r.splice_donor && !r.splice_acceptor);
        r = duckvep_splice_classify(&stx, &sex, 0, 206u); /* ie -> acceptor             */
        ASSERT(r.splice_acceptor && !r.splice_donor);
        r = duckvep_splice_classify(&stx, &sex, 0, 204u); /* is+3=ie-2 -> intronic interior */
        ASSERT(r.intronic);
    }
    PASS();
}

/* ===================================================================== *
 * Insertion (interbase) splice classification — VEP _intron_effects rule.
 *
 * VEP models a pure insertion as vf->start = P+1, vf->end = P (start > end), where
 * P is the anchor base to its 5' side. A zone [lo,hi] is touched only when BOTH
 * flanking bases fall inside (P in [lo, hi-1]); it also fires exact-edge special
 * cases at the exon/intron and donor/acceptor rim. Every anchor below is the
 * insertion analogue of a point anchor in splice_classify_known_scene, and the
 * expected value is hand-derived from VEP 116 BaseTranscriptVariationAllele.pm.
 * Scene: 2-exon transcript, intron = [1201,1400] (is=1201, ie=1400); tx0 '+',
 * tx1 '-'. These pin the ClinVar 10:53827390 T>TA family fix (an exon/intron
 * boundary insertion is splice_region, NOT splice_donor_region).
 * ===================================================================== */
TEST splice_classify_insertion_interbase_scene(void) {
    static uint16_t chrom[2]  = {0u, 0u};
    static uint32_t tstart[2] = {1000u, 1000u};
    static uint32_t tend[2]   = {2000u, 2000u};
    static int8_t   strand[2] = {(int8_t)1, (int8_t)-1};
    static uint64_t flags[2]  = {0u, 0u};
    static uint32_t exoff[2]  = {0u, 0u};
    static uint16_t excnt[2]  = {2u, 2u};
    static uint32_t cds_s[2]  = {1100u, 1100u};
    static uint32_t cds_e[2]  = {1900u, 1900u};
    static uint32_t es[2]     = {1000u, 1401u};
    static uint32_t ee[2]     = {1200u, 2000u};
    static uint32_t z2[2]     = {0u, 0u};
    static int8_t   zp2[2]    = {0, 0};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_splice_state_t p; /* point */
    duckvep_splice_state_t i; /* insertion (interbase) */

    memset(&tx, 0, sizeof tx); memset(&ex, 0, sizeof ex);
    tx.chrom_id = chrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = strand;
    tx.flags = flags; tx.exon_offset = exoff; tx.exon_count = excnt;
    tx.cds_start1 = cds_s; tx.cds_end1 = cds_e; tx.transcript_count = 2u;
    ex.start1 = es; ex.end1 = ee; ex.cdna_start1 = z2; ex.cdna_end1 = z2;
    ex.phase = zp2; ex.end_phase = zp2; ex.exon_count = 2u;

    /* --- '+' strand, intron [1201,1400] --- */
    /* exon-end/intron-start boundary insertion (P=1200, between exon base 1200 and
     * is=1201): VEP fires splice_region via the rs==is special edge, NOT any donor
     * term. The point at 1200 is a plain exon-side splice_region too, but the
     * insertion reaches it only through the special-edge, so both must agree here. */
    i = duckvep_splice_classify_span(&tx, &ex, 0, 1200u, 1200u, 1u);
    ASSERT(i.splice_region && !i.splice_donor && !i.splice_donor_region && !i.intronic);

    /* essential donor: only the insertion straddling is..is+1 (P=1201) is donor */
    i = duckvep_splice_classify_span(&tx, &ex, 0, 1201u, 1201u, 1u);
    ASSERT(i.splice_donor && !i.splice_acceptor && !i.splice_region);

    /* is+4 (1205): POINT is the donor 5th base; the INSERTION is donor_region
     * (fifth-base single-point zone can never contain both insertion flanks). */
    p = duckvep_splice_classify(&tx, &ex, 0, 1205u);
    ASSERT(p.splice_donor_5th && !p.splice_donor_region);
    i = duckvep_splice_classify_span(&tx, &ex, 0, 1205u, 1205u, 1u);
    ASSERT(i.splice_donor_region && !i.splice_donor_5th);

    /* is+5 (1206): POINT is donor_region (top edge); the INSERTION shrinks off that
     * edge (needs P<=is+4) and becomes plain splice_region. */
    p = duckvep_splice_classify(&tx, &ex, 0, 1206u);
    ASSERT(p.splice_donor_region);
    i = duckvep_splice_classify_span(&tx, &ex, 0, 1206u, 1206u, 1u);
    ASSERT(i.splice_region && !i.splice_donor_region);

    /* --- '-' strand: the exact ClinVar 10:53827390 T>TA pattern. On '-' the donor
     * region is at the ie end. A POINT at ie-2 (1398) is splice_donor_region; the
     * INSERTION anchored there is splice_region (rim shrink + re==ie-2 edge). --- */
    p = duckvep_splice_classify(&tx, &ex, 1, 1398u);
    ASSERT(p.splice_donor_region);
    i = duckvep_splice_classify_span(&tx, &ex, 1, 1398u, 1398u, 1u);
    ASSERT(i.splice_region && !i.splice_donor_region && !i.splice_donor && i.intronic);
    PASS();
}

/* ===================================================================== *
 * Polypyrimidine-tract exon gate. VEP's splice_polypyrimidine_tract_variant
 * OverlapConsequence carries include => { exon => 0, intron => 1 }: the tract
 * term is emitted only when the variant does NOT overlap an exon. A deletion
 * running from the tract across the acceptor into the coding exon keeps
 * splice_acceptor but drops polypyrimidine (ClinVar X:134393492 family).
 * Scene: 2 exons [100,200]+[301,500], intron [201,300] (is=201, ie=300);
 * '+' strand acceptor at ie, tract [ie-16,ie-2] = [284,298].
 * ===================================================================== */
TEST splice_ppt_exon_gate_scene(void) {
    static uint16_t chrom[1]  = {0u};
    static uint32_t tstart[1] = {100u};
    static uint32_t tend[1]   = {500u};
    static int8_t   strand[1] = {(int8_t)1};
    static uint64_t flags[1]  = {0u};
    static uint32_t exoff[1]  = {0u};
    static uint16_t excnt[1]  = {2u};
    static uint32_t cds_s[1]  = {150u};
    static uint32_t cds_e[1]  = {450u};
    static uint32_t es[2]     = {100u, 301u};
    static uint32_t ee[2]     = {200u, 500u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_effect_ctx_t ctx;

    memset(&tx, 0, sizeof tx); memset(&ex, 0, sizeof ex);
    tx.chrom_id = chrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = strand;
    tx.flags = flags; tx.exon_offset = exoff; tx.exon_count = excnt;
    tx.cds_start1 = cds_s; tx.cds_end1 = cds_e; tx.transcript_count = 1u;
    ex.start1 = es; ex.end1 = ee; ex.exon_count = 2u;

    /* deletion inside the tract, wholly intronic: PPT fires, no exon overlap */
    duckvep_effect_ctx_fill(&tx, &ex, 0u, 0u, 285u, 290u, 0u, 3u, 8u, &ctx);
    ASSERT(!ctx.region_state.overlaps_exon);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_PPT)) != 0u);

    /* deletion from the tract across the acceptor into the exon: exon overlap
     * suppresses PPT, but the essential acceptor site still fires */
    duckvep_effect_ctx_fill(&tx, &ex, 0u, 0u, 285u, 305u, 0u, 3u, 8u, &ctx);
    ASSERT(ctx.region_state.overlaps_exon);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_PPT)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_ACCEPTOR)) != 0u);
    PASS();
}

/* ===================================================================== *
 * The VEP-shaped consequence rule table (duckvep_effect_eval).
 *
 * Pins the pre-bits -> SO mapping DIRECTLY, independent of region_mask /
 * splice_classify: feed hand-built pre-bit sets and assert the emitted SO bitset.
 * This is the decision layer that replaced the old structural_consequence_mask
 * if/else + the splice OR block + the codon SO mapping — all three now live as
 * rows in one table, so this anchor guards the whole consequence assembly. The key
 * interactions: the finalized CODING_UNKNOWN predicate emits
 * coding_sequence_variant unless a sequence delta refined it; the specific codon
 * term then comes from its fact bit. Complete transcript overlap instead selects
 * the transcript-level coding/non-coding fallback predicates.
 * ===================================================================== */
TEST effect_rule_table_known_pre_bits(void) {
    /* upstream / downstream are exclusive single terms */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_UPSTREAM_GENE),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_UPSTREAM) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_CODING)));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_DOWNSTREAM_GENE),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_DOWNSTREAM)));

    /* Finalized CDS, no delta -> generic coding_sequence_variant. */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_CDS) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_CODING) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_CODING_UNKNOWN)));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_CODING_TRANSCRIPT),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_CODING) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_CODING_TRANSCRIPT)));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_NONCODING) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_NONCODING_GENE)));
    /* CDS + a refined delta -> the specific codon term, coding_sequence SUPPRESSED */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_CDS) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_DELTA) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_MISSENSE)));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_CDS) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_DELTA) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_STOP_RETAINED)));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_START_LOST) | DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_CDS) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_DELTA) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_START_LOST) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_STOP_GAINED)));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_CDS) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_DELTA) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_FRAMESHIFT)));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_CDS) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_DELTA) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_INFRAME_DELETION)));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_CDS) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_DELTA) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_INFRAME_INSERTION)));
    /* VEP coding_unknown does not exclude inframe_insertion. A terminal insertion
     * whose local peptide ends in X therefore emits both terms. */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
              DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_CDS) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_DELTA) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_CODING_UNKNOWN) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_INFRAME_INSERTION)));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_CDS) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_DELTA) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_PROTEIN_ALTERING)));

    /* UTR sides */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_UTR5)));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_UTR3)));

    /* non-coding EXON -> the exon term ONLY (NOT non_coding_transcript_variant):
     * VEP within_non_coding_gene excludes non_coding_exon_variant — the two terms are
     * mutually exclusive (VariationEffect.pm:495; regression for the rule-120 fix). */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT_EXON),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_EXON) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_NONCODING) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_NONCODING_EXON)));

    /* intron_variant requires VEP within_intron (PRE_WITHIN_INTRON = _intron_effects->{intronic}),
     * NOT mere placement. Deep intron: within_intron true. non_coding_transcript uses PLACEMENT. */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INTRON),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_INTRON) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_INTRON) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_CODING)));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INTRON) | DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_INTRON) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_INTRON) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_NONCODING) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_NONCODING_GENE)));

    /* Essential splice donor/acceptor: intron PLACEMENT but within_intron FALSE (the
     * dinucleotides are start/end_splice_site, not intronic) -> intron_variant does NOT fire,
     * so the splice term emits alone. */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_SPLICE_DONOR),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_INTRON) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_CODING) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_DONOR)));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_SPLICE_ACCEPTOR),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_INTRON) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_CODING) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_ACCEPTOR)));
    /* 5th base / donor region: within_intron TRUE there -> intron_variant co-emits. */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INTRON) | DUCKVEP_SO(DUCKVEP_SO_SPLICE_DONOR_REGION),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_INTRON) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_INTRON) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_CODING) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_DONOR_REGION)));
    /* Non-coding transcript at an essential splice site: the finalized
     * within_non_coding_gene predicate co-emits with the splice term, while
     * within_intron remains false so intron_variant does not fire. */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_SPLICE_DONOR) | DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_INTRON) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_NONCODING) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_NONCODING_GENE) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_DONOR)));
    /* GENERIC-SPAN readiness (the point of the formalization): a future deletion crossing both
     * an essential dinucleotide AND deep intron sets the splice bit AND within_intron -> BOTH
     * terms emit, which the SNV-point forbidden-mask could NOT express. Pinned at the fact level. */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INTRON) | DUCKVEP_SO(DUCKVEP_SO_SPLICE_DONOR),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_INTRON) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_INTRON) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_CODING) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_DONOR)));
    PASS();
}

/* VEP's NMD_transcript_variant means "inside a transcript whose curated
 * biotype is nonsense_mediated_decay". It is independent of the ordinary
 * region consequence and does not extend into the upstream/downstream halo. */
TEST nmd_transcript_predicate_known_scene(void) {
    static uint16_t chrom[1] = {0u};
    static uint32_t tstart[1] = {100u};
    static uint32_t tend[1] = {300u};
    static int8_t strand[1] = {1};
    static uint64_t flags[1] = {(uint64_t)DUCKVEP_TX_BIOTYPE_NMD};
    static uint32_t exoff[1] = {0u};
    static uint16_t excnt[1] = {2u};
    static uint32_t cds_start[1] = {110u};
    static uint32_t cds_end[1] = {290u};
    static uint32_t exon_start[2] = {100u, 250u};
    static uint32_t exon_end[2] = {150u, 300u};
    static uint32_t cdna_start[2] = {1u, 52u};
    static uint32_t cdna_end[2] = {51u, 102u};
    static int8_t phase[2] = {0, 0};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_effect_ctx_t ctx;
    uint64_t mask;

    memset(&tx, 0, sizeof tx);
    memset(&ex, 0, sizeof ex);
    tx.chrom_id = chrom;
    tx.start1 = tstart;
    tx.end1 = tend;
    tx.strand = strand;
    tx.flags = flags;
    tx.exon_offset = exoff;
    tx.exon_count = excnt;
    tx.cds_start1 = cds_start;
    tx.cds_end1 = cds_end;
    tx.transcript_count = 1u;
    ex.start1 = exon_start;
    ex.end1 = exon_end;
    ex.cdna_start1 = cdna_start;
    ex.cdna_end1 = cdna_end;
    ex.phase = phase;
    ex.end_phase = phase;
    ex.exon_count = 2u;

    duckvep_effect_ctx_fill(&tx, &ex, 0u, 0u, 200u, 200u, 0u, 3u, 8u,
                            &ctx);
    duckvep_effect_ctx_finalize(&ctx);
    mask = duckvep_effect_eval(ctx.pre_bits);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INTRON) |
              DUCKVEP_SO(DUCKVEP_SO_NMD_TRANSCRIPT), mask);

    duckvep_effect_ctx_fill(&tx, &ex, 0u, 0u, 90u, 90u, 0u, 3u, 8u,
                            &ctx);
    ASSERT((ctx.pre_bits &
            DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_NMD_TRANSCRIPT)) == 0u);
    PASS();
}

/* VEP Plugins release/116 NMD.pm is a separate prediction from the curated
 * NMD transcript biotype. Pin its executable coordinate thresholds on both
 * strands, including the source's inclusive 51-base penultimate-exon offset. */
TEST variant_induced_nmd_prediction_known_scene(void) {
    static uint16_t chrom[1] = {0u};
    static uint32_t tstart[1] = {100u};
    static uint32_t tend[1] = {599u};
    static int8_t strand[1] = {1};
    static uint64_t flags[1] = {0u};
    static uint32_t exoff[1] = {0u};
    static uint16_t excnt[1] = {3u};
    static uint32_t cds_start[1] = {100u};
    static uint32_t cds_end[1] = {599u};
    static uint32_t exon_start[3] = {100u, 300u, 500u};
    static uint32_t exon_end[3] = {199u, 399u, 599u};
    static uint32_t cdna_start[3] = {1u, 101u, 201u};
    static uint32_t cdna_end[3] = {100u, 200u, 300u};
    static int8_t phase[3] = {0, 0, 0};
    static uint32_t reverse_exon_start[3] = {500u, 300u, 100u};
    static uint32_t reverse_exon_end[3] = {599u, 399u, 199u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_event_t event;
    duckvep_nmd_result_t nmd;
    uint64_t stop_gained = DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED);

    memset(&tx, 0, sizeof tx);
    memset(&ex, 0, sizeof ex);
    memset(&event, 0, sizeof event);
    tx.chrom_id = chrom;
    tx.start1 = tstart;
    tx.end1 = tend;
    tx.strand = strand;
    tx.flags = flags;
    tx.exon_offset = exoff;
    tx.exon_count = excnt;
    tx.cds_start1 = cds_start;
    tx.cds_end1 = cds_end;
    tx.transcript_count = 1u;
    ex.start1 = exon_start;
    ex.end1 = exon_end;
    ex.cdna_start1 = cdna_start;
    ex.cdna_end1 = cdna_end;
    ex.phase = phase;
    ex.end_phase = phase;
    ex.exon_count = 3u;
    event.ref_diff_length = 1u;

    event.start1 = event.end1 = 320u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_PREDICTED_TRIGGERING, nmd.prediction);
    ASSERT_EQ(0u, nmd.escape_reasons);

    event.start1 = event.end1 = 300u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_PREDICTED_ESCAPING, nmd.prediction);
    ASSERT_EQ(DUCKVEP_NMD_ESCAPE_EARLY_CDS, nmd.escape_reasons);
    event.start1 = event.end1 = 301u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_PREDICTED_TRIGGERING, nmd.prediction);

    event.start1 = event.end1 = 348u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_ESCAPE_PENULTIMATE_EXON_END,
              nmd.escape_reasons);
    event.start1 = event.end1 = 347u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_PREDICTED_TRIGGERING, nmd.prediction);

    event.start1 = event.end1 = 500u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_ESCAPE_LAST_EXON, nmd.escape_reasons);

    excnt[0] = 1u;
    exon_end[0] = 599u;
    cdna_end[0] = 500u;
    event.start1 = event.end1 = 300u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_ESCAPE_INTRONLESS |
              DUCKVEP_NMD_ESCAPE_LAST_EXON, nmd.escape_reasons);
    excnt[0] = 3u;
    exon_end[0] = 199u;
    cdna_end[0] = 100u;

    event.start1 = event.end1 = 320u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event,
                        DUCKVEP_SO(DUCKVEP_SO_MISSENSE), &nmd);
    ASSERT_EQ(DUCKVEP_NMD_NOT_APPLICABLE, nmd.prediction);

    tx.exon_offset = NULL;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_UNRESOLVED, nmd.prediction);
    tx.exon_offset = exoff;

    event.start1 = event.end1 = 250u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event,
                        DUCKVEP_SO(DUCKVEP_SO_SPLICE_DONOR), &nmd);
    ASSERT_EQ(DUCKVEP_NMD_UNRESOLVED, nmd.prediction);

    strand[0] = -1;
    ex.start1 = reverse_exon_start;
    ex.end1 = reverse_exon_end;
    event.start1 = event.end1 = 351u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_ESCAPE_PENULTIMATE_EXON_END,
              nmd.escape_reasons);
    event.start1 = event.end1 = 352u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_PREDICTED_TRIGGERING, nmd.prediction);
    PASS();
}

/* VEP's tier machine is not severity ranking: a tier-1/2 match suppresses
 * later tiers, while every matching rule in the assigned tier co-emits. */
TEST effect_rule_tiers_suppress_only_later_tiers(void) {
    static const duckvep_consequence_rule_t rules[] = {
        { DUCKVEP_PRE(DUCKVEP_PRE_UPSTREAM),   0u, DUCKVEP_SO(DUCKVEP_SO_UPSTREAM_GENE),   1u, 40u, DUCKVEP_IMPACT_MODIFIER },
        { DUCKVEP_PRE(DUCKVEP_PRE_DOWNSTREAM), 0u, DUCKVEP_SO(DUCKVEP_SO_DOWNSTREAM_GENE), 2u,  1u, DUCKVEP_IMPACT_MODIFIER },
        { DUCKVEP_PRE(DUCKVEP_PRE_INTRON),     0u, DUCKVEP_SO(DUCKVEP_SO_INTRON),          2u, 99u, DUCKVEP_IMPACT_MODIFIER },
        { DUCKVEP_PRE(DUCKVEP_PRE_UTR5),       0u, DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR),     3u,  2u, DUCKVEP_IMPACT_MODIFIER },
        { DUCKVEP_PRE(DUCKVEP_PRE_UTR3),       0u, DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR),     3u,  3u, DUCKVEP_IMPACT_MODIFIER }
    };
    uint64_t pre;

    /* No tier-1 match; both matching tier-2 rules emit, regardless of rank, and
     * tier 3 is suppressed. */
    pre = DUCKVEP_PRE(DUCKVEP_PRE_DOWNSTREAM) |
          DUCKVEP_PRE(DUCKVEP_PRE_INTRON) |
          DUCKVEP_PRE(DUCKVEP_PRE_UTR5);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_DOWNSTREAM_GENE) | DUCKVEP_SO(DUCKVEP_SO_INTRON),
              duckvep_effect_eval_rules(pre, rules, sizeof rules / sizeof rules[0]));

    /* A tier-1 hit suppresses both tier 2 and tier 3. */
    pre |= DUCKVEP_PRE(DUCKVEP_PRE_UPSTREAM);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_UPSTREAM_GENE),
              duckvep_effect_eval_rules(pre, rules, sizeof rules / sizeof rules[0]));

    /* Without a suppressing tier, matching tier-3 rules co-emit. */
    pre = DUCKVEP_PRE(DUCKVEP_PRE_UTR5) | DUCKVEP_PRE(DUCKVEP_PRE_UTR3);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR) | DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR),
              duckvep_effect_eval_rules(pre, rules, sizeof rules / sizeof rules[0]));
    PASS();
}


/* One compatibility matrix is shared by the public borrowed-view boundary and
 * adapter tile. Compound operation/direction semantics require a richer event
 * type; reject them rather than synthesizing simultaneous gain/loss predicates. */
TEST sv_metadata_validity_matrix_known(void) {
    ASSERT(duckvep_sv_metadata_valid(DUCKVEP_SV_UNKNOWN,
                                     DUCKVEP_COPY_CHANGE_UNKNOWN));
    ASSERT(duckvep_sv_metadata_valid(DUCKVEP_SV_UNKNOWN,
                                     DUCKVEP_COPY_CHANGE_GAIN));
    ASSERT(duckvep_sv_metadata_valid(DUCKVEP_SV_CNV,
                                     DUCKVEP_COPY_CHANGE_LOSS));
    ASSERT(duckvep_sv_metadata_valid(DUCKVEP_SV_CNV,
                                     DUCKVEP_COPY_CHANGE_NEUTRAL));
    ASSERT(duckvep_sv_metadata_valid(DUCKVEP_SV_CNV,
                                     DUCKVEP_COPY_CHANGE_GAIN));
    ASSERT(duckvep_sv_metadata_valid(DUCKVEP_SV_INSERTION,
                                     DUCKVEP_COPY_CHANGE_UNKNOWN));
    ASSERT(duckvep_sv_metadata_valid(DUCKVEP_SV_INSERTION,
                                     DUCKVEP_COPY_CHANGE_NEUTRAL));
    ASSERT(duckvep_sv_metadata_valid(DUCKVEP_SV_INSERTION,
                                     DUCKVEP_COPY_CHANGE_GAIN));
    ASSERT(duckvep_sv_metadata_valid(DUCKVEP_SV_DELETION,
                                     DUCKVEP_COPY_CHANGE_LOSS));
    ASSERT(duckvep_sv_metadata_valid(DUCKVEP_SV_DUPLICATION,
                                     DUCKVEP_COPY_CHANGE_GAIN));
    ASSERT(duckvep_sv_metadata_valid(DUCKVEP_SV_INVERSION,
                                     DUCKVEP_COPY_CHANGE_NEUTRAL));
    ASSERT(duckvep_sv_metadata_valid(DUCKVEP_SV_BREAKEND,
                                     DUCKVEP_COPY_CHANGE_UNKNOWN));

    ASSERT(!duckvep_sv_metadata_valid(DUCKVEP_SV_NONE,
                                      DUCKVEP_COPY_CHANGE_UNKNOWN));
    ASSERT(!duckvep_sv_metadata_valid(DUCKVEP_SV_DELETION,
                                      DUCKVEP_COPY_CHANGE_GAIN));
    ASSERT(!duckvep_sv_metadata_valid(DUCKVEP_SV_DUPLICATION,
                                      DUCKVEP_COPY_CHANGE_LOSS));
    ASSERT(!duckvep_sv_metadata_valid(DUCKVEP_SV_INSERTION,
                                      DUCKVEP_COPY_CHANGE_LOSS));
    ASSERT(!duckvep_sv_metadata_valid(DUCKVEP_SV_INVERSION,
                                      DUCKVEP_COPY_CHANGE_GAIN));
    ASSERT(!duckvep_sv_metadata_valid(DUCKVEP_SV_BREAKEND,
                                      DUCKVEP_COPY_CHANGE_LOSS));
    ASSERT(!duckvep_sv_metadata_valid((duckvep_sv_type_t)255,
                                      DUCKVEP_COPY_CHANGE_UNKNOWN));
    ASSERT(!duckvep_sv_metadata_valid(DUCKVEP_SV_CNV,
                                      (duckvep_copy_change_t)255));
    PASS();
}

/* Structural operation, copy state, and geometry are separate facts. This pins
 * the VEP-shaped producer and proves that applying it preserves the operation
 * class bits needed by future include/predicate separation. */
TEST sv_predicate_facts_known(void) {
    duckvep_event_t event;
    duckvep_region_state_t region;
    duckvep_sv_effect_t sv;
    duckvep_effect_ctx_t ctx;

    memset(&event, 0, sizeof event);
    memset(&region, 0, sizeof region);
    memset(&ctx, 0, sizeof ctx);
    event.kind = (uint8_t)DUCKVEP_KIND_SV;
    event.sv_type = (uint8_t)DUCKVEP_SV_CNV;
    event.copy_change = (uint8_t)DUCKVEP_COPY_CHANGE_GAIN;
    region.complete_overlap_feature = 1u;
    sv = duckvep_sv_effect_fill(&event, &region);
    ASSERT(sv.copy_number_gain);
    ASSERT(sv.insertion);
    ASSERT(sv.feature_amplification);
    ASSERT(!sv.deletion);
    ASSERT(!sv.feature_ablation);
    duckvep_effect_ctx_apply_sv(&ctx, &sv);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INSERTION)) != 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_FEATURE_AMPLIFICATION)) != 0u);

    memset(&region, 0, sizeof region);
    memset(&ctx, 0, sizeof ctx);
    event.copy_change = (uint8_t)DUCKVEP_COPY_CHANGE_LOSS;
    region.within_cdna = 1u;
    region.partial_overlap_feature = 1u;
    sv = duckvep_sv_effect_fill(&event, &region);
    ASSERT(sv.copy_number_loss);
    ASSERT(sv.deletion);
    ASSERT(sv.feature_truncation);
    ASSERT(!sv.insertion);
    duckvep_effect_ctx_apply_sv(&ctx, &sv);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_DELETION)) != 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_FEATURE_TRUNCATION)) != 0u);

    memset(&region, 0, sizeof region);
    event.copy_change = (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN;
    region.complete_overlap_feature = 1u;
    sv = duckvep_sv_effect_fill(&event, &region);
    ASSERT(!sv.copy_number_gain);
    ASSERT(!sv.copy_number_loss);
    ASSERT(!sv.feature_amplification);
    ASSERT(!sv.feature_ablation);
    PASS();
}

/* Layer-1 rendering primitive (duckvep_so_render / duckvep_impact_name): the
 * adapter builds CSQ strings from the kernel's emitted mask via these. Pins the
 * '&'-joined VEP severity-rank render, the snprintf-style truncation return, and
 * the impact labels — testable with no DuckDB. */
TEST so_render_and_impact_name_known(void) {
    char buf[128];
    size_t n;

    /* empty mask -> "" , returns 0 */
    n = duckvep_so_render(0u, '&', buf, sizeof buf);
    ASSERT_EQ((size_t)0u, n);
    ASSERT_EQ(0, strcmp(buf, ""));

    /* single term */
    n = duckvep_so_render(DUCKVEP_SO(DUCKVEP_SO_INTRON), '&', buf, sizeof buf);
    ASSERT_EQ(0, strcmp(buf, "intron_variant"));
    ASSERT_EQ((size_t)14u, n);

    /* Stable bit indices do not encode severity: rendering follows VEP rank, so
     * splice_donor (rank 3) precedes intron (rank 28). */
    n = duckvep_so_render(DUCKVEP_SO(DUCKVEP_SO_SPLICE_DONOR) | DUCKVEP_SO(DUCKVEP_SO_INTRON),
                          '&', buf, sizeof buf);
    ASSERT_EQ(0, strcmp(buf, "splice_donor_variant&intron_variant"));

    /* truncation: NUL-terminated, never overruns, returns the FULL needed length */
    {
        char small[8];
        size_t need = duckvep_so_render(DUCKVEP_SO(DUCKVEP_SO_INTRON), '&', small, sizeof small);
        ASSERT_EQ((size_t)14u, need);          /* full "intron_variant" length */
        ASSERT_EQ(0, strcmp(small, "intron_"));/* 7 chars + NUL fit in 8 */
    }

    ASSERT_EQ(3u, duckvep_so_rank(DUCKVEP_SO_SPLICE_DONOR));
    ASSERT_EQ(28u, duckvep_so_rank(DUCKVEP_SO_INTRON));
    ASSERT_EQ(2u, duckvep_so_tier(DUCKVEP_SO_MATURE_MIRNA));
    ASSERT_EQ(3u, duckvep_so_tier(DUCKVEP_SO_INTRON));

    ASSERT_EQ(0, strcmp("HIGH", duckvep_impact_name(DUCKVEP_IMPACT_HIGH)));
    ASSERT_EQ(0, strcmp("MODERATE", duckvep_impact_name(DUCKVEP_IMPACT_MODERATE)));
    ASSERT_EQ(0, strcmp("LOW", duckvep_impact_name(DUCKVEP_IMPACT_LOW)));
    ASSERT_EQ(0, strcmp("MODIFIER", duckvep_impact_name(DUCKVEP_IMPACT_MODIFIER)));
    PASS();
}

/* ===================================================================== *
 * Coding coordinate projection (genomic -> cDNA -> CDS -> peptide).
 *
 * Oracle: a brute-force transcript-order base walk over small random exon
 * models. This deliberately does NOT use the model's cdna_start1/cdna_end1
 * arrays, while the implementation does, so it catches adapter/order/strand
 * mistakes instead of restating the implementation. The phase anchor covers the
 * positive Ensembl start-exon phase convention audited from VEP/Ensembl source:
 * a positive phase shifts CDS numbering by that many bases.
 * ===================================================================== */

#define KPROP_MAX_PROJ_EXONS 5u

struct kprop_proj_scene {
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t       ex;
    uint16_t chrom; uint32_t tstart; uint32_t tend; int8_t strand; uint64_t flags;
    uint32_t exoff; uint16_t excnt; uint32_t cds_s; uint32_t cds_e;
    uint32_t es[KPROP_MAX_PROJ_EXONS];
    uint32_t ee[KPROP_MAX_PROJ_EXONS];
    uint32_t cs[KPROP_MAX_PROJ_EXONS];
    uint32_t ce[KPROP_MAX_PROJ_EXONS];
    int8_t phase[KPROP_MAX_PROJ_EXONS];
    int8_t end_phase[KPROP_MAX_PROJ_EXONS];
};

static uint32_t proj_brute_cdna_len(const struct kprop_proj_scene *s) {
    uint32_t n = 0u;
    uint32_t i;
    for (i = 0u; i < (uint32_t)s->excnt; i++) n += (uint32_t)(s->ee[i] - s->es[i] + 1u);
    return n;
}

static int proj_brute_genomic_to_cdna(const struct kprop_proj_scene *s, uint32_t genomic,
                                      uint32_t *cdna, uint32_t *exon_idx) {
    uint32_t c = 1u;
    uint32_t i;
    int fwd = s->strand >= 0;
    for (i = 0u; i < (uint32_t)s->excnt; i++) {
        uint32_t len = (uint32_t)(s->ee[i] - s->es[i] + 1u);
        uint32_t j;
        for (j = 0u; j < len; j++, c++) {
            uint32_t g = fwd ? (uint32_t)(s->es[i] + j) : (uint32_t)(s->ee[i] - j);
            if (g == genomic) {
                *cdna = c;
                if (exon_idx != NULL) *exon_idx = i;
                return 1;
            }
        }
    }
    return 0;
}

static int proj_brute_cdna_to_genomic(const struct kprop_proj_scene *s, uint32_t cdna,
                                      uint32_t *genomic, uint32_t *exon_idx) {
    uint32_t c = 1u;
    uint32_t i;
    int fwd = s->strand >= 0;
    for (i = 0u; i < (uint32_t)s->excnt; i++) {
        uint32_t len = (uint32_t)(s->ee[i] - s->es[i] + 1u);
        uint32_t j;
        for (j = 0u; j < len; j++, c++) {
            if (c == cdna) {
                *genomic = fwd ? (uint32_t)(s->es[i] + j) : (uint32_t)(s->ee[i] - j);
                if (exon_idx != NULL) *exon_idx = i;
                return 1;
            }
        }
    }
    return 0;
}

static int proj_brute_coding(const struct kprop_proj_scene *s, uint32_t genomic,
                             duckvep_coding_projection_t *out) {
    duckvep_coding_projection_t r;
    uint32_t cdna = 0u, exon_idx = 0u;
    uint32_t start_cdna = 0u, end_cdna = 0u, start_exon = 0u, dummy = 0u;
    uint32_t cds;
    uint32_t start_g;
    uint32_t end_g;
    uint8_t phase;
    memset(&r, 0, sizeof r);
    if (s->cds_s == 0u || s->cds_e == 0u) return 0;
    if (!proj_brute_genomic_to_cdna(s, genomic, &cdna, &exon_idx)) return 0;
    start_g = s->strand >= 0 ? s->cds_s : s->cds_e;
    end_g   = s->strand >= 0 ? s->cds_e : s->cds_s;
    if (!proj_brute_genomic_to_cdna(s, start_g, &start_cdna, &start_exon)) return 0;
    if (!proj_brute_genomic_to_cdna(s, end_g, &end_cdna, &dummy)) return 0;
    if (end_cdna < start_cdna || cdna < start_cdna || cdna > end_cdna) return 0;
    phase = s->phase[start_exon] > 0 ? (uint8_t)s->phase[start_exon] : 0u;
    cds = (uint32_t)(cdna - start_cdna + 1u + (uint32_t)phase);
    r.cdna_pos = cdna;
    r.cds_pos = cds;
    r.protein_pos = (uint32_t)((cds - 1u) / 3u + 1u);
    r.codon_offset = (uint8_t)((cds - 1u) % 3u);
    r.codon_start_cds = (uint32_t)(cds - (uint32_t)r.codon_offset);
    r.exon_idx = exon_idx;
    r.phase_offset = phase;
    *out = r;
    return 1;
}

static void kprop_proj_scene_finish(struct kprop_proj_scene *s) {
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e;
    s->tx.transcript_count = 1u;
    s->ex.start1 = s->es; s->ex.end1 = s->ee; s->ex.cdna_start1 = s->cs; s->ex.cdna_end1 = s->ce;
    s->ex.phase = s->phase; s->ex.end_phase = s->end_phase; s->ex.exon_count = (size_t)s->excnt;
}

static enum theft_alloc_res kprop_proj_alloc(struct theft *t, void *env, void **instance) {
    struct kprop_proj_scene *s = (struct kprop_proj_scene *)calloc(1u, sizeof *s);
    uint32_t nex = (uint32_t)kprop_bounded(t, KPROP_MAX_PROJ_EXONS) + 1u;
    uint32_t asc_s[KPROP_MAX_PROJ_EXONS] = {0u};
    uint32_t asc_e[KPROP_MAX_PROJ_EXONS] = {0u};
    uint32_t base = (uint32_t)kprop_bounded(t, 10000u) + 100u;
    uint32_t cursor = base;
    uint32_t cdna = 1u;
    uint32_t i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;

    for (i = 0u; i < nex; i++) {
        uint32_t len = (uint32_t)kprop_bounded(t, 12u) + 3u;
        uint32_t gap = (uint32_t)kprop_bounded(t, 12u) + 1u;
        asc_s[i] = cursor;
        asc_e[i] = cursor + len - 1u;
        cursor = asc_e[i] + gap + 1u;
    }

    s->chrom = 0u;
    s->tstart = asc_s[0];
    s->tend = asc_e[nex - 1u];
    s->strand = (kprop_bounded(t, 2u) == 0u) ? (int8_t)1 : (int8_t)-1;
    s->exoff = 0u;
    s->excnt = (uint16_t)nex;
    for (i = 0u; i < nex; i++) {
        uint32_t src = s->strand >= 0 ? i : (uint32_t)(nex - 1u - i);
        uint32_t len = (uint32_t)(asc_e[src] - asc_s[src] + 1u);
        s->es[i] = asc_s[src];
        s->ee[i] = asc_e[src];
        s->cs[i] = cdna;
        s->ce[i] = cdna + len - 1u;
        s->phase[i] = (int8_t)-1;
        s->end_phase[i] = (int8_t)-1;
        cdna += len;
    }
    kprop_proj_scene_finish(s);

    if (kprop_bounded(t, 4u) != 0u) { /* mostly coding, but keep non-coding cases */
        uint32_t total = proj_brute_cdna_len(s);
        uint32_t c_start = (uint32_t)kprop_bounded(t, total) + 1u;
        uint32_t c_len = (uint32_t)kprop_bounded(t, (uint64_t)(total - c_start + 1u)) + 1u;
        uint32_t c_end = c_start + c_len - 1u;
        uint32_t g_start = 0u, g_end = 0u, start_exon = 0u;
        (void)proj_brute_cdna_to_genomic(s, c_start, &g_start, &start_exon);
        (void)proj_brute_cdna_to_genomic(s, c_end, &g_end, NULL);
        s->cds_s = g_start < g_end ? g_start : g_end;
        s->cds_e = g_start < g_end ? g_end : g_start;
        s->phase[start_exon] = (int8_t)kprop_bounded(t, 3u); /* Ensembl phase 0/1/2 */
    }

    *instance = s;
    return THEFT_ALLOC_OK;
}

static void kprop_proj_free(void *instance, void *env) { (void)env; free(instance); }
static struct theft_type_info kprop_proj_info = {
    .alloc = kprop_proj_alloc,
    .free  = kprop_proj_free,
};

static enum theft_trial_res prop_projection_matches_bruteforce(struct theft *t, void *arg1) {
    const struct kprop_proj_scene *s = (const struct kprop_proj_scene *)arg1;
    uint32_t total = proj_brute_cdna_len(s);
    uint32_t cdna;
    uint32_t outside = s->tstart > 1u ? (uint32_t)(s->tstart - 1u) : (uint32_t)(s->tend + 1u);
    (void)t;

    for (cdna = 1u; cdna <= total; cdna++) {
        uint32_t want_g = 0u, want_exon = 0u;
        uint32_t got_g = 0u, got_exon = 0u;
        uint32_t got_cdna = 0u;
        duckvep_coding_projection_t want_cp, got_cp;
        int want_coding;
        int got_coding;
        if (!proj_brute_cdna_to_genomic(s, cdna, &want_g, &want_exon)) return THEFT_TRIAL_FAIL;
        if (!duckvep_project_cdna_to_genomic(&s->tx, &s->ex, 0u, cdna, &got_g, &got_exon)) return THEFT_TRIAL_FAIL;
        if (got_g != want_g || got_exon != want_exon) return THEFT_TRIAL_FAIL;
        if (!duckvep_project_genomic_to_cdna(&s->tx, &s->ex, 0u, want_g, &got_cdna, &got_exon)) return THEFT_TRIAL_FAIL;
        if (got_cdna != cdna || got_exon != want_exon) return THEFT_TRIAL_FAIL;

        want_coding = proj_brute_coding(s, want_g, &want_cp);
        got_coding = duckvep_project_coding_base(&s->tx, &s->ex, 0u, want_g, &got_cp);
        if (want_coding != got_coding) return THEFT_TRIAL_FAIL;
        if (want_coding) {
            if (got_cp.cdna_pos != want_cp.cdna_pos || got_cp.cds_pos != want_cp.cds_pos ||
                got_cp.protein_pos != want_cp.protein_pos ||
                got_cp.codon_start_cds != want_cp.codon_start_cds ||
                got_cp.exon_idx != want_cp.exon_idx ||
                got_cp.codon_offset != want_cp.codon_offset ||
                got_cp.phase_offset != want_cp.phase_offset) {
                return THEFT_TRIAL_FAIL;
            }
        }
    }

    if (duckvep_project_genomic_to_cdna(&s->tx, &s->ex, 0u, outside, &cdna, NULL)) return THEFT_TRIAL_FAIL;
    if (duckvep_project_coding_base(&s->tx, &s->ex, 0u, outside, &(duckvep_coding_projection_t){0})) return THEFT_TRIAL_FAIL;
    return THEFT_TRIAL_PASS;
}

TEST projection_known_forward_reverse_and_phase(void) {
    duckvep_coding_projection_t p;
    uint32_t cdna = 0u, genomic = 0u, exon_idx = 0u;

    /* Forward: exons [100,109]=>cDNA 1..10, [200,209]=>11..20; CDS 103..205. */
    {
        struct kprop_proj_scene s;
        memset(&s, 0, sizeof s);
        s.chrom = 0u; s.tstart = 100u; s.tend = 209u; s.strand = (int8_t)1;
        s.excnt = 2u; s.cds_s = 103u; s.cds_e = 205u;
        s.es[0] = 100u; s.ee[0] = 109u; s.cs[0] = 1u;  s.ce[0] = 10u; s.phase[0] = 0;
        s.es[1] = 200u; s.ee[1] = 209u; s.cs[1] = 11u; s.ce[1] = 20u; s.phase[1] = 0;
        kprop_proj_scene_finish(&s);
        ASSERT_EQ(1, duckvep_project_genomic_to_cdna(&s.tx, &s.ex, 0u, 100u, &cdna, &exon_idx));
        ASSERT_EQ(1u, cdna); ASSERT_EQ(0u, exon_idx);
        ASSERT_EQ(1, duckvep_project_cdna_to_genomic(&s.tx, &s.ex, 0u, 11u, &genomic, &exon_idx));
        ASSERT_EQ(200u, genomic); ASSERT_EQ(1u, exon_idx);
        ASSERT_EQ(1, duckvep_project_coding_base(&s.tx, &s.ex, 0u, 105u, &p));
        ASSERT_EQ(6u, p.cdna_pos); ASSERT_EQ(3u, p.cds_pos); ASSERT_EQ(1u, p.protein_pos);
        ASSERT_EQ(1u, p.codon_start_cds); ASSERT_EQ(2u, (uint32_t)p.codon_offset);
        ASSERT_EQ(1, duckvep_project_coding_base(&s.tx, &s.ex, 0u, 200u, &p));
        ASSERT_EQ(11u, p.cdna_pos); ASSERT_EQ(8u, p.cds_pos); ASSERT_EQ(3u, p.protein_pos);
        ASSERT_EQ(0, duckvep_project_genomic_to_cdna(&s.tx, &s.ex, 0u, 150u, &cdna, NULL));
    }

    /* Reverse: transcript order is high-to-low. CDS genomic 104..205 starts at 205. */
    {
        struct kprop_proj_scene s;
        memset(&s, 0, sizeof s);
        s.chrom = 0u; s.tstart = 100u; s.tend = 209u; s.strand = (int8_t)-1;
        s.excnt = 2u; s.cds_s = 104u; s.cds_e = 205u;
        s.es[0] = 200u; s.ee[0] = 209u; s.cs[0] = 1u;  s.ce[0] = 10u; s.phase[0] = 0;
        s.es[1] = 100u; s.ee[1] = 109u; s.cs[1] = 11u; s.ce[1] = 20u; s.phase[1] = 0;
        kprop_proj_scene_finish(&s);
        ASSERT_EQ(1, duckvep_project_genomic_to_cdna(&s.tx, &s.ex, 0u, 205u, &cdna, &exon_idx));
        ASSERT_EQ(5u, cdna); ASSERT_EQ(0u, exon_idx);
        ASSERT_EQ(1, duckvep_project_coding_base(&s.tx, &s.ex, 0u, 203u, &p));
        ASSERT_EQ(7u, p.cdna_pos); ASSERT_EQ(3u, p.cds_pos); ASSERT_EQ(1u, p.protein_pos);
        ASSERT_EQ(1, duckvep_project_coding_base(&s.tx, &s.ex, 0u, 104u, &p));
        ASSERT_EQ(16u, p.cdna_pos); ASSERT_EQ(12u, p.cds_pos); ASSERT_EQ(4u, p.protein_pos);
    }

    /* Positive Ensembl phase convention: phase 2 means first coding base is CDS position 3. */
    {
        struct kprop_proj_scene s;
        memset(&s, 0, sizeof s);
        s.chrom = 0u; s.tstart = 100u; s.tend = 110u; s.strand = (int8_t)1;
        s.excnt = 1u; s.cds_s = 100u; s.cds_e = 110u;
        s.es[0] = 100u; s.ee[0] = 110u; s.cs[0] = 1u; s.ce[0] = 11u; s.phase[0] = 2;
        kprop_proj_scene_finish(&s);
        ASSERT_EQ(1, duckvep_project_coding_base(&s.tx, &s.ex, 0u, 100u, &p));
        ASSERT_EQ(1u, p.cdna_pos); ASSERT_EQ(3u, p.cds_pos); ASSERT_EQ(1u, p.protein_pos);
        ASSERT_EQ(1u, p.codon_start_cds); ASSERT_EQ(2u, (uint32_t)p.codon_offset);
        ASSERT_EQ(2u, (uint32_t)p.phase_offset);
        ASSERT_EQ(1, duckvep_project_coding_base(&s.tx, &s.ex, 0u, 101u, &p));
        ASSERT_EQ(4u, p.cds_pos); ASSERT_EQ(2u, p.protein_pos); ASSERT_EQ(4u, p.codon_start_cds);
    }

    /* Non-zero exon_offset / tx_idx>0: returned exon_idx is absolute in exon model. */
    {
        static uint16_t chrom[2] = {0u, 0u};
        static uint32_t tstart[2] = {1u, 100u};
        static uint32_t tend[2] = {10u, 209u};
        static int8_t strand[2] = {(int8_t)1, (int8_t)1};
        static uint64_t flags[2] = {0u, 0u};
        static uint32_t exoff[2] = {0u, 1u};
        static uint16_t excnt[2] = {1u, 2u};
        static uint32_t cds_s[2] = {0u, 103u};
        static uint32_t cds_e[2] = {0u, 205u};
        static uint32_t es[3] = {1u, 100u, 200u};
        static uint32_t ee[3] = {10u, 109u, 209u};
        static uint32_t cs[3] = {1u, 1u, 11u};
        static uint32_t ce[3] = {10u, 10u, 20u};
        static int8_t phase[3] = {0, 0, 0};
        duckvep_transcript_model_t tx;
        duckvep_exon_model_t ex;
        memset(&tx, 0, sizeof tx); memset(&ex, 0, sizeof ex);
        tx.chrom_id = chrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = strand;
        tx.flags = flags; tx.exon_offset = exoff; tx.exon_count = excnt;
        tx.cds_start1 = cds_s; tx.cds_end1 = cds_e; tx.transcript_count = 2u;
        ex.start1 = es; ex.end1 = ee; ex.cdna_start1 = cs; ex.cdna_end1 = ce;
        ex.phase = phase; ex.end_phase = phase; ex.exon_count = 3u;
        ASSERT_EQ(1, duckvep_project_genomic_to_cdna(&tx, &ex, 1u, 200u, &cdna, &exon_idx));
        ASSERT_EQ(11u, cdna); ASSERT_EQ(2u, exon_idx);
        ASSERT_EQ(1, duckvep_project_coding_base(&tx, &ex, 1u, 200u, &p));
        ASSERT_EQ(8u, p.cds_pos); ASSERT_EQ(2u, p.exon_idx);
    }
    PASS();
}

TEST projection_matches_bruteforce_for_any_small_transcript(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "coordinate projection == brute-force transcript-order base walk";
    cfg.prop1 = prop_projection_matches_bruteforce;
    cfg.type_info[0] = &kprop_proj_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

/* ===================================================================== *
 * Codon translation + coding-change classification, incl. the mitochondrial
 * table. The canonical genetic code is the independent oracle.
 * ===================================================================== */

#define STD  DUCKVEP_CODON_TABLE_STANDARD
#define MITO DUCKVEP_CODON_TABLE_VERT_MITO

TEST codon_translate_known(void) {
    /* standard (transl_table=1) */
    ASSERT_EQ('M', duckvep_translate_codon("ATG", STD));
    ASSERT_EQ('W', duckvep_translate_codon("TGG", STD));
    ASSERT_EQ('*', duckvep_translate_codon("TAA", STD));
    ASSERT_EQ('*', duckvep_translate_codon("TAG", STD));
    ASSERT_EQ('*', duckvep_translate_codon("TGA", STD));
    ASSERT_EQ('I', duckvep_translate_codon("ATA", STD));
    ASSERT_EQ('R', duckvep_translate_codon("AGA", STD));
    ASSERT_EQ('R', duckvep_translate_codon("AGG", STD));
    ASSERT_EQ('G', duckvep_translate_codon("ggg", STD)); /* case-insensitive */
    ASSERT_EQ('X', duckvep_translate_codon("ANG", STD)); /* invalid base */
    /* vertebrate mitochondrial (transl_table=2) — the four edits */
    ASSERT_EQ('W', duckvep_translate_codon("TGA", MITO));
    ASSERT_EQ('M', duckvep_translate_codon("ATA", MITO));
    ASSERT_EQ('*', duckvep_translate_codon("AGA", MITO));
    ASSERT_EQ('*', duckvep_translate_codon("AGG", MITO));
    ASSERT_EQ('M', duckvep_translate_codon("ATG", MITO)); /* start unchanged */
    PASS();
}

TEST codon_change_known(void) {
    duckvep_codon_result_t r;
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_SYNONYMOUS, duckvep_codon_change("AAA", "AAA", STD).change);
    r = duckvep_codon_change("AAA", "AAC", STD); /* K -> N */
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_MISSENSE, r.change);
    ASSERT_EQ('K', r.aa_ref); ASSERT_EQ('N', r.aa_alt);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_STOP_GAINED, duckvep_codon_change("TAC", "TAA", STD).change); /* Y->* */
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_STOP_LOST,   duckvep_codon_change("TAA", "TAC", STD).change); /* *->Y */
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_INVALID,     duckvep_codon_change("ANG", "AAA", STD).change);
    /* MT matters: TGG->TGA is stop_gained nuclear but synonymous in mito (TGA=W). */
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_STOP_GAINED, duckvep_codon_change("TGG", "TGA", STD).change);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_SYNONYMOUS,  duckvep_codon_change("TGG", "TGA", MITO).change);
    /* CGA->AGA: synonymous nuclear (R->R) but stop_gained in mito (AGA=*). */
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_SYNONYMOUS,  duckvep_codon_change("CGA", "AGA", STD).change);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_STOP_GAINED, duckvep_codon_change("CGA", "AGA", MITO).change);
    PASS();
}

struct kprop_codon { char ref[4]; char alt[4]; duckvep_codon_table_t table; };

static enum theft_alloc_res kprop_codon_alloc(struct theft *t, void *env, void **instance) {
    static const char bases5[5] = {'A', 'C', 'G', 'T', 'N'}; /* N -> invalid path */
    struct kprop_codon *c = (struct kprop_codon *)calloc(1u, sizeof *c);
    int i;
    (void)env;
    if (c == NULL) return THEFT_ALLOC_ERROR;
    for (i = 0; i < 3; i++) {
        c->ref[i] = bases5[kprop_bounded(t, 5u)];
        c->alt[i] = bases5[kprop_bounded(t, 5u)];
    }
    c->table = (kprop_bounded(t, 2u) == 0u) ? STD : MITO;
    *instance = c;
    return THEFT_ALLOC_OK;
}
static void kprop_codon_free(void *instance, void *env) { (void)env; free(instance); }
static struct theft_type_info kprop_codon_info = {
    .alloc = kprop_codon_alloc,
    .free  = kprop_codon_free,
};

/* The classification must be consistent with translation: exactly one change
 * bit, and each bit's defining relationship between aa_ref/aa_alt holds. */
static enum theft_trial_res prop_codon_change_consistent(struct theft *t, void *arg1) {
    const struct kprop_codon *c = (const struct kprop_codon *)arg1;
    char ar = duckvep_translate_codon(c->ref, c->table);
    char aa = duckvep_translate_codon(c->alt, c->table);
    duckvep_codon_result_t r = duckvep_codon_change(c->ref, c->alt, c->table);
    (void)t;
    if (r.aa_ref != ar || r.aa_alt != aa) return THEFT_TRIAL_FAIL;
    if (ar == 'X' || aa == 'X') {
        return (r.change == (uint32_t)DUCKVEP_CODON_INVALID) ? THEFT_TRIAL_PASS : THEFT_TRIAL_FAIL;
    }
    if (popcount_u32(r.change) != 1) return THEFT_TRIAL_FAIL;
    switch (r.change) {
        case DUCKVEP_CODON_SYNONYMOUS:  if (ar != aa) return THEFT_TRIAL_FAIL; break;
        case DUCKVEP_CODON_MISSENSE:    if (ar == aa || ar == '*' || aa == '*') return THEFT_TRIAL_FAIL; break;
        case DUCKVEP_CODON_STOP_GAINED: if (!(aa == '*' && ar != '*')) return THEFT_TRIAL_FAIL; break;
        case DUCKVEP_CODON_STOP_LOST:   if (!(ar == '*' && aa != '*')) return THEFT_TRIAL_FAIL; break;
        default: return THEFT_TRIAL_FAIL;
    }
    return THEFT_TRIAL_PASS;
}

TEST codon_change_consistent_with_translation(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "codon change classification consistent with translation";
    cfg.prop1 = prop_codon_change_consistent;
    cfg.type_info[0] = &kprop_codon_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

/* ===================================================================== *
 * Sequence-backed SNV codon edit application.
 *
 * Oracle: a codon-slice edit derived from the generated CDS position — copy the
 * reference codon from the generated CDS and replace exactly the projected base
 * after strand orientation — then classify via the already independently-tested
 * genetic-code table. This pins the hot-path contract between projection and the
 * codon classifier without claiming indel/haplotype support yet.
 * ===================================================================== */

static char coding_test_comp(char b) {
    switch (b) {
        case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
        default: return 'N';
    }
}
static char coding_test_genomic_from_tx(char tx_base, int8_t strand) {
    return strand < 0 ? coding_test_comp(tx_base) : tx_base;
}

TEST coding_snv_from_cds_known_cases(void) {
    duckvep_coding_projection_t p;
    duckvep_coding_snv_result_t r;

    memset(&p, 0, sizeof p);
    p.cds_pos = 4u; p.codon_start_cds = 4u; p.protein_pos = 2u; p.codon_offset = 0u;
    ASSERT_EQ(DUCKVEP_CODING_SNV_OK,
              duckvep_coding_snv_from_cds((const uint8_t *)"ATGGCTTGG", 9u, &p,
                                          'G', 'A', (int8_t)1, STD, &r));
    ASSERT_STR_EQ("GCT", r.ref_codon);
    ASSERT_STR_EQ("ACT", r.alt_codon);
    ASSERT_EQ('A', r.aa_ref); ASSERT_EQ('T', r.aa_alt);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_MISSENSE, r.change);
    ASSERT_EQ(4u, r.cds_pos); ASSERT_EQ(2u, r.protein_pos);

    /* Negative strand: genomic C->A orients to transcript G->T, TGG->TGT. */
    memset(&p, 0, sizeof p);
    p.cds_pos = 3u; p.codon_start_cds = 1u; p.protein_pos = 1u; p.codon_offset = 2u;
    ASSERT_EQ(DUCKVEP_CODING_SNV_OK,
              duckvep_coding_snv_from_cds((const uint8_t *)"TGG", 3u, &p,
                                          'C', 'A', (int8_t)-1, STD, &r));
    ASSERT_STR_EQ("TGG", r.ref_codon);
    ASSERT_STR_EQ("TGT", r.alt_codon);
    ASSERT_EQ('G', r.ref_base_tx); ASSERT_EQ('T', r.alt_base_tx);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_MISSENSE, r.change);

    /* Codon table still matters after sequence-backed edit: TGG->TGA. */
    ASSERT_EQ(DUCKVEP_CODING_SNV_OK,
              duckvep_coding_snv_from_cds((const uint8_t *)"TGG", 3u, &p,
                                          'G', 'A', (int8_t)1, STD, &r));
    ASSERT_STR_EQ("TGA", r.alt_codon);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_STOP_GAINED, r.change);
    ASSERT_EQ(DUCKVEP_CODING_SNV_OK,
              duckvep_coding_snv_from_cds((const uint8_t *)"TGG", 3u, &p,
                                          'G', 'A', (int8_t)1, MITO, &r));
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_SYNONYMOUS, r.change);

    ASSERT_EQ(DUCKVEP_CODING_SNV_REF_MISMATCH,
              duckvep_coding_snv_from_cds((const uint8_t *)"AAA", 3u, &p,
                                          'C', 'G', (int8_t)1, STD, &r));
    ASSERT_EQ('\0', r.ref_codon[0]);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_INVALID, r.change);
    ASSERT_EQ('X', r.aa_ref);

    p.codon_start_cds = 2u; p.cds_pos = 2u; p.codon_offset = 0u; p.protein_pos = 1u;
    ASSERT_EQ(DUCKVEP_CODING_SNV_INVALID_ARG,
              duckvep_coding_snv_from_cds((const uint8_t *)"AAAAAA", 6u, &p,
                                          'A', 'G', (int8_t)1, STD, &r));
    ASSERT_EQ('\0', r.ref_codon[0]);

    p.codon_start_cds = 4u; p.cds_pos = 4u; p.codon_offset = 0u; p.protein_pos = 2u;
    ASSERT_EQ(DUCKVEP_CODING_SNV_CODON_OUT_OF_RANGE,
              duckvep_coding_snv_from_cds((const uint8_t *)"AAAAA", 5u, &p,
                                          'A', 'G', (int8_t)1, STD, &r));
    ASSERT_EQ('\0', r.ref_codon[0]);

    p.codon_start_cds = 1u; p.cds_pos = 3u; p.codon_offset = 2u; p.protein_pos = 1u;
    ASSERT_EQ(DUCKVEP_CODING_SNV_INVALID_BASE,
              duckvep_coding_snv_from_cds((const uint8_t *)"AAA", 3u, &p,
                                          'A', 'N', (int8_t)1, STD, &r));
    ASSERT_EQ('\0', r.ref_codon[0]);

    ASSERT_EQ(DUCKVEP_CODING_SNV_INVALID_ARG,
              duckvep_coding_snv_from_cds((const uint8_t *)"AAA", 3u, &p,
                                          'A', 'G', (int8_t)0, STD, &r));
    p.protein_pos = 99u;
    ASSERT_EQ(DUCKVEP_CODING_SNV_INVALID_ARG,
              duckvep_coding_snv_from_cds((const uint8_t *)"AAA", 3u, &p,
                                          'A', 'G', (int8_t)1, STD, &r));
    p.protein_pos = 1u;

    /* N in non-edited phase-padding positions is allowed but classifies invalid. */
    ASSERT_EQ(DUCKVEP_CODING_SNV_OK,
              duckvep_coding_snv_from_cds((const uint8_t *)"NNA", 3u, &p,
                                          'A', 'G', (int8_t)1, STD, &r));
    ASSERT_STR_EQ("NNA", r.ref_codon);
    ASSERT_STR_EQ("NNG", r.alt_codon);
    ASSERT_EQ((uint32_t)DUCKVEP_CODON_INVALID, r.change);
    PASS();
}

#define KPROP_MAX_CODING_CDS 90u

struct kprop_coding_snv {
    uint8_t seq[KPROP_MAX_CODING_CDS];
    size_t len;
    duckvep_coding_projection_t p;
    char genomic_ref, genomic_alt;
    char ref_tx, alt_tx;
    char ref_codon[4], alt_codon[4];
    duckvep_codon_table_t table;
    int8_t strand;
};

static enum theft_alloc_res kprop_coding_snv_alloc(struct theft *t, void *env, void **instance) {
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    struct kprop_coding_snv *s = (struct kprop_coding_snv *)calloc(1u, sizeof *s);
    uint32_t codons = (uint32_t)kprop_bounded(t, 20u) + 1u;
    uint32_t codon_idx = (uint32_t)kprop_bounded(t, codons);
    uint32_t offset = (uint32_t)kprop_bounded(t, 3u);
    size_t i;
    uint32_t alt_pick;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->len = (size_t)codons * 3u;
    if (s->len > KPROP_MAX_CODING_CDS) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < s->len; i++) s->seq[i] = (uint8_t)bases[kprop_bounded(t, 4u)];
    s->strand = (kprop_bounded(t, 2u) == 0u) ? (int8_t)1 : (int8_t)-1;
    s->table = (kprop_bounded(t, 2u) == 0u) ? STD : MITO;
    s->p.codon_start_cds = codon_idx * 3u + 1u;
    s->p.codon_offset = (uint8_t)offset;
    s->p.cds_pos = s->p.codon_start_cds + offset;
    s->p.protein_pos = codon_idx + 1u;
    s->ref_tx = (char)s->seq[s->p.cds_pos - 1u];
    alt_pick = (uint32_t)kprop_bounded(t, 3u) + 1u;
    s->alt_tx = bases[((s->ref_tx == 'A' ? 0u : s->ref_tx == 'C' ? 1u : s->ref_tx == 'G' ? 2u : 3u) + alt_pick) % 4u];
    s->genomic_ref = coding_test_genomic_from_tx(s->ref_tx, s->strand);
    s->genomic_alt = coding_test_genomic_from_tx(s->alt_tx, s->strand);
    for (i = 0u; i < 3u; i++) {
        s->ref_codon[i] = (char)s->seq[s->p.codon_start_cds - 1u + (uint32_t)i];
        s->alt_codon[i] = s->ref_codon[i];
    }
    s->ref_codon[3] = '\0'; s->alt_codon[3] = '\0';
    s->alt_codon[offset] = s->alt_tx;
    *instance = s;
    return THEFT_ALLOC_OK;
}
static void kprop_coding_snv_free(void *instance, void *env) { (void)env; free(instance); }
static struct theft_type_info kprop_coding_snv_info = {
    .alloc = kprop_coding_snv_alloc,
    .free  = kprop_coding_snv_free,
};

static enum theft_trial_res prop_coding_snv_matches_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding_snv *s = (const struct kprop_coding_snv *)arg1;
    duckvep_coding_snv_result_t r;
    duckvep_codon_result_t cr;
    (void)t;
    if (duckvep_coding_snv_from_cds(s->seq, s->len, &s->p, s->genomic_ref, s->genomic_alt,
                                    s->strand, s->table, &r) != DUCKVEP_CODING_SNV_OK) {
        return THEFT_TRIAL_FAIL;
    }
    cr = duckvep_codon_change(s->ref_codon, s->alt_codon, s->table);
    if (strcmp(r.ref_codon, s->ref_codon) != 0) return THEFT_TRIAL_FAIL;
    if (strcmp(r.alt_codon, s->alt_codon) != 0) return THEFT_TRIAL_FAIL;
    if (r.ref_base_tx != s->ref_tx || r.alt_base_tx != s->alt_tx) return THEFT_TRIAL_FAIL;
    if (r.aa_ref != cr.aa_ref || r.aa_alt != cr.aa_alt || r.change != cr.change) return THEFT_TRIAL_FAIL;
    if (r.cds_pos != s->p.cds_pos || r.codon_start_cds != s->p.codon_start_cds ||
        r.protein_pos != s->p.protein_pos || r.codon_offset != s->p.codon_offset) {
        return THEFT_TRIAL_FAIL;
    }
    return THEFT_TRIAL_PASS;
}

TEST coding_snv_from_cds_matches_oracle_for_any_valid_snv(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "sequence-backed SNV codon edit == codon-slice edit oracle";
    cfg.prop1 = prop_coding_snv_matches_oracle;
    cfg.type_info[0] = &kprop_coding_snv_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

/* ===================================================================== *
 * Haplosaurus-shaped multi-edit CDS mutation + translation.
 *
 * Oracle: rebuild the observed CDS from left to right against original
 * coordinates, appending untouched reference segments and independently
 * oriented alternate alleles. The implementation applies edits in-place from
 * right to left; the oracle never memmoves the working sequence. These tests
 * pin the source-audited Haplosaurus behavior without claiming sample/phase-set
 * aggregation or VEP differential parity yet.
 * ===================================================================== */

#define KPROP_HAPLO_MAX_EDITS 5u
#define KPROP_HAPLO_CDS_LEN   36u
#define KPROP_HAPLO_CAP       80u
#define KPROP_HAPLO_ALLELE_MAX 4u

static char haplo_test_oriented_base(const uint8_t *seq, uint32_t len, uint32_t idx,
                                     int reverse_complement) {
    char b;
    if (seq == NULL || idx >= len) return '\0';
    b = (char)seq[reverse_complement ? (len - 1u - idx) : idx];
    if (b >= 'a' && b <= 'z') b = (char)(b - ('a' - 'A'));
    if (b == 'U') b = 'T';
    if (b != 'A' && b != 'C' && b != 'G' && b != 'T') return '\0';
    return reverse_complement ? coding_test_comp(b) : b;
}

static uint8_t haplo_test_variant_from_tx_base(char b, int reverse_complement) {
    return (uint8_t)(reverse_complement ? coding_test_comp(b) : b);
}

static int haplo_oracle_rebuild(const uint8_t *ref, size_t ref_len,
                                const duckvep_haplotype_edit_t *edits, size_t edit_count,
                                int8_t transcript_strand,
                                uint8_t *out, size_t out_cap, size_t *out_len,
                                int64_t *length_diff, uint32_t *flags) {
    size_t cursor = 0u;
    size_t n = 0u;
    size_t ei;
    int saw_frameshift = 0;
    int64_t diff_total = 0;
    uint32_t f = 0u;

    for (ei = edit_count; ei > 0u; ei--) {
        const duckvep_haplotype_edit_t *e = &edits[ei - 1u];
        size_t start0 = (size_t)e->cds_start - 1u;
        uint32_t j;
        int reverse = e->variant_strand != transcript_strand;
        int64_t d = (int64_t)e->alt_len - (int64_t)e->ref_len;
        if (e->cds_start == 0u || start0 < cursor || start0 > ref_len ||
            (size_t)e->ref_len > ref_len - start0) return 0;
        if (n + (start0 - cursor) > out_cap) return 0;
        memcpy(out + n, ref + cursor, start0 - cursor);
        n += start0 - cursor;
        for (j = 0u; j < e->ref_len; j++) {
            char expected = haplo_test_oriented_base(e->ref, e->ref_len, j, reverse);
            if (expected == '\0' || expected != (char)ref[start0 + (size_t)j]) return 0;
        }
        if (n + (size_t)e->alt_len > out_cap) return 0;
        for (j = 0u; j < e->alt_len; j++) {
            char b = haplo_test_oriented_base(e->alt, e->alt_len, j, reverse);
            if (b == '\0') return 0;
            out[n++] = (uint8_t)b;
        }
        cursor = start0 + (size_t)e->ref_len;
        diff_total += d;
        if (d != 0) f |= DUCKVEP_HAPLOTYPE_FLAG_INDEL;
        if ((d % 3) != 0) saw_frameshift = 1;
    }
    if (n + (ref_len - cursor) > out_cap) return 0;
    memcpy(out + n, ref + cursor, ref_len - cursor);
    n += ref_len - cursor;
    if (n < out_cap) out[n] = (uint8_t)'\0';
    if (saw_frameshift) {
        if ((diff_total % 3) == 0) f |= DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT;
        else f |= DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT;
    }
    *out_len = n;
    *length_diff = diff_total;
    *flags = f;
    return 1;
}

TEST haplotype_apply_and_translate_known_cases(void) {
    duckvep_haplotype_edit_t edits[2];
    duckvep_haplotype_result_t r;
    uint8_t cds[64];
    uint8_t protein[32];
    size_t cds_len, protein_len;

    memset(edits, 0, sizeof edits);
    edits[0].cds_start = 2u; edits[0].ref_len = 1u; edits[0].ref = (const uint8_t *)"T";
    edits[0].alt_len = 1u; edits[0].alt = (const uint8_t *)"G"; edits[0].variant_strand = (int8_t)1;
    edits[1].cds_start = 1u; edits[1].ref_len = 1u; edits[1].ref = (const uint8_t *)"A";
    edits[1].alt_len = 1u; edits[1].alt = (const uint8_t *)"T"; edits[1].variant_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATGAAATAA", 9u,
                                                edits, 2u, (int8_t)1, cds, sizeof cds, &cds_len, &r));
    ASSERT_EQ(9u, cds_len);
    ASSERT_STR_EQ("TGGAAATAA", (const char *)cds);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_translate_cds(cds, cds_len, STD, protein, sizeof protein, &protein_len, &r));
    ASSERT_STR_EQ("WK*", (const char *)protein);
    ASSERT_EQ(3u, protein_len);

    /* Negative transcript strand: genomic C>A orients to transcript G>T. */
    memset(edits, 0, sizeof edits);
    edits[0].cds_start = 4u; edits[0].ref_len = 1u; edits[0].ref = (const uint8_t *)"C";
    edits[0].alt_len = 1u; edits[0].alt = (const uint8_t *)"A"; edits[0].variant_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATGGCTTAA", 9u,
                                                edits, 1u, (int8_t)-1, cds, sizeof cds, &cds_len, &r));
    ASSERT_STR_EQ("ATGTCTTAA", (const char *)cds);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_translate_cds(cds, cds_len, STD, protein, sizeof protein, &protein_len, &r));
    ASSERT_STR_EQ("MS*", (const char *)protein);

    /* Two non-triplet edits with net length diff 0 are a resolved frameshift. */
    memset(edits, 0, sizeof edits);
    edits[0].cds_start = 10u; edits[0].ref_len = 1u; edits[0].ref = (const uint8_t *)"G";
    edits[0].alt_len = 0u; edits[0].alt = NULL; edits[0].variant_strand = (int8_t)1;
    edits[1].cds_start = 4u; edits[1].ref_len = 0u; edits[1].ref = NULL;
    edits[1].alt_len = 1u; edits[1].alt = (const uint8_t *)"C"; edits[1].variant_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATGAAACCCGGGTAA", 15u,
                                                edits, 2u, (int8_t)1, cds, sizeof cds, &cds_len, &r));
    ASSERT_STR_EQ("ATGCAAACCCGGTAA", (const char *)cds);
    ASSERT_EQ(0, r.length_diff);
    ASSERT((r.flags & DUCKVEP_HAPLOTYPE_FLAG_INDEL) != 0u);
    ASSERT((r.flags & DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT) != 0u);
    ASSERT((r.flags & DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT) == 0u);

    /* Protein translation truncates after the first stop and is standalone: it
     * must not read pre-initialized apply flags from result. */
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATGAAATAACCC", 12u,
                                                NULL, 0u, (int8_t)1, cds, sizeof cds, &cds_len, &r));
    {
        duckvep_haplotype_result_t fresh_result;
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                  duckvep_haplotype_translate_cds(cds, cds_len, STD, protein, sizeof protein,
                                                  &protein_len, &fresh_result));
        ASSERT_STR_EQ("MK*", (const char *)protein);
        ASSERT_EQ(12u, fresh_result.cds_len);
        ASSERT((fresh_result.flags & DUCKVEP_HAPLOTYPE_FLAG_STOP_TRUNCATED) != 0u);
    }

    /* Contract checks: descending order and reference validation are explicit. */
    edits[0].cds_start = 1u; edits[0].ref_len = 1u; edits[0].ref = (const uint8_t *)"A";
    edits[0].alt_len = 1u; edits[0].alt = (const uint8_t *)"C"; edits[0].variant_strand = (int8_t)1;
    edits[1].cds_start = 2u; edits[1].ref_len = 1u; edits[1].ref = (const uint8_t *)"T";
    edits[1].alt_len = 1u; edits[1].alt = (const uint8_t *)"G"; edits[1].variant_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_EDIT_ORDER,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATG", 3u, edits, 2u,
                                                (int8_t)1, cds, sizeof cds, &cds_len, &r));
    edits[0].cds_start = 1u; edits[0].ref = (const uint8_t *)"C";
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_REF_MISMATCH,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATG", 3u, edits, 1u,
                                                (int8_t)1, cds, sizeof cds, &cds_len, &r));
    edits[0].cds_start = 4u; edits[0].ref_len = 0u; edits[0].ref = NULL;
    edits[0].alt_len = 1u; edits[0].alt = (const uint8_t *)"A"; edits[0].variant_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATG", 3u, edits, 1u,
                                                (int8_t)1, cds, sizeof cds, &cds_len, &r));
    ASSERT_STR_EQ("ATGA", (const char *)cds);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATG", 3u, edits, 1u,
                                                (int8_t)1, cds, 3u, &cds_len, &r));
    edits[0].cds_start = 2u; edits[0].ref_len = 2u; edits[0].ref = (const uint8_t *)"TG";
    edits[0].alt_len = 1u; edits[0].alt = (const uint8_t *)"C"; edits[0].variant_strand = (int8_t)1;
    edits[1].cds_start = 3u; edits[1].ref_len = 1u; edits[1].ref = (const uint8_t *)"G";
    edits[1].alt_len = 1u; edits[1].alt = (const uint8_t *)"A"; edits[1].variant_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_EDIT_ORDER,
              duckvep_haplotype_apply_cds_edits((const uint8_t *)"ATG", 3u, edits, 2u,
                                                (int8_t)1, cds, sizeof cds, &cds_len, &r));
    PASS();
}

TEST haplotype_partition_known_cases(void) {
    duckvep_haplotype_edit_t edits[2];
    duckvep_haplotype_block_t blocks[2];
    size_t required = 0u;

    memset(edits, 0, sizeof edits);
    edits[0].cds_start = 1u; edits[0].ref_len = 1u; edits[0].ref = (const uint8_t *)"A";
    edits[0].alt_len = 1u; edits[0].alt = (const uint8_t *)"C"; edits[0].variant_strand = (int8_t)1;
    edits[1].cds_start = 3u; edits[1].ref_len = 1u; edits[1].ref = (const uint8_t *)"G";
    edits[1].alt_len = 1u; edits[1].alt = (const uint8_t *)"T"; edits[1].variant_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_partition(edits, 2u, blocks, 2u, &required));
    ASSERT_EQ(1u, required);
    ASSERT_EQ(0u, blocks[0].edit_begin);
    ASSERT_EQ(2u, blocks[0].edit_count);
    ASSERT_EQ(1u, blocks[0].cds_start);
    ASSERT_EQ(3u, blocks[0].cds_end);
    ASSERT_EQ(0, blocks[0].length_diff);
    ASSERT_EQ(0u, blocks[0].flags);

    edits[1].cds_start = 7u;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_partition(edits, 2u, blocks, 2u, &required));
    ASSERT_EQ(2u, required);
    ASSERT_EQ(1u, blocks[0].edit_count);
    ASSERT_EQ(1u, blocks[1].edit_count);

    /* A downstream deletion cannot be flushed while the insertion has left the
     * frame displaced.  Together they restore it. */
    memset(edits, 0, sizeof edits);
    edits[0].cds_start = 4u; edits[0].alt_len = 1u;
    edits[0].alt = (const uint8_t *)"A"; edits[0].variant_strand = (int8_t)1;
    edits[1].cds_start = 10u; edits[1].ref_len = 1u;
    edits[1].ref = (const uint8_t *)"G"; edits[1].variant_strand = (int8_t)1;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_partition(edits, 2u, blocks, 2u, &required));
    ASSERT_EQ(1u, required);
    ASSERT_EQ(0, blocks[0].length_diff);
    ASSERT((blocks[0].flags & DUCKVEP_HAPLOTYPE_FLAG_INDEL) != 0u);
    ASSERT((blocks[0].flags & DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT) != 0u);
    ASSERT((blocks[0].flags & DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT) == 0u);

    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_partition(edits, 1u, blocks, 2u, &required));
    ASSERT_EQ(1u, required);
    ASSERT_EQ(1, blocks[0].length_diff);
    ASSERT((blocks[0].flags & DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT) != 0u);
    ASSERT((blocks[0].flags & DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT) == 0u);

    /* The sizing pass publishes the required count but never a partial array. */
    edits[0].cds_start = 1u; edits[0].ref_len = 1u;
    edits[0].ref = (const uint8_t *)"A"; edits[0].alt_len = 1u;
    edits[0].alt = (const uint8_t *)"C";
    edits[1].cds_start = 7u; edits[1].ref_len = 1u;
    edits[1].ref = (const uint8_t *)"G"; edits[1].alt_len = 1u;
    edits[1].alt = (const uint8_t *)"T";
    memset(blocks, 0xa5, sizeof blocks);
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL,
              duckvep_haplotype_partition(edits, 2u, blocks, 1u, &required));
    ASSERT_EQ(2u, required);
    {
        const unsigned char *raw = (const unsigned char *)blocks;
        size_t i;
        for (i = 0u; i < sizeof blocks; i++) ASSERT_EQ(0xa5u, raw[i]);
    }

    edits[1].cds_start = 1u;
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_EDIT_ORDER,
              duckvep_haplotype_partition(edits, 2u, blocks, 2u, &required));
    edits[1].cds_start = 2u;
    edits[0].ref_len = 2u; edits[0].ref = (const uint8_t *)"AT";
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_EDIT_ORDER,
              duckvep_haplotype_partition(edits, 2u, blocks, 2u, &required));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_INVALID_ARG,
              duckvep_haplotype_partition(NULL, 1u, blocks, 2u, &required));
    PASS();
}

struct kprop_haplo_case {
    uint8_t ref[KPROP_HAPLO_CDS_LEN + 1u];
    duckvep_haplotype_edit_t edits[KPROP_HAPLO_MAX_EDITS];
    uint8_t ref_alleles[KPROP_HAPLO_MAX_EDITS][3u];
    uint8_t alt_alleles[KPROP_HAPLO_MAX_EDITS][KPROP_HAPLO_ALLELE_MAX];
    size_t edit_count;
    int8_t transcript_strand;
};

static enum theft_alloc_res kprop_haplo_alloc(struct theft *t, void *env, void **instance) {
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    struct kprop_haplo_case *c = (struct kprop_haplo_case *)calloc(1u, sizeof *c);
    uint32_t starts_tmp[KPROP_HAPLO_MAX_EDITS];
    uint32_t ref_len_tmp[KPROP_HAPLO_MAX_EDITS];
    uint32_t alt_len_tmp[KPROP_HAPLO_MAX_EDITS];
    uint32_t cursor = 1u;
    size_t target;
    size_t n = 0u;
    size_t i;
    (void)env;
    if (c == NULL) return THEFT_ALLOC_ERROR;
    for (i = 0u; i < KPROP_HAPLO_CDS_LEN; i++) c->ref[i] = (uint8_t)bases[kprop_bounded(t, 4u)];
    c->ref[KPROP_HAPLO_CDS_LEN] = (uint8_t)'\0';
    target = (size_t)kprop_bounded(t, (uint64_t)KPROP_HAPLO_MAX_EDITS + 1u);
    c->transcript_strand = kprop_bounded(t, 2u) == 0u ? (int8_t)1 : (int8_t)-1;

    /* Build a valid non-overlapping edit set in ascending original CDS order,
     * allowing adjacency, insertions at the CDS end, and variable gaps; reverse
     * it below because the kernel contract is descending application order. */
    while (n < target && cursor <= KPROP_HAPLO_CDS_LEN + 1u) {
        uint32_t gap = (uint32_t)kprop_bounded(t, 3u); /* 0..2: adjacent is common. */
        uint32_t start = cursor + gap;
        uint32_t max_ref;
        uint32_t ref_len;
        if (start > KPROP_HAPLO_CDS_LEN + 1u) break;
        max_ref = start <= KPROP_HAPLO_CDS_LEN ? KPROP_HAPLO_CDS_LEN - start + 1u : 0u;
        if (max_ref > 3u) max_ref = 3u;
        ref_len = max_ref == 0u ? 0u : (uint32_t)kprop_bounded(t, (uint64_t)max_ref + 1u);
        starts_tmp[n] = start;
        ref_len_tmp[n] = ref_len;
        alt_len_tmp[n] = (uint32_t)kprop_bounded(t, KPROP_HAPLO_ALLELE_MAX + 1u); /* 0..4 */
        n++;
        cursor = start + (ref_len > 0u ? ref_len : 1u) + (uint32_t)kprop_bounded(t, 2u);
    }
    c->edit_count = n;

    for (i = 0u; i < c->edit_count; i++) {
        size_t src = c->edit_count - 1u - i;
        duckvep_haplotype_edit_t *e = &c->edits[i];
        uint32_t ref_len = ref_len_tmp[src];
        uint32_t alt_len = alt_len_tmp[src];
        uint32_t j;
        int reverse;
        e->cds_start = starts_tmp[src];
        e->ref_len = ref_len;
        e->alt_len = alt_len;
        e->variant_strand = kprop_bounded(t, 2u) == 0u ? (int8_t)1 : (int8_t)-1;
        reverse = e->variant_strand != c->transcript_strand;
        for (j = 0u; j < ref_len; j++) {
            char tx = (char)c->ref[(size_t)e->cds_start - 1u + (size_t)j];
            c->ref_alleles[i][reverse ? (ref_len - 1u - j) : j] = haplo_test_variant_from_tx_base(tx, reverse);
        }
        for (j = 0u; j < alt_len; j++) {
            char tx = bases[kprop_bounded(t, 4u)];
            c->alt_alleles[i][reverse ? (alt_len - 1u - j) : j] = haplo_test_variant_from_tx_base(tx, reverse);
        }
        e->ref = ref_len ? c->ref_alleles[i] : NULL;
        e->alt = alt_len ? c->alt_alleles[i] : NULL;
    }
    *instance = c;
    return THEFT_ALLOC_OK;
}
static void kprop_haplo_free(void *instance, void *env) { (void)env; free(instance); }
static struct theft_type_info kprop_haplo_info = {
    .alloc = kprop_haplo_alloc,
    .free  = kprop_haplo_free,
};

static enum theft_trial_res prop_haplotype_apply_matches_rebuild_oracle(struct theft *t, void *arg1) {
    const struct kprop_haplo_case *c = (const struct kprop_haplo_case *)arg1;
    uint8_t got[KPROP_HAPLO_CAP];
    uint8_t want[KPROP_HAPLO_CAP];
    size_t got_len = 0u, want_len = 0u;
    int64_t want_diff = 0;
    uint32_t want_flags = 0u;
    duckvep_haplotype_result_t r;
    (void)t;
    if (!haplo_oracle_rebuild(c->ref, KPROP_HAPLO_CDS_LEN, c->edits, c->edit_count,
                              c->transcript_strand, want, sizeof want, &want_len,
                              &want_diff, &want_flags)) {
        return THEFT_TRIAL_ERROR;
    }
    if (duckvep_haplotype_apply_cds_edits(c->ref, KPROP_HAPLO_CDS_LEN, c->edits, c->edit_count,
                                          c->transcript_strand, got, sizeof got, &got_len, &r) !=
        DUCKVEP_HAPLOTYPE_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (got_len != want_len || memcmp(got, want, got_len) != 0) return THEFT_TRIAL_FAIL;
    if (r.cds_len != want_len || r.length_diff != want_diff || r.flags != want_flags ||
        r.applied_edits != c->edit_count) {
        return THEFT_TRIAL_FAIL;
    }
    return THEFT_TRIAL_PASS;
}

TEST haplotype_apply_matches_rebuild_oracle_for_any_valid_edit_set(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "multi-edit CDS haplotype apply == left-to-right rebuild oracle";
    cfg.prop1 = prop_haplotype_apply_matches_rebuild_oracle;
    cfg.type_info[0] = &kprop_haplo_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

static enum theft_trial_res prop_haplotype_partition_preserves_interactions(
    struct theft *t, void *arg1) {
    const struct kprop_haplo_case *c = (const struct kprop_haplo_case *)arg1;
    duckvep_haplotype_edit_t ascending[KPROP_HAPLO_MAX_EDITS];
    duckvep_haplotype_block_t blocks[KPROP_HAPLO_MAX_EDITS];
    size_t block_count = 0u;
    size_t covered = 0u;
    size_t bi;
    (void)t;

    for (bi = 0u; bi < c->edit_count; bi++) {
        ascending[bi] = c->edits[c->edit_count - 1u - bi];
    }
    if (duckvep_haplotype_partition(ascending, c->edit_count, blocks,
                                    KPROP_HAPLO_MAX_EDITS, &block_count) !=
        DUCKVEP_HAPLOTYPE_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (block_count > c->edit_count ||
        (c->edit_count == 0u && block_count != 0u)) {
        return THEFT_TRIAL_FAIL;
    }

    for (bi = 0u; bi < block_count; bi++) {
        const duckvep_haplotype_block_t *block = &blocks[bi];
        int64_t difference = 0;
        int64_t shift_before = 0;
        uint32_t expected_flags = 0u;
        uint32_t expected_end = 0u;
        int saw_frameshift = 0;
        size_t j;

        if (block->edit_begin != covered || block->edit_count == 0u ||
            block->edit_begin + block->edit_count > c->edit_count) {
            return THEFT_TRIAL_FAIL;
        }
        for (j = block->edit_begin;
             j < block->edit_begin + block->edit_count; j++) {
            const duckvep_haplotype_edit_t *edit = &ascending[j];
            int64_t edit_difference = (int64_t)edit->alt_len - (int64_t)edit->ref_len;
            uint32_t edit_end = edit->cds_start +
                (edit->ref_len == 0u ? 0u : edit->ref_len - 1u);
            int should_stay_open = 0;

            if (edit_end > expected_end) expected_end = edit_end;
            difference += edit_difference;
            if (edit_difference != 0) {
                expected_flags |= DUCKVEP_HAPLOTYPE_FLAG_INDEL;
                if (edit_difference % 3 != 0) saw_frameshift = 1;
            }
            if (j + 1u < c->edit_count) {
                int64_t current_start = (int64_t)edit->cds_start - 1 + shift_before;
                int64_t current_end = current_start +
                    (edit->alt_len == 0u ? 0 : (int64_t)edit->alt_len - 1);
                int64_t next_start = (int64_t)ascending[j + 1u].cds_start - 1 + difference;
                should_stay_open = difference % 3 != 0 ||
                    current_end / 3 == next_start / 3;
            }
            shift_before += edit_difference;

            if (j + 1u < block->edit_begin + block->edit_count) {
                if (!should_stay_open) return THEFT_TRIAL_FAIL;
            } else if (j + 1u < c->edit_count && should_stay_open) {
                return THEFT_TRIAL_FAIL;
            }
        }
        if (saw_frameshift) {
            expected_flags |= difference % 3 == 0
                ? DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT
                : DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT;
        }
        if (block->cds_start != ascending[block->edit_begin].cds_start ||
            block->cds_end != expected_end ||
            block->length_diff != difference ||
            block->flags != expected_flags) {
            return THEFT_TRIAL_FAIL;
        }
        covered += block->edit_count;
    }
    return covered == c->edit_count ? THEFT_TRIAL_PASS : THEFT_TRIAL_FAIL;
}

TEST haplotype_partition_preserves_interactions_for_any_valid_edit_set(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "haplotype blocks preserve every frame and same-codon interaction";
    cfg.prop1 = prop_haplotype_partition_preserves_interactions;
    cfg.type_info[0] = &kprop_haplo_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

/* ===================================================================== *
 * annotate_tile FUSION: the composed engine == the sum of its tested parts.
 *
 * annotate_tile fuses the sweep (candidate pairs) with the classifier (region
 * mask) and the SO mapping. The strongest hermetic check: for every generated
 * scene, the rows annotate_tile emits must equal the INDEPENDENT composition —
 * the same sweep set, each pair mapped to its structural SO terms by a from-spec
 * map written here (not the kernel's mapper). kprop_scene transcripts carry no
 * exons and no CDS, so they exercise the upstream/downstream/intron branches at
 * scale + the sweep<->builder<->filter wiring; the deterministic scene below
 * pins the exon/UTR/CDS/splice branches exactly.
 * ===================================================================== */

static enum theft_trial_res prop_annotate_matches_composition(struct theft *t, void *arg1) {
    const struct kprop_scene *s = (const struct kprop_scene *)arg1;
    duckvep_exon_model_t exons;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_options_init_t init;
    duckvep_error_t err;
    static duckvep_consequence_t rows[KPROP_MAX_PAIRS];
    static uint64_t sweep_buf[KPROP_MAX_PAIRS];
    duckvep_result_builder_t rb;
    uint32_t active[KPROP_MAX_TX];
    uint32_t candidates[KPROP_MAX_TX];
    struct pair_collector col;
    duckvep_status_t st = DUCKVEP_OK;
    enum theft_trial_res res = THEFT_TRIAL_PASS;
    size_t n_sweep, i;
    const uint32_t HALO = KPROP_SWEEP_HALO; /* fixed >0 so options.halo is unambiguous */
    (void)t;

    memset(&exons, 0, sizeof exons); /* no exons -> non_coding_transcript (placement) / up / down */
    memset(&err, 0, sizeof err);

    /* Scene transcripts are sorted by (chrom_id, start1) by construction, so the
     * model must open. */
    if (duckvep_model_open(&s->tx, &exons, NULL, &model, &err) != DUCKVEP_OK) {
        return THEFT_TRIAL_FAIL;
    }
    memset(&init, 0, sizeof init);
    /* up == down == halo: every sweep-admitted pair is within the directional
     * reach, so the filter is inert, and options.halo (now clamped to
     * max(halo,up,down)) stays exactly HALO. */
    init.upstream_dist = HALO;
    init.downstream_dist = HALO;
    init.halo = HALO;
    if (duckvep_options_open(&init, &opts, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }

    duckvep_result_builder_init(&rb, rows, KPROP_MAX_PAIRS);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }

    /* Independent composition: same sweep, then a from-spec structural map. */
    col.buf = sweep_buf; col.n = 0u; col.cap = KPROP_MAX_PAIRS;
    n_sweep = duckvep_sweep_candidates(&s->v, &s->tx, HALO,
                                       active, KPROP_MAX_TX,
                                       candidates, KPROP_MAX_TX,
                                       pair_sink, &col, &st);
    if (st != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }

    /* Every in-window pair maps to a non-empty mask (non_coding_transcript/up/down) and the
     * filter is inert, so annotate emits exactly the sweep set in sweep order. */
    if (duckvep_result_builder_count(&rb) != n_sweep) { res = THEFT_TRIAL_FAIL; goto done; }

    for (i = 0u; i < n_sweep; i++) {
        uint32_t vi = (uint32_t)(col.buf[i] >> 32);
        uint32_t ti = (uint32_t)(col.buf[i] & 0xffffffffu);
        uint32_t start = s->v.pos1[vi];
        uint32_t end = s->v.variant_kind[vi] == (uint8_t)DUCKVEP_KIND_SV
                         ? s->v.end1[vi]
                         : start;
        uint32_t ts = s->tx.start1[ti];
        uint32_t te = s->tx.end1[ti];
        int fwd = s->tx.strand[ti] >= 0;
        uint32_t exp_region;
        uint64_t exp_mask;
        const duckvep_consequence_t *row = &rows[i];

        if (end < ts)        exp_region = fwd ? (uint32_t)DUCKVEP_REGION_UPSTREAM
                                              : (uint32_t)DUCKVEP_REGION_DOWNSTREAM;
        else if (start > te) exp_region = fwd ? (uint32_t)DUCKVEP_REGION_DOWNSTREAM
                                              : (uint32_t)DUCKVEP_REGION_UPSTREAM;
        else                 exp_region = (uint32_t)DUCKVEP_REGION_INTRON;

        if (exp_region == (uint32_t)DUCKVEP_REGION_UPSTREAM)
            exp_mask = DUCKVEP_SO(DUCKVEP_SO_UPSTREAM_GENE);
        else if (exp_region == (uint32_t)DUCKVEP_REGION_DOWNSTREAM)
            exp_mask = DUCKVEP_SO(DUCKVEP_SO_DOWNSTREAM_GENE);
        else
            /* zero-exon non-coding transcript: no exons -> no introns -> within_intron is
             * FALSE, so NO intron_variant; the in-span position is non_coding_transcript_variant
             * (placement). This is VEP-correct (the old INTRON co-emission was an over-call). */
            exp_mask = DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT);

        if (row->variant_idx != vi || row->tx_idx != ti) { res = THEFT_TRIAL_FAIL; goto done; }
        if (row->region_mask != exp_region) { res = THEFT_TRIAL_FAIL; goto done; }
        if (row->consequence_mask != exp_mask) { res = THEFT_TRIAL_FAIL; goto done; }
        if (row->impact != (uint8_t)DUCKVEP_IMPACT_MODIFIER) { res = THEFT_TRIAL_FAIL; goto done; }
        if (row->cdna_pos != -1 || row->cds_pos != -1 || row->protein_pos != -1) {
            res = THEFT_TRIAL_FAIL; goto done;
        }
    }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return res;
}

TEST annotate_matches_composition_for_any_scene(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile == sweep + classify + structural-SO composition";
    cfg.prop1 = prop_annotate_matches_composition;
    cfg.type_info[0] = &kprop_scene_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

static int consequence_rows_equal(const duckvep_consequence_t *a,
                                  const duckvep_consequence_t *b) {
    return a->variant_idx == b->variant_idx &&
           a->tx_idx == b->tx_idx &&
           a->gene_idx == b->gene_idx &&
           a->consequence_mask == b->consequence_mask &&
           a->region_mask == b->region_mask &&
           a->flags == b->flags &&
           a->impact == b->impact &&
           a->sequence_status == b->sequence_status &&
           a->nmd_prediction == b->nmd_prediction &&
           a->nmd_escape_reasons == b->nmd_escape_reasons &&
           a->cdna_pos == b->cdna_pos &&
           a->cds_pos == b->cds_pos &&
           a->protein_pos == b->protein_pos &&
           a->aa_ref == b->aa_ref &&
           a->aa_alt == b->aa_alt;
}

TEST annotate_cursor_resumes_known_scene(void) {
    static const uint16_t tchrom[2] = {0u, 0u};
    static const uint32_t tstart[2] = {100u, 100u};
    static const uint32_t tend[2]   = {200u, 200u};
    static const int8_t   strand[2] = {1, 1};
    static const uint64_t flags[2]  = {0u, 0u};
    static const uint32_t zero32[2] = {0u, 0u};
    static const uint16_t zero16[2] = {0u, 0u};
    uint16_t vchrom[1] = {0u};
    uint32_t vpos[1] = {150u};
    uint32_t vend[1] = {150u};
    uint8_t vkind[1] = {(uint8_t)DUCKVEP_KIND_SNV};
    duckvep_variant_batch_t v;
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_annotate_cursor_t *cur = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rowbuf[1];
    duckvep_consequence_t got[2];
    duckvep_result_builder_t rb;
    size_t got_n = 0u;

    memset(&v, 0, sizeof v); memset(&tx, 0, sizeof tx); memset(&ex, 0, sizeof ex); memset(&err, 0, sizeof err);
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind; v.count = 1u;
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = strand;
    tx.flags = flags; tx.exon_offset = zero32; tx.exon_count = zero16;
    tx.cds_start1 = zero32; tx.cds_end1 = zero32; tx.transcript_count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_cursor_open(model, &v, opts, ws, &cur, &err));

    while (!duckvep_annotate_cursor_done(cur)) {
        duckvep_status_t st;
        duckvep_result_builder_init(&rb, rowbuf, 1u);
        st = duckvep_annotate_cursor_fill(cur, &rb, &err);
        ASSERT(st == DUCKVEP_OK || st == DUCKVEP_ERR_RESULT_FULL);
        if (duckvep_result_builder_count(&rb) == 1u) {
            ASSERT(got_n < 2u);
            got[got_n++] = rowbuf[0];
        }
    }
    ASSERT_EQ(2u, got_n);
    ASSERT_EQ(0u, got[0].variant_idx); ASSERT_EQ(0u, got[0].tx_idx);
    ASSERT_EQ(0u, got[1].variant_idx); ASSERT_EQ(1u, got[1].tx_idx);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT), got[0].consequence_mask);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT), got[1].consequence_mask);

    duckvep_annotate_cursor_close(cur);
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST sorted_point_cursor_survives_tiles_and_resets_on_rewind(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1] = {310u};
    static const int8_t strand[1] = {1};
    static const uint64_t flags[1] = {0u};
    static const uint32_t exoff[1] = {0u};
    static const uint16_t excnt[1] = {3u};
    static const uint32_t zero[1] = {0u};
    static const uint32_t es[3] = {100u, 200u, 300u};
    static const uint32_t ee[3] = {120u, 220u, 310u};
    uint16_t vchrom[3] = {0u, 0u, 0u};
    uint32_t vpos[3];
    uint32_t vend[3];
    uint8_t vkind[3] = {
        (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_model_t *other_model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t rows[3];
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx);
    memset(&ex, 0, sizeof ex);
    memset(&v, 0, sizeof v);
    memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = strand; tx.flags = flags; tx.exon_offset = exoff;
    tx.exon_count = excnt; tx.cds_start1 = zero; tx.cds_end1 = zero;
    tx.transcript_count = 1u;
    ex.start1 = es; ex.end1 = ee; ex.exon_count = 3u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend;
    v.variant_kind = vkind;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_model_open(&tx, &ex, NULL, &other_model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));

    vpos[0] = vend[0] = 105u;
    vpos[1] = vend[1] = 150u;
    v.count = 2u;
    duckvep_result_builder_init(&rb, rows, 3u);
    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
              duckvep_annotate_tile(other_model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(65u, err.where_code);
    duckvep_result_builder_reset(&rb);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(2u, rb.count);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_EXON, rows[0].region_mask);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_INTRON, rows[1].region_mask);

    vpos[0] = vend[0] = 205u;
    vpos[1] = vend[1] = 250u;
    vpos[2] = vend[2] = 305u;
    v.count = 3u;
    duckvep_result_builder_init(&rb, rows, 3u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, rb.count);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_EXON, rows[0].region_mask);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_INTRON, rows[1].region_mask);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_EXON, rows[2].region_mask);

    /* A new run may rewind. The workspace detects that boundary and invalidates
     * its per-transcript ranks lazily, rather than carrying exon 3 into exon 1. */
    vpos[0] = vend[0] = 110u;
    v.count = 1u;
    duckvep_result_builder_init(&rb, rows, 3u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, rb.count);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_EXON, rows[0].region_mask);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(other_model);
    duckvep_model_close(model);
    PASS();
}

TEST padded_snv_rewind_uses_vep_feature_span(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1] = {220u};
    static const int8_t strand[1] = {1};
    static const uint64_t flags[1] = {0u};
    static const uint32_t exoff[1] = {0u};
    static const uint16_t excnt[1] = {2u};
    static const uint32_t zero[1] = {0u};
    static const uint32_t es[2] = {100u, 200u};
    static const uint32_t ee[2] = {120u, 220u};
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2] = {100u, 101u};
    static const uint32_t vend[2] = {210u, 105u};
    static const uint8_t vkind[2] = {
        (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV
    };
    static const uint32_t roff[2] = {0u, 222u};
    static const uint32_t aoff[2] = {111u, 227u};
    static const uint16_t rlen[2] = {111u, 5u};
    static const uint16_t alen[2] = {111u, 5u};
    uint8_t alleles[232];
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(alleles, 'A', sizeof alleles);
    alleles[aoff[0] + alen[0] - 1u] = 'C'; /* effect position 210 */
    alleles[aoff[1] + alen[1] - 1u] = 'C'; /* effect position 105 */
    memset(&tx, 0, sizeof tx); memset(&ex, 0, sizeof ex);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = strand; tx.flags = flags; tx.exon_offset = exoff;
    tx.exon_count = excnt; tx.cds_start1 = zero; tx.cds_end1 = zero;
    tx.transcript_count = 1u;
    ex.start1 = es; ex.end1 = ee; ex.exon_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend;
    v.variant_kind = vkind; v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen;
    v.allele_bytes = alleles; v.allele_bytes_len = sizeof alleles;
    v.count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(2u, rb.count);
    ASSERT_EQ((uint32_t)(DUCKVEP_REGION_EXON | DUCKVEP_REGION_INTRON),
              rows[0].region_mask);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_EXON, rows[1].region_mask);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* The resumable cursor is a transport primitive: arbitrary output-buffer splits
 * must equal one full annotate_tile call over the same tile. This property uses
 * the same random zero-exon scenes as the composition oracle and varies the chunk
 * cap from the generated dimensions, so it proves cursor state across variant,
 * active-set, and output-buffer boundaries without restating the sweep. */
static enum theft_trial_res prop_annotate_cursor_matches_tile(struct theft *t, void *arg1) {
    const struct kprop_scene *s = (const struct kprop_scene *)arg1;
    duckvep_exon_model_t exons;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_annotate_cursor_t *cur = NULL;
    duckvep_options_init_t init;
    duckvep_error_t err;
    static duckvep_consequence_t full[KPROP_MAX_PAIRS];
    static duckvep_consequence_t got[KPROP_MAX_PAIRS];
    duckvep_consequence_t chunk[9];
    duckvep_result_builder_t rb;
    duckvep_status_t st;
    enum theft_trial_res res = THEFT_TRIAL_PASS;
    size_t full_n, got_n = 0u, i;
    size_t chunk_cap;
    const uint32_t HALO = KPROP_SWEEP_HALO;
    (void)t;

    memset(&exons, 0, sizeof exons);
    memset(&err, 0, sizeof err);
    memset(&init, 0, sizeof init);
    init.upstream_dist = HALO;
    init.downstream_dist = HALO;
    init.halo = HALO;
    chunk_cap = 1u + ((s->v.count + s->tx.transcript_count + (size_t)s->halo) % 9u);

    if (duckvep_model_open(&s->tx, &exons, NULL, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(&init, &opts, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }

    duckvep_result_builder_init(&rb, full, KPROP_MAX_PAIRS);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    full_n = duckvep_result_builder_count(&rb);

    if (duckvep_annotate_cursor_open(model, &s->v, opts, ws, &cur, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    while (!duckvep_annotate_cursor_done(cur)) {
        duckvep_result_builder_init(&rb, chunk, chunk_cap);
        st = duckvep_annotate_cursor_fill(cur, &rb, &err);
        if (st != DUCKVEP_OK && st != DUCKVEP_ERR_RESULT_FULL) { res = THEFT_TRIAL_FAIL; goto done; }
        if (st == DUCKVEP_ERR_RESULT_FULL && duckvep_result_builder_count(&rb) == 0u) {
            res = THEFT_TRIAL_FAIL; goto done;
        }
        for (i = 0u; i < duckvep_result_builder_count(&rb); i++) {
            if (got_n >= KPROP_MAX_PAIRS) { res = THEFT_TRIAL_FAIL; goto done; }
            got[got_n++] = chunk[i];
        }
    }

    if (got_n != full_n) { res = THEFT_TRIAL_FAIL; goto done; }
    for (i = 0u; i < full_n; i++) {
        if (!consequence_rows_equal(&got[i], &full[i])) { res = THEFT_TRIAL_FAIL; goto done; }
    }

done:
    duckvep_annotate_cursor_close(cur);
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return res;
}

TEST annotate_cursor_matches_tile_for_any_output_split(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate cursor output splits == one annotate_tile";
    cfg.prop1 = prop_annotate_cursor_matches_tile;
    cfg.type_info[0] = &kprop_scene_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

TEST event_length_delta_pre_bits_follow_trimmed_alleles(void) {
    duckvep_effect_ctx_t ctx;
    duckvep_event_t event;

    memset(&ctx, 0, sizeof ctx);
    memset(&event, 0, sizeof event);
    event.kind = (uint8_t)DUCKVEP_KIND_INDEL;
    event.ref_diff_length = 1u;
    event.alt_diff_length = 4u;
    duckvep_effect_ctx_apply_event(&ctx, &event);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INSERTION)) != 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_DELETION)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SV)) == 0u);

    memset(&ctx, 0, sizeof ctx);
    memset(&event, 0, sizeof event);
    event.kind = (uint8_t)DUCKVEP_KIND_INDEL;
    event.ref_diff_length = 4u;
    event.alt_diff_length = 1u;
    duckvep_effect_ctx_apply_event(&ctx, &event);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_DELETION)) != 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INSERTION)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SV)) == 0u);

    memset(&ctx, 0, sizeof ctx);
    memset(&event, 0, sizeof event);
    event.kind = (uint8_t)DUCKVEP_KIND_MNV;
    event.ref_diff_length = 2u;
    event.alt_diff_length = 2u;
    duckvep_effect_ctx_apply_event(&ctx, &event);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INSERTION)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_DELETION)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SV)) == 0u);

    memset(&ctx, 0, sizeof ctx);
    memset(&event, 0, sizeof event);
    event.kind = (uint8_t)DUCKVEP_KIND_SV;
    duckvep_effect_ctx_apply_event(&ctx, &event);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SV)) != 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INSERTION)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_DELETION)) == 0u);

    PASS();
}

/* Deterministic SO-mapping anchor: a 2-exon coding transcript on chr0 (exons
 * [1000,1300] & [1700,2000], CDS [1100,1900]) and a single-exon non-coding
 * transcript on chr1 ([1000,1300]). Each variant lands in exactly one structural
 * bucket; the expected consequence_mask + impact are hand-computed from the SO
 * spec, so this cannot pass vacuously and pins every branch of the mapping. */
TEST annotate_structural_known_scene(void) {
    static const uint16_t tchrom[2]  = {0u, 1u};
    static const uint32_t tstart[2]  = {1000u, 1000u};
    static const uint32_t tend[2]    = {2000u, 1300u};
    static const int8_t   tstrand[2] = {1, 1};
    static const uint64_t tflags[2]  = {0u, 0u};
    static const uint32_t texoff[2]  = {0u, 2u};
    static const uint16_t texcnt[2]  = {2u, 1u};
    static const uint32_t tcds_s[2]  = {1100u, 0u};
    static const uint32_t tcds_e[2]  = {1900u, 0u};
    static const uint32_t estart[3]  = {1000u, 1700u, 1000u};
    static const uint32_t eend[3]    = {1300u, 2000u, 1300u};

    static const uint16_t vchrom[9] = {0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 1u};
    static const uint32_t vpos[9]   = {900u, 1050u, 1200u, 1301u, 1500u, 1750u, 1950u, 2100u, 1150u};
    static const uint8_t  vkind[9]  = {0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u, 0u};

    static const uint64_t exp_mask[9] = {
        DUCKVEP_SO(DUCKVEP_SO_UPSTREAM_GENE),                                       /* 900  up   */
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR),                                         /* 1050 5'utr*/
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE),                                     /* 1200 cds  */
        DUCKVEP_SO(DUCKVEP_SO_SPLICE_DONOR),                                       /* 1301 = 1st intronic base (essential donor); splice_donor ALONE (intron_variant suppressed at the essential site, VEP-faithful), HIGH impact flows through annotate_tile */
        DUCKVEP_SO(DUCKVEP_SO_INTRON),                                              /* 1500 intr */
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE),                                     /* 1750 cds  */
        DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR),                                         /* 1950 3'utr*/
        DUCKVEP_SO(DUCKVEP_SO_DOWNSTREAM_GENE),                                     /* 2100 down */
        DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT_EXON)                          /* 1150 nc exon -> exon term only */
    };
    static const uint8_t exp_impact[9] = {
        DUCKVEP_IMPACT_MODIFIER, DUCKVEP_IMPACT_MODIFIER, DUCKVEP_IMPACT_MODIFIER,
        DUCKVEP_IMPACT_HIGH,     DUCKVEP_IMPACT_MODIFIER, DUCKVEP_IMPACT_MODIFIER,
        DUCKVEP_IMPACT_MODIFIER, DUCKVEP_IMPACT_MODIFIER, DUCKVEP_IMPACT_MODIFIER
    };

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[16];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v);
    memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend; exons.exon_count = 3u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind; v.count = 9u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err)); /* defaults */
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 16u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(9u, duckvep_result_builder_count(&rb));

    for (i = 0u; i < 9u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(exp_impact[i], rows[i].impact);
        ASSERT_EQ(-1, rows[i].cds_pos);
        if (i == 2u || i == 5u) {
            ASSERT((rows[i].flags &
                    (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED) != 0u);
            ASSERT_EQ(DUCKVEP_SEQUENCE_MISSING, rows[i].sequence_status);
        } else {
            ASSERT_EQ(0u, rows[i].flags &
                         (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED);
            ASSERT_EQ(DUCKVEP_SEQUENCE_NOT_APPLICABLE,
                      rows[i].sequence_status);
        }
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_padded_small_variants_use_differing_region_topology(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {100u};
    static const uint32_t tend[1]    = {300u};
    static const int8_t   tstrand[1] = {1};
    static const uint64_t tflags[1]  = {0u};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {2u};
    static const uint32_t tcds_s[1]  = {120u};
    static const uint32_t tcds_e[1]  = {220u};
    static const uint32_t estart[2]  = {100u, 200u};
    static const uint32_t eend[2]    = {150u, 300u};

    static const uint16_t vchrom[3] = {0u, 0u, 0u};
    static const uint32_t vpos[3]   = {150u, 199u, 220u};
    static const uint32_t vend[3]   = {151u, 200u, 221u};
    static const uint8_t  vkind[3]  = {
        (uint8_t)DUCKVEP_KIND_DEL,
        (uint8_t)DUCKVEP_KIND_DEL,
        (uint8_t)DUCKVEP_KIND_DEL
    };
    static const uint8_t  bytes[9] = {
        'C','G', 'C', /* changed base 151 = first intronic donor base */
        'G','T', 'G', /* changed base 200 = exon start, not acceptor dinucleotide */
        'T','A', 'T'  /* changed base 221 = 3'UTR, not CDS stop anchor */
    };
    static const uint32_t roff[3] = {0u, 3u, 6u};
    static const uint32_t aoff[3] = {2u, 5u, 8u};
    static const uint16_t rlen[3] = {2u, 2u, 2u};
    static const uint16_t alen[3] = {1u, 1u, 1u};
    static const uint64_t exp_mask[3] = {
        DUCKVEP_SO(DUCKVEP_SO_SPLICE_DONOR),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) | DUCKVEP_SO(DUCKVEP_SO_SPLICE_REGION),
        DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR)
    };

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend; exons.exon_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = bytes; v.allele_bytes_len = sizeof bytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 3u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 3u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(-1, rows[i].cds_pos);
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* End-to-end structural/CNV anchor. Copy direction is explicit: an undirected
 * CNV does not become a gain or loss merely from its interval. Tier-1 complete
 * feature events suppress the ordinary placement terms; contained events retain
 * both their structural predicate and the applicable topology consequence. */
TEST annotate_sv_cnv_known_scene(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1] = {300u};
    static const int8_t strand[1] = {1};
    static const uint64_t flags[1] = {0u};
    static const uint32_t exoff[1] = {0u};
    static const uint16_t excnt[1] = {2u};
    static const uint32_t cds_s[1] = {120u};
    static const uint32_t cds_e[1] = {280u};
    static const uint32_t es[2] = {100u, 250u};
    static const uint32_t ee[2] = {150u, 300u};

    static const uint16_t vchrom[6] = {0u, 0u, 0u, 0u, 0u, 0u};
    static const uint32_t vstart[6] = {50u, 60u, 125u, 135u, 180u, 260u};
    static const uint32_t vend[6] = {350u, 340u, 130u, 140u, 220u, 265u};
    static const uint8_t vkind[6] = {
        DUCKVEP_KIND_SV, DUCKVEP_KIND_SV, DUCKVEP_KIND_SV,
        DUCKVEP_KIND_SV, DUCKVEP_KIND_SV, DUCKVEP_KIND_SV
    };
    static const uint8_t sv_type[6] = {
        DUCKVEP_SV_CNV, DUCKVEP_SV_DELETION, DUCKVEP_SV_CNV,
        DUCKVEP_SV_CNV, DUCKVEP_SV_INVERSION, DUCKVEP_SV_CNV
    };
    static const uint8_t copy_change[6] = {
        DUCKVEP_COPY_CHANGE_GAIN, DUCKVEP_COPY_CHANGE_LOSS,
        DUCKVEP_COPY_CHANGE_GAIN, DUCKVEP_COPY_CHANGE_LOSS,
        DUCKVEP_COPY_CHANGE_UNKNOWN, DUCKVEP_COPY_CHANGE_NEUTRAL
    };
    static const uint64_t expected[6] = {
        DUCKVEP_SO(DUCKVEP_SO_TRANSCRIPT_AMPLIFICATION),
        DUCKVEP_SO(DUCKVEP_SO_TRANSCRIPT_ABLATION),
        DUCKVEP_SO(DUCKVEP_SO_FEATURE_ELONGATION) |
            DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE),
        DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION) |
            DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE),
        DUCKVEP_SO(DUCKVEP_SO_INTRON),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE)
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;
    duckvep_error_t err;
    size_t i;

    memset(&tx, 0, sizeof tx);
    memset(&ex, 0, sizeof ex);
    memset(&v, 0, sizeof v);
    memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = strand;
    tx.flags = flags; tx.exon_offset = exoff; tx.exon_count = excnt;
    tx.cds_start1 = cds_s; tx.cds_end1 = cds_e; tx.transcript_count = 1u;
    ex.start1 = es; ex.end1 = ee; ex.exon_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vstart; v.end1 = vend;
    v.variant_kind = vkind; v.sv_type = sv_type; v.copy_change = copy_change;
    v.count = 6u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(6u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 6u; i++) {
        ASSERT_EQ((uint32_t)i, rows[i].variant_idx);
        ASSERT_EQ(expected[i], rows[i].consequence_mask);
        ASSERT_EQ(i < 4u ? (uint8_t)DUCKVEP_IMPACT_HIGH
                         : (uint8_t)DUCKVEP_IMPACT_MODIFIER,
                  rows[i].impact);
        ASSERT_EQ(-1, rows[i].cdna_pos);
        ASSERT_EQ(-1, rows[i].cds_pos);
        ASSERT_EQ(-1, rows[i].protein_pos);
    }

    /* Direct kernel callers bypass the adapter tile, so the borrowed-view
     * validator must also reject contradictory structural metadata. */
    {
        static const uint8_t bad_type[1] = {DUCKVEP_SV_DELETION};
        static const uint8_t bad_copy[1] = {DUCKVEP_COPY_CHANGE_GAIN};
        duckvep_variant_batch_t bad = v;
        bad.count = 1u;
        bad.sv_type = bad_type;
        bad.copy_change = bad_copy;
        duckvep_result_builder_reset(&rb);
        ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
                  duckvep_annotate_tile(model, &bad, opts, ws, &rb, &err));
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* Complete neutral/undirected structural spans use VEP's transcript-level
 * fallback predicates. They must not leak the point-era coding_sequence or
 * non_coding_transcript_exon fallback merely because the span crosses CDS/exon. */
TEST annotate_complete_neutral_sv_uses_transcript_fallbacks(void) {
    static const uint16_t tchrom[2] = {0u, 1u};
    static const uint32_t tstart[2] = {100u, 100u};
    static const uint32_t tend[2] = {300u, 300u};
    static const int8_t strand[2] = {1, 1};
    static const uint64_t flags[2] = {0u, 0u};
    static const uint32_t exoff[2] = {0u, 1u};
    static const uint16_t excnt[2] = {1u, 1u};
    static const uint32_t cds_s[2] = {100u, 0u};
    static const uint32_t cds_e[2] = {300u, 0u};
    static const uint32_t es[2] = {100u, 100u};
    static const uint32_t ee[2] = {300u, 300u};
    static const uint16_t vchrom[2] = {0u, 1u};
    static const uint32_t vstart[2] = {50u, 50u};
    static const uint32_t vend[2] = {350u, 350u};
    static const uint8_t vkind[2] = {DUCKVEP_KIND_SV, DUCKVEP_KIND_SV};
    static const uint8_t sv_type[2] = {DUCKVEP_SV_CNV, DUCKVEP_SV_CNV};
    static const uint8_t copy_change[2] = {
        DUCKVEP_COPY_CHANGE_NEUTRAL, DUCKVEP_COPY_CHANGE_NEUTRAL
    };
    static const uint64_t expected[2] = {
        DUCKVEP_SO(DUCKVEP_SO_CODING_TRANSCRIPT),
        DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT)
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    duckvep_error_t err;
    size_t i;

    memset(&tx, 0, sizeof tx);
    memset(&ex, 0, sizeof ex);
    memset(&v, 0, sizeof v);
    memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = strand;
    tx.flags = flags; tx.exon_offset = exoff; tx.exon_count = excnt;
    tx.cds_start1 = cds_s; tx.cds_end1 = cds_e; tx.transcript_count = 2u;
    ex.start1 = es; ex.end1 = ee; ex.exon_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vstart; v.end1 = vend;
    v.variant_kind = vkind; v.sv_type = sv_type; v.copy_change = copy_change;
    v.count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(2u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 2u; i++) {
        ASSERT_EQ((uint32_t)i, rows[i].variant_idx);
        ASSERT_EQ(expected[i], rows[i].consequence_mask);
        ASSERT_EQ((uint8_t)DUCKVEP_IMPACT_MODIFIER, rows[i].impact);
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* Regression: the emitted row->region_mask is reconciled to AGREE with emission.
 * duckvep_region_mask (the coarse pre-classifier) sets REGION_SPLICE at EVERY exon
 * boundary including the transcript's outer 5'/3' ends, which are NOT real splice
 * sites — so it over-calls. annotate_pair must emit the authoritative splice fact
 * (splice_classify), never the coarse second opinion. This pins: (a) a position in
 * the coarse splice reach of the OUTER 5' end where region_mask must have NO SPLICE
 * bit and emit only 5_prime_utr; (b) a real donor site where SPLICE is retained. */
TEST annotate_region_mask_truthful_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {2000u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {2u};
    static const uint32_t tcds_s[1]  = {1100u};
    static const uint32_t tcds_e[1]  = {1900u};
    static const uint32_t estart[2]  = {1000u, 1700u};
    static const uint32_t eend[2]    = {1300u, 2000u};

    /* 1002 = 3rd base of the FIRST exon: inside the coarse exonic splice reach of the
     * transcript's outer 5' end (1000), but no intron there -> not a splice site.
     * 1301 = first intronic base of intron [1301,1699] -> a real essential donor. */
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2]   = {1002u, 1301u};
    static const uint8_t  vkind[2]  = {0u, 0u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v);
    memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend; exons.exon_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind; v.count = 2u;

    /* Document the coarse over-call we are suppressing: the pre-classifier DOES flag
     * SPLICE at 1002 (outer-end reach). The emitted row below must not. */
    ASSERT((duckvep_region_mask(&tx, &exons, 0, 1002u, 3u, 8u) &
            (uint32_t)DUCKVEP_REGION_SPLICE) != 0u);

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(2u, duckvep_result_builder_count(&rb));

    /* (a) outer-end over-call suppressed: pure 5'UTR, region_mask has NO SPLICE bit. */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR), rows[0].consequence_mask);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_UTR, rows[0].region_mask);
    ASSERT_EQ(0u, rows[0].region_mask & (uint32_t)DUCKVEP_REGION_SPLICE);

    /* (b) real essential donor: emitted as splice_donor ALONE because the dinucleotide is a
     * splice site, not within_intron. The region_mask still records intron placement plus
     * splice proximity; it is a structural diagnostic, not the emitted consequence set. */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_SPLICE_DONOR), rows[1].consequence_mask);
    ASSERT((rows[1].region_mask & (uint32_t)DUCKVEP_REGION_SPLICE) != 0u);
    ASSERT((rows[1].region_mask & (uint32_t)DUCKVEP_REGION_INTRON) != 0u);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* Reverse-strand SO-mapping anchor: same coding transcript geometry as above but
 * on the '-' strand. Strand must flip upstream<->downstream and 5'<->3' UTR
 * relative to genomic coordinates; CDS and intron stay strand-independent. */
TEST annotate_reverse_strand_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {2000u};
    static const int8_t   tstrand[1] = {-1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {2u};
    static const uint32_t tcds_s[1]  = {1100u};
    static const uint32_t tcds_e[1]  = {1900u};
    static const uint32_t estart[2]  = {1700u, 1000u};
    static const uint32_t eend[2]    = {2000u, 1300u};

    static const uint16_t vchrom[6] = {0u, 0u, 0u, 0u, 0u, 0u};
    static const uint32_t vpos[6]   = {900u, 1050u, 1200u, 1500u, 1950u, 2100u};
    static const uint8_t  vkind[6]  = {0u, 0u, 0u, 0u, 0u, 0u};

    static const uint64_t exp_mask[6] = {
        DUCKVEP_SO(DUCKVEP_SO_DOWNSTREAM_GENE),    /* 900  genomic-upstream -> downstream on '-' */
        DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR),        /* 1050 before CDS start -> 3' on '-'         */
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE),    /* 1200 cds                                    */
        DUCKVEP_SO(DUCKVEP_SO_INTRON),             /* 1500 intron                                 */
        DUCKVEP_SO(DUCKVEP_SO_5_PRIME_UTR),        /* 1950 after CDS end -> 5' on '-'             */
        DUCKVEP_SO(DUCKVEP_SO_UPSTREAM_GENE)       /* 2100 genomic-downstream -> upstream on '-'  */
    };

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend; exons.exon_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind; v.count = 6u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(6u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 6u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* Codon-bucket fusion: a single-exon coding transcript on '+' strand with a CDS
 * sequence pool. annotate_tile must refine the generic coding_sequence_variant
 * into the specific codon consequence via the projection + coding-SNV kernels.
 * CDS = ATG AAA TTT TGG CCC (M K F W P). Expected, hand-computed from the genetic
 * code (the oracle): a synonymous, a missense, and a stop_gained SNV, with the
 * right cds_pos / protein_pos / aa. */
TEST annotate_codon_snv_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1014u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {15u};
    static const int8_t   ephase[1]  = {0};

    /* genomic == transcript (+ strand); cds_bytes[0] is CDS position 1. */
    static const uint8_t cds_bytes[15] = {
        'A','T','G', 'A','A','A', 'T','T','T', 'T','G','G', 'C','C','C'
    };
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {15u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    /* 3 SNVs: cds6 A>G (AAA->AAG, K->K syn), cds8 T>A (TTT->TAT, F->Y mis),
     * cds11 G>A (TGG->TAG, W->* stop_gained). */
    static const uint16_t vchrom[3] = {0u, 0u, 0u};
    static const uint32_t vpos[3]   = {1005u, 1007u, 1010u};
    static const uint8_t  vkind[3]  = {0u, 0u, 0u}; /* SNV */
    static const uint8_t  abytes[6] = {'A','G', 'T','A', 'G','A'};
    static const uint32_t roff[3]   = {0u, 2u, 4u};
    static const uint32_t aoff[3]   = {1u, 3u, 5u};
    static const uint16_t rlen[3]   = {1u, 1u, 1u};
    static const uint16_t alen[3]   = {1u, 1u, 1u};

    static const uint64_t exp_mask[3] = {
        DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED)
    };
    static const int32_t exp_cds[3]     = {6, 8, 11};
    static const int32_t exp_protein[3] = {2, 3, 4};
    static const char    exp_aa_ref[3]  = {'K', 'F', 'W'};
    static const char    exp_aa_alt[3]  = {'K', 'Y', '*'};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 15u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 3u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, duckvep_result_builder_count(&rb));

    for (i = 0u; i < 3u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(exp_cds[i], rows[i].cds_pos);
        ASSERT_EQ(exp_protein[i], rows[i].protein_pos);
        ASSERT_EQ((uint8_t)exp_aa_ref[i], rows[i].aa_ref);
        ASSERT_EQ((uint8_t)exp_aa_alt[i], rows[i].aa_alt);
        ASSERT_EQ(DUCKVEP_SEQUENCE_RESOLVED, rows[i].sequence_status);
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* Reverse-strand codon fusion: single-exon '-' transcript over genomic
 * [1000,1014]; the transcript CDS is the reverse-complement of the genome, so
 * cdna = 1015 - genomic_pos. Transcript CDS = ATG AAA TTT TGG CCC (M K F W P), so
 * the GENOMIC bases are the complement of the transcript bases. Same three edits
 * as the '+' scene, expressed in genomic orientation; annotate_tile must reverse-
 * complement them (via coding_snv_from_cds) and land the identical aa/cds/protein.
 * Variants are listed in ascending genomic order (the kernel's sort precondition);
 * each is hand-derived genomic->transcript via complement (audited before run):
 *   1004 C>T genomic -> transcript G>A at cds11 -> TGG->TAG = W->* stop_gained
 *   1007 A>T genomic -> transcript T>A at cds8  -> TTT->TAT = F->Y missense
 *   1009 T>C genomic -> transcript A>G at cds6  -> AAA->AAG = K->K synonymous
 */
TEST annotate_codon_reverse_strand_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {-1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1014u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};  /* cdna 1 at the high genomic end (1014) */
    static const uint32_t ecdna_e[1] = {15u};
    static const int8_t   ephase[1]  = {0};

    /* TRANSCRIPT-oriented CDS (what the kernel translates). */
    static const uint8_t cds_bytes[15] = {
        'A','T','G', 'A','A','A', 'T','T','T', 'T','G','G', 'C','C','C'
    };
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {15u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    /* GENOMIC ref/alt (complement of the transcript bases). */
    static const uint16_t vchrom[3] = {0u, 0u, 0u};
    static const uint32_t vpos[3]   = {1004u, 1007u, 1009u};
    static const uint8_t  vkind[3]  = {0u, 0u, 0u};
    static const uint8_t  abytes[6] = {'C','T', 'A','T', 'T','C'}; /* sg: C>T, mis: A>T, syn: T>C */
    static const uint32_t roff[3]   = {0u, 2u, 4u};
    static const uint32_t aoff[3]   = {1u, 3u, 5u};
    static const uint16_t rlen[3]   = {1u, 1u, 1u};
    static const uint16_t alen[3]   = {1u, 1u, 1u};

    static const uint64_t exp_mask[3] = {
        DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS)
    };
    static const int32_t exp_cds[3]     = {11, 8, 6};
    static const int32_t exp_protein[3] = {4, 3, 2};
    static const char    exp_aa_ref[3]  = {'W', 'F', 'K'};
    static const char    exp_aa_alt[3]  = {'*', 'Y', 'K'};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 15u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 3u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 3u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(exp_cds[i], rows[i].cds_pos);
        ASSERT_EQ(exp_protein[i], rows[i].protein_pos);
        ASSERT_EQ((uint8_t)exp_aa_ref[i], rows[i].aa_ref);
        ASSERT_EQ((uint8_t)exp_aa_alt[i], rows[i].aa_alt);
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* Ref-mismatch fallback (BLACK-BOX behavior, not a branch proof): a coding SNV
 * whose stated genomic REF does not match the reference base must NOT be
 * force-classified into a codon term — the row stays the generic
 * coding_sequence_variant with no codon coordinates (cds_pos == -1). (The
 * REF_MISMATCH status itself is asserted directly by the coding_snv_from_cds unit
 * tests; here we only pin annotate_tile's observable fallback.) */
TEST annotate_codon_ref_mismatch_falls_back(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1008u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1008u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {9u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'A','A','A', 'T','T','T'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    /* pos 1003 = cds 4 (ref 'A'); we LIE and say REF is 'C'. */
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {1003u};
    static const uint8_t  vkind[1]  = {0u};
    static const uint8_t  abytes[2] = {'C', 'G'}; /* wrong ref C, alt G */
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {1u};
    static const uint16_t rlen[1]   = {1u};
    static const uint16_t alen[1]   = {1u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE), rows[0].consequence_mask);
    ASSERT_EQ(-1, rows[0].cds_pos);      /* not refined */
    ASSERT_EQ(0, (int)rows[0].aa_ref);
    ASSERT_EQ(DUCKVEP_SEQUENCE_REFERENCE_MISMATCH,
              rows[0].sequence_status);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* stop_retained anchor: a synonymous change AT the stop codon (stop -> stop) is
 * stop_retained_variant, not synonymous_variant. CDS = ATG AAA TAA (M K *), with a
 * 3' UTR after it (exon to 1014, CDS to 1008) so the terminal stop is not within
 * splice reach of the exon end. SNV at cds8 A>G: TAA -> TGA (still a stop). */
TEST annotate_codon_stop_retained_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u}; /* CDS ends before the exon (3' UTR follows) */
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {15u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'A','A','A', 'T','A','A'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {1007u}; /* cds 8, 2nd base of the stop codon */
    static const uint8_t  vkind[1]  = {0u};
    static const uint8_t  abytes[2] = {'A', 'G'}; /* TAA -> TGA */
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {1u};
    static const uint16_t rlen[1]   = {1u};
    static const uint16_t alen[1]   = {1u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED), rows[0].consequence_mask);
    ASSERT_EQ(8, rows[0].cds_pos);
    ASSERT_EQ(3, rows[0].protein_pos);
    ASSERT_EQ((uint8_t)'*', rows[0].aa_ref);
    ASSERT_EQ((uint8_t)'*', rows[0].aa_alt);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* start_lost anchor: VEP's SNV start_lost is not just "missense at protein 1".
 * It requires overlap with the annotated start codon plus translation_start==1;
 * then missense is suppressed, while stop_gained/synonymous can co-emit. The
 * non-ATG TAC start codon mirrors the in-tree VEP --gff witness fixture. */
TEST annotate_codon_start_lost_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {900u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {900u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {115u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'T','A','C', 'A','A','A', 'T','A','A'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};
    static const uint64_t nf_flags[1] = {(uint64_t)DUCKVEP_TX_CDS_START_NF};

    static const uint16_t vchrom[3] = {0u, 0u, 0u};
    static const uint32_t vpos[3]   = {1000u, 1002u, 1002u};
    static const uint8_t  vkind[3]  = {0u, 0u, 0u};
    static const uint8_t  abytes[6] = {'T','A', 'C','A', 'C','T'};
    static const uint32_t roff[3]   = {0u, 2u, 4u};
    static const uint32_t aoff[3]   = {1u, 3u, 5u};
    static const uint16_t rlen[3]   = {1u, 1u, 1u};
    static const uint16_t alen[3]   = {1u, 1u, 1u};
    static const uint64_t exp_mask[3] = {
        DUCKVEP_SO(DUCKVEP_SO_START_LOST),
        DUCKVEP_SO(DUCKVEP_SO_START_LOST) | DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
        DUCKVEP_SO(DUCKVEP_SO_START_LOST) | DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS)
    };
    static const char exp_alt[3] = {'N', '*', 'Y'};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 3u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 3u; i++) {
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(1, rows[i].protein_pos);
        ASSERT_EQ((uint8_t)'Y', rows[i].aa_ref);
        ASSERT_EQ((uint8_t)exp_alt[i], rows[i].aa_alt);
    }

    duckvep_workspace_close(ws); ws = NULL;
    duckvep_options_close(opts); opts = NULL;
    duckvep_model_close(model); model = NULL;

    /* VEP _overlaps_start_codon returns false for cds_start_NF transcripts, so the
     * same first-codon SNV falls back to ordinary missense rather than start_lost. */
    tx.flags = nf_flags;
    v.count = 1u;
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_MISSENSE), rows[0].consequence_mask);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* Multi-exon codon stitching: a 2-exon '+' transcript (exon0 [1000,1010] cdna
 * 1-11, exon1 [2000,2010] cdna 12-22), CDS cdna 1..21 = ATG AAA AAA CGT AAA AAA AAA
 * (M K K R K K K). codon4 (CGT) STRADDLES the junction: cds10,11 in exon0
 * (genomic 1009,1010), cds12 in exon1 (genomic 2000). Three uploaded alleles:
 *   - genomic 1004 (cds5, mid-exon0, no splice): AAA->ATA = K->I missense.
 *   - genomic 1009 CGA>CTA spans the exon/intron boundary but differs only at
 *     genomic 1010. VEP maps the whole equal-length feature, cannot form peptide
 *     alleles, and emits coding_sequence_variant plus splice_region_variant.
 *   - genomic 2000 (cds12, the split base, exon1 ACCEPTOR): CGT->CGC = R->R
 *     synonymous AND splice_region (a real internal boundary). Correctly emitting
 *     synonymous here proves the codon was assembled across the junction: a wrong
 *     contiguous-genomic read of 2000..2002 would be TAA (a stop), not Arg, so an
 *     R->R synonymous result is only reachable via the spliced CDS (1009,1010,2000). */
TEST annotate_codon_multi_exon_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {2010u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {2u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {2009u};
    static const uint32_t estart[2]  = {1000u, 2000u};
    static const uint32_t eend[2]    = {1010u, 2010u};
    static const uint32_t ecdna_s[2] = {1u, 12u};
    static const uint32_t ecdna_e[2] = {11u, 22u};
    static const int8_t   ephase[2]  = {0, 0};
    /* spliced transcript CDS (21 nt). */
    static const uint8_t cds_bytes[21] = {
        'A','T','G', 'A','A','A', 'A','A','A', 'C','G','T',
        'A','A','A', 'A','A','A', 'A','A','A'
    };
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {21u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    static const uint16_t vchrom[3] = {0u, 0u, 0u};
    static const uint32_t vpos[3]   = {1004u, 1009u, 2000u};
    static const uint32_t vend[3]   = {1004u, 1011u, 2000u};
    static const uint8_t  vkind[3]  = {0u, 0u, 0u};
    static const uint8_t  abytes[10] = {
        'A','T', 'C','G','A', 'C','T','A', 'T','C'
    }; /* 1004 A>T, 1009 CGA>CTA, 2000 T>C */
    static const uint32_t roff[3]   = {0u, 2u, 8u};
    static const uint32_t aoff[3]   = {1u, 5u, 9u};
    static const uint16_t rlen[3]   = {1u, 3u, 1u};
    static const uint16_t alen[3]   = {1u, 3u, 1u};

    static const uint64_t exp_mask[3] = {
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
            DUCKVEP_SO(DUCKVEP_SO_SPLICE_REGION),
        DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS) | DUCKVEP_SO(DUCKVEP_SO_SPLICE_REGION)
    };
    static const int32_t exp_cds[3]     = {5, -1, 12};
    static const int32_t exp_protein[3] = {2, -1, 4};
    static const char    exp_aa_ref[3]  = {'K', '\0', 'R'};
    static const char    exp_aa_alt[3]  = {'I', '\0', 'R'};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 21u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 3u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 3u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(exp_cds[i], rows[i].cds_pos);
        ASSERT_EQ(exp_protein[i], rows[i].protein_pos);
        ASSERT_EQ((uint8_t)exp_aa_ref[i], rows[i].aa_ref);
        ASSERT_EQ((uint8_t)exp_aa_alt[i], rows[i].aa_alt);
    }
    ASSERT_EQ(DUCKVEP_SEQUENCE_NOT_APPLICABLE, rows[1].sequence_status);
    ASSERT_EQ(0u, rows[1].flags &
                   (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* Same-codon MNV anchor: a multi-base substitution contained within one codon uses the
 * same codon-change oracle as SNVs. Two-codon non-terminal body MNVs can emit the narrow
 * missense slice below; broader codon/exon/CDS-boundary MNVs still fall back to
 * coding_sequence_variant until the general edit-set path lands. */
TEST annotate_codon_mnv_same_codon_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {15u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'A','A','A', 'T','T','T'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    /* 2bp MNV at genomic 1003-1004 (cds 4-5), AA>GG. */
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {1003u};
    static const uint32_t vend[1]   = {1004u};
    static const uint8_t  vkind[1]  = {(uint8_t)DUCKVEP_KIND_MNV};
    static const uint8_t  abytes[4] = {'A','A', 'G','G'};
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {2u};
    static const uint16_t rlen[1]   = {2u};
    static const uint16_t alen[1]   = {2u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_MISSENSE), rows[0].consequence_mask);
    ASSERT_EQ(-1, rows[0].cdna_pos);
    ASSERT_EQ(-1, rows[0].cds_pos);
    ASSERT_EQ(2, rows[0].protein_pos);
    ASSERT_EQ((uint8_t)'K', rows[0].aa_ref);
    ASSERT_EQ((uint8_t)'G', rows[0].aa_alt);
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    ASSERT_EQ(1u, stats->substitution_context);
    ASSERT_EQ(0u, stats->mnv_direct_fallback);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_mnv_start_lost_route_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1008u};
    static const int8_t   tstrand[1] = {1};
    static const uint64_t tflags[1]  = {0u};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1008u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {9u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'G','A','A', 'T','A','A'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {1000u};
    static const uint32_t vend[1]   = {1002u};
    static const uint8_t  vkind[1]  = {(uint8_t)DUCKVEP_KIND_MNV};
    static const uint8_t  abytes[6] = {'A','T','G', 'G','C','A'};
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {3u};
    static const uint16_t rlen[1]   = {3u};
    static const uint16_t alen[1]   = {3u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_START_LOST), rows[0].consequence_mask);
    ASSERT_EQ(1, rows[0].protein_pos);
    ASSERT_EQ((uint8_t)'M', rows[0].aa_ref);
    ASSERT_EQ((uint8_t)'A', rows[0].aa_alt);
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    ASSERT_EQ(1u, stats->substitution_context);
    ASSERT_EQ(0u, stats->mnv_direct_fallback);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_cursor_mnv_start_lost_route_matches_tile_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1008u};
    static const int8_t   tstrand[1] = {1};
    static const uint64_t tflags[1]  = {0u};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1008u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {9u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'G','A','A', 'T','A','A'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {1000u};
    static const uint32_t vend[1]   = {1002u};
    static const uint8_t  vkind[1]  = {(uint8_t)DUCKVEP_KIND_MNV};
    static const uint8_t  abytes[6] = {'A','T','G', 'G','C','A'};
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {3u};
    static const uint16_t rlen[1]   = {3u};
    static const uint16_t alen[1]   = {3u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_annotate_cursor_t *cur = NULL;
    duckvep_error_t err;
    duckvep_consequence_t tile_rows[4];
    duckvep_consequence_t cursor_rows[4];
    duckvep_consequence_t chunk[1];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    duckvep_workspace_delta_route_stats_t tile_stats;
    duckvep_workspace_delta_route_stats_t cursor_stats;
    size_t cursor_n = 0u;
    int saw_full = 0;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    memset(&tile_stats, 0, sizeof tile_stats); memset(&cursor_stats, 0, sizeof cursor_stats);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));

    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, tile_rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    tile_stats = *stats;

    duckvep_workspace_delta_route_stats_reset(ws);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_cursor_open(model, &v, opts, ws, &cur, &err));
    while (!duckvep_annotate_cursor_done(cur)) {
        duckvep_status_t st;
        size_t i;
        duckvep_result_builder_init(&rb, chunk, 1u);
        st = duckvep_annotate_cursor_fill(cur, &rb, &err);
        ASSERT(st == DUCKVEP_OK || st == DUCKVEP_ERR_RESULT_FULL);
        if (st == DUCKVEP_ERR_RESULT_FULL) saw_full = 1;
        for (i = 0u; i < duckvep_result_builder_count(&rb); i++) {
            ASSERT(cursor_n < 4u);
            cursor_rows[cursor_n++] = chunk[i];
        }
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    cursor_stats = *stats;

    ASSERT(saw_full);
    ASSERT_EQ(1u, cursor_n);
    ASSERT(consequence_rows_equal(&tile_rows[0], &cursor_rows[0]));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_START_LOST), cursor_rows[0].consequence_mask);
    ASSERT_EQ(1u, tile_stats.substitution_context);
    ASSERT_EQ(0u, tile_stats.mnv_direct_fallback);
    ASSERT_EQ(tile_stats.substitution_context, cursor_stats.substitution_context);
    ASSERT_EQ(tile_stats.mnv_direct_fallback, cursor_stats.mnv_direct_fallback);

    duckvep_annotate_cursor_close(cur);
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_mnv_len3_and_cross_codon_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {15u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'A','A','A', 'T','T','T'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    /* Variant 0 replaces all of codon 2 (AAA>GGG = K->G). Variant 1 spans codons
     * 2/3 (cds 6-7) but touches the terminal codon in this 9nt CDS, so it stays the
     * generic CDS bucket until terminal-boundary edit-set support. */
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2]   = {1003u, 1005u};
    static const uint32_t vend[2]   = {1005u, 1006u};
    static const uint8_t  vkind[2]  = {(uint8_t)DUCKVEP_KIND_MNV, (uint8_t)DUCKVEP_KIND_MNV};
    static const uint8_t  abytes[10] = {'A','A','A', 'G','G','G', 'A','T', 'G','C'};
    static const uint32_t roff[2]   = {0u, 6u};
    static const uint32_t aoff[2]   = {3u, 8u};
    static const uint16_t rlen[2]   = {3u, 2u};
    static const uint16_t alen[2]   = {3u, 2u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(2u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_MISSENSE), rows[0].consequence_mask);
    ASSERT_EQ(-1, rows[0].cdna_pos);
    ASSERT_EQ(-1, rows[0].cds_pos);
    ASSERT_EQ(2, rows[0].protein_pos);
    ASSERT_EQ((uint8_t)'K', rows[0].aa_ref);
    ASSERT_EQ((uint8_t)'G', rows[0].aa_alt);
    /* Cross-codon MNV (codon 2 synonymous K->K, codon 3 missense F->L): the generalized
     * window classifier resolves the whole window as a coarse missense_variant, where the
     * old two-codon slice bailed to coding_sequence_variant. Both variants now route through
     * the authoritative interpreter (accepted), so there is no fallback. */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_MISSENSE), rows[1].consequence_mask);
    ASSERT_EQ(-1, rows[1].cdna_pos);
    ASSERT_EQ(-1, rows[1].cds_pos);
    ASSERT_EQ(-1, rows[1].protein_pos);
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    ASSERT_EQ(2u, stats->substitution_context);
    ASSERT_EQ(0u, stats->mnv_direct_fallback);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_mnv_cross_codon_missense_known_scene(void) {
    static const uint16_t tchrom[2]  = {0u, 1u};
    static const uint32_t tstart[2]  = {1000u, 2000u};
    static const uint32_t tend[2]    = {1014u, 2014u};
    static const int8_t   tstrand[2] = {1, -1};
    static const uint64_t tflags[2]  = {0u, 0u};
    static const uint32_t texoff[2]  = {0u, 1u};
    static const uint16_t texcnt[2]  = {1u, 1u};
    static const uint32_t tcds_s[2]  = {1000u, 2000u};
    static const uint32_t tcds_e[2]  = {1014u, 2014u};
    static const uint32_t estart[2]  = {1000u, 2000u};
    static const uint32_t eend[2]    = {1014u, 2014u};
    static const uint32_t ecdna_s[2] = {1u, 1u};
    static const uint32_t ecdna_e[2] = {15u, 15u};
    static const int8_t   ephase[2]  = {0, 0};
    static const uint8_t  cds_bytes[30] = {
        'A','T','G', 'A','A','A', 'T','T','A', 'G','G','G', 'T','T','T',
        'A','T','G', 'A','A','A', 'T','T','A', 'G','G','G', 'T','T','T'
    };
    static const uint64_t cds_off[2]  = {0u, 15u};
    static const uint32_t cds_lenv[2] = {15u, 15u};
    static const uint8_t  cds_tab[2]  = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };

    static const uint16_t vchrom[5] = {0u, 0u, 0u, 0u, 1u};
    static const uint32_t vpos[5]   = {1002u, 1005u, 1005u, 1011u, 2008u};
    static const uint32_t vend[5]   = {1003u, 1006u, 1006u, 1012u, 2009u};
    static const uint8_t  vkind[5]  = {
        (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV,
        (uint8_t)DUCKVEP_KIND_MNV
    };
    static const uint8_t  abytes[20] = {
        'G','A', 'T','T',  /* +: ATG/AAA -> ATT/TAA = start_lost (M->I) & stop_gained (K->*) */
        'A','T', 'G','G',  /* +: AAA/TTA -> AAG/GTA = one syn, one missense => missense */
        'A','T', 'G','C',  /* +: AAA/TTA -> AAG/CTA = both synonymous => synonymous */
        'G','T', 'A','A',  /* +: GGG/TTT -> GGA/ATT = one syn, one missense => missense */
        'A','T', 'C','C'   /* -: genomic AT>CC -> transcript AT>GG => missense */
    };
    static const uint32_t roff[5] = {0u, 4u, 8u, 12u, 16u};
    static const uint32_t aoff[5] = {2u, 6u, 10u, 14u, 18u};
    static const uint16_t rlen[5] = {2u, 2u, 2u, 2u, 2u};
    static const uint16_t alen[5] = {2u, 2u, 2u, 2u, 2u};
    /* The generalized window classifier resolves every two-codon window (including the
     * start&stop composite and both-synonymous), so none fall back to coding_sequence. */
    static const uint64_t exp_mask[5] = {
        DUCKVEP_SO(DUCKVEP_SO_START_LOST) | DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE)
    };

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 5u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(5u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 5u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(-1, rows[i].cdna_pos);
        ASSERT_EQ(-1, rows[i].cds_pos);
        ASSERT_EQ(-1, rows[i].protein_pos);
        ASSERT_EQ((uint8_t)0u, rows[i].aa_ref);
        ASSERT_EQ((uint8_t)0u, rows[i].aa_alt);
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    ASSERT_EQ(5u, stats->substitution_context);
    ASSERT_EQ(0u, stats->mnv_direct_fallback);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_mnv_reverse_strand_same_codon_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1008u};
    static const int8_t   tstrand[1] = {-1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1008u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {9u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'A','A','A', 'T','T','T'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    /* Genomic 1004-1005 is transcript CDS 5-4 on the reverse strand. TT>CC genomic
     * reverse-complements to AA>GG in transcript codon 2: AAA -> GGA = K->G. */
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {1004u};
    static const uint32_t vend[1]   = {1005u};
    static const uint8_t  vkind[1]  = {(uint8_t)DUCKVEP_KIND_MNV};
    static const uint8_t  abytes[4] = {'T','T', 'C','C'};
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {2u};
    static const uint16_t rlen[1]   = {2u};
    static const uint16_t alen[1]   = {2u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_MISSENSE), rows[0].consequence_mask);
    ASSERT_EQ(-1, rows[0].cdna_pos);
    ASSERT_EQ(-1, rows[0].cds_pos);
    ASSERT_EQ(2, rows[0].protein_pos);
    ASSERT_EQ((uint8_t)'K', rows[0].aa_ref);
    ASSERT_EQ((uint8_t)'G', rows[0].aa_alt);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* VEP keeps an equal-length uploaded feature intact for CDS mapping even when
 * semantic trimming leaves one changed coding base. If the uploaded span also
 * reaches the transcript's 3' UTR, one CDS endpoint is a mapper gap and VEP
 * cannot form peptide alleles: coding_unknown wins. The adjacent CDS-only
 * controls prove this is a mapping boundary, not blanket stop-loss suppression. */
TEST annotate_equal_length_cds_to_utr3_mapping_gap_known_scene(void) {
    static const uint16_t tchrom[2] = {0u, 1u};
    static const uint32_t tstart[2] = {1000u, 2000u};
    static const uint32_t tend[2] = {1011u, 2011u};
    static const int8_t tstrand[2] = {1, -1};
    static const uint64_t tflags[2] = {0u, 0u};
    static const uint32_t texoff[2] = {0u, 1u};
    static const uint16_t texcnt[2] = {1u, 1u};
    static const uint32_t tcds_s[2] = {1000u, 2003u};
    static const uint32_t tcds_e[2] = {1008u, 2011u};
    static const uint32_t estart[2] = {1000u, 2000u};
    static const uint32_t eend[2] = {1011u, 2011u};
    static const uint32_t ecdna_s[2] = {1u, 1u};
    static const uint32_t ecdna_e[2] = {12u, 12u};
    static const int8_t ephase[2] = {0, 0};
    static const uint8_t cds_bytes[18] = {
        'A','T','G', 'A','A','A', 'T','A','A',
        'A','T','G', 'A','A','A', 'T','A','A'
    };
    static const uint64_t cds_off[2] = {0u, 9u};
    static const uint32_t cds_len[2] = {9u, 9u};
    static const uint8_t cds_tab[2] = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };
    static const uint8_t post_cds[6] = {'A','A','A', 'A','A','A'};

    static const uint16_t vchrom[4] = {0u, 0u, 1u, 1u};
    static const uint32_t vpos[4] = {1007u, 1008u, 2002u, 2003u};
    static const uint32_t vend[4] = {1008u, 1009u, 2003u, 2004u};
    static const uint8_t vkind[4] = {
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV
    };
    static const uint8_t abytes[16] = {
        'A','A', 'C','A', /* + CDS-only: TAA -> TCA, stop_lost */
        'A','A', 'T','A', /* + CDS/3'UTR feature: coding_unknown */
        'T','T', 'T','G', /* - 3'UTR/CDS feature: coding_unknown */
        'T','T', 'G','T'  /* - CDS-only: TAA -> TAC, stop_lost */
    };
    static const uint32_t roff[4] = {0u, 4u, 8u, 12u};
    static const uint32_t aoff[4] = {2u, 6u, 10u, 14u};
    static const uint16_t rlen[4] = {2u, 2u, 2u, 2u};
    static const uint16_t alen[4] = {2u, 2u, 2u, 2u};
    static const uint64_t expected[4] = {
        DUCKVEP_SO(DUCKVEP_SO_STOP_LOST),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
            DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
            DUCKVEP_SO(DUCKVEP_SO_3_PRIME_UTR),
        DUCKVEP_SO(DUCKVEP_SO_STOP_LOST)
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_len; seq.codon_table = cds_tab;
    seq.post_cds_bases = post_cds; seq.transcript_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 4u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(4u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 4u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(expected[i], rows[i].consequence_mask);
        if (i == 1u || i == 2u) {
            ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_NOT_APPLICABLE,
                      rows[i].sequence_status);
            ASSERT_EQ(0u, rows[i].flags &
                          (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED);
        }
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_padded_small_variant_delta_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1014u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {15u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[15] = {
        'A','T','G', 'A','A','A', 'C','C','C', 'G','G','G', 'T','T','T'
    };
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {15u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    static const uint16_t vchrom[3] = {0u, 0u, 0u};
    static const uint32_t vpos[3]   = {1002u, 1002u, 1005u};
    static const uint32_t vend[3]   = {1005u, 1006u, 1006u};
    static const uint8_t  vkind[3]  = {
        (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_DEL,
        (uint8_t)DUCKVEP_KIND_INS
    };
    static const uint8_t  abytes[24] = {
        /* Complete feature spans codons 1-2 while only CDS5 changes. VEP sees
         * MK>MR, so start_lost and start_retained_variant are both true. */
        'G','A','A','A', 'G','A','G','A',
        'G','A','A','A','C', 'G','C',       /* GAAAC>GC deletes codon2 AAA */
        'A','C', 'A','G','C','C','C'        /* AC>AGCCC inserts GCC after codon2 */
    };
    static const uint32_t roff[3] = {0u, 8u, 15u};
    static const uint32_t aoff[3] = {4u, 13u, 17u};
    static const uint16_t rlen[3] = {4u, 5u, 2u};
    static const uint16_t alen[3] = {4u, 2u, 5u};
    static const uint64_t exp_mask[3] = {
        DUCKVEP_SO(DUCKVEP_SO_START_LOST) |
            DUCKVEP_SO(DUCKVEP_SO_START_RETAINED),
        DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION),
        DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION)
    };
    static const int32_t exp_protein[3] = {-1, 2, 3};
    static const int32_t exp_cdna[3] = {-1, -1, -1};
    static const int32_t exp_cds[3] = {-1, -1, -1};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 3u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 3u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(exp_cdna[i], rows[i].cdna_pos);
        ASSERT_EQ(exp_cds[i], rows[i].cds_pos);
        ASSERT_EQ(exp_protein[i], rows[i].protein_pos);
    }
    ASSERT_EQ((uint8_t)0u, rows[0].aa_ref);
    ASSERT_EQ((uint8_t)0u, rows[0].aa_alt);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_equal_length_feature_window_both_strands_known_scene(void) {
    static const uint16_t tchrom[2] = {0u, 1u};
    static const uint32_t tstart[2] = {1000u, 2000u};
    static const uint32_t tend[2] = {1014u, 2014u};
    static const int8_t tstrand[2] = {1, -1};
    static const uint64_t tflags[2] = {0u, 0u};
    static const uint32_t texoff[2] = {0u, 1u};
    static const uint16_t texcnt[2] = {1u, 1u};
    static const uint32_t tcds_s[2] = {1000u, 2000u};
    static const uint32_t tcds_e[2] = {1014u, 2014u};
    static const uint32_t estart[2] = {1000u, 2000u};
    static const uint32_t eend[2] = {1014u, 2014u};
    static const uint32_t ecdna_s[2] = {1u, 1u};
    static const uint32_t ecdna_e[2] = {15u, 15u};
    static const int8_t ephase[2] = {0, 0};
    static const uint8_t cds_bytes[30] = {
        'A','T','G', 'G','A','A', 'C','C','C', 'T','G','G', 'T','A','A',
        'A','T','G', 'G','A','A', 'C','C','C', 'T','G','G', 'T','A','A'
    };
    static const uint64_t cds_off[2] = {0u, 15u};
    static const uint32_t cds_len[2] = {15u, 15u};
    static const uint8_t cds_tab[2] = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };

    /* Each two-base feature differs at one base. Its paired one-base record has
     * the same semantic CDS edit but does not select the adjacent start/stop codon. */
    static const uint16_t vchrom[12] = {
        0u, 0u, 0u, 0u, 0u, 0u, 1u, 1u, 1u, 1u, 1u, 1u
    };
    static const uint32_t vpos[12] = {
        1002u, 1003u, 1011u, 1011u, 1011u, 1011u,
        2002u, 2002u, 2003u, 2003u, 2011u, 2011u
    };
    static const uint32_t vend[12] = {
        1003u, 1003u, 1011u, 1011u, 1012u, 1012u,
        2003u, 2003u, 2003u, 2003u, 2012u, 2011u
    };
    static const uint8_t vkind[12] = {
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV, (uint8_t)DUCKVEP_KIND_SNV
    };
    static const uint8_t abytes[36] = {
        'G','G', 'G','A', 'G','A', 'G','A', 'G','T',
        'G','T', 'A','T', 'G','T', 'T','T',
        'A','C', 'A','A', 'A','C', 'A','T', 'C','A', 'C','T',
        'C','C', 'T','C', 'C','T'
    };
    static const uint32_t roff[12] = {
        0u, 4u, 6u, 8u, 10u, 14u, 18u, 22u, 26u, 28u, 30u, 34u
    };
    static const uint32_t aoff[12] = {
        2u, 5u, 7u, 9u, 12u, 16u, 20u, 24u, 27u, 29u, 32u, 35u
    };
    static const uint16_t rlen[12] = {
        2u, 1u, 1u, 1u, 2u, 2u, 2u, 2u, 1u, 1u, 2u, 1u
    };
    static const uint16_t alen[12] = {
        2u, 1u, 1u, 1u, 2u, 2u, 2u, 2u, 1u, 1u, 2u, 1u
    };
    static const uint64_t expected[12] = {
        DUCKVEP_SO(DUCKVEP_SO_START_LOST) |
            DUCKVEP_SO(DUCKVEP_SO_START_RETAINED),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED),
        DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE),
        DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
        DUCKVEP_SO(DUCKVEP_SO_START_LOST) |
            DUCKVEP_SO(DUCKVEP_SO_START_RETAINED),
        DUCKVEP_SO(DUCKVEP_SO_MISSENSE)
    };
    static const int32_t expected_protein[12] = {
        -1, 2, 4, 4, -1, -1, -1, -1, 4, 4, -1, 2
    };
    static const uint8_t expected_ref_aa[12] = {
        0u, 'E', 'W', 'W', 0u, 0u, 0u, 0u, 'W', 'W', 0u, 'E'
    };
    static const uint8_t expected_alt_aa[12] = {
        0u, 'K', '*', 'C', 0u, 0u, 0u, 0u, 'C', '*', 0u, 'K'
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[12];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_len; seq.codon_table = cds_tab;
    seq.transcript_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 12u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, rows, 12u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(12u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 12u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(expected[i], rows[i].consequence_mask);
        ASSERT_EQ(expected_protein[i], rows[i].protein_pos);
        ASSERT_EQ(expected_ref_aa[i], rows[i].aa_ref);
        ASSERT_EQ(expected_alt_aa[i], rows[i].aa_alt);
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    ASSERT(stats != NULL);
    ASSERT_EQ(6u, stats->substitution_context);
    ASSERT_EQ(0u, stats->mnv_direct_fallback);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_indel_frameshift_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1008u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1008u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {9u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'A','A','A', 'T','T','T'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    static const uint16_t vchrom[7] = {0u, 0u, 0u, 0u, 0u, 0u, 0u};
    static const uint32_t vpos[7]   = {1000u, 1002u, 1003u, 1004u, 1006u, 1006u, 1006u};
    static const uint32_t vend[7]   = {1000u, 1005u, 1003u, 1005u, 1006u, 1007u, 1006u};
    static const uint8_t  vkind[7]  = {(uint8_t)DUCKVEP_KIND_INS, (uint8_t)DUCKVEP_KIND_DEL, (uint8_t)DUCKVEP_KIND_INS, (uint8_t)DUCKVEP_KIND_DEL, (uint8_t)DUCKVEP_KIND_INS, (uint8_t)DUCKVEP_KIND_DEL, (uint8_t)DUCKVEP_KIND_INS};
    static const uint8_t  abytes[23] = {'A', 'A','T', 'G','A','A','A', 'G', 'A', 'A','T', 'A','A', 'A', 'A', 'A','T', 'A','A', 'A', 'T', 'T','A'};
    static const uint32_t roff[7]   = {0u, 3u, 8u, 11u, 14u, 17u, 20u};
    static const uint32_t aoff[7]   = {1u, 7u, 9u, 13u, 15u, 19u, 21u};
    static const uint16_t rlen[7]   = {1u, 4u, 1u, 2u, 1u, 2u, 1u};
    static const uint16_t alen[7]   = {2u, 1u, 2u, 1u, 2u, 1u, 2u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 7u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(7u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT) |
                  DUCKVEP_SO(DUCKVEP_SO_START_LOST),
              rows[0].consequence_mask);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION), rows[1].consequence_mask);
    ASSERT_EQ(-1, rows[1].cds_pos);
    ASSERT_EQ(2, rows[1].protein_pos);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT), rows[2].consequence_mask);
    ASSERT_EQ(-1, rows[2].cds_pos);
    ASSERT_EQ(2, rows[2].protein_pos);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT), rows[3].consequence_mask);
    ASSERT_EQ(-1, rows[3].cds_pos);
    ASSERT_EQ(2, rows[3].protein_pos);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE), rows[4].consequence_mask); /* wrong INS REF */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE), rows[5].consequence_mask); /* wrong DEL REF */
    /* Final-codon +1 INS with the ATG start intact: the general CodingContext resolves
     * the frameshift the direct body-only restriction rejected at the terminal codon. */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT), rows[6].consequence_mask);
    ASSERT_EQ(-1, rows[6].cds_pos);
    ASSERT_EQ(3, rows[6].protein_pos);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_terminal_stop_frame_change_uses_transcript_tail(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {1000u};
    static const uint32_t tend[1] = {1011u};
    static const int8_t tstrand[1] = {1};
    static const uint32_t texoff[1] = {0u};
    static const uint16_t texcnt[1] = {1u};
    static const uint32_t tcds_s[1] = {1000u};
    static const uint32_t tcds_e[1] = {1008u};
    static const uint32_t estart[1] = {1000u};
    static const uint32_t eend[1] = {1011u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {12u};
    static const int8_t ephase[1] = {0};
    static const uint8_t cds_bytes[9] = {
        'A','T','G', 'A','A','A', 'T','A','A'
    };
    static const uint8_t post_cds[DUCKVEP_POST_CDS_BASE_COUNT] = {
        'A','C','G'
    };
    static const uint64_t cds_off[1] = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t cds_tab[1] = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2] = {1006u, 1006u};
    static const uint32_t vend[2] = {1007u, 1006u};
    static const uint8_t vkind[2] = {
        (uint8_t)DUCKVEP_KIND_DEL,
        (uint8_t)DUCKVEP_KIND_INS
    };
    static const uint8_t abytes[6] = {
        'T','A', 'T',
        'T', 'T','T'
    };
    static const uint32_t roff[2] = {0u, 3u};
    static const uint32_t aoff[2] = {2u, 4u};
    static const uint16_t rlen[2] = {2u, 1u};
    static const uint16_t alen[2] = {1u, 2u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq);
    memset(&v, 0, sizeof v);
    memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.post_cds_bases = post_cds; seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(2u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED), rows[0].consequence_mask);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED, rows[0].sequence_status);
    ASSERT_EQ((uint8_t)'*', rows[0].aa_ref);
    ASSERT_EQ((uint8_t)'*', rows[0].aa_alt);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_STOP_LOST), rows[1].consequence_mask);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED, rows[1].sequence_status);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);

    seq.post_cds_bases = NULL;
    model = NULL; opts = NULL; ws = NULL;
    memset(&err, 0, sizeof err);
    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 2u);
    v.count = 1u;
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE), rows[0].consequence_mask);
    ASSERT((rows[0].flags &
            (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED) != 0u);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_MISSING_TRANSCRIPT_TAIL,
              rows[0].sequence_status);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_delins_frameshift_known_scene(void) {
    static const uint16_t tchrom[2]  = {0u, 1u};
    static const uint32_t tstart[2]  = {1000u, 2000u};
    static const uint32_t tend[2]    = {1014u, 2014u};
    static const int8_t   tstrand[2] = {1, -1};
    static const uint64_t tflags[2]  = {0u, 0u};
    static const uint32_t texoff[2]  = {0u, 1u};
    static const uint16_t texcnt[2]  = {1u, 1u};
    static const uint32_t tcds_s[2]  = {1000u, 2000u};
    static const uint32_t tcds_e[2]  = {1014u, 2014u};
    static const uint32_t estart[2]  = {1000u, 2000u};
    static const uint32_t eend[2]    = {1014u, 2014u};
    static const uint32_t ecdna_s[2] = {1u, 1u};
    static const uint32_t ecdna_e[2] = {15u, 15u};
    static const int8_t   ephase[2]  = {0, 0};
    static const uint8_t  cds_bytes[30] = {
        'A','T','G', 'A','A','A', 'C','C','C', 'G','G','G', 'T','T','T',
        'A','T','G', 'A','A','A', 'C','C','C', 'G','G','G', 'T','T','T'
    };
    static const uint64_t cds_off[2]  = {0u, 15u};
    static const uint32_t cds_lenv[2] = {15u, 15u};
    static const uint8_t  cds_tab[2]  = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };

    static const uint16_t vchrom[4] = {0u, 0u, 0u, 1u};
    static const uint32_t vpos[4]   = {1003u, 1004u, 1004u, 2009u};
    static const uint32_t vend[4]   = {1005u, 1005u, 1005u, 2010u};
    static const uint8_t  vkind[4]  = {
        (uint8_t)DUCKVEP_KIND_INDEL,
        (uint8_t)DUCKVEP_KIND_INDEL,
        (uint8_t)DUCKVEP_KIND_INDEL,
        (uint8_t)DUCKVEP_KIND_INDEL
    };
    static const uint8_t  abytes[18] = {
        'A','A','A', 'A','G',          /* + padded: raw AAA>AG trims to AA>G */
        'A','A', 'G',                  /* +: cds5-6 AA>G, net -1 frameshift */
        'A','A', 'G','G','G','G','G',  /* +: net +3 delins, edit-set backlog */
        'T','T', 'C'                   /* -: genomic TT>C == transcript AA>G */
    };
    static const uint32_t roff[4] = {0u, 5u, 8u, 15u};
    static const uint32_t aoff[4] = {3u, 7u, 10u, 17u};
    static const uint16_t rlen[4] = {3u, 2u, 2u, 2u};
    static const uint16_t alen[4] = {2u, 1u, 5u, 1u};
    static const uint64_t exp_mask[4] = {
        DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT),
        DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT),
        DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING),
        DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT)
    };
    static const int32_t exp_protein[4] = {2, 2, 2, 2};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[5];
    duckvep_result_builder_t rb;
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 4u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 5u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(4u, duckvep_result_builder_count(&rb));
    for (i = 0u; i < 4u; i++) {
        ASSERT_EQ_FMT((uint32_t)i, rows[i].variant_idx, "%u");
        ASSERT_EQ(exp_mask[i], rows[i].consequence_mask);
        ASSERT_EQ(-1, rows[i].cdna_pos);
        ASSERT_EQ(-1, rows[i].cds_pos);
        ASSERT_EQ(exp_protein[i], rows[i].protein_pos);
    }

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_inframe_insertion_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1014u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {15u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[15] = {'A','T','G', 'A','A','A', 'C','C','C', 'G','G','G', 'T','T','T'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {15u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2]   = {1005u, 1005u};
    static const uint32_t vend[2]   = {1005u, 1005u};
    static const uint8_t  vkind[2]  = {(uint8_t)DUCKVEP_KIND_INS, (uint8_t)DUCKVEP_KIND_INS};
    static const uint8_t  abytes[10] = {'A', 'A','G','C','C', 'A', 'A','T','A','A'};
    static const uint32_t roff[2]   = {0u, 5u};
    static const uint32_t aoff[2]   = {1u, 6u};
    static const uint16_t rlen[2]   = {1u, 1u};
    static const uint16_t alen[2]   = {4u, 4u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 15u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(2u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION), rows[0].consequence_mask);
    ASSERT_EQ(-1, rows[0].cds_pos);
    ASSERT_EQ(3, rows[0].protein_pos);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION) |
                  DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
              rows[1].consequence_mask);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_inframe_insertion_reverse_known_scene(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1014u};
    static const int8_t   tstrand[1] = {-1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1014u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1014u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {15u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[15] = {'A','T','G', 'A','A','A', 'C','C','C', 'G','G','G', 'T','T','T'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {15u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {1008u}; /* transcript CDS 7 on the '-' strand */
    static const uint32_t vend[1]   = {1008u};
    static const uint8_t  vkind[1]  = {(uint8_t)DUCKVEP_KIND_INS};
    static const uint8_t  abytes[5] = {'G', 'G','G','G','C'}; /* ref G, inserted transcript GCC */
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {1u};
    static const uint16_t rlen[1]   = {1u};
    static const uint16_t alen[1]   = {4u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 15u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION), rows[0].consequence_mask);
    ASSERT_EQ(3, rows[0].protein_pos);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_protein_altering_insertion_known_scene(void) {
    static const uint16_t tchrom[2]  = {0u, 1u};
    static const uint32_t tstart[2]  = {1000u, 2000u};
    static const uint32_t tend[2]    = {1014u, 2014u};
    static const int8_t   tstrand[2] = {1, -1};
    static const uint64_t tflags[2]  = {0u, 0u};
    static const uint32_t texoff[2]  = {0u, 1u};
    static const uint16_t texcnt[2]  = {1u, 1u};
    static const uint32_t tcds_s[2]  = {1000u, 2000u};
    static const uint32_t tcds_e[2]  = {1014u, 2014u};
    static const uint32_t estart[2]  = {1000u, 2000u};
    static const uint32_t eend[2]    = {1014u, 2014u};
    static const uint32_t ecdna_s[2] = {1u, 1u};
    static const uint32_t ecdna_e[2] = {15u, 15u};
    static const int8_t   ephase[2]  = {0, 0};
    static const uint8_t  cds_bytes[30] = {
        'A','T','G', 'A','A','A', 'C','C','C', 'G','G','G', 'T','T','T',
        'A','T','G', 'A','A','A', 'C','C','C', 'G','G','G', 'T','T','T'
    };
    static const uint64_t cds_off[2]  = {0u, 15u};
    static const uint32_t cds_lenv[2] = {15u, 15u};
    static const uint8_t  cds_tab[2]  = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };

    static const uint16_t vchrom[3] = {0u, 0u, 1u};
    static const uint32_t vpos[3]   = {1003u, 1003u, 2010u};
    static const uint32_t vend[3]   = {1003u, 1003u, 2010u};
    static const uint8_t  vkind[3]  = {
        (uint8_t)DUCKVEP_KIND_INS,
        (uint8_t)DUCKVEP_KIND_INS,
        (uint8_t)DUCKVEP_KIND_INS
    };
    static const uint8_t  abytes[15] = {
        'A', 'A','G','C','C',  /* + strand: A + transcript GCC -> protein_altering */
        'A', 'A','C','C','T',  /* + strand: protein-altering junction + TAA stop */
        'T', 'T','G','G','C'   /* - strand: transcript GCC after CDS 4 */
    };
    static const uint32_t roff[3] = {0u, 5u, 10u};
    static const uint32_t aoff[3] = {1u, 6u, 11u};
    static const uint16_t rlen[3] = {1u, 1u, 1u};
    static const uint16_t alen[3] = {4u, 4u, 4u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[6];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 2u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 2u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 3u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 6u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(3u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING), rows[0].consequence_mask);
    ASSERT_EQ(2, rows[0].protein_pos);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING) |
                  DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
              rows[1].consequence_mask);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING), rows[2].consequence_mask);
    ASSERT_EQ(2, rows[2].protein_pos);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_codon_indel_reverse_ref_mismatch_falls_back(void) {
    static const uint16_t tchrom[1]  = {0u};
    static const uint32_t tstart[1]  = {1000u};
    static const uint32_t tend[1]    = {1008u};
    static const int8_t   tstrand[1] = {-1};
    static const uint32_t texoff[1]  = {0u};
    static const uint16_t texcnt[1]  = {1u};
    static const uint32_t tcds_s[1]  = {1000u};
    static const uint32_t tcds_e[1]  = {1008u};
    static const uint32_t estart[1]  = {1000u};
    static const uint32_t eend[1]    = {1008u};
    static const uint32_t ecdna_s[1] = {1u};
    static const uint32_t ecdna_e[1] = {9u};
    static const int8_t   ephase[1]  = {0};
    static const uint8_t  cds_bytes[9] = {'A','T','G', 'A','A','A', 'T','T','T'};
    static const uint64_t cds_off[1]  = {0u};
    static const uint32_t cds_lenv[1] = {9u};
    static const uint8_t  cds_tab[1]  = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};

    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {1004u};
    static const uint32_t vend[1]   = {1004u};
    static const uint8_t  vkind[1]  = {(uint8_t)DUCKVEP_KIND_INS};
    static const uint8_t  abytes[3] = {'A', 'A','T'}; /* wrong genomic REF at 1004; actual is T */
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {1u};
    static const uint16_t rlen[1]   = {1u};
    static const uint16_t alen[1]   = {2u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = texoff; tx.exon_count = texcnt;
    tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
    exons.start1 = estart; exons.end1 = eend;
    exons.cdna_start1 = ecdna_s; exons.cdna_end1 = ecdna_e;
    exons.phase = ephase; exons.end_phase = ephase; exons.exon_count = 1u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = 9u;
    seq.cds_offset = cds_off; seq.cds_length = cds_lenv; seq.codon_table = cds_tab;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = abytes; v.allele_bytes_len = sizeof abytes;
    v.ref_offset = roff; v.alt_offset = aoff;
    v.ref_length = rlen; v.alt_length = alen; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE), rows[0].consequence_mask);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_rejects_unsorted_variant_batch(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1]   = {200u};
    static const int8_t tstrand[1]  = {1};
    static const uint32_t tzero[1]  = {0u};
    static const uint16_t ezero[1]  = {0u};
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2]   = {160u, 150u};
    static const uint8_t vkind[2]   = {0u, 0u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t rows[4];
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = tzero; tx.exon_count = ezero;
    tx.cds_start1 = tzero; tx.cds_end1 = tzero; tx.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind; v.count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 4u);
    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(0u, duckvep_result_builder_count(&rb));

    duckvep_workspace_close(ws); duckvep_options_close(opts); duckvep_model_close(model);
    PASS();
}

TEST annotate_rejects_missing_alleles_for_nonpoint_small_variant(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1]   = {200u};
    static const int8_t tstrand[1]  = {1};
    static const uint64_t flags[1]  = {0u};
    static const uint32_t tzero[1]  = {0u};
    static const uint16_t ezero[1]  = {0u};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {150u};
    static const uint32_t vend[1]   = {151u};
    static const uint8_t vkind[1]   = {(uint8_t)DUCKVEP_KIND_DEL};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t row;
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = flags; tx.exon_offset = tzero; tx.exon_count = ezero;
    tx.cds_start1 = tzero; tx.cds_end1 = tzero; tx.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, &row, 1u);
    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(0u, duckvep_result_builder_count(&rb));

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_accepts_right_anchored_pure_insertion(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1]   = {200u};
    static const int8_t tstrand[1]  = {1};
    static const uint64_t flags[1]  = {0u};
    static const uint32_t tzero[1]  = {0u};
    static const uint16_t ezero[1]  = {0u};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {150u};
    static const uint32_t vend[1]   = {150u};
    static const uint8_t vkind[1]   = {(uint8_t)DUCKVEP_KIND_INS};
    static const uint8_t bytes[3]   = {'A', 'T', 'A'}; /* A>TA inserts before POS */
    static const uint32_t roff[1]   = {0u};
    static const uint16_t rlen[1]   = {1u};
    static const uint32_t aoff[1]   = {1u};
    static const uint16_t alen[1]   = {2u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t row;
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = flags; tx.exon_offset = tzero; tx.exon_count = ezero;
    tx.cds_start1 = tzero; tx.cds_end1 = tzero; tx.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = bytes; v.allele_bytes_len = sizeof bytes;
    v.ref_offset = roff; v.ref_length = rlen; v.alt_offset = aoff; v.alt_length = alen;
    v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, &row, 1u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb));

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_rejects_kind_allele_shape_mismatch(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1]   = {200u};
    static const int8_t tstrand[1]  = {1};
    static const uint64_t flags[1]  = {0u};
    static const uint32_t tzero[1]  = {0u};
    static const uint16_t ezero[1]  = {0u};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {150u};
    static const uint32_t vend[1]   = {150u};
    static const uint8_t vkind[1]   = {(uint8_t)DUCKVEP_KIND_SNV};
    static const uint8_t bytes[3]   = {'T', 'A', 'T'};
    static const uint32_t roff[1]   = {0u};
    static const uint16_t rlen[1]   = {2u};
    static const uint32_t aoff[1]   = {2u};
    static const uint16_t alen[1]   = {1u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t row;
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = flags; tx.exon_offset = tzero; tx.exon_count = ezero;
    tx.cds_start1 = tzero; tx.cds_end1 = tzero; tx.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.variant_kind = vkind;
    v.allele_bytes = bytes; v.allele_bytes_len = sizeof bytes;
    v.ref_offset = roff; v.ref_length = rlen; v.alt_offset = aoff; v.alt_length = alen;
    v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, &row, 1u);
    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(0u, duckvep_result_builder_count(&rb));

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_rejects_allele_slice_outside_pool(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1]   = {200u};
    static const int8_t tstrand[1]  = {1};
    static const uint32_t tzero[1]  = {0u};
    static const uint16_t ezero[1]  = {0u};
    static const uint64_t cds_off[1] = {0u};
    static const uint32_t cds_len[1] = {0u};
    static const uint8_t codon_table[1] = {(uint8_t)DUCKVEP_CODON_TABLE_STANDARD};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1]   = {150u};
    static const uint8_t vkind[1]   = {(uint8_t)DUCKVEP_KIND_SNV};
    static const uint8_t alleles[1] = {'A'};
    static const uint32_t roff[1]   = {0u};
    static const uint32_t aoff[1]   = {1u}; /* one past the pool */
    static const uint16_t one[1]    = {1u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t row;
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = tzero; tx.exon_count = ezero;
    tx.cds_start1 = tzero; tx.cds_end1 = tzero; tx.transcript_count = 1u;
    seq.cds_offset = cds_off; seq.cds_length = cds_len; seq.codon_table = codon_table;
    seq.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind;
    v.ref_offset = roff; v.ref_length = one; v.alt_offset = aoff; v.alt_length = one;
    v.allele_bytes = alleles; v.allele_bytes_len = sizeof alleles; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, &row, 1u);
    ASSERT_EQ(DUCKVEP_ERR_OUT_OF_RANGE,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(0u, duckvep_result_builder_count(&rb));

    duckvep_workspace_close(ws); duckvep_options_close(opts); duckvep_model_close(model);
    PASS();
}

/* A full result builder reports DUCKVEP_ERR_RESULT_FULL and never truncates a
 * partial row: two intronic pairs, capacity one. */
TEST annotate_result_full_is_reported(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1]   = {200u};
    static const int8_t   tstrand[1]= {1};
    static const uint32_t tcds[1]   = {0u};
    static const uint32_t teoff[1]  = {0u};
    static const uint16_t tecnt[1]  = {0u};
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2]   = {150u, 160u};
    static const uint8_t  vkind[2]  = {0u, 0u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[1];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = teoff; tx.exon_count = tecnt;
    tx.cds_start1 = tcds; tx.cds_end1 = tcds; tx.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind; v.count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 1u);
    ASSERT_EQ(DUCKVEP_ERR_RESULT_FULL,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(1u, duckvep_result_builder_count(&rb)); /* exactly the first row, no partial */

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

TEST annotate_rejects_result_count_past_capacity(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u}, tend[1] = {200u}, zero32[1] = {0u};
    static const uint16_t zero16[1] = {0u};
    static const int8_t tstrand[1] = {1};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1] = {150u};
    static const uint8_t vkind[1] = {(uint8_t)DUCKVEP_KIND_SNV};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_consequence_t row;
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = zero32; tx.exon_count = zero16;
    tx.cds_start1 = zero32; tx.cds_end1 = zero32; tx.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind; v.count = 1u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, &row, 1u);
    rb.count = 2u;
    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
              duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(2u, rb.count);

    duckvep_workspace_close(ws); duckvep_options_close(opts); duckvep_model_close(model);
    PASS();
}

/* Distances are uint32 coordinates. A signed narrowing would turn this 4-billion
 * base upstream distance negative and bypass the 3-billion-base directional cut. */
TEST annotate_directional_distance_uses_u32(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {UINT32_C(4000000001)};
    static const uint32_t tend[1] = {UINT32_C(4000000010)};
    static const uint32_t zero32[1] = {0u};
    static const uint16_t zero16[1] = {0u};
    static const int8_t tstrand[1] = {1};
    static const uint16_t vchrom[1] = {0u};
    static const uint32_t vpos[1] = {1u};
    static const uint8_t vkind[1] = {(uint8_t)DUCKVEP_KIND_SNV};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_options_init_t init;
    duckvep_consequence_t row;
    duckvep_result_builder_t rb;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&init, 0, sizeof init); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = zero32; tx.exon_count = zero16;
    tx.cds_start1 = zero32; tx.cds_end1 = zero32; tx.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind; v.count = 1u;
    init.upstream_dist = UINT32_C(3000000000);
    init.downstream_dist = 1u;
    init.halo = UINT32_MAX;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(&init, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, &row, 1u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));
    ASSERT_EQ(0u, rb.count);

    duckvep_workspace_close(ws); duckvep_options_close(opts); duckvep_model_close(model);
    PASS();
}

/* The symmetric sweep halo over-admits when up/downstream distances differ; the
 * directional filter must drop pairs beyond the per-direction window. halo 5000,
 * up = down = 100: variants 200bp away are dropped, 50bp away are kept. */
TEST annotate_directional_distance_filter(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {1000u};
    static const uint32_t tend[1]   = {2000u};
    static const int8_t   tstrand[1]= {1};
    static const uint32_t tcds[1]   = {0u};
    static const uint32_t teoff[1]  = {0u};
    static const uint16_t tecnt[1]  = {0u};
    static const uint16_t vchrom[4] = {0u, 0u, 0u, 0u};
    static const uint32_t vpos[4]   = {800u, 950u, 2050u, 2300u}; /* up200, up50, down50, down300 */
    static const uint8_t  vkind[4]  = {0u, 0u, 0u, 0u};

    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t v;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_options_init_t init;
    duckvep_error_t err;
    duckvep_consequence_t rows[8];
    duckvep_result_builder_t rb;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&v, 0, sizeof v); memset(&err, 0, sizeof err);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
    tx.flags = k_zero_flags; tx.exon_offset = teoff; tx.exon_count = tecnt;
    tx.cds_start1 = tcds; tx.cds_end1 = tcds; tx.transcript_count = 1u;
    v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vpos; v.variant_kind = vkind; v.count = 4u;

    memset(&init, 0, sizeof init);
    init.upstream_dist = 100u; init.downstream_dist = 100u; init.halo = 5000u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(&init, &opts, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    duckvep_result_builder_init(&rb, rows, 8u);
    ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &v, opts, ws, &rb, &err));

    ASSERT_EQ(2u, duckvep_result_builder_count(&rb));
    ASSERT_EQ_FMT(1u, rows[0].variant_idx, "%u"); /* up50 kept */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_UPSTREAM_GENE), rows[0].consequence_mask);
    ASSERT_EQ_FMT(2u, rows[1].variant_idx, "%u"); /* down50 kept */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_DOWNSTREAM_GENE), rows[1].consequence_mask);

    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    PASS();
}

/* model_open validates once and rejects malformed models with STABLE where_codes
 * (golden anchors: 13 exon-range, 14 cds-range, 15 unsorted — see duckvep_kernel.c). */
TEST model_open_rejects_invalid_models(void) {
    duckvep_model_t *model = NULL;
    duckvep_error_t err;
    duckvep_exon_model_t exons;

#if SIZE_MAX > UINT32_MAX
    /* tx_idx is uint32 throughout the ABI; reject an oversized borrowed view
     * before touching any of its one-element sentinel columns. */
    {
        static const uint16_t c[1] = {0u};
        static const uint32_t z32[1] = {0u};
        static const int8_t st[1] = {1};
        static const uint16_t z16[1] = {0u};
        duckvep_transcript_model_t tx;
        memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons); memset(&err, 0, sizeof err);
        tx.chrom_id = c; tx.start1 = z32; tx.end1 = z32; tx.strand = st;
        tx.flags = k_zero_flags; tx.exon_offset = z32; tx.exon_count = z16;
        tx.cds_start1 = z32; tx.cds_end1 = z32;
        tx.transcript_count = (size_t)UINT32_MAX + 1u;
        ASSERT_EQ(DUCKVEP_ERR_OUT_OF_RANGE,
                  duckvep_model_open(&tx, &exons, NULL, &model, &err));
        ASSERT_EQ(18u, err.where_code);
        ASSERT(model == NULL);
    }
#endif

    /* exon slice out of range */
    {
        static const uint16_t c[1] = {0u};
        static const uint32_t s1[1] = {100u}, e1[1] = {200u};
        static const int8_t st[1] = {1};
        static const uint32_t eoff[1] = {0u}; static const uint16_t ecnt[1] = {5u};
        static const uint32_t cz[1] = {0u};
        static const uint32_t es[2] = {100u, 150u}, ee[2] = {120u, 200u};
        duckvep_transcript_model_t tx;
        memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons); memset(&err, 0, sizeof err);
        tx.chrom_id = c; tx.start1 = s1; tx.end1 = e1; tx.strand = st;
        tx.flags = k_zero_flags; tx.exon_offset = eoff; tx.exon_count = ecnt; tx.cds_start1 = cz; tx.cds_end1 = cz;
        tx.transcript_count = 1u;
        exons.start1 = es; exons.end1 = ee; exons.exon_count = 2u; /* slice [0,5) exceeds 2 */
        ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID, duckvep_model_open(&tx, &exons, NULL, &model, &err));
        ASSERT_EQ(13u, err.where_code);
        ASSERT_EQ(NULL, model);
    }
    /* cds outside the transcript span */
    {
        static const uint16_t c[1] = {0u};
        static const uint32_t s1[1] = {100u}, e1[1] = {200u};
        static const int8_t st[1] = {1};
        static const uint32_t eoff[1] = {0u}; static const uint16_t ecnt[1] = {0u};
        static const uint32_t cs[1] = {50u}, ce[1] = {150u}; /* 50 < start1 100 */
        duckvep_transcript_model_t tx;
        memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons); memset(&err, 0, sizeof err);
        tx.chrom_id = c; tx.start1 = s1; tx.end1 = e1; tx.strand = st;
        tx.flags = k_zero_flags; tx.exon_offset = eoff; tx.exon_count = ecnt; tx.cds_start1 = cs; tx.cds_end1 = ce;
        tx.transcript_count = 1u;
        ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID, duckvep_model_open(&tx, &exons, NULL, &model, &err));
        ASSERT_EQ(14u, err.where_code);
    }
    /* transcripts not sorted by (chrom_id, start1) */
    {
        static const uint16_t c[2] = {0u, 0u};
        static const uint32_t s1[2] = {200u, 100u}; /* descending start on same chrom */
        static const uint32_t e1[2] = {300u, 150u};
        static const int8_t st[2] = {1, 1};
        static const uint32_t eoff[2] = {0u, 0u}; static const uint16_t ecnt[2] = {0u, 0u};
        static const uint32_t cz[2] = {0u, 0u};
        duckvep_transcript_model_t tx;
        memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons); memset(&err, 0, sizeof err);
        tx.chrom_id = c; tx.start1 = s1; tx.end1 = e1; tx.strand = st;
        tx.flags = k_zero_flags; tx.exon_offset = eoff; tx.exon_count = ecnt; tx.cds_start1 = cz; tx.cds_end1 = cz;
        tx.transcript_count = 2u;
        ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID, duckvep_model_open(&tx, &exons, NULL, &model, &err));
        ASSERT_EQ(15u, err.where_code);
    }
    /* NULL flags with transcript_count>0 is rejected (distilled tx_flags are mandatory;
     * the kernel will read them for biotype terms — pin the invariant now). */
    {
        static const uint16_t c[1] = {0u};
        static const uint32_t s1[1] = {100u}, e1[1] = {200u};
        static const int8_t st[1] = {1};
        static const uint32_t eoff[1] = {0u}; static const uint16_t ecnt[1] = {0u};
        static const uint32_t cz[1] = {0u};
        duckvep_transcript_model_t tx;
        memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons); memset(&err, 0, sizeof err);
        tx.chrom_id = c; tx.start1 = s1; tx.end1 = e1; tx.strand = st;
        tx.flags = NULL; /* the invariant under test */
        tx.exon_offset = eoff; tx.exon_count = ecnt; tx.cds_start1 = cz; tx.cds_end1 = cz;
        tx.transcript_count = 1u;
        ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG, duckvep_model_open(&tx, &exons, NULL, &model, &err));
        ASSERT_EQ(11u, err.where_code); /* DVW_MODEL_NULL_VIEW */
        ASSERT(model == NULL);
    }
    PASS();
}

/* ===================================================================== *
 * Codon-bucket FUSION vs the coding-SNV kernel oracle.
 *
 * For a generated single-exon '+' coding transcript (random CDS sequence) and a
 * random single-base SNV at a non-boundary CDS position, annotate_tile's emitted
 * row must EQUAL the result of calling the tested coding-SNV kernel directly:
 * same SO term, same cds_pos / protein_pos / aa. The SNV is kept >=3 bases inside
 * each exon end so no splice_region bit is involved (that interaction is anchored
 * separately). This proves the fused codon path == sweep . classify . coding-SNV
 * composition across all codon-change classes (syn/missense/stop_gained/lost) at
 * scale, not just the deterministic anchor.
 * ===================================================================== */

#define KPROP_MAX_CODONS 8u

struct kprop_coding {
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t       ex;
    duckvep_sequence_pool_t    seq;
    duckvep_variant_batch_t    v;
    uint16_t chrom; uint32_t tstart; uint32_t tend; int8_t strand; uint64_t flags;
    uint32_t exoff; uint16_t excnt; uint32_t cds_s; uint32_t cds_e;
    uint32_t es; uint32_t ee; uint32_t ecds; uint32_t ecde; int8_t eph; int8_t eeph;
    uint8_t *cds; uint64_t cds_off0; uint32_t cds_lenv; uint8_t ctab;
    uint16_t vchrom; uint32_t vpos; uint32_t vend; uint8_t vkind;
    uint8_t abytes[16]; uint32_t roff; uint32_t aoff; uint16_t rlen; uint16_t alen;
    uint8_t expect_cds[64]; uint32_t expect_len; uint8_t expect_shape; uint8_t expect_region;
    uint8_t expect_protein_altering;
};

static void kprop_coding_free(void *instance, void *env) {
    struct kprop_coding *s = (struct kprop_coding *)instance;
    (void)env;
    if (s == NULL) return;
    free(s->cds);
    free(s);
}

static enum theft_alloc_res kprop_coding_alloc(struct theft *t, void *env, void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 2u) + 3u; /* 3..8 */
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t ci, ri, ai, i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];

    /* SNV at a non-boundary CDS index: [3 .. cds_len-4] (cds_len>=9). */
    ci = (uint32_t)kprop_bounded(t, cds_len - 6u) + 3u;
    ri = 0u; while (BASES[ri] != (char)s->cds[ci]) ri++;          /* ref base index */
    ai = (ri + 1u + (uint32_t)kprop_bounded(t, 3u)) % 4u;         /* a different base */

    s->chrom = 0u; s->strand = 1;
    s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u;

    s->vchrom = 0u; s->vpos = base + ci; s->vend = s->vpos; s->vkind = (uint8_t)DUCKVEP_KIND_SNV;
    s->abytes[0] = s->cds[ci]; s->abytes[1] = (uint8_t)BASES[ai];
    s->roff = 0u; s->aoff = 1u; s->rlen = 1u; s->alen = 1u;

    s->cds_off0 = 0u; s->cds_lenv = cds_len; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;

    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = sizeof s->abytes;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_coding_info = {
    .alloc = kprop_coding_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_start_codon_alloc(struct theft *t, void *env, void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 2u) + 3u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t ci, ri, ai, i, tries;
    char first[4];
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }

    first[3] = '\0';
    for (tries = 0u; tries < 32u; tries++) {
        first[0] = BASES[kprop_bounded(t, 4u)];
        first[1] = BASES[kprop_bounded(t, 4u)];
        first[2] = BASES[kprop_bounded(t, 4u)];
        if (duckvep_translate_codon(first, DUCKVEP_CODON_TABLE_STANDARD) != '*') break;
    }
    if (tries == 32u) { first[0] = 'A'; first[1] = 'T'; first[2] = 'G'; }
    s->cds[0] = (uint8_t)first[0]; s->cds[1] = (uint8_t)first[1]; s->cds[2] = (uint8_t)first[2];
    for (i = 3u; i < cds_len; i++) s->cds[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];

    ci = (uint32_t)kprop_bounded(t, 3u); /* a base inside the annotated start codon */
    ri = 0u; while (BASES[ri] != (char)s->cds[ci]) ri++;
    ai = (ri + 1u + (uint32_t)kprop_bounded(t, 3u)) % 4u;

    s->chrom = 0u; s->strand = 1; s->flags = 0u;
    s->tstart = base - 20u; s->tend = base + cds_len + 19u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = s->tstart; s->ee = s->tend;
    s->ecds = 1u; s->ecde = cds_len + 40u; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u;

    s->vchrom = 0u; s->vpos = base + ci; s->vend = s->vpos; s->vkind = (uint8_t)DUCKVEP_KIND_SNV;
    s->abytes[0] = s->cds[ci]; s->abytes[1] = (uint8_t)BASES[ai];
    s->roff = 0u; s->aoff = 1u; s->rlen = 1u; s->alen = 1u;

    s->cds_off0 = 0u; s->cds_lenv = cds_len; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = sizeof s->abytes;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_start_codon_info = {
    .alloc = kprop_start_codon_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_mnv_alloc(struct theft *t, void *env, void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 2u) + 3u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t codon = (uint32_t)kprop_bounded(t, ncodons);
    uint32_t mnv_len = (uint32_t)kprop_bounded(t, 2u) + 2u;
    uint32_t codon_off = mnv_len == 3u ? 0u : (uint32_t)kprop_bounded(t, 2u);
    uint32_t ci = codon * 3u + codon_off;
    uint32_t i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];

    s->chrom = 0u; s->strand = 1; s->flags = 0u;
    s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u;

    s->vchrom = 0u; s->vpos = base + ci; s->vend = s->vpos + mnv_len - 1u;
    s->vkind = (uint8_t)DUCKVEP_KIND_MNV;
    for (i = 0u; i < mnv_len; i++) s->abytes[i] = s->cds[ci + i];
    for (i = 0u; i < mnv_len; i++) {
        uint32_t ri = 0u;
        uint32_t ai;
        while (BASES[ri] != (char)s->abytes[i]) ri++;
        ai = (ri + 1u + (uint32_t)kprop_bounded(t, 3u)) % 4u;
        s->abytes[mnv_len + i] = (uint8_t)BASES[ai];
    }
    s->roff = 0u; s->aoff = mnv_len; s->rlen = (uint16_t)mnv_len; s->alen = (uint16_t)mnv_len;

    s->cds_off0 = 0u; s->cds_lenv = cds_len; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = (size_t)mnv_len * 2u;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_mnv_info = {
    .alloc = kprop_mnv_alloc,
    .free  = kprop_coding_free,
};

static char kprop_complement_base(char b) {
    switch (b) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

static uint32_t kprop_cds_pos_for_genomic(const struct kprop_coding *s, uint32_t pos1) {
    uint32_t off0 = pos1 - s->es;
    return s->strand > 0 ? off0 + 1u : s->cds_lenv - off0;
}

static uint32_t kprop_genomic_pos_for_cds(const struct kprop_coding *s, uint32_t cds_pos) {
    uint32_t off0 = s->strand > 0 ? cds_pos - 1u : s->cds_lenv - cds_pos;
    return s->es + off0;
}

static char kprop_genomic_base_at(const struct kprop_coding *s, uint32_t pos1) {
    uint32_t cds_pos = kprop_cds_pos_for_genomic(s, pos1);
    char b = (char)s->cds[cds_pos - 1u];
    return s->strand > 0 ? b : kprop_complement_base(b);
}

static uint8_t kprop_base_not(uint8_t a, uint8_t b) {
    static const uint8_t BASES[4] = {'A', 'C', 'G', 'T'};
    uint32_t i;
    for (i = 0u; i < 4u; i++) {
        if (BASES[i] != a && BASES[i] != b) return BASES[i];
    }
    return (uint8_t)'N';
}

static uint8_t kprop_base_not3(uint8_t a, uint8_t b, uint8_t c) {
    static const uint8_t BASES[4] = {'A', 'C', 'G', 'T'};
    uint32_t i;
    for (i = 0u; i < 4u; i++) {
        if (BASES[i] != a && BASES[i] != b && BASES[i] != c) return BASES[i];
    }
    return (uint8_t)'N';
}

enum {
    KPROP_CDS_EDIT_SNV = 0u,
    KPROP_CDS_EDIT_MNV = 1u,
    KPROP_CDS_EDIT_INS = 2u,
    KPROP_CDS_EDIT_DEL = 3u,
    KPROP_CDS_EDIT_INDEL = 4u
};

enum {
    KPROP_CDS_EDIT_START = 0u,
    KPROP_CDS_EDIT_BODY = 1u,
    KPROP_CDS_EDIT_STOP = 2u
};

enum {
    KPROP_CONTEXT_DELTA_SYNONYMOUS = 0u,
    KPROP_CONTEXT_DELTA_MISSENSE = 1u,
    KPROP_CONTEXT_DELTA_STOP_GAINED = 2u,
    KPROP_CONTEXT_DELTA_STOP_LOST = 3u,
    KPROP_CONTEXT_DELTA_STOP_RETAINED = 4u
};

static uint32_t kprop_pick_cds_start(struct theft *t, uint32_t cds_len,
                                     uint32_t ref_len, uint32_t region) {
    uint32_t max_start = cds_len - ref_len + 1u;
    if (region == KPROP_CDS_EDIT_START) return 1u;
    if (region == KPROP_CDS_EDIT_STOP) return max_start;
    return (uint32_t)kprop_bounded(t, max_start - 6u) + 4u;
}

static void kprop_fill_expected_cds(struct kprop_coding *s, uint32_t cds_start,
                                    uint32_t ref_len, const uint8_t *alt_tx,
                                    uint32_t alt_len) {
    uint32_t prefix = cds_start - 1u;
    uint32_t suffix_start = prefix + ref_len;
    uint32_t suffix_len = s->cds_lenv - suffix_start;
    memcpy(s->expect_cds, s->cds, prefix);
    if (alt_len > 0u) memcpy(s->expect_cds + prefix, alt_tx, alt_len);
    memcpy(s->expect_cds + prefix + alt_len, s->cds + suffix_start, suffix_len);
    s->expect_len = prefix + alt_len + suffix_len;
}

static void kprop_fill_variant_alt_from_tx(struct kprop_coding *s, uint32_t off,
                                           const uint8_t *alt_tx, uint32_t alt_len) {
    uint32_t i;
    for (i = 0u; i < alt_len; i++) {
        uint8_t b = s->strand > 0 ? alt_tx[i]
            : (uint8_t)kprop_complement_base((char)alt_tx[alt_len - 1u - i]);
        s->abytes[off + i] = b;
    }
}

static void kprop_wire_coding_scene(struct kprop_coding *s, uint32_t cds_len) {
    s->cds_off0 = 0u; s->cds_lenv = cds_len; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = sizeof s->abytes;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;
}

static enum theft_alloc_res kprop_cds_edit_builder_alloc(struct theft *t, void *env,
                                                         void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 4u) + 5u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t shape = (uint32_t)kprop_bounded(t, 5u);
    uint32_t region = (uint32_t)kprop_bounded(t, 3u);
    uint32_t ref_len = 1u;
    uint32_t alt_len = 1u;
    uint32_t cds_start = 1u;
    uint8_t alt_tx[8];
    uint32_t i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? (int8_t)1 : (int8_t)-1;
    s->flags = 0u; s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u; s->vchrom = 0u;
    kprop_wire_coding_scene(s, cds_len);

    s->expect_shape = (uint8_t)shape;
    s->expect_region = (uint8_t)region;
    if (shape == KPROP_CDS_EDIT_SNV) {
        ref_len = 1u; alt_len = 1u;
        cds_start = kprop_pick_cds_start(t, cds_len, ref_len, region);
        alt_tx[0] = kprop_base_not(s->cds[cds_start - 1u], (uint8_t)'N');
        s->vkind = (uint8_t)DUCKVEP_KIND_SNV;
        s->vpos = kprop_genomic_pos_for_cds(s, cds_start); s->vend = s->vpos;
        s->abytes[0] = (uint8_t)kprop_genomic_base_at(s, s->vpos);
        kprop_fill_variant_alt_from_tx(s, 1u, alt_tx, alt_len);
        s->roff = 0u; s->aoff = 1u; s->rlen = 1u; s->alen = 1u;
        kprop_fill_expected_cds(s, cds_start, ref_len, alt_tx, alt_len);
    } else if (shape == KPROP_CDS_EDIT_MNV || shape == KPROP_CDS_EDIT_INDEL) {
        ref_len = (uint32_t)kprop_bounded(t, 3u) + 1u;
        alt_len = shape == KPROP_CDS_EDIT_MNV
            ? ref_len
            : (uint32_t)kprop_bounded(t, 4u) + 1u;
        if (shape == KPROP_CDS_EDIT_INDEL && alt_len == ref_len) alt_len = alt_len == 4u ? 1u : alt_len + 1u;
        cds_start = kprop_pick_cds_start(t, cds_len, ref_len, region);
        {
            uint8_t b0 = s->cds[cds_start - 1u];
            uint8_t b1 = ref_len > 1u ? s->cds[cds_start] : b0;
            uint8_t b2 = ref_len > 2u ? s->cds[cds_start + 1u] : b1;
            uint8_t alt_base = kprop_base_not3(b0, b1, b2);
            for (i = 0u; i < alt_len; i++) alt_tx[i] = alt_base;
        }
        s->vkind = shape == KPROP_CDS_EDIT_MNV ? (uint8_t)DUCKVEP_KIND_MNV : (uint8_t)DUCKVEP_KIND_INDEL;
        s->vpos = s->strand > 0 ? kprop_genomic_pos_for_cds(s, cds_start)
                                : kprop_genomic_pos_for_cds(s, cds_start + ref_len - 1u);
        s->vend = s->vpos + ref_len - 1u;
        for (i = 0u; i < ref_len; i++) s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, s->vpos + i);
        kprop_fill_variant_alt_from_tx(s, ref_len, alt_tx, alt_len);
        s->roff = 0u; s->aoff = ref_len; s->rlen = (uint16_t)ref_len; s->alen = (uint16_t)alt_len;
        kprop_fill_expected_cds(s, cds_start, ref_len, alt_tx, alt_len);
    } else if (shape == KPROP_CDS_EDIT_INS) {
        uint32_t anchor_cds = region == KPROP_CDS_EDIT_START ? 1u
            : (region == KPROP_CDS_EDIT_STOP ? cds_len
                                             : (uint32_t)kprop_bounded(t, cds_len - 6u) + 4u);
        uint32_t insert_cds = s->strand > 0 ? anchor_cds + 1u : anchor_cds;
        alt_len = (uint32_t)kprop_bounded(t, 3u) + 1u;
        for (i = 0u; i < alt_len; i++) alt_tx[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];
        s->vkind = (uint8_t)DUCKVEP_KIND_INS;
        s->vpos = kprop_genomic_pos_for_cds(s, anchor_cds); s->vend = s->vpos;
        s->abytes[0] = (uint8_t)kprop_genomic_base_at(s, s->vpos);
        s->abytes[1] = s->abytes[0];
        kprop_fill_variant_alt_from_tx(s, 2u, alt_tx, alt_len);
        s->roff = 0u; s->aoff = 1u; s->rlen = 1u; s->alen = (uint16_t)(alt_len + 1u);
        kprop_fill_expected_cds(s, insert_cds, 0u, alt_tx, alt_len);
    } else {
        uint32_t anchor_cds;
        ref_len = (uint32_t)kprop_bounded(t, 3u) + 1u;
        cds_start = (uint32_t)kprop_bounded(t, cds_len - ref_len - 1u) + 2u;
        if (s->strand > 0) {
            anchor_cds = cds_start - 1u;
        } else {
            if (cds_start + ref_len > cds_len) cds_start = cds_len - ref_len;
            anchor_cds = cds_start + ref_len;
        }
        s->expect_region = cds_start <= 3u ? (uint8_t)KPROP_CDS_EDIT_START
            : (cds_start + ref_len - 1u > cds_len - 3u ? (uint8_t)KPROP_CDS_EDIT_STOP
                                                        : (uint8_t)KPROP_CDS_EDIT_BODY);
        s->vkind = (uint8_t)DUCKVEP_KIND_DEL;
        s->vpos = kprop_genomic_pos_for_cds(s, anchor_cds);
        s->vend = s->vpos + ref_len;
        for (i = 0u; i <= ref_len; i++) s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, s->vpos + i);
        s->abytes[ref_len + 1u] = s->abytes[0];
        s->roff = 0u; s->aoff = ref_len + 1u; s->rlen = (uint16_t)(ref_len + 1u); s->alen = 1u;
        kprop_fill_expected_cds(s, cds_start, ref_len, alt_tx, 0u);
    }

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_cds_edit_builder_info = {
    .alloc = kprop_cds_edit_builder_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_cds_edit_set_mnv_alloc(struct theft *t, void *env,
                                                         void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 4u) + 5u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t len = (uint32_t)kprop_bounded(t, 5u) + 3u;
    uint32_t region = (uint32_t)kprop_bounded(t, 3u);
    uint32_t cds_start;
    uint8_t alt_tx[8];
    uint32_t i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? (int8_t)1 : (int8_t)-1;
    s->flags = 0u; s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u; s->vchrom = 0u;
    kprop_wire_coding_scene(s, cds_len);

    cds_start = kprop_pick_cds_start(t, cds_len, len, region);
    s->expect_shape = KPROP_CDS_EDIT_MNV;
    s->expect_region = (uint8_t)region;
    for (i = 0u; i < len; i++) {
        uint8_t ref = s->cds[cds_start - 1u + i];
        int diff = (i == 0u || i + 1u == len) ? 1 : (kprop_bounded(t, 2u) == 0u);
        if (i == 1u) diff = 0; /* force at least one retained internal base */
        alt_tx[i] = diff ? kprop_base_not(ref, (uint8_t)'N') : ref;
    }
    s->vkind = (uint8_t)DUCKVEP_KIND_MNV;
    s->vpos = s->strand > 0 ? kprop_genomic_pos_for_cds(s, cds_start)
                            : kprop_genomic_pos_for_cds(s, cds_start + len - 1u);
    s->vend = s->vpos + len - 1u;
    for (i = 0u; i < len; i++) s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, s->vpos + i);
    kprop_fill_variant_alt_from_tx(s, len, alt_tx, len);
    s->roff = 0u; s->aoff = len; s->rlen = (uint16_t)len; s->alen = (uint16_t)len;
    kprop_fill_expected_cds(s, cds_start, len, alt_tx, len);
    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_cds_edit_set_mnv_info = {
    .alloc = kprop_cds_edit_set_mnv_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_context_delta_alloc(struct theft *t, void *env,
                                                      void **instance) {
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 4u) + 5u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t category = (uint32_t)kprop_bounded(t, 5u);
    uint32_t codon_idx = (uint32_t)kprop_bounded(t, ncodons - 1u) + 1u;
    uint32_t cds_start = codon_idx * 3u + 1u;
    const uint8_t *ref_tx;
    const uint8_t *alt_tx;
    uint32_t i;
    static const uint8_t syn_ref[3] = {'G','A','A'};
    static const uint8_t syn_alt[3] = {'G','A','G'};
    static const uint8_t mis_ref[3] = {'G','A','A'};
    static const uint8_t mis_alt[3] = {'G','A','C'};
    static const uint8_t sg_ref[3] = {'T','G','G'};
    static const uint8_t sg_alt[3] = {'T','G','A'};
    static const uint8_t sl_ref[3] = {'T','A','A'};
    static const uint8_t sl_alt[3] = {'C','A','A'};
    static const uint8_t sr_ref[3] = {'T','A','A'};
    static const uint8_t sr_alt[3] = {'T','A','G'};
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    s->cds[0] = 'A'; s->cds[1] = 'T'; s->cds[2] = 'G';
    for (i = 3u; i < cds_len; i += 3u) {
        s->cds[i] = 'G'; s->cds[i + 1u] = 'A'; s->cds[i + 2u] = 'A';
    }
    if (category == KPROP_CONTEXT_DELTA_SYNONYMOUS) {
        ref_tx = syn_ref; alt_tx = syn_alt;
    } else if (category == KPROP_CONTEXT_DELTA_MISSENSE) {
        ref_tx = mis_ref; alt_tx = mis_alt;
    } else if (category == KPROP_CONTEXT_DELTA_STOP_GAINED) {
        ref_tx = sg_ref; alt_tx = sg_alt;
    } else if (category == KPROP_CONTEXT_DELTA_STOP_LOST) {
        ref_tx = sl_ref; alt_tx = sl_alt;
    } else {
        ref_tx = sr_ref; alt_tx = sr_alt;
    }
    memcpy(s->cds + (size_t)cds_start - 1u, ref_tx, 3u);

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? (int8_t)1 : (int8_t)-1;
    s->flags = 0u; s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u; s->vchrom = 0u;
    kprop_wire_coding_scene(s, cds_len);

    s->expect_shape = KPROP_CDS_EDIT_MNV;
    s->expect_region = (uint8_t)category;
    s->vkind = (uint8_t)DUCKVEP_KIND_MNV;
    s->vpos = s->strand > 0 ? kprop_genomic_pos_for_cds(s, cds_start)
                            : kprop_genomic_pos_for_cds(s, cds_start + 2u);
    s->vend = s->vpos + 2u;
    for (i = 0u; i < 3u; i++) s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, s->vpos + i);
    kprop_fill_variant_alt_from_tx(s, 3u, alt_tx, 3u);
    s->roff = 0u; s->aoff = 3u; s->rlen = 3u; s->alen = 3u;
    kprop_fill_expected_cds(s, cds_start, 3u, alt_tx, 3u);
    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_context_delta_info = {
    .alloc = kprop_context_delta_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_cross_codon_mnv_alloc(struct theft *t, void *env, void **instance) {
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 4u) + 5u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t codon_idx;
    uint32_t codon_start;
    uint32_t first_cds;
    uint32_t mode;
    uint32_t len;
    char alt_tx[3];
    uint32_t i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)'A';

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? 1 : -1; s->flags = 0u;
    s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u;
    s->cds_lenv = cds_len;

    codon_idx = (uint32_t)kprop_bounded(t, ncodons - 3u) + 1u; /* body codon, plus next body codon */
    codon_start = codon_idx * 3u + 1u;
    mode = (uint32_t)kprop_bounded(t, 4u);
    if (mode < 2u) {
        /* AAA/TTA, editing last base of codon1 + first of codon2.
         * mode 0: one synonymous codon and one missense codon.
         * mode 1: both codons synonymous, so this narrow slice must fall back. */
        s->cds[(size_t)codon_start - 1u] = (uint8_t)'A';
        s->cds[(size_t)codon_start] = (uint8_t)'A';
        s->cds[(size_t)codon_start + 1u] = (uint8_t)'A';
        s->cds[(size_t)codon_start + 2u] = (uint8_t)'T';
        s->cds[(size_t)codon_start + 3u] = (uint8_t)'T';
        s->cds[(size_t)codon_start + 4u] = (uint8_t)'A';
        first_cds = codon_start + 2u;
        len = 2u;
        alt_tx[0] = 'G';
        alt_tx[1] = mode == 0u ? 'G' : 'C';
    } else {
        /* TCA/AAA, editing last two bases of codon1 + first of codon2.
         * mode 2 introduces TAG stop and must fall back; mode 3 is non-stop missense. */
        s->cds[(size_t)codon_start - 1u] = (uint8_t)'T';
        s->cds[(size_t)codon_start] = (uint8_t)'C';
        s->cds[(size_t)codon_start + 1u] = (uint8_t)'A';
        s->cds[(size_t)codon_start + 2u] = (uint8_t)'A';
        s->cds[(size_t)codon_start + 3u] = (uint8_t)'A';
        s->cds[(size_t)codon_start + 4u] = (uint8_t)'A';
        first_cds = codon_start + 1u;
        len = 3u;
        if (mode == 2u) {
            alt_tx[0] = 'A'; alt_tx[1] = 'G'; alt_tx[2] = 'C';
        } else {
            alt_tx[0] = 'G'; alt_tx[1] = 'G'; alt_tx[2] = 'C';
        }
    }

    s->vchrom = 0u;
    s->vpos = s->strand > 0 ? kprop_genomic_pos_for_cds(s, first_cds)
                            : kprop_genomic_pos_for_cds(s, first_cds + len - 1u);
    s->vend = s->vpos + len - 1u;
    s->vkind = (uint8_t)DUCKVEP_KIND_MNV;
    for (i = 0u; i < len; i++) {
        uint32_t gpos = s->vpos + i;
        uint32_t cds_pos = kprop_cds_pos_for_genomic(s, gpos);
        uint32_t tx_i = cds_pos - first_cds;
        char alt = alt_tx[tx_i];
        s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, gpos);
        s->abytes[len + i] = (uint8_t)(s->strand > 0 ? alt : kprop_complement_base(alt));
    }
    s->roff = 0u; s->aoff = len; s->rlen = (uint16_t)len; s->alen = (uint16_t)len;

    s->cds_off0 = 0u; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = (size_t)len * 2u;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_cross_codon_mnv_info = {
    .alloc = kprop_cross_codon_mnv_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_frameshift_indel_alloc(struct theft *t, void *env, void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    static const uint8_t STOPS[3][3] = {
        {'T','A','A'}, {'T','A','G'}, {'T','G','A'}
    };
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 2u) + 3u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t g_off;
    uint32_t i;
    uint32_t extra;
    uint32_t cds_pos;
    uint32_t mode;
    int force_terminal;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];
    force_terminal = kprop_bounded(t, 8u) == 0u;
    if (force_terminal) {
        uint32_t stop = (uint32_t)kprop_bounded(t, 3u);
        memcpy(s->cds + cds_len - 3u, STOPS[stop], 3u);
    }

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? 1 : -1; s->flags = 0u;
    s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u;
    s->cds_lenv = cds_len;

    mode = (uint32_t)kprop_bounded(t, 3u);
    if (mode == 0u) {
        extra = (uint32_t)kprop_bounded(t, 2u) + 1u; /* net +1 or +2, both frameshift */
        if (force_terminal) {
            cds_pos = s->strand > 0 ? cds_len - 2u : cds_len - 1u;
        } else {
            cds_pos = (uint32_t)kprop_bounded(t, cds_len - 6u) + 4u;
        }
        s->vpos = kprop_genomic_pos_for_cds(s, cds_pos);
        s->vend = s->vpos; s->vkind = (uint8_t)DUCKVEP_KIND_INS;
        s->abytes[0] = (uint8_t)kprop_genomic_base_at(s, s->vpos);
        s->abytes[1] = s->abytes[0];
        for (i = 0u; i < extra; i++) s->abytes[2u + i] = (uint8_t)BASES[kprop_bounded(t, 4u)];
        s->roff = 0u; s->aoff = 1u; s->rlen = 1u; s->alen = (uint16_t)(1u + extra);
    } else if (mode == 1u) {
        extra = (uint32_t)kprop_bounded(t, 2u) + 1u; /* delete 1 or 2 bases after anchor */
        cds_pos = (uint32_t)kprop_bounded(t, cds_len - extra - 5u) + 4u;
        if (s->strand > 0) {
            s->vpos = kprop_genomic_pos_for_cds(s, cds_pos - 1u);
        } else {
            s->vpos = kprop_genomic_pos_for_cds(s, cds_pos + extra);
        }
        g_off = s->vpos - base;
        s->vend = s->vpos + extra; s->vkind = (uint8_t)DUCKVEP_KIND_DEL;
        for (i = 0u; i <= extra; i++) s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, base + g_off + i);
        s->abytes[extra + 1u] = s->abytes[0];
        s->roff = 0u; s->aoff = extra + 1u; s->rlen = (uint16_t)(extra + 1u); s->alen = 1u;
    } else {
        uint32_t delins_mode = (uint32_t)kprop_bounded(t, 4u);
        uint32_t ref_len = delins_mode < 2u ? 1u : delins_mode;
        uint32_t alt_len = delins_mode < 2u ? delins_mode + 2u : 1u;
        uint32_t first_cds = (uint32_t)kprop_bounded(t, cds_len - ref_len - 5u) + 4u;
        s->vpos = s->strand > 0
            ? kprop_genomic_pos_for_cds(s, first_cds)
            : kprop_genomic_pos_for_cds(s, first_cds + ref_len - 1u);
        s->vend = s->vpos + ref_len - 1u;
        s->vkind = (uint8_t)DUCKVEP_KIND_INDEL;
        for (i = 0u; i < ref_len; i++) {
            s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, s->vpos + i);
        }
        for (i = 0u; i < alt_len; i++) {
            s->abytes[ref_len + i] = kprop_base_not(s->abytes[0], s->abytes[ref_len - 1u]);
        }
        s->roff = 0u; s->aoff = ref_len;
        s->rlen = (uint16_t)ref_len; s->alen = (uint16_t)alt_len;
    }

    s->vchrom = 0u;
    s->cds_off0 = 0u; s->cds_lenv = cds_len; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = sizeof s->abytes;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_frameshift_indel_info = {
    .alloc = kprop_frameshift_indel_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_inframe_deletion_alloc(struct theft *t, void *env, void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 3u) + 4u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t codon_idx;
    uint32_t del_start_cds;
    uint32_t del_end_cds;
    uint32_t anchor_cds;
    uint32_t i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? 1 : -1; s->flags = 0u;
    s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u;
    s->cds_lenv = cds_len;

    codon_idx = (uint32_t)kprop_bounded(t, ncodons - 2u) + 1u; /* exclude first/last codon */
    del_start_cds = codon_idx * 3u + 1u;
    del_end_cds = del_start_cds + 2u;
    /* Keep the deleted codon non-stop: real coding transcripts have no internal stop, and
     * removing an internal stop is stop_lost (not a clean in-frame deletion) — the classifier
     * defers those, so exclude them here to keep this generator to clean in-frame deletions.
     * TAA/TAG/TGA -> T[C]A/T[C]G/T[C]A (Ser), a minimal non-stop rewrite of the middle base. */
    {
        uint8_t *dc = &s->cds[(size_t)del_start_cds - 1u];
        if (dc[0] == (uint8_t)'T' &&
            ((dc[1] == (uint8_t)'A' && (dc[2] == (uint8_t)'A' || dc[2] == (uint8_t)'G')) ||
             (dc[1] == (uint8_t)'G' && dc[2] == (uint8_t)'A'))) {
            dc[1] = (uint8_t)'C';
        }
    }
    anchor_cds = s->strand > 0 ? del_start_cds - 1u : del_end_cds + 1u;
    s->vpos = kprop_genomic_pos_for_cds(s, anchor_cds);
    s->vend = s->vpos + 3u;
    s->vkind = (uint8_t)DUCKVEP_KIND_DEL;
    for (i = 0u; i < 4u; i++) s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, s->vpos + i);
    s->abytes[4] = s->abytes[0];
    s->roff = 0u; s->aoff = 4u; s->rlen = 4u; s->alen = 1u;

    s->vchrom = 0u;
    s->cds_off0 = 0u; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = sizeof s->abytes;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_inframe_deletion_info = {
    .alloc = kprop_inframe_deletion_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_inframe_insertion_alloc(struct theft *t, void *env, void **instance) {
    static const char BASES[4] = {'A', 'C', 'G', 'T'};
    static const char INSERTED[3] = {'G', 'C', 'C'}; /* alanine, not a stop codon */
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 4u) + 5u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t codon_idx;
    uint32_t before_cds;
    uint32_t anchor_cds;
    uint32_t i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)BASES[kprop_bounded(t, 4u)];

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? 1 : -1; s->flags = 0u;
    s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u;
    s->cds_lenv = cds_len;

    codon_idx = (uint32_t)kprop_bounded(t, ncodons - 4u) + 2u; /* after codon 2..n-2 */
    before_cds = codon_idx * 3u;
    anchor_cds = s->strand > 0 ? before_cds : before_cds + 1u;
    s->vpos = kprop_genomic_pos_for_cds(s, anchor_cds);
    s->vend = s->vpos;
    s->vkind = (uint8_t)DUCKVEP_KIND_INS;
    s->abytes[0] = (uint8_t)kprop_genomic_base_at(s, s->vpos);
    s->abytes[1] = s->abytes[0];
    for (i = 0u; i < 3u; i++) {
        char b = s->strand > 0 ? INSERTED[i] : kprop_complement_base(INSERTED[2u - i]);
        s->abytes[2u + i] = (uint8_t)b;
    }
    s->roff = 0u; s->aoff = 1u; s->rlen = 1u; s->alen = 4u;

    s->vchrom = 0u;
    s->cds_off0 = 0u; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = sizeof s->abytes;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_inframe_insertion_info = {
    .alloc = kprop_inframe_insertion_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_delins_shape_alloc(struct theft *t, void *env, void **instance) {
    static const uint8_t ALT3_PRESERVED[3] = { 'A', 'A', 'A' };
    static const uint8_t ALT3_ALTERED[3] = { 'G', 'C', 'C' };
    static const uint8_t ALT6_PRESERVED[6] = { 'A', 'A', 'G', 'G', 'C', 'C' };
    static const uint8_t ALT6_ALTERED[6] = { 'G', 'C', 'C', 'G', 'C', 'T' };
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 4u) + 5u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t lengthening = (uint32_t)kprop_bounded(t, 2u) == 0u;
    uint32_t preserved = (uint32_t)kprop_bounded(t, 2u) == 0u;
    uint32_t codon_idx;
    uint32_t cds_start;
    uint32_t ref_len;
    uint32_t alt_len;
    const uint8_t *alt_tx;
    uint32_t i;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)'A';

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? 1 : -1; s->flags = 0u;
    s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u; s->vchrom = 0u;
    kprop_wire_coding_scene(s, cds_len);

    ref_len = lengthening ? 3u : 6u;
    alt_len = lengthening ? 6u : 3u;
    alt_tx = lengthening
        ? (preserved ? ALT6_PRESERVED : ALT6_ALTERED)
        : (preserved ? ALT3_PRESERVED : ALT3_ALTERED);
    codon_idx = (uint32_t)kprop_bounded(t, ncodons - 3u) + 1u;
    cds_start = codon_idx * 3u + 1u;
    s->expect_shape = lengthening ? (uint8_t)KPROP_CDS_EDIT_INS
                                  : (uint8_t)KPROP_CDS_EDIT_DEL;
    s->expect_region = (uint8_t)KPROP_CDS_EDIT_BODY;
    s->expect_protein_altering = preserved ? 0u : 1u;
    s->vkind = (uint8_t)DUCKVEP_KIND_INDEL;
    s->vpos = s->strand > 0 ? kprop_genomic_pos_for_cds(s, cds_start)
                            : kprop_genomic_pos_for_cds(s, cds_start + ref_len - 1u);
    s->vend = s->vpos + ref_len - 1u;
    for (i = 0u; i < ref_len; i++) {
        s->abytes[i] = (uint8_t)kprop_genomic_base_at(s, s->vpos + i);
    }
    kprop_fill_variant_alt_from_tx(s, ref_len, alt_tx, alt_len);
    s->roff = 0u; s->aoff = ref_len;
    s->rlen = (uint16_t)ref_len; s->alen = (uint16_t)alt_len;
    kprop_fill_expected_cds(s, cds_start, ref_len, alt_tx, alt_len);

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_delins_shape_info = {
    .alloc = kprop_delins_shape_alloc,
    .free  = kprop_coding_free,
};

static enum theft_alloc_res kprop_protein_altering_insertion_alloc(struct theft *t, void *env, void **instance) {
    static const char INSERTED[3] = {'G', 'C', 'C'}; /* AAA + inserted GCC stays non-stop */
    struct kprop_coding *s = (struct kprop_coding *)calloc(1u, sizeof *s);
    uint32_t ncodons = (uint32_t)kprop_bounded(t, KPROP_MAX_CODONS - 4u) + 5u;
    uint32_t cds_len = 3u * ncodons;
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFF00000u) + 1000u;
    uint32_t codon_idx;
    uint32_t before_cds;
    uint32_t anchor_cds;
    uint32_t i;
    uint32_t off;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->cds = (uint8_t *)malloc(cds_len);
    if (s->cds == NULL) { free(s); return THEFT_ALLOC_ERROR; }
    for (i = 0u; i < cds_len; i++) s->cds[i] = (uint8_t)'A';

    s->chrom = 0u; s->strand = kprop_bounded(t, 2u) == 0u ? 1 : -1; s->flags = 0u;
    s->tstart = base; s->tend = base + cds_len - 1u;
    s->cds_s = base; s->cds_e = base + cds_len - 1u;
    s->es = base; s->ee = base + cds_len - 1u;
    s->ecds = 1u; s->ecde = cds_len; s->eph = 0; s->eeph = 0;
    s->exoff = 0u; s->excnt = 1u;
    s->cds_lenv = cds_len;

    codon_idx = (uint32_t)kprop_bounded(t, ncodons - 2u) + 1u; /* codon 2..n-1 */
    off = (uint32_t)kprop_bounded(t, 2u) + 1u;                 /* after base 1 or 2 */
    before_cds = codon_idx * 3u + off;
    anchor_cds = s->strand > 0 ? before_cds : before_cds + 1u;
    s->vpos = kprop_genomic_pos_for_cds(s, anchor_cds);
    s->vend = s->vpos;
    s->vkind = (uint8_t)DUCKVEP_KIND_INS;
    s->abytes[0] = (uint8_t)kprop_genomic_base_at(s, s->vpos);
    s->abytes[1] = s->abytes[0];
    for (i = 0u; i < 3u; i++) {
        char b = s->strand > 0 ? INSERTED[i] : kprop_complement_base(INSERTED[2u - i]);
        s->abytes[2u + i] = (uint8_t)b;
    }
    s->roff = 0u; s->aoff = 1u; s->rlen = 1u; s->alen = 4u;

    s->vchrom = 0u;
    s->cds_off0 = 0u; s->ctab = (uint8_t)DUCKVEP_CODON_TABLE_STANDARD;
    s->ex.start1 = &s->es; s->ex.end1 = &s->ee;
    s->ex.cdna_start1 = &s->ecds; s->ex.cdna_end1 = &s->ecde;
    s->ex.phase = &s->eph; s->ex.end_phase = &s->eeph; s->ex.exon_count = 1u;
    s->tx.chrom_id = &s->chrom; s->tx.start1 = &s->tstart; s->tx.end1 = &s->tend;
    s->tx.strand = &s->strand; s->tx.flags = &s->flags;
    s->tx.exon_offset = &s->exoff; s->tx.exon_count = &s->excnt;
    s->tx.cds_start1 = &s->cds_s; s->tx.cds_end1 = &s->cds_e; s->tx.transcript_count = 1u;
    s->seq.cds_bytes = s->cds; s->seq.cds_bytes_len = cds_len;
    s->seq.cds_offset = &s->cds_off0; s->seq.cds_length = &s->cds_lenv;
    s->seq.codon_table = &s->ctab; s->seq.transcript_count = 1u;
    s->v.chrom_id = &s->vchrom; s->v.pos1 = &s->vpos; s->v.end1 = &s->vend;
    s->v.ref_offset = &s->roff; s->v.ref_length = &s->rlen;
    s->v.alt_offset = &s->aoff; s->v.alt_length = &s->alen;
    s->v.allele_bytes = s->abytes; s->v.allele_bytes_len = sizeof s->abytes;
    s->v.variant_kind = &s->vkind; s->v.count = 1u;

    *instance = s;
    return THEFT_ALLOC_OK;
}

static struct theft_type_info kprop_protein_altering_insertion_info = {
    .alloc = kprop_protein_altering_insertion_alloc,
    .free  = kprop_coding_free,
};

static uint64_t codon_change_to_so(uint32_t change, char aa_ref) {
    if (change & DUCKVEP_CODON_STOP_GAINED)     return DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED);
    if (change & DUCKVEP_CODON_STOP_LOST)       return DUCKVEP_SO(DUCKVEP_SO_STOP_LOST);
    if (change & DUCKVEP_CODON_MISSENSE)        return DUCKVEP_SO(DUCKVEP_SO_MISSENSE);
    if (change & DUCKVEP_CODON_SYNONYMOUS) {
        return (aa_ref == '*') ? DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED)
                               : DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS);
    }
    return 0u;
}

static uint64_t start_codon_snv_to_so(uint32_t change, char aa_ref, char aa_alt) {
    uint64_t mask;
    if (aa_alt == 'M') return codon_change_to_so(change, aa_ref);
    mask = DUCKVEP_SO(DUCKVEP_SO_START_LOST);
    if (change & DUCKVEP_CODON_STOP_GAINED) mask |= DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED);
    else if (change & DUCKVEP_CODON_STOP_LOST) mask |= DUCKVEP_SO(DUCKVEP_SO_STOP_LOST);
    else if (change & DUCKVEP_CODON_SYNONYMOUS) {
        mask |= (aa_ref == '*') ? DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED)
                                : DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS);
    }
    return mask;
}

/* Independent VEP-116 oracle for one codon selected by an equal-length uploaded
 * feature. Unlike a minimized SNV, the start predicates inspect the complete
 * feature: unequal local peptides make start_lost true, while the rebuilt ATG
 * independently makes start_retained_variant true. */
static uint64_t vep_feature_codon_to_so(
    char ref_aa, char alt_aa, const char alt_codon[4], int overlaps_start) {

    uint64_t mask = 0u;
    int start_lost = overlaps_start && ref_aa != alt_aa;
    int start_retained = overlaps_start &&
        alt_codon[0] == 'A' && alt_codon[1] == 'T' && alt_codon[2] == 'G';

    if (start_lost) mask |= DUCKVEP_SO(DUCKVEP_SO_START_LOST);
    if (start_retained) mask |= DUCKVEP_SO(DUCKVEP_SO_START_RETAINED);
    if (ref_aa == '*' && alt_aa != '*') {
        mask |= DUCKVEP_SO(DUCKVEP_SO_STOP_LOST);
    } else if (ref_aa == '*' && alt_aa == '*') {
        mask |= DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED);
    } else if (ref_aa != '*' && alt_aa == '*') {
        mask |= DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED);
    } else if (ref_aa == alt_aa) {
        if (!start_retained) mask |= DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS);
    } else if (!start_lost) {
        mask |= DUCKVEP_SO(DUCKVEP_SO_MISSENSE);
    }
    return mask;
}

/* Coverage accounting so "across codon-change classes" is a checked claim, not a
 * hope: the run must actually observe each class at least once (the framework's
 * coverage signal, applied locally). */
static struct { uint32_t syn; uint32_t mis; uint32_t sg; uint32_t sl; uint32_t sr; } g_codon_cov;

static enum theft_trial_res prop_annotate_codon_matches_kernel(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    duckvep_coding_projection_t proj;
    duckvep_coding_snv_result_t res;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;
    memset(&err, 0, sizeof err);

    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }

    duckvep_result_builder_init(&rb, rows, 2u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_result_builder_count(&rb) != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }

    /* Oracle: the tested kernels invoked directly. */
    if (!duckvep_project_coding_base(&s->tx, &s->ex, 0u, s->vpos, &proj)) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_coding_snv_from_cds(s->cds, (size_t)s->cds_lenv, &proj,
                                    (char)s->abytes[0], (char)s->abytes[1], 1,
                                    DUCKVEP_CODON_TABLE_STANDARD, &res) != DUCKVEP_CODING_SNV_OK) {
        tr = THEFT_TRIAL_FAIL; goto done;
    }

    if (res.change & DUCKVEP_CODON_STOP_GAINED)     g_codon_cov.sg++;
    else if (res.change & DUCKVEP_CODON_STOP_LOST)  g_codon_cov.sl++;
    else if (res.change & DUCKVEP_CODON_MISSENSE)   g_codon_cov.mis++;
    else if (res.change & DUCKVEP_CODON_SYNONYMOUS) {
        if (res.aa_ref == '*') g_codon_cov.sr++; else g_codon_cov.syn++;
    }

    if (rows[0].consequence_mask != codon_change_to_so(res.change, res.aa_ref)) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].cds_pos != (int32_t)res.cds_pos) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].protein_pos != (int32_t)res.protein_pos) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].aa_ref != (uint8_t)res.aa_ref) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].aa_alt != (uint8_t)res.aa_alt) { tr = THEFT_TRIAL_FAIL; goto done; }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return tr;
}

static struct { uint32_t sl; uint32_t sl_sg; uint32_t sl_syn; uint32_t alt_m; } g_start_cov;

static enum theft_trial_res prop_annotate_start_lost_matches_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    duckvep_coding_projection_t proj;
    duckvep_coding_snv_result_t res;
    uint64_t want;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;
    memset(&err, 0, sizeof err);

    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    duckvep_result_builder_init(&rb, rows, 2u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_result_builder_count(&rb) != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }

    if (!duckvep_project_coding_base(&s->tx, &s->ex, 0u, s->vpos, &proj)) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (proj.codon_start_cds != 1u || proj.protein_pos != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_coding_snv_from_cds(s->cds, (size_t)s->cds_lenv, &proj,
                                    (char)s->abytes[0], (char)s->abytes[1], 1,
                                    DUCKVEP_CODON_TABLE_STANDARD, &res) != DUCKVEP_CODING_SNV_OK) {
        tr = THEFT_TRIAL_FAIL; goto done;
    }
    want = start_codon_snv_to_so(res.change, res.aa_ref, res.aa_alt);
    if (res.aa_alt == 'M') g_start_cov.alt_m++;
    if (want & DUCKVEP_SO(DUCKVEP_SO_START_LOST)) {
        g_start_cov.sl++;
        if (want & DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED)) g_start_cov.sl_sg++;
        if (want & DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS))  g_start_cov.sl_syn++;
    }

    if (rows[0].consequence_mask != want) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].protein_pos != 1) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].aa_ref != (uint8_t)res.aa_ref) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].aa_alt != (uint8_t)res.aa_alt) { tr = THEFT_TRIAL_FAIL; goto done; }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return tr;
}

TEST annotate_codon_matches_kernel_for_any_cds_snv(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile codon refinement == coding-SNV kernel oracle";
    cfg.prop1 = prop_annotate_codon_matches_kernel;
    cfg.type_info[0] = &kprop_coding_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_codon_cov, 0, sizeof g_codon_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    /* The "across codon-change classes" claim is only honest if the run actually
     * hit each class. (Deterministic under the fixed seed.) */
    ASSERT(g_codon_cov.syn > 0u);
    ASSERT(g_codon_cov.mis > 0u);
    ASSERT(g_codon_cov.sg > 0u);
    ASSERT(g_codon_cov.sl > 0u);
    ASSERT(g_codon_cov.sr > 0u); /* stop->stop synonymous reaches stop_retained */
    fprintf(stderr, "[codon coverage] syn=%u mis=%u stop_gained=%u stop_lost=%u stop_retained=%u\n",
            g_codon_cov.syn, g_codon_cov.mis, g_codon_cov.sg, g_codon_cov.sl, g_codon_cov.sr);
    PASS();
}

static struct { uint32_t len2; uint32_t len3; } g_mnv_cov;

static enum theft_trial_res prop_annotate_mnv_matches_codon_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    char ref_codon[4];
    char alt_codon[4];
    duckvep_codon_result_t cr;
    uint32_t cds0;
    uint32_t codon_start;
    uint32_t protein_pos;
    uint8_t off0;
    uint32_t i;
    uint64_t want;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;
    memset(&err, 0, sizeof err);

    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }

    duckvep_result_builder_init(&rb, rows, 2u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_result_builder_count(&rb) != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }

    /* Independent single-exon '+' oracle: compute the CDS/codon coordinates directly,
     * edit the codon slice, then classify via the canonical genetic-code helper. */
    cds0 = s->vpos - s->tstart + 1u;
    codon_start = ((cds0 - 1u) / 3u) * 3u + 1u;
    protein_pos = ((codon_start - 1u) / 3u) + 1u;
    off0 = (uint8_t)(cds0 - codon_start);
    ref_codon[0] = (char)s->cds[codon_start - 1u];
    ref_codon[1] = (char)s->cds[codon_start];
    ref_codon[2] = (char)s->cds[codon_start + 1u];
    ref_codon[3] = '\0';
    alt_codon[0] = ref_codon[0]; alt_codon[1] = ref_codon[1]; alt_codon[2] = ref_codon[2]; alt_codon[3] = '\0';
    for (i = 0u; i < (uint32_t)s->alen; i++) {
        alt_codon[off0 + i] = (char)s->abytes[s->aoff + i];
    }
    cr = duckvep_codon_change(ref_codon, alt_codon, DUCKVEP_CODON_TABLE_STANDARD);
    want = (codon_start == 1u && protein_pos == 1u)
        ? vep_feature_codon_to_so(cr.aa_ref, cr.aa_alt, alt_codon, 1)
        : codon_change_to_so(cr.change, cr.aa_ref);

    if (s->alen == 2u) g_mnv_cov.len2++; else if (s->alen == 3u) g_mnv_cov.len3++;
    if (rows[0].consequence_mask != want) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].cdna_pos != -1) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].cds_pos != -1) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].protein_pos != (int32_t)protein_pos) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].aa_ref != (uint8_t)cr.aa_ref) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].aa_alt != (uint8_t)cr.aa_alt) { tr = THEFT_TRIAL_FAIL; goto done; }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return tr;
}

TEST annotate_mnv_same_codon_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile same-codon MNV == codon oracle";
    cfg.prop1 = prop_annotate_mnv_matches_codon_oracle;
    cfg.type_info[0] = &kprop_mnv_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_mnv_cov, 0, sizeof g_mnv_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_mnv_cov.len2 > 0u);
    ASSERT(g_mnv_cov.len3 > 0u);
    fprintf(stderr, "[mnv coverage] len2=%u len3=%u\n", g_mnv_cov.len2, g_mnv_cov.len3);
    PASS();
}

static int kprop_cross_codon_mnv_oracle_missense(const struct kprop_coding *s) {
    uint32_t codon_start[2] = {0u, 0u};
    uint8_t seen[2] = {0u, 0u};
    char ref_codon[2][4];
    char alt_codon[2][4];
    uint32_t codon_count = 0u;
    uint32_t j;
    int changed = 0;

    for (j = 0u; j < (uint32_t)s->rlen; j++) {
        uint32_t gpos = s->vpos + j;
        uint32_t cds_pos = kprop_cds_pos_for_genomic(s, gpos);
        uint32_t start = ((cds_pos - 1u) / 3u) * 3u + 1u;
        uint32_t off = cds_pos - start;
        uint32_t k;
        uint32_t slot = UINT32_MAX;
        uint8_t bit;
        char ref_tx;
        char alt_tx;

        if (start <= 1u || start + 2u > s->cds_lenv - 3u) return 0;
        for (k = 0u; k < codon_count; k++) {
            if (codon_start[k] == start) { slot = k; break; }
        }
        if (slot == UINT32_MAX) {
            if (codon_count >= 2u) return 0;
            slot = codon_count++;
            codon_start[slot] = start;
            ref_codon[slot][0] = (char)s->cds[start - 1u];
            ref_codon[slot][1] = (char)s->cds[start];
            ref_codon[slot][2] = (char)s->cds[start + 1u];
            ref_codon[slot][3] = '\0';
            memcpy(alt_codon[slot], ref_codon[slot], sizeof alt_codon[slot]);
        }
        bit = (uint8_t)(1u << off);
        if (seen[slot] & bit) return 0;
        seen[slot] |= bit;
        ref_tx = s->strand > 0 ? (char)s->abytes[s->roff + j]
                               : kprop_complement_base((char)s->abytes[s->roff + j]);
        alt_tx = s->strand > 0 ? (char)s->abytes[s->aoff + j]
                               : kprop_complement_base((char)s->abytes[s->aoff + j]);
        if (ref_codon[slot][off] != ref_tx) return 0;
        alt_codon[slot][off] = alt_tx;
    }
    if (codon_count != 2u) return 0;
    if (codon_start[0] > codon_start[1]) {
        uint32_t tmp = codon_start[0];
        codon_start[0] = codon_start[1];
        codon_start[1] = tmp;
    }
    if (codon_start[1] - codon_start[0] != 3u) return 0;

    for (j = 0u; j < 2u; j++) {
        duckvep_codon_result_t cr = duckvep_codon_change(
            ref_codon[j], alt_codon[j], DUCKVEP_CODON_TABLE_STANDARD);
        if (cr.change & DUCKVEP_CODON_INVALID) return 0;
        if (cr.aa_ref == '*' || cr.aa_alt == '*') return 0;
        if (cr.aa_ref != cr.aa_alt) changed = 1;
    }
    return changed;
}

/* Full independent two-codon-window oracle for a body MNV: re-collects the (<=2) affected
 * body codons exactly like kprop_cross_codon_mnv_oracle_missense, then classifies the whole
 * window with VEP's precedence (stop_lost > stop_retained > stop_gained > synonymous/
 * missense) from the genetic code. `.valid == 0` means genuinely unsupported (N base, invalid
 * codon, or not two adjacent body codons) — where the kernel falls back to
 * coding_sequence_variant. Since the generator never touches codon 0, start_lost never
 * applies here. This is the oracle the generalized window classifier is graded against. */
struct kprop_cc_facts {
    int valid;
    int synonymous, missense, stop_gained, stop_lost, stop_retained;
};

static struct kprop_cc_facts kprop_cross_codon_mnv_oracle_facts(const struct kprop_coding *s) {
    struct kprop_cc_facts f;
    uint32_t codon_start[2] = {0u, 0u};
    uint8_t seen[2] = {0u, 0u};
    char ref_codon[2][4];
    char alt_codon[2][4];
    uint32_t codon_count = 0u;
    uint32_t j;
    int has_ref_stop = 0;
    int has_alt_stop = 0;
    int ref_stop_idx = -1;
    int alt_stop_idx = -1;
    int win_equal = 1;

    memset(&f, 0, sizeof f);
    for (j = 0u; j < (uint32_t)s->rlen; j++) {
        uint32_t gpos = s->vpos + j;
        uint32_t cds_pos = kprop_cds_pos_for_genomic(s, gpos);
        uint32_t start = ((cds_pos - 1u) / 3u) * 3u + 1u;
        uint32_t off = cds_pos - start;
        uint32_t k;
        uint32_t slot = UINT32_MAX;
        uint8_t bit;
        char ref_tx;
        char alt_tx;

        if (start <= 1u || start + 2u > s->cds_lenv - 3u) return f;
        for (k = 0u; k < codon_count; k++) {
            if (codon_start[k] == start) { slot = k; break; }
        }
        if (slot == UINT32_MAX) {
            if (codon_count >= 2u) return f;
            slot = codon_count++;
            codon_start[slot] = start;
            ref_codon[slot][0] = (char)s->cds[start - 1u];
            ref_codon[slot][1] = (char)s->cds[start];
            ref_codon[slot][2] = (char)s->cds[start + 1u];
            ref_codon[slot][3] = '\0';
            memcpy(alt_codon[slot], ref_codon[slot], sizeof alt_codon[slot]);
        }
        bit = (uint8_t)(1u << off);
        if (seen[slot] & bit) return f;
        seen[slot] |= bit;
        ref_tx = s->strand > 0 ? (char)s->abytes[s->roff + j]
                               : kprop_complement_base((char)s->abytes[s->roff + j]);
        alt_tx = s->strand > 0 ? (char)s->abytes[s->aoff + j]
                               : kprop_complement_base((char)s->abytes[s->aoff + j]);
        if (ref_codon[slot][off] != ref_tx) return f;
        alt_codon[slot][off] = alt_tx;
    }
    if (codon_count != 2u) return f;
    if (codon_start[0] > codon_start[1]) {
        uint32_t tmp = codon_start[0];
        codon_start[0] = codon_start[1];
        codon_start[1] = tmp;
    }
    if (codon_start[1] - codon_start[0] != 3u) return f;

    for (j = 0u; j < 2u; j++) {
        duckvep_codon_result_t cr = duckvep_codon_change(
            ref_codon[j], alt_codon[j], DUCKVEP_CODON_TABLE_STANDARD);
        if (cr.change & DUCKVEP_CODON_INVALID) return f;
        if (cr.aa_ref == 'X' || cr.aa_alt == 'X') return f;
        if (cr.aa_ref != cr.aa_alt) win_equal = 0;
        if (cr.aa_ref == '*' && !has_ref_stop) { ref_stop_idx = (int)j; has_ref_stop = 1; }
        if (cr.aa_alt == '*' && !has_alt_stop) { alt_stop_idx = (int)j; has_alt_stop = 1; }
    }
    f.valid = 1;
    if (has_ref_stop && !has_alt_stop) f.stop_lost = 1;
    else if (has_ref_stop && has_alt_stop && ref_stop_idx == alt_stop_idx) f.stop_retained = 1;
    else if (has_alt_stop && !has_ref_stop) f.stop_gained = 1;
    else if (win_equal) f.synonymous = 1;
    else f.missense = 1;
    return f;
}

/* SO consequence mask that annotate_tile emits for a resolved two-codon window (coarse:
 * protein_pos -1, no AA pair), or coding_sequence_variant when the oracle is unsupported. */
static uint64_t kprop_cross_codon_facts_to_so(const struct kprop_cc_facts *f) {
    if (!f->valid) return DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE);
    if (f->stop_gained) return DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED);
    if (f->stop_lost) return DUCKVEP_SO(DUCKVEP_SO_STOP_LOST);
    if (f->stop_retained) return DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED);
    if (f->synonymous) return DUCKVEP_SO(DUCKVEP_SO_SYNONYMOUS);
    return DUCKVEP_SO(DUCKVEP_SO_MISSENSE);
}

static struct { uint32_t missense; uint32_t synonymous; uint32_t stop_gained; uint32_t fwd; uint32_t rev; uint32_t len2; uint32_t len3; } g_cross_mnv_cov;

static enum theft_trial_res prop_annotate_cross_codon_mnv_matches_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    struct kprop_cc_facts f;
    uint64_t want;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;
    memset(&err, 0, sizeof err);

    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }

    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, rows, 2u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_result_builder_count(&rb) != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }

    f = kprop_cross_codon_mnv_oracle_facts(s);
    want = kprop_cross_codon_facts_to_so(&f);
    if (f.missense) g_cross_mnv_cov.missense++;
    else if (f.synonymous) g_cross_mnv_cov.synonymous++;
    else if (f.stop_gained) g_cross_mnv_cov.stop_gained++;
    if (s->strand > 0) g_cross_mnv_cov.fwd++; else g_cross_mnv_cov.rev++;
    if (s->rlen == 2u) g_cross_mnv_cov.len2++; else if (s->rlen == 3u) g_cross_mnv_cov.len3++;

    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].consequence_mask != want) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (f.valid) {
        if (stats->substitution_context != 1u ||
            stats->mnv_direct_fallback != 0u) {
            tr = THEFT_TRIAL_FAIL; goto done;
        }
    } else {
        if (stats->substitution_context != 0u ||
            stats->mnv_direct_fallback != 1u) {
            tr = THEFT_TRIAL_FAIL; goto done;
        }
    }
    if (rows[0].cdna_pos != -1 || rows[0].cds_pos != -1 || rows[0].protein_pos != -1) {
        tr = THEFT_TRIAL_FAIL;
        goto done;
    }
    if (rows[0].aa_ref != (uint8_t)0u || rows[0].aa_alt != (uint8_t)0u) {
        tr = THEFT_TRIAL_FAIL;
        goto done;
    }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return tr;
}

TEST annotate_cross_codon_mnv_missense_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile two-codon body MNV missense == codon-window oracle";
    cfg.prop1 = prop_annotate_cross_codon_mnv_matches_oracle;
    cfg.type_info[0] = &kprop_cross_codon_mnv_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_cross_mnv_cov, 0, sizeof g_cross_mnv_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_cross_mnv_cov.missense > 0u);
    ASSERT(g_cross_mnv_cov.synonymous > 0u);
    ASSERT(g_cross_mnv_cov.stop_gained > 0u);
    ASSERT(g_cross_mnv_cov.fwd > 0u);
    ASSERT(g_cross_mnv_cov.rev > 0u);
    ASSERT(g_cross_mnv_cov.len2 > 0u);
    ASSERT(g_cross_mnv_cov.len3 > 0u);
    fprintf(stderr, "[cross-mnv coverage] missense=%u synonymous=%u stop_gained=%u fwd=%u rev=%u len2=%u len3=%u\n",
            g_cross_mnv_cov.missense, g_cross_mnv_cov.synonymous, g_cross_mnv_cov.stop_gained,
            g_cross_mnv_cov.fwd, g_cross_mnv_cov.rev,
            g_cross_mnv_cov.len2, g_cross_mnv_cov.len3);
    PASS();
}

static int kprop_cds_edits_equal(const duckvep_haplotype_edit_t *a,
                                  const duckvep_haplotype_edit_t *b) {
    return a != NULL && b != NULL &&
           a->cds_start == b->cds_start &&
           a->ref_len == b->ref_len &&
           a->ref == b->ref &&
           a->alt_len == b->alt_len &&
           a->alt == b->alt &&
           a->variant_strand == b->variant_strand;
}

static int kprop_context_zero(const duckvep_coding_context_t *ctx) {
    duckvep_coding_context_t z;
    memset(&z, 0, sizeof z);
    return ctx != NULL && memcmp(ctx, &z, sizeof z) == 0;
}

static char kprop_context_norm_base(uint8_t c) {
    switch ((char)c) {
        case 'A': case 'a': return 'A';
        case 'C': case 'c': return 'C';
        case 'G': case 'g': return 'G';
        case 'T': case 't': case 'U': case 'u': return 'T';
        case 'N': case 'n': return 'N';
        default: return '\0';
    }
}

static int kprop_translate_full_oracle(const uint8_t *cds, size_t cds_len,
                                       duckvep_codon_table_t table,
                                       uint8_t *pep, size_t *pep_len) {
    size_t codons = cds_len / 3u;
    size_t i;
    if (cds == NULL || pep == NULL || pep_len == NULL) return 0;
    for (i = 0u; i < codons; i++) {
        char codon[4];
        uint32_t j;
        int has_n = 0;
        for (j = 0u; j < 3u; j++) {
            char b = kprop_context_norm_base(cds[i * 3u + (size_t)j]);
            if (b == '\0') return 0;
            if (b == 'N') has_n = 1;
            codon[j] = b;
        }
        codon[3] = '\0';
        pep[i] = has_n ? (uint8_t)'X'
                       : (uint8_t)duckvep_translate_codon(codon, table);
    }
    pep[codons] = (uint8_t)'\0';
    *pep_len = codons;
    return 1;
}

static int kprop_cds_changed_oracle(const uint8_t *ref_cds, size_t ref_len,
                                    const uint8_t *alt_cds, size_t alt_len) {
    size_t i;
    if (ref_len != alt_len) return 1;
    for (i = 0u; i < ref_len; i++) {
        char rb = kprop_context_norm_base(ref_cds[i]);
        char ab = kprop_context_norm_base(alt_cds[i]);
        if (rb == '\0' || ab == '\0') return -1;
        if (rb != ab) return 1;
    }
    return 0;
}

static void kprop_peptide_window_oracle(const uint8_t *ref_pep, size_t ref_len,
                                        const uint8_t *alt_pep, size_t alt_len,
                                        uint32_t *rf, uint32_t *rl,
                                        uint32_t *af, uint32_t *al) {
    size_t prefix = 0u;
    size_t suffix = 0u;
    size_t ref_mid;
    size_t alt_mid;
    *rf = 0u; *rl = 0u; *af = 0u; *al = 0u;
    while (prefix < ref_len && prefix < alt_len && ref_pep[prefix] == alt_pep[prefix]) {
        prefix++;
    }
    if (prefix == ref_len && prefix == alt_len) return;
    while (suffix < ref_len - prefix && suffix < alt_len - prefix &&
           ref_pep[ref_len - 1u - suffix] == alt_pep[alt_len - 1u - suffix]) {
        suffix++;
    }
    ref_mid = ref_len - prefix - suffix;
    alt_mid = alt_len - prefix - suffix;
    if (ref_mid > 0u) { *rf = (uint32_t)prefix + 1u; *rl = (uint32_t)(ref_len - suffix); }
    if (alt_mid > 0u) { *af = (uint32_t)prefix + 1u; *al = (uint32_t)(alt_len - suffix); }
}

static uint32_t kprop_context_flags_oracle(const duckvep_edit_set_t *edit_set) {
    uint32_t flags = 0u;
    int saw_frameshift = 0;
    int64_t total_diff = 0;
    size_t i;
    for (i = 0u; i < edit_set->count; i++) {
        int64_t d = (int64_t)edit_set->edits[i].alt_len -
                    (int64_t)edit_set->edits[i].ref_len;
        total_diff += d;
        if (d != 0) flags |= DUCKVEP_HAPLOTYPE_FLAG_INDEL;
        if ((d % 3) != 0) saw_frameshift = 1;
    }
    if (saw_frameshift) {
        if ((total_diff % 3) == 0) flags |= DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT;
        else flags |= DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT;
    }
    return flags;
}

static uint32_t kprop_single_variant_flags_oracle(uint32_t ref_len, uint32_t alt_len) {
    int64_t d = (int64_t)alt_len - (int64_t)ref_len;
    uint32_t flags = 0u;
    if (d != 0) flags |= DUCKVEP_HAPLOTYPE_FLAG_INDEL;
    if ((d % 3) != 0) flags |= DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT;
    return flags;
}

TEST cds_edit_builder_projects_exon_boundary_insertions_on_both_strands(void) {
    static const uint16_t chrom[2] = {0u, 1u};
    static const uint32_t tx_start[2] = {100u, 100u};
    static const uint32_t tx_end[2] = {250u, 250u};
    static const int8_t strand[2] = {1, -1};
    static const uint32_t exon_offset[2] = {0u, 2u};
    static const uint16_t exon_count[2] = {2u, 2u};
    static const uint32_t cds_start[2] = {120u, 110u};
    static const uint32_t cds_end[2] = {240u, 230u};
    static const uint32_t exon_start[4] = {100u, 200u, 200u, 100u};
    static const uint32_t exon_end[4] = {150u, 250u, 250u, 150u};
    static const uint32_t cdna_start[4] = {1u, 52u, 1u, 52u};
    static const uint32_t cdna_end[4] = {51u, 102u, 51u, 102u};
    static const int8_t phase[4] = {0, 0, 0, 0};
    static const uint64_t cds_offset[2] = {0u, 72u};
    static const uint32_t cds_length[2] = {72u, 72u};
    static const uint8_t codon_table[2] = {
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD,
        (uint8_t)DUCKVEP_CODON_TABLE_STANDARD
    };
    static const uint32_t pos[2] = {199u, 199u};
    static const uint8_t kind[2] = {
        (uint8_t)DUCKVEP_KIND_INS, (uint8_t)DUCKVEP_KIND_INS
    };
    static const uint8_t alleles[6] = {'G','G','T', 'G','G','T'};
    static const uint32_t ref_offset[2] = {0u, 3u};
    static const uint32_t alt_offset[2] = {1u, 4u};
    static const uint16_t ref_length[2] = {1u, 1u};
    static const uint16_t alt_length[2] = {2u, 2u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_sequence_pool_t seq;
    duckvep_variant_batch_t variants;
    duckvep_haplotype_edit_t edit;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t results;
    duckvep_error_t error;
    uint8_t cds[144];
    size_t i;

    memset(&tx, 0, sizeof tx); memset(&exons, 0, sizeof exons);
    memset(&seq, 0, sizeof seq); memset(&variants, 0, sizeof variants);
    memset(&error, 0, sizeof error);
    memset(cds, 'A', sizeof cds);
    cds[0] = cds[72] = 'A'; cds[1] = cds[73] = 'T'; cds[2] = cds[74] = 'G';
    cds[69] = cds[141] = 'T'; cds[70] = cds[142] = 'A';
    cds[71] = cds[143] = 'A';
    tx.chrom_id = chrom; tx.start1 = tx_start; tx.end1 = tx_end; tx.strand = strand;
    tx.flags = k_zero_flags; tx.exon_offset = exon_offset; tx.exon_count = exon_count;
    tx.cds_start1 = cds_start; tx.cds_end1 = cds_end; tx.transcript_count = 2u;
    exons.start1 = exon_start; exons.end1 = exon_end;
    exons.cdna_start1 = cdna_start; exons.cdna_end1 = cdna_end;
    exons.phase = phase; exons.end_phase = phase; exons.exon_count = 4u;
    seq.cds_bytes = cds; seq.cds_bytes_len = sizeof cds;
    seq.cds_offset = cds_offset; seq.cds_length = cds_length;
    seq.codon_table = codon_table; seq.transcript_count = 2u;
    variants.chrom_id = chrom; variants.pos1 = pos; variants.end1 = pos;
    variants.variant_kind = kind; variants.allele_bytes = alleles;
    variants.allele_bytes_len = sizeof alleles; variants.ref_offset = ref_offset;
    variants.alt_offset = alt_offset; variants.ref_length = ref_length;
    variants.alt_length = alt_length; variants.count = 2u;

    for (i = 0u; i < 2u; i++) {
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
                  duckvep_variant_cds_edit_build(&tx, &exons, &seq, &variants,
                                                 (uint32_t)i, i, strand[i], &edit));
        ASSERT_EQ(32u, edit.cds_start);
        ASSERT_EQ(0u, edit.ref_len);
        ASSERT_EQ(1u, edit.alt_len);
    }

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, &model, &error));
    ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &options, &error));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &workspace, &error));
    duckvep_result_builder_init(&results, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &variants, options, workspace,
                                    &results, &error));
    ASSERT_EQ(2u, duckvep_result_builder_count(&results));
    for (i = 0u; i < 2u; i++) {
        ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT) |
                  DUCKVEP_SO(DUCKVEP_SO_SPLICE_REGION),
                  rows[i].consequence_mask);
        ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_RESOLVED, rows[i].sequence_status);
    }

    duckvep_workspace_close(workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
    PASS();
}

TEST cds_edit_builder_known_scene(void) {
    static uint8_t plus_cds[12] = {'A','T','G','A','A','A','C','C','C','T','A','A'};
    static uint8_t minus_cds[12] = {'A','T','G','A','A','A','C','C','C','T','A','A'};
    struct kprop_coding s;
    duckvep_haplotype_edit_t edit;
    duckvep_haplotype_edit_t scratch[3];
    duckvep_edit_set_t edit_set;
    duckvep_haplotype_result_t result;
    uint8_t mutated[32];
    size_t mutated_len;
    uint8_t payload[2] = {'A','G'};

    memset(&s, 0, sizeof s);
    s.cds = plus_cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1011u; s.cds_s = 1000u; s.cds_e = 1011u;
    s.es = 1000u; s.ee = 1011u; s.ecds = 1u; s.ecde = 12u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, 12u);
    s.vpos = 1000u; s.vend = 1000u; s.vkind = (uint8_t)DUCKVEP_KIND_SNV;
    s.abytes[0] = 'A'; s.abytes[1] = 'C';
    s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 1u;
    ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
              duckvep_variant_cds_edit_build(&s.tx, &s.ex, &s.seq, &s.v, 0u, 0u,
                                             s.strand, &edit));
    ASSERT_EQ(1u, edit.cds_start);
    ASSERT_EQ(1u, edit.ref_len);
    ASSERT_EQ(1u, edit.alt_len);
    ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
              duckvep_variant_cds_edit_set_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                 0u, 0u, s.strand, scratch, 1u,
                                                 &edit_set));
    ASSERT_EQ(1u, edit_set.count);
    ASSERT(edit_set.edits == scratch);
    ASSERT(kprop_cds_edits_equal(&edit, &scratch[0]));
    edit_set.edits = scratch;
    edit_set.count = 99u;
    ASSERT_EQ(DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL,
              duckvep_variant_cds_edit_set_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                 0u, 0u, s.strand, NULL, 0u,
                                                 &edit_set));
    ASSERT(edit_set.edits == NULL);
    ASSERT_EQ(0u, edit_set.count);

    memset(&s, 0, sizeof s);
    s.cds = plus_cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1011u; s.cds_s = 1000u; s.cds_e = 1011u;
    s.es = 1000u; s.ee = 1011u; s.ecds = 1u; s.ecde = 12u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, 12u);
    {
        uint8_t alt_mnv[5] = {'T','A','T','C','G'}; /* diff islands at CDS 4,6,8 */
        uint32_t i;
        s.vpos = kprop_genomic_pos_for_cds(&s, 4u); s.vend = s.vpos + 4u;
        s.vkind = (uint8_t)DUCKVEP_KIND_MNV;
        for (i = 0u; i < 5u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
        kprop_fill_variant_alt_from_tx(&s, 5u, alt_mnv, 5u);
        s.roff = 0u; s.aoff = 5u; s.rlen = 5u; s.alen = 5u;
        kprop_fill_expected_cds(&s, 4u, 5u, alt_mnv, 5u);
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
                  duckvep_variant_cds_edit_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                 0u, 0u, s.strand, &edit));
        ASSERT_EQ(5u, edit.ref_len);
        ASSERT_EQ(DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL,
                  duckvep_variant_cds_edit_set_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                     0u, 0u, s.strand, scratch, 2u,
                                                     &edit_set));
        ASSERT(edit_set.edits == NULL);
        ASSERT_EQ(0u, edit_set.count);
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
                  duckvep_variant_cds_edit_set_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                     0u, 0u, s.strand, scratch, 3u,
                                                     &edit_set));
        ASSERT_EQ(3u, edit_set.count);
        ASSERT_EQ(8u, edit_set.edits[0].cds_start);
        ASSERT_EQ(6u, edit_set.edits[1].cds_start);
        ASSERT_EQ(4u, edit_set.edits[2].cds_start);
        ASSERT_EQ(1u, edit_set.edits[0].ref_len);
        ASSERT_EQ(1u, edit_set.edits[1].ref_len);
        ASSERT_EQ(1u, edit_set.edits[2].ref_len);
        ASSERT(edit_set.edits[0].ref == s.abytes + 4u);
        ASSERT(edit_set.edits[1].ref == s.abytes + 2u);
        ASSERT(edit_set.edits[2].ref == s.abytes);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                  duckvep_haplotype_apply_cds_edits(s.cds, s.cds_lenv,
                                                    edit_set.edits, edit_set.count,
                                                    s.strand, mutated, sizeof mutated,
                                                    &mutated_len, &result));
        ASSERT_EQ((size_t)s.expect_len, mutated_len);
        ASSERT(memcmp(mutated, s.expect_cds, mutated_len) == 0);
        s.abytes[1] = kprop_base_not(s.abytes[1], (uint8_t)'N');
        s.abytes[6] = s.abytes[1]; /* retained internal mismatch must still validate */
        ASSERT_EQ(DUCKVEP_CDS_EDIT_REF_MISMATCH,
                  duckvep_variant_cds_edit_set_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                     0u, 0u, s.strand, scratch, 3u,
                                                     &edit_set));
        ASSERT(edit_set.edits == NULL);
        ASSERT_EQ(0u, edit_set.count);
    }

    memset(&s, 0, sizeof s);
    s.cds = minus_cds; s.chrom = 1u; s.strand = -1; s.flags = 0u;
    s.tstart = 2000u; s.tend = 2011u; s.cds_s = 2000u; s.cds_e = 2011u;
    s.es = 2000u; s.ee = 2011u; s.ecds = 1u; s.ecde = 12u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 1u;
    kprop_wire_coding_scene(&s, 12u);
    s.vpos = kprop_genomic_pos_for_cds(&s, 5u); s.vend = s.vpos;
    s.vkind = (uint8_t)DUCKVEP_KIND_INS;
    s.abytes[0] = (uint8_t)kprop_genomic_base_at(&s, s.vpos);
    s.abytes[1] = s.abytes[0];
    kprop_fill_variant_alt_from_tx(&s, 2u, payload, 2u);
    s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 3u;
    kprop_fill_expected_cds(&s, 5u, 0u, payload, 2u);
    ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
              duckvep_variant_cds_edit_build(&s.tx, &s.ex, &s.seq, &s.v, 0u, 0u,
                                             s.strand, &edit));
    ASSERT_EQ(5u, edit.cds_start);
    ASSERT_EQ(0u, edit.ref_len);
    ASSERT_EQ(2u, edit.alt_len);
    ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
              duckvep_variant_cds_edit_set_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                 0u, 0u, s.strand, scratch, 1u,
                                                 &edit_set));
    ASSERT_EQ(1u, edit_set.count);
    ASSERT(kprop_cds_edits_equal(&edit, &scratch[0]));
    ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
              duckvep_haplotype_apply_cds_edits(s.cds, s.cds_lenv,
                                                edit_set.edits, edit_set.count,
                                                s.strand, mutated, sizeof mutated,
                                                &mutated_len, &result));
    ASSERT_EQ((size_t)s.expect_len, mutated_len);
    ASSERT(memcmp(mutated, s.expect_cds, mutated_len) == 0);
    s.abytes[0] = kprop_base_not(s.abytes[0], (uint8_t)'N');
    s.abytes[1] = s.abytes[0];
    ASSERT_EQ(DUCKVEP_CDS_EDIT_REF_MISMATCH,
              duckvep_variant_cds_edit_build(&s.tx, &s.ex, &s.seq, &s.v, 0u, 0u,
                                             s.strand, &edit));
    edit_set.edits = scratch;
    edit_set.count = 99u;
    ASSERT_EQ(DUCKVEP_CDS_EDIT_REF_MISMATCH,
              duckvep_variant_cds_edit_set_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                 0u, 0u, s.strand, scratch, 1u,
                                                 &edit_set));
    ASSERT(edit_set.edits == NULL);
    ASSERT_EQ(0u, edit_set.count);

    memset(&s, 0, sizeof s);
    s.cds = minus_cds; s.chrom = 1u; s.strand = -1; s.flags = 0u;
    s.tstart = 2000u; s.tend = 2011u; s.cds_s = 2000u; s.cds_e = 2011u;
    s.es = 2000u; s.ee = 2011u; s.ecds = 1u; s.ecde = 12u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 1u;
    kprop_wire_coding_scene(&s, 12u);
    {
        uint8_t alt_mnv[5] = {'T','A','T','C','G'}; /* transcript-order diff CDS 4,6,8 */
        uint32_t i;
        s.vpos = kprop_genomic_pos_for_cds(&s, 8u); s.vend = s.vpos + 4u;
        s.vkind = (uint8_t)DUCKVEP_KIND_MNV;
        for (i = 0u; i < 5u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
        kprop_fill_variant_alt_from_tx(&s, 5u, alt_mnv, 5u);
        s.roff = 0u; s.aoff = 5u; s.rlen = 5u; s.alen = 5u;
        kprop_fill_expected_cds(&s, 4u, 5u, alt_mnv, 5u);
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OK,
                  duckvep_variant_cds_edit_set_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                     0u, 0u, s.strand, scratch, 3u,
                                                     &edit_set));
        ASSERT_EQ(3u, edit_set.count);
        ASSERT_EQ(8u, edit_set.edits[0].cds_start);
        ASSERT_EQ(6u, edit_set.edits[1].cds_start);
        ASSERT_EQ(4u, edit_set.edits[2].cds_start);
        ASSERT(edit_set.edits[0].ref == s.abytes);
        ASSERT(edit_set.edits[1].ref == s.abytes + 2u);
        ASSERT(edit_set.edits[2].ref == s.abytes + 4u);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                  duckvep_haplotype_apply_cds_edits(s.cds, s.cds_lenv,
                                                    edit_set.edits, edit_set.count,
                                                    s.strand, mutated, sizeof mutated,
                                                    &mutated_len, &result));
        ASSERT_EQ((size_t)s.expect_len, mutated_len);
        ASSERT(memcmp(mutated, s.expect_cds, mutated_len) == 0);
    }

    {
        static const uint16_t tchrom[1] = {2u};
        static const uint32_t tstart[1] = {3000u};
        static const uint32_t tend[1] = {3014u};
        static const int8_t tstrand[1] = {1};
        static const uint64_t tflags[1] = {0u};
        static const uint32_t texoff[1] = {0u};
        static const uint16_t texcnt[1] = {2u};
        static const uint32_t tcds_s[1] = {3000u};
        static const uint32_t tcds_e[1] = {3014u};
        static const uint32_t estart[2] = {3000u, 3010u};
        static const uint32_t eend[2] = {3002u, 3014u};
        static const uint32_t ecdna_s[2] = {1u, 4u};
        static const uint32_t ecdna_e[2] = {3u, 8u};
        static const int8_t ephase[2] = {0, 0};
        static const uint8_t cds[8] = {'A','T','G','A','A','A','C','C'};
        static const uint64_t cds_off[1] = {0u};
        static const uint32_t cds_len[1] = {8u};
        static const uint8_t abytes[18] = {
            'G','A','A','A','A','A','A','A','T',
            'C','C','C','C','C','C','C','C','C'
        };
        static const uint16_t vchrom[1] = {2u};
        static const uint32_t vpos[1] = {3002u};
        static const uint32_t vend[1] = {3010u};
        static const uint32_t roff[1] = {0u};
        static const uint32_t aoff[1] = {9u};
        static const uint16_t rlen[1] = {9u};
        static const uint16_t alen[1] = {9u};
        static const uint8_t vkind[1] = {(uint8_t)DUCKVEP_KIND_MNV};
        duckvep_transcript_model_t tx;
        duckvep_exon_model_t ex;
        duckvep_sequence_pool_t seq;
        duckvep_variant_batch_t v;
        memset(&tx, 0, sizeof tx); memset(&ex, 0, sizeof ex);
        memset(&seq, 0, sizeof seq); memset(&v, 0, sizeof v);
        tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = tstrand;
        tx.flags = tflags; tx.exon_offset = texoff; tx.exon_count = texcnt;
        tx.cds_start1 = tcds_s; tx.cds_end1 = tcds_e; tx.transcript_count = 1u;
        ex.start1 = estart; ex.end1 = eend; ex.cdna_start1 = ecdna_s; ex.cdna_end1 = ecdna_e;
        ex.phase = ephase; ex.end_phase = ephase; ex.exon_count = 2u;
        seq.cds_bytes = cds; seq.cds_bytes_len = sizeof cds;
        seq.cds_offset = cds_off; seq.cds_length = cds_len; seq.transcript_count = 1u;
        v.chrom_id = vchrom; v.pos1 = vpos; v.end1 = vend; v.ref_offset = roff; v.alt_offset = aoff;
        v.ref_length = rlen; v.alt_length = alen; v.allele_bytes = abytes;
        v.allele_bytes_len = sizeof abytes; v.variant_kind = vkind; v.count = 1u;
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OUT_OF_CDS,
                  duckvep_variant_cds_edit_build(&tx, &ex, &seq, &v, 0u, 0u, 1, &edit));
        edit_set.edits = scratch;
        edit_set.count = 99u;
        ASSERT_EQ(DUCKVEP_CDS_EDIT_OUT_OF_CDS,
                  duckvep_variant_cds_edit_set_build(&tx, &ex, &seq, &v, 0u, 0u,
                                                     1, scratch, 1u, &edit_set));
        ASSERT(edit_set.edits == NULL);
        ASSERT_EQ(0u, edit_set.count);
    }

    {
        duckvep_haplotype_edit_t cap_scratch[1];
        edit_set.edits = cap_scratch;
        edit_set.count = 99u;
        ASSERT_EQ(DUCKVEP_CDS_EDIT_INVALID_ARG,
                  duckvep_variant_cds_edit_set_build(NULL, &s.ex, &s.seq, &s.v,
                                                     0u, 0u, s.strand,
                                                     cap_scratch, 1u, &edit_set));
        ASSERT(edit_set.edits == NULL);
        ASSERT_EQ(0u, edit_set.count);
    }

    PASS();
}

TEST coding_context_known_scene(void) {
    {
        static const uint8_t ref_cds[9] = {'A','T','G','T','A','A','G','A','A'};
        static const uint8_t ref_base[1] = {'A'};
        static const uint8_t alt_base[1] = {'C'};
        duckvep_haplotype_edit_t edit;
        duckvep_edit_set_t edit_set;
        duckvep_coding_context_t ctx;
        duckvep_haplotype_result_t hres;
        uint8_t alt_cds[16];
        uint8_t ref_pep[8];
        uint8_t alt_pep[8];
        size_t trunc_len = 0u;

        edit.cds_start = 8u; edit.ref_len = 1u; edit.ref = ref_base;
        edit.alt_len = 1u; edit.alt = alt_base; edit.variant_strand = 1;
        edit_set.edits = &edit; edit_set.count = 1u;
        memset(ref_pep, 0xA5, sizeof ref_pep);
        memset(alt_pep, 0x5A, sizeof alt_pep);
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(ref_cds, sizeof ref_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ((size_t)9u, ctx.ref_cds_len);
        ASSERT_EQ((size_t)9u, ctx.alt_cds_len);
        ASSERT_EQ((size_t)3u, ctx.ref_peptide_len);
        ASSERT_EQ((size_t)3u, ctx.alt_peptide_len);
        ASSERT(memcmp(ctx.ref_peptide, "M*E", 3u) == 0);
        ASSERT(memcmp(ctx.alt_peptide, "M*A", 3u) == 0);
        ASSERT_EQ((uint8_t)'\0', ctx.ref_peptide[ctx.ref_peptide_len]);
        ASSERT_EQ((uint8_t)'\0', ctx.alt_peptide[ctx.alt_peptide_len]);
        ASSERT_EQ(1u, ctx.cds_changed);
        ASSERT_EQ((size_t)1u, ctx.applied_edits);
        ASSERT_EQ(0, ctx.length_diff);
        ASSERT_EQ(3u, ctx.ref_first_changed_codon);
        ASSERT_EQ(3u, ctx.ref_last_changed_codon);
        ASSERT_EQ(3u, ctx.alt_first_changed_codon);
        ASSERT_EQ(3u, ctx.alt_last_changed_codon);
        ASSERT_EQ(DUCKVEP_HAPLOTYPE_OK,
                  duckvep_haplotype_translate_cds(ctx.alt_cds, ctx.alt_cds_len,
                                                  DUCKVEP_CODON_TABLE_STANDARD,
                                                  alt_pep, sizeof alt_pep,
                                                  &trunc_len, &hres));
        ASSERT_EQ((size_t)2u, trunc_len);
        ASSERT(memcmp(alt_pep, "M*", 2u) == 0);
        ASSERT((hres.flags & DUCKVEP_HAPLOTYPE_FLAG_STOP_TRUNCATED) != 0u);

        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_ALT_CDS_BUFFER_TOO_SMALL,
                  duckvep_coding_context_build(ref_cds, sizeof ref_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, 8u, ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT(kprop_context_zero(&ctx));
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_REF_PEPTIDE_BUFFER_TOO_SMALL,
                  duckvep_coding_context_build(ref_cds, sizeof ref_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep, 3u,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT(kprop_context_zero(&ctx));
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_ALT_PEPTIDE_BUFFER_TOO_SMALL,
                  duckvep_coding_context_build(ref_cds, sizeof ref_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep, alt_pep, 3u,
                                               &ctx));
        ASSERT(kprop_context_zero(&ctx));
    }

    {
        static const uint8_t ref_cds[6] = {'G','C','T','G','A','A'};
        static const uint8_t ref_base[1] = {'T'};
        static const uint8_t alt_base[1] = {'C'};
        duckvep_haplotype_edit_t edit;
        duckvep_edit_set_t edit_set;
        duckvep_coding_context_t ctx;
        uint8_t alt_cds[8];
        uint8_t ref_pep[8];
        uint8_t alt_pep[8];

        edit.cds_start = 3u; edit.ref_len = 1u; edit.ref = ref_base;
        edit.alt_len = 1u; edit.alt = alt_base; edit.variant_strand = 1;
        edit_set.edits = &edit; edit_set.count = 1u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(ref_cds, sizeof ref_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(1u, ctx.cds_changed);
        ASSERT_EQ((size_t)2u, ctx.ref_peptide_len);
        ASSERT_EQ((size_t)2u, ctx.alt_peptide_len);
        ASSERT(memcmp(ctx.ref_peptide, ctx.alt_peptide, 2u) == 0);
        ASSERT_EQ(0u, ctx.ref_first_changed_codon);
        ASSERT_EQ(0u, ctx.ref_last_changed_codon);
        ASSERT_EQ(0u, ctx.alt_first_changed_codon);
        ASSERT_EQ(0u, ctx.alt_last_changed_codon);
    }

    {
        static const uint8_t ref_cds[8] = {'A','T','G','G','A','A','C','T'};
        static const uint8_t ref_base[1] = {'T'};
        static const uint8_t alt_base[1] = {'A'};
        duckvep_haplotype_edit_t edit;
        duckvep_edit_set_t edit_set;
        duckvep_coding_context_t ctx;
        uint8_t alt_cds[8];
        uint8_t ref_pep[8];
        uint8_t alt_pep[8];

        edit.cds_start = 8u; edit.ref_len = 1u; edit.ref = ref_base;
        edit.alt_len = 1u; edit.alt = alt_base; edit.variant_strand = 1;
        edit_set.edits = &edit; edit_set.count = 1u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(ref_cds, sizeof ref_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ((size_t)8u, ctx.alt_cds_len);
        ASSERT_EQ((size_t)2u, ctx.ref_peptide_len);
        ASSERT_EQ((size_t)2u, ctx.alt_peptide_len);
        ASSERT(memcmp(ctx.ref_peptide, "ME", 2u) == 0);
        ASSERT(memcmp(ctx.alt_peptide, "ME", 2u) == 0);
        ASSERT_EQ(1u, ctx.cds_changed);
        ASSERT_EQ(0u, ctx.ref_first_changed_codon);
        ASSERT_EQ(0u, ctx.alt_first_changed_codon);
    }

    {
        static const uint8_t ref_cds[6] = {'A','T','N','G','A','A'};
        duckvep_edit_set_t edit_set;
        duckvep_coding_context_t ctx;
        uint8_t alt_cds[8];
        uint8_t ref_pep[8];
        uint8_t alt_pep[8];
        edit_set.edits = NULL; edit_set.count = 0u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(ref_cds, sizeof ref_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(0u, ctx.cds_changed);
        ASSERT(memcmp(ctx.ref_peptide, "XE", 2u) == 0);
        ASSERT(memcmp(ctx.alt_peptide, "XE", 2u) == 0);
    }

    PASS();
}

TEST variant_coding_context_known_scene(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'T','A','A',  'G','A','A',  'C','C','C',  'T','T','T'
    };
    struct kprop_coding s;
    duckvep_haplotype_edit_t edit_scratch[3];
    duckvep_coding_context_t ctx;
    uint8_t alt_cds[24];
    uint8_t ref_pep[16];
    uint8_t alt_pep[16];
    uint8_t alt_mnv[5] = {'A','A','C','C','G'};
    uint32_t i;

    memset(&s, 0, sizeof s);
    s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
    s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, 15u);
    s.vpos = kprop_genomic_pos_for_cds(&s, 7u); s.vend = s.vpos + 4u;
    s.vkind = (uint8_t)DUCKVEP_KIND_MNV;
    for (i = 0u; i < 5u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
    kprop_fill_variant_alt_from_tx(&s, 5u, alt_mnv, 5u);
    s.roff = 0u; s.aoff = 5u; s.rlen = 5u; s.alen = 5u;
    kprop_fill_expected_cds(&s, 7u, 5u, alt_mnv, 5u);

    ctx.ref_cds = cds;
    ctx.ref_cds_len = 99u;
    ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_EDIT_BUFFER_TOO_SMALL,
              duckvep_variant_coding_context_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                   0u, 0u, s.strand,
                                                   edit_scratch, 2u,
                                                   alt_cds, sizeof alt_cds,
                                                   ref_pep, sizeof ref_pep,
                                                   alt_pep, sizeof alt_pep,
                                                   &ctx));
    ASSERT(kprop_context_zero(&ctx));

    memset(ref_pep, 0xA5, sizeof ref_pep);
    memset(alt_pep, 0x5A, sizeof alt_pep);
    ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_OK,
              duckvep_variant_coding_context_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                   0u, 0u, s.strand,
                                                   edit_scratch, 3u,
                                                   alt_cds, sizeof alt_cds,
                                                   ref_pep, sizeof ref_pep,
                                                   alt_pep, sizeof alt_pep,
                                                   &ctx));
    ASSERT(ctx.ref_cds == s.cds);
    ASSERT(ctx.alt_cds == alt_cds);
    ASSERT_EQ((size_t)15u, ctx.ref_cds_len);
    ASSERT_EQ((size_t)15u, ctx.alt_cds_len);
    ASSERT(memcmp(ctx.alt_cds, s.expect_cds, (size_t)s.expect_len) == 0);
    ASSERT_EQ((size_t)5u, ctx.ref_peptide_len);
    ASSERT_EQ((size_t)5u, ctx.alt_peptide_len);
    ASSERT(memcmp(ctx.ref_peptide, "M*EPF", 5u) == 0);
    ASSERT(memcmp(ctx.alt_peptide, "M*NRF", 5u) == 0);
    ASSERT_EQ((uint8_t)'\0', ctx.ref_peptide[ctx.ref_peptide_len]);
    ASSERT_EQ((uint8_t)'\0', ctx.alt_peptide[ctx.alt_peptide_len]);
    ASSERT_EQ((size_t)3u, ctx.applied_edits);
    ASSERT_EQ(1u, ctx.cds_changed);
    ASSERT_EQ(3u, ctx.ref_first_changed_codon);
    ASSERT_EQ(4u, ctx.ref_last_changed_codon);
    ASSERT_EQ(3u, ctx.alt_first_changed_codon);
    ASSERT_EQ(4u, ctx.alt_last_changed_codon);

    {
        static uint8_t mito_cds[6] = {'A','T','A','G','A','A'};
        uint8_t ref_base[1] = {'A'};
        uint8_t alt_base[1] = {'G'};
        memset(&s, 0, sizeof s);
        s.cds = mito_cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
        s.tstart = 2000u; s.tend = 2005u; s.cds_s = 2000u; s.cds_e = 2005u;
        s.es = 2000u; s.ee = 2005u; s.ecds = 1u; s.ecde = 6u; s.eph = 0; s.eeph = 0;
        s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
        kprop_wire_coding_scene(&s, 6u);
        s.ctab = (uint8_t)DUCKVEP_CODON_TABLE_VERT_MITO;
        s.vpos = kprop_genomic_pos_for_cds(&s, 6u); s.vend = s.vpos;
        s.vkind = (uint8_t)DUCKVEP_KIND_SNV;
        s.abytes[0] = ref_base[0]; s.abytes[1] = alt_base[0];
        s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 1u;
        ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_OK,
                  duckvep_variant_coding_context_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                       0u, 0u, s.strand,
                                                       edit_scratch, 1u,
                                                       alt_cds, sizeof alt_cds,
                                                       ref_pep, sizeof ref_pep,
                                                       alt_pep, sizeof alt_pep,
                                                       &ctx));
        ASSERT(memcmp(ctx.ref_peptide, "ME", 2u) == 0);
    }

    PASS();
}

TEST coding_context_delta_known_scene(void) {
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;

    {
        static const uint8_t ref_cds[6] = {'A','T','G','G','A','A'};
        static const uint8_t alt_cds[6] = {'A','T','G','G','A','G'};
        static const uint8_t pep[3] = {'M','E','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 6u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 6u;
        ctx.ref_peptide = pep; ctx.ref_peptide_len = 2u;
        ctx.alt_peptide = pep; ctx.alt_peptide_len = 2u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.synonymous && !delta.start_lost);
        ASSERT_EQ(2, delta.protein_pos);
        ASSERT_EQ((uint8_t)'E', delta.ref_aa);
        ASSERT_EQ((uint8_t)'E', delta.alt_aa);
    }

    {
        static const uint8_t ref_cds[6] = {'A','T','G','G','A','A'};
        static const uint8_t alt_cds[6] = {'A','T','G','G','A','C'};
        static const uint8_t ref_pep[3] = {'M','E','\0'};
        static const uint8_t alt_pep[3] = {'M','D','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 6u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 6u;
        ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 2u;
        ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 2u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.missense && !delta.start_lost);
        ASSERT_EQ((uint8_t)'E', delta.ref_aa);
        ASSERT_EQ((uint8_t)'D', delta.alt_aa);
    }

    {
        static const uint8_t ref_cds[6] = {'A','T','G','T','G','G'};
        static const uint8_t alt_cds[6] = {'A','T','G','T','G','A'};
        static const uint8_t ref_pep[3] = {'M','W','\0'};
        static const uint8_t alt_pep[3] = {'M','*','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 6u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 6u;
        ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 2u;
        ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 2u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_gained && !delta.start_lost);
    }

    {
        static const uint8_t ref_cds[6] = {'A','T','G','T','A','A'};
        static const uint8_t alt_cds[6] = {'A','T','G','C','A','A'};
        static const uint8_t ref_pep[3] = {'M','*','\0'};
        static const uint8_t alt_pep[3] = {'M','Q','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 6u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 6u;
        ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 2u;
        ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 2u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_lost && !delta.start_lost);
    }

    {
        static const uint8_t ref_cds[6] = {'A','T','G','T','A','A'};
        static const uint8_t alt_cds[6] = {'A','T','G','T','A','G'};
        static const uint8_t pep[3] = {'M','*','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 6u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 6u;
        ctx.ref_peptide = pep; ctx.ref_peptide_len = 2u;
        ctx.alt_peptide = pep; ctx.alt_peptide_len = 2u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_retained && !delta.start_lost);
    }

    {
        static const uint8_t ref_cds[3] = {'A','T','G'};
        static const uint8_t alt_cds[3] = {'G','T','G'};
        static const uint8_t ref_pep[2] = {'M','\0'};
        static const uint8_t alt_pep[2] = {'V','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 3u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 3u;
        ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 1u;
        ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 1u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.start_lost && !delta.missense);
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx,
                                                    (uint64_t)DUCKVEP_TX_CDS_START_NF,
                                                    &delta));
        ASSERT(delta.valid && !delta.start_lost && delta.missense);
    }

    {
        /* Two-codon body substitution (codons 2-3, past the start): the generalized window
         * classifier resolves it as missense_variant — EE -> DD — where the old two-codon
         * slice bailed to coding_sequence_variant. Coarse window: protein_pos -1, no AA pair. */
        static const uint8_t ref_cds[9] = {'A','T','G','G','A','A','G','A','A'};
        static const uint8_t alt_cds[9] = {'A','T','G','G','A','C','G','A','C'};
        static const uint8_t ref_pep[4] = {'M','E','E','\0'};
        static const uint8_t alt_pep[4] = {'M','D','D','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 9u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 9u;
        ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 3u;
        ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 3u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.missense && !delta.synonymous && !delta.stop_gained &&
               !delta.stop_lost && !delta.stop_retained && !delta.start_lost);
        ASSERT(delta.protein_pos == -1 && delta.ref_aa == (uint8_t)0u &&
               delta.alt_aa == (uint8_t)0u);
    }

    {
        static const uint8_t ref_cds[3] = {'G','A','N'};
        static const uint8_t alt_cds[3] = {'G','A','C'};
        static const uint8_t ref_pep[2] = {'X','\0'};
        static const uint8_t alt_pep[2] = {'D','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 3u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 3u;
        ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 1u;
        ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 1u;
        ctx.cds_changed = 1u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(!delta.valid);
        /* Do not rely on peptide X alone: inconsistent callers with non-X peptide bytes
         * still cannot classify an ambiguous changed codon. */
        {
            static const uint8_t fake_ref_pep[2] = {'E','\0'};
            memset(&ctx, 0, sizeof ctx);
            ctx.ref_cds = ref_cds; ctx.ref_cds_len = 3u;
            ctx.alt_cds = alt_cds; ctx.alt_cds_len = 3u;
            ctx.ref_peptide = fake_ref_pep; ctx.ref_peptide_len = 1u;
            ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 1u;
            ctx.cds_changed = 1u;
        }
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(!delta.valid);
    }

    {
        static const uint8_t ref_cds[3] = {'G','A','A'};
        static const uint8_t alt_cds[4] = {'G','A','A','A'};
        static const uint8_t pep[2] = {'E','\0'};
        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ref_cds; ctx.ref_cds_len = 3u;
        ctx.alt_cds = alt_cds; ctx.alt_cds_len = 4u;
        ctx.ref_peptide = pep; ctx.ref_peptide_len = 1u;
        ctx.alt_peptide = pep; ctx.alt_peptide_len = 1u;
        ctx.cds_changed = 1u;
        ctx.length_diff = 1;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(!delta.valid);
    }

    PASS();
}

static struct {
    uint32_t snv;
    uint32_t mnv;
    uint32_t ins;
    uint32_t del;
    uint32_t indel;
    uint32_t fwd;
    uint32_t rev;
    uint32_t start;
    uint32_t body;
    uint32_t stop;
} g_cds_edit_builder_cov;

static enum theft_trial_res prop_cds_edit_builder_matches_splice_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edit;
    duckvep_haplotype_result_t result;
    duckvep_event_t event;
    uint32_t projected_cds_start;
    uint32_t projected_cds_end;
    uint8_t mutated[64];
    size_t mutated_len = 0u;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;

    {
        duckvep_cds_edit_status_t st = duckvep_variant_cds_edit_build(
            &s->tx, &s->ex, &s->seq, &s->v, 0u, 0u, s->strand, &edit);
        if (st != DUCKVEP_CDS_EDIT_OK) {
            fprintf(stderr, "\n[cds-edit-builder fail] build status=%d shape=%u strand=%d pos=%u end=%u ref_len=%u alt_len=%u\n",
                    (int)st, (unsigned)s->expect_shape, (int)s->strand,
                    s->vpos, s->vend, (unsigned)s->rlen, (unsigned)s->alen);
            return THEFT_TRIAL_FAIL;
        }
    }
    if (!duckvep_event_prepare_small(
            s->vpos, s->v.allele_bytes + s->roff, s->rlen,
            s->v.allele_bytes + s->aoff, s->alen, &event) ||
        !duckvep_project_event_to_cds(&s->tx, &s->ex, 0u, &event,
                                      &projected_cds_start,
                                      &projected_cds_end) ||
        projected_cds_start != edit.cds_start ||
        projected_cds_end != edit.cds_start +
            (edit.ref_len == 0u ? 0u : edit.ref_len - 1u)) {
        fprintf(stderr,
                "\n[cds-edit-builder fail] shared event projection disagrees shape=%u strand=%d edit=%u+%u projected=%u-%u\n",
                (unsigned)s->expect_shape, (int)s->strand,
                edit.cds_start, edit.ref_len, projected_cds_start,
                projected_cds_end);
        return THEFT_TRIAL_FAIL;
    }
    {
        duckvep_haplotype_status_t hst = duckvep_haplotype_apply_cds_edits(
            s->cds, s->cds_lenv, &edit, 1u, s->strand, mutated, sizeof mutated,
            &mutated_len, &result);
        if (hst != DUCKVEP_HAPLOTYPE_OK) {
            fprintf(stderr, "\n[cds-edit-builder fail] haplo status=%d shape=%u strand=%d cds_start=%u ref_len=%u alt_len=%u\n",
                    (int)hst, (unsigned)s->expect_shape, (int)s->strand,
                    edit.cds_start, edit.ref_len, edit.alt_len);
            return THEFT_TRIAL_FAIL;
        }
    }
    if (mutated_len != (size_t)s->expect_len ||
        memcmp(mutated, s->expect_cds, (size_t)s->expect_len) != 0) {
        fprintf(stderr,
                "\n[cds-edit-builder fail] mismatch shape=%u strand=%d vpos=%u end=%u edit_start=%u ref_len=%u alt_len=%u got_len=%zu want_len=%u\n",
                (unsigned)s->expect_shape, (int)s->strand, s->vpos, s->vend,
                edit.cds_start, edit.ref_len, edit.alt_len, mutated_len,
                (unsigned)s->expect_len);
        tr = THEFT_TRIAL_FAIL;
        goto done;
    }
    if (s->expect_shape == KPROP_CDS_EDIT_SNV) g_cds_edit_builder_cov.snv++;
    else if (s->expect_shape == KPROP_CDS_EDIT_MNV) g_cds_edit_builder_cov.mnv++;
    else if (s->expect_shape == KPROP_CDS_EDIT_INS) g_cds_edit_builder_cov.ins++;
    else if (s->expect_shape == KPROP_CDS_EDIT_DEL) g_cds_edit_builder_cov.del++;
    else if (s->expect_shape == KPROP_CDS_EDIT_INDEL) g_cds_edit_builder_cov.indel++;
    else tr = THEFT_TRIAL_FAIL;
    if (s->strand > 0) g_cds_edit_builder_cov.fwd++; else g_cds_edit_builder_cov.rev++;
    if (s->expect_region == KPROP_CDS_EDIT_START) g_cds_edit_builder_cov.start++;
    else if (s->expect_region == KPROP_CDS_EDIT_BODY) g_cds_edit_builder_cov.body++;
    else if (s->expect_region == KPROP_CDS_EDIT_STOP) g_cds_edit_builder_cov.stop++;
    else tr = THEFT_TRIAL_FAIL;
done:
    return tr;
}

TEST cds_edit_builder_matches_direct_splice_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "variant CDS edit builder == direct CDS splice oracle";
    cfg.prop1 = prop_cds_edit_builder_matches_splice_oracle;
    cfg.type_info[0] = &kprop_cds_edit_builder_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_cds_edit_builder_cov, 0, sizeof g_cds_edit_builder_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_cds_edit_builder_cov.snv > 0u);
    ASSERT(g_cds_edit_builder_cov.mnv > 0u);
    ASSERT(g_cds_edit_builder_cov.ins > 0u);
    ASSERT(g_cds_edit_builder_cov.del > 0u);
    ASSERT(g_cds_edit_builder_cov.indel > 0u);
    ASSERT(g_cds_edit_builder_cov.fwd > 0u);
    ASSERT(g_cds_edit_builder_cov.rev > 0u);
    ASSERT(g_cds_edit_builder_cov.start > 0u);
    ASSERT(g_cds_edit_builder_cov.body > 0u);
    ASSERT(g_cds_edit_builder_cov.stop > 0u);
    fprintf(stderr,
            "[cds-edit-builder coverage] snv=%u mnv=%u ins=%u del=%u indel=%u fwd=%u rev=%u start=%u body=%u stop=%u\n",
            g_cds_edit_builder_cov.snv, g_cds_edit_builder_cov.mnv,
            g_cds_edit_builder_cov.ins, g_cds_edit_builder_cov.del,
            g_cds_edit_builder_cov.indel, g_cds_edit_builder_cov.fwd,
            g_cds_edit_builder_cov.rev, g_cds_edit_builder_cov.start,
            g_cds_edit_builder_cov.body, g_cds_edit_builder_cov.stop);
    PASS();
}

static struct {
    uint32_t snv;
    uint32_t mnv;
    uint32_t ins;
    uint32_t del;
    uint32_t indel;
    uint32_t fwd;
    uint32_t rev;
    uint32_t start;
    uint32_t body;
    uint32_t stop;
    uint32_t cap0;
} g_cds_edit_set_cov;

static enum theft_trial_res prop_cds_edit_set_matches_single_builder(struct theft *t,
                                                                     void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edit;
    duckvep_haplotype_edit_t scratch[1];
    duckvep_edit_set_t edit_set;
    duckvep_haplotype_result_t result;
    uint8_t mutated[64];
    size_t mutated_len = 0u;
    duckvep_cds_edit_status_t single_st;
    duckvep_cds_edit_status_t set_st;
    (void)t;

    single_st = duckvep_variant_cds_edit_build(&s->tx, &s->ex, &s->seq, &s->v,
                                               0u, 0u, s->strand, &edit);
    edit_set.edits = scratch;
    edit_set.count = 99u;
    set_st = duckvep_variant_cds_edit_set_build(&s->tx, &s->ex, &s->seq, &s->v,
                                                0u, 0u, s->strand, scratch, 1u,
                                                &edit_set);
    if (single_st != set_st) return THEFT_TRIAL_FAIL;
    if (set_st != DUCKVEP_CDS_EDIT_OK) {
        return edit_set.edits == NULL && edit_set.count == 0u
            ? THEFT_TRIAL_PASS : THEFT_TRIAL_FAIL;
    }
    if (edit_set.edits != scratch || edit_set.count != 1u ||
        !kprop_cds_edits_equal(&edit, &scratch[0])) {
        return THEFT_TRIAL_FAIL;
    }
    {
        duckvep_edit_set_t too_small;
        too_small.edits = scratch;
        too_small.count = 99u;
        if (duckvep_variant_cds_edit_set_build(&s->tx, &s->ex, &s->seq, &s->v,
                                               0u, 0u, s->strand, NULL, 0u,
                                               &too_small) != DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (too_small.edits != NULL || too_small.count != 0u) return THEFT_TRIAL_FAIL;
        g_cds_edit_set_cov.cap0++;
    }
    if (duckvep_haplotype_apply_cds_edits(s->cds, s->cds_lenv, edit_set.edits,
                                          edit_set.count, s->strand, mutated,
                                          sizeof mutated, &mutated_len,
                                          &result) != DUCKVEP_HAPLOTYPE_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (mutated_len != (size_t)s->expect_len ||
        memcmp(mutated, s->expect_cds, (size_t)s->expect_len) != 0) {
        return THEFT_TRIAL_FAIL;
    }
    if (s->expect_shape == KPROP_CDS_EDIT_SNV) g_cds_edit_set_cov.snv++;
    else if (s->expect_shape == KPROP_CDS_EDIT_MNV) g_cds_edit_set_cov.mnv++;
    else if (s->expect_shape == KPROP_CDS_EDIT_INS) g_cds_edit_set_cov.ins++;
    else if (s->expect_shape == KPROP_CDS_EDIT_DEL) g_cds_edit_set_cov.del++;
    else if (s->expect_shape == KPROP_CDS_EDIT_INDEL) g_cds_edit_set_cov.indel++;
    else return THEFT_TRIAL_FAIL;
    if (s->strand > 0) g_cds_edit_set_cov.fwd++; else g_cds_edit_set_cov.rev++;
    if (s->expect_region == KPROP_CDS_EDIT_START) g_cds_edit_set_cov.start++;
    else if (s->expect_region == KPROP_CDS_EDIT_BODY) g_cds_edit_set_cov.body++;
    else if (s->expect_region == KPROP_CDS_EDIT_STOP) g_cds_edit_set_cov.stop++;
    else return THEFT_TRIAL_FAIL;
    return THEFT_TRIAL_PASS;
}

TEST cds_edit_set_builder_matches_single_edit_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "variant CDS edit-set builder == single-edit splice oracle";
    cfg.prop1 = prop_cds_edit_set_matches_single_builder;
    cfg.type_info[0] = &kprop_cds_edit_builder_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_cds_edit_set_cov, 0, sizeof g_cds_edit_set_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_cds_edit_set_cov.snv > 0u);
    ASSERT(g_cds_edit_set_cov.mnv > 0u);
    ASSERT(g_cds_edit_set_cov.ins > 0u);
    ASSERT(g_cds_edit_set_cov.del > 0u);
    ASSERT(g_cds_edit_set_cov.indel > 0u);
    ASSERT(g_cds_edit_set_cov.fwd > 0u);
    ASSERT(g_cds_edit_set_cov.rev > 0u);
    ASSERT(g_cds_edit_set_cov.start > 0u);
    ASSERT(g_cds_edit_set_cov.body > 0u);
    ASSERT(g_cds_edit_set_cov.stop > 0u);
    ASSERT(g_cds_edit_set_cov.cap0 > 0u);
    fprintf(stderr,
            "[cds-edit-set coverage] snv=%u mnv=%u ins=%u del=%u indel=%u fwd=%u rev=%u start=%u body=%u stop=%u cap0=%u\n",
            g_cds_edit_set_cov.snv, g_cds_edit_set_cov.mnv,
            g_cds_edit_set_cov.ins, g_cds_edit_set_cov.del,
            g_cds_edit_set_cov.indel, g_cds_edit_set_cov.fwd,
            g_cds_edit_set_cov.rev, g_cds_edit_set_cov.start,
            g_cds_edit_set_cov.body, g_cds_edit_set_cov.stop,
            g_cds_edit_set_cov.cap0);
    PASS();
}

static struct {
    uint32_t fwd;
    uint32_t rev;
    uint32_t start;
    uint32_t body;
    uint32_t stop;
    uint32_t multi;
    uint32_t capfail;
} g_cds_edit_set_mnv_cov;

static enum theft_trial_res prop_cds_edit_set_splits_mnv_islands(struct theft *t,
                                                                 void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t scratch[8];
    duckvep_edit_set_t edit_set;
    duckvep_haplotype_result_t result;
    uint8_t mutated[64];
    uint8_t covered[64];
    uint32_t exp_start[64];
    uint32_t exp_len[64];
    size_t exp_n = 0u;
    size_t mutated_len = 0u;
    size_t i;
    (void)t;

    memset(covered, 0, sizeof covered);
    if (duckvep_variant_cds_edit_set_build(&s->tx, &s->ex, &s->seq, &s->v,
                                           0u, 0u, s->strand, scratch, 8u,
                                           &edit_set) != DUCKVEP_CDS_EDIT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (edit_set.count < 2u) return THEFT_TRIAL_FAIL;
    i = 0u;
    while (i < s->cds_lenv) {
        uint32_t start;
        uint32_t len;
        while (i < s->cds_lenv && s->cds[i] == s->expect_cds[i]) i++;
        if (i >= s->cds_lenv) break;
        start = (uint32_t)i;
        len = 0u;
        while (i < s->cds_lenv && s->cds[i] != s->expect_cds[i]) {
            i++;
            len++;
        }
        exp_start[exp_n] = start;
        exp_len[exp_n] = len;
        exp_n++;
    }
    if (edit_set.count != exp_n) return THEFT_TRIAL_FAIL;
    for (i = 0u; i < edit_set.count; i++) {
        uint32_t j;
        size_t exp_i = exp_n - 1u - i;
        if (edit_set.edits[i].cds_start != exp_start[exp_i] + 1u ||
            edit_set.edits[i].ref_len != exp_len[exp_i] ||
            edit_set.edits[i].alt_len != exp_len[exp_i]) {
            return THEFT_TRIAL_FAIL;
        }
        if (i + 1u < edit_set.count &&
            edit_set.edits[i].cds_start <= edit_set.edits[i + 1u].cds_start) {
            return THEFT_TRIAL_FAIL;
        }
        for (j = 0u; j < edit_set.edits[i].ref_len; j++) {
            uint32_t cds_pos = edit_set.edits[i].cds_start + j;
            if (cds_pos == 0u || cds_pos > s->cds_lenv) return THEFT_TRIAL_FAIL;
            covered[cds_pos - 1u] = 1u;
        }
    }
    for (i = 0u; i < s->cds_lenv; i++) {
        int diff = s->cds[i] != s->expect_cds[i];
        if ((covered[i] != 0u) != diff) return THEFT_TRIAL_FAIL;
    }
    if (edit_set.count > 1u) {
        duckvep_edit_set_t fail_set;
        fail_set.edits = scratch;
        fail_set.count = 99u;
        if (duckvep_variant_cds_edit_set_build(&s->tx, &s->ex, &s->seq, &s->v,
                                               0u, 0u, s->strand, scratch,
                                               edit_set.count - 1u,
                                               &fail_set) != DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (fail_set.edits != NULL || fail_set.count != 0u) return THEFT_TRIAL_FAIL;
        g_cds_edit_set_mnv_cov.capfail++;
    }
    if (duckvep_haplotype_apply_cds_edits(s->cds, s->cds_lenv, edit_set.edits,
                                          edit_set.count, s->strand, mutated,
                                          sizeof mutated, &mutated_len,
                                          &result) != DUCKVEP_HAPLOTYPE_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (mutated_len != (size_t)s->expect_len ||
        memcmp(mutated, s->expect_cds, (size_t)s->expect_len) != 0) {
        return THEFT_TRIAL_FAIL;
    }
    if (s->strand > 0) g_cds_edit_set_mnv_cov.fwd++; else g_cds_edit_set_mnv_cov.rev++;
    if (s->expect_region == KPROP_CDS_EDIT_START) g_cds_edit_set_mnv_cov.start++;
    else if (s->expect_region == KPROP_CDS_EDIT_BODY) g_cds_edit_set_mnv_cov.body++;
    else if (s->expect_region == KPROP_CDS_EDIT_STOP) g_cds_edit_set_mnv_cov.stop++;
    else return THEFT_TRIAL_FAIL;
    g_cds_edit_set_mnv_cov.multi++;
    return THEFT_TRIAL_PASS;
}

TEST cds_edit_set_builder_splits_mnv_diff_islands(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "variant CDS edit-set builder splits MNV diff islands";
    cfg.prop1 = prop_cds_edit_set_splits_mnv_islands;
    cfg.type_info[0] = &kprop_cds_edit_set_mnv_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_cds_edit_set_mnv_cov, 0, sizeof g_cds_edit_set_mnv_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_cds_edit_set_mnv_cov.fwd > 0u);
    ASSERT(g_cds_edit_set_mnv_cov.rev > 0u);
    ASSERT(g_cds_edit_set_mnv_cov.start > 0u);
    ASSERT(g_cds_edit_set_mnv_cov.body > 0u);
    ASSERT(g_cds_edit_set_mnv_cov.stop > 0u);
    ASSERT(g_cds_edit_set_mnv_cov.multi > 0u);
    ASSERT(g_cds_edit_set_mnv_cov.capfail > 0u);
    fprintf(stderr,
            "[cds-edit-set-mnv coverage] fwd=%u rev=%u start=%u body=%u stop=%u multi=%u capfail=%u\n",
            g_cds_edit_set_mnv_cov.fwd, g_cds_edit_set_mnv_cov.rev,
            g_cds_edit_set_mnv_cov.start, g_cds_edit_set_mnv_cov.body,
            g_cds_edit_set_mnv_cov.stop, g_cds_edit_set_mnv_cov.multi,
            g_cds_edit_set_mnv_cov.capfail);
    PASS();
}

static struct {
    uint32_t snv;
    uint32_t mnv;
    uint32_t ins;
    uint32_t del;
    uint32_t indel;
    uint32_t fwd;
    uint32_t rev;
    uint32_t pep_same;
    uint32_t pep_diff;
    uint32_t capfail;
} g_coding_context_cov;

static enum theft_trial_res prop_coding_context_matches_direct_oracles(struct theft *t,
                                                                       void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t scratch[8];
    duckvep_edit_set_t edit_set;
    duckvep_coding_context_t ctx;
    uint8_t alt_cds[80];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    uint8_t ref_want[32];
    uint8_t alt_want[32];
    size_t ref_want_len = 0u;
    size_t alt_want_len = 0u;
    int cds_changed;
    uint32_t rf, rl, af, al;
    (void)t;

    if (duckvep_variant_cds_edit_set_build(&s->tx, &s->ex, &s->seq, &s->v,
                                           0u, 0u, s->strand, scratch, 8u,
                                           &edit_set) != DUCKVEP_CDS_EDIT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    memset(ref_pep, 0xA5, sizeof ref_pep);
    memset(alt_pep, 0x5A, sizeof alt_pep);
    if (duckvep_coding_context_build(s->cds, s->cds_lenv, &edit_set, s->strand,
                                     (duckvep_codon_table_t)s->ctab,
                                     alt_cds, sizeof alt_cds,
                                     ref_pep, sizeof ref_pep,
                                     alt_pep, sizeof alt_pep,
                                     &ctx) != DUCKVEP_CODING_CONTEXT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.alt_cds != alt_cds || ctx.ref_cds != s->cds) return THEFT_TRIAL_FAIL;
    if (ctx.alt_cds_len != (size_t)s->expect_len ||
        memcmp(ctx.alt_cds, s->expect_cds, (size_t)s->expect_len) != 0) {
        return THEFT_TRIAL_FAIL;
    }
    if (!kprop_translate_full_oracle(s->cds, s->cds_lenv,
                                     (duckvep_codon_table_t)s->ctab,
                                     ref_want, &ref_want_len)) {
        return THEFT_TRIAL_FAIL;
    }
    if (!kprop_translate_full_oracle(s->expect_cds, (size_t)s->expect_len,
                                     (duckvep_codon_table_t)s->ctab,
                                     alt_want, &alt_want_len)) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.ref_peptide_len != ref_want_len || ctx.alt_peptide_len != alt_want_len) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.ref_peptide[ctx.ref_peptide_len] != (uint8_t)'\0' ||
        ctx.alt_peptide[ctx.alt_peptide_len] != (uint8_t)'\0') {
        return THEFT_TRIAL_FAIL;
    }
    if (memcmp(ctx.ref_peptide, ref_want, ref_want_len) != 0 ||
        memcmp(ctx.alt_peptide, alt_want, alt_want_len) != 0) {
        return THEFT_TRIAL_FAIL;
    }
    cds_changed = kprop_cds_changed_oracle(s->cds, s->cds_lenv,
                                           s->expect_cds, (size_t)s->expect_len);
    if (cds_changed < 0 || ctx.cds_changed != (uint8_t)cds_changed) return THEFT_TRIAL_FAIL;
    if (ctx.applied_edits != edit_set.count) return THEFT_TRIAL_FAIL;
    if (ctx.length_diff != (int64_t)s->expect_len - (int64_t)s->cds_lenv) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.flags != kprop_context_flags_oracle(&edit_set)) return THEFT_TRIAL_FAIL;
    kprop_peptide_window_oracle(ref_want, ref_want_len, alt_want, alt_want_len,
                                &rf, &rl, &af, &al);
    if (ctx.ref_first_changed_codon != rf || ctx.ref_last_changed_codon != rl ||
        ctx.alt_first_changed_codon != af || ctx.alt_last_changed_codon != al) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.alt_cds_len > 0u) {
        duckvep_coding_context_t fail_ctx;
        if (duckvep_coding_context_build(s->cds, s->cds_lenv, &edit_set, s->strand,
                                         (duckvep_codon_table_t)s->ctab,
                                         alt_cds, ctx.alt_cds_len - 1u,
                                         ref_pep, sizeof ref_pep,
                                         alt_pep, sizeof alt_pep,
                                         &fail_ctx) !=
            DUCKVEP_CODING_CONTEXT_ALT_CDS_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (!kprop_context_zero(&fail_ctx)) return THEFT_TRIAL_FAIL;
        g_coding_context_cov.capfail++;
    }
    if (ctx.ref_peptide_len > 0u) {
        duckvep_coding_context_t fail_ctx;
        if (duckvep_coding_context_build(s->cds, s->cds_lenv, &edit_set, s->strand,
                                         (duckvep_codon_table_t)s->ctab,
                                         alt_cds, sizeof alt_cds,
                                         ref_pep, ctx.ref_peptide_len,
                                         alt_pep, sizeof alt_pep,
                                         &fail_ctx) !=
            DUCKVEP_CODING_CONTEXT_REF_PEPTIDE_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (!kprop_context_zero(&fail_ctx)) return THEFT_TRIAL_FAIL;
        g_coding_context_cov.capfail++;
    }
    if (ctx.alt_peptide_len > 0u) {
        duckvep_coding_context_t fail_ctx;
        if (duckvep_coding_context_build(s->cds, s->cds_lenv, &edit_set, s->strand,
                                         (duckvep_codon_table_t)s->ctab,
                                         alt_cds, sizeof alt_cds,
                                         ref_pep, sizeof ref_pep,
                                         alt_pep, ctx.alt_peptide_len,
                                         &fail_ctx) !=
            DUCKVEP_CODING_CONTEXT_ALT_PEPTIDE_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (!kprop_context_zero(&fail_ctx)) return THEFT_TRIAL_FAIL;
        g_coding_context_cov.capfail++;
    }

    if (s->expect_shape == KPROP_CDS_EDIT_SNV) g_coding_context_cov.snv++;
    else if (s->expect_shape == KPROP_CDS_EDIT_MNV) g_coding_context_cov.mnv++;
    else if (s->expect_shape == KPROP_CDS_EDIT_INS) g_coding_context_cov.ins++;
    else if (s->expect_shape == KPROP_CDS_EDIT_DEL) g_coding_context_cov.del++;
    else if (s->expect_shape == KPROP_CDS_EDIT_INDEL) g_coding_context_cov.indel++;
    else return THEFT_TRIAL_FAIL;
    if (s->strand > 0) g_coding_context_cov.fwd++; else g_coding_context_cov.rev++;
    if (rf == 0u && af == 0u) g_coding_context_cov.pep_same++;
    else g_coding_context_cov.pep_diff++;
    return THEFT_TRIAL_PASS;
}

TEST coding_context_matches_direct_oracles(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "coding context == direct CDS splice + full peptide oracles";
    cfg.prop1 = prop_coding_context_matches_direct_oracles;
    cfg.type_info[0] = &kprop_cds_edit_builder_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_coding_context_cov, 0, sizeof g_coding_context_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_coding_context_cov.snv > 0u);
    ASSERT(g_coding_context_cov.mnv > 0u);
    ASSERT(g_coding_context_cov.ins > 0u);
    ASSERT(g_coding_context_cov.del > 0u);
    ASSERT(g_coding_context_cov.indel > 0u);
    ASSERT(g_coding_context_cov.fwd > 0u);
    ASSERT(g_coding_context_cov.rev > 0u);
    ASSERT(g_coding_context_cov.pep_diff > 0u);
    ASSERT(g_coding_context_cov.capfail > 0u);
    fprintf(stderr,
            "[coding-context coverage] snv=%u mnv=%u ins=%u del=%u indel=%u fwd=%u rev=%u pep_same=%u pep_diff=%u capfail=%u\n",
            g_coding_context_cov.snv, g_coding_context_cov.mnv,
            g_coding_context_cov.ins, g_coding_context_cov.del,
            g_coding_context_cov.indel, g_coding_context_cov.fwd,
            g_coding_context_cov.rev, g_coding_context_cov.pep_same,
            g_coding_context_cov.pep_diff, g_coding_context_cov.capfail);
    PASS();
}

static struct {
    uint32_t snv;
    uint32_t mnv;
    uint32_t ins;
    uint32_t del;
    uint32_t indel;
    uint32_t fwd;
    uint32_t rev;
    uint32_t pep_same;
    uint32_t pep_diff;
    uint32_t capfail;
} g_variant_coding_context_cov;

static enum theft_trial_res prop_variant_coding_context_matches_oracles(struct theft *t,
                                                                        void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edit_scratch[8];
    duckvep_haplotype_edit_t expected_edits[8];
    duckvep_edit_set_t expected_set;
    duckvep_coding_context_t ctx;
    uint8_t alt_cds[80];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    uint8_t ref_want[32];
    uint8_t alt_want[32];
    size_t ref_want_len = 0u;
    size_t alt_want_len = 0u;
    int cds_changed;
    uint32_t rf, rl, af, al;
    (void)t;

    if (duckvep_variant_cds_edit_set_build(&s->tx, &s->ex, &s->seq, &s->v,
                                           0u, 0u, s->strand, expected_edits, 8u,
                                           &expected_set) != DUCKVEP_CDS_EDIT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (duckvep_variant_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->strand, edit_scratch, 8u,
                                             alt_cds, sizeof alt_cds,
                                             ref_pep, sizeof ref_pep,
                                             alt_pep, sizeof alt_pep,
                                             &ctx) != DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.ref_cds != s->cds || ctx.alt_cds != alt_cds) return THEFT_TRIAL_FAIL;
    if (ctx.alt_cds_len != (size_t)s->expect_len ||
        memcmp(ctx.alt_cds, s->expect_cds, (size_t)s->expect_len) != 0) {
        return THEFT_TRIAL_FAIL;
    }
    if (!kprop_translate_full_oracle(s->cds, s->cds_lenv,
                                     (duckvep_codon_table_t)s->ctab,
                                     ref_want, &ref_want_len)) {
        return THEFT_TRIAL_FAIL;
    }
    if (!kprop_translate_full_oracle(s->expect_cds, (size_t)s->expect_len,
                                     (duckvep_codon_table_t)s->ctab,
                                     alt_want, &alt_want_len)) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.ref_peptide_len != ref_want_len || ctx.alt_peptide_len != alt_want_len) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.ref_peptide[ctx.ref_peptide_len] != (uint8_t)'\0' ||
        ctx.alt_peptide[ctx.alt_peptide_len] != (uint8_t)'\0') {
        return THEFT_TRIAL_FAIL;
    }
    if (memcmp(ctx.ref_peptide, ref_want, ref_want_len) != 0 ||
        memcmp(ctx.alt_peptide, alt_want, alt_want_len) != 0) {
        return THEFT_TRIAL_FAIL;
    }
    cds_changed = kprop_cds_changed_oracle(s->cds, s->cds_lenv,
                                           s->expect_cds, (size_t)s->expect_len);
    if (cds_changed < 0 || ctx.cds_changed != (uint8_t)cds_changed) return THEFT_TRIAL_FAIL;
    if (ctx.applied_edits != expected_set.count) return THEFT_TRIAL_FAIL;
    if (ctx.length_diff != (int64_t)s->expect_len - (int64_t)s->cds_lenv) {
        return THEFT_TRIAL_FAIL;
    }
    if (ctx.flags != kprop_context_flags_oracle(&expected_set)) return THEFT_TRIAL_FAIL;
    if (ctx.flags != kprop_single_variant_flags_oracle(s->cds_lenv, s->expect_len)) {
        return THEFT_TRIAL_FAIL;
    }
    kprop_peptide_window_oracle(ref_want, ref_want_len, alt_want, alt_want_len,
                                &rf, &rl, &af, &al);
    if (ctx.ref_first_changed_codon != rf || ctx.ref_last_changed_codon != rl ||
        ctx.alt_first_changed_codon != af || ctx.alt_last_changed_codon != al) {
        return THEFT_TRIAL_FAIL;
    }

    {
        duckvep_coding_context_t fail_ctx;
        if (duckvep_variant_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                                 0u, 0u, s->strand, NULL, 0u,
                                                 alt_cds, sizeof alt_cds,
                                                 ref_pep, sizeof ref_pep,
                                                 alt_pep, sizeof alt_pep,
                                                 &fail_ctx) !=
            DUCKVEP_VARIANT_CODING_CONTEXT_EDIT_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (!kprop_context_zero(&fail_ctx)) return THEFT_TRIAL_FAIL;
        g_variant_coding_context_cov.capfail++;
    }
    if (ctx.alt_cds_len > 0u) {
        duckvep_coding_context_t fail_ctx;
        if (duckvep_variant_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                                 0u, 0u, s->strand, edit_scratch, 8u,
                                                 alt_cds, ctx.alt_cds_len - 1u,
                                                 ref_pep, sizeof ref_pep,
                                                 alt_pep, sizeof alt_pep,
                                                 &fail_ctx) !=
            DUCKVEP_VARIANT_CODING_CONTEXT_ALT_CDS_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (!kprop_context_zero(&fail_ctx)) return THEFT_TRIAL_FAIL;
        g_variant_coding_context_cov.capfail++;
    }
    if (ctx.ref_peptide_len > 0u) {
        duckvep_coding_context_t fail_ctx;
        if (duckvep_variant_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                                 0u, 0u, s->strand, edit_scratch, 8u,
                                                 alt_cds, sizeof alt_cds,
                                                 ref_pep, ctx.ref_peptide_len,
                                                 alt_pep, sizeof alt_pep,
                                                 &fail_ctx) !=
            DUCKVEP_VARIANT_CODING_CONTEXT_REF_PEPTIDE_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (!kprop_context_zero(&fail_ctx)) return THEFT_TRIAL_FAIL;
        g_variant_coding_context_cov.capfail++;
    }
    if (ctx.alt_peptide_len > 0u) {
        duckvep_coding_context_t fail_ctx;
        if (duckvep_variant_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                                 0u, 0u, s->strand, edit_scratch, 8u,
                                                 alt_cds, sizeof alt_cds,
                                                 ref_pep, sizeof ref_pep,
                                                 alt_pep, ctx.alt_peptide_len,
                                                 &fail_ctx) !=
            DUCKVEP_VARIANT_CODING_CONTEXT_ALT_PEPTIDE_BUFFER_TOO_SMALL) {
            return THEFT_TRIAL_FAIL;
        }
        if (!kprop_context_zero(&fail_ctx)) return THEFT_TRIAL_FAIL;
        g_variant_coding_context_cov.capfail++;
    }

    if (s->expect_shape == KPROP_CDS_EDIT_SNV) g_variant_coding_context_cov.snv++;
    else if (s->expect_shape == KPROP_CDS_EDIT_MNV) g_variant_coding_context_cov.mnv++;
    else if (s->expect_shape == KPROP_CDS_EDIT_INS) g_variant_coding_context_cov.ins++;
    else if (s->expect_shape == KPROP_CDS_EDIT_DEL) g_variant_coding_context_cov.del++;
    else if (s->expect_shape == KPROP_CDS_EDIT_INDEL) g_variant_coding_context_cov.indel++;
    else return THEFT_TRIAL_FAIL;
    if (s->strand > 0) g_variant_coding_context_cov.fwd++;
    else g_variant_coding_context_cov.rev++;
    if (rf == 0u && af == 0u) g_variant_coding_context_cov.pep_same++;
    else g_variant_coding_context_cov.pep_diff++;
    return THEFT_TRIAL_PASS;
}

TEST variant_coding_context_matches_oracles(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "variant coding context == direct CDS splice + full peptide oracles";
    cfg.prop1 = prop_variant_coding_context_matches_oracles;
    cfg.type_info[0] = &kprop_cds_edit_builder_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_variant_coding_context_cov, 0, sizeof g_variant_coding_context_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_variant_coding_context_cov.snv > 0u);
    ASSERT(g_variant_coding_context_cov.mnv > 0u);
    ASSERT(g_variant_coding_context_cov.ins > 0u);
    ASSERT(g_variant_coding_context_cov.del > 0u);
    ASSERT(g_variant_coding_context_cov.indel > 0u);
    ASSERT(g_variant_coding_context_cov.fwd > 0u);
    ASSERT(g_variant_coding_context_cov.rev > 0u);
    ASSERT(g_variant_coding_context_cov.pep_diff > 0u);
    ASSERT(g_variant_coding_context_cov.capfail > 0u);
    fprintf(stderr,
            "[variant-coding-context coverage] snv=%u mnv=%u ins=%u del=%u indel=%u fwd=%u rev=%u pep_same=%u pep_diff=%u capfail=%u\n",
            g_variant_coding_context_cov.snv, g_variant_coding_context_cov.mnv,
            g_variant_coding_context_cov.ins, g_variant_coding_context_cov.del,
            g_variant_coding_context_cov.indel, g_variant_coding_context_cov.fwd,
            g_variant_coding_context_cov.rev, g_variant_coding_context_cov.pep_same,
            g_variant_coding_context_cov.pep_diff,
            g_variant_coding_context_cov.capfail);
    PASS();
}

static struct {
    uint32_t syn;
    uint32_t mis;
    uint32_t stop_gained;
    uint32_t stop_lost;
    uint32_t stop_retained;
    uint32_t fwd;
    uint32_t rev;
} g_context_delta_cov;

static enum theft_trial_res prop_context_delta_matches_codon_oracle(struct theft *t,
                                                                    void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edit_scratch[4];
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;
    uint8_t alt_cds[64];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    size_t first_diff = 0u;
    size_t codon_start;
    size_t codon_idx;
    char ref_codon[4];
    char alt_codon[4];
    duckvep_codon_result_t cr;
    uint32_t j;
    (void)t;

    if (duckvep_variant_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->strand,
                                             edit_scratch, 4u,
                                             alt_cds, sizeof alt_cds,
                                             ref_pep, sizeof ref_pep,
                                             alt_pep, sizeof alt_pep,
                                             &ctx) != DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (duckvep_coding_context_delta_fill(&ctx, s->flags, &delta) !=
        DUCKVEP_CONTEXT_DELTA_OK) {
        return THEFT_TRIAL_FAIL;
    }
    while (first_diff < s->cds_lenv && s->cds[first_diff] == s->expect_cds[first_diff]) {
        first_diff++;
    }
    if (first_diff >= s->cds_lenv) return THEFT_TRIAL_FAIL;
    codon_start = first_diff - (first_diff % 3u);
    codon_idx = codon_start / 3u;
    for (j = 0u; j < 3u; j++) {
        ref_codon[j] = (char)s->cds[codon_start + (size_t)j];
        alt_codon[j] = (char)s->expect_cds[codon_start + (size_t)j];
    }
    ref_codon[3] = '\0';
    alt_codon[3] = '\0';
    cr = duckvep_codon_change(ref_codon, alt_codon, DUCKVEP_CODON_TABLE_STANDARD);
    if (cr.change & DUCKVEP_CODON_INVALID) return THEFT_TRIAL_FAIL;
    if (!delta.valid || delta.cdna_pos != -1 || delta.cds_pos != -1 ||
        delta.protein_pos != (int32_t)codon_idx + 1 ||
        delta.ref_aa != (uint8_t)cr.aa_ref || delta.alt_aa != (uint8_t)cr.aa_alt) {
        return THEFT_TRIAL_FAIL;
    }
    if (s->expect_region == KPROP_CONTEXT_DELTA_SYNONYMOUS) {
        if (!delta.synonymous || delta.missense || delta.stop_gained ||
            delta.stop_lost || delta.stop_retained) return THEFT_TRIAL_FAIL;
        g_context_delta_cov.syn++;
    } else if (s->expect_region == KPROP_CONTEXT_DELTA_MISSENSE) {
        if (!delta.missense || delta.synonymous || delta.stop_gained ||
            delta.stop_lost || delta.stop_retained) return THEFT_TRIAL_FAIL;
        g_context_delta_cov.mis++;
    } else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_GAINED) {
        if (!delta.stop_gained || delta.synonymous || delta.missense ||
            delta.stop_lost || delta.stop_retained) return THEFT_TRIAL_FAIL;
        g_context_delta_cov.stop_gained++;
    } else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_LOST) {
        if (!delta.stop_lost || delta.synonymous || delta.missense ||
            delta.stop_gained || delta.stop_retained) return THEFT_TRIAL_FAIL;
        g_context_delta_cov.stop_lost++;
    } else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_RETAINED) {
        if (!delta.stop_retained || delta.synonymous || delta.missense ||
            delta.stop_gained || delta.stop_lost) return THEFT_TRIAL_FAIL;
        g_context_delta_cov.stop_retained++;
    } else return THEFT_TRIAL_FAIL;
    if (s->strand > 0) g_context_delta_cov.fwd++;
    else g_context_delta_cov.rev++;
    return THEFT_TRIAL_PASS;
}

TEST coding_context_delta_matches_codon_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "coding context delta == single-codon oracle";
    cfg.prop1 = prop_context_delta_matches_codon_oracle;
    cfg.type_info[0] = &kprop_context_delta_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_context_delta_cov, 0, sizeof g_context_delta_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_context_delta_cov.syn > 0u);
    ASSERT(g_context_delta_cov.mis > 0u);
    ASSERT(g_context_delta_cov.stop_gained > 0u);
    ASSERT(g_context_delta_cov.stop_lost > 0u);
    ASSERT(g_context_delta_cov.stop_retained > 0u);
    ASSERT(g_context_delta_cov.fwd > 0u);
    ASSERT(g_context_delta_cov.rev > 0u);
    fprintf(stderr,
            "[context-delta coverage] syn=%u mis=%u stop_gained=%u stop_lost=%u stop_retained=%u fwd=%u rev=%u\n",
            g_context_delta_cov.syn, g_context_delta_cov.mis,
            g_context_delta_cov.stop_gained, g_context_delta_cov.stop_lost,
            g_context_delta_cov.stop_retained, g_context_delta_cov.fwd,
            g_context_delta_cov.rev);
    PASS();
}

static int kprop_delta_is_frameshift_at(const duckvep_sequence_delta_t *d,
                                        int32_t protein_pos) {
    return d != NULL && d->valid && d->frameshift &&
           !d->synonymous && !d->missense && !d->stop_gained && !d->stop_lost &&
           !d->stop_retained && !d->start_lost && !d->start_retained &&
           !d->inframe_deletion && !d->inframe_insertion && !d->protein_altering &&
           !d->coding_unknown &&
           d->cdna_pos == -1 && d->cds_pos == -1 && d->protein_pos == protein_pos &&
           d->ref_aa == (uint8_t)0u && d->alt_aa == (uint8_t)0u;
}

/* Independent oracle for VEP's LOCAL-window frameshift stop_gained. Re-derives the
 * codon window from the raw CDS bytes (substr(alt_cds, codon_cds_start-1,
 * codon_len + net_delta)) and tests whether a WHOLE codon there is a standard stop the
 * reference window lacked. Deliberately re-translates from bytes (stops enumerated by
 * hand for the standard table) rather than reading the kernel's pre-translated peptides,
 * so it does not restate the implementation. Only valid for the standard codon table. */
static int kprop_is_standard_stop(const uint8_t *b) {
    return b[0] == (uint8_t)'T' &&
           ((b[1] == (uint8_t)'A' && (b[2] == (uint8_t)'A' || b[2] == (uint8_t)'G')) ||
            (b[1] == (uint8_t)'G' && b[2] == (uint8_t)'A'));
}

static int kprop_frameshift_local_stop_oracle(
    const uint8_t *ref_cds, size_t ref_cds_len,
    const uint8_t *alt_cds, size_t alt_cds_len,
    uint32_t cds_start, uint32_t ref_len, int64_t length_diff) {

    uint64_t first = cds_start;
    uint64_t last, tv_s, tv_e, codon_cds_start, codon_len, off, i, whole;
    int64_t win;
    int alt_stop = 0;
    int ref_stop = 0;

    if (first == 0u) return 0;
    tv_s = ((first - 1u) / 3u) + 1u;
    codon_cds_start = tv_s * 3u - 2u;
    off = codon_cds_start - 1u;
    /* Pure insertion before the first CDS base spans no reference codon; never form
     * last == 0 (whose (last - 1) would underflow), matching the kernel helper. */
    if (ref_len > 0u) {
        last = first + ref_len - 1u;
        tv_e = ((last - 1u) / 3u) + 1u;
        codon_len = (tv_e >= tv_s) ? (tv_e - tv_s + 1u) * 3u : 0u;
    } else if (first > 1u) {
        last = first - 1u;
        tv_e = ((last - 1u) / 3u) + 1u;
        codon_len = (tv_e >= tv_s) ? (tv_e - tv_s + 1u) * 3u : 0u;
    } else {
        codon_len = 0u;
    }

    win = (int64_t)codon_len + length_diff;
    if (win > 0 && off <= alt_cds_len) {
        uint64_t avail = (uint64_t)alt_cds_len - off;
        if ((uint64_t)win > avail) win = (int64_t)avail;
        whole = (uint64_t)win / 3u;
        for (i = 0u; i < whole; i++) {
            if (kprop_is_standard_stop(alt_cds + off + i * 3u)) { alt_stop = 1; break; }
        }
    }
    if (codon_len > 0u && off <= ref_cds_len) {
        uint64_t avail = (uint64_t)ref_cds_len - off;
        uint64_t rw = codon_len;
        if (rw > avail) rw = avail;
        whole = rw / 3u;
        for (i = 0u; i < whole; i++) {
            if (kprop_is_standard_stop(ref_cds + off + i * 3u)) { ref_stop = 1; break; }
        }
    }
    return alt_stop && !ref_stop;
}

/* Deterministic anchor for VEP's local-window frameshift stop_gained: the recomputed
 * codon AT THE JUNCTION being a stop composites frameshift_variant&stop_gained, while a
 * frameshift whose junction codon is not a stop stays a bare frameshift EVEN WHEN the
 * shifted downstream frame contains a stop — the discriminator proving the engine uses
 * VEP's local window, not a downstream retranslation. */
TEST coding_context_delta_frameshift_stop_gained_scene(void) {
    duckvep_haplotype_edit_t edit;
    duckvep_edit_set_t edit_set;
    uint8_t alt_cds[40];
    uint8_t ref_pep[20];
    uint8_t alt_pep[20];
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;

    memset(&edit, 0, sizeof edit);
    edit.variant_strand = 1;
    edit_set.edits = &edit; edit_set.count = 1u;

    /* (A) +1 frameshift whose recomputed codon 3 becomes TAG (stop). Replace the 'C' at
     *     CDS pos 7 with "TA": ATG AAA CGT GAA -> ATG AAA TAG TGA A. */
    {
        static uint8_t cdsA[12] = {
            'A','T','G',  'A','A','A',  'C','G','T',  'G','A','A'
        };
        static const uint8_t insTA[2] = { 'T','A' };
        edit.cds_start = 7u; edit.ref = cdsA + 6u; edit.ref_len = 1u;
        edit.alt = insTA; edit.alt_len = 2u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cdsA, sizeof cdsA, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep,
                                               sizeof ref_pep, alt_pep, sizeof alt_pep,
                                               &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.frameshift && delta.stop_gained);
        ASSERT_EQ(3, delta.protein_pos);
        /* oracle agreement */
        ASSERT_EQ(1, kprop_frameshift_local_stop_oracle(cdsA, sizeof cdsA, ctx.alt_cds,
                                                         ctx.alt_cds_len, 7u, 1u,
                                                         ctx.length_diff));
    }

    /* (B) +1 frameshift whose recomputed codon 3 is CAA (not a stop) but the shifted
     *     frame carries TGA one codon downstream. stop_gained MUST be 0. Replace the 'C'
     *     at CDS pos 7 with "CA": ATG AAA CAT GAA -> ATG AAA CAA TGA A. */
    {
        static uint8_t cdsB[12] = {
            'A','T','G',  'A','A','A',  'C','A','T',  'G','A','A'
        };
        static const uint8_t insCA[2] = { 'C','A' };
        int downstream_stop = 0;
        size_t i;
        edit.cds_start = 7u; edit.ref = cdsB + 6u; edit.ref_len = 1u;
        edit.alt = insCA; edit.alt_len = 2u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cdsB, sizeof cdsB, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep,
                                               sizeof ref_pep, alt_pep, sizeof alt_pep,
                                               &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.frameshift && !delta.stop_gained);
        ASSERT_EQ(3, delta.protein_pos);
        /* non-vacuous: a downstream stop really is present in the shifted alt peptide,
         * proving the local window (not full retranslation) drove the 0 verdict. */
        for (i = 3u; i < ctx.alt_peptide_len; i++) {
            if (ctx.alt_peptide[i] == (uint8_t)'*') { downstream_stop = 1; break; }
        }
        ASSERT(downstream_stop);
        ASSERT_EQ(0, kprop_frameshift_local_stop_oracle(cdsB, sizeof cdsB, ctx.alt_cds,
                                                         ctx.alt_cds_len, 7u, 1u,
                                                         ctx.length_diff));
    }

    /* (C) PURE insertion (ref_len 0) — the dominant real ClinVar frameshift&stop_gained
     *     shape (e.g. C>CTTTAA). Insert "TAAG" before CDS pos 7 (codon boundary,
     *     codon_len 0): ATG AAA GGG CCC -> ATG AAA TAA GGG GCC C, recomputed codon 3 = TAA
     *     (stop). Exercises the ref_len==0 window branch and the first_cds>1 path. */
    {
        static uint8_t cdsC[12] = {
            'A','T','G',  'A','A','A',  'G','G','G',  'C','C','C'
        };
        static const uint8_t insTAAG[4] = { 'T','A','A','G' };
        edit.cds_start = 7u; edit.ref = NULL; edit.ref_len = 0u;
        edit.alt = insTAAG; edit.alt_len = 4u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cdsC, sizeof cdsC, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep,
                                               sizeof ref_pep, alt_pep, sizeof alt_pep,
                                               &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.frameshift && delta.stop_gained);
        ASSERT_EQ(3, delta.protein_pos);
        ASSERT_EQ(1, kprop_frameshift_local_stop_oracle(cdsC, sizeof cdsC, ctx.alt_cds,
                                                        ctx.alt_cds_len, 7u, 0u,
                                                        ctx.length_diff));
    }

    /* (D) Codon-table awareness: the recomputed junction codon is TGA, a stop under the
     *     standard table but tryptophan under the vertebrate-mitochondrial table. The
     *     helper reads the peptide translated with the build's table, so stop_gained must
     *     flip. Replace the 'C' at CDS pos 7 with "TG": ATG AAA CAG GAA -> ATG AAA TGA GGA A. */
    {
        static uint8_t cdsD[12] = {
            'A','T','G',  'A','A','A',  'C','A','G',  'G','A','A'
        };
        static const uint8_t insTG[2] = { 'T','G' };
        edit.cds_start = 7u; edit.ref = cdsD + 6u; edit.ref_len = 1u;
        edit.alt = insTG; edit.alt_len = 2u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cdsD, sizeof cdsD, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep,
                                               sizeof ref_pep, alt_pep, sizeof alt_pep,
                                               &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.frameshift && delta.stop_gained);

        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cdsD, sizeof cdsD, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_VERT_MITO,
                                               alt_cds, sizeof alt_cds, ref_pep,
                                               sizeof ref_pep, alt_pep, sizeof alt_pep,
                                               &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.frameshift && !delta.stop_gained);
    }

    /* (E) Regression guard for the pure-insertion first_cds==1 window (formerly a
     *     (last_cds-1) underflow). Insert "TAAC" before CDS pos 1 on a CDS_start_NF
     *     transcript (so the start-codon guard is bypassed): recomputed codon 1 = TAA
     *     (stop). Must resolve cleanly with no out-of-bounds read. */
    {
        static uint8_t cdsE[12] = {
            'A','T','G',  'A','A','A',  'G','G','G',  'C','C','C'
        };
        static const uint8_t insTAAC[4] = { 'T','A','A','C' };
        edit.cds_start = 1u; edit.ref = NULL; edit.ref_len = 0u;
        edit.alt = insTAAC; edit.alt_len = 4u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cdsE, sizeof cdsE, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep,
                                               sizeof ref_pep, alt_pep, sizeof alt_pep,
                                               &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx,
                                                    (uint64_t)DUCKVEP_TX_CDS_START_NF,
                                                    &delta));
        ASSERT(delta.valid && delta.frameshift && delta.stop_gained);
        ASSERT_EQ(1, kprop_frameshift_local_stop_oracle(cdsE, sizeof cdsE, ctx.alt_cds,
                                                        ctx.alt_cds_len, 1u, 0u,
                                                        ctx.length_diff));
    }

    PASS();
}

/* Deterministic anchor for the general CodingContext frameshift fact. A net CDS length
 * change not divisible by three resolves to frameshift; a start-codon edit also carries
 * VEP's start_lost fact unless CDS_START_NF suppresses that predicate. */
TEST coding_context_delta_frameshift_known_scene(void) {
    /* M K P G F *  (start codon intact unless the edit lands in codon 1). */
    static uint8_t cds[18] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','T',  'T','A','A'
    };
    static const uint8_t ins1[2]  = { 'C','A' };  /* replace 1 base with 2 -> +1 */
    static const uint8_t del2a[1] = { 'C' };      /* replace 2 bases with 1 -> -1 */
    duckvep_haplotype_edit_t edit;
    duckvep_edit_set_t edit_set;
    uint8_t alt_cds[40];
    uint8_t ref_pep[20];
    uint8_t alt_pep[20];
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;

    memset(&edit, 0, sizeof edit);
    edit.variant_strand = 1;
    edit_set.edits = &edit; edit_set.count = 1u;

    /* (1) +1 frameshift at codon 3 (body), ATG intact -> frameshift at protein pos 3. */
    edit.cds_start = 7u; edit.ref = cds + 6u; edit.ref_len = 1u;
    edit.alt = ins1; edit.alt_len = 2u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(cds, sizeof cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(kprop_delta_is_frameshift_at(&delta, 3));

    /* (2) -1 frameshift at codon 4 (body) -> frameshift at protein pos 4. */
    edit.cds_start = 10u; edit.ref = cds + 9u; edit.ref_len = 2u;
    edit.alt = del2a; edit.alt_len = 1u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(cds, sizeof cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(kprop_delta_is_frameshift_at(&delta, 4));

    /* (3) +1 frameshift disrupting the start codon (cds 1). */
    edit.cds_start = 1u; edit.ref = cds; edit.ref_len = 1u;
    edit.alt = ins1; edit.alt_len = 2u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(cds, sizeof cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(delta.valid && delta.frameshift && delta.start_lost &&
           !delta.start_retained);

    /* (4) same start-disrupting edit but CDS_start_NF set: no start_lost to add, so the
     *     frameshift fact is emitted. */
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx,
                                                (uint64_t)DUCKVEP_TX_CDS_START_NF, &delta));
    ASSERT(delta.valid && delta.frameshift && !delta.start_lost &&
           !delta.start_retained);

    PASS();
}

/* Randomized frame-arithmetic oracle: over random valid CDS + single body edits, the
 * general CodingContext emits the frameshift fact exactly when the net length change is
 * not divisible by three, and never flags an in-frame length change as a frameshift. The
 * oracle is the frame arithmetic itself, independent of the classifier implementation. */
TEST coding_context_delta_frameshift_matches_length_oracle(void) {
    static const char bases[4] = { 'A', 'C', 'G', 'T' };
    uint64_t rng = kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED) ^ UINT64_C(0xF00DFACE);
    uint64_t trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    uint64_t t;
    unsigned fs_seen = 0u;
    unsigned inframe_seen = 0u;
    unsigned fs_stop_seen = 0u;

    for (t = 0u; t < trials; t++) {
        uint8_t cds[36];
        uint8_t alt_bytes[4];
        uint8_t alt_cds[64];
        uint8_t ref_pep[32];
        uint8_t alt_pep[32];
        duckvep_haplotype_edit_t edit;
        duckvep_edit_set_t edit_set;
        duckvep_coding_context_t ctx;
        duckvep_sequence_delta_t delta;
        size_t ncodon, cds_len, i;
        uint32_t codon_idx, cds_start;
        unsigned ref_len, alt_len;
        int net;

        rng = rng * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
        ncodon = 7u + (size_t)((rng >> 33) % 6u); /* 7..12 codons */
        cds_len = ncodon * 3u;
        cds[0] = 'A'; cds[1] = 'T'; cds[2] = 'G';
        for (i = 3u; i < cds_len; i++) {
            rng = rng * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
            cds[i] = (uint8_t)bases[(rng >> 33) % 4u];
        }
        /* Edit in a body codon: not codon 1, and leaving >= 2 codons of tail. */
        rng = rng * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
        codon_idx = 1u + (uint32_t)((rng >> 33) % (ncodon - 3u)); /* 1..ncodon-3 */
        cds_start = codon_idx * 3u + 1u;
        rng = rng * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
        ref_len = 1u + (unsigned)((rng >> 33) % 4u); /* 1..4 */
        rng = rng * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
        alt_len = 1u + (unsigned)((rng >> 33) % 4u); /* 1..4 */
        for (i = 0u; i < alt_len; i++) {
            rng = rng * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
            alt_bytes[i] = (uint8_t)bases[(rng >> 33) % 4u];
        }

        memset(&edit, 0, sizeof edit);
        edit.variant_strand = 1; edit.cds_start = cds_start;
        edit.ref = cds + (cds_start - 1u); edit.ref_len = (uint16_t)ref_len;
        edit.alt = alt_bytes; edit.alt_len = (uint16_t)alt_len;
        edit_set.edits = &edit; edit_set.count = 1u;

        if (duckvep_coding_context_build(cds, cds_len, &edit_set, 1,
                                         DUCKVEP_CODON_TABLE_STANDARD,
                                         alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                                         alt_pep, sizeof alt_pep, &ctx) !=
            DUCKVEP_CODING_CONTEXT_OK) {
            continue;
        }
        (void)duckvep_coding_context_delta_fill(&ctx, 0u, &delta);
        net = (int)alt_len - (int)ref_len;
        if (net != 0 && (net % 3) != 0) {
            int stop_oracle;
            ASSERT(delta.valid && delta.frameshift);
            ASSERT(!delta.inframe_insertion && !delta.inframe_deletion &&
                   !delta.missense && !delta.synonymous && !delta.start_lost);
            /* stop_gained must equal VEP's independent local-window oracle, never the
             * downstream retranslation the full alt peptide would suggest. */
            stop_oracle = kprop_frameshift_local_stop_oracle(cds, cds_len, ctx.alt_cds,
                                                             ctx.alt_cds_len, cds_start,
                                                             ref_len, ctx.length_diff);
            ASSERT_EQ(stop_oracle, delta.stop_gained ? 1 : 0);
            if (delta.stop_gained) fs_stop_seen++;
            fs_seen++;
        } else if (net != 0) {
            if (delta.valid) ASSERT(!delta.frameshift);
            inframe_seen++;
        }
    }
    ASSERT(fs_seen > 0u);
    ASSERT(inframe_seen > 0u);
    fprintf(stderr,
            "[frameshift length-oracle coverage] frameshift=%u inframe_len=%u stop_gained=%u\n",
            fs_seen, inframe_seen, fs_stop_seen);
    PASS();
}

static int kprop_sequence_delta_equal(const duckvep_sequence_delta_t *a,
                                      const duckvep_sequence_delta_t *b) {
    return a != NULL && b != NULL &&
           a->cdna_pos == b->cdna_pos &&
           a->cds_pos == b->cds_pos &&
           a->protein_pos == b->protein_pos &&
           a->ref_aa == b->ref_aa &&
           a->alt_aa == b->alt_aa &&
           a->synonymous == b->synonymous &&
           a->missense == b->missense &&
           a->stop_gained == b->stop_gained &&
           a->stop_lost == b->stop_lost &&
           a->stop_retained == b->stop_retained &&
           a->start_lost == b->start_lost &&
           a->start_retained == b->start_retained &&
           a->frameshift == b->frameshift &&
           a->inframe_deletion == b->inframe_deletion &&
           a->inframe_insertion == b->inframe_insertion &&
           a->protein_altering == b->protein_altering &&
           a->coding_unknown == b->coding_unknown &&
           a->valid == b->valid;
}

static int kprop_delta_is_coarse_cross_codon_missense(const duckvep_sequence_delta_t *d) {
    return d != NULL && d->valid && d->missense && !d->synonymous &&
           !d->stop_gained && !d->stop_lost && !d->stop_retained &&
           !d->start_lost && !d->start_retained && !d->frameshift &&
           !d->inframe_deletion &&
           !d->inframe_insertion && !d->protein_altering && !d->coding_unknown &&
           d->cdna_pos == -1 && d->cds_pos == -1 && d->protein_pos == -1 &&
           d->ref_aa == (uint8_t)0u && d->alt_aa == (uint8_t)0u;
}

static int kprop_delta_is_inframe_deletion_at(const duckvep_sequence_delta_t *d,
                                               int32_t protein_pos) {
    return d != NULL && d->valid && d->inframe_deletion &&
           !d->synonymous && !d->missense && !d->stop_gained && !d->stop_lost &&
           !d->stop_retained && !d->start_lost && !d->start_retained &&
           !d->frameshift &&
           !d->inframe_insertion && !d->protein_altering && !d->coding_unknown &&
           d->cdna_pos == -1 && d->cds_pos == -1 && d->protein_pos == protein_pos &&
           d->ref_aa == (uint8_t)0u && d->alt_aa == (uint8_t)0u;
}

static int kprop_delta_is_inframe_insertion_at(const duckvep_sequence_delta_t *d,
                                                int32_t protein_pos) {
    return d != NULL && d->valid && d->inframe_insertion &&
           !d->synonymous && !d->missense && !d->stop_gained && !d->stop_lost &&
           !d->stop_retained && !d->start_lost && !d->start_retained &&
           !d->frameshift &&
           !d->inframe_deletion && !d->protein_altering && !d->coding_unknown &&
           d->cdna_pos == -1 && d->cds_pos == -1 && d->protein_pos == protein_pos &&
           d->ref_aa == (uint8_t)0u && d->alt_aa == (uint8_t)0u;
}

static int kprop_delta_is_protein_altering_at(const duckvep_sequence_delta_t *d,
                                               int32_t protein_pos) {
    return d != NULL && d->valid && d->protein_altering &&
           !d->synonymous && !d->missense && !d->stop_gained && !d->stop_lost &&
           !d->stop_retained && !d->start_lost && !d->start_retained &&
           !d->frameshift && !d->inframe_deletion && !d->inframe_insertion &&
           !d->coding_unknown &&
           d->cdna_pos == -1 && d->cds_pos == -1 && d->protein_pos == protein_pos &&
           d->ref_aa == (uint8_t)0u && d->alt_aa == (uint8_t)0u;
}

TEST coding_context_delta_inframe_deletion_known_scene(void) {
    static const uint8_t short_cds[9] = {
        'A','T','G',  'A','A','A',  'T','T','T'
    };
    static const uint8_t mixed_short_cds[9] = {
        'a','U','g',  'a','a','a',  'u','U','U'
    };
    static const uint8_t repeat_cds[15] = {
        'A','T','G',  'A','A','A',  'A','A','A',  'C','C','C',  'T','T','T'
    };
    static const uint8_t ambiguous_cds[9] = {
        'A','T','G',  'N','N','N',  'T','T','T'
    };
    static const uint8_t short_alt_cds[6] = {
        'A','T','G',  'T','T','T'
    };
    static const uint8_t fake_ref_pep[4] = { 'M', 'K', 'F', '\0' };
    static const uint8_t fake_alt_pep[3] = { 'M', 'F', '\0' };
    static const uint8_t x_ref_pep[4] = { 'M', 'X', 'F', '\0' };
    static const uint8_t multi_cds[18] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','A',  'T','T','T'
    };
    duckvep_haplotype_edit_t edit;
    duckvep_haplotype_edit_t edits[2];
    duckvep_edit_set_t edit_set;
    uint8_t alt_cds[32];
    uint8_t ref_pep[16];
    uint8_t alt_pep[16];
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;

    memset(&edit, 0, sizeof edit);
    edit.cds_start = 4u; edit.ref_len = 3u; edit.ref = short_cds + 3u;
    edit.alt_len = 0u; edit.alt = NULL; edit.variant_strand = 1;
    edit_set.edits = &edit; edit_set.count = 1u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(short_cds, sizeof short_cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds,
                                           ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(kprop_delta_is_inframe_deletion_at(&delta, 2));

    memset(&edit, 0, sizeof edit);
    edit.cds_start = 4u; edit.ref_len = 3u; edit.ref = mixed_short_cds + 3u;
    edit.alt_len = 0u; edit.alt = NULL; edit.variant_strand = 1;
    edit_set.edits = &edit; edit_set.count = 1u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(mixed_short_cds, sizeof mixed_short_cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds,
                                           ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(kprop_delta_is_inframe_deletion_at(&delta, 2));

    memset(&edit, 0, sizeof edit);
    edit.cds_start = 4u; edit.ref_len = 3u; edit.ref = repeat_cds + 3u;
    edit.alt_len = 0u; edit.alt = NULL; edit.variant_strand = 1;
    edit_set.edits = &edit; edit_set.count = 1u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(repeat_cds, sizeof repeat_cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds,
                                           ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT(ctx.has_single_edit);
    ASSERT_EQ(4u, ctx.single_edit_cds_start);
    ASSERT_EQ(3u, ctx.single_edit_ref_len);
    ASSERT_EQ(0u, ctx.single_edit_alt_len);
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(kprop_delta_is_inframe_deletion_at(&delta, 2));

    memset(&edit, 0, sizeof edit);
    edit.cds_start = 4u; edit.ref_len = 6u; edit.ref = multi_cds + 3u;
    edit.alt_len = 0u; edit.alt = NULL; edit.variant_strand = 1;
    edit_set.edits = &edit; edit_set.count = 1u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(multi_cds, sizeof multi_cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds,
                                           ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(kprop_delta_is_inframe_deletion_at(&delta, 2));

    /* Non-codon-aligned pure in-frame deletion (cds_start 5, mid-codon): removes CDS 5-7 from
     * M K P G L F, merging codons 2-3 into M [T] G L F — a clean one-residue in-frame deletion
     * with a changed junction codon. VEP's codon-allele trim (AAACCC vs ACC -> alt empties,
     * ref len 3 % 3 == 0) calls this inframe_deletion regardless of alignment; the generalized
     * classifier now resolves it (the old codon-aligned-only slice bailed to unsupported).
     * protein_pos is the first affected codon, ((5-1)/3)+1 = 2. */
    edit.cds_start = 5u; edit.ref_len = 3u; edit.ref = multi_cds + 4u;
    edit.alt_len = 0u; edit.alt = NULL; edit.variant_strand = 1;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(multi_cds, sizeof multi_cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds,
                                           ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(kprop_delta_is_inframe_deletion_at(&delta, 2));

    /* Non-aligned deletion whose merged junction codon is a premature stop: ATG TCA GAA GGG
     * TTT AAA (M S E G F K), delete CDS 5-7 -> ATG TAA GGG ... = M * ...  The alt diff window
     * carries a stop the ref window lacked, so this is a stop-composite case (VEP stop_gained),
     * DEFERRED here rather than mislabelled inframe_deletion — exercising the alt-window nonstop
     * guard that the (now clean-only) generator never produces. */
    {
        static const uint8_t junction_stop_cds[18] = {
            'A','T','G',  'T','C','A',  'G','A','A',  'G','G','G',  'T','T','T',  'A','A','A'
        };
        edit.cds_start = 5u; edit.ref_len = 3u; edit.ref = junction_stop_cds + 4u;
        edit.alt_len = 0u; edit.alt = NULL; edit.variant_strand = 1;
        edit_set.edits = &edit; edit_set.count = 1u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(junction_stop_cds, sizeof junction_stop_cds,
                                               &edit_set, 1, DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(!delta.valid);
    }

    /* Removing the terminal stop codon: the ref diff window IS the stop, so this is stop_lost,
     * deferred here rather than mislabelled inframe_deletion (exercises the ref-window guard). */
    {
        static const uint8_t stop_del_cds[12] = {
            'A','T','G',  'A','A','A',  'C','C','C',  'T','A','A'
        };
        edit.cds_start = 10u; edit.ref_len = 3u; edit.ref = stop_del_cds + 9u;
        edit.alt_len = 0u; edit.alt = NULL; edit.variant_strand = 1;
        edit_set.edits = &edit; edit_set.count = 1u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(stop_del_cds, sizeof stop_del_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(!delta.valid);
    }

    /* Deleting the complete start codon is both inframe_deletion and start_lost. */
    {
        static const uint8_t start_del_cds[12] = {
            'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G'
        };
        edit.cds_start = 1u; edit.ref_len = 3u; edit.ref = start_del_cds;
        edit.alt_len = 0u; edit.alt = NULL; edit.variant_strand = 1;
        edit_set.edits = &edit; edit_set.count = 1u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(start_del_cds, sizeof start_del_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_deletion && delta.start_lost &&
               !delta.start_retained);
    }

    edits[0].cds_start = 10u; edits[0].ref_len = 3u; edits[0].ref = multi_cds + 9u;
    edits[0].alt_len = 0u; edits[0].alt = NULL; edits[0].variant_strand = 1;
    edits[1].cds_start = 4u; edits[1].ref_len = 3u; edits[1].ref = multi_cds + 3u;
    edits[1].alt_len = 0u; edits[1].alt = NULL; edits[1].variant_strand = 1;
    edit_set.edits = edits; edit_set.count = 2u;
    ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
              duckvep_coding_context_build(multi_cds, sizeof multi_cds, &edit_set, 1,
                                           DUCKVEP_CODON_TABLE_STANDARD,
                                           alt_cds, sizeof alt_cds,
                                           ref_pep, sizeof ref_pep,
                                           alt_pep, sizeof alt_pep, &ctx));
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(!delta.valid);

    memset(&ctx, 0, sizeof ctx);
    ctx.ref_cds = multi_cds; ctx.ref_cds_len = sizeof multi_cds;
    ctx.alt_cds = multi_cds; ctx.alt_cds_len = sizeof multi_cds;
    ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 6u;
    ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 6u;
    ctx.length_diff = -3;
    ctx.cds_changed = 1u;
    ctx.applied_edits = 1u; ctx.has_single_edit = 1u;
    ctx.single_edit_cds_start = 4u; ctx.single_edit_ref_len = 3u;
    ctx.single_edit_alt_len = 0u;
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(!delta.valid);

    memset(&ctx, 0, sizeof ctx);
    ctx.ref_cds = ambiguous_cds; ctx.ref_cds_len = sizeof ambiguous_cds;
    ctx.alt_cds = short_alt_cds; ctx.alt_cds_len = sizeof short_alt_cds;
    ctx.ref_peptide = fake_ref_pep; ctx.ref_peptide_len = 3u;
    ctx.alt_peptide = fake_alt_pep; ctx.alt_peptide_len = 2u;
    ctx.length_diff = -3;
    ctx.cds_changed = 1u;
    ctx.applied_edits = 1u; ctx.has_single_edit = 1u;
    ctx.single_edit_cds_start = 4u; ctx.single_edit_ref_len = 3u;
    ctx.single_edit_alt_len = 0u;
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(!delta.valid);

    memset(&ctx, 0, sizeof ctx);
    ctx.ref_cds = short_cds; ctx.ref_cds_len = sizeof short_cds;
    ctx.alt_cds = short_alt_cds; ctx.alt_cds_len = sizeof short_alt_cds;
    ctx.ref_peptide = x_ref_pep; ctx.ref_peptide_len = 3u;
    ctx.alt_peptide = fake_alt_pep; ctx.alt_peptide_len = 2u;
    ctx.length_diff = -3;
    ctx.cds_changed = 1u;
    ctx.applied_edits = 1u; ctx.has_single_edit = 1u;
    ctx.single_edit_cds_start = 4u; ctx.single_edit_ref_len = 3u;
    ctx.single_edit_alt_len = 0u;
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(!delta.valid);

    memset(&ctx, 0, sizeof ctx);
    ctx.ref_cds = multi_cds; ctx.ref_cds_len = sizeof multi_cds;
    ctx.alt_cds = multi_cds; ctx.alt_cds_len = sizeof multi_cds;
    ctx.ref_peptide = ref_pep; ctx.ref_peptide_len = 6u;
    ctx.alt_peptide = alt_pep; ctx.alt_peptide_len = 6u;
    ctx.length_diff = INT64_MIN;
    ctx.cds_changed = 1u;
    ctx.applied_edits = 1u; ctx.has_single_edit = 1u;
    ctx.single_edit_cds_start = 4u; ctx.single_edit_ref_len = 3u;
    ctx.single_edit_alt_len = 0u;
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(!delta.valid);
    PASS();
}

TEST coding_context_delta_inframe_insertion_known_scene(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','T'
    };
    static const uint8_t direct_cds[12] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'T','T','T'
    };
    static const uint8_t terminal_cds[12] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'T','A','A'
    };
    static const uint8_t internal_stop_cds[15] = {
        'A','T','G',  'T','A','A',  'C','C','C',  'G','G','G',  'T','T','T'
    };
    static const uint8_t stop_window_cds[12] = {
        'A','T','G',  'A','A','A',  'T','A','C',  'G','G','G'
    };
    static const uint8_t insert_gcc[3] = { 'G','C','C' };
    static const uint8_t insert_gcc2[6] = { 'G','C','C', 'G','C','C' };
    static const uint8_t insert_stop[3] = { 'T','A','A' };
    static const uint8_t insert_stop_before_ref_flank[3] = { 'A','T','A' };
    static const uint8_t insert_ref_flank_before_stop[6] = {
        'A','C','T', 'A','A','G'
    };
    static const uint8_t insert_atg[3] = { 'A','T','G' };
    static const uint8_t insert_start_retained[9] = {
        'G','C','G', 'T','T','G', 'G','C','A'
    };
    static const uint8_t insert_start_retained_stop[9] = {
        'G','T','C', 'A','T','C', 'C','T','A'
    };
    static const uint8_t terminal_insert_agc[3] = { 'A','G','C' };
    static const uint8_t terminal_insert_aggt[4] = { 'A','G','G','T' };
    static const uint8_t terminal_insert_stop_lost[11] = {
        'C','G','A','T','G','T','T','A','T','G','A'
    };
    static const uint8_t terminal_before_taaa[4] = { 'T','A','A','A' };
    static const uint8_t terminal_before_c[1] = { 'C' };
    static const uint8_t terminal_before_stop_gained[6] = {
        'A','A','A','T','A','A'
    };
    static const int8_t strands[2] = { 1, -1 };
    size_t case_idx;

    for (case_idx = 0u; case_idx < 2u; case_idx++) {
        struct kprop_coding s;
        duckvep_haplotype_edit_t edits[4];
        uint8_t alt_cds[32];
        uint8_t ref_pep[16];
        uint8_t alt_pep[16];
        duckvep_coding_context_t ctx;
        duckvep_sequence_delta_t delta;
        uint32_t anchor_cds;
        uint32_t i;

        memset(&s, 0, sizeof s);
        s.cds = cds; s.chrom = 0u; s.strand = strands[case_idx]; s.flags = 0u;
        s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
        s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u; s.eph = 0; s.eeph = 0;
        s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
        kprop_wire_coding_scene(&s, 15u);
        anchor_cds = s.strand > 0 ? 6u : 7u;
        s.vpos = kprop_genomic_pos_for_cds(&s, anchor_cds); s.vend = s.vpos;
        s.vkind = (uint8_t)DUCKVEP_KIND_INS;
        s.abytes[0] = (uint8_t)kprop_genomic_base_at(&s, s.vpos);
        s.abytes[1] = s.abytes[0];
        for (i = 0u; i < 3u; i++) {
            char b = s.strand > 0 ? (char)insert_gcc[i]
                                  : kprop_complement_base((char)insert_gcc[2u - i]);
            s.abytes[2u + i] = (uint8_t)b;
        }
        s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 4u;
        ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_OK,
                  duckvep_variant_coding_context_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                       0u, 0u, s.strand,
                                                       edits, 4u, alt_cds, sizeof alt_cds,
                                                       ref_pep, sizeof ref_pep,
                                                       alt_pep, sizeof alt_pep,
                                                       &ctx));
        ASSERT(ctx.has_single_edit);
        ASSERT_EQ(7u, ctx.single_edit_cds_start);
        ASSERT_EQ(0u, ctx.single_edit_ref_len);
        ASSERT_EQ(3u, ctx.single_edit_alt_len);
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_inframe_insertion_at(&delta, 3));
    }

    {
        duckvep_haplotype_edit_t edit;
        duckvep_edit_set_t edit_set;
        uint8_t alt_cds[32];
        uint8_t ref_pep[16];
        uint8_t alt_pep[16];
        duckvep_coding_context_t ctx;
        duckvep_sequence_delta_t delta;

        memset(&edit, 0, sizeof edit);
        edit.cds_start = 7u; edit.ref_len = 0u; edit.ref = NULL;
        edit.alt_len = 6u; edit.alt = insert_gcc2; edit.variant_strand = 1;
        edit_set.edits = &edit; edit_set.count = 1u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(direct_cds, sizeof direct_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_inframe_insertion_at(&delta, 3));

        edit.alt_len = 3u; edit.alt = insert_stop;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(direct_cds, sizeof direct_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_insertion && delta.stop_gained &&
               !delta.stop_lost && !delta.stop_retained);

        /* The two VEP predicates deliberately see different peptide values. Inserting
         * ATA before the last base of TAC makes the local peptide Y -> *Y. The in-frame
         * predicate truncates the alternate to '*', while protein_altering sees raw '*Y';
         * neither shape term applies, leaving stop_gained alone. */
        edit.cds_start = 9u;
        edit.alt_len = sizeof insert_stop_before_ref_flank;
        edit.alt = insert_stop_before_ref_flank;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      stop_window_cds, sizeof stop_window_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_gained && !delta.inframe_insertion &&
               !delta.protein_altering && !delta.frameshift);

        /* Nearby positive control: inserting ACTAAG before the second base of TAC makes
         * Y -> Y*D. The preserved Y occurs before the new stop, so VEP emits both terms. */
        edit.cds_start = 8u;
        edit.alt_len = sizeof insert_ref_flank_before_stop;
        edit.alt = insert_ref_flank_before_stop;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      stop_window_cds, sizeof stop_window_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_gained && delta.inframe_insertion &&
               !delta.protein_altering && !delta.frameshift);

        /* VEP truncates *both* sides to the first stop for the exact-stop
         * insertion guard. Inserting GCC inside an internal TAA makes the local
         * peptide * -> *P: stop_retained applies, but inframe_insertion does not. */
        edit.cds_start = 6u;
        edit.alt_len = sizeof insert_gcc;
        edit.alt = insert_gcc;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      internal_stop_cds, sizeof internal_stop_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_retained && !delta.inframe_insertion &&
               !delta.protein_altering && !delta.stop_gained && !delta.stop_lost);

        /* Moving the same ATG insertion across the terminal TAA distinguishes
         * stop_lost from inframe_insertion&stop_retained_variant. */
        edit.alt = insert_atg;
        edit.alt_len = sizeof insert_atg;
        edit.cds_start = 11u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_lost && !delta.stop_gained &&
               !delta.stop_retained && !delta.inframe_insertion);

        edit.cds_start = 12u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.stop_retained && delta.inframe_insertion &&
               !delta.stop_gained && !delta.stop_lost);

        /* VEP's terminal-stop insertion predicates are deliberately not reducible
         * to length modulo three. These are executable witnesses from the pinned
         * VEP-116 state machine; see design/duckvep_errata.md. */
        edit.cds_start = 11u;
        edit.alt_len = sizeof terminal_insert_agc;
        edit.alt = terminal_insert_agc;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_insertion && delta.stop_retained &&
               !delta.frameshift && !delta.stop_lost && !delta.coding_unknown);

        edit.alt_len = sizeof terminal_insert_aggt;
        edit.alt = terminal_insert_aggt;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_insertion && delta.coding_unknown &&
               !delta.frameshift && !delta.stop_lost && !delta.stop_retained);

        edit.alt_len = sizeof terminal_insert_stop_lost;
        edit.alt = terminal_insert_stop_lost;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_insertion && delta.stop_lost &&
               !delta.frameshift && !delta.stop_retained && !delta.coding_unknown);

        edit.cds_start = 10u;
        edit.alt_len = sizeof terminal_before_taaa;
        edit.alt = terminal_before_taaa;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_insertion && delta.stop_retained &&
               !delta.frameshift && !delta.stop_gained);

        edit.alt_len = sizeof terminal_before_c;
        edit.alt = terminal_before_c;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.frameshift && !delta.inframe_insertion &&
               !delta.stop_retained && !delta.stop_gained);

        edit.alt_len = sizeof terminal_before_stop_gained;
        edit.alt = terminal_before_stop_gained;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, 1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_insertion && delta.stop_gained &&
               !delta.frameshift && !delta.stop_retained);

        /* The coding context is transcript-oriented; the same predicate state
         * must be independent of genomic strand. */
        edit.cds_start = 11u;
        edit.alt_len = sizeof terminal_insert_aggt;
        edit.alt = terminal_insert_aggt;
        edit.variant_strand = -1;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(
                      terminal_cds, sizeof terminal_cds, &edit_set, -1,
                      DUCKVEP_CODON_TABLE_STANDARD, alt_cds, sizeof alt_cds,
                      ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.inframe_insertion && delta.coding_unknown &&
               !delta.frameshift && !delta.stop_lost && !delta.stop_retained);
        edit.variant_strand = 1;

        /* Mid-codon insertion that PRESERVES the flanking residue is inframe_insertion, not
         * protein_altering and not a bail: GCC after CDS 5 makes codon 2 AAA->AAG (still Lys)
         * and inserts Pro, so the ref window is empty. VEP calls this inframe_insertion; the old
         * codon-boundary-only classifier bailed. protein_pos = (5/3)+1 = 2. */
        edit.alt_len = sizeof insert_gcc;
        edit.alt = insert_gcc;
        edit.cds_start = 6u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(direct_cds, sizeof direct_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_inframe_insertion_at(&delta, 2));

        /* Insertion immediately after the complete start codon (before_cds == 3): the start Met
         * is preserved, Ala inserted after it -> inframe_insertion at protein_pos 2. The old
         * classifier required before_cds > 3 and bailed. */
        edit.cds_start = 4u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(direct_cds, sizeof direct_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_inframe_insertion_at(&delta, 2));

        /* VEP 116 evaluates the start and insertion predicates independently.
         * Inserting nine bases before CDS base 3 can preserve ATG while adding
         * residues: both start_retained_variant and inframe_insertion apply. */
        edit.cds_start = 3u;
        edit.alt_len = sizeof insert_start_retained;
        edit.alt = insert_start_retained;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(direct_cds, sizeof direct_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.start_retained && delta.inframe_insertion &&
               !delta.start_lost && !delta.stop_gained);

        /* The same state can gain a stop: AT + GTCATCCTA + the original G
         * translates MSS*. VEP emits all three facts. */
        edit.alt = insert_start_retained_stop;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(direct_cds, sizeof direct_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.start_retained && delta.inframe_insertion &&
               delta.stop_gained && !delta.start_lost);

        /* Insertion between the penultimate and last codon (before_cds == 9, one short of the
         * old before_cds < cds_len-3 bound): flank-preserving -> inframe_insertion at pp 4. */
        edit.cds_start = 10u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(direct_cds, sizeof direct_cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_inframe_insertion_at(&delta, 4));
    }
    PASS();
}

TEST coding_context_delta_delins_known_scene(void) {
    static uint8_t cds[18] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','T',  'G','A','A'
    };
    static const uint8_t alt3[3] = { 'G','C','C' };
    static const uint8_t alt6[6] = { 'G','C','C',  'G','C','T' };
    static const uint8_t stop_alt6[6] = { 'G','C','C',  'T','A','A' };
    static const uint8_t ambiguous_cds[18] = {
        'A','T','G',  'A','A','A',  'N','N','N',  'G','G','G',  'T','T','T',  'G','A','A'
    };
    static const uint8_t ambiguous_alt_cds[21] = {
        'A','T','G',  'A','A','A',  'G','C','C',  'G','C','T',
        'G','G','G',  'T','T','T',  'G','A','A'
    };
    static const uint8_t fake_ref_pep[6] = { 'M', 'K', 'X', 'G', 'F', 'E' };
    static const uint8_t fake_alt_pep[7] = { 'M', 'K', 'A', 'A', 'G', 'F', 'E' };
    static const int8_t strands[2] = { 1, -1 };
    size_t case_idx;

    for (case_idx = 0u; case_idx < 2u; case_idx++) {
        struct kprop_coding s;
        duckvep_haplotype_edit_t edits[4];
        uint8_t alt_cds[40];
        uint8_t ref_pep[20];
        uint8_t alt_pep[20];
        duckvep_coding_context_t ctx;
        duckvep_sequence_delta_t delta;
        uint32_t cds_start;
        uint32_t i;

        memset(&s, 0, sizeof s);
        s.cds = cds; s.chrom = 0u; s.strand = strands[case_idx]; s.flags = 0u;
        s.tstart = 1000u; s.tend = 1017u; s.cds_s = 1000u; s.cds_e = 1017u;
        s.es = 1000u; s.ee = 1017u; s.ecds = 1u; s.ecde = 18u; s.eph = 0; s.eeph = 0;
        s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
        kprop_wire_coding_scene(&s, 18u);

        cds_start = 7u;
        s.vpos = s.strand > 0 ? kprop_genomic_pos_for_cds(&s, cds_start)
                              : kprop_genomic_pos_for_cds(&s, cds_start + 2u);
        s.vend = s.vpos + 2u; s.vkind = (uint8_t)DUCKVEP_KIND_INDEL;
        for (i = 0u; i < 3u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
        kprop_fill_variant_alt_from_tx(&s, 3u, alt6, 6u);
        s.roff = 0u; s.aoff = 3u; s.rlen = 3u; s.alen = 6u;
        ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_OK,
                  duckvep_variant_coding_context_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                       0u, 0u, s.strand,
                                                       edits, 4u, alt_cds, sizeof alt_cds,
                                                       ref_pep, sizeof ref_pep,
                                                       alt_pep, sizeof alt_pep,
                                                       &ctx));
        ASSERT(ctx.has_single_edit);
        ASSERT_EQ(7u, ctx.single_edit_cds_start);
        ASSERT_EQ(3u, ctx.single_edit_ref_len);
        ASSERT_EQ(6u, ctx.single_edit_alt_len);
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_protein_altering_at(&delta, 3));

        s.vpos = s.strand > 0 ? kprop_genomic_pos_for_cds(&s, cds_start)
                              : kprop_genomic_pos_for_cds(&s, cds_start + 5u);
        s.vend = s.vpos + 5u;
        for (i = 0u; i < 6u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
        kprop_fill_variant_alt_from_tx(&s, 6u, alt3, 3u);
        s.roff = 0u; s.aoff = 6u; s.rlen = 6u; s.alen = 3u;
        ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_OK,
                  duckvep_variant_coding_context_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                       0u, 0u, s.strand,
                                                       edits, 4u, alt_cds, sizeof alt_cds,
                                                       ref_pep, sizeof ref_pep,
                                                       alt_pep, sizeof alt_pep,
                                                       &ctx));
        ASSERT(ctx.has_single_edit);
        ASSERT_EQ(7u, ctx.single_edit_cds_start);
        ASSERT_EQ(6u, ctx.single_edit_ref_len);
        ASSERT_EQ(3u, ctx.single_edit_alt_len);
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_protein_altering_at(&delta, 3));
    }

    {
        duckvep_haplotype_edit_t edit;
        duckvep_haplotype_edit_t edits[2];
        duckvep_edit_set_t edit_set;
        uint8_t alt_cds[40];
        uint8_t ref_pep[20];
        uint8_t alt_pep[20];
        duckvep_coding_context_t ctx;
        duckvep_sequence_delta_t delta;

        memset(&edit, 0, sizeof edit);
        edit.ref_len = 3u; edit.ref = cds; edit.alt_len = 6u; edit.alt = alt6;
        edit.variant_strand = 1; edit_set.edits = &edit; edit_set.count = 1u;

        edit.cds_start = 1u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cds, sizeof cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.start_lost && !delta.inframe_insertion &&
               !delta.protein_altering);

        edit.cds_start = 16u; edit.ref = cds + 15u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cds, sizeof cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_protein_altering_at(&delta, 6));

        edit.cds_start = 8u; edit.ref = cds + 7u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cds, sizeof cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(kprop_delta_is_protein_altering_at(&delta, 3));

        edit.cds_start = 7u; edit.ref = cds + 6u; edit.alt = stop_alt6;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cds, sizeof cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(delta.valid && delta.protein_altering && delta.stop_gained &&
               !delta.inframe_insertion && !delta.inframe_deletion);

        edits[0].cds_start = 13u; edits[0].ref_len = 3u; edits[0].ref = cds + 12u;
        edits[0].alt_len = 3u; edits[0].alt = alt3; edits[0].variant_strand = 1;
        edits[1].cds_start = 7u; edits[1].ref_len = 3u; edits[1].ref = cds + 6u;
        edits[1].alt_len = 6u; edits[1].alt = alt6; edits[1].variant_strand = 1;
        edit_set.edits = edits; edit_set.count = 2u;
        ASSERT_EQ(DUCKVEP_CODING_CONTEXT_OK,
                  duckvep_coding_context_build(cds, sizeof cds, &edit_set, 1,
                                               DUCKVEP_CODON_TABLE_STANDARD,
                                               alt_cds, sizeof alt_cds,
                                               ref_pep, sizeof ref_pep,
                                               alt_pep, sizeof alt_pep, &ctx));
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(!delta.valid);

        memset(&ctx, 0, sizeof ctx);
        ctx.ref_cds = ambiguous_cds; ctx.ref_cds_len = sizeof ambiguous_cds;
        ctx.alt_cds = ambiguous_alt_cds; ctx.alt_cds_len = sizeof ambiguous_alt_cds;
        ctx.ref_peptide = fake_ref_pep; ctx.ref_peptide_len = sizeof fake_ref_pep;
        ctx.alt_peptide = fake_alt_pep; ctx.alt_peptide_len = sizeof fake_alt_pep;
        ctx.length_diff = 3; ctx.cds_changed = 1u;
        ctx.applied_edits = 1u; ctx.has_single_edit = 1u;
        ctx.single_edit_cds_start = 7u; ctx.single_edit_ref_len = 3u;
        ctx.single_edit_alt_len = 6u;
        ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
                  duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
        ASSERT(!delta.valid);
    }

    {
        struct kprop_coding s;
        duckvep_model_t *model = NULL;
        duckvep_options_t *opts = NULL;
        duckvep_workspace_t *ws = NULL;
        const duckvep_workspace_delta_route_stats_t *stats;
        duckvep_error_t err;
        duckvep_consequence_t rows[2];
        duckvep_result_builder_t rb;
        uint32_t i;

        memset(&s, 0, sizeof s); memset(&err, 0, sizeof err);
        s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
        s.tstart = 1000u; s.tend = 1017u; s.cds_s = 1000u; s.cds_e = 1017u;
        s.es = 1000u; s.ee = 1017u; s.ecds = 1u; s.ecde = 18u; s.eph = 0; s.eeph = 0;
        s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
        kprop_wire_coding_scene(&s, 18u);
        s.vpos = kprop_genomic_pos_for_cds(&s, 7u); s.vend = s.vpos + 2u;
        s.vkind = (uint8_t)DUCKVEP_KIND_INDEL;
        for (i = 0u; i < 3u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
        kprop_fill_variant_alt_from_tx(&s, 3u, alt6, 6u);
        s.roff = 0u; s.aoff = 3u; s.rlen = 3u; s.alen = 6u;
        ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&s.tx, &s.ex, &s.seq, &model, &err));
        ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
        ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
        duckvep_workspace_delta_route_stats_reset(ws);
        duckvep_result_builder_init(&rb, rows, 2u);
        ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &s.v, opts, ws, &rb, &err));
        ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
        /* Replacing P with AA increases the CDS by one codon but preserves neither
         * peptide edge. VEP calls the local shape protein_altering_variant, not an
         * in-frame insertion inferred from net length alone. */
        ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING), rows[0].consequence_mask);
        stats = duckvep_workspace_delta_route_stats(ws);
        ASSERT(stats != NULL);
        ASSERT_EQ(0u, stats->substitution_context);
        ASSERT_EQ(0u, stats->mnv_direct_fallback);
        ASSERT_EQ(0u, stats->del_context);
        ASSERT_EQ(0u, stats->del_direct_fallback);
        ASSERT_EQ(0u, stats->ins_context);
        ASSERT_EQ(0u, stats->ins_direct_fallback);
        ASSERT_EQ(1u, stats->indel_context);
        ASSERT_EQ(0u, stats->indel_direct_fallback);
        duckvep_workspace_close(ws);
        duckvep_options_close(opts);
        duckvep_model_close(model);
    }
    PASS();
}

TEST sequence_delta_with_scratch_indel_known_scene(void) {
    static uint8_t cds[18] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','T',  'G','A','A'
    };
    static const uint8_t alt3[3] = { 'G','C','C' };
    static const uint8_t alt6[6] = { 'G','C','C',  'G','C','T' };
    static const uint8_t alt2[2] = { 'G','T' };
    static const int8_t strands[2] = { 1, -1 };
    size_t case_idx;

    for (case_idx = 0u; case_idx < 2u; case_idx++) {
        struct kprop_coding s;
        duckvep_haplotype_edit_t edits[4];
        uint8_t alt_cds[40];
        uint8_t ref_pep[20];
        uint8_t alt_pep[20];
        duckvep_delta_scratch_t scratch;
        duckvep_sequence_delta_t delta;
        uint32_t cds_start = 7u;
        uint32_t i;

        memset(&s, 0, sizeof s);
        s.cds = cds; s.chrom = 0u; s.strand = strands[case_idx]; s.flags = 0u;
        s.tstart = 1000u; s.tend = 1017u; s.cds_s = 1000u; s.cds_e = 1017u;
        s.es = 1000u; s.ee = 1017u; s.ecds = 1u; s.ecde = 18u; s.eph = 0; s.eeph = 0;
        s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
        kprop_wire_coding_scene(&s, 18u);
        memset(&scratch, 0, sizeof scratch);
        scratch.edits = edits; scratch.edits_cap = 4u;
        scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
        scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
        scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;

        s.vpos = s.strand > 0 ? kprop_genomic_pos_for_cds(&s, cds_start)
                              : kprop_genomic_pos_for_cds(&s, cds_start + 2u);
        s.vend = s.vpos + 2u; s.vkind = (uint8_t)DUCKVEP_KIND_INDEL;
        for (i = 0u; i < 3u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
        kprop_fill_variant_alt_from_tx(&s, 3u, alt6, 6u);
        s.roff = 0u; s.aoff = 3u; s.rlen = 3u; s.alen = 6u;
        duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_INDEL, &s.tx, &s.ex,
                                                 &s.seq, &s.v, 0u, 0u, s.vpos,
                                                 s.strand, &scratch, &delta);
        ASSERT(kprop_delta_is_protein_altering_at(&delta, 3));

        s.vpos = s.strand > 0 ? kprop_genomic_pos_for_cds(&s, cds_start)
                              : kprop_genomic_pos_for_cds(&s, cds_start + 5u);
        s.vend = s.vpos + 5u;
        for (i = 0u; i < 6u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
        kprop_fill_variant_alt_from_tx(&s, 6u, alt3, 3u);
        s.roff = 0u; s.aoff = 6u; s.rlen = 6u; s.alen = 3u;
        duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_INDEL, &s.tx, &s.ex,
                                                 &s.seq, &s.v, 0u, 0u, s.vpos,
                                                 s.strand, &scratch, &delta);
        ASSERT(kprop_delta_is_protein_altering_at(&delta, 3));

        s.vpos = s.strand > 0 ? kprop_genomic_pos_for_cds(&s, cds_start)
                              : kprop_genomic_pos_for_cds(&s, cds_start);
        s.vend = s.vpos; s.vkind = (uint8_t)DUCKVEP_KIND_INDEL;
        s.abytes[0] = (uint8_t)kprop_genomic_base_at(&s, s.vpos);
        kprop_fill_variant_alt_from_tx(&s, 1u, alt2, 2u);
        s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 2u;
        duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_INDEL, &s.tx, &s.ex,
                                                 &s.seq, &s.v, 0u, 0u, s.vpos,
                                                 s.strand, &scratch, &delta);
        /* +1 delins at codon 3 (body) with the ATG start intact: the general
         * CodingContext now resolves the frameshift the direct body-only path skipped. */
        ASSERT(kprop_delta_is_frameshift_at(&delta, 3));
        duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_INDEL, &s.tx, &s.ex,
                                                 &s.seq, &s.v, 0u, 0u, s.vpos,
                                                 s.strand, NULL, &delta);
        ASSERT(delta.valid);
        ASSERT(delta.frameshift);
        ASSERT(!delta.inframe_insertion);
        ASSERT(!delta.inframe_deletion);
    }
    PASS();
}

TEST annotate_delins_boundary_no_route_known_scene(void) {
    static uint8_t cds[18] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','T',  'G','A','A'
    };
    static const uint8_t alt3[3] = { 'G','C','C' };
    static const uint8_t alt6[6] = { 'G','C','C',  'G','C','T' };
    static const uint32_t starts[4] = { 1u, 16u, 13u, 8u };
    static const uint32_t ref_lens[4] = { 3u, 3u, 6u, 3u };
    static const uint32_t alt_lens[4] = { 6u, 6u, 3u, 6u };
    static const uint64_t expected_masks[4] = {
        DUCKVEP_SO(DUCKVEP_SO_START_LOST),
        DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING),
        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE),
        DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING)
    };
    static const int8_t strands[2] = { 1, -1 };
    size_t case_idx;
    size_t strand_idx;

    for (strand_idx = 0u; strand_idx < 2u; strand_idx++) {
        for (case_idx = 0u; case_idx < 4u; case_idx++) {
            struct kprop_coding s;
            duckvep_model_t *model = NULL;
            duckvep_options_t *opts = NULL;
            duckvep_workspace_t *ws = NULL;
            const duckvep_workspace_delta_route_stats_t *stats;
            duckvep_error_t err;
            duckvep_consequence_t rows[2];
            duckvep_result_builder_t rb;
            const uint8_t *alt_tx = alt_lens[case_idx] == 3u ? alt3 : alt6;
            uint32_t cds_start = starts[case_idx];
            uint32_t ref_len = ref_lens[case_idx];
            uint32_t alt_len = alt_lens[case_idx];
            uint32_t i;

            memset(&s, 0, sizeof s); memset(&err, 0, sizeof err);
            s.cds = cds; s.chrom = 0u; s.strand = strands[strand_idx]; s.flags = 0u;
            s.tstart = 1000u; s.tend = 1017u; s.cds_s = 1000u; s.cds_e = 1017u;
            s.es = 1000u; s.ee = 1017u; s.ecds = 1u; s.ecde = 18u; s.eph = 0; s.eeph = 0;
            s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
            kprop_wire_coding_scene(&s, 18u);
            s.vpos = s.strand > 0 ? kprop_genomic_pos_for_cds(&s, cds_start)
                                  : kprop_genomic_pos_for_cds(&s, cds_start + ref_len - 1u);
            s.vend = s.vpos + ref_len - 1u; s.vkind = (uint8_t)DUCKVEP_KIND_INDEL;
            for (i = 0u; i < ref_len; i++) {
                s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
            }
            kprop_fill_variant_alt_from_tx(&s, ref_len, alt_tx, alt_len);
            s.roff = 0u; s.aoff = ref_len;
            s.rlen = (uint16_t)ref_len; s.alen = (uint16_t)alt_len;

            ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&s.tx, &s.ex, &s.seq, &model, &err));
            ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
            ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
            duckvep_workspace_delta_route_stats_reset(ws);
            duckvep_result_builder_init(&rb, rows, 2u);
            ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &s.v, opts, ws, &rb, &err));
            ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
            ASSERT_EQ(expected_masks[case_idx], rows[0].consequence_mask);
            stats = duckvep_workspace_delta_route_stats(ws);
            ASSERT(stats != NULL);
            ASSERT_EQ(0u, stats->substitution_context);
            ASSERT_EQ(0u, stats->mnv_direct_fallback);
            ASSERT_EQ(0u, stats->del_context);
            ASSERT_EQ(0u, stats->del_direct_fallback);
            ASSERT_EQ(0u, stats->ins_context);
            ASSERT_EQ(0u, stats->ins_direct_fallback);

            duckvep_workspace_close(ws);
            duckvep_options_close(opts);
            duckvep_model_close(model);
        }
    }
    PASS();
}

TEST annotate_inframe_insertion_route_known_scene(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','T'
    };
    static const int8_t strands[2] = { 1, -1 };
    size_t case_idx;

    for (case_idx = 0u; case_idx < 2u; case_idx++) {
        struct kprop_coding s;
        duckvep_model_t *model = NULL;
        duckvep_options_t *opts = NULL;
        duckvep_workspace_t *ws = NULL;
        const duckvep_workspace_delta_route_stats_t *stats;
        duckvep_error_t err;
        duckvep_consequence_t rows[2];
        duckvep_result_builder_t rb;
        uint8_t alt_tx[3] = { 'G','C','C' };
        uint32_t anchor_cds;
        uint32_t i;

        memset(&s, 0, sizeof s); memset(&err, 0, sizeof err);
        s.cds = cds; s.chrom = 0u; s.strand = strands[case_idx]; s.flags = 0u;
        s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
        s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u; s.eph = 0; s.eeph = 0;
        s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
        kprop_wire_coding_scene(&s, 15u);
        anchor_cds = s.strand > 0 ? 6u : 7u;
        s.vpos = kprop_genomic_pos_for_cds(&s, anchor_cds); s.vend = s.vpos;
        s.vkind = (uint8_t)DUCKVEP_KIND_INS;
        s.abytes[0] = (uint8_t)kprop_genomic_base_at(&s, s.vpos);
        s.abytes[1] = s.abytes[0];
        for (i = 0u; i < 3u; i++) {
            char b = s.strand > 0 ? (char)alt_tx[i]
                                  : kprop_complement_base((char)alt_tx[2u - i]);
            s.abytes[2u + i] = (uint8_t)b;
        }
        s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 4u;

        ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&s.tx, &s.ex, &s.seq, &model, &err));
        ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
        ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
        duckvep_workspace_delta_route_stats_reset(ws);
        duckvep_result_builder_init(&rb, rows, 2u);
        ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &s.v, opts, ws, &rb, &err));
        ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
        ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION), rows[0].consequence_mask);
        ASSERT_EQ(-1, rows[0].cdna_pos);
        ASSERT_EQ(-1, rows[0].cds_pos);
        ASSERT_EQ(3, rows[0].protein_pos);
        stats = duckvep_workspace_delta_route_stats(ws);
        ASSERT(stats != NULL);
        ASSERT_EQ(0u, stats->substitution_context);
        ASSERT_EQ(0u, stats->mnv_direct_fallback);
        ASSERT_EQ(0u, stats->del_context);
        ASSERT_EQ(0u, stats->del_direct_fallback);
        ASSERT_EQ(1u, stats->ins_context);
        ASSERT_EQ(0u, stats->ins_direct_fallback);

        duckvep_workspace_close(ws);
        duckvep_options_close(opts);
        duckvep_model_close(model);
    }

    {
        struct kprop_coding s;
        duckvep_model_t *model = NULL;
        duckvep_options_t *opts = NULL;
        duckvep_workspace_t *ws = NULL;
        const duckvep_workspace_delta_route_stats_t *stats;
        duckvep_error_t err;
        duckvep_consequence_t rows[2];
        duckvep_result_builder_t rb;
        uint8_t alt_tx[3] = { 'T','A','A' };
        uint32_t i;

        memset(&s, 0, sizeof s); memset(&err, 0, sizeof err);
        s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
        s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
        s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u; s.eph = 0; s.eeph = 0;
        s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
        kprop_wire_coding_scene(&s, 15u);
        s.vpos = kprop_genomic_pos_for_cds(&s, 6u); s.vend = s.vpos;
        s.vkind = (uint8_t)DUCKVEP_KIND_INS;
        s.abytes[0] = (uint8_t)kprop_genomic_base_at(&s, s.vpos);
        s.abytes[1] = s.abytes[0];
        for (i = 0u; i < 3u; i++) s.abytes[2u + i] = alt_tx[i];
        s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 4u;

        ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&s.tx, &s.ex, &s.seq, &model, &err));
        ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
        ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
        duckvep_workspace_delta_route_stats_reset(ws);
        duckvep_result_builder_init(&rb, rows, 2u);
        ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &s.v, opts, ws, &rb, &err));
        ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
        ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION) |
                      DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED),
                  rows[0].consequence_mask);
        stats = duckvep_workspace_delta_route_stats(ws);
        ASSERT(stats != NULL);
        ASSERT_EQ(1u, stats->ins_context);
        ASSERT_EQ(0u, stats->ins_direct_fallback);

        duckvep_workspace_close(ws);
        duckvep_options_close(opts);
        duckvep_model_close(model);
    }
    PASS();
}

TEST annotate_inframe_deletion_route_known_scene(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','T'
    };
    static const int8_t strands[2] = { 1, -1 };
    size_t case_idx;

    for (case_idx = 0u; case_idx < 2u; case_idx++) {
        struct kprop_coding s;
        duckvep_model_t *model = NULL;
        duckvep_options_t *opts = NULL;
        duckvep_workspace_t *ws = NULL;
        const duckvep_workspace_delta_route_stats_t *stats;
        duckvep_error_t err;
        duckvep_consequence_t rows[2];
        duckvep_result_builder_t rb;
        uint32_t anchor_cds;
        uint32_t i;

        memset(&s, 0, sizeof s); memset(&err, 0, sizeof err);
        s.cds = cds; s.chrom = 0u; s.strand = strands[case_idx]; s.flags = 0u;
        s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
        s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u; s.eph = 0; s.eeph = 0;
        s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
        kprop_wire_coding_scene(&s, 15u);
        anchor_cds = s.strand > 0 ? 3u : 7u;
        s.vpos = kprop_genomic_pos_for_cds(&s, anchor_cds); s.vend = s.vpos + 3u;
        s.vkind = (uint8_t)DUCKVEP_KIND_DEL;
        for (i = 0u; i < 4u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
        s.abytes[4] = s.abytes[0];
        s.roff = 0u; s.aoff = 4u; s.rlen = 4u; s.alen = 1u;

        ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&s.tx, &s.ex, &s.seq, &model, &err));
        ASSERT_EQ(DUCKVEP_OK, duckvep_options_open(NULL, &opts, &err));
        ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
        duckvep_workspace_delta_route_stats_reset(ws);
        duckvep_result_builder_init(&rb, rows, 2u);
        ASSERT_EQ(DUCKVEP_OK, duckvep_annotate_tile(model, &s.v, opts, ws, &rb, &err));
        ASSERT_EQ(1u, duckvep_result_builder_count(&rb));
        ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION), rows[0].consequence_mask);
        ASSERT_EQ(-1, rows[0].cdna_pos);
        ASSERT_EQ(-1, rows[0].cds_pos);
        ASSERT_EQ(2, rows[0].protein_pos);
        stats = duckvep_workspace_delta_route_stats(ws);
        ASSERT(stats != NULL);
        ASSERT_EQ(0u, stats->substitution_context);
        ASSERT_EQ(0u, stats->mnv_direct_fallback);
        ASSERT_EQ(1u, stats->del_context);
        ASSERT_EQ(0u, stats->del_direct_fallback);

        duckvep_workspace_close(ws);
        duckvep_options_close(opts);
        duckvep_model_close(model);
    }
    PASS();
}

TEST workspace_delta_scratch_caps_known(void) {
    static const uint16_t chrom[3] = {0u, 0u, 0u};
    static const uint32_t start1[3] = {100u, 200u, 300u};
    static const uint32_t end1[3] = {100u, 208u, 313u};
    static const int8_t strand[3] = {1, 1, -1};
    static const uint64_t flags[3] = {0u, 0u, 0u};
    static const uint32_t exon_off[3] = {0u, 0u, 1u};
    static const uint16_t exon_cnt[3] = {0u, 1u, 1u};
    static const uint32_t cds_s[3] = {0u, 200u, 300u};
    static const uint32_t cds_e[3] = {0u, 208u, 313u};
    static const uint64_t cds_off[3] = {0u, 0u, 9u};
    static const uint32_t cds_len[3] = {0u, 9u, 14u};
    static const uint8_t table[3] = {0u, 1u, 1u};
    static const uint32_t exon_start[2] = {200u, 300u};
    static const uint32_t exon_end[2] = {208u, 313u};
    static const uint32_t exon_cdna_start[2] = {1u, 1u};
    static const uint32_t exon_cdna_end[2] = {9u, 14u};
    static const int8_t exon_phase[2] = {0, 0};
    static const uint8_t cds_bytes[23] = {
        'A','T','G','G','A','A','T','A','A',
        'A','T','G','C','C','C','G','G','G','T','T','T','A','A'
    };
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_sequence_pool_t seq;
    duckvep_model_t *model = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_delta_scratch_t *scratch;
    duckvep_error_t err;

    memset(&tx, 0, sizeof tx); memset(&ex, 0, sizeof ex); memset(&seq, 0, sizeof seq);
    memset(&err, 0, sizeof err);
    tx.chrom_id = chrom; tx.start1 = start1; tx.end1 = end1; tx.strand = strand;
    tx.flags = flags; tx.exon_offset = exon_off; tx.exon_count = exon_cnt;
    tx.cds_start1 = cds_s; tx.cds_end1 = cds_e; tx.transcript_count = 3u;
    ex.start1 = exon_start; ex.end1 = exon_end;
    ex.cdna_start1 = exon_cdna_start; ex.cdna_end1 = exon_cdna_end;
    ex.phase = exon_phase; ex.end_phase = exon_phase; ex.exon_count = 2u;
    seq.cds_bytes = cds_bytes; seq.cds_bytes_len = sizeof cds_bytes;
    seq.cds_offset = cds_off; seq.cds_length = cds_len; seq.codon_table = table;
    seq.transcript_count = 3u;
    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, &seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    scratch = duckvep_workspace_delta_scratch(ws);
    ASSERT(scratch != NULL);
    ASSERT_EQ((size_t)14u + (size_t)UINT16_MAX, scratch->alt_cds_cap);
    ASSERT_EQ(5u, scratch->ref_peptide_cap);
    ASSERT_EQ(((size_t)14u + (size_t)UINT16_MAX) / 3u + 1u, scratch->alt_peptide_cap);
    ASSERT_EQ(7u, scratch->edits_cap);
    ASSERT(scratch->alt_cds != NULL);
    ASSERT(scratch->ref_peptide != NULL);
    ASSERT(scratch->alt_peptide != NULL);
    ASSERT(scratch->edits != NULL);
    duckvep_workspace_close(ws); ws = NULL;
    duckvep_model_close(model); model = NULL;

    memset(&seq, 0, sizeof seq); memset(&err, 0, sizeof err);
    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, NULL, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    scratch = duckvep_workspace_delta_scratch(ws);
    ASSERT(scratch != NULL);
    ASSERT_EQ(0u, scratch->alt_cds_cap);
    ASSERT_EQ(0u, scratch->ref_peptide_cap);
    ASSERT_EQ(0u, scratch->alt_peptide_cap);
    ASSERT_EQ(0u, scratch->edits_cap);
    ASSERT(scratch->alt_cds == NULL);
    ASSERT(scratch->ref_peptide == NULL);
    ASSERT(scratch->alt_peptide == NULL);
    ASSERT(scratch->edits == NULL);
    duckvep_workspace_close(ws); ws = NULL;
    duckvep_model_close(model); model = NULL;

    memset(&seq, 0, sizeof seq); memset(&err, 0, sizeof err);
    seq.cds_offset = cds_off; seq.cds_length = cds_len; seq.codon_table = table;
    seq.cds_bytes = NULL; seq.cds_bytes_len = 0u; seq.transcript_count = 3u;
    {
        static const uint64_t zero_off[3] = {0u, 0u, 0u};
        static const uint32_t zero_len[3] = {0u, 0u, 0u};
        seq.cds_offset = zero_off;
        seq.cds_length = zero_len;
        ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &ex, &seq, &model, &err));
    }
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    scratch = duckvep_workspace_delta_scratch(ws);
    ASSERT(scratch != NULL);
    ASSERT_EQ(0u, scratch->alt_cds_cap);
    ASSERT_EQ(0u, scratch->ref_peptide_cap);
    ASSERT_EQ(0u, scratch->alt_peptide_cap);
    ASSERT_EQ(0u, scratch->edits_cap);
    ASSERT(scratch->alt_cds == NULL);
    ASSERT(scratch->ref_peptide == NULL);
    ASSERT(scratch->alt_peptide == NULL);
    ASSERT(scratch->edits == NULL);
    duckvep_workspace_close(ws);
    duckvep_model_close(model);
    PASS();
}

TEST workspace_delta_scratch_builds_lengthening_context(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'G','A','A',  'C','C','C',  'G','G','G',  'T','T','T'
    };
    struct kprop_coding s;
    duckvep_model_t *model = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_delta_scratch_t *ws_scratch;
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;
    duckvep_error_t err;
    uint8_t alt_tx[3] = {'G','C','C'};
    uint32_t i;

    memset(&s, 0, sizeof s); memset(&err, 0, sizeof err);
    s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
    s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, 15u);
    s.vpos = kprop_genomic_pos_for_cds(&s, 6u); s.vend = s.vpos;
    s.vkind = (uint8_t)DUCKVEP_KIND_INS;
    s.abytes[0] = (uint8_t)kprop_genomic_base_at(&s, s.vpos);
    s.abytes[1] = s.abytes[0];
    for (i = 0u; i < 3u; i++) s.abytes[2u + i] = alt_tx[i];
    s.roff = 0u; s.aoff = 1u; s.rlen = 1u; s.alen = 4u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&s.tx, &s.ex, &s.seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    ws_scratch = duckvep_workspace_delta_scratch(ws);
    ASSERT(ws_scratch != NULL);
    ASSERT(ws_scratch->alt_cds_cap >= (size_t)s.cds_lenv + (size_t)UINT16_MAX);
    ASSERT(ws_scratch->alt_peptide_cap >= ws_scratch->alt_cds_cap / 3u + 1u);
    ASSERT_EQ(DUCKVEP_VARIANT_CODING_CONTEXT_OK,
              duckvep_variant_coding_context_build(&s.tx, &s.ex, &s.seq, &s.v,
                                                   0u, 0u, s.strand,
                                                   ws_scratch->edits,
                                                   ws_scratch->edits_cap,
                                                   ws_scratch->alt_cds,
                                                   ws_scratch->alt_cds_cap,
                                                   ws_scratch->ref_peptide,
                                                   ws_scratch->ref_peptide_cap,
                                                   ws_scratch->alt_peptide,
                                                   ws_scratch->alt_peptide_cap,
                                                   &ctx));
    ASSERT_EQ(3, ctx.length_diff);
    ASSERT_EQ((size_t)18u, ctx.alt_cds_len);
    ASSERT_EQ((size_t)5u, ctx.ref_peptide_len);
    ASSERT_EQ((size_t)6u, ctx.alt_peptide_len);
    ASSERT_EQ(DUCKVEP_CONTEXT_DELTA_OK,
              duckvep_coding_context_delta_fill(&ctx, 0u, &delta));
    ASSERT(kprop_delta_is_inframe_insertion_at(&delta, 3));
    duckvep_workspace_close(ws);
    duckvep_model_close(model);
    PASS();
}

TEST workspace_delta_scratch_usable_for_mnv(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'G','A','A',  'G','A','A',  'G','A','A',  'T','T','T'
    };
    struct kprop_coding s;
    duckvep_model_t *model = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_delta_scratch_t *ws_scratch;
    duckvep_haplotype_edit_t stack_edits[4];
    uint8_t stack_alt_cds[32];
    uint8_t stack_ref_pep[16];
    uint8_t stack_alt_pep[16];
    duckvep_delta_scratch_t stack_scratch;
    duckvep_sequence_delta_t got;
    duckvep_sequence_delta_t want;
    duckvep_error_t err;
    uint8_t alt_tx[2] = {'C','C'};
    uint32_t i;

    memset(&s, 0, sizeof s); memset(&err, 0, sizeof err);
    s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
    s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, 15u);
    s.vpos = kprop_genomic_pos_for_cds(&s, 4u); s.vend = s.vpos + 1u;
    s.vkind = (uint8_t)DUCKVEP_KIND_MNV;
    for (i = 0u; i < 2u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
    kprop_fill_variant_alt_from_tx(&s, 2u, alt_tx, 2u);
    s.roff = 0u; s.aoff = 2u; s.rlen = 2u; s.alen = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&s.tx, &s.ex, &s.seq, &model, &err));
    ASSERT_EQ(DUCKVEP_OK, duckvep_workspace_open(model, &ws, &err));
    ws_scratch = duckvep_workspace_delta_scratch(ws);
    ASSERT(ws_scratch != NULL);
    ASSERT(ws_scratch->edits_cap >= 1u);
    ASSERT(ws_scratch->alt_cds_cap >= s.cds_lenv);
    ASSERT(ws_scratch->ref_peptide_cap >= s.cds_lenv / 3u + 1u);
    ASSERT(ws_scratch->alt_peptide_cap >= s.cds_lenv / 3u + 1u);

    memset(&stack_scratch, 0, sizeof stack_scratch);
    stack_scratch.edits = stack_edits; stack_scratch.edits_cap = 4u;
    stack_scratch.alt_cds = stack_alt_cds; stack_scratch.alt_cds_cap = sizeof stack_alt_cds;
    stack_scratch.ref_peptide = stack_ref_pep; stack_scratch.ref_peptide_cap = sizeof stack_ref_pep;
    stack_scratch.alt_peptide = stack_alt_pep; stack_scratch.alt_peptide_cap = sizeof stack_alt_pep;

    duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_MNV, &s.tx, &s.ex, &s.seq,
                                             &s.v, 0u, 0u, s.vpos, s.strand,
                                             &stack_scratch, &want);
    duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_MNV, &s.tx, &s.ex, &s.seq,
                                             &s.v, 0u, 0u, s.vpos, s.strand,
                                             ws_scratch, &got);
    ASSERT(want.valid);
    ASSERT(got.valid);
    ASSERT_EQ(want.synonymous, got.synonymous);
    ASSERT_EQ(want.missense, got.missense);
    ASSERT_EQ(want.stop_gained, got.stop_gained);
    ASSERT_EQ(want.stop_lost, got.stop_lost);
    ASSERT_EQ(want.stop_retained, got.stop_retained);
    ASSERT_EQ(want.start_lost, got.start_lost);
    ASSERT_EQ(want.start_retained, got.start_retained);
    ASSERT_EQ(want.frameshift, got.frameshift);
    ASSERT_EQ(want.inframe_deletion, got.inframe_deletion);
    ASSERT_EQ(want.inframe_insertion, got.inframe_insertion);
    ASSERT_EQ(want.protein_altering, got.protein_altering);
    ASSERT_EQ(want.protein_pos, got.protein_pos);
    ASSERT_EQ(want.cdna_pos, got.cdna_pos);
    ASSERT_EQ(want.cds_pos, got.cds_pos);
    ASSERT_EQ(want.ref_aa, got.ref_aa);
    ASSERT_EQ(want.alt_aa, got.alt_aa);
    duckvep_workspace_close(ws);
    duckvep_model_close(model);
    PASS();
}

TEST sequence_delta_scratch_rejects_unequal_mnv_kind(void) {
    static uint8_t cds[9] = {
        'A','T','G',  'A','A','A',  'T','T','T'
    };
    struct kprop_coding s;
    duckvep_haplotype_edit_t stack_edits[4];
    uint8_t stack_alt_cds[32];
    uint8_t stack_ref_pep[16];
    uint8_t stack_alt_pep[16];
    duckvep_delta_scratch_t scratch;
    duckvep_sequence_delta_t delta;
    uint8_t alt_tx[1] = { 'G' };
    uint32_t i;

    memset(&s, 0, sizeof s);
    s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1008u; s.cds_s = 1000u; s.cds_e = 1008u;
    s.es = 1000u; s.ee = 1008u; s.ecds = 1u; s.ecde = 9u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, 9u);
    s.vpos = kprop_genomic_pos_for_cds(&s, 3u); s.vend = s.vpos + 3u;
    s.vkind = (uint8_t)DUCKVEP_KIND_MNV;
    for (i = 0u; i < 4u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
    kprop_fill_variant_alt_from_tx(&s, 4u, alt_tx, 1u);
    s.roff = 0u; s.aoff = 4u; s.rlen = 4u; s.alen = 1u;

    memset(&scratch, 0, sizeof scratch);
    scratch.edits = stack_edits; scratch.edits_cap = 4u;
    scratch.alt_cds = stack_alt_cds; scratch.alt_cds_cap = sizeof stack_alt_cds;
    scratch.ref_peptide = stack_ref_pep; scratch.ref_peptide_cap = sizeof stack_ref_pep;
    scratch.alt_peptide = stack_alt_pep; scratch.alt_peptide_cap = sizeof stack_alt_pep;

    duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_MNV, &s.tx, &s.ex, &s.seq,
                                             &s.v, 0u, 0u, s.vpos, s.strand,
                                             &scratch, &delta);
    ASSERT(!delta.valid);
    PASS();
}

TEST sequence_delta_annotation_wrapper_del_insufficient_scratch_known(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'A','A','A',  'C','C','C',  'G','G','G',  'T','T','T'
    };
    struct kprop_coding s;
    duckvep_haplotype_edit_t stack_edits[1];
    uint8_t stack_alt_cds[32];
    uint8_t stack_ref_pep[16];
    uint8_t stack_alt_pep[16];
    duckvep_delta_scratch_t scratch;
    duckvep_sequence_delta_t shape;
    duckvep_sequence_delta_t routed;
    duckvep_sequence_delta_route_t route;
    uint32_t i;

    memset(&s, 0, sizeof s);
    s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
    s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, 15u);
    s.vpos = kprop_genomic_pos_for_cds(&s, 3u); s.vend = s.vpos + 3u;
    s.vkind = (uint8_t)DUCKVEP_KIND_DEL;
    for (i = 0u; i < 4u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
    s.abytes[4] = s.abytes[0];
    s.roff = 0u; s.aoff = 4u; s.rlen = 4u; s.alen = 1u;

    memset(&scratch, 0, sizeof scratch);
    scratch.edits = stack_edits; scratch.edits_cap = 0u;
    scratch.alt_cds = stack_alt_cds; scratch.alt_cds_cap = sizeof stack_alt_cds;
    scratch.ref_peptide = stack_ref_pep; scratch.ref_peptide_cap = sizeof stack_ref_pep;
    scratch.alt_peptide = stack_alt_pep; scratch.alt_peptide_cap = sizeof stack_alt_pep;

    duckvep_sequence_delta_fill(DUCKVEP_KIND_DEL, &s.tx, &s.ex, &s.seq,
                                &s.v, 0u, 0u, s.vpos, s.strand, &shape);
    duckvep_sequence_delta_fill_for_annotation_trace(DUCKVEP_KIND_DEL, &s.tx, &s.ex,
                                                     &s.seq, &s.v, 0u, 0u, s.vpos,
                                                     s.strand, &scratch, NULL, &route,
                                                     &routed);
    ASSERT(kprop_delta_is_inframe_deletion_at(&shape, 2));
    ASSERT(kprop_sequence_delta_equal(&shape, &routed));
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_DEL_DIRECT_FALLBACK, route);
    PASS();
}

TEST sequence_delta_annotation_wrapper_start_lost_mnv(void) {
    static uint8_t cds[9] = {'A','T','G', 'G','A','A', 'T','A','A'};
    struct kprop_coding s;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[32];
    uint8_t ref_pep[16];
    uint8_t alt_pep[16];
    duckvep_delta_scratch_t scratch;
    duckvep_sequence_delta_t shape;
    duckvep_sequence_delta_t direct;
    duckvep_sequence_delta_t routed;
    uint8_t alt_tx[2] = {'C','C'};
    uint32_t i;

    memset(&s, 0, sizeof s);
    s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1008u; s.cds_s = 1000u; s.cds_e = 1008u;
    s.es = 1000u; s.ee = 1008u; s.ecds = 1u; s.ecde = 9u; s.eph = 0; s.eeph = 0;
    s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, 9u);
    s.vpos = kprop_genomic_pos_for_cds(&s, 1u); s.vend = s.vpos + 1u;
    s.vkind = (uint8_t)DUCKVEP_KIND_MNV;
    for (i = 0u; i < 2u; i++) s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
    kprop_fill_variant_alt_from_tx(&s, 2u, alt_tx, 2u);
    s.roff = 0u; s.aoff = 2u; s.rlen = 2u; s.alen = 2u;

    scratch.edits = edits; scratch.edits_cap = 4u;
    scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
    scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
    scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;

    duckvep_sequence_delta_fill(DUCKVEP_KIND_MNV, &s.tx, &s.ex, &s.seq, &s.v,
                                0u, 0u, s.vpos, s.strand, &shape);
    duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_MNV, &s.tx, &s.ex, &s.seq,
                                             &s.v, 0u, 0u, s.vpos, s.strand,
                                             &scratch, &direct);
    duckvep_sequence_delta_fill_for_annotation(DUCKVEP_KIND_MNV, &s.tx, &s.ex,
                                               &s.seq, &s.v, 0u, 0u, s.vpos,
                                               s.strand, &scratch, NULL, &routed);
    ASSERT(shape.valid);
    ASSERT(direct.valid);
    ASSERT(routed.valid);
    ASSERT(kprop_sequence_delta_equal(&shape, &direct));
    ASSERT(kprop_sequence_delta_equal(&shape, &routed));
    ASSERT(routed.start_lost);
    ASSERT(!routed.missense);
    ASSERT_EQ((uint8_t)'M', routed.ref_aa);
    ASSERT_EQ(1, routed.protein_pos);

    s.flags = (uint64_t)DUCKVEP_TX_CDS_START_NF;
    kprop_wire_coding_scene(&s, 9u);
    duckvep_sequence_delta_fill_for_annotation(DUCKVEP_KIND_MNV, &s.tx, &s.ex,
                                               &s.seq, &s.v, 0u, 0u, s.vpos,
                                               s.strand, &scratch, NULL, &routed);
    ASSERT(routed.valid);
    ASSERT(!routed.start_lost);
    ASSERT(routed.missense);
    PASS();
}

static int kprop_cross_codon_scene_deltas(
    int8_t                         strand,
    uint8_t                       *cds,
    uint32_t                       cds_len,
    uint32_t                       first_cds,
    uint32_t                       ref_len,
    const uint8_t                 *alt_tx,
    uint32_t                       alt_len,
    duckvep_sequence_delta_t      *direct,
    duckvep_sequence_delta_t      *shape,
    duckvep_sequence_delta_t      *routed,
    duckvep_sequence_delta_route_t *route) {

    struct kprop_coding s;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[64];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    duckvep_delta_scratch_t scratch;
    uint32_t i;

    if (cds == NULL || alt_tx == NULL || direct == NULL || shape == NULL || routed == NULL ||
        route == NULL || ref_len == 0u || ref_len + alt_len > sizeof s.abytes) {
        return 0;
    }
    memset(&s, 0, sizeof s);
    s.cds = cds; s.chrom = 0u; s.strand = strand; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1000u + cds_len - 1u;
    s.cds_s = s.tstart; s.cds_e = s.tend;
    s.es = s.tstart; s.ee = s.tend; s.ecds = 1u; s.ecde = cds_len;
    s.eph = 0; s.eeph = 0; s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    kprop_wire_coding_scene(&s, cds_len);
    s.vpos = strand > 0 ? kprop_genomic_pos_for_cds(&s, first_cds)
                         : kprop_genomic_pos_for_cds(&s, first_cds + ref_len - 1u);
    s.vend = s.vpos + ref_len - 1u;
    s.vkind = (uint8_t)DUCKVEP_KIND_MNV;
    for (i = 0u; i < ref_len; i++) {
        s.abytes[i] = (uint8_t)kprop_genomic_base_at(&s, s.vpos + i);
    }
    kprop_fill_variant_alt_from_tx(&s, ref_len, alt_tx, alt_len);
    s.roff = 0u; s.aoff = ref_len; s.rlen = (uint16_t)ref_len; s.alen = (uint16_t)alt_len;

    scratch.edits = edits; scratch.edits_cap = 4u;
    scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
    scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
    scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;

    duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_MNV, &s.tx, &s.ex, &s.seq,
                                             &s.v, 0u, 0u, s.vpos, s.strand,
                                             &scratch, direct);
    duckvep_sequence_delta_fill(DUCKVEP_KIND_MNV, &s.tx, &s.ex, &s.seq, &s.v,
                                0u, 0u, s.vpos, s.strand, shape);
    duckvep_sequence_delta_fill_for_annotation_trace(DUCKVEP_KIND_MNV, &s.tx, &s.ex,
                                                     &s.seq, &s.v, 0u, 0u, s.vpos,
                                                     s.strand, &scratch, NULL, route,
                                                     routed);
    return 1;
}

TEST sequence_delta_with_scratch_cross_codon_known_scene(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'G','A','A',  'G','A','A',  'G','A','A',  'T','T','T'
    };
    duckvep_sequence_delta_t direct;
    duckvep_sequence_delta_t shape;
    duckvep_sequence_delta_t routed;
    duckvep_sequence_delta_route_t route;
    uint8_t alt_tx[2] = {'C','C'};

    ASSERT(kprop_cross_codon_scene_deltas(1, cds, 15u, 6u, 2u, alt_tx, 2u,
                                          &direct, &shape, &routed, &route));
    ASSERT(kprop_delta_is_coarse_cross_codon_missense(&direct));
    ASSERT(kprop_delta_is_coarse_cross_codon_missense(&shape));
    ASSERT(kprop_sequence_delta_equal(&shape, &direct));
    ASSERT(kprop_sequence_delta_equal(&shape, &routed));
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);
    PASS();
}

TEST sequence_delta_with_scratch_cross_codon_reverse_known_scene(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'G','A','A',  'G','A','A',  'G','A','A',  'T','T','T'
    };
    duckvep_sequence_delta_t direct;
    duckvep_sequence_delta_t shape;
    duckvep_sequence_delta_t routed;
    duckvep_sequence_delta_route_t route;
    uint8_t alt_tx[2] = {'C','C'};

    ASSERT(kprop_cross_codon_scene_deltas(-1, cds, 15u, 6u, 2u, alt_tx, 2u,
                                          &direct, &shape, &routed, &route));
    ASSERT(kprop_delta_is_coarse_cross_codon_missense(&direct));
    ASSERT(kprop_delta_is_coarse_cross_codon_missense(&shape));
    ASSERT(kprop_sequence_delta_equal(&shape, &direct));
    ASSERT(kprop_sequence_delta_equal(&shape, &routed));
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);
    PASS();
}

TEST sequence_delta_with_scratch_cross_codon_negative_scenes(void) {
    static uint8_t syn_cds[15] = {
        'A','T','G',  'A','A','A',  'T','T','A',  'G','A','A',  'T','T','T'
    };
    static uint8_t stop_cds[15] = {
        'A','T','G',  'T','C','A',  'A','A','A',  'G','A','A',  'T','T','T'
    };
    static uint8_t terminal_cds[9] = {
        'A','T','G',  'G','A','A',  'G','A','A'
    };
    static uint8_t wide_cds[15] = {
        'A','T','G',  'G','A','A',  'G','A','A',  'G','A','A',  'T','T','T'
    };
    static uint8_t length_cds[15] = {
        'A','T','G',  'G','A','A',  'G','A','A',  'G','A','A',  'T','T','T'
    };
    duckvep_sequence_delta_t direct;
    duckvep_sequence_delta_t shape;
    duckvep_sequence_delta_t routed;
    duckvep_sequence_delta_route_t route;
    uint8_t syn_alt[2] = {'G','C'};
    uint8_t stop_alt[3] = {'A','G','C'};
    uint8_t terminal_alt[2] = {'C','C'};
    uint8_t wide_alt[4] = {'C','C','C','C'};
    uint8_t length_alt[1] = {'C'};

    /* These two-codon windows were unsupported by the old two-codon-missense-only slice and
     * fell back to coding_sequence_variant. The generalized window classifier now resolves
     * each one authoritatively (synonymous / stop_gained / missense, incl. the terminal-codon
     * and >3-base-wide windows). The direct shape-specific filler still only knows the missense
     * cross-codon window, so it stays invalid on the synonymous and stop-gained cases — which
     * is exactly the capability the interpreter adds. The router echoes the interpreter. */

    /* both codons synonymous -> synonymous_variant */
    ASSERT(kprop_cross_codon_scene_deltas(1, syn_cds, 15u, 6u, 2u, syn_alt, 2u,
                                          &direct, &shape, &routed, &route));
    ASSERT(direct.valid && direct.synonymous && direct.protein_pos == -1 &&
           !direct.missense && !direct.stop_gained && !direct.stop_lost &&
           !direct.stop_retained && !direct.start_lost && !direct.start_retained);
    ASSERT(!shape.valid);
    ASSERT(kprop_sequence_delta_equal(&direct, &routed));
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);

    /* junction codon becomes a stop -> stop_gained */
    ASSERT(kprop_cross_codon_scene_deltas(1, stop_cds, 15u, 5u, 3u, stop_alt, 3u,
                                          &direct, &shape, &routed, &route));
    ASSERT(direct.valid && direct.stop_gained && direct.protein_pos == -1 &&
           !direct.missense && !direct.synonymous && !direct.stop_lost &&
           !direct.stop_retained && !direct.start_lost && !direct.start_retained);
    ASSERT(!shape.valid);
    ASSERT(kprop_sequence_delta_equal(&direct, &routed));
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);

    /* window includes the terminal (last) codon -> missense */
    ASSERT(kprop_cross_codon_scene_deltas(1, terminal_cds, 9u, 6u, 2u, terminal_alt, 2u,
                                          &direct, &shape, &routed, &route));
    ASSERT(direct.valid && direct.missense && direct.protein_pos == -1);
    ASSERT(!shape.valid);
    ASSERT(kprop_sequence_delta_equal(&direct, &routed));
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);

    /* four-base window spanning two codons -> missense */
    ASSERT(kprop_cross_codon_scene_deltas(1, wide_cds, 15u, 5u, 4u, wide_alt, 4u,
                                          &direct, &shape, &routed, &route));
    ASSERT(direct.valid && direct.missense && direct.protein_pos == -1);
    ASSERT(!shape.valid);
    ASSERT(kprop_sequence_delta_equal(&direct, &routed));
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);

    /* unequal ref/alt length on a KIND_MNV is a disguised delins: genuinely unsupported. */
    ASSERT(kprop_cross_codon_scene_deltas(1, length_cds, 15u, 6u, 2u, length_alt, 1u,
                                          &direct, &shape, &routed, &route));
    ASSERT(!direct.valid);
    ASSERT(!shape.valid);
    ASSERT(!routed.valid);
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_MNV_DIRECT_FALLBACK, route);
    PASS();
}

TEST feature_substitution_window_fails_closed_known_scene(void) {
    static uint8_t cds[15] = {
        'A','T','G',  'N','A','A',  'C','C','C',  'T','G','G',  'T','A','A'
    };
    struct kprop_coding s;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[64];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    duckvep_delta_scratch_t scratch;
    duckvep_event_t event;
    duckvep_sequence_delta_t delta;
    duckvep_sequence_delta_route_t route;

    memset(&s, 0, sizeof s);
    s.cds = cds; s.chrom = 0u; s.strand = 1; s.flags = 0u;
    s.tstart = 1000u; s.tend = 1014u; s.cds_s = 1000u; s.cds_e = 1014u;
    s.es = 1000u; s.ee = 1014u; s.ecds = 1u; s.ecde = 15u;
    s.eph = 0; s.eeph = 0; s.exoff = 0u; s.excnt = 1u; s.vchrom = 0u;
    s.vkind = (uint8_t)DUCKVEP_KIND_SNV;
    s.roff = 0u; s.aoff = 2u; s.rlen = 2u; s.alen = 2u;
    kprop_wire_coding_scene(&s, sizeof cds);

    scratch.edits = edits; scratch.edits_cap = 4u;
    scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
    scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
    scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;

    /* The semantic edit is G>A at CDS position 3. The complete uploaded feature
     * also retains an N in the next codon, so VEP cannot use the widened peptide
     * window. Do not retry the apparently valid one-base edit. */
    s.vpos = 1002u; s.vend = 1003u;
    s.abytes[0] = 'G'; s.abytes[1] = 'N';
    s.abytes[2] = 'A'; s.abytes[3] = 'N';
    duckvep_event_load(&s.v, 0u, &event);
    duckvep_sequence_delta_fill_for_annotation_trace(
        DUCKVEP_KIND_SNV, &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u,
        event.start1, s.strand, &scratch, &event, &route, &delta);
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);
    ASSERT(!delta.valid);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_AMBIGUOUS, delta.sequence_status);

    /* The retained A is not the transcript's G at CDS position 12. The changed
     * T>C base itself is valid, so validating only the trimmed edit would miss
     * this reference mismatch and emit a supported terminal-stop consequence. */
    s.vpos = 1011u; s.vend = 1012u;
    s.abytes[0] = 'A'; s.abytes[1] = 'T';
    s.abytes[2] = 'A'; s.abytes[3] = 'C';
    duckvep_event_load(&s.v, 0u, &event);
    duckvep_sequence_delta_fill_for_annotation_trace(
        DUCKVEP_KIND_SNV, &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u,
        event.start1, s.strand, &scratch, &event, &route, &delta);
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);
    ASSERT(!delta.valid);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_REFERENCE_MISMATCH,
              delta.sequence_status);

    /* On the reverse strand genomic 1002-1003 maps to CDS positions 13-12.
     * The same checks must follow genomic byte order while comparing against
     * the reverse-complemented transcript bases. */
    s.strand = -1;
    kprop_wire_coding_scene(&s, sizeof cds);
    s.vpos = 1002u; s.vend = 1003u;
    s.abytes[0] = 'G'; s.abytes[1] = 'C';
    s.abytes[2] = 'G'; s.abytes[3] = 'T';
    duckvep_event_load(&s.v, 0u, &event);
    duckvep_sequence_delta_fill_for_annotation_trace(
        DUCKVEP_KIND_SNV, &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u,
        event.start1, s.strand, &scratch, &event, &route, &delta);
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);
    ASSERT(!delta.valid);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_REFERENCE_MISMATCH,
              delta.sequence_status);

    s.abytes[0] = 'N'; s.abytes[1] = 'C';
    s.abytes[2] = 'N'; s.abytes[3] = 'T';
    duckvep_event_load(&s.v, 0u, &event);
    duckvep_sequence_delta_fill_for_annotation_trace(
        DUCKVEP_KIND_SNV, &s.tx, &s.ex, &s.seq, &s.v, 0u, 0u,
        event.start1, s.strand, &scratch, &event, &route, &delta);
    ASSERT_EQ(DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT, route);
    ASSERT(!delta.valid);
    ASSERT_EQ((uint8_t)DUCKVEP_SEQUENCE_AMBIGUOUS, delta.sequence_status);
    PASS();
}

static struct {
    uint32_t syn;
    uint32_t mis;
    uint32_t stop_gained;
    uint32_t stop_lost;
    uint32_t stop_retained;
    uint32_t fwd;
    uint32_t rev;
} g_delta_wrapper_cov;

static enum theft_trial_res prop_delta_annotation_wrapper_matches_direct_shape(struct theft *t,
                                                                         void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[64];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    duckvep_delta_scratch_t scratch;
    duckvep_sequence_delta_t shape;
    duckvep_sequence_delta_t routed;
    (void)t;

    scratch.edits = edits; scratch.edits_cap = 4u;
    scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
    scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
    scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;

    duckvep_sequence_delta_fill((duckvep_variant_kind_t)s->vkind, &s->tx, &s->ex,
                                &s->seq, &s->v, 0u, 0u, s->vpos, s->strand,
                                &shape);
    duckvep_sequence_delta_fill_for_annotation((duckvep_variant_kind_t)s->vkind,
                                               &s->tx, &s->ex, &s->seq, &s->v,
                                               0u, 0u, s->vpos, s->strand,
                                               &scratch, NULL, &routed);
    if (!shape.valid || !routed.valid) return THEFT_TRIAL_FAIL;
    if (!kprop_sequence_delta_equal(&shape, &routed)) return THEFT_TRIAL_FAIL;

    if (s->expect_region == KPROP_CONTEXT_DELTA_SYNONYMOUS) g_delta_wrapper_cov.syn++;
    else if (s->expect_region == KPROP_CONTEXT_DELTA_MISSENSE) g_delta_wrapper_cov.mis++;
    else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_GAINED) g_delta_wrapper_cov.stop_gained++;
    else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_LOST) g_delta_wrapper_cov.stop_lost++;
    else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_RETAINED) g_delta_wrapper_cov.stop_retained++;
    else return THEFT_TRIAL_FAIL;
    if (s->strand > 0) g_delta_wrapper_cov.fwd++;
    else g_delta_wrapper_cov.rev++;
    return THEFT_TRIAL_PASS;
}

TEST sequence_delta_annotation_wrapper_matches_direct_shape(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "sequence delta annotation wrapper MNV == direct shape";
    cfg.prop1 = prop_delta_annotation_wrapper_matches_direct_shape;
    cfg.type_info[0] = &kprop_context_delta_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_delta_wrapper_cov, 0, sizeof g_delta_wrapper_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_delta_wrapper_cov.syn > 0u);
    ASSERT(g_delta_wrapper_cov.mis > 0u);
    ASSERT(g_delta_wrapper_cov.stop_gained > 0u);
    ASSERT(g_delta_wrapper_cov.stop_lost > 0u);
    ASSERT(g_delta_wrapper_cov.stop_retained > 0u);
    ASSERT(g_delta_wrapper_cov.fwd > 0u);
    ASSERT(g_delta_wrapper_cov.rev > 0u);
    fprintf(stderr,
            "[delta-wrapper coverage] syn=%u mis=%u stop_gained=%u stop_lost=%u stop_retained=%u fwd=%u rev=%u\n",
            g_delta_wrapper_cov.syn, g_delta_wrapper_cov.mis,
            g_delta_wrapper_cov.stop_gained, g_delta_wrapper_cov.stop_lost,
            g_delta_wrapper_cov.stop_retained, g_delta_wrapper_cov.fwd,
            g_delta_wrapper_cov.rev);
    PASS();
}

static struct {
    uint32_t syn;
    uint32_t mis;
    uint32_t stop_gained;
    uint32_t stop_lost;
    uint32_t stop_retained;
    uint32_t fwd;
    uint32_t rev;
    uint32_t full;
} g_cursor_route_cov;

static enum theft_trial_res prop_cursor_padded_snv_matches_tile(struct theft *t,
                                                                void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_variant_batch_t variants = s->v;
    const uint8_t semantic_kind = (uint8_t)DUCKVEP_KIND_SNV;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_annotate_cursor_t *cur = NULL;
    duckvep_error_t err;
    duckvep_consequence_t tile_rows[4];
    duckvep_consequence_t cursor_rows[4];
    duckvep_consequence_t chunk[1];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    duckvep_workspace_delta_route_stats_t tile_stats;
    duckvep_workspace_delta_route_stats_t cursor_stats;
    size_t tile_n;
    size_t cursor_n = 0u;
    int saw_full = 0;
    enum theft_trial_res res = THEFT_TRIAL_PASS;
    (void)t;

    memset(&err, 0, sizeof err);
    memset(&tile_stats, 0, sizeof tile_stats);
    memset(&cursor_stats, 0, sizeof cursor_stats);
    /* The corpus stores whole codons whose differing region is one base. At the
     * public boundary these are SNVs, regardless of their padded representation. */
    variants.variant_kind = &semantic_kind;
    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, &model, &err) != DUCKVEP_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }

    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, tile_rows, 4u);
    if (duckvep_annotate_tile(model, &variants, opts, ws, &rb, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    tile_n = duckvep_result_builder_count(&rb);
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { res = THEFT_TRIAL_FAIL; goto done; }
    tile_stats = *stats;

    duckvep_workspace_delta_route_stats_reset(ws);
    if (duckvep_annotate_cursor_open(model, &variants, opts, ws, &cur, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    while (!duckvep_annotate_cursor_done(cur)) {
        duckvep_status_t st;
        size_t i;
        duckvep_result_builder_init(&rb, chunk, 1u);
        st = duckvep_annotate_cursor_fill(cur, &rb, &err);
        if (st != DUCKVEP_OK && st != DUCKVEP_ERR_RESULT_FULL) { res = THEFT_TRIAL_FAIL; goto done; }
        if (st == DUCKVEP_ERR_RESULT_FULL) saw_full = 1;
        for (i = 0u; i < duckvep_result_builder_count(&rb); i++) {
            if (cursor_n >= 4u) { res = THEFT_TRIAL_FAIL; goto done; }
            cursor_rows[cursor_n++] = chunk[i];
        }
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { res = THEFT_TRIAL_FAIL; goto done; }
    cursor_stats = *stats;

    if (!saw_full || tile_n != cursor_n || tile_n == 0u) { res = THEFT_TRIAL_FAIL; goto done; }
    if (tile_stats.substitution_context != 1u ||
        tile_stats.mnv_direct_fallback != 0u ||
        tile_stats.substitution_context != cursor_stats.substitution_context ||
        tile_stats.mnv_direct_fallback != cursor_stats.mnv_direct_fallback) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    {
        size_t i;
        for (i = 0u; i < tile_n; i++) {
            if (!consequence_rows_equal(&tile_rows[i], &cursor_rows[i])) {
                res = THEFT_TRIAL_FAIL; goto done;
            }
        }
    }

    if (s->expect_region == KPROP_CONTEXT_DELTA_SYNONYMOUS) g_cursor_route_cov.syn++;
    else if (s->expect_region == KPROP_CONTEXT_DELTA_MISSENSE) g_cursor_route_cov.mis++;
    else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_GAINED) g_cursor_route_cov.stop_gained++;
    else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_LOST) g_cursor_route_cov.stop_lost++;
    else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_RETAINED) g_cursor_route_cov.stop_retained++;
    else { res = THEFT_TRIAL_FAIL; goto done; }
    if (s->strand > 0) g_cursor_route_cov.fwd++;
    else g_cursor_route_cov.rev++;
    if (saw_full) g_cursor_route_cov.full++;

done:
    duckvep_annotate_cursor_close(cur);
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return res;
}

TEST annotate_cursor_padded_snv_matches_tile_for_any_output_split(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate cursor padded SNV == tile under output splits";
    cfg.prop1 = prop_cursor_padded_snv_matches_tile;
    cfg.type_info[0] = &kprop_context_delta_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_cursor_route_cov, 0, sizeof g_cursor_route_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_cursor_route_cov.syn > 0u);
    ASSERT(g_cursor_route_cov.mis > 0u);
    ASSERT(g_cursor_route_cov.stop_gained > 0u);
    ASSERT(g_cursor_route_cov.stop_lost > 0u);
    ASSERT(g_cursor_route_cov.stop_retained > 0u);
    ASSERT(g_cursor_route_cov.fwd > 0u);
    ASSERT(g_cursor_route_cov.rev > 0u);
    ASSERT(g_cursor_route_cov.full > 0u);
    fprintf(stderr,
            "[cursor-route coverage] syn=%u mis=%u stop_gained=%u stop_lost=%u stop_retained=%u fwd=%u rev=%u full=%u\n",
            g_cursor_route_cov.syn, g_cursor_route_cov.mis,
            g_cursor_route_cov.stop_gained, g_cursor_route_cov.stop_lost,
            g_cursor_route_cov.stop_retained, g_cursor_route_cov.fwd,
            g_cursor_route_cov.rev, g_cursor_route_cov.full);
    PASS();
}

static struct {
    uint32_t context;
    uint32_t fwd;
    uint32_t rev;
    uint32_t len2;
    uint32_t len3;
} g_cursor_cross_route_cov;

static enum theft_trial_res prop_cursor_cross_codon_mnv_route_matches_tile(
    struct theft *t, void *arg1) {

    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_annotate_cursor_t *cur = NULL;
    duckvep_error_t err;
    duckvep_consequence_t tile_rows[4];
    duckvep_consequence_t cursor_rows[4];
    duckvep_consequence_t chunk[1];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    duckvep_workspace_delta_route_stats_t tile_stats;
    duckvep_workspace_delta_route_stats_t cursor_stats;
    size_t tile_n;
    size_t cursor_n = 0u;
    enum theft_trial_res res = THEFT_TRIAL_PASS;
    (void)t;

    memset(&err, 0, sizeof err);
    memset(&tile_stats, 0, sizeof tile_stats);
    memset(&cursor_stats, 0, sizeof cursor_stats);
    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, &model, &err) != DUCKVEP_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }

    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, tile_rows, 4u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    tile_n = duckvep_result_builder_count(&rb);
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { res = THEFT_TRIAL_FAIL; goto done; }
    tile_stats = *stats;

    duckvep_workspace_delta_route_stats_reset(ws);
    if (duckvep_annotate_cursor_open(model, &s->v, opts, ws, &cur, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    while (!duckvep_annotate_cursor_done(cur)) {
        duckvep_status_t st;
        size_t i;
        duckvep_result_builder_init(&rb, chunk, 1u);
        st = duckvep_annotate_cursor_fill(cur, &rb, &err);
        if (st != DUCKVEP_OK && st != DUCKVEP_ERR_RESULT_FULL) { res = THEFT_TRIAL_FAIL; goto done; }
        for (i = 0u; i < duckvep_result_builder_count(&rb); i++) {
            if (cursor_n >= 4u) { res = THEFT_TRIAL_FAIL; goto done; }
            cursor_rows[cursor_n++] = chunk[i];
        }
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { res = THEFT_TRIAL_FAIL; goto done; }
    cursor_stats = *stats;

    if (tile_n != cursor_n || tile_n != 1u) { res = THEFT_TRIAL_FAIL; goto done; }
    if (!consequence_rows_equal(&tile_rows[0], &cursor_rows[0])) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    if (tile_stats.substitution_context != cursor_stats.substitution_context ||
        tile_stats.mnv_direct_fallback != cursor_stats.mnv_direct_fallback) {
        res = THEFT_TRIAL_FAIL; goto done;
    }

    /* Every generator mode (missense / synonymous / stop-gained cross-codon window)
     * resolves through the context interpreter. The
     * property under test is that a chunked cursor and a single tile agree on both the row
     * and the route stats for any output split — verified above. */
    if (tile_stats.substitution_context != 1u ||
        tile_stats.mnv_direct_fallback != 0u) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    g_cursor_cross_route_cov.context++;
    if (s->strand > 0) g_cursor_cross_route_cov.fwd++; else g_cursor_cross_route_cov.rev++;
    if (s->rlen == 2u) g_cursor_cross_route_cov.len2++;
    else if (s->rlen == 3u) g_cursor_cross_route_cov.len3++;
    else { res = THEFT_TRIAL_FAIL; goto done; }

done:
    duckvep_annotate_cursor_close(cur);
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return res;
}

TEST annotate_cursor_cross_codon_mnv_route_matches_tile_for_any_output_split(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate cursor cross-codon MNV route == tile";
    cfg.prop1 = prop_cursor_cross_codon_mnv_route_matches_tile;
    cfg.type_info[0] = &kprop_cross_codon_mnv_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_cursor_cross_route_cov, 0, sizeof g_cursor_cross_route_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_cursor_cross_route_cov.context > 0u);
    ASSERT(g_cursor_cross_route_cov.fwd > 0u);
    ASSERT(g_cursor_cross_route_cov.rev > 0u);
    ASSERT(g_cursor_cross_route_cov.len2 > 0u);
    ASSERT(g_cursor_cross_route_cov.len3 > 0u);
    fprintf(stderr,
            "[cursor-cross-route coverage] context=%u fwd=%u rev=%u len2=%u len3=%u\n",
            g_cursor_cross_route_cov.context,
            g_cursor_cross_route_cov.fwd, g_cursor_cross_route_cov.rev,
            g_cursor_cross_route_cov.len2, g_cursor_cross_route_cov.len3);
    PASS();
}

static struct {
    uint32_t syn;
    uint32_t mis;
    uint32_t stop_gained;
    uint32_t stop_lost;
    uint32_t stop_retained;
    uint32_t fwd;
    uint32_t rev;
    uint32_t capfail;
} g_delta_scratch_cov;

static enum theft_trial_res prop_delta_scratch_mnv_matches_codon_oracle(struct theft *t,
                                                                        void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[64];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    duckvep_delta_scratch_t scratch;
    duckvep_delta_scratch_t too_small;
    duckvep_sequence_delta_t delta;
    size_t first_diff = 0u;
    size_t codon_start;
    size_t codon_idx;
    char ref_codon[4];
    char alt_codon[4];
    duckvep_codon_result_t cr;
    uint32_t j;
    (void)t;

    scratch.edits = edits; scratch.edits_cap = 4u;
    scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
    scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
    scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;

    too_small = scratch;
    too_small.edits_cap = 0u;
    duckvep_sequence_delta_fill_with_scratch((duckvep_variant_kind_t)s->vkind,
                                             &s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->vpos, s->strand,
                                             &too_small, &delta);
    if (delta.valid) return THEFT_TRIAL_FAIL;
    g_delta_scratch_cov.capfail++;

    duckvep_sequence_delta_fill_with_scratch((duckvep_variant_kind_t)s->vkind,
                                             &s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->vpos, s->strand,
                                             &scratch, &delta);
    if (!delta.valid) return THEFT_TRIAL_FAIL;

    while (first_diff < s->cds_lenv && s->cds[first_diff] == s->expect_cds[first_diff]) {
        first_diff++;
    }
    if (first_diff >= s->cds_lenv) return THEFT_TRIAL_FAIL;
    codon_start = first_diff - (first_diff % 3u);
    codon_idx = codon_start / 3u;
    for (j = 0u; j < 3u; j++) {
        ref_codon[j] = (char)s->cds[codon_start + (size_t)j];
        alt_codon[j] = (char)s->expect_cds[codon_start + (size_t)j];
    }
    ref_codon[3] = '\0';
    alt_codon[3] = '\0';
    cr = duckvep_codon_change(ref_codon, alt_codon, DUCKVEP_CODON_TABLE_STANDARD);
    if (cr.change & DUCKVEP_CODON_INVALID) return THEFT_TRIAL_FAIL;
    if (delta.protein_pos != (int32_t)(codon_idx + 1u) ||
        delta.ref_aa != (uint8_t)cr.aa_ref || delta.alt_aa != (uint8_t)cr.aa_alt ||
        delta.cdna_pos != -1 || delta.cds_pos != -1) {
        return THEFT_TRIAL_FAIL;
    }
    if (delta.start_lost || delta.start_retained || delta.frameshift ||
        delta.inframe_deletion ||
        delta.inframe_insertion || delta.protein_altering || delta.coding_unknown) {
        return THEFT_TRIAL_FAIL;
    }
    if (s->expect_region == KPROP_CONTEXT_DELTA_SYNONYMOUS) {
        if (!delta.synonymous || delta.missense || delta.stop_gained ||
            delta.stop_lost || delta.stop_retained) return THEFT_TRIAL_FAIL;
        g_delta_scratch_cov.syn++;
    } else if (s->expect_region == KPROP_CONTEXT_DELTA_MISSENSE) {
        if (!delta.missense || delta.synonymous || delta.stop_gained ||
            delta.stop_lost || delta.stop_retained) return THEFT_TRIAL_FAIL;
        g_delta_scratch_cov.mis++;
    } else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_GAINED) {
        if (!delta.stop_gained || delta.synonymous || delta.missense ||
            delta.stop_lost || delta.stop_retained) return THEFT_TRIAL_FAIL;
        g_delta_scratch_cov.stop_gained++;
    } else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_LOST) {
        if (!delta.stop_lost || delta.synonymous || delta.missense ||
            delta.stop_gained || delta.stop_retained) return THEFT_TRIAL_FAIL;
        g_delta_scratch_cov.stop_lost++;
    } else if (s->expect_region == KPROP_CONTEXT_DELTA_STOP_RETAINED) {
        if (!delta.stop_retained || delta.synonymous || delta.missense ||
            delta.stop_gained || delta.stop_lost) return THEFT_TRIAL_FAIL;
        g_delta_scratch_cov.stop_retained++;
    } else return THEFT_TRIAL_FAIL;
    if (s->strand > 0) g_delta_scratch_cov.fwd++;
    else g_delta_scratch_cov.rev++;
    return THEFT_TRIAL_PASS;
}

TEST sequence_delta_with_scratch_mnv_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "sequence delta scratch MNV == single-codon oracle";
    cfg.prop1 = prop_delta_scratch_mnv_matches_codon_oracle;
    cfg.type_info[0] = &kprop_context_delta_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_delta_scratch_cov, 0, sizeof g_delta_scratch_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_delta_scratch_cov.syn > 0u);
    ASSERT(g_delta_scratch_cov.mis > 0u);
    ASSERT(g_delta_scratch_cov.stop_gained > 0u);
    ASSERT(g_delta_scratch_cov.stop_lost > 0u);
    ASSERT(g_delta_scratch_cov.stop_retained > 0u);
    ASSERT(g_delta_scratch_cov.fwd > 0u);
    ASSERT(g_delta_scratch_cov.rev > 0u);
    ASSERT(g_delta_scratch_cov.capfail > 0u);
    fprintf(stderr,
            "[delta-scratch coverage] syn=%u mis=%u stop_gained=%u stop_lost=%u stop_retained=%u fwd=%u rev=%u capfail=%u\n",
            g_delta_scratch_cov.syn, g_delta_scratch_cov.mis,
            g_delta_scratch_cov.stop_gained, g_delta_scratch_cov.stop_lost,
            g_delta_scratch_cov.stop_retained, g_delta_scratch_cov.fwd,
            g_delta_scratch_cov.rev, g_delta_scratch_cov.capfail);
    PASS();
}

static struct { uint32_t missense; uint32_t synonymous; uint32_t stop_gained; uint32_t fwd; uint32_t rev; uint32_t len2; uint32_t len3; } g_delta_cross_scratch_cov;

static enum theft_trial_res prop_delta_scratch_cross_codon_mnv_matches_oracle(
    struct theft *t, void *arg1) {

    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[64];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    duckvep_delta_scratch_t scratch;
    duckvep_sequence_delta_t direct;
    duckvep_sequence_delta_t shape;
    duckvep_sequence_delta_t routed;
    duckvep_sequence_delta_route_t route;
    struct kprop_cc_facts f;
    (void)t;

    scratch.edits = edits; scratch.edits_cap = 4u;
    scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
    scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
    scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;

    duckvep_sequence_delta_fill_with_scratch((duckvep_variant_kind_t)s->vkind,
                                             &s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->vpos, s->strand,
                                             &scratch, &direct);
    duckvep_sequence_delta_fill((duckvep_variant_kind_t)s->vkind,
                                &s->tx, &s->ex, &s->seq, &s->v,
                                0u, 0u, s->vpos, s->strand, &shape);
    duckvep_sequence_delta_fill_for_annotation_trace((duckvep_variant_kind_t)s->vkind,
                                                     &s->tx, &s->ex, &s->seq, &s->v,
                                                     0u, 0u, s->vpos, s->strand,
                                                     &scratch, NULL, &route, &routed);

    /* The generator only builds two adjacent body codons, so the window oracle always
     * resolves; the authoritative interpreter emits it directly and the router echoes it. */
    f = kprop_cross_codon_mnv_oracle_facts(s);
    if (!f.valid) return THEFT_TRIAL_FAIL;
    if (route != DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT) return THEFT_TRIAL_FAIL;
    if (!direct.valid || !routed.valid) return THEFT_TRIAL_FAIL;
    if (!kprop_sequence_delta_equal(&direct, &routed)) return THEFT_TRIAL_FAIL;
    /* Coarse multi-codon window: protein_pos -1, no AA pair. */
    if (direct.protein_pos != -1 || direct.cdna_pos != -1 || direct.cds_pos != -1 ||
        direct.ref_aa != (uint8_t)0u || direct.alt_aa != (uint8_t)0u) {
        return THEFT_TRIAL_FAIL;
    }
    /* SO facts equal the independent window oracle. */
    if (direct.synonymous != (uint8_t)f.synonymous ||
        direct.missense != (uint8_t)f.missense ||
        direct.stop_gained != (uint8_t)f.stop_gained ||
        direct.stop_lost != (uint8_t)f.stop_lost ||
        direct.stop_retained != (uint8_t)f.stop_retained ||
        direct.start_lost || direct.start_retained || direct.frameshift ||
        direct.inframe_deletion ||
        direct.inframe_insertion || direct.protein_altering || direct.coding_unknown) {
        return THEFT_TRIAL_FAIL;
    }
    /* The direct shape-specific filler only resolved the missense cross-codon window; the
     * synonymous and stop-gained windows are exactly what the authoritative interpreter adds. */
    if (f.missense) {
        if (!shape.valid || !kprop_delta_is_coarse_cross_codon_missense(&shape)) {
            return THEFT_TRIAL_FAIL;
        }
        g_delta_cross_scratch_cov.missense++;
    } else if (f.synonymous) {
        if (shape.valid) return THEFT_TRIAL_FAIL;
        g_delta_cross_scratch_cov.synonymous++;
    } else if (f.stop_gained) {
        if (shape.valid) return THEFT_TRIAL_FAIL;
        g_delta_cross_scratch_cov.stop_gained++;
    } else {
        return THEFT_TRIAL_FAIL;
    }
    if (s->strand > 0) g_delta_cross_scratch_cov.fwd++; else g_delta_cross_scratch_cov.rev++;
    if (s->rlen == 2u) g_delta_cross_scratch_cov.len2++;
    else if (s->rlen == 3u) g_delta_cross_scratch_cov.len3++;
    else return THEFT_TRIAL_FAIL;
    return THEFT_TRIAL_PASS;
}

TEST sequence_delta_with_scratch_cross_codon_mnv_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "sequence delta scratch two-codon MNV window == codon-window oracle";
    cfg.prop1 = prop_delta_scratch_cross_codon_mnv_matches_oracle;
    cfg.type_info[0] = &kprop_cross_codon_mnv_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_delta_cross_scratch_cov, 0, sizeof g_delta_cross_scratch_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_delta_cross_scratch_cov.missense > 0u);
    ASSERT(g_delta_cross_scratch_cov.synonymous > 0u);
    ASSERT(g_delta_cross_scratch_cov.stop_gained > 0u);
    ASSERT(g_delta_cross_scratch_cov.fwd > 0u);
    ASSERT(g_delta_cross_scratch_cov.rev > 0u);
    ASSERT(g_delta_cross_scratch_cov.len2 > 0u);
    ASSERT(g_delta_cross_scratch_cov.len3 > 0u);
    fprintf(stderr,
            "[delta-cross-scratch coverage] missense=%u synonymous=%u stop_gained=%u fwd=%u rev=%u len2=%u len3=%u\n",
            g_delta_cross_scratch_cov.missense, g_delta_cross_scratch_cov.synonymous,
            g_delta_cross_scratch_cov.stop_gained, g_delta_cross_scratch_cov.fwd,
            g_delta_cross_scratch_cov.rev,
            g_delta_cross_scratch_cov.len2, g_delta_cross_scratch_cov.len3);
    PASS();
}

static struct {
    uint32_t ins;
    uint32_t del;
    uint32_t delins;
    uint32_t delins_plus1;
    uint32_t delins_plus2;
    uint32_t delins_minus1;
    uint32_t delins_minus2;
    uint32_t rev;
    uint32_t stop_gained;
    uint32_t terminal_stop;
    uint32_t terminal_reverse;
    uint32_t terminal_missing_tail;
} g_frameshift_cov;

static enum theft_trial_res prop_annotate_frameshift_indel_matches_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    duckvep_haplotype_edit_t edit_scratch[8];
    duckvep_coding_context_t fsctx;
    uint8_t alt_cds[80];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    uint32_t protein_cds;
    uint32_t j;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;
    memset(&err, 0, sizeof err);

    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }

    duckvep_result_builder_init(&rb, rows, 2u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_result_builder_count(&rb) != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }

    if (s->vkind == (uint8_t)DUCKVEP_KIND_INS) {
        protein_cds = kprop_cds_pos_for_genomic(s, s->vpos);
        g_frameshift_cov.ins++;
    } else if (s->vkind == (uint8_t)DUCKVEP_KIND_DEL) {
        protein_cds = UINT32_MAX;
        for (j = 1u; j < (uint32_t)s->rlen; j++) {
            uint32_t cds_pos = kprop_cds_pos_for_genomic(s, s->vpos + j);
            if (cds_pos < protein_cds) protein_cds = cds_pos;
        }
        if (protein_cds == UINT32_MAX) { tr = THEFT_TRIAL_FAIL; goto done; }
        g_frameshift_cov.del++;
    } else if (s->vkind == (uint8_t)DUCKVEP_KIND_INDEL) {
        protein_cds = UINT32_MAX;
        for (j = 0u; j < (uint32_t)s->rlen; j++) {
            uint32_t cds_pos = kprop_cds_pos_for_genomic(s, s->vpos + j);
            if (cds_pos < protein_cds) protein_cds = cds_pos;
        }
        if (protein_cds == UINT32_MAX) { tr = THEFT_TRIAL_FAIL; goto done; }
        g_frameshift_cov.delins++;
        if (s->alen == s->rlen + 1u) g_frameshift_cov.delins_plus1++;
        else if (s->alen == s->rlen + 2u) g_frameshift_cov.delins_plus2++;
        else if (s->rlen == s->alen + 1u) g_frameshift_cov.delins_minus1++;
        else if (s->rlen == s->alen + 2u) g_frameshift_cov.delins_minus2++;
    } else { tr = THEFT_TRIAL_FAIL; goto done; }
    if (s->strand < 0) g_frameshift_cov.rev++;

    /* Frameshift fact plus VEP's local-window stop composite. Rebuild the same coding
     * context the engine uses, then independently re-derive stop_gained from the CDS bytes
     * (kprop_frameshift_local_stop_oracle) for the standard table. For the mitochondrial
     * table the standard-stop enumeration does not apply, so require the frameshift bit and
     * allow at most the stop composite. */
    {
        int have_ctx = duckvep_variant_coding_context_build(
                           &s->tx, &s->ex, &s->seq, &s->v, 0u, 0u, s->strand, edit_scratch, 8u,
                           alt_cds, sizeof alt_cds, ref_pep, sizeof ref_pep,
                           alt_pep, sizeof alt_pep, &fsctx) ==
                       DUCKVEP_VARIANT_CODING_CONTEXT_OK;
        int overlaps_terminal_stop = 0;
        uint32_t terminal_start = 0u;
        if (have_ctx && fsctx.ref_cds_len >= 3u) {
            terminal_start = (uint32_t)fsctx.ref_cds_len - 2u;
        }
        if (have_ctx && terminal_start > 0u &&
            fsctx.ref_peptide_len > 0u &&
            fsctx.ref_peptide[fsctx.ref_peptide_len - 1u] == (uint8_t)'*') {
            if (fsctx.single_edit_ref_len == 0u) {
                overlaps_terminal_stop =
                    fsctx.single_edit_cds_start > terminal_start &&
                    fsctx.single_edit_cds_start <= fsctx.ref_cds_len;
            } else {
                uint64_t edit_end = (uint64_t)fsctx.single_edit_cds_start +
                    (uint64_t)fsctx.single_edit_ref_len - 1u;
                overlaps_terminal_stop =
                    fsctx.single_edit_cds_start <= fsctx.ref_cds_len &&
                    edit_end >= (uint64_t)terminal_start;
            }
        }
        if (overlaps_terminal_stop) {
            if (fsctx.alt_cds_len >= fsctx.ref_cds_len) {
                char codon[3];
                char alt_aa;
                uint64_t want;
                memcpy(codon, fsctx.alt_cds + fsctx.ref_cds_len - 3u,
                       sizeof codon);
                alt_aa = duckvep_translate_codon(
                    codon, DUCKVEP_CODON_TABLE_STANDARD);
                if (alt_aa == '*') {
                    want = DUCKVEP_SO(DUCKVEP_SO_STOP_RETAINED);
                } else if (fsctx.single_edit_cds_start >= terminal_start) {
                    want = DUCKVEP_SO(DUCKVEP_SO_STOP_LOST);
                } else {
                    want = DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT) |
                           DUCKVEP_SO(DUCKVEP_SO_STOP_LOST);
                }
                if (rows[0].consequence_mask != want ||
                    rows[0].sequence_status !=
                        (uint8_t)DUCKVEP_SEQUENCE_RESOLVED) {
                    tr = THEFT_TRIAL_FAIL;
                } else if ((want & DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT)) == 0u &&
                           (rows[0].protein_pos != (int32_t)fsctx.ref_peptide_len ||
                            rows[0].aa_ref != (uint8_t)'*' ||
                            rows[0].aa_alt != (uint8_t)alt_aa)) {
                    tr = THEFT_TRIAL_FAIL;
                } else {
                    g_frameshift_cov.terminal_stop++;
                    if (s->strand < 0) g_frameshift_cov.terminal_reverse++;
                }
            } else {
                if (rows[0].consequence_mask !=
                        DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) ||
                    rows[0].sequence_status !=
                        (uint8_t)DUCKVEP_SEQUENCE_MISSING_TRANSCRIPT_TAIL ||
                    rows[0].protein_pos != -1) {
                    tr = THEFT_TRIAL_FAIL;
                } else {
                    g_frameshift_cov.terminal_missing_tail++;
                }
            }
            goto done;
        }
        if (have_ctx && (duckvep_codon_table_t)s->ctab == DUCKVEP_CODON_TABLE_STANDARD) {
            int expected_stop = kprop_frameshift_local_stop_oracle(
                fsctx.ref_cds, fsctx.ref_cds_len, fsctx.alt_cds, fsctx.alt_cds_len,
                fsctx.single_edit_cds_start, fsctx.single_edit_ref_len, fsctx.length_diff);
            uint64_t want = DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT) |
                            (expected_stop ? DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED) : 0u);
            if (rows[0].consequence_mask != want) { tr = THEFT_TRIAL_FAIL; goto done; }
            if (expected_stop) g_frameshift_cov.stop_gained++;
        } else {
            uint64_t allowed = DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT) |
                               DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED);
            if ((rows[0].consequence_mask & DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT)) == 0u ||
                (rows[0].consequence_mask & ~allowed) != 0u) { tr = THEFT_TRIAL_FAIL; goto done; }
        }
        if (rows[0].cdna_pos != -1 || rows[0].cds_pos != -1) { tr = THEFT_TRIAL_FAIL; goto done; }
        /* Frameshift protein_position = the first affected codon. The authoritative path
         * anchors it to the edit's CDS start (single_edit_cds_start) — the first shifted base,
         * which is VEP's frameshift protein position and, for an insertion at a codon boundary,
         * differs from the anchor base's codon. The edit's cds_start is itself validated
         * against an independent oracle by the cds-edit-builder property tests; when the
         * standalone context build is unavailable we fall back to the genomic-derived codon. */
        {
            int32_t want_pp = have_ctx
                ? (int32_t)(((fsctx.single_edit_cds_start - 1u) / 3u) + 1u)
                : (int32_t)(((protein_cds - 1u) / 3u) + 1u);
            if (rows[0].protein_pos != want_pp) { tr = THEFT_TRIAL_FAIL; goto done; }
        }
    }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return tr;
}

TEST annotate_frameshift_indel_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile simple frameshift indel == CDS-position oracle";
    cfg.prop1 = prop_annotate_frameshift_indel_matches_oracle;
    cfg.type_info[0] = &kprop_frameshift_indel_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_frameshift_cov, 0, sizeof g_frameshift_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_frameshift_cov.ins > 0u);
    ASSERT(g_frameshift_cov.del > 0u);
    ASSERT(g_frameshift_cov.delins > 0u);
    ASSERT(g_frameshift_cov.delins_plus1 > 0u);
    ASSERT(g_frameshift_cov.delins_plus2 > 0u);
    ASSERT(g_frameshift_cov.delins_minus1 > 0u);
    ASSERT(g_frameshift_cov.delins_minus2 > 0u);
    ASSERT(g_frameshift_cov.rev > 0u);
    ASSERT(g_frameshift_cov.terminal_stop > 0u);
    ASSERT(g_frameshift_cov.terminal_reverse > 0u);
    fprintf(stderr,
            "[frameshift coverage] ins=%u del=%u delins=%u(+1=%u +2=%u -1=%u -2=%u) reverse=%u stop_gained=%u terminal_stop=%u terminal_reverse=%u terminal_missing_tail=%u\n",
            g_frameshift_cov.ins, g_frameshift_cov.del, g_frameshift_cov.delins,
            g_frameshift_cov.delins_plus1, g_frameshift_cov.delins_plus2,
            g_frameshift_cov.delins_minus1, g_frameshift_cov.delins_minus2,
            g_frameshift_cov.rev, g_frameshift_cov.stop_gained,
            g_frameshift_cov.terminal_stop,
            g_frameshift_cov.terminal_reverse,
            g_frameshift_cov.terminal_missing_tail);
    PASS();
}

static struct { uint32_t fwd; uint32_t rev; } g_inframe_deletion_cov;

static enum theft_trial_res prop_annotate_inframe_deletion_matches_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    uint32_t min_cds = UINT32_MAX;
    uint32_t j;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;
    memset(&err, 0, sizeof err);

    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    duckvep_workspace_delta_route_stats_reset(ws);

    duckvep_result_builder_init(&rb, rows, 2u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_result_builder_count(&rb) != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }

    for (j = 1u; j < (uint32_t)s->rlen; j++) {
        uint32_t cds_pos = kprop_cds_pos_for_genomic(s, s->vpos + j);
        if (cds_pos < min_cds) min_cds = cds_pos;
    }
    if (min_cds == UINT32_MAX || ((min_cds - 1u) % 3u) != 0u) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (s->strand > 0) g_inframe_deletion_cov.fwd++; else g_inframe_deletion_cov.rev++;

    if (rows[0].consequence_mask != DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION)) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].cdna_pos != -1 || rows[0].cds_pos != -1) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].protein_pos != (int32_t)(((min_cds - 1u) / 3u) + 1u)) { tr = THEFT_TRIAL_FAIL; goto done; }
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL || stats->del_context != 1u ||
        stats->del_direct_fallback != 0u ||
        stats->substitution_context != 0u ||
        stats->mnv_direct_fallback != 0u) {
        tr = THEFT_TRIAL_FAIL; goto done;
    }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return tr;
}

TEST annotate_inframe_deletion_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile codon-aligned in-frame deletion == CDS-position oracle";
    cfg.prop1 = prop_annotate_inframe_deletion_matches_oracle;
    cfg.type_info[0] = &kprop_inframe_deletion_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_inframe_deletion_cov, 0, sizeof g_inframe_deletion_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_inframe_deletion_cov.fwd > 0u);
    ASSERT(g_inframe_deletion_cov.rev > 0u);
    fprintf(stderr, "[inframe_deletion coverage] forward=%u reverse=%u\n",
            g_inframe_deletion_cov.fwd, g_inframe_deletion_cov.rev);
    PASS();
}

static struct { uint32_t fwd; uint32_t rev; uint32_t full; } g_cursor_del_route_cov;

static enum theft_trial_res prop_cursor_del_route_matches_tile(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_annotate_cursor_t *cur = NULL;
    duckvep_error_t err;
    duckvep_consequence_t tile_rows[4];
    duckvep_consequence_t cursor_rows[4];
    duckvep_consequence_t chunk[1];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    duckvep_workspace_delta_route_stats_t tile_stats;
    duckvep_workspace_delta_route_stats_t cursor_stats;
    size_t tile_n;
    size_t cursor_n = 0u;
    uint32_t min_cds = UINT32_MAX;
    uint32_t j;
    int saw_full = 0;
    enum theft_trial_res res = THEFT_TRIAL_PASS;
    (void)t;

    memset(&err, 0, sizeof err);
    memset(&tile_stats, 0, sizeof tile_stats);
    memset(&cursor_stats, 0, sizeof cursor_stats);
    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, &model, &err) != DUCKVEP_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }

    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, tile_rows, 4u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    tile_n = duckvep_result_builder_count(&rb);
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { res = THEFT_TRIAL_FAIL; goto done; }
    tile_stats = *stats;

    duckvep_workspace_delta_route_stats_reset(ws);
    if (duckvep_annotate_cursor_open(model, &s->v, opts, ws, &cur, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    while (!duckvep_annotate_cursor_done(cur)) {
        duckvep_status_t st;
        size_t i;
        duckvep_result_builder_init(&rb, chunk, 1u);
        st = duckvep_annotate_cursor_fill(cur, &rb, &err);
        if (st != DUCKVEP_OK && st != DUCKVEP_ERR_RESULT_FULL) { res = THEFT_TRIAL_FAIL; goto done; }
        if (st == DUCKVEP_ERR_RESULT_FULL) saw_full = 1;
        for (i = 0u; i < duckvep_result_builder_count(&rb); i++) {
            if (cursor_n >= 4u) { res = THEFT_TRIAL_FAIL; goto done; }
            cursor_rows[cursor_n++] = chunk[i];
        }
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { res = THEFT_TRIAL_FAIL; goto done; }
    cursor_stats = *stats;

    if (!saw_full || tile_n != cursor_n || tile_n != 1u) { res = THEFT_TRIAL_FAIL; goto done; }
    if (tile_stats.del_context != 1u ||
        tile_stats.del_direct_fallback != 0u ||
        tile_stats.substitution_context != 0u ||
        tile_stats.mnv_direct_fallback != 0u ||
        tile_stats.del_context != cursor_stats.del_context ||
        tile_stats.del_direct_fallback != cursor_stats.del_direct_fallback) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    if (!consequence_rows_equal(&tile_rows[0], &cursor_rows[0])) { res = THEFT_TRIAL_FAIL; goto done; }
    for (j = 1u; j < (uint32_t)s->rlen; j++) {
        uint32_t cds_pos = kprop_cds_pos_for_genomic(s, s->vpos + j);
        if (cds_pos < min_cds) min_cds = cds_pos;
    }
    if (min_cds == UINT32_MAX || tile_rows[0].consequence_mask != DUCKVEP_SO(DUCKVEP_SO_INFRAME_DELETION) ||
        tile_rows[0].protein_pos != (int32_t)(((min_cds - 1u) / 3u) + 1u)) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    if (s->strand > 0) g_cursor_del_route_cov.fwd++;
    else g_cursor_del_route_cov.rev++;
    if (saw_full) g_cursor_del_route_cov.full++;

done:
    duckvep_annotate_cursor_close(cur);
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return res;
}

TEST annotate_cursor_del_route_matches_tile_for_any_output_split(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate cursor DEL route == tile under output splits";
    cfg.prop1 = prop_cursor_del_route_matches_tile;
    cfg.type_info[0] = &kprop_inframe_deletion_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_cursor_del_route_cov, 0, sizeof g_cursor_del_route_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_cursor_del_route_cov.fwd > 0u);
    ASSERT(g_cursor_del_route_cov.rev > 0u);
    ASSERT(g_cursor_del_route_cov.full > 0u);
    fprintf(stderr, "[cursor-del-route coverage] forward=%u reverse=%u full=%u\n",
            g_cursor_del_route_cov.fwd, g_cursor_del_route_cov.rev,
            g_cursor_del_route_cov.full);
    PASS();
}

static struct { uint32_t fwd; uint32_t rev; } g_context_inframe_deletion_cov;

static enum theft_trial_res prop_context_inframe_deletion_matches_oracle(struct theft *t,
                                                                        void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[64];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;
    uint32_t min_cds = UINT32_MAX;
    uint32_t j;
    (void)t;

    if (duckvep_variant_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->strand,
                                             edits, 4u, alt_cds, sizeof alt_cds,
                                             ref_pep, sizeof ref_pep,
                                             alt_pep, sizeof alt_pep,
                                             &ctx) != DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (!ctx.has_single_edit || ctx.single_edit_alt_len != 0u ||
        ctx.single_edit_ref_len != 3u || ctx.applied_edits != 1u) {
        return THEFT_TRIAL_FAIL;
    }
    for (j = 1u; j < (uint32_t)s->rlen; j++) {
        uint32_t cds_pos = kprop_cds_pos_for_genomic(s, s->vpos + j);
        if (cds_pos < min_cds) min_cds = cds_pos;
    }
    if (min_cds == UINT32_MAX) return THEFT_TRIAL_FAIL;
    if (ctx.single_edit_cds_start != min_cds) return THEFT_TRIAL_FAIL;
    if (duckvep_coding_context_delta_fill(&ctx, 0u, &delta) != DUCKVEP_CONTEXT_DELTA_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (!kprop_delta_is_inframe_deletion_at(
            &delta, (int32_t)(((min_cds - 1u) / 3u) + 1u))) {
        return THEFT_TRIAL_FAIL;
    }
    if (s->strand > 0) g_context_inframe_deletion_cov.fwd++;
    else g_context_inframe_deletion_cov.rev++;
    return THEFT_TRIAL_PASS;
}

TEST coding_context_delta_inframe_deletion_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "coding context delta in-frame deletion == edit-origin oracle";
    cfg.prop1 = prop_context_inframe_deletion_matches_oracle;
    cfg.type_info[0] = &kprop_inframe_deletion_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_context_inframe_deletion_cov, 0, sizeof g_context_inframe_deletion_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_context_inframe_deletion_cov.fwd > 0u);
    ASSERT(g_context_inframe_deletion_cov.rev > 0u);
    fprintf(stderr, "[context-inframe-deletion coverage] forward=%u reverse=%u\n",
            g_context_inframe_deletion_cov.fwd, g_context_inframe_deletion_cov.rev);
    PASS();
}

static struct { uint32_t fwd; uint32_t rev; } g_context_inframe_insertion_cov;

static enum theft_trial_res prop_context_inframe_insertion_matches_oracle(struct theft *t,
                                                                         void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[64];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;
    uint32_t anchor_cds;
    uint32_t before_cds;
    (void)t;

    if (duckvep_variant_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->strand,
                                             edits, 4u, alt_cds, sizeof alt_cds,
                                             ref_pep, sizeof ref_pep,
                                             alt_pep, sizeof alt_pep,
                                             &ctx) != DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    anchor_cds = kprop_cds_pos_for_genomic(s, s->vpos);
    before_cds = s->strand > 0 ? anchor_cds : anchor_cds - 1u;
    if ((before_cds % 3u) != 0u || before_cds <= 3u || before_cds >= s->cds_lenv - 3u) {
        return THEFT_TRIAL_FAIL;
    }
    if (!ctx.has_single_edit || ctx.single_edit_ref_len != 0u ||
        ctx.single_edit_alt_len != 3u || ctx.applied_edits != 1u ||
        ctx.single_edit_cds_start != before_cds + 1u) {
        return THEFT_TRIAL_FAIL;
    }
    if (duckvep_coding_context_delta_fill(&ctx, 0u, &delta) != DUCKVEP_CONTEXT_DELTA_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (!kprop_delta_is_inframe_insertion_at(
            &delta, (int32_t)(before_cds / 3u + 1u))) {
        return THEFT_TRIAL_FAIL;
    }
    if (s->strand > 0) g_context_inframe_insertion_cov.fwd++;
    else g_context_inframe_insertion_cov.rev++;
    return THEFT_TRIAL_PASS;
}

TEST coding_context_delta_inframe_insertion_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "coding context delta in-frame insertion == edit-origin oracle";
    cfg.prop1 = prop_context_inframe_insertion_matches_oracle;
    cfg.type_info[0] = &kprop_inframe_insertion_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_context_inframe_insertion_cov, 0, sizeof g_context_inframe_insertion_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_context_inframe_insertion_cov.fwd > 0u);
    ASSERT(g_context_inframe_insertion_cov.rev > 0u);
    fprintf(stderr, "[context-inframe-insertion coverage] forward=%u reverse=%u\n",
            g_context_inframe_insertion_cov.fwd, g_context_inframe_insertion_cov.rev);
    PASS();
}

static struct {
    uint32_t fwd;
    uint32_t rev;
    uint32_t lengthen;
    uint32_t shorten;
    uint32_t inframe;
    uint32_t protein_altering;
} g_context_delins_shape_cov;

static enum theft_trial_res prop_context_delins_shape_matches_oracle(struct theft *t,
                                                                    void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[96];
    uint8_t ref_pep[40];
    uint8_t alt_pep[40];
    duckvep_coding_context_t ctx;
    duckvep_sequence_delta_t delta;
    int32_t protein_pos;
    (void)t;

    if (duckvep_variant_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->strand,
                                             edits, 4u, alt_cds, sizeof alt_cds,
                                             ref_pep, sizeof ref_pep,
                                             alt_pep, sizeof alt_pep,
                                             &ctx) != DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (!ctx.has_single_edit || ctx.applied_edits != 1u ||
        (ctx.single_edit_ref_len == 0u && ctx.single_edit_alt_len == 0u) ||
        ctx.length_diff == 0 || (ctx.length_diff % 3) != 0 ||
        ctx.single_edit_cds_start <= 3u) {
        return THEFT_TRIAL_FAIL;
    }
    if ((uint64_t)ctx.single_edit_cds_start +
            (uint64_t)ctx.single_edit_ref_len - 1u >
        (uint64_t)ctx.ref_cds_len - 3u) {
        return THEFT_TRIAL_FAIL;
    }
    protein_pos = (int32_t)(((ctx.single_edit_cds_start - 1u) / 3u) + 1u);
    if (duckvep_coding_context_delta_fill(&ctx, 0u, &delta) != DUCKVEP_CONTEXT_DELTA_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (s->expect_protein_altering) {
        if (!kprop_delta_is_protein_altering_at(&delta, protein_pos)) {
            return THEFT_TRIAL_FAIL;
        }
        g_context_delins_shape_cov.protein_altering++;
    } else if (ctx.length_diff > 0) {
        if (!kprop_delta_is_inframe_insertion_at(&delta, protein_pos)) {
            return THEFT_TRIAL_FAIL;
        }
        g_context_delins_shape_cov.inframe++;
    } else {
        if (!kprop_delta_is_inframe_deletion_at(&delta, protein_pos)) {
            return THEFT_TRIAL_FAIL;
        }
        g_context_delins_shape_cov.inframe++;
    }
    if (ctx.length_diff > 0) g_context_delins_shape_cov.lengthen++;
    else g_context_delins_shape_cov.shorten++;
    if (s->strand > 0) g_context_delins_shape_cov.fwd++;
    else g_context_delins_shape_cov.rev++;
    return THEFT_TRIAL_PASS;
}

TEST coding_context_delta_delins_shape_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "coding context delins shape == local-edge oracle";
    cfg.prop1 = prop_context_delins_shape_matches_oracle;
    cfg.type_info[0] = &kprop_delins_shape_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_context_delins_shape_cov, 0, sizeof g_context_delins_shape_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_context_delins_shape_cov.fwd > 0u);
    ASSERT(g_context_delins_shape_cov.rev > 0u);
    ASSERT(g_context_delins_shape_cov.lengthen > 0u);
    ASSERT(g_context_delins_shape_cov.shorten > 0u);
    ASSERT(g_context_delins_shape_cov.inframe > 0u);
    ASSERT(g_context_delins_shape_cov.protein_altering > 0u);
    fprintf(stderr,
            "[context-delins-shape coverage] forward=%u reverse=%u lengthen=%u "
            "shorten=%u inframe=%u protein_altering=%u\n",
            g_context_delins_shape_cov.fwd, g_context_delins_shape_cov.rev,
            g_context_delins_shape_cov.lengthen, g_context_delins_shape_cov.shorten,
            g_context_delins_shape_cov.inframe,
            g_context_delins_shape_cov.protein_altering);
    PASS();
}

static struct { uint32_t fwd; uint32_t rev; uint32_t lengthen; uint32_t shorten; }
    g_delta_scratch_indel_cov;

static enum theft_trial_res prop_delta_scratch_indel_matches_oracle(struct theft *t,
                                                                   void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_haplotype_edit_t edits[4];
    uint8_t alt_cds[96];
    uint8_t ref_pep[40];
    uint8_t alt_pep[40];
    duckvep_delta_scratch_t scratch;
    duckvep_sequence_delta_t delta;
    duckvep_coding_context_t ctx;
    int32_t protein_pos;
    (void)t;

    memset(&scratch, 0, sizeof scratch);
    scratch.edits = edits; scratch.edits_cap = 4u;
    scratch.alt_cds = alt_cds; scratch.alt_cds_cap = sizeof alt_cds;
    scratch.ref_peptide = ref_pep; scratch.ref_peptide_cap = sizeof ref_pep;
    scratch.alt_peptide = alt_pep; scratch.alt_peptide_cap = sizeof alt_pep;

    if (duckvep_variant_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v,
                                             0u, 0u, s->strand,
                                             edits, 4u, alt_cds, sizeof alt_cds,
                                             ref_pep, sizeof ref_pep,
                                             alt_pep, sizeof alt_pep,
                                             &ctx) != DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        return THEFT_TRIAL_FAIL;
    }
    protein_pos = (int32_t)(((ctx.single_edit_cds_start - 1u) / 3u) + 1u);
    duckvep_sequence_delta_fill_with_scratch(DUCKVEP_KIND_INDEL, &s->tx, &s->ex,
                                             &s->seq, &s->v, 0u, 0u, s->vpos,
                                             s->strand, &scratch, &delta);
    if (s->expect_protein_altering) {
        if (!kprop_delta_is_protein_altering_at(&delta, protein_pos)) {
            return THEFT_TRIAL_FAIL;
        }
    } else if (ctx.length_diff > 0) {
        if (!kprop_delta_is_inframe_insertion_at(&delta, protein_pos)) {
            return THEFT_TRIAL_FAIL;
        }
    } else if (ctx.length_diff < 0) {
        if (!kprop_delta_is_inframe_deletion_at(&delta, protein_pos)) {
            return THEFT_TRIAL_FAIL;
        }
    } else return THEFT_TRIAL_FAIL;
    if (ctx.length_diff > 0) g_delta_scratch_indel_cov.lengthen++;
    else g_delta_scratch_indel_cov.shorten++;
    if (s->strand > 0) g_delta_scratch_indel_cov.fwd++;
    else g_delta_scratch_indel_cov.rev++;
    return THEFT_TRIAL_PASS;
}

TEST sequence_delta_with_scratch_indel_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "sequence delta scratch INDEL == local delins-shape oracle";
    cfg.prop1 = prop_delta_scratch_indel_matches_oracle;
    cfg.type_info[0] = &kprop_delins_shape_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_delta_scratch_indel_cov, 0, sizeof g_delta_scratch_indel_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_delta_scratch_indel_cov.fwd > 0u);
    ASSERT(g_delta_scratch_indel_cov.rev > 0u);
    ASSERT(g_delta_scratch_indel_cov.lengthen > 0u);
    ASSERT(g_delta_scratch_indel_cov.shorten > 0u);
    fprintf(stderr,
            "[delta-scratch-indel coverage] forward=%u reverse=%u lengthen=%u shorten=%u\n",
            g_delta_scratch_indel_cov.fwd, g_delta_scratch_indel_cov.rev,
            g_delta_scratch_indel_cov.lengthen, g_delta_scratch_indel_cov.shorten);
    PASS();
}

static struct { uint32_t fwd; uint32_t rev; } g_inframe_insertion_cov;

static enum theft_trial_res prop_annotate_inframe_insertion_matches_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    uint32_t anchor_cds;
    uint32_t before_cds;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;
    memset(&err, 0, sizeof err);

    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    duckvep_workspace_delta_route_stats_reset(ws);

    duckvep_result_builder_init(&rb, rows, 2u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_result_builder_count(&rb) != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }

    anchor_cds = kprop_cds_pos_for_genomic(s, s->vpos);
    before_cds = s->strand > 0 ? anchor_cds : anchor_cds - 1u;
    if ((before_cds % 3u) != 0u || before_cds <= 3u || before_cds >= s->cds_lenv - 3u) {
        tr = THEFT_TRIAL_FAIL; goto done;
    }
    if (s->strand > 0) g_inframe_insertion_cov.fwd++; else g_inframe_insertion_cov.rev++;

    if (rows[0].consequence_mask != DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION)) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].cdna_pos != -1 || rows[0].cds_pos != -1) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].protein_pos != (int32_t)((before_cds / 3u) + 1u)) { tr = THEFT_TRIAL_FAIL; goto done; }
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL || stats->substitution_context != 0u ||
        stats->mnv_direct_fallback != 0u ||
        stats->del_context != 0u ||
        stats->del_direct_fallback != 0u ||
        stats->ins_context != 1u ||
        stats->ins_direct_fallback != 0u) {
        tr = THEFT_TRIAL_FAIL; goto done;
    }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return tr;
}

TEST annotate_inframe_insertion_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile codon-boundary in-frame insertion == CDS-position oracle";
    cfg.prop1 = prop_annotate_inframe_insertion_matches_oracle;
    cfg.type_info[0] = &kprop_inframe_insertion_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_inframe_insertion_cov, 0, sizeof g_inframe_insertion_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_inframe_insertion_cov.fwd > 0u);
    ASSERT(g_inframe_insertion_cov.rev > 0u);
    fprintf(stderr, "[inframe_insertion coverage] forward=%u reverse=%u\n",
            g_inframe_insertion_cov.fwd, g_inframe_insertion_cov.rev);
    PASS();
}

static struct { uint32_t fwd; uint32_t rev; uint32_t full; } g_cursor_ins_route_cov;

static enum theft_trial_res prop_cursor_ins_route_matches_tile(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_annotate_cursor_t *cur = NULL;
    duckvep_error_t err;
    duckvep_consequence_t tile_rows[4];
    duckvep_consequence_t cursor_rows[4];
    duckvep_consequence_t chunk[1];
    duckvep_result_builder_t rb;
    const duckvep_workspace_delta_route_stats_t *stats;
    duckvep_workspace_delta_route_stats_t tile_stats;
    duckvep_workspace_delta_route_stats_t cursor_stats;
    size_t tile_n;
    size_t cursor_n = 0u;
    uint32_t anchor_cds;
    uint32_t before_cds;
    int saw_full = 0;
    enum theft_trial_res res = THEFT_TRIAL_PASS;
    (void)t;

    memset(&err, 0, sizeof err);
    memset(&tile_stats, 0, sizeof tile_stats);
    memset(&cursor_stats, 0, sizeof cursor_stats);
    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, &model, &err) != DUCKVEP_OK) {
        return THEFT_TRIAL_FAIL;
    }
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { res = THEFT_TRIAL_FAIL; goto done; }

    duckvep_workspace_delta_route_stats_reset(ws);
    duckvep_result_builder_init(&rb, tile_rows, 4u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    tile_n = duckvep_result_builder_count(&rb);
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { res = THEFT_TRIAL_FAIL; goto done; }
    tile_stats = *stats;

    duckvep_workspace_delta_route_stats_reset(ws);
    if (duckvep_annotate_cursor_open(model, &s->v, opts, ws, &cur, &err) != DUCKVEP_OK) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    while (!duckvep_annotate_cursor_done(cur)) {
        duckvep_status_t st;
        size_t i;
        duckvep_result_builder_init(&rb, chunk, 1u);
        st = duckvep_annotate_cursor_fill(cur, &rb, &err);
        if (st != DUCKVEP_OK && st != DUCKVEP_ERR_RESULT_FULL) { res = THEFT_TRIAL_FAIL; goto done; }
        if (st == DUCKVEP_ERR_RESULT_FULL) saw_full = 1;
        for (i = 0u; i < duckvep_result_builder_count(&rb); i++) {
            if (cursor_n >= 4u) { res = THEFT_TRIAL_FAIL; goto done; }
            cursor_rows[cursor_n++] = chunk[i];
        }
    }
    stats = duckvep_workspace_delta_route_stats(ws);
    if (stats == NULL) { res = THEFT_TRIAL_FAIL; goto done; }
    cursor_stats = *stats;

    if (!saw_full || tile_n != cursor_n || tile_n != 1u) { res = THEFT_TRIAL_FAIL; goto done; }
    if (tile_stats.ins_context != 1u ||
        tile_stats.ins_direct_fallback != 0u ||
        tile_stats.substitution_context != 0u ||
        tile_stats.mnv_direct_fallback != 0u ||
        tile_stats.del_context != 0u ||
        tile_stats.del_direct_fallback != 0u ||
        tile_stats.ins_context != cursor_stats.ins_context ||
        tile_stats.ins_direct_fallback != cursor_stats.ins_direct_fallback) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    if (!consequence_rows_equal(&tile_rows[0], &cursor_rows[0])) { res = THEFT_TRIAL_FAIL; goto done; }
    anchor_cds = kprop_cds_pos_for_genomic(s, s->vpos);
    before_cds = s->strand > 0 ? anchor_cds : anchor_cds - 1u;
    if ((before_cds % 3u) != 0u ||
        tile_rows[0].consequence_mask != DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION) ||
        tile_rows[0].protein_pos != (int32_t)((before_cds / 3u) + 1u)) {
        res = THEFT_TRIAL_FAIL; goto done;
    }
    if (s->strand > 0) g_cursor_ins_route_cov.fwd++;
    else g_cursor_ins_route_cov.rev++;
    if (saw_full) g_cursor_ins_route_cov.full++;

done:
    duckvep_annotate_cursor_close(cur);
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return res;
}

TEST annotate_cursor_ins_route_matches_tile_for_any_output_split(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate cursor INS route == tile under output splits";
    cfg.prop1 = prop_cursor_ins_route_matches_tile;
    cfg.type_info[0] = &kprop_inframe_insertion_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_cursor_ins_route_cov, 0, sizeof g_cursor_ins_route_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_cursor_ins_route_cov.fwd > 0u);
    ASSERT(g_cursor_ins_route_cov.rev > 0u);
    ASSERT(g_cursor_ins_route_cov.full > 0u);
    fprintf(stderr, "[cursor-ins-route coverage] forward=%u reverse=%u full=%u\n",
            g_cursor_ins_route_cov.fwd, g_cursor_ins_route_cov.rev,
            g_cursor_ins_route_cov.full);
    PASS();
}

/* Independent VEP inframe_insertion test: is the ref peptide window a prefix OR suffix of the
 * alt window (an empty ref window trivially yes)? Mirrors VariationEffect::inframe_insertion,
 * re-derived here so the test does not restate the kernel helper. */
static int kprop_pep_window_prefix_or_suffix(
    const uint8_t *outer, uint32_t outer_off, uint32_t outer_len,
    const uint8_t *inner, uint32_t inner_off, uint32_t inner_len) {
    uint32_t i;
    int ok;
    if (inner_len == 0u) return 1;
    if (inner_len > outer_len) return 0;
    ok = 1;
    for (i = 0u; i < inner_len; i++) {
        if (outer[outer_off + i] != inner[inner_off + i]) { ok = 0; break; }
    }
    if (ok) return 1;
    ok = 1;
    for (i = 0u; i < inner_len; i++) {
        if (outer[outer_off + outer_len - inner_len + i] != inner[inner_off + i]) { ok = 0; break; }
    }
    return ok;
}

static struct { uint32_t fwd; uint32_t rev; uint32_t inframe; uint32_t altering; } g_protein_altering_cov;

static enum theft_trial_res prop_annotate_protein_altering_insertion_matches_oracle(struct theft *t, void *arg1) {
    const struct kprop_coding *s = (const struct kprop_coding *)arg1;
    duckvep_model_t *model = NULL;
    duckvep_options_t *opts = NULL;
    duckvep_workspace_t *ws = NULL;
    duckvep_error_t err;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t rb;
    duckvep_haplotype_edit_t edit_scratch[4];
    duckvep_coding_context_t ctx;
    uint8_t alt_cds[80];
    uint8_t ref_pep[32];
    uint8_t alt_pep[32];
    uint32_t anchor_cds;
    uint32_t before_cds;
    uint32_t rf, rl, af, al;
    int inframe;
    uint64_t want;
    enum theft_trial_res tr = THEFT_TRIAL_PASS;
    (void)t;
    memset(&err, 0, sizeof err);

    if (duckvep_model_open(&s->tx, &s->ex, &s->seq, &model, &err) != DUCKVEP_OK) return THEFT_TRIAL_FAIL;
    if (duckvep_options_open(NULL, &opts, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_workspace_open(model, &ws, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }

    duckvep_result_builder_init(&rb, rows, 2u);
    if (duckvep_annotate_tile(model, &s->v, opts, ws, &rb, &err) != DUCKVEP_OK) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (duckvep_result_builder_count(&rb) != 1u) { tr = THEFT_TRIAL_FAIL; goto done; }

    anchor_cds = kprop_cds_pos_for_genomic(s, s->vpos);
    before_cds = s->strand > 0 ? anchor_cds : anchor_cds - 1u;
    if ((before_cds % 3u) == 0u || before_cds <= 3u || before_cds >= s->cds_lenv - 3u) {
        tr = THEFT_TRIAL_FAIL; goto done;
    }

    /* Independent expected term: build the coding context (the alt CDS construction is itself
     * validated by the cds-edit-builder oracle), diff the peptides, and apply VEP's prefix/
     * suffix rule. A non-boundary insertion is inframe_insertion when the flanking residues are
     * preserved (empty ref window, or ref window prefix/suffix of alt) and protein_altering when
     * the junction residue also changes. */
    if (duckvep_variant_coding_context_build(&s->tx, &s->ex, &s->seq, &s->v, 0u, 0u, s->strand,
                                             edit_scratch, 4u, alt_cds, sizeof alt_cds,
                                             ref_pep, sizeof ref_pep, alt_pep, sizeof alt_pep,
                                             &ctx) != DUCKVEP_VARIANT_CODING_CONTEXT_OK) {
        tr = THEFT_TRIAL_FAIL; goto done;
    }
    kprop_peptide_window_oracle(ctx.ref_peptide, ctx.ref_peptide_len,
                                ctx.alt_peptide, ctx.alt_peptide_len, &rf, &rl, &af, &al);
    if (af == 0u) { tr = THEFT_TRIAL_FAIL; goto done; }
    inframe = (rf == 0u) ||
              kprop_pep_window_prefix_or_suffix(ctx.alt_peptide, af - 1u, al - af + 1u,
                                                ctx.ref_peptide, rf - 1u, rl - rf + 1u);
    want = inframe ? DUCKVEP_SO(DUCKVEP_SO_INFRAME_INSERTION)
                   : DUCKVEP_SO(DUCKVEP_SO_PROTEIN_ALTERING);
    if (s->strand > 0) g_protein_altering_cov.fwd++; else g_protein_altering_cov.rev++;
    if (inframe) g_protein_altering_cov.inframe++; else g_protein_altering_cov.altering++;

    if (rows[0].consequence_mask != want) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].cdna_pos != -1 || rows[0].cds_pos != -1) { tr = THEFT_TRIAL_FAIL; goto done; }
    if (rows[0].protein_pos != (int32_t)((before_cds / 3u) + 1u)) { tr = THEFT_TRIAL_FAIL; goto done; }

done:
    duckvep_workspace_close(ws);
    duckvep_options_close(opts);
    duckvep_model_close(model);
    return tr;
}

TEST annotate_protein_altering_insertion_matches_oracle(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile non-boundary in-frame insertion == peptide-window oracle";
    cfg.prop1 = prop_annotate_protein_altering_insertion_matches_oracle;
    cfg.type_info[0] = &kprop_protein_altering_insertion_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_protein_altering_cov, 0, sizeof g_protein_altering_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_protein_altering_cov.fwd > 0u);
    ASSERT(g_protein_altering_cov.rev > 0u);
    ASSERT(g_protein_altering_cov.inframe > 0u);
    ASSERT(g_protein_altering_cov.altering > 0u);
    fprintf(stderr, "[non-boundary insertion coverage] forward=%u reverse=%u inframe_insertion=%u protein_altering=%u\n",
            g_protein_altering_cov.fwd, g_protein_altering_cov.rev,
            g_protein_altering_cov.inframe, g_protein_altering_cov.altering);
    PASS();
}

TEST annotate_start_lost_matches_oracle_for_any_start_codon_snv(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "annotate_tile start_lost SNV == start-codon oracle";
    cfg.prop1 = prop_annotate_start_lost_matches_oracle;
    cfg.type_info[0] = &kprop_start_codon_info;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_start_cov, 0, sizeof g_start_cov);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    ASSERT(g_start_cov.sl > 0u);
    ASSERT(g_start_cov.sl_sg > 0u);
    ASSERT(g_start_cov.sl_syn > 0u);
    ASSERT(g_start_cov.alt_m > 0u);
    fprintf(stderr, "[start_lost coverage] start_lost=%u co_stop_gained=%u co_syn=%u alt_m_no_start=%u\n",
            g_start_cov.sl, g_start_cov.sl_sg, g_start_cov.sl_syn, g_start_cov.alt_m);
    PASS();
}

/* ========================================================================
 * Variant-tile storage and sorted-run validation.
 * Pure C, no DuckDB. Small-variant plus single-interval SV/CNV transport. The tile
 * requires contiguous sequence-region runs, nondecreasing positions within a run,
 * and no sequence region reappearing after its run closes. Run state persists across
 * reset_rows.
 * ====================================================================== */

/* Append one SNV with the common-case shape (end1==pos1, ref/alt single base). */
static duckvep_tile_status_t tile_put_snv(duckvep_variant_tile_t *tile, uint16_t c,
                                          uint32_t p, char ref, char alt,
                                          const char *id) {
    return duckvep_variant_tile_append(tile, c, p, p, (uint8_t)DUCKVEP_KIND_SNV,
                                       &ref, 1u, &alt, 1u, id, strlen(id));
}

static duckvep_tile_status_t tile_put_mnv(duckvep_variant_tile_t *tile, uint16_t c,
                                          uint32_t p, const char *ref, const char *alt,
                                          size_t len, const char *id) {
    return duckvep_variant_tile_append(tile, c, p, p + (uint32_t)len - 1u,
                                       (uint8_t)DUCKVEP_KIND_MNV, ref, len, alt, len,
                                       id, strlen(id));
}

static duckvep_tile_status_t tile_put_ins(duckvep_variant_tile_t *tile, uint16_t c,
                                          uint32_t p, const char *ref, const char *alt,
                                          size_t alt_len, const char *id) {
    return duckvep_variant_tile_append(tile, c, p, p, (uint8_t)DUCKVEP_KIND_INS,
                                       ref, 1u, alt, alt_len, id, strlen(id));
}

static duckvep_tile_status_t tile_put_del(duckvep_variant_tile_t *tile, uint16_t c,
                                          uint32_t p, const char *ref, size_t ref_len,
                                          const char *alt, const char *id) {
    return duckvep_variant_tile_append(tile, c, p, p + (uint32_t)ref_len - 1u,
                                       (uint8_t)DUCKVEP_KIND_DEL, ref, ref_len, alt, 1u,
                                       id, strlen(id));
}

static duckvep_tile_status_t tile_put_indel(duckvep_variant_tile_t *tile, uint16_t c,
                                            uint32_t p, const char *ref, size_t ref_len,
                                            const char *alt, size_t alt_len,
                                            const char *id) {
    return duckvep_variant_tile_append(tile, c, p, p + (uint32_t)ref_len - 1u,
                                       (uint8_t)DUCKVEP_KIND_INDEL, ref, ref_len,
                                       alt, alt_len, id, strlen(id));
}

TEST tile_rejects_breakends_and_malformed_small_variants(void) {
    duckvep_variant_tile_t *tile = duckvep_variant_tile_create(8u, 4u);
    char r = 'A', a = 'C';
    ASSERT(tile != NULL);

    /* breakends require a two-locus event representation. */
    ASSERT_EQ(DUCKVEP_TILE_INVALID,
              duckvep_variant_tile_append_event(
                  tile, 0u, 10u, 10u, (uint8_t)DUCKVEP_KIND_SV,
                  (uint8_t)DUCKVEP_SV_BREAKEND,
                  (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN,
                  &r, 1u, &a, 1u, "v", 1u));
    /* NONE and contradictory direction/subtype pairs are invalid. */
    ASSERT_EQ(DUCKVEP_TILE_INVALID,
              duckvep_variant_tile_append_event(
                  tile, 0u, 10u, 20u, (uint8_t)DUCKVEP_KIND_SV,
                  (uint8_t)DUCKVEP_SV_NONE,
                  (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN,
                  &r, 1u, &a, 1u, "v", 1u));
    ASSERT_EQ(DUCKVEP_TILE_INVALID,
              duckvep_variant_tile_append_event(
                  tile, 0u, 10u, 20u, (uint8_t)DUCKVEP_KIND_SV,
                  (uint8_t)DUCKVEP_SV_DELETION,
                  (uint8_t)DUCKVEP_COPY_CHANGE_GAIN,
                  &r, 1u, &a, 1u, "v", 1u));
    ASSERT_EQ(DUCKVEP_TILE_INVALID,
              duckvep_variant_tile_append_event(
                  tile, 0u, 10u, 20u, (uint8_t)DUCKVEP_KIND_SV,
                  (uint8_t)DUCKVEP_SV_DUPLICATION,
                  (uint8_t)DUCKVEP_COPY_CHANGE_LOSS,
                  &r, 1u, &a, 1u, "v", 1u));
    ASSERT_EQ(DUCKVEP_TILE_INVALID,
              duckvep_variant_tile_append_event(
                  tile, 0u, 10u, 20u, (uint8_t)DUCKVEP_KIND_SV,
                  (uint8_t)DUCKVEP_SV_INSERTION,
                  (uint8_t)DUCKVEP_COPY_CHANGE_LOSS,
                  &r, 1u, &a, 1u, "v", 1u));
    ASSERT_EQ(DUCKVEP_TILE_INVALID,
              duckvep_variant_tile_append_event(
                  tile, 0u, 10u, 20u, (uint8_t)DUCKVEP_KIND_SV,
                  (uint8_t)DUCKVEP_SV_INVERSION,
                  (uint8_t)DUCKVEP_COPY_CHANGE_GAIN,
                  &r, 1u, &a, 1u, "v", 1u));
    /* ref_len != 1 */
    ASSERT_EQ(DUCKVEP_TILE_INVALID,
              duckvep_variant_tile_append(tile, 0u, 10u, 10u, (uint8_t)DUCKVEP_KIND_SNV,
                                          "AC", 2u, &a, 1u, "v", 1u));
    /* end1 != pos1 */
    ASSERT_EQ(DUCKVEP_TILE_INVALID,
              duckvep_variant_tile_append(tile, 0u, 10u, 11u, (uint8_t)DUCKVEP_KIND_SNV,
                                          &r, 1u, &a, 1u, "v", 1u));
    /* malformed MNV: unequal lengths and bad end coordinate */
    ASSERT_EQ(DUCKVEP_TILE_INVALID,
              duckvep_variant_tile_append(tile, 0u, 10u, 11u, (uint8_t)DUCKVEP_KIND_MNV,
                                          "AC", 2u, "T", 1u, "v", 1u));
    ASSERT_EQ(DUCKVEP_TILE_INVALID,
              duckvep_variant_tile_append(tile, 0u, 10u, 10u, (uint8_t)DUCKVEP_KIND_MNV,
                                          "AC", 2u, "GT", 2u, "v", 1u));
    /* malformed INS/DEL/INDEL anchor shapes */
    ASSERT_EQ(DUCKVEP_TILE_INVALID,
              duckvep_variant_tile_append(tile, 0u, 10u, 11u, (uint8_t)DUCKVEP_KIND_INS,
                                          "A", 1u, "AT", 2u, "v", 1u));
    ASSERT_EQ(DUCKVEP_TILE_INVALID,
              duckvep_variant_tile_append(tile, 0u, 10u, 11u, (uint8_t)DUCKVEP_KIND_DEL,
                                          "AC", 2u, "GT", 2u, "v", 1u));
    ASSERT_EQ(DUCKVEP_TILE_INVALID,
              duckvep_variant_tile_append(tile, 0u, 10u, 11u, (uint8_t)DUCKVEP_KIND_INDEL,
                                          "AC", 2u, "GT", 2u, "v", 1u));
    ASSERT_EQ(DUCKVEP_TILE_INVALID,
              duckvep_variant_tile_append(tile, 0u, 10u, 10u, (uint8_t)DUCKVEP_KIND_INDEL,
                                          "AC", 2u, "G", 1u, "v", 1u));
    /* null required field */
    ASSERT_EQ(DUCKVEP_TILE_INVALID,
              duckvep_variant_tile_append(tile, 0u, 10u, 10u, (uint8_t)DUCKVEP_KIND_SNV,
                                          NULL, 1u, &a, 1u, "v", 1u));
    /* chrom_id outside interned space */
    ASSERT_EQ(DUCKVEP_TILE_INVALID, tile_put_snv(tile, 4u, 10u, 'A', 'C', "v"));
    ASSERT_EQ(0u, duckvep_variant_tile_count(tile)); /* nothing stored */
    duckvep_variant_tile_destroy(tile);
    PASS();
}

TEST tile_appends_and_batch_view_roundtrips(void) {
    duckvep_variant_tile_t *tile = duckvep_variant_tile_create(8u, 4u);
    duckvep_variant_batch_t v;
    size_t idlen = 0u;
    const char *id;
    ASSERT(tile != NULL);

    ASSERT_EQ(DUCKVEP_TILE_APPENDED, tile_put_snv(tile, 1u, 100u, 'A', 'G', "rs1"));
    ASSERT_EQ(DUCKVEP_TILE_APPENDED, tile_put_snv(tile, 1u, 100u, 'C', 'T', "rs2")); /* equal pos OK */
    ASSERT_EQ(DUCKVEP_TILE_APPENDED, tile_put_mnv(tile, 1u, 205u, "AA", "GG", 2u, "rs3"));
    ASSERT_EQ(DUCKVEP_TILE_APPENDED, tile_put_ins(tile, 1u, 210u, "A", "AT", 2u, "rs4"));
    ASSERT_EQ(DUCKVEP_TILE_APPENDED, tile_put_del(tile, 1u, 220u, "AC", 2u, "A", "rs5"));
    ASSERT_EQ(DUCKVEP_TILE_APPENDED, tile_put_indel(tile, 1u, 230u, "AC", 2u, "G", 1u, "rs6"));
    ASSERT_EQ(6u, duckvep_variant_tile_count(tile));

    duckvep_variant_tile_batch(tile, &v);
    ASSERT_EQ(6u, v.count);
    ASSERT_EQ(1u, v.chrom_id[0]);
    ASSERT_EQ(100u, v.pos1[0]);
    ASSERT_EQ(205u, v.pos1[2]);
    ASSERT_EQ(206u, v.end1[2]);
    ASSERT_EQ(210u, v.pos1[3]);
    ASSERT_EQ(210u, v.end1[3]);
    ASSERT_EQ(221u, v.end1[4]);
    ASSERT_EQ((uint8_t)DUCKVEP_KIND_SNV, v.variant_kind[1]);
    ASSERT_EQ((uint8_t)DUCKVEP_KIND_MNV, v.variant_kind[2]);
    ASSERT_EQ((uint8_t)DUCKVEP_KIND_INS, v.variant_kind[3]);
    ASSERT_EQ((uint8_t)DUCKVEP_KIND_DEL, v.variant_kind[4]);
    ASSERT_EQ((uint8_t)DUCKVEP_KIND_INDEL, v.variant_kind[5]);
    /* allele pool: ref then alt per row */
    ASSERT_EQ('A', (char)v.allele_bytes[v.ref_offset[0]]);
    ASSERT_EQ('G', (char)v.allele_bytes[v.alt_offset[0]]);
    ASSERT_EQ('T', (char)v.allele_bytes[v.alt_offset[1]]);
    ASSERT_EQ('A', (char)v.allele_bytes[v.ref_offset[2]]);
    ASSERT_EQ('G', (char)v.allele_bytes[v.alt_offset[2]]);
    ASSERT_EQ(2u, v.ref_length[2]);
    ASSERT_EQ(2u, v.alt_length[2]);
    ASSERT_EQ(1u, v.ref_length[3]);
    ASSERT_EQ(2u, v.alt_length[3]);
    ASSERT_EQ(2u, v.ref_length[4]);
    ASSERT_EQ(1u, v.alt_length[4]);
    ASSERT_EQ(2u, v.ref_length[5]);
    ASSERT_EQ(1u, v.alt_length[5]);
    /* variant_id metadata indexed by variant_idx */
    id = duckvep_variant_tile_variant_id(tile, 5u, &idlen);
    ASSERT_EQ(3u, idlen);
    ASSERT(id != NULL && memcmp(id, "rs6", 3u) == 0);
    ASSERT(duckvep_variant_tile_variant_id(tile, 6u, &idlen) == NULL); /* out of range */
    duckvep_variant_tile_destroy(tile);
    PASS();
}

TEST tile_preserves_structural_operation_and_copy_direction(void) {
    duckvep_variant_tile_t *tile = duckvep_variant_tile_create(4u, 2u);
    duckvep_variant_batch_t v;
    ASSERT(tile != NULL);
    ASSERT_EQ(
        DUCKVEP_TILE_APPENDED,
        duckvep_variant_tile_append_event(
            tile, 0u, 100u, 400u, (uint8_t)DUCKVEP_KIND_SV,
            (uint8_t)DUCKVEP_SV_CNV, (uint8_t)DUCKVEP_COPY_CHANGE_GAIN,
            "N", 1u, "<CNV>", 5u, "cnv1", 4u));
    ASSERT_EQ(
        DUCKVEP_TILE_APPENDED,
        duckvep_variant_tile_append_event(
            tile, 0u, 500u, 700u, (uint8_t)DUCKVEP_KIND_SV,
            (uint8_t)DUCKVEP_SV_INVERSION,
            (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN,
            "N", 1u, "<INV>", 5u, "inv1", 4u));
    duckvep_variant_tile_batch(tile, &v);
    ASSERT_EQ(2u, v.count);
    ASSERT(v.sv_type != NULL);
    ASSERT(v.copy_change != NULL);
    ASSERT_EQ((uint8_t)DUCKVEP_SV_CNV, v.sv_type[0]);
    ASSERT_EQ((uint8_t)DUCKVEP_COPY_CHANGE_GAIN, v.copy_change[0]);
    ASSERT_EQ((uint8_t)DUCKVEP_SV_INVERSION, v.sv_type[1]);
    ASSERT_EQ((uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN, v.copy_change[1]);
    duckvep_variant_tile_destroy(tile);
    PASS();
}

TEST tile_rejects_pos_decrease_within_run(void) {
    duckvep_variant_tile_t *tile = duckvep_variant_tile_create(8u, 4u);
    ASSERT(tile != NULL);
    ASSERT_EQ(DUCKVEP_TILE_APPENDED, tile_put_snv(tile, 2u, 500u, 'A', 'C', "a"));
    ASSERT_EQ(DUCKVEP_TILE_INVALID, tile_put_snv(tile, 2u, 499u, 'A', 'C', "b"));
    ASSERT(duckvep_variant_tile_error(tile) != NULL);
    ASSERT_EQ(1u, duckvep_variant_tile_count(tile)); /* the bad row not stored */
    duckvep_variant_tile_destroy(tile);
    PASS();
}

TEST tile_chrom_change_flushes_closes_and_blocks_reappearance(void) {
    duckvep_variant_tile_t *tile = duckvep_variant_tile_create(8u, 4u);
    ASSERT(tile != NULL);
    ASSERT_EQ(DUCKVEP_TILE_APPENDED, tile_put_snv(tile, 1u, 10u, 'A', 'C', "a"));

    /* moving to a new chrom signals FULL_CHROM (flush current tile first) */
    ASSERT_EQ(DUCKVEP_TILE_FULL_CHROM, tile_put_snv(tile, 2u, 5u, 'A', 'C', "b"));
    duckvep_variant_tile_close_run(tile);          /* closes chrom 1 */
    ASSERT_EQ(0u, duckvep_variant_tile_count(tile));
    ASSERT_EQ(DUCKVEP_TILE_APPENDED, tile_put_snv(tile, 2u, 5u, 'A', 'C', "b")); /* retry opens chrom 2 */

    /* chrom 1 reappearing after its run closed is caught immediately as a hard
     * error, even mid-run for chrom 2 — no flush/retry dance needed. */
    ASSERT_EQ(DUCKVEP_TILE_INVALID, tile_put_snv(tile, 1u, 99u, 'A', 'C', "c"));
    ASSERT(duckvep_variant_tile_error(tile) != NULL);
    duckvep_variant_tile_destroy(tile);
    PASS();
}

TEST tile_capacity_full_then_reset_preserves_run(void) {
    duckvep_variant_tile_t *tile = duckvep_variant_tile_create(2u, 4u);
    ASSERT(tile != NULL);
    ASSERT_EQ(DUCKVEP_TILE_APPENDED, tile_put_snv(tile, 0u, 10u, 'A', 'C', "a"));
    ASSERT_EQ(DUCKVEP_TILE_APPENDED, tile_put_snv(tile, 0u, 20u, 'A', 'C', "b"));
    /* third row on a cap-2 tile: flush signal, same chrom continues */
    ASSERT_EQ(DUCKVEP_TILE_FULL_CAPACITY, tile_put_snv(tile, 0u, 30u, 'A', 'C', "c"));
    duckvep_variant_tile_reset_rows(tile);
    ASSERT_EQ(0u, duckvep_variant_tile_count(tile));
    ASSERT_EQ(DUCKVEP_TILE_APPENDED, tile_put_snv(tile, 0u, 30u, 'A', 'C', "c"));
    /* monotonicity is enforced ACROSS the tile boundary (last_pos1 preserved) */
    ASSERT_EQ(DUCKVEP_TILE_INVALID, tile_put_snv(tile, 0u, 25u, 'A', 'C', "d"));
    duckvep_variant_tile_destroy(tile);
    PASS();
}

/* ---- property: a sorted stream survives the controller protocol intact ---- */
#define TILEPROP_MAX_EVENTS 64u
#define TILEPROP_MAX_RUNS   8u
struct tileprop_stream {
    size_t   n;
    uint16_t chrom[TILEPROP_MAX_EVENTS];
    uint32_t pos[TILEPROP_MAX_EVENTS];
    char     ref[TILEPROP_MAX_EVENTS];
    char     alt[TILEPROP_MAX_EVENTS];
    uint16_t nruns;
    uint32_t cap; /* tile capacity to use */
};

static enum theft_alloc_res tileprop_alloc(struct theft *t, void *env, void **instance) {
    struct tileprop_stream *s = (struct tileprop_stream *)calloc(1u, sizeof *s);
    static const char bases[4] = {'A', 'C', 'G', 'T'};
    size_t i;
    uint16_t run = 0u;
    uint32_t pos = 0u;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->n = (size_t)kprop_bounded(t, (uint64_t)TILEPROP_MAX_EVENTS + 1u);
    s->cap = (uint32_t)kprop_bounded(t, 16u) + 1u;
    /* Build a SORTED stream: contiguous chrom runs with strictly increasing chrom
     * ids (so a contig never reappears), pos RESETS only when the chrom changes and
     * is otherwise nondecreasing within the run. */
    for (i = 0; i < s->n; i++) {
        uint8_t ref_idx;
        uint8_t alt_step;
        if (i > 0u && run + 1u < TILEPROP_MAX_RUNS && theft_random_bits(t, 2) == 0u) {
            run = (uint16_t)(run + 1u); /* start a new contig */
            pos = 1u;
        } else {
            pos += (uint32_t)kprop_bounded(t, 8u); /* nondecreasing within the run */
            if (pos == 0u) pos = 1u;               /* keep pos1 >= 1 */
        }
        s->chrom[i] = run;
        s->pos[i] = pos;
        ref_idx = (uint8_t)theft_random_bits(t, 2);
        alt_step = (uint8_t)kprop_bounded(t, 3u) + 1u;
        s->ref[i] = bases[ref_idx];
        /* The tile rejects REF==ALT because it is not a variant. Choose one of
         * the other three bases while retaining a uniformly varied SNV stream. */
        s->alt[i] = bases[(ref_idx + alt_step) & 3u];
    }
    s->nruns = (uint16_t)(run + 1u);
    *instance = s;
    return THEFT_ALLOC_OK;
}

static void tileprop_free(void *instance, void *env) {
    (void)env;
    free(instance);
}

static struct theft_type_info tileprop_info = {
    .alloc = tileprop_alloc,
    .free  = tileprop_free,
};

/* Drain the tile into the collected[] arrays, asserting each flushed batch is a
 * single chrom with nondecreasing pos. Returns 0 on an invariant violation. */
static int tileprop_drain(const duckvep_variant_tile_t *tile, size_t *out_n,
                          uint16_t *cc, uint32_t *cp) {
    duckvep_variant_batch_t v;
    size_t j;
    duckvep_variant_tile_batch(tile, &v);
    for (j = 0; j < v.count; j++) {
        if (j > 0u) {
            if (v.chrom_id[j] != v.chrom_id[0]) return 0;       /* one chrom per tile */
            if (v.pos1[j] < v.pos1[j - 1u]) return 0;           /* nondecreasing */
        }
        cc[*out_n] = v.chrom_id[j];
        cp[*out_n] = v.pos1[j];
        (*out_n)++;
    }
    return 1;
}

static enum theft_trial_res prop_tile_controller_preserves_stream(struct theft *t, void *arg1) {
    const struct tileprop_stream *s = (const struct tileprop_stream *)arg1;
    duckvep_variant_tile_t *tile;
    uint16_t cc[TILEPROP_MAX_EVENTS];
    uint32_t cp[TILEPROP_MAX_EVENTS];
    size_t collected = 0u;
    size_t i;
    (void)t;

    tile = duckvep_variant_tile_create((size_t)s->cap, (size_t)s->nruns + 1u);
    if (tile == NULL) return THEFT_TRIAL_ERROR;

    for (i = 0; i < s->n; i++) {
        for (;;) { /* retry after a flush */
            duckvep_tile_status_t st = duckvep_variant_tile_append(
                tile, s->chrom[i], s->pos[i], s->pos[i], (uint8_t)DUCKVEP_KIND_SNV,
                &s->ref[i], 1u, &s->alt[i], 1u, "x", 1u);
            if (st == DUCKVEP_TILE_APPENDED) break;
            if (st == DUCKVEP_TILE_INVALID) { /* a sorted stream must never be invalid */
                duckvep_variant_tile_destroy(tile);
                return THEFT_TRIAL_FAIL;
            }
            if (!tileprop_drain(tile, &collected, cc, cp)) {
                duckvep_variant_tile_destroy(tile);
                return THEFT_TRIAL_FAIL;
            }
            if (st == DUCKVEP_TILE_FULL_CHROM) duckvep_variant_tile_close_run(tile);
            else duckvep_variant_tile_reset_rows(tile);
        }
    }
    if (!tileprop_drain(tile, &collected, cc, cp)) {
        duckvep_variant_tile_destroy(tile);
        return THEFT_TRIAL_FAIL;
    }
    duckvep_variant_tile_destroy(tile);

    /* every input row emitted exactly once, in order, unchanged */
    if (collected != s->n) return THEFT_TRIAL_FAIL;
    for (i = 0; i < s->n; i++) {
        if (cc[i] != s->chrom[i] || cp[i] != s->pos[i]) return THEFT_TRIAL_FAIL;
    }
    return THEFT_TRIAL_PASS;
}

TEST tile_controller_preserves_sorted_stream_for_any_input(void) {
    struct theft_run_config cfg;
    static struct theft_type_info *arg_infos[1];
    arg_infos[0] = &tileprop_info;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "tile_controller_preserves_sorted_stream";
    cfg.prop1 = prop_tile_controller_preserves_stream;
    cfg.type_info[0] = arg_infos[0];
    cfg.trials = 500;
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

GREATEST_MAIN_DEFS();

int main(int argc, char **argv) {
    GREATEST_MAIN_BEGIN();
    RUN_TEST(kernel_version_is_well_formed);
    RUN_TEST(null_args_are_rejected);
    RUN_TEST(model_open_rejects_projection_and_sequence_invariant_mutations);
    RUN_TEST(annotate_tile_rejects_null_model_for_any_batch);
    RUN_TEST(event_load_without_variant_kind_uses_supplied_interval);
    RUN_TEST(event_load_trims_small_variant_differing_region);
    RUN_TEST(event_prepare_small_preserves_anchor_side_and_semantic_kind);
    RUN_TEST(event_normalization_matches_trim_oracle);
    RUN_TEST(sweep_vep_feature_span_candidates_match_oracle);
    RUN_TEST(sweep_known_scene_exact_pairs);
    RUN_TEST(sweep_span_tail_does_not_poison_point_frontier);
    RUN_TEST(sweep_uses_full_span_only_for_structural_events);
    RUN_TEST(sweep_small_variant_differing_region_tail_is_not_persistent);
    RUN_TEST(sweep_sink_stop_is_immediate);
    RUN_TEST(sweep_saturation_near_uint32_max);
    RUN_TEST(sweep_active_cap_overflow_is_reported);
    RUN_TEST(sweep_matches_bruteforce_for_any_scene);
    RUN_TEST(seeded_sweep_matches_bruteforce_for_any_scene);
    RUN_TEST(splice_classify_known_scene);
    RUN_TEST(splice_classify_insertion_interbase_scene);
    RUN_TEST(splice_ppt_exon_gate_scene);
    RUN_TEST(effect_rule_table_known_pre_bits);
    RUN_TEST(nmd_transcript_predicate_known_scene);
    RUN_TEST(variant_induced_nmd_prediction_known_scene);
    RUN_TEST(effect_rule_tiers_suppress_only_later_tiers);
    RUN_TEST(event_length_delta_pre_bits_follow_trimmed_alleles);
    RUN_TEST(annotate_structural_known_scene);
    RUN_TEST(annotate_padded_small_variants_use_differing_region_topology);
    RUN_TEST(annotate_sv_cnv_known_scene);
    RUN_TEST(annotate_complete_neutral_sv_uses_transcript_fallbacks);
    RUN_TEST(annotate_region_mask_truthful_known_scene);
    RUN_TEST(sv_metadata_validity_matrix_known);
    RUN_TEST(sv_predicate_facts_known);
    RUN_TEST(so_render_and_impact_name_known);
    RUN_TEST(annotate_reverse_strand_known_scene);
    RUN_TEST(annotate_codon_snv_known_scene);
    RUN_TEST(annotate_codon_reverse_strand_known_scene);
    RUN_TEST(annotate_codon_ref_mismatch_falls_back);
    RUN_TEST(annotate_codon_stop_retained_known_scene);
    RUN_TEST(annotate_codon_start_lost_known_scene);
    RUN_TEST(annotate_codon_multi_exon_known_scene);
    RUN_TEST(annotate_codon_mnv_same_codon_known_scene);
    RUN_TEST(annotate_codon_mnv_start_lost_route_known_scene);
    RUN_TEST(annotate_cursor_mnv_start_lost_route_matches_tile_known_scene);
    RUN_TEST(annotate_codon_mnv_len3_and_cross_codon_known_scene);
    RUN_TEST(annotate_codon_mnv_cross_codon_missense_known_scene);
    RUN_TEST(annotate_codon_mnv_reverse_strand_same_codon_known_scene);
    RUN_TEST(annotate_equal_length_cds_to_utr3_mapping_gap_known_scene);
    RUN_TEST(annotate_codon_padded_small_variant_delta_known_scene);
    RUN_TEST(annotate_equal_length_feature_window_both_strands_known_scene);
    RUN_TEST(annotate_codon_indel_frameshift_known_scene);
    RUN_TEST(annotate_terminal_stop_frame_change_uses_transcript_tail);
    RUN_TEST(annotate_codon_delins_frameshift_known_scene);
    RUN_TEST(annotate_codon_inframe_insertion_known_scene);
    RUN_TEST(annotate_codon_inframe_insertion_reverse_known_scene);
    RUN_TEST(annotate_codon_protein_altering_insertion_known_scene);
    RUN_TEST(annotate_codon_indel_reverse_ref_mismatch_falls_back);
    RUN_TEST(annotate_codon_matches_kernel_for_any_cds_snv);
    RUN_TEST(annotate_mnv_same_codon_matches_oracle);
    RUN_TEST(annotate_cross_codon_mnv_missense_matches_oracle);
    RUN_TEST(cds_edit_builder_projects_exon_boundary_insertions_on_both_strands);
    RUN_TEST(cds_edit_builder_known_scene);
    RUN_TEST(coding_context_known_scene);
    RUN_TEST(variant_coding_context_known_scene);
    RUN_TEST(coding_context_delta_known_scene);
    RUN_TEST(coding_context_delta_inframe_deletion_known_scene);
    RUN_TEST(coding_context_delta_inframe_insertion_known_scene);
    RUN_TEST(coding_context_delta_delins_known_scene);
    RUN_TEST(coding_context_delta_frameshift_known_scene);
    RUN_TEST(coding_context_delta_frameshift_stop_gained_scene);
    RUN_TEST(sequence_delta_with_scratch_indel_known_scene);
    RUN_TEST(annotate_delins_boundary_no_route_known_scene);
    RUN_TEST(annotate_inframe_insertion_route_known_scene);
    RUN_TEST(annotate_inframe_deletion_route_known_scene);
    RUN_TEST(cds_edit_builder_matches_direct_splice_oracle);
    RUN_TEST(cds_edit_set_builder_matches_single_edit_oracle);
    RUN_TEST(cds_edit_set_builder_splits_mnv_diff_islands);
    RUN_TEST(coding_context_matches_direct_oracles);
    RUN_TEST(variant_coding_context_matches_oracles);
    RUN_TEST(coding_context_delta_matches_codon_oracle);
    RUN_TEST(coding_context_delta_frameshift_matches_length_oracle);
    RUN_TEST(workspace_delta_scratch_caps_known);
    RUN_TEST(workspace_delta_scratch_builds_lengthening_context);
    RUN_TEST(workspace_delta_scratch_usable_for_mnv);
    RUN_TEST(sequence_delta_scratch_rejects_unequal_mnv_kind);
    RUN_TEST(sequence_delta_annotation_wrapper_del_insufficient_scratch_known);
    RUN_TEST(sequence_delta_annotation_wrapper_start_lost_mnv);
    RUN_TEST(sequence_delta_with_scratch_cross_codon_known_scene);
    RUN_TEST(sequence_delta_with_scratch_cross_codon_reverse_known_scene);
    RUN_TEST(sequence_delta_with_scratch_cross_codon_negative_scenes);
    RUN_TEST(feature_substitution_window_fails_closed_known_scene);
    RUN_TEST(sequence_delta_annotation_wrapper_matches_direct_shape);
    RUN_TEST(annotate_cursor_padded_snv_matches_tile_for_any_output_split);
    RUN_TEST(annotate_cursor_cross_codon_mnv_route_matches_tile_for_any_output_split);
    RUN_TEST(sequence_delta_with_scratch_mnv_matches_oracle);
    RUN_TEST(sequence_delta_with_scratch_cross_codon_mnv_matches_oracle);
    RUN_TEST(annotate_frameshift_indel_matches_oracle);
    RUN_TEST(annotate_inframe_deletion_matches_oracle);
    RUN_TEST(annotate_cursor_del_route_matches_tile_for_any_output_split);
    RUN_TEST(coding_context_delta_inframe_deletion_matches_oracle);
    RUN_TEST(coding_context_delta_inframe_insertion_matches_oracle);
    RUN_TEST(coding_context_delta_delins_shape_matches_oracle);
    RUN_TEST(sequence_delta_with_scratch_indel_matches_oracle);
    RUN_TEST(annotate_inframe_insertion_matches_oracle);
    RUN_TEST(annotate_cursor_ins_route_matches_tile_for_any_output_split);
    RUN_TEST(annotate_protein_altering_insertion_matches_oracle);
    RUN_TEST(annotate_start_lost_matches_oracle_for_any_start_codon_snv);
    RUN_TEST(annotate_rejects_unsorted_variant_batch);
    RUN_TEST(annotate_rejects_missing_alleles_for_nonpoint_small_variant);
    RUN_TEST(annotate_accepts_right_anchored_pure_insertion);
    RUN_TEST(annotate_rejects_kind_allele_shape_mismatch);
    RUN_TEST(annotate_rejects_allele_slice_outside_pool);
    RUN_TEST(annotate_result_full_is_reported);
    RUN_TEST(annotate_rejects_result_count_past_capacity);
    RUN_TEST(annotate_directional_distance_uses_u32);
    RUN_TEST(annotate_directional_distance_filter);
    RUN_TEST(model_open_rejects_invalid_models);
    RUN_TEST(annotate_matches_composition_for_any_scene);
    RUN_TEST(annotate_cursor_resumes_known_scene);
    RUN_TEST(sorted_point_cursor_survives_tiles_and_resets_on_rewind);
    RUN_TEST(padded_snv_rewind_uses_vep_feature_span);
    RUN_TEST(annotate_cursor_matches_tile_for_any_output_split);
    RUN_TEST(region_mask_known_scene);
    RUN_TEST(region_mask_short_exon_splice_reach_clamped);
    RUN_TEST(region_span_can_cross_cds_intron_and_splice_windows);
    RUN_TEST(region_mask_invariants_hold);
    RUN_TEST(sorted_point_classifier_matches_exhaustive_for_any_transcript);
    RUN_TEST(projection_known_forward_reverse_and_phase);
    RUN_TEST(projection_matches_bruteforce_for_any_small_transcript);
    RUN_TEST(codon_translate_known);
    RUN_TEST(codon_change_known);
    RUN_TEST(codon_change_consistent_with_translation);
    RUN_TEST(coding_snv_from_cds_known_cases);
    RUN_TEST(coding_snv_from_cds_matches_oracle_for_any_valid_snv);
    RUN_TEST(haplotype_partition_known_cases);
    RUN_TEST(haplotype_partition_preserves_interactions_for_any_valid_edit_set);
    RUN_TEST(haplotype_apply_and_translate_known_cases);
    RUN_TEST(haplotype_apply_matches_rebuild_oracle_for_any_valid_edit_set);
    RUN_TEST(tile_rejects_breakends_and_malformed_small_variants);
    RUN_TEST(tile_appends_and_batch_view_roundtrips);
    RUN_TEST(tile_preserves_structural_operation_and_copy_direction);
    RUN_TEST(tile_rejects_pos_decrease_within_run);
    RUN_TEST(tile_chrom_change_flushes_closes_and_blocks_reappearance);
    RUN_TEST(tile_capacity_full_then_reset_preserves_run);
    RUN_TEST(tile_controller_preserves_sorted_stream_for_any_input);
    GREATEST_MAIN_END();
}
