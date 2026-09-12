#include "duckvep_property.h"

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
    uint32_t first_stop_position1[1] = {3u};
    uint32_t peptide_edit_offset[2] = {0u, 1u};
    uint32_t peptide_edit_position1[1] = {2u};
    uint8_t peptide_edit_alt[1] = {(uint8_t)'U'};
    uint8_t invalid_tail[2] = {'A','X'};
    uint64_t pre_cds_offset[1] = {0u};
    uint32_t pre_cds_length[1] = {0u};
    uint64_t post_cds_offset[1] = {0u};
    uint32_t post_cds_length[1] = {2u};
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

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    duckvep_model_close(model); model = NULL;

    seq.first_stop_position1 = first_stop_position1;
    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    duckvep_model_close(model); model = NULL;

    first_stop_position1[0] = 2u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(71u, err.where_code);
    first_stop_position1[0] = 3u;

    seq.peptide_edit_offset = peptide_edit_offset;
    seq.peptide_edit_position1 = peptide_edit_position1;
    seq.peptide_edit_alt = peptide_edit_alt;
    seq.peptide_edit_count = 1u;
    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    duckvep_model_close(model); model = NULL;

    peptide_edit_position1[0] = 4u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(73u, err.where_code);
    peptide_edit_position1[0] = 2u;

    peptide_edit_alt[0] = (uint8_t)'?';
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(73u, err.where_code);
    peptide_edit_alt[0] = (uint8_t)'U';

    strand[0] = 0;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(66u, err.where_code); strand[0] = 1;

    tx_start[0] = 99u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(67u, err.where_code); tx_start[0] = 100u;

    cdna_start[1] = 7u; cdna_end[1] = 10u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(68u, err.where_code); cdna_start[1] = 6u; cdna_end[1] = 9u;

    phase[0] = 3;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(69u, err.where_code); phase[0] = 0;

    cds_start[0] = 105u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(70u, err.where_code); cds_start[0] = 100u;

    sequence_length[0] = 8u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(71u, err.where_code); sequence_length[0] = 9u;

    codon_table[0] = 8u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(71u, err.where_code); codon_table[0] = 1u;

    cds_bytes[4] = 'X';
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(71u, err.where_code); cds_bytes[4] = 'A';

    seq.flank_bytes = invalid_tail; seq.flank_bytes_len = sizeof invalid_tail;
    seq.pre_cds_offset = pre_cds_offset; seq.pre_cds_length = pre_cds_length;
    seq.post_cds_offset = post_cds_offset; seq.post_cds_length = post_cds_length;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
    ASSERT_EQ(71u, err.where_code);
    seq.flank_bytes = NULL; seq.flank_bytes_len = 0u;
    seq.pre_cds_offset = NULL; seq.pre_cds_length = NULL;
    seq.post_cds_offset = NULL; seq.post_cds_length = NULL;

    exons.phase = NULL; exons.end_phase = NULL;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, &seq, NULL, &model, &err));
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
    static const uint32_t halos[11] = {
        0u, 1u, 50u, 100u, 4999u, 5000u, 5001u, 10000u, 50000u,
        65535u, UINT32_MAX
    };
    uint32_t base = (uint32_t)kprop_bounded(t, 0xFFFF0000u) + 1u;
    uint32_t span = (uint32_t)kprop_bounded(t, 20000u) + 1u;
    (void)env;
    if (s == NULL) return THEFT_ALLOC_ERROR;
    s->halo = halos[kprop_bounded(t, 11u)];

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

struct theft_type_info kprop_scene_info = {
    .alloc = kprop_scene_alloc,
    .free  = kprop_scene_free,
};

int pair_sink(uint32_t vi, uint32_t ti, void *ctx) {
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

    /* Candidate ordinals are part of the streaming contract: vector/tile
     * restarts must not permute one variant's result list. */
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
TEST breakend_parser_checks_shapes_and_limits(void) {
    const char *invalid[] = {"", "[chr:1]A", "A[chr:1]", "A[chr:1[A",
        "[chr:1[", "A[:1[", "A[chr:[", "A[chr:-1[", "A[chr:+1[",
        "A[chr:1x[", "A[chr:1[[", "A[chr 1:2[", "A[chr,1:2[", "A[",
        "A..", "..A", "A.T.", "U[chr:1[", "[chr:1[R"};
    duckvep_breakend_t out, zero = {0};
    for (size_t i = 0u; i < sizeof invalid / sizeof invalid[0]; i++) {
        memset(&out, 0xa5, sizeof out);
        ASSERT_EQ(DUCKVEP_BREAKEND_INVALID, duckvep_breakend_parse(
            (const uint8_t *)invalid[i], strlen(invalid[i]), &out));
        ASSERT_MEM_EQ(&zero, &out, sizeof out);
    }
    const char *ordinary[] = {".", "*", "<DUP>", "<ctg.1>", "ACGTN"};
    for (size_t i = 0u; i < sizeof ordinary / sizeof ordinary[0]; i++) {
        ASSERT_EQ(DUCKVEP_BREAKEND_NOT_BREAKEND, duckvep_breakend_parse(
            (const uint8_t *)ordinary[i], strlen(ordinary[i]), &out));
        ASSERT_MEM_EQ(&zero, &out, sizeof out);
    }
    const char maximum[] = "a[chr:with:semicolon;:18446744073709551615[";
    ASSERT_EQ(DUCKVEP_BREAKEND_OK, duckvep_breakend_parse(
        (const uint8_t *)maximum, sizeof maximum - 1u, &out));
    ASSERT_EQ(UINT64_MAX, out.mate_position);
    ASSERT_EQ(19u, out.mate_chrom_length);
    ASSERT_MEM_EQ("chr:with:semicolon;", out.mate_chrom, 19u);
    const char overflow[] = "a[chr:18446744073709551616[";
    ASSERT_EQ(DUCKVEP_BREAKEND_POSITION_OVERFLOW, duckvep_breakend_parse(
        (const uint8_t *)overflow, sizeof overflow - 1u, &out));
    ASSERT_MEM_EQ(&zero, &out, sizeof out);
    const uint8_t embedded_nul[] = {'A','[','c',0,':','1','['};
    ASSERT_EQ(DUCKVEP_BREAKEND_INVALID, duckvep_breakend_parse(
        embedded_nul, sizeof embedded_nul, &out));
    const char telomere[] = "]<ctg:1>:000]tN";
    ASSERT_EQ(DUCKVEP_BREAKEND_OK, duckvep_breakend_parse(
        (const uint8_t *)telomere, sizeof telomere - 1u, &out));
    ASSERT_EQ(0u, out.mate_position);
    ASSERT_EQ(7u, out.mate_chrom_length);
    ASSERT_MEM_EQ("<ctg:1>", out.mate_chrom, 7u);
    ASSERT_EQ(0u, out.local_join_after);
    ASSERT_EQ(0u, out.mate_extends_right);
    ASSERT_EQ(2u, out.replacement_length);
    ASSERT_MEM_EQ("tN", out.replacement, 2u);
    ASSERT_EQ(DUCKVEP_BREAKEND_INVALID, duckvep_breakend_parse(NULL, 1u, &out));
    ASSERT_EQ(DUCKVEP_BREAKEND_INVALID, duckvep_breakend_parse(embedded_nul, 1u, NULL));
    PASS();
}

struct kprop_breakend {
    char chrom[40], replacement[40];
    uint64_t position;
    unsigned form;
};

static enum theft_alloc_res kprop_breakend_alloc(struct theft *t, void *env, void **instance) {
    static const char bases[] = "ACGTNacgtn";
    static const char name_bytes[] = "ACGT0123456789:;._-";
    struct kprop_breakend *b = calloc(1u, sizeof *b);
    (void)env;
    if (b == NULL) return THEFT_ALLOC_ERROR;
    size_t n = (size_t)kprop_bounded(t, 32u) + 1u;
    for (size_t i = 0u; i < n; i++)
        b->replacement[i] = bases[kprop_bounded(t, sizeof bases - 1u)];
    n = (size_t)kprop_bounded(t, 32u) + 1u;
    for (size_t i = 0u; i < n; i++)
        b->chrom[i] = name_bytes[kprop_bounded(t, sizeof name_bytes - 1u)];
    b->position = theft_random_bits(t, 64u);
    b->form = (unsigned)kprop_bounded(t, 6u);
    *instance = b;
    return THEFT_ALLOC_OK;
}

static int kprop_breakend_render(const struct kprop_breakend *b, char *rendered, size_t capacity) {
    int after = b->form < 2u || b->form == 4u;
    int right = b->form == 0u || b->form == 3u;
    if (b->form >= 4u)
        return snprintf(rendered, capacity, after ? "%s." : ".%s", b->replacement);
    char bracket = right ? '[' : ']';
    return after
        ? snprintf(rendered, capacity, "%s%c%s:%" PRIu64 "%c",
            b->replacement, bracket, b->chrom, b->position, bracket)
        : snprintf(rendered, capacity, "%c%s:%" PRIu64 "%c%s",
            bracket, b->chrom, b->position, bracket, b->replacement);
}

static enum theft_trial_res prop_breakend_recovers_constructed_components(struct theft *t, void *arg) {
    const struct kprop_breakend *b = arg;
    char rendered[128];
    int length = kprop_breakend_render(b, rendered, sizeof rendered);
    int after = b->form < 2u || b->form == 4u;
    int right = b->form == 0u || b->form == 3u;
    int paired = b->form < 4u;
    (void)t;
    if (length < 1 || (size_t)length >= sizeof rendered) return THEFT_TRIAL_ERROR;
    /* Exact, unterminated allocation catches accidental strlen/strchr reads. */
    uint8_t *bytes = malloc((size_t)length);
    if (bytes == NULL) return THEFT_TRIAL_ERROR;
    memcpy(bytes, rendered, (size_t)length);
    duckvep_breakend_t parsed;
    int ok = duckvep_breakend_parse(bytes, (size_t)length, &parsed) == DUCKVEP_BREAKEND_OK;
    if (ok) {
        ok = parsed.local_join_after == after && parsed.has_mate == paired &&
            parsed.replacement_length == strlen(b->replacement) &&
            memcmp(parsed.replacement, b->replacement, parsed.replacement_length) == 0;
        if (paired) ok = ok && parsed.mate_extends_right == right &&
            parsed.mate_position == b->position && parsed.mate_chrom_length == strlen(b->chrom) &&
            memcmp(parsed.mate_chrom, b->chrom, parsed.mate_chrom_length) == 0;
        else ok = ok && parsed.mate_chrom == NULL && parsed.mate_chrom_length == 0u;
    }
    free(bytes);
    return ok ? THEFT_TRIAL_PASS : THEFT_TRIAL_FAIL;
}

TEST breakend_parser_recovers_constructed_components(void) {
    struct theft_type_info type = {.alloc = kprop_breakend_alloc, .free = theft_generic_free_cb};
    struct theft_run_config cfg = {0};
    cfg.name = "breakend_parser_recovers_constructed_components";
    cfg.prop1 = prop_breakend_recovers_constructed_components;
    cfg.type_info[0] = &type;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&cfg));
    PASS();
}

static struct {
    uint64_t generated, mutations, passed, cells[6][16];
} g_breakend_mutation_cov;

static enum theft_trial_res prop_breakend_rejects_mutated_components(struct theft *t, void *arg) {
    struct kprop_breakend scene = *(const struct kprop_breakend *)arg;
    g_breakend_mutation_cov.generated++;
    for (scene.form = 0u; scene.form < 6u; scene.form++) {
        char valid[128], mutant[128];
        int n = kprop_breakend_render(&scene, valid, sizeof valid);
        if (n < 1 || (size_t)n >= sizeof valid) return THEFT_TRIAL_ERROR;
        /* A parser that rejects everything must fail before testing corruption. */
        enum theft_trial_res positive = prop_breakend_recovers_constructed_components(t, &scene);
        if (positive != THEFT_TRIAL_PASS) return positive;
        size_t length = (size_t)n;
        int after = scene.form < 2u || scene.form == 4u;
        int paired = scene.form < 4u;
        size_t sequence = after ? 0u : paired ? length - strlen(scene.replacement) : 1u;
        size_t first = after ? strlen(scene.replacement) : 0u;
        size_t colon = first + 1u + strlen(scene.chrom);
        size_t second = after ? length - 1u : sequence - 1u;
        for (unsigned fault = 0u; fault < (paired ? 16u : 4u); fault++) {
            size_t used = length;
            duckvep_breakend_status_t expected = DUCKVEP_BREAKEND_INVALID;
            memcpy(mutant, valid, length);
            switch (fault) {
                case 0: mutant[sequence] = 'R'; break;
                case 1: mutant[sequence] = '\0'; break;
                case 2: mutant[sequence] = (char)0x80; break;
                case 3: mutant[sequence] = '.'; break;
                case 4: mutant[second] = valid[first] == '[' ? ']' : '['; break;
                case 5: mutant[used++] = '['; break;
                case 6: mutant[first + 1u] = ' '; break;
                case 7: mutant[first + 1u] = '\0'; break;
                case 8: mutant[first + 1u] = ','; break;
                case 9:
                    memmove(mutant + first + 1u, mutant + colon, length - colon);
                    used -= colon - first - 1u;
                    break;
                case 10:
                    memmove(mutant + colon + 1u, mutant + second, length - second);
                    used -= second - colon - 1u;
                    break;
                case 11: mutant[colon + 1u] = 'x'; break;
                case 12: mutant[colon + 1u] = '-'; break;
                case 13: {
                    static const char overflow[] = "18446744073709551616";
                    size_t digits = sizeof overflow - 1u;
                    used = colon + 1u + digits + length - second;
                    if (used > sizeof mutant) return THEFT_TRIAL_ERROR;
                    memmove(mutant + colon + 1u + digits, mutant + second, length - second);
                    memcpy(mutant + colon + 1u, overflow, digits);
                    expected = DUCKVEP_BREAKEND_POSITION_OVERFLOW;
                    break;
                }
                case 14:
                    memmove(mutant + second, mutant + second + 1u, length - second - 1u);
                    used--;
                    break;
                case 15: mutant[colon + 1u] = '\0'; break;
            }
            uint8_t *bytes = malloc(used);
            if (bytes == NULL) return THEFT_TRIAL_ERROR;
            memcpy(bytes, mutant, used);
            struct {
                uint64_t before;
                duckvep_breakend_t parsed;
                uint64_t after;
            } guarded;
            memset(&guarded, 0xa5, sizeof guarded);
            duckvep_breakend_t zero;
            memset(&zero, 0, sizeof zero);
            const uint64_t canary = UINT64_C(0xa5a5a5a5a5a5a5a5);
            g_breakend_mutation_cov.mutations++;
            g_breakend_mutation_cov.cells[scene.form][fault]++;
            duckvep_breakend_status_t status = duckvep_breakend_parse(bytes, used, &guarded.parsed);
            int ok = status == expected && !memcmp(&guarded.parsed, &zero, sizeof zero) &&
                guarded.before == canary && guarded.after == canary && !memcmp(bytes, mutant, used);
            free(bytes);
            if (!ok) {
                fprintf(stderr, "[breakend mutation failure] form=%u fault=%u status=%u expected=%u\n",
                    scene.form, fault, (unsigned)status, (unsigned)expected);
                return THEFT_TRIAL_FAIL;
            }
            /* Failure must not poison the same destination for a valid record. */
            duckvep_breakend_t reference;
            if (duckvep_breakend_parse((const uint8_t *)valid, length, &reference) != DUCKVEP_BREAKEND_OK ||
                duckvep_breakend_parse((const uint8_t *)valid, length, &guarded.parsed) != DUCKVEP_BREAKEND_OK ||
                reference.mate_chrom != guarded.parsed.mate_chrom ||
                reference.mate_chrom_length != guarded.parsed.mate_chrom_length ||
                reference.mate_position != guarded.parsed.mate_position ||
                reference.replacement != guarded.parsed.replacement ||
                reference.replacement_length != guarded.parsed.replacement_length ||
                reference.local_join_after != guarded.parsed.local_join_after ||
                reference.mate_extends_right != guarded.parsed.mate_extends_right ||
                reference.has_mate != guarded.parsed.has_mate ||
                guarded.before != canary || guarded.after != canary) return THEFT_TRIAL_FAIL;
            g_breakend_mutation_cov.passed++;
        }
    }
    return THEFT_TRIAL_PASS;
}

TEST breakend_parser_rejects_mutated_components(void) {
    struct theft_type_info type = {.alloc = kprop_breakend_alloc, .free = theft_generic_free_cb};
    struct theft_run_config cfg = {0};
    cfg.name = "breakend_parser_rejects_mutated_components";
    cfg.prop1 = prop_breakend_rejects_mutated_components;
    cfg.type_info[0] = &type;
    cfg.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    cfg.seed = (theft_seed)kprop_env_u64("DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    memset(&g_breakend_mutation_cov, 0, sizeof g_breakend_mutation_cov);
    enum theft_run_res result = theft_run(&cfg);
    uint64_t minimum = UINT64_MAX;
    for (unsigned form = 0u; form < 6u; form++)
        for (unsigned fault = 0u; fault < (form < 4u ? 16u : 4u); fault++)
            if (g_breakend_mutation_cov.cells[form][fault] < minimum)
                minimum = g_breakend_mutation_cov.cells[form][fault];
    fprintf(stderr, "[breakend mutation coverage] generated=%" PRIu64 " cases=%" PRIu64
        " passed=%" PRIu64 " cells=72 min_per_cell=%" PRIu64 "\n",
        g_breakend_mutation_cov.generated, g_breakend_mutation_cov.mutations,
        g_breakend_mutation_cov.passed, minimum);
    ASSERT_EQ(THEFT_RUN_PASS, result);
    ASSERT_EQ(g_breakend_mutation_cov.mutations, g_breakend_mutation_cov.passed);
    ASSERT_EQ(g_breakend_mutation_cov.generated * 72u, g_breakend_mutation_cov.mutations);
    ASSERT_EQ(g_breakend_mutation_cov.generated, minimum);
    PASS();
}

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
    ASSERT_EQ(1u, event.feature_end1);
    ASSERT_EQ(1u, event.start1);
    ASSERT_EQ(1u, event.end1);
    ASSERT_EQ(0u, event.ref_diff_offset);
    ASSERT_EQ(1u, event.ref_diff_length);
    ASSERT_EQ(0u, event.alt_diff_length);

    ASSERT(duckvep_event_prepare_small(1u, insert_first_ref, 1u,
                                      insert_first_alt, 2u, &event));
    ASSERT_EQ((uint8_t)DUCKVEP_KIND_INS, event.kind);
    ASSERT_EQ(1u, event.feature_start1);
    ASSERT_EQ(0u, event.feature_end1);
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

typedef struct vep116_feature_oracle {
    uint32_t start1;
    uint32_t end1;
    uint16_t allele_offset;
    uint16_t prefix;
    uint16_t suffix;
    uint16_t ref_diff_length;
    uint16_t alt_diff_length;
} vep116_feature_oracle_t;

/* Independent transcription of Parser::post_process_vfs + minimise_alleles in
 * VEP 116. Length-changing biallelic pairs are fully prefix/suffix minimized;
 * equal-length substitutions keep the complete uploaded VariationFeature. */
static int vep116_feature_oracle(
    uint32_t                  pos1,
    const uint8_t            *ref,
    uint16_t                  ref_length,
    const uint8_t            *alt,
    uint16_t                  alt_length,
    vep116_feature_oracle_t   *out) {

    uint16_t prefix = 0u;
    uint16_t suffix = 0u;
    uint16_t ref_rem;
    uint16_t alt_rem;

    if (pos1 == 0u || ref == NULL || alt == NULL || out == NULL ||
        ref_length == 0u || alt_length == 0u) {
        return 0;
    }
    while (prefix < ref_length && prefix < alt_length &&
           ref[prefix] == alt[prefix]) {
        prefix++;
    }
    ref_rem = (uint16_t)(ref_length - prefix);
    alt_rem = (uint16_t)(alt_length - prefix);
    while (suffix < ref_rem && suffix < alt_rem &&
           ref[(uint16_t)(ref_length - 1u - suffix)] ==
           alt[(uint16_t)(alt_length - 1u - suffix)]) {
        suffix++;
    }

    out->prefix = prefix;
    out->suffix = suffix;
    out->ref_diff_length = (uint16_t)(ref_length - prefix - suffix);
    out->alt_diff_length = (uint16_t)(alt_length - prefix - suffix);
    if (out->ref_diff_length == 0u && out->alt_diff_length == 0u) return 0;

    if (ref_length != alt_length) {
        out->start1 = pos1 + (uint32_t)prefix;
        out->end1 = pos1 + (uint32_t)ref_length - 1u - (uint32_t)suffix;
        out->allele_offset = prefix;
    } else {
        out->start1 = pos1;
        out->end1 = pos1 + (uint32_t)ref_length - 1u;
        out->allele_offset = 0u;
    }
    return 1;
}

static enum theft_trial_res prop_event_normalization_matches_trim_oracle(
    struct theft *t,
    void         *arg1) {
    const struct kprop_event_norm *s = (const struct kprop_event_norm *)arg1;
    const uint8_t *ref = s->bytes + s->roff;
    const uint8_t *alt = s->bytes + s->aoff;
    const uint8_t *feature_ref;
    const uint8_t *feature_alt;
    uint16_t feature_ref_length;
    uint16_t feature_alt_length;
    vep116_feature_oracle_t feature;
    uint16_t prefix;
    uint16_t suffix;
    uint16_t ref_diff_len;
    uint16_t alt_diff_len;
    uint32_t start;
    uint32_t end;
    duckvep_event_t event;
    (void)t;

    if (!vep116_feature_oracle(s->pos, ref, s->rlen, alt, s->alen, &feature)) {
        return THEFT_TRIAL_FAIL;
    }
    prefix = feature.prefix;
    suffix = feature.suffix;
    ref_diff_len = feature.ref_diff_length;
    alt_diff_len = feature.alt_diff_length;
    if (ref_diff_len == 0u && alt_diff_len > 0u) {
        start = prefix > 0u ? s->pos + (uint32_t)prefix - 1u : s->pos;
        end = start;
    } else {
        start = s->pos + (uint32_t)prefix;
        end = ref_diff_len > 0u ? start + (uint32_t)ref_diff_len - 1u : start;
    }

    duckvep_event_load(&s->v, 0u, &event);
    if (event.raw_start1 != s->pos || event.raw_end1 != s->end) return THEFT_TRIAL_FAIL;
    if (event.feature_start1 != feature.start1 ||
        event.feature_end1 != feature.end1 ||
        event.feature_allele_offset != feature.allele_offset) {
        return THEFT_TRIAL_FAIL;
    }
    if (!duckvep_event_feature_alleles(
            &s->v, 0u, &event, &feature_ref, &feature_ref_length,
            &feature_alt, &feature_alt_length)) {
        return THEFT_TRIAL_FAIL;
    }
    if (s->rlen != s->alen) {
        if (feature_ref != ref + prefix || feature_alt != alt + prefix ||
            feature_ref_length != ref_diff_len ||
            feature_alt_length != alt_diff_len) {
            return THEFT_TRIAL_FAIL;
        }
    } else if (feature_ref != ref || feature_alt != alt ||
               feature_ref_length != s->rlen ||
               feature_alt_length != s->alen) {
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

struct theft_type_info kprop_allele_sweep_scene_info = {
    .alloc = kprop_allele_sweep_scene_alloc,
    .free = kprop_allele_sweep_scene_free,
};

static uint32_t allele_sweep_oracle_end(const duckvep_variant_batch_t *v, size_t i) {
    const uint8_t *ref = v->allele_bytes + v->ref_offset[i];
    const uint8_t *alt = v->allele_bytes + v->alt_offset[i];
    uint16_t ref_len = v->ref_length[i];
    uint16_t alt_len = v->alt_length[i];
    vep116_feature_oracle_t feature;

    if (!vep116_feature_oracle(
            v->pos1[i], ref, ref_len, alt, alt_len, &feature)) {
        return v->pos1[i];
    }

    return feature.start1 > feature.end1 ? feature.start1 : feature.end1;
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

TEST sweep_rejects_null_transcript_model(void) {
    duckvep_variant_batch_t variants;
    duckvep_sweep_cursor_t cursor;
    uint32_t active[1];
    uint32_t candidates[1];

    memset(&variants, 0, sizeof variants);
    duckvep_sweep_cursor_init(&cursor, &variants, NULL, 0u,
                              active, 1u, candidates, 1u);
    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG, cursor.status);
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
