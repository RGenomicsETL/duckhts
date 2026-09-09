#include "duckvep_property.h"

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

int popcount_u32(uint32_t x) {
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
    int in_exon = 0, frameshift_intron = 0, splice_ref = 0;
    size_t e;
    (void)t;

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
        if (e + 1u < s->ex.exon_count &&
            ee < s->es[e + 1u] && ee != UINT32_MAX) {
            uint32_t gap_start = ee + 1u;
            uint32_t gap_end = s->es[e + 1u] - 1u;

            if (gap_end - gap_start <= 12u &&
                s->pos >= gap_start && s->pos <= gap_end) {
                frameshift_intron = 1;
            }
        }
    }
    if (frameshift_intron) {
        splice_ref = 0;
        if (s->cds_s != 0u && primary == (uint32_t)DUCKVEP_REGION_INTRON)
            return THEFT_TRIAL_FAIL;
        if (s->cds_s == 0u) {
            /* The 12-base exon stretch is only a candidate lookup. VEP's
             * non_coding_exon_variant predicate rechecks the real exon. */
            if (primary != (uint32_t)DUCKVEP_REGION_INTRON)
                return THEFT_TRIAL_FAIL;
        } else if (s->pos >= s->cds_s && s->pos <= s->cds_e) {
            if (primary != (uint32_t)DUCKVEP_REGION_CDS)
                return THEFT_TRIAL_FAIL;
        } else {
            if (primary != (uint32_t)DUCKVEP_REGION_UTR)
                return THEFT_TRIAL_FAIL;
        }
    } else {
        if (popcount_u32(primary) != 1) return THEFT_TRIAL_FAIL;
        if (primary == (uint32_t)DUCKVEP_REGION_INTRON && in_exon)
            return THEFT_TRIAL_FAIL;
        if (primary != (uint32_t)DUCKVEP_REGION_INTRON && !in_exon)
            return THEFT_TRIAL_FAIL;
    }
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
           a->within_frameshift_intron == b->within_frameshift_intron &&
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
           a->intronic == b->intronic &&
           a->within_frameshift_intron == b->within_frameshift_intron &&
           a->any == b->any;
}

static enum theft_trial_res prop_sorted_point_classifier_matches_exhaustive(
    struct theft *t, void *arg1) {

    const struct kprop_point_scene *s =
        (const struct kprop_point_scene *)arg1;
    uint16_t rank = UINT16_MAX;
    uint16_t monotone_rank = UINT16_MAX;
    uint32_t first = s->tstart > 24u ? s->tstart - 24u : 1u;
    uint32_t last = s->tend + 24u;
    uint32_t pos;
    uint64_t width = (uint64_t)last - (uint64_t)first;
    uint64_t step;
    (void)t;

    /* Walk forward, then rewind the same live cursor through every coordinate.
     * Uploaded records remain position-sorted, but allele minimization can
     * move an execution coordinate backwards between adjacent rows. */
    for (step = 0u; step <= width * 2u; step++) {
        duckvep_region_state_t slow_region;
        duckvep_region_state_t fast_region;
        duckvep_splice_state_t slow_splice;
        duckvep_splice_state_t fast_splice;

        pos = step <= width
            ? first + (uint32_t)step
            : last - (uint32_t)(step - width);
        slow_region = duckvep_region_classify_span(
            &s->tx, &s->ex, 0u, pos, pos, 0u, 0u);
        slow_splice = duckvep_splice_classify_span_with_windows(
            &s->tx, &s->ex, 0u, pos, pos, 0u,
            s->splice_exonic, s->splice_intronic);
        duckvep_classify_point_sorted(
            &s->tx, &s->ex, 0u, pos,
            s->splice_exonic, s->splice_intronic,
            1u,
            &rank, &fast_region, &fast_splice);
        if (!region_states_equal(&slow_region, &fast_region) ||
            !splice_states_equal(&slow_splice, &fast_splice))
            return THEFT_TRIAL_FAIL;
        if (step <= width) {
            duckvep_classify_point_sorted(
                &s->tx, &s->ex, 0u, pos,
                s->splice_exonic, s->splice_intronic,
                0u,
                &monotone_rank, &fast_region, &fast_splice);
            if (!region_states_equal(&slow_region, &fast_region) ||
                !splice_states_equal(&slow_splice, &fast_splice))
                return THEFT_TRIAL_FAIL;
        }
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

static enum theft_trial_res prop_sorted_span_classifier_matches_exhaustive(
    struct theft *t, void *arg1) {

    const struct kprop_point_scene *s =
        (const struct kprop_point_scene *)arg1;
    uint16_t rank = UINT16_MAX;
    uint16_t monotone_rank = UINT16_MAX;
    uint32_t first = s->tstart > 24u ? s->tstart - 24u : 1u;
    uint32_t last = s->tend + 24u;
    uint32_t pos;
    uint64_t width = (uint64_t)last - (uint64_t)first;
    uint64_t step;
    uint8_t ref[16];
    uint8_t alt[16];
    (void)t;

    for (step = 0u; step <= width * 2u; step++) {
        duckvep_region_state_t slow_region;
        duckvep_region_state_t fast_region;
        duckvep_splice_state_t slow_splice;
        duckvep_splice_state_t fast_splice;
        uint32_t span_length;
        uint32_t end1;
        uint16_t ref_length;
        uint16_t alt_length;
        uint16_t i;

        pos = step <= width
            ? first + (uint32_t)step
            : last - (uint32_t)(step - width);
        span_length = 1u + pos % 17u;
        end1 = pos + span_length - 1u;
        switch (pos % 5u) {
        case 0u:
            ref_length = 0u;
            alt_length = (uint16_t)(1u + pos % 9u);
            break;
        case 1u:
            ref_length = (uint16_t)(1u + pos % 9u);
            alt_length = 0u;
            break;
        case 2u:
            ref_length = (uint16_t)(2u + pos % 8u);
            alt_length = ref_length;
            break;
        case 3u:
            ref_length = (uint16_t)(2u + pos % 8u);
            alt_length = (uint16_t)(1u + (pos / 3u) % 9u);
            break;
        default:
            ref_length = 1u;
            alt_length = 1u;
            break;
        }
        for (i = 0u; i < 16u; i++) {
            static const uint8_t bases[4] = {'A', 'C', 'G', 'T'};
            ref[i] = bases[(pos + i) & 3u];
            alt[i] = (i % 3u) == 1u
                ? bases[(pos + i + 1u) & 3u] : ref[i];
        }

        slow_region = duckvep_region_classify_span(
            &s->tx, &s->ex, 0u, pos, end1, 0u, 0u);
        fast_region = duckvep_region_classify_span_sorted(
            &s->tx, &s->ex, 0u, pos, end1, 0u, 0u, 1u, &rank);
        slow_splice = duckvep_splice_classify_differing_regions_with_windows(
            &s->tx, &s->ex, 0u, pos,
            ref, ref_length, alt, alt_length,
            s->splice_exonic, s->splice_intronic);
        fast_splice =
            duckvep_splice_classify_differing_regions_sorted_with_windows(
                &s->tx, &s->ex, 0u, pos,
                ref, ref_length, alt, alt_length,
                s->splice_exonic, s->splice_intronic, rank);
        if (!region_states_equal(&slow_region, &fast_region) ||
            !splice_states_equal(&slow_splice, &fast_splice)) {
            uint32_t current_start = 0u;
            uint32_t current_end = 0u;

            if ((uint32_t)rank < (uint32_t)s->excnt) {
                uint32_t rank_index = s->strand >= 0
                    ? (uint32_t)rank
                    : (uint32_t)s->excnt - (uint32_t)rank - 1u;
                current_start = s->es[rank_index];
                current_end = s->ee[rank_index];
            }
            fprintf(stderr,
                "\n[sorted-span mismatch] pos=%u end=%u rank=%u strand=%d "
                "tx=%u-%u cds=%u-%u exons=%u current=%u-%u "
                "ref=%u alt=%u win=%u,%u "
                "region=%u/%u intron=%u/%u exon=%u/%u cds=%u/%u utr=%u,%u/%u,%u "
                "splice=%u,%u,%u,%u,%u,%u,%u/%u,%u,%u,%u,%u,%u,%u\n",
                pos, end1, (unsigned)rank, (int)s->strand,
                s->tstart, s->tend, s->cds_s, s->cds_e,
                (unsigned)s->excnt, current_start, current_end,
                (unsigned)ref_length,
                (unsigned)alt_length, s->splice_exonic, s->splice_intronic,
                slow_region.region_mask, fast_region.region_mask,
                slow_region.overlaps_intron, fast_region.overlaps_intron,
                slow_region.overlaps_exon, fast_region.overlaps_exon,
                slow_region.overlaps_cds, fast_region.overlaps_cds,
                slow_region.overlaps_utr5, slow_region.overlaps_utr3,
                fast_region.overlaps_utr5, fast_region.overlaps_utr3,
                slow_splice.splice_donor, slow_splice.splice_acceptor,
                slow_splice.splice_donor_5th,
                slow_splice.splice_donor_region,
                slow_splice.splice_polypyrimidine,
                slow_splice.splice_region, slow_splice.intronic,
                fast_splice.splice_donor, fast_splice.splice_acceptor,
                fast_splice.splice_donor_5th,
                fast_splice.splice_donor_region,
                fast_splice.splice_polypyrimidine,
                fast_splice.splice_region, fast_splice.intronic);
            return THEFT_TRIAL_FAIL;
        }
        if (step <= width) {
            fast_region = duckvep_region_classify_span_sorted(
                &s->tx, &s->ex, 0u, pos, end1, 0u, 0u, 0u,
                &monotone_rank);
            if (!region_states_equal(&slow_region, &fast_region))
                return THEFT_TRIAL_FAIL;
        }
    }
    return THEFT_TRIAL_PASS;
}

TEST sorted_span_classifier_matches_exhaustive_for_any_transcript(void) {
    struct theft_run_config cfg;
    memset(&cfg, 0, sizeof cfg);
    cfg.name = "sorted span cursor classifier == exhaustive exon/gap scans";
    cfg.prop1 = prop_sorted_span_classifier_matches_exhaustive;
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

/* VEP 116 calls overlap() on the nominal UTR interval before checking that the
 * event has a cDNA mapping. If CDS and transcript share an endpoint, that UTR
 * interval is inverted by one base; a deletion spanning the endpoint still
 * satisfies VEP's comparison and emits both CDS and the strand-correct UTR. */
TEST region_span_crossing_empty_utr_keeps_vep_predicate(void) {
    static const uint16_t chrom[2] = {0u, 0u};
    static const uint32_t tstart[2] = {100u, 100u};
    static const uint32_t tend[2] = {200u, 200u};
    static const int8_t strand[2] = {1, -1};
    static const uint64_t flags[2] = {0u, 0u};
    static const uint32_t exoff[2] = {0u, 0u};
    static const uint16_t excnt[2] = {1u, 1u};
    static const uint32_t cds_s[2] = {100u, 100u};
    static const uint32_t cds_e[2] = {200u, 200u};
    static const uint32_t es[1] = {100u};
    static const uint32_t ee[1] = {200u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t ex;
    duckvep_region_state_t low_plus;
    duckvep_region_state_t high_plus;
    duckvep_region_state_t low_minus;
    duckvep_region_state_t high_minus;

    memset(&tx, 0, sizeof tx);
    memset(&ex, 0, sizeof ex);
    tx.chrom_id = chrom; tx.start1 = tstart; tx.end1 = tend; tx.strand = strand;
    tx.flags = flags; tx.exon_offset = exoff; tx.exon_count = excnt;
    tx.cds_start1 = cds_s; tx.cds_end1 = cds_e; tx.transcript_count = 2u;
    ex.start1 = es; ex.end1 = ee; ex.exon_count = 1u;

    low_plus = duckvep_region_classify_span(&tx, &ex, 0u, 98u, 102u, 0u, 0u);
    high_plus = duckvep_region_classify_span(&tx, &ex, 0u, 198u, 202u, 0u, 0u);
    low_minus = duckvep_region_classify_span(&tx, &ex, 1u, 98u, 102u, 0u, 0u);
    high_minus = duckvep_region_classify_span(&tx, &ex, 1u, 198u, 202u, 0u, 0u);

    ASSERT(low_plus.overlaps_cds && low_plus.overlaps_utr5 &&
           !low_plus.overlaps_utr3);
    ASSERT(high_plus.overlaps_cds && high_plus.overlaps_utr3 &&
           !high_plus.overlaps_utr5);
    ASSERT(low_minus.overlaps_cds && low_minus.overlaps_utr3 &&
           !low_minus.overlaps_utr5);
    ASSERT(high_minus.overlaps_cds && high_minus.overlaps_utr5 &&
           !high_minus.overlaps_utr3);
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

    /* The sorted-span shortcut must keep the inclusive edge of VEP's splice
     * cache and PPT windows.  With a one-base intron, the reverse-strand PPT
     * window extends into the next exon.  An insertion at start+16 has a
     * reversed mismatch interval whose preceding base is exactly on that
     * edge.  Treating the 16-base margin as half-open drops the term. */
    {
        static uint16_t tchrom[1] = {0u};
        static uint32_t tts[1]    = {1000u};
        static uint32_t tte[1]    = {1200u};
        static int8_t   tstr[1]   = {(int8_t)-1};
        static uint64_t tfl[1]    = {0u};
        static uint32_t teo[1]    = {0u};
        static uint16_t tec[1]    = {2u};
        static uint32_t tcs[1]    = {0u};
        static uint32_t tce[1]    = {0u};
        static uint32_t tes[2]    = {1105u, 1000u};
        static uint32_t tee[2]    = {1200u, 1103u};
        static uint32_t tz[2]     = {0u, 0u};
        static int8_t   tp[2]     = {0, 0};
        static const uint8_t alt[2] = {'A', 'A'};
        duckvep_transcript_model_t rtx;
        duckvep_exon_model_t rex;
        duckvep_region_state_t rr;
        duckvep_splice_state_t slow;
        duckvep_splice_state_t fast;
        uint16_t rank = UINT16_MAX;

        memset(&rtx, 0, sizeof rtx); memset(&rex, 0, sizeof rex);
        rtx.chrom_id = tchrom; rtx.start1 = tts; rtx.end1 = tte;
        rtx.strand = tstr; rtx.flags = tfl; rtx.exon_offset = teo;
        rtx.exon_count = tec; rtx.cds_start1 = tcs; rtx.cds_end1 = tce;
        rtx.transcript_count = 1u;
        rex.start1 = tes; rex.end1 = tee; rex.cdna_start1 = tz;
        rex.cdna_end1 = tz; rex.phase = tp; rex.end_phase = tp;
        rex.exon_count = 2u;

        rr = duckvep_region_classify_span_sorted(
            &rtx, &rex, 0u, 1121u, 1121u, 0u, 0u, 0u, &rank);
        slow = duckvep_splice_classify_differing_regions_with_windows(
            &rtx, &rex, 0u, 1121u, NULL, 0u, alt, 2u, 16u, 11u);
        fast = duckvep_splice_classify_differing_regions_sorted_with_windows(
            &rtx, &rex, 0u, 1121u, NULL, 0u, alt, 2u, 16u, 11u, rank);
        ASSERT(rr.overlaps_exon && rank == 1u);
        ASSERT(slow.splice_polypyrimidine && !slow.intronic);
        ASSERT(fast.splice_polypyrimidine == slow.splice_polypyrimidine);
        ASSERT(fast.any == slow.any);
    }

    /* --- VEP frameshift intron: BaseTranscriptVariation marks any intron with
     * abs(end-start) <= 12, then _intron_effects skips every ordinary intron and
     * splice predicate when a differing region overlaps it. within_cds treats
     * the gap as coding, but peptide projection remains unavailable. --- */
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
        duckvep_effect_ctx_t ctx;
        uint64_t mask;

        memset(&stx, 0, sizeof stx); memset(&sex, 0, sizeof sex);
        stx.chrom_id = schrom; stx.start1 = sts; stx.end1 = ste; stx.strand = sstr;
        stx.flags = sfl; stx.exon_offset = seo; stx.exon_count = sec;
        stx.cds_start1 = scs; stx.cds_end1 = sce; stx.transcript_count = 1u;
        sex.start1 = ses; sex.end1 = see; sex.cdna_start1 = sz2; sex.cdna_end1 = sz2;
        sex.phase = szp2; sex.end_phase = szp2; sex.exon_count = 2u;

        r = duckvep_splice_classify(&stx, &sex, 0, 201u);
        ASSERT(r.within_frameshift_intron && !r.any && !r.intronic);
        r = duckvep_splice_classify(&stx, &sex, 0, 206u);
        ASSERT(r.within_frameshift_intron && !r.any && !r.intronic);
        r = duckvep_splice_classify(&stx, &sex, 0, 204u);
        ASSERT(r.within_frameshift_intron && !r.any && !r.intronic);
        {
            static const uint8_t ref[2] = {'T', 'T'};
            static const uint8_t alt[2] = {'A', 'A'};

            r = duckvep_splice_classify_differing_regions(
                &stx, &sex, 0u, 200u, ref, 2u, alt, 2u);
            ASSERT(r.within_frameshift_intron && !r.any && !r.intronic);
        }

        duckvep_effect_ctx_fill(
            &stx, &sex, 0u, 0u, 204u, 204u, 0u,
            DUCKVEP_DEFAULT_SPLICE_REGION_EXONIC,
            DUCKVEP_DEFAULT_SPLICE_REGION_INTRONIC, &ctx);
        ASSERT(ctx.region_state.within_frameshift_intron);
        ASSERT(ctx.region_state.overlaps_cds);
        ASSERT(!ctx.region_state.overlaps_intron);
        ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_DONOR)) == 0u);
        ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_ACCEPTOR)) == 0u);
        duckvep_effect_ctx_finalize(&ctx);
        mask = duckvep_effect_eval(ctx.pre_bits);
        ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE), mask);
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
    duckvep_effect_ctx_fill(
        &tx, &ex, 0u, 0u, 285u, 290u, 0u,
        DUCKVEP_DEFAULT_SPLICE_REGION_EXONIC,
        DUCKVEP_DEFAULT_SPLICE_REGION_INTRONIC, &ctx);
    ASSERT(!ctx.region_state.overlaps_exon);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_PPT)) != 0u);

    /* deletion from the tract across the acceptor into the exon: exon overlap
     * suppresses PPT, but the essential acceptor site still fires */
    duckvep_effect_ctx_fill(
        &tx, &ex, 0u, 0u, 285u, 305u, 0u,
        DUCKVEP_DEFAULT_SPLICE_REGION_EXONIC,
        DUCKVEP_DEFAULT_SPLICE_REGION_INTRONIC, &ctx);
    ASSERT(ctx.region_state.overlaps_exon);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_PPT)) == 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_ACCEPTOR)) != 0u);
    PASS();
}

/* VEP sets one transcript-wide `_has_frameshift_intron` cache bit when any
 * intron is at most 13 bases long. Its coarse exon predicate then stretches
 * every exon by 12 bases, including exons beside unrelated ordinary introns.
 * The exact region remains intronic, but the coarse exon=0 gate suppresses PPT.
 * This is the GRCh37 ENST00000262952/ENST00000398240/ENST00000543616 state. */
TEST annotate_remote_frameshift_intron_suppresses_ppt(void) {
    static const uint16_t tchrom[1] = {0u};
    static const uint32_t tstart[1] = {100u};
    static const uint32_t tend[1] = {500u};
    static const int8_t strand[1] = {(int8_t)1};
    static const uint64_t flags[1] = {0u};
    static const uint32_t exon_offset[1] = {0u};
    static const uint16_t exon_count[1] = {3u};
    static const uint32_t no_cds[1] = {0u};
    static const uint32_t exon_start[3] = {100u, 301u, 357u};
    static const uint32_t exon_end[3] = {200u, 350u, 500u};
    static const uint32_t cdna_start[3] = {1u, 102u, 152u};
    static const uint32_t cdna_end[3] = {101u, 151u, 295u};
    static const int8_t phase[3] = {-1, -1, -1};
    static const uint16_t vchrom[2] = {0u, 0u};
    static const uint32_t vpos[2] = {289u, 353u};
    static const uint8_t vkind[2] = {
        (uint8_t)DUCKVEP_KIND_SNV,
        (uint8_t)DUCKVEP_KIND_SNV
    };
    static const uint8_t alleles[4] = {'C', 'G', 'A', 'T'};
    static const uint32_t ref_offset[2] = {0u, 2u};
    static const uint32_t alt_offset[2] = {1u, 3u};
    static const uint16_t allele_length[2] = {1u, 1u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t variants;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_consequence_t rows[2];
    duckvep_result_builder_t builder;
    duckvep_error_t error;
    uint64_t expected;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&variants, 0, sizeof variants);
    memset(&error, 0, sizeof error);
    tx.chrom_id = tchrom; tx.start1 = tstart; tx.end1 = tend;
    tx.strand = strand; tx.flags = flags; tx.exon_offset = exon_offset;
    tx.exon_count = exon_count; tx.cds_start1 = no_cds;
    tx.cds_end1 = no_cds; tx.transcript_count = 1u;
    exons.start1 = exon_start; exons.end1 = exon_end;
    exons.cdna_start1 = cdna_start; exons.cdna_end1 = cdna_end;
    exons.phase = phase; exons.end_phase = phase; exons.exon_count = 3u;
    variants.chrom_id = vchrom; variants.pos1 = vpos; variants.end1 = vpos;
    variants.ref_offset = ref_offset; variants.ref_length = allele_length;
    variants.alt_offset = alt_offset; variants.alt_length = allele_length;
    variants.allele_bytes = alleles;
    variants.allele_bytes_len = sizeof alleles;
    variants.variant_kind = vkind; variants.count = 2u;

    ASSERT_EQ(DUCKVEP_OK,
              duckvep_model_open(&tx, &exons, NULL, NULL, &model, &error));
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_options_open(NULL, &options, &error));
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_workspace_open(model, &workspace, &error));
    duckvep_result_builder_init(&builder, rows, 2u);
    ASSERT_EQ(DUCKVEP_OK,
              duckvep_annotate_tile(model, &variants, options, workspace,
                                    &builder, &error));
    ASSERT_EQ(2u, duckvep_result_builder_count(&builder));
    expected = DUCKVEP_SO(DUCKVEP_SO_INTRON) |
               DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT);
    ASSERT_EQ(expected, rows[0].consequence_mask);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_INTRON, rows[0].region_mask);

    /* The second point is inside the six-base frameshift intron. The same
     * transcript-wide stretch finds an exon candidate, but VEP's
     * non_coding_exon_variant double-check uses the real exon coordinates;
     * _intron_effects also suppresses intron_variant in this gap. */
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT),
              rows[1].consequence_mask);
    ASSERT_EQ((uint32_t)DUCKVEP_REGION_INTRON, rows[1].region_mask);

    duckvep_workspace_close(workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
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
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
              DUCKVEP_SO(DUCKVEP_SO_INCOMPLETE_TERMINAL_CODON),
              duckvep_effect_eval(DUCKVEP_PRE(DUCKVEP_PRE_CDS) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_DELTA) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_CODING_UNKNOWN) |
                                  DUCKVEP_PRE(DUCKVEP_PRE_PARTIAL_CODON)));
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

TEST consequence_predicate_flags_survive_so_mask_omission(void) {
    duckvep_sequence_delta_t source;
    duckvep_sequence_delta_t restored;
    uint32_t flags;
    uint64_t mask;

    memset(&source, 0, sizeof source);
    memset(&restored, 0, sizeof restored);
    source.valid = 1u;
    source.sequence_status = (uint8_t)DUCKVEP_SEQUENCE_RESOLVED;
    source.frameshift = 1u;
    flags = duckvep_sequence_delta_consequence_flags(&source, 1);
    ASSERT((flags &
            (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_PREDICATES_VALID) != 0u);
    ASSERT((flags &
            (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_FRAMESHIFT) != 0u);

    /* Executable VEP 116 witnesses can emit this term set while HGVSp still
     * follows the raw frameshift predicate. The compact SO mask therefore is
     * deliberately not the storage authority for peptide predicates. */
    mask = DUCKVEP_SO(DUCKVEP_SO_CODING_SEQUENCE) |
           DUCKVEP_SO(DUCKVEP_SO_SPLICE_ACCEPTOR);
    ASSERT((mask & DUCKVEP_SO(DUCKVEP_SO_SPLICE_ACCEPTOR)) != 0u);
    ASSERT((mask & DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT)) == 0u);
    ASSERT_EQ(1, duckvep_sequence_delta_apply_consequence_flags(
                     flags, &restored));
    ASSERT(restored.frameshift != 0u);

    /* Frameshift and stop_retained bits are positive evidence. start_lost and
     * stop_lost are cached before HGVS placement, so a valid sidecar must also
     * preserve their false state. */
    memset(&source, 0, sizeof source);
    source.valid = 1u;
    flags = duckvep_sequence_delta_consequence_flags(&source, 1);
    restored.frameshift = 1u;
    restored.start_lost = 1u;
    restored.stop_lost = 1u;
    ASSERT_EQ(1, duckvep_sequence_delta_apply_consequence_flags(
                     flags, &restored));
    ASSERT(restored.frameshift != 0u);
    ASSERT(!restored.start_lost);
    ASSERT(!restored.stop_lost);

    source.start_lost = 1u;
    flags = duckvep_sequence_delta_consequence_flags(&source, 1);
    ASSERT_EQ(1, duckvep_sequence_delta_apply_consequence_flags(
                     flags, &restored));
    ASSERT(restored.start_lost != 0u);

    source.stop_lost = 1u;
    flags = duckvep_sequence_delta_consequence_flags(&source, 1);
    ASSERT_EQ(1, duckvep_sequence_delta_apply_consequence_flags(
                     flags, &restored));
    ASSERT(restored.stop_lost != 0u);

    memset(&source, 0, sizeof source);
    flags = duckvep_sequence_delta_consequence_flags(&source, 1);
    ASSERT((flags &
            (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_UNRESOLVED) != 0u);
    ASSERT_EQ(0, duckvep_sequence_delta_apply_consequence_flags(
                     flags, &restored));
    PASS();
}

/* Regulatory and motif features are independent resident interval arrays. This
 * pins the VEP tier-2 predicates evaluated for one event/feature pair,
 * including the reversed insertion interval at a feature edge. */
TEST interval_feature_consequences_known_scene(void) {
    duckvep_event_t event;
    static const uint8_t insertion_ref[1] = {'A'};
    static const uint8_t insertion_alt[2] = {'A', 'T'};
    uint64_t expected;

    memset(&event, 0, sizeof event);
    event.kind = (uint8_t)DUCKVEP_KIND_SNV;
    event.feature_start1 = 110u;
    event.feature_end1 = 110u;
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_TF_BINDING_SITE),
              duckvep_effect_eval_interval_feature(
                  DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE,
                  &event, 0u, 100u, 120u));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_REGULATORY_REGION),
              duckvep_effect_eval_interval_feature(
                  DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION,
                  &event, 0u, 100u, 120u));

    event.kind = (uint8_t)DUCKVEP_KIND_SV;
    event.feature_start1 = 90u;
    event.feature_end1 = 130u;
    event.sv_type = (uint8_t)DUCKVEP_SV_DELETION;
    event.copy_change = (uint8_t)DUCKVEP_COPY_CHANGE_LOSS;
    expected = DUCKVEP_SO(DUCKVEP_SO_TFBS_ABLATION) |
               DUCKVEP_SO(DUCKVEP_SO_TF_BINDING_SITE);
    ASSERT_EQ(expected,
              duckvep_effect_eval_interval_feature(
                  DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE,
                  &event, 0u, 100u, 120u));
    expected = DUCKVEP_SO(DUCKVEP_SO_REGULATORY_REGION_ABLATION) |
               DUCKVEP_SO(DUCKVEP_SO_REGULATORY_REGION);
    ASSERT_EQ(expected,
              duckvep_effect_eval_interval_feature(
                  DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION,
                  &event, 0u, 100u, 120u));

    event.sv_type = (uint8_t)DUCKVEP_SV_DUPLICATION;
    event.copy_change = (uint8_t)DUCKVEP_COPY_CHANGE_GAIN;
    expected = DUCKVEP_SO(DUCKVEP_SO_TFBS_AMPLIFICATION) |
               DUCKVEP_SO(DUCKVEP_SO_TF_BINDING_SITE);
    ASSERT_EQ(expected,
              duckvep_effect_eval_interval_feature(
                  DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE,
                  &event, 0u, 100u, 120u));
    expected = DUCKVEP_SO(DUCKVEP_SO_REGULATORY_REGION_AMPLIFICATION) |
               DUCKVEP_SO(DUCKVEP_SO_REGULATORY_REGION);
    ASSERT_EQ(expected,
              duckvep_effect_eval_interval_feature(
                  DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION,
                  &event, 0u, 100u, 120u));

    /* VEP retains an oversized/unexpanded CNV:TR as a distinct structural
     * tandem-repeat class, but its regulatory consequence predicates are the
     * same gain predicates as a tandem duplication. */
    event.sv_type = (uint8_t)DUCKVEP_SV_TANDEM_REPEAT;
    event.copy_change = (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN;
    ASSERT_EQ(expected,
              duckvep_effect_eval_interval_feature(
                  DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION,
                  &event, 0u, 100u, 120u));

    /* A partial deletion overlaps but does not contain the feature, so VEP's
     * feature_ablation predicate remains false. */
    event.feature_start1 = 110u;
    event.feature_end1 = 130u;
    event.sv_type = (uint8_t)DUCKVEP_SV_DELETION;
    event.copy_change = (uint8_t)DUCKVEP_COPY_CHANGE_LOSS;
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_TF_BINDING_SITE),
              duckvep_effect_eval_interval_feature(
                  DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE,
                  &event, 0u, 100u, 120u));

    ASSERT(duckvep_event_prepare_small(
        120u, insertion_ref, 1u, insertion_alt, 2u, &event));
    ASSERT_EQ(0u, duckvep_effect_eval_interval_feature(
        DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION,
        &event, 0u, 100u, 120u));
    ASSERT(duckvep_event_prepare_small(
        110u, insertion_ref, 1u, insertion_alt, 2u, &event));
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_REGULATORY_REGION),
              duckvep_effect_eval_interval_feature(
                  DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION,
                  &event, 0u, 100u, 120u));

    event.feature_start1 = 200u;
    event.feature_end1 = 200u;
    ASSERT_EQ(0u, duckvep_effect_eval_interval_feature(
        DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE,
        &event, 0u, 100u, 120u));
    ASSERT_EQ(0u, duckvep_effect_eval_interval_feature(
        (duckvep_interval_feature_kind_t)0,
        &event, 0u, 100u, 120u));
    PASS();
}

TEST model_open_rejects_invalid_interval_feature_models(void) {
    duckvep_transcript_model_t transcripts;
    duckvep_exon_model_t exons;
    duckvep_interval_feature_model_t features;
    duckvep_model_t *model = NULL;
    duckvep_error_t error;
    uint16_t chrom[2] = {0u, 0u};
    uint32_t start[2] = {100u, 200u};
    uint32_t end[2] = {120u, 220u};
    uint8_t kind[2] = {
        (uint8_t)DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION,
        (uint8_t)DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE
    };

    memset(&transcripts, 0, sizeof transcripts);
    memset(&exons, 0, sizeof exons);
    memset(&features, 0, sizeof features);
    memset(&error, 0, sizeof error);
    features.chrom_id = chrom;
    features.start1 = start;
    features.end1 = end;
    features.kind = kind;
    features.feature_count = 2u;

    ASSERT_EQ(DUCKVEP_OK, duckvep_model_open(
        &transcripts, &exons, NULL, &features, &model, &error));
    duckvep_model_close(model);
    model = NULL;

    features.kind = NULL;
    ASSERT_EQ(DUCKVEP_ERR_INVALID_ARG,
              duckvep_model_open(
                  &transcripts, &exons, NULL, &features, &model, &error));
    features.kind = kind;

    start[0] = 0u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(
                  &transcripts, &exons, NULL, &features, &model, &error));
    start[0] = 100u;

    kind[0] = 0u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(
                  &transcripts, &exons, NULL, &features, &model, &error));
    kind[0] = (uint8_t)DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION;

    start[0] = 201u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(
                  &transcripts, &exons, NULL, &features, &model, &error));
    PASS();
}

struct interval_feature_observation {
    uint32_t variant_idx;
    uint32_t feature_idx;
    uint64_t consequence_mask;
};

static int interval_feature_observation_cmp(const void *left,
                                            const void *right) {
    const struct interval_feature_observation *a =
        (const struct interval_feature_observation *)left;
    const struct interval_feature_observation *b =
        (const struct interval_feature_observation *)right;

    if (a->variant_idx != b->variant_idx)
        return a->variant_idx < b->variant_idx ? -1 : 1;
    if (a->feature_idx != b->feature_idx)
        return a->feature_idx < b->feature_idx ? -1 : 1;
    if (a->consequence_mask != b->consequence_mask)
        return a->consequence_mask < b->consequence_mask ? -1 : 1;
    return 0;
}

static int interval_feature_kernel_observations(
    const duckvep_variant_batch_t *variants,
    const duckvep_interval_feature_model_t *features,
    struct interval_feature_observation *observed,
    size_t capacity,
    size_t *observed_count) {
    duckvep_transcript_model_t transcripts;
    duckvep_exon_model_t exons;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_annotate_cursor_t *cursor = NULL;
    duckvep_error_t error;
    uint32_t seed[KPROP_MAX_TX];
    size_t seed_count = 0u;
    size_t count = 0u;
    size_t feature;
    int ok = 0;

    memset(&transcripts, 0, sizeof transcripts);
    memset(&exons, 0, sizeof exons);
    memset(&error, 0, sizeof error);
    if (duckvep_model_open(
            &transcripts, &exons, NULL, features, &model, &error) !=
        DUCKVEP_OK) goto done;
    if (duckvep_options_open(NULL, &options, &error) != DUCKVEP_OK)
        goto done;
    if (duckvep_workspace_open(model, &workspace, &error) != DUCKVEP_OK)
        goto done;
    if (duckvep_annotate_cursor_open(model, variants, options, workspace,
                                     &cursor, &error) != DUCKVEP_OK)
        goto done;

    /* This is the exact seed contract used by the cgranges adapter: only
     * intervals overlapping the first raw point enter the persistent active
     * set; a longer differing region is discovered from the forward frontier. */
    for (feature = 0u; feature < features->feature_count; feature++) {
        if (features->chrom_id[feature] == variants->chrom_id[0] &&
            features->start1[feature] <= variants->pos1[0] &&
            features->end1[feature] >= variants->pos1[0]) {
            if (seed_count >= KPROP_MAX_TX) goto done;
            seed[seed_count++] = (uint32_t)feature;
        }
    }
    if (duckvep_annotate_cursor_seed_interval_features(
            cursor, seed, seed_count, &error) != DUCKVEP_OK) goto done;

    while (!duckvep_annotate_cursor_done(cursor)) {
        duckvep_consequence_t row;
        duckvep_result_builder_t builder;
        duckvep_status_t status;

        duckvep_result_builder_init(&builder, &row, 1u);
        status = duckvep_annotate_cursor_fill(cursor, &builder, &error);
        if (status != DUCKVEP_OK && status != DUCKVEP_ERR_RESULT_FULL)
            goto done;
        if (builder.count > 1u) goto done;
        if (builder.count == 1u) {
            if (count >= capacity ||
                row.overlap_object_kind ==
                    (uint8_t)DUCKVEP_OVERLAP_OBJECT_TRANSCRIPT ||
                row.interval_feature_idx >= features->feature_count)
                goto done;
            observed[count].variant_idx = row.variant_idx;
            observed[count].feature_idx = row.interval_feature_idx;
            observed[count].consequence_mask = row.consequence_mask;
            count++;
        }
    }
    *observed_count = count;
    ok = 1;
done:
    duckvep_annotate_cursor_close(cursor);
    duckvep_workspace_close(workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
    return ok;
}

static int interval_feature_bruteforce_observations(
    const duckvep_variant_batch_t *variants,
    const duckvep_interval_feature_model_t *features,
    struct interval_feature_observation *observed,
    size_t capacity,
    size_t *observed_count) {
    size_t count = 0u;
    size_t variant;

    for (variant = 0u; variant < variants->count; variant++) {
        duckvep_event_t event;
        size_t feature;

        duckvep_event_load(variants, variant, &event);
        for (feature = 0u; feature < features->feature_count; feature++) {
            uint64_t mask;

            if (variants->chrom_id[variant] != features->chrom_id[feature])
                continue;
            mask = duckvep_effect_eval_interval_feature(
                (duckvep_interval_feature_kind_t)features->kind[feature],
                &event, features->chrom_id[feature],
                features->start1[feature], features->end1[feature]);
            if (mask == 0u) continue;
            if (count >= capacity) return 0;
            observed[count].variant_idx = (uint32_t)variant;
            observed[count].feature_idx = (uint32_t)feature;
            observed[count].consequence_mask = mask;
            count++;
        }
    }
    *observed_count = count;
    return 1;
}

static int interval_feature_run_matches_bruteforce(
    const duckvep_variant_batch_t *variants,
    const duckvep_interval_feature_model_t *features) {
    struct interval_feature_observation kernel[KPROP_MAX_PAIRS];
    struct interval_feature_observation brute[KPROP_MAX_PAIRS];
    size_t kernel_count = 0u;
    size_t brute_count = 0u;
    size_t index;

    if (!interval_feature_kernel_observations(
            variants, features, kernel, KPROP_MAX_PAIRS, &kernel_count) ||
        !interval_feature_bruteforce_observations(
            variants, features, brute, KPROP_MAX_PAIRS, &brute_count))
        return 0;
    if (kernel_count != brute_count) return 0;
    qsort(kernel, kernel_count, sizeof kernel[0],
          interval_feature_observation_cmp);
    qsort(brute, brute_count, sizeof brute[0],
          interval_feature_observation_cmp);
    for (index = 0u; index < kernel_count; index++) {
        if (interval_feature_observation_cmp(&kernel[index],
                                             &brute[index]) != 0)
            return 0;
    }
    return 1;
}

static int interval_feature_breakend_pairs_match_bruteforce(
    const duckvep_variant_batch_t *source,
    const duckvep_interval_feature_model_t *features) {
    duckvep_transcript_model_t transcripts;
    duckvep_exon_model_t exons;
    duckvep_variant_batch_t variants;
    duckvep_interval_feature_pairs_t pairs;
    duckvep_model_t *model = NULL;
    duckvep_options_t *options = NULL;
    duckvep_workspace_t *workspace = NULL;
    duckvep_consequence_t rows[KPROP_MAX_PAIRS];
    struct interval_feature_observation expected[KPROP_MAX_PAIRS];
    duckvep_result_builder_t results;
    duckvep_error_t error;
    uint32_t ends[KPROP_MAX_VARIANTS];
    uint16_t mate_chrom[KPROP_MAX_VARIANTS];
    uint32_t mate_pos[KPROP_MAX_VARIANTS];
    uint8_t kinds[KPROP_MAX_VARIANTS];
    uint8_t sv_types[KPROP_MAX_VARIANTS];
    uint8_t copy_changes[KPROP_MAX_VARIANTS];
    uint32_t pair_variant[KPROP_MAX_PAIRS];
    uint32_t pair_feature[KPROP_MAX_PAIRS];
    size_t pair_count = 0u;
    size_t expected_count = 0u;
    size_t variant;
    int ok = 0;

    if (source->count > KPROP_MAX_VARIANTS ||
        features->feature_count > KPROP_MAX_TX)
        return 0;
    memset(&transcripts, 0, sizeof transcripts);
    memset(&exons, 0, sizeof exons);
    memset(&variants, 0, sizeof variants);
    memset(&error, 0, sizeof error);
    variants = *source;
    variants.end1 = ends;
    variants.mate_chrom_id = mate_chrom;
    variants.mate_pos1 = mate_pos;
    variants.ref_offset = NULL;
    variants.ref_length = NULL;
    variants.alt_offset = NULL;
    variants.alt_length = NULL;
    variants.allele_bytes = NULL;
    variants.allele_bytes_len = 0u;
    variants.variant_kind = kinds;
    variants.sv_type = sv_types;
    variants.copy_change = copy_changes;

    for (variant = 0u; variant < variants.count; variant++) {
        size_t mate_feature = features->feature_count == 0u
            ? 0u : (variant * 5u + 1u) % features->feature_count;
        size_t feature;

        if (variants.pos1[variant] == UINT32_MAX) goto done;
        ends[variant] = variants.pos1[variant];
        kinds[variant] = (uint8_t)DUCKVEP_KIND_SV;
        sv_types[variant] = (uint8_t)DUCKVEP_SV_BREAKEND;
        copy_changes[variant] = (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN;
        if (features->feature_count != 0u) {
            mate_chrom[variant] = features->chrom_id[mate_feature];
            mate_pos[variant] = variant % 3u == 0u
                ? features->start1[mate_feature]
                : (variant % 3u == 1u
                    ? features->end1[mate_feature]
                    : features->end1[mate_feature] == UINT32_MAX
                        ? features->end1[mate_feature]
                        : features->end1[mate_feature] + 1u);
        } else {
            mate_chrom[variant] = variants.chrom_id[variant];
            mate_pos[variant] = variants.pos1[variant];
        }
        for (feature = 0u; feature < features->feature_count; feature++) {
            uint32_t local_point = variants.pos1[variant] + 1u;
            int local_hit = features->chrom_id[feature] ==
                    variants.chrom_id[variant] &&
                local_point >= features->start1[feature] &&
                local_point <= features->end1[feature];
            int mate_hit = features->chrom_id[feature] == mate_chrom[variant] &&
                mate_pos[variant] >= features->start1[feature] &&
                mate_pos[variant] <= features->end1[feature];
            int local_close = features->chrom_id[feature] ==
                    variants.chrom_id[variant] &&
                (local_point < features->start1[feature]
                    ? features->start1[feature] - local_point <=
                        DUCKVEP_BREAKEND_ALLELE_DISTANCE
                    : local_point > features->end1[feature]
                        ? local_point - features->end1[feature] <=
                            DUCKVEP_BREAKEND_ALLELE_DISTANCE
                        : 1);

            if (pair_count >= KPROP_MAX_PAIRS) goto done;
            pair_variant[pair_count] = (uint32_t)variant;
            pair_feature[pair_count] = (uint32_t)feature;
            pair_count++;
            if (local_hit || mate_hit) {
                if (expected_count >= KPROP_MAX_PAIRS) goto done;
                expected[expected_count].variant_idx = (uint32_t)variant;
                expected[expected_count].feature_idx = (uint32_t)feature;
                expected[expected_count].consequence_mask = local_hit
                    ? (features->kind[feature] ==
                            (uint8_t)DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION
                        ? DUCKVEP_SO(DUCKVEP_SO_REGULATORY_REGION)
                        : DUCKVEP_SO(DUCKVEP_SO_TF_BINDING_SITE))
                    : DUCKVEP_SO(DUCKVEP_SO_FEATURE_TRUNCATION) |
                        (local_close
                            ? DUCKVEP_SO(DUCKVEP_SO_INTERGENIC)
                            : UINT64_C(0));
                expected_count++;
            }
        }
    }

    pairs.variant_idx = pair_variant;
    pairs.feature_idx = pair_feature;
    pairs.count = pair_count;
    if (duckvep_model_open(
            &transcripts, &exons, NULL, features, &model, &error) !=
        DUCKVEP_OK) goto done;
    if (duckvep_options_open(NULL, &options, &error) != DUCKVEP_OK)
        goto done;
    if (duckvep_workspace_open(model, &workspace, &error) != DUCKVEP_OK)
        goto done;
    duckvep_result_builder_init(&results, rows, KPROP_MAX_PAIRS);
    if (duckvep_annotate_interval_feature_pairs(
            model, &variants, &pairs, options, workspace, &results, &error) !=
        DUCKVEP_OK || results.count != expected_count)
        goto done;
    for (variant = 0u; variant < expected_count; variant++) {
        if (rows[variant].variant_idx != expected[variant].variant_idx ||
            rows[variant].interval_feature_idx !=
                expected[variant].feature_idx ||
            rows[variant].consequence_mask !=
                expected[variant].consequence_mask)
            goto done;
    }
    ok = 1;
done:
    duckvep_workspace_close(workspace);
    duckvep_options_close(options);
    duckvep_model_close(model);
    return ok;
}

static enum theft_trial_res prop_interval_feature_cursor_matches_bruteforce(
    struct theft *theft, void *argument) {
    const struct kprop_allele_sweep_scene *scene =
        (const struct kprop_allele_sweep_scene *)argument;
    duckvep_interval_feature_model_t features;
    duckvep_variant_batch_t small;
    duckvep_variant_batch_t structural;
    uint8_t kinds[KPROP_MAX_TX];
    uint8_t small_kinds[KPROP_MAX_VARIANTS];
    uint8_t structural_kinds[KPROP_MAX_VARIANTS];
    uint8_t structural_types[KPROP_MAX_VARIANTS];
    uint8_t copy_changes[KPROP_MAX_VARIANTS];
    size_t index;
    (void)theft;

    memset(&features, 0, sizeof features);
    for (index = 0u; index < scene->tx.transcript_count; index++) {
        kinds[index] = (uint8_t)(index % 2u == 0u
            ? DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION
            : DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE);
    }
    features.chrom_id = scene->tchrom;
    features.start1 = scene->tstart;
    features.end1 = scene->tend;
    features.kind = kinds;
    features.feature_count = scene->tx.transcript_count;

    /* The shared sweep generator deliberately varies the uploaded allele
     * representation. Derive the semantic event kind from its differing
     * region, exactly as the SQL adapter does once per ALT, so this property
     * exercises the kernel contract rather than feeding contradictory kinds. */
    small = scene->v;
    for (index = 0u; index < small.count; index++) {
        duckvep_event_t event;
        const uint8_t *ref = small.allele_bytes + small.ref_offset[index];
        const uint8_t *alt = small.allele_bytes + small.alt_offset[index];

        if (!duckvep_event_prepare_small(
                small.pos1[index], ref, small.ref_length[index],
                alt, small.alt_length[index], &event))
            return THEFT_TRIAL_ERROR;
        small_kinds[index] = event.kind;
    }
    small.variant_kind = small_kinds;
    if (!interval_feature_run_matches_bruteforce(&small, &features))
        return THEFT_TRIAL_FAIL;

    memset(&structural, 0, sizeof structural);
    structural = small;
    structural.ref_offset = NULL;
    structural.ref_length = NULL;
    structural.alt_offset = NULL;
    structural.alt_length = NULL;
    structural.allele_bytes = NULL;
    structural.allele_bytes_len = 0u;
    structural.variant_kind = structural_kinds;
    structural.sv_type = structural_types;
    structural.copy_change = copy_changes;
    for (index = 0u; index < structural.count; index++) {
        structural_kinds[index] = (uint8_t)DUCKVEP_KIND_SV;
        if (index % 4u == 0u) {
            structural_types[index] = (uint8_t)DUCKVEP_SV_DELETION;
            copy_changes[index] = (uint8_t)DUCKVEP_COPY_CHANGE_LOSS;
        } else if (index % 4u == 1u) {
            structural_types[index] = (uint8_t)DUCKVEP_SV_DUPLICATION;
            copy_changes[index] = (uint8_t)DUCKVEP_COPY_CHANGE_GAIN;
        } else if (index % 4u == 2u) {
            structural_types[index] = (uint8_t)DUCKVEP_SV_TANDEM_REPEAT;
            copy_changes[index] = (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN;
        } else {
            structural_types[index] = (uint8_t)DUCKVEP_SV_INVERSION;
            copy_changes[index] = (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN;
        }
    }
    if (!interval_feature_run_matches_bruteforce(&structural, &features))
        return THEFT_TRIAL_FAIL;
    if (!interval_feature_breakend_pairs_match_bruteforce(&small, &features))
        return THEFT_TRIAL_FAIL;
    return THEFT_TRIAL_PASS;
}

TEST interval_feature_cursor_matches_bruteforce_for_random_events(void) {
    struct theft_run_config config;

    memset(&config, 0, sizeof config);
    config.name = "regulation sweep/BND pairs == independent feature oracles";
    config.prop1 = prop_interval_feature_cursor_matches_bruteforce;
    config.type_info[0] = &kprop_allele_sweep_scene_info;
    config.trials = kprop_env_u64("DUCKVEP_PROP_TRIALS", KPROP_DEFAULT_TRIALS);
    config.seed = (theft_seed)kprop_env_u64(
        "DUCKVEP_PROP_SEED", KPROP_DEFAULT_SEED);
    ASSERT_EQ(THEFT_RUN_PASS, theft_run(&config));
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

    duckvep_effect_ctx_fill(
        &tx, &ex, 0u, 0u, 200u, 200u, 0u,
        DUCKVEP_DEFAULT_SPLICE_REGION_EXONIC,
        DUCKVEP_DEFAULT_SPLICE_REGION_INTRONIC, &ctx);
    duckvep_effect_ctx_finalize(&ctx);
    mask = duckvep_effect_eval(ctx.pre_bits);
    ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_INTRON) |
              DUCKVEP_SO(DUCKVEP_SO_NMD_TRANSCRIPT), mask);

    duckvep_effect_ctx_fill(
        &tx, &ex, 0u, 0u, 90u, 90u, 0u,
        DUCKVEP_DEFAULT_SPLICE_REGION_EXONIC,
        DUCKVEP_DEFAULT_SPLICE_REGION_INTRONIC, &ctx);
    ASSERT((ctx.pre_bits &
            DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_NMD_TRANSCRIPT)) == 0u);
    PASS();
}

/* Ensembl's repeated miRNA attributes are projected once into small genomic
 * exon segments. The tier-2 mature term must replace, not accompany, the
 * generic noncoding-exon consequence on either transcript strand. */
TEST mature_mirna_predicate_known_scene(void) {
    uint16_t chrom[1] = {0u};
    uint32_t transcript_start[1] = {100u};
    uint32_t transcript_end[1] = {120u};
    int8_t strand[1] = {1};
    uint64_t flags[1] = {(uint64_t)DUCKVEP_TX_BIOTYPE_MIRNA};
    uint32_t exon_offset[1] = {0u};
    uint16_t exon_count[1] = {1u};
    uint32_t cds_start[1] = {0u};
    uint32_t cds_end[1] = {0u};
    uint32_t exon_start[1] = {100u};
    uint32_t exon_end[1] = {120u};
    uint32_t cdna_start[1] = {1u};
    uint32_t cdna_end[1] = {21u};
    int8_t phase[1] = {-1};
    uint32_t mature_offset[2] = {0u, 1u};
    uint32_t mature_start[1] = {105u};
    uint32_t mature_end[1] = {110u};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_effect_ctx_t ctx;
    duckvep_model_t *model = NULL;
    duckvep_error_t error;
    duckvep_event_t event;
    static const uint8_t insertion_ref[1] = {'A'};
    static const uint8_t insertion_alt[2] = {'A', 'C'};
    size_t orientation;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&error, 0, sizeof error);
    tx.chrom_id = chrom; tx.start1 = transcript_start;
    tx.end1 = transcript_end; tx.strand = strand; tx.flags = flags;
    tx.exon_offset = exon_offset; tx.exon_count = exon_count;
    tx.cds_start1 = cds_start; tx.cds_end1 = cds_end;
    tx.transcript_count = 1u;
    tx.mature_mirna_offset = mature_offset;
    tx.mature_mirna_start1 = mature_start;
    tx.mature_mirna_end1 = mature_end;
    tx.mature_mirna_count = 1u;
    exons.start1 = exon_start; exons.end1 = exon_end;
    exons.cdna_start1 = cdna_start; exons.cdna_end1 = cdna_end;
    exons.phase = phase; exons.end_phase = phase; exons.exon_count = 1u;

    ASSERT_EQ(DUCKVEP_OK,
              duckvep_model_open(&tx, &exons, NULL, NULL, &model, &error));
    duckvep_model_close(model);
    model = NULL;

    for (orientation = 0u; orientation < 2u; orientation++) {
        uint64_t mask;

        strand[0] = orientation == 0u ? 1 : -1;
        duckvep_effect_ctx_fill(
            &tx, &exons, 0u, 0u, 106u, 106u, 0u,
            DUCKVEP_DEFAULT_SPLICE_REGION_EXONIC,
            DUCKVEP_DEFAULT_SPLICE_REGION_INTRONIC, &ctx);
        duckvep_effect_ctx_finalize(&ctx);
        mask = duckvep_effect_eval(ctx.pre_bits);
        ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_MATURE_MIRNA), mask);

        duckvep_effect_ctx_fill(
            &tx, &exons, 0u, 0u, 102u, 102u, 0u,
            DUCKVEP_DEFAULT_SPLICE_REGION_EXONIC,
            DUCKVEP_DEFAULT_SPLICE_REGION_INTRONIC, &ctx);
        duckvep_effect_ctx_finalize(&ctx);
        mask = duckvep_effect_eval(ctx.pre_bits);
        ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT_EXON),
                  mask);

        /* VEP represents an insertion after P as the reversed feature span
         * (P+1,P). At the mature segment's final base, it therefore remains a
         * generic non-coding exon event rather than a mature-miRNA event. */
        ASSERT(duckvep_event_prepare_small(
            mature_end[0], insertion_ref, 1u, insertion_alt, 2u, &event));
        duckvep_effect_ctx_fill(
            &tx, &exons, 0u, 0u, mature_end[0], mature_end[0], 1u,
            DUCKVEP_DEFAULT_SPLICE_REGION_EXONIC,
            DUCKVEP_DEFAULT_SPLICE_REGION_INTRONIC, &ctx);
        duckvep_effect_ctx_apply_event(&tx, &ctx, &event);
        duckvep_effect_ctx_finalize(&ctx);
        mask = duckvep_effect_eval(ctx.pre_bits);
        ASSERT_EQ(DUCKVEP_SO(DUCKVEP_SO_NON_CODING_TRANSCRIPT_EXON),
                  mask);
    }

    mature_start[0] = 99u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, NULL, NULL, &model, &error));
    mature_start[0] = 105u;
    mature_offset[1] = 2u;
    ASSERT_EQ(DUCKVEP_ERR_MODEL_INVALID,
              duckvep_model_open(&tx, &exons, NULL, NULL, &model, &error));
    PASS();
}

/* VEP Plugins release/116 NMD.pm is a separate prediction from the curated
 * NMD transcript biotype. Pin its executable coordinate thresholds on both
 * strands, including the source's inclusive 51-base penultimate-exon offset. */
static void nmd_point_event(duckvep_event_t *event, uint32_t pos1) {
    memset(event, 0, sizeof *event);
    event->raw_start1 = pos1;
    event->raw_end1 = pos1;
    event->feature_start1 = pos1;
    event->feature_end1 = pos1;
    event->start1 = pos1;
    event->end1 = pos1;
    event->ref_diff_length = 1u;
}

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
    nmd_point_event(&event, 320u);
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_PREDICTED_TRIGGERING, nmd.prediction);
    ASSERT_EQ(0u, nmd.escape_reasons);

    nmd_point_event(&event, 300u);
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_PREDICTED_ESCAPING, nmd.prediction);
    ASSERT_EQ(DUCKVEP_NMD_ESCAPE_EARLY_CDS, nmd.escape_reasons);
    nmd_point_event(&event, 301u);
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_PREDICTED_TRIGGERING, nmd.prediction);

    /* NMD.pm reads TranscriptVariation::cds_end. A pure insertion keeps a
     * reversed CDS interval, so insertion after CDS 101 still has cds_end 101
     * even though the allele JSON renderer reports an expanded 101..102. */
    nmd_point_event(&event, 301u);
    event.feature_start1 = 301u;
    event.feature_end1 = 300u;
    event.insertion_boundary0 = 300u;
    event.start1 = 300u;
    event.interbase = 1u;
    event.anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT;
    event.ref_diff_length = 0u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_ESCAPE_EARLY_CDS, nmd.escape_reasons);
    nmd_point_event(&event, 302u);
    event.feature_start1 = 302u;
    event.feature_end1 = 301u;
    event.insertion_boundary0 = 301u;
    event.start1 = 301u;
    event.interbase = 1u;
    event.anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT;
    event.ref_diff_length = 0u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_PREDICTED_TRIGGERING, nmd.prediction);

    /* Definedness, not nonzero coordinates, gates NMD.pm. An insertion before
     * CDS base 1 projects to TranscriptVariation CDS 1..0 and therefore escapes. */
    nmd_point_event(&event, 100u);
    event.feature_start1 = 100u;
    event.feature_end1 = 99u;
    event.insertion_boundary0 = 99u;
    event.start1 = 99u;
    event.interbase = 1u;
    event.anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT;
    event.ref_diff_length = 0u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_PREDICTED_ESCAPING, nmd.prediction);
    ASSERT_EQ(DUCKVEP_NMD_ESCAPE_EARLY_CDS, nmd.escape_reasons);

    nmd_point_event(&event, 348u);
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_ESCAPE_PENULTIMATE_EXON_END,
              nmd.escape_reasons);
    nmd_point_event(&event, 347u);
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_PREDICTED_TRIGGERING, nmd.prediction);

    nmd_point_event(&event, 500u);
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_ESCAPE_LAST_EXON, nmd.escape_reasons);

    excnt[0] = 1u;
    exon_end[0] = 599u;
    cdna_end[0] = 500u;
    nmd_point_event(&event, 300u);
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_ESCAPE_INTRONLESS |
              DUCKVEP_NMD_ESCAPE_LAST_EXON, nmd.escape_reasons);
    excnt[0] = 3u;
    exon_end[0] = 199u;
    cdna_end[0] = 100u;

    nmd_point_event(&event, 320u);
    duckvep_nmd_predict(&tx, &ex, 0u, &event,
                        DUCKVEP_SO(DUCKVEP_SO_MISSENSE), NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_NOT_APPLICABLE, nmd.prediction);

    tx.exon_offset = NULL;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_UNRESOLVED, nmd.prediction);
    tx.exon_offset = exoff;

    cds_start[0] = 0u;
    cds_end[0] = 0u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event,
                        DUCKVEP_SO(DUCKVEP_SO_SPLICE_DONOR), NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_UNRESOLVED, nmd.prediction);
    cds_start[0] = 100u;
    cds_end[0] = 599u;

    nmd_point_event(&event, 250u);
    duckvep_nmd_predict(&tx, &ex, 0u, &event,
                        DUCKVEP_SO(DUCKVEP_SO_SPLICE_DONOR), NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_UNRESOLVED, nmd.prediction);

    strand[0] = -1;
    ex.start1 = reverse_exon_start;
    ex.end1 = reverse_exon_end;
    nmd_point_event(&event, 351u);
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_ESCAPE_PENULTIMATE_EXON_END,
              nmd.escape_reasons);
    nmd_point_event(&event, 352u);
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_PREDICTED_TRIGGERING, nmd.prediction);

    nmd_point_event(&event, 399u);
    event.feature_start1 = 399u;
    event.feature_end1 = 398u;
    event.insertion_boundary0 = 398u;
    event.start1 = 398u;
    event.interbase = 1u;
    event.anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT;
    event.ref_diff_length = 0u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_ESCAPE_EARLY_CDS, nmd.escape_reasons);
    nmd_point_event(&event, 398u);
    event.feature_start1 = 398u;
    event.feature_end1 = 397u;
    event.insertion_boundary0 = 397u;
    event.start1 = 397u;
    event.interbase = 1u;
    event.anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT;
    event.ref_diff_length = 0u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_PREDICTED_TRIGGERING, nmd.prediction);

    nmd_point_event(&event, 600u);
    event.feature_start1 = 600u;
    event.feature_end1 = 599u;
    event.insertion_boundary0 = 599u;
    event.start1 = 599u;
    event.interbase = 1u;
    event.anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT;
    event.ref_diff_length = 0u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_PREDICTED_ESCAPING, nmd.prediction);
    ASSERT_EQ(DUCKVEP_NMD_ESCAPE_EARLY_CDS, nmd.escape_reasons);

    /* NMD.pm projects the complete VariationFeature and reads its genomic end,
     * not the smaller mismatch island used to edit the CDS. These equal-length
     * padded features pin both coordinate-authority differences on each strand. */
    strand[0] = 1;
    ex.start1 = exon_start;
    ex.end1 = exon_end;
    nmd_point_event(&event, 300u);
    event.raw_end1 = 301u;
    event.feature_end1 = 301u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_PREDICTED_TRIGGERING, nmd.prediction);
    ASSERT_EQ(0u, nmd.escape_reasons);

    nmd_point_event(&event, 347u);
    event.raw_end1 = 348u;
    event.feature_end1 = 348u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_ESCAPE_PENULTIMATE_EXON_END,
              nmd.escape_reasons);

    strand[0] = -1;
    ex.start1 = reverse_exon_start;
    ex.end1 = reverse_exon_end;
    nmd_point_event(&event, 351u);
    event.raw_end1 = 352u;
    event.feature_end1 = 352u;
    duckvep_nmd_predict(&tx, &ex, 0u, &event, stop_gained, NULL, &nmd);
    ASSERT_EQ(DUCKVEP_NMD_PREDICTED_TRIGGERING, nmd.prediction);
    ASSERT_EQ(0u, nmd.escape_reasons);
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

TEST generated_effect_lookup_matches_rule_interpreter(void) {
    unsigned a;
    unsigned b;
    uint64_t all = 0u;

    /* The generator requires each consequence rule to consume one unique bit,
     * forbids compound/negative masks, and partitions rules into three suppression
     * groups. Within one group evaluation is an OR homomorphism; between groups a
     * single pair exercises the complete suppression relation. Consequently zero,
     * all singletons, all pairs, and the all-bits case are a complete basis for any
     * input subset, not merely sampled masks. Apply the same basis after the
     * structural VariationFeature-class mask. */
    ASSERT_EQ(duckvep_effect_eval_reference(0u), duckvep_effect_eval(0u));
    ASSERT_EQ(duckvep_effect_eval_structural_reference(0u),
              duckvep_effect_eval_structural(0u));
    for (a = 0u; a < (unsigned)DUCKVEP_PRE_BIT_COUNT; a++) {
        uint64_t one = DUCKVEP_PRE(a);
        all |= one;
        ASSERT_EQ(duckvep_effect_eval_reference(one),
                  duckvep_effect_eval(one));
        ASSERT_EQ(duckvep_effect_eval_structural_reference(one),
                  duckvep_effect_eval_structural(one));
        for (b = 0u; b < (unsigned)DUCKVEP_PRE_BIT_COUNT; b++) {
            uint64_t pair = one | DUCKVEP_PRE(b);
            ASSERT_EQ(duckvep_effect_eval_reference(pair),
                      duckvep_effect_eval(pair));
            ASSERT_EQ(duckvep_effect_eval_structural_reference(pair),
                      duckvep_effect_eval_structural(pair));
        }
    }
    ASSERT_EQ(duckvep_effect_eval_reference(all),
              duckvep_effect_eval(all));
    ASSERT_EQ(duckvep_effect_eval_structural_reference(all),
              duckvep_effect_eval_structural(all));
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
    ASSERT(duckvep_sv_metadata_valid(DUCKVEP_SV_TANDEM_REPEAT,
                                     DUCKVEP_COPY_CHANGE_UNKNOWN));
    ASSERT(duckvep_sv_metadata_valid(DUCKVEP_SV_TANDEM_REPEAT,
                                     DUCKVEP_COPY_CHANGE_GAIN));

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
    ASSERT(!duckvep_sv_metadata_valid(DUCKVEP_SV_TANDEM_REPEAT,
                                      DUCKVEP_COPY_CHANGE_LOSS));
    ASSERT(!duckvep_sv_metadata_valid((duckvep_sv_type_t)255,
                                      DUCKVEP_COPY_CHANGE_UNKNOWN));
    ASSERT(!duckvep_sv_metadata_valid(DUCKVEP_SV_CNV,
                                      (duckvep_copy_change_t)255));

    ASSERT(duckvep_sv_geometry_valid(DUCKVEP_SV_DELETION, 100u, 200u));
    ASSERT(duckvep_sv_geometry_valid(DUCKVEP_SV_INSERTION, 100u, 100u));
    ASSERT(duckvep_sv_geometry_valid(DUCKVEP_SV_BREAKEND, 100u, 100u));
    ASSERT(!duckvep_sv_geometry_valid(DUCKVEP_SV_INSERTION, 100u, 101u));
    ASSERT(!duckvep_sv_geometry_valid(DUCKVEP_SV_BREAKEND, 100u, 101u));
    ASSERT(!duckvep_sv_geometry_valid(
        DUCKVEP_SV_INSERTION, UINT32_MAX, UINT32_MAX));
    ASSERT(!duckvep_sv_geometry_valid(DUCKVEP_SV_DELETION, 0u, 1u));
    ASSERT(!duckvep_sv_geometry_valid(DUCKVEP_SV_DELETION, 2u, 1u));
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
    sv = duckvep_sv_effect_fill(NULL, NULL, 0u, 0, &event, &region);
    ASSERT(sv.copy_number_gain);
    ASSERT(sv.insertion);
    ASSERT(sv.feature_amplification);
    ASSERT(!sv.deletion);
    ASSERT(!sv.feature_ablation);
    duckvep_effect_ctx_apply_sv(&ctx, &sv);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INSERTION)) != 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_FEATURE_AMPLIFICATION)) != 0u);

    /* Structural tandem repeats remain distinguishable in the event ABI but
     * VEP maps them to the tandem-duplication gain/insertion facts. */
    memset(&region, 0, sizeof region);
    memset(&ctx, 0, sizeof ctx);
    event.sv_type = (uint8_t)DUCKVEP_SV_TANDEM_REPEAT;
    event.copy_change = (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN;
    region.complete_overlap_feature = 1u;
    sv = duckvep_sv_effect_fill(NULL, NULL, 0u, 0, &event, &region);
    ASSERT(sv.copy_number_gain);
    ASSERT(sv.insertion);
    ASSERT(sv.feature_amplification);
    ASSERT(!sv.deletion);
    duckvep_effect_ctx_apply_sv(&ctx, &sv);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INSERTION)) != 0u);
    ASSERT((ctx.pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_FEATURE_AMPLIFICATION)) != 0u);

    memset(&region, 0, sizeof region);
    memset(&ctx, 0, sizeof ctx);
    event.sv_type = (uint8_t)DUCKVEP_SV_CNV;
    event.copy_change = (uint8_t)DUCKVEP_COPY_CHANGE_LOSS;
    region.within_cdna = 1u;
    region.partial_overlap_feature = 1u;
    sv = duckvep_sv_effect_fill(NULL, NULL, 0u, 0, &event, &region);
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
    sv = duckvep_sv_effect_fill(NULL, NULL, 0u, 0, &event, &region);
    ASSERT(!sv.copy_number_gain);
    ASSERT(!sv.copy_number_loss);
    ASSERT(!sv.feature_amplification);
    ASSERT(!sv.feature_ablation);
    PASS();
}

/* VEP's structural start predicates are not a plain genomic overlap. Its
 * shared _overlaps_start_codon guard first requires both feature endpoints to
 * map to cDNA. Keep that two-endpoint mapping rule explicit: the same span emits the
 * start-lost/start-retained pair when its far endpoint is exonic, but not when
 * that endpoint falls in the intron between the two exons. */
TEST sv_start_codon_requires_two_exonic_endpoints(void) {
    static const uint16_t chrom[1] = {0u};
    static const uint32_t tx_start[1] = {100u};
    static const uint32_t tx_end[1] = {300u};
    static const int8_t strand[1] = {1};
    static const uint64_t flags[1] = {0u};
    static const uint32_t exon_offset[1] = {0u};
    static const uint16_t exon_count[1] = {2u};
    static const uint32_t cds_start[1] = {120u};
    static const uint32_t cds_end[1] = {280u};
    static const uint32_t exon_start[2] = {100u, 250u};
    static const uint32_t exon_end[2] = {150u, 300u};
    static const uint32_t cdna_start[2] = {1u, 52u};
    static const uint32_t cdna_end[2] = {51u, 102u};
    static const int8_t phase[2] = {0, 0};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_region_state_t region;
    duckvep_event_t event;
    duckvep_sv_effect_t effect;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&region, 0, sizeof region);
    memset(&event, 0, sizeof event);
    tx.chrom_id = chrom; tx.start1 = tx_start; tx.end1 = tx_end;
    tx.strand = strand; tx.flags = flags; tx.exon_offset = exon_offset;
    tx.exon_count = exon_count; tx.cds_start1 = cds_start;
    tx.cds_end1 = cds_end; tx.transcript_count = 1u;
    exons.start1 = exon_start; exons.end1 = exon_end;
    exons.cdna_start1 = cdna_start; exons.cdna_end1 = cdna_end;
    exons.phase = phase; exons.end_phase = phase; exons.exon_count = 2u;
    event.kind = (uint8_t)DUCKVEP_KIND_SV;
    event.sv_type = (uint8_t)DUCKVEP_SV_INVERSION;
    event.copy_change = (uint8_t)DUCKVEP_COPY_CHANGE_NEUTRAL;
    event.start1 = 110u; event.end1 = 260u;
    event.feature_start1 = 110u; event.feature_end1 = 260u;
    region.within_cdna = 1u;
    region.overlaps_cds = 1u;

    effect = duckvep_sv_effect_fill(&tx, &exons, 0u, 0, &event, &region);
    ASSERT(effect.start_lost);
    ASSERT(effect.start_retained);

    event.end1 = 200u;
    event.feature_end1 = 200u;
    effect = duckvep_sv_effect_fill(&tx, &exons, 0u, 0, &event, &region);
    ASSERT(!effect.start_lost);
    ASSERT(!effect.start_retained);
    PASS();
}

/* VEP's structural start-retained predicate uses projected cDNA, but its
 * start-lost branch adds a contiguous genomic coding-region test. A codon split
 * across exons can therefore be retained without being lost. Reversed insertion
 * geometry at that exon edge is clamped to its one exonic flank. */
TEST sv_split_start_codon_and_edge_insertion(void) {
    static const uint16_t chrom[1] = {0u};
    static const uint32_t tx_start[1] = {100u};
    static const uint32_t tx_end[1] = {250u};
    static const int8_t strand[1] = {1};
    static const uint64_t flags[1] = {0u};
    static const uint32_t exon_offset[1] = {0u};
    static const uint16_t exon_count[1] = {2u};
    static const uint32_t cds_start[1] = {120u};
    static const uint32_t cds_end[1] = {240u};
    static const uint32_t exon_start[2] = {100u, 200u};
    static const uint32_t exon_end[2] = {121u, 250u};
    static const uint32_t cdna_start[2] = {1u, 23u};
    static const uint32_t cdna_end[2] = {22u, 73u};
    static const int8_t phase[2] = {0, 1};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_region_state_t region;
    duckvep_event_t event;
    duckvep_sv_effect_t effect;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&region, 0, sizeof region);
    memset(&event, 0, sizeof event);
    tx.chrom_id = chrom; tx.start1 = tx_start; tx.end1 = tx_end;
    tx.strand = strand; tx.flags = flags; tx.exon_offset = exon_offset;
    tx.exon_count = exon_count; tx.cds_start1 = cds_start;
    tx.cds_end1 = cds_end; tx.transcript_count = 1u;
    exons.start1 = exon_start; exons.end1 = exon_end;
    exons.cdna_start1 = cdna_start; exons.cdna_end1 = cdna_end;
    exons.phase = phase; exons.end_phase = phase; exons.exon_count = 2u;
    region.within_cdna = 1u;
    region.overlaps_cds = 1u;
    event.kind = (uint8_t)DUCKVEP_KIND_SV;
    event.sv_type = (uint8_t)DUCKVEP_SV_INVERSION;
    event.copy_change = (uint8_t)DUCKVEP_COPY_CHANGE_NEUTRAL;

    event.start1 = 200u; event.end1 = 205u;
    event.feature_start1 = 200u; event.feature_end1 = 205u;
    effect = duckvep_sv_effect_fill(&tx, &exons, 0u, 0, &event, &region);
    ASSERT(effect.start_retained);
    ASSERT(!effect.start_lost);

    event.sv_type = (uint8_t)DUCKVEP_SV_INSERTION;
    event.start1 = 121u; event.end1 = 121u;
    event.feature_start1 = 122u; event.feature_end1 = 121u;
    event.insertion_boundary0 = 121u;
    event.interbase = 1u;
    event.anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT;
    effect = duckvep_sv_effect_fill(&tx, &exons, 0u, 0, &event, &region);
    ASSERT(effect.start_retained);
    ASSERT(effect.start_lost);
    PASS();
}

/* Mapper::map_insert exposes an insertion as a reversed cDNA interval on both
 * strands. On a reverse transcript, the genomic interval itself would project
 * in ascending cDNA order; sorting those two flank projections would therefore
 * invent a start-codon overlap immediately before translation begins. */
TEST sv_insertion_start_uses_reversed_cdna_interval(void) {
    static const uint16_t chrom[1] = {0u};
    static const uint32_t tx_start[1] = {100u};
    static const uint32_t tx_end[1] = {200u};
    static const int8_t strand[1] = {-1};
    static const uint64_t flags[1] = {0u};
    static const uint32_t exon_offset[1] = {0u};
    static const uint16_t exon_count[1] = {1u};
    static const uint32_t cds_start[1] = {120u};
    static const uint32_t cds_end[1] = {180u};
    static const uint32_t exon_start[1] = {100u};
    static const uint32_t exon_end[1] = {200u};
    static const uint32_t cdna_start[1] = {1u};
    static const uint32_t cdna_end[1] = {101u};
    static const int8_t phase[1] = {0};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_region_state_t region;
    duckvep_event_t event;
    duckvep_sv_effect_t effect;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&region, 0, sizeof region);
    memset(&event, 0, sizeof event);
    tx.chrom_id = chrom; tx.start1 = tx_start; tx.end1 = tx_end;
    tx.strand = strand; tx.flags = flags; tx.exon_offset = exon_offset;
    tx.exon_count = exon_count; tx.cds_start1 = cds_start;
    tx.cds_end1 = cds_end; tx.transcript_count = 1u;
    exons.start1 = exon_start; exons.end1 = exon_end;
    exons.cdna_start1 = cdna_start; exons.cdna_end1 = cdna_end;
    exons.phase = phase; exons.end_phase = phase; exons.exon_count = 1u;
    region.within_cdna = 1u;
    event.kind = (uint8_t)DUCKVEP_KIND_SV;
    event.sv_type = (uint8_t)DUCKVEP_SV_INSERTION;
    event.copy_change = (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN;
    event.interbase = 1u;
    event.anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT;

    /* Between genomic 180 and 181: immediately before the reverse-strand CDS. */
    event.start1 = 180u; event.end1 = 180u;
    event.feature_start1 = 181u; event.feature_end1 = 180u;
    event.insertion_boundary0 = 180u;
    effect = duckvep_sv_effect_fill(&tx, &exons, 0u, 0u, &event, &region);
    ASSERT(!effect.start_lost);
    ASSERT(!effect.start_retained);

    /* Between 179 and 180: after the first coding base in transcript order. */
    event.start1 = 179u; event.end1 = 179u;
    event.feature_start1 = 180u; event.feature_end1 = 179u;
    event.insertion_boundary0 = 179u;
    effect = duckvep_sv_effect_fill(&tx, &exons, 0u, 0u, &event, &region);
    ASSERT(effect.start_lost);
    ASSERT(effect.start_retained);
    PASS();
}

/* VEP 116's structural frameshift and inframe-deletion predicates do not use
 * the same mapper guard. A deletion wholly inside an exon may be a frameshift
 * while crossing from CDS into UTR; the divisible-by-three case requires one
 * CDS Coordinate with no mapper Gap and therefore stays generic coding. */
TEST sv_inframe_deletion_rejects_cds_to_utr_span(void) {
    static const uint16_t chrom[1] = {0u};
    static const uint32_t tx_start[1] = {100u};
    static const uint32_t tx_end[1] = {300u};
    static const int8_t strand[1] = {1};
    static const uint64_t flags[1] = {0u};
    static const uint32_t exon_offset[1] = {0u};
    static const uint16_t exon_count[1] = {1u};
    static const uint32_t cds_start[1] = {120u};
    static const uint32_t cds_end[1] = {279u};
    static const uint32_t exon_start[1] = {100u};
    static const uint32_t exon_end[1] = {300u};
    static const uint32_t cdna_start[1] = {1u};
    static const uint32_t cdna_end[1] = {201u};
    static const int8_t phase[1] = {0};
    duckvep_transcript_model_t tx;
    duckvep_exon_model_t exons;
    duckvep_region_state_t region;
    duckvep_event_t event;
    duckvep_sv_effect_t effect;

    memset(&tx, 0, sizeof tx);
    memset(&exons, 0, sizeof exons);
    memset(&region, 0, sizeof region);
    memset(&event, 0, sizeof event);
    tx.chrom_id = chrom; tx.start1 = tx_start; tx.end1 = tx_end;
    tx.strand = strand; tx.flags = flags; tx.exon_offset = exon_offset;
    tx.exon_count = exon_count; tx.cds_start1 = cds_start;
    tx.cds_end1 = cds_end; tx.transcript_count = 1u;
    exons.start1 = exon_start; exons.end1 = exon_end;
    exons.cdna_start1 = cdna_start; exons.cdna_end1 = cdna_end;
    exons.phase = phase; exons.end_phase = phase; exons.exon_count = 1u;
    event.kind = (uint8_t)DUCKVEP_KIND_SV;
    event.sv_type = (uint8_t)DUCKVEP_SV_DELETION;
    event.copy_change = (uint8_t)DUCKVEP_COPY_CHANGE_LOSS;
    region.overlaps_cds = 1u;

    event.start1 = 260u; event.end1 = 289u; /* 30 bp, CDS + UTR. */
    effect = duckvep_sv_effect_fill(&tx, &exons, 0u, 0, &event, &region);
    ASSERT(!effect.frameshift);
    ASSERT(!effect.inframe_deletion);

    event.end1 = 290u; /* 31 bp, same mapper shape. */
    effect = duckvep_sv_effect_fill(&tx, &exons, 0u, 0, &event, &region);
    ASSERT(effect.frameshift);
    ASSERT(!effect.inframe_deletion);

    event.start1 = 240u; event.end1 = 269u; /* 30 bp, wholly CDS. */
    effect = duckvep_sv_effect_fill(&tx, &exons, 0u, 0, &event, &region);
    ASSERT(!effect.frameshift);
    ASSERT(effect.inframe_deletion);

    /* partial_codon is event-relative. This transcript's 160-base CDS ends in
     * one incomplete base. A deletion beginning in the preceding complete
     * codon still loses the stop; one beginning at CDS base 160 does not. */
    event.start1 = 277u; event.end1 = 279u;
    event.feature_start1 = 277u; event.feature_end1 = 279u;
    effect = duckvep_sv_effect_fill(&tx, &exons, 0u, 160u, &event, &region);
    ASSERT(effect.stop_lost);
    event.start1 = 279u; event.end1 = 279u;
    event.feature_start1 = 279u; event.feature_end1 = 279u;
    effect = duckvep_sv_effect_fill(&tx, &exons, 0u, 160u, &event, &region);
    ASSERT(!effect.stop_lost);
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

static duckvep_impact_t so_impact_metadata_oracle(uint64_t mask) {
    duckvep_impact_t best = DUCKVEP_IMPACT_MODIFIER;
    unsigned bit;

    for (bit = 0u; bit < 64u; bit++) {
        duckvep_impact_t impact;

        if ((mask & (UINT64_C(1) << bit)) == 0u) continue;
        impact = duckvep_so_bit_impact((duckvep_so_bit_t)bit);
        if (impact > best) best = impact;
    }
    return best;
}

/* One generated mask lookup replaces the former set-bit scan in every compact
 * result row.  Exhaustive singleton and pair coverage proves the maximum-impact
 * reduction against the independently indexed generated metadata; additional
 * bits cannot change a maximum except through one of these pair relations. */
TEST so_impact_masks_match_metadata_oracle(void) {
    unsigned first;
    unsigned second;

    ASSERT_EQ(so_impact_metadata_oracle(0u), duckvep_so_impact(0u));
    for (first = 0u; first < 64u; first++) {
        uint64_t first_mask = UINT64_C(1) << first;

        ASSERT_EQ(so_impact_metadata_oracle(first_mask),
                  duckvep_so_impact(first_mask));
        for (second = first; second < 64u; second++) {
            uint64_t mask = first_mask | (UINT64_C(1) << second);

            ASSERT_EQ(so_impact_metadata_oracle(mask),
                      duckvep_so_impact(mask));
        }
    }
    PASS();
}
