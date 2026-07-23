/*
 * duckvep_classify.c — span-aware structural and splice classification.
 * No DuckDB/htslib/Arrow; no allocation. Pure per-pair arithmetic.
 */
#include "duckvep_classify.h"

#include <string.h>

#if defined(__GNUC__) || defined(__clang__)
#define DUCKVEP_CLASSIFY_INLINE __attribute__((always_inline)) inline
#else
#define DUCKVEP_CLASSIFY_INLINE inline
#endif

static int span_overlaps_u32(uint32_t as, uint32_t ae,
                             uint32_t bs, uint32_t be) {
    return ae >= bs && as <= be;
}

static int span_overlaps_i64(int64_t start1, int64_t end1,
                             int64_t lo, int64_t hi) {
    return hi >= lo && end1 >= lo && start1 <= hi;
}

static uint32_t min_u32(uint32_t a, uint32_t b) { return a < b ? a : b; }
static uint32_t max_u32(uint32_t a, uint32_t b) { return a > b ? a : b; }

static int gap_between_exons(const duckvep_exon_model_t *exons,
                             size_t a_idx,
                             size_t b_idx,
                             uint32_t *gap_start,
                             uint32_t *gap_end) {
    uint32_t a_s = exons->start1[a_idx];
    uint32_t a_e = exons->end1[a_idx];
    uint32_t b_s = exons->start1[b_idx];
    uint32_t b_e = exons->end1[b_idx];

    if (a_e < b_s && a_e != UINT32_MAX && b_s > 0u) {
        *gap_start = a_e + 1u;
        *gap_end = b_s - 1u;
    } else if (b_e < a_s && b_e != UINT32_MAX && a_s > 0u) {
        *gap_start = b_e + 1u;
        *gap_end = a_s - 1u;
    } else {
        return 0;
    }
    return *gap_start <= *gap_end;
}

/* BaseTranscriptVariation::_create_intron_trees marks an intron as a
 * frameshift intron when abs(end - start) <= 12. Normal prepared exon gaps have
 * start <= end, so this is an inclusive intron length of at most 13 bases. */
static int gap_is_frameshift_intron(uint32_t gap_start, uint32_t gap_end) {
    return gap_end >= gap_start &&
        gap_end - gap_start <= DUCKVEP_VEP_FRAMESHIFT_INTRON_MAX_SPAN;
}

enum duckvep_splice_accum_bit {
    DUCKVEP_SPLICE_ACCUM_START_SS          = 1u << 0,
    DUCKVEP_SPLICE_ACCUM_END_SS            = 1u << 1,
    DUCKVEP_SPLICE_ACCUM_FIFTH             = 1u << 2,
    DUCKVEP_SPLICE_ACCUM_FIFTH_REV         = 1u << 3,
    DUCKVEP_SPLICE_ACCUM_DREG              = 1u << 4,
    DUCKVEP_SPLICE_ACCUM_DREG_REV          = 1u << 5,
    DUCKVEP_SPLICE_ACCUM_PPT               = 1u << 6,
    DUCKVEP_SPLICE_ACCUM_PPT_REV           = 1u << 7,
    DUCKVEP_SPLICE_ACCUM_INTRONIC          = 1u << 8,
    DUCKVEP_SPLICE_ACCUM_FRAMESHIFT_INTRON = 1u << 9,
    DUCKVEP_SPLICE_ACCUM_REGION            = 1u << 10
};

/* These are independent Boolean facts accumulated across the introns selected
 * by VEP's cache geometry.  Keep them as one hot bitset rather than eleven
 * `int`s: every point/span classifier creates this object, while the final
 * strand-relative projection below remains the sole public-state authority. */
struct duckvep_splice_accum {
    uint16_t bits;
};

static int splice_region_overlaps_gap(int64_t start1, int64_t end1,
                                      int64_t intron_start,
                                      int64_t intron_end,
                                      uint8_t interbase,
                                      uint32_t splice_exonic,
                                      uint32_t splice_intronic) {
    int overlaps = 0;

    /* The first two intronic bases have the essential donor/acceptor terms.
     * splice_region_variant starts at base three and extends through the
     * caller-selected intronic reach. */
    if (splice_intronic > 2u) {
        int64_t reach = (int64_t)splice_intronic - 1;
        overlaps = span_overlaps_i64(start1, end1,
                                     intron_start + 2,
                                     intron_start + reach) ||
                   span_overlaps_i64(start1, end1,
                                     intron_end - reach,
                                     intron_end - 2);
    }
    if (!overlaps && splice_exonic != 0u) {
        overlaps = span_overlaps_i64(start1, end1,
                                     intron_start - (int64_t)splice_exonic,
                                     intron_start - 1) ||
                   span_overlaps_i64(start1, end1,
                                     intron_end + 1,
                                     intron_end + (int64_t)splice_exonic);
    }
    /* VEP treats insertions at these zone edges specially because its feature
     * interval is reversed. Keep those source-compatible edge cases while the
     * ordinary windows above honor the configured reaches. */
    if (!overlaps && interbase &&
        ((splice_exonic != 0u &&
          (start1 == intron_start || end1 == intron_end)) ||
         (splice_intronic > 2u &&
          (start1 == intron_start + 2 || end1 == intron_end - 2)))) {
        overlaps = 1;
    }
    return overlaps;
}

/* Accumulate VEP's splice predicates for one already-resolved intron. Span,
 * mismatch-island, and sorted-point traversal deliberately share this code:
 * each caller may choose fewer gaps, but none owns a second definition of the
 * biology. */
static DUCKVEP_CLASSIFY_INLINE int splice_accum_add_gap_bounds(
    uint32_t gap_start,
    uint32_t gap_end,
    int64_t irs,
    int64_t ire,
    int64_t lo_s,
    int64_t hi_s,
    uint8_t interbase,
    int intron_cached,
    int boundary_cached,
    uint32_t splice_exonic,
    uint32_t splice_intronic,
    struct duckvep_splice_accum *a,
    int64_t *reach_lo,
    int64_t *reach_hi) {
    int64_t is;
    int64_t ie;

    is = (int64_t)gap_start;
    ie = (int64_t)gap_end;
    if (reach_lo != NULL || reach_hi != NULL) {
        int64_t intronic_reach = splice_intronic > 1u
            ? (int64_t)splice_intronic - 1 : 0;
        int64_t internal_reach = intronic_reach > 16
            ? intronic_reach : 16;
        int64_t lo = is - (int64_t)splice_exonic;
        int64_t hi = ie + (int64_t)splice_exonic;
        int64_t other_lo = ie - internal_reach;
        int64_t other_hi = is + internal_reach;

        if (reach_lo != NULL) *reach_lo = lo < other_lo ? lo : other_lo;
        if (reach_hi != NULL) *reach_hi = hi > other_hi ? hi : other_hi;
    }

    /* VEP checks this before every ordinary intron/splice predicate. A
     * differing region that overlaps a short assembly-gap intron is treated as
     * coding-but-unprojectable, not as donor, acceptor, splice region, PPT, or
     * intron_variant. */
    if (gap_is_frameshift_intron(gap_start, gap_end) &&
        span_overlaps_i64(irs, ire, is, ie)) {
        a->bits |= DUCKVEP_SPLICE_ACCUM_FRAMESHIFT_INTRON;
        return 1;
    }

    if (boundary_cached) {
        if (span_overlaps_i64(irs, ire, is, is + 1))
            a->bits |= DUCKVEP_SPLICE_ACCUM_START_SS;
        if (span_overlaps_i64(irs, ire, ie - 1, ie))
            a->bits |= DUCKVEP_SPLICE_ACCUM_END_SS;
        if (span_overlaps_i64(irs, ire, is + 4, is + 4))
            a->bits |= DUCKVEP_SPLICE_ACCUM_FIFTH;
        if (span_overlaps_i64(irs, ire, ie - 4, ie - 4))
            a->bits |= DUCKVEP_SPLICE_ACCUM_FIFTH_REV;
        if (span_overlaps_i64(irs, ire, is + 2, is + 5))
            a->bits |= DUCKVEP_SPLICE_ACCUM_DREG;
        if (span_overlaps_i64(irs, ire, ie - 5, ie - 2))
            a->bits |= DUCKVEP_SPLICE_ACCUM_DREG_REV;
        if (splice_region_overlaps_gap(irs, ire, is, ie, interbase,
                                       splice_exonic, splice_intronic))
            a->bits |= DUCKVEP_SPLICE_ACCUM_REGION;
    }
    if (intron_cached || boundary_cached) {
        if (span_overlaps_i64(lo_s, hi_s, ie - 16, ie - 2))
            a->bits |= DUCKVEP_SPLICE_ACCUM_PPT;
        if (span_overlaps_i64(lo_s, hi_s, is + 2, is + 16))
            a->bits |= DUCKVEP_SPLICE_ACCUM_PPT_REV;
    }
    if (intron_cached &&
        (span_overlaps_i64(irs, ire, is + 2, ie - 2) ||
         (interbase && (irs == is + 2 || ire == ie - 2)))) {
        a->bits |= DUCKVEP_SPLICE_ACCUM_INTRONIC;
    }
    return 1;
}

static int splice_accum_add_gap(const duckvep_exon_model_t *exons,
                                size_t left_idx, size_t right_idx,
                                int64_t irs, int64_t ire,
                                int64_t lo_s, int64_t hi_s,
                                uint8_t interbase,
                                int intron_cached,
                                int boundary_cached,
                                uint32_t splice_exonic,
                                uint32_t splice_intronic,
                                struct duckvep_splice_accum *a,
                                int64_t *reach_lo, int64_t *reach_hi) {
    uint32_t gap_start;
    uint32_t gap_end;

    if (!gap_between_exons(exons, left_idx, right_idx,
                           &gap_start, &gap_end)) {
        return 0;
    }
    return splice_accum_add_gap_bounds(
        gap_start, gap_end, irs, ire, lo_s, hi_s, interbase,
        intron_cached, boundary_cached, splice_exonic, splice_intronic,
        a, reach_lo, reach_hi);
}

static duckvep_splice_state_t splice_accum_finish(
    const struct duckvep_splice_accum *a, int fwd) {
    duckvep_splice_state_t st;

    memset(&st, 0, sizeof st);
    st.splice_donor = (uint8_t)((a->bits & (fwd
        ? DUCKVEP_SPLICE_ACCUM_START_SS
        : DUCKVEP_SPLICE_ACCUM_END_SS)) != 0u);
    st.splice_acceptor = (uint8_t)((a->bits & (fwd
        ? DUCKVEP_SPLICE_ACCUM_END_SS
        : DUCKVEP_SPLICE_ACCUM_START_SS)) != 0u);
    st.splice_donor_5th = (uint8_t)((a->bits & (fwd
        ? DUCKVEP_SPLICE_ACCUM_FIFTH
        : DUCKVEP_SPLICE_ACCUM_FIFTH_REV)) != 0u);
    st.splice_donor_region =
        (uint8_t)((a->bits & (fwd ? DUCKVEP_SPLICE_ACCUM_DREG
                                  : DUCKVEP_SPLICE_ACCUM_DREG_REV)) != 0u &&
                  !st.splice_donor_5th);
    st.splice_polypyrimidine = (uint8_t)((a->bits & (fwd
        ? DUCKVEP_SPLICE_ACCUM_PPT
        : DUCKVEP_SPLICE_ACCUM_PPT_REV)) != 0u);
    st.intronic = (uint8_t)((a->bits &
        DUCKVEP_SPLICE_ACCUM_INTRONIC) != 0u);
    st.splice_region =
        (uint8_t)((a->bits & DUCKVEP_SPLICE_ACCUM_REGION) != 0u &&
                  !st.splice_donor && !st.splice_acceptor &&
                  !st.splice_donor_region && !st.splice_donor_5th);
    st.any = (uint8_t)(st.splice_donor || st.splice_acceptor ||
                       st.splice_donor_5th || st.splice_donor_region ||
                       st.splice_polypyrimidine || st.splice_region);
    st.within_frameshift_intron =
        (uint8_t)((a->bits & DUCKVEP_SPLICE_ACCUM_FRAMESHIFT_INTRON) != 0u);
    return st;
}

static void splice_accum_add_region(
    const duckvep_exon_model_t *exons,
    size_t                      exon_offset,
    size_t                      exon_count,
    int64_t                     region_start1,
    int64_t                     region_end1,
    int64_t                     feature_min1,
    int64_t                     feature_max1,
    uint32_t                    splice_exonic,
    uint32_t                    splice_intronic,
    struct duckvep_splice_accum *acc) {

    uint8_t interbase = (uint8_t)(region_start1 == region_end1 + 1);
    int64_t lo = region_start1 < region_end1 ? region_start1 : region_end1;
    int64_t hi = region_start1 > region_end1 ? region_start1 : region_end1;
    int64_t cache_exonic = splice_exonic > 3u
        ? (int64_t)splice_exonic : 3;
    int64_t cache_intronic = splice_intronic > 1u
        ? (int64_t)splice_intronic - 1 : 0;
    size_t k;

    /* This is only a superset cache gate. Keep every fixed VEP predicate
     * reachable (the old -3/+7 window), and widen it when the caller asks the
     * generic splice-region predicate to see farther. Predicate coordinates
     * themselves remain exact in the shared intron predicate helper. */
    if (cache_intronic < 7) cache_intronic = 7;

    for (k = 0u; k + 1u < exon_count; k++) {
        uint32_t gap_start;
        uint32_t gap_end;
        int intron_cached;
        int boundary_cached;
        int frameshift_overlap;

        if (!gap_between_exons(exons, exon_offset + k, exon_offset + k + 1u,
                               &gap_start, &gap_end)) {
            continue;
        }
        boundary_cached =
            span_overlaps_i64(feature_min1, feature_max1,
                              (int64_t)gap_start - cache_exonic,
                              (int64_t)gap_start + cache_intronic) ||
            span_overlaps_i64(feature_min1, feature_max1,
                              (int64_t)gap_end - cache_intronic,
                              (int64_t)gap_end + cache_exonic);
        /* With Set::IntervalTree enabled, VEP 116 preselects introns over a
         * three-base flank on either side of the actual intron. The selected
         * list is then cached before _intron_effects visits mismatch islands.
         * A lengthening allele whose REF feature stops just inside the exon can
         * therefore expose an ALT-only island in the intron. Preserve that
         * prefilter here; the later predicate still requires the island itself
         * to overlap the intronic interior [is+2,ie-2]. */
        intron_cached = span_overlaps_i64(
            feature_min1, feature_max1,
            (int64_t)gap_start - 3, (int64_t)gap_end + 3);
        frameshift_overlap =
            gap_is_frameshift_intron(gap_start, gap_end) &&
            span_overlaps_i64(region_start1, region_end1,
                              (int64_t)gap_start, (int64_t)gap_end);
        if (!intron_cached && !boundary_cached && !frameshift_overlap) {
            continue;
        }
        (void)splice_accum_add_gap_bounds(
            gap_start, gap_end,
            region_start1, region_end1, lo, hi, interbase,
            intron_cached, boundary_cached,
            splice_exonic, splice_intronic,
            acc, NULL, NULL);
        if (boundary_cached && !frameshift_overlap) {
            /* VEP 116 assigns, rather than ORs, splice_region for each
             * differing island/intron pair. Later islands can therefore clear
             * an earlier hit; the other splice facts remain accumulators. */
            if (splice_region_overlaps_gap(
                    region_start1, region_end1,
                    (int64_t)gap_start, (int64_t)gap_end, interbase,
                    splice_exonic, splice_intronic)) {
                acc->bits |= DUCKVEP_SPLICE_ACCUM_REGION;
            } else {
                acc->bits &= (uint16_t)~DUCKVEP_SPLICE_ACCUM_REGION;
            }
        }
    }
}

static int overlaps_coarse_exon_reach(uint32_t start1, uint32_t end1,
                                      uint32_t es, uint32_t ee,
                                      uint32_t exonic,
                                      uint32_t intronic) {
    if (exonic > 0u) {
        int64_t left_hi = (int64_t)es + (int64_t)exonic - 1;
        int64_t right_lo = (int64_t)ee - (int64_t)exonic + 1;
        /* The exonic splice reach covers only bases INSIDE the exon. Clamp both zones to
         * [es, ee] so an exon shorter than `exonic` does not spill its reach past the
         * opposite boundary and wrongly flag a flanking intron/region base as splice. */
        if (left_hi > (int64_t)ee) left_hi = (int64_t)ee;
        if (right_lo < (int64_t)es) right_lo = (int64_t)es;
        if (span_overlaps_i64(start1, end1, (int64_t)es, left_hi) ||
            span_overlaps_i64(start1, end1, right_lo, (int64_t)ee)) {
            return 1;
        }
    }
    if (intronic > 0u) {
        if (span_overlaps_i64(start1, end1,
                              (int64_t)es - (int64_t)intronic,
                              (int64_t)es - 1) ||
            span_overlaps_i64(start1, end1,
                              (int64_t)ee + 1,
                              (int64_t)ee + (int64_t)intronic)) {
            return 1;
        }
    }
    return 0;
}

/* Region rules are shared by the exhaustive scan and the sorted contained-span
 * proof. Candidate traversal may change; the biological fact definitions do
 * not. */
static void region_add_exon_overlap(duckvep_region_state_t *st,
                                    uint32_t start1, uint32_t end1,
                                    uint32_t es, uint32_t ee,
                                    uint32_t cds_s, uint32_t cds_e,
                                    int fwd) {
    uint32_t ov_s;
    uint32_t ov_e;

    if (st == NULL || !span_overlaps_u32(start1, end1, es, ee)) return;
    st->within_cdna = 1u;
    st->overlaps_exon = 1u;
    if (cds_s == 0u) return;

    ov_s = max_u32(start1, es);
    ov_e = min_u32(end1, ee);
    if (span_overlaps_u32(ov_s, ov_e, cds_s, cds_e)) {
        st->overlaps_cds = 1u;
    }
    if (ov_s < cds_s) {
        uint32_t utr_e = min_u32(ov_e, cds_s - 1u);
        if (ov_s <= utr_e) {
            if (fwd) st->overlaps_utr5 = 1u;
            else     st->overlaps_utr3 = 1u;
        }
    }
    if (ov_e > cds_e && cds_e != UINT32_MAX) {
        uint32_t utr_s = max_u32(ov_s, cds_e + 1u);
        if (utr_s <= ov_e) {
            if (fwd) st->overlaps_utr3 = 1u;
            else     st->overlaps_utr5 = 1u;
        }
    }
}

static void region_add_frameshift_intron(duckvep_region_state_t *st,
                                         uint32_t start1,
                                         uint32_t end1,
                                         uint32_t cds_s,
                                         uint32_t cds_e,
                                         int fwd) {
    if (st == NULL) return;
    st->within_frameshift_intron = 1u;
    st->within_cdna = 1u;
    if (cds_s == 0u) {
        /* _overlapped_exons uses the transcript-wide 12-base stretch to find
         * candidates, but non_coding_exon_variant then checks the unextended
         * exon coordinates. A point in the short intron is therefore within
         * the transcript cDNA mapper without satisfying the exon predicate.
         * Keep its physical region intronic; _intron_effects independently
         * suppresses the intron_variant predicate for a frameshift intron. */
        st->overlaps_intron = 1u;
        return;
    }
    if (span_overlaps_u32(start1, end1, cds_s, cds_e))
        st->overlaps_cds = 1u;
    if (start1 < cds_s) {
        if (fwd) st->overlaps_utr5 = 1u;
        else     st->overlaps_utr3 = 1u;
    }
    if (end1 > cds_e) {
        if (fwd) st->overlaps_utr3 = 1u;
        else     st->overlaps_utr5 = 1u;
    }
}

static void region_add_vep_endpoint_utr(duckvep_region_state_t *st,
                                        uint32_t start1, uint32_t end1,
                                        uint32_t ts, uint32_t te,
                                        uint32_t cds_s, uint32_t cds_e,
                                        int fwd) {
    int before_coding;
    int after_coding;

    if (st == NULL || cds_s == 0u || !st->within_cdna) {
        return;
    }

    /* VEP applies its four-comparison overlap test even when a UTR interval is
     * inverted at a CDS/transcript endpoint. Keep that source-compatible state
     * in one helper for both traversal paths. Complete transcript overlap does
     * not suppress these predicates: an equal-length uploaded span containing
     * a transcript whose CDS shares both endpoints can therefore report both
     * otherwise-empty UTRs. */
    before_coding = cds_s > 0u && end1 >= ts && start1 <= cds_s - 1u;
    after_coding = cds_e < UINT32_MAX && end1 >= cds_e + 1u && start1 <= te;
    if (fwd) {
        if (before_coding) st->overlaps_utr5 = 1u;
        if (after_coding)  st->overlaps_utr3 = 1u;
    } else {
        if (before_coding) st->overlaps_utr3 = 1u;
        if (after_coding)  st->overlaps_utr5 = 1u;
    }
}

static void region_finish_mask(duckvep_region_state_t *st,
                               uint32_t cds_s, int coarse_splice) {
    if (st == NULL) return;
    if (cds_s == 0u && st->overlaps_exon)
        st->region_mask |= (uint32_t)DUCKVEP_REGION_EXON;
    if (st->overlaps_cds)
        st->region_mask |= (uint32_t)DUCKVEP_REGION_CDS;
    if (st->overlaps_utr5 || st->overlaps_utr3)
        st->region_mask |= (uint32_t)DUCKVEP_REGION_UTR;
    if (st->overlaps_intron)
        st->region_mask |= (uint32_t)DUCKVEP_REGION_INTRON;
    if (coarse_splice)
        st->region_mask |= (uint32_t)DUCKVEP_REGION_SPLICE;
}

duckvep_region_state_t duckvep_region_classify_span(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          start1,
    uint32_t                          end1,
    uint32_t                          splice_exonic,
    uint32_t                          splice_intronic) {

    duckvep_region_state_t st;
    uint32_t ts = transcripts->start1[tx_idx];
    uint32_t te = transcripts->end1[tx_idx];
    int fwd = transcripts->strand[tx_idx] >= 0;
    uint32_t cds_s = transcripts->cds_start1[tx_idx];
    uint32_t cds_e = transcripts->cds_end1[tx_idx];
    size_t off = (size_t)transcripts->exon_offset[tx_idx];
    size_t cnt = (size_t)transcripts->exon_count[tx_idx];
    uint32_t min_exon_start = UINT32_MAX;
    uint32_t max_exon_end = 0u;
    int coarse_splice = 0;
    size_t e;

    memset(&st, 0, sizeof st);
    if (end1 < ts) {
        st.region_mask = fwd ? (uint32_t)DUCKVEP_REGION_UPSTREAM
                             : (uint32_t)DUCKVEP_REGION_DOWNSTREAM;
        return st;
    }
    if (start1 > te) {
        st.region_mask = fwd ? (uint32_t)DUCKVEP_REGION_DOWNSTREAM
                             : (uint32_t)DUCKVEP_REGION_UPSTREAM;
        return st;
    }

    st.within_feature = 1u;
    st.complete_overlap_feature = (uint8_t)(start1 <= ts && end1 >= te);
    st.complete_within_feature = (uint8_t)(start1 >= ts && end1 <= te);
    st.partial_overlap_feature =
        (uint8_t)(!st.complete_overlap_feature && !st.complete_within_feature);

    for (e = 0u; e < cnt; e++) {
        size_t ei = off + e;
        uint32_t es = exons->start1[ei];
        uint32_t ee = exons->end1[ei];

        min_exon_start = min_u32(min_exon_start, es);
        max_exon_end = max_u32(max_exon_end, ee);

        region_add_exon_overlap(&st, start1, end1, es, ee,
                                cds_s, cds_e, fwd);
        if ((splice_exonic | splice_intronic) != 0u &&
            overlaps_coarse_exon_reach(start1, end1, es, ee,
                                       splice_exonic, splice_intronic)) {
            coarse_splice = 1;
        }
    }

    for (e = 0u; e + 1u < cnt; e++) {
        uint32_t is;
        uint32_t ie;
        if (gap_between_exons(exons, off + e, off + e + 1u, &is, &ie) &&
            span_overlaps_u32(start1, end1, is, ie)) {
            if (gap_is_frameshift_intron(is, ie)) {
                region_add_frameshift_intron(
                    &st, start1, end1, cds_s, cds_e, fwd);
            } else {
                st.overlaps_intron = 1u;
            }
        }
    }
    /* Prepared transcript spans normally coincide with their outer exons, but
     * malformed/imported models need not. Preserve point-event exclusivity and
     * make span semantics total by treating any in-feature base outside the
     * exon envelope as non-exonic transcript sequence. Internal exon gaps were
     * detected above. */
    if (cnt == 0u ||
        (!st.overlaps_exon && !st.within_frameshift_intron) ||
        max_u32(start1, ts) < min_exon_start ||
        min_u32(end1, te) > max_exon_end) {
        st.overlaps_intron = 1u;
    }

    /* VEP's UTR predicates do not ask whether an annotated UTR base exists.
     * They call overlap() on [transcript_start, cds_start - 1] or
     * [cds_end + 1, transcript_end], then separately require a cDNA mapping.
     * When transcript and CDS share an endpoint that UTR interval is empty, but
     * VEP's four-comparison overlap still returns true for an event spanning the
     * inverted interval. Preserve that observable state: a deletion crossing a
     * CDS endpoint may carry both coding_sequence_variant and a UTR term even
     * when the transcript has no UTR bases on that side. */
    region_add_vep_endpoint_utr(&st, start1, end1, ts, te,
                                cds_s, cds_e, fwd);
    if (st.within_frameshift_intron && !st.overlaps_exon)
        coarse_splice = 0;
    region_finish_mask(&st, cds_s, coarse_splice);
    return st;
}

uint32_t duckvep_region_mask(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          pos,
    uint32_t                          splice_exonic,
    uint32_t                          splice_intronic) {
    return duckvep_region_classify_span(transcripts, exons, tx_idx, pos, pos,
                                        splice_exonic, splice_intronic).region_mask;
}

duckvep_splice_state_t duckvep_splice_classify_span_with_windows(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          start1,
    uint32_t                          end1,
    uint8_t                           interbase,
    uint32_t                          splice_exonic,
    uint32_t                          splice_intronic) {

    struct duckvep_splice_accum acc;
    int fwd = transcripts->strand[tx_idx] >= 0;
    size_t off = (size_t)transcripts->exon_offset[tx_idx];
    size_t cnt = (size_t)transcripts->exon_count[tx_idx];
    /* VEP models a pure insertion as vf->start = P+1, vf->end = P (start > end),
     * where P == end1 is the anchor base to its 5' side. Feeding (P+1, P) into the
     * same overlap predicate reproduces VEP's "both flanking bases inside the zone"
     * rule for insertions and leaves substitutions/deletions unchanged. Computed in
     * int64 so P+1 cannot wrap at the uint32 ceiling. irs/ire also drive the
     * exact-edge insertion special cases below. */
    int64_t irs = interbase ? (int64_t)end1 + 1 : (int64_t)start1;
    int64_t ire = (int64_t)end1;
    /* VEP computes the polypyrimidine tract in its intron-INTERIOR loop with the
     * SORTED (min,max) span — an insertion touches the tract when EITHER flank is
     * inside, not both (its other zones live in the boundary loop and use the
     * interbase rule). Use the ascending span for the two PPT windows only. */
    int64_t lo_s = interbase ? ire : (int64_t)start1;
    int64_t hi_s = interbase ? irs : (int64_t)end1;

    memset(&acc, 0, sizeof acc);

    for (size_t k = 0u; k + 1u < cnt; k++) {
        (void)splice_accum_add_gap(exons, off + k, off + k + 1u,
                                   irs, ire, lo_s, hi_s, interbase,
                                   1, 1, splice_exonic, splice_intronic,
                                   &acc, NULL, NULL);
    }
    return splice_accum_finish(&acc, fwd);
}

duckvep_splice_state_t duckvep_splice_classify_span(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          start1,
    uint32_t                          end1,
    uint8_t                           interbase) {

    return duckvep_splice_classify_span_with_windows(
        transcripts, exons, tx_idx, start1, end1, interbase,
        DUCKVEP_DEFAULT_SPLICE_REGION_EXONIC,
        DUCKVEP_DEFAULT_SPLICE_REGION_INTRONIC);
}

duckvep_splice_state_t duckvep_splice_classify_differing_regions_with_windows(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          feature_start1,
    const uint8_t                    *feature_ref,
    uint16_t                          feature_ref_length,
    const uint8_t                    *feature_alt,
    uint16_t                          feature_alt_length,
    uint32_t                          splice_exonic,
    uint32_t                          splice_intronic) {

    struct duckvep_splice_accum acc;
    size_t off;
    size_t cnt;
    int fwd;
    int64_t feature_end1;
    int64_t feature_min1;
    int64_t feature_max1;

    memset(&acc, 0, sizeof acc);
    if (transcripts == NULL || exons == NULL ||
        tx_idx >= transcripts->transcript_count ||
        (feature_ref_length > 0u && feature_ref == NULL) ||
        (feature_alt_length > 0u && feature_alt == NULL)) {
        return splice_accum_finish(&acc, 1);
    }

    off = (size_t)transcripts->exon_offset[tx_idx];
    cnt = (size_t)transcripts->exon_count[tx_idx];
    fwd = transcripts->strand[tx_idx] >= 0;
    feature_end1 = (int64_t)feature_start1 +
                   (int64_t)feature_ref_length - 1;
    feature_min1 = feature_end1 < (int64_t)feature_start1
        ? feature_end1 : (int64_t)feature_start1;
    feature_max1 = feature_end1 > (int64_t)feature_start1
        ? feature_end1 : (int64_t)feature_start1;

    if (feature_ref_length > 1u && feature_alt_length > 1u) {
        uint16_t length = feature_ref_length > feature_alt_length
            ? feature_ref_length : feature_alt_length;
        uint16_t i = 0u;

        /* Perl's string XOR pads the shorter operand with zero bytes. DNA
         * alleles cannot contain zero, so every overhanging byte is a mismatch. */
        while (i < length) {
            uint8_t ref_base = i < feature_ref_length ? feature_ref[i] : 0u;
            uint8_t alt_base = i < feature_alt_length ? feature_alt[i] : 0u;
            uint16_t first;
            uint16_t last;

            if (ref_base == alt_base) {
                i++;
                continue;
            }
            first = i;
            last = i;
            i++;
            while (i < length) {
                ref_base = i < feature_ref_length ? feature_ref[i] : 0u;
                alt_base = i < feature_alt_length ? feature_alt[i] : 0u;
                if (ref_base == alt_base) break;
                last = i;
                i++;
            }
            splice_accum_add_region(
                exons, off, cnt,
                (int64_t)feature_start1 + (int64_t)first,
                (int64_t)feature_start1 + (int64_t)last,
                feature_min1, feature_max1,
                splice_exonic, splice_intronic,
                &acc);
        }
    } else {
        /* This is VEP's fallback {s => 0, e => ref_length - 1}. A zero-length
         * reference therefore becomes the reversed insertion interval P+1,P. */
        splice_accum_add_region(
            exons, off, cnt, (int64_t)feature_start1,
            (int64_t)feature_start1 + (int64_t)feature_ref_length - 1,
            feature_min1, feature_max1,
            splice_exonic, splice_intronic,
            &acc);
    }
    return splice_accum_finish(&acc, fwd);
}

duckvep_splice_state_t duckvep_splice_classify_differing_regions(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          feature_start1,
    const uint8_t                    *feature_ref,
    uint16_t                          feature_ref_length,
    const uint8_t                    *feature_alt,
    uint16_t                          feature_alt_length) {

    return duckvep_splice_classify_differing_regions_with_windows(
        transcripts, exons, tx_idx, feature_start1,
        feature_ref, feature_ref_length, feature_alt, feature_alt_length,
        DUCKVEP_DEFAULT_SPLICE_REGION_EXONIC,
        DUCKVEP_DEFAULT_SPLICE_REGION_INTRONIC);
}

static size_t point_exon_index(const duckvep_transcript_model_t *transcripts,
                               size_t tx_idx, uint32_t rank) {
    size_t off = (size_t)transcripts->exon_offset[tx_idx];
    uint32_t cnt = (uint32_t)transcripts->exon_count[tx_idx];

    if (transcripts->strand[tx_idx] >= 0)
        return off + (size_t)rank;
    return off + (size_t)(cnt - rank - 1u);
}

static uint32_t point_exon_seek(const duckvep_transcript_model_t *transcripts,
                                const duckvep_exon_model_t *exons,
                                size_t tx_idx, uint32_t pos, uint32_t rank) {
    uint32_t cnt = (uint32_t)transcripts->exon_count[tx_idx];

    if (rank == UINT16_MAX || rank > cnt) {
        uint32_t lo = 0u;
        uint32_t hi = cnt;

        /* This search runs once when a transcript first enters a sorted run.
         * Later variants only advance the cursor. */
        while (lo < hi) {
            uint32_t mid = lo + (hi - lo) / 2u;
            size_t ei = point_exon_index(transcripts, tx_idx, mid);

            if (exons->end1[ei] < pos)
                lo = mid + 1u;
            else
                hi = mid;
        }
        rank = lo;
    } else {
        while (rank < cnt &&
               exons->end1[point_exon_index(transcripts, tx_idx, rank)] < pos)
            rank++;
    }
    return rank;
}

static uint32_t point_exon_seek_with_rewinds(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          pos,
    uint32_t                          rank) {

    uint32_t cnt = (uint32_t)transcripts->exon_count[tx_idx];

    /* A padded or long allele can move its semantic/feature coordinate behind
     * the preceding uploaded record. When the previous exon can again contain
     * or follow the coordinate, force the ordinary seek to restart with one
     * binary search. Monotone vectors call point_exon_seek() directly and pay
     * none of this repair check. */
    if (rank > 0u && rank <= cnt &&
        exons->end1[point_exon_index(
            transcripts, tx_idx, rank - 1u)] >= pos) {
        rank = UINT16_MAX;
    }
    return point_exon_seek(transcripts, exons, tx_idx, pos, rank);
}

int duckvep_point_outside_splice_reach(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          pos,
    uint32_t                          splice_exonic,
    uint32_t                          splice_intronic) {

    uint32_t cnt = (uint32_t)transcripts->exon_count[tx_idx];
    uint32_t first_gap_start;
    uint32_t last_gap_end;
    uint32_t ignored;
    uint32_t margin = splice_exonic > splice_intronic
        ? splice_exonic : splice_intronic;

    if (cnt < 2u) return 1;
    if (margin < DUCKVEP_VEP_SPLICE_FIXED_REACH)
        margin = DUCKVEP_VEP_SPLICE_FIXED_REACH;
    if (!gap_between_exons(
            exons,
            point_exon_index(transcripts, tx_idx, 0u),
            point_exon_index(transcripts, tx_idx, 1u),
            &first_gap_start, &ignored) ||
        !gap_between_exons(
            exons,
            point_exon_index(transcripts, tx_idx, cnt - 2u),
            point_exon_index(transcripts, tx_idx, cnt - 1u),
            &ignored, &last_gap_end)) {
        return 0;
    }
    return (int64_t)pos < (int64_t)first_gap_start - (int64_t)margin ||
           (int64_t)pos > (int64_t)last_gap_end + (int64_t)margin;
}

duckvep_region_state_t duckvep_region_classify_span_sorted(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          start1,
    uint32_t                          end1,
    uint32_t                          splice_exonic,
    uint32_t                          splice_intronic,
    uint8_t                           repair_rewinds,
    uint16_t                         *exon_rank_io) {

    duckvep_region_state_t st;
    uint32_t ts;
    uint32_t te;
    uint32_t cds_s;
    uint32_t cds_e;
    uint32_t cnt;
    uint32_t rank;
    int fwd;

    memset(&st, 0, sizeof st);
    if (transcripts == NULL || exons == NULL || exon_rank_io == NULL ||
        tx_idx >= transcripts->transcript_count) {
        return st;
    }

    cnt = (uint32_t)transcripts->exon_count[tx_idx];
    rank = repair_rewinds
        ? point_exon_seek_with_rewinds(
            transcripts, exons, tx_idx, start1, (uint32_t)*exon_rank_io)
        : point_exon_seek(
            transcripts, exons, tx_idx, start1, (uint32_t)*exon_rank_io);
    *exon_rank_io = (uint16_t)rank;
    ts = transcripts->start1[tx_idx];
    te = transcripts->end1[tx_idx];
    fwd = transcripts->strand[tx_idx] >= 0;

    if (end1 < ts) {
        st.region_mask = fwd ? (uint32_t)DUCKVEP_REGION_UPSTREAM
                             : (uint32_t)DUCKVEP_REGION_DOWNSTREAM;
        return st;
    }
    if (start1 > te) {
        st.region_mask = fwd ? (uint32_t)DUCKVEP_REGION_DOWNSTREAM
                             : (uint32_t)DUCKVEP_REGION_UPSTREAM;
        return st;
    }

    /* Coarse splice-mask callers and any span crossing a transcript boundary
     * retain the exhaustive definition. Production consequence annotation
     * passes zero reaches here because its exact splice facts are separate. */
    if (start1 > end1 || start1 < ts || end1 > te ||
        (splice_exonic | splice_intronic) != 0u) {
        return duckvep_region_classify_span(
            transcripts, exons, tx_idx, start1, end1,
            splice_exonic, splice_intronic);
    }

    st.within_feature = 1u;
    st.complete_overlap_feature = (uint8_t)(start1 <= ts && end1 >= te);
    st.complete_within_feature = 1u;
    cds_s = transcripts->cds_start1[tx_idx];
    cds_e = transcripts->cds_end1[tx_idx];

    if (rank < cnt) {
        size_t ei = point_exon_index(transcripts, tx_idx, rank);
        uint32_t es = exons->start1[ei];
        uint32_t ee = exons->end1[ei];

        if (start1 >= es && end1 <= ee) {
            region_add_exon_overlap(&st, start1, end1, es, ee,
                                    cds_s, cds_e, fwd);
        } else if (end1 < es) {
            /* The cursor identifies the next exon. The preceding exon ended
             * before start1, so the complete span lies in one gap. */
            uint32_t gap_start;
            uint32_t gap_end;

            if (rank > 0u &&
                gap_between_exons(
                    exons,
                    point_exon_index(transcripts, tx_idx, rank - 1u),
                    point_exon_index(transcripts, tx_idx, rank),
                    &gap_start, &gap_end) &&
                gap_is_frameshift_intron(gap_start, gap_end)) {
                region_add_frameshift_intron(
                    &st, start1, end1, cds_s, cds_e, fwd);
            } else {
                st.overlaps_intron = 1u;
            }
        } else {
            return duckvep_region_classify_span(
                transcripts, exons, tx_idx, start1, end1, 0u, 0u);
        }
    } else {
        /* Valid prepared models normally end with the outer exon. Preserve the
         * exhaustive classifier's total semantics for any in-feature tail. */
        st.overlaps_intron = 1u;
    }

    region_add_vep_endpoint_utr(&st, start1, end1, ts, te,
                                cds_s, cds_e, fwd);
    region_finish_mask(&st, cds_s, 0);
    return st;
}

duckvep_splice_state_t
duckvep_splice_classify_differing_regions_sorted_with_windows(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          feature_start1,
    const uint8_t                    *feature_ref,
    uint16_t                          feature_ref_length,
    const uint8_t                    *feature_alt,
    uint16_t                          feature_alt_length,
    uint32_t                          splice_exonic,
    uint32_t                          splice_intronic,
    uint16_t                          exon_rank) {

    duckvep_splice_state_t empty;
    uint32_t cnt;
    uint32_t mismatch_span_length;
    uint32_t empty_margin;

    memset(&empty, 0, sizeof empty);
    if (transcripts == NULL || exons == NULL ||
        tx_idx >= transcripts->transcript_count ||
        (feature_ref_length > 0u && feature_ref == NULL) ||
        (feature_alt_length > 0u && feature_alt == NULL)) {
        return empty;
    }

    cnt = (uint32_t)transcripts->exon_count[tx_idx];
    if (cnt < 2u || exon_rank == UINT16_MAX || (uint32_t)exon_rank > cnt)
        return cnt < 2u ? empty
                        : duckvep_splice_classify_differing_regions_with_windows(
                            transcripts, exons, tx_idx, feature_start1,
                            feature_ref, feature_ref_length,
                            feature_alt, feature_alt_length,
                            splice_exonic, splice_intronic);

    /* A 5 kb consequence halo creates many upstream/downstream transcript
     * candidates. None can have a splice consequence when even the largest
     * raw/mismatch span lies beyond the outermost intron plus the exon-side
     * cache reach. Prove that from the first/last genomic gap before entering
     * the exhaustive mismatch-island loop. The ALT length is included because
     * Perl's padded XOR can expose an ALT-only island past the REF feature. */
    {
        uint32_t first_gap_start;
        uint32_t last_gap_end;
        uint32_t ignored_gap_endpoint;
        uint32_t outer_margin = splice_exonic > splice_intronic
            ? splice_exonic : splice_intronic;
        uint32_t span_length = feature_ref_length > feature_alt_length
            ? (uint32_t)feature_ref_length : (uint32_t)feature_alt_length;
        int64_t query_lo = (int64_t)feature_start1;
        int64_t query_hi = span_length == 0u
            ? query_lo - 1 : query_lo + (int64_t)span_length - 1;

        /* VEP leaves its short-intron PPT and generic splice windows
         * unclamped, so an intron-side window may extend into an outer exon. */
        if (outer_margin < DUCKVEP_VEP_SPLICE_FIXED_REACH)
            outer_margin = DUCKVEP_VEP_SPLICE_FIXED_REACH;
        if (feature_ref_length == 0u) query_lo--;
        if (gap_between_exons(
                exons,
                point_exon_index(transcripts, tx_idx, 0u),
                point_exon_index(transcripts, tx_idx, 1u),
                &first_gap_start, &ignored_gap_endpoint) &&
            gap_between_exons(
                exons,
                point_exon_index(transcripts, tx_idx, cnt - 2u),
                point_exon_index(transcripts, tx_idx, cnt - 1u),
                &ignored_gap_endpoint, &last_gap_end) &&
            (query_hi < (int64_t)first_gap_start - (int64_t)outer_margin ||
             query_lo > (int64_t)last_gap_end + (int64_t)outer_margin)) {
            return empty;
        }
    }

    /* Most coding edits are far from a junction. The 16-base floor covers VEP's
     * fixed PPT window; the configurable reaches cover its deliberately
     * unclamped short-intron windows. */
    mismatch_span_length = feature_ref_length > feature_alt_length
        ? (uint32_t)feature_ref_length : (uint32_t)feature_alt_length;
    empty_margin = splice_exonic > splice_intronic
        ? splice_exonic : splice_intronic;
    if (empty_margin < DUCKVEP_VEP_SPLICE_FIXED_REACH)
        empty_margin = DUCKVEP_VEP_SPLICE_FIXED_REACH;
    if ((uint32_t)exon_rank < cnt && mismatch_span_length != 0u) {
        size_t ei = point_exon_index(
            transcripts, tx_idx, (uint32_t)exon_rank);
        uint64_t mismatch_end1 = (uint64_t)feature_start1 +
                                 (uint64_t)mismatch_span_length - 1u;

        /* The fixed PPT window is inclusive.  In a short intron it can extend
         * exactly 16 bases into the neighboring exon, and an insertion's
         * reversed interval also reaches the base immediately before
         * feature_start1.  Only skip the exhaustive predicates when both
         * sides are strictly beyond that inclusive reach. */
        if ((uint64_t)feature_start1 >
                (uint64_t)exons->start1[ei] + empty_margin &&
            mismatch_end1 + empty_margin < (uint64_t)exons->end1[ei]) {
            return empty;
        }
    }

    return duckvep_splice_classify_differing_regions_with_windows(
        transcripts, exons, tx_idx, feature_start1,
        feature_ref, feature_ref_length,
        feature_alt, feature_alt_length,
        splice_exonic, splice_intronic);
}

void duckvep_classify_point_sorted(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          pos,
    uint32_t                          splice_exonic,
    uint32_t                          splice_intronic,
    uint8_t                           repair_rewinds,
    uint16_t                         *exon_rank_io,
    duckvep_region_state_t           *region_out,
    duckvep_splice_state_t           *splice_out) {

    duckvep_region_state_t region;
    duckvep_splice_state_t splice;
    struct duckvep_splice_accum acc;
    uint32_t ts;
    uint32_t te;
    uint32_t cds_s;
    uint32_t cds_e;
    uint32_t cnt;
    uint32_t rank;
    uint32_t current_exon_start = 0u;
    uint32_t current_exon_end = 0u;
    int fwd;
    int in_exon = 0;

    memset(&region, 0, sizeof region);
    memset(&splice, 0, sizeof splice);
    memset(&acc, 0, sizeof acc);
    if (region_out == NULL || splice_out == NULL || exon_rank_io == NULL ||
        transcripts == NULL || exons == NULL ||
        tx_idx >= transcripts->transcript_count) {
        if (region_out != NULL) *region_out = region;
        if (splice_out != NULL) *splice_out = splice;
        return;
    }

    ts = transcripts->start1[tx_idx];
    te = transcripts->end1[tx_idx];
    fwd = transcripts->strand[tx_idx] >= 0;
    if (pos < ts) {
        region.region_mask = fwd ? (uint32_t)DUCKVEP_REGION_UPSTREAM
                                 : (uint32_t)DUCKVEP_REGION_DOWNSTREAM;
        if (duckvep_point_outside_splice_reach(
                transcripts, exons, tx_idx, pos,
                splice_exonic, splice_intronic)) {
            *region_out = region;
            *splice_out = splice;
            return;
        }
        cnt = (uint32_t)transcripts->exon_count[tx_idx];
        rank = repair_rewinds
            ? point_exon_seek_with_rewinds(
                transcripts, exons, tx_idx, pos, *exon_rank_io)
            : point_exon_seek(
                transcripts, exons, tx_idx, pos, *exon_rank_io);
        *exon_rank_io = (uint16_t)rank;
        goto exact_splice;
    }
    if (pos > te) {
        region.region_mask = fwd ? (uint32_t)DUCKVEP_REGION_DOWNSTREAM
                                 : (uint32_t)DUCKVEP_REGION_UPSTREAM;
        if (duckvep_point_outside_splice_reach(
                transcripts, exons, tx_idx, pos,
                splice_exonic, splice_intronic)) {
            *region_out = region;
            *splice_out = splice;
            return;
        }
        cnt = (uint32_t)transcripts->exon_count[tx_idx];
        rank = repair_rewinds
            ? point_exon_seek_with_rewinds(
                transcripts, exons, tx_idx, pos, *exon_rank_io)
            : point_exon_seek(
                transcripts, exons, tx_idx, pos, *exon_rank_io);
        *exon_rank_io = (uint16_t)rank;
        goto exact_splice;
    }

    region.within_feature = 1u;
    region.complete_overlap_feature = (uint8_t)(ts == te && pos == ts);
    region.complete_within_feature = 1u;
    cnt = (uint32_t)transcripts->exon_count[tx_idx];
    rank = repair_rewinds
        ? point_exon_seek_with_rewinds(
            transcripts, exons, tx_idx, pos, *exon_rank_io)
        : point_exon_seek(
            transcripts, exons, tx_idx, pos, *exon_rank_io);
    *exon_rank_io = (uint16_t)rank;

    if (rank < cnt) {
        size_t ei = point_exon_index(transcripts, tx_idx, rank);
        uint32_t es = exons->start1[ei];

        if (pos >= es) {
            in_exon = 1;
            current_exon_start = es;
            current_exon_end = exons->end1[ei];
            region.within_cdna = 1u;
            region.overlaps_exon = 1u;
            cds_s = transcripts->cds_start1[tx_idx];
            cds_e = transcripts->cds_end1[tx_idx];
            if (cds_s != 0u) {
                if (pos >= cds_s && pos <= cds_e) {
                    region.overlaps_cds = 1u;
                } else if (pos < cds_s) {
                    if (fwd) region.overlaps_utr5 = 1u;
                    else     region.overlaps_utr3 = 1u;
                } else {
                    if (fwd) region.overlaps_utr3 = 1u;
                    else     region.overlaps_utr5 = 1u;
                }
            }
        }
    }
    if (!in_exon) {
        uint32_t gap_start;
        uint32_t gap_end;

        if (rank > 0u && rank < cnt &&
            gap_between_exons(
                exons,
                point_exon_index(transcripts, tx_idx, rank - 1u),
                point_exon_index(transcripts, tx_idx, rank),
                &gap_start, &gap_end) &&
            gap_is_frameshift_intron(gap_start, gap_end)) {
            cds_s = transcripts->cds_start1[tx_idx];
            cds_e = transcripts->cds_end1[tx_idx];
            region_add_frameshift_intron(
                &region, pos, pos, cds_s, cds_e, fwd);
        } else {
            region.overlaps_intron = 1u;
        }
    }

    if (transcripts->cds_start1[tx_idx] == 0u && region.overlaps_exon)
        region.region_mask |= (uint32_t)DUCKVEP_REGION_EXON;
    if (region.overlaps_cds)
        region.region_mask |= (uint32_t)DUCKVEP_REGION_CDS;
    if (region.overlaps_utr5 || region.overlaps_utr3)
        region.region_mask |= (uint32_t)DUCKVEP_REGION_UTR;
    if (region.overlaps_intron)
        region.region_mask |= (uint32_t)DUCKVEP_REGION_INTRON;

    if (in_exon) {
        uint32_t margin = splice_exonic > splice_intronic
            ? splice_exonic : splice_intronic;

        /* VEP's fixed polypyrimidine windows reach 16 intronic bases and may
         * extend into a short neighboring exon. Configured generic splice
         * windows can be wider. A point strictly farther than the larger
         * reach from both ends of its current exon cannot satisfy any splice
         * predicate, so avoid visiting either neighboring intron. */
        if (margin < DUCKVEP_VEP_SPLICE_FIXED_REACH)
            margin = DUCKVEP_VEP_SPLICE_FIXED_REACH;
        if ((int64_t)pos >
                (int64_t)current_exon_start + (int64_t)margin &&
            (int64_t)pos <
                (int64_t)current_exon_end - (int64_t)margin) {
            *region_out = region;
            *splice_out = splice;
            return;
        }
    }

    /* The polypyrimidine predicate reaches 16 intronic bases; a caller may
     * request a wider generic splice-region window. Return the already-complete
     * deep-intron state only beyond both reaches. */
    if (!in_exon && rank > 0u && rank < cnt &&
        splice_exonic <= DUCKVEP_DEFAULT_SPLICE_REGION_EXONIC &&
        splice_intronic <= DUCKVEP_DEFAULT_SPLICE_REGION_INTRONIC) {
        uint32_t gap_start;
        uint32_t gap_end;
        int64_t intronic_reach;

        intronic_reach = splice_intronic > 1u
            ? (int64_t)splice_intronic - 1 : 0;
        if (intronic_reach < 16) intronic_reach = 16;

        if (gap_between_exons(exons,
            point_exon_index(transcripts, tx_idx, rank - 1u),
            point_exon_index(transcripts, tx_idx, rank),
            &gap_start, &gap_end) &&
            (int64_t)pos > (int64_t)gap_start + intronic_reach &&
            (int64_t)pos < (int64_t)gap_end - intronic_reach) {
            splice.intronic = 1u;
            *region_out = region;
            *splice_out = splice;
            return;
        }
    }

exact_splice:
    if (cnt >= 2u) {
        uint32_t pivot = rank == 0u ? 0u : rank - 1u;
        uint32_t g;
        int64_t p = (int64_t)pos;

        if (pivot >= cnt - 1u) pivot = cnt - 2u;
        g = pivot;
        for (;;) {
            int64_t reach_hi;
            int found = splice_accum_add_gap(exons,
                point_exon_index(transcripts, tx_idx, g),
                point_exon_index(transcripts, tx_idx, g + 1u),
                p, p, p, p, 0u, 1, 1,
                splice_exonic, splice_intronic,
                &acc, NULL, &reach_hi);
            if ((found && p > reach_hi) || g == 0u) break;
            g--;
        }
        for (g = pivot + 1u; g + 1u < cnt; g++) {
            int64_t reach_lo;
            int found = splice_accum_add_gap(exons,
                point_exon_index(transcripts, tx_idx, g),
                point_exon_index(transcripts, tx_idx, g + 1u),
                p, p, p, p, 0u, 1, 1,
                splice_exonic, splice_intronic,
                &acc, &reach_lo, NULL);
            if (found && p < reach_lo) break;
        }
    }
    splice = splice_accum_finish(&acc, fwd);

    *region_out = region;
    *splice_out = splice;
}

duckvep_splice_state_t duckvep_splice_classify(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          pos) {
    return duckvep_splice_classify_span(transcripts, exons, tx_idx, pos, pos, 0u);
}
