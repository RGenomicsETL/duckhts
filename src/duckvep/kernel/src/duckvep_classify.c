/*
 * duckvep_classify.c — span-aware structural and splice classification.
 * No DuckDB/htslib/Arrow; no allocation. Pure per-pair arithmetic.
 */
#include "duckvep_classify.h"

#include <string.h>

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

        if (span_overlaps_u32(start1, end1, es, ee)) {
            st.within_cdna = 1u;
            st.overlaps_exon = 1u;
            if (cds_s == 0u) {
                /* Non-coding exon; no CDS/UTR split. */
            } else {
                uint32_t ov_s = max_u32(start1, es);
                uint32_t ov_e = min_u32(end1, ee);
                if (span_overlaps_u32(ov_s, ov_e, cds_s, cds_e)) {
                    st.overlaps_cds = 1u;
                }
                if (ov_s < cds_s) {
                    uint32_t utr_e = min_u32(ov_e, cds_s - 1u);
                    if (ov_s <= utr_e) {
                        if (fwd) st.overlaps_utr5 = 1u;
                        else     st.overlaps_utr3 = 1u;
                    }
                }
                if (ov_e > cds_e && cds_e != UINT32_MAX) {
                    uint32_t utr_s = max_u32(ov_s, cds_e + 1u);
                    if (utr_s <= ov_e) {
                        if (fwd) st.overlaps_utr3 = 1u;
                        else     st.overlaps_utr5 = 1u;
                    }
                }
            }
        }
        if (overlaps_coarse_exon_reach(start1, end1, es, ee,
                                        splice_exonic, splice_intronic)) {
            coarse_splice = 1;
        }
    }

    for (e = 0u; e + 1u < cnt; e++) {
        uint32_t is;
        uint32_t ie;
        if (gap_between_exons(exons, off + e, off + e + 1u, &is, &ie) &&
            span_overlaps_u32(start1, end1, is, ie)) {
            st.overlaps_intron = 1u;
        }
    }
    /* Prepared transcript spans normally coincide with their outer exons, but
     * malformed/imported models need not. Preserve point-event exclusivity and
     * make span semantics total by treating any in-feature base outside the
     * exon envelope as non-exonic transcript sequence. Internal exon gaps were
     * detected above. */
    if (cnt == 0u || !st.overlaps_exon ||
        max_u32(start1, ts) < min_exon_start ||
        min_u32(end1, te) > max_exon_end) {
        st.overlaps_intron = 1u;
    }

    if (cds_s == 0u && st.overlaps_exon) {
        st.region_mask |= (uint32_t)DUCKVEP_REGION_EXON;
    }
    if (st.overlaps_cds) {
        st.region_mask |= (uint32_t)DUCKVEP_REGION_CDS;
    }
    if (st.overlaps_utr5 || st.overlaps_utr3) {
        st.region_mask |= (uint32_t)DUCKVEP_REGION_UTR;
    }
    if (st.overlaps_intron) {
        st.region_mask |= (uint32_t)DUCKVEP_REGION_INTRON;
    }
    if (coarse_splice) {
        st.region_mask |= (uint32_t)DUCKVEP_REGION_SPLICE;
    }
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

duckvep_splice_state_t duckvep_splice_classify_span(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          start1,
    uint32_t                          end1,
    uint8_t                           interbase) {

    duckvep_splice_state_t st;
    int fwd = transcripts->strand[tx_idx] >= 0;
    size_t off = (size_t)transcripts->exon_offset[tx_idx];
    size_t cnt = (size_t)transcripts->exon_count[tx_idx];
    int start_ss = 0, end_ss = 0;
    int fifth = 0, fifth_rev = 0;
    int dreg = 0, dreg_rev = 0;
    int ppt = 0, ppt_rev = 0;
    int intronic = 0;
    int region = 0;
    size_t k;
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

    memset(&st, 0, sizeof st);
    if (cnt < 2u) return st;

    for (k = 0u; k + 1u < cnt; k++) {
        uint32_t gap_start;
        uint32_t gap_end;
        int64_t is;
        int64_t ie;
        if (!gap_between_exons(exons, off + k, off + k + 1u,
                               &gap_start, &gap_end)) {
            continue;
        }
        is = (int64_t)gap_start;
        ie = (int64_t)gap_end;

        if (span_overlaps_i64(irs, ire, is, is + 1)) start_ss = 1;
        if (span_overlaps_i64(irs, ire, ie - 1, ie)) end_ss = 1;
        if (span_overlaps_i64(irs, ire, is + 4, is + 4)) fifth = 1;
        if (span_overlaps_i64(irs, ire, ie - 4, ie - 4)) fifth_rev = 1;
        if (span_overlaps_i64(irs, ire, is + 2, is + 5)) dreg = 1;
        if (span_overlaps_i64(irs, ire, ie - 5, ie - 2)) dreg_rev = 1;
        if (span_overlaps_i64(lo_s, hi_s, ie - 16, ie - 2)) ppt = 1;
        if (span_overlaps_i64(lo_s, hi_s, is + 2, is + 16)) ppt_rev = 1;
        if (span_overlaps_i64(irs, ire, is + 2, ie - 2) ||
            (interbase && (irs == is + 2 || ire == ie - 2))) {
            intronic = 1;
        }
        /* VEP _intron_overlap: intron interior 3-8, exon rim 1-3, plus insertion
         * exact-edge cases at the exon/intron and donor/acceptor boundaries. */
        if (span_overlaps_i64(irs, ire, is + 2, is + 7) ||
            span_overlaps_i64(irs, ire, ie - 7, ie - 2) ||
            span_overlaps_i64(irs, ire, is - 3, is - 1) ||
            span_overlaps_i64(irs, ire, ie + 1, ie + 3) ||
            (interbase && (irs == is || ire == ie ||
                           irs == is + 2 || ire == ie - 2))) {
            region = 1;
        }
    }

    st.splice_donor = (uint8_t)(fwd ? start_ss : end_ss);
    st.splice_acceptor = (uint8_t)(fwd ? end_ss : start_ss);
    st.splice_donor_5th = (uint8_t)(fwd ? fifth : fifth_rev);
    st.splice_donor_region =
        (uint8_t)((fwd ? dreg : dreg_rev) && !st.splice_donor_5th);
    st.splice_polypyrimidine = (uint8_t)(fwd ? ppt : ppt_rev);
    st.intronic = (uint8_t)intronic;
    st.splice_region =
        (uint8_t)(region && !st.splice_donor && !st.splice_acceptor &&
                  !st.splice_donor_region && !st.splice_donor_5th);
    st.any = (uint8_t)(st.splice_donor || st.splice_acceptor ||
                       st.splice_donor_5th || st.splice_donor_region ||
                       st.splice_polypyrimidine || st.splice_region);
    return st;
}

duckvep_splice_state_t duckvep_splice_classify(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          pos) {
    return duckvep_splice_classify_span(transcripts, exons, tx_idx, pos, pos, 0u);
}
