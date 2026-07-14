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

struct duckvep_splice_accum {
    int start_ss;
    int end_ss;
    int fifth;
    int fifth_rev;
    int dreg;
    int dreg_rev;
    int ppt;
    int ppt_rev;
    int intronic;
    int region;
};

static int splice_region_overlaps_gap(int64_t start1, int64_t end1,
                                      int64_t intron_start,
                                      int64_t intron_end,
                                      uint8_t interbase) {
    return span_overlaps_i64(start1, end1,
                             intron_start + 2, intron_start + 7) ||
           span_overlaps_i64(start1, end1,
                             intron_end - 7, intron_end - 2) ||
           span_overlaps_i64(start1, end1,
                             intron_start - 3, intron_start - 1) ||
           span_overlaps_i64(start1, end1,
                             intron_end + 1, intron_end + 3) ||
           (interbase && (start1 == intron_start || end1 == intron_end ||
                          start1 == intron_start + 2 ||
                          end1 == intron_end - 2));
}

/* Accumulate VEP's splice predicates for one intron. Span and sorted-point
 * traversal deliberately share this code: the optimization chooses fewer gaps,
 * never a second definition of the biology. */
static int splice_accum_add_gap(const duckvep_exon_model_t *exons,
                                size_t left_idx, size_t right_idx,
                                int64_t irs, int64_t ire,
                                int64_t lo_s, int64_t hi_s,
                                uint8_t interbase,
                                int intron_cached,
                                int boundary_cached,
                                struct duckvep_splice_accum *a,
                                int64_t *reach_lo, int64_t *reach_hi) {
    uint32_t gap_start;
    uint32_t gap_end;
    int64_t is;
    int64_t ie;

    if (!gap_between_exons(exons, left_idx, right_idx,
                           &gap_start, &gap_end))
        return 0;
    is = (int64_t)gap_start;
    ie = (int64_t)gap_end;
    if (reach_lo != NULL)
        *reach_lo = is - 3 < ie - 16 ? is - 3 : ie - 16;
    if (reach_hi != NULL)
        *reach_hi = is + 16 > ie + 3 ? is + 16 : ie + 3;

    if (boundary_cached) {
        if (span_overlaps_i64(irs, ire, is, is + 1)) a->start_ss = 1;
        if (span_overlaps_i64(irs, ire, ie - 1, ie)) a->end_ss = 1;
        if (span_overlaps_i64(irs, ire, is + 4, is + 4)) a->fifth = 1;
        if (span_overlaps_i64(irs, ire, ie - 4, ie - 4)) a->fifth_rev = 1;
        if (span_overlaps_i64(irs, ire, is + 2, is + 5)) a->dreg = 1;
        if (span_overlaps_i64(irs, ire, ie - 5, ie - 2)) a->dreg_rev = 1;
        if (splice_region_overlaps_gap(irs, ire, is, ie, interbase))
            a->region = 1;
    }
    if (intron_cached || boundary_cached) {
        if (span_overlaps_i64(lo_s, hi_s, ie - 16, ie - 2)) a->ppt = 1;
        if (span_overlaps_i64(lo_s, hi_s, is + 2, is + 16)) a->ppt_rev = 1;
    }
    if (intron_cached &&
        (span_overlaps_i64(irs, ire, is + 2, ie - 2) ||
         (interbase && (irs == is + 2 || ire == ie - 2)))) {
        a->intronic = 1;
    }
    return 1;
}

static duckvep_splice_state_t splice_accum_finish(
    const struct duckvep_splice_accum *a, int fwd) {
    duckvep_splice_state_t st;

    memset(&st, 0, sizeof st);
    st.splice_donor = (uint8_t)(fwd ? a->start_ss : a->end_ss);
    st.splice_acceptor = (uint8_t)(fwd ? a->end_ss : a->start_ss);
    st.splice_donor_5th = (uint8_t)(fwd ? a->fifth : a->fifth_rev);
    st.splice_donor_region =
        (uint8_t)((fwd ? a->dreg : a->dreg_rev) && !st.splice_donor_5th);
    st.splice_polypyrimidine = (uint8_t)(fwd ? a->ppt : a->ppt_rev);
    st.intronic = (uint8_t)a->intronic;
    st.splice_region =
        (uint8_t)(a->region && !st.splice_donor && !st.splice_acceptor &&
                  !st.splice_donor_region && !st.splice_donor_5th);
    st.any = (uint8_t)(st.splice_donor || st.splice_acceptor ||
                       st.splice_donor_5th || st.splice_donor_region ||
                       st.splice_polypyrimidine || st.splice_region);
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
    struct duckvep_splice_accum *acc) {

    uint8_t interbase = (uint8_t)(region_start1 == region_end1 + 1);
    int64_t lo = region_start1 < region_end1 ? region_start1 : region_end1;
    int64_t hi = region_start1 > region_end1 ? region_start1 : region_end1;
    size_t k;

    for (k = 0u; k + 1u < exon_count; k++) {
        uint32_t gap_start;
        uint32_t gap_end;
        int intron_cached;
        int boundary_cached;

        if (!gap_between_exons(exons, exon_offset + k, exon_offset + k + 1u,
                               &gap_start, &gap_end)) {
            continue;
        }
        boundary_cached =
            span_overlaps_i64(feature_min1, feature_max1,
                              (int64_t)gap_start - 3,
                              (int64_t)gap_start + 7) ||
            span_overlaps_i64(feature_min1, feature_max1,
                              (int64_t)gap_end - 7,
                              (int64_t)gap_end + 3);
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
        (void)splice_accum_add_gap(
            exons, exon_offset + k, exon_offset + k + 1u,
            region_start1, region_end1, lo, hi, interbase,
            intron_cached, boundary_cached,
            acc, NULL, NULL);
        if (boundary_cached) {
            /* VEP 116 assigns, rather than ORs, splice_region for each
             * differing island/intron pair. Later islands can therefore clear
             * an earlier hit; the other splice facts remain accumulators. */
            acc->region = splice_region_overlaps_gap(
                region_start1, region_end1,
                (int64_t)gap_start, (int64_t)gap_end, interbase);
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
                                   1, 1,
                                   &acc, NULL, NULL);
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
                &acc);
        }
    } else {
        /* This is VEP's fallback {s => 0, e => ref_length - 1}. A zero-length
         * reference therefore becomes the reversed insertion interval P+1,P. */
        splice_accum_add_region(
            exons, off, cnt, (int64_t)feature_start1,
            (int64_t)feature_start1 + (int64_t)feature_ref_length - 1,
            feature_min1, feature_max1,
            &acc);
    }
    return splice_accum_finish(&acc, fwd);
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

void duckvep_classify_point_sorted(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          pos,
    uint32_t                          splice_exonic,
    uint32_t                          splice_intronic,
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
        cnt = (uint32_t)transcripts->exon_count[tx_idx];
        rank = point_exon_seek(transcripts, exons, tx_idx, pos, *exon_rank_io);
        *exon_rank_io = (uint16_t)rank;
        goto exact_splice;
    }
    if (pos > te) {
        region.region_mask = fwd ? (uint32_t)DUCKVEP_REGION_DOWNSTREAM
                                 : (uint32_t)DUCKVEP_REGION_UPSTREAM;
        cnt = (uint32_t)transcripts->exon_count[tx_idx];
        rank = point_exon_seek(transcripts, exons, tx_idx, pos, *exon_rank_io);
        *exon_rank_io = (uint16_t)rank;
        goto exact_splice;
    }

    region.within_feature = 1u;
    region.complete_overlap_feature = (uint8_t)(ts == te && pos == ts);
    region.complete_within_feature = 1u;
    cnt = (uint32_t)transcripts->exon_count[tx_idx];
    rank = point_exon_seek(transcripts, exons, tx_idx, pos, *exon_rank_io);
    *exon_rank_io = (uint16_t)rank;

    if (rank < cnt) {
        size_t ei = point_exon_index(transcripts, tx_idx, rank);
        uint32_t es = exons->start1[ei];

        if (pos >= es) {
            in_exon = 1;
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
    if (!in_exon)
        region.overlaps_intron = 1u;

    if ((splice_exonic | splice_intronic) != 0u) {
        if (in_exon) {
            size_t ei = point_exon_index(transcripts, tx_idx, rank);
            if (overlaps_coarse_exon_reach(pos, pos,
                                           exons->start1[ei], exons->end1[ei],
                                           splice_exonic, splice_intronic))
                region.region_mask |= (uint32_t)DUCKVEP_REGION_SPLICE;
            if (rank > 0u) {
                size_t pi = point_exon_index(transcripts, tx_idx, rank - 1u);
                if (overlaps_coarse_exon_reach(pos, pos,
                                               exons->start1[pi], exons->end1[pi],
                                               splice_exonic, splice_intronic))
                    region.region_mask |= (uint32_t)DUCKVEP_REGION_SPLICE;
            }
            if (rank + 1u < cnt) {
                size_t ni = point_exon_index(transcripts, tx_idx, rank + 1u);
                if (overlaps_coarse_exon_reach(pos, pos,
                                               exons->start1[ni], exons->end1[ni],
                                               splice_exonic, splice_intronic))
                    region.region_mask |= (uint32_t)DUCKVEP_REGION_SPLICE;
            }
        } else {
            if (rank < cnt) {
                size_t ni = point_exon_index(transcripts, tx_idx, rank);
                if (overlaps_coarse_exon_reach(pos, pos,
                                               exons->start1[ni], exons->end1[ni],
                                               splice_exonic, splice_intronic))
                    region.region_mask |= (uint32_t)DUCKVEP_REGION_SPLICE;
            }
            if (rank > 0u) {
                size_t pi = point_exon_index(transcripts, tx_idx, rank - 1u);
                if (overlaps_coarse_exon_reach(pos, pos,
                                               exons->start1[pi], exons->end1[pi],
                                               splice_exonic, splice_intronic))
                    region.region_mask |= (uint32_t)DUCKVEP_REGION_SPLICE;
            }
        }
    }

    if (transcripts->cds_start1[tx_idx] == 0u && in_exon)
        region.region_mask |= (uint32_t)DUCKVEP_REGION_EXON;
    if (region.overlaps_cds)
        region.region_mask |= (uint32_t)DUCKVEP_REGION_CDS;
    if (region.overlaps_utr5 || region.overlaps_utr3)
        region.region_mask |= (uint32_t)DUCKVEP_REGION_UTR;
    if (region.overlaps_intron)
        region.region_mask |= (uint32_t)DUCKVEP_REGION_INTRON;

    /* Exact VEP splice predicates reach at most 16 intronic bases from an
     * internal boundary. Return the already-complete deep-intron state
     * without running every splice predicate. */
    if (!in_exon && rank > 0u && rank < cnt) {
        uint32_t gap_start;
        uint32_t gap_end;

        if (gap_between_exons(exons,
            point_exon_index(transcripts, tx_idx, rank - 1u),
            point_exon_index(transcripts, tx_idx, rank),
            &gap_start, &gap_end) &&
            (int64_t)pos > (int64_t)gap_start + 16 &&
            (int64_t)pos < (int64_t)gap_end - 16) {
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
                p, p, p, p, 0u, 1, 1, &acc, NULL, &reach_hi);
            if ((found && p > reach_hi) || g == 0u) break;
            g--;
        }
        for (g = pivot + 1u; g + 1u < cnt; g++) {
            int64_t reach_lo;
            int found = splice_accum_add_gap(exons,
                point_exon_index(transcripts, tx_idx, g),
                point_exon_index(transcripts, tx_idx, g + 1u),
                p, p, p, p, 0u, 1, 1, &acc, &reach_lo, NULL);
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
