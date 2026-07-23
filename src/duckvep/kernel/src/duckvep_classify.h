/*
 * duckvep_classify.h — span-aware transcript topology facts (INTERNAL).
 *
 * One event may overlap several structural regions. The classifier therefore
 * returns an explicit fact record rather than forcing point-event mutual
 * exclusivity. `region_mask` remains the compact public diagnostic summary.
 */
#ifndef DUCKVEP_CLASSIFY_H
#define DUCKVEP_CLASSIFY_H

#include "duckvep_kernel.h"

#include <stddef.h>
#include <stdint.h>

/* Fixed VEP 116 predicate reaches. Caller-configurable generic splice-region
 * distances remain in duckvep_options_init_t. */
#define DUCKVEP_VEP_FRAMESHIFT_INTRON_MAX_SPAN 12u
#define DUCKVEP_VEP_SPLICE_FIXED_REACH          16u

typedef struct duckvep_region_state {
    uint32_t region_mask;
    uint8_t  within_feature;
    uint8_t  complete_overlap_feature; /* event contains the transcript span */
    uint8_t  complete_within_feature;  /* event is contained by transcript span */
    uint8_t  partial_overlap_feature;  /* crosses exactly one transcript boundary */
    uint8_t  within_cdna;              /* overlaps at least one exon */
    uint8_t  overlaps_exon;
    uint8_t  overlaps_intron;
    uint8_t  within_frameshift_intron; /* gap span <= fixed VEP threshold */
    uint8_t  overlaps_cds;
    uint8_t  overlaps_utr5;
    uint8_t  overlaps_utr3;
} duckvep_region_state_t;

/* Classify one 1-based inclusive event span against one transcript. The coarse
 * splice reaches are retained only for the diagnostic region mask; consequence
 * emission uses duckvep_splice_classify_span below. */
duckvep_region_state_t duckvep_region_classify_span(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          start1,
    uint32_t                          end1,
    uint32_t                          splice_exonic,
    uint32_t                          splice_intronic);

/* Sorted small-span fast path. `exon_rank_io` is the next genomic-order exon
 * whose end is not before start1. Spans contained in one exon or one intron are
 * classified in constant time; exon-intron-junction-crossing spans use the
 * exhaustive classifier above. Set `repair_rewinds` only when normalized
 * feature coordinates can move backwards within the current input vector. */
duckvep_region_state_t duckvep_region_classify_span_sorted(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          start1,
    uint32_t                          end1,
    uint32_t                          splice_exonic,
    uint32_t                          splice_intronic,
    uint8_t                           repair_rewinds,
    uint16_t                         *exon_rank_io);

/* Point compatibility wrapper. */
uint32_t duckvep_region_mask(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          pos,
    uint32_t                          splice_exonic,
    uint32_t                          splice_intronic);

/* Return nonzero only when a point is provably outside every splice predicate
 * window of the transcript. The 16-base floor retains VEP's fixed
 * polypyrimidine windows; caller-selected generic windows may widen it. */
int duckvep_point_outside_splice_reach(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          pos,
    uint32_t                          splice_exonic,
    uint32_t                          splice_intronic);

/* VEP-source-grounded per-intron splice predicate facts. A span can overlap more
 * than one window and therefore legitimately set several facts. Outer transcript
 * ends are excluded because only gaps between consecutive exons are inspected. */
typedef struct duckvep_splice_state {
    uint8_t splice_donor;
    uint8_t splice_acceptor;
    uint8_t splice_donor_5th;
    uint8_t splice_donor_region;
    uint8_t splice_polypyrimidine;
    uint8_t splice_region;
    uint8_t intronic;
    uint8_t within_frameshift_intron;
    uint8_t any;
} duckvep_splice_state_t;

/* Fused point-event classifier for the sorted hot path. `exon_rank_io` is the
 * genomic-order rank of the first exon whose end is not before `pos`; initialize
 * it to UINT16_MAX for a new transcript run. Successive calls for the same
 * transcript must have non-decreasing positions unless `repair_rewinds` is set.
 * In repair mode a backwards normalized coordinate performs one binary seek.
 * The model already limits exon counts to uint16_t. Exons remain stored in
 * transcript order, so genomic rank
 * runs forward on both strands while the implementation maps minus-strand ranks
 * back to the borrowed exon slice. `splice_exonic` and `splice_intronic`
 * configure the exact generic splice-region predicate; `region_out` omits the
 * older coarse proximity bit because `splice_out` is authoritative. */
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
    duckvep_splice_state_t           *splice_out);

/* `interbase` marks a pure insertion: it lands BETWEEN two reference bases rather
 * than replacing any. VEP models it as vf->start = end1+1, vf->end = end1 (start >
 * end) so a zone [lo,hi] is only touched when BOTH flanking bases fall inside it
 * (P in [lo, hi-1]); it also fires several exact-edge special cases at the intron
 * rim. For a substitution/deletion pass interbase = 0 and the plain span applies. */
duckvep_splice_state_t duckvep_splice_classify_span(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          start1,
    uint32_t                          end1,
    uint8_t                           interbase);

/* Annotation uses the resolved per-call windows. The unsuffixed helper above
 * remains the default-window oracle used by focused predicate tests. */
duckvep_splice_state_t duckvep_splice_classify_span_with_windows(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          start1,
    uint32_t                          end1,
    uint8_t                           interbase,
    uint32_t                          splice_exonic,
    uint32_t                          splice_intronic);

/* VEP does not feed one minimized interval to its splice predicates. For
 * feature alleles longer than one base it XORs REF and ALT, groups contiguous
 * non-zero bytes, and accumulates every resulting mismatch island before
 * applying donor/acceptor/region suppression. Missing bytes in the shorter
 * allele compare as zero, matching Perl's string XOR. `feature_start1` and the
 * allele slices are VEP's complete uploaded equal-length pair or its fully
 * minimized length-changing pair. */
duckvep_splice_state_t duckvep_splice_classify_differing_regions(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          feature_start1,
    const uint8_t                    *feature_ref,
    uint16_t                          feature_ref_length,
    const uint8_t                    *feature_alt,
    uint16_t                          feature_alt_length);

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
    uint32_t                          splice_intronic);

/* Return the same differing-region splice state while proving the common
 * deep-inside-one-exon case empty from the sorted exon cursor. */
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
    uint16_t                          exon_rank);

/* Point compatibility wrapper. */
duckvep_splice_state_t duckvep_splice_classify(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          pos);

#endif /* DUCKVEP_CLASSIFY_H */
