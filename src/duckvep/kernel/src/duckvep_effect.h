/*
 * duckvep_effect.h — VEP-shaped consequence decision program (INTERNAL).
 *
 * Cheap event/topology facts and cached predicate results gate a generated,
 * tier-ordered rule table. Rank remains severity metadata; tier alone controls
 * suppression. Sequence deltas are populated lazily for supported coding edits.
 */
#ifndef DUCKVEP_EFFECT_H
#define DUCKVEP_EFFECT_H

#include "duckvep_classify.h"
#include "duckvep_delta.h"
#include "duckvep_event.h"
#include "duckvep_kernel.h"
#include "duckvep_sv.h"

#include <stddef.h>
#include <stdint.h>

typedef enum duckvep_pre_bit {
    DUCKVEP_PRE_UPSTREAM             = 0,
    DUCKVEP_PRE_DOWNSTREAM           = 1,
    DUCKVEP_PRE_INTRON               = 2,
    DUCKVEP_PRE_EXON                 = 3,
    DUCKVEP_PRE_CDS                  = 4,
    DUCKVEP_PRE_UTR5                 = 5,
    DUCKVEP_PRE_UTR3                 = 6,
    DUCKVEP_PRE_CODING               = 7,
    DUCKVEP_PRE_NONCODING            = 8,
    DUCKVEP_PRE_SPLICE_DONOR         = 9,
    DUCKVEP_PRE_SPLICE_ACCEPTOR      = 10,
    DUCKVEP_PRE_SPLICE_DONOR_5TH     = 11,
    DUCKVEP_PRE_SPLICE_DONOR_REGION  = 12,
    DUCKVEP_PRE_SPLICE_PPT           = 13,
    DUCKVEP_PRE_SPLICE_REGION        = 14,
    DUCKVEP_PRE_DELTA                = 15,
    DUCKVEP_PRE_SYNONYMOUS           = 16,
    DUCKVEP_PRE_MISSENSE             = 17,
    DUCKVEP_PRE_STOP_GAINED          = 18,
    DUCKVEP_PRE_STOP_LOST            = 19,
    DUCKVEP_PRE_STOP_RETAINED        = 20,
    DUCKVEP_PRE_INSERTION            = 21,
    DUCKVEP_PRE_DELETION             = 22,
    DUCKVEP_PRE_SV                   = 23,
    DUCKVEP_PRE_WITHIN_INTRON        = 24,
    DUCKVEP_PRE_START_LOST           = 25,
    DUCKVEP_PRE_FRAMESHIFT           = 26,
    DUCKVEP_PRE_INFRAME_DELETION     = 27,
    /* Cached structural/CNV predicate truth. */
    DUCKVEP_PRE_FEATURE_ABLATION      = 28,
    DUCKVEP_PRE_FEATURE_AMPLIFICATION = 29,
    DUCKVEP_PRE_FEATURE_ELONGATION    = 30,
    DUCKVEP_PRE_FEATURE_TRUNCATION    = 31,
    /* VEP predicate truth that is not equivalent to raw placement. These
     * distinctions first matter for whole-transcript spans: complete overlap
     * suppresses non_coding_exon_variant/coding_unknown and instead enables the
     * transcript-level fallback predicates. */
    DUCKVEP_PRE_NONCODING_EXON         = 32,
    DUCKVEP_PRE_WITHIN_NONCODING_GENE  = 33,
    DUCKVEP_PRE_CODING_UNKNOWN         = 34,
    DUCKVEP_PRE_CODING_TRANSCRIPT      = 35,
    DUCKVEP_PRE_INFRAME_INSERTION      = 36,
    DUCKVEP_PRE_PROTEIN_ALTERING       = 37,
    DUCKVEP_PRE_WITHIN_NMD_TRANSCRIPT  = 38,
    DUCKVEP_PRE_START_RETAINED          = 39,
    DUCKVEP_PRE_BIT_COUNT               = 40
} duckvep_pre_bit_t;

#define DUCKVEP_PRE(b) (UINT64_C(1) << (b))

typedef struct duckvep_effect_ctx {
    uint32_t               variant_idx;
    uint32_t               tx_idx;
    uint32_t               start1;
    uint32_t               end1;
    int8_t                 strand;
    uint32_t               region;
    uint64_t               pre_bits;
    duckvep_region_state_t region_state;
    duckvep_splice_state_t splice;
} duckvep_effect_ctx_t;

/* `interbase` marks a pure insertion (it lands between two reference bases, so
 * event.start1 == event.end1 == the 5' anchor). It is forwarded to the splice
 * classifier so insertion splice-tier facts follow VEP's interbase overlap rule. */
void duckvep_effect_ctx_fill(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          start1,
    uint32_t                          end1,
    uint8_t                           interbase,
    uint32_t                          splice_region_exonic,
    uint32_t                          splice_region_intronic,
    duckvep_effect_ctx_t             *out);

/* Region placement and splice overlap are distinct for a pure insertion at an
 * exon/CDS boundary. The mapper may place the event on the right flank while
 * VEP's splice predicates retain the reversed interbase interval (P+1,P). */
void duckvep_effect_ctx_fill_geometry(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          region_start1,
    uint32_t                          region_end1,
    uint32_t                          splice_start1,
    uint32_t                          splice_end1,
    uint8_t                           interbase,
    uint32_t                          splice_region_exonic,
    uint32_t                          splice_region_intronic,
    duckvep_effect_ctx_t             *out);

/* Convert already-classified region and splice facts into the one predicate
 * bitset consumed by the generated consequence program. This is the join point
 * for VEP's feature-span region geometry and its multi-island splice geometry. */
void duckvep_effect_ctx_fill_classified(
    const duckvep_transcript_model_t *transcripts,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          region_start1,
    uint32_t                          region_end1,
    const duckvep_region_state_t     *region_state,
    const duckvep_splice_state_t     *splice_state,
    duckvep_effect_ctx_t             *out);

/* Sorted SNV hot path. The cursor is per transcript and survives adjacent
 * annotation tiles; UINT16_MAX means that this transcript has not yet been
 * visited in the current monotone run. */
void duckvep_effect_ctx_fill_point_sorted(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          pos,
    uint32_t                          splice_region_exonic,
    uint32_t                          splice_region_intronic,
    uint16_t                         *exon_rank_io,
    duckvep_effect_ctx_t             *out);

/* Apply variant-class facts from canonical event topology, not the caller's
 * broad transport kind. VEP's insertion/deletion predicates are allele-length
 * predicates after REF/ALT normalization/trimming; `kind` only chooses the broad
 * classification path. */
void duckvep_effect_ctx_apply_event(duckvep_effect_ctx_t *ctx,
                                    const duckvep_event_t *event);

void duckvep_effect_ctx_apply_delta(
    duckvep_effect_ctx_t           *ctx,
    const duckvep_sequence_delta_t *delta);

void duckvep_effect_ctx_apply_sv(
    duckvep_effect_ctx_t      *ctx,
    const duckvep_sv_effect_t *sv);

/* Derive VEP predicates that depend on the complete fact set. Call exactly once
 * after variant-class, structural, and optional sequence-delta facts are applied
 * and before evaluating the generated consequence program. */
void duckvep_effect_ctx_finalize(duckvep_effect_ctx_t *ctx);

typedef struct duckvep_nmd_result {
    uint8_t prediction;
    uint8_t escape_reasons;
} duckvep_nmd_result_t;

/* Apply the pinned VEP Plugins release/116 NMD.pm location rules after the SO
 * consequence set is known. Eligible splice consequences without a contiguous
 * coding projection are unresolved. */
void duckvep_nmd_predict(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    uint64_t                          consequence_mask,
    duckvep_nmd_result_t             *out);

typedef struct duckvep_consequence_rule {
    uint64_t required;
    uint64_t forbidden;
    uint64_t so_mask;
    uint8_t  tier;
    uint8_t  rank;
    uint8_t  impact;
} duckvep_consequence_rule_t;

uint64_t duckvep_effect_eval_rules(
    uint64_t                          pre_bits,
    const duckvep_consequence_rule_t *rules,
    size_t                            rule_count);

uint64_t duckvep_effect_eval(uint64_t pre_bits);

#endif /* DUCKVEP_EFFECT_H */
