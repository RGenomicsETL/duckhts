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
#include "duckvep_so.h"
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
    DUCKVEP_PRE_WITHIN_MATURE_MIRNA     = 40,
    DUCKVEP_PRE_PARTIAL_CODON            = 41,
    /* Facts owned by the separate regulatory/motif interval-feature lane. */
    DUCKVEP_PRE_TFBS_ABLATION             = 42,
    DUCKVEP_PRE_TFBS_AMPLIFICATION        = 43,
    DUCKVEP_PRE_WITHIN_TFBS               = 44,
    DUCKVEP_PRE_REGULATORY_ABLATION       = 45,
    DUCKVEP_PRE_REGULATORY_AMPLIFICATION  = 46,
    DUCKVEP_PRE_WITHIN_REGULATORY_REGION  = 47,
    /* VEP's raw equal-feature-length pre-predicate. It has no consequence of
     * its own, but filters frameshift and in-frame terms after transcript
     * projection. */
    DUCKVEP_PRE_SNP                        = 48,
    DUCKVEP_PRE_BIT_COUNT                  = 49
} duckvep_pre_bit_t;

#define DUCKVEP_PRE(b) (UINT64_C(1) << (b))

typedef struct duckvep_effect_ctx {
    uint32_t               variant_idx;
    uint32_t               tx_idx;
    uint32_t               start1;
    uint32_t               end1;
    int8_t                 strand;
    uint8_t                protein_coding;
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
    uint8_t                           repair_rewinds,
    uint16_t                         *exon_rank_io,
    duckvep_effect_ctx_t             *out);

/* Apply VEP's coarse OverlapConsequence `include => { exon => 0 }` gate.
 * `predicate_overlaps_exon` is deliberately distinct from exact genomic exon
 * overlap: VEP stretches every exon by 12 bases when the transcript contains
 * any frameshift intron. */
void duckvep_effect_ctx_apply_exon_gate(
    duckvep_effect_ctx_t *ctx,
    int                   predicate_overlaps_exon);

/* Apply variant-class facts from canonical event topology, not the caller's
 * broad transport kind. VEP's ordinary insertion/deletion/SNP predicates use
 * complete uploaded feature REF/ALT lengths before transcript projection.
 * The semantic differing-region lengths are the fallback only when that raw
 * relation is unavailable. */
void duckvep_effect_ctx_apply_event(
    const duckvep_transcript_model_t *transcripts,
    duckvep_effect_ctx_t             *ctx,
    const duckvep_event_t            *event);

void duckvep_effect_ctx_apply_delta(
    duckvep_effect_ctx_t           *ctx,
    const duckvep_sequence_delta_t *delta);

void duckvep_effect_ctx_apply_sv(
    duckvep_effect_ctx_t      *ctx,
    const duckvep_sv_effect_t *sv);

/* VEP's regulatory and motif consequences use the same event geometry and
 * structural operation facts as transcripts, but their features are ordinary
 * intervals. The resident feature sweep discovers candidates; this helper
 * evaluates one event/feature pair and returns only VEP SO bits.
 * `feature_start1` and `feature_end1` are one-based inclusive. */
uint64_t duckvep_effect_eval_interval_feature(
    duckvep_interval_feature_kind_t feature_kind,
    const duckvep_event_t          *event,
    uint16_t                        feature_chrom_id,
    uint32_t                        feature_start1,
    uint32_t                        feature_end1);

/* Derive VEP predicates that depend on the complete fact set. Call exactly once
 * after variant-class, structural, and optional sequence-delta facts are applied
 * and before evaluating the generated consequence program. */
void duckvep_effect_ctx_finalize(duckvep_effect_ctx_t *ctx);

typedef struct duckvep_nmd_result {
    uint8_t prediction;
    uint8_t escape_reasons;
} duckvep_nmd_result_t;

/* Executable thresholds in VEP Plugins release/116 NMD.pm. The names retain
 * the observed inclusive coordinate semantics instead of paraphrasing them as
 * 100/50-base prose. */
#define DUCKVEP_NMD_PENULTIMATE_EXON_OFFSET   51u

#define DUCKVEP_NMD_ELIGIBLE_SO_MASK ( \
    DUCKVEP_SO(DUCKVEP_SO_STOP_GAINED) | \
    DUCKVEP_SO(DUCKVEP_SO_FRAMESHIFT) | \
    DUCKVEP_SO(DUCKVEP_SO_SPLICE_DONOR) | \
    DUCKVEP_SO(DUCKVEP_SO_SPLICE_ACCEPTOR))

/* Apply the pinned VEP Plugins release/116 NMD.pm location rules after the SO
 * consequence set is known. `sequence_delta` may carry the first-101-CDS-base
 * comparison already resolved by coding classification; otherwise NMD performs
 * the exhaustive projection needed by splice and unresolved sequence cases. */
void duckvep_nmd_predict(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    uint64_t                          consequence_mask,
    const duckvep_sequence_delta_t   *sequence_delta,
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

/* Reference interpreter over the generated VEP rule program. Tests compare
 * the compact production lookup against this independent execution path. */
uint64_t duckvep_effect_eval_reference(uint64_t pre_bits);

/* StructuralVariationFeature inherits only VEP rules whose pinned metadata
 * names BaseVariationFeature. This is not equivalent to evaluating every
 * transcript predicate over a wide interval: splice, missense, synonymous and
 * several peptide predicates are ordinary-VariationFeature-only. */
uint64_t duckvep_effect_eval_structural_reference(uint64_t pre_bits);

uint64_t duckvep_effect_eval(uint64_t pre_bits);
uint64_t duckvep_effect_eval_structural(uint64_t pre_bits);

#endif /* DUCKVEP_EFFECT_H */
