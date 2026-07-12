/* duckvep_effect.c — VEP-shaped consequence evaluator. */
#include "duckvep_effect.h"

#include "duckvep_so.h"

static void duckvep_effect_ctx_set_topology(
    const duckvep_transcript_model_t *transcripts,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          start1,
    uint32_t                          end1,
    const duckvep_region_state_t     *region_state,
    const duckvep_splice_state_t     *splice_state,
    duckvep_effect_ctx_t             *out) {

    int coding = transcripts->cds_start1[tx_idx] != 0u;
    uint32_t region;
    uint64_t pre = 0u;

    out->variant_idx = variant_idx;
    out->tx_idx = (uint32_t)tx_idx;
    out->start1 = start1;
    out->end1 = end1;
    out->strand = transcripts->strand[tx_idx];

    out->region_state = *region_state;
    region = out->region_state.region_mask;

    /* Span placement is not mutually exclusive: one event can overlap CDS, UTR,
     * exon and intron compartments. */
    if (region & (uint32_t)DUCKVEP_REGION_UPSTREAM)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_UPSTREAM);
    if (region & (uint32_t)DUCKVEP_REGION_DOWNSTREAM)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_DOWNSTREAM);
    if (out->region_state.overlaps_cds)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_CDS);
    if (out->region_state.overlaps_utr5)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_UTR5);
    if (out->region_state.overlaps_utr3)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_UTR3);
    if (!coding && out->region_state.overlaps_exon)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_EXON);
    if (out->region_state.overlaps_intron)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_INTRON);

    pre |= coding ? DUCKVEP_PRE(DUCKVEP_PRE_CODING)
                  : DUCKVEP_PRE(DUCKVEP_PRE_NONCODING);

    /* VariationEffect::within_nmd_transcript is exactly the conjunction of
     * within_feature and the curated nonsense_mediated_decay biotype. It is
     * not a prediction that this allele creates a new NMD substrate. */
    if (out->region_state.within_feature &&
        (transcripts->flags[tx_idx] &
         (uint64_t)DUCKVEP_TX_BIOTYPE_NMD) != 0u) {
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_NMD_TRANSCRIPT);
    }

    out->splice = *splice_state;
    if (out->splice.splice_donor)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_DONOR);
    if (out->splice.splice_acceptor)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_ACCEPTOR);
    if (out->splice.splice_donor_5th)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_DONOR_5TH);
    if (out->splice.splice_donor_region)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_DONOR_REGION);
    /* VEP's splice_polypyrimidine_tract_variant OverlapConsequence carries
     * include => { exon => 0, intron => 1 }: the tract is reported only when the
     * variant does NOT overlap an exon. A deletion that runs from the tract across
     * the acceptor into the coding exon therefore keeps splice_acceptor but drops
     * the polypyrimidine term. (The essential donor/acceptor and the region terms
     * are gated on intron_boundary, not exon, so they are unaffected.) */
    if (out->splice.splice_polypyrimidine && !out->region_state.overlaps_exon)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_PPT);
    if (out->splice.splice_region)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_REGION);
    if (out->splice.intronic)
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_INTRON);

    region &= ~(uint32_t)DUCKVEP_REGION_SPLICE;
    if (out->splice.any) region |= (uint32_t)DUCKVEP_REGION_SPLICE;
    out->region = region;
    out->pre_bits = pre;
}

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
    duckvep_effect_ctx_t             *out) {

    duckvep_region_state_t region;
    duckvep_splice_state_t splice;

    region = duckvep_region_classify_span(transcripts, exons, tx_idx,
                                           start1, end1,
                                           splice_region_exonic,
                                           splice_region_intronic);
    splice = duckvep_splice_classify_span(transcripts, exons, tx_idx,
                                           start1, end1, interbase);
    duckvep_effect_ctx_set_topology(transcripts, variant_idx, tx_idx,
                                    start1, end1, &region, &splice, out);
}

void duckvep_effect_ctx_fill_point_sorted(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          pos,
    uint32_t                          splice_region_exonic,
    uint32_t                          splice_region_intronic,
    uint16_t                         *exon_rank_io,
    duckvep_effect_ctx_t             *out) {

    duckvep_region_state_t region;
    duckvep_splice_state_t splice;

    duckvep_classify_point_sorted(transcripts, exons, tx_idx, pos,
                                  splice_region_exonic,
                                  splice_region_intronic, exon_rank_io,
                                  &region, &splice);
    duckvep_effect_ctx_set_topology(transcripts, variant_idx, tx_idx,
                                    pos, pos, &region, &splice, out);
}

void duckvep_effect_ctx_apply_event(duckvep_effect_ctx_t *ctx,
                                    const duckvep_event_t *event) {

    if (ctx == NULL || event == NULL) return;
    if (event->kind == (uint8_t)DUCKVEP_KIND_SV) {
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_SV);
        return;
    }

    /* Mirrors Ensembl VariationEffect insertion/deletion predicates: inspect the
     * normalized allele length delta, not the broad declared kind. A lengthening
     * delins is an insertion predicate, a shortening delins is a deletion
     * predicate, and an equal-length substitution is neither. */
    if (event->ref_diff_length == 0u && event->alt_diff_length > 0u) {
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_INSERTION);
    } else if (event->ref_diff_length > 0u && event->alt_diff_length == 0u) {
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_DELETION);
    } else if (event->ref_diff_length > 0u && event->alt_diff_length > 0u) {
        if (event->alt_diff_length > event->ref_diff_length) {
            ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_INSERTION);
        } else if (event->ref_diff_length > event->alt_diff_length) {
            ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_DELETION);
        }
    }
}

void duckvep_effect_ctx_apply_delta(
    duckvep_effect_ctx_t           *ctx,
    const duckvep_sequence_delta_t *delta) {

    if (ctx == NULL || delta == NULL || !delta->valid) return;
    ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_DELTA);
    if (delta->synonymous) ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_SYNONYMOUS);
    if (delta->missense) ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_MISSENSE);
    if (delta->stop_gained) ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_STOP_GAINED);
    if (delta->stop_lost) ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_STOP_LOST);
    if (delta->stop_retained) ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_STOP_RETAINED);
    if (delta->start_lost) ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_START_LOST);
    if (delta->frameshift) ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_FRAMESHIFT);
    if (delta->inframe_deletion)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_INFRAME_DELETION);
    if (delta->inframe_insertion)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_INFRAME_INSERTION);
    if (delta->protein_altering)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_PROTEIN_ALTERING);
}

void duckvep_effect_ctx_apply_sv(
    duckvep_effect_ctx_t      *ctx,
    const duckvep_sv_effect_t *sv) {

    if (ctx == NULL || sv == NULL) return;
    if (sv->deletion)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_DELETION);
    if (sv->insertion)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_INSERTION);
    if (sv->feature_ablation)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_FEATURE_ABLATION);
    if (sv->feature_amplification)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_FEATURE_AMPLIFICATION);
    if (sv->feature_elongation)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_FEATURE_ELONGATION);
    if (sv->feature_truncation)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_FEATURE_TRUNCATION);
}

void duckvep_effect_ctx_finalize(duckvep_effect_ctx_t *ctx) {
    uint64_t pre;
    int coding;
    int noncoding;
    int noncoding_exon;

    if (ctx == NULL) return;
    pre = ctx->pre_bits;
    coding = (pre & DUCKVEP_PRE(DUCKVEP_PRE_CODING)) != 0u;
    noncoding = (pre & DUCKVEP_PRE(DUCKVEP_PRE_NONCODING)) != 0u;

    /* VariationEffect::non_coding_exon_variant is not just exon placement: a
     * complete feature overlap is explicitly excluded. This distinction is
     * observable for neutral CNVs, inversions, and unknown structural events. */
    noncoding_exon = noncoding && ctx->region_state.overlaps_exon &&
                     !ctx->region_state.complete_overlap_feature;
    if (noncoding_exon) {
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_NONCODING_EXON);
    }
    if (noncoding && ctx->region_state.within_feature && !noncoding_exon) {
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_NONCODING_GENE);
    }

    /* VariationEffect::coding_unknown returns false for complete feature
     * overlap. In that case coding_transcript_variant is the transcript-level
     * fallback. For the current supported slices, DELTA means a more specific
     * sequence predicate was resolved; no DELTA retains coding_unknown. */
    if (coding && ctx->region_state.complete_overlap_feature) {
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_CODING_TRANSCRIPT);
    } else if ((pre & DUCKVEP_PRE(DUCKVEP_PRE_CDS)) != 0u &&
               (pre & DUCKVEP_PRE(DUCKVEP_PRE_DELTA)) == 0u) {
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_CODING_UNKNOWN);
    }

    ctx->pre_bits = pre;
}

#include "duckvep_effect_rules.inc"

uint64_t duckvep_effect_eval_rules(
    uint64_t                          pre_bits,
    const duckvep_consequence_rule_t *rules,
    size_t                            rule_count) {

    uint64_t mask = 0u;
    uint8_t assigned_tier = 0u;
    size_t i;

    if (rules == NULL && rule_count > 0u) return 0u;
    for (i = 0u; i < rule_count; i++) {
        const duckvep_consequence_rule_t *r = &rules[i];
        if (assigned_tier != 0u && r->tier > assigned_tier) break;
        if ((pre_bits & r->required) == r->required &&
            (pre_bits & r->forbidden) == 0u) {
            mask |= r->so_mask;
            if (r->tier <= 2u && assigned_tier == 0u) assigned_tier = r->tier;
        }
    }
    return mask;
}

uint64_t duckvep_effect_eval(uint64_t pre_bits) {
    return duckvep_effect_eval_rules(pre_bits, k_rules,
                                     sizeof k_rules / sizeof k_rules[0]);
}
