/* duckvep_effect.c — VEP-shaped consequence evaluator. */
#include "duckvep_effect.h"

#include "duckvep_projection.h"
#include "duckvep_so.h"

void duckvep_nmd_predict(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    uint64_t                          consequence_mask,
    const duckvep_sequence_delta_t   *sequence_delta,
    duckvep_nmd_result_t             *out) {

    uint32_t exon_offset;
    uint16_t exon_count;
    uint32_t variant_end1;
    uint32_t variant_cds_start;
    uint32_t variant_cds_end;
    int early_cds_escape;
    uint8_t reasons = 0u;

    if (out == NULL) return;
    out->prediction = (uint8_t)DUCKVEP_NMD_NOT_APPLICABLE;
    out->escape_reasons = 0u;
    if (transcripts == NULL || exons == NULL || event == NULL ||
        (consequence_mask & DUCKVEP_NMD_ELIGIBLE_SO_MASK) == 0u) {
        return;
    }
    if (tx_idx >= transcripts->transcript_count ||
        transcripts->strand == NULL ||
        transcripts->cds_start1 == NULL ||
        transcripts->cds_end1 == NULL ||
        transcripts->exon_offset == NULL ||
        transcripts->exon_count == NULL ||
        exons->start1 == NULL || exons->end1 == NULL) {
        out->prediction = (uint8_t)DUCKVEP_NMD_UNRESOLVED;
        return;
    }
    if (transcripts->cds_start1[tx_idx] == 0u ||
        transcripts->cds_end1[tx_idx] == 0u) {
        out->prediction = (uint8_t)DUCKVEP_NMD_UNRESOLVED;
        return;
    }

    if (sequence_delta != NULL && sequence_delta->valid != 0u &&
        sequence_delta->nmd_early_cds_fact !=
            (uint8_t)DUCKVEP_NMD_EARLY_CDS_UNKNOWN) {
        early_cds_escape = sequence_delta->nmd_early_cds_fact ==
            (uint8_t)DUCKVEP_NMD_EARLY_CDS_ENDS_THROUGH_101;
    } else {
        /* NMD.pm returns no prediction unless TranscriptVariation has both CDS
         * coordinates. Keep that absence explicit for splice-crossing events. */
        if (!duckvep_project_feature_to_cds(
                transcripts, exons, tx_idx, event,
                &variant_cds_start, &variant_cds_end)) {
            out->prediction = (uint8_t)DUCKVEP_NMD_UNRESOLVED;
            return;
        }
        /* Projector success is the definedness test used by NMD.pm. A pure
         * insertion immediately before CDS base 1 is the valid reversed range
         * 1..0, and the plugin classifies its zero end as an early-CDS escape. */
        early_cds_escape =
            variant_cds_end <= DUCKVEP_NMD_EARLY_CDS_MAX_END1;
    }
    if (early_cds_escape) {
        reasons |= (uint8_t)DUCKVEP_NMD_ESCAPE_EARLY_CDS;
    }

    exon_offset = transcripts->exon_offset[tx_idx];
    exon_count = transcripts->exon_count[tx_idx];
    if (exon_count == 0u || exon_offset > exons->exon_count ||
        (size_t)exon_count > exons->exon_count - exon_offset) {
        out->prediction = (uint8_t)DUCKVEP_NMD_UNRESOLVED;
        return;
    }
    if (exon_count == 1u) {
        reasons |= (uint8_t)DUCKVEP_NMD_ESCAPE_INTRONLESS;
    }

    /* NMD.pm reads VariationFeature::seq_region_end directly. Equal-length
     * alleles retain unchanged uploaded bases, so the semantic edit endpoint
     * is not interchangeable with this coordinate. */
    variant_end1 = event->feature_end1;
    {
        size_t last = (size_t)exon_offset + (size_t)exon_count - 1u;
        if (variant_end1 >= exons->start1[last] &&
            variant_end1 <= exons->end1[last]) {
            reasons |= (uint8_t)DUCKVEP_NMD_ESCAPE_LAST_EXON;
        }
    }
    if (exon_count >= 2u) {
        size_t penultimate = (size_t)exon_offset + (size_t)exon_count - 2u;
        uint32_t lo;
        uint32_t hi;

        if (transcripts->strand[tx_idx] < 0) {
            lo = exons->start1[penultimate];
            hi = lo > UINT32_MAX - DUCKVEP_NMD_PENULTIMATE_EXON_OFFSET
                ? UINT32_MAX
                : lo + DUCKVEP_NMD_PENULTIMATE_EXON_OFFSET;
        } else {
            hi = exons->end1[penultimate];
            lo = hi > DUCKVEP_NMD_PENULTIMATE_EXON_OFFSET
                ? hi - DUCKVEP_NMD_PENULTIMATE_EXON_OFFSET : 0u;
        }
        if (variant_end1 >= lo && variant_end1 <= hi) {
            reasons |= (uint8_t)DUCKVEP_NMD_ESCAPE_PENULTIMATE_EXON_END;
        }
    }

    out->escape_reasons = reasons;
    out->prediction = reasons != 0u
        ? (uint8_t)DUCKVEP_NMD_PREDICTED_ESCAPING
        : (uint8_t)DUCKVEP_NMD_PREDICTED_TRIGGERING;
}

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
    out->protein_coding = (uint8_t)(
        (transcripts->flags[tx_idx] &
         (uint64_t)DUCKVEP_TX_BIOTYPE_PROTEIN_CODING) != 0u);

    out->region_state = *region_state;
    if (splice_state->within_frameshift_intron) {
        out->region_state.within_frameshift_intron = 1u;
        if (coding && end1 >= transcripts->cds_start1[tx_idx] &&
            start1 <= transcripts->cds_end1[tx_idx]) {
            out->region_state.overlaps_cds = 1u;
            out->region_state.region_mask |= (uint32_t)DUCKVEP_REGION_CDS;
        }
    }
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

    /* VariationEffect::within_mature_miRNA maps each curated cDNA range to
     * genomic exon segments, then applies the same inclusive overlap predicate
     * to the uploaded feature. The prepared model stores those mapped segments,
     * so the common non-miRNA lane pays one flag test and no projection cost. */
    if (out->region_state.within_feature &&
        (transcripts->flags[tx_idx] &
         (uint64_t)DUCKVEP_TX_BIOTYPE_MIRNA) != 0u &&
        transcripts->mature_mirna_offset != NULL) {
        size_t segment;
        size_t begin = transcripts->mature_mirna_offset[tx_idx];
        size_t finish = transcripts->mature_mirna_offset[tx_idx + 1u];

        for (segment = begin; segment < finish; segment++) {
            if (end1 >= transcripts->mature_mirna_start1[segment] &&
                start1 <= transcripts->mature_mirna_end1[segment]) {
                pre |= DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_MATURE_MIRNA);
                break;
            }
        }
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
    duckvep_effect_ctx_apply_exon_gate(
        out, out->region_state.overlaps_exon != 0u);
}

void duckvep_effect_ctx_apply_exon_gate(
    duckvep_effect_ctx_t *ctx,
    int                   predicate_overlaps_exon) {

    if (ctx == NULL || !predicate_overlaps_exon ||
        !ctx->splice.splice_polypyrimidine) {
        return;
    }

    ctx->splice.splice_polypyrimidine = 0u;
    ctx->pre_bits &= ~DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_PPT);
    ctx->splice.any = (uint8_t)(
        ctx->splice.splice_donor || ctx->splice.splice_acceptor ||
        ctx->splice.splice_donor_5th || ctx->splice.splice_donor_region ||
        ctx->splice.splice_region);
    if (ctx->splice.any) {
        ctx->region |= (uint32_t)DUCKVEP_REGION_SPLICE;
        ctx->region_state.region_mask |= (uint32_t)DUCKVEP_REGION_SPLICE;
    } else {
        ctx->region &= ~(uint32_t)DUCKVEP_REGION_SPLICE;
        ctx->region_state.region_mask &= ~(uint32_t)DUCKVEP_REGION_SPLICE;
    }
}

void duckvep_effect_ctx_fill_classified(
    const duckvep_transcript_model_t *transcripts,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          region_start1,
    uint32_t                          region_end1,
    const duckvep_region_state_t     *region_state,
    const duckvep_splice_state_t     *splice_state,
    duckvep_effect_ctx_t             *out) {

    duckvep_effect_ctx_set_topology(
        transcripts, variant_idx, tx_idx, region_start1, region_end1,
        region_state, splice_state, out);
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

    duckvep_effect_ctx_fill_geometry(
        transcripts, exons, variant_idx, tx_idx,
        start1, end1, start1, end1, interbase,
        splice_region_exonic, splice_region_intronic, out);
}

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
    duckvep_effect_ctx_t             *out) {

    duckvep_region_state_t region;
    duckvep_splice_state_t splice;

    region = duckvep_region_classify_span(transcripts, exons, tx_idx,
                                           region_start1, region_end1, 0u, 0u);
    splice = duckvep_splice_classify_span_with_windows(
        transcripts, exons, tx_idx, splice_start1, splice_end1, interbase,
        splice_region_exonic, splice_region_intronic);
    duckvep_effect_ctx_fill_classified(
        transcripts, variant_idx, tx_idx, region_start1, region_end1,
        &region, &splice, out);
}

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
    duckvep_effect_ctx_t             *out) {

    duckvep_region_state_t region;
    duckvep_splice_state_t splice;

    /* The exact VEP predicate consumes the configured generic-region reaches;
     * the coarse diagnostic pass is still unnecessary. */
    duckvep_classify_point_sorted(transcripts, exons, tx_idx, pos,
                                  splice_region_exonic,
                                  splice_region_intronic,
                                  repair_rewinds,
                                  exon_rank_io, &region, &splice);
    duckvep_effect_ctx_set_topology(transcripts, variant_idx, tx_idx,
                                    pos, pos, &region, &splice, out);
}

void duckvep_effect_ctx_apply_event(
    const duckvep_transcript_model_t *transcripts,
    duckvep_effect_ctx_t             *ctx,
    const duckvep_event_t            *event) {

    size_t segment;

    if (ctx == NULL || event == NULL) return;

    /* Most topology predicates consume the transcript-relative placement point
     * selected by the mapper. within_mature_miRNA is different: VEP compares
     * the raw VariationFeature coordinates with each projected miRNA segment.
     * A minimized insertion has the reversed interval (P+1,P), so an insertion
     * after the final mature base must not inherit that base's mature term. */
    if (transcripts != NULL &&
        ctx->tx_idx < transcripts->transcript_count &&
        (transcripts->flags[ctx->tx_idx] &
         (uint64_t)DUCKVEP_TX_BIOTYPE_MIRNA) != 0u &&
        transcripts->mature_mirna_offset != NULL) {
        size_t begin = transcripts->mature_mirna_offset[ctx->tx_idx];
        size_t finish = transcripts->mature_mirna_offset[ctx->tx_idx + 1u];

        ctx->pre_bits &= ~DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_MATURE_MIRNA);
        if (ctx->region_state.within_feature) {
            for (segment = begin; segment < finish; segment++) {
                if (event->feature_end1 >=
                        transcripts->mature_mirna_start1[segment] &&
                    event->feature_start1 <=
                        transcripts->mature_mirna_end1[segment]) {
                    ctx->pre_bits |=
                        DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_MATURE_MIRNA);
                    break;
                }
            }
        }
    }
    if (event->kind == (uint8_t)DUCKVEP_KIND_SV) {
        if (ctx->region_state.complete_overlap_feature) {
            /* BaseVariationFeatureOverlapAllele::_bvfo_preds returns early for
             * a structural event that contains the complete transcript. It
             * retains complete-overlap, within-feature and biotype facts, but
             * deliberately does not populate CDS, UTR, exon, intron or splice
             * preconditions. Transcript ablation/amplification or the coding/
             * non-coding transcript fallback therefore remain the only paths. */
            ctx->pre_bits &= ~(
                DUCKVEP_PRE(DUCKVEP_PRE_CDS) |
                DUCKVEP_PRE(DUCKVEP_PRE_UTR5) |
                DUCKVEP_PRE(DUCKVEP_PRE_UTR3) |
                DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_DONOR) |
                DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_ACCEPTOR) |
                DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_DONOR_5TH) |
                DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_DONOR_REGION) |
                DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_PPT) |
                DUCKVEP_PRE(DUCKVEP_PRE_SPLICE_REGION) |
                DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_INTRON));
        }
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_SV);
        return;
    }

    /* BaseVariationFeatureOverlapAllele classifies the VEP feature allele
     * before transcript projection. An equal-length genomic allele spanning
     * an intron remains SNP-class even when its outer mapped CDS replacement
     * changes length. */
    if (event->feature_length_relation ==
        (uint8_t)DUCKVEP_FEATURE_LENGTH_EQUAL) {
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_SNP);
    } else if (event->feature_length_relation ==
               (uint8_t)DUCKVEP_FEATURE_LENGTH_INCREASE) {
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_INSERTION);
    } else if (event->feature_length_relation ==
               (uint8_t)DUCKVEP_FEATURE_LENGTH_DECREASE) {
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_DELETION);
    } else if (event->ref_diff_length == 0u &&
               event->alt_diff_length > 0u) {
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_INSERTION);
    } else if (event->ref_diff_length > 0u &&
               event->alt_diff_length == 0u) {
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_DELETION);
    } else if (event->alt_diff_length > event->ref_diff_length) {
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_INSERTION);
    } else if (event->ref_diff_length > event->alt_diff_length) {
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_DELETION);
    }
    /* VariationEffect::feature_ablation is shared by ordinary
     * VariationFeatureOverlapAlleles and structural overlap alleles. The
     * predicate is exactly complete transcript overlap plus the normalized
     * deletion fact; it is not restricted to symbolic/SV input. This matters
     * for long literal alleles from pangenome call sets whose shortened ALT
     * contains an entire transcript. */
    if (ctx->region_state.complete_overlap_feature &&
        (ctx->pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_DELETION)) != 0u) {
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_FEATURE_ABLATION);
    }
}

void duckvep_effect_ctx_apply_delta(
    duckvep_effect_ctx_t           *ctx,
    const duckvep_sequence_delta_t *delta) {

    if (ctx == NULL || delta == NULL || !delta->valid) return;
    ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_DELTA);
    if (delta->synonymous) ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_SYNONYMOUS);
    if (delta->missense &&
        (ctx->pre_bits &
         (DUCKVEP_PRE(DUCKVEP_PRE_INSERTION) |
          DUCKVEP_PRE(DUCKVEP_PRE_DELETION))) == 0u) {
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_MISSENSE);
    }
    if (delta->stop_gained) ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_STOP_GAINED);
    if (delta->stop_lost) ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_STOP_LOST);
    if (delta->stop_retained) ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_STOP_RETAINED);
    if (delta->start_lost) ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_START_LOST);
    if (delta->start_retained) ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_START_RETAINED);
    if (delta->frameshift &&
        (ctx->pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_SNP)) == 0u) {
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_FRAMESHIFT);
    }
    if (delta->inframe_deletion &&
        (ctx->pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_DELETION)) != 0u) {
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_INFRAME_DELETION);
    }
    if (delta->inframe_insertion &&
        (ctx->pre_bits & DUCKVEP_PRE(DUCKVEP_PRE_INSERTION)) != 0u) {
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_INFRAME_INSERTION);
    }
    if (delta->protein_altering)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_PROTEIN_ALTERING);
    if (delta->coding_unknown)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_CODING_UNKNOWN);
    if (delta->partial_codon)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_PARTIAL_CODON);
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
    if (sv->start_lost)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_START_LOST);
    if (sv->start_retained)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_START_RETAINED);
    if (sv->stop_lost)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_STOP_LOST);
    if (sv->frameshift)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_FRAMESHIFT);
    if (sv->inframe_deletion)
        ctx->pre_bits |= DUCKVEP_PRE(DUCKVEP_PRE_INFRAME_DELETION);
}

static int duckvep_breakend_point_close_to_interval(
    uint16_t point_chrom_id,
    uint32_t point_pos1,
    uint16_t feature_chrom_id,
    uint32_t feature_start1,
    uint32_t feature_end1) {

    if (point_chrom_id != feature_chrom_id || point_pos1 == 0u) return 0;
    if (point_pos1 < feature_start1) {
        return feature_start1 - point_pos1 <=
            DUCKVEP_BREAKEND_ALLELE_DISTANCE;
    }
    if (point_pos1 > feature_end1) {
        return point_pos1 - feature_end1 <=
            DUCKVEP_BREAKEND_ALLELE_DISTANCE;
    }
    return 1;
}

uint64_t duckvep_effect_eval_interval_feature(
    duckvep_interval_feature_kind_t feature_kind,
    const duckvep_event_t          *event,
    uint16_t                        feature_chrom_id,
    uint32_t                        feature_start1,
    uint32_t                        feature_end1) {

    uint64_t pre = 0u;
    int within;
    int complete_overlap;
    int deletion;
    int copy_number_gain;

    if (event == NULL || feature_start1 == 0u ||
        feature_end1 < feature_start1 ||
        (feature_kind != DUCKVEP_INTERVAL_FEATURE_REGULATORY_REGION &&
         feature_kind != DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE)) {
        return 0u;
    }

    /* RegFeat::annotate_InputBuffer discovers a feature when either exact BND
     * point overlaps it, then StructuralVariationOverlap admits every endpoint
     * no farther than MAX_DISTANCE_FROM_TRANSCRIPT (5000 bases) from that same
     * feature. The local point has already undergone BaseVCF4 anchor removal;
     * the mate coordinate remains verbatim.
     *
     * VariationEffect's ordinary within-feature predicate reads the local
     * StructuralVariationFeature for every overlap allele. Therefore a local
     * exact hit wins even when the mate also hits. With a mate-only exact hit,
     * the mate allele contributes feature_truncation; a local endpoint in the
     * 5000-base admission halo but outside the feature contributes VEP's
     * intergenic fallback to the same object-level consequence set. */
    if (event->kind == (uint8_t)DUCKVEP_KIND_SV &&
        event->sv_type == (uint8_t)DUCKVEP_SV_BREAKEND) {
        int local_within = feature_chrom_id == event->chrom_id &&
            event->feature_start1 >= feature_start1 &&
            event->feature_start1 <= feature_end1;
        int mate_within = event->has_mate &&
            feature_chrom_id == event->mate_chrom_id &&
            event->mate_pos1 >= feature_start1 &&
            event->mate_pos1 <= feature_end1;

        if (local_within) {
            pre = feature_kind == DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE
                ? DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_TFBS)
                : DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_REGULATORY_REGION);
            return duckvep_effect_eval(pre);
        }
        if (mate_within) {
            uint64_t mask = duckvep_effect_eval(
                DUCKVEP_PRE(DUCKVEP_PRE_FEATURE_TRUNCATION));

            if (duckvep_breakend_point_close_to_interval(
                    event->chrom_id, event->feature_start1,
                    feature_chrom_id, feature_start1, feature_end1)) {
                mask |= DUCKVEP_SO(DUCKVEP_SO_INTERGENIC);
            }
            return mask;
        }
        return 0u;
    } else {
        if (feature_chrom_id != event->chrom_id) return 0u;
        within = event->feature_end1 >= feature_start1 &&
                 event->feature_start1 <= feature_end1;
        complete_overlap = event->feature_start1 <= feature_start1 &&
                           event->feature_end1 >= feature_end1;
    }

    /* These are the executable VariationEffect::within_feature and
     * complete_overlap_feature comparisons. Keeping the insertion's reversed
     * feature interval is intentional: an insertion exactly outside a feature
     * must not acquire its consequence from an arbitrary anchor base. */
    if (!within) return 0u;

    if (event->kind == (uint8_t)DUCKVEP_KIND_SV) {
        deletion = event->sv_type == (uint8_t)DUCKVEP_SV_DELETION ||
                   event->copy_change == (uint8_t)DUCKVEP_COPY_CHANGE_LOSS;
        copy_number_gain =
            event->sv_type == (uint8_t)DUCKVEP_SV_DUPLICATION ||
            event->sv_type == (uint8_t)DUCKVEP_SV_TANDEM_DUPLICATION ||
            event->sv_type == (uint8_t)DUCKVEP_SV_TANDEM_REPEAT ||
            event->copy_change == (uint8_t)DUCKVEP_COPY_CHANGE_GAIN;
    } else {
        deletion = event->ref_diff_length > event->alt_diff_length;
        /* VEP does not call every lengthening allele a copy-number gain.
         * Small tandem-repeat/duplication recognition needs allele content and
         * stays out until that typed event fact is available. */
        copy_number_gain = 0;
    }

    if (feature_kind == DUCKVEP_INTERVAL_FEATURE_TF_BINDING_SITE) {
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_TFBS);
        if (complete_overlap && deletion)
            pre |= DUCKVEP_PRE(DUCKVEP_PRE_TFBS_ABLATION);
        if (complete_overlap && copy_number_gain)
            pre |= DUCKVEP_PRE(DUCKVEP_PRE_TFBS_AMPLIFICATION);
    } else {
        pre |= DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_REGULATORY_REGION);
        if (complete_overlap && deletion)
            pre |= DUCKVEP_PRE(DUCKVEP_PRE_REGULATORY_ABLATION);
        if (complete_overlap && copy_number_gain)
            pre |= DUCKVEP_PRE(DUCKVEP_PRE_REGULATORY_AMPLIFICATION);
    }
    return duckvep_effect_eval(pre);
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
     * fallback. Usually a resolved DELTA suppresses coding_unknown, but VEP can
     * set both explicitly (for example an in-frame insertion whose local peptide
     * ends in X); duckvep_effect_ctx_apply_delta preserves that fact above. */
    if (ctx->region_state.complete_overlap_feature) {
        /* This predicate is false before VEP inspects any peptide state. A
         * sequence delta may independently contain X/coding-unknown facts, but
         * it must not reintroduce the term after complete overlap has been
         * established. */
        pre &= ~DUCKVEP_PRE(DUCKVEP_PRE_CODING_UNKNOWN);
        if (coding && ctx->protein_coding &&
            (pre & DUCKVEP_PRE(DUCKVEP_PRE_WITHIN_NMD_TRANSCRIPT)) == 0u) {
            pre |= DUCKVEP_PRE(DUCKVEP_PRE_CODING_TRANSCRIPT);
        }
    } else if ((pre & DUCKVEP_PRE(DUCKVEP_PRE_CDS)) != 0u &&
               (pre & DUCKVEP_PRE(DUCKVEP_PRE_DELTA)) == 0u &&
               ((pre & DUCKVEP_PRE(DUCKVEP_PRE_SV)) == 0u ||
                (pre & (DUCKVEP_PRE(DUCKVEP_PRE_FRAMESHIFT) |
                        DUCKVEP_PRE(DUCKVEP_PRE_INFRAME_INSERTION) |
                        DUCKVEP_PRE(DUCKVEP_PRE_INFRAME_DELETION))) == 0u)) {
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

uint64_t duckvep_effect_eval_reference(uint64_t pre_bits) {
    return duckvep_effect_eval_rules(pre_bits, k_rules,
                                     sizeof k_rules / sizeof k_rules[0]);
}

uint64_t duckvep_effect_eval_structural_reference(uint64_t pre_bits) {
    return duckvep_effect_eval_rules(
        pre_bits & k_rule_structural_pre_mask,
        k_rules, sizeof k_rules / sizeof k_rules[0]);
}

static unsigned duckvep_first_set_pre_bit(uint64_t mask) {
#if defined(__GNUC__) || defined(__clang__)
    return (unsigned)__builtin_ctzll(mask);
#else
    unsigned bit = 0u;
    while ((mask & UINT64_C(1)) == 0u) {
        mask >>= 1;
        bit++;
    }
    return bit;
#endif
}

uint64_t duckvep_effect_eval(uint64_t pre_bits) {
    uint64_t bits;
    uint64_t mask = 0u;

    bits = pre_bits & k_rule_tier1_pre_mask;
    if (bits == 0u) bits = pre_bits & k_rule_tier2_pre_mask;
    if (bits == 0u) bits = pre_bits & k_rule_fallback_pre_mask;
    while (bits != 0u) {
        mask |= k_rule_so_by_pre_bit[duckvep_first_set_pre_bit(bits)];
        bits &= bits - UINT64_C(1);
    }
    return mask;
}

uint64_t duckvep_effect_eval_structural(uint64_t pre_bits) {
    return duckvep_effect_eval(pre_bits & k_rule_structural_pre_mask);
}
