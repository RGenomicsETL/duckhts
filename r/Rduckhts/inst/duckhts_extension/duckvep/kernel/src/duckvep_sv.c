/* duckvep_sv.c — VEP-shaped structural/CNV predicates. */
#include "duckvep_sv.h"
#include "duckvep_projection.h"

#include <string.h>

static int sv_span_overlaps(uint32_t start1, uint32_t end1,
                            uint32_t other_start1, uint32_t other_end1) {
    return end1 >= other_start1 && start1 <= other_end1;
}

static int sv_feature_overlaps_start_codon(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event) {

    uint32_t cdna_start;
    uint32_t cdna_end;
    uint32_t coding_start_cdna;
    uint32_t coding_start_genomic;
    int start_mapped;
    int end_mapped;

    if (transcripts == NULL || exons == NULL || event == NULL ||
        tx_idx >= transcripts->transcript_count ||
        event->feature_start1 == 0u || event->feature_end1 == 0u) {
        return 0;
    }
    if (transcripts->cds_cdna_start1 != NULL &&
        transcripts->cds_cdna_start1[tx_idx] != 0u) {
        coding_start_cdna = transcripts->cds_cdna_start1[tx_idx];
    } else {
        coding_start_genomic = transcripts->strand[tx_idx] >= 0
            ? transcripts->cds_start1[tx_idx]
            : transcripts->cds_end1[tx_idx];
        if (coding_start_genomic == 0u ||
            !duckvep_project_genomic_to_cdna(
                transcripts, exons, tx_idx, coding_start_genomic,
                &coding_start_cdna, NULL)) {
            return 0;
        }
    }

    start_mapped = duckvep_project_genomic_to_cdna(
        transcripts, exons, tx_idx, event->feature_start1,
        &cdna_start, NULL);
    end_mapped = duckvep_project_genomic_to_cdna(
        transcripts, exons, tx_idx, event->feature_end1,
        &cdna_end, NULL);
    if (event->interbase) {
        if (!start_mapped && !end_mapped) return 0;
        if (start_mapped && end_mapped) {
            uint32_t lo = cdna_start < cdna_end ? cdna_start : cdna_end;
            uint32_t hi = cdna_start > cdna_end ? cdna_start : cdna_end;
            cdna_start = hi;
            cdna_end = lo;
        } else {
            uint32_t mapped = start_mapped ? cdna_start : cdna_end;
            int insertion_after_mapped =
                (end_mapped && transcripts->strand[tx_idx] >= 0) ||
                (start_mapped && transcripts->strand[tx_idx] < 0);

            /* Mapper::map_insert keeps start=end+1 even when only one flank
             * maps. Which side of that base owns the insertion depends on the
             * transcript strand; collapsing this to a point invents start
             * overlap at CDS and exon entrances. */
            if (insertion_after_mapped) {
                if (mapped == UINT32_MAX) return 0;
                cdna_start = mapped + 1u;
                cdna_end = mapped;
            } else {
                cdna_start = mapped;
                cdna_end = mapped - 1u;
            }
        }
    } else if (!start_mapped || !end_mapped) {
        return 0;
    } else if (cdna_start > cdna_end) {
        uint32_t tmp = cdna_start;
        cdna_start = cdna_end;
        cdna_end = tmp;
    }
    return cdna_end >= coding_start_cdna &&
           cdna_start <= coding_start_cdna + 2u;
}

static int sv_event_translation_start_is_partial_codon(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          translateable_cds_length,
    const duckvep_event_t            *event) {

    duckvep_coding_projection_t projection;

    if (translateable_cds_length == 0u ||
        !duckvep_project_feature_translation_start(
            transcripts, exons, tx_idx, event, &projection)) {
        return 0;
    }
    return duckvep_cds_position_is_partial_codon(
        (size_t)translateable_cds_length, projection.cds_pos);
}

static int sv_cds_overlap_within_one_exon(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event) {

    size_t offset;
    size_t count;
    size_t i;

    if (transcripts == NULL || exons == NULL || event == NULL ||
        tx_idx >= transcripts->transcript_count ||
        transcripts->cds_start1[tx_idx] == 0u ||
        event->end1 < transcripts->cds_start1[tx_idx] ||
        event->start1 > transcripts->cds_end1[tx_idx]) {
        return 0;
    }
    offset = transcripts->exon_offset[tx_idx];
    count = transcripts->exon_count[tx_idx];
    if (offset > exons->exon_count || count > exons->exon_count - offset) {
        return 0;
    }
    for (i = 0u; i < count; i++) {
        size_t exon_idx = offset + i;
        if (event->start1 >= exons->start1[exon_idx] &&
            event->end1 <= exons->end1[exon_idx]) {
            return 1;
        }
    }
    return 0;
}

int duckvep_sv_metadata_valid(duckvep_sv_type_t sv_type,
                              duckvep_copy_change_t copy_change) {
    if (sv_type <= DUCKVEP_SV_NONE ||
        sv_type > DUCKVEP_SV_TANDEM_REPEAT ||
        copy_change < DUCKVEP_COPY_CHANGE_UNKNOWN ||
        copy_change > DUCKVEP_COPY_CHANGE_GAIN) {
        return 0;
    }

    switch (sv_type) {
    case DUCKVEP_SV_DELETION:
        return copy_change == DUCKVEP_COPY_CHANGE_UNKNOWN ||
               copy_change == DUCKVEP_COPY_CHANGE_LOSS;
    case DUCKVEP_SV_DUPLICATION:
    case DUCKVEP_SV_TANDEM_DUPLICATION:
    case DUCKVEP_SV_TANDEM_REPEAT:
        return copy_change == DUCKVEP_COPY_CHANGE_UNKNOWN ||
               copy_change == DUCKVEP_COPY_CHANGE_GAIN;
    case DUCKVEP_SV_INSERTION:
        return copy_change != DUCKVEP_COPY_CHANGE_LOSS;
    case DUCKVEP_SV_INVERSION:
    case DUCKVEP_SV_BREAKEND:
        return copy_change == DUCKVEP_COPY_CHANGE_UNKNOWN ||
               copy_change == DUCKVEP_COPY_CHANGE_NEUTRAL;
    case DUCKVEP_SV_UNKNOWN:
    case DUCKVEP_SV_CNV:
        return 1;
    case DUCKVEP_SV_NONE:
    default:
        return 0;
    }
}

int duckvep_sv_geometry_valid(duckvep_sv_type_t sv_type,
                              uint32_t start1, uint32_t end1) {
    if (start1 == 0u || end1 < start1) return 0;
    if (sv_type == DUCKVEP_SV_INSERTION) {
        return start1 == end1 && start1 != UINT32_MAX;
    }
    if (sv_type == DUCKVEP_SV_BREAKEND) {
        return start1 == end1;
    }
    return 1;
}

duckvep_sv_effect_t duckvep_sv_effect_fill(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          translateable_cds_length,
    const duckvep_event_t            *event,
    const duckvep_region_state_t     *region) {

    duckvep_sv_effect_t st;
    memset(&st, 0, sizeof st);
    if (event == NULL || region == NULL ||
        event->kind != (uint8_t)DUCKVEP_KIND_SV) {
        return st;
    }

    st.copy_number_gain =
        (uint8_t)(event->copy_change == (uint8_t)DUCKVEP_COPY_CHANGE_GAIN ||
                  event->sv_type == (uint8_t)DUCKVEP_SV_DUPLICATION ||
                  event->sv_type == (uint8_t)DUCKVEP_SV_TANDEM_DUPLICATION ||
                  event->sv_type == (uint8_t)DUCKVEP_SV_TANDEM_REPEAT);
    st.copy_number_loss =
        (uint8_t)(event->copy_change == (uint8_t)DUCKVEP_COPY_CHANGE_LOSS ||
                  event->sv_type == (uint8_t)DUCKVEP_SV_DELETION);
    st.deletion =
        (uint8_t)(event->sv_type == (uint8_t)DUCKVEP_SV_DELETION ||
                  st.copy_number_loss);
    st.insertion =
        (uint8_t)(event->sv_type == (uint8_t)DUCKVEP_SV_INSERTION ||
                  st.copy_number_gain);
    st.chromosome_breakpoint =
        (uint8_t)(event->sv_type == (uint8_t)DUCKVEP_SV_BREAKEND);

    if (st.chromosome_breakpoint) {
        int local_inside = event->chrom_id == transcripts->chrom_id[tx_idx] &&
            event->feature_start1 >= transcripts->start1[tx_idx] &&
            event->feature_start1 <= transcripts->end1[tx_idx];
        int mate_inside = event->has_mate &&
            event->mate_chrom_id == transcripts->chrom_id[tx_idx] &&
            event->mate_pos1 >= transcripts->start1[tx_idx] &&
            event->mate_pos1 <= transcripts->end1[tx_idx];

        /* StructuralVariationOverlap builds one alternate allele for every
         * close endpoint. feature_truncation is the exceptional predicate that
         * inspects that allele's endpoint instead of the local feature. The
         * compact transcript row is their consequence-set union. */
        st.feature_truncation = (uint8_t)(local_inside || mate_inside);
        if (event->chrom_id != transcripts->chrom_id[tx_idx]) return st;
    }

    st.feature_ablation =
        (uint8_t)(region->complete_overlap_feature && st.deletion);
    st.feature_amplification =
        (uint8_t)(region->complete_overlap_feature && st.copy_number_gain);
    st.feature_elongation =
        (uint8_t)(region->within_cdna && region->complete_within_feature &&
                  (st.copy_number_gain || st.insertion));

    /* Other losses/deletions must overlap cDNA and be partial-overlap or wholly
     * within the feature. Breakpoint truncation was resolved from both loci
     * above and must not be overwritten by the local region state. */
    if (!st.chromosome_breakpoint) {
        st.feature_truncation =
            (uint8_t)(region->within_cdna &&
                      (region->partial_overlap_feature ||
                       region->complete_within_feature) &&
                      (st.copy_number_loss || st.deletion));
    }

    if (!region->complete_overlap_feature &&
        transcripts != NULL && tx_idx < transcripts->transcript_count &&
        transcripts->cds_start1[tx_idx] != 0u) {
        uint32_t cds_start1 = transcripts->cds_start1[tx_idx];
        uint32_t cds_end1 = transcripts->cds_end1[tx_idx];
        uint32_t start_codon_start1;
        uint32_t start_codon_end1;
        uint32_t stop_codon_start1;
        uint32_t stop_codon_end1;
        uint64_t flags = transcripts->flags[tx_idx];
        uint32_t event_length = event->end1 - event->start1 + 1u;
        int overlaps_start_codon = sv_feature_overlaps_start_codon(
            transcripts, exons, tx_idx, event);
        int cds_overlap_within_one_exon = sv_cds_overlap_within_one_exon(
            transcripts, exons, tx_idx, event);

        if (transcripts->strand[tx_idx] >= 0) {
            start_codon_start1 = cds_start1;
            start_codon_end1 = cds_end1 - cds_start1 >= 2u
                ? cds_start1 + 2u : cds_end1;
            stop_codon_start1 = cds_end1 - cds_start1 >= 2u
                ? cds_end1 - 2u : cds_start1;
            stop_codon_end1 = cds_end1;
        } else {
            start_codon_start1 = cds_end1 - cds_start1 >= 2u
                ? cds_end1 - 2u : cds_start1;
            start_codon_end1 = cds_end1;
            stop_codon_start1 = cds_start1;
            stop_codon_end1 = cds_end1 - cds_start1 >= 2u
                ? cds_start1 + 2u : cds_end1;
        }

        /* start_retained uses the cDNA-mapped codon, while structural
         * start_lost adds a separate contiguous genomic-codon overlap. A start
         * codon split across exons can therefore retain without being lost. */
        if ((flags & (uint64_t)DUCKVEP_TX_CDS_START_NF) == 0u &&
            overlaps_start_codon) {
            st.start_retained = 1u;
            if (sv_span_overlaps(event->start1, event->end1,
                                 start_codon_start1, start_codon_end1)) {
                st.start_lost = 1u;
            }
        }
        if (st.deletion && region->overlaps_cds &&
            !sv_event_translation_start_is_partial_codon(
                transcripts, exons, tx_idx, translateable_cds_length, event) &&
            sv_span_overlaps(event->start1, event->end1,
                             stop_codon_start1, stop_codon_end1)) {
            st.stop_lost = 1u;
        }

        /* The two VEP 116 structural predicates are intentionally asymmetric.
         * frameshift() only requires a deletion wholly inside one exon, with
         * some CDS overlap supplied by the consequence-rule precondition.
         * inframe_deletion() additionally requires cds_coords to contain one
         * Coordinate and no Gap, which means the whole span is inside the CDS.
         * A length-divisible deletion crossing from CDS into UTR is therefore
         * coding_sequence_variant, not inframe_deletion. */
        if (st.deletion && cds_overlap_within_one_exon) {
            st.frameshift = (uint8_t)(event_length % 3u != 0u);
            st.inframe_deletion = (uint8_t)(
                event_length % 3u == 0u &&
                event->start1 >= cds_start1 && event->end1 <= cds_end1);
        }
    }
    return st;
}
