/* duckvep_hgvs.c — typed HGVS facts over the shared transcript edit. */
#include "duckvep_hgvs.h"

#include "duckvep_codon.h"
#include "duckvep_dna.h"
#include "duckvep_projection.h"

#include <string.h>

static duckvep_hgvs_status_t hgvs_oriented_allele_base(
    const uint8_t *allele,
    size_t         length,
    int8_t         strand,
    size_t         rotation,
    size_t         index,
    uint8_t       *base_out) {

    size_t oriented_index;
    char base;

    if (base_out == NULL || allele == NULL || length == 0u ||
        index >= length ||
        (strand != (int8_t)1 && strand != (int8_t)-1)) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    index = (index + rotation % length) % length;
    oriented_index = strand > 0 ? index : length - 1u - index;
    base = duckvep_dna_normalize((char)allele[oriented_index], 0);
    if (base == '\0') return DUCKVEP_HGVS_INVALID_ALLELE;
    if (strand < 0) base = duckvep_dna_complement(base);
    *base_out = (uint8_t)base;
    return DUCKVEP_HGVS_OK;
}

static duckvep_hgvs_status_t hgvs_feature_is_inversion(
    const duckvep_transcript_edit_t *edit,
    int                              *is_inversion) {

    size_t i;

    if (is_inversion == NULL) return DUCKVEP_HGVS_INVALID_ARG;
    *is_inversion = 0;
    if (edit == NULL || edit->feature_ref_length <= 1u ||
        edit->feature_ref_length != edit->feature_alt_length ||
        edit->feature_ref == NULL || edit->feature_alt == NULL) {
        return DUCKVEP_HGVS_OK;
    }
    for (i = 0u; i < (size_t)edit->feature_ref_length; i++) {
        uint8_t alt_base;
        uint8_t ref_base;
        duckvep_hgvs_status_t status = hgvs_oriented_allele_base(
            edit->feature_alt, edit->feature_alt_length,
            edit->transcript_strand, 0u, i, &alt_base);
        if (status != DUCKVEP_HGVS_OK) return status;
        status = hgvs_oriented_allele_base(
            edit->feature_ref, edit->feature_ref_length,
            edit->transcript_strand, 0u,
            (size_t)edit->feature_ref_length - 1u - i, &ref_base);
        if (status != DUCKVEP_HGVS_OK) return status;
        if (alt_base != (uint8_t)duckvep_dna_complement((char)ref_base)) {
            return DUCKVEP_HGVS_OK;
        }
    }
    *is_inversion = 1;
    return DUCKVEP_HGVS_OK;
}

static duckvep_hgvs_status_t hgvs_coordinate_from_cdna_bounds(
    const duckvep_transcript_coordinate_t *coordinate,
    int                                     coding,
    uint32_t                                coding_start_cdna,
    uint32_t                                coding_end_cdna,
    duckvep_hgvs_coordinate_t              *out) {

    duckvep_hgvs_coordinate_t result;

    if (coordinate == NULL || out == NULL ||
        coordinate->cdna_anchor1 == 0u ||
        (coordinate->exonic != 0u && coordinate->intron_offset != 0)) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    memset(&result, 0, sizeof result);
    result.intron_offset = coordinate->intron_offset;
    if (!coding) {
        result.kind = (uint8_t)DUCKVEP_HGVS_COORDINATE_N;
        result.base = (int64_t)coordinate->cdna_anchor1;
        *out = result;
        return DUCKVEP_HGVS_OK;
    }
    if (coding_start_cdna == 0u || coding_end_cdna < coding_start_cdna) {
        return DUCKVEP_HGVS_INVALID_PROJECTION;
    }
    if (coordinate->cdna_anchor1 > coding_end_cdna) {
        result.kind = (uint8_t)DUCKVEP_HGVS_COORDINATE_C_STAR;
        result.base = (int64_t)coordinate->cdna_anchor1 -
                      (int64_t)coding_end_cdna;
    } else if (coordinate->cdna_anchor1 == coding_end_cdna &&
               coordinate->intron_offset != 0) {
        result.kind = (uint8_t)DUCKVEP_HGVS_COORDINATE_C_STAR;
        result.base = 0;
    } else {
        result.kind = (uint8_t)DUCKVEP_HGVS_COORDINATE_C;
        result.base = coordinate->cdna_anchor1 >= coding_start_cdna
            ? (int64_t)coordinate->cdna_anchor1 -
              (int64_t)coding_start_cdna + 1
            : (int64_t)coordinate->cdna_anchor1 -
              (int64_t)coding_start_cdna;
        if (result.base == 0) return DUCKVEP_HGVS_INVALID_PROJECTION;
    }
    *out = result;
    return DUCKVEP_HGVS_OK;
}

duckvep_hgvs_status_t duckvep_hgvs_coordinate_from_transcript(
    const duckvep_transcript_model_t      *transcripts,
    const duckvep_exon_model_t            *exons,
    size_t                                 tx_idx,
    const duckvep_transcript_coordinate_t *coordinate,
    duckvep_hgvs_coordinate_t             *out) {

    uint32_t coding_start_cdna;
    uint32_t coding_end_cdna;
    int coding;

    if (out == NULL) return DUCKVEP_HGVS_INVALID_ARG;
    memset(out, 0, sizeof *out);
    if (transcripts == NULL || exons == NULL || coordinate == NULL ||
        tx_idx >= transcripts->transcript_count ||
        transcripts->cds_start1 == NULL || transcripts->cds_end1 == NULL ||
        coordinate->cdna_anchor1 == 0u ||
        (coordinate->exonic != 0u && coordinate->intron_offset != 0)) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    coding = transcripts->cds_start1[tx_idx] != 0u ||
             transcripts->cds_end1[tx_idx] != 0u;
    if (transcripts->cds_start1[tx_idx] == 0u ||
        transcripts->cds_end1[tx_idx] == 0u) {
        if (coding) return DUCKVEP_HGVS_INVALID_PROJECTION;
        coding_start_cdna = 0u;
        coding_end_cdna = 0u;
    } else if (!duckvep_project_coding_cdna_bounds(
            transcripts, exons, tx_idx, &coding_start_cdna,
            &coding_end_cdna, NULL, NULL)) {
        return DUCKVEP_HGVS_INVALID_PROJECTION;
    }
    return hgvs_coordinate_from_cdna_bounds(
        coordinate, coding, coding_start_cdna, coding_end_cdna, out);
}

duckvep_hgvs_status_t duckvep_hgvs_dna_fact_build(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_transcript_edit_t  *edit,
    duckvep_hgvs_dna_fact_t          *out) {

    duckvep_hgvs_dna_fact_t result;
    duckvep_hgvs_status_t status;
    uint32_t coding_start_cdna;
    uint32_t coding_end_cdna;
    uint32_t cds_pos;
    uint8_t phase_offset;
    int coding;

    if (out == NULL) return DUCKVEP_HGVS_INVALID_ARG;
    memset(out, 0, sizeof *out);
    if (transcripts == NULL || exons == NULL || edit == NULL ||
        transcripts->strand == NULL ||
        (size_t)edit->tx_idx >= transcripts->transcript_count ||
        edit->transcript_strand != transcripts->strand[edit->tx_idx] ||
        (edit->transcript_strand != (int8_t)1 &&
         edit->transcript_strand != (int8_t)-1) ||
        (edit->ref_length != 0u && edit->ref == NULL) ||
        (edit->alt_length != 0u && edit->alt == NULL) ||
        (edit->ref_length == 0u && edit->alt_length == 0u)) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    memset(&result, 0, sizeof result);
    coding = transcripts->cds_start1[edit->tx_idx] != 0u ||
             transcripts->cds_end1[edit->tx_idx] != 0u;
    coding_start_cdna = 0u;
    coding_end_cdna = 0u;
    phase_offset = 0u;
    if (coding &&
        (transcripts->cds_start1[edit->tx_idx] == 0u ||
         transcripts->cds_end1[edit->tx_idx] == 0u ||
         !duckvep_project_coding_cdna_bounds(
             transcripts, exons, edit->tx_idx, &coding_start_cdna,
             &coding_end_cdna, NULL, &phase_offset))) {
        return DUCKVEP_HGVS_INVALID_PROJECTION;
    }
    status = hgvs_coordinate_from_cdna_bounds(
        &edit->first, coding, coding_start_cdna, coding_end_cdna,
        &result.first);
    if (status != DUCKVEP_HGVS_OK) return status;
    status = hgvs_coordinate_from_cdna_bounds(
        &edit->last, coding, coding_start_cdna, coding_end_cdna,
        &result.last);
    if (status != DUCKVEP_HGVS_OK) return status;

    /* TranscriptVariationAllele::hgvs_transcript bypasses
     * _get_cDNA_position only for a literal exonic VariationFeature SNP. Its
     * TranscriptVariation::cds_start already includes the positive starting
     * phase, so preserve that representation-specific fast path instead of
     * shifting generic indel, MNV, or intronic coordinates. */
    if (edit->feature_ref_length == 1u &&
        edit->feature_alt_length == 1u &&
        edit->feature_first.exonic != 0u &&
        edit->feature_last.exonic != 0u &&
        edit->feature_first.genomic_pos1 != 0u &&
        edit->feature_first.genomic_pos1 ==
            edit->feature_last.genomic_pos1 && coding &&
        edit->feature_first.cdna_anchor1 >= coding_start_cdna &&
        edit->feature_first.cdna_anchor1 <= coding_end_cdna) {
        cds_pos = edit->feature_first.cdna_anchor1 - coding_start_cdna + 1u +
            (uint32_t)phase_offset;
        result.first.kind = (uint8_t)DUCKVEP_HGVS_COORDINATE_C;
        result.first.base = (int64_t)cds_pos;
        result.first.intron_offset = 0;
        result.last = result.first;
    }

    result.numbering = result.first.kind ==
            (uint8_t)DUCKVEP_HGVS_COORDINATE_N
        ? (uint8_t)DUCKVEP_HGVS_NUMBERING_N
        : (uint8_t)DUCKVEP_HGVS_NUMBERING_C;
    if ((result.numbering == (uint8_t)DUCKVEP_HGVS_NUMBERING_N) !=
        (result.last.kind == (uint8_t)DUCKVEP_HGVS_COORDINATE_N)) {
        return DUCKVEP_HGVS_INVALID_PROJECTION;
    }
    if (edit->ref_length == 1u && edit->alt_length == 1u) {
        result.shape = (uint8_t)DUCKVEP_HGVS_DNA_SUBSTITUTION;
    } else if (edit->alt_length == 0u) {
        result.shape = (uint8_t)DUCKVEP_HGVS_DNA_DELETION;
    } else if (edit->ref_length == 0u) {
        result.shape = (uint8_t)DUCKVEP_HGVS_DNA_INSERTION;
    } else {
        result.shape = (uint8_t)DUCKVEP_HGVS_DNA_REPLACEMENT;
    }
    result.ref = edit->ref;
    result.alt = edit->alt;
    result.ref_length = edit->ref_length;
    result.alt_length = edit->alt_length;
    result.placed_start1 = edit->event.start1;
    result.placed_end1 = edit->event.end1;
    result.placed_insertion_boundary0 = edit->event.interbase != 0u
        ? edit->event.insertion_boundary0 : 0u;
    result.transcript_strand = edit->transcript_strand;
    *out = result;
    return DUCKVEP_HGVS_OK;
}

static duckvep_hgvs_status_t hgvs_reference_base(
    const duckvep_hgvs_reference_window_t *reference,
    uint32_t                                position1,
    uint8_t                                *base_out) {

    uint64_t offset;
    char base;

    if (base_out == NULL) return DUCKVEP_HGVS_INVALID_ARG;
    *base_out = 0u;
    if (reference == NULL || reference->bases == NULL ||
        reference->length == 0u || reference->start1 == 0u ||
        position1 < reference->start1) {
        return DUCKVEP_HGVS_MISSING_REFERENCE;
    }
    offset = (uint64_t)position1 - (uint64_t)reference->start1;
    if (offset >= (uint64_t)reference->length) {
        return DUCKVEP_HGVS_MISSING_REFERENCE;
    }
    base = duckvep_dna_normalize((char)reference->bases[(size_t)offset], 0);
    if (base == '\0' || base == 'N') return DUCKVEP_HGVS_INVALID_ALLELE;
    *base_out = (uint8_t)base;
    return DUCKVEP_HGVS_OK;
}

duckvep_hgvs_status_t duckvep_hgvs_genomic_search_interval(
    const duckvep_event_t *event,
    uint32_t               sequence_length,
    uint32_t              *start1_out,
    uint32_t              *end1_out) {

    uint64_t variation_start1;
    uint64_t variation_end1;
    uint64_t fetch_start1;
    uint64_t fetch_end1;
    int shiftable;

    if (start1_out == NULL || end1_out == NULL) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    *start1_out = 0u;
    *end1_out = 0u;
    if (event == NULL || sequence_length == 0u ||
        ((event->raw_start1 != 0u || event->raw_end1 != 0u) &&
         (event->raw_start1 == 0u || event->raw_end1 < event->raw_start1 ||
          event->raw_end1 > sequence_length))) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    shiftable = 0;
    if (event->kind == (uint8_t)DUCKVEP_KIND_INS) {
        if (event->interbase == 0u ||
            event->insertion_boundary0 > sequence_length) {
            return DUCKVEP_HGVS_INVALID_PROJECTION;
        }
        variation_start1 = (uint64_t)event->insertion_boundary0 + 1u;
        variation_end1 = event->insertion_boundary0;
        shiftable = 1;
    } else if (event->kind == (uint8_t)DUCKVEP_KIND_DEL) {
        if (event->interbase != 0u || event->start1 == 0u ||
            event->end1 < event->start1 ||
            event->end1 > sequence_length) {
            return DUCKVEP_HGVS_INVALID_PROJECTION;
        }
        variation_start1 = event->start1;
        variation_end1 = event->end1;
        shiftable = 1;
    } else if ((event->kind == (uint8_t)DUCKVEP_KIND_SNV ||
                event->kind == (uint8_t)DUCKVEP_KIND_MNV ||
                event->kind == (uint8_t)DUCKVEP_KIND_INDEL) &&
               event->interbase == 0u && event->start1 != 0u &&
               event->end1 >= event->start1 &&
               event->end1 <= sequence_length &&
               event->ref_diff_length != 0u &&
               (uint64_t)event->end1 - (uint64_t)event->start1 + 1u ==
                   (uint64_t)event->ref_diff_length) {
        /* Substitutions and replacements do not enter VEP's 3-prime shift
         * loop, but hgvs_transcript still derives their REF from the
         * reference slice. Fetch exactly the semantic REF span so the adapter
         * can prove that assertion without paying for a +/-1000 window. */
        variation_start1 = event->start1;
        variation_end1 = event->end1;
    } else {
        return DUCKVEP_HGVS_UNSUPPORTED_EDIT;
    }

    if (shiftable) {
        fetch_start1 = variation_start1 > DUCKVEP_HGVS_SHIFT_LIMIT
            ? variation_start1 - DUCKVEP_HGVS_SHIFT_LIMIT : 1u;
        fetch_end1 = variation_end1 + DUCKVEP_HGVS_SHIFT_LIMIT;
        if (fetch_end1 > sequence_length) fetch_end1 = sequence_length;
    } else {
        fetch_start1 = variation_start1;
        fetch_end1 = variation_end1;
    }
    if (fetch_start1 == 0u || fetch_start1 > fetch_end1 ||
        fetch_end1 > UINT32_MAX) {
        return DUCKVEP_HGVS_OUT_OF_RANGE;
    }
    *start1_out = (uint32_t)fetch_start1;
    *end1_out = (uint32_t)fetch_end1;
    return DUCKVEP_HGVS_OK;
}

duckvep_hgvs_status_t duckvep_hgvs_reference_fetch_interval(
    const duckvep_event_t *event,
    uint32_t               sequence_length,
    uint32_t              *start1_out,
    uint32_t              *end1_out) {

    duckvep_hgvs_status_t status;
    uint32_t shift_start1;
    uint32_t shift_end1;
    uint64_t fetch_start1;
    uint64_t fetch_end1;
    int shiftable;

    if (start1_out == NULL || end1_out == NULL) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    *start1_out = 0u;
    *end1_out = 0u;
    status = duckvep_hgvs_genomic_search_interval(
        event, sequence_length, &shift_start1, &shift_end1);
    if (status != DUCKVEP_HGVS_OK) return status;

    fetch_start1 = shift_start1;
    fetch_end1 = shift_end1;
    if (event->raw_start1 != 0u && event->raw_start1 < fetch_start1) {
        fetch_start1 = event->raw_start1;
    }
    if (event->raw_end1 != 0u && event->raw_end1 > fetch_end1) {
        fetch_end1 = event->raw_end1;
    }

    shiftable = event->kind == (uint8_t)DUCKVEP_KIND_INS ||
        event->kind == (uint8_t)DUCKVEP_KIND_DEL;
    if (event->feature_start1 != 0u &&
        event->feature_end1 >= event->feature_start1) {
        uint64_t feature_start1 = event->feature_start1;
        uint64_t feature_end1 = event->feature_end1;

        if (feature_end1 > sequence_length) {
            return DUCKVEP_HGVS_INVALID_ARG;
        }
        if (shiftable) {
            feature_start1 = feature_start1 > DUCKVEP_HGVS_SHIFT_LIMIT
                ? feature_start1 - DUCKVEP_HGVS_SHIFT_LIMIT : 1u;
            feature_end1 += DUCKVEP_HGVS_SHIFT_LIMIT;
            if (feature_end1 > sequence_length) {
                feature_end1 = sequence_length;
            }
        }
        if (feature_start1 < fetch_start1) fetch_start1 = feature_start1;
        if (feature_end1 > fetch_end1) fetch_end1 = feature_end1;
    }

    if (event->kind == (uint8_t)DUCKVEP_KIND_INS) {
        uint64_t source_start1;
        uint64_t source_end1;
        uint64_t pattern_length = event->alt_diff_length;

        if (pattern_length == 0u) {
            return DUCKVEP_HGVS_INVALID_PROJECTION;
        }
        source_start1 = shift_start1 > pattern_length
            ? (uint64_t)shift_start1 - pattern_length : 1u;
        source_end1 = (uint64_t)shift_end1 + pattern_length;
        if (source_end1 > sequence_length) source_end1 = sequence_length;
        if (source_start1 < fetch_start1) fetch_start1 = source_start1;
        if (source_end1 > fetch_end1) fetch_end1 = source_end1;
    }

    if (fetch_start1 == 0u || fetch_start1 > fetch_end1 ||
        fetch_end1 > sequence_length || fetch_end1 > UINT32_MAX) {
        return DUCKVEP_HGVS_OUT_OF_RANGE;
    }
    *start1_out = (uint32_t)fetch_start1;
    *end1_out = (uint32_t)fetch_end1;
    return DUCKVEP_HGVS_OK;
}

duckvep_hgvs_status_t duckvep_hgvs_uploaded_reference_validate(
    const duckvep_hgvs_reference_window_t *reference,
    const duckvep_event_t                 *event,
    const uint8_t                         *raw_ref,
    size_t                                 raw_ref_length) {

    size_t i;

    if (reference == NULL || event == NULL || raw_ref == NULL ||
        raw_ref_length == 0u || event->raw_start1 == 0u ||
        event->raw_end1 < event->raw_start1 ||
        (uint64_t)event->raw_end1 - (uint64_t)event->raw_start1 + 1u !=
            (uint64_t)raw_ref_length ||
        reference->chrom_id != event->chrom_id) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    for (i = 0u; i < raw_ref_length; i++) {
        uint8_t uploaded;
        uint8_t reference_base;
        duckvep_hgvs_status_t status = hgvs_oriented_allele_base(
            raw_ref, raw_ref_length, (int8_t)1, 0u, i, &uploaded);
        if (status != DUCKVEP_HGVS_OK) return status;
        status = hgvs_reference_base(
            reference, event->raw_start1 + (uint32_t)i, &reference_base);
        if (status != DUCKVEP_HGVS_OK) return status;
        if (uploaded != reference_base) {
            return DUCKVEP_HGVS_REFERENCE_MISMATCH;
        }
    }
    return DUCKVEP_HGVS_OK;
}

static duckvep_hgvs_status_t hgvs_project_genomic_pair(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                             tx_idx,
    uint32_t                           genomic_low1,
    uint32_t                           genomic_high1,
    int8_t                             strand,
    duckvep_hgvs_coordinate_t        *first_out,
    duckvep_hgvs_coordinate_t        *last_out) {

    duckvep_transcript_coordinate_t low;
    duckvep_transcript_coordinate_t high;
    duckvep_hgvs_status_t status;

    if (first_out == NULL || last_out == NULL || genomic_low1 == 0u ||
        genomic_high1 < genomic_low1) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    if (duckvep_project_transcript_coordinate(
            transcripts, exons, tx_idx, genomic_low1, &low) !=
            DUCKVEP_TRANSCRIPT_EDIT_OK ||
        duckvep_project_transcript_coordinate(
            transcripts, exons, tx_idx, genomic_high1, &high) !=
            DUCKVEP_TRANSCRIPT_EDIT_OK) {
        return DUCKVEP_HGVS_OUT_OF_RANGE;
    }
    status = duckvep_hgvs_coordinate_from_transcript(
        transcripts, exons, tx_idx, strand > 0 ? &low : &high, first_out);
    if (status != DUCKVEP_HGVS_OK) return status;
    return duckvep_hgvs_coordinate_from_transcript(
        transcripts, exons, tx_idx, strand > 0 ? &high : &low, last_out);
}

/* VEP 116 clamps a VariationFeature that overlaps only part of a transcript,
 * derives REF from that clamped transcript slice, then clips common REF/ALT
 * bases. This is observably different from clipping the uploaded alleles
 * first: a terminal CG>GC can become c.*10_*9insG. */
static duckvep_hgvs_status_t hgvs_dna_fact_build_clamped_feature(
    const duckvep_transcript_model_t      *transcripts,
    const duckvep_exon_model_t            *exons,
    const duckvep_hgvs_reference_window_t *reference,
    const duckvep_transcript_edit_t       *edit,
    uint32_t                               shift_offset,
    duckvep_hgvs_dna_fact_t               *out) {

    duckvep_hgvs_dna_fact_t result;
    duckvep_transcript_coordinate_t first_transcript;
    duckvep_transcript_coordinate_t last_transcript;
    duckvep_hgvs_status_t status;
    const uint8_t *reference_allele;
    const uint8_t *alternate_allele;
    size_t reference_length;
    size_t alternate_length;
    size_t prefix = 0u;
    size_t suffix = 0u;
    size_t reference_offset;
    size_t alternate_offset;
    uint32_t genomic_low1;
    uint32_t genomic_high1;
    uint32_t genomic_first1;
    uint32_t genomic_last1;
    uint32_t repeat_count = 0u;
    int8_t strand;

    if (out == NULL) return DUCKVEP_HGVS_INVALID_ARG;
    memset(out, 0, sizeof *out);
    if (transcripts == NULL || exons == NULL || reference == NULL ||
        edit == NULL || (size_t)edit->tx_idx >= transcripts->transcript_count ||
        edit->event.interbase != 0u || edit->feature_alt == NULL ||
        edit->feature_first.genomic_pos1 == 0u ||
        edit->feature_last.genomic_pos1 == 0u ||
        reference->chrom_id != edit->event.chrom_id) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    strand = edit->transcript_strand;
    genomic_low1 = edit->feature_first.genomic_pos1 <
            edit->feature_last.genomic_pos1
        ? edit->feature_first.genomic_pos1
        : edit->feature_last.genomic_pos1;
    genomic_high1 = edit->feature_first.genomic_pos1 >
            edit->feature_last.genomic_pos1
        ? edit->feature_first.genomic_pos1
        : edit->feature_last.genomic_pos1;
    /* _var2transcript_slice_coords clamps the uploaded feature first, while
     * hgvs_transcript later adds the shift computed from the complete genomic
     * deletion. Keep that ordering: shift this fixed-width clamped slice in
     * transcript 3-prime direction, then derive its displayed REF. */
    if (strand > 0) {
        if (genomic_low1 > UINT32_MAX - shift_offset ||
            genomic_high1 > UINT32_MAX - shift_offset) {
            return DUCKVEP_HGVS_OUT_OF_RANGE;
        }
        genomic_low1 += shift_offset;
        genomic_high1 += shift_offset;
    } else {
        if (genomic_low1 <= shift_offset || genomic_high1 <= shift_offset) {
            return DUCKVEP_HGVS_NOT_APPLICABLE;
        }
        genomic_low1 -= shift_offset;
        genomic_high1 -= shift_offset;
    }
    /* hgvs_transcript() adds the shift to the already clamped transcript
     * slice.  If that carries either displayed endpoint outside the
     * transcript, VEP returns undef before it attempts to read the reference
     * sequence.  Report that ordinary absence instead of mislabelling a
     * complete indexed reference as missing. */
    if (genomic_low1 < transcripts->start1[edit->tx_idx] ||
        genomic_high1 > transcripts->end1[edit->tx_idx]) {
        return DUCKVEP_HGVS_NOT_APPLICABLE;
    }
    if (genomic_low1 < reference->start1 ||
        (uint64_t)genomic_high1 - (uint64_t)reference->start1 >=
            (uint64_t)reference->length) {
        return DUCKVEP_HGVS_MISSING_REFERENCE;
    }
    reference_length = (size_t)((uint64_t)genomic_high1 -
        (uint64_t)genomic_low1 + 1u);
    alternate_length = (size_t)edit->feature_alt_length;
    if (reference_length > UINT16_MAX || alternate_length > UINT16_MAX) {
        return DUCKVEP_HGVS_OUT_OF_RANGE;
    }
    reference_allele = reference->bases +
        (size_t)(genomic_low1 - reference->start1);
    alternate_allele = edit->feature_alt;

    /* hgvs_variant_notation() tests multiplication before _clip_alleles().
     * That order matters for a feature clamped at a transcript endpoint:
     * terminal CG>CC is a direct C>CC duplication, while CGT>CAC first
     * becomes an insertion only after clipping and then has no projectable
     * second coordinate.  Do not promote the latter from insAC to a
     * duplication merely because AC happens to copy terminal sequence. */
    if (reference_length != 0u &&
        alternate_length > reference_length &&
        alternate_length % reference_length == 0u) {
        size_t i;
        int repeated = 1;

        repeat_count = (uint32_t)(alternate_length / reference_length);
        for (i = 0u; i < alternate_length; i++) {
            uint8_t reference_base;
            uint8_t alternate_base;

            status = hgvs_oriented_allele_base(
                reference_allele, reference_length, strand, 0u,
                i % reference_length, &reference_base);
            if (status != DUCKVEP_HGVS_OK) return status;
            status = hgvs_oriented_allele_base(
                alternate_allele, alternate_length, strand, 0u, i,
                &alternate_base);
            if (status != DUCKVEP_HGVS_OK) return status;
            if (reference_base != alternate_base) {
                repeated = 0;
                break;
            }
        }
        if (!repeated) repeat_count = 0u;
    }
    if (repeat_count >= 2u) {
        status = hgvs_project_genomic_pair(
            transcripts, exons, edit->tx_idx, genomic_low1, genomic_high1,
            strand, &result.first, &result.last);
        if (status != DUCKVEP_HGVS_OK) return status;
        result.numbering = result.first.kind ==
                (uint8_t)DUCKVEP_HGVS_COORDINATE_N
            ? (uint8_t)DUCKVEP_HGVS_NUMBERING_N
            : (uint8_t)DUCKVEP_HGVS_NUMBERING_C;
        if ((result.numbering == (uint8_t)DUCKVEP_HGVS_NUMBERING_N) !=
            (result.last.kind == (uint8_t)DUCKVEP_HGVS_COORDINATE_N)) {
            return DUCKVEP_HGVS_INVALID_PROJECTION;
        }
        result.ref = reference_allele;
        result.alt = alternate_allele;
        result.ref_length = (uint16_t)reference_length;
        result.alt_length = (uint16_t)alternate_length;
        result.repeat_count = repeat_count;
        result.shape = repeat_count == 2u
            ? (uint8_t)DUCKVEP_HGVS_DNA_DUPLICATION
            : (uint8_t)DUCKVEP_HGVS_DNA_REPEAT;
        result.transcript_strand = strand;
        result.shift_offset = (int32_t)shift_offset;
        result.placed_start1 = genomic_low1;
        result.placed_end1 = genomic_high1;
        *out = result;
        return DUCKVEP_HGVS_OK;
    }

    while (prefix < reference_length && prefix < alternate_length) {
        uint8_t ref_base;
        uint8_t alt_base;
        status = hgvs_oriented_allele_base(
            reference_allele, reference_length, strand, 0u, prefix,
            &ref_base);
        if (status != DUCKVEP_HGVS_OK) return status;
        status = hgvs_oriented_allele_base(
            alternate_allele, alternate_length, strand, 0u, prefix,
            &alt_base);
        if (status != DUCKVEP_HGVS_OK) return status;
        if (ref_base != alt_base) break;
        prefix++;
    }
    while (suffix < reference_length - prefix &&
           suffix < alternate_length - prefix) {
        uint8_t ref_base;
        uint8_t alt_base;
        status = hgvs_oriented_allele_base(
            reference_allele, reference_length, strand, 0u,
            reference_length - 1u - suffix, &ref_base);
        if (status != DUCKVEP_HGVS_OK) return status;
        status = hgvs_oriented_allele_base(
            alternate_allele, alternate_length, strand, 0u,
            alternate_length - 1u - suffix, &alt_base);
        if (status != DUCKVEP_HGVS_OK) return status;
        if (ref_base != alt_base) break;
        suffix++;
    }
    reference_length -= prefix + suffix;
    alternate_length -= prefix + suffix;
    if (reference_length == 0u && alternate_length == 0u) {
        return DUCKVEP_HGVS_NOT_APPLICABLE;
    }

    if (strand > 0) {
        genomic_first1 = genomic_low1 + (uint32_t)prefix;
        genomic_last1 = genomic_high1 - (uint32_t)suffix;
        reference_offset = prefix;
        alternate_offset = prefix;
    } else {
        genomic_first1 = genomic_high1 - (uint32_t)prefix;
        genomic_last1 = genomic_low1 + (uint32_t)suffix;
        reference_offset = suffix;
        alternate_offset = suffix;
    }
    if (duckvep_project_transcript_coordinate(
            transcripts, exons, edit->tx_idx, genomic_first1,
            &first_transcript) != DUCKVEP_TRANSCRIPT_EDIT_OK ||
        duckvep_project_transcript_coordinate(
            transcripts, exons, edit->tx_idx, genomic_last1,
            &last_transcript) != DUCKVEP_TRANSCRIPT_EDIT_OK) {
        /* look_for_slice_start clamps a partially overlapping feature, but
         * hgvs_transcript still returns undef when clipping or the added
         * shift carries the printable edit outside that transcript slice. */
        return DUCKVEP_HGVS_NOT_APPLICABLE;
    } else {
        status = duckvep_hgvs_coordinate_from_transcript(
            transcripts, exons, edit->tx_idx, &first_transcript,
            &result.first);
        if (status != DUCKVEP_HGVS_OK) return status;
        status = duckvep_hgvs_coordinate_from_transcript(
            transcripts, exons, edit->tx_idx, &last_transcript,
            &result.last);
        if (status != DUCKVEP_HGVS_OK) return status;
    }
    result.numbering = result.first.kind ==
            (uint8_t)DUCKVEP_HGVS_COORDINATE_N
        ? (uint8_t)DUCKVEP_HGVS_NUMBERING_N
        : (uint8_t)DUCKVEP_HGVS_NUMBERING_C;
    if ((result.numbering == (uint8_t)DUCKVEP_HGVS_NUMBERING_N) !=
        (result.last.kind == (uint8_t)DUCKVEP_HGVS_COORDINATE_N)) {
        return DUCKVEP_HGVS_INVALID_PROJECTION;
    }
    result.ref = reference_allele + reference_offset;
    result.alt = alternate_allele + alternate_offset;
    result.ref_length = (uint16_t)reference_length;
    result.alt_length = (uint16_t)alternate_length;
    if (reference_length == 1u && alternate_length == 1u) {
        result.shape = (uint8_t)DUCKVEP_HGVS_DNA_SUBSTITUTION;
    } else if (alternate_length == 0u) {
        result.shape = (uint8_t)DUCKVEP_HGVS_DNA_DELETION;
    } else if (reference_length == 0u) {
        result.shape = (uint8_t)DUCKVEP_HGVS_DNA_INSERTION;
    } else {
        result.shape = (uint8_t)DUCKVEP_HGVS_DNA_REPLACEMENT;
    }
    result.transcript_strand = strand;
    result.shift_offset = (int32_t)shift_offset;
    if (reference_length == 0u) {
        result.placed_insertion_boundary0 = genomic_first1 < genomic_last1
            ? genomic_first1 : genomic_last1;
        result.placed_start1 = genomic_first1;
        result.placed_end1 = genomic_first1;
    } else {
        result.placed_start1 = genomic_first1 < genomic_last1
            ? genomic_first1 : genomic_last1;
        result.placed_end1 = genomic_first1 > genomic_last1
            ? genomic_first1 : genomic_last1;
    }
    *out = result;
    return DUCKVEP_HGVS_OK;
}

duckvep_hgvs_status_t
duckvep_hgvs_dna_fact_build_genomic_shifted_with_lookup(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_hgvs_reference_window_t *shift_reference,
    const duckvep_hgvs_reference_window_t *lookup_reference,
    const duckvep_transcript_edit_t  *edit,
    duckvep_hgvs_dna_fact_t          *out) {

    duckvep_hgvs_dna_fact_t result;
    duckvep_hgvs_status_t status;
    size_t pattern_length;
    uint32_t event_low1;
    uint32_t event_high1;
    uint32_t shifted_low1;
    uint32_t shifted_high1;
    uint32_t first_flank1;
    uint32_t flank_available;
    uint32_t pre_length;
    uint32_t post_length;
    uint32_t post_offset;
    uint32_t limit;
    uint32_t shift = 0u;
    int alternate;
    int is_inversion;
    int clamped_feature;

    if (out == NULL) return DUCKVEP_HGVS_INVALID_ARG;
    memset(out, 0, sizeof *out);
    status = duckvep_hgvs_dna_fact_build(transcripts, exons, edit, &result);
    if (status != DUCKVEP_HGVS_OK) return status;
    if (result.ref_length != 0u && lookup_reference != NULL) {
        size_t i;

        if (lookup_reference->chrom_id != edit->event.chrom_id ||
            edit->event.interbase != 0u || edit->event.start1 == 0u ||
            edit->event.end1 < edit->event.start1 ||
            (uint64_t)edit->event.end1 - (uint64_t)edit->event.start1 + 1u !=
                (uint64_t)result.ref_length) {
            return DUCKVEP_HGVS_INVALID_PROJECTION;
        }
        for (i = 0u; i < (size_t)result.ref_length; i++) {
            uint8_t uploaded;
            uint8_t reference_base;

            status = hgvs_oriented_allele_base(
                result.ref, result.ref_length, (int8_t)1, 0u, i,
                &uploaded);
            if (status != DUCKVEP_HGVS_OK) return status;
            status = hgvs_reference_base(
                lookup_reference, edit->event.start1 + (uint32_t)i,
                &reference_base);
            if (status != DUCKVEP_HGVS_OK) return status;
            if (uploaded != reference_base) {
                return DUCKVEP_HGVS_REFERENCE_MISMATCH;
            }
        }
    }
    clamped_feature = edit->event.interbase == 0u &&
        edit->event.feature_start1 != 0u &&
        edit->event.feature_end1 >= edit->event.feature_start1 &&
        (edit->event.feature_start1 <
             transcripts->start1[edit->tx_idx] ||
         edit->event.feature_end1 > transcripts->end1[edit->tx_idx]);
    if (clamped_feature &&
        result.shape != (uint8_t)DUCKVEP_HGVS_DNA_DELETION) {
        return hgvs_dna_fact_build_clamped_feature(
            transcripts, exons, lookup_reference, edit, 0u, out);
    }
    status = hgvs_feature_is_inversion(edit, &is_inversion);
    if (status != DUCKVEP_HGVS_OK) return status;
    if (is_inversion && result.ref_length > 1u &&
        result.ref_length == result.alt_length) {
        result.shape = (uint8_t)DUCKVEP_HGVS_DNA_INVERSION;
    }
    if (result.shape != (uint8_t)DUCKVEP_HGVS_DNA_DELETION &&
        result.shape != (uint8_t)DUCKVEP_HGVS_DNA_INSERTION) {
        *out = result;
        return DUCKVEP_HGVS_OK;
    }
    alternate = result.shape == (uint8_t)DUCKVEP_HGVS_DNA_INSERTION;
    pattern_length = alternate ? (size_t)result.alt_length
                               : (size_t)result.ref_length;
    if (pattern_length == 0u || shift_reference == NULL ||
        lookup_reference == NULL ||
        shift_reference->chrom_id != edit->event.chrom_id ||
        lookup_reference->chrom_id != edit->event.chrom_id ||
        shift_reference->bases == NULL || shift_reference->length == 0u ||
        lookup_reference->bases == NULL || lookup_reference->length == 0u ||
        shift_reference->start1 == 0u || lookup_reference->start1 == 0u ||
        shift_reference->length > UINT32_MAX ||
        lookup_reference->length > UINT32_MAX ||
        (uint64_t)shift_reference->start1 +
                (uint64_t)shift_reference->length - 1u > UINT32_MAX ||
        (uint64_t)lookup_reference->start1 +
                (uint64_t)lookup_reference->length - 1u >
            UINT32_MAX) {
        return shift_reference == NULL || lookup_reference == NULL ||
               shift_reference->bases == NULL ||
               lookup_reference->bases == NULL ||
               shift_reference->length == 0u ||
               lookup_reference->length == 0u
            ? DUCKVEP_HGVS_MISSING_REFERENCE
            : DUCKVEP_HGVS_INVALID_ARG;
    }
    if (alternate) {
        if (edit->event.interbase == 0u ||
            edit->event.insertion_boundary0 == UINT32_MAX ||
            edit->event.insertion_boundary0 == 0u) {
            return DUCKVEP_HGVS_INVALID_PROJECTION;
        }
        event_low1 = edit->event.insertion_boundary0;
        event_high1 = event_low1 + 1u;
    } else {
        if (edit->event.interbase != 0u || edit->event.start1 == 0u ||
            edit->event.end1 < edit->event.start1 ||
            (uint64_t)edit->event.end1 - (uint64_t)edit->event.start1 + 1u !=
                (uint64_t)pattern_length) {
            return DUCKVEP_HGVS_INVALID_PROJECTION;
        }
        event_low1 = edit->event.start1;
        event_high1 = edit->event.end1;
    }
    /* _genomic_shift does not slice pre/post relative to the event after the
     * +/-1000 region has been constrained to the sequence region. It takes
     * the first 1000 bytes as pre_seq and the last 1000 as post_seq. Those
     * strings overlap on short sequence regions and near either endpoint;
     * reproducing that observable VEP 116 behavior is required for exact
     * HGVS compatibility rather than a biologically cleaner byte walk. */
    pre_length = shift_reference->length < (size_t)DUCKVEP_HGVS_SHIFT_LIMIT
        ? (uint32_t)shift_reference->length
        : (uint32_t)DUCKVEP_HGVS_SHIFT_LIMIT;
    post_length = pre_length;
    post_offset = (uint32_t)shift_reference->length - post_length;
    if (edit->transcript_strand > 0) {
        if ((uint64_t)shift_reference->start1 + (uint64_t)post_offset >
            UINT32_MAX) {
            return DUCKVEP_HGVS_OUT_OF_RANGE;
        }
        first_flank1 = shift_reference->start1 + post_offset;
        flank_available = post_length;
    } else {
        if (pre_length == 0u ||
            (uint64_t)shift_reference->start1 +
                (uint64_t)pre_length - 1u >
                UINT32_MAX) {
            return DUCKVEP_HGVS_OUT_OF_RANGE;
        }
        first_flank1 = shift_reference->start1 + pre_length - 1u;
        flank_available = pre_length;
    }
    /* TranscriptVariationAllele::perform_shift leaves pattern_length - 1
     * reference bases in the fetched flank when the flank is at least as long
     * as the allele. This slightly odd loop limit is part of VEP 116 output. */
    limit = flank_available < (uint32_t)pattern_length
        ? flank_available
        : flank_available - (uint32_t)pattern_length + 1u;

    while (shift < limit) {
        uint8_t pattern_base;
        uint8_t reference_base;
        size_t pattern_index = edit->transcript_strand > 0
            ? (size_t)(shift % (uint32_t)pattern_length)
            : pattern_length - 1u -
                (size_t)(shift % (uint32_t)pattern_length);
        uint32_t reference_position = edit->transcript_strand > 0
            ? first_flank1 + shift : first_flank1 - shift;
        const uint8_t *pattern = alternate ? result.alt : result.ref;
        status = hgvs_oriented_allele_base(
            pattern, pattern_length, (int8_t)1, 0u, pattern_index,
            &pattern_base);
        if (status != DUCKVEP_HGVS_OK) return status;
        status = hgvs_reference_base(
            shift_reference, reference_position, &reference_base);
        if (status != DUCKVEP_HGVS_OK) return status;
        if (pattern_base != reference_base) break;
        shift++;
    }
    if (edit->transcript_strand > 0) {
        if (event_low1 > UINT32_MAX - shift ||
            event_high1 > UINT32_MAX - shift) {
            return DUCKVEP_HGVS_OUT_OF_RANGE;
        }
        shifted_low1 = event_low1 + shift;
        shifted_high1 = event_high1 + shift;
    } else {
        if (event_low1 <= shift || event_high1 <= shift) {
            return DUCKVEP_HGVS_OUT_OF_RANGE;
        }
        shifted_low1 = event_low1 - shift;
        shifted_high1 = event_high1 - shift;
    }
    if (clamped_feature) {
        return hgvs_dna_fact_build_clamped_feature(
            transcripts, exons, lookup_reference, edit, shift, out);
    }
    result.shift_offset = (int32_t)shift;
    if (alternate) {
        uint8_t anchor_side = edit->event.anchor_side;
        if (anchor_side != (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT &&
            anchor_side != (uint8_t)DUCKVEP_EVENT_ANCHOR_RIGHT) {
            return DUCKVEP_HGVS_INVALID_PROJECTION;
        }
        result.placed_insertion_boundary0 = shifted_low1;
        result.placed_start1 = anchor_side ==
                (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT
            ? shifted_low1 : shifted_high1;
        result.placed_end1 = result.placed_start1;
    } else {
        result.placed_start1 = shifted_low1;
        result.placed_end1 = shifted_high1;
        result.placed_insertion_boundary0 = 0u;
    }

    if (alternate) {
        int direction;
        int duplication_found = 0;

        /* hgvs_variant_notation() checks the transcript 3-prime sequence
         * first and the 5-prime sequence second.  This order matters on the
         * short-region overlap path where _genomic_shift can leave a tandem
         * insertion unshifted but either adjacent copy can still name a dup. */
        for (direction = 0; direction < 2; direction++) {
            uint32_t source_low1;
            uint32_t source_high1;
            size_t i;
            int duplicated = 1;
            int higher_genomic_side =
                (edit->transcript_strand > 0) == (direction == 0);

            if (higher_genomic_side) {
                if ((uint64_t)shifted_low1 + (uint64_t)pattern_length >
                        UINT32_MAX) {
                    duplicated = 0;
                    source_low1 = source_high1 = 0u;
                } else {
                    source_low1 = shifted_low1 + 1u;
                    source_high1 =
                        shifted_low1 + (uint32_t)pattern_length;
                }
            } else if ((uint64_t)shifted_low1 + 1u <
                           (uint64_t)pattern_length) {
                duplicated = 0;
                source_low1 = source_high1 = 0u;
            } else {
                source_high1 = shifted_low1;
                source_low1 = source_high1 -
                    (uint32_t)pattern_length + 1u;
            }
            for (i = 0u; duplicated && i < pattern_length; i++) {
                uint8_t alt_base;
                uint8_t reference_base;
                uint32_t reference_position = edit->transcript_strand > 0
                    ? source_low1 + (uint32_t)i
                    : source_high1 - (uint32_t)i;
                status = duckvep_hgvs_dna_base(&result, 1, i, &alt_base);
                if (status != DUCKVEP_HGVS_OK) return status;
                status = hgvs_reference_base(
                    lookup_reference, reference_position, &reference_base);
                if (status != DUCKVEP_HGVS_OK) {
                    duplicated = 0;
                    break;
                }
                if (edit->transcript_strand < 0) {
                    reference_base = (uint8_t)duckvep_dna_complement(
                        (char)reference_base);
                }
                if (alt_base != reference_base) duplicated = 0;
            }
            if (duplicated) {
                status = hgvs_project_genomic_pair(
                    transcripts, exons, edit->tx_idx, source_low1,
                    source_high1, edit->transcript_strand,
                    &result.first, &result.last);
                if (status == DUCKVEP_HGVS_OUT_OF_RANGE) continue;
                if (status != DUCKVEP_HGVS_OK) return status;
                result.shape = (uint8_t)DUCKVEP_HGVS_DNA_DUPLICATION;
                duplication_found = 1;
                break;
            }
        }
        /* The shared transcript edit already projected the original two
         * insertion flanks. Re-project only when the HGVS-only shift moved
         * them; a duplication names its copied source span above. */
        if (!duplication_found && shift != 0u) {
            status = hgvs_project_genomic_pair(
                transcripts, exons, edit->tx_idx, shifted_low1,
                shifted_high1, edit->transcript_strand,
                &result.first, &result.last);
            if (status != DUCKVEP_HGVS_OK) {
                /* hgvs_transcript returns undef when its shift offset carries
                 * either insertion flank beyond the transcript slice. This
                 * is an ordinary absent annotation, not an unsupported HGVS
                 * operation. A duplication is different: VEP names its
                 * copied source span, so a source wholly inside the
                 * transcript remains printable at a transcript endpoint. */
                return status == DUCKVEP_HGVS_OUT_OF_RANGE
                    ? DUCKVEP_HGVS_NOT_APPLICABLE : status;
            }
        }
    } else if (shift != 0u) {
        /* Zero shift retains the coordinates built from the shared transcript
         * edit. Only moved deletion endpoints need another exon lookup. */
        status = hgvs_project_genomic_pair(
            transcripts, exons, edit->tx_idx, shifted_low1, shifted_high1,
            edit->transcript_strand, &result.first, &result.last);
        if (status != DUCKVEP_HGVS_OK) {
            /* hgvs_transcript returns undef when its shift offset carries
             * either deletion endpoint beyond the transcript slice. */
            return status == DUCKVEP_HGVS_OUT_OF_RANGE
                ? DUCKVEP_HGVS_NOT_APPLICABLE : status;
        }
    }
    *out = result;
    return DUCKVEP_HGVS_OK;
}

duckvep_hgvs_status_t duckvep_hgvs_dna_fact_build_genomic_shifted(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_hgvs_reference_window_t *reference,
    const duckvep_transcript_edit_t  *edit,
    duckvep_hgvs_dna_fact_t          *out) {

    return duckvep_hgvs_dna_fact_build_genomic_shifted_with_lookup(
        transcripts, exons, reference, reference, edit, out);
}

static duckvep_hgvs_status_t hgvs_materialize_genomic_allele(
    const duckvep_hgvs_dna_fact_t *fact,
    int                             alternate,
    uint8_t                        *scratch,
    size_t                          length) {

    size_t i;

    if (fact == NULL || scratch == NULL || length == 0u ||
        (fact->transcript_strand != (int8_t)1 &&
         fact->transcript_strand != (int8_t)-1)) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    for (i = 0u; i < length; i++) {
        size_t oriented_index = fact->transcript_strand > 0
            ? i : length - 1u - i;
        uint8_t base;
        duckvep_hgvs_status_t status = duckvep_hgvs_dna_base(
            fact, alternate, oriented_index, &base);
        if (status != DUCKVEP_HGVS_OK) return status;
        scratch[i] = fact->transcript_strand > 0
            ? base : (uint8_t)duckvep_dna_complement((char)base);
    }
    return DUCKVEP_HGVS_OK;
}

static duckvep_hgvs_status_t hgvs_status_from_cds_edit(
    duckvep_cds_edit_status_t status) {

    switch (status) {
        case DUCKVEP_CDS_EDIT_OK:
            return DUCKVEP_HGVS_OK;
        case DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL:
            return DUCKVEP_HGVS_BUFFER_TOO_SMALL;
        case DUCKVEP_CDS_EDIT_INVALID_ALLELE:
            return DUCKVEP_HGVS_INVALID_ALLELE;
        case DUCKVEP_CDS_EDIT_REF_MISMATCH:
            return DUCKVEP_HGVS_REFERENCE_MISMATCH;
        case DUCKVEP_CDS_EDIT_UNSUPPORTED_KIND:
            return DUCKVEP_HGVS_UNSUPPORTED_EDIT;
        case DUCKVEP_CDS_EDIT_OUT_OF_CDS:
        case DUCKVEP_CDS_EDIT_NON_CONTIGUOUS:
            return DUCKVEP_HGVS_INVALID_PROJECTION;
        case DUCKVEP_CDS_EDIT_INVALID_ARG:
        case DUCKVEP_CDS_EDIT_INVALID_EVENT:
        default:
            return DUCKVEP_HGVS_INVALID_ARG;
    }
}

static duckvep_hgvs_status_t hgvs_shifted_event_build(
    const duckvep_transcript_edit_t *edit,
    const duckvep_hgvs_dna_fact_t   *fact,
    duckvep_event_t                  *out) {

    duckvep_event_t event;
    int insertion;

    if (out == NULL) return DUCKVEP_HGVS_INVALID_ARG;
    memset(out, 0, sizeof *out);
    if (edit == NULL || fact == NULL || fact->shift_offset < 0 ||
        fact->transcript_strand != edit->transcript_strand) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    insertion = fact->ref_length == 0u && fact->alt_length != 0u;
    if (!insertion && !(fact->ref_length != 0u && fact->alt_length == 0u)) {
        return DUCKVEP_HGVS_UNSUPPORTED_EDIT;
    }

    event = edit->event;
    event.raw_start1 = fact->placed_start1;
    event.raw_end1 = fact->placed_end1;
    event.ref_diff_offset = 0u;
    event.alt_diff_offset = 0u;
    event.anchor_ref_offset = 0u;
    event.feature_allele_offset = 0u;
    event.has_mate = 0u;
    event.mate_pos1 = 0u;
    event.mate_chrom_id = 0u;
    event.sv_type = 0u;
    event.copy_change = 0u;

    if (insertion) {
        uint32_t right_flank1;
        uint32_t anchor_position1;

        if ((fact->shape != (uint8_t)DUCKVEP_HGVS_DNA_INSERTION &&
             fact->shape != (uint8_t)DUCKVEP_HGVS_DNA_DUPLICATION) ||
            fact->placed_insertion_boundary0 == 0u ||
            fact->placed_insertion_boundary0 == UINT32_MAX ||
            (event.anchor_side != (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT &&
             event.anchor_side != (uint8_t)DUCKVEP_EVENT_ANCHOR_RIGHT)) {
            return DUCKVEP_HGVS_INVALID_PROJECTION;
        }
        right_flank1 = fact->placed_insertion_boundary0 + 1u;
        anchor_position1 = event.anchor_side ==
                (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT
            ? fact->placed_insertion_boundary0 : right_flank1;
        event.kind = (uint8_t)DUCKVEP_KIND_INS;
        event.interbase = 1u;
        event.insertion_boundary0 = fact->placed_insertion_boundary0;
        event.start1 = anchor_position1;
        event.end1 = anchor_position1;
        event.feature_start1 = right_flank1;
        event.feature_end1 = fact->placed_insertion_boundary0;
        event.ref_diff_length = 0u;
        event.alt_diff_length = fact->alt_length;
    } else {
        if (fact->shape != (uint8_t)DUCKVEP_HGVS_DNA_DELETION ||
            fact->placed_start1 == 0u ||
            fact->placed_end1 < fact->placed_start1 ||
            (uint64_t)fact->placed_end1 - (uint64_t)fact->placed_start1 + 1u !=
                (uint64_t)fact->ref_length) {
            return DUCKVEP_HGVS_INVALID_PROJECTION;
        }
        event.kind = (uint8_t)DUCKVEP_KIND_DEL;
        event.interbase = 0u;
        event.anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_NONE;
        event.insertion_boundary0 = 0u;
        event.start1 = fact->placed_start1;
        event.end1 = fact->placed_end1;
        event.feature_start1 = fact->placed_start1;
        event.feature_end1 = fact->placed_end1;
        event.ref_diff_length = fact->ref_length;
        event.alt_diff_length = 0u;
    }
    *out = event;
    return DUCKVEP_HGVS_OK;
}

duckvep_hgvs_status_t duckvep_hgvs_protein_coordinates_defined(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_transcript_edit_t  *edit,
    const duckvep_hgvs_dna_fact_t    *fact,
    int                              *defined_out) {

    uint32_t first_cds;
    uint32_t last_cds;

    if (defined_out == NULL) return DUCKVEP_HGVS_INVALID_ARG;
    *defined_out = 0;
    if (transcripts == NULL || exons == NULL || edit == NULL || fact == NULL ||
        fact->shift_offset < 0 ||
        (size_t)edit->tx_idx >= transcripts->transcript_count) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    /* hgvs_protein first reuses the original pre-consequence `coding` fact.
     * It then asks genomic2pep() for the complete uploaded VariationFeature
     * after HGVS-only placement, not for the minimized semantic edit printed
     * by HGVSc. A leading/trailing UTR, intron, or transcript-external mapper
     * Gap consequently suppresses HGVSp even when the shifted edit itself can
     * be replayed against CDS. */
    if (!duckvep_project_feature_has_coding_precondition_unshifted(
            transcripts, exons, edit->tx_idx, &edit->event)) {
        return DUCKVEP_HGVS_OK;
    }
    *defined_out = duckvep_project_complete_feature_translation_bounds(
        transcripts, exons, edit->tx_idx, &edit->event,
        (uint32_t)fact->shift_offset, &first_cds, &last_cds);
    return DUCKVEP_HGVS_OK;
}

duckvep_hgvs_status_t duckvep_hgvs_shifted_cds_edit_build(
    const duckvep_transcript_model_t       *transcripts,
    const duckvep_exon_model_t             *exons,
    const duckvep_sequence_pool_t          *seq,
    const duckvep_hgvs_reference_window_t  *reference,
    const duckvep_transcript_edit_t        *edit,
    const duckvep_hgvs_dna_fact_t          *fact,
    uint8_t                                *allele_scratch,
    size_t                                  allele_scratch_cap,
    size_t                                 *required_out,
    duckvep_event_t                        *event_out,
    duckvep_haplotype_edit_t               *out) {

    duckvep_event_t event;
    duckvep_prepared_cds_allele_t allele;
    duckvep_hgvs_status_t status;
    duckvep_haplotype_edit_t result_edit;
    duckvep_cds_edit_status_t edit_status;
    size_t allele_length;
    size_t required;
    int insertion;

    if (required_out != NULL) *required_out = 0u;
    if (event_out != NULL) memset(event_out, 0, sizeof *event_out);
    if (out != NULL) memset(out, 0, sizeof *out);
    if (transcripts == NULL || exons == NULL || seq == NULL || edit == NULL ||
        fact == NULL || required_out == NULL || event_out == NULL || out == NULL ||
        edit->tx_idx >= transcripts->transcript_count ||
        fact->transcript_strand != edit->transcript_strand ||
        fact->shift_offset < 0) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    insertion = fact->ref_length == 0u && fact->alt_length != 0u;
    if (!insertion && !(fact->ref_length != 0u && fact->alt_length == 0u)) {
        return DUCKVEP_HGVS_UNSUPPORTED_EDIT;
    }
    if (insertion) {
        if (fact->shape != (uint8_t)DUCKVEP_HGVS_DNA_INSERTION &&
            fact->shape != (uint8_t)DUCKVEP_HGVS_DNA_DUPLICATION) {
            return DUCKVEP_HGVS_UNSUPPORTED_EDIT;
        }
        allele_length = (size_t)fact->alt_length;
    } else {
        if (fact->shape != (uint8_t)DUCKVEP_HGVS_DNA_DELETION) {
            return DUCKVEP_HGVS_UNSUPPORTED_EDIT;
        }
        allele_length = (size_t)fact->ref_length;
    }
    required = allele_length + (insertion ? 1u : 0u);
    *required_out = required;
    if (allele_scratch == NULL || allele_scratch_cap < required) {
        return DUCKVEP_HGVS_BUFFER_TOO_SMALL;
    }
    status = hgvs_materialize_genomic_allele(
        fact, insertion, allele_scratch, allele_length);
    if (status != DUCKVEP_HGVS_OK) return status;

    status = hgvs_shifted_event_build(edit, fact, &event);
    if (status != DUCKVEP_HGVS_OK) return status;

    memset(&allele, 0, sizeof allele);
    allele.event = &event;
    allele.variant_strand = (int8_t)1;
    if (insertion) {
        uint32_t anchor_position1;

        anchor_position1 = event.anchor_side ==
                (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT
            ? fact->placed_insertion_boundary0 :
              fact->placed_insertion_boundary0 + 1u;
        status = hgvs_reference_base(
            reference, anchor_position1, allele_scratch + allele_length);
        if (status != DUCKVEP_HGVS_OK) return status;
        allele.alt = allele_scratch;
        allele.alt_length = fact->alt_length;
        allele.anchor_ref = allele_scratch + allele_length;
    } else {
        allele.ref = allele_scratch;
        allele.ref_length = fact->ref_length;
    }
    memset(&result_edit, 0, sizeof result_edit);
    edit_status = duckvep_cds_edit_build_prepared_allele(
        transcripts, exons, seq, edit->tx_idx, edit->transcript_strand,
        &allele, UINT32_MAX, &result_edit);
    if (!insertion &&
        (edit_status == DUCKVEP_CDS_EDIT_OUT_OF_CDS ||
         edit_status == DUCKVEP_CDS_EDIT_NON_CONTIGUOUS)) {
        /*
         * HGVS-only 3-prime placement can move both ends of a deletion into
         * CDS while leaving one or more introns inside its genomic span.
         * TranscriptVariationAllele::_get_alternate_cds then replaces the
         * outer mapped CDS range, exactly like the independent consequence
         * path. Reuse that compatibility helper so HGVSp does not invent a
         * second mapper-gap interpretation.
         */
        edit_status = duckvep_compat_vep116_outer_cds_edit_build(
            transcripts, exons, seq, edit->tx_idx,
            edit->transcript_strand, &event, NULL, 0u,
            edit->transcript_strand, &result_edit);
    }
    if (!insertion && edit_status == DUCKVEP_CDS_EDIT_REF_MISMATCH) {
        uint32_t cds_start1;
        uint32_t cds_end1;
        uint64_t cds_offset;
        size_t cds_length;

        /* HGVS-only 3-prime placement mutates the TVA coordinates before
         * _get_alternate_cds applies a deletion. On VEP's overlapping short-
         * region reference slices, the rotated original deletion string need
         * not equal the prepared CDS at that nonlocal copy; VEP still removes
         * the transcript bases at the shifted coordinates. The uploaded REF
         * was validated at its original genomic span before this helper. */
        if (seq->cds_bytes == NULL || seq->cds_offset == NULL ||
            seq->cds_length == NULL ||
            edit->tx_idx >= seq->transcript_count ||
            !duckvep_project_event_to_cds(
                transcripts, exons, edit->tx_idx, &event,
                &cds_start1, &cds_end1) ||
            cds_end1 < cds_start1 ||
            (uint64_t)cds_end1 - (uint64_t)cds_start1 + 1u !=
                (uint64_t)fact->ref_length) {
            return DUCKVEP_HGVS_REFERENCE_MISMATCH;
        }
        cds_offset = seq->cds_offset[edit->tx_idx];
        cds_length = (size_t)seq->cds_length[edit->tx_idx];
        if (cds_offset > (uint64_t)seq->cds_bytes_len ||
            (uint64_t)cds_length >
                (uint64_t)seq->cds_bytes_len - cds_offset ||
            cds_start1 == 0u || (uint64_t)cds_end1 > cds_length) {
            return DUCKVEP_HGVS_INVALID_PROJECTION;
        }
        allele.ref = seq->cds_bytes + (size_t)cds_offset +
            (size_t)cds_start1 - 1u;
        allele.variant_strand = edit->transcript_strand;
        edit_status = duckvep_cds_edit_build_prepared_allele(
            transcripts, exons, seq, edit->tx_idx, edit->transcript_strand,
            &allele, UINT32_MAX, &result_edit);
    }
    status = hgvs_status_from_cds_edit(edit_status);
    if (status != DUCKVEP_HGVS_OK) return status;
    *event_out = event;
    *out = result_edit;
    return DUCKVEP_HGVS_OK;
}

duckvep_hgvs_status_t duckvep_hgvs_dna_base(
    const duckvep_hgvs_dna_fact_t *fact,
    int                             alternate,
    size_t                          index,
    uint8_t                        *base_out) {

    const uint8_t *allele;
    size_t length;
    if (base_out == NULL) return DUCKVEP_HGVS_INVALID_ARG;
    *base_out = 0u;
    if (fact == NULL ||
        (fact->transcript_strand != (int8_t)1 &&
         fact->transcript_strand != (int8_t)-1)) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    allele = alternate ? fact->alt : fact->ref;
    length = alternate ? (size_t)fact->alt_length
                       : (size_t)fact->ref_length;
    if (allele == NULL || index >= length) return DUCKVEP_HGVS_OUT_OF_RANGE;
    if (fact->shift_offset < 0) return DUCKVEP_HGVS_INVALID_ARG;
    return hgvs_oriented_allele_base(
        allele, length, fact->transcript_strand,
        (size_t)fact->shift_offset, index, base_out);
}

typedef struct hgvs_writer {
    char *buffer;
    size_t capacity;
    size_t required;
} hgvs_writer_t;

static int hgvs_writer_char(hgvs_writer_t *writer, char value) {
    if (writer == NULL || writer->required == SIZE_MAX) return 0;
    if (writer->buffer != NULL && writer->capacity > 0u &&
        writer->required < writer->capacity - 1u) {
        writer->buffer[writer->required] = value;
    }
    writer->required++;
    return 1;
}

static int hgvs_writer_literal(hgvs_writer_t *writer, const char *value) {
    size_t i;
    if (value == NULL) return 0;
    for (i = 0u; value[i] != '\0'; i++) {
        if (!hgvs_writer_char(writer, value[i])) return 0;
    }
    return 1;
}

static int hgvs_writer_uint64(hgvs_writer_t *writer, uint64_t value) {
    char digits[20];
    size_t count = 0u;
    do {
        digits[count++] = (char)('0' + (value % 10u));
        value /= 10u;
    } while (value != 0u);
    while (count > 0u) {
        if (!hgvs_writer_char(writer, digits[--count])) return 0;
    }
    return 1;
}

static int hgvs_writer_int64(hgvs_writer_t *writer, int64_t value,
                             int positive_sign) {
    uint64_t magnitude;
    if (value < 0) {
        if (!hgvs_writer_char(writer, '-')) return 0;
        magnitude = (uint64_t)(-(value + 1)) + 1u;
    } else {
        if (positive_sign && !hgvs_writer_char(writer, '+')) return 0;
        magnitude = (uint64_t)value;
    }
    return hgvs_writer_uint64(writer, magnitude);
}

static int hgvs_writer_coordinate(
    hgvs_writer_t                    *writer,
    const duckvep_hgvs_coordinate_t *coordinate) {

    int offset_consumed = 0;

    if (writer == NULL || coordinate == NULL) return 0;
    switch ((duckvep_hgvs_coordinate_kind_t)coordinate->kind) {
        case DUCKVEP_HGVS_COORDINATE_N:
            if (coordinate->base <= 0 ||
                !hgvs_writer_uint64(writer, (uint64_t)coordinate->base)) {
                return 0;
            }
            break;
        case DUCKVEP_HGVS_COORDINATE_C:
            if (coordinate->base == 0 ||
                !hgvs_writer_int64(writer, coordinate->base, 0)) {
                return 0;
            }
            break;
        case DUCKVEP_HGVS_COORDINATE_C_STAR:
            if (coordinate->base < 0 || !hgvs_writer_char(writer, '*')) {
                return 0;
            }
            if (coordinate->base > 0) {
                if (!hgvs_writer_uint64(writer, (uint64_t)coordinate->base)) {
                    return 0;
                }
            } else {
                if (coordinate->intron_offset == 0 ||
                    !hgvs_writer_int64(
                        writer, (int64_t)coordinate->intron_offset, 0)) {
                    return 0;
                }
                offset_consumed = 1;
            }
            break;
        default:
            return 0;
    }
    if (!offset_consumed && coordinate->intron_offset != 0 &&
        !hgvs_writer_int64(writer, (int64_t)coordinate->intron_offset, 1)) {
        return 0;
    }
    return 1;
}

static duckvep_hgvs_status_t hgvs_writer_finish(
    hgvs_writer_t *writer,
    size_t        *required_out) {

    if (writer == NULL || required_out == NULL) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    if (writer->buffer != NULL && writer->capacity > 0u) {
        size_t terminator = writer->required < writer->capacity
            ? writer->required : writer->capacity - 1u;
        writer->buffer[terminator] = '\0';
    }
    *required_out = writer->required;
    if (writer->buffer == NULL || writer->required >= writer->capacity) {
        return DUCKVEP_HGVS_BUFFER_TOO_SMALL;
    }
    return DUCKVEP_HGVS_OK;
}

static int hgvs_coordinate_equal(
    const duckvep_hgvs_coordinate_t *left,
    const duckvep_hgvs_coordinate_t *right) {

    return left != NULL && right != NULL && left->kind == right->kind &&
           left->base == right->base &&
           left->intron_offset == right->intron_offset;
}

static duckvep_hgvs_status_t hgvs_writer_allele(
    hgvs_writer_t                  *writer,
    const duckvep_hgvs_dna_fact_t *fact,
    int                             alternate) {

    size_t length = alternate ? (size_t)fact->alt_length
                              : (size_t)fact->ref_length;
    size_t i;
    for (i = 0u; i < length; i++) {
        uint8_t base;
        duckvep_hgvs_status_t status = duckvep_hgvs_dna_base(
            fact, alternate, i, &base);
        if (status != DUCKVEP_HGVS_OK) return status;
        if (!hgvs_writer_char(writer, (char)base)) {
            return DUCKVEP_HGVS_OUT_OF_RANGE;
        }
    }
    return DUCKVEP_HGVS_OK;
}

duckvep_hgvs_status_t duckvep_hgvs_dna_render_basic(
    const duckvep_hgvs_dna_fact_t *fact,
    char                           *buffer,
    size_t                          capacity,
    size_t                         *required_out) {

    hgvs_writer_t writer;
    duckvep_hgvs_status_t status;
    int same_coordinate;

    if (required_out != NULL) *required_out = 0u;
    if (buffer != NULL && capacity > 0u) buffer[0] = '\0';
    if (fact == NULL || required_out == NULL ||
        (buffer == NULL && capacity != 0u) || fact->shift_offset < 0 ||
        (fact->numbering != (uint8_t)DUCKVEP_HGVS_NUMBERING_C &&
         fact->numbering != (uint8_t)DUCKVEP_HGVS_NUMBERING_N)) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    memset(&writer, 0, sizeof writer);
    writer.buffer = buffer;
    writer.capacity = capacity;
    if (!hgvs_writer_char(&writer,
            fact->numbering == (uint8_t)DUCKVEP_HGVS_NUMBERING_C ? 'c' : 'n') ||
        !hgvs_writer_char(&writer, '.')) {
        return DUCKVEP_HGVS_OUT_OF_RANGE;
    }
    same_coordinate = hgvs_coordinate_equal(&fact->first, &fact->last);

    switch ((duckvep_hgvs_dna_shape_t)fact->shape) {
        case DUCKVEP_HGVS_DNA_SUBSTITUTION:
            if (!same_coordinate || fact->ref_length != 1u ||
                fact->alt_length != 1u ||
                !hgvs_writer_coordinate(&writer, &fact->first)) {
                return DUCKVEP_HGVS_INVALID_PROJECTION;
            }
            status = hgvs_writer_allele(&writer, fact, 0);
            if (status != DUCKVEP_HGVS_OK) return status;
            if (!hgvs_writer_char(&writer, '>')) {
                return DUCKVEP_HGVS_OUT_OF_RANGE;
            }
            status = hgvs_writer_allele(&writer, fact, 1);
            if (status != DUCKVEP_HGVS_OK) return status;
            break;
        case DUCKVEP_HGVS_DNA_DELETION:
        case DUCKVEP_HGVS_DNA_REPLACEMENT:
            if (!hgvs_writer_coordinate(&writer, &fact->first) ||
                (!same_coordinate &&
                 (!hgvs_writer_char(&writer, '_') ||
                  !hgvs_writer_coordinate(&writer, &fact->last)))) {
                return DUCKVEP_HGVS_INVALID_PROJECTION;
            }
            if (fact->shape == (uint8_t)DUCKVEP_HGVS_DNA_DELETION) {
                if (fact->ref_length == 0u || fact->alt_length != 0u ||
                    !hgvs_writer_literal(&writer, "del")) {
                    return DUCKVEP_HGVS_INVALID_ARG;
                }
            } else {
                if (fact->ref_length == 0u || fact->alt_length == 0u ||
                    !hgvs_writer_literal(&writer, "delins")) {
                    return DUCKVEP_HGVS_INVALID_ARG;
                }
                status = hgvs_writer_allele(&writer, fact, 1);
                if (status != DUCKVEP_HGVS_OK) return status;
            }
            break;
        case DUCKVEP_HGVS_DNA_INSERTION:
            if (same_coordinate || fact->ref_length != 0u ||
                fact->alt_length == 0u ||
                !hgvs_writer_coordinate(&writer, &fact->first) ||
                !hgvs_writer_char(&writer, '_') ||
                !hgvs_writer_coordinate(&writer, &fact->last) ||
                !hgvs_writer_literal(&writer, "ins")) {
                return DUCKVEP_HGVS_INVALID_PROJECTION;
            }
            status = hgvs_writer_allele(&writer, fact, 1);
            if (status != DUCKVEP_HGVS_OK) return status;
            break;
        case DUCKVEP_HGVS_DNA_INVERSION:
        case DUCKVEP_HGVS_DNA_DUPLICATION:
            if (!hgvs_writer_coordinate(&writer, &fact->first) ||
                (!same_coordinate &&
                 (!hgvs_writer_char(&writer, '_') ||
                  !hgvs_writer_coordinate(&writer, &fact->last))) ||
                !hgvs_writer_literal(
                    &writer,
                    fact->shape == (uint8_t)DUCKVEP_HGVS_DNA_INVERSION
                        ? "inv" : "dup")) {
                return DUCKVEP_HGVS_INVALID_PROJECTION;
            }
            break;
        case DUCKVEP_HGVS_DNA_REPEAT:
            if (fact->repeat_count < 2u ||
                !hgvs_writer_coordinate(&writer, &fact->first) ||
                (!same_coordinate &&
                 (!hgvs_writer_char(&writer, '_') ||
                  !hgvs_writer_coordinate(&writer, &fact->last))) ||
                !hgvs_writer_char(&writer, '[') ||
                !hgvs_writer_uint64(&writer, fact->repeat_count) ||
                !hgvs_writer_char(&writer, ']')) {
                return DUCKVEP_HGVS_INVALID_PROJECTION;
            }
            break;
        default:
            return DUCKVEP_HGVS_INVALID_ARG;
    }
    return hgvs_writer_finish(&writer, required_out);
}

static int hgvs_protein_full_window(
    const duckvep_coding_context_t  *context,
    duckvep_coding_peptide_window_t *window) {

    if (context == NULL || window == NULL || context->ref_cds == NULL ||
        context->ref_peptide_len > UINT32_MAX ||
        context->alt_peptide_len > UINT32_MAX ||
        (!context->virtual_single_edit &&
         (context->ref_peptide == NULL || context->alt_peptide == NULL))) {
        return 0;
    }
    memset(window, 0, sizeof *window);
    window->ref_length = context->ref_peptide_len;
    window->alt_length = context->alt_peptide_len;
    window->ref_whole_length = context->ref_peptide_len;
    window->alt_whole_length = context->alt_peptide_len;
    window->reference_span_length = context->ref_peptide_len;
    return 1;
}

static int hgvs_protein_terminal_partial_insertion_window(
    const duckvep_coding_context_t  *context,
    duckvep_coding_peptide_window_t *window) {

    size_t start0;
    size_t nt_offset;

    if (window == NULL || context == NULL ||
        !duckvep_compat_enabled(
            (duckvep_compat_profile_t)context->compatibility_profile,
            DUCKVEP_COMPAT_HGVS_TERMINAL_PARTIAL_INSERTION) ||
        !duckvep_coding_context_is_terminal_partial_insertion(context)) {
        return 0;
    }
    start0 = (size_t)context->single_edit_cds_start - 1u;
    nt_offset = (start0 / 3u) * 3u;
    if (nt_offset > context->ref_cds_len ||
        (size_t)context->single_edit_alt_len >
            SIZE_MAX - (context->ref_cds_len - nt_offset)) {
        return 0;
    }
    memset(window, 0, sizeof *window);
    window->peptide_offset = nt_offset / 3u;
    window->ref_nt_length = context->ref_cds_len - nt_offset;
    window->alt_nt_length =
        window->ref_nt_length + (size_t)context->single_edit_alt_len;
    window->ref_whole_length = window->ref_nt_length / 3u;
    window->alt_whole_length = window->alt_nt_length / 3u;
    window->ref_partial_x =
        (uint8_t)((window->ref_nt_length % 3u) != 0u);
    window->alt_partial_x =
        (uint8_t)((window->alt_nt_length % 3u) != 0u);
    if (window->alt_partial_x && window->alt_whole_length == 1u &&
        duckvep_coding_context_peptide_base(
            context, 1, window->peptide_offset) == (uint8_t)'*') {
        window->alt_partial_x = 0u;
    }
    window->ref_length =
        window->ref_whole_length + (size_t)window->ref_partial_x;
    window->alt_length =
        window->alt_whole_length + (size_t)window->alt_partial_x;
    window->reference_span_length = window->ref_length;
    return 1;
}

static int hgvs_protein_uses_terminal_partial_insertion_view(
    const duckvep_coding_context_t *context) {

    return context != NULL &&
        duckvep_compat_enabled(
            (duckvep_compat_profile_t)context->compatibility_profile,
            DUCKVEP_COMPAT_HGVS_TERMINAL_PARTIAL_INSERTION) &&
        duckvep_coding_context_is_terminal_partial_insertion(context);
}

static uint8_t hgvs_protein_window_base(
    const duckvep_coding_context_t        *context,
    const duckvep_coding_peptide_window_t *window,
    int                                    alternate,
    size_t                                 index) {

    if (hgvs_protein_uses_terminal_partial_insertion_view(context) &&
        window != NULL && window->ref_nt_length != 0u) {
        size_t whole_length = alternate
            ? window->alt_whole_length : window->ref_whole_length;
        size_t length = alternate ? window->alt_length : window->ref_length;
        uint8_t partial_x = alternate
            ? window->alt_partial_x : window->ref_partial_x;

        if (index >= length) return 0u;
        if (index == whole_length && partial_x) return (uint8_t)'X';
        return duckvep_coding_context_peptide_base(
            context, alternate, window->peptide_offset + index);
    }
    return duckvep_coding_context_peptide_window_base(
        context, window, alternate, index);
}

static size_t hgvs_protein_reference_length(
    const duckvep_coding_context_t *context) {

    size_t length;
    if (context == NULL) return 0u;
    length = context->ref_peptide_len;
    if (length != 0u && duckvep_coding_context_peptide_base(
            context, 0, length - 1u) == (uint8_t)'*') {
        length--;
    }
    return length;
}

/* VEP 116's _trim_incomplete_codon() assigns keep_length in its equality
 * guard.  An alternate CDS of at least three bases therefore survives in
 * full, including a one- or two-base terminal remainder; a one- or two-base
 * alternate CDS is trimmed to empty.  Apply that executable state before the
 * transcript 3-prime UTR is appended. */
static size_t hgvs_protein_alt_cds_length(
    const duckvep_coding_context_t *context) {

    if (context == NULL) return 0u;
    if (duckvep_compat_enabled(
            (duckvep_compat_profile_t)context->compatibility_profile,
            DUCKVEP_COMPAT_HGVS_INCOMPLETE_CODON_ASSIGNMENT)) {
        return context->alt_cds_len < 3u ? 0u : context->alt_cds_len;
    }
    return context->alt_cds_len - (context->alt_cds_len % 3u);
}

static duckvep_hgvs_status_t hgvs_protein_extended_alt_base(
    const duckvep_coding_context_t *context,
    size_t                          position0,
    uint8_t                        *residue_out) {

    size_t edited_cds_length;
    size_t nucleotide_position;
    size_t total_length;
    char codon[4];
    size_t i;
    int has_n = 0;

    if (residue_out == NULL) return DUCKVEP_HGVS_INVALID_ARG;
    *residue_out = 0u;
    if (context == NULL || !duckvep_codon_table_supported(
            (duckvep_codon_table_t)context->codon_table)) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    edited_cds_length = hgvs_protein_alt_cds_length(context);
    if (context->post_cds_complete == 0u ||
        (context->post_cds_length != 0u &&
         context->post_cds_bases == NULL) ||
        edited_cds_length > SIZE_MAX - context->post_cds_length ||
        position0 > SIZE_MAX / 3u) {
        return DUCKVEP_HGVS_MISSING_PEPTIDE;
    }
    total_length = edited_cds_length + context->post_cds_length;
    nucleotide_position = position0 * 3u;
    if (nucleotide_position > total_length ||
        total_length - nucleotide_position < 3u) {
        return DUCKVEP_HGVS_OUT_OF_RANGE;
    }
    for (i = 0u; i < 3u; i++) {
        size_t index = nucleotide_position + i;
        char base = index < edited_cds_length
            ? duckvep_coding_context_cds_base(context, 1, index)
            : duckvep_dna_normalize(
                  (char)context->post_cds_bases[
                      index - edited_cds_length], 1);
        if (base == '\0') return DUCKVEP_HGVS_INVALID_ALLELE;
        if (base == 'N') has_n = 1;
        codon[i] = base;
    }
    codon[3] = '\0';
    *residue_out = has_n ? (uint8_t)'X'
        : (uint8_t)duckvep_translate_codon(
              codon, (duckvep_codon_table_t)context->codon_table);
    return DUCKVEP_HGVS_OK;
}

static size_t hgvs_protein_extended_alt_length(
    const duckvep_coding_context_t *context) {

    size_t edited_cds_length;
    if (context == NULL) return 0u;
    edited_cds_length = hgvs_protein_alt_cds_length(context);
    if (context->post_cds_complete == 0u ||
        edited_cds_length > SIZE_MAX - context->post_cds_length) {
        return context->alt_cds_len / 3u;
    }
    return (edited_cds_length + context->post_cds_length) / 3u;
}

static void hgvs_protein_stop_distance(
    const duckvep_coding_context_t *context,
    uint32_t                        variant_position1,
    int                             frameshift,
    uint32_t                       *distance_out,
    uint8_t                        *known_out) {

    size_t edited_cds_length;
    size_t stop_position0;
    size_t reference_length;
    duckvep_coding_context_t vep_context;
    int found;

    if (distance_out != NULL) *distance_out = 0u;
    if (known_out != NULL) *known_out = 0u;
    if (context == NULL || distance_out == NULL || known_out == NULL ||
        variant_position1 == 0u) {
        return;
    }
    edited_cds_length = hgvs_protein_alt_cds_length(context);
    /* VEP 116's _stop_loss_extra_AA() calls $alt_cds->translate() without
     * passing the transcript's codon table. BioPerl therefore uses standard
     * table 1 for this late stop search even when the ordinary peptide and
     * consequence path used (for example) vertebrate mitochondrial table 2.
     * Preserve that executable inconsistency only in this formatter query. */
    vep_context = *context;
    vep_context.codon_table = (uint8_t)duckvep_compat_late_stop_codon_table(
        (duckvep_compat_profile_t)context->compatibility_profile,
        (duckvep_codon_table_t)context->codon_table);
    /* This shortcut was proved with context->codon_table. Reusing its first
     * stop fact after switching tables can skip an earlier mitochondrial TGA
     * that VEP's full table-1 translation sees before the edited codon. */
    vep_context.ref_first_stop_known = 0u;
    vep_context.ref_first_stop_position1 = 0u;
    if (context->post_cds_complete == 0u ||
        (context->post_cds_length != 0u &&
         context->post_cds_bases == NULL) ||
        duckvep_coding_context_first_alt_stop(
            &vep_context, edited_cds_length, context->post_cds_bases,
            context->post_cds_length, &stop_position0, &found) !=
                DUCKVEP_CODING_CONTEXT_OK ||
        !found) {
        return;
    }
    reference_length = hgvs_protein_reference_length(context);
    {
        int64_t distance = frameshift
            ? (int64_t)stop_position0 + 1 -
                  ((int64_t)variant_position1 - 1)
            : (int64_t)stop_position0 - (int64_t)reference_length;
        if (distance > 0 && (uint64_t)distance <= UINT32_MAX) {
            *distance_out = (uint32_t)distance;
            *known_out = 1u;
        }
    }
}

duckvep_hgvs_status_t
duckvep_hgvs_protein_frameshift_termination_replay(
    const duckvep_coding_context_t *late_context,
    duckvep_hgvs_protein_fact_t    *fact) {

    if (late_context == NULL || fact == NULL ||
        !duckvep_compat_profile_valid(
            (duckvep_compat_profile_t)
                late_context->compatibility_profile) ||
        !duckvep_compat_profile_valid(
            (duckvep_compat_profile_t)fact->compatibility_profile) ||
        fact->compatibility_profile !=
            late_context->compatibility_profile ||
        (fact->context != NULL &&
         (!duckvep_compat_profile_valid(
              (duckvep_compat_profile_t)
                  fact->context->compatibility_profile) ||
          fact->context->compatibility_profile !=
              fact->compatibility_profile)) ||
        fact->shape != (uint8_t)DUCKVEP_HGVS_PROTEIN_FRAMESHIFT ||
        fact->first_position1 == 0u) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    fact->termination_distance = 0u;
    fact->termination_known = 0u;
    /* _get_hgvs_protein_format passes start - 1 to
     * _stop_loss_extra_AA(), whose first statement rejects zero. */
    if (fact->first_position1 == 1u) return DUCKVEP_HGVS_OK;
    hgvs_protein_stop_distance(
        late_context, fact->first_position1, 1,
        &fact->termination_distance, &fact->termination_known);
    return DUCKVEP_HGVS_OK;
}

static duckvep_hgvs_status_t hgvs_protein_fact_set_residues(
    duckvep_hgvs_protein_fact_t *fact) {

    uint8_t residue;
    duckvep_hgvs_status_t status;

    if (fact == NULL ||
        ((fact->first_position1 == 0u || fact->last_position1 == 0u) &&
         !((duckvep_hgvs_protein_shape_t)fact->shape ==
               DUCKVEP_HGVS_PROTEIN_INSERTION &&
           fact->first_position1 == 1u && fact->last_position1 == 0u))) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    if (fact->ref_length != 0u) {
        status = duckvep_hgvs_protein_base(fact, 0, 0u, &residue);
        if (status != DUCKVEP_HGVS_OK) return status;
        fact->reference_first = residue;
        status = duckvep_hgvs_protein_base(
            fact, 0, fact->ref_length - 1u, &residue);
        if (status != DUCKVEP_HGVS_OK) return status;
        fact->reference_last = residue;
    }
    if (fact->alt_length != 0u) {
        status = duckvep_hgvs_protein_base(fact, 1, 0u, &residue);
        if (status != DUCKVEP_HGVS_OK) return status;
        fact->alternate_first = residue;
    }
    return DUCKVEP_HGVS_OK;
}

duckvep_hgvs_status_t duckvep_hgvs_protein_base(
    const duckvep_hgvs_protein_fact_t *fact,
    int                                alternate,
    size_t                             index,
    uint8_t                           *residue_out) {

    size_t offset;
    size_t length;
    size_t rotation;
    size_t local_index;

    if (residue_out == NULL) return DUCKVEP_HGVS_INVALID_ARG;
    *residue_out = 0u;
    if (fact == NULL || fact->context == NULL ||
        !duckvep_compat_profile_valid(
            (duckvep_compat_profile_t)fact->compatibility_profile) ||
        !duckvep_compat_profile_valid(
            (duckvep_compat_profile_t)
                fact->context->compatibility_profile) ||
        fact->compatibility_profile !=
            fact->context->compatibility_profile) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    offset = alternate ? fact->alt_offset : fact->ref_offset;
    length = alternate ? fact->alt_length : fact->ref_length;
    rotation = alternate ? fact->alt_rotation : fact->ref_rotation;
    if (index >= length) return DUCKVEP_HGVS_OUT_OF_RANGE;
    local_index = offset + (index + rotation % length) % length;
    *residue_out = hgvs_protein_window_base(
        fact->context, &fact->window, alternate, local_index);
    return *residue_out == 0u ? DUCKVEP_HGVS_MISSING_PEPTIDE
                              : DUCKVEP_HGVS_OK;
}

static duckvep_hgvs_status_t hgvs_protein_frameshift_fact(
    const duckvep_coding_context_t *context,
    uint32_t                        start_position1,
    duckvep_hgvs_protein_fact_t    *fact) {

    size_t ref_length;
    size_t alt_length;
    size_t position0;
    uint8_t reference;
    uint8_t alternate;

    if (context == NULL || fact == NULL || start_position1 == 0u) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    ref_length = hgvs_protein_reference_length(context);
    alt_length = hgvs_protein_extended_alt_length(context);
    position0 = (size_t)start_position1 - 1u;
    if (position0 >= alt_length) {
        if (position0 > ref_length) return DUCKVEP_HGVS_OUT_OF_RANGE;
        reference = position0 < ref_length
            ? duckvep_coding_context_peptide_base(context, 0, position0)
            : (uint8_t)'*';
        if (reference == 0u) return DUCKVEP_HGVS_MISSING_PEPTIDE;
        fact->shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_DELETION;
        fact->first_position1 = start_position1;
        fact->last_position1 = start_position1;
        fact->reference_first = reference;
        fact->reference_last = reference;
        return DUCKVEP_HGVS_OK;
    }
    while (position0 < alt_length) {
        duckvep_hgvs_status_t status;
        reference = position0 < ref_length
            ? duckvep_coding_context_peptide_base(context, 0, position0)
            : position0 == ref_length ? (uint8_t)'*' : 0u;
        status = hgvs_protein_extended_alt_base(
            context, position0, &alternate);
        if (status != DUCKVEP_HGVS_OK) return status;
        if (reference == (uint8_t)'*' && alternate == (uint8_t)'*') {
            fact->shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_EQUAL;
            fact->first_position1 = (uint32_t)position0 + 1u;
            fact->last_position1 = fact->first_position1;
            fact->reference_first = reference;
            fact->reference_last = reference;
            fact->alternate_first = alternate;
            return DUCKVEP_HGVS_OK;
        }
        if (reference == 0u || reference != alternate) break;
        position0++;
        if (position0 >= UINT32_MAX) return DUCKVEP_HGVS_OUT_OF_RANGE;
    }
    if (position0 >= alt_length || reference == 0u) {
        return DUCKVEP_HGVS_UNSUPPORTED_PROTEIN;
    }
    fact->first_position1 = (uint32_t)position0 + 1u;
    fact->last_position1 = fact->first_position1;
    fact->reference_first = reference;
    fact->reference_last = reference;
    fact->alternate_first = alternate;
    fact->shape = alternate == (uint8_t)'*'
        ? (uint8_t)DUCKVEP_HGVS_PROTEIN_SUBSTITUTION
        : (uint8_t)DUCKVEP_HGVS_PROTEIN_FRAMESHIFT;
    if (fact->shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_FRAMESHIFT) {
        hgvs_protein_stop_distance(
            context, fact->first_position1, 1,
            &fact->termination_distance, &fact->termination_known);
    }
    return DUCKVEP_HGVS_OK;
}

static duckvep_hgvs_status_t hgvs_protein_shift_simple(
    duckvep_hgvs_protein_fact_t *fact) {

    const duckvep_coding_context_t *context;
    size_t pattern_length;
    size_t post_start0;
    size_t available;
    size_t limit;
    size_t shift = 0u;
    int alternate;

    if (fact == NULL || fact->context == NULL) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    alternate = fact->shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_INSERTION;
    if (!alternate &&
        fact->shape != (uint8_t)DUCKVEP_HGVS_PROTEIN_DELETION) {
        return DUCKVEP_HGVS_OK;
    }
    pattern_length = alternate ? fact->alt_length : fact->ref_length;
    if (pattern_length == 0u) {
        return DUCKVEP_HGVS_INVALID_PROJECTION;
    }
    if (fact->last_position1 == 0u && !alternate) {
        return DUCKVEP_HGVS_INVALID_PROJECTION;
    }
    context = fact->context;
    post_start0 = (size_t)fact->last_position1;
    /* TranscriptVariationAllele::_get_surrounding_peptides() returns undef
     * when length(_peptide) <= post_pos.  Since post_pos is the one-based
     * residue after the changed peptide, VEP refuses even the otherwise
     * matching final residue as a 3-prime shift source.  Preserve that
     * executable endpoint test as well as excluding the synthetic stop. */
    available =
        post_start0 < SIZE_MAX &&
        hgvs_protein_reference_length(context) > post_start0 + 1u
            ? hgvs_protein_reference_length(context) - post_start0 : 0u;
    /* _shift_3prime() iterates only through
     * length(post_seq) - length(changed_peptide). When the complete changed
     * peptide is longer than the remaining reference peptide, Perl performs
     * no rotation at all even if their first residues match. */
    limit = available < pattern_length
        ? 0u : available - pattern_length + 1u;
    while (shift < limit) {
        uint8_t pattern;
        uint8_t following;
        duckvep_hgvs_status_t status = duckvep_hgvs_protein_base(
            fact, alternate, 0u, &pattern);
        if (status != DUCKVEP_HGVS_OK) return status;
        following = duckvep_coding_context_peptide_base(
            context, 0, post_start0 + shift);
        if (following == 0u) return DUCKVEP_HGVS_MISSING_PEPTIDE;
        if (pattern != following) break;
        if (fact->first_position1 == UINT32_MAX ||
            fact->last_position1 == UINT32_MAX) {
            return DUCKVEP_HGVS_OUT_OF_RANGE;
        }
        fact->first_position1++;
        fact->last_position1++;
        if (alternate) fact->alt_rotation++;
        else fact->ref_rotation++;
        shift++;
    }
    return DUCKVEP_HGVS_OK;
}

static duckvep_hgvs_status_t hgvs_protein_insertion_finish(
    duckvep_hgvs_protein_fact_t *fact) {

    const duckvep_coding_context_t *context;
    uint32_t source_start1;
    uint32_t source_end1;
    uint32_t low;
    size_t i;
    int duplicated = 1;

    if (fact == NULL || fact->context == NULL || fact->alt_length == 0u ||
        (fact->first_position1 == 0u && fact->last_position1 == 0u)) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    context = fact->context;
    if ((uint64_t)fact->first_position1 <= (uint64_t)fact->alt_length) {
        duplicated = 0;
        source_start1 = source_end1 = 0u;
    } else {
        source_end1 = fact->first_position1 - 1u;
        source_start1 = source_end1 - (uint32_t)fact->alt_length + 1u;
    }
    if (duplicated &&
        (size_t)source_end1 > hgvs_protein_reference_length(context)) {
        /* TranscriptVariationAllele::_peptide excludes the terminal stop.
         * Its duplication check therefore cannot use that synthetic residue
         * as the source of a protein duplication. */
        duplicated = 0;
    }
    for (i = 0u; duplicated && i < fact->alt_length; i++) {
        uint8_t inserted;
        uint8_t reference;
        duckvep_hgvs_status_t status = duckvep_hgvs_protein_base(
            fact, 1, i, &inserted);
        if (status != DUCKVEP_HGVS_OK) return status;
        reference = duckvep_coding_context_peptide_base(
            context, 0, (size_t)source_start1 - 1u + i);
        if (reference == 0u || inserted != reference) duplicated = 0;
    }
    if (duplicated) {
        fact->shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_DUPLICATION;
        fact->first_position1 = source_start1;
        fact->last_position1 = source_end1;
        fact->reference_first = duckvep_coding_context_peptide_base(
            context, 0, (size_t)source_start1 - 1u);
        fact->reference_last = duckvep_coding_context_peptide_base(
            context, 0, (size_t)source_end1 - 1u);
        return fact->reference_first == 0u || fact->reference_last == 0u
            ? DUCKVEP_HGVS_MISSING_PEPTIDE : DUCKVEP_HGVS_OK;
    }

    low = fact->first_position1 < fact->last_position1
        ? fact->first_position1 : fact->last_position1;
    if (low == 0u) {
        size_t reference_length = hgvs_protein_reference_length(context);
        uint8_t terminal_reference;

        /* Perl substr(_peptide, -1, 2) is reached for an insertion clipped to
         * start=1,end=0. It returns the final translated residue and VEP then
         * formats that same residue at positions zero and one. Preserve this
         * executable VEP-116 compatibility rule rather than imposing valid
         * HGVS coordinates. */
        if (!duckvep_compat_enabled(
                (duckvep_compat_profile_t)context->compatibility_profile,
                DUCKVEP_COMPAT_HGVS_NEGATIVE_SUBSTR) ||
            reference_length == 0u) {
            return DUCKVEP_HGVS_NOT_APPLICABLE;
        }
        terminal_reference = duckvep_coding_context_peptide_base(
            context, 0, reference_length - 1u);
        if (terminal_reference == 0u) return DUCKVEP_HGVS_MISSING_PEPTIDE;
        fact->first_position1 = 0u;
        fact->last_position1 = 1u;
        fact->reference_first = terminal_reference;
        fact->reference_last = terminal_reference;
        return DUCKVEP_HGVS_OK;
    }
    if ((size_t)low > hgvs_protein_reference_length(context)) {
        /* _get_surrounding_peptides() may append original_ref='*' and expose
         * the synthetic terminal stop as the second insertion flank. */
        return DUCKVEP_HGVS_NOT_APPLICABLE;
    }
    fact->first_position1 = low;
    fact->last_position1 = low + 1u;
    fact->reference_first = duckvep_coding_context_peptide_base(
        context, 0, (size_t)low - 1u);
    fact->reference_last = duckvep_coding_context_peptide_base(
        context, 0, (size_t)low);
    if (fact->reference_first == 0u) {
        return DUCKVEP_HGVS_MISSING_PEPTIDE;
    }
    /* At the peptide end, a real reference stop is retained in the coding
     * context but excluded from hgvs_protein_reference_length(). VEP exposes
     * it as the second insertion flank. Without that stop the insertion is
     * absent, not a peptide-data failure. */
    if ((size_t)low == hgvs_protein_reference_length(context)) {
        uint8_t original_reference_first =
            fact->window.ref_length == 0u ? 0u :
            hgvs_protein_window_base(
                context, &fact->window, 0, 0u);

        /* _clip_alleles caches the complete local reference peptide before
         * clipping. _get_surrounding_peptides() appends that cached peptide
         * to TranscriptVariation::_peptide only when it begins with '*'.
         * This is what makes the terminal stop available as the second flank,
         * including a nucleotide delins whose peptide-level edit clips to an
         * insertion. Length or nucleotide edit kind is not the authority. */
        if (fact->reference_last == 0u ||
            original_reference_first != (uint8_t)'*') {
            return DUCKVEP_HGVS_NOT_APPLICABLE;
        }
    }
    return fact->reference_last == 0u
        ? DUCKVEP_HGVS_MISSING_PEPTIDE : DUCKVEP_HGVS_OK;
}

duckvep_hgvs_status_t duckvep_hgvs_protein_fact_build(
    const duckvep_coding_context_t *context,
    const duckvep_sequence_delta_t *delta,
    duckvep_hgvs_protein_fact_t    *out) {

    duckvep_hgvs_protein_fact_t fact;
    size_t ref_length;
    size_t alt_length;
    size_t prefix = 0u;
    size_t suffix = 0u;
    uint64_t first64;
    uint64_t last64;
    duckvep_hgvs_status_t status;
    int stop_pair_early = 0;
    int xaa_as_ter;

    if (out == NULL) return DUCKVEP_HGVS_INVALID_ARG;
    memset(out, 0, sizeof *out);
    if (context == NULL || delta == NULL || delta->valid == 0u ||
        !duckvep_compat_profile_valid(
            (duckvep_compat_profile_t)context->compatibility_profile) ||
        context->ref_cds == NULL || context->ref_peptide_len == 0u ||
        (!context->virtual_single_edit &&
         (context->ref_peptide == NULL || context->alt_peptide == NULL))) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    if (context->ref_peptide_len > UINT32_MAX ||
        context->alt_peptide_len > UINT32_MAX) {
        return DUCKVEP_HGVS_OUT_OF_RANGE;
    }
    memset(&fact, 0, sizeof fact);
    fact.context = context;
    fact.compatibility_profile = context->compatibility_profile;
    xaa_as_ter = duckvep_compat_enabled(
        (duckvep_compat_profile_t)fact.compatibility_profile,
        DUCKVEP_COMPAT_HGVS_XAA_AS_TER);
    if (!duckvep_coding_context_peptide_window_open(
            context, &fact.window) &&
        !hgvs_protein_full_window(context, &fact.window)) {
        return DUCKVEP_HGVS_MISSING_PEPTIDE;
    }
    if (hgvs_protein_uses_terminal_partial_insertion_view(context) &&
        !hgvs_protein_terminal_partial_insertion_window(
            context, &fact.window)) {
        return DUCKVEP_HGVS_MISSING_PEPTIDE;
    }
    if (fact.window.ref_length == 0u && fact.window.alt_length == 0u) {
        /* VEP has no peptide allele to format (notably a terminal insertion
         * whose translated stop is retained), so hgvs_protein returns undef. */
        return DUCKVEP_HGVS_NOT_APPLICABLE;
    }
    first64 = (uint64_t)fact.window.peptide_offset + 1u;
    last64 = (uint64_t)fact.window.peptide_offset +
             (uint64_t)fact.window.ref_length;
    if (first64 > UINT32_MAX || last64 > UINT32_MAX) {
        return DUCKVEP_HGVS_OUT_OF_RANGE;
    }
    fact.first_position1 = (uint32_t)first64;
    fact.last_position1 = (uint32_t)last64;
    ref_length = fact.window.ref_length;
    alt_length = fact.window.alt_length;

    if (ref_length == alt_length && ref_length != 0u) {
        size_t i;
        int identical = 1;
        for (i = 0u; i < ref_length; i++) {
            uint8_t reference = hgvs_protein_window_base(
                context, &fact.window, 0, i);
            uint8_t alternate = hgvs_protein_window_base(
                context, &fact.window, 1, i);
            if (reference == 0u || alternate == 0u) {
                return DUCKVEP_HGVS_MISSING_PEPTIDE;
            }
            if (reference != alternate) {
                identical = 0;
                break;
            }
        }
        /* hgvs_protein skips _clip_alleles when the complete local peptides
         * are equal, so equality keeps translation_start rather than moving
         * one position past the matched window. */
        if (identical) {
            fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_EQUAL;
            fact.ref_length = ref_length;
            fact.alt_length = alt_length;
            fact.reference_first = hgvs_protein_window_base(
                context, &fact.window, 0, 0u);
            fact.reference_last = fact.reference_first;
            fact.alternate_first = fact.reference_first;
            fact.last_position1 = fact.first_position1;
            *out = fact;
            return DUCKVEP_HGVS_OK;
        }
    }

    while (prefix < ref_length && prefix < alt_length) {
        uint8_t reference = hgvs_protein_window_base(
            context, &fact.window, 0, prefix);
        uint8_t alternate = hgvs_protein_window_base(
            context, &fact.window, 1, prefix);
        if (reference == 0u || alternate == 0u) {
            return DUCKVEP_HGVS_MISSING_PEPTIDE;
        }
        /* _clip_alleles returns equality immediately when both local peptide
         * strings encounter a stop, even if later bytes differ. */
        if (reference == (uint8_t)'*' && alternate == (uint8_t)'*') {
            /* _clip_alleles returns before publishing its local strings and
             * coordinates, but _get_hgvs_protein_type still runs afterward.
             * Preserve the complete window here; a longer alternate peptide
             * is consequently a delins, not synonymous merely because both
             * strings begin with a stop. */
            stop_pair_early = 1;
            break;
        }
        if (reference != alternate) break;
        prefix++;
    }
    while (!stop_pair_early && suffix < ref_length - prefix &&
           suffix < alt_length - prefix) {
        uint8_t reference = hgvs_protein_window_base(
            context, &fact.window, 0, ref_length - 1u - suffix);
        uint8_t alternate = hgvs_protein_window_base(
            context, &fact.window, 1, alt_length - 1u - suffix);
        if (reference == 0u || alternate == 0u) {
            return DUCKVEP_HGVS_MISSING_PEPTIDE;
        }
        if (reference != alternate) break;
        suffix++;
    }
    if (stop_pair_early) {
        fact.ref_offset = 0u;
        fact.alt_offset = 0u;
        fact.ref_length = ref_length;
        fact.alt_length = alt_length;
    } else {
        fact.ref_offset = prefix;
        fact.alt_offset = prefix;
        fact.ref_length = ref_length - prefix - suffix;
        fact.alt_length = alt_length - prefix - suffix;
        if ((uint64_t)fact.first_position1 + prefix > UINT32_MAX ||
            suffix > (size_t)fact.last_position1) {
            return DUCKVEP_HGVS_OUT_OF_RANGE;
        }
        fact.first_position1 += (uint32_t)prefix;
        fact.last_position1 -= (uint32_t)suffix;
    }

    if (delta->frameshift) {
        uint32_t mapper_end = fact.last_position1;
        status = hgvs_protein_frameshift_fact(
            context, (uint32_t)first64, &fact);
        if (status != DUCKVEP_HGVS_OK) return status;
        if (delta->start_lost) {
            /* _get_fs_peptides replaces only start/ref/alt.  The mapper's
             * original translation_end remains in the HGVS hash before the
             * later start_lost override clears type and changes alt to '?'. */
            fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_START_LOST;
            fact.last_position1 = mapper_end;
            fact.ref_length = 1u;
            fact.reference_last = fact.reference_first;
            fact.start_lost_flanking = 1u;
        } else if (delta->stop_lost &&
                   fact.shape ==
                       (uint8_t)DUCKVEP_HGVS_PROTEIN_DELETION) {
            /* _get_hgvs_protein_format tests cached stop_lost plus a del/>
             * peptide type before its later fs branch. A frame-changing edit
             * that removes the terminal peptide is therefore rendered as a
             * stop extension, not as an ordinary deletion. */
            fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_EXTENSION;
            hgvs_protein_stop_distance(
                context, fact.first_position1, 0,
                &fact.termination_distance, &fact.termination_known);
        }
        *out = fact;
        return DUCKVEP_HGVS_OK;
    }

    if (delta->start_lost) {
        uint8_t reference;
        uint8_t reference_last;

        if (fact.ref_length == 0u) {
            /* _get_hgvs_peptides() treats this as a peptide insertion before
             * its final start_lost override. Rotate the changed peptide across
             * matching following reference residues exactly as VEP's
             * _check_peptides_post_var()/_shift_3prime() do. */
            fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_INSERTION;
            status = hgvs_protein_shift_simple(&fact);
            if (status != DUCKVEP_HGVS_OK) return status;
            uint32_t low = fact.first_position1 < fact.last_position1
                ? fact.first_position1 : fact.last_position1;
            size_t translation_length =
                hgvs_protein_reference_length(context);
            if (fact.alt_length == 0u || translation_length == 0u ||
                (low != 0u && (size_t)low >= translation_length)) {
                return DUCKVEP_HGVS_NOT_APPLICABLE;
            }
            if (low == 0u) {
                if (!duckvep_compat_enabled(
                        (duckvep_compat_profile_t)
                            context->compatibility_profile,
                        DUCKVEP_COMPAT_HGVS_NEGATIVE_SUBSTR)) {
                    return DUCKVEP_HGVS_NOT_APPLICABLE;
                }
                reference = duckvep_coding_context_peptide_base(
                    context, 0, translation_length - 1u);
                if (reference == 0u) {
                    return DUCKVEP_HGVS_MISSING_PEPTIDE;
                }
                fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_START_LOST;
                fact.reference_first = reference;
                fact.reference_last = reference;
                fact.ref_length = 1u;
                fact.start_lost_flanking = 1u;
                *out = fact;
                return DUCKVEP_HGVS_OK;
            }
            reference = duckvep_coding_context_peptide_base(
                context, 0, (size_t)low - 1u);
            reference_last = duckvep_coding_context_peptide_base(
                context, 0, (size_t)low);
            if (reference == 0u || reference_last == 0u) {
                return DUCKVEP_HGVS_MISSING_PEPTIDE;
            }
            fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_START_LOST;
            fact.reference_first = reference;
            fact.reference_last = reference_last;
            fact.ref_length = 2u;
            fact.start_lost_flanking = 1u;
            *out = fact;
            return DUCKVEP_HGVS_OK;
        }
        reference = hgvs_protein_window_base(
            context, &fact.window, 0, fact.ref_offset);
        reference_last = hgvs_protein_window_base(
            context, &fact.window, 0,
            fact.ref_offset + fact.ref_length - 1u);
        if (reference == 0u || reference_last == 0u) {
            return DUCKVEP_HGVS_MISSING_PEPTIDE;
        }
        fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_START_LOST;
        fact.reference_first = reference;
        fact.reference_last = reference_last;
        *out = fact;
        return DUCKVEP_HGVS_OK;
    }

    if (fact.ref_length == 0u && fact.alt_length == 0u) {
        uint8_t reference;
        size_t anchor0 = fact.first_position1 == 0u
            ? 0u : (size_t)fact.first_position1 - 1u;
        if (anchor0 >= context->ref_peptide_len && anchor0 != 0u) anchor0--;
        reference = duckvep_coding_context_peptide_base(
            context, 0, anchor0);
        if (reference == 0u) return DUCKVEP_HGVS_MISSING_PEPTIDE;
        fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_EQUAL;
        fact.first_position1 = (uint32_t)anchor0 + 1u;
        fact.last_position1 = fact.first_position1;
        fact.reference_first = reference;
        fact.reference_last = reference;
        fact.alternate_first = reference;
        *out = fact;
        return DUCKVEP_HGVS_OK;
    }
    if (fact.ref_length == fact.alt_length && fact.ref_length != 0u) {
        size_t i;
        int normalized_equal = 1;
        for (i = 0u; i < fact.ref_length; i++) {
            uint8_t reference = hgvs_protein_window_base(
                context, &fact.window, 0, fact.ref_offset + i);
            uint8_t alternate = hgvs_protein_window_base(
                context, &fact.window, 1, fact.alt_offset + i);
            if (reference == 0u || alternate == 0u) {
                return DUCKVEP_HGVS_MISSING_PEPTIDE;
            }
            if (xaa_as_ter && reference == (uint8_t)'*')
                reference = (uint8_t)'X';
            if (xaa_as_ter && alternate == (uint8_t)'*')
                alternate = (uint8_t)'X';
            if (reference != alternate) {
                normalized_equal = 0;
                break;
            }
        }
        if (normalized_equal) {
            fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_EQUAL;
            fact.last_position1 = fact.first_position1;
            status = hgvs_protein_fact_set_residues(&fact);
            if (status != DUCKVEP_HGVS_OK) return status;
            *out = fact;
            return DUCKVEP_HGVS_OK;
        }
    }
    if (fact.ref_length == 0u) {
        fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_INSERTION;
    } else if (fact.alt_length == 0u) {
        fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_DELETION;
    } else if (fact.ref_length == 1u && fact.alt_length == 1u) {
        fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_SUBSTITUTION;
    } else {
        fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_DELINS;
    }

    if (fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_INSERTION ||
        fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_DELETION) {
        status = hgvs_protein_shift_simple(&fact);
        if (status != DUCKVEP_HGVS_OK) return status;
    }
    status = hgvs_protein_fact_set_residues(&fact);
    if (status != DUCKVEP_HGVS_OK) return status;
    if (fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_INSERTION) {
        status = hgvs_protein_insertion_finish(&fact);
        if (status != DUCKVEP_HGVS_OK) return status;
    }

    if (delta->stop_lost &&
        (fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_SUBSTITUTION ||
         fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_DELETION)) {
        fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_EXTENSION;
        hgvs_protein_stop_distance(
            context, fact.first_position1, 0,
            &fact.termination_distance, &fact.termination_known);
    } else if (delta->stop_lost &&
               (fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_INSERTION ||
                fact.shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_DELINS)) {
        size_t i;
        int alternate_contains_stop = 0;

        for (i = 0u; i < fact.alt_length; i++) {
            uint8_t residue;
            status = duckvep_hgvs_protein_base(&fact, 1, i, &residue);
            if (status != DUCKVEP_HGVS_OK) return status;
            if (residue == (uint8_t)'*' ||
                (xaa_as_ter && residue == (uint8_t)'X')) {
                alternate_contains_stop = 1;
                break;
            }
        }
        hgvs_protein_stop_distance(
            context, fact.first_position1, 0,
            &fact.termination_distance, &fact.termination_known);
        /* VEP truncates delins/ins peptide text at its first Ter before it
         * considers extTer. An already represented alternate stop therefore
         * ends the operation; it is not also an extension. */
        fact.extends_stop = fact.termination_known &&
            !alternate_contains_stop;
    }
    *out = fact;
    return DUCKVEP_HGVS_OK;
}

duckvep_hgvs_status_t duckvep_hgvs_protein_fact_build_single_residue(
    uint32_t                     position1,
    uint8_t                      reference,
    uint8_t                      alternate,
    uint32_t                     consequence_flags,
    duckvep_compat_profile_t     compatibility_profile,
    duckvep_hgvs_protein_fact_t *out) {

    duckvep_hgvs_protein_fact_t fact;

    if (out == NULL) return DUCKVEP_HGVS_INVALID_ARG;
    memset(out, 0, sizeof *out);
    if (!duckvep_compat_profile_valid(compatibility_profile)) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    if (position1 == 0u || reference == 0u || alternate == 0u ||
        (consequence_flags & (uint32_t)
             DUCKVEP_CONSEQUENCE_FLAG_SEQUENCE_PREDICATES_VALID) == 0u ||
        (consequence_flags & ((uint32_t)DUCKVEP_CONSEQUENCE_FLAG_FRAMESHIFT |
                              (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_STOP_LOST)) !=
            0u) {
        return DUCKVEP_HGVS_NOT_APPLICABLE;
    }
    memset(&fact, 0, sizeof fact);
    fact.compatibility_profile = (uint8_t)compatibility_profile;
    fact.first_position1 = position1;
    fact.last_position1 = position1;
    fact.reference_first = reference;
    fact.reference_last = reference;
    fact.alternate_first = alternate;
    if ((consequence_flags &
         (uint32_t)DUCKVEP_CONSEQUENCE_FLAG_START_LOST) != 0u) {
        fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_START_LOST;
        fact.ref_length = 1u;
        fact.start_lost_flanking = 1u;
    } else if (reference == alternate) {
        /* The renderer's zero-length equality form uses the cached residue and
         * does not require a peptide-window pointer. */
        fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_EQUAL;
    } else {
        fact.shape = (uint8_t)DUCKVEP_HGVS_PROTEIN_SUBSTITUTION;
        fact.ref_length = 1u;
        fact.alt_length = 1u;
    }
    *out = fact;
    return DUCKVEP_HGVS_OK;
}

static const char *hgvs_protein_residue_name(
    uint8_t                   residue,
    duckvep_compat_profile_t  profile) {
    switch ((char)residue) {
        case 'A': return "Ala";
        case 'R': return "Arg";
        case 'N': return "Asn";
        case 'D': return "Asp";
        case 'C': return "Cys";
        case 'Q': return "Gln";
        case 'E': return "Glu";
        case 'G': return "Gly";
        case 'H': return "His";
        case 'I': return "Ile";
        case 'L': return "Leu";
        case 'K': return "Lys";
        case 'M': return "Met";
        case 'F': return "Phe";
        case 'P': return "Pro";
        case 'S': return "Ser";
        case 'T': return "Thr";
        case 'W': return "Trp";
        case 'Y': return "Tyr";
        case 'V': return "Val";
        case 'B': return "Asx";
        case 'Z': return "Glx";
        case 'J': return "Xle";
        case 'U': return "Sec";
        case 'O': return "Pyl";
        /* _get_hgvs_protein_type turns '*' into X; Bio::SeqUtils then emits
         * Xaa and VEP replaces Xaa with Ter. Preserve that observable VEP 116
         * conflation for both an explicit stop and its synthetic X residue. */
        case '*': return "Ter";
        case 'X': return duckvep_compat_enabled(
            profile, DUCKVEP_COMPAT_HGVS_XAA_AS_TER) ? "Ter" : "Xaa";
        default: return NULL;
    }
}

static int hgvs_writer_protein_residue(
    hgvs_writer_t *writer,
    uint8_t        residue,
    duckvep_compat_profile_t profile) {

    const char *name = hgvs_protein_residue_name(residue, profile);
    return name != NULL && hgvs_writer_literal(writer, name);
}

static duckvep_hgvs_status_t hgvs_writer_protein_sequence(
    hgvs_writer_t                       *writer,
    const duckvep_hgvs_protein_fact_t   *fact,
    int                                  alternate) {

    size_t length = alternate ? fact->alt_length : fact->ref_length;
    size_t i;
    duckvep_compat_profile_t profile =
        (duckvep_compat_profile_t)fact->compatibility_profile;
    for (i = 0u; i < length; i++) {
        uint8_t residue;
        duckvep_hgvs_status_t status = duckvep_hgvs_protein_base(
            fact, alternate, i, &residue);
        if (status != DUCKVEP_HGVS_OK) return status;
        if (!hgvs_writer_protein_residue(writer, residue, profile)) {
            return DUCKVEP_HGVS_INVALID_ALLELE;
        }
        /* _get_hgvs_protein_format truncates insertion/delins peptide text at
         * the first Ter before it writes the operation. */
        if (alternate &&
            (residue == (uint8_t)'*' ||
             (duckvep_compat_enabled(
                  profile, DUCKVEP_COMPAT_HGVS_XAA_AS_TER) &&
              residue == (uint8_t)'X'))) {
            break;
        }
    }
    return DUCKVEP_HGVS_OK;
}

static int hgvs_writer_protein_locus(
    hgvs_writer_t *writer,
    uint8_t        residue,
    uint32_t       position1,
    duckvep_compat_profile_t profile) {

    return hgvs_writer_protein_residue(writer, residue, profile) &&
           hgvs_writer_uint64(writer, position1);
}

static int hgvs_writer_termination(
    hgvs_writer_t                       *writer,
    const duckvep_hgvs_protein_fact_t   *fact) {

    return hgvs_writer_literal(writer, "Ter") &&
           (fact->termination_known
                ? hgvs_writer_uint64(writer, fact->termination_distance)
                : hgvs_writer_char(writer, '?'));
}

duckvep_hgvs_status_t duckvep_hgvs_protein_render(
    const duckvep_hgvs_protein_fact_t *fact,
    int                                predicted,
    char                              *buffer,
    size_t                             capacity,
    size_t                            *required_out) {

    hgvs_writer_t writer;
    duckvep_hgvs_status_t status;
    duckvep_compat_profile_t profile;

    if (required_out != NULL) *required_out = 0u;
    if (buffer != NULL && capacity > 0u) buffer[0] = '\0';
    if (fact == NULL || required_out == NULL ||
        (buffer == NULL && capacity != 0u) || predicted < 0 || predicted > 1 ||
        ((fact->first_position1 == 0u || fact->last_position1 == 0u) &&
         !((fact->shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_INSERTION &&
            fact->first_position1 == 0u && fact->last_position1 == 1u) ||
           (fact->shape == (uint8_t)DUCKVEP_HGVS_PROTEIN_START_LOST &&
            fact->first_position1 == 1u && fact->last_position1 == 0u)))) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    memset(&writer, 0, sizeof writer);
    profile = (duckvep_compat_profile_t)fact->compatibility_profile;
    if (!duckvep_compat_profile_valid(profile) ||
        (fact->context != NULL &&
         (!duckvep_compat_profile_valid(
              (duckvep_compat_profile_t)
                  fact->context->compatibility_profile) ||
          fact->context->compatibility_profile !=
              fact->compatibility_profile))) {
        return DUCKVEP_HGVS_INVALID_ARG;
    }
    if ((fact->first_position1 == 0u || fact->last_position1 == 0u) &&
        !duckvep_compat_enabled(
            profile, DUCKVEP_COMPAT_HGVS_NEGATIVE_SUBSTR)) {
        return DUCKVEP_HGVS_NOT_APPLICABLE;
    }
    writer.buffer = buffer;
    writer.capacity = capacity;
    if (!hgvs_writer_literal(&writer, "p.") ||
        (predicted && !hgvs_writer_char(&writer, '('))) {
        return DUCKVEP_HGVS_OUT_OF_RANGE;
    }

    switch ((duckvep_hgvs_protein_shape_t)fact->shape) {
        case DUCKVEP_HGVS_PROTEIN_EQUAL:
            /* _get_fs_peptides() can discover equality only at the synthetic
             * reference stop appended for its comparison.  That residue is
             * cached in the fact but is outside the local peptide window. */
            if (fact->ref_length == 0u && fact->reference_first != 0u) {
                if (!hgvs_writer_protein_residue(
                        &writer, fact->reference_first, profile)) {
                    return DUCKVEP_HGVS_INVALID_PROJECTION;
                }
            } else {
                status = hgvs_writer_protein_sequence(&writer, fact, 0);
                if (status != DUCKVEP_HGVS_OK) return status;
            }
            if (!hgvs_writer_uint64(&writer, fact->first_position1) ||
                !hgvs_writer_char(&writer, '=')) {
                return DUCKVEP_HGVS_INVALID_PROJECTION;
            }
            break;
        case DUCKVEP_HGVS_PROTEIN_START_LOST:
            if (fact->start_lost_flanking) {
                if (!hgvs_writer_protein_residue(
                        &writer, fact->reference_first, profile) ||
                    (fact->ref_length > 1u &&
                     !hgvs_writer_protein_residue(
                         &writer, fact->reference_last, profile))) {
                    return DUCKVEP_HGVS_INVALID_PROJECTION;
                }
            } else {
                status = hgvs_writer_protein_sequence(&writer, fact, 0);
                if (status != DUCKVEP_HGVS_OK) return status;
            }
            if (!hgvs_writer_uint64(&writer, fact->first_position1) ||
                (fact->first_position1 != fact->last_position1 &&
                 (!hgvs_writer_char(&writer, '_') ||
                  !hgvs_writer_char(&writer, '?') ||
                  !hgvs_writer_uint64(&writer, fact->last_position1))) ||
                (fact->first_position1 == fact->last_position1 &&
                 !hgvs_writer_char(&writer, '?'))) {
                return DUCKVEP_HGVS_INVALID_PROJECTION;
            }
            break;
        case DUCKVEP_HGVS_PROTEIN_SUBSTITUTION:
            if (!hgvs_writer_protein_locus(
                    &writer, fact->reference_first,
                    fact->first_position1, profile) ||
                !hgvs_writer_protein_residue(
                    &writer, fact->alternate_first, profile)) {
                return DUCKVEP_HGVS_INVALID_PROJECTION;
            }
            break;
        case DUCKVEP_HGVS_PROTEIN_DELETION:
            if (!hgvs_writer_protein_locus(
                    &writer, fact->reference_first,
                    fact->first_position1, profile) ||
                (fact->first_position1 != fact->last_position1 &&
                 (!hgvs_writer_char(&writer, '_') ||
                  !hgvs_writer_protein_locus(
                      &writer, fact->reference_last,
                      fact->last_position1, profile))) ||
                !hgvs_writer_literal(&writer, "del")) {
                return DUCKVEP_HGVS_INVALID_PROJECTION;
            }
            break;
        case DUCKVEP_HGVS_PROTEIN_INSERTION:
            if (!hgvs_writer_protein_locus(
                    &writer, fact->reference_first,
                    fact->first_position1, profile) ||
                !hgvs_writer_char(&writer, '_') ||
                !hgvs_writer_protein_locus(
                    &writer, fact->reference_last,
                    fact->last_position1, profile) ||
                !hgvs_writer_literal(&writer, "ins")) {
                return DUCKVEP_HGVS_INVALID_PROJECTION;
            }
            status = hgvs_writer_protein_sequence(&writer, fact, 1);
            if (status != DUCKVEP_HGVS_OK) return status;
            if (fact->extends_stop &&
                (!hgvs_writer_literal(&writer, "ext") ||
                 !hgvs_writer_termination(&writer, fact))) {
                return DUCKVEP_HGVS_OUT_OF_RANGE;
            }
            break;
        case DUCKVEP_HGVS_PROTEIN_DELINS:
            if (!hgvs_writer_protein_locus(
                    &writer, fact->reference_first,
                    fact->first_position1, profile) ||
                (fact->first_position1 != fact->last_position1 &&
                 (!hgvs_writer_char(&writer, '_') ||
                  !hgvs_writer_protein_locus(
                      &writer, fact->reference_last,
                      fact->last_position1, profile))) ||
                !hgvs_writer_literal(&writer, "delins")) {
                return DUCKVEP_HGVS_INVALID_PROJECTION;
            }
            status = hgvs_writer_protein_sequence(&writer, fact, 1);
            if (status != DUCKVEP_HGVS_OK) return status;
            if (fact->extends_stop &&
                (!hgvs_writer_literal(&writer, "ext") ||
                 !hgvs_writer_termination(&writer, fact))) {
                return DUCKVEP_HGVS_OUT_OF_RANGE;
            }
            break;
        case DUCKVEP_HGVS_PROTEIN_DUPLICATION:
            if (!hgvs_writer_protein_locus(
                    &writer, fact->reference_first,
                    fact->first_position1, profile) ||
                (fact->first_position1 != fact->last_position1 &&
                 (!hgvs_writer_char(&writer, '_') ||
                  !hgvs_writer_protein_locus(
                      &writer, fact->reference_last,
                      fact->last_position1, profile))) ||
                !hgvs_writer_literal(&writer, "dup")) {
                return DUCKVEP_HGVS_INVALID_PROJECTION;
            }
            break;
        case DUCKVEP_HGVS_PROTEIN_FRAMESHIFT:
            if (!hgvs_writer_protein_locus(
                    &writer, fact->reference_first,
                    fact->first_position1, profile) ||
                !hgvs_writer_protein_residue(
                    &writer, fact->alternate_first, profile) ||
                !hgvs_writer_literal(&writer, "fs") ||
                !hgvs_writer_termination(&writer, fact)) {
                return DUCKVEP_HGVS_INVALID_PROJECTION;
            }
            break;
        case DUCKVEP_HGVS_PROTEIN_EXTENSION:
            if (!hgvs_writer_protein_locus(
                    &writer, fact->reference_first,
                    fact->first_position1, profile) ||
                (fact->first_position1 != fact->last_position1 &&
                 (!hgvs_writer_char(&writer, '_') ||
                  !hgvs_writer_protein_locus(
                      &writer, fact->reference_last,
                      fact->last_position1, profile)))) {
                return DUCKVEP_HGVS_INVALID_PROJECTION;
            }
            if (fact->alternate_first != 0u) {
                if (!hgvs_writer_protein_residue(
                        &writer, fact->alternate_first, profile)) {
                    return DUCKVEP_HGVS_INVALID_ALLELE;
                }
            } else if (!hgvs_writer_literal(&writer, "del")) {
                return DUCKVEP_HGVS_OUT_OF_RANGE;
            }
            if (!hgvs_writer_literal(&writer, "ext") ||
                !hgvs_writer_termination(&writer, fact)) {
                return DUCKVEP_HGVS_OUT_OF_RANGE;
            }
            break;
        default:
            return DUCKVEP_HGVS_INVALID_ARG;
    }
    if (predicted && !hgvs_writer_char(&writer, ')')) {
        return DUCKVEP_HGVS_OUT_OF_RANGE;
    }
    return hgvs_writer_finish(&writer, required_out);
}
