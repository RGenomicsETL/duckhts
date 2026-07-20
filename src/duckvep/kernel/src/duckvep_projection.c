/*
 * duckvep_projection.c — transcript coordinate projection.
 * See duckvep_projection.h. No allocation; no sequence access.
 */
#include "duckvep_projection.h"

#include <string.h>

static int valid_tx_exon_slice(const duckvep_transcript_model_t *tx,
                               const duckvep_exon_model_t *ex,
                               size_t tx_idx,
                               size_t *off_out,
                               size_t *cnt_out) {
    size_t off;
    size_t cnt;

    if (tx == NULL || ex == NULL || off_out == NULL || cnt_out == NULL) return 0;
    if (tx_idx >= tx->transcript_count) return 0;
    if (tx->exon_offset == NULL || tx->exon_count == NULL || tx->strand == NULL) {
        return 0;
    }
    if (ex->start1 == NULL || ex->end1 == NULL ||
        ex->cdna_start1 == NULL || ex->cdna_end1 == NULL) {
        return 0;
    }

    off = (size_t)tx->exon_offset[tx_idx];
    cnt = (size_t)tx->exon_count[tx_idx];
    if (cnt == 0u || off > ex->exon_count || cnt > ex->exon_count - off) return 0;

    *off_out = off;
    *cnt_out = cnt;
    return 1;
}

static uint8_t coding_phase_offset(const duckvep_exon_model_t *ex, size_t exon_idx) {
    int8_t ph;
    if (ex == NULL || ex->phase == NULL || exon_idx >= ex->exon_count) return 0u;
    ph = ex->phase[exon_idx];
    return ph > 0 ? (uint8_t)ph : 0u;
}

int duckvep_project_genomic_to_cdna(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          genomic_pos,
    uint32_t                         *cdna_pos_out,
    uint32_t                         *exon_idx_out) {

    size_t off;
    size_t cnt;
    size_t lo;
    size_t hi;
    size_t i;
    int fwd;

    if (cdna_pos_out == NULL) return 0;
    *cdna_pos_out = 0u;
    if (exon_idx_out != NULL) *exon_idx_out = 0u;
    if (!valid_tx_exon_slice(transcripts, exons, tx_idx, &off, &cnt)) return 0;

    fwd = transcripts->strand[tx_idx] >= 0;
    lo = 0u;
    hi = cnt;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2u;
        size_t ei = off + mid;

        if (fwd ? exons->end1[ei] < genomic_pos
                : exons->start1[ei] > genomic_pos) {
            lo = mid + 1u;
        } else {
            hi = mid;
        }
    }
    if (lo >= cnt) return 0;
    i = off + lo;
    {
        uint32_t es = exons->start1[i];
        uint32_t ee = exons->end1[i];
        uint32_t cs = exons->cdna_start1[i];
        uint32_t ce = exons->cdna_end1[i];
        uint32_t cdna;

        if (es > ee || cs == 0u || ce < cs ||
            genomic_pos < es || genomic_pos > ee) {
            return 0;
        }
        cdna = fwd ? (uint32_t)(cs + (genomic_pos - es))
                   : (uint32_t)(cs + (ee - genomic_pos));
        if (cdna > ce) return 0;
        *cdna_pos_out = cdna;
        if (exon_idx_out != NULL) *exon_idx_out = (uint32_t)i;
        return 1;
    }
}

int duckvep_project_cdna_to_genomic(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          cdna_pos,
    uint32_t                         *genomic_pos_out,
    uint32_t                         *exon_idx_out) {

    size_t off;
    size_t cnt;
    size_t lo;
    size_t hi;
    size_t i;
    int fwd;

    if (genomic_pos_out == NULL || cdna_pos == 0u) return 0;
    *genomic_pos_out = 0u;
    if (exon_idx_out != NULL) *exon_idx_out = 0u;
    if (!valid_tx_exon_slice(transcripts, exons, tx_idx, &off, &cnt)) return 0;

    fwd = transcripts->strand[tx_idx] >= 0;
    lo = 0u;
    hi = cnt;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2u;
        size_t ei = off + mid;

        if (exons->cdna_end1[ei] < cdna_pos)
            lo = mid + 1u;
        else
            hi = mid;
    }
    if (lo >= cnt) return 0;
    i = off + lo;
    {
        uint32_t es = exons->start1[i];
        uint32_t ee = exons->end1[i];
        uint32_t cs = exons->cdna_start1[i];
        uint32_t ce = exons->cdna_end1[i];
        uint32_t genomic;

        if (es > ee || cs == 0u || ce < cs ||
            cdna_pos < cs || cdna_pos > ce) {
            return 0;
        }
        genomic = fwd ? (uint32_t)(es + (cdna_pos - cs))
                      : (uint32_t)(ee - (cdna_pos - cs));
        if (genomic < es || genomic > ee) return 0;
        *genomic_pos_out = genomic;
        if (exon_idx_out != NULL) *exon_idx_out = (uint32_t)i;
        return 1;
    }
}

int duckvep_project_coding_cdna_bounds(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                         *coding_start_cdna_out,
    uint32_t                         *coding_end_cdna_out,
    uint32_t                         *coding_start_exon_idx_out,
    uint8_t                          *phase_offset_out) {

    uint32_t start_genomic;
    uint32_t end_genomic;
    uint32_t coding_start_cdna;
    uint32_t coding_end_cdna;
    uint32_t coding_start_exon_idx;
    uint32_t ignored_exon_idx;
    uint8_t phase;

    if (coding_start_cdna_out != NULL) *coding_start_cdna_out = 0u;
    if (coding_end_cdna_out != NULL) *coding_end_cdna_out = 0u;
    if (coding_start_exon_idx_out != NULL) {
        *coding_start_exon_idx_out = 0u;
    }
    if (phase_offset_out != NULL) *phase_offset_out = 0u;
    if (transcripts == NULL || exons == NULL ||
        coding_start_cdna_out == NULL || coding_end_cdna_out == NULL ||
        tx_idx >= transcripts->transcript_count ||
        transcripts->strand == NULL || transcripts->cds_start1 == NULL ||
        transcripts->cds_end1 == NULL ||
        transcripts->cds_start1[tx_idx] == 0u ||
        transcripts->cds_end1[tx_idx] < transcripts->cds_start1[tx_idx]) {
        return 0;
    }

    if (transcripts->cds_cdna_start1 != NULL &&
        transcripts->cds_cdna_end1 != NULL &&
        transcripts->cds_start_exon_index != NULL &&
        transcripts->cds_phase_offset != NULL &&
        transcripts->cds_cdna_start1[tx_idx] != 0u) {
        coding_start_cdna = transcripts->cds_cdna_start1[tx_idx];
        coding_end_cdna = transcripts->cds_cdna_end1[tx_idx];
        coding_start_exon_idx =
            transcripts->cds_start_exon_index[tx_idx];
        phase = transcripts->cds_phase_offset[tx_idx];
        if (coding_end_cdna < coding_start_cdna ||
            coding_start_exon_idx >= exons->exon_count || phase > 2u) {
            return 0;
        }
    } else {
        start_genomic = transcripts->strand[tx_idx] > 0
            ? transcripts->cds_start1[tx_idx]
            : transcripts->cds_end1[tx_idx];
        end_genomic = transcripts->strand[tx_idx] > 0
            ? transcripts->cds_end1[tx_idx]
            : transcripts->cds_start1[tx_idx];
        if (!duckvep_project_genomic_to_cdna(
                transcripts, exons, tx_idx, start_genomic,
                &coding_start_cdna, &coding_start_exon_idx) ||
            !duckvep_project_genomic_to_cdna(
                transcripts, exons, tx_idx, end_genomic,
                &coding_end_cdna, &ignored_exon_idx) ||
            coding_end_cdna < coding_start_cdna) {
            return 0;
        }
        phase = coding_phase_offset(exons, (size_t)coding_start_exon_idx);
    }

    *coding_start_cdna_out = coding_start_cdna;
    *coding_end_cdna_out = coding_end_cdna;
    if (coding_start_exon_idx_out != NULL) {
        *coding_start_exon_idx_out = coding_start_exon_idx;
    }
    if (phase_offset_out != NULL) *phase_offset_out = phase;
    return 1;
}

static int project_mapper_insert_cdna_bounds(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    uint32_t                         *cdna_start_out,
    uint32_t                         *cdna_end_out) {

    uint32_t after_genomic;
    uint32_t before_genomic;
    uint32_t after_cdna = 0u;
    uint32_t before_cdna = 0u;
    int after_exonic;
    int before_exonic;

    if (transcripts == NULL || exons == NULL || event == NULL ||
        cdna_start_out == NULL || cdna_end_out == NULL ||
        tx_idx >= transcripts->transcript_count ||
        transcripts->strand == NULL || event->feature_start1 == 0u ||
        event->feature_end1 == UINT32_MAX ||
        event->feature_start1 != event->feature_end1 + 1u) {
        return 0;
    }
    /* In transcript order, map_insert's returned start is the flank after the
     * insertion and its returned end is the flank before it. When one flank
     * maps to a Gap, Mapper::map_insert drops that item and moves the surviving
     * coordinate one base inward to retain start=end+1 geometry. */
    after_genomic = transcripts->strand[tx_idx] > 0
        ? event->feature_start1 : event->feature_end1;
    before_genomic = transcripts->strand[tx_idx] > 0
        ? event->feature_end1 : event->feature_start1;
    after_exonic = after_genomic != 0u &&
        duckvep_project_genomic_to_cdna(
            transcripts, exons, tx_idx, after_genomic, &after_cdna, NULL);
    before_exonic = before_genomic != 0u &&
        duckvep_project_genomic_to_cdna(
            transcripts, exons, tx_idx, before_genomic, &before_cdna, NULL);
    if (!after_exonic && !before_exonic) return 0;
    if (after_exonic && before_exonic) {
        *cdna_start_out = after_cdna;
        *cdna_end_out = before_cdna;
    } else if (before_exonic) {
        if (before_cdna == UINT32_MAX) return 0;
        *cdna_start_out = before_cdna + 1u;
        *cdna_end_out = before_cdna;
    } else {
        *cdna_start_out = after_cdna;
        *cdna_end_out = after_cdna - 1u;
    }
    return 1;
}

int duckvep_project_feature_overlaps_start_codon_unshifted(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event) {

    uint32_t cdna_start;
    uint32_t cdna_end;
    uint32_t coding_start;
    uint32_t coding_end;

    if (transcripts == NULL || exons == NULL || event == NULL ||
        tx_idx >= transcripts->transcript_count ||
        transcripts->strand == NULL || transcripts->flags == NULL ||
        event->feature_start1 == 0u ||
        (transcripts->flags[tx_idx] &
         (uint64_t)DUCKVEP_TX_CDS_START_NF) != 0u) {
        return 0;
    }
    if (event->feature_end1 != UINT32_MAX &&
        event->feature_start1 == event->feature_end1 + 1u) {
        if (!project_mapper_insert_cdna_bounds(
                transcripts, exons, tx_idx, event,
                &cdna_start, &cdna_end)) {
            return 0;
        }
    } else {
        uint32_t first_genomic = transcripts->strand[tx_idx] > 0
            ? event->feature_start1 : event->feature_end1;
        uint32_t last_genomic = transcripts->strand[tx_idx] > 0
            ? event->feature_end1 : event->feature_start1;
        if (event->feature_end1 == 0u ||
            !duckvep_project_genomic_to_cdna(
                transcripts, exons, tx_idx, first_genomic,
                &cdna_start, NULL) ||
            !duckvep_project_genomic_to_cdna(
                transcripts, exons, tx_idx, last_genomic,
                &cdna_end, NULL)) {
            return 0;
        }
    }
    if (cdna_start == 0u || cdna_end == 0u ||
        !duckvep_project_coding_cdna_bounds(
            transcripts, exons, tx_idx, &coding_start, &coding_end,
            NULL, NULL) ||
        coding_start > UINT32_MAX - 2u) {
        return 0;
    }
    (void)coding_end;
    return cdna_end >= coding_start && cdna_start <= coding_start + 2u;
}

int duckvep_project_feature_has_coding_precondition_unshifted(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event) {

    size_t offset;
    size_t count;
    size_t i;
    uint32_t low;
    uint32_t high;
    uint32_t cds_low;
    uint32_t cds_high;

    if (!valid_tx_exon_slice(transcripts, exons, tx_idx, &offset, &count) ||
        event == NULL || transcripts->start1 == NULL ||
        transcripts->end1 == NULL || transcripts->cds_start1 == NULL ||
        transcripts->cds_end1 == NULL || event->feature_start1 == 0u ||
        event->feature_end1 == 0u ||
        transcripts->cds_start1[tx_idx] == 0u ||
        transcripts->cds_end1[tx_idx] < transcripts->cds_start1[tx_idx]) {
        return 0;
    }
    /* BaseVariationFeatureOverlapAllele::_bvfo_preds applies overlap() to
     * reversed insertion coordinates before normalizing them for exon trees. */
    if (event->feature_end1 < transcripts->start1[tx_idx] ||
        event->feature_start1 > transcripts->end1[tx_idx]) {
        return 0;
    }
    low = event->feature_start1 < event->feature_end1
        ? event->feature_start1 : event->feature_end1;
    high = event->feature_start1 > event->feature_end1
        ? event->feature_start1 : event->feature_end1;
    for (i = offset; i < offset + count; i++) {
        cds_low = exons->start1[i] > transcripts->cds_start1[tx_idx]
            ? exons->start1[i] : transcripts->cds_start1[tx_idx];
        cds_high = exons->end1[i] < transcripts->cds_end1[tx_idx]
            ? exons->end1[i] : transcripts->cds_end1[tx_idx];
        if (cds_low <= cds_high && high >= cds_low && low <= cds_high) {
            return 1;
        }
    }
    return 0;
}

int duckvep_project_complete_feature_translation_bounds(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    uint32_t                          genomic_shift,
    uint32_t                         *cds_start_out,
    uint32_t                         *cds_end_out) {

    duckvep_coding_projection_t first;
    duckvep_coding_projection_t last;
    uint32_t first_cdna;
    uint32_t last_cdna;
    uint32_t first_exon;
    uint32_t last_exon;
    uint32_t feature_start;
    uint32_t feature_end;
    uint32_t first_genomic;
    uint32_t last_genomic;
    int first_exonic;
    int last_exonic;
    int first_coding;
    int last_coding;
    int8_t strand;

    if (cds_start_out == NULL || cds_end_out == NULL) return 0;
    *cds_start_out = 0u;
    *cds_end_out = 0u;
    if (transcripts == NULL || exons == NULL || event == NULL ||
        tx_idx >= transcripts->transcript_count ||
        transcripts->strand == NULL || event->feature_start1 == 0u ||
        event->feature_end1 == 0u) {
        return 0;
    }
    strand = transcripts->strand[tx_idx];
    feature_start = event->feature_start1;
    feature_end = event->feature_end1;
    if (strand > 0) {
        if (feature_start > UINT32_MAX - genomic_shift ||
            feature_end > UINT32_MAX - genomic_shift) {
            return 0;
        }
        feature_start += genomic_shift;
        feature_end += genomic_shift;
    } else {
        if (feature_start <= genomic_shift || feature_end <= genomic_shift) {
            return 0;
        }
        feature_start -= genomic_shift;
        feature_end -= genomic_shift;
    }
    first_genomic = strand > 0 ? feature_start : feature_end;
    last_genomic = strand > 0 ? feature_end : feature_start;
    if (feature_end != UINT32_MAX && feature_start == feature_end + 1u) {
        /* Mapper::map_insert first maps both flanks as a two-base interval.
         * Intronic/transcript-external Gap items are then removed, while an
         * exonic UTR flank survives genomic2cdna() and becomes a Gap only in
         * genomic2cds(). Thus one coding exon flank is sufficient when the
         * other flank is not exonic, but two exonic flanks must both be in
         * CDS. This is observably different from projecting both raw genomic
         * endpoints independently. */
        first_exonic = duckvep_project_genomic_to_cdna(
            transcripts, exons, tx_idx, first_genomic,
            &first_cdna, &first_exon);
        last_exonic = duckvep_project_genomic_to_cdna(
            transcripts, exons, tx_idx, last_genomic,
            &last_cdna, &last_exon);
        (void)first_cdna;
        (void)last_cdna;
        (void)first_exon;
        (void)last_exon;
        first_coding = duckvep_project_coding_base(
            transcripts, exons, tx_idx, first_genomic, &first);
        last_coding = duckvep_project_coding_base(
            transcripts, exons, tx_idx, last_genomic, &last);
        if ((first_exonic && last_exonic &&
             (!first_coding || !last_coding)) ||
            (!first_exonic && !last_exonic) ||
            (first_exonic && !last_exonic && !first_coding) ||
            (!first_exonic && last_exonic && !last_coding)) {
            return 0;
        }
        if (first_coding && last_coding) {
            *cds_start_out = first.cds_pos;
            *cds_end_out = last.cds_pos;
        } else if (first_coding) {
            *cds_start_out = first.cds_pos;
            *cds_end_out = first.cds_pos;
        } else {
            *cds_start_out = last.cds_pos;
            *cds_end_out = last.cds_pos;
        }
        return 1;
    }
    if (!duckvep_project_coding_base(
            transcripts, exons, tx_idx, first_genomic, &first) ||
        !duckvep_project_coding_base(
            transcripts, exons, tx_idx, last_genomic, &last)) {
        return 0;
    }
    *cds_start_out = first.cds_pos;
    *cds_end_out = last.cds_pos;
    return 1;
}

int duckvep_project_coding_base(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          genomic_pos,
    duckvep_coding_projection_t      *out) {

    duckvep_coding_projection_t r;
    uint32_t cdna = 0u;
    uint32_t exon_idx = 0u;
    uint32_t coding_start_cdna = 0u;
    uint32_t coding_end_cdna = 0u;
    uint32_t cds_pos;
    uint8_t phase;

    if (out == NULL || transcripts == NULL) return 0;
    memset(&r, 0, sizeof r);
    *out = r;
    if (tx_idx >= transcripts->transcript_count) return 0;
    if (transcripts->strand == NULL || transcripts->cds_start1 == NULL ||
        transcripts->cds_end1 == NULL) {
        return 0;
    }
    if (transcripts->cds_start1[tx_idx] == 0u || transcripts->cds_end1[tx_idx] == 0u ||
        transcripts->cds_end1[tx_idx] < transcripts->cds_start1[tx_idx]) {
        return 0;
    }

    if (!duckvep_project_genomic_to_cdna(transcripts, exons, tx_idx, genomic_pos,
                                         &cdna, &exon_idx)) {
        return 0;
    }

    if (!duckvep_project_coding_cdna_bounds(
            transcripts, exons, tx_idx, &coding_start_cdna,
            &coding_end_cdna, NULL, &phase)) {
        return 0;
    }
    if (cdna < coding_start_cdna || cdna > coding_end_cdna) return 0;
    cds_pos = (uint32_t)(cdna - coding_start_cdna + 1u + (uint32_t)phase);

    r.cdna_pos = cdna;
    r.cds_pos = cds_pos;
    r.protein_pos = (uint32_t)((cds_pos - 1u) / 3u + 1u);
    r.codon_offset = (uint8_t)((cds_pos - 1u) % 3u);
    r.codon_start_cds = (uint32_t)(cds_pos - (uint32_t)r.codon_offset);
    r.exon_idx = exon_idx;
    r.phase_offset = phase;
    *out = r;
    return 1;
}

int duckvep_project_feature_translation_start(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    duckvep_coding_projection_t      *out) {

    uint32_t first_genomic1;

    if (out == NULL) return 0;
    memset(out, 0, sizeof *out);
    if (transcripts == NULL || event == NULL || event->interbase ||
        tx_idx >= transcripts->transcript_count ||
        event->feature_start1 == 0u ||
        event->feature_end1 < event->feature_start1) {
        return 0;
    }

    /* TranscriptMapper reverses mapper pieces for a reverse-strand
     * transcript. Its first item is therefore the high genomic endpoint. */
    first_genomic1 = transcripts->strand[tx_idx] >= 0
        ? event->feature_start1 : event->feature_end1;
    return duckvep_project_coding_base(
        transcripts, exons, tx_idx, first_genomic1, out);
}

int duckvep_cds_position_is_partial_codon(
    size_t   cds_length,
    uint32_t cds_position1) {

    size_t codon_start0;
    size_t last_codon_length;

    if (cds_position1 == 0u || (size_t)cds_position1 > cds_length ||
        (cds_length % 3u) == 0u) {
        return 0;
    }
    codon_start0 = ((size_t)cds_position1 - 1u) / 3u * 3u;
    last_codon_length = cds_length - codon_start0;
    return last_codon_length > 0u && last_codon_length < 3u;
}

int duckvep_project_event_to_cds(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    uint32_t                         *cds_start_out,
    uint32_t                         *cds_end_out) {

    duckvep_coding_projection_t first_projection;
    duckvep_coding_projection_t last_projection;
    uint32_t min_cds = UINT32_MAX;
    uint32_t max_cds = 0u;

    if (cds_start_out == NULL || cds_end_out == NULL) return 0;
    *cds_start_out = 0u;
    *cds_end_out = 0u;
    if (transcripts == NULL || exons == NULL || event == NULL ||
        tx_idx >= transcripts->transcript_count) {
        return 0;
    }
    if (event->interbase) {
        uint32_t boundary = event->insertion_boundary0;
        uint32_t right = duckvep_event_right_flank1(event);
        uint32_t projected_pos = event->start1;
        uint32_t cds_boundary;

        if (!duckvep_project_coding_base(transcripts, exons, tx_idx,
                                         projected_pos, &first_projection)) {
            projected_pos = event->anchor_side ==
                                (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT
                              ? right
                              : boundary;
            if (projected_pos == 0u ||
                !duckvep_project_coding_base(transcripts, exons, tx_idx,
                                             projected_pos, &first_projection)) {
                return 0;
            }
        }
        if (projected_pos == boundary) {
            if (transcripts->strand[tx_idx] > 0) {
                if (first_projection.cds_pos == UINT32_MAX) return 0;
                cds_boundary = first_projection.cds_pos + 1u;
            } else {
                cds_boundary = first_projection.cds_pos;
            }
        } else if (projected_pos == right) {
            if (transcripts->strand[tx_idx] > 0) {
                cds_boundary = first_projection.cds_pos;
            } else {
                if (first_projection.cds_pos == UINT32_MAX) return 0;
                cds_boundary = first_projection.cds_pos + 1u;
            }
        } else {
            return 0;
        }
        if (cds_boundary == 0u) return 0;
        *cds_start_out = cds_boundary;
        *cds_end_out = cds_boundary;
        return 1;
    }

    if (event->ref_diff_length == 0u || event->start1 == 0u ||
        (uint64_t)event->ref_diff_length - 1u >
            (uint64_t)UINT32_MAX - (uint64_t)event->start1) {
        return 0;
    }
    if (!duckvep_project_coding_base(
            transcripts, exons, tx_idx, event->start1, &first_projection) ||
        !duckvep_project_coding_base(
            transcripts, exons, tx_idx,
            event->start1 + (uint32_t)event->ref_diff_length - 1u,
            &last_projection)) {
        return 0;
    }
    min_cds = first_projection.cds_pos < last_projection.cds_pos
        ? first_projection.cds_pos : last_projection.cds_pos;
    max_cds = first_projection.cds_pos > last_projection.cds_pos
        ? first_projection.cds_pos : last_projection.cds_pos;
    if (min_cds == UINT32_MAX || max_cds < min_cds ||
        max_cds - min_cds + 1u != (uint32_t)event->ref_diff_length) {
        return 0;
    }
    *cds_start_out = min_cds;
    *cds_end_out = max_cds;
    return 1;
}

int duckvep_project_feature_to_cds(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    uint32_t                         *cds_start_out,
    uint32_t                         *cds_end_out) {

    duckvep_coding_projection_t first;
    duckvep_coding_projection_t last;
    uint32_t feature_start1;
    uint32_t feature_end1;

    if (cds_start_out == NULL || cds_end_out == NULL) return 0;
    *cds_start_out = 0u;
    *cds_end_out = 0u;
    if (transcripts == NULL || exons == NULL || event == NULL ||
        tx_idx >= transcripts->transcript_count) {
        return 0;
    }

    feature_start1 = event->feature_start1;
    feature_end1 = event->feature_end1;
    if (event->interbase) {
        uint32_t cds_boundary;
        uint32_t ignored_end;

        if (feature_start1 == 0u || feature_end1 == UINT32_MAX ||
            feature_start1 != feature_end1 + 1u) {
            return 0;
        }
        /* VEP's mapper projects an empty insertion interval from whichever
         * genomic flank maps, including an exon edge. The ordinary insertion
         * projector implements that rule and returns the transcript boundary.
         * TranscriptVariation then exposes the boundary as a reversed range. */
        if (!duckvep_project_event_to_cds(
                transcripts, exons, tx_idx, event,
                &cds_boundary, &ignored_end) || cds_boundary == 0u) {
            return 0;
        }
        *cds_start_out = cds_boundary;
        *cds_end_out = cds_boundary - 1u;
        return 1;
    }
    if (feature_start1 == 0u || feature_end1 == 0u ||
        feature_end1 < feature_start1) {
        return 0;
    }
    if (!duckvep_project_coding_base(
            transcripts, exons, tx_idx, feature_start1, &first) ||
        !duckvep_project_coding_base(
            transcripts, exons, tx_idx, feature_end1, &last)) {
        return 0;
    }
    *cds_start_out = first.cds_pos < last.cds_pos
        ? first.cds_pos : last.cds_pos;
    *cds_end_out = first.cds_pos > last.cds_pos
        ? first.cds_pos : last.cds_pos;
    return *cds_start_out != 0u && *cds_end_out != 0u;
}
