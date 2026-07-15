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

int duckvep_project_coding_base(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          genomic_pos,
    duckvep_coding_projection_t      *out) {

    duckvep_coding_projection_t r;
    uint32_t cdna = 0u;
    uint32_t exon_idx = 0u;
    uint32_t coding_start_genomic;
    uint32_t coding_end_genomic;
    uint32_t coding_start_cdna = 0u;
    uint32_t coding_end_cdna = 0u;
    uint32_t coding_start_exon_idx = 0u;
    uint32_t dummy_exon_idx = 0u;
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
        if (coding_start_exon_idx >= exons->exon_count || phase > 2u) {
            return 0;
        }
    } else {
        if (transcripts->strand[tx_idx] >= 0) {
            coding_start_genomic = transcripts->cds_start1[tx_idx];
            coding_end_genomic = transcripts->cds_end1[tx_idx];
        } else {
            coding_start_genomic = transcripts->cds_end1[tx_idx];
            coding_end_genomic = transcripts->cds_start1[tx_idx];
        }
        if (!duckvep_project_genomic_to_cdna(transcripts, exons, tx_idx,
                                             coding_start_genomic,
                                             &coding_start_cdna,
                                             &coding_start_exon_idx) ||
            !duckvep_project_genomic_to_cdna(transcripts, exons, tx_idx,
                                             coding_end_genomic,
                                             &coding_end_cdna,
                                             &dummy_exon_idx)) {
            return 0;
        }
        phase = coding_phase_offset(exons, (size_t)coding_start_exon_idx);
    }
    if (coding_end_cdna < coding_start_cdna) return 0;
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
