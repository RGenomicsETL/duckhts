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
    size_t i;
    int fwd;

    if (cdna_pos_out == NULL) return 0;
    *cdna_pos_out = 0u;
    if (exon_idx_out != NULL) *exon_idx_out = 0u;
    if (!valid_tx_exon_slice(transcripts, exons, tx_idx, &off, &cnt)) return 0;

    fwd = transcripts->strand[tx_idx] >= 0;
    for (i = off; i < off + cnt; i++) {
        uint32_t es = exons->start1[i];
        uint32_t ee = exons->end1[i];
        uint32_t cs = exons->cdna_start1[i];
        uint32_t ce = exons->cdna_end1[i];
        uint32_t cdna;

        if (es > ee || cs == 0u || ce < cs) return 0;
        if (genomic_pos < es || genomic_pos > ee) continue;

        cdna = fwd ? (uint32_t)(cs + (genomic_pos - es))
                   : (uint32_t)(cs + (ee - genomic_pos));
        if (cdna > ce) return 0;
        *cdna_pos_out = cdna;
        if (exon_idx_out != NULL) *exon_idx_out = (uint32_t)i;
        return 1;
    }
    return 0;
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
    size_t i;
    int fwd;

    if (genomic_pos_out == NULL || cdna_pos == 0u) return 0;
    *genomic_pos_out = 0u;
    if (exon_idx_out != NULL) *exon_idx_out = 0u;
    if (!valid_tx_exon_slice(transcripts, exons, tx_idx, &off, &cnt)) return 0;

    fwd = transcripts->strand[tx_idx] >= 0;
    for (i = off; i < off + cnt; i++) {
        uint32_t es = exons->start1[i];
        uint32_t ee = exons->end1[i];
        uint32_t cs = exons->cdna_start1[i];
        uint32_t ce = exons->cdna_end1[i];
        uint32_t genomic;

        if (es > ee || cs == 0u || ce < cs) return 0;
        if (cdna_pos < cs || cdna_pos > ce) continue;

        genomic = fwd ? (uint32_t)(es + (cdna_pos - cs))
                      : (uint32_t)(ee - (cdna_pos - cs));
        if (genomic < es || genomic > ee) return 0;
        *genomic_pos_out = genomic;
        if (exon_idx_out != NULL) *exon_idx_out = (uint32_t)i;
        return 1;
    }
    return 0;
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
                                         &coding_start_exon_idx)) {
        return 0;
    }
    if (!duckvep_project_genomic_to_cdna(transcripts, exons, tx_idx,
                                         coding_end_genomic,
                                         &coding_end_cdna,
                                         &dummy_exon_idx)) {
        return 0;
    }
    if (coding_end_cdna < coding_start_cdna) return 0;
    if (cdna < coding_start_cdna || cdna > coding_end_cdna) return 0;

    phase = coding_phase_offset(exons, (size_t)coding_start_exon_idx);
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
