/*
 * duckvep_projection.h — genomic -> cDNA -> CDS -> peptide coordinate
 * projection (INTERNAL).
 *
 * Pure coordinate arithmetic over the borrowed transcript/exon SoA model: no
 * DuckDB, htslib, sequence access, or allocation. This is the bridge between
 * structural classification ("this candidate is in CDS") and codon-level
 * consequence classification ("which codon/protein residue is affected?").
 *
 * Implements the coordinate convention audited from Ensembl/VEP's
 * TranscriptMapper path (`BaseTranscriptVariation.pm`) and VEP's GFF/GTF loader
 * (`BaseGXF.pm`). This unit is still a structural arithmetic kernel, not proof
 * of full VEP differential conformance. In particular, CDS coordinates include
 * the first coding exon's positive Ensembl phase offset (cds_start_NF /
 * incomplete start handling): phase 0 starts at CDS 1, phase 1 starts at CDS 2,
 * phase 2 starts at CDS 3.
 *
 * Exons are expected in TRANSCRIPT ORDER in the model slice, with cDNA starts
 * and ends precomputed by the adapter. For '-' transcripts this means genomic
 * high-to-low order, while exon start1/end1 remain normal genomic coordinates.
 */
#ifndef DUCKVEP_PROJECTION_H
#define DUCKVEP_PROJECTION_H

#include "duckvep_kernel.h"
#include "duckvep_event.h"

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct duckvep_coding_projection {
    uint32_t cdna_pos;        /* 1-based cDNA coordinate */
    uint32_t cds_pos;         /* 1-based CDS coordinate, phase-adjusted */
    uint32_t protein_pos;     /* 1-based peptide coordinate */
    uint32_t codon_start_cds; /* 1-based CDS coordinate of codon start */
    uint32_t exon_idx;        /* absolute exon index in duckvep_exon_model_t */
    uint8_t  codon_offset;    /* 0,1,2 within codon using phase-adjusted CDS */
    uint8_t  phase_offset;    /* positive Ensembl phase on the translation-start exon */
} duckvep_coding_projection_t;

/* Map a genomic base to cDNA. Returns 1 if the position is exonic for tx_idx,
 * else 0. `exon_idx_out` is optional. */
int duckvep_project_genomic_to_cdna(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          genomic_pos,
    uint32_t                         *cdna_pos_out,
    uint32_t                         *exon_idx_out);

/* Inverse of genomic_to_cdna for exonic positions. Returns 1 if `cdna_pos` is
 * inside an exon for tx_idx, else 0. */
int duckvep_project_cdna_to_genomic(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          cdna_pos,
    uint32_t                         *genomic_pos_out,
    uint32_t                         *exon_idx_out);

/* Full coding projection for a single genomic base. Returns 1 only when the
 * position is exonic AND inside the transcript's coding interval; otherwise 0.
 * The result is ready to identify the affected codon/protein residue. */
int duckvep_project_coding_base(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          genomic_pos,
    duckvep_coding_projection_t      *out);

/* Return VEP's BaseTranscriptVariation::translation_start projection for a
 * non-insertion feature. The first mapper item is transcript-oriented: a
 * leading UTR/intron Gap leaves translation_start undefined even when a later
 * part of the uploaded span overlaps CDS. */
int duckvep_project_feature_translation_start(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    duckvep_coding_projection_t      *out);

/* VEP partial_codon asks whether the first affected CDS coordinate belongs to
 * the incomplete codon at the end of the translateable CDS. */
int duckvep_cds_position_is_partial_codon(
    size_t   cds_length,
    uint32_t cds_position1);

/* Project one prepared small-variant event to its inclusive CDS-coordinate
 * bounds. Every changed reference base must map contiguously; splice-crossing
 * spans therefore fail instead of inventing a coding range. A pure insertion
 * returns the CDS boundary before which its payload is inserted. */
int duckvep_project_event_to_cds(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    uint32_t                         *cds_start_out,
    uint32_t                         *cds_end_out);

/* Project VEP's VariationFeature endpoints to CDS coordinates. Unlike the
 * semantic edit projector above, this intentionally retains unchanged bases
 * in equal-length uploaded alleles and permits an intron between two mapped
 * endpoints. BaseTranscriptVariation::cds_coords makes the first or last
 * mapper Gap undefined; this helper therefore requires both endpoints to map.
 * A pure insertion reuses the one-flank insertion-boundary projector, including
 * exon edges, and returns VEP's reversed TranscriptVariation CDS range
 * (`cds_start > cds_end`). NMD.pm reads that lower `cds_end`; it is not the
 * expanded allele range rendered in VEP JSON. */
int duckvep_project_feature_to_cds(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    uint32_t                         *cds_start_out,
    uint32_t                         *cds_end_out);

#ifdef __cplusplus
}
#endif

#endif /* DUCKVEP_PROJECTION_H */
