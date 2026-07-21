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

/* Resolve the spliced-transcript cDNA coordinates of the first and last CDS
 * bases. Prepared models provide these values directly; borrowed-view callers
 * may omit the cache and use the same exon projection fallback. The optional
 * exon/phase outputs identify the translation-start exon and its positive
 * Ensembl phase padding. This is the shared coordinate authority for coding
 * projection, HGVS c. numbering, and complete-transcript sequence views. */
int duckvep_project_coding_cdna_bounds(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                         *coding_start_cdna_out,
    uint32_t                         *coding_end_cdna_out,
    uint32_t                         *coding_start_exon_idx_out,
    uint8_t                          *phase_offset_out);

/* Reproduce VariationEffect::_overlaps_start_codon for the unshifted complete
 * VariationFeature. The mapper preserves insertion geometry as start=end+1;
 * this is observably different from asking whether a normalized CDS edit
 * touches CDS position one. Mapper::map_insert drops an intronic Gap flank,
 * while an ordinary feature requires both endpoint coordinates. A
 * cds_start_NF transcript always returns false. */
int duckvep_project_feature_overlaps_start_codon_unshifted(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event);

/* Reproduce the cached `_pre_consequence_predicates->{coding}` admission for
 * a complete small-variant feature. This is broader than requiring both
 * feature endpoints in CDS: one mapped CDS segment is sufficient, including
 * one flank of an insertion at a UTR/CDS edge. */
int duckvep_project_feature_has_coding_precondition_unshifted(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event);

/* Reproduce BaseTranscriptVariation::translation_coords for the complete
 * VariationFeature after a signed genomic HGVS shift. Both the first and last
 * mapper items must be coding coordinates; internal introns are permitted.
 * Mapper::map_insert removes an intronic/transcript-external flank before the
 * CDS conversion, but retains an exonic UTR flank as a terminal Gap. */
int duckvep_project_complete_feature_translation_bounds(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    uint32_t                          genomic_shift,
    uint32_t                         *cds_start_out,
    uint32_t                         *cds_end_out);

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
