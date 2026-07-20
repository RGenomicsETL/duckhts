/*
 * duckvep_transcript_edit.h — shared transcript edit projection (INTERNAL).
 *
 * Independent HGVS rendering and later phased haplotype execution must
 * consume the same prepared-event and CDS-edit authorities as consequence
 * classification.  The existing duckvep_haplotype_edit_t is intentionally
 * CDS-only; this wider HGVS-facing IR retains transcript coordinates for
 * exonic, intronic, UTR, and non-coding events while borrowing that same
 * optional CDS edit set.  Production consequence annotation still consumes
 * the prepared event and CDS edit helpers directly; it does not construct
 * this wider carrier in its hot loop.
 *
 * No strings are rendered here.  Alleles remain borrowed in genomic-forward
 * orientation and carry transcript_strand explicitly.  A renderer or mutation
 * consumer reverse-complements while reading when transcript_strand is -1.
 */
#ifndef DUCKVEP_TRANSCRIPT_EDIT_H
#define DUCKVEP_TRANSCRIPT_EDIT_H

#include "duckvep_delta.h"

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum duckvep_transcript_edit_status {
    DUCKVEP_TRANSCRIPT_EDIT_OK = 0,
    DUCKVEP_TRANSCRIPT_EDIT_INVALID_ARG,
    DUCKVEP_TRANSCRIPT_EDIT_UNSUPPORTED_KIND,
    DUCKVEP_TRANSCRIPT_EDIT_INVALID_EVENT,
    DUCKVEP_TRANSCRIPT_EDIT_OUTSIDE_TRANSCRIPT,
    DUCKVEP_TRANSCRIPT_EDIT_INVALID_PROJECTION,
    DUCKVEP_TRANSCRIPT_EDIT_OUT_OF_RANGE
} duckvep_transcript_edit_status_t;

/* One genomic base expressed relative to the spliced transcript.  Exonic
 * positions have intron_offset == 0 and cdna_anchor1 is the base itself.
 * Intronic positions retain the closest exon base as cdna_anchor1 and a
 * signed HGVS offset: +N follows that exon base in transcript direction and
 * -N precedes it.  VEP resolves an exactly equidistant intronic base toward
 * the lower-genomic exon on '+' and the higher-genomic exon on '-'; the
 * projector reproduces that strand-dependent tie rule.
 */
typedef struct duckvep_transcript_coordinate {
    uint32_t genomic_pos1;
    uint32_t cdna_anchor1;
    uint32_t exon_idx;       /* absolute model exon index anchoring the coordinate */
    int32_t  intron_offset;  /* zero exonic; otherwise signed transcript offset   */
    uint8_t  exonic;
} duckvep_transcript_coordinate_t;

/* One independent allele projected to one transcript.  `event` retains the
 * unchanged semantic genomic edit. `first`/`last` are the VEP-compatible
 * transcript slice after an endpoint-crossing feature is clipped to the
 * transcript span, in transcript 5'-to-3' order even on the reverse strand.
 * A pure insertion uses the two genomic flanks around its interbase coordinate.
 * `cds_edits` is empty when the event is outside CDS or cannot be represented
 * as a contiguous CDS edit; `cds_status` preserves that reason rather than
 * making transcript projection fail.
 */
typedef struct duckvep_transcript_edit {
    uint32_t variant_idx;
    uint32_t tx_idx;
    duckvep_event_t event;
    /* Transcript-slice coordinates for the minimized edit after VEP's
     * endpoint clipping. The unchanged semantic coordinates remain in event. */
    duckvep_transcript_coordinate_t first;
    duckvep_transcript_coordinate_t last;
    /* Transcript-slice coordinates for VEP's VariationFeature span after the
     * same endpoint clipping. Equal-length uploaded substitutions retain the
     * complete feature span; length-changing alleles use the minimized span.
     * Pure insertions retain their two flanking transcript bases. */
    duckvep_transcript_coordinate_t feature_first;
    duckvep_transcript_coordinate_t feature_last;
    const uint8_t *raw_ref;
    const uint8_t *raw_alt;
    uint16_t raw_ref_length;
    uint16_t raw_alt_length;
    const uint8_t *feature_ref;
    const uint8_t *feature_alt;
    uint16_t feature_ref_length;
    uint16_t feature_alt_length;
    const uint8_t *ref;
    const uint8_t *alt;
    uint16_t ref_length;
    uint16_t alt_length;
    int8_t transcript_strand;
    duckvep_edit_set_t cds_edits;
    uint8_t cds_status; /* duckvep_cds_edit_status_t */
} duckvep_transcript_edit_t;

/* Project an arbitrary genomic base inside a transcript span to the coordinate
 * representation above.  This includes introns; use
 * duckvep_project_genomic_to_cdna() when only exonic bases are admissible.
 */
DUCKVEP_INTERNAL_API duckvep_transcript_edit_status_t
duckvep_project_transcript_coordinate(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    size_t                            tx_idx,
    uint32_t                          genomic_pos1,
    duckvep_transcript_coordinate_t  *out);

/* Build the shared transcript edit for one small allele.  `cds_scratch` is
 * caller-owned and is borrowed through out->cds_edits until reused.  Passing a
 * zero-capacity scratch is valid: transcript coordinates are still returned
 * and cds_status reports whether a CDS edit required capacity.
 */
DUCKVEP_INTERNAL_API duckvep_transcript_edit_status_t
duckvep_transcript_edit_build(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *variants,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    duckvep_haplotype_edit_t         *cds_scratch,
    size_t                            cds_scratch_cap,
    duckvep_transcript_edit_t        *out);

/* Prepared-event form of duckvep_transcript_edit_build(). Annotation cursors
 * already own one validated event per input ALT and must pass that same value
 * to consequence, HGVS, and later haplotype consumers. This entry point never
 * reloads or re-trims REF/ALT. */
DUCKVEP_INTERNAL_API duckvep_transcript_edit_status_t
duckvep_transcript_edit_build_prepared(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *variants,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    const duckvep_event_t            *event,
    duckvep_haplotype_edit_t         *cds_scratch,
    size_t                            cds_scratch_cap,
    duckvep_transcript_edit_t        *out);

#ifdef __cplusplus
}
#endif

#endif /* DUCKVEP_TRANSCRIPT_EDIT_H */
