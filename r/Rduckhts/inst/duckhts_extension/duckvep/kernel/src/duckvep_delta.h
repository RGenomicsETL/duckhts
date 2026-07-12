/*
 * duckvep_delta.h — the sequence delta: the FACTS a coding edit produces, and the
 * functions that produce them (INTERNAL).
 *
 * This is VEP's codon/peptide layer (`BaseTranscriptVariationAllele` codon +
 * peptide predicates, and the Haplo `_mutate_sequences` path). It is deliberately
 * its OWN translation unit, separate from the orchestration in duckvep_kernel.c and
 * the rule table in duckvep_effect.c: the consequence engine reads the delta as
 * FACTS (never SO labels — duckvep_effect_eval maps them). All coding edit shapes
 * enter through this module instead of growing the annotation sweep.
 *
 * Borrowed-view discipline: the producers take the SoA views + the borrowed CDS
 * sequence pool directly, NOT the opaque duckvep_model — so this module never
 * depends on the engine's private model layout. No DuckDB/htslib/Arrow; no malloc.
 */
#ifndef DUCKVEP_DELTA_H
#define DUCKVEP_DELTA_H

#include "duckvep_haplotype.h" /* duckvep_haplotype_edit_t (the CDS-edit element) */
#include "duckvep_kernel.h"    /* SoA views, duckvep_sequence_pool_t, variant batch */
#include "duckvep_event.h"

#include <stddef.h>
#include <stdint.h>

#ifndef DUCKVEP_INTERNAL_API
# if defined(__GNUC__) || defined(__clang__)
#  define DUCKVEP_INTERNAL_API __attribute__((visibility("hidden")))
# else
#  define DUCKVEP_INTERNAL_API
# endif
#endif

/* The sequence delta: the FACTS a coding edit produces (== VEP's codon/peptide
 * predicates). The rule table — not this struct — turns these into SO terms, so the
 * codon path never spells SO labels. Unsupported contexts leave `valid` clear and
 * render as an auditable coding_sequence_variant fallback. */
typedef struct duckvep_sequence_delta {
    int32_t cdna_pos, cds_pos, protein_pos; /* 1-based; -1 when absent */
    uint8_t ref_aa, alt_aa;                 /* 1-letter AA; '*' = stop  */
    uint8_t synonymous, missense, stop_gained, stop_lost, stop_retained;
    uint8_t start_lost, frameshift, inframe_deletion, inframe_insertion;
    uint8_t protein_altering;
    uint8_t valid;                          /* 1 when the delta is filled */
} duckvep_sequence_delta_t;

/* One or more projected CDS edits in translation orientation. A single small
 * variant may produce several edits; for example, an MNV with retained internal
 * bases splits into multiple differing islands. A phased haplotype is another
 * N-edit set grouped by sample x phase_set x haplotype x transcript. Both are
 * applied and translated once by the oracle-tested CDS mutation core
 * (duckvep_haplotype_apply_cds_edits):
 * "haplotypes are MNVs at the coding-delta layer; they differ in grouping/flushing,
 * not in consequence logic". The cds-edit element is duckvep_haplotype_edit_t (the C
 * equivalent of the Rust oracle's CdsEdit, predictor.rs CdsEdit).
 *
 * EDIT CONTRACT (matches duckvep_haplotype_edit_t, NOT a new one): alleles are in
 * variant_strand orientation with a per-edit `variant_strand`; the apply helper
 * reverse-complements when variant_strand != transcript_strand (do NOT pre-orient
 * here). cds_start is 1-based; a pure insertion is ref_len==0 inserted BEFORE
 * cds_start. The Rust oracle's variant_to_cds_edit places an insertion AFTER cds_lo
 * (predictor.rs:248-253), so comparisons must account for that +1 convention.
 *
 * Apply/translate split: mutate via duckvep_haplotype_apply_cds_edits, then translate.
 * duckvep_haplotype_translate_cds truncates after the first stop for haplotype protein
 * output. VEP coding predicates instead use the full codon-window translation in
 * duckvep_coding_context_build (predictor.rs build_coding_context:1138). */
typedef struct duckvep_edit_set {
    const duckvep_haplotype_edit_t *edits; /* borrowed; variant_strand orientation     */
    size_t                          count; /* N edits; shared by one allele or haplotype */
} duckvep_edit_set_t;

/* Caller-owned scratch for the CodingContext path. This is internal POD only: no
 * allocation, no ownership transfer, no hidden fallback capacity. */
typedef struct duckvep_delta_scratch {
    duckvep_haplotype_edit_t *edits;
    size_t                    edits_cap;
    uint8_t                  *alt_cds;
    size_t                    alt_cds_cap;
    uint8_t                  *ref_peptide;
    size_t                    ref_peptide_cap;
    uint8_t                  *alt_peptide;
    size_t                    alt_peptide_cap;
} duckvep_delta_scratch_t;

typedef enum duckvep_cds_edit_status {
    DUCKVEP_CDS_EDIT_OK = 0,
    DUCKVEP_CDS_EDIT_INVALID_ARG,
    DUCKVEP_CDS_EDIT_UNSUPPORTED_KIND,
    DUCKVEP_CDS_EDIT_INVALID_EVENT,
    DUCKVEP_CDS_EDIT_OUT_OF_CDS,
    DUCKVEP_CDS_EDIT_NON_CONTIGUOUS,
    DUCKVEP_CDS_EDIT_BUFFER_TOO_SMALL,
    DUCKVEP_CDS_EDIT_INVALID_ALLELE,
    DUCKVEP_CDS_EDIT_REF_MISMATCH
} duckvep_cds_edit_status_t;

typedef enum duckvep_coding_context_status {
    DUCKVEP_CODING_CONTEXT_OK = 0,
    DUCKVEP_CODING_CONTEXT_INVALID_ARG,
    DUCKVEP_CODING_CONTEXT_ALT_CDS_BUFFER_TOO_SMALL,
    DUCKVEP_CODING_CONTEXT_REF_PEPTIDE_BUFFER_TOO_SMALL,
    DUCKVEP_CODING_CONTEXT_ALT_PEPTIDE_BUFFER_TOO_SMALL,
    DUCKVEP_CODING_CONTEXT_OUT_OF_RANGE,
    DUCKVEP_CODING_CONTEXT_INVALID_BASE,
    DUCKVEP_CODING_CONTEXT_REF_MISMATCH,
    DUCKVEP_CODING_CONTEXT_EDIT_ORDER
} duckvep_coding_context_status_t;

/* Consequence-layer coding context over an already projected edit set. Peptides are
 * translated over every complete CDS codon and do NOT truncate after internal stops;
 * '*' is an ordinary peptide byte here. `N` codons translate to `X`, while non-ACGTUN
 * CDS bases are invalid. Changed-codon spans are peptide-diff windows after common
 * prefix/suffix trimming: 0/0 denotes an empty side (for example, a pure peptide
 * insertion on the reference side) or no peptide difference on both sides. These spans
 * are peptide coordinates, not genomic/CDS coordinates. */
typedef struct duckvep_coding_context {
    const uint8_t *ref_cds;     size_t ref_cds_len;
    const uint8_t *alt_cds;     size_t alt_cds_len;
    const uint8_t *ref_peptide; size_t ref_peptide_len;
    const uint8_t *alt_peptide; size_t alt_peptide_len;
    int64_t length_diff;
    uint32_t flags;
    size_t applied_edits;
    uint32_t single_edit_cds_start, single_edit_ref_len, single_edit_alt_len;
    uint8_t has_single_edit;
    uint8_t cds_changed;
    uint32_t ref_first_changed_codon, ref_last_changed_codon;
    uint32_t alt_first_changed_codon, alt_last_changed_codon;
} duckvep_coding_context_t;

typedef enum duckvep_variant_coding_context_status {
    DUCKVEP_VARIANT_CODING_CONTEXT_OK = 0,
    DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ARG,
    DUCKVEP_VARIANT_CODING_CONTEXT_UNSUPPORTED_KIND,
    DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_EVENT,
    DUCKVEP_VARIANT_CODING_CONTEXT_OUT_OF_CDS,
    DUCKVEP_VARIANT_CODING_CONTEXT_NON_CONTIGUOUS,
    DUCKVEP_VARIANT_CODING_CONTEXT_EDIT_BUFFER_TOO_SMALL,
    DUCKVEP_VARIANT_CODING_CONTEXT_ALT_CDS_BUFFER_TOO_SMALL,
    DUCKVEP_VARIANT_CODING_CONTEXT_REF_PEPTIDE_BUFFER_TOO_SMALL,
    DUCKVEP_VARIANT_CODING_CONTEXT_ALT_PEPTIDE_BUFFER_TOO_SMALL,
    DUCKVEP_VARIANT_CODING_CONTEXT_OUT_OF_RANGE,
    DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_ALLELE,
    DUCKVEP_VARIANT_CODING_CONTEXT_INVALID_BASE,
    DUCKVEP_VARIANT_CODING_CONTEXT_REF_MISMATCH,
    DUCKVEP_VARIANT_CODING_CONTEXT_EDIT_ORDER
} duckvep_variant_coding_context_status_t;

typedef enum duckvep_context_delta_status {
    DUCKVEP_CONTEXT_DELTA_OK = 0,
    DUCKVEP_CONTEXT_DELTA_INVALID_ARG,
    DUCKVEP_CONTEXT_DELTA_UNSUPPORTED
} duckvep_context_delta_status_t;

typedef enum duckvep_sequence_delta_route {
    DUCKVEP_DELTA_ROUTE_DIRECT = 0,
    DUCKVEP_DELTA_ROUTE_MNV_CONTEXT,
    DUCKVEP_DELTA_ROUTE_MNV_DIRECT_FALLBACK,
    DUCKVEP_DELTA_ROUTE_DEL_CONTEXT,
    DUCKVEP_DELTA_ROUTE_DEL_DIRECT_FALLBACK,
    DUCKVEP_DELTA_ROUTE_INS_CONTEXT,
    DUCKVEP_DELTA_ROUTE_INS_DIRECT_FALLBACK,
    DUCKVEP_DELTA_ROUTE_INDEL_CONTEXT,
    DUCKVEP_DELTA_ROUTE_INDEL_DIRECT_FALLBACK
} duckvep_sequence_delta_route_t;

/* Project one allele-trimmed small variant into the edit element consumed by the
 * haplotype/CDS mutation helper. This is an internal projection primitive, not a
 * consequence emitter: callers still decide which edit sets are supported.
 * The returned edit keeps allele bytes in genomic/variant orientation and sets
 * variant_strand=+1; duckvep_haplotype_apply_cds_edits() performs reverse-complementing
 * when transcript_strand is -1. Only edits whose differing REF span is wholly CDS and
 * contiguous in transcript CDS coordinates are accepted. */
DUCKVEP_INTERNAL_API duckvep_cds_edit_status_t duckvep_variant_cds_edit_build(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    int8_t                            transcript_strand,
    duckvep_haplotype_edit_t         *edit);

/* Project one small allele into the edit set consumed by the CodingContext and
 * haplotype kernels. Non-MNV shapes emit the one edit produced by
 * duckvep_variant_cds_edit_build(), while fully-CDS
 * equal-length MNVs are split into maximal internal differing islands so retained bases
 * are not edited. Edits are in descending original-CDS order for
 * duckvep_haplotype_apply_cds_edits(). The caller owns `scratch`; `out` borrows it only
 * until the caller reuses the buffer. On every failure, `out` is reset to { NULL, 0 }. */
DUCKVEP_INTERNAL_API duckvep_cds_edit_status_t duckvep_variant_cds_edit_set_build(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    int8_t                            transcript_strand,
    duckvep_haplotype_edit_t         *scratch,
    size_t                            scratch_cap,
    duckvep_edit_set_t               *out);

/* Build a CodingContext over a caller-owned edit set without emitting
 * any SO terms. `ctx` is zeroed on entry and remains all-zero on failure. `alt_cds` is
 * an explicit-length byte buffer; peptide buffers require room for NUL termination and
 * are NUL-terminated on success. `ctx`, `ref_cds`, and all scratch buffers must be
 * non-overlapping; scratch contents are unspecified after failure. */
DUCKVEP_INTERNAL_API duckvep_coding_context_status_t duckvep_coding_context_build(
    const uint8_t               *ref_cds,
    size_t                       ref_cds_len,
    const duckvep_edit_set_t    *edit_set,
    int8_t                       transcript_strand,
    duckvep_codon_table_t        table,
    uint8_t                     *alt_cds_scratch,
    size_t                       alt_cds_cap,
    uint8_t                     *ref_peptide_scratch,
    size_t                       ref_peptide_cap,
    uint8_t                     *alt_peptide_scratch,
    size_t                       alt_peptide_cap,
    duckvep_coding_context_t    *ctx);

/* Compose one small variant directly into a CodingContext: build the
 * edit set, borrow the transcript CDS from `seq`, and call duckvep_coding_context_build.
 * This is not an SO emitter; it is reached by explicit scratch-aware direct calls
 * and by the annotation wrapper for supported MNV/DEL/INS/INDEL contexts.
 * `ctx` is zeroed on entry and remains all-zero on failure. On success, `ctx->ref_cds`
 * borrows from `seq`; alternate CDS and peptides borrow caller scratch. `ctx`, sequence
 * storage, and all scratch buffers must be non-overlapping. */
DUCKVEP_INTERNAL_API duckvep_variant_coding_context_status_t duckvep_variant_coding_context_build(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    int8_t                            transcript_strand,
    duckvep_haplotype_edit_t         *edit_scratch,
    size_t                            edit_scratch_cap,
    uint8_t                          *alt_cds_scratch,
    size_t                            alt_cds_cap,
    uint8_t                          *ref_peptide_scratch,
    size_t                            ref_peptide_cap,
    uint8_t                          *alt_peptide_scratch,
    size_t                            alt_peptide_cap,
    duckvep_coding_context_t         *ctx);

/* Classify a CodingContext into sequence-delta facts. The classifier handles
 * length-preserving substitutions across complete codon windows and guarded
 * single-edit frameshift, in-frame insertion, deletion, and delins contexts.
 * Ambiguous bases, incomplete codons, multi-edit contexts, and cases requiring an
 * incomplete compound consequence return UNSUPPORTED with `delta` invalid. The raw
 * dispatcher has no shape-specific fallback. */
DUCKVEP_INTERNAL_API duckvep_context_delta_status_t duckvep_coding_context_delta_fill(
    const duckvep_coding_context_t *ctx,
    uint64_t                        tx_flags,
    duckvep_sequence_delta_t       *delta);

/* Raw scratch-aware sequence-delta dispatcher for one (variant, transcript) CDS-bucket
 * candidate. Passing scratch exercises the internal CodingContext path for supported
 * MNV and indel shapes. This entry point has no shape-specific fallback. Dispatches on the
 * variant SHAPE (duckvep_variant_kind_t):
 *   - SNV  -> the oracle-tested genomic-codon fast path (project the base to the
 *             coding frame, orient + apply, classify the codon change). Sets the
 *             codon FACTS (never SO labels — duckvep_effect_eval maps them).
 *   - MNV  -> with scratch, build variant CodingContext then classify narrow
 *             one-codon facts plus the two-adjacent-body-codon missense seam; without
 *             scratch, use the direct same-codon path plus the narrow
 *             two-codon non-terminal missense path.
 *   - INS / DEL -> allele-trimmed single-edit frameshift, codon-boundary in-frame
 *             insertion, non-boundary protein-altering insertion, plus codon-aligned
 *             in-frame deletion, only when the affected bases are in non-terminal CDS
 *             body codons.
 *   - INDEL -> allele-trimmed frameshift delins when the replaced REF span is contiguous
 *             non-terminal CDS body and the net length change is not divisible by three;
 *             in-frame/protein-altering delins still fall back to the CDS bucket.
 *   - SV    -> currently leaves the delta INVALID. Structural edit interpretation
 *             belongs here; adapter code must not bypass this dispatcher by calling
 *             the internal CodingContext classifier directly.
 * `delta` is reset (valid = 0) on entry. With the delta invalid, the rule table keeps
 * coding_sequence_variant for the CDS bucket. `seq` is the borrowed CDS pool, or NULL
 * when the model has no sequence (delta always invalid). */
DUCKVEP_INTERNAL_API void duckvep_sequence_delta_fill_with_scratch(
    duckvep_variant_kind_t            kind,
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          pos,
    int8_t                            strand,
    duckvep_delta_scratch_t          *scratch,
    duckvep_sequence_delta_t         *delta);

DUCKVEP_INTERNAL_API void duckvep_sequence_delta_fill(
    duckvep_variant_kind_t            kind,
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          pos,
    int8_t                            strand,
    duckvep_sequence_delta_t         *delta);

/* Production annotation wrapper. It uses caller-owned workspace scratch for the
 * CodingContext path and, while that classifier remains incomplete, uses the direct
 * shape-specific path only when the context leaves `delta` invalid. Route reporting is
 * test instrumentation for making that fallback visible. */
DUCKVEP_INTERNAL_API void duckvep_sequence_delta_fill_for_annotation(
    duckvep_variant_kind_t            kind,
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          pos,
    int8_t                            strand,
    duckvep_delta_scratch_t          *scratch,
    const duckvep_event_t            *prepared_event,
    duckvep_sequence_delta_t         *delta);

DUCKVEP_INTERNAL_API void duckvep_sequence_delta_fill_for_annotation_trace(
    duckvep_variant_kind_t            kind,
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    uint32_t                          pos,
    int8_t                            strand,
    duckvep_delta_scratch_t          *scratch,
    const duckvep_event_t            *prepared_event,
    duckvep_sequence_delta_route_t   *route,
    duckvep_sequence_delta_t         *delta);

#endif /* DUCKVEP_DELTA_H */
