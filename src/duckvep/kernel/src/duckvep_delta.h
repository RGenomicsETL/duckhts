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
#include "duckvep_compat.h"
#include "duckvep_event.h"

#include <stddef.h>
#include <stdint.h>

#ifndef DUCKVEP_INTERNAL_API
# if (defined(__GNUC__) || defined(__clang__)) && !defined(_WIN32)
#  define DUCKVEP_INTERNAL_API __attribute__((visibility("hidden")))
# else
#  define DUCKVEP_INTERNAL_API
# endif
#endif

/* The sequence delta: the FACTS a coding edit produces (== VEP's codon/peptide
 * predicates). The rule table — not this struct — turns these into SO terms, so the
 * codon path never spells SO labels. Unsupported contexts leave `valid` clear and
 * render as an auditable coding_sequence_variant fallback. */
typedef enum duckvep_nmd_early_cds_fact {
    DUCKVEP_NMD_EARLY_CDS_UNKNOWN = 0,
    DUCKVEP_NMD_EARLY_CDS_ENDS_THROUGH_101,
    DUCKVEP_NMD_EARLY_CDS_ENDS_AFTER_101
} duckvep_nmd_early_cds_fact_t;

/* Executable inclusive threshold in VEP Plugins release/116 NMD.pm. */
#define DUCKVEP_NMD_EARLY_CDS_MAX_END1 101u

typedef struct duckvep_sequence_delta {
    int32_t cdna_pos, cds_pos, protein_pos; /* 1-based; -1 when absent */
    uint8_t ref_aa, alt_aa;                 /* 1-letter AA; '*' = stop  */
    uint8_t synonymous, missense, stop_gained, stop_lost, stop_retained;
    uint8_t start_lost, start_retained, frameshift, inframe_deletion, inframe_insertion;
    uint8_t protein_altering, coding_unknown, partial_codon;
    uint8_t valid;                          /* 1 when the delta is filled */
    uint8_t sequence_status;                /* duckvep_sequence_status_t  */
    /* The one projected-CDS comparison consumed by VEP 116 NMD.pm. This fits
     * in the struct's existing padding; UNKNOWN makes NMD project the event. */
    uint8_t nmd_early_cds_fact;             /* duckvep_nmd_early_cds_fact_t */
} duckvep_sequence_delta_t;

#ifdef __cplusplus
static_assert(sizeof(duckvep_sequence_delta_t) == 32u,
              "duckvep sequence delta no longer fits its hot-path budget");
#else
_Static_assert(sizeof(duckvep_sequence_delta_t) == 32u,
               "duckvep sequence delta no longer fits its hot-path budget");
#endif

/* Encode/apply positive raw sequence predicates carried beside a compact SO
 * mask. Executable VEP witnesses contain emitted term sets that omit a raw
 * predicate still consumed by hgvs_protein, so downstream HGVS must never
 * reconstruct these facts by inspecting consequence_mask. Absence of a bit is
 * not a closed-world false assertion for generalized splice-overlapping edits;
 * apply augments facts independently proved by the HGVS CDS replay. */
DUCKVEP_INTERNAL_API uint32_t duckvep_sequence_delta_consequence_flags(
    const duckvep_sequence_delta_t *delta,
    int                              attempted);

DUCKVEP_INTERNAL_API int duckvep_sequence_delta_apply_consequence_flags(
    uint32_t                   flags,
    duckvep_sequence_delta_t *delta);

/* True when the prepared transcript CDS ends in a one- or two-base codon.
 * VEP's partial_codon predicate derives this from `_translateable_seq` length;
 * it is not reliably equivalent to the cds_end_NF transcript attribute (notably
 * for mitochondrial transcripts). */
DUCKVEP_INTERNAL_API int duckvep_transcript_has_partial_terminal_codon(
    const duckvep_sequence_pool_t *seq,
    size_t                         tx_idx);

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

/* Prepared semantic allele borrowed by CDS projection after upload parsing or
 * HGVS 3-prime placement has already chosen the event coordinates. Alleles
 * are expressed in `variant_strand` orientation. `anchor_ref` is required only
 * for an interbase insertion and names the retained flanking reference base;
 * it is kept separate because the semantic REF length is zero. */
typedef struct duckvep_prepared_cds_allele {
    const duckvep_event_t *event;
    const uint8_t         *ref;
    const uint8_t         *alt;
    const uint8_t         *anchor_ref;
    uint16_t               ref_length;
    uint16_t               alt_length;
    int8_t                 variant_strand;
} duckvep_prepared_cds_allele_t;

/* Borrowed complete spliced-transcript sequence view for a sequence-backed
 * coding transcript. The prepared CDS may contain positive start-phase
 * padding; `phase_offset` hides those synthetic bytes from cDNA addressing.
 * Non-coding transcript sequence storage is not yet part of the model and is
 * reported as MISSING rather than synthesized. */
typedef enum duckvep_transcript_sequence_status {
    DUCKVEP_TRANSCRIPT_SEQUENCE_OK = 0,
    DUCKVEP_TRANSCRIPT_SEQUENCE_MISSING,
    DUCKVEP_TRANSCRIPT_SEQUENCE_INVALID
} duckvep_transcript_sequence_status_t;

typedef struct duckvep_transcript_sequence_view {
    const uint8_t *pre_cds;
    const uint8_t *cds;
    const uint8_t *post_cds;
    size_t pre_cds_len;
    size_t cds_len;
    size_t post_cds_len;
    uint32_t coding_start_cdna;
    uint32_t coding_end_cdna;
    uint8_t phase_offset;
} duckvep_transcript_sequence_view_t;

DUCKVEP_INTERNAL_API duckvep_transcript_sequence_status_t
duckvep_transcript_sequence_open(
    const duckvep_transcript_model_t      *transcripts,
    const duckvep_exon_model_t            *exons,
    const duckvep_sequence_pool_t         *seq,
    size_t                                 tx_idx,
    duckvep_transcript_sequence_view_t    *out);

/* Return one normalized transcript-oriented cDNA byte, or NUL when the
 * coordinate/view is invalid. Ambiguous N is retained so each consumer can
 * apply its own documented policy. */
DUCKVEP_INTERNAL_API char duckvep_transcript_sequence_base(
    const duckvep_transcript_sequence_view_t *sequence,
    uint32_t                                  cdna_pos);

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
    /* A model-backed single edit can expose the alternate CDS as a borrowed
     * view instead of copying and translating the complete transcript. The
     * accessors in duckvep_delta.c preserve the same coding predicates; the
     * materialized fields above remain the compound-edit/Haplosaurus oracle. */
    const uint8_t *single_edit_alt;
    const uint8_t *local_alt_cds;
    const uint8_t *local_ref_peptide;
    const uint8_t *local_alt_peptide;
    size_t local_cds_offset, local_ref_cds_len, local_alt_cds_len;
    size_t local_peptide_offset;
    size_t local_ref_peptide_len, local_alt_peptide_len;
    int8_t transcript_strand;
    int8_t single_edit_variant_strand;
    uint8_t virtual_single_edit;
    /* VEP's genomic _overlaps_stop_codon_cil result for a minimized pure
     * insertion. CDS distance cannot derive this: an intron may lie between
     * the insertion and terminal codon. */
    uint8_t insertion_length_reaches_terminal_stop;
    uint8_t local_ref_unambiguous;
    uint8_t local_alt_unambiguous;
    /* Sparse Ensembl Translation SeqEdits. Positions are one-based and sorted.
     * They are an overlay on the reference peptide only: VEP deliberately
     * leaves the alternate peptide as the raw codon translation. */
    const uint32_t *ref_peptide_edit_position1;
    const uint8_t *ref_peptide_edit_alt;
    size_t ref_peptide_edit_count;
    /* Borrowed complete transcript-oriented 5-prime sequence. A zero length
     * with pre_cds_complete set means the CDS begins at cDNA base one. */
    const uint8_t *pre_cds_bases;
    size_t pre_cds_length;
    const uint8_t *post_cds_bases;
    size_t post_cds_length;
    uint8_t codon_table;
    int64_t length_diff;
    uint32_t flags;
    size_t applied_edits;
    uint32_t single_edit_cds_start, single_edit_ref_len, single_edit_alt_len;
    uint8_t has_single_edit;
    uint8_t cds_changed;
    uint8_t pre_cds_complete;
    uint8_t post_cds_complete;
    uint8_t ref_first_stop_known;
    uint8_t compatibility_profile; /* duckvep_compat_profile_t */
    uint8_t feature_length_relation; /* duckvep_feature_length_relation_t */
    uint32_t ref_first_stop_position1;
    uint32_t ref_first_changed_codon, ref_last_changed_codon;
    uint32_t alt_first_changed_codon, alt_last_changed_codon;
} duckvep_coding_context_t;

/* One authority for the non-obvious single-edit state consumed by both the
 * consequence peptide view and the HGVS peptide view.  The two consumers may
 * deliberately render different VEP-116 views, but they must not rediscover
 * the state with duplicate predicates. */
DUCKVEP_INTERNAL_API int
duckvep_coding_context_is_terminal_partial_insertion(
    const duckvep_coding_context_t *context);

/* The compact consequence sidecar is a closed-world authority for cached
 * start/stop-lost predicates, but frameshift is positive evidence only.  A
 * length-changing splice-overlapping edit may acquire a frameshift during the
 * complete HGVS replay even when the consequence-side flag is absent. */
DUCKVEP_INTERNAL_API int
duckvep_sequence_delta_consequence_flags_complete_for_hgvs(
    const duckvep_coding_context_t *context,
    uint32_t                        flags);

/* VEP's codon-rounded peptide strings for one projected CDS edit. The window
 * borrows `duckvep_coding_context_t`; consumers read residues through
 * duckvep_coding_context_peptide_window_base() so virtual single-edit
 * contexts, materialized haplotypes, terminal partial codons, and Ensembl
 * Translation SeqEdits keep one authority. `peptide_offset` is zero-based in
 * the complete protein; the remaining lengths count one-letter residues. */
typedef struct duckvep_coding_peptide_window {
    size_t peptide_offset;
    size_t reference_span_length;
    size_t ref_nt_length;
    size_t alt_nt_length;
    size_t ref_whole_length;
    size_t alt_whole_length;
    size_t ref_length;
    size_t alt_length;
    uint8_t ref_partial_x;
    uint8_t alt_partial_x;
} duckvep_coding_peptide_window_t;

/* Open the exact codon-rounded local peptide window consumed by VEP 116 for a
 * supported single edit. Consequence predicates currently use it for
 * length-changing alleles; HGVS also uses it for equal-length edits. */
DUCKVEP_INTERNAL_API int duckvep_coding_context_peptide_window_open(
    const duckvep_coding_context_t  *ctx,
    duckvep_coding_peptide_window_t *window);

/* Read one local residue from an opened window. `alternate` selects the
 * altered peptide when nonzero. Returns NUL for invalid or out-of-range input.
 * Reference reads include sparse Ensembl Translation SeqEdits; alternate
 * reads intentionally retain VEP's raw codon translation. */
DUCKVEP_INTERNAL_API uint8_t duckvep_coding_context_peptide_window_base(
    const duckvep_coding_context_t        *ctx,
    const duckvep_coding_peptide_window_t *window,
    int                                    alternate,
    size_t                                 index);

/* Read one residue from a complete coding context with the same reference
 * edit and virtual-single-edit semantics as the local window above. */
DUCKVEP_INTERNAL_API uint8_t duckvep_coding_context_peptide_base(
    const duckvep_coding_context_t *ctx,
    int                             alternate,
    size_t                          position0);

/* Read one normalized nucleotide from the reference or alternate CDS view.
 * This preserves virtual single-edit storage and is the nucleotide authority
 * for consumers, such as VEP-compatible frameshift HGVS, that translate
 * across the CDS-to-3-prime-UTR join without materializing another sequence. */
DUCKVEP_INTERNAL_API char duckvep_coding_context_cds_base(
    const duckvep_coding_context_t *ctx,
    int                             alternate,
    size_t                          position0);

/* Scan one alternate CDS prefix followed by a caller-owned transcript suffix
 * and report the first translated stop. This is the sequential authority for
 * consumers that need a stop position without materializing the complete
 * alternate peptide. Ambiguous N codons translate to X and do not stop the
 * scan. `stop_position0` is a zero-based peptide coordinate when `found` is
 * set. */
DUCKVEP_INTERNAL_API duckvep_coding_context_status_t
duckvep_coding_context_first_alt_stop(
    const duckvep_coding_context_t *ctx,
    size_t                          cds_prefix_length,
    const uint8_t                  *suffix,
    size_t                          suffix_length,
    size_t                         *stop_position0,
    int                            *found);

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
    DUCKVEP_CONTEXT_DELTA_UNSUPPORTED,
    DUCKVEP_CONTEXT_DELTA_MISSING_TRANSCRIPT_TAIL,
    DUCKVEP_CONTEXT_DELTA_MISSING_TRANSCRIPT_FLANK
} duckvep_context_delta_status_t;

typedef enum duckvep_sequence_delta_route {
    DUCKVEP_DELTA_ROUTE_DIRECT = 0,
    DUCKVEP_DELTA_ROUTE_SIMPLE_INDEL,
    DUCKVEP_DELTA_ROUTE_SUBSTITUTION_CONTEXT,
    DUCKVEP_DELTA_ROUTE_DEL_CONTEXT,
    DUCKVEP_DELTA_ROUTE_INS_CONTEXT,
    DUCKVEP_DELTA_ROUTE_INDEL_CONTEXT,
    DUCKVEP_DELTA_ROUTE_BOUNDARY_CONTEXT,
    DUCKVEP_DELTA_ROUTE_UPLOADED_FEATURE_CONTEXT
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

/* Project one already prepared semantic allele into the common CDS edit IR.
 * This is the shared entry point for uploaded variants, HGVS-shifted alleles,
 * and later phased transcript edit sets. It validates the shifted REF against
 * the prepared CDS and does not trim, normalize, or move the event. */
DUCKVEP_INTERNAL_API duckvep_cds_edit_status_t
duckvep_cds_edit_build_prepared_allele(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    size_t                            tx_idx,
    int8_t                            transcript_strand,
    const duckvep_prepared_cds_allele_t *allele,
    uint32_t                          exon_hint,
    duckvep_haplotype_edit_t         *edit);

/* Reproduce VEP 116's independent-event outer-CDS replacement for a literal
 * feature whose genomic span contains one or more introns but whose two
 * endpoints map to CDS. VEP replaces the contiguous CDS range between the
 * mapped endpoints with the complete feature ALT. This is intentionally not
 * the phased edit-set projector: a phased caller must preserve every exon and
 * intron segment instead of collapsing the feature to one CDS edit. */
DUCKVEP_INTERNAL_API duckvep_cds_edit_status_t
duckvep_compat_vep116_outer_cds_edit_build(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    size_t                            tx_idx,
    int8_t                            transcript_strand,
    const duckvep_event_t            *event,
    const uint8_t                    *alternate,
    uint32_t                          alternate_length,
    int8_t                            alternate_strand,
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

/* Prepared-event form of duckvep_variant_cds_edit_set_build(). The caller has
 * already interpreted the uploaded REF/ALT pair into `event`; this function
 * must consume that exact geometry instead of trimming the raw alleles again.
 * `exon_hint` is an absolute model exon index, or UINT32_MAX when unavailable.
 * The output ownership and failure contract are identical to the wrapper
 * above. */
DUCKVEP_INTERNAL_API duckvep_cds_edit_status_t
duckvep_variant_cds_edit_set_build_prepared(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    int8_t                            transcript_strand,
    const duckvep_event_t            *event,
    uint32_t                          exon_hint,
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

/* Build and enrich a CodingContext from an already projected edit set and the
 * immutable model. One length-changing edit may retain the virtual local
 * sequence representation; multi-edit/haplotype sets materialize once. The
 * resulting context borrows model sequence, edit alleles, and caller scratch.
 * `event` is optional for phased sets; an independent insertion supplies it so
 * the VEP terminal-stop guard remains available. */
DUCKVEP_INTERNAL_API duckvep_variant_coding_context_status_t
duckvep_model_coding_context_build(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    size_t                            tx_idx,
    int8_t                            transcript_strand,
    const duckvep_event_t            *event,
    const duckvep_edit_set_t         *edit_set,
    uint8_t                          *alt_cds_scratch,
    size_t                            alt_cds_cap,
    uint8_t                          *ref_peptide_scratch,
    size_t                            ref_peptide_cap,
    uint8_t                          *alt_peptide_scratch,
    size_t                            alt_peptide_cap,
    duckvep_coding_context_t         *ctx);

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

/* Prepared-event form used when the annotation cursor has already interpreted
 * REF/ALT once. `exon_hint` is an absolute model exon index or UINT32_MAX.
 * The helper may retain the model-backed virtual single-edit representation;
 * all consumers must use the coding-context peptide accessors rather than
 * assuming complete alternate buffers were materialized. */
DUCKVEP_INTERNAL_API duckvep_variant_coding_context_status_t
duckvep_variant_coding_context_build_prepared(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *v,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    int8_t                            transcript_strand,
    const duckvep_event_t            *event,
    uint32_t                          exon_hint,
    duckvep_haplotype_edit_t         *edit_scratch,
    size_t                            edit_scratch_cap,
    uint8_t                          *alt_cds_scratch,
    size_t                            alt_cds_cap,
    uint8_t                          *ref_peptide_scratch,
    size_t                            ref_peptide_cap,
    uint8_t                          *alt_peptide_scratch,
    size_t                            alt_peptide_cap,
    duckvep_coding_context_t         *ctx);

/* VEP maps an equal-length uploaded multi-base feature as one peptide window,
 * even when semantic allele trimming leaves a smaller substitution.  This is
 * the shared representation-sensitive path used by both consequence facts and
 * independent-event HGVSp.  NOT_APPLICABLE means the caller should use the
 * ordinary semantic edit set.  DELTA_ONLY means the uploaded feature was
 * authoritative but no reusable peptide context exists (for example, a
 * partial terminal codon or an explicit sequence/projection failure).
 * CONTEXT_READY returns both the exact consequence delta and its borrowed
 * coding context; the context remains valid until caller scratch is reused. */
typedef enum duckvep_feature_substitution_result {
    DUCKVEP_FEATURE_SUBSTITUTION_NOT_APPLICABLE = 0,
    DUCKVEP_FEATURE_SUBSTITUTION_DELTA_ONLY = 1,
    DUCKVEP_FEATURE_SUBSTITUTION_CONTEXT_READY = 2
} duckvep_feature_substitution_result_t;

DUCKVEP_INTERNAL_API duckvep_feature_substitution_result_t
duckvep_feature_substitution_context_fill(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_sequence_pool_t    *seq,
    const duckvep_variant_batch_t    *variants,
    uint32_t                          variant_idx,
    size_t                            tx_idx,
    int8_t                            transcript_strand,
    duckvep_delta_scratch_t          *scratch,
    const duckvep_event_t            *prepared_event,
    uint32_t                          exon_hint,
    duckvep_coding_context_t         *context_out,
    duckvep_sequence_delta_t         *delta_out);

/* Classify a CodingContext into sequence-delta facts. The classifier handles
 * length-preserving substitutions across complete codon windows and guarded
 * single-edit frameshift, in-frame insertion, deletion, and delins contexts.
 * Ambiguous bases, incomplete codons, length-changing multi-edit contexts, and cases
 * requiring an incomplete compound consequence return UNSUPPORTED with `delta`
 * invalid. Length-preserving multi-edit substitutions use the same complete-CDS
 * comparison as one uploaded MNV. The raw dispatcher has no shape-specific fallback. */
DUCKVEP_INTERNAL_API duckvep_context_delta_status_t duckvep_coding_context_delta_fill(
    const duckvep_coding_context_t *ctx,
    uint64_t                        tx_flags,
    duckvep_sequence_delta_t       *delta);

/* Pure-C reference dispatcher for one (variant, transcript) CDS-bucket candidate. The
 * property suite calls this directly; production annotation calls
 * duckvep_sequence_delta_fill_for_annotation instead. Passing scratch exercises the shared
 * CodingContext interpreter. Omitting scratch selects deliberately narrow direct reference
 * classifiers so properties can compare two implementations; a failed production context
 * is never retried here. Dispatches on the variant SHAPE (duckvep_variant_kind_t):
 *   - SNV  -> the oracle-tested genomic-codon fast path (project the base to the
 *             coding frame, orient + apply, classify the codon change). Sets the
 *             codon FACTS (never SO labels — duckvep_effect_eval maps them).
 *   - MNV  -> with scratch, build and classify the CodingContext; without scratch,
 *             use the direct same-codon path plus a narrow two-codon reference path.
 *   - INS / DEL -> allele-trimmed single-edit frameshift, codon-boundary in-frame
 *             insertion, non-boundary protein-altering insertion, plus codon-aligned
 *             in-frame deletion, only when the affected bases are in non-terminal CDS
 *             body codons.
 *   - INDEL -> the no-scratch reference accepts allele-trimmed frameshift delins only when
 *             the replaced REF span is contiguous non-terminal CDS body and the net length
 *             change is not divisible by three.
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
 * CodingContext path. Projection or classification failure remains reason-coded;
 * annotation never retries through a shape-specific consequence classifier. Route
 * reporting is internal test instrumentation. */
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
    uint32_t                          classified_region_mask,
    uint32_t                          exon_hint,
    duckvep_sequence_delta_route_t   *route,
    duckvep_sequence_delta_t         *delta);

/* Synchronous traced form used by consumers that must reuse the exact coding
 * interpretation before worker scratch is recycled. `context_out` borrows the
 * supplied scratch and is valid only until the next delta build on that
 * workspace. A successful context is reported independently of whether the
 * consequence delta is valid, because HGVS and later phased edit-set consumers
 * may still need the translated state. Ordinary annotation calls the wrapper
 * above and pays no trace-copy cost. */
DUCKVEP_INTERNAL_API void duckvep_sequence_delta_fill_for_annotation_observed(
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
    uint32_t                          classified_region_mask,
    uint32_t                          exon_hint,
    duckvep_sequence_delta_route_t   *route,
    duckvep_sequence_delta_t         *delta,
    duckvep_coding_context_t         *context_out,
    duckvep_variant_coding_context_status_t *context_status_out);

#endif /* DUCKVEP_DELTA_H */
