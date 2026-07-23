/*
 * duckvep_hgvs.h — typed HGVS facts over the shared transcript edit (INTERNAL).
 *
 * This layer owns HGVS coordinate numbering and transcript-oriented allele
 * access. It does not allocate or render identifiers/strings. Sequence-aware
 * operation refinement (dup/inv/repeat), 3-prime shifting, and protein facts
 * build on these numeric facts instead of re-projecting the uploaded allele.
 */
#ifndef DUCKVEP_HGVS_H
#define DUCKVEP_HGVS_H

#include "duckvep_transcript_edit.h"

#include <stddef.h>
#include <stdint.h>

/* VEP 116 TranscriptVariationAllele::_genomic_shift fetches at most this many
 * genomic bases on either side. Keep one authority for interval construction,
 * reference lookup, and transcript-oriented shifting. */
#define DUCKVEP_HGVS_SHIFT_LIMIT 1000u

#ifdef __cplusplus
extern "C" {
#endif

typedef enum duckvep_hgvs_status {
    DUCKVEP_HGVS_OK = 0,
    DUCKVEP_HGVS_INVALID_ARG,
    DUCKVEP_HGVS_INVALID_PROJECTION,
    DUCKVEP_HGVS_OUT_OF_RANGE,
    DUCKVEP_HGVS_INVALID_ALLELE,
    DUCKVEP_HGVS_MISSING_REFERENCE,
    DUCKVEP_HGVS_REFERENCE_MISMATCH,
    DUCKVEP_HGVS_MISSING_PEPTIDE,
    /* VEP intentionally returns undef for this otherwise valid topology. */
    DUCKVEP_HGVS_NOT_APPLICABLE,
    DUCKVEP_HGVS_UNSUPPORTED_PROTEIN,
    DUCKVEP_HGVS_UNSUPPORTED_EDIT,
    DUCKVEP_HGVS_BUFFER_TOO_SMALL
} duckvep_hgvs_status_t;

typedef enum duckvep_hgvs_numbering {
    DUCKVEP_HGVS_NUMBERING_C = 1,
    DUCKVEP_HGVS_NUMBERING_N = 2
} duckvep_hgvs_numbering_t;

/* HGVS c./n. coordinate atom before text materialization.
 *
 *   N:      base is the positive transcript coordinate.
 *   C:      negative bases are 5-prime UTR; positive bases are CDS-relative.
 *   C_STAR: base is the positive 3-prime UTR offset. A zero base preserves
 *           VEP's special stop-codon-adjacent intron form: c.* plus the
 *           signed intron_offset (for example c.*1).
 */
typedef enum duckvep_hgvs_coordinate_kind {
    DUCKVEP_HGVS_COORDINATE_N = 1,
    DUCKVEP_HGVS_COORDINATE_C = 2,
    DUCKVEP_HGVS_COORDINATE_C_STAR = 3
} duckvep_hgvs_coordinate_kind_t;

typedef struct duckvep_hgvs_coordinate {
    int64_t base;
    int32_t intron_offset;
    uint8_t kind; /* duckvep_hgvs_coordinate_kind_t */
} duckvep_hgvs_coordinate_t;

/* Initial edit shape. REPLACEMENT is the VEP delins default for any non-empty
 * multi-base or unequal-length REF/ALT pair. Sequence-aware refinement may
 * later prove inversion, duplication, or repeat syntax without changing the
 * underlying transcript edit. */
typedef enum duckvep_hgvs_dna_shape {
    DUCKVEP_HGVS_DNA_SUBSTITUTION = 1,
    DUCKVEP_HGVS_DNA_DELETION = 2,
    DUCKVEP_HGVS_DNA_INSERTION = 3,
    DUCKVEP_HGVS_DNA_REPLACEMENT = 4,
    DUCKVEP_HGVS_DNA_INVERSION = 5,
    DUCKVEP_HGVS_DNA_DUPLICATION = 6,
    DUCKVEP_HGVS_DNA_REPEAT = 7
} duckvep_hgvs_dna_shape_t;

typedef struct duckvep_hgvs_dna_fact {
    duckvep_hgvs_coordinate_t first;
    duckvep_hgvs_coordinate_t last;
    /* Semantic genomic placement after optional VEP 3-prime shifting.
     * These fields remain the inserted/deleted edit location when `first` and
     * `last` are later rewritten to the copied source span for `dup` syntax. */
    uint32_t placed_start1;
    uint32_t placed_end1;
    uint32_t placed_insertion_boundary0;
    const uint8_t *ref;
    const uint8_t *alt;
    uint16_t ref_length;
    uint16_t alt_length;
    uint32_t repeat_count;
    int32_t shift_offset;
    int8_t transcript_strand;
    uint8_t numbering; /* duckvep_hgvs_numbering_t */
    uint8_t shape;     /* duckvep_hgvs_dna_shape_t */
} duckvep_hgvs_dna_fact_t;

/* Borrowed forward-reference window. `start1` is the one-based genomic
 * coordinate of bases[0]. This keeps faidx ownership outside the host-neutral
 * kernel. HGVS shifting and reference lookup deliberately use separate views:
 * VEP's exact constrained +/-1000 genomic slice controls the shift, while a
 * wider lookup view may validate retained VCF REF padding and identify a
 * copied duplication source without changing that shift slice. */
typedef struct duckvep_hgvs_reference_window {
    const uint8_t *bases;
    size_t length;
    uint32_t start1;
    uint16_t chrom_id;
} duckvep_hgvs_reference_window_t;

/* Derive the exact genomic slice fetched by VEP 116 _genomic_shift.  A
 * deletion uses [start-1000, end+1000].  A pure insertion retains VEP's
 * reversed VariationFeature geometry (start = interbase+1, end = interbase),
 * so its slice is [interbase+1-1000, interbase+1000].  Both ends are clipped
 * to the declared sequence-region length.  Keeping this arithmetic in the
 * host-neutral kernel prevents the faidx adapter and shift oracle from
 * disagreeing by one reference base. */
DUCKVEP_INTERNAL_API duckvep_hgvs_status_t
duckvep_hgvs_genomic_search_interval(
    const duckvep_event_t *event,
    uint32_t               sequence_length,
    uint32_t              *start1_out,
    uint32_t              *end1_out);

/* Derive the bounded genomic interval needed for reference assertions and
 * operation refinement around the exact VEP shift slice. The lookup interval
 * also covers the complete uploaded REF, a shifted partially-overlapping
 * feature, and both possible adjacent duplication sources for a literal
 * insertion. It must not be substituted for the exact shift interval above. */
DUCKVEP_INTERNAL_API duckvep_hgvs_status_t
duckvep_hgvs_reference_fetch_interval(
    const duckvep_event_t *event,
    uint32_t               sequence_length,
    uint32_t              *start1_out,
    uint32_t              *end1_out);

/* Validate the complete uploaded VCF REF span against a fetched genomic
 * reference window. This is deliberately separate from the minimized edit:
 * insertion anchors and retained padding are part of the uploaded assertion
 * even though they are absent from the semantic differing region. Adapters
 * should call this once per ALT before reusing the window across transcripts. */
DUCKVEP_INTERNAL_API duckvep_hgvs_status_t
duckvep_hgvs_uploaded_reference_validate(
    const duckvep_hgvs_reference_window_t *reference,
    const duckvep_event_t                 *event,
    const uint8_t                         *raw_ref,
    size_t                                 raw_ref_length);

/* Typed protein edit after VEP-compatible peptide clipping and 3-prime
 * placement.  The fact borrows `context`; `window` is a value snapshot whose
 * residues remain accessible only while the context and its caller-owned
 * scratch buffers remain alive.  ref/alt offsets are relative to `window`.
 * `termination_distance` is the HGVS Ter distance for frameshift/extension
 * syntax when termination_known is set; otherwise the renderer writes `?`.
 */
typedef enum duckvep_hgvs_protein_shape {
    DUCKVEP_HGVS_PROTEIN_EQUAL = 1,
    DUCKVEP_HGVS_PROTEIN_SUBSTITUTION = 2,
    DUCKVEP_HGVS_PROTEIN_DELETION = 3,
    DUCKVEP_HGVS_PROTEIN_INSERTION = 4,
    DUCKVEP_HGVS_PROTEIN_DELINS = 5,
    DUCKVEP_HGVS_PROTEIN_DUPLICATION = 6,
    DUCKVEP_HGVS_PROTEIN_FRAMESHIFT = 7,
    DUCKVEP_HGVS_PROTEIN_START_LOST = 8,
    DUCKVEP_HGVS_PROTEIN_EXTENSION = 9
} duckvep_hgvs_protein_shape_t;

typedef struct duckvep_hgvs_protein_fact {
    const duckvep_coding_context_t *context;
    duckvep_coding_peptide_window_t window;
    size_t ref_offset;
    size_t ref_length;
    size_t alt_offset;
    size_t alt_length;
    size_t ref_rotation;
    size_t alt_rotation;
    uint32_t first_position1;
    uint32_t last_position1;
    uint32_t termination_distance;
    uint8_t reference_first;
    uint8_t reference_last;
    uint8_t alternate_first;
    uint8_t shape; /* duckvep_hgvs_protein_shape_t */
    uint8_t termination_known;
    uint8_t extends_stop;
    uint8_t compatibility_profile; /* duckvep_compat_profile_t */
    /* start_lost after an insertion uses VEP's two surrounding reference
     * residues, which need not be part of the local changed-peptide window. */
    uint8_t start_lost_flanking;
} duckvep_hgvs_protein_fact_t;

/* Convert one transcript coordinate into VEP 116 c./n. numbering. */
DUCKVEP_INTERNAL_API duckvep_hgvs_status_t
duckvep_hgvs_coordinate_from_transcript(
    const duckvep_transcript_model_t      *transcripts,
    const duckvep_exon_model_t            *exons,
    size_t                                 tx_idx,
    const duckvep_transcript_coordinate_t *coordinate,
    duckvep_hgvs_coordinate_t             *out);

/* Build unshifted independent-event DNA facts from one shared transcript edit.
 * No raw allele interpretation, genomic projection, or sequence lookup occurs
 * here. */
DUCKVEP_INTERNAL_API duckvep_hgvs_status_t duckvep_hgvs_dna_fact_build(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_transcript_edit_t  *edit,
    duckvep_hgvs_dna_fact_t          *out);

/* VEP 116 --shift_hgvs semantics used for ordinary Ensembl transcripts.
 * TranscriptVariationAllele::_return_3prime first scans a genomic reference
 * window in transcript-strand direction, with the hard-coded 1000-base cap,
 * then projects the shifted event back to c./n. coordinates. The input edit
 * is not mutated; the fact retains the shifted semantic genomic placement so
 * HGVSp can rebuild the same rotated edit. Other edit shapes are returned
 * unchanged without requiring a reference window. RefSeq RNA-edited
 * transcripts require VEP's separate transcript-sequence fallback and are
 * not handled by this entry point. */
DUCKVEP_INTERNAL_API duckvep_hgvs_status_t
duckvep_hgvs_dna_fact_build_genomic_shifted(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_hgvs_reference_window_t *reference,
    const duckvep_transcript_edit_t  *edit,
    duckvep_hgvs_dna_fact_t          *out);

/* Variant of the shift builder with explicit reference authorities. The
 * shift view must be exactly duckvep_hgvs_genomic_search_interval(); the
 * lookup view may be wider and is used for REF validation, clamped feature
 * materialization, and adjacent duplication-source comparison. */
DUCKVEP_INTERNAL_API duckvep_hgvs_status_t
duckvep_hgvs_dna_fact_build_genomic_shifted_with_lookup(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_hgvs_reference_window_t *shift_reference,
    const duckvep_hgvs_reference_window_t *lookup_reference,
    const duckvep_transcript_edit_t  *edit,
    duckvep_hgvs_dna_fact_t          *out);

/* Reproduce hgvs_protein's mapper precondition after HGVS 3-prime placement.
 * VEP requires both the first and last genomic2pep mapper items to be peptide
 * coordinates; a leading or trailing Gap makes protein HGVS non-applicable
 * even when some part of the feature overlaps CDS. `defined_out` is zero for
 * that ordinary state and one when protein reconstruction may proceed. */
DUCKVEP_INTERNAL_API duckvep_hgvs_status_t
duckvep_hgvs_protein_coordinates_defined(
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t       *exons,
    const duckvep_transcript_edit_t  *edit,
    const duckvep_hgvs_dna_fact_t    *fact,
    int                              *defined_out);

/* Project the semantic insertion/deletion retained by a shifted DNA fact into
 * the shared CDS edit IR. The rotated allele is materialized once in genomic-
 * forward orientation in caller-owned `allele_scratch`; an insertion also
 * uses one trailing byte for its retained VCF reference anchor. The returned
 * edit borrows `allele_scratch` until it is reused. `required_out` receives
 * the required byte count even when capacity is insufficient.
 *
 * This is the composition point that makes c./n. placement and protein
 * mutation consume one VEP-shifted edit. It does not left-align, minimize, or
 * reinterpret the uploaded allele. */
DUCKVEP_INTERNAL_API duckvep_hgvs_status_t
duckvep_hgvs_shifted_cds_edit_build(
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
    duckvep_haplotype_edit_t               *out);

/* Return one normalized A/C/G/T allele byte in transcript 5'-to-3' order.
 * `alternate` selects ALT when nonzero and REF otherwise. */
DUCKVEP_INTERNAL_API duckvep_hgvs_status_t duckvep_hgvs_dna_base(
    const duckvep_hgvs_dna_fact_t *fact,
    int                             alternate,
    size_t                          index,
    uint8_t                        *base_out);

/* Render the unshifted c./n. DNA suffix (for example c.76A>G or
 * n.12_13insAC) into a caller-owned buffer. The transcript stable identifier
 * is deliberately a relational projection outside the kernel. `required_out`
 * receives the byte count excluding NUL, including when capacity is zero or
 * truncation occurs. A missing or undersized output buffer returns
 * DUCKVEP_HGVS_BUFFER_TOO_SMALL; callers may reserve required_out + 1 bytes
 * and retry. */
DUCKVEP_INTERNAL_API duckvep_hgvs_status_t duckvep_hgvs_dna_render_basic(
    const duckvep_hgvs_dna_fact_t *fact,
    char                           *buffer,
    size_t                          capacity,
    size_t                         *required_out);

/* Build an independent-event HGVSp fact from the same coding context and
 * predicate facts used by consequence classification.  No SO term is used as
 * an input.  Materialized multi-edit contexts are accepted, which is the
 * later Haplosaurus reuse point; virtual and materialized single-edit
 * contexts share the codon-rounded peptide-window helper above. */
DUCKVEP_INTERNAL_API duckvep_hgvs_status_t duckvep_hgvs_protein_fact_build(
    const duckvep_coding_context_t *context,
    const duckvep_sequence_delta_t *delta,
    duckvep_hgvs_protein_fact_t    *out);

/* VEP deletes its HGVS shift hash before the late frameshift stop search.
 * Adapters that replay the restored original allele at the original CDS
 * coordinate use this helper to replace the provisional Ter distance without
 * rebuilding the already cached frameshift residues. */
DUCKVEP_INTERNAL_API duckvep_hgvs_status_t
duckvep_hgvs_protein_frameshift_termination_replay(
    const duckvep_coding_context_t *late_context,
    duckvep_hgvs_protein_fact_t    *fact);

/* Build the subset of HGVSp facts completely represented by the consequence
 * result's proven single-residue sidecar. This avoids reconstructing a coding
 * context after consequence annotation has already established the peptide
 * position and both residues. Stop-loss/extension, frameshift, and any
 * multi-residue operation remain on duckvep_hgvs_protein_fact_build(). */
DUCKVEP_INTERNAL_API duckvep_hgvs_status_t
duckvep_hgvs_protein_fact_build_single_residue(
    uint32_t                         position1,
    uint8_t                          reference,
    uint8_t                          alternate,
    uint32_t                         consequence_flags,
    duckvep_compat_profile_t         compatibility_profile,
    duckvep_hgvs_protein_fact_t     *out);

/* Return one one-letter peptide residue from the fact's clipped REF/ALT
 * sequence after any protein-level 3-prime rotation. */
DUCKVEP_INTERNAL_API duckvep_hgvs_status_t duckvep_hgvs_protein_base(
    const duckvep_hgvs_protein_fact_t *fact,
    int                                alternate,
    size_t                             index,
    uint8_t                           *residue_out);

/* Render the allocation-free VEP-style protein suffix (for example
 * p.(Arg97ProfsTer23)).  Ensembl's default three-letter amino-acid form is
 * intentional; protein identifiers remain a relational projection outside
 * the kernel. Facts built by the general interpreter borrow their coding
 * context; proven single-residue facts are self-contained. A missing or
 * undersized output buffer returns DUCKVEP_HGVS_BUFFER_TOO_SMALL and reports
 * the required byte count excluding NUL. */
DUCKVEP_INTERNAL_API duckvep_hgvs_status_t duckvep_hgvs_protein_render(
    const duckvep_hgvs_protein_fact_t *fact,
    int                                predicted,
    char                              *buffer,
    size_t                             capacity,
    size_t                            *required_out);

#ifdef __cplusplus
}
#endif

#endif /* DUCKVEP_HGVS_H */
