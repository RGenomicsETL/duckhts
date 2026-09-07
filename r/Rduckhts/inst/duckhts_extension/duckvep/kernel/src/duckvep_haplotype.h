/*
 * duckvep_haplotype.h — multi-edit CDS haplotype mutation helpers (INTERNAL).
 *
 * Callers group and project phased variants, then pass a
 * transcript-oriented CDS plus per-haplotype edits. The kernel rebuilds a
 * caller-owned scratch sequence in one reverse CDS-coordinate pass, translates
 * once, and emits aggregate indel/frameshift flags. No allocation, no DuckDB/htslib.
 *
 * This header owns sequence mutation and translation only. Stream grouping and
 * DuckDB materialization stay outside the kernel.
 */
#ifndef DUCKVEP_HAPLOTYPE_H
#define DUCKVEP_HAPLOTYPE_H

#include "duckvep_codon.h"

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum duckvep_haplotype_status {
    DUCKVEP_HAPLOTYPE_OK = 0,
    DUCKVEP_HAPLOTYPE_INVALID_ARG,
    DUCKVEP_HAPLOTYPE_OUT_OF_RANGE,
    DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL,
    DUCKVEP_HAPLOTYPE_INVALID_BASE,
    DUCKVEP_HAPLOTYPE_REF_MISMATCH,
    DUCKVEP_HAPLOTYPE_EDIT_ORDER,
    DUCKVEP_HAPLOTYPE_INPUT_INCOMPLETE
} duckvep_haplotype_status_t;

enum {
    DUCKVEP_HAPLOTYPE_FLAG_INDEL               = 1u << 0,
    DUCKVEP_HAPLOTYPE_FLAG_FRAMESHIFT          = 1u << 1,
    DUCKVEP_HAPLOTYPE_FLAG_RESOLVED_FRAMESHIFT = 1u << 2,
    DUCKVEP_HAPLOTYPE_FLAG_STOP_TRUNCATED      = 1u << 3
};

typedef struct duckvep_haplotype_edit {
    uint32_t       cds_start;      /* 1-based position in the original CDS. */
    uint32_t       ref_len;        /* Number of CDS bases replaced; 0 = insertion before cds_start. */
    const uint8_t *ref;            /* Allele/reference bases in variant_strand orientation when ref_len > 0. */
    uint32_t       alt_len;        /* Alternate allele length after removing VCF deletion '-' characters. */
    const uint8_t *alt;            /* Allele bases in variant_strand orientation when alt_len > 0. */
    int8_t         variant_strand; /* +1/-1; allele is reverse-complemented if != transcript_strand. */
} duckvep_haplotype_edit_t;

typedef struct duckvep_haplotype_result {
    size_t  cds_len;
    int64_t length_diff;
    uint32_t flags;
    size_t  applied_edits;
} duckvep_haplotype_result_t;

/* One interaction block on the reference and replayed CDS axes. Edits remain in the same
 * block while their cumulative length change displaces the reading frame, or
 * while the next edit still touches the same alternate-sequence codon.
 * The spans include retained bases between edits, not an alignment or HGVS
 * normalization. A pure insertion has ref_len=0; a pure deletion has alt_len=0. */
typedef struct duckvep_haplotype_block {
    size_t   edit_begin;
    size_t   edit_count;
    uint32_t cds_start; /* 1-based reference CDS; insertion is before this base. */
    uint32_t ref_len;
    size_t   alt_start0; /* 0-based offset in the complete replayed CDS. */
    size_t   alt_len;
    int64_t  length_diff;
    uint32_t flags;
} duckvep_haplotype_block_t;

/* Construct the pinned VEP-116 Haplosaurus reference protein from a borrowed
 * transcript-oriented CDS and optional sorted single-residue Ensembl SeqEdits.
 * Uses consensus translation, removes the last complete translated stop,
 * forces a legitimate initial start to M, then applies peptide edits. Finally
 * the container's exact uppercase raw-CDS suffix TAA/TAG/TGA appends '*', even
 * for nonstandard tables or partial CDS. Internal stops are not truncated.
 * Peptide positions are one-based, strictly increasing and <= cds_length/3;
 * alternates are A-Z or '*'. A terminal SeqEdit may append a removed residue.
 * Capacity must be at least cds_length/3 + 2 (extra stop plus NUL). Shorter CDS
 * returns INPUT_INCOMPLETE, not an invented empty reference protein.
 * Input spans and length output must not overlap peptide[0..capacity).
 * No allocation. Length is zero on failure; invalid bases may leave partial bytes. */
duckvep_haplotype_status_t duckvep_haplotype_reference_protein(
    const uint8_t *cds, size_t cds_length, duckvep_codon_table_t table,
    const uint32_t *edit_positions1, const uint8_t *edit_alternates, size_t edit_count,
    uint8_t *peptide, size_t capacity, size_t *length);

/* Partition edits sorted by ascending original CDS coordinate. The caller has
 * already grouped them by model, transcript, sample, phase set, and haplotype.
 * Call with blocks=NULL/capacity=0 to obtain required_blocks. No partial block
 * array is published when capacity is insufficient. */
duckvep_haplotype_status_t duckvep_haplotype_partition(
    const duckvep_haplotype_edit_t *edits,
    size_t                          edit_count,
    duckvep_haplotype_block_t      *blocks,
    size_t                          block_capacity,
    size_t                         *required_blocks);

/* Test whether a half-open span on the rebuilt CDS intersects this block's
 * frame displacement. The edit array is ascending and the block is one of its
 * actual partitions. A displacement starts at the first frame-changing edit
 * and ends after the ALT bases of the restoring edit; a final open frame
 * continues downstream. Inserted/replacement bases belong to the displaced
 * span when either the entering or leaving frame is displaced. Earlier closed
 * blocks may shift ALT by whole codons. This is geometry, not a consequence or
 * an assertion that a restored frame rescues a protein. In particular, a caller
 * can query the three bases of an already translated stop. No allocation;
 * failure leaves intersects zero, and an empty query never intersects. */
duckvep_haplotype_status_t duckvep_haplotype_block_frame_intersects(
    const duckvep_haplotype_edit_t  *edits,
    size_t                          edit_count,
    const duckvep_haplotype_block_t *block,
    size_t                          alt_start0,
    size_t                          alt_length,
    int                            *intersects);

/* Apply edits to `ref_cds`, writing the mutated CDS to `cds_out` and its length
 * to `cds_len_out`. Edits must be sorted by descending original CDS coordinate
 * and must not overlap in original CDS space, including two insertions at the
 * same interbase site. Such edits require prior conflict resolution with their
 * source provenance; input order must not choose the inserted sequence.
 * This mirrors Ensembl's reverse
 * mapping order and permits one linear rebuild without allocation or sorting.
 * The reference CDS, edit descriptors and allele byte spans are borrowed and
 * must not overlap the output span `cds_out[0..cds_cap)`. An overlap returns
 * INVALID_ARG before writing sequence bytes. `cds_cap` need only cover the final
 * CDS; no intermediate mutated sequence is constructed. Output length/result
 * storage must be separate from the input and sequence buffers.
 *
 * `ref`/`alt` alleles are oriented from variant_strand to transcript_strand
 * before validation/application. Bases must be A/C/G/T (case-insensitive; U is
 * accepted as T). The source CDS may contain N, but edited/reference-validated
 * positions must match A/C/G/T after orientation.
 */
duckvep_haplotype_status_t duckvep_haplotype_apply_cds_edits(
    const uint8_t                    *ref_cds,
    size_t                            ref_cds_len,
    const duckvep_haplotype_edit_t    *edits,
    size_t                            edit_count,
    int8_t                            transcript_strand,
    uint8_t                          *cds_out,
    size_t                            cds_cap,
    size_t                           *cds_len_out,
    duckvep_haplotype_result_t       *result);


#ifdef __cplusplus
}
#endif

#endif /* DUCKVEP_HAPLOTYPE_H */
