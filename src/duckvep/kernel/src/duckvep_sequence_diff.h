/* Aligned differences over borrowed sequence bytes (INTERNAL).
 * Ensembl Variation release/116, 2fb834b987ede3824e200197a838ce11e91aeb4b:
 * Utils::Sequence::align_seqs and TranscriptHaplotype::_get_raw_diffs.
 * The supported alignment is the pure-Perl NW path, not Bio::Ext::Align.
 */
#ifndef DUCKVEP_SEQUENCE_DIFF_H
#define DUCKVEP_SEQUENCE_DIFF_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    size_t ref_start0, alt_start0, ref_length, alt_length;
    size_t alignment_start0;
} duckvep_sequence_difference_t;

typedef struct {
    uint64_t *scores; /* Two rows; at least 2*(alt_length+1) elements for NW. */
    size_t score_capacity;
    uint8_t *trace;
    size_t trace_capacity; /* Cells, one byte per traceback direction. */
} duckvep_sequence_diff_scratch_t;

typedef enum {
    DUCKVEP_SEQUENCE_DIFF_OK,
    DUCKVEP_SEQUENCE_DIFF_INVALID_ARG,
    DUCKVEP_SEQUENCE_DIFF_SCORE_FULL,
    DUCKVEP_SEQUENCE_DIFF_TRACE_FULL,
    DUCKVEP_SEQUENCE_DIFF_OUTPUT_FULL
} duckvep_sequence_diff_status_t;

typedef struct {
    size_t count, alignment_length, trace_cells;
} duckvep_sequence_diff_result_t;

/* Positional comparison when align_indels == 0 (pad the shorter sequence with
 * gaps); otherwise global NW, match +1 / mismatch -1 / gap -1, with traceback
 * ties preferring reference deletion, then alternate insertion, then diagonal.
 * Adjacent differences join only when both sides have the same gap/non-gap type.
 * Results borrow ungapped slices of the input; empty spans denote insertions or
 * deletions. Coordinates are zero-based on the two sequences and the alignment.
 *
 * No allocation. All mutable spans and result storage must be distinct from one
 * another and from the inputs. Inputs may alias each other. Bytes must be nonzero
 * ASCII without '-' (the alignment gap). count reports the required difference
 * capacity on OUTPUT_FULL; trace_cells reports required cells on TRACE_FULL.
 * No output difference is written until its complete capacity has been checked.
 */
duckvep_sequence_diff_status_t duckvep_sequence_differences(
    const uint8_t *reference, size_t ref_length, const uint8_t *alternate, size_t alt_length,
    int align_indels, const duckvep_sequence_diff_scratch_t *scratch,
    duckvep_sequence_difference_t *differences, size_t difference_capacity,
    duckvep_sequence_diff_result_t *result);

#endif
