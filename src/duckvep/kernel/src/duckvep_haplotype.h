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
    DUCKVEP_HAPLOTYPE_EDIT_ORDER
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
    size_t  protein_len;
    int64_t length_diff;
    uint32_t flags;
    size_t  applied_edits;
} duckvep_haplotype_result_t;

/* One interaction block on the original CDS axis. Edits remain in the same
 * block while their cumulative length change displaces the reading frame, or
 * while the next edit still touches the same alternate-sequence codon. */
typedef struct duckvep_haplotype_block {
    size_t   edit_begin;
    size_t   edit_count;
    uint32_t cds_start;
    uint32_t cds_end;
    int64_t  length_diff;
    uint32_t flags;
} duckvep_haplotype_block_t;

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

/* Apply edits to `ref_cds`, writing the mutated CDS to `cds_out` and its length
 * to `cds_len_out`. Edits must be sorted by descending original CDS coordinate
 * and must not overlap in original CDS space; this mirrors Ensembl's reverse
 * mapping order and permits one linear rebuild without allocation or sorting.
 * `ref_cds == cds_out` remains supported for compatibility, but distinct input
 * and output buffers are the linear streaming path. Exact aliasing requires
 * `cds_cap` to cover the largest intermediate sequence produced while applying
 * the descending edits, not only the final sequence. Partially overlapping
 * buffers are not supported.
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

/* Translate a mutated CDS once and truncate after the first stop, matching the
 * Haplosaurus container's post-translation behavior that keeps the first '*'
 * and drops later amino acids. `protein_out` is NUL-terminated on success.
 * `result` is overwritten with translation-side lengths/flags; callers that
 * need apply+translate aggregate flags should OR the two result flag fields. */
duckvep_haplotype_status_t duckvep_haplotype_translate_cds(
    const uint8_t                    *cds,
    size_t                            cds_len,
    duckvep_codon_table_t             table,
    uint8_t                          *protein_out,
    size_t                            protein_cap,
    size_t                           *protein_len_out,
    duckvep_haplotype_result_t       *result);

#ifdef __cplusplus
}
#endif

#endif /* DUCKVEP_HAPLOTYPE_H */
