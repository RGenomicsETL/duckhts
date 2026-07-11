/*
 * duckvep_coding.h — sequence-backed coding SNV codon edit application
 * (INTERNAL).
 *
 * This is the narrow bridge from coordinate projection to codon consequence:
 * given a transcript-oriented translateable/CDS sequence and a projected coding
 * base, extract the reference codon, orient a genomic SNV allele onto the
 * transcript strand, apply it to the codon, and delegate translation/classifying
 * to duckvep_codon_change().
 *
 * Scope is deliberately small: SNVs/MNV-at-one-base only. In-frame indels,
 * frameshifts, start/stop-retained subtleties, and haplotype-combined edits land
 * in later kernels. The caller provides `cds_seq` such that cds_seq[0]
 * corresponds to CDS position 1; if the transcript has a positive Ensembl phase,
 * the adapter may include the phase-padding bases at the front (often `N`).
 */
#ifndef DUCKVEP_CODING_H
#define DUCKVEP_CODING_H

#include "duckvep_codon.h"
#include "duckvep_projection.h"

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum duckvep_coding_snv_status {
    DUCKVEP_CODING_SNV_OK = 0,
    DUCKVEP_CODING_SNV_INVALID_ARG,
    DUCKVEP_CODING_SNV_CODON_OUT_OF_RANGE,
    DUCKVEP_CODING_SNV_INVALID_BASE,
    DUCKVEP_CODING_SNV_REF_MISMATCH
} duckvep_coding_snv_status_t;

typedef struct duckvep_coding_snv_result {
    char     ref_codon[4];
    char     alt_codon[4];
    char     aa_ref;
    char     aa_alt;
    uint32_t change;          /* DUCKVEP_CODON_* */
    uint32_t cds_pos;
    uint32_t codon_start_cds;
    uint32_t protein_pos;
    uint8_t  codon_offset;
    char     ref_base_tx;     /* reference base after transcript-strand orientation */
    char     alt_base_tx;     /* alternate base after transcript-strand orientation */
} duckvep_coding_snv_result_t;

/* Apply a single genomic SNV allele to the projected transcript codon.
 *
 * `transcript_strand` is +1/-1; negative-strand genomic REF/ALT are reverse-
 * complemented before codon edit. Returns a status rather than silently producing
 * a consequence on ref mismatch or codon boundary overflow. On non-OK, `out` is
 * still zeroed and `change` is DUCKVEP_CODON_INVALID when `out != NULL`.
 */
duckvep_coding_snv_status_t duckvep_coding_snv_from_cds(
    const uint8_t                      *cds_seq,
    size_t                              cds_len,
    const duckvep_coding_projection_t  *projection,
    char                                genomic_ref,
    char                                genomic_alt,
    int8_t                              transcript_strand,
    duckvep_codon_table_t               table,
    duckvep_coding_snv_result_t        *out);

#ifdef __cplusplus
}
#endif

#endif /* DUCKVEP_CODING_H */
