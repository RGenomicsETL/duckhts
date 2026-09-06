/* Test-only R .C bridge to the unmodified, host-neutral mutation kernel.
 * R owns every input/output buffer; this bridge owns only edit descriptors.
 * This is not a SQL adapter or an implementation of carrier stream grouping. */
#include "duckvep_haplotype.h"
#include <limits.h>
#include <stdlib.h>
#include <string.h>

void duckhts_test_haplotype(
    uint8_t *reference, int *reference_length, int *strand,
    int *starts, char **refs, char **alts, int *edit_count,
    uint8_t *cds, int *cds_capacity, uint8_t *protein, int *protein_capacity,
    int *status, int *cds_length, int *protein_length, int *flags) {
    duckvep_haplotype_edit_t *edits = NULL;
    duckvep_haplotype_result_t applied, translated;
    size_t cds_len = 0, protein_len = 0;
    *status = DUCKVEP_HAPLOTYPE_INVALID_ARG;
    *cds_length = *protein_length = *flags = 0;
    if (*reference_length < 0 || *edit_count < 0 || *edit_count > 1000000 ||
        *cds_capacity < 1 || *protein_capacity < 1 ||
        (*strand != -1 && *strand != 1)) return;
    if (*edit_count) {
        edits = calloc((size_t)*edit_count, sizeof(*edits));
        if (!edits) return;
    }
    for (int i = 0; i < *edit_count; i++) {
        size_t ref_len = strlen(refs[i]), alt_len = strlen(alts[i]);
        if (starts[i] < 1 || ref_len > UINT32_MAX || alt_len > UINT32_MAX) goto cleanup;
        edits[i].cds_start = (uint32_t)starts[i];
        edits[i].ref_len = (uint32_t)ref_len;
        edits[i].alt_len = (uint32_t)alt_len;
        edits[i].ref = (const uint8_t *)refs[i];
        edits[i].alt = (const uint8_t *)alts[i];
        edits[i].variant_strand = 1;
    }
    *status = duckvep_haplotype_apply_cds_edits(reference, (size_t)*reference_length,
        edits, (size_t)*edit_count, (int8_t)*strand, cds, (size_t)*cds_capacity,
        &cds_len, &applied);
    if (*status != DUCKVEP_HAPLOTYPE_OK) goto cleanup;
    *status = duckvep_haplotype_translate_cds(cds, cds_len,
        DUCKVEP_CODON_TABLE_STANDARD, protein, (size_t)*protein_capacity,
        &protein_len, &translated);
    if (*status != DUCKVEP_HAPLOTYPE_OK) goto cleanup;
    if (cds_len > INT_MAX || protein_len > INT_MAX) {
        *status = DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
        goto cleanup;
    }
    *cds_length = (int)cds_len;
    *protein_length = (int)protein_len;
    *flags = (int)(applied.flags | translated.flags);
cleanup:
    free(edits);
}
