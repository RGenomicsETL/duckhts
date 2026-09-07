/* Observe the existing coding-context evaluator on actual edit sets. This test
 * bridge owns scratch, not biological rules, and never flattens a compound edit. */
#include "duckvep_delta.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

void duckhts_test_compound_coding(
    char **reference, int *strand, int *starts, char **refs, char **alts,
    int *count, int *capacity, uint8_t *cds, uint8_t *protein,
    int *statuses, int *facts, int *lengths) {
    duckvep_haplotype_edit_t *edits = NULL;
    uint8_t *ref_protein = NULL;
    duckvep_coding_context_t context;
    duckvep_sequence_delta_t delta;
    statuses[0] = DUCKVEP_CODING_CONTEXT_INVALID_ARG;
    statuses[1] = DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    memset(facts, 0, 17u * sizeof(*facts));
    memset(lengths, 0, 2u * sizeof(*lengths));
    if (*count < 0 || *count > 1000000 || *capacity < 1 ||
        (*strand != -1 && *strand != 1)) return;
    edits = calloc(*count ? (size_t)*count : 1u, sizeof(*edits));
    ref_protein = malloc((size_t)*capacity);
    if (!edits || !ref_protein) goto cleanup;
    for (int i = 0; i < *count; i++) {
        size_t r = strlen(refs[i]), a = strlen(alts[i]);
        if (starts[i] < 1 || r > UINT32_MAX || a > UINT32_MAX) goto cleanup;
        edits[i] = (duckvep_haplotype_edit_t){(uint32_t)starts[i], (uint32_t)r,
            (const uint8_t *)refs[i], (uint32_t)a, (const uint8_t *)alts[i], 1};
    }
    duckvep_edit_set_t set = {edits, (size_t)*count};
    statuses[0] = duckvep_coding_context_build((const uint8_t *)*reference,
        strlen(*reference), &set, (int8_t)*strand, DUCKVEP_CODON_TABLE_STANDARD,
        cds, (size_t)*capacity, ref_protein, (size_t)*capacity,
        protein, (size_t)*capacity, &context);
    if (statuses[0] != DUCKVEP_CODING_CONTEXT_OK) goto cleanup;
    if (context.alt_cds_len > INT_MAX || context.alt_peptide_len > INT_MAX) {
        statuses[0] = DUCKVEP_CODING_CONTEXT_OUT_OF_RANGE;
        goto cleanup;
    }
    lengths[0] = (int)context.alt_cds_len;
    lengths[1] = (int)context.alt_peptide_len;
    statuses[1] = duckvep_coding_context_delta_fill(&context, 0u, &delta);
    int observed[] = {delta.valid, delta.sequence_status, delta.synonymous, delta.missense,
        delta.stop_gained, delta.stop_lost, delta.stop_retained, delta.start_lost,
        delta.start_retained, delta.frameshift, delta.inframe_deletion,
        delta.inframe_insertion, delta.protein_altering, delta.coding_unknown,
        delta.partial_codon, (int)context.flags, (int)context.applied_edits};
    memcpy(facts, observed, sizeof observed);
cleanup:
    free(ref_protein);
    free(edits);
}
