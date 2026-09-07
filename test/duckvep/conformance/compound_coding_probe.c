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
    duckvep_haplotype_block_t *blocks = NULL;
    uint8_t *ref_protein = NULL;
    duckvep_coding_context_t context;
    duckvep_sequence_delta_t delta;
    statuses[0] = DUCKVEP_CODING_CONTEXT_INVALID_ARG;
    statuses[1] = DUCKVEP_CONTEXT_DELTA_INVALID_ARG;
    memset(facts, 0, 22u * sizeof(*facts));
    memset(lengths, 0, 2u * sizeof(*lengths));
    if (*count < 0 || *count > 1000000 || *capacity < 1 ||
        (*strand != -1 && *strand != 1)) return;
    edits = calloc(*count ? (size_t)*count : 1u, sizeof(*edits));
    blocks = calloc(*count ? (size_t)*count : 1u, sizeof(*blocks));
    ref_protein = malloc((size_t)*capacity);
    if (!edits || !blocks || !ref_protein) goto cleanup;
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
    /* The actual complete path owns both peptide axes. Opening its interaction
     * windows must not reapply edits or substitute the path's total length change. */
    for (int i = 0; i < *count / 2; i++) {
        duckvep_haplotype_edit_t edit = edits[i];
        edits[i] = edits[*count - 1 - i]; edits[*count - 1 - i] = edit;
    }
    size_t block_count = 0u;
    if (duckvep_haplotype_partition(edits, (size_t)*count, blocks, (size_t)*count,
            &block_count) != DUCKVEP_HAPLOTYPE_OK) {
        facts[21] = -1;
        goto cleanup;
    }
    facts[17] = (int)block_count;
    for (size_t b = 0u; b < block_count; b++) {
        duckvep_coding_peptide_window_t view;
        if (!duckvep_coding_context_block_window_open(&context, blocks + b, &view)) {
            if (!facts[21]) facts[21] = (int)b + 1;
            continue;
        }
        facts[18]++;
        facts[19] += view.ref_peptide_offset != view.alt_peptide_offset;
        for (int side = 0; side < 2; side++) {
            size_t offset = side ? view.alt_peptide_offset : view.ref_peptide_offset;
            size_t whole = side ? view.alt_whole_length : view.ref_whole_length;
            size_t length = side ? view.alt_length : view.ref_length;
            for (size_t i = 0u; i < length; i++) {
                uint8_t expected = i == whole ? 'X' : (side ? protein : ref_protein)[offset + i];
                facts[20] += expected != duckvep_coding_context_peptide_window_base(
                    &context, &view, side, i);
            }
        }
    }
cleanup:
    free(ref_protein);
    free(blocks);
    free(edits);
}
