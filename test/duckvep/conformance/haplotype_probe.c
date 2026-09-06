/* Test-only R .C bridges to the host-neutral mutation and carrier-path kernels.
 * R owns every input/output buffer; the bridges own bounded test scratch.
 * These are not SQL adapters and do not decode GT or choose a phase policy. */
#include "duckvep_haplotype.h"
#include "duckvep_carriers.h"
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

static int descending_edit(const void *a, const void *b) {
    uint32_t x = ((const duckvep_haplotype_edit_t *)a)->cds_start;
    uint32_t y = ((const duckvep_haplotype_edit_t *)b)->cds_start;
    return x < y ? 1 : x > y ? -1 : 0;
}

/* Exercise the real sparse index with the same three diploid fixture samples.
 * R supplies the allele-slot assignments that generated the VCF, not grouped
 * haplotypes. The bridge owns bounded test scratch and adapts each distinct
 * completed event path to the existing mutation/translation kernel once.
 * No genomic projection or phase-policy inference is implemented here. */
void duckhts_test_carrier_haplotypes(
    uint8_t *reference, int *reference_length, int *strand, int *positions,
    int *starts, char **refs, char **alts, int *edit_count, int *order, int *lanes,
    int *capacity, uint8_t *cds, uint8_t *protein, int *cds_lengths,
    int *protein_lengths, int *flags, int *contributors, int *contributor_counts,
    double *metrics, int *status) {
    duckvep_carrier_buffers_t b = {0};
    duckvep_carriers_t stream;
    duckvep_haplotype_edit_t *all_edits = NULL, *leaf_edits = NULL;
    uint64_t *events = NULL;
    uint8_t *cds_scratch = NULL, *protein_scratch = NULL;
    uint16_t chrom = 0u;
    uint32_t end1 = UINT32_MAX, completed, active;
    duckvep_transcript_model_t model = {0};
    duckvep_carrier_transcript_t transcript;
    duckvep_carrier_bucket_t transcript_index[2];
    duckvep_carriers_status_t stream_status;
    *status = DUCKVEP_HAPLOTYPE_INVALID_ARG;
    memset(metrics, 0, 6u * sizeof(*metrics));
    if (*reference_length < 0 || *edit_count < 1 || *edit_count > 1000000 ||
        *capacity < 1 || *capacity > INT_MAX / 6 || (*strand != -1 && *strand != 1)) return;
    size_t count = (size_t)*edit_count;
    all_edits = calloc(count, sizeof(*all_edits));
    leaf_edits = calloc(count, sizeof(*leaf_edits));
    events = calloc(count, sizeof(*events));
    cds_scratch = malloc((size_t)*capacity);
    protein_scratch = malloc((size_t)*capacity);
    b.transcripts = &transcript;
    b.active_transcripts = &active;
    b.transcript_index = transcript_index;
    b.transcript_capacity = 1u; b.transcript_buckets = 2u;
    b.call_capacity = 6u; b.call_buckets = 16u;
    b.prefix_capacity = 3u * (uint32_t)count;
    b.prefix_buckets = 1u;
    while (b.prefix_buckets < 2u * b.prefix_capacity) b.prefix_buckets *= 2u;
    b.calls = calloc(b.call_capacity, sizeof(*b.calls));
    b.call_index = calloc(b.call_buckets, sizeof(*b.call_index));
    b.prefixes = calloc(b.prefix_capacity, sizeof(*b.prefixes));
    b.prefix_index = calloc(b.prefix_buckets, sizeof(*b.prefix_index));
    if (!all_edits || !leaf_edits || !events || !cds_scratch || !protein_scratch ||
        !b.calls || !b.call_index || !b.prefixes || !b.prefix_index) goto cleanup;
    for (size_t i = 0u; i < count; i++) {
        size_t ref_len = strlen(refs[i]), alt_len = strlen(alts[i]);
        if (starts[i] < 1 || positions[i] < 1 || ref_len > UINT32_MAX || alt_len > UINT32_MAX ||
            order[i] < 0 || (size_t)order[i] >= count) goto cleanup;
        all_edits[i] = (duckvep_haplotype_edit_t){(uint32_t)starts[i], (uint32_t)ref_len,
            (const uint8_t *)refs[i], (uint32_t)alt_len, (const uint8_t *)alts[i], 1};
    }
    model.transcript_count = 1u; model.chrom_id = &chrom; model.end1 = &end1;
    stream_status = duckvep_carriers_init(&stream, &model, &b);
    if (stream_status != DUCKVEP_CARRIERS_OK) goto stream_error;
    for (size_t i = 0u; i < count; i++) {
        size_t event = (size_t)order[i];
        for (uint32_t sample = 0u; sample < 3u; sample++) {
            /* Resume even between carriers of the same input event, simulating
             * a source batch ending at every possible carrier row. */
            stream_status = duckvep_carriers_advance(&stream, 0u,
                (uint32_t)positions[event], event + 1u, &completed);
            if (stream_status != DUCKVEP_CARRIERS_OK) goto stream_error;
            int lane = lanes[event * 3u + sample];
            if (lane < 1 || lane > 2) goto cleanup;
            duckvep_carrier_key_t key = {sample, 10, (uint16_t)lane, 2u, 1u};
            stream_status = duckvep_carriers_push(&stream, 0u, &key);
            if (stream_status != DUCKVEP_CARRIERS_OK) goto stream_error;
        }
    }
    /* Reference lanes are implicit in the sparse index. The fixture adapter
     * materializes them from one shared reference translation for comparison
     * with Haplosaurus, which reports all six lanes. */
    size_t cds_len, protein_len;
    duckvep_haplotype_result_t applied, translated;
    *status = duckvep_haplotype_apply_cds_edits(reference, (size_t)*reference_length,
        NULL, 0u, (int8_t)*strand, cds_scratch, (size_t)*capacity, &cds_len, &applied);
    if (*status != DUCKVEP_HAPLOTYPE_OK) goto cleanup;
    *status = duckvep_haplotype_translate_cds(cds_scratch, cds_len, DUCKVEP_CODON_TABLE_STANDARD,
        protein_scratch, (size_t)*capacity, &protein_len, &translated);
    if (*status != DUCKVEP_HAPLOTYPE_OK) goto cleanup;
    for (size_t slot = 0u; slot < 6u; slot++) {
        memcpy(cds + slot * (size_t)*capacity, cds_scratch, cds_len);
        memcpy(protein + slot * (size_t)*capacity, protein_scratch, protein_len);
        cds_lengths[slot] = (int)cds_len; protein_lengths[slot] = (int)protein_len;
        flags[slot] = (int)(applied.flags | translated.flags);
        contributor_counts[slot] = 0;
    }
    stream_status = duckvep_carriers_finish(&stream, &completed);
    if (stream_status != DUCKVEP_CARRIERS_TRANSCRIPT_READY) goto stream_error;
    duckvep_carrier_leaf_t leaf;
    while ((stream_status = duckvep_carriers_next_leaf(&stream, &leaf)) == DUCKVEP_CARRIERS_OK) {
        size_t required;
        stream_status = duckvep_carriers_leaf_events(&stream, leaf.id, events, count, &required);
        if (stream_status != DUCKVEP_CARRIERS_OK) goto stream_error;
        for (size_t i = 0u; i < required; i++) {
            if (!events[i] || events[i] > count) goto invalid;
            leaf_edits[i] = all_edits[events[i] - 1u];
        }
        qsort(leaf_edits, required, sizeof(*leaf_edits), descending_edit);
        *status = duckvep_haplotype_apply_cds_edits(reference, (size_t)*reference_length,
            leaf_edits, required, (int8_t)*strand, cds_scratch, (size_t)*capacity, &cds_len, &applied);
        if (*status != DUCKVEP_HAPLOTYPE_OK) goto cleanup;
        *status = duckvep_haplotype_translate_cds(cds_scratch, cds_len, DUCKVEP_CODON_TABLE_STANDARD,
            protein_scratch, (size_t)*capacity, &protein_len, &translated);
        if (*status != DUCKVEP_HAPLOTYPE_OK) goto cleanup;
        metrics[3] += 1.0;
        metrics[4] += (double)cds_len;
        for (uint32_t id = leaf.first_call; id;) {
            const duckvep_carrier_call_t *call = duckvep_carriers_call(&stream, id);
            if (!call || call->key.sample_index >= 3u || call->key.lane > 2u) goto invalid;
            size_t slot = (size_t)call->key.sample_index * 2u + call->key.lane - 1u;
            memcpy(cds + slot * (size_t)*capacity, cds_scratch, cds_len);
            memcpy(protein + slot * (size_t)*capacity, protein_scratch, protein_len);
            cds_lengths[slot] = (int)cds_len; protein_lengths[slot] = (int)protein_len;
            flags[slot] = (int)(applied.flags | translated.flags);
            contributor_counts[slot] = (int)required;
            for (size_t i = 0u; i < required; i++) contributors[slot * count + i] = (int)events[i] - 1;
            metrics[5] += 1.0;
            id = call->next_leaf;
        }
    }
    if (stream_status != DUCKVEP_CARRIERS_DONE) goto stream_error;
    stream_status = duckvep_carriers_release(&stream);
    if (stream_status != DUCKVEP_CARRIERS_OK) goto stream_error;
    stream_status = duckvep_carriers_finish(&stream, &completed);
    if (stream_status != DUCKVEP_CARRIERS_DONE) goto stream_error;
    metrics[0] = stream.peak_transcripts;
    metrics[1] = stream.peak_calls;
    metrics[2] = stream.peak_prefixes;
    *status = DUCKVEP_HAPLOTYPE_OK;
    goto cleanup;
invalid:
    *status = DUCKVEP_HAPLOTYPE_INVALID_ARG;
    goto cleanup;
stream_error:
    *status = 100 + (int)stream_status;
cleanup:
    free(all_edits); free(leaf_edits); free(events);
    free(cds_scratch); free(protein_scratch);
    free(b.calls); free(b.call_index); free(b.prefixes); free(b.prefix_index);
}
