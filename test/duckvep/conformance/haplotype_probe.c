/* Test-only R .C bridges to the host-neutral mutation and carrier-path kernels.
 * R owns every input/output buffer; the bridges own bounded test scratch.
 * These are not SQL adapters and do not decode GT or choose a phase policy. */
#include "duckvep_haplotype.h"
#include "duckvep_carriers.h"
#include "duckvep_delta.h"
#include "duckvep_haplotype_stream.h"
#include <limits.h>
#include <stdlib.h>
#include <string.h>

/* Keep the R bridge's existing failure categories during the native API change. */
static duckvep_haplotype_status_t translation_status(duckvep_translation_status_t status) {
    switch (status) {
    case DUCKVEP_TRANSLATION_OK: return DUCKVEP_HAPLOTYPE_OK;
    case DUCKVEP_TRANSLATION_BUFFER_TOO_SMALL: return DUCKVEP_HAPLOTYPE_BUFFER_TOO_SMALL;
    case DUCKVEP_TRANSLATION_INVALID_BASE: return DUCKVEP_HAPLOTYPE_INVALID_BASE;
    default: return DUCKVEP_HAPLOTYPE_INVALID_ARG;
    }
}

void duckhts_test_haplotype(
    uint8_t *reference, int *reference_length, int *strand,
    int *starts, char **refs, char **alts, int *edit_count,
    uint8_t *cds, int *cds_capacity, uint8_t *protein, int *protein_capacity,
    int *status, int *cds_length, int *protein_length, int *flags) {
    duckvep_haplotype_edit_t *edits = NULL;
    duckvep_haplotype_result_t applied;
    duckvep_translation_t translated;
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
    duckvep_translation_status_t tst = duckvep_translate_cds(cds, cds_len,
        DUCKVEP_CODON_TABLE_STANDARD, protein, (size_t)*protein_capacity, &translated);
    *status = translation_status(tst);
    if (*status != DUCKVEP_HAPLOTYPE_OK) goto cleanup;
    protein_len = translated.first_stop_position1 ? translated.first_stop_position1 : translated.length;
    protein[protein_len] = 0u;
    if (cds_len > INT_MAX || protein_len > INT_MAX) {
        *status = DUCKVEP_HAPLOTYPE_OUT_OF_RANGE;
        goto cleanup;
    }
    *cds_length = (int)cds_len;
    *protein_length = (int)protein_len;
    *flags = (int)applied.flags;
    if (protein_len < translated.length) *flags |= DUCKVEP_HAPLOTYPE_FLAG_STOP_TRUNCATED;
cleanup:
    free(edits);
}

/* Exercise the native replay stream with the same three diploid fixture samples.
 * R supplies genomic VCF alleles, ranked all-coding exons, and the allele-slot
 * assignments that generated the VCF, not projected CDS edit coordinates or
 * grouped haplotypes. Existing event/CDS projection authorities run once per
 * event, then each distinct completed path is rebuilt and translated once.
 * The bridge owns bounded test scratch; it does not infer a phase policy. */
void duckhts_test_carrier_haplotypes(
    uint8_t *reference, int *reference_length, int *strand,
    int *exon_starts, int *exon_ends, int *exon_count, int *positions,
    char **refs, char **alts, int *edit_count, int *order, int *lanes,
    int *capacity, uint8_t *cds, uint8_t *protein, int *cds_lengths,
    int *protein_lengths, int *flags, int *contributors, int *contributor_counts,
    double *metrics, int *status) {
    duckvep_carrier_buffers_t b = {0};
    duckvep_haplotype_stream_t stream;
    duckvep_haplotype_stream_buffers_t storage = {0};
    uint32_t *exon_storage = NULL;
    uint8_t *cds_scratch = NULL, *protein_scratch = NULL;
    uint16_t chrom = 0u;
    uint32_t start1 = UINT32_MAX, end1 = 0u, active, zero = 0u;
    uint16_t model_exon_count;
    uint64_t cds_offset = 0u;
    uint32_t cds_length;
    int8_t model_strand;
    duckvep_transcript_model_t model = {0};
    duckvep_exon_model_t exons = {0};
    duckvep_sequence_pool_t sequences = {0};
    duckvep_carrier_transcript_t transcript;
    duckvep_carrier_bucket_t transcript_index[2];
    duckvep_haplotype_stream_status_t stream_status;
    *status = DUCKVEP_HAPLOTYPE_INVALID_ARG;
    memset(metrics, 0, 6u * sizeof(*metrics));
    if (*reference_length < 1 || *edit_count < 1 || *edit_count > 1000000 ||
        *exon_count < 1 || *exon_count > UINT16_MAX ||
        *capacity < 1 || *capacity > INT_MAX / 6 || (*strand != -1 && *strand != 1)) return;
    size_t count = (size_t)*edit_count;
    storage.event_capacity = storage.projection_capacity = (uint32_t)count;
    storage.leaf_capacity = count;
    for (size_t i = 0u; i < count; i++) {
        size_t r = strlen(refs[i]), a = strlen(alts[i]);
        if (r > UINT16_MAX || a > UINT16_MAX ||
            r + a > SIZE_MAX - storage.allele_capacity) goto cleanup;
        storage.allele_capacity += r + a;
        size_t edits = r == a && r > 1u ? (r + 1u) / 2u : 1u;
        if (edits > SIZE_MAX - storage.edit_capacity) goto cleanup;
        storage.edit_capacity += edits;
    }
    storage.events = calloc(count, sizeof(*storage.events));
    storage.projections = calloc(count, sizeof(*storage.projections));
    storage.leaf_events = calloc(count, sizeof(*storage.leaf_events));
    storage.contributors = calloc(count, sizeof(*storage.contributors));
    storage.edits = calloc(storage.edit_capacity, sizeof(*storage.edits));
    storage.blocks = calloc(storage.edit_capacity, sizeof(*storage.blocks));
    storage.alleles = malloc(storage.allele_capacity);
    exon_storage = calloc(4u * (size_t)*exon_count, sizeof(*exon_storage));
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
    if (!storage.events || !storage.projections || !storage.leaf_events ||
        !storage.contributors || !storage.edits || !storage.blocks || !storage.alleles ||
        !exon_storage || !cds_scratch || !protein_scratch ||
        !b.calls || !b.call_index || !b.prefixes || !b.prefix_index) goto cleanup;
    model_exon_count = (uint16_t)*exon_count;
    model_strand = (int8_t)*strand;
    cds_length = (uint32_t)*reference_length;
    uint32_t cdna = 0u;
    for (size_t i = 0u; i < model_exon_count; i++) {
        if (exon_starts[i] < 1 || exon_ends[i] < exon_starts[i] ||
            (i && (*strand > 0 ? exon_starts[i] <= exon_ends[i - 1u]
                               : exon_ends[i] >= exon_starts[i - 1u]))) goto cleanup;
        uint32_t es = (uint32_t)exon_starts[i], ee = (uint32_t)exon_ends[i];
        uint32_t length = ee - es + 1u;
        if (length > UINT32_MAX - cdna) goto cleanup;
        exon_storage[i] = es;
        exon_storage[model_exon_count + i] = ee;
        exon_storage[2u * model_exon_count + i] = cdna + 1u;
        cdna += length;
        exon_storage[3u * model_exon_count + i] = cdna;
        if (es < start1) start1 = es;
        if (ee > end1) end1 = ee;
    }
    if (cdna != cds_length) goto cleanup;
    model.transcript_count = 1u; model.chrom_id = &chrom;
    model.start1 = model.cds_start1 = &start1;
    model.end1 = model.cds_end1 = &end1;
    model.strand = &model_strand;
    model.exon_offset = &zero; model.exon_count = &model_exon_count;
    exons.exon_count = model_exon_count; exons.start1 = exon_storage;
    exons.end1 = exon_storage + model_exon_count;
    exons.cdna_start1 = exon_storage + 2u * model_exon_count;
    exons.cdna_end1 = exon_storage + 3u * model_exon_count;
    sequences.cds_bytes = reference; sequences.cds_bytes_len = cds_length;
    sequences.cds_offset = &cds_offset; sequences.cds_length = &cds_length;
    sequences.transcript_count = 1u;
    storage.carriers = b;
    storage.cds = cds_scratch; storage.cds_capacity = (size_t)*capacity;
    storage.protein = protein_scratch; storage.protein_capacity = (size_t)*capacity;
    stream_status = duckvep_haplotype_stream_init(&stream, &model, &exons, &sequences, &storage);
    if (stream_status != DUCKVEP_HAPLOTYPE_STREAM_OK) goto stream_error;
    for (size_t i = 0u; i < count; i++) {
        if (order[i] < 0 || (size_t)order[i] >= count) goto cleanup;
        size_t event = (size_t)order[i];
        if (positions[event] < 1) goto cleanup;
        duckvep_haplotype_source_t source = {event + 1u,
            (const uint8_t *)refs[event], (const uint8_t *)alts[event],
            (uint32_t)positions[event], 0u,
            (uint16_t)strlen(refs[event]), (uint16_t)strlen(alts[event])};
        stream_status = duckvep_haplotype_stream_begin(&stream, &source);
        if (stream_status != DUCKVEP_HAPLOTYPE_STREAM_OK) goto stream_error;
        stream_status = duckvep_haplotype_stream_project(&stream, zero);
        if (stream_status != DUCKVEP_HAPLOTYPE_STREAM_OK) goto stream_error;
        for (uint32_t sample = 0u; sample < 3u; sample++) {
            /* Each push may be separated by an input batch or vector edge. */
            int lane = lanes[event * 3u + sample];
            if (lane < 1 || lane > 2) goto cleanup;
            int32_t gt[2] = {lane == 1, lane == 2};
            const uint8_t phase[2] = {1u, 1u};
            duckvep_haplotype_phase_set_t set = {10, 1u};
            duckvep_haplotype_call_t call = {gt, phase, sample, 1u, 2u, set, DUCKVEP_PHASE_STRICT};
            stream_status = duckvep_haplotype_stream_push_call(&stream, zero, &call, &set, 1u);
            if (stream_status != DUCKVEP_HAPLOTYPE_STREAM_OK) goto stream_error;
        }
    }
    /* Reference lanes are implicit in the sparse index. The fixture adapter
     * materializes them from one shared reference translation for comparison
     * with Haplosaurus, which reports all six lanes. */
    size_t cds_len, protein_len;
    duckvep_haplotype_result_t applied;
    duckvep_translation_t translated;
    *status = duckvep_haplotype_apply_cds_edits(reference, (size_t)*reference_length,
        NULL, 0u, (int8_t)*strand, cds_scratch, (size_t)*capacity, &cds_len, &applied);
    if (*status != DUCKVEP_HAPLOTYPE_OK) goto cleanup;
    duckvep_translation_status_t tst = duckvep_translate_cds(cds_scratch, cds_len,
        DUCKVEP_CODON_TABLE_STANDARD, protein_scratch, (size_t)*capacity, &translated);
    *status = translation_status(tst);
    if (*status != DUCKVEP_HAPLOTYPE_OK) goto cleanup;
    protein_len = translated.first_stop_position1 ? translated.first_stop_position1 : translated.length;
    for (size_t slot = 0u; slot < 6u; slot++) {
        memcpy(cds + slot * (size_t)*capacity, cds_scratch, cds_len);
        memcpy(protein + slot * (size_t)*capacity, protein_scratch, protein_len);
        cds_lengths[slot] = (int)cds_len; protein_lengths[slot] = (int)protein_len;
        flags[slot] = (int)applied.flags;
        if (protein_len < translated.length) flags[slot] |= DUCKVEP_HAPLOTYPE_FLAG_STOP_TRUNCATED;
        contributor_counts[slot] = 0;
    }
    stream_status = duckvep_haplotype_stream_finish(&stream);
    if (stream_status != DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY) goto stream_error;
    duckvep_haplotype_leaf_t leaf;
    while ((stream_status = duckvep_haplotype_stream_next(&stream, &leaf)) ==
           DUCKVEP_HAPLOTYPE_STREAM_OK) {
        if (leaf.projection_status != DUCKVEP_CDS_EDIT_OK) {
            *status = 200 + (int)leaf.projection_status;
            goto cleanup;
        }
        if (leaf.sequence_status != DUCKVEP_HAPLOTYPE_OK) {
            *status = (int)leaf.sequence_status;
            goto cleanup;
        }
        for (uint32_t id = leaf.carriers.first_call; id;) {
            const duckvep_carrier_call_t *call = duckvep_carriers_call(&stream.carriers, id);
            if (!call || call->key.sample_index >= 3u || call->key.lane > 2u) goto invalid;
            size_t slot = (size_t)call->key.sample_index * 2u + call->key.lane - 1u;
            memcpy(cds + slot * (size_t)*capacity, leaf.cds, leaf.cds_length);
            memcpy(protein + slot * (size_t)*capacity, leaf.protein, leaf.protein_length);
            cds_lengths[slot] = (int)leaf.cds_length;
            protein_lengths[slot] = (int)leaf.protein_length;
            flags[slot] = (int)leaf.flags;
            contributor_counts[slot] = (int)leaf.contributor_count;
            for (size_t i = 0u; i < leaf.contributor_count; i++) {
                uint64_t event = leaf.contributors[i].source.event_id;
                if (!event || event > count) goto invalid;
                contributors[slot * count + i] = (int)event - 1;
            }
            metrics[5] += 1.0;
            id = call->next_leaf;
        }
    }
    if (stream_status != DUCKVEP_HAPLOTYPE_STREAM_DONE) goto stream_error;
    stream_status = duckvep_haplotype_stream_finish(&stream);
    if (stream_status != DUCKVEP_HAPLOTYPE_STREAM_DONE) goto stream_error;
    metrics[0] = stream.carriers.peak_transcripts;
    metrics[1] = stream.carriers.peak_calls;
    metrics[2] = stream.carriers.peak_prefixes;
    metrics[3] = (double)stream.completed_leaves;
    metrics[4] = (double)stream.translated_bases;
    *status = DUCKVEP_HAPLOTYPE_OK;
    goto cleanup;
invalid:
    *status = DUCKVEP_HAPLOTYPE_INVALID_ARG;
    goto cleanup;
stream_error:
    *status = stream_status == DUCKVEP_HAPLOTYPE_STREAM_INPUT_ORDER
        ? 100 + DUCKVEP_CARRIERS_INPUT_ORDER
        : stream_status == DUCKVEP_HAPLOTYPE_STREAM_CARRIER_ERROR
        ? 100 + (int)stream.carrier_error : 300 + (int)stream_status;
cleanup:
    free(storage.events); free(storage.projections); free(storage.leaf_events);
    free(storage.contributors); free(storage.edits); free(storage.blocks);
    free(storage.alleles); free(exon_storage);
    free(cds_scratch); free(protein_scratch);
    free(b.calls); free(b.call_index); free(b.prefixes); free(b.prefix_index);
}
