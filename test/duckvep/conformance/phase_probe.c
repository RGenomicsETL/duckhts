/* Test-only R .C bridge. The production reducer owns raw-GT interpretation;
 * this bridge only converts its fixed result into R-owned integer arrays. */
#include "duckvep_phase.h"
#include "duckvep_haplotype_stream.h"
#include <string.h>

void duckhts_test_vep116_raw_gt(char **gt, int *alt_count, int *count,
    int *status, int *disposition, int *ploidy, int *missing, int *slots,
    int *first, int *second) {
    for (int i = 0; i < *count; i++) {
        duckvep_raw_gt_t parsed;
        status[i] = duckvep_phase_parse_vep116_raw((const uint8_t *)gt[i], strlen(gt[i]),
            (uint32_t)alt_count[i], &parsed);
        disposition[i] = parsed.disposition;
        ploidy[i] = parsed.source_ploidy;
        missing[i] = parsed.source_has_missing;
        slots[i] = (int)parsed.parsed_slots;
        first[i] = parsed.allele_index[0] == UINT32_MAX ? -1 : (int)parsed.allele_index[0];
        second[i] = parsed.allele_index[1] == UINT32_MAX ? -1 : (int)parsed.allele_index[1];
    }
}

/* The finite phase audit has two source records and one all-coding exon per
 * profile. This bridge only binds that fixture to the production record/stream
 * APIs. Output includes observations as well as physical-edit source IDs, so
 * a conditional REF no-op cannot masquerade as an applied upstream variant. */
void duckhts_test_raw_phase_haplotypes(char **reference, int *genomic_start,
    char **gt, int *positions, char **refs, char **alts, int *alt_counts, int *count,
    int *capacity, uint8_t *cds, uint8_t *protein, int *cds_lengths, int *protein_lengths,
    int *sequence_status, int *evidence, int *edit_masks, int *source_indices,
    int *source_evidence, int *errors) {
    size_t length = strlen(reference[0]);
    if (!length || length > 512u || *genomic_start < 1 || *capacity < 512 ||
        *count < 1 || alt_counts[0] < 1 || alt_counts[0] > 2 ||
        alt_counts[1] < 1 || alt_counts[1] > 2) {
        for (int i = 0; i < *count; i++) errors[i] = 1;
        return;
    }
    uint16_t chrom = 0u, exon_count = 1u;
    uint32_t start = (uint32_t)*genomic_start, end = start + (uint32_t)length - 1u;
    uint32_t zero = 0u, one = 1u, cds_length = (uint32_t)length;
    uint64_t offset = 0u;
    int8_t strand = 1;
    duckvep_transcript_model_t model = {0};
    model.transcript_count = 1u; model.chrom_id = &chrom;
    model.start1 = model.cds_start1 = &start; model.end1 = model.cds_end1 = &end;
    model.strand = &strand; model.exon_offset = &zero; model.exon_count = &exon_count;
    duckvep_exon_model_t exons = {0};
    exons.exon_count = 1u; exons.start1 = &start; exons.end1 = &end;
    exons.cdna_start1 = &one; exons.cdna_end1 = &cds_length;
    duckvep_sequence_pool_t sequences = {0};
    sequences.transcript_count = 1u; sequences.cds_offset = &offset;
    sequences.cds_length = &cds_length; sequences.cds_bytes = (const uint8_t *)reference[0];
    sequences.cds_bytes_len = length;
    duckvep_carrier_transcript_t transcript;
    duckvep_carrier_call_t calls[2];
    duckvep_carrier_prefix_t prefixes[16];
    duckvep_carrier_bucket_t tx_index[2], call_index[4], prefix_index[32];
    uint32_t active;
    duckvep_haplotype_stored_event_t events[8];
    duckvep_haplotype_projection_t projections[8];
    duckvep_carrier_event_t leaf_events[8];
    duckvep_haplotype_contributor_t contributors[8];
    duckvep_haplotype_edit_t edits[16];
    duckvep_haplotype_block_t blocks[16];
    uint64_t edit_ids[16];
    uint8_t alleles[128], cds_scratch[512], protein_scratch[512];
    duckvep_haplotype_stream_buffers_t buffers = {
        .carriers = {.transcripts = &transcript, .calls = calls, .prefixes = prefixes,
            .active_transcripts = &active, .transcript_index = tx_index, .call_index = call_index,
            .prefix_index = prefix_index, .transcript_capacity = 1u, .call_capacity = 2u,
            .prefix_capacity = 16u, .transcript_buckets = 2u, .call_buckets = 4u, .prefix_buckets = 32u},
        .events = events, .projections = projections, .alleles = alleles,
        .event_capacity = 8u, .projection_capacity = 8u, .allele_capacity = sizeof(alleles),
        .leaf_events = leaf_events, .contributors = contributors, .edits = edits,
        .edit_event_ids = edit_ids, .blocks = blocks, .leaf_capacity = 8u, .edit_capacity = 16u,
        .cds = cds_scratch, .protein = protein_scratch,
        .cds_capacity = sizeof(cds_scratch), .protein_capacity = sizeof(protein_scratch)};
    for (int profile = 0; profile < *count; profile++) {
        duckvep_haplotype_stream_t stream;
        duckvep_haplotype_stream_status_t status = duckvep_haplotype_stream_init(
            &stream, &model, &exons, &sequences, &buffers);
        for (size_t i = 0u; i < 4u; i++) source_indices[(size_t)profile * 4u + i] = -2;
        size_t alt_base = 0u;
        for (int record = 0; record < 2 && status == DUCKVEP_HAPLOTYPE_STREAM_OK; record++) {
            duckvep_raw_gt_t parsed;
            const char *raw = gt[(size_t)profile * 2u + (size_t)record];
            if (duckvep_phase_parse_vep116_raw((const uint8_t *)raw, strlen(raw),
                    (uint32_t)alt_counts[record], &parsed) != DUCKVEP_RAW_GT_OK) {
                status = DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
                break;
            }
            for (int ordinal = 0; ordinal <= alt_counts[record] + 1; ordinal++) {
                uint32_t index = ordinal > alt_counts[record] ? UINT32_MAX : (uint32_t)ordinal;
                const char *alt = index == UINT32_MAX ? "" : !index ? refs[record]
                    : alts[alt_base + index - 1u];
                size_t ref_len = strlen(refs[record]), alt_len = strlen(alt);
                if (!ref_len || ref_len > UINT16_MAX || alt_len > UINT16_MAX || positions[record] < 1) {
                    status = DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG;
                    break;
                }
                duckvep_haplotype_source_t source = {(uint64_t)record + 1u,
                    (const uint8_t *)refs[record], (const uint8_t *)alt,
                    (uint32_t)positions[record], 0u, (uint16_t)ref_len, (uint16_t)alt_len, index, 1u, 0u};
                status = duckvep_haplotype_stream_begin(&stream, &source);
                if (status == DUCKVEP_HAPLOTYPE_STREAM_OK)
                    status = duckvep_haplotype_stream_project(&stream, 0u);
                if (status == DUCKVEP_HAPLOTYPE_STREAM_OK)
                    status = duckvep_haplotype_stream_push_raw_call(&stream, 0u, 0u, &parsed, 1u);
                if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) break;
            }
            alt_base += (size_t)alt_counts[record];
        }
        if (status == DUCKVEP_HAPLOTYPE_STREAM_OK) status = duckvep_haplotype_stream_finish(&stream);
        if (status != DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY) {
            errors[profile] = 100 + (int)status;
            continue;
        }
        duckvep_haplotype_leaf_t leaf;
        unsigned seen = 0u;
        while ((status = duckvep_haplotype_stream_next(&stream, &leaf)) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
            if (leaf.projection_status != DUCKVEP_CDS_EDIT_OK ||
                (leaf.sequence_status != DUCKVEP_HAPLOTYPE_OK && leaf.sequence_status != DUCKVEP_HAPLOTYPE_CONDITIONAL)) {
                errors[profile] = 1000 + (int)leaf.sequence_status + 100 * (int)leaf.projection_status;
                break;
            }
            for (uint32_t id = leaf.carriers.first_call; id;) {
                const duckvep_carrier_call_t *call = duckvep_carriers_call(&stream.carriers, id);
                if (!call || call->key.lane < 1u || call->key.lane > 2u ||
                    call->key.ploidy != 2u || call->key.sample_index) { errors[profile] = 2; break; }
                unsigned bit = 1u << (call->key.lane - 1u);
                if (seen & bit) { errors[profile] = 2; break; }
                seen |= bit;
                size_t row = (size_t)profile * 2u + call->key.lane - 1u;
                memcpy(cds + row * (size_t)*capacity, leaf.cds, leaf.cds_length);
                memcpy(protein + row * (size_t)*capacity, leaf.protein, leaf.protein_length);
                cds_lengths[row] = (int)leaf.cds_length; protein_lengths[row] = (int)leaf.protein_length;
                sequence_status[row] = (int)leaf.sequence_status; evidence[row] = leaf.evidence_flags;
                for (size_t e = 0u; e < leaf.edit_count; e++) {
                    if (leaf.edit_event_ids[e] < 1u || leaf.edit_event_ids[e] > 2u) {
                        errors[profile] = 3; break;
                    }
                    edit_masks[row] |= 1 << (leaf.edit_event_ids[e] - 1u);
                }
                for (size_t e = 0u; e < leaf.contributor_count; e++) {
                    size_t record = (size_t)leaf.contributors[e].source.event_id - 1u;
                    uint32_t index = leaf.contributors[e].source.allele_index;
                    if (record > 1u || source_indices[row * 2u + record] != -2) {
                        errors[profile] = 4; break;
                    }
                    source_indices[row * 2u + record] = index == UINT32_MAX ? -1 : (int)index;
                    source_evidence[row * 2u + record] = leaf.contributors[e].evidence_flags;
                }
                id = call->next_leaf;
            }
        }
        if (status == DUCKVEP_HAPLOTYPE_STREAM_DONE) status = duckvep_haplotype_stream_finish(&stream);
        if (!errors[profile] && status != DUCKVEP_HAPLOTYPE_STREAM_DONE) errors[profile] = 100 + (int)status;
        if (!errors[profile] && seen != 3u) errors[profile] = 5;
    }
}
