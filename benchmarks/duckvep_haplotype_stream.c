/* Benchmark adapter over the production stream, not another replay engine.
 * The registered four-event fixture is repeated across transcript clusters.
 * Allocation/model setup and the complete validation pass are outside timing. */
#define _POSIX_C_SOURCE 200809L
#include "duckvep_haplotype_stream.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

enum { SECONDS, WORKSPACE, MODEL_BYTES, RECORDS, PROJECTIONS, CALLS,
    PEAK_TRANSCRIPTS, PEAK_CARRIERS, PEAK_PREFIXES, PREFIXES_CREATED, LEAVES,
    TRANSLATED_BASES, OUTPUT_CARRIERS, CDS_BYTES, PROTEIN_BYTES, EDITS,
    CONTRIBUTORS, BLOCKS, PEAK_EVENTS, PEAK_PROJECTIONS, PEAK_ALLELES, METRICS };

static uint32_t buckets(uint32_t n) {
    uint32_t b = 1u;
    while (b < n * 2u) b *= 2u;
    return b;
}

static double seconds(void) {
    struct timespec t;
    return clock_gettime(CLOCK_MONOTONIC, &t) ? -1.0 : t.tv_sec + t.tv_nsec / 1e9;
}

static int drain(duckvep_haplotype_stream_t *stream, uint32_t samples, uint32_t overlap,
    const int *genotypes, char **expected_cds, char **expected_protein,
    uint8_t *seen_transcripts, uint8_t *seen_carriers, double *metrics, int verify) {
    duckvep_haplotype_leaf_t leaf;
    duckvep_haplotype_stream_status_t status;
    while ((status = duckvep_haplotype_stream_next(stream, &leaf)) == DUCKVEP_HAPLOTYPE_STREAM_OK) {
        if (leaf.projection_status || leaf.sequence_status || !leaf.cds || !leaf.protein) return 0;
        metrics[OUTPUT_CARRIERS] += leaf.carriers.call_count;
        metrics[CDS_BYTES] += leaf.cds_length;
        metrics[PROTEIN_BYTES] += leaf.protein_length;
        metrics[EDITS] += leaf.edit_count;
        metrics[CONTRIBUTORS] += leaf.contributor_count;
        metrics[BLOCKS] += leaf.block_count;
        if (!verify) continue;
        uint32_t tx = leaf.carriers.transcript_index, cluster = tx / overlap;
        if (tx >= stream->carriers.model->transcript_count) return 0;
        unsigned mask = 0u;
        for (size_t i = 0u; i < leaf.contributor_count; i++) {
            uint64_t id = leaf.contributors[i].source.event_id;
            if (id < (uint64_t)cluster * 4u + 1u || id > (uint64_t)cluster * 4u + 4u) return 0;
            unsigned bit = 1u << (id - (uint64_t)cluster * 4u - 1u);
            if (mask & bit || leaf.contributors[i].evidence_flags != DUCKVEP_CARRIER_CALLED) return 0;
            mask |= bit;
        }
        int path = mask == 3u ? 0 : mask == 12u ? 1 : mask == 15u ? 2 : -1;
        if (path < 0 || seen_transcripts[tx] & (1u << path)) return 0;
        seen_transcripts[tx] |= (uint8_t)(1u << path);
        if (leaf.cds_length != strlen(expected_cds[path]) ||
            memcmp(leaf.cds, expected_cds[path], leaf.cds_length) ||
            leaf.protein_length != strlen(expected_protein[path]) ||
            memcmp(leaf.protein, expected_protein[path], leaf.protein_length)) return 0;
        memset(seen_carriers, 0, samples * 2u);
        uint32_t count = 0u;
        for (uint32_t id = leaf.carriers.first_call; id;) {
            const duckvep_carrier_call_t *call = duckvep_carriers_call(&stream->carriers, id);
            if (!call || call->key.sample_index >= samples || call->key.lane < 1u ||
                call->key.lane > 2u || call->key.ploidy != 2u ||
                !call->key.phase_set_present || call->key.phase_set != 10) return 0;
            uint32_t slot = call->key.sample_index * 2u + call->key.lane - 1u;
            if (seen_carriers[slot]) return 0;
            seen_carriers[slot] = 1u;
            id = call->next_leaf;
            count++;
        }
        if (count != leaf.carriers.call_count) return 0;
        for (uint32_t sample = 0u; sample < samples; sample++) {
            for (unsigned lane = 0u; lane < 2u; lane++) {
                unsigned expected = 0u;
                for (unsigned event = 0u; event < 4u; event++)
                    if (genotypes[event * 8u + (sample % 4u) * 2u + lane]) expected |= 1u << event;
                if (seen_carriers[sample * 2u + lane] != (expected == mask)) return 0;
            }
        }
    }
    return status == DUCKVEP_HAPLOTYPE_STREAM_DONE;
}

void duckhts_bench_haplotype_stream(char **reference, int *transcripts, int *samples,
    int *overlap, int *positions, char **refs, char **alts, int *genotypes,
    char **expected_cds, char **expected_protein, double *metrics, int *result) {
    *result = 1;
    memset(metrics, 0, METRICS * sizeof(*metrics));
    if (*transcripts < 1 || *transcripts > 100000 || *samples < 4 || *samples > 4096 ||
        *samples % 4 || *overlap < 1 || *overlap > 512 || *transcripts % *overlap ||
        strlen(*reference) != 180u) return;
    uint32_t n = (uint32_t)*transcripts, c = (uint32_t)*samples, o = (uint32_t)*overlap;
    for (unsigned i = 0u; i < 4u; i++) {
        if (positions[i] < 1 || positions[i] > 178 || (i && positions[i] <= positions[i - 1]) ||
            strlen(refs[i]) > 2u || strlen(alts[i]) > 2u) return;
        for (unsigned j = 0u; j < 8u; j++)
            if (genotypes[i * 8u + j] != 0 && genotypes[i * 8u + j] != 1) return;
    }
    duckvep_transcript_model_t model = {0};
    duckvep_exon_model_t exons = {0};
    duckvep_sequence_pool_t sequences = {0};
    duckvep_haplotype_stream_buffers_t b = {0};
    duckvep_haplotype_stream_t stream;
    uint16_t *chrom = NULL, *exon_count = NULL;
    int8_t *strand = NULL;
    uint32_t *starts = NULL, *ends = NULL, *offsets = NULL, *lengths = NULL, *cdna_starts = NULL;
    uint64_t *cds_offsets = NULL;
    uint8_t *seen_transcripts = NULL, *seen_carriers = NULL;
    size_t model_bytes = sizeof model + sizeof exons + sizeof sequences;
    size_t workspace = sizeof stream + sizeof b;
    b.carriers.transcript_capacity = o; b.carriers.transcript_buckets = buckets(o);
    b.carriers.call_capacity = 2u * c * o; b.carriers.call_buckets = buckets(2u * c * o);
    b.carriers.prefix_capacity = 16u * o; b.carriers.prefix_buckets = buckets(16u * o);
    b.event_capacity = 5u; b.projection_capacity = 4u * o + 1u;
    b.allele_capacity = 64u; b.leaf_capacity = b.edit_capacity = 4u;
    b.cds_capacity = 184u; b.protein_capacity = 64u;
#define MODEL_ARRAYS(X) \
    X(chrom, n) X(exon_count, n) X(strand, n) X(starts, n) X(ends, n) \
    X(offsets, n) X(lengths, n) X(cdna_starts, n) X(cds_offsets, n)
#define WORK_ARRAYS(X) \
    X(b.carriers.transcripts, o) X(b.carriers.active_transcripts, o) \
    X(b.carriers.transcript_index, b.carriers.transcript_buckets) \
    X(b.carriers.calls, b.carriers.call_capacity) X(b.carriers.call_index, b.carriers.call_buckets) \
    X(b.carriers.prefixes, b.carriers.prefix_capacity) X(b.carriers.prefix_index, b.carriers.prefix_buckets) \
    X(b.events, b.event_capacity) X(b.projections, b.projection_capacity) X(b.alleles, b.allele_capacity) \
    X(b.leaf_events, 4u) X(b.contributors, 4u) X(b.edits, 4u) X(b.edit_event_ids, 4u) \
    X(b.blocks, 4u) X(b.cds, b.cds_capacity) X(b.protein, b.protein_capacity)
#define COUNT_MODEL(p, count) model_bytes += (count) * sizeof(*(p));
#define COUNT_WORK(p, count) workspace += (count) * sizeof(*(p));
    MODEL_ARRAYS(COUNT_MODEL)
    WORK_ARRAYS(COUNT_WORK)
#define ALLOCATE(p, count) if (!((p) = calloc((count), sizeof(*(p))))) goto cleanup;
    MODEL_ARRAYS(ALLOCATE)
    WORK_ARRAYS(ALLOCATE)
    ALLOCATE(seen_transcripts, n)
    ALLOCATE(seen_carriers, 2u * c)
    for (uint32_t i = 0u; i < n; i++) {
        starts[i] = (i / o) * 1000u + 100u; ends[i] = starts[i] + 179u;
        offsets[i] = i; lengths[i] = 180u; cdna_starts[i] = 1u;
        exon_count[i] = 1u; strand[i] = 1;
    }
    model.transcript_count = n; model.chrom_id = chrom; model.strand = strand;
    model.start1 = model.cds_start1 = starts; model.end1 = model.cds_end1 = ends;
    model.exon_offset = offsets; model.exon_count = exon_count;
    exons.exon_count = n; exons.start1 = starts; exons.end1 = ends;
    exons.cdna_start1 = cdna_starts; exons.cdna_end1 = lengths;
    sequences.transcript_count = n; sequences.cds_bytes = (const uint8_t *)*reference;
    sequences.cds_bytes_len = 180u; sequences.cds_offset = cds_offsets; sequences.cds_length = lengths;
    /* First pass validates every occupied path and every carrier against the
     * separately rebuilt fixture. Second pass times only initialized execution
     * and a count/length sink, not SQL alignment or nested-row materialization. */
    for (int pass = 0; pass < 2; pass++) {
        memset(metrics, 0, METRICS * sizeof(*metrics));
        if (duckvep_haplotype_stream_init(&stream, &model, &exons, &sequences, &b)) goto cleanup;
        double start = seconds();
        if (start < 0) goto cleanup;
        for (uint32_t cluster = 0u; cluster < n / o; cluster++) {
            for (unsigned event = 0u; event < 4u; event++) {
                duckvep_haplotype_source_t source = {(uint64_t)cluster * 4u + event + 1u,
                    (const uint8_t *)refs[event], (const uint8_t *)alts[event],
                    cluster * 1000u + 99u + (uint32_t)positions[event], 0u,
                    (uint16_t)strlen(refs[event]), (uint16_t)strlen(alts[event]), 0u, 0u};
                duckvep_haplotype_stream_status_t status;
                while ((status = duckvep_haplotype_stream_begin(&stream, &source)) ==
                        DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY)
                    if (!drain(&stream, c, o, genotypes, expected_cds, expected_protein,
                            seen_transcripts, seen_carriers, metrics, !pass)) goto cleanup;
                if (status != DUCKVEP_HAPLOTYPE_STREAM_OK) goto cleanup;
                for (uint32_t tx = cluster * o; tx < (cluster + 1u) * o; tx++) {
                    if (duckvep_haplotype_stream_project(&stream, tx)) goto cleanup;
                    for (uint32_t sample = 0u; sample < c; sample++) {
                        const int *gt = genotypes + event * 8u + (sample % 4u) * 2u;
                        int32_t alleles[] = {gt[0], gt[1]};
                        const uint8_t phase[] = {1u, 1u};
                        duckvep_haplotype_phase_set_t ps = {10, 1u};
                        duckvep_haplotype_call_t call = {alleles, phase, sample, 1u, 2u, ps,
                            DUCKVEP_PHASE_STRICT};
                        uint32_t before = stream.carriers.prefix_count;
                        if (duckvep_haplotype_stream_push_call(&stream, tx, &call, &ps, 1u)) goto cleanup;
                        metrics[PREFIXES_CREATED] += stream.carriers.prefix_count - before;
                        metrics[CALLS]++;
                    }
                }
            }
        }
        duckvep_haplotype_stream_status_t status;
        while ((status = duckvep_haplotype_stream_finish(&stream)) == DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY)
            if (!drain(&stream, c, o, genotypes, expected_cds, expected_protein,
                    seen_transcripts, seen_carriers, metrics, !pass)) goto cleanup;
        if (status != DUCKVEP_HAPLOTYPE_STREAM_DONE) goto cleanup;
        metrics[SECONDS] = seconds() - start;
        if (!pass) for (uint32_t i = 0u; i < n; i++) if (seen_transcripts[i] != 7u) goto cleanup;
    }
    metrics[WORKSPACE] = (double)workspace; metrics[MODEL_BYTES] = (double)model_bytes;
    metrics[RECORDS] = (double)stream.input_events; metrics[PROJECTIONS] = (double)stream.projected_events;
    metrics[PEAK_TRANSCRIPTS] = stream.carriers.peak_transcripts;
    metrics[PEAK_CARRIERS] = stream.carriers.peak_calls; metrics[PEAK_PREFIXES] = stream.carriers.peak_prefixes;
    metrics[LEAVES] = (double)stream.completed_leaves; metrics[TRANSLATED_BASES] = (double)stream.translated_bases;
    metrics[PEAK_EVENTS] = stream.peak_events; metrics[PEAK_PROJECTIONS] = stream.peak_projections;
    metrics[PEAK_ALLELES] = (double)stream.peak_alleles;
    if (metrics[SECONDS] < 0) goto cleanup;
    *result = 0;
cleanup:
#define FREE(p, count) free(p);
    MODEL_ARRAYS(FREE)
    WORK_ARRAYS(FREE)
    free(seen_transcripts); free(seen_carriers);
}
