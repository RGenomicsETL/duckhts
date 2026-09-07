/* Model-scoped literal-event replay over sparse carrier paths (INTERNAL).
 * All storage is caller-owned. begin copies REF/ALT once; project adds one
 * candidate transcript; push adds prepared carrier evidence, while push_call
 * interprets a decoded GT. Candidate and phase-domain discovery and output
 * materialization belong to the host query plan.
 *
 * Input is sorted by (chrom_id, pos1, event_id). begin may report a transcript
 * ready before consuming the input: drain next until DONE, then retry begin with
 * the same input. finish uses the same drain protocol. Each next result borrows
 * scratch until the next next/begin/finish call; its carrier list stays valid
 * throughout that transcript's drain. The owner may pause between any two calls.
 *
 * Event/projection/allele rings retain only the oldest still-active genomic
 * window, including intervening events pinned behind a longer transcript.
 * Limits are explicit; exhaustion is an error, never silent loss or growth.
 * Errors latch until reinitialization. Per-leaf projection/edit conflicts are
 * results with complete provenance, not errors that drop an occupied path.
 */
#ifndef DUCKVEP_HAPLOTYPE_STREAM_H
#define DUCKVEP_HAPLOTYPE_STREAM_H

#include "duckvep_carriers.h"
#include "duckvep_delta.h"
#include "duckvep_phase.h"

typedef enum {
    DUCKVEP_HAPLOTYPE_STREAM_OK,
    DUCKVEP_HAPLOTYPE_STREAM_TRANSCRIPT_READY,
    DUCKVEP_HAPLOTYPE_STREAM_DONE,
    DUCKVEP_HAPLOTYPE_STREAM_INVALID_ARG,
    DUCKVEP_HAPLOTYPE_STREAM_INPUT_ORDER,
    DUCKVEP_HAPLOTYPE_STREAM_EVENT_FULL,
    DUCKVEP_HAPLOTYPE_STREAM_PROJECTION_FULL,
    DUCKVEP_HAPLOTYPE_STREAM_ALLELE_FULL,
    DUCKVEP_HAPLOTYPE_STREAM_LEAF_FULL,
    DUCKVEP_HAPLOTYPE_STREAM_EDIT_FULL,
    DUCKVEP_HAPLOTYPE_STREAM_SEQUENCE_FULL,
    DUCKVEP_HAPLOTYPE_STREAM_CARRIER_ERROR,
    DUCKVEP_HAPLOTYPE_STREAM_INTERNAL_ERROR
} duckvep_haplotype_stream_status_t;

typedef struct {
    uint64_t event_id;
    const uint8_t *ref, *alt;
    uint32_t pos1;
    uint16_t chrom_id, ref_len, alt_len;
} duckvep_haplotype_source_t;

typedef struct {
    duckvep_haplotype_source_t source;
    duckvep_event_t prepared;
    uint64_t serial;
    size_t allele_consumed;
    uint32_t projection_begin, projection_count, last_end1;
} duckvep_haplotype_stored_event_t;

typedef struct {
    duckvep_haplotype_edit_t edit;
    uint32_t transcript_index;
    duckvep_cds_edit_status_t status;
} duckvep_haplotype_projection_t;

typedef struct {
    duckvep_haplotype_source_t source;
    duckvep_cds_edit_status_t projection_status;
    uint8_t evidence_flags;
} duckvep_haplotype_contributor_t;

typedef struct {
    int64_t value;
    uint8_t present;
} duckvep_haplotype_phase_set_t;

typedef struct {
    const int32_t *alleles;
    const uint8_t *phase_before; /* NULL means no phase information. */
    uint32_t sample_index, alt_index;
    uint16_t ploidy;
    duckvep_haplotype_phase_set_t phase_set;
    duckvep_phase_policy_t policy;
} duckvep_haplotype_call_t;

typedef struct {
    duckvep_carrier_buffers_t carriers;
    duckvep_haplotype_stored_event_t *events;
    duckvep_haplotype_projection_t *projections;
    uint8_t *alleles;
    uint32_t event_capacity, projection_capacity;
    size_t allele_capacity;
    /* Scratch for one distinct occupied path, reused across its carriers. */
    duckvep_carrier_event_t *leaf_events;
    duckvep_haplotype_contributor_t *contributors;
    duckvep_haplotype_edit_t *edits;
    uint64_t *edit_event_ids; /* Parallel to edits, including after sorting/reversal. */
    duckvep_haplotype_block_t *blocks; /* At most edit_capacity interaction blocks. */
    size_t leaf_capacity, edit_capacity;
    uint8_t *cds, *protein;
    size_t cds_capacity, protein_capacity;
} duckvep_haplotype_stream_buffers_t;

typedef struct {
    duckvep_carrier_leaf_t carriers;
    const duckvep_haplotype_contributor_t *contributors;
    size_t contributor_count;
    size_t edit_count; /* Physical differing islands, not source-event count. */
    /* One source identity per physical edit in ascending reference CDS order.
     * A source may occur more than once or in several blocks. Borrowed only for
     * known sequences; block.edit_begin/edit_count select the corresponding IDs. */
    const uint64_t *edit_event_ids;
    const duckvep_haplotype_block_t *blocks; /* Ascending reference CDS order. */
    size_t block_count;
    const uint8_t *reference_cds; /* Model-owned; block spans borrow this and cds. */
    const uint8_t *cds, *protein; /* Length-delimited views; protein excludes residues after first stop. */
    duckvep_translation_t translation; /* Full translation remains in protein storage for coding facts. */
    size_t cds_length, protein_length;
    uint32_t flags;
    uint8_t evidence_flags; /* OR of contributor evidence, distinct from sequence flags. */
    uint8_t stop_in_displaced_frame; /* First translated stop intersects a frame excursion. */
    /* First failed projection, or edit/rebuild status when projection is OK.
     * Failed paths have no CDS/protein/blocks; all contributors/carriers remain. */
    duckvep_cds_edit_status_t projection_status;
    duckvep_haplotype_status_t sequence_status;
} duckvep_haplotype_leaf_t;

typedef struct {
    duckvep_carriers_t carriers;
    const duckvep_exon_model_t *exons;
    const duckvep_sequence_pool_t *sequences;
    duckvep_haplotype_stream_buffers_t buffers;
    uint32_t event_begin, event_count, projection_begin, projection_count;
    uint32_t current_event, closing;
    size_t allele_begin, allele_count;
    uint64_t serial, last_event_id;
    uint32_t last_pos1;
    uint16_t last_chrom;
    uint8_t have_input, have_current, initialized;
    uint8_t have_phase_policy;
    duckvep_phase_policy_t phase_policy;
    duckvep_haplotype_stream_status_t error;
    duckvep_carriers_status_t carrier_error;
    uint64_t input_events, projected_events, completed_leaves, translated_bases;
    uint32_t peak_events, peak_projections;
    size_t peak_alleles;
} duckvep_haplotype_stream_t;

/* Distinct buffers; immutable model views remain pinned until teardown. */
duckvep_haplotype_stream_status_t duckvep_haplotype_stream_init(
    duckvep_haplotype_stream_t *stream,
    const duckvep_transcript_model_t *transcripts,
    const duckvep_exon_model_t *exons,
    const duckvep_sequence_pool_t *sequences,
    const duckvep_haplotype_stream_buffers_t *buffers);

/* The owner may discard input alleles after OK; input must not alias workspace
 * storage. Event IDs need only
 * be unique, with increasing IDs used to order events at the same coordinate. */
duckvep_haplotype_stream_status_t duckvep_haplotype_stream_begin(
    duckvep_haplotype_stream_t *stream, const duckvep_haplotype_source_t *event);

/* Add candidates in strictly increasing model-local ordinal order. Calls for
 * an earlier candidate may be consumed before projecting the next candidate;
 * an event's entire cohort need not be buffered. Each candidate must overlap
 * the genomic span (an insertion may touch the transcript end). */
duckvep_haplotype_stream_status_t duckvep_haplotype_stream_project(
    duckvep_haplotype_stream_t *stream, uint32_t transcript_index);

/* Add this event to one lane in every candidate transcript, using the prepared
 * projection. No call matrix and no repeated projection across samples. */
duckvep_haplotype_stream_status_t duckvep_haplotype_stream_push(
    duckvep_haplotype_stream_t *stream, const duckvep_carrier_key_t *key,
    uint8_t evidence_flags);

/* Interpret the complete decoded GT for this source ALT and one candidate
 * transcript. The caller's query plan supplies ALL phase sets for this sample
 * and transcript, including those first observed in later records. Arrays are
 * borrowed for the call, sorted absent-first then by signed value, and unique.
 * An empty domain means only the absent/default set. Domains must stay logically
 * identical across the transcript; discovery/broadcast planning belongs to the
 * host, not a second first-party store of the whole query's phase-set catalogue.
 * Homozygous/haploid and wholly unphased evidence is broadcast across the domain.
 * Partial phase affects only its declared set and the unresolved slots within it.
 * Compatibility mode requires the absent-only domain, compacts called slots,
 * and retains missing-call evidence on every lane. Uncertain paths return
 * INPUT_INCOMPLETE with no CDS/protein, including in compatibility mode; this
 * does not certify VEP's conditional sequence output for missing genotypes.
 * Source ALT ordinals must be positive and refer to the current begin event.
 * Missing GT with unknown ploidy is an error; known-ploidy missing slots use -1.
 * One phase policy applies to the stream; changing it is an error.
 * A failure latches even if a broadcast already updated an earlier carrier.
 */
duckvep_haplotype_stream_status_t duckvep_haplotype_stream_push_call(
    duckvep_haplotype_stream_t *stream, uint32_t transcript_index,
    const duckvep_haplotype_call_t *call,
    const duckvep_haplotype_phase_set_t *phase_sets, size_t phase_set_count);

duckvep_haplotype_stream_status_t duckvep_haplotype_stream_finish(
    duckvep_haplotype_stream_t *stream);

/* Each occupied edit prefix is rebuilt/translated once. DONE releases the
 * completed transcript; retry begin/finish to accept input or close the next.
 * Iterate leaf.carriers.first_call using duckvep_carriers_call(&stream->carriers,
 * id), following next_leaf, before calling next again. */
duckvep_haplotype_stream_status_t duckvep_haplotype_stream_next(
    duckvep_haplotype_stream_t *stream, duckvep_haplotype_leaf_t *leaf);

#endif
