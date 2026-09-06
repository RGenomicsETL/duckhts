/* Sparse, model-scoped carrier paths over projected event identities (INTERNAL).
 * No allocation, DuckDB, HTSlib, phase-policy inference or allele decoding.
 *
 * The owner pins one immutable transcript model and owns all buffers. Projected
 * alleles/provenance stay in its event store until a transcript is drained.
 * This index stores only event IDs: identical prefixes share nodes across
 * samples/phase sets/lanes; reference paths require no entry. Input is ordered
 * by (chrom_id, pos1, event_id), with a unique event ID per source event. Calls
 * for one event may span arbitrarily many input batches.
 *
 * advance/finish report a completed transcript before accepting further input.
 * Drain its distinct leaves, materialize each leaf and its carrier list, then
 * release it. All returned IDs/views expire at release. Releasing returns the
 * transcript's slots to the pools; memory does not grow with past transcripts.
 */
#ifndef DUCKVEP_CARRIERS_H
#define DUCKVEP_CARRIERS_H

#include "duckvep_kernel.h"

typedef enum {
    DUCKVEP_CARRIERS_OK,
    DUCKVEP_CARRIERS_TRANSCRIPT_READY,
    DUCKVEP_CARRIERS_DONE,
    DUCKVEP_CARRIERS_INVALID_ARG,
    DUCKVEP_CARRIERS_INPUT_ORDER,
    DUCKVEP_CARRIERS_TRANSCRIPT_FULL,
    DUCKVEP_CARRIERS_CALL_FULL,
    DUCKVEP_CARRIERS_PREFIX_FULL,
    DUCKVEP_CARRIERS_OUTPUT_FULL,
    DUCKVEP_CARRIERS_DUPLICATE_CALL,
    DUCKVEP_CARRIERS_PLOIDY_CONFLICT
} duckvep_carriers_status_t;

typedef struct {
    uint32_t sample_index;
    int64_t phase_set;
    uint16_t lane;                  /* One-based, at most ploidy. */
    uint16_t ploidy;
    uint8_t phase_set_present;      /* Absent and a present PS=0 are distinct. */
} duckvep_carrier_key_t;

/* Slots and hash buckets are caller-owned implementation storage. A zero ID
 * means no slot; live IDs are one-based offsets into their respective arrays. */
typedef struct {
    uint64_t hash;
    uint32_t id;
} duckvep_carrier_bucket_t;

typedef struct {
    uint32_t transcript_index;
    uint32_t first_call, first_prefix, next_free;
} duckvep_carrier_transcript_t;

typedef struct {
    duckvep_carrier_key_t key;
    uint32_t transcript, leaf, next_transcript, next_leaf, previous_leaf, next_free;
} duckvep_carrier_call_t;

typedef struct {
    uint64_t event_id;
    uint32_t transcript, parent, depth, first_call, call_count, next_transcript, next_free;
} duckvep_carrier_prefix_t;

typedef struct {
    duckvep_carrier_transcript_t *transcripts;
    duckvep_carrier_call_t *calls;
    duckvep_carrier_prefix_t *prefixes;
    uint32_t *active_transcripts;   /* Min-heap by model chromosome, end, ordinal. */
    duckvep_carrier_bucket_t *transcript_index, *call_index, *prefix_index;
    uint32_t transcript_capacity, call_capacity, prefix_capacity;
    uint32_t transcript_buckets, call_buckets, prefix_buckets;
} duckvep_carrier_buffers_t;

typedef struct {
    const duckvep_transcript_model_t *model;
    duckvep_carrier_buffers_t buffers;
    uint32_t free_transcript, free_call, free_prefix;
    uint32_t transcript_count, call_count, prefix_count;
    uint32_t peak_transcripts, peak_calls, peak_prefixes;
    uint32_t closing, next_leaf;
    uint32_t pos1, pending_pos1;
    uint64_t event_id, pending_event_id;
    uint16_t chrom, pending_chrom;
    uint8_t initialized, have_event, pending, pending_eof, finished, drained;
} duckvep_carriers_t;

typedef struct {
    uint32_t id, transcript_index, event_count, first_call, call_count;
} duckvep_carrier_leaf_t;

/* Buffer arrays must be distinct. Each hash size is a power of two at least
 * twice its pool capacity, or both sizes are zero. Init clears all storage. */
duckvep_carriers_status_t duckvep_carriers_init(
    duckvep_carriers_t *stream, const duckvep_transcript_model_t *model,
    const duckvep_carrier_buffers_t *buffers);

/* Retry the same tuple after draining each TRANSCRIPT_READY. The completed
 * model-local ordinal is written to transcript_index. Equal tuples continue
 * the same event across batches; decreasing tuples fail without mutation. */
duckvep_carriers_status_t duckvep_carriers_advance(
    duckvep_carriers_t *stream, uint16_t chrom, uint32_t pos1, uint64_t event_id,
    uint32_t *transcript_index);
duckvep_carriers_status_t duckvep_carriers_finish(
    duckvep_carriers_t *stream, uint32_t *transcript_index);

/* Append the current non-reference event to one carrier. Capacity failures,
 * duplicate calls and changed ploidy for an existing carrier key leave all
 * path/pool state unchanged.
 * Phase interpretation belongs to the consumer; this preserves explicit keys. */
duckvep_carriers_status_t duckvep_carriers_push(
    duckvep_carriers_t *stream, uint32_t transcript_index,
    const duckvep_carrier_key_t *key);

/* Each distinct occupied prefix is emitted once, including a prefix that is
 * both a completed carrier path and an ancestor of a longer carrier path. */
duckvep_carriers_status_t duckvep_carriers_next_leaf(
    duckvep_carriers_t *stream, duckvep_carrier_leaf_t *leaf);
/* Event IDs are returned in arrival order. No output is written on OUTPUT_FULL. */
duckvep_carriers_status_t duckvep_carriers_leaf_events(
    const duckvep_carriers_t *stream, uint32_t leaf, uint64_t *events,
    size_t capacity, size_t *required);
/* Iterate from leaf.first_call through next_leaf. Borrowed until release. */
const duckvep_carrier_call_t *duckvep_carriers_call(
    const duckvep_carriers_t *stream, uint32_t call);
/* Requires next_leaf to have returned DONE. No other transcript is released. */
duckvep_carriers_status_t duckvep_carriers_release(duckvep_carriers_t *stream);

#endif
