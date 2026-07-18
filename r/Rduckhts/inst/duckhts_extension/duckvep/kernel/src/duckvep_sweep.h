/*
 * duckvep_sweep.h — candidate (event x interval) generation via a sorted,
 * resumable sweep. INTERNAL: exposed for white-box property tests.
 *
 * Preconditions:
 *   - events sorted ascending by (chrom_id, pos1/start1);
 *   - intervals sorted ascending by (chrom_id, start1).
 * An interval I is a candidate for event E iff, on the same chromosome,
 *
 *   I.start1 <= E.effective_end1 + halo && I.end1 + halo >= E.start1.
 *
 * For production annotation, effective_end1 is the supplied end1 for the
 * single-interval SV/CNV events and the allele-normalized differing-region end for
 * small variants. Raw pos1 remains the monotone sweep key; shared VCF padding never
 * becomes the topology start. Geometry-only white-box callers may omit variant_kind
 * to sweep an explicitly supplied interval.
 *
 * Event ends need not be monotone. The cursor therefore keeps a persistent point
 * frontier at raw E.start1+halo, then appends the sorted interval tail whose starts
 * lie in (raw E.start1+halo, E.effective_end1+halo] only for the current span. That tail is not
 * retained in the active set, so one chromosome-wide CNV cannot poison the hot
 * path for later point variants. Point events retain the zero-copy active slice;
 * proper spans use caller-owned `candidates` scratch for active + tail.
 *
 * Both scratch arrays are caller-owned and allocation-free. They must be distinct.
 */
#ifndef DUCKVEP_SWEEP_H
#define DUCKVEP_SWEEP_H

#include "duckvep_kernel.h"
#include "duckvep_event.h"

#include <stddef.h>
#include <stdint.h>

typedef struct duckvep_sweep_cursor {
    const duckvep_variant_batch_t    *variants;
    const duckvep_event_t            *events; /* optional once-prepared geometry */
    const uint16_t                    *interval_chrom_id;
    const uint32_t                    *interval_start1;
    const uint32_t                    *interval_end1;
    size_t                             interval_count;
    uint32_t                          halo;
    uint32_t                         *active;
    size_t                            active_cap;
    uint32_t                         *candidates;
    size_t                            candidate_cap;

    size_t                            vi;
    size_t                            interval_pos;
    size_t                            interval_hi;
    size_t                            add;
    size_t                            nact;
    uint16_t                          chrom;
    uint8_t                           have_chrom;
    duckvep_status_t                  status;
} duckvep_sweep_cursor_t;

void duckvep_sweep_cursor_init(
    duckvep_sweep_cursor_t            *cursor,
    const duckvep_variant_batch_t     *variants,
    const duckvep_transcript_model_t  *transcripts,
    uint32_t                           halo,
    uint32_t                          *active,
    size_t                             active_cap,
    uint32_t                          *candidates,
    size_t                             candidate_cap);

/* The sweep only needs sorted contig/start/end columns. This constructor lets
 * the regulation lane reuse the exact candidate algorithm without fabricating
 * transcript records or duplicating interval traversal logic. */
void duckvep_interval_sweep_cursor_init(
    duckvep_sweep_cursor_t        *cursor,
    const duckvep_variant_batch_t *variants,
    const uint16_t                *interval_chrom_id,
    const uint32_t                *interval_start1,
    const uint32_t                *interval_end1,
    size_t                         interval_count,
    uint32_t                       halo,
    uint32_t                      *active,
    size_t                         active_cap,
    uint32_t                      *candidates,
    size_t                         candidate_cap);

/* Start a fresh cursor at the first event without scanning every earlier
 * interval on the chromosome. `seed` is the unique set that overlaps the
 * first event's point window; it may come from cgranges or another interval
 * index. Subsequent events use the ordinary forward sweep. */
int duckvep_sweep_cursor_seed(
    duckvep_sweep_cursor_t *cursor,
    const uint32_t         *seed,
    size_t                  seed_count);

/* Advance to the next event and expose its exact candidate interval slice.
 * Returns 1 and fills outputs, 0 at end, -1 on error (cursor->status). The slice
 * is borrowed and valid only until the next cursor call. Empty slices are valid. */
int duckvep_sweep_cursor_next(
    duckvep_sweep_cursor_t *cursor,
    uint32_t                *variant_idx,
    const uint32_t         **interval_indices,
    size_t                  *interval_count);

typedef int (*duckvep_pair_sink)(uint32_t variant_idx, uint32_t tx_idx, void *ctx);

size_t duckvep_sweep_candidates(
    const duckvep_variant_batch_t    *variants,
    const duckvep_transcript_model_t *transcripts,
    uint32_t                          halo,
    uint32_t                         *active,
    size_t                            active_cap,
    uint32_t                         *candidates,
    size_t                            candidate_cap,
    duckvep_pair_sink                 sink,
    void                             *ctx,
    duckvep_status_t                 *status);

#endif /* DUCKVEP_SWEEP_H */
