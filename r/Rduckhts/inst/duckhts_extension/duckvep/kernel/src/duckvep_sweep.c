/*
 * duckvep_sweep.c — span-aware sorted candidate sweep. See duckvep_sweep.h.
 * No DuckDB/htslib/Arrow; no allocation.
 */
#include "duckvep_sweep.h"

#include "duckvep_event.h"

#include <string.h>

static uint32_t sat_add_u32(uint32_t a, uint32_t b) {
    uint32_t s = (uint32_t)(a + b);
    return s < a ? UINT32_MAX : s;
}

void duckvep_sweep_cursor_init(
    duckvep_sweep_cursor_t            *cursor,
    const duckvep_variant_batch_t     *variants,
    const duckvep_transcript_model_t  *transcripts,
    uint32_t                           halo,
    uint32_t                          *active,
    size_t                             active_cap,
    uint32_t                          *candidates,
    size_t                             candidate_cap) {

    if (cursor == NULL) return;
    memset(cursor, 0, sizeof *cursor);
    cursor->variants = variants;
    cursor->transcripts = transcripts;
    cursor->halo = halo;
    cursor->active = active;
    cursor->active_cap = active_cap;
    cursor->candidates = candidates;
    cursor->candidate_cap = candidate_cap;
    cursor->status = DUCKVEP_OK;

    if (variants == NULL || transcripts == NULL ||
        (active == NULL && active_cap > 0u) ||
        (candidates == NULL && candidate_cap > 0u) ||
        (active != NULL && candidates != NULL && active == candidates)) {
        cursor->status = DUCKVEP_ERR_INVALID_ARG;
    }
}

int duckvep_sweep_cursor_seed(duckvep_sweep_cursor_t *cursor,
                              const uint32_t *seed,
                              size_t seed_count) {
    const duckvep_transcript_model_t *transcripts;
    uint16_t chrom;
    uint32_t event_start;
    uint32_t point_hi;
    size_t chrom_begin;
    size_t chrom_end;
    size_t add;
    size_t i;

    if (cursor == NULL || cursor->status != DUCKVEP_OK || cursor->vi != 0u ||
        cursor->have_chrom || cursor->variants == NULL ||
        cursor->transcripts == NULL || cursor->variants->count == 0u ||
        (seed_count != 0u && seed == NULL)) {
        if (cursor != NULL) cursor->status = DUCKVEP_ERR_INVALID_ARG;
        return 0;
    }
    if (seed_count > cursor->active_cap) {
        cursor->status = DUCKVEP_ERR_RESULT_FULL;
        return 0;
    }

    transcripts = cursor->transcripts;
    chrom = cursor->variants->chrom_id[0];
    event_start = cursor->variants->pos1[0];
    point_hi = sat_add_u32(event_start, cursor->halo);

    chrom_begin = 0u;
    chrom_end = transcripts->transcript_count;
    while (chrom_begin < chrom_end) {
        size_t middle = chrom_begin + (chrom_end - chrom_begin) / 2u;
        if (transcripts->chrom_id[middle] < chrom) chrom_begin = middle + 1u;
        else chrom_end = middle;
    }
    chrom_end = chrom_begin;
    while (chrom_end < transcripts->transcript_count &&
           transcripts->chrom_id[chrom_end] == chrom) {
        chrom_end++;
    }
    add = chrom_begin;
    while (add < chrom_end && transcripts->start1[add] <= point_hi) add++;

    for (i = 0u; i < seed_count; i++) {
        uint32_t tx = seed[i];
        if ((size_t)tx < chrom_begin || (size_t)tx >= add ||
            transcripts->chrom_id[tx] != chrom ||
            sat_add_u32(transcripts->end1[tx], cursor->halo) < event_start) {
            cursor->status = DUCKVEP_ERR_INVALID_ARG;
            return 0;
        }
        cursor->active[i] = tx;
    }

    cursor->ti = chrom_begin;
    cursor->tx_hi = chrom_end;
    cursor->add = add;
    cursor->nact = seed_count;
    cursor->chrom = chrom;
    cursor->have_chrom = 1u;
    return 1;
}

int duckvep_sweep_cursor_next(
    duckvep_sweep_cursor_t *cursor,
    uint32_t                *variant_idx,
    const uint32_t         **tx_indices,
    size_t                  *tx_count) {

    const duckvep_variant_batch_t *variants;
    const duckvep_transcript_model_t *transcripts;
    uint16_t chrom;
    uint32_t event_start;
    uint32_t point_hi;
    uint32_t span_hi;
    size_t k;

    if (variant_idx != NULL) *variant_idx = 0u;
    if (tx_indices != NULL) *tx_indices = NULL;
    if (tx_count != NULL) *tx_count = 0u;
    if (cursor == NULL || variant_idx == NULL || tx_indices == NULL || tx_count == NULL) {
        if (cursor != NULL) cursor->status = DUCKVEP_ERR_INVALID_ARG;
        return -1;
    }
    if (cursor->status != DUCKVEP_OK) return -1;

    variants = cursor->variants;
    transcripts = cursor->transcripts;
    if (cursor->vi >= variants->count) return 0;

    chrom = variants->chrom_id[cursor->vi];
    if (!cursor->have_chrom || chrom != cursor->chrom) {
        if (cursor->have_chrom) cursor->ti = cursor->tx_hi;
        cursor->chrom = chrom;
        cursor->have_chrom = 1u;
        cursor->nact = 0u;

        while (cursor->ti < transcripts->transcript_count &&
               transcripts->chrom_id[cursor->ti] < chrom) {
            cursor->ti++;
        }
        cursor->tx_hi = cursor->ti;
        while (cursor->tx_hi < transcripts->transcript_count &&
               transcripts->chrom_id[cursor->tx_hi] == chrom) {
            cursor->tx_hi++;
        }
        cursor->add = cursor->ti;
    }

    event_start = variants->pos1[cursor->vi];
    point_hi = sat_add_u32(event_start, cursor->halo);
    span_hi = sat_add_u32(
        cursor->events != NULL
            ? duckvep_event_feature_max1(&cursor->events[cursor->vi])
            : duckvep_event_effective_end1_at(variants, cursor->vi),
        cursor->halo);

    /* Persistent point frontier. Event starts are monotone, so each transcript is
     * admitted once and an expired transcript can never overlap a later event. */
    while (cursor->add < cursor->tx_hi &&
           transcripts->start1[cursor->add] <= point_hi) {
        if (cursor->nact >= cursor->active_cap) {
            cursor->status = DUCKVEP_ERR_RESULT_FULL;
            return -1;
        }
        cursor->active[cursor->nact++] = (uint32_t)cursor->add;
        cursor->add++;
    }

    k = 0u;
    while (k < cursor->nact) {
        uint32_t end1 = transcripts->end1[cursor->active[k]];
        if (sat_add_u32(end1, cursor->halo) < event_start) {
            cursor->active[k] = cursor->active[--cursor->nact];
        } else {
            k++;
        }
    }

    *variant_idx = (uint32_t)cursor->vi;
    if (span_hi == point_hi) {
        *tx_indices = cursor->active;
        *tx_count = cursor->nact;
    } else {
        size_t future_hi = cursor->add;
        size_t ncand = cursor->nact;

        while (future_hi < cursor->tx_hi &&
               transcripts->start1[future_hi] <= span_hi) {
            future_hi++;
        }
        if (ncand > cursor->candidate_cap ||
            future_hi - cursor->add > cursor->candidate_cap - ncand) {
            cursor->status = DUCKVEP_ERR_RESULT_FULL;
            return -1;
        }
        if (ncand > 0u) {
            memcpy(cursor->candidates, cursor->active,
                   ncand * sizeof *cursor->candidates);
        }
        for (k = cursor->add; k < future_hi; k++) {
            cursor->candidates[ncand++] = (uint32_t)k;
        }
        *tx_indices = cursor->candidates;
        *tx_count = ncand;
    }
    cursor->vi++;
    return 1;
}

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
    duckvep_status_t                 *status) {

    duckvep_sweep_cursor_t cursor;
    size_t pairs = 0u;
    int rc;

    if (status != NULL) *status = DUCKVEP_OK;
    duckvep_sweep_cursor_init(&cursor, variants, transcripts, halo,
                              active, active_cap, candidates, candidate_cap);
    if (cursor.status != DUCKVEP_OK) {
        if (status != NULL) *status = cursor.status;
        return 0u;
    }

    for (;;) {
        uint32_t vi;
        const uint32_t *tx_indices_local;
        size_t tx_count_local;
        size_t j;

        rc = duckvep_sweep_cursor_next(&cursor, &vi, &tx_indices_local,
                                       &tx_count_local);
        if (rc <= 0) break;
        if (sink == NULL) {
            pairs += tx_count_local;
            continue;
        }
        for (j = 0u; j < tx_count_local; j++) {
            if (!sink(vi, tx_indices_local[j], ctx)) {
                if (status != NULL) *status = cursor.status;
                return pairs;
            }
            pairs++;
        }
    }

    if (status != NULL) *status = cursor.status;
    return pairs;
}
