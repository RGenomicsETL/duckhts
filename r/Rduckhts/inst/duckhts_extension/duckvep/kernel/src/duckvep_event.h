/*
 * duckvep_event.h — execution event geometry and structural semantics.
 * INTERNAL, header-only: the public batch remains a borrowed SoA, while the hot
 * path loads one compact event value per ALT.
 *
 * Single-interval SV/CNV rows use their supplied inclusive [pos1,end1] span.
 * Small variants carry TWO geometries: the raw VCF anchor/span, and the VEP-style
 * differing region derived by trimming shared REF/ALT prefix + suffix. The sweep
 * keeps raw pos1 as the monotone key, but topology/splice predicates consume the
 * differing region so padded VCF anchors do not create fake CDS/splice overlap.
 */
#ifndef DUCKVEP_EVENT_H
#define DUCKVEP_EVENT_H

#include "duckvep_kernel.h"

#include <stddef.h>
#include <stdint.h>

typedef struct duckvep_event {
    uint16_t chrom_id;
    uint32_t raw_start1;
    uint32_t raw_end1;
    uint32_t start1;          /* differing-region topology coordinate */
    uint32_t end1;            /* inclusive; pure insertions use a point */
    uint16_t ref_diff_offset; /* offset inside REF where differing region starts */
    uint16_t alt_diff_offset; /* offset inside ALT where differing region starts */
    uint16_t ref_diff_length;
    uint16_t alt_diff_length;
    uint8_t  interbase;       /* pure insertion after trimming */
    uint8_t  kind;
    uint8_t  sv_type;
    uint8_t  copy_change;
} duckvep_event_t;

static inline uint32_t duckvep_event_sat_add_u32_u16(uint32_t x, uint16_t y) {
    return x > UINT32_MAX - (uint32_t)y ? UINT32_MAX : x + (uint32_t)y;
}

static inline int duckvep_event_allele_slices_ok(
    const duckvep_variant_batch_t *batch,
    size_t                         idx) {

    uint64_t pool;
    uint64_t roff, rlen, aoff, alen;
    if (batch == NULL || batch->allele_bytes == NULL || batch->ref_offset == NULL ||
        batch->alt_offset == NULL || batch->ref_length == NULL ||
        batch->alt_length == NULL) {
        return 0;
    }
    pool = (uint64_t)batch->allele_bytes_len;
    roff = (uint64_t)batch->ref_offset[idx];
    rlen = (uint64_t)batch->ref_length[idx];
    aoff = (uint64_t)batch->alt_offset[idx];
    alen = (uint64_t)batch->alt_length[idx];
    return roff <= pool && rlen <= pool - roff &&
           aoff <= pool && alen <= pool - aoff &&
           rlen <= UINT16_MAX && alen <= UINT16_MAX;
}

static inline void duckvep_event_load_raw_interval(
    const duckvep_variant_batch_t *batch,
    size_t                         idx,
    duckvep_event_t               *event) {

    event->raw_start1 = batch->pos1[idx];
    event->raw_end1 = batch->end1[idx];
    event->start1 = batch->pos1[idx];
    event->end1 = batch->end1[idx];
    event->ref_diff_offset = 0u;
    event->alt_diff_offset = 0u;
    event->ref_diff_length = 0u;
    event->alt_diff_length = 0u;
    event->interbase = 0u;
}

static inline void duckvep_event_load_small_differing_region(
    const duckvep_variant_batch_t *batch,
    size_t                         idx,
    duckvep_event_t               *event) {

    const uint8_t *ref;
    const uint8_t *alt;
    uint16_t ref_len;
    uint16_t alt_len;
    uint16_t prefix = 0u;
    uint16_t suffix = 0u;
    uint16_t ref_rem;
    uint16_t alt_rem;
    uint16_t ref_diff_len;
    uint16_t alt_diff_len;
    uint32_t diff_start;

    duckvep_event_load_raw_interval(batch, idx, event);
    event->end1 = event->start1; /* small-variant fallback when alleles are absent */
    if (!duckvep_event_allele_slices_ok(batch, idx)) return;

    ref_len = batch->ref_length[idx];
    alt_len = batch->alt_length[idx];
    ref = batch->allele_bytes + batch->ref_offset[idx];
    alt = batch->allele_bytes + batch->alt_offset[idx];

    while (prefix < ref_len && prefix < alt_len && ref[prefix] == alt[prefix]) {
        prefix++;
    }
    ref_rem = (uint16_t)(ref_len - prefix);
    alt_rem = (uint16_t)(alt_len - prefix);
    while (suffix < ref_rem && suffix < alt_rem &&
           ref[(uint16_t)(ref_len - 1u - suffix)] ==
           alt[(uint16_t)(alt_len - 1u - suffix)]) {
        suffix++;
    }

    ref_diff_len = (uint16_t)(ref_len - prefix - suffix);
    alt_diff_len = (uint16_t)(alt_len - prefix - suffix);
    diff_start = duckvep_event_sat_add_u32_u16(batch->pos1[idx], prefix);

    event->ref_diff_offset = prefix;
    event->alt_diff_offset = prefix;
    event->ref_diff_length = ref_diff_len;
    event->alt_diff_length = alt_diff_len;
    event->interbase = (uint8_t)(ref_diff_len == 0u && alt_diff_len > 0u);
    if (event->interbase) {
        /* A pure insertion has no REF bases to classify. Keep the interbase edit facts
         * but use VCF's anchor-side point for topology; VEP treats exon-end insertions
         * as exonic splice-region/coding edits, not as substitutions of the first
         * intronic donor base to the right. */
        event->start1 = prefix > 0u
            ? duckvep_event_sat_add_u32_u16(batch->pos1[idx], (uint16_t)(prefix - 1u))
            : batch->pos1[idx];
        event->end1 = event->start1;
    } else {
        event->start1 = diff_start;
        if (ref_diff_len > 0u) {
            event->end1 = duckvep_event_sat_add_u32_u16(diff_start,
                                                        (uint16_t)(ref_diff_len - 1u));
        } else {
            event->end1 = diff_start;
        }
    }
}

static inline uint32_t duckvep_event_effective_end1_at(
    const duckvep_variant_batch_t *batch, size_t idx) {
    duckvep_event_t event;
    /* Geometry-only white-box sweep callers may omit variant_kind and thereby
     * request the supplied interval directly. Production annotation validates
     * variant_kind and uses full intervals only for structural events. */
    if (batch->variant_kind == NULL ||
        batch->variant_kind[idx] == (uint8_t)DUCKVEP_KIND_SV) {
        return batch->end1[idx];
    }
    duckvep_event_load_small_differing_region(batch, idx, &event);
    return event.end1;
}

static inline uint8_t duckvep_event_sv_type_at(
    const duckvep_variant_batch_t *batch, size_t idx) {
    if (batch->variant_kind == NULL) {
        return (uint8_t)DUCKVEP_SV_UNKNOWN;
    }
    if (batch->variant_kind[idx] != (uint8_t)DUCKVEP_KIND_SV) {
        return (uint8_t)DUCKVEP_SV_NONE;
    }
    return batch->sv_type != NULL ? batch->sv_type[idx]
                                  : (uint8_t)DUCKVEP_SV_UNKNOWN;
}

static inline uint8_t duckvep_event_copy_change_at(
    const duckvep_variant_batch_t *batch, size_t idx) {
    uint8_t sv_type;
    if (batch->variant_kind == NULL ||
        batch->variant_kind[idx] != (uint8_t)DUCKVEP_KIND_SV) {
        return (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN;
    }
    sv_type = duckvep_event_sv_type_at(batch, idx);
    if (sv_type == (uint8_t)DUCKVEP_SV_DELETION) {
        return (uint8_t)DUCKVEP_COPY_CHANGE_LOSS;
    }
    if (sv_type == (uint8_t)DUCKVEP_SV_DUPLICATION ||
        sv_type == (uint8_t)DUCKVEP_SV_TANDEM_DUPLICATION) {
        return (uint8_t)DUCKVEP_COPY_CHANGE_GAIN;
    }
    return batch->copy_change != NULL ? batch->copy_change[idx]
                                      : (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN;
}

static inline void duckvep_event_load(
    const duckvep_variant_batch_t *batch, size_t idx, duckvep_event_t *event) {
    event->chrom_id = batch->chrom_id[idx];
    event->kind = batch->variant_kind != NULL ? batch->variant_kind[idx]
                                               : (uint8_t)DUCKVEP_KIND_SV;
    event->sv_type = duckvep_event_sv_type_at(batch, idx);
    event->copy_change = duckvep_event_copy_change_at(batch, idx);
    if (batch->variant_kind == NULL || event->kind == (uint8_t)DUCKVEP_KIND_SV) {
        duckvep_event_load_raw_interval(batch, idx, event);
    } else {
        duckvep_event_load_small_differing_region(batch, idx, event);
    }
}

#endif /* DUCKVEP_EVENT_H */
