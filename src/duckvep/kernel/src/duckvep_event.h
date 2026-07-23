/*
 * duckvep_event.h — execution event geometry and structural semantics.
 * INTERNAL, header-only: the public batch remains a borrowed SoA, while the hot
 * path loads one compact event value per ALT.
 *
 * Single-interval SV/CNV rows use their supplied inclusive [pos1,end1] span.
 * Small variants carry three geometries:
 *
 *   raw_*      the uploaded VCF REF span;
 *   feature_*  VEP's VariationFeature span. VEP 116 fully minimizes every
 *              length-changing biallelic allele during parser post-processing,
 *              even without --minimal; equal-length substitutions retain their
 *              uploaded span;
 *   start/end  the fully trimmed semantic edit used for CDS projection.
 *
 * Region predicates consume feature_*. Splice predicates inspect the separate
 * mismatch islands inside the feature alleles. Sequence projection consumes the
 * semantic edit. Keeping these explicit prevents one convenient interval from
 * becoming three subtly different biological authorities.
 */
#ifndef DUCKVEP_EVENT_H
#define DUCKVEP_EVENT_H

#include "duckvep_kernel.h"

#include <stddef.h>
#include <stdint.h>

/* VEP's ordinary-allele pre-predicate is computed from feature REF/ALT
 * lengths before transcript projection.  This can differ from a mapped CDS
 * edit when the uploaded feature spans an intron. */
typedef enum duckvep_feature_length_relation {
    DUCKVEP_FEATURE_LENGTH_UNKNOWN = 0,
    DUCKVEP_FEATURE_LENGTH_EQUAL,
    DUCKVEP_FEATURE_LENGTH_INCREASE,
    DUCKVEP_FEATURE_LENGTH_DECREASE
} duckvep_feature_length_relation_t;

typedef struct duckvep_event {
    uint16_t chrom_id;
    uint16_t mate_chrom_id;
    uint32_t raw_start1;
    uint32_t raw_end1;
    uint32_t mate_pos1;
    uint32_t feature_start1;  /* VEP VariationFeature start; may exceed end */
    uint32_t feature_end1;
    uint32_t start1;          /* fully trimmed semantic edit coordinate */
    uint32_t end1;            /* inclusive; pure insertions use a point */
    uint32_t insertion_boundary0; /* genomic interbase boundary; 0 = before base 1 */
    uint16_t ref_diff_offset; /* offset inside REF where differing region starts */
    uint16_t alt_diff_offset; /* offset inside ALT where differing region starts */
    uint16_t anchor_ref_offset; /* REF base used to validate a pure insertion */
    uint16_t feature_allele_offset; /* VEP feature allele offset inside REF/ALT */
    uint16_t ref_diff_length;
    uint16_t alt_diff_length;
    uint8_t  interbase;       /* pure insertion after trimming */
    uint8_t  anchor_side;     /* duckvep_event_anchor_t */
    uint8_t  kind;
    uint8_t  sv_type;
    uint8_t  copy_change;
    uint8_t  has_mate;
    uint8_t  feature_length_relation; /* duckvep_feature_length_relation_t */
} duckvep_event_t;

typedef enum duckvep_event_anchor {
    DUCKVEP_EVENT_ANCHOR_NONE = 0,
    DUCKVEP_EVENT_ANCHOR_LEFT,
    DUCKVEP_EVENT_ANCHOR_RIGHT
} duckvep_event_anchor_t;

static inline uint8_t duckvep_event_ascii_upper(uint8_t byte) {
    return byte >= (uint8_t)'a' && byte <= (uint8_t)'z'
        ? (uint8_t)(byte - ((uint8_t)'a' - (uint8_t)'A')) : byte;
}

static inline uint32_t duckvep_event_sat_add_u32_u16(uint32_t x, uint16_t y) {
    return x > UINT32_MAX - (uint32_t)y ? UINT32_MAX : x + (uint32_t)y;
}

/* VEP places a pure insertion between boundary P and P+1. Projection validates
 * one retained anchor base; predicates that cross a right-hand transcript,
 * exon, or CDS boundary also need P+1 explicitly. */
static inline uint32_t duckvep_event_right_flank1(
    const duckvep_event_t *event) {

    if (event == NULL || !event->interbase) return event != NULL ? event->start1 : 0u;
    return event->insertion_boundary0 == UINT32_MAX
        ? UINT32_MAX
        : event->insertion_boundary0 + 1u;
}

static inline uint32_t duckvep_event_feature_max1(
    const duckvep_event_t *event) {

    if (event == NULL) return 0u;
    return event->feature_start1 > event->feature_end1
        ? event->feature_start1
        : event->feature_end1;
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

/* Borrow VEP's feature alleles. Length-changing pairs use the fully minimized
 * differing slices; equal-length substitutions keep the complete uploaded
 * alleles. Region, splice, and sequence predicates must all use this helper so
 * feature geometry has one authority. */
static inline int duckvep_event_feature_alleles(
    const duckvep_variant_batch_t *batch,
    size_t                         idx,
    const duckvep_event_t         *event,
    const uint8_t                **ref,
    uint16_t                      *ref_length,
    const uint8_t                **alt,
    uint16_t                      *alt_length) {

    uint16_t offset;
    uint16_t raw_ref_length;
    uint16_t raw_alt_length;

    if (batch == NULL || event == NULL || ref == NULL || ref_length == NULL ||
        alt == NULL || alt_length == NULL || idx >= batch->count ||
        !duckvep_event_allele_slices_ok(batch, idx)) {
        return 0;
    }
    offset = event->feature_allele_offset;
    raw_ref_length = batch->ref_length[idx];
    raw_alt_length = batch->alt_length[idx];
    if (offset > raw_ref_length || offset > raw_alt_length) return 0;

    *ref = batch->allele_bytes + batch->ref_offset[idx] + offset;
    *alt = batch->allele_bytes + batch->alt_offset[idx] + offset;
    if (raw_ref_length != raw_alt_length) {
        *ref_length = event->ref_diff_length;
        *alt_length = event->alt_diff_length;
    } else {
        *ref_length = raw_ref_length;
        *alt_length = raw_alt_length;
    }
    return 1;
}

static inline void duckvep_event_load_raw_interval(
    const duckvep_variant_batch_t *batch,
    size_t                         idx,
    duckvep_event_t               *event) {

    event->raw_start1 = batch->pos1[idx];
    event->raw_end1 = batch->end1[idx];
    event->feature_start1 = batch->pos1[idx];
    event->feature_end1 = batch->end1[idx];
    event->start1 = batch->pos1[idx];
    event->end1 = batch->end1[idx];
    event->insertion_boundary0 = 0u;
    event->ref_diff_offset = 0u;
    event->alt_diff_offset = 0u;
    event->anchor_ref_offset = 0u;
    event->feature_allele_offset = 0u;
    event->ref_diff_length = 0u;
    event->alt_diff_length = 0u;
    event->interbase = 0u;
    event->anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_NONE;
    event->feature_length_relation =
        (uint8_t)DUCKVEP_FEATURE_LENGTH_UNKNOWN;
}

/* Decode one ordinary REF/ALT pair into its lossless effect geometry. This is
 * representation interpretation, not left alignment or canonical rewriting.
 * In particular, VCF permits a right padding base at contig position 1. The
 * explicit boundary and anchor side preserve that case without inventing a
 * genomic base zero. */
static inline int duckvep_event_prepare_small(
    uint32_t         pos1,
    const uint8_t   *ref,
    uint16_t         ref_len,
    const uint8_t   *alt,
    uint16_t         alt_len,
    duckvep_event_t *event) {

    uint16_t prefix = 0u;
    uint16_t suffix = 0u;
    uint16_t ref_rem;
    uint16_t alt_rem;
    uint16_t ref_diff_len;
    uint16_t alt_diff_len;
    uint32_t diff_start;

    if (event == NULL || ref == NULL || alt == NULL || pos1 == 0u ||
        ref_len == 0u || alt_len == 0u ||
        (uint32_t)(ref_len - 1u) > UINT32_MAX - pos1) {
        return 0;
    }

    event->raw_start1 = pos1;
    event->raw_end1 = pos1 + (uint32_t)ref_len - 1u;
    event->feature_start1 = pos1;
    event->feature_end1 = event->raw_end1;
    event->start1 = pos1;
    event->end1 = pos1;
    event->insertion_boundary0 = 0u;
    event->ref_diff_offset = 0u;
    event->alt_diff_offset = 0u;
    event->anchor_ref_offset = 0u;
    event->feature_allele_offset = 0u;
    event->ref_diff_length = 0u;
    event->alt_diff_length = 0u;
    event->interbase = 0u;
    event->anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_NONE;
    event->feature_length_relation =
        ref_len == alt_len
        ? (uint8_t)DUCKVEP_FEATURE_LENGTH_EQUAL
        : alt_len > ref_len
            ? (uint8_t)DUCKVEP_FEATURE_LENGTH_INCREASE
            : (uint8_t)DUCKVEP_FEATURE_LENGTH_DECREASE;

    while (prefix < ref_len && prefix < alt_len &&
           duckvep_event_ascii_upper(ref[prefix]) ==
               duckvep_event_ascii_upper(alt[prefix])) {
        prefix++;
    }
    ref_rem = (uint16_t)(ref_len - prefix);
    alt_rem = (uint16_t)(alt_len - prefix);
    while (suffix < ref_rem && suffix < alt_rem &&
           duckvep_event_ascii_upper(
               ref[(uint16_t)(ref_len - 1u - suffix)]) ==
           duckvep_event_ascii_upper(
               alt[(uint16_t)(alt_len - 1u - suffix)])) {
        suffix++;
    }

    ref_diff_len = (uint16_t)(ref_len - prefix - suffix);
    alt_diff_len = (uint16_t)(alt_len - prefix - suffix);
    if (ref_diff_len == 0u && alt_diff_len == 0u) return 0;

    event->ref_diff_offset = prefix;
    event->alt_diff_offset = prefix;
    event->ref_diff_length = ref_diff_len;
    event->alt_diff_length = alt_diff_len;
    if (ref_diff_len == 0u) {
        event->kind = (uint8_t)DUCKVEP_KIND_INS;
    } else if (alt_diff_len == 0u) {
        event->kind = (uint8_t)DUCKVEP_KIND_DEL;
    } else if (ref_diff_len == alt_diff_len) {
        event->kind = ref_diff_len == 1u
            ? (uint8_t)DUCKVEP_KIND_SNV
            : (uint8_t)DUCKVEP_KIND_MNV;
    } else {
        event->kind = (uint8_t)DUCKVEP_KIND_INDEL;
    }

    diff_start = duckvep_event_sat_add_u32_u16(pos1, prefix);
    if (ref_len != alt_len) {
        /* Parser::post_process_vfs calls minimise_alleles for every remaining
         * non-minimal indel. trim_sequences adjusts both ends and represents an
         * insertion with start=end+1. */
        event->feature_allele_offset = prefix;
        event->feature_start1 = diff_start;
        event->feature_end1 = event->raw_end1 - (uint32_t)suffix;
    }
    event->interbase = (uint8_t)(ref_diff_len == 0u);
    if (event->interbase) {
        event->insertion_boundary0 = (pos1 - 1u) + (uint32_t)prefix;
        if (prefix > 0u) {
            event->anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT;
            event->anchor_ref_offset = (uint16_t)(prefix - 1u);
            event->start1 = event->insertion_boundary0;
        } else {
            /* With no shared prefix, the retained VCF padding is on the right.
             * The first REF base is the validation anchor immediately after the
             * insertion boundary. */
            event->anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_RIGHT;
            event->anchor_ref_offset = 0u;
            event->start1 = event->insertion_boundary0 + 1u;
        }
        event->end1 = event->start1;
    } else {
        event->start1 = diff_start;
        event->end1 = duckvep_event_sat_add_u32_u16(
            diff_start, (uint16_t)(ref_diff_len - 1u));
    }
    return 1;
}

static inline void duckvep_event_load_small_differing_region(
    const duckvep_variant_batch_t *batch,
    size_t                         idx,
    duckvep_event_t               *event) {

    duckvep_event_load_raw_interval(batch, idx, event);
    event->kind = batch->variant_kind != NULL
        ? batch->variant_kind[idx]
        : (uint8_t)DUCKVEP_KIND_SV;
    event->end1 = event->start1; /* small-variant fallback when alleles are absent */
    event->feature_end1 = event->feature_start1;
    if (!duckvep_event_allele_slices_ok(batch, idx)) return;
    (void)duckvep_event_prepare_small(
        batch->pos1[idx], batch->allele_bytes + batch->ref_offset[idx],
        batch->ref_length[idx], batch->allele_bytes + batch->alt_offset[idx],
        batch->alt_length[idx], event);
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
    return duckvep_event_feature_max1(&event);
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
        sv_type == (uint8_t)DUCKVEP_SV_TANDEM_DUPLICATION ||
        sv_type == (uint8_t)DUCKVEP_SV_TANDEM_REPEAT) {
        return (uint8_t)DUCKVEP_COPY_CHANGE_GAIN;
    }
    return batch->copy_change != NULL ? batch->copy_change[idx]
                                      : (uint8_t)DUCKVEP_COPY_CHANGE_UNKNOWN;
}

static inline void duckvep_event_load(
    const duckvep_variant_batch_t *batch, size_t idx, duckvep_event_t *event) {
    event->chrom_id = batch->chrom_id[idx];
    event->mate_chrom_id = 0u;
    event->mate_pos1 = 0u;
    event->has_mate = 0u;
    event->kind = batch->variant_kind != NULL ? batch->variant_kind[idx]
                                               : (uint8_t)DUCKVEP_KIND_SV;
    event->sv_type = duckvep_event_sv_type_at(batch, idx);
    event->copy_change = duckvep_event_copy_change_at(batch, idx);
    if (batch->variant_kind == NULL || event->kind == (uint8_t)DUCKVEP_KIND_SV) {
        duckvep_event_load_raw_interval(batch, idx, event);
        if (event->kind == (uint8_t)DUCKVEP_KIND_SV &&
            event->sv_type == (uint8_t)DUCKVEP_SV_INSERTION) {
            /* VEP represents a symbolic insertion at boundary P as the
             * reversed feature interval P+1,P after removing the VCF anchor.
             * The public structural adapter stores P as a valid point interval
             * for sorting; reconstruct the interbase geometry once here. */
            event->feature_start1 = event->raw_start1 == UINT32_MAX
                ? UINT32_MAX : event->raw_start1 + 1u;
            event->feature_end1 = event->raw_start1;
            event->start1 = event->raw_start1;
            event->end1 = event->raw_start1;
            event->insertion_boundary0 = event->raw_start1;
            event->interbase = 1u;
            event->anchor_side = (uint8_t)DUCKVEP_EVENT_ANCHOR_LEFT;
        } else if (event->kind == (uint8_t)DUCKVEP_KIND_SV &&
                   event->sv_type == (uint8_t)DUCKVEP_SV_BREAKEND) {
            /* BaseVCF4::get_start removes the local VCF anchor by advancing
             * POS once. StructuralVariationFeature::_parse_breakends retains
             * the mate coordinate verbatim. This asymmetry is observable in
             * VEP 116 consequence predicates. */
            event->feature_start1 = event->raw_start1 == UINT32_MAX
                ? UINT32_MAX : event->raw_start1 + 1u;
            event->feature_end1 = event->feature_start1;
            event->start1 = event->feature_start1;
            event->end1 = event->feature_start1;
            if (batch->mate_chrom_id != NULL && batch->mate_pos1 != NULL) {
                event->mate_chrom_id = batch->mate_chrom_id[idx];
                event->mate_pos1 = batch->mate_pos1[idx];
                event->has_mate = 1u;
            }
        }
    } else {
        duckvep_event_load_small_differing_region(batch, idx, event);
    }
}

#endif /* DUCKVEP_EVENT_H */
