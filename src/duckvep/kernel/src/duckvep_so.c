/*
 * duckvep_so.c — SO term names + VEP class-model metadata. See duckvep_so.h.
 * No DuckDB/htslib/Arrow; no allocation.
 */
#include "duckvep_so.h"

#include <stddef.h>

/* Stable bit assignments live in duckvep_so.h; VEP metadata and the rank-sorted
 * render order are generated from the pinned class model plus an explicit
 * term->bit binding ledger. */
struct duckvep_so_metadata {
    const char       *name;
    duckvep_impact_t  impact;
    uint8_t           rank;
    uint8_t           tier;
};

#include "duckvep_so_metadata.inc"

typedef char duckvep_so_binding_count_must_match_public_enum[
    (int)DUCKVEP_GENERATED_SO_BINDING_COUNT == (int)DUCKVEP_SO_BIT_COUNT ? 1 : -1];

const char *duckvep_so_name(duckvep_so_bit_t bit) {
    if ((int)bit < 0 || bit >= DUCKVEP_SO_BIT_COUNT) return NULL;
    return k_so[bit].name;
}

uint8_t duckvep_so_rank(duckvep_so_bit_t bit) {
    if ((int)bit < 0 || bit >= DUCKVEP_SO_BIT_COUNT) return UINT8_MAX;
    return k_so[bit].rank;
}

uint8_t duckvep_so_tier(duckvep_so_bit_t bit) {
    if ((int)bit < 0 || bit >= DUCKVEP_SO_BIT_COUNT) return UINT8_MAX;
    return k_so[bit].tier;
}

duckvep_impact_t duckvep_so_bit_impact(duckvep_so_bit_t bit) {
    if ((int)bit < 0 || bit >= DUCKVEP_SO_BIT_COUNT) return DUCKVEP_IMPACT_MODIFIER;
    return k_so[bit].impact;
}

duckvep_impact_t duckvep_so_impact(uint64_t mask) {
    duckvep_impact_t best = DUCKVEP_IMPACT_MODIFIER;
    int bit;
    for (bit = 0; bit < DUCKVEP_SO_BIT_COUNT; bit++) {
        if ((mask & DUCKVEP_SO((duckvep_so_bit_t)bit)) != 0u) {
            duckvep_impact_t imp = k_so[bit].impact;
            if (imp > best) best = imp;
        }
    }
    return best;
}

const char *duckvep_impact_name(duckvep_impact_t impact) {
    switch (impact) {
    case DUCKVEP_IMPACT_HIGH:     return "HIGH";
    case DUCKVEP_IMPACT_MODERATE: return "MODERATE";
    case DUCKVEP_IMPACT_LOW:      return "LOW";
    case DUCKVEP_IMPACT_MODIFIER:
    default:                      return "MODIFIER";
    }
}

size_t duckvep_so_render(uint64_t mask, char sep, char *buf, size_t buflen) {
    size_t need = 0u;
    size_t have = 0u;
    size_t oi;
    int first = 1;

    for (oi = 0u; oi < DUCKVEP_SO_BIT_COUNT; oi++) {
        duckvep_so_bit_t bit = k_so_render_order[oi];
        const char *name;
        size_t i;
        if ((mask & DUCKVEP_SO(bit)) == 0u) continue;
        name = k_so[bit].name;
        if (name == NULL) continue;

        if (!first) {
            if (have + 1u < buflen) buf[have++] = sep;
            need++;
        }
        first = 0;
        for (i = 0u; name[i] != '\0'; i++) {
            if (have + 1u < buflen) buf[have++] = name[i];
            need++;
        }
    }
    if (buflen > 0u) buf[have] = '\0';
    return need;
}
