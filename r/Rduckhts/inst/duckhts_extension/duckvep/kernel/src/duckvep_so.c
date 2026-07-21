/*
 * duckvep_so.c — SO term names + VEP class-model metadata. See duckvep_so.h.
 * No DuckDB/htslib/Arrow; no allocation.
 */
#include "duckvep_so.h"

#include <stddef.h>
#include <string.h>

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

static unsigned duckvep_first_set_bit(uint64_t mask) {
#if defined(__GNUC__) || defined(__clang__)
    return (unsigned)__builtin_ctzll(mask);
#else
    unsigned bit = 0u;
    while ((mask & UINT64_C(1)) == 0u) {
        mask >>= 1;
        bit++;
    }
    return bit;
#endif
}

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
    if ((mask & k_so_impact_high_mask) != 0u)
        return DUCKVEP_IMPACT_HIGH;
    if ((mask & k_so_impact_moderate_mask) != 0u)
        return DUCKVEP_IMPACT_MODERATE;
    if ((mask & k_so_impact_low_mask) != 0u)
        return DUCKVEP_IMPACT_LOW;
    return DUCKVEP_IMPACT_MODIFIER;
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

    if (mask == 0u) {
        if (buflen > 0u) buf[0] = '\0';
        return 0u;
    }
    if ((mask & (mask - UINT64_C(1))) == 0u) {
        unsigned bit = duckvep_first_set_bit(mask);
        const char *name = bit < (unsigned)DUCKVEP_SO_BIT_COUNT
                             ? k_so[bit].name : NULL;
        size_t len;
        size_t copy;

        if (name == NULL) {
            if (buflen > 0u) buf[0] = '\0';
            return 0u;
        }
        len = strlen(name);
        copy = buflen > 0u && len >= buflen ? buflen - 1u : len;
        if (buflen > 0u) {
            memcpy(buf, name, copy);
            buf[copy] = '\0';
        }
        return len;
    }

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
