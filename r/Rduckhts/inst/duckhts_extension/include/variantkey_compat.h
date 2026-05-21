#ifndef DUCKHTS_VARIANTKEY_COMPAT_H
#define DUCKHTS_VARIANTKEY_COMPAT_H

#include "variantkey/variantkey.h"
#include "variantkey/regionkey.h"

#define DUCKHTS_VARIANTKEY_POS0_MAX 0x0FFFFFFFU
#define DUCKHTS_VARIANTKEY_POS1_MAX ((uint64_t)DUCKHTS_VARIANTKEY_POS0_MAX + 1ULL)

static inline int duckhts_variantkey_pos1_to_pos0(int64_t pos1, uint32_t *out) {
    if (!out) return 0;
    if (pos1 < 1 || (uint64_t)pos1 > DUCKHTS_VARIANTKEY_POS1_MAX) return 0;
    *out = (uint32_t)(pos1 - 1);
    return 1;
}

static inline int duckhts_regionkey_pos0_valid(int64_t pos0) {
    return (pos0 >= 0 && (uint64_t)pos0 <= (uint64_t)RK_MAX_POS);
}

static inline int duckhts_variantkey_hex_is_valid(const char *s, size_t len) {
    size_t i;
    if (!s || len != 16) return 0;
    for (i = 0; i < 16; i++) {
        unsigned char c = (unsigned char)s[i];
        if (!((c >= '0' && c <= '9') ||
              (c >= 'a' && c <= 'f') ||
              (c >= 'A' && c <= 'F'))) {
            return 0;
        }
    }
    return 1;
}

#endif
