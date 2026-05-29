#include <stdint.h>

#include "duckhts_simd_internal.h"

static void duckhts_base_counts_scalar(const char *seq, size_t len,
                                       duckhts_simd_base_counts_t *out) {
    uint64_t gc = 0;
    uint64_t called = 0;

    out->gc = 0;
    out->called = 0;
    out->invalid = 0;

    for (size_t i = 0; i < len; i++) {
        unsigned char c = (unsigned char)seq[i] & 0xDFu;
        switch (c) {
        case 'G':
        case 'C':
            gc++;
            called++;
            break;
        case 'A':
        case 'T':
            called++;
            break;
        case 'N':
            break;
        default:
            out->invalid = 1;
            out->gc = gc;
            out->called = called;
            return;
        }
    }

    out->gc = gc;
    out->called = called;
}

void duckhts_simd_scalar_register(duckhts_simd_builder_t *builder) {
    duckhts_simd_builder_consider_base_counts(builder,
                                              DUCKHTS_SIMD_CAP_SCALAR,
                                              "scalar",
                                              1000,
                                              duckhts_base_counts_scalar);
}
