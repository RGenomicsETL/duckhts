#include <stdint.h>

#if defined(__wasm_simd128__)
#include <wasm_simd128.h>
#endif

#include "duckhts_simd_internal.h"

#if defined(__wasm_simd128__)
#define DUCKHTS_SIMD_COMPILED_WASM_SIMD128 1
#if defined(__GNUC__) || defined(__clang__)
#define DUCKHTS_POPCOUNT32(x) __builtin_popcount((unsigned int)(x))
#else
static unsigned int DUCKHTS_POPCOUNT32(uint32_t x) {
    unsigned int n = 0;
    while (x) {
        n += (unsigned int)(x & 1u);
        x >>= 1;
    }
    return n;
}
#endif

static void duckhts_base_counts_wasm_simd128(const char *seq, size_t len,
                                             duckhts_simd_base_counts_t *out) {
    const v128_t case_mask = wasm_i8x16_splat((char)0xDF);
    const v128_t v_a = wasm_i8x16_splat('A');
    const v128_t v_c = wasm_i8x16_splat('C');
    const v128_t v_g = wasm_i8x16_splat('G');
    const v128_t v_t = wasm_i8x16_splat('T');
    const v128_t v_n = wasm_i8x16_splat('N');
    uint64_t gc = 0;
    uint64_t called = 0;
    size_t i = 0;

    out->gc = 0;
    out->called = 0;
    out->invalid = 0;

    for (; i + 16 <= len; i += 16) {
        v128_t bytes = wasm_v128_load((const void *)(seq + i));
        v128_t up = wasm_v128_and(bytes, case_mask);
        uint32_t a_mask = (uint32_t)wasm_i8x16_bitmask(wasm_i8x16_eq(up, v_a));
        uint32_t c_mask = (uint32_t)wasm_i8x16_bitmask(wasm_i8x16_eq(up, v_c));
        uint32_t g_mask = (uint32_t)wasm_i8x16_bitmask(wasm_i8x16_eq(up, v_g));
        uint32_t t_mask = (uint32_t)wasm_i8x16_bitmask(wasm_i8x16_eq(up, v_t));
        uint32_t n_mask = (uint32_t)wasm_i8x16_bitmask(wasm_i8x16_eq(up, v_n));
        uint32_t gc_mask = c_mask | g_mask;
        uint32_t called_mask = a_mask | c_mask | g_mask | t_mask;
        uint32_t valid_mask = called_mask | n_mask;

        if (valid_mask != 0xFFFFu) {
            out->invalid = 1;
            out->gc = gc;
            out->called = called;
            return;
        }
        gc += (uint64_t)DUCKHTS_POPCOUNT32(gc_mask);
        called += (uint64_t)DUCKHTS_POPCOUNT32(called_mask);
    }

    for (; i < len; i++) {
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

void duckhts_simd_wasm_simd128_register(duckhts_simd_builder_t *builder) {
    duckhts_simd_builder_consider_base_counts(builder,
                                              DUCKHTS_SIMD_CAP_WASM_SIMD128,
                                              "wasm_simd128",
                                              60,
                                              duckhts_base_counts_wasm_simd128);
}

int duckhts_simd_wasm_simd128_compiled(void) { return 1; }
int duckhts_simd_wasm_simd128_cpu_supported(void) { return 1; }
int duckhts_simd_wasm_simd128_available(void) { return 1; }

#else
void duckhts_simd_wasm_simd128_register(duckhts_simd_builder_t *builder) { (void)builder; }
int duckhts_simd_wasm_simd128_compiled(void) { return 0; }
int duckhts_simd_wasm_simd128_cpu_supported(void) { return 0; }
int duckhts_simd_wasm_simd128_available(void) { return 0; }
#endif
