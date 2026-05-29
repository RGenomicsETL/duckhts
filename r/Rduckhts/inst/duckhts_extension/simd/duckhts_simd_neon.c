#include <stdint.h>

#if defined(__aarch64__) && !defined(DUCKDB_WASM_EXTENSION)
#include <arm_neon.h>
#endif

#include "duckhts_simd_internal.h"

#if defined(__aarch64__) && !defined(DUCKDB_WASM_EXTENSION)
#define DUCKHTS_SIMD_COMPILED_NEON 1

static inline uint64_t duckhts_neon_sum_u8(uint8x16_t x) {
    return (uint64_t)vaddvq_u8(x);
}

static void duckhts_base_counts_neon(const char *seq, size_t len,
                                     duckhts_simd_base_counts_t *out) {
    const uint8x16_t case_mask = vdupq_n_u8(0xDFu);
    const uint8x16_t v_a = vdupq_n_u8((uint8_t)'A');
    const uint8x16_t v_c = vdupq_n_u8((uint8_t)'C');
    const uint8x16_t v_g = vdupq_n_u8((uint8_t)'G');
    const uint8x16_t v_t = vdupq_n_u8((uint8_t)'T');
    const uint8x16_t v_n = vdupq_n_u8((uint8_t)'N');
    uint64_t gc = 0;
    uint64_t called = 0;
    size_t i = 0;

    out->gc = 0;
    out->called = 0;
    out->invalid = 0;

    for (; i + 16 <= len; i += 16) {
        uint8x16_t bytes = vld1q_u8((const uint8_t *)(const void *)(seq + i));
        uint8x16_t up = vandq_u8(bytes, case_mask);
        uint8x16_t a_mask = vceqq_u8(up, v_a);
        uint8x16_t c_mask = vceqq_u8(up, v_c);
        uint8x16_t g_mask = vceqq_u8(up, v_g);
        uint8x16_t t_mask = vceqq_u8(up, v_t);
        uint8x16_t n_mask = vceqq_u8(up, v_n);
        uint8x16_t gc_mask = vorrq_u8(c_mask, g_mask);
        uint8x16_t called_mask = vorrq_u8(vorrq_u8(a_mask, c_mask), vorrq_u8(g_mask, t_mask));
        uint8x16_t valid_mask = vorrq_u8(called_mask, n_mask);

        if (duckhts_neon_sum_u8(vshrq_n_u8(valid_mask, 7)) != 16u) {
            out->invalid = 1;
            out->gc = gc;
            out->called = called;
            return;
        }
        gc += duckhts_neon_sum_u8(vshrq_n_u8(gc_mask, 7));
        called += duckhts_neon_sum_u8(vshrq_n_u8(called_mask, 7));
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

void duckhts_simd_neon_register(duckhts_simd_builder_t *builder) {
    duckhts_simd_builder_consider_base_counts(builder,
                                              DUCKHTS_SIMD_CAP_NEON,
                                              "neon",
                                              50,
                                              duckhts_base_counts_neon);
}

int duckhts_simd_neon_compiled(void) { return 1; }
int duckhts_simd_neon_cpu_supported(void) { return 1; }
int duckhts_simd_neon_available(void) { return 1; }

#else
void duckhts_simd_neon_register(duckhts_simd_builder_t *builder) { (void)builder; }
int duckhts_simd_neon_compiled(void) { return 0; }
int duckhts_simd_neon_cpu_supported(void) { return 0; }
int duckhts_simd_neon_available(void) { return 0; }
#endif
