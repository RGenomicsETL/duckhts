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

    for (; len - i >= 16; i += 16) {
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

/* BAM nt16 stores two codes per byte, with the even-indexed base in the high
 * nibble.  Each vector therefore accounts for 32 bases.  Codes other than
 * A/C/G/T/N are the IUPAC/'=' remainder, exactly as in the scalar oracle. */
static void duckhts_bam_nt16_counts_wasm_simd128(
    const uint8_t *packed_seq, int32_t n_bases,
    duckhts_simd_bam_nt16_counts_t *out) {
    const v128_t lo_mask = wasm_i8x16_splat(0x0F);
    const v128_t v_a = wasm_i8x16_splat(1);
    const v128_t v_c = wasm_i8x16_splat(2);
    const v128_t v_g = wasm_i8x16_splat(4);
    const v128_t v_t = wasm_i8x16_splat(8);
    const v128_t v_n = wasm_i8x16_splat(15);
    uint64_t a = 0;
    uint64_t c = 0;
    uint64_t g = 0;
    uint64_t t = 0;
    uint64_t n = 0;
    size_t i = 0;
    size_t total;

    out->a = out->c = out->g = out->t = 0;
    out->n = out->iupac = out->gc = out->called = 0;
    if (!packed_seq || n_bases <= 0) return;
    total = (size_t)n_bases;

    for (; total - i >= 32; i += 32) {
        v128_t bytes = wasm_v128_load((const void *)(packed_seq + (i >> 1)));
        v128_t hi = wasm_v128_and(wasm_u8x16_shr(bytes, 4), lo_mask);
        v128_t lo = wasm_v128_and(bytes, lo_mask);

#define DUCKHTS_COUNT_NT16(code, accumulator)                                      \
        do {                                                                        \
            uint32_t hi_mask = (uint32_t)wasm_i8x16_bitmask(                       \
                wasm_i8x16_eq(hi, (code)));                                         \
            uint32_t lo_bits = (uint32_t)wasm_i8x16_bitmask(                       \
                wasm_i8x16_eq(lo, (code)));                                         \
            (accumulator) += (uint64_t)DUCKHTS_POPCOUNT32(hi_mask);                 \
            (accumulator) += (uint64_t)DUCKHTS_POPCOUNT32(lo_bits);                 \
        } while (0)

        DUCKHTS_COUNT_NT16(v_a, a);
        DUCKHTS_COUNT_NT16(v_c, c);
        DUCKHTS_COUNT_NT16(v_g, g);
        DUCKHTS_COUNT_NT16(v_t, t);
        DUCKHTS_COUNT_NT16(v_n, n);
#undef DUCKHTS_COUNT_NT16
    }

    out->iupac = (uint64_t)i - (a + c + g + t + n);
    for (; i < total; i++) {
        uint8_t byte = packed_seq[i >> 1];
        unsigned int code = (i & 1u) != 0u ? (unsigned int)(byte & 0x0Fu)
                                           : (unsigned int)(byte >> 4);
        switch (code) {
        case 1: a++; break;
        case 2: c++; break;
        case 4: g++; break;
        case 8: t++; break;
        case 15: n++; break;
        default: out->iupac++; break;
        }
    }

    out->a = a;
    out->c = c;
    out->g = g;
    out->t = t;
    out->n = n;
    out->gc = c + g;
    out->called = a + c + g + t;
}

/* One unpacked nt16 code per byte.  Invalid input must leave gc/called at
 * zero, even when the bad code follows valid vector blocks; this is the scalar
 * contract consumed by seq_gc_content's NULL path. */
static void duckhts_nt16_gc_counts_wasm_simd128(
    const uint8_t *codes, size_t n, duckhts_simd_base_counts_t *out) {
    const v128_t v_a = wasm_i8x16_splat(1);
    const v128_t v_c = wasm_i8x16_splat(2);
    const v128_t v_g = wasm_i8x16_splat(4);
    const v128_t v_t = wasm_i8x16_splat(8);
    const v128_t v_n = wasm_i8x16_splat(15);
    uint64_t gc = 0;
    uint64_t called = 0;
    size_t i = 0;

    out->gc = 0;
    out->called = 0;
    out->invalid = 0;
    if (!codes) return;

    for (; n - i >= 16; i += 16) {
        v128_t x = wasm_v128_load((const void *)(codes + i));
        uint32_t ma = (uint32_t)wasm_i8x16_bitmask(wasm_i8x16_eq(x, v_a));
        uint32_t mc = (uint32_t)wasm_i8x16_bitmask(wasm_i8x16_eq(x, v_c));
        uint32_t mg = (uint32_t)wasm_i8x16_bitmask(wasm_i8x16_eq(x, v_g));
        uint32_t mt = (uint32_t)wasm_i8x16_bitmask(wasm_i8x16_eq(x, v_t));
        uint32_t mn = (uint32_t)wasm_i8x16_bitmask(wasm_i8x16_eq(x, v_n));

        if ((ma | mc | mg | mt | mn) != 0xFFFFu) {
            out->invalid = 1;
            return;
        }
        gc += (uint64_t)DUCKHTS_POPCOUNT32(mc | mg);
        called += (uint64_t)DUCKHTS_POPCOUNT32(ma | mc | mg | mt);
    }

    for (; i < n; i++) {
        uint8_t code = codes[i];
        switch (code) {
        case 2:
        case 4:
            gc++;
            called++;
            break;
        case 1:
        case 8:
            called++;
            break;
        case 15:
            break;
        default:
            out->invalid = 1;
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
    duckhts_simd_builder_consider_bam_nt16_counts(
        builder, DUCKHTS_SIMD_CAP_WASM_SIMD128, "wasm_simd128", 60,
        duckhts_bam_nt16_counts_wasm_simd128);
    duckhts_simd_builder_consider_nt16_gc_counts(
        builder, DUCKHTS_SIMD_CAP_WASM_SIMD128, "wasm_simd128", 60,
        duckhts_nt16_gc_counts_wasm_simd128);
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
