#include <config.h>

#include <stdint.h>

#if defined(__x86_64__) && defined(HAVE_X86INTRIN_H) && HAVE_X86INTRIN_H
#include <x86intrin.h>
#endif

#include "duckhts_simd_internal.h"

#if defined(__x86_64__) && defined(HAVE_X86INTRIN_H) && HAVE_X86INTRIN_H && \
    defined(HAVE_ATTRIBUTE_TARGET_SSSE3) && HAVE_ATTRIBUTE_TARGET_SSSE3 && \
    defined(HAVE_BUILTIN_CPU_SUPPORT_SSSE3) && HAVE_BUILTIN_CPU_SUPPORT_SSSE3 && \
    defined(HAVE_AVX512) && HAVE_AVX512 && defined(HAVE_POPCNT) && HAVE_POPCNT && \
    !defined(DUCKDB_WASM_EXTENSION)
#define DUCKHTS_SIMD_COMPILED_AVX512 1
#if defined(__GNUC__) || defined(__clang__)
#define DUCKHTS_POPCOUNT64(x) __builtin_popcountll((unsigned long long)(x))
#else
static unsigned int DUCKHTS_POPCOUNT64(uint64_t x) {
    unsigned int n = 0;
    while (x) {
        n += (unsigned int)(x & UINT64_C(1));
        x >>= 1;
    }
    return n;
}
#endif

/* Compile only this function for AVX512F+AVX512BW+POPCNT.  The rest of the
 * translation unit stays on baseline CFLAGS so runtime availability probes can
 * execute safely on machines without AVX-512. */
__attribute__((target("avx512f,avx512bw,popcnt")))
static void duckhts_base_counts_avx512(const char *seq, size_t len,
                                       duckhts_simd_base_counts_t *out) {
    const __m512i case_mask = _mm512_set1_epi8((char)0xDF);
    const __m512i v_a = _mm512_set1_epi8('A');
    const __m512i v_c = _mm512_set1_epi8('C');
    const __m512i v_g = _mm512_set1_epi8('G');
    const __m512i v_t = _mm512_set1_epi8('T');
    const __m512i v_n = _mm512_set1_epi8('N');
    uint64_t gc = 0;
    uint64_t called = 0;
    size_t i = 0;

    out->gc = 0;
    out->called = 0;
    out->invalid = 0;

    for (; i + 64 <= len; i += 64) {
        __m512i bytes = _mm512_loadu_si512((const void *)(seq + i));
        __m512i up = _mm512_and_si512(bytes, case_mask);
        uint64_t a_mask = (uint64_t)_mm512_cmpeq_epi8_mask(up, v_a);
        uint64_t c_mask = (uint64_t)_mm512_cmpeq_epi8_mask(up, v_c);
        uint64_t g_mask = (uint64_t)_mm512_cmpeq_epi8_mask(up, v_g);
        uint64_t t_mask = (uint64_t)_mm512_cmpeq_epi8_mask(up, v_t);
        uint64_t n_mask = (uint64_t)_mm512_cmpeq_epi8_mask(up, v_n);
        uint64_t gc_mask = c_mask | g_mask;
        uint64_t called_mask = a_mask | c_mask | g_mask | t_mask;
        uint64_t valid_mask = called_mask | n_mask;

        if (valid_mask != UINT64_MAX) {
            out->invalid = 1;
            out->gc = gc;
            out->called = called;
            return;
        }
        gc += (uint64_t)DUCKHTS_POPCOUNT64(gc_mask);
        called += (uint64_t)DUCKHTS_POPCOUNT64(called_mask);
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

static int duckhts_cpu_has_avx512(void) {
#if defined(__GNUC__) || defined(__clang__)
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx512f") != 0 &&
           __builtin_cpu_supports("avx512bw") != 0 &&
           __builtin_cpu_supports("popcnt") != 0;
#else
    return 0;
#endif
}

static const duckhts_simd_ops_t duckhts_simd_ops_avx512 = {
    "avx512",
    duckhts_base_counts_avx512
};

int duckhts_simd_avx512_compiled(void) { return 1; }

int duckhts_simd_avx512_cpu_supported(void) {
    return duckhts_cpu_has_avx512();
}

int duckhts_simd_avx512_available(void) {
    return duckhts_simd_avx512_compiled() && duckhts_simd_avx512_cpu_supported();
}

const duckhts_simd_ops_t *duckhts_simd_avx512_ops_if_available(void) {
    return duckhts_simd_avx512_available() ? &duckhts_simd_ops_avx512 : (const duckhts_simd_ops_t *)0;
}
#else
int duckhts_simd_avx512_compiled(void) { return 0; }
int duckhts_simd_avx512_cpu_supported(void) { return 0; }
int duckhts_simd_avx512_available(void) { return 0; }

const duckhts_simd_ops_t *duckhts_simd_avx512_ops_if_available(void) {
    return (const duckhts_simd_ops_t *)0;
}
#endif
