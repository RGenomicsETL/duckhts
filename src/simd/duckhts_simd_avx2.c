#include <stdint.h>

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && \
    !defined(DUCKDB_WASM_EXTENSION)
#include <x86intrin.h>
#endif

#include "duckhts_simd_internal.h"

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && \
    !defined(DUCKDB_WASM_EXTENSION)
#define DUCKHTS_SIMD_COMPILED_AVX2 1
#define DUCKHTS_POPCOUNT32(x) __builtin_popcount((unsigned int)(x))

/* Compile only this function for AVX2+POPCNT.  The translation unit itself
 * stays on the baseline CFLAGS so availability/probe code remains safe on
 * non-AVX2 hosts. */
__attribute__((target("avx2,popcnt")))
static void duckhts_base_counts_avx2(const char *seq, size_t len,
                                     duckhts_simd_base_counts_t *out) {
    const __m256i case_mask = _mm256_set1_epi8((char)0xDF);
    const __m256i v_a = _mm256_set1_epi8('A');
    const __m256i v_c = _mm256_set1_epi8('C');
    const __m256i v_g = _mm256_set1_epi8('G');
    const __m256i v_t = _mm256_set1_epi8('T');
    const __m256i v_n = _mm256_set1_epi8('N');
    uint64_t gc = 0;
    uint64_t called = 0;
    size_t i = 0;

    out->gc = 0;
    out->called = 0;
    out->invalid = 0;

    for (; i + 32 <= len; i += 32) {
        __m256i bytes = _mm256_loadu_si256((const __m256i *)(const void *)(seq + i));
        __m256i up = _mm256_and_si256(bytes, case_mask);
        uint32_t a_mask = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(up, v_a));
        uint32_t c_mask = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(up, v_c));
        uint32_t g_mask = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(up, v_g));
        uint32_t t_mask = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(up, v_t));
        uint32_t n_mask = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(up, v_n));
        uint32_t gc_mask = c_mask | g_mask;
        uint32_t called_mask = a_mask | c_mask | g_mask | t_mask;
        uint32_t valid_mask = called_mask | n_mask;

        if (valid_mask != UINT32_MAX) {
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

static int duckhts_cpu_has_avx2(void) {
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx2") != 0 && __builtin_cpu_supports("popcnt") != 0;
}

void duckhts_simd_avx2_register(duckhts_simd_builder_t *builder) {
    duckhts_simd_builder_consider_base_counts(builder,
                                              DUCKHTS_SIMD_CAP_AVX2,
                                              "avx2",
                                              20,
                                              duckhts_base_counts_avx2);
}

int duckhts_simd_avx2_compiled(void) { return 1; }

int duckhts_simd_avx2_cpu_supported(void) {
    return duckhts_cpu_has_avx2();
}

int duckhts_simd_avx2_available(void) {
    return duckhts_simd_avx2_compiled() && duckhts_simd_avx2_cpu_supported();
}

#else
void duckhts_simd_avx2_register(duckhts_simd_builder_t *builder) { (void)builder; }
int duckhts_simd_avx2_compiled(void) { return 0; }
int duckhts_simd_avx2_cpu_supported(void) { return 0; }
int duckhts_simd_avx2_available(void) { return 0; }
#endif
