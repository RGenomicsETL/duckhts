#include <stdint.h>
#include <string.h>

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

    for (; len - i >= 32; i += 32) {
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

/* BAM nt16 packed sequence: two 4-bit codes per byte.  Vectorize the body in
 * 64-base (32-byte) blocks where both nibbles of every byte are valid bases,
 * then clean up the tail (and any odd final base) in scalar.  Split each byte
 * into its high and low nibble, compare against the five ACGTN codes, and
 * popcount the compare masks; IUPAC/`=` is the processed remainder. */
__attribute__((target("avx2,popcnt")))
static void duckhts_bam_nt16_counts_avx2(const uint8_t *packed_seq,
                                         int32_t n_bases,
                                         duckhts_simd_bam_nt16_counts_t *out) {
    const __m256i lo_mask = _mm256_set1_epi8(0x0F);
    const __m256i v_a = _mm256_set1_epi8(1);
    const __m256i v_c = _mm256_set1_epi8(2);
    const __m256i v_g = _mm256_set1_epi8(4);
    const __m256i v_t = _mm256_set1_epi8(8);
    const __m256i v_n = _mm256_set1_epi8(15);
    uint64_t a = 0, c = 0, g = 0, t = 0, n = 0;
    size_t n_total;
    size_t i = 0;

    out->a = out->c = out->g = out->t = 0;
    out->n = out->iupac = out->gc = out->called = 0;
    if (!packed_seq || n_bases <= 0) return;
    n_total = (size_t)n_bases;

    for (; n_total - i >= 64; i += 64) {
        __m256i bytes = _mm256_loadu_si256((const __m256i *)(const void *)(packed_seq + (i >> 1)));
        __m256i hi = _mm256_and_si256(_mm256_srli_epi16(bytes, 4), lo_mask);
        __m256i lo = _mm256_and_si256(bytes, lo_mask);
        a += (uint64_t)DUCKHTS_POPCOUNT32((uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(hi, v_a)));
        a += (uint64_t)DUCKHTS_POPCOUNT32((uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(lo, v_a)));
        c += (uint64_t)DUCKHTS_POPCOUNT32((uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(hi, v_c)));
        c += (uint64_t)DUCKHTS_POPCOUNT32((uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(lo, v_c)));
        g += (uint64_t)DUCKHTS_POPCOUNT32((uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(hi, v_g)));
        g += (uint64_t)DUCKHTS_POPCOUNT32((uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(lo, v_g)));
        t += (uint64_t)DUCKHTS_POPCOUNT32((uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(hi, v_t)));
        t += (uint64_t)DUCKHTS_POPCOUNT32((uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(lo, v_t)));
        n += (uint64_t)DUCKHTS_POPCOUNT32((uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(hi, v_n)));
        n += (uint64_t)DUCKHTS_POPCOUNT32((uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(lo, v_n)));
    }

    /* IUPAC/`=` for the vectorized region is the processed remainder. */
    out->iupac = (uint64_t)i - (a + c + g + t + n);

    for (; i < n_total; i++) {
        uint8_t byte = packed_seq[i >> 1];
        unsigned code = (i & 1u) ? (unsigned)(byte & 0x0Fu)
                                 : (unsigned)(byte >> 4);
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

static int duckhts_cpu_has_avx2(void) {
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx2") != 0 && __builtin_cpu_supports("popcnt") != 0;
}

/* Unpacked nt16 GC classify. One nt16 code per byte, so no nibble unpack:
   compare each byte against A/C/G/T/N. If any 32-byte block contains a byte
   matching none of them (an invalid `=`/IUPAC code), set invalid and return,
   matching the scalar reference and the text path (NULL on non-ACGTN). */
__attribute__((target("avx2,popcnt")))
static void duckhts_nt16_gc_counts_avx2(const uint8_t *codes, size_t n,
                                        duckhts_simd_base_counts_t *out) {
    const __m256i v_a = _mm256_set1_epi8(1);
    const __m256i v_c = _mm256_set1_epi8(2);
    const __m256i v_g = _mm256_set1_epi8(4);
    const __m256i v_t = _mm256_set1_epi8(8);
    const __m256i v_n = _mm256_set1_epi8(15);
    uint64_t gc = 0;
    uint64_t called = 0;
    size_t i = 0;

    out->gc = 0;
    out->called = 0;
    out->invalid = 0;
    if (!codes) return;

    for (; n - i >= 32; i += 32) {
        __m256i x = _mm256_loadu_si256((const __m256i *)(const void *)(codes + i));
        uint32_t ma = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(x, v_a));
        uint32_t mc = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(x, v_c));
        uint32_t mg = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(x, v_g));
        uint32_t mt = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(x, v_t));
        uint32_t mn = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(x, v_n));
        if ((ma | mc | mg | mt | mn) != 0xFFFFFFFFu) {
            out->invalid = 1;
            return;
        }
        gc += (uint64_t)DUCKHTS_POPCOUNT32(mc | mg);
        called += (uint64_t)DUCKHTS_POPCOUNT32(ma | mc | mg | mt);
    }
    {
        static const uint8_t nt16_class[16] = {0, 2, 3, 0, 3, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1};
        for (; i < n; i++) {
            uint8_t code = codes[i];
            uint8_t klass = (code < 16) ? nt16_class[code] : 0;
            if (klass == 0) {
                out->invalid = 1;
                return;
            }
            if (klass >= 2) {
                called++;
                if (klass == 3) gc++;
            }
        }
    }
    out->gc = gc;
    out->called = called;
}

static inline void duckhts_fastq_qc_cycle_add(duckhts_simd_fastq_cycle_t *cycle,
                                               uint64_t phred,
                                               uint32_t bit,
                                               uint32_t a_mask,
                                               uint32_t c_mask,
                                               uint32_t g_mask,
                                               uint32_t t_mask,
                                               uint32_t n_mask) {
    uint64_t is_a = (a_mask & bit) != 0;
    uint64_t is_c = (c_mask & bit) != 0;
    uint64_t is_g = (g_mask & bit) != 0;
    uint64_t is_t = (t_mask & bit) != 0;
    uint64_t is_n = (n_mask & bit) != 0;
    uint64_t is_other = (is_a | is_c | is_g | is_t | is_n) ^ 1u;

    cycle->quality_sum += phred;
    cycle->a += is_a;
    cycle->c += is_c;
    cycle->g += is_g;
    cycle->t += is_t;
    cycle->n += is_n;
    cycle->other += is_other;
    cycle->a_quality_sum += phred * is_a;
    cycle->c_quality_sum += phred * is_c;
    cycle->g_quality_sum += phred * is_g;
    cycle->t_quality_sum += phred * is_t;
}

__attribute__((target("avx2,popcnt")))
static void duckhts_fastq_qc_avx2(const char *sequence, const char *quality,
                                  size_t len,
                                  duckhts_simd_fastq_cycle_t *cycles,
                                  duckhts_simd_fastq_read_t *out) {
    const __m256i case_mask = _mm256_set1_epi8((char)0xDF);
    const __m256i v_a = _mm256_set1_epi8('A');
    const __m256i v_c = _mm256_set1_epi8('C');
    const __m256i v_g = _mm256_set1_epi8('G');
    const __m256i v_t = _mm256_set1_epi8('T');
    const __m256i v_n = _mm256_set1_epi8('N');
    const __m256i v_bang = _mm256_set1_epi8('!');
    const __m256i v_tilde = _mm256_set1_epi8('~');
    const __m256i v_q20_minus_one = _mm256_set1_epi8('5' - 1);
    const __m256i v_q30_minus_one = _mm256_set1_epi8('?' - 1);
    const __m256i v_q40_minus_one = _mm256_set1_epi8('I' - 1);
    size_t i = 0;

    memset(out, 0, sizeof(*out));
    if (!sequence || !quality || !cycles) return;

    for (; len - i >= 32; i += 32) {
        __m256i seq = _mm256_loadu_si256((const __m256i *)(const void *)(sequence + i));
        __m256i qual = _mm256_loadu_si256((const __m256i *)(const void *)(quality + i));
        __m256i upper = _mm256_and_si256(seq, case_mask);
        uint32_t a_mask = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(upper, v_a));
        uint32_t c_mask = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(upper, v_c));
        uint32_t g_mask = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(upper, v_g));
        uint32_t t_mask = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(upper, v_t));
        uint32_t n_mask = (uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(upper, v_n));
        uint32_t valid_base_mask = a_mask | c_mask | g_mask | t_mask | n_mask;
        uint32_t invalid_quality_mask = (uint32_t)_mm256_movemask_epi8(
            _mm256_or_si256(_mm256_cmpgt_epi8(v_bang, qual),
                            _mm256_cmpgt_epi8(qual, v_tilde)));

        if (invalid_quality_mask != 0) {
            out->invalid_quality = 1;
            return;
        }

        out->q20 += (uint64_t)DUCKHTS_POPCOUNT32((uint32_t)_mm256_movemask_epi8(
            _mm256_cmpgt_epi8(qual, v_q20_minus_one)));
        out->q30 += (uint64_t)DUCKHTS_POPCOUNT32((uint32_t)_mm256_movemask_epi8(
            _mm256_cmpgt_epi8(qual, v_q30_minus_one)));
        out->q40 += (uint64_t)DUCKHTS_POPCOUNT32((uint32_t)_mm256_movemask_epi8(
            _mm256_cmpgt_epi8(qual, v_q40_minus_one)));
        out->a += (uint64_t)DUCKHTS_POPCOUNT32(a_mask);
        out->c += (uint64_t)DUCKHTS_POPCOUNT32(c_mask);
        out->g += (uint64_t)DUCKHTS_POPCOUNT32(g_mask);
        out->t += (uint64_t)DUCKHTS_POPCOUNT32(t_mask);
        out->n += (uint64_t)DUCKHTS_POPCOUNT32(n_mask);
        out->other += (uint64_t)DUCKHTS_POPCOUNT32(~valid_base_mask);

        for (unsigned lane = 0; lane < 32; lane++) {
            uint64_t phred = (uint64_t)((unsigned char)quality[i + lane] - (unsigned char)'!');
            uint32_t bit = (uint32_t)1u << lane;
            out->quality_sum += phred;
            duckhts_fastq_qc_cycle_add(&cycles[i + lane], phred, bit,
                                       a_mask, c_mask, g_mask, t_mask, n_mask);
        }
    }

    for (; i < len; i++) {
        unsigned char base = (unsigned char)sequence[i] & 0xDFu;
        unsigned char ascii_quality = (unsigned char)quality[i];
        uint64_t is_a = base == (unsigned char)'A';
        uint64_t is_c = base == (unsigned char)'C';
        uint64_t is_g = base == (unsigned char)'G';
        uint64_t is_t = base == (unsigned char)'T';
        uint64_t is_n = base == (unsigned char)'N';
        uint64_t is_other = (is_a | is_c | is_g | is_t | is_n) ^ 1u;
        uint64_t phred;

        if (ascii_quality < (unsigned char)'!' || ascii_quality > (unsigned char)'~') {
            out->invalid_quality = 1;
            return;
        }
        phred = (uint64_t)(ascii_quality - (unsigned char)'!');
        out->q20 += ascii_quality >= (unsigned char)'5';
        out->q30 += ascii_quality >= (unsigned char)'?';
        out->q40 += ascii_quality >= (unsigned char)'I';
        out->quality_sum += phred;
        out->a += is_a;
        out->c += is_c;
        out->g += is_g;
        out->t += is_t;
        out->n += is_n;
        out->other += is_other;
        cycles[i].quality_sum += phred;
        cycles[i].a += is_a;
        cycles[i].c += is_c;
        cycles[i].g += is_g;
        cycles[i].t += is_t;
        cycles[i].n += is_n;
        cycles[i].other += is_other;
        cycles[i].a_quality_sum += phred * is_a;
        cycles[i].c_quality_sum += phred * is_c;
        cycles[i].g_quality_sum += phred * is_g;
        cycles[i].t_quality_sum += phred * is_t;
    }
}

void duckhts_simd_avx2_register(duckhts_simd_builder_t *builder) {
    duckhts_simd_builder_consider_base_counts(builder,
                                              DUCKHTS_SIMD_CAP_AVX2,
                                              "avx2",
                                              20,
                                              duckhts_base_counts_avx2);
    duckhts_simd_builder_consider_bam_nt16_counts(builder,
                                                  DUCKHTS_SIMD_CAP_AVX2,
                                                  "avx2",
                                                  20,
                                                  duckhts_bam_nt16_counts_avx2);
    duckhts_simd_builder_consider_nt16_gc_counts(builder,
                                                 DUCKHTS_SIMD_CAP_AVX2,
                                                 "avx2",
                                                 20,
                                                 duckhts_nt16_gc_counts_avx2);
    duckhts_simd_builder_consider_fastq_qc(builder,
                                           DUCKHTS_SIMD_CAP_AVX2,
                                           "avx2",
                                           20,
                                           duckhts_fastq_qc_avx2);
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
