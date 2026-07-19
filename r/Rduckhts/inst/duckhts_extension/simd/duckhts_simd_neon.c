#include <stdint.h>
#include <string.h>

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

    for (; len - i >= 16; i += 16) {
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

/* BAM nt16 packed sequence: even bases occupy the high nibble, odd bases the
 * low nibble.  Process only complete 32-base blocks here so an odd final base
 * never observes the unused low nibble of its byte. */
static void duckhts_bam_nt16_counts_neon(const uint8_t *packed_seq,
                                         int32_t n_bases,
                                         duckhts_simd_bam_nt16_counts_t *out) {
    const uint8x16_t lo_mask = vdupq_n_u8(0x0Fu);
    const uint8x16_t v_a = vdupq_n_u8(1);
    const uint8x16_t v_c = vdupq_n_u8(2);
    const uint8x16_t v_g = vdupq_n_u8(4);
    const uint8x16_t v_t = vdupq_n_u8(8);
    const uint8x16_t v_n = vdupq_n_u8(15);
    uint64_t a = 0;
    uint64_t c = 0;
    uint64_t g = 0;
    uint64_t t = 0;
    uint64_t n = 0;
    size_t n_total;
    size_t i = 0;

    out->a = 0;
    out->c = 0;
    out->g = 0;
    out->t = 0;
    out->n = 0;
    out->iupac = 0;
    out->gc = 0;
    out->called = 0;
    if (!packed_seq || n_bases <= 0) return;
    n_total = (size_t)n_bases;

    for (; n_total - i >= 32; i += 32) {
        uint8x16_t bytes = vld1q_u8(packed_seq + (i >> 1));
        uint8x16_t hi = vshrq_n_u8(bytes, 4);
        uint8x16_t lo = vandq_u8(bytes, lo_mask);

        a += duckhts_neon_sum_u8(vshrq_n_u8(vceqq_u8(hi, v_a), 7));
        a += duckhts_neon_sum_u8(vshrq_n_u8(vceqq_u8(lo, v_a), 7));
        c += duckhts_neon_sum_u8(vshrq_n_u8(vceqq_u8(hi, v_c), 7));
        c += duckhts_neon_sum_u8(vshrq_n_u8(vceqq_u8(lo, v_c), 7));
        g += duckhts_neon_sum_u8(vshrq_n_u8(vceqq_u8(hi, v_g), 7));
        g += duckhts_neon_sum_u8(vshrq_n_u8(vceqq_u8(lo, v_g), 7));
        t += duckhts_neon_sum_u8(vshrq_n_u8(vceqq_u8(hi, v_t), 7));
        t += duckhts_neon_sum_u8(vshrq_n_u8(vceqq_u8(lo, v_t), 7));
        n += duckhts_neon_sum_u8(vshrq_n_u8(vceqq_u8(hi, v_n), 7));
        n += duckhts_neon_sum_u8(vshrq_n_u8(vceqq_u8(lo, v_n), 7));
    }

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

/* One unpacked nt16 code per byte.  A/C/G/T/N are the only valid codes; any
 * other byte invalidates the whole result, leaving gc/called at zero. */
static void duckhts_nt16_gc_counts_neon(const uint8_t *codes, size_t n,
                                        duckhts_simd_base_counts_t *out) {
    const uint8x16_t v_a = vdupq_n_u8(1);
    const uint8x16_t v_c = vdupq_n_u8(2);
    const uint8x16_t v_g = vdupq_n_u8(4);
    const uint8x16_t v_t = vdupq_n_u8(8);
    const uint8x16_t v_n = vdupq_n_u8(15);
    uint64_t gc = 0;
    uint64_t called = 0;
    size_t i = 0;

    out->gc = 0;
    out->called = 0;
    out->invalid = 0;
    if (!codes) return;

    for (; n - i >= 16; i += 16) {
        uint8x16_t x = vld1q_u8(codes + i);
        uint8x16_t a_mask = vceqq_u8(x, v_a);
        uint8x16_t c_mask = vceqq_u8(x, v_c);
        uint8x16_t g_mask = vceqq_u8(x, v_g);
        uint8x16_t t_mask = vceqq_u8(x, v_t);
        uint8x16_t n_mask = vceqq_u8(x, v_n);
        uint8x16_t called_mask = vorrq_u8(vorrq_u8(a_mask, c_mask),
                                          vorrq_u8(g_mask, t_mask));
        uint8x16_t valid_mask = vorrq_u8(called_mask, n_mask);

        if (duckhts_neon_sum_u8(vshrq_n_u8(valid_mask, 7)) != 16u) {
            out->invalid = 1;
            return;
        }
        gc += duckhts_neon_sum_u8(vshrq_n_u8(vorrq_u8(c_mask, g_mask), 7));
        called += duckhts_neon_sum_u8(vshrq_n_u8(called_mask, 7));
    }

    for (; i < n; i++) {
        switch (codes[i]) {
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

static inline void duckhts_fastq_qc_cycle_add_neon(
    duckhts_simd_fastq_cycle_t *cycle, unsigned char base, uint64_t phred) {
    uint64_t is_a = base == (unsigned char)'A';
    uint64_t is_c = base == (unsigned char)'C';
    uint64_t is_g = base == (unsigned char)'G';
    uint64_t is_t = base == (unsigned char)'T';
    uint64_t is_n = base == (unsigned char)'N';
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

static void duckhts_fastq_qc_neon(const char *sequence, const char *quality,
                                  size_t len,
                                  duckhts_simd_fastq_cycle_t *cycles,
                                  duckhts_simd_fastq_read_t *out) {
    const uint8x16_t case_mask = vdupq_n_u8(0xDFu);
    const uint8x16_t v_a = vdupq_n_u8((uint8_t)'A');
    const uint8x16_t v_c = vdupq_n_u8((uint8_t)'C');
    const uint8x16_t v_g = vdupq_n_u8((uint8_t)'G');
    const uint8x16_t v_t = vdupq_n_u8((uint8_t)'T');
    const uint8x16_t v_n = vdupq_n_u8((uint8_t)'N');
    const uint8x16_t v_bang = vdupq_n_u8((uint8_t)'!');
    const uint8x16_t v_tilde = vdupq_n_u8((uint8_t)'~');
    const uint8x16_t v_q20 = vdupq_n_u8((uint8_t)'5');
    const uint8x16_t v_q30 = vdupq_n_u8((uint8_t)'?');
    const uint8x16_t v_q40 = vdupq_n_u8((uint8_t)'I');
    size_t i = 0;

    memset(out, 0, sizeof(*out));
    if (!sequence || !quality || !cycles) return;

    for (; len - i >= 16; i += 16) {
        uint8x16_t seq = vld1q_u8((const uint8_t *)(const void *)(sequence + i));
        uint8x16_t qual = vld1q_u8((const uint8_t *)(const void *)(quality + i));
        uint8x16_t upper = vandq_u8(seq, case_mask);
        uint8x16_t a_mask = vceqq_u8(upper, v_a);
        uint8x16_t c_mask = vceqq_u8(upper, v_c);
        uint8x16_t g_mask = vceqq_u8(upper, v_g);
        uint8x16_t t_mask = vceqq_u8(upper, v_t);
        uint8x16_t n_mask = vceqq_u8(upper, v_n);
        uint8x16_t valid_base_mask = vorrq_u8(
            vorrq_u8(a_mask, c_mask), vorrq_u8(vorrq_u8(g_mask, t_mask), n_mask));
        uint8x16_t invalid_quality_mask = vorrq_u8(
            vcltq_u8(qual, v_bang), vcgtq_u8(qual, v_tilde));

        if (duckhts_neon_sum_u8(vshrq_n_u8(invalid_quality_mask, 7)) != 0) {
            out->invalid_quality = 1;
            return;
        }

        out->q20 += duckhts_neon_sum_u8(vshrq_n_u8(vcgeq_u8(qual, v_q20), 7));
        out->q30 += duckhts_neon_sum_u8(vshrq_n_u8(vcgeq_u8(qual, v_q30), 7));
        out->q40 += duckhts_neon_sum_u8(vshrq_n_u8(vcgeq_u8(qual, v_q40), 7));
        out->a += duckhts_neon_sum_u8(vshrq_n_u8(a_mask, 7));
        out->c += duckhts_neon_sum_u8(vshrq_n_u8(c_mask, 7));
        out->g += duckhts_neon_sum_u8(vshrq_n_u8(g_mask, 7));
        out->t += duckhts_neon_sum_u8(vshrq_n_u8(t_mask, 7));
        out->n += duckhts_neon_sum_u8(vshrq_n_u8(n_mask, 7));
        out->other += 16u - duckhts_neon_sum_u8(vshrq_n_u8(valid_base_mask, 7));

        for (unsigned lane = 0; lane < 16; lane++) {
            unsigned char base = (unsigned char)sequence[i + lane] & 0xDFu;
            uint64_t phred = (uint64_t)((unsigned char)quality[i + lane] - (unsigned char)'!');
            out->quality_sum += phred;
            duckhts_fastq_qc_cycle_add_neon(&cycles[i + lane], base, phred);
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
        duckhts_fastq_qc_cycle_add_neon(&cycles[i], base, phred);
    }
}

void duckhts_simd_neon_register(duckhts_simd_builder_t *builder) {
    duckhts_simd_builder_consider_base_counts(builder,
                                              DUCKHTS_SIMD_CAP_NEON,
                                              "neon",
                                              50,
                                              duckhts_base_counts_neon);
    duckhts_simd_builder_consider_bam_nt16_counts(builder,
                                                  DUCKHTS_SIMD_CAP_NEON,
                                                  "neon",
                                                  50,
                                                  duckhts_bam_nt16_counts_neon);
    duckhts_simd_builder_consider_nt16_gc_counts(builder,
                                                 DUCKHTS_SIMD_CAP_NEON,
                                                 "neon",
                                                 50,
                                                 duckhts_nt16_gc_counts_neon);
    duckhts_simd_builder_consider_fastq_qc(builder,
                                           DUCKHTS_SIMD_CAP_NEON,
                                           "neon",
                                           50,
                                           duckhts_fastq_qc_neon);
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
