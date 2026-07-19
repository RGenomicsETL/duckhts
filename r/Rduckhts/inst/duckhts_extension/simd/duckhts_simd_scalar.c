#include <stdint.h>
#include <string.h>

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

/* BAM nt16 packed sequence: two bases per byte, even index in the high nibble
 * (matching htslib `bam_seqi`).  Classes are the BAM nt16 codes (A=1, C=2,
 * G=4, T=8, N=15); every other code (0 `=` and the IUPAC ambiguity codes)
 * lands in `iupac`, reconciling to htslib `seq_nt16_int != 4` for `called`. */
static void duckhts_bam_nt16_counts_scalar(const uint8_t *packed_seq,
                                           int32_t n_bases,
                                           duckhts_simd_bam_nt16_counts_t *out) {
    /* nt16 code (0..15) -> class bucket: 0=A 1=C 2=G 3=T 4=N 5=IUPAC/'='.
     * Branchless histogram (no data-dependent branch) so the oracle reflects a
     * real classifier rather than a mispredict-bound switch. */
    static const uint8_t cls[16] = {
        5, 0, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 4
    };
    uint64_t h[6] = {0, 0, 0, 0, 0, 0};

    out->a = out->c = out->g = out->t = 0;
    out->n = out->iupac = out->gc = out->called = 0;
    if (!packed_seq || n_bases <= 0) return;

    for (int32_t i = 0; i < n_bases; i++) {
        uint8_t byte = packed_seq[i >> 1];
        unsigned code = (i & 1) ? (byte & 0x0Fu) : (unsigned)(byte >> 4);
        h[cls[code]]++;
    }

    out->a = h[0];
    out->c = h[1];
    out->g = h[2];
    out->t = h[3];
    out->n = h[4];
    out->iupac = h[5];
    out->gc = h[1] + h[2];
    out->called = h[0] + h[1] + h[2] + h[3];
}

/* Unpacked nt16 GC classify (reference oracle). One nt16 code per byte.
   class LUT: 0 = invalid (=, IUPAC), 1 = N, 2 = A/T (called), 3 = C/G
   (called + gc). Any invalid code sets out->invalid and returns early, so
   the result matches the text seq_gc_content path (which returns NULL). */
static void duckhts_nt16_gc_counts_scalar(const uint8_t *codes, size_t n,
                                          duckhts_simd_base_counts_t *out) {
    static const uint8_t nt16_class[16] = {0, 2, 3, 0, 3, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1};
    uint64_t gc = 0;
    uint64_t called = 0;
    size_t i;

    out->gc = 0;
    out->called = 0;
    out->invalid = 0;
    if (!codes) return;

    for (i = 0; i < n; i++) {
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
    out->gc = gc;
    out->called = called;
}

static void duckhts_fastq_qc_scalar(const char *sequence, const char *quality,
                                    size_t len,
                                    duckhts_simd_fastq_cycle_t *cycles,
                                    duckhts_simd_fastq_read_t *out) {
    memset(out, 0, sizeof(*out));
    if (!sequence || !quality || !cycles) return;

    for (size_t i = 0; i < len; i++) {
        duckhts_simd_fastq_cycle_t *cycle = &cycles[i];
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
}

void duckhts_simd_scalar_register(duckhts_simd_builder_t *builder) {
    duckhts_simd_builder_consider_base_counts(builder,
                                              DUCKHTS_SIMD_CAP_SCALAR,
                                              "scalar",
                                              1000,
                                              duckhts_base_counts_scalar);
    duckhts_simd_builder_consider_bam_nt16_counts(builder,
                                                  DUCKHTS_SIMD_CAP_SCALAR,
                                                  "scalar",
                                                  1000,
                                                  duckhts_bam_nt16_counts_scalar);
    duckhts_simd_builder_consider_nt16_gc_counts(builder,
                                                 DUCKHTS_SIMD_CAP_SCALAR,
                                                 "scalar",
                                                 1000,
                                                 duckhts_nt16_gc_counts_scalar);
    duckhts_simd_builder_consider_fastq_qc(builder,
                                           DUCKHTS_SIMD_CAP_SCALAR,
                                           "scalar",
                                           1000,
                                           duckhts_fastq_qc_scalar);
}
