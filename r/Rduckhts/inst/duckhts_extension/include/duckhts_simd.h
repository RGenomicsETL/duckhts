#ifndef DUCKHTS_SIMD_H
#define DUCKHTS_SIMD_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    uint64_t gc;
    uint64_t called;
    int invalid;
} duckhts_simd_base_counts_t;

/* BAM packed-sequence (nt16) class counts.  Input is the htslib nibble-packed
 * sequence (`bam_get_seq`), two bases per byte.  Classes follow the BAM nt16
 * encoding (A=1, C=2, G=4, T=8, N=15) and reconcile to htslib `seq_nt16_int`
 * buckets: `called` is exactly A+C+G+T, everything else (N, `=`, IUPAC) is
 * not-ACGT, matching what samtools/htslib would count. */
typedef struct {
    uint64_t a;
    uint64_t c;
    uint64_t g;
    uint64_t t;
    uint64_t n;       /* nt16 15: true N and htslib-unrepresentable bytes */
    uint64_t iupac;   /* nt16 0 (=) and ambiguity codes: not ACGT, not N */
    uint64_t gc;      /* c + g */
    uint64_t called;  /* a + c + g + t */
} duckhts_simd_bam_nt16_counts_t;

/* Per-cycle FASTQ sufficient statistics. Callers allocate one entry per
 * possible read cycle; the kernel only mutates entries [0, len). */
typedef struct {
    uint64_t quality_sum;
    uint64_t a;
    uint64_t c;
    uint64_t g;
    uint64_t t;
    uint64_t n;
    uint64_t other;
    uint64_t a_quality_sum;
    uint64_t c_quality_sum;
    uint64_t g_quality_sum;
    uint64_t t_quality_sum;
} duckhts_simd_fastq_cycle_t;

/* One-read reductions produced while the cycle state is updated. Quality is
 * canonical Phred+33 text, as emitted by read_fastq(...,
 * quality_representation := 'string'). */
typedef struct {
    uint64_t q20;
    uint64_t q30;
    uint64_t q40;
    uint64_t quality_sum;
    uint64_t a;
    uint64_t c;
    uint64_t g;
    uint64_t t;
    uint64_t n;
    uint64_t other;
    int invalid_quality;
} duckhts_simd_fastq_read_t;

void duckhts_simd_init(void);
const char *duckhts_simd_selected_backend(void);
const char *duckhts_simd_requested_backend(void);
int duckhts_simd_backend_compiled(const char *backend);
int duckhts_simd_backend_cpu_supported(const char *backend);
int duckhts_simd_backend_available(const char *backend);
int duckhts_simd_set_backend(const char *backend, char *err, size_t err_len);
void duckhts_simd_base_counts(const char *seq, size_t len,
                              duckhts_simd_base_counts_t *out);
void duckhts_simd_bam_nt16_counts(const uint8_t *packed_seq, int32_t n_bases,
                                  duckhts_simd_bam_nt16_counts_t *out);
void duckhts_simd_fastq_qc_read(const char *sequence, const char *quality,
                                size_t len,
                                duckhts_simd_fastq_cycle_t *cycles,
                                duckhts_simd_fastq_read_t *out);

#endif /* DUCKHTS_SIMD_H */
