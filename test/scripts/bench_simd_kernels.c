/* Standalone C microbenchmark for DuckHTS SIMD logical kernels.
 *
 * This harness isolates the byte-oriented kernels from DuckDB/R so it measures
 * the kernel itself, not query or transport overhead.  It compiles the kernel
 * translation units directly (textual include) so it can call each backend's
 * concrete implementation without the DuckDB-coupled dispatch.c, then for every
 * compiled+cpu-supported backend it:
 *   1. verifies the backend's counts are bit-identical to the scalar oracle
 *      (correctness gate -- a divergence fails the run), and
 *   2. times throughput and reports the scalar-relative speedup.
 *
 * It is a developer tool, not part of any product build.  Build and run with:
 *     make bench-simd-kernels
 *
 * Per the dispatch-matrix contract, scalar is the correctness oracle and any
 * backend whose result differs from scalar is a bug, not a faster path.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* The backend register functions call duckhts_simd_builder_consider_*, which
 * live in the DuckDB-coupled dispatch.c.  We never invoke the registrars here
 * (we call the kernels directly), so provide no-op stubs to satisfy the link. */
#include "duckhts_simd_internal.h"
void duckhts_simd_builder_consider_base_counts(duckhts_simd_builder_t *b, duckhts_simd_cap_t c,
                                               const char *n, int p, duckhts_base_counts_fn f) {
    (void)b; (void)c; (void)n; (void)p; (void)f;
}
void duckhts_simd_builder_consider_bam_nt16_counts(duckhts_simd_builder_t *b, duckhts_simd_cap_t c,
                                                   const char *n, int p, duckhts_bam_nt16_counts_fn f) {
    (void)b; (void)c; (void)n; (void)p; (void)f;
}

/* Pull in the concrete kernels.  On non-x86 hosts the AVX2 TU compiles to its
 * unavailable stubs and we simply skip that backend below. */
#include "../../src/simd/duckhts_simd_scalar.c"
#include "../../src/simd/duckhts_simd_avx2.c"

/* Forward decls for the static kernels we benchmark (defined in the includes). */
static void duckhts_bam_nt16_counts_scalar(const uint8_t *, int32_t, duckhts_simd_bam_nt16_counts_t *);
#if defined(DUCKHTS_SIMD_COMPILED_AVX2)
__attribute__((target("avx2,popcnt")))
static void duckhts_bam_nt16_counts_avx2(const uint8_t *, int32_t, duckhts_simd_bam_nt16_counts_t *);
#endif

typedef void (*nt16_fn)(const uint8_t *, int32_t, duckhts_simd_bam_nt16_counts_t *);

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Pack n_bases of a realistic nt16 distribution (~mostly ACGT, ~1% N, a few
 * IUPAC) into ceil(n_bases/2) bytes, even base in the high nibble (bam_seqi). */
static uint8_t *make_packed_nt16(int32_t n_bases) {
    static const uint8_t codes[16] = {1,2,4,8, 1,2,4,8, 1,2,4,8, 1,2,15,6}; /* mostly ACGT, one N, one S */
    size_t nbytes = (size_t)((n_bases + 1) / 2);
    uint8_t *buf = (uint8_t *)malloc(nbytes ? nbytes : 1);
    uint64_t rng = 0x9E3779B97F4A7C15ull;
    if (!buf) return NULL;
    memset(buf, 0, nbytes);
    for (int32_t i = 0; i < n_bases; i++) {
        rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
        uint8_t code = codes[rng & 0x0F];
        if (i & 1) buf[i >> 1] |= code;
        else       buf[i >> 1] |= (uint8_t)(code << 4);
    }
    return buf;
}

static int counts_equal(const duckhts_simd_bam_nt16_counts_t *a,
                        const duckhts_simd_bam_nt16_counts_t *b) {
    return a->a == b->a && a->c == b->c && a->g == b->g && a->t == b->t &&
           a->n == b->n && a->iupac == b->iupac && a->gc == b->gc && a->called == b->called;
}

static double bench_backend(nt16_fn fn, const uint8_t *buf, int32_t n_bases,
                           int iters, duckhts_simd_bam_nt16_counts_t *out) {
    volatile uint64_t sink = 0;
    double best = 1e30;
    for (int it = 0; it < iters; it++) {
        duckhts_simd_bam_nt16_counts_t c;
        double t0 = now_sec();
        fn(buf, n_bases, &c);
        double dt = now_sec() - t0;
        sink += c.gc;            /* defeat dead-code elimination */
        if (dt < best) best = dt;
        if (it == 0) *out = c;
    }
    (void)sink;
    return best;
}

int main(int argc, char **argv) {
    int32_t n_bases = (argc > 1) ? atoi(argv[1]) : 200000000; /* 200 Mbase ~ 100 MB packed */
    int iters = (argc > 2) ? atoi(argv[2]) : 7;
    uint8_t *buf = make_packed_nt16(n_bases);
    duckhts_simd_bam_nt16_counts_t ref, got;
    double gbase = (double)n_bases / 1e9;
    double scalar_t, t;
    int rc = 0;

    if (!buf) { fprintf(stderr, "bench: allocation failed\n"); return 2; }

    printf("DuckHTS SIMD kernel microbenchmark\n");
    printf("kernel=bam_nt16_counts  bases=%d (%.1f MB packed)  iters=%d (best-of)\n\n",
           n_bases, (double)((n_bases + 1) / 2) / 1e6, iters);
    printf("%-10s  %-9s  %12s  %10s  %s\n", "backend", "available", "Gbase/s", "speedup", "correct");
    printf("%-10s  %-9s  %12s  %10s  %s\n", "----------", "---------", "------------", "----------", "-------");

    scalar_t = bench_backend(duckhts_bam_nt16_counts_scalar, buf, n_bases, iters, &ref);
    printf("%-10s  %-9s  %12.3f  %9.2fx  %s\n", "scalar", "yes",
           gbase / scalar_t, 1.0, "oracle");

#if defined(DUCKHTS_SIMD_COMPILED_AVX2)
    if (duckhts_simd_avx2_cpu_supported()) {
        t = bench_backend(duckhts_bam_nt16_counts_avx2, buf, n_bases, iters, &got);
        int ok = counts_equal(&ref, &got);
        printf("%-10s  %-9s  %12.3f  %9.2fx  %s\n", "avx2", "yes",
               gbase / t, scalar_t / t, ok ? "OK" : "MISMATCH");
        if (!ok) {
            fprintf(stderr, "bench: avx2 diverged from scalar oracle\n");
            rc = 1;
        }
    } else {
        printf("%-10s  %-9s  %12s  %10s  %s\n", "avx2", "no (cpu)", "-", "-", "-");
    }
#else
    printf("%-10s  %-9s  %12s  %10s  %s\n", "avx2", "no (build)", "-", "-", "-");
#endif

    printf("\noracle counts: A=%llu C=%llu G=%llu T=%llu N=%llu IUPAC=%llu  gc=%llu called=%llu\n",
           (unsigned long long)ref.a, (unsigned long long)ref.c, (unsigned long long)ref.g,
           (unsigned long long)ref.t, (unsigned long long)ref.n, (unsigned long long)ref.iupac,
           (unsigned long long)ref.gc, (unsigned long long)ref.called);
    free(buf);
    return rc;
}
