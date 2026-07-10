#define _POSIX_C_SOURCE 200809L

/* Standalone contract tests and microbenchmarks for the byte-oriented SIMD
 * kernels.  Backend translation units are linked normally.  Their registrar
 * callbacks are the only supported way this program obtains private kernel
 * function pointers.
 */
#include <inttypes.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "duckhts_simd_internal.h"

#define ARRAY_LEN(a) (sizeof(a) / sizeof((a)[0]))
#define TEST_SEED UINT64_C(0xd1b54a32d192ed03)

typedef struct backend {
    const char *name;
    int (*compiled)(void);
    int (*cpu_supported)(void);
    duckhts_base_counts_fn base_counts;
    duckhts_bam_nt16_counts_fn bam_nt16_counts;
    duckhts_nt16_gc_counts_fn nt16_gc_counts;
} backend_t;

static int always_true(void) { return 1; }

static backend_t backends[] = {
    {"scalar", always_true, always_true, NULL, NULL, NULL},
    {"avx2", duckhts_simd_avx2_compiled, duckhts_simd_avx2_cpu_supported, NULL, NULL, NULL},
    {"neon", duckhts_simd_neon_compiled, duckhts_simd_neon_cpu_supported, NULL, NULL, NULL},
    {"wasm_simd128", duckhts_simd_wasm_simd128_compiled,
     duckhts_simd_wasm_simd128_cpu_supported, NULL, NULL, NULL}
};

static int registration_errors;

static backend_t *find_backend(const char *name) {
    size_t i;

    for (i = 0; i < ARRAY_LEN(backends); i++) {
        if (strcmp(backends[i].name, name) == 0) return &backends[i];
    }
    fprintf(stderr, "simd-contract: registrar supplied unknown backend=%s\n", name);
    registration_errors++;
    return NULL;
}

void duckhts_simd_builder_consider_base_counts(duckhts_simd_builder_t *builder,
                                               duckhts_simd_cap_t cap,
                                               const char *backend,
                                               int priority,
                                               duckhts_base_counts_fn fn) {
    backend_t *slot = find_backend(backend);
    (void)builder;
    (void)cap;
    (void)priority;
    if (slot == NULL) return;
    if (slot->base_counts != NULL) registration_errors++;
    slot->base_counts = fn;
}

void duckhts_simd_builder_consider_bam_nt16_counts(duckhts_simd_builder_t *builder,
                                                   duckhts_simd_cap_t cap,
                                                   const char *backend,
                                                   int priority,
                                                   duckhts_bam_nt16_counts_fn fn) {
    backend_t *slot = find_backend(backend);
    (void)builder;
    (void)cap;
    (void)priority;
    if (slot == NULL) return;
    if (slot->bam_nt16_counts != NULL) registration_errors++;
    slot->bam_nt16_counts = fn;
}

void duckhts_simd_builder_consider_nt16_gc_counts(duckhts_simd_builder_t *builder,
                                                  duckhts_simd_cap_t cap,
                                                  const char *backend,
                                                  int priority,
                                                  duckhts_nt16_gc_counts_fn fn) {
    backend_t *slot = find_backend(backend);
    (void)builder;
    (void)cap;
    (void)priority;
    if (slot == NULL) return;
    if (slot->nt16_gc_counts != NULL) registration_errors++;
    slot->nt16_gc_counts = fn;
}

static void capture_backends(void) {
    duckhts_simd_builder_t unused;

    memset(&unused, 0, sizeof(unused));
    duckhts_simd_scalar_register(&unused);
    duckhts_simd_avx2_register(&unused);
    duckhts_simd_neon_register(&unused);
    duckhts_simd_wasm_simd128_register(&unused);
}

static int backend_runnable(const backend_t *backend) {
    return backend->compiled() != 0 && backend->cpu_supported() != 0;
}

static uint64_t rng_next(uint64_t *state) {
    uint64_t x = *state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *state = x;
    return x * UINT64_C(2685821657736338717);
}

static void packed_oracle(const uint8_t *packed, int32_t n_bases,
                          duckhts_simd_bam_nt16_counts_t *out) {
    int32_t i;

    memset(out, 0, sizeof(*out));
    if (packed == NULL || n_bases <= 0) return;
    for (i = 0; i < n_bases; i++) {
        uint8_t byte = packed[(size_t)i >> 1];
        uint8_t code = ((i & 1) != 0) ? (uint8_t)(byte & UINT8_C(15)) :
                                       (uint8_t)(byte >> 4);
        switch (code) {
        case 1: out->a++; break;
        case 2: out->c++; break;
        case 4: out->g++; break;
        case 8: out->t++; break;
        case 15: out->n++; break;
        default: out->iupac++; break;
        }
    }
    out->gc = out->c + out->g;
    out->called = out->a + out->c + out->g + out->t;
}

static void unpacked_oracle(const uint8_t *codes, size_t n,
                            duckhts_simd_base_counts_t *out) {
    size_t i;
    uint64_t gc = 0;
    uint64_t called = 0;

    memset(out, 0, sizeof(*out));
    if (codes == NULL) return;
    for (i = 0; i < n; i++) {
        switch (codes[i]) {
        case 1: called++; break;
        case 2: gc++; called++; break;
        case 4: gc++; called++; break;
        case 8: called++; break;
        case 15: break;
        default:
            out->invalid = 1;
            return;
        }
    }
    out->gc = gc;
    out->called = called;
}

static int packed_equal(const duckhts_simd_bam_nt16_counts_t *a,
                        const duckhts_simd_bam_nt16_counts_t *b) {
    return a->a == b->a && a->c == b->c && a->g == b->g && a->t == b->t &&
           a->n == b->n && a->iupac == b->iupac && a->gc == b->gc &&
           a->called == b->called;
}

static int gc_equal(const duckhts_simd_base_counts_t *a,
                    const duckhts_simd_base_counts_t *b) {
    return a->gc == b->gc && a->called == b->called && a->invalid == b->invalid;
}

static void dump_bytes(const uint8_t *p, size_t n) {
    size_t i;
    size_t shown = n < 32 ? n : 32;

    for (i = 0; i < shown; i++) fprintf(stderr, "%02x", (unsigned)p[i]);
    if (shown < n) fprintf(stderr, "...");
}

static int check_packed_case(const uint8_t *packed, int32_t n_bases,
                             size_t offset, uint64_t case_no, const char *shape) {
    duckhts_simd_bam_nt16_counts_t want;
    size_t i;

    packed_oracle(packed, n_bases, &want);
    if (want.gc != want.c + want.g ||
        want.called != want.a + want.c + want.g + want.t ||
        want.a + want.c + want.g + want.t + want.n + want.iupac != (uint64_t)n_bases) {
        fprintf(stderr, "simd-contract: oracle structural failure kernel=bam_nt16_counts "
                "seed=%" PRIx64 " case=%" PRIu64 " shape=%s len=%" PRId32 "\n",
                TEST_SEED, case_no, shape, n_bases);
        return 1;
    }
    for (i = 0; i < ARRAY_LEN(backends); i++) {
        duckhts_simd_bam_nt16_counts_t got;
        backend_t *backend = &backends[i];
        if (!backend_runnable(backend)) continue;
        if (backend->bam_nt16_counts == NULL) {
            fprintf(stderr, "simd-contract: missing kernel=bam_nt16_counts backend=%s\n",
                    backend->name);
            return 1;
        }
        memset(&got, 0xa5, sizeof(got));
        backend->bam_nt16_counts(packed, n_bases, &got);
        if (!packed_equal(&want, &got) || got.gc != got.c + got.g ||
            got.called != got.a + got.c + got.g + got.t) {
            fprintf(stderr, "simd-contract: mismatch kernel=bam_nt16_counts backend=%s "
                    "seed=%" PRIx64 " case=%" PRIu64 " shape=%s len=%" PRId32
                    " offset=%zu bytes=", backend->name, TEST_SEED, case_no, shape,
                    n_bases, offset);
            dump_bytes(packed, ((size_t)n_bases + 1u) / 2u);
            fprintf(stderr, "\n  want={a=%" PRIu64 ",c=%" PRIu64 ",g=%" PRIu64
                    ",t=%" PRIu64 ",n=%" PRIu64 ",iupac=%" PRIu64
                    ",gc=%" PRIu64 ",called=%" PRIu64 "}\n"
                    "  got ={a=%" PRIu64 ",c=%" PRIu64 ",g=%" PRIu64
                    ",t=%" PRIu64 ",n=%" PRIu64 ",iupac=%" PRIu64
                    ",gc=%" PRIu64 ",called=%" PRIu64 "}\n",
                    want.a, want.c, want.g, want.t, want.n, want.iupac, want.gc, want.called,
                    got.a, got.c, got.g, got.t, got.n, got.iupac, got.gc, got.called);
            return 1;
        }
    }
    return 0;
}

static int check_unpacked_case(const uint8_t *codes, size_t n, size_t offset,
                               uint64_t case_no, const char *shape) {
    duckhts_simd_base_counts_t want;
    size_t i;

    unpacked_oracle(codes, n, &want);
    if ((!want.invalid && (want.gc > want.called || want.called > n)) ||
        (want.invalid && (want.gc != 0 || want.called != 0))) {
        fprintf(stderr, "simd-contract: oracle structural failure kernel=nt16_gc_counts "
                "seed=%" PRIx64 " case=%" PRIu64 " shape=%s len=%zu\n",
                TEST_SEED, case_no, shape, n);
        return 1;
    }
    for (i = 0; i < ARRAY_LEN(backends); i++) {
        duckhts_simd_base_counts_t got;
        backend_t *backend = &backends[i];
        if (!backend_runnable(backend)) continue;
        if (backend->nt16_gc_counts == NULL) {
            fprintf(stderr, "simd-contract: missing kernel=nt16_gc_counts backend=%s\n",
                    backend->name);
            return 1;
        }
        memset(&got, 0xa5, sizeof(got));
        backend->nt16_gc_counts(codes, n, &got);
        if (!gc_equal(&want, &got) || (!got.invalid && got.gc > got.called) ||
            (got.invalid && (got.gc != 0 || got.called != 0))) {
            fprintf(stderr, "simd-contract: mismatch kernel=nt16_gc_counts backend=%s "
                    "seed=%" PRIx64 " case=%" PRIu64 " shape=%s len=%zu offset=%zu codes=",
                    backend->name, TEST_SEED, case_no, shape, n, offset);
            dump_bytes(codes, n);
            fprintf(stderr, "\n  want={gc=%" PRIu64 ",called=%" PRIu64 ",invalid=%d}"
                    " got={gc=%" PRIu64 ",called=%" PRIu64 ",invalid=%d}\n",
                    want.gc, want.called, want.invalid, got.gc, got.called, got.invalid);
            return 1;
        }
    }
    return 0;
}

static int check_null_and_empty(void) {
    uint8_t byte = UINT8_C(0xff);
    size_t i;

    for (i = 0; i < ARRAY_LEN(backends); i++) {
        duckhts_simd_bam_nt16_counts_t packed;
        duckhts_simd_base_counts_t unpacked;
        backend_t *backend = &backends[i];
        if (!backend_runnable(backend)) continue;
        if (backend->bam_nt16_counts == NULL || backend->nt16_gc_counts == NULL) {
            fprintf(stderr, "simd-contract: incomplete nt16 registrar backend=%s\n", backend->name);
            return 1;
        }
        memset(&packed, 0xa5, sizeof(packed));
        backend->bam_nt16_counts(NULL, 17, &packed);
        {
            duckhts_simd_bam_nt16_counts_t zero;
            memset(&zero, 0, sizeof(zero));
            if (!packed_equal(&packed, &zero)) {
                fprintf(stderr, "simd-contract: nonzero NULL result kernel=bam_nt16_counts backend=%s\n",
                        backend->name);
                return 1;
            }
        }
        memset(&packed, 0xa5, sizeof(packed));
        backend->bam_nt16_counts(&byte, -1, &packed);
        {
            duckhts_simd_bam_nt16_counts_t zero;
            memset(&zero, 0, sizeof(zero));
            if (!packed_equal(&packed, &zero)) {
                fprintf(stderr, "simd-contract: nonzero negative-length result "
                        "kernel=bam_nt16_counts backend=%s\n", backend->name);
                return 1;
            }
        }
        memset(&unpacked, 0xa5, sizeof(unpacked));
        backend->nt16_gc_counts(NULL, 17u, &unpacked);
        {
            duckhts_simd_base_counts_t zero;
            memset(&zero, 0, sizeof(zero));
            if (!gc_equal(&unpacked, &zero)) {
                fprintf(stderr, "simd-contract: nonzero NULL result kernel=nt16_gc_counts backend=%s\n",
                        backend->name);
                return 1;
            }
        }
    }
    return 0;
}

static int run_contract_tests(void) {
    static const size_t lengths[] = {
        0, 1, 2, 3, 7, 15, 16, 17, 31, 32, 33, 47, 63, 64, 65,
        79, 95, 96, 97, 127, 128, 129, 191, 192, 193, 255, 256, 257
    };
    static const size_t offsets[] = {0, 1, 3, 7, 15, 31};
    uint8_t storage[8192 + 64];
    uint64_t state = TEST_SEED;
    uint64_t case_no = 0;
    size_t li;
    size_t oi;
    size_t i;

    if (registration_errors != 0) {
        fprintf(stderr, "simd-contract: %d registrar error(s)\n", registration_errors);
        return 1;
    }
    if (check_null_and_empty() != 0) return 1;

    /* Exhaust every packed nt16 code, every useful vector tail, unaligned
     * inputs, and poison the unused low nibble of odd-length records. */
    for (li = 0; li < ARRAY_LEN(lengths); li++) {
        size_t n = lengths[li];
        for (oi = 0; oi < ARRAY_LEN(offsets); oi++) {
            uint8_t *p = storage + offsets[oi];
            size_t bytes = (n + 1u) / 2u;
            memset(storage, 0xcd, sizeof(storage));
            for (i = 0; i < bytes; i++) {
                uint8_t hi = (uint8_t)((2u * i) & 15u);
                uint8_t lo = (uint8_t)((2u * i + 1u) & 15u);
                p[i] = (uint8_t)((uint8_t)(hi << 4) | lo);
            }
            if ((n & 1u) != 0) p[bytes - 1u] = (uint8_t)((p[bytes - 1u] & 0xf0u) | 0x0eu);
            if (check_packed_case(p, (int32_t)n, offsets[oi], case_no++, "boundary") != 0)
                return 1;
        }
    }

    /* Exhaust all byte values for the unpacked classifier, including nt16
     * ambiguity codes and bytes outside the nibble range. */
    for (i = 0; i <= UINT8_MAX; i++) {
        uint8_t code = (uint8_t)i;
        if (check_unpacked_case(&code, 1u, 0u, case_no++, "all-byte-values") != 0)
            return 1;
    }
    for (li = 0; li < ARRAY_LEN(lengths); li++) {
        size_t n = lengths[li];
        for (oi = 0; oi < ARRAY_LEN(offsets); oi++) {
            uint8_t *p = storage + offsets[oi];
            static const uint8_t legal[] = {1, 2, 4, 8, 15};
            memset(storage, 0xcd, sizeof(storage));
            for (i = 0; i < n; i++) p[i] = legal[i % ARRAY_LEN(legal)];
            if (check_unpacked_case(p, n, offsets[oi], case_no++, "boundary-valid") != 0)
                return 1;
            if (n != 0) {
                static const size_t positions[] = {0, 1, 15, 16, 31, 32, 63, 64};
                size_t pi;
                for (pi = 0; pi < ARRAY_LEN(positions); pi++) {
                    size_t pos = positions[pi] < n ? positions[pi] : n - 1u;
                    uint8_t saved = p[pos];
                    p[pos] = (pi & 1u) != 0 ? UINT8_C(3) : UINT8_C(255);
                    if (check_unpacked_case(p, n, offsets[oi], case_no++,
                                            "boundary-invalid") != 0) return 1;
                    p[pos] = saved;
                }
            }
        }
    }

    /* Deterministic fuzzing biases lengths around vector boundaries but also
     * covers multi-vector records. */
    for (i = 0; i < 3000u; i++) {
        size_t n = (size_t)(rng_next(&state) % 4097u);
        size_t offset = (size_t)(rng_next(&state) & 31u);
        uint8_t *p = storage + offset;
        size_t bytes = (n + 1u) / 2u;
        size_t j;
        for (j = 0; j < bytes; j++) p[j] = (uint8_t)rng_next(&state);
        if ((n & 1u) != 0) p[bytes - 1u] = (uint8_t)((p[bytes - 1u] & 0xf0u) |
                                                      (uint8_t)(rng_next(&state) & 15u));
        if (check_packed_case(p, (int32_t)n, offset, case_no++, "fuzz") != 0) return 1;

        for (j = 0; j < n; j++) {
            static const uint8_t pool[] = {1, 2, 4, 8, 15, 0, 3, 5, 14, 16, 255};
            p[j] = pool[rng_next(&state) % ARRAY_LEN(pool)];
        }
        if (check_unpacked_case(p, n, offset, case_no++, "fuzz") != 0) return 1;
    }

    printf("simd kernel contracts: OK (%" PRIu64 " cases, seed=%" PRIx64 ")\n",
           case_no, TEST_SEED);
    for (i = 0; i < ARRAY_LEN(backends); i++) {
        printf("  backend=%-12s compiled=%s cpu=%s tested=%s\n", backends[i].name,
               backends[i].compiled() ? "yes" : "no",
               backends[i].cpu_supported() ? "yes" : "no",
               backend_runnable(&backends[i]) ? "yes" : "no");
    }
    return 0;
}

static double now_seconds(void) {
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) return 0.0;
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static double bench_packed(duckhts_bam_nt16_counts_fn fn, const uint8_t *p,
                           int32_t n, int iterations,
                           duckhts_simd_bam_nt16_counts_t *result) {
    double best = 1e100;
    volatile uint64_t sink = 0;
    int i;
    for (i = 0; i < iterations; i++) {
        duckhts_simd_bam_nt16_counts_t out;
        double begin = now_seconds();
        fn(p, n, &out);
        {
            double elapsed = now_seconds() - begin;
            if (elapsed < best) best = elapsed;
        }
        sink += out.called;
        if (i == 0) *result = out;
    }
    (void)sink;
    return best;
}

static double bench_unpacked(duckhts_nt16_gc_counts_fn fn, const uint8_t *p,
                             size_t n, int iterations,
                             duckhts_simd_base_counts_t *result) {
    double best = 1e100;
    volatile uint64_t sink = 0;
    int i;
    for (i = 0; i < iterations; i++) {
        duckhts_simd_base_counts_t out;
        double begin = now_seconds();
        fn(p, n, &out);
        {
            double elapsed = now_seconds() - begin;
            if (elapsed < best) best = elapsed;
        }
        sink += out.called;
        if (i == 0) *result = out;
    }
    (void)sink;
    return best;
}

static int run_benchmarks(int32_t n_bases, int iterations) {
    static const uint8_t legal[] = {1, 2, 4, 8, 15};
    uint8_t *packed;
    uint8_t *unpacked;
    uint64_t state = TEST_SEED;
    size_t packed_bytes = ((size_t)n_bases + 1u) / 2u;
    size_t i;
    double scalar_packed = 0.0;
    double scalar_unpacked = 0.0;

    packed = (uint8_t *)malloc(packed_bytes == 0 ? 1u : packed_bytes);
    unpacked = (uint8_t *)malloc((size_t)n_bases == 0 ? 1u : (size_t)n_bases);
    if (packed == NULL || unpacked == NULL) {
        fprintf(stderr, "simd-bench: allocation failed\n");
        free(packed);
        free(unpacked);
        return 2;
    }
    memset(packed, 0, packed_bytes);
    for (i = 0; i < (size_t)n_bases; i++) {
        uint8_t code = legal[rng_next(&state) % ARRAY_LEN(legal)];
        unpacked[i] = code;
        if ((i & 1u) != 0) packed[i >> 1] |= code;
        else packed[i >> 1] = (uint8_t)(code << 4);
    }

    printf("\nSIMD kernel microbenchmark: bases=%" PRId32 " iterations=%d (best)\n",
           n_bases, iterations);
    printf("%-18s %-14s %12s %10s\n", "kernel", "backend", "Gbase/s", "speedup");
    for (i = 0; i < ARRAY_LEN(backends); i++) {
        duckhts_simd_bam_nt16_counts_t out;
        double elapsed;
        backend_t *backend = &backends[i];
        if (!backend_runnable(backend)) continue;
        elapsed = bench_packed(backend->bam_nt16_counts, packed, n_bases, iterations, &out);
        if (i == 0) scalar_packed = elapsed;
        printf("%-18s %-14s %12.3f %9.2fx\n", "bam_nt16_counts", backend->name,
               ((double)n_bases / 1e9) / elapsed, scalar_packed / elapsed);
    }
    for (i = 0; i < ARRAY_LEN(backends); i++) {
        duckhts_simd_base_counts_t out;
        double elapsed;
        backend_t *backend = &backends[i];
        if (!backend_runnable(backend)) continue;
        elapsed = bench_unpacked(backend->nt16_gc_counts, unpacked, (size_t)n_bases,
                                 iterations, &out);
        if (i == 0) scalar_unpacked = elapsed;
        printf("%-18s %-14s %12.3f %9.2fx\n", "nt16_gc_counts", backend->name,
               ((double)n_bases / 1e9) / elapsed, scalar_unpacked / elapsed);
    }
    free(packed);
    free(unpacked);
    return 0;
}

static int parse_positive(const char *text, long maximum, long *value) {
    char *end = NULL;
    long parsed = strtol(text, &end, 10);
    if (text[0] == '\0' || end == NULL || *end != '\0' || parsed <= 0 || parsed > maximum)
        return 0;
    *value = parsed;
    return 1;
}

int main(int argc, char **argv) {
    capture_backends();
    if (argc == 1 || (argc == 2 && strcmp(argv[1], "--check") == 0))
        return run_contract_tests();
    if (strcmp(argv[1], "--bench") == 0) {
        long bases = 50000000;
        long iterations = 5;
        if (argc > 4 || (argc >= 3 && !parse_positive(argv[2], INT32_MAX, &bases)) ||
            (argc == 4 && !parse_positive(argv[3], INT_MAX, &iterations))) {
            fprintf(stderr, "usage: %s [--check | --bench [bases [iterations]]]\n", argv[0]);
            return 2;
        }
        return run_benchmarks((int32_t)bases, (int)iterations);
    }
    fprintf(stderr, "usage: %s [--check | --bench [bases [iterations]]]\n", argv[0]);
    return 2;
}
