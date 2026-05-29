DuckHTS SIMD seq_gc_content Benchmark
================

<!-- benchmark_simd_seq_gc.md is generated from benchmark_simd_seq_gc.Rmd. -->

# Benchmark

This report measures the first DuckHTS runtime-SIMD dispatch path:
`seq_gc_content(...)` uses the shared byte-oriented base-count vtable
with a scalar fallback and optional platform-specific SIMD
implementations such as AVX2. The benchmark requests backends explicitly
through `duckhts_simd_set_backend(...)`; unavailable requests are
reported as skipped so the same report remains valid on x86, ARM, and
wasm builds.

The goal is to decide whether a SIMD path is worth keeping before
expanding the same dispatch scaffold to quality-score scans, BAM
packed-sequence GC counts, depth reductions, and variant normalization
helpers.

# Run

``` sh
make bench-simd-seq-gc
```

Useful overrides:

- `SIMD_SEQ_GC_ROWS`: number of rows, default `200000`
- `SIMD_SEQ_GC_LEN`: sequence length per row, default `512`
- `SIMD_SEQ_GC_RUNS`: timed repeats per backend, default `7`
- `SIMD_SEQ_GC_BACKENDS`: comma-separated backend requests, default
  `scalar,auto,avx2`
- `SIMD_SEQ_GC_PATTERN`: repeated sequence pattern, default
  `ACGTNNacgtnn`
- `SIMD_SEQ_GC_FORCE=1`: rerun even if cached JSON exists
- `DUCKHTS_EXTENSION`: extension path override

# Configuration

| setting          | value                                  |
|:-----------------|:---------------------------------------|
| rows             | 200,000                                |
| sequence_length  | 512                                    |
| total_bases      | 102,400,000                            |
| iterations       | 7                                      |
| backend_requests | scalar,auto,avx2                       |
| pattern          | ACGTNNacgtnn                           |
| extension        | build/release/duckhts.duckdb_extension |

# Results

| backend_request | requested_backend | selected_backend | status | median_sec | min_sec | mbases_per_sec_median | speedup_vs_scalar | skip_reason |
|:----------------|:------------------|:-----------------|:-------|-----------:|--------:|----------------------:|------------------:|:------------|
| scalar          | scalar            | scalar           | ok     |      0.033 |   0.033 |              3063.023 |             1.000 | NA          |
| auto            | auto              | avx2             | ok     |      0.007 |   0.007 |             13753.426 |             4.490 | NA          |
| avx2            | avx2              | avx2             | ok     |      0.007 |   0.007 |             13865.699 |             4.527 | NA          |

# Interpretation

The fastest measured backend was `avx2` selected from request `avx2`, at
13865.7 Mbases/s.

That is 4.53x the scalar baseline for this synthetic workload.

The benchmark is intentionally synthetic: it isolates the byte-scanning
kernel and does not include FASTA/VCF/BAM I/O. Real end-to-end wins will
be smaller when decompression, SQL materialization, and joins dominate.

# Raw timing vectors

| backend_request | selected_backend | iteration |  seconds |
|:----------------|:-----------------|----------:|---------:|
| scalar          | scalar           |         1 | 0.032968 |
| scalar          | scalar           |         2 | 0.033207 |
| scalar          | scalar           |         3 | 0.033448 |
| scalar          | scalar           |         4 | 0.033431 |
| scalar          | scalar           |         5 | 0.033472 |
| scalar          | scalar           |         6 | 0.033747 |
| scalar          | scalar           |         7 | 0.033162 |
| auto            | avx2             |         1 | 0.007593 |
| auto            | avx2             |         2 | 0.007727 |
| auto            | avx2             |         3 | 0.007348 |
| auto            | avx2             |         4 | 0.007426 |
| auto            | avx2             |         5 | 0.007245 |
| auto            | avx2             |         6 | 0.007445 |
| auto            | avx2             |         7 | 0.007709 |
| avx2            | avx2             |         1 | 0.007555 |
| avx2            | avx2             |         2 | 0.007414 |
| avx2            | avx2             |         3 | 0.007475 |
| avx2            | avx2             |         4 | 0.007385 |
| avx2            | avx2             |         5 | 0.007228 |
| avx2            | avx2             |         6 | 0.007283 |
| avx2            | avx2             |         7 | 0.007237 |
