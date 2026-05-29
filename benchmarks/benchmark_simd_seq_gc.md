DuckHTS SIMD seq_gc_content Benchmark
================

<!-- benchmark_simd_seq_gc.md is generated from benchmark_simd_seq_gc.Rmd. -->

# Benchmark

This report measures the first DuckHTS runtime-SIMD dispatch path:
`seq_gc_content(...)` uses the immutable byte-oriented dispatch table
and the logical `seq_base_counts` kernel. The `auto` policy resolves
kernels from the compiled-and-CPU-supported capability mask, while
concrete backend requests measure that backend’s kernel implementation
with scalar fallback for missing slots. Unavailable requests are
reported as skipped so the same report remains valid on x86, ARM, and
wasm builds.

The goal is to decide whether a SIMD path is worth keeping before
expanding the same dispatch scaffold to quality-score scans, delimiter
scans, BAM packed-sequence GC counts, depth reductions, and variant
normalization helpers.

# Run

``` sh
make bench-simd-seq-gc
```

Useful overrides:

- `SIMD_SEQ_GC_ROWS`: number of rows, default `200000`
- `SIMD_SEQ_GC_LEN`: sequence length per row, default `512`
- `SIMD_SEQ_GC_RUNS`: timed repeats per backend, default `7`
- `SIMD_SEQ_GC_BACKENDS`: comma-separated backend requests, default
  `scalar,auto,avx2,avx512`
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
| backend_requests | scalar,auto,avx2,avx512                |
| pattern          | ACGTNNacgtnn                           |
| extension        | build/release/duckhts.duckdb_extension |

# Results

| backend_request | requested_backend | selected_backend | kernel_backend | kernel_capability | kernel_scalar_fallback | status  | median_sec | min_sec | mbases_per_sec_median | speedup_vs_scalar | skip_reason                                      |
|:----------------|:------------------|:-----------------|:---------------|:------------------|:-----------------------|:--------|-----------:|--------:|----------------------:|------------------:|:-------------------------------------------------|
| scalar          | scalar            | scalar           | scalar         | scalar            | FALSE                  | ok      |      0.033 |   0.033 |              3101.752 |             1.000 | NA                                               |
| auto            | auto              | avx2             | avx2           | avx2              | FALSE                  | ok      |      0.007 |   0.007 |             13882.232 |             4.476 | NA                                               |
| avx2            | avx2              | avx2             | avx2           | avx2              | FALSE                  | ok      |      0.007 |   0.007 |             13706.436 |             4.419 | NA                                               |
| avx512          | NA                | NA               | NA             | NA                | NA                     | skipped |         NA |      NA |                    NA |                NA | backend request is not available in this process |

# Interpretation

The fastest measured `seq_base_counts` kernel backend was `avx2`
selected from request `auto`, at 13882.2 Mbases/s.

That is 4.48x the scalar baseline for this synthetic workload.

The benchmark is intentionally synthetic: it isolates the byte-scanning
kernel and does not include FASTA/VCF/BAM I/O. Real end-to-end wins will
be smaller when decompression, SQL materialization, and joins dominate.

# Raw timing vectors

| backend_request | selected_backend | kernel_backend | iteration |  seconds |
|:----------------|:-----------------|:---------------|----------:|---------:|
| scalar          | scalar           | scalar         |         1 | 0.033094 |
| scalar          | scalar           | scalar         |         2 | 0.033038 |
| scalar          | scalar           | scalar         |         3 | 0.032769 |
| scalar          | scalar           | scalar         |         4 | 0.033014 |
| scalar          | scalar           | scalar         |         5 | 0.032885 |
| scalar          | scalar           | scalar         |         6 | 0.033200 |
| scalar          | scalar           | scalar         |         7 | 0.032910 |
| auto            | avx2             | avx2           |         1 | 0.007531 |
| auto            | avx2             | avx2           |         2 | 0.007415 |
| auto            | avx2             | avx2           |         3 | 0.007320 |
| auto            | avx2             | avx2           |         4 | 0.007272 |
| auto            | avx2             | avx2           |         5 | 0.007376 |
| auto            | avx2             | avx2           |         6 | 0.007291 |
| auto            | avx2             | avx2           |         7 | 0.007691 |
| avx2            | avx2             | avx2           |         1 | 0.007427 |
| avx2            | avx2             | avx2           |         2 | 0.007484 |
| avx2            | avx2             | avx2           |         3 | 0.007468 |
| avx2            | avx2             | avx2           |         4 | 0.007610 |
| avx2            | avx2             | avx2           |         5 | 0.007794 |
| avx2            | avx2             | avx2           |         6 | 0.007471 |
| avx2            | avx2             | avx2           |         7 | 0.007346 |
