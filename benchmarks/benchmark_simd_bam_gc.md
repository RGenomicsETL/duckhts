DuckHTS SIMD BAM GC-content Benchmark
================

<!-- benchmark_simd_bam_gc.md is generated from benchmark_simd_bam_gc.Rmd. -->

# Benchmark

This report measures `seq_gc_content(SEQ)` on real BAM read sequences.
It is intended to complement the synthetic `seq_gc_content(...)`
benchmark: the input strings come from `read_bam(...)`, so read lengths,
lowercase/uppercase handling, `N` bases, DuckDB vector flow, and BAM
decoding overhead are closer to a real HTS workload.

Two workloads are reported for each SIMD backend request:

- `bam_scan`: scan the BAM and compute `seq_gc_content(SEQ)` in the same
  query. This is the end-to-end path and includes BAM decoding,
  decompression, DuckDB materialization, and GC-content computation.
- `materialized_seq`: materialize BAM `SEQ` strings once, then time only
  the real BAM sequence strings flowing through `seq_gc_content(...)`.
  This better isolates the SIMD byte-scanning kernel on realistic
  sequence lengths.

The default local input is
`HG00106.chrom11.ILLUMINA.bwa.GBR.exome.20130415.bam` when present,
falling back to `test/data/nanopore.bam` for smaller checkouts.

# Run

``` sh
make bench-simd-bam-gc
```

Useful overrides:

- `SIMD_BAM_GC_BAM`: BAM path; default is the local HG00106 exome BAM if
  present
- `SIMD_BAM_GC_MAX_READS`: read limit before aggregation, default `0`
  for all reads
- `SIMD_BAM_GC_RUNS`: timed repeats per backend/workload, default `5`
- `SIMD_BAM_GC_THREADS`: DuckDB threads, default `1` to isolate SIMD
  effects
- `SIMD_BAM_GC_BACKENDS`: comma-separated backend requests, default
  `scalar,auto,avx2,avx512`
- `SIMD_BAM_GC_MODES`: comma-separated workloads, default
  `bam_scan,materialized_seq`
- `SIMD_BAM_GC_FORCE=1`: rerun even if cached JSON exists
- `DUCKHTS_EXTENSION`: extension path override

# Configuration

| setting          | value                                               |
|:-----------------|:----------------------------------------------------|
| bam              | HG00106.chrom11.ILLUMINA.bwa.GBR.exome.20130415.bam |
| max_reads        | all                                                 |
| iterations       | 5                                                   |
| threads          | 1                                                   |
| backend_requests | scalar,auto,avx2,avx512                             |
| workloads        | bam_scan,materialized_seq                           |
| extension        | build/release/duckhts.duckdb_extension              |

# Results

| benchmark        | backend_request | requested_backend | selected_backend | kernel_backend | status  |   reads | total_bases | load_sec | median_sec | min_sec | mbases_per_sec_median | speedup_vs_scalar | skip_reason                                      |
|:-----------------|:----------------|:------------------|:-----------------|:---------------|:--------|--------:|------------:|---------:|-----------:|--------:|----------------------:|------------------:|:-------------------------------------------------|
| bam_scan         | scalar          | scalar            | scalar           | scalar         | ok      | 3245905 |   327836405 |       NA |      1.148 |   1.136 |               285.545 |             1.000 | NA                                               |
| materialized_seq | scalar          | scalar            | scalar           | scalar         | ok      | 3245905 |   327836405 |    0.610 |      0.577 |   0.575 |               568.385 |             1.000 | NA                                               |
| bam_scan         | auto            | auto              | avx2             | avx2           | ok      | 3245905 |   327836405 |       NA |      0.601 |   0.595 |               545.205 |             1.909 | NA                                               |
| materialized_seq | auto            | auto              | avx2             | avx2           | ok      | 3245905 |   327836405 |    0.603 |      0.107 |   0.107 |              3065.541 |             5.393 | NA                                               |
| bam_scan         | avx2            | avx2              | avx2             | avx2           | ok      | 3245905 |   327836405 |       NA |      0.610 |   0.601 |               537.468 |             1.882 | NA                                               |
| materialized_seq | avx2            | avx2              | avx2             | avx2           | ok      | 3245905 |   327836405 |    0.599 |      0.111 |   0.110 |              2949.096 |             5.189 | NA                                               |
| bam_scan         | avx512          | NA                | NA               | NA             | skipped |      NA |          NA |       NA |         NA |      NA |                    NA |                NA | backend request is not available in this process |
| materialized_seq | avx512          | NA                | NA               | NA             | skipped |      NA |          NA |       NA |         NA |      NA |                    NA |                NA | backend request is not available in this process |

# Interpretation

For `bam_scan`, the fastest measured `seq_base_counts` kernel backend
was `avx2` selected from request `auto`, at 545.2 Mbases/s (1.91x
scalar).

For `materialized_seq`, the fastest measured `seq_base_counts` kernel
backend was `avx2` selected from request `auto`, at 3065.5 Mbases/s
(5.39x scalar).

`bam_scan` is the practical end-to-end view. `materialized_seq` removes
repeated BAM I/O and shows the SIMD kernel headroom on realistic read
strings.

# Raw timing vectors

| benchmark        | backend_request | selected_backend | kernel_backend | iteration |  seconds |
|:-----------------|:----------------|:-----------------|:---------------|----------:|---------:|
| bam_scan         | scalar          | scalar           | scalar         |         1 | 1.148109 |
| bam_scan         | scalar          | scalar           | scalar         |         2 | 1.143463 |
| bam_scan         | scalar          | scalar           | scalar         |         3 | 1.156094 |
| bam_scan         | scalar          | scalar           | scalar         |         4 | 1.149882 |
| bam_scan         | scalar          | scalar           | scalar         |         5 | 1.136231 |
| materialized_seq | scalar          | scalar           | scalar         |         1 | 0.577441 |
| materialized_seq | scalar          | scalar           | scalar         |         2 | 0.576786 |
| materialized_seq | scalar          | scalar           | scalar         |         3 | 0.579256 |
| materialized_seq | scalar          | scalar           | scalar         |         4 | 0.575871 |
| materialized_seq | scalar          | scalar           | scalar         |         5 | 0.574703 |
| bam_scan         | auto            | avx2             | avx2           |         1 | 0.609519 |
| bam_scan         | auto            | avx2             | avx2           |         2 | 0.594682 |
| bam_scan         | auto            | avx2             | avx2           |         3 | 0.601308 |
| bam_scan         | auto            | avx2             | avx2           |         4 | 0.597413 |
| bam_scan         | auto            | avx2             | avx2           |         5 | 0.606577 |
| materialized_seq | auto            | avx2             | avx2           |         1 | 0.106942 |
| materialized_seq | auto            | avx2             | avx2           |         2 | 0.106646 |
| materialized_seq | auto            | avx2             | avx2           |         3 | 0.107784 |
| materialized_seq | auto            | avx2             | avx2           |         4 | 0.106550 |
| materialized_seq | auto            | avx2             | avx2           |         5 | 0.106989 |
| bam_scan         | avx2            | avx2             | avx2           |         1 | 0.609964 |
| bam_scan         | avx2            | avx2             | avx2           |         2 | 0.616871 |
| bam_scan         | avx2            | avx2             | avx2           |         3 | 0.613045 |
| bam_scan         | avx2            | avx2             | avx2           |         4 | 0.608880 |
| bam_scan         | avx2            | avx2             | avx2           |         5 | 0.601385 |
| materialized_seq | avx2            | avx2             | avx2           |         1 | 0.112241 |
| materialized_seq | avx2            | avx2             | avx2           |         2 | 0.110804 |
| materialized_seq | avx2            | avx2             | avx2           |         3 | 0.111500 |
| materialized_seq | avx2            | avx2             | avx2           |         4 | 0.110436 |
| materialized_seq | avx2            | avx2             | avx2           |         5 | 0.111165 |
