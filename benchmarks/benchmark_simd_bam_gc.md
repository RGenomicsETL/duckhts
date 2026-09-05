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

The default input is the committed `test/data/nanopore.bam` fixture. Set
`SIMD_BAM_GC_BAM` to benchmark a separately staged BAM under the DuckHTS
cache. The shared R driver runs each backend in a fresh R/DuckDB
process, checks every repeat and cross-backend aggregate, and measures
elapsed `Sys.time()` including DBI result retrieval. Warm-up and
connection creation are outside the timer. Historical Python-host
measurements are not directly comparable to this R host.

# Run

``` sh
make bench-simd-bam-gc
```

Useful overrides:

- `SIMD_BAM_GC_BAM`: BAM path; default is the committed nanopore fixture
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

| setting          | value                                    |
|:-----------------|:-----------------------------------------|
| bam              | test/data/nanopore.bam                   |
| max_reads        | all                                      |
| iterations       | 5                                        |
| threads          | 1                                        |
| backend_requests | scalar,auto,avx2,avx512                  |
| workloads        | bam_scan,materialized_seq                |
| extension        | build/release/duckhts.duckdb_extension   |
| source revision  | a1509946e9f7d425ce9af755b3dad4d41f03dba3 |
| source tree      | 346d710310f00739fa0361f0119a23cc69fc40e8 |
| input MD5        | 850dc34ada8d7023ee7146c7953da90b         |
| extension MD5    | cb5a8d25a505cb50ed5fadf496040190         |
| R                | R version 4.6.0 (2026-04-24)             |
| DuckDB R package | 1.5.3                                    |

# Results

| benchmark        | backend_request | requested_backend | selected_backend | kernel_backend | status  | reads | total_bases | load_sec | median_sec | min_sec | mbases_per_sec_median | speedup_vs_scalar | skip_reason                                      |
|:-----------------|:----------------|:------------------|:-----------------|:---------------|:--------|------:|------------:|---------:|-----------:|--------:|----------------------:|------------------:|:-------------------------------------------------|
| bam_scan         | scalar          | scalar            | scalar           | scalar         | ok      |   186 |      249110 |       NA |      0.003 |   0.002 |                97.458 |             1.000 | NA                                               |
| materialized_seq | scalar          | scalar            | scalar           | scalar         | ok      |   186 |      249110 |    0.002 |      0.002 |   0.001 |               164.906 |             1.000 | NA                                               |
| bam_scan         | auto            | auto              | avx2             | avx2           | ok      |   186 |      249110 |       NA |      0.002 |   0.002 |               142.953 |             1.467 | NA                                               |
| materialized_seq | auto            | auto              | avx2             | avx2           | ok      |   186 |      249110 |    0.002 |      0.001 |   0.001 |               371.698 |             2.254 | NA                                               |
| bam_scan         | avx2            | avx2              | avx2             | avx2           | ok      |   186 |      249110 |       NA |      0.002 |   0.002 |               139.722 |             1.434 | NA                                               |
| materialized_seq | avx2            | avx2              | avx2             | avx2           | ok      |   186 |      249110 |    0.002 |      0.001 |   0.001 |               393.538 |             2.386 | NA                                               |
| bam_scan         | avx512          | NA                | NA               | NA             | skipped |    NA |          NA |       NA |         NA |      NA |                    NA |                NA | backend request is not available in this process |
| materialized_seq | avx512          | NA                | NA               | NA             | skipped |    NA |          NA |       NA |         NA |      NA |                    NA |                NA | backend request is not available in this process |

# Interpretation

For `bam_scan`, the fastest measured `seq_base_counts` kernel backend
was `avx2` selected from request `auto`, at 143.0 Mbases/s (1.47x
scalar).

For `materialized_seq`, the fastest measured `seq_base_counts` kernel
backend was `avx2` selected from request `avx2`, at 393.5 Mbases/s
(2.39x scalar).

`bam_scan` is the practical end-to-end view. `materialized_seq` removes
repeated BAM I/O and shows the SIMD kernel headroom on realistic read
strings.

# Raw timing vectors

| benchmark        | backend_request | selected_backend | kernel_backend | iteration |  seconds |
|:-----------------|:----------------|:-----------------|:---------------|----------:|---------:|
| bam_scan         | scalar          | scalar           | scalar         |         1 | 0.002684 |
| bam_scan         | scalar          | scalar           | scalar         |         2 | 0.002556 |
| bam_scan         | scalar          | scalar           | scalar         |         3 | 0.002474 |
| bam_scan         | scalar          | scalar           | scalar         |         4 | 0.002501 |
| bam_scan         | scalar          | scalar           | scalar         |         5 | 0.002772 |
| materialized_seq | scalar          | scalar           | scalar         |         1 | 0.001462 |
| materialized_seq | scalar          | scalar           | scalar         |         2 | 0.001511 |
| materialized_seq | scalar          | scalar           | scalar         |         3 | 0.001477 |
| materialized_seq | scalar          | scalar           | scalar         |         4 | 0.001597 |
| materialized_seq | scalar          | scalar           | scalar         |         5 | 0.001565 |
| bam_scan         | auto            | avx2             | avx2           |         1 | 0.002084 |
| bam_scan         | auto            | avx2             | avx2           |         2 | 0.001743 |
| bam_scan         | auto            | avx2             | avx2           |         3 | 0.001725 |
| bam_scan         | auto            | avx2             | avx2           |         4 | 0.001720 |
| bam_scan         | auto            | avx2             | avx2           |         5 | 0.001798 |
| materialized_seq | auto            | avx2             | avx2           |         1 | 0.000642 |
| materialized_seq | auto            | avx2             | avx2           |         2 | 0.000670 |
| materialized_seq | auto            | avx2             | avx2           |         3 | 0.000646 |
| materialized_seq | auto            | avx2             | avx2           |         4 | 0.000715 |
| materialized_seq | auto            | avx2             | avx2           |         5 | 0.000676 |
| bam_scan         | avx2            | avx2             | avx2           |         1 | 0.001786 |
| bam_scan         | avx2            | avx2             | avx2           |         2 | 0.001643 |
| bam_scan         | avx2            | avx2             | avx2           |         3 | 0.001567 |
| bam_scan         | avx2            | avx2             | avx2           |         4 | 0.001874 |
| bam_scan         | avx2            | avx2             | avx2           |         5 | 0.001783 |
| materialized_seq | avx2            | avx2             | avx2           |         1 | 0.000622 |
| materialized_seq | avx2            | avx2             | avx2           |         2 | 0.000624 |
| materialized_seq | avx2            | avx2             | avx2           |         3 | 0.000641 |
| materialized_seq | avx2            | avx2             | avx2           |         4 | 0.000657 |
| materialized_seq | avx2            | avx2             | avx2           |         5 | 0.000633 |
