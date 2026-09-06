Indexed BCF/VCF scan initialization: small-fixture cost
================

This paired measurement exercises automatic contig scanning with a
usable worker index. Each committed fixture contains two records on one
occupied contig and 256 empty header contigs. It measures setup and
materialization latency, not cohort throughput. The separate SQL/R
regressions corrupt or remove a private index after preparing a query
and require an error instead of zero rows.

The nearest larger full-column report is
[`benchmark_reader_allocations.md`](benchmark_reader_allocations.md):
4,048,342 GIAB HG002 records in wide/tidy form, one sequential scan
handle and zero HTSlib decompression workers. It does not exercise
automatic parallel index reloads and is not an identical-workload
baseline. This first small-fixture scan-init report measures the pre-fix
binary alongside the candidate.

## Reproduction and identity

Set `BCF_SCAN_BASELINE_EXTENSION` and `BCF_SCAN_BASELINE_REVISION` to a
separately built baseline. `DUCKHTS_EXTENSION` selects the candidate.
From the repository root:

``` sh
taskset -c 0,2,4,6 Rscript -e 'rmarkdown::render("benchmarks/benchmark_bcf_scan_init.Rmd")'
```

| property             | value                                                                       |
|:---------------------|:----------------------------------------------------------------------------|
| baseline revision    | 536cf59425fe5be9aa4e2518feaadd226901146b                                    |
| candidate revision   | f06064bade5bd6dadcd6dceffb92669155a6d913                                    |
| candidate src tree   | e906cdc916567d4c995a52f2bf7d14267c0f2e48                                    |
| baseline binary MD5  | 5296739913e78217f12229647f5299a5                                            |
| candidate binary MD5 | 0d03c255a6609e1b67a36672c0b4168c                                            |
| R                    | R version 4.6.0 (2026-04-24)                                                |
| DuckDB R package     | 1.5.3                                                                       |
| HTSlib               | 1.24                                                                        |
| CPU affinity         | pid 970927’s current affinity list: 0,2,4,6                                 |
| repetitions          | 9 timed calls after warm-up; fresh R process per build/input/thread setting |

| file                                        | bytes | md5                              |
|:--------------------------------------------|------:|:---------------------------------|
| test/data/parallel_empty_contigs.bcf        |  1351 | e7576013a2a4cfc843130eff617ea527 |
| test/data/parallel_empty_contigs.vcf.gz     |   372 | 8583a37cb25fe73d7dc8925652804942 |
| test/data/parallel_empty_contigs.bcf.csi    |   117 | 1cb71253d787009c2867430c349f2769 |
| test/data/parallel_empty_contigs.vcf.gz.tbi |   106 | 9f978fe03cf1a2d047683b8c2e578eae |

| build     | format     | duckdb_threads | htslib_workers_per_handle | header_contigs | input_records | materialized_rows | materialized_columns | median_ms | min_ms | max_ms |
|:----------|:-----------|---------------:|--------------------------:|---------------:|--------------:|------------------:|---------------------:|----------:|-------:|-------:|
| baseline  | BCF/CSI    |              1 |                         0 |            257 |             2 |                 2 |                    8 |     1.203 |  1.162 |  1.338 |
| candidate | BCF/CSI    |              1 |                         0 |            257 |             2 |                 2 |                    8 |     1.161 |  1.138 |  1.296 |
| baseline  | BCF/CSI    |              4 |                         0 |            257 |             2 |                 2 |                    8 |     1.254 |  1.217 |  1.486 |
| candidate | BCF/CSI    |              4 |                         0 |            257 |             2 |                 2 |                    8 |     1.217 |  1.197 |  1.329 |
| baseline  | VCF.gz/TBI |              1 |                         0 |            257 |             2 |                 2 |                    8 |     0.952 |  0.884 |  1.076 |
| candidate | VCF.gz/TBI |              1 |                         0 |            257 |             2 |                 2 |                    8 |     0.902 |  0.890 |  0.978 |
| baseline  | VCF.gz/TBI |              4 |                         0 |            257 |             2 |                 2 |                    8 |     0.925 |  0.885 |  0.986 |
| candidate | VCF.gz/TBI |              4 |                         0 |            257 |             2 |                 2 |                    8 |     0.902 |  0.864 |  0.995 |

Each timed `CREATE OR REPLACE TEMP TABLE AS SELECT *` produces two
complete materialized rows, not merely a count. Elapsed `Sys.time()`
includes the DBI call; connection creation, warm-up, and complete
typed-row verification are outside the timer. All repetitions, formats,
thread settings and builds agree exactly on those rows, including GT and
ALT, without asserting scan arrival order.

This does not measure the latency of an index failure, remote I/O,
large-cohort throughput, RSS, or maximum worker concurrency. Thread
settings are DuckDB’s requested limits. Small timings can vary with
scheduling; no general speedup or no-regression claim follows from them.
