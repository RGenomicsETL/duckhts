Indexed BCF/VCF scan initialization: small-fixture cost
================

This paired measurement exercises automatic full-file scan planning and
complete materialization on committed BCF/VCF fixtures staged through
duckhtsbench’s registry. The registry pins encoded file bytes; staging
is a network-free verified copy, not fixture reconstruction. The
original two-record inputs have one occupied contig and 256 empty header
contigs. A second five-record corpus has two occupied contigs, two
samples, equal-position distinct alleles, and a physical duplicate; its
VCF header declares all, some, or none of the contigs. It uses
non-colocated indexes. Timings measure setup and small-result
materialization, not cohort throughput.

Both builds plan VCF work using the Tabix dictionary. The
one-occupied-contig VCF therefore selects a single full-file stream
instead of 257 indexed contig tasks; the five-record VCFs use two contig
tasks irrespective of header completeness. BCF retains its header
dictionary. This is a deliberate scheduling change, not a promise that
the requested DuckDB thread limit equals the number of active workers.
The candidate retains the parsed index at bind and removes worker index
reloads; the baseline reloads an independent index in each indexed scan
worker. Data and indexes remain unchanged during this timing workload.

The nearest larger full-column report is
[`benchmark_reader_allocations.md`](benchmark_reader_allocations.md):
4,048,342 GIAB HG002 records in wide/tidy form, one sequential scan
handle and zero HTSlib decompression workers. It does not exercise
automatic parallel index reloads and is not an identical-workload
baseline. The previous identical two/five-record workload is retained at
[`40c2b52:benchmark_bcf_scan_init.md`](https://github.com/RGenomicsETL/duckhts/blob/40c2b52793740064956709e44b825553efe76cd3/benchmarks/benchmark_bcf_scan_init.md).
This report measures the pre-fix binary alongside the candidate on the
same host.

## Reproduction and identity

Install duckhtsbench from this checkout. Set
`BCF_SCAN_BASELINE_EXTENSION` and `BCF_SCAN_BASELINE_REVISION` to a
separately built baseline. `DUCKHTS_EXTENSION` selects the candidate.
From the repository root:

``` sh
taskset -c 0,2,4,6 Rscript -e 'rmarkdown::render("benchmarks/benchmark_bcf_scan_init.Rmd")'
```

| property             | value                                                                       |
|:---------------------|:----------------------------------------------------------------------------|
| baseline revision    | d36b2a591cab42e9a000f9e71cfef058392c959d                                    |
| candidate revision   | c0e11568bd658cbb10e35ee60d2a5770afbde668                                    |
| candidate src tree   | f447adfb5ed8be98e33f85c2c9335b1e4f4e1e0e                                    |
| baseline binary MD5  | 9b96299db0b2f7b4625a6704cca982a9                                            |
| candidate binary MD5 | 37996465549069014fa3ae36bdb7589c                                            |
| R                    | R version 4.6.0 (2026-04-24)                                                |
| DuckDB R package     | 1.5.3                                                                       |
| HTSlib               | 1.24                                                                        |
| CPU affinity         | pid 1143669’s current affinity list: 0,2,4,6                                |
| repetitions          | 9 timed calls after warm-up; fresh R process per build/input/thread setting |

| file                                      | bytes | md5                              |
|:------------------------------------------|------:|:---------------------------------|
| parallel_empty_contigs.bcf                |  1351 | e7576013a2a4cfc843130eff617ea527 |
| parallel_empty_contigs.vcf.gz             |   372 | 8583a37cb25fe73d7dc8925652804942 |
| bcf_scan_contigs.full.vcf.gz              |   424 | abd4cc644200e4812df5ccbf185e3a4f |
| bcf_scan_contigs.partial.vcf.gz           |   421 | 4dcaa2047c179520b1259d7671c8a0e9 |
| bcf_scan_contigs.none.vcf.gz              |   404 | 01c9b6bbc25b859aebf67f74d9e265c3 |
| parallel_empty_contigs.bcf.csi            |   117 | 1cb71253d787009c2867430c349f2769 |
| parallel_empty_contigs.vcf.gz.tbi         |   106 | 9f978fe03cf1a2d047683b8c2e578eae |
| bcf_scan_contigs.full.vcf.gz.index.tbi    |   126 | fbb48c2c7cd524e891c35b7bc11f6401 |
| bcf_scan_contigs.partial.vcf.gz.index.tbi |   126 | 953769f59bf8df49d7ea3ee4de9c4bb7 |
| bcf_scan_contigs.none.vcf.gz.index.tbi    |   126 | 698f40d722cfee51932a386cb0b2f937 |
| bcf_scan_contigs.bcf                      |   470 | 8e006ec9582009fc56a1616e5dc67010 |

| build     | format     | header      | duckdb_threads | htslib_workers_per_handle | header_contigs | input_records | materialized_rows | materialized_columns | exact_rows | median_ms | min_ms | max_ms |
|:----------|:-----------|:------------|---------------:|--------------------------:|---------------:|--------------:|------------------:|---------------------:|:-----------|----------:|-------:|-------:|
| baseline  | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                 2 |                    8 | TRUE       |     1.179 |  1.137 |  1.336 |
| candidate | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                 2 |                    8 | TRUE       |     1.148 |  1.114 |  1.396 |
| baseline  | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                 2 |                    8 | TRUE       |     1.267 |  1.202 |  1.622 |
| candidate | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                 2 |                    8 | TRUE       |     1.303 |  1.162 |  1.487 |
| baseline  | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                 2 |                    8 | TRUE       |     0.904 |  0.858 |  0.951 |
| candidate | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                 2 |                    8 | TRUE       |     0.897 |  0.873 |  0.971 |
| baseline  | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                 2 |                    8 | TRUE       |     0.902 |  0.868 |  1.038 |
| candidate | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                 2 |                    8 | TRUE       |     0.884 |  0.850 |  0.940 |
| baseline  | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                 5 |                   10 | TRUE       |     0.901 |  0.879 |  1.003 |
| candidate | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                 5 |                   10 | TRUE       |     0.923 |  0.892 |  0.992 |
| baseline  | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                 5 |                   10 | TRUE       |     0.962 |  0.935 |  1.021 |
| candidate | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                 5 |                   10 | TRUE       |     0.957 |  0.919 |  0.999 |
| baseline  | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                 5 |                   10 | TRUE       |     0.980 |  0.931 |  1.044 |
| candidate | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                 5 |                   10 | TRUE       |     0.959 |  0.924 |  1.060 |
| baseline  | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                 5 |                   10 | TRUE       |     0.990 |  0.956 |  1.196 |
| candidate | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                 5 |                   10 | TRUE       |     1.012 |  0.973 |  1.061 |
| baseline  | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                 5 |                   10 | TRUE       |     0.980 |  0.969 |  1.066 |
| candidate | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                 5 |                   10 | TRUE       |     1.000 |  0.944 |  1.132 |
| baseline  | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                 5 |                   10 | TRUE       |     1.040 |  0.948 |  1.092 |
| candidate | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                 5 |                   10 | TRUE       |     1.005 |  0.969 |  1.095 |

Each timed `CREATE OR REPLACE TEMP TABLE AS SELECT *` materializes all
columns, not merely a count. The output denominator is two or five
complete rows on both builds, including partial-header scans fixed by
the preceding change. Elapsed `Sys.time()` includes the DBI call;
connection creation, warm-up, and complete typed-row verification are
outside the timer. Every candidate result and every valid baseline
agrees exactly with the complete BCF sequential rows, including GT, ALT
and physical duplicate multiplicity, without asserting scan arrival
order. The BCF oracle also agrees across builds.

This does not measure the latency of an index failure, remote I/O,
large-cohort throughput, RSS, or maximum worker concurrency. Thread
settings are DuckDB’s requested limits. Small timings can vary with
scheduling; no general speedup or no-regression claim follows from them.
