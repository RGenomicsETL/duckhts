Indexed BCF/VCF scan initialization: small-fixture cost
================

This paired measurement exercises automatic full-file scan planning and
complete materialization on committed BCF/VCF fixtures. The original
two-record inputs have one occupied contig and 256 empty header contigs.
A second five-record corpus has two occupied contigs, two samples,
equal-position distinct alleles, and a physical duplicate; its VCF
header declares all, some, or none of the contigs. It uses non-colocated
indexes. Timings measure setup and small-result materialization, not
cohort throughput.

VCF work planning now uses the Tabix dictionary. The one-occupied-contig
VCF therefore selects a single full-file stream instead of 257 indexed
contig tasks; the five-record VCFs use two contig tasks irrespective of
header completeness. BCF retains its header dictionary. This is a
deliberate scheduling change, not a promise that the requested DuckDB
thread limit equals the number of active workers.

The nearest larger full-column report is
[`benchmark_reader_allocations.md`](benchmark_reader_allocations.md):
4,048,342 GIAB HG002 records in wide/tidy form, one sequential scan
handle and zero HTSlib decompression workers. It does not exercise
automatic parallel index reloads and is not an identical-workload
baseline. The previous identical two-record workload is retained at
[`6d6b5d2:benchmark_bcf_scan_init.md`](https://github.com/RGenomicsETL/duckhts/blob/6d6b5d2ad1ce7ef8640fd2fda152b930a8ff14a0/benchmarks/benchmark_bcf_scan_init.md).
This report measures the pre-fix binary alongside the candidate on the
same host.

## Reproduction and identity

Set `BCF_SCAN_BASELINE_EXTENSION` and `BCF_SCAN_BASELINE_REVISION` to a
separately built baseline. `DUCKHTS_EXTENSION` selects the candidate.
From the repository root:

``` sh
taskset -c 0,2,4,6 Rscript -e 'rmarkdown::render("benchmarks/benchmark_bcf_scan_init.Rmd")'
```

| property             | value                                                                       |
|:---------------------|:----------------------------------------------------------------------------|
| baseline revision    | 76eba5b1ed3a48207d1aadb8ab4c8dc723d382cb                                    |
| candidate revision   | 99a8f380260c73692e929f48c932e34dfc624462                                    |
| candidate src tree   | 74f4dc2fe246e595b7976452e17f17bbcafdac5a                                    |
| baseline binary MD5  | 0d03c255a6609e1b67a36672c0b4168c                                            |
| candidate binary MD5 | 9b96299db0b2f7b4625a6704cca982a9                                            |
| R                    | R version 4.6.0 (2026-04-24)                                                |
| DuckDB R package     | 1.5.3                                                                       |
| HTSlib               | 1.24                                                                        |
| CPU affinity         | pid 1052870’s current affinity list: 0,2,4,6                                |
| repetitions          | 9 timed calls after warm-up; fresh R process per build/input/thread setting |

| file                                                | bytes | md5                              |
|:----------------------------------------------------|------:|:---------------------------------|
| test/data/parallel_empty_contigs.bcf                |  1351 | e7576013a2a4cfc843130eff617ea527 |
| test/data/parallel_empty_contigs.vcf.gz             |   372 | 8583a37cb25fe73d7dc8925652804942 |
| test/data/bcf_scan_contigs.full.vcf.gz              |   424 | abd4cc644200e4812df5ccbf185e3a4f |
| test/data/bcf_scan_contigs.partial.vcf.gz           |   421 | 4dcaa2047c179520b1259d7671c8a0e9 |
| test/data/bcf_scan_contigs.none.vcf.gz              |   404 | 01c9b6bbc25b859aebf67f74d9e265c3 |
| test/data/parallel_empty_contigs.bcf.csi            |   117 | 1cb71253d787009c2867430c349f2769 |
| test/data/parallel_empty_contigs.vcf.gz.tbi         |   106 | 9f978fe03cf1a2d047683b8c2e578eae |
| test/data/bcf_scan_contigs.full.vcf.gz.index.tbi    |   126 | fbb48c2c7cd524e891c35b7bc11f6401 |
| test/data/bcf_scan_contigs.partial.vcf.gz.index.tbi |   126 | 953769f59bf8df49d7ea3ee4de9c4bb7 |
| test/data/bcf_scan_contigs.none.vcf.gz.index.tbi    |   126 | 698f40d722cfee51932a386cb0b2f937 |
| test/data/bcf_scan_contigs.bcf                      |   470 | 8e006ec9582009fc56a1616e5dc67010 |

| build     | format     | header      | duckdb_threads | htslib_workers_per_handle | header_contigs | input_records | materialized_rows | materialized_columns | exact_rows | median_ms | min_ms | max_ms |
|:----------|:-----------|:------------|---------------:|--------------------------:|---------------:|--------------:|------------------:|---------------------:|:-----------|----------:|-------:|-------:|
| baseline  | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                 2 |                    8 | TRUE       |     1.173 |  1.132 |  1.470 |
| candidate | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                 2 |                    8 | TRUE       |     1.197 |  1.158 |  1.397 |
| baseline  | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                 2 |                    8 | TRUE       |     1.202 |  1.178 |  1.485 |
| candidate | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                 2 |                    8 | TRUE       |     1.320 |  1.231 |  1.480 |
| baseline  | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                 2 |                    8 | TRUE       |     0.891 |  0.835 |  0.970 |
| candidate | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                 2 |                    8 | TRUE       |     0.885 |  0.866 |  0.921 |
| baseline  | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                 2 |                    8 | TRUE       |     0.912 |  0.882 |  0.942 |
| candidate | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                 2 |                    8 | TRUE       |     0.887 |  0.835 |  0.965 |
| baseline  | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                 5 |                   10 | TRUE       |     0.919 |  0.873 |  0.973 |
| candidate | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                 5 |                   10 | TRUE       |     0.931 |  0.905 |  1.034 |
| baseline  | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                 5 |                   10 | TRUE       |     0.959 |  0.910 |  1.059 |
| candidate | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                 5 |                   10 | TRUE       |     0.931 |  0.911 |  1.031 |
| baseline  | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                 1 |                   10 | FALSE      |        NA |     NA |     NA |
| candidate | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                 5 |                   10 | TRUE       |     0.952 |  0.919 |  0.983 |
| baseline  | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                 1 |                   10 | FALSE      |        NA |     NA |     NA |
| candidate | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                 5 |                   10 | TRUE       |     0.959 |  0.933 |  1.074 |
| baseline  | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                 5 |                   10 | TRUE       |     0.944 |  0.916 |  0.995 |
| candidate | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                 5 |                   10 | TRUE       |     0.956 |  0.923 |  0.979 |
| baseline  | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                 5 |                   10 | TRUE       |     1.055 |  0.963 |  1.479 |
| candidate | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                 5 |                   10 | TRUE       |     0.971 |  0.945 |  1.168 |

Each timed `CREATE OR REPLACE TEMP TABLE AS SELECT *` materializes all
columns, not merely a count. The output denominator is two or five
complete rows, except the defective baseline partial-header scan, which
emits only one of five rows. Its missing timings are intentional: an
incomplete result is not comparable. Elapsed `Sys.time()` includes the
DBI call; connection creation, warm-up, and complete typed-row
verification are outside the timer. Every candidate result and every
valid baseline agrees exactly with the complete BCF sequential rows,
including GT, ALT and physical duplicate multiplicity, without asserting
scan arrival order. The BCF oracle also agrees across builds.

This does not measure the latency of an index failure, remote I/O,
large-cohort throughput, RSS, or maximum worker concurrency. Thread
settings are DuckDB’s requested limits. Small timings can vary with
scheduling; no general speedup or no-regression claim follows from them.
