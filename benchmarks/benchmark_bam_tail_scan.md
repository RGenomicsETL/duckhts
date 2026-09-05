Automatic BAM full-file scan: small-fixture cost
================

This compares automatic BAM scans at two named reader revisions, using
the shared `benchmark_simd_bam_gc.R` driver on its committed nanopore
fixture, with and without projecting `FILE_OFFSET` via an `IS NOT NULL`
filter. It also times full materialization with standard/AUX tags and
text or numeric CIGAR/SEQ/QUAL, covering formatter scratch and DuckDB
list writes. It exercises many empty contig claims and full `SEQ`
decoding, but this 186-record fixture is dominated by setup cost: it is
not a WGS throughput benchmark. There are no unplaced reads in this
performance input, so both builds do equal work. The separate SQL/R
regression corpus tests retained no-coordinate records, physical
duplicates, BAI/CSI, CRAM, and output chunks.

The nearest checked-in SEQ/GC workload is this report’s previous
rendering at `e87b9f3ee5e27068964e4e57e107bacd7b5def3e`: source
`8b057efe41f865831b00b217ede522c61137fdc0`, 186 records / 249,110 bases,
1 or 4 DuckDB workers and 2 HTSlib workers per handle. Its candidate
medians without `FILE_OFFSET` were 1.296 ms and 1.857 ms respectively;
with `FILE_OFFSET`, 1.314 ms and 1.857 ms. That run used Python/DuckDB;
its frontend timings are not directly comparable to this R/DBI run. Both
binaries below are measured with the same R driver. The unchanged
baseline binary below contains that same reader source. This is the
first rendered all-column text/numeric materialization measurement; its
same-run baseline uses the pre-fix binary, not an invented historical
workload.

The older larger-input report is [`benchmark_simd_bam_gc.md` at the
baseline
revision](https://github.com/RGenomicsETL/duckhts/blob/e87b9f3ee5e27068964e4e57e107bacd7b5def3e/benchmarks/benchmark_simd_bam_gc.md),
whose historical run used 3,245,905 HG00106 exome reads / 327,836,405
bases and one DuckDB worker. That input is not staged here; its timings
are not a comparable regression baseline. The same-run comparison below
uses the current driver’s committed default input. Remote I/O,
large-input throughput, and CRAM performance remain unmeasured.

## Reproduction and identity

Build each named source revision separately. Set
`BAM_SCAN_BASELINE_EXTENSION` to the baseline binary,
`BAM_SCAN_BASELINE_REVISION` to its revision, then render this report
from the repository root. `DUCKHTS_EXTENSION` selects the candidate
binary. The R driver uses `callr`, `duckdb`, `DBI`, and `digest`.

| property             | value                                                                 |
|:---------------------|:----------------------------------------------------------------------|
| baseline revision    | e87b9f3ee5e27068964e4e57e107bacd7b5def3e                              |
| candidate revision   | a1509946e9f7d425ce9af755b3dad4d41f03dba3                              |
| candidate src tree   | 346d710310f00739fa0361f0119a23cc69fc40e8                              |
| baseline binary MD5  | 0eed0f04a038542f0c97b130c0fc928c                                      |
| candidate binary MD5 | cb5a8d25a505cb50ed5fadf496040190                                      |
| input                | test/data/nanopore.bam                                                |
| input MD5            | 850dc34ada8d7023ee7146c7953da90b                                      |
| input bytes          | 283081                                                                |
| HTSlib               | 1.24                                                                  |
| host                 | Linux                                                                 |
| R                    | R version 4.6.0 (2026-04-24)                                          |
| DuckDB R package     | 1.5.3                                                                 |
| CPU affinity         | pid 847456’s current affinity list: 0,2,4,6                           |
| repetitions          | 9 per build/thread setting after warm-up; fresh R process per setting |

| build     | workload                | duckdb_workers | htslib_workers_per_handle | input_records | decoded_sequence_rows | input_bases | sql_output_rows |   gc_sum | backend | median_seconds | min_seconds | max_seconds |
|:----------|:------------------------|---------------:|--------------------------:|--------------:|----------------------:|------------:|----------------:|---------:|:--------|---------------:|------------:|------------:|
| baseline  | bam_scan                |              1 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001746 |    0.001666 |    0.001945 |
| baseline  | bam_scan_offset         |              1 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001832 |    0.001709 |    0.001913 |
| baseline  | bam_materialize_text    |              1 |                         2 |           186 |                   186 |      249110 |             186 | 75.87703 | avx2    |       0.003503 |    0.003134 |    0.003731 |
| baseline  | bam_materialize_numeric |              1 |                         2 |           186 |                   186 |      249110 |             186 | 75.87703 | avx2    |       0.004378 |    0.003674 |    0.006364 |
| candidate | bam_scan                |              1 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001779 |    0.001629 |    0.001972 |
| candidate | bam_scan_offset         |              1 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001712 |    0.001671 |    0.001857 |
| candidate | bam_materialize_text    |              1 |                         2 |           186 |                   186 |      249110 |             186 | 75.87703 | avx2    |       0.003431 |    0.003262 |    0.004255 |
| candidate | bam_materialize_numeric |              1 |                         2 |           186 |                   186 |      249110 |             186 | 75.87703 | avx2    |       0.003886 |    0.003611 |    0.004093 |
| baseline  | bam_scan                |              4 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.002269 |    0.002139 |    0.002527 |
| baseline  | bam_scan_offset         |              4 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.002182 |    0.002016 |    0.005234 |
| baseline  | bam_materialize_text    |              4 |                         2 |           186 |                   186 |      249110 |             186 | 75.87703 | avx2    |       0.003378 |    0.003242 |    0.003917 |
| baseline  | bam_materialize_numeric |              4 |                         2 |           186 |                   186 |      249110 |             186 | 75.87703 | avx2    |       0.003926 |    0.003495 |    0.004133 |
| candidate | bam_scan                |              4 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.002418 |    0.002175 |    0.003792 |
| candidate | bam_scan_offset         |              4 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.002343 |    0.002171 |    0.003561 |
| candidate | bam_materialize_text    |              4 |                         2 |           186 |                   186 |      249110 |             186 | 75.87703 | avx2    |       0.003508 |    0.003175 |    0.007134 |
| candidate | bam_materialize_numeric |              4 |                         2 |           186 |                   186 |      249110 |             186 | 75.87703 | avx2    |       0.004000 |    0.003764 |    0.004395 |

|     | workload                | row_multiset_sha256                                              |
|:----|:------------------------|:-----------------------------------------------------------------|
| 3   | bam_materialize_text    | 476467e0e9b89afb79c0622121e283462fd3c50b6049f191b2b82bb179ae74a8 |
| 4   | bam_materialize_numeric | f8f2594075c3ba5412bd9fa67f7544f941f4842422fa92b37e35ff252d7cba08 |

Each SEQ/GC query scans 186 BAM records, materializes their sequences
into DuckDB vectors, and returns one aggregate row containing the read
count, total bases, and sum of GC fractions. All nine repetitions
validate those values; equality here is not a substitute for the
independent SAM-field and row-multiset regression tests. The
`bam_scan_offset` query additionally materializes and filters
`FILE_OFFSET`; all 186 BGZF BAM records must still reach the aggregate.
Exact offset values are tested separately against physical BGZF/BAM
record ends in SQL and installed R. The full-materialization queries
time `CREATE OR REPLACE TEMP TABLE AS SELECT *` with standard/AUX tags
enabled. Their output denominator is 186 materialized table rows, not
the one-row CTAS completion response. Both representations include
249,110 sequence and quality elements; numeric mode also writes packed
CIGAR lists. Outside each timed call, the driver checks the aggregate
and a SHA-256 digest of the complete typed row multiset (including every
projected tag), sorted by SQL `ORDER BY ALL` and serialized by R. Each
digest agrees across revisions and thread settings; text and numeric
schemas intentionally differ. This is a small-input latency measurement,
not a materialized-byte or RSS benchmark. The adversarial
formatter/fault fixtures are separate correctness tests. No claim of
general performance parity or speedup follows from these small timings.
Elapsed wall time uses R’s `Sys.time()` and includes the DBI call/result
overhead; connection creation, warm-up, and post-materialization
verification are not timed.
