Automatic BAM tail ownership: small-fixture scan cost
================

This compares the automatic BAM scan before and after the no-coordinate
tail fix, using the existing `benchmark_simd_bam_gc.py` workload on its
committed nanopore fixture. It exercises many empty contig claims and
full `SEQ` decoding, but this 186-record fixture is dominated by setup
cost: it is not a WGS throughput benchmark. There are no unplaced reads
in this performance input, so both builds do equal work. The separate
SQL/R regression corpus tests retained no-coordinate records, physical
duplicates, BAI/CSI, CRAM, and output chunks.

The nearest checked-in report is
[`benchmark_simd_bam_gc.md`](benchmark_simd_bam_gc.md), whose historical
run used 3,245,905 HG00106 exome reads / 327,836,405 bases and one
DuckDB worker. That input is not staged here; its timings are not a
comparable regression baseline. The same-run comparison below uses the
current driver’s committed default input. Remote I/O, large-input
throughput, and CRAM performance remain unmeasured.

## Reproduction and identity

Build each named source revision separately. Set
`BAM_SCAN_BASELINE_EXTENSION` to the pre-fix binary,
`BAM_SCAN_BASELINE_REVISION` to its revision, then render this report
from the repository root. `DUCKHTS_EXTENSION` selects the candidate
binary; `PYTHON` selects the Python environment with DuckDB installed.

| property             | value                                                                      |
|:---------------------|:---------------------------------------------------------------------------|
| baseline revision    | 3878b976380bd4d3816d4fc416a78b59eaf0580f                                   |
| candidate revision   | 67db30ac9dc40d0de4f6bb8ae55c6cfa3ea48011                                   |
| candidate src tree   | 352c4c79b2eb84b5ef7847e316420b8f4a0a238c                                   |
| baseline binary MD5  | 92aef5253e9920fe0b59df156412991d                                           |
| candidate binary MD5 | df999b32c0759ea8a9fa50b4ae0cecd7                                           |
| input                | test/data/nanopore.bam                                                     |
| input MD5            | 850dc34ada8d7023ee7146c7953da90b                                           |
| input bytes          | 283081                                                                     |
| HTSlib               | 1.24                                                                       |
| host                 | Linux                                                                      |
| CPU affinity         | pid 366859’s current affinity list: 0,2,4,6                                |
| repetitions          | 9 per build/thread setting after warm-up; fresh Python process per setting |

| build     | duckdb_workers | htslib_workers_per_handle | input_records | decoded_sequence_rows | input_bases | sql_output_rows |   gc_sum | backend | median_seconds | min_seconds | max_seconds |
|:----------|---------------:|--------------------------:|--------------:|----------------------:|------------:|----------------:|---------:|:--------|---------------:|------------:|------------:|
| baseline  |              1 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001314 |    0.001199 |    0.001399 |
| candidate |              1 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001340 |    0.001274 |    0.001407 |
| baseline  |              4 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001739 |    0.001664 |    0.002049 |
| candidate |              4 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001810 |    0.001654 |    0.005563 |

Each query scans 186 BAM records, materializes their sequences into
DuckDB vectors, and returns one aggregate row containing the read count,
total bases, and sum of GC fractions. All nine repetitions validate
those values; equality here is not a substitute for the independent
SAM-field and row-multiset regression tests. No claim of general
performance parity or speedup follows from these small timings.
