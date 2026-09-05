Automatic BAM full-file scan: small-fixture cost
================

This compares automatic BAM scans at two named reader revisions, using
the existing `benchmark_simd_bam_gc.py` workload on its committed
nanopore fixture. It exercises many empty contig claims and full `SEQ`
decoding, but this 186-record fixture is dominated by setup cost: it is
not a WGS throughput benchmark. There are no unplaced reads in this
performance input, so both builds do equal work. The separate SQL/R
regression corpus tests retained no-coordinate records, physical
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
`BAM_SCAN_BASELINE_EXTENSION` to the baseline binary,
`BAM_SCAN_BASELINE_REVISION` to its revision, then render this report
from the repository root. `DUCKHTS_EXTENSION` selects the candidate
binary; `PYTHON` selects the Python environment with DuckDB installed.

| property             | value                                                                      |
|:---------------------|:---------------------------------------------------------------------------|
| baseline revision    | da5daa2cdd75aa9a920bf6a71f18174661b15198                                   |
| candidate revision   | 05fcdfdfcb32eb94a987d946e65d377040644c1a                                   |
| candidate src tree   | b7e5c3bbd984e0e339780d71d6324098622d2fa5                                   |
| baseline binary MD5  | 5bfcfbbd4d206288713427cb9d402806                                           |
| candidate binary MD5 | ec990d55a894db4fd760d84b7c66026b                                           |
| input                | test/data/nanopore.bam                                                     |
| input MD5            | 850dc34ada8d7023ee7146c7953da90b                                           |
| input bytes          | 283081                                                                     |
| HTSlib               | 1.24                                                                       |
| host                 | Linux                                                                      |
| CPU affinity         | pid 595985’s current affinity list: 0,2,4,6                                |
| repetitions          | 9 per build/thread setting after warm-up; fresh Python process per setting |

| build     | duckdb_workers | htslib_workers_per_handle | input_records | decoded_sequence_rows | input_bases | sql_output_rows |   gc_sum | backend | median_seconds | min_seconds | max_seconds |
|:----------|---------------:|--------------------------:|--------------:|----------------------:|------------:|----------------:|---------:|:--------|---------------:|------------:|------------:|
| baseline  |              1 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001684 |    0.001316 |    0.012549 |
| candidate |              1 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001561 |    0.001486 |    0.001701 |
| baseline  |              4 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.002069 |    0.001731 |    0.002444 |
| candidate |              4 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.002239 |    0.001794 |    0.004775 |

Each query scans 186 BAM records, materializes their sequences into
DuckDB vectors, and returns one aggregate row containing the read count,
total bases, and sum of GC fractions. All nine repetitions validate
those values; equality here is not a substitute for the independent
SAM-field and row-multiset regression tests. No claim of general
performance parity or speedup follows from these small timings.
