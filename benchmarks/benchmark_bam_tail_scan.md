Automatic BAM full-file scan: small-fixture cost
================

This compares automatic BAM scans at two named reader revisions, using
the existing `benchmark_simd_bam_gc.py` workload on its committed
nanopore fixture, with and without projecting `FILE_OFFSET` via an
`IS NOT NULL` filter. It exercises many empty contig claims and full
`SEQ` decoding, but this 186-record fixture is dominated by setup cost:
it is not a WGS throughput benchmark. There are no unplaced reads in
this performance input, so both builds do equal work. The separate SQL/R
regression corpus tests retained no-coordinate records, physical
duplicates, BAI/CSI, CRAM, and output chunks.

The nearest identical checked-in workload is this report’s previous
rendering at `8e743aceb4109b5894ff23cf8f5c8e1b9a1607a1`: source
`05fcdfdfcb32eb94a987d946e65d377040644c1a`, 186 records / 249,110 bases,
1 or 4 DuckDB workers and 2 HTSlib workers per handle. Its candidate
medians without `FILE_OFFSET` were 1.561 ms and 2.239 ms respectively.
The unchanged baseline binary below contains that same reader source.
This is the first rendered measurement requiring `FILE_OFFSET`; its
same-run baseline uses the pre-fix binary, not an invented historical
offset workload.

The older larger-input report is
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
| baseline revision    | 8e743aceb4109b5894ff23cf8f5c8e1b9a1607a1                                   |
| candidate revision   | 8b057efe41f865831b00b217ede522c61137fdc0                                   |
| candidate src tree   | 13ea540bb2cf2e03664b0faef6aa84a910f89b93                                   |
| baseline binary MD5  | ec990d55a894db4fd760d84b7c66026b                                           |
| candidate binary MD5 | 0eed0f04a038542f0c97b130c0fc928c                                           |
| input                | test/data/nanopore.bam                                                     |
| input MD5            | 850dc34ada8d7023ee7146c7953da90b                                           |
| input bytes          | 283081                                                                     |
| HTSlib               | 1.24                                                                       |
| host                 | Linux                                                                      |
| CPU affinity         | pid 715047’s current affinity list: 0,2,4,6                                |
| repetitions          | 9 per build/thread setting after warm-up; fresh Python process per setting |

| build     | workload        | duckdb_workers | htslib_workers_per_handle | input_records | decoded_sequence_rows | input_bases | sql_output_rows |   gc_sum | backend | median_seconds | min_seconds | max_seconds |
|:----------|:----------------|---------------:|--------------------------:|--------------:|----------------------:|------------:|----------------:|---------:|:--------|---------------:|------------:|------------:|
| baseline  | bam_scan        |              1 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001272 |    0.001248 |    0.001407 |
| baseline  | bam_scan_offset |              1 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001364 |    0.001293 |    0.001572 |
| candidate | bam_scan        |              1 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001296 |    0.001241 |    0.001437 |
| candidate | bam_scan_offset |              1 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001314 |    0.001257 |    0.001385 |
| baseline  | bam_scan        |              4 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001929 |    0.001759 |    0.005324 |
| baseline  | bam_scan_offset |              4 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001785 |    0.001526 |    0.002132 |
| candidate | bam_scan        |              4 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001857 |    0.001605 |    0.005024 |
| candidate | bam_scan_offset |              4 |                         2 |           186 |                   186 |      249110 |               1 | 75.87703 | avx2    |       0.001857 |    0.001625 |    0.002067 |

Each query scans 186 BAM records, materializes their sequences into
DuckDB vectors, and returns one aggregate row containing the read count,
total bases, and sum of GC fractions. All nine repetitions validate
those values; equality here is not a substitute for the independent
SAM-field and row-multiset regression tests. The `bam_scan_offset` query
additionally materializes and filters `FILE_OFFSET`; all 186 BGZF BAM
records must still reach the aggregate. Exact offset values are tested
separately against physical BGZF/BAM record ends in SQL and installed R.
No claim of general performance parity or speedup follows from these
small timings.
