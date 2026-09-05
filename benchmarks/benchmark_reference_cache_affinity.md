DuckHTS worker-local reference cache
================

<!-- benchmark_reference_cache.md is generated from benchmark_reference_cache.Rmd. -->

This paired benchmark exercises `bcftools_norm_row` and
`bcftools_munge_row` against the full, uncompressed Broad hg38 v0
reference from the `duckhtsbench` registry (`liftover_grch38_fasta`).
The registered 3,249,912,778-byte size and MD5 are validated before use.
Input rows are synthetic SNVs at the first 1,000,000 A/C/G/T positions
starting at chr22:20,000,000; each ALT is C for reference A, otherwise
A. Each row requests one reference base and must return one unchanged,
correctly oriented row. Every timed run checks positions, alleles,
status/filter, row count, and the input-ID sum, not just throughput.

Run from the repo root after building both revisions with the same
toolchain:

``` sh
REFERENCE_CACHE_BASELINE_EXTENSION=/path/to/baseline/duckhts.duckdb_extension \
REFERENCE_CACHE_BASELINE_REV=<commit> \
Rscript -e 'rmarkdown::render("benchmarks/benchmark_reference_cache.Rmd", knit_root_dir=getwd())'
```

The candidate defaults to `build/release/duckhts.duckdb_extension`.
Stage the registered reference explicitly with
`duckhts_bench_fetch("liftover_grch38_fasta")` before rendering. No
benchmark download or alternate cache layout is embedded here. Linux is
required for peak RSS and the optional htslib interposition probe.

## Revisions and environment

    ##   revision                                   commit
    ##   baseline 47424232570311f09785b096f75e0fb083244734
    ##  candidate ae7c5242db2376ce3f958a269234aad4594cc9e0
    ##                     extension_md5
    ##  bcf15693b036abeccc30941db47821b1
    ##  e83dbba4a843ebd84125a71cee449b94

    ## R: R version 4.6.0 (2026-04-24)
    ## DuckDB: 1.5.3

    ## Linux 6.8.0-78-generic x86_64 GNU/Linux

    ## Model name:                           13th Gen Intel(R) Core(TM) i5-13500 BIOS Model name:                      13th Gen Intel(R) Core(TM) i5-13500 To Be Filled By O.E.M. CPU @ 2.4GHz

    ## Rows per query: 1000000 ; measured repetitions: 3

    ## Ubuntu clang version 18.1.3 (1ubuntu1)

    ## pid 3930354's current affinity list: 0,2,4,6

## Timings

Each function/thread/revision combination uses a fresh R process;
repetitions within that process reuse worker caches when DuckDB retains
workers. Table construction, FASTA staging, indexing, and R startup are
excluded from elapsed and CPU time. CPU is user + system process CPU
during the query. Peak RSS is the **whole-process high-water mark**,
including setup, not cache-only memory. The operating-system page cache
is warm; this is not a cold-storage benchmark.

    ##   revision operator threads median_elapsed_s median_cpu_s median_peak_rss_kib
    ##   baseline    munge       1            0.284        0.284              231080
    ##  candidate    munge       1            0.247        0.247              231116
    ##   baseline     norm       1            0.316        0.316              231364
    ##  candidate     norm       1            0.318        0.318              230984
    ##   baseline    munge       4            0.235        0.845              231684
    ##  candidate    munge       4            0.067        0.255              230900
    ##   baseline     norm       4            0.087        0.330              231196
    ##  candidate     norm       4            0.088        0.334              230860

    ##   revision operator threads run elapsed   cpu
    ##   baseline     norm       1   1   0.316 0.318
    ##   baseline     norm       1   2   0.316 0.315
    ##   baseline     norm       1   3   0.316 0.316
    ##  candidate     norm       1   1   0.314 0.313
    ##  candidate     norm       1   2   0.318 0.318
    ##  candidate     norm       1   3   0.321 0.320
    ##   baseline     norm       4   1   0.087 0.330
    ##   baseline     norm       4   2   0.086 0.327
    ##   baseline     norm       4   3   0.088 0.334
    ##  candidate     norm       4   1   0.089 0.339
    ##  candidate     norm       4   2   0.088 0.333
    ##  candidate     norm       4   3   0.088 0.334
    ##   baseline    munge       1   1   0.284 0.284
    ##   baseline    munge       1   2   0.286 0.287
    ##   baseline    munge       1   3   0.283 0.282
    ##  candidate    munge       1   1   0.247 0.246
    ##  candidate    munge       1   2   0.249 0.249
    ##  candidate    munge       1   3   0.247 0.247
    ##   baseline    munge       4   1   0.235 0.845
    ##   baseline    munge       4   2   0.236 0.848
    ##   baseline    munge       4   3   0.220 0.795
    ##  candidate    munge       4   1   0.068 0.259
    ##  candidate    munge       4   2   0.067 0.255
    ##  candidate    munge       4   3   0.067 0.253

## Reference work, counted separately

These are separate untimed runs with a benchmark-only interposition
library; counter overhead is absent from the timing table. All rows take
one reference lookup. Thus a cache hit is `input rows - htslib fetches`.
Fetched bases include read-ahead. One FASTA key is used per worker, so
there are no capacity evictions; handle destruction can occur when
DuckDB tears down workers. The independent native test, rather than this
locality workload, exercises more than eight references and verifies
destruction, retained-window limits, and exact bytes. The candidate also
runs at eight requested workers on the same full plain FASTA as a
regression stress case; handle-load counts expose actual workers.

    ##   revision operator threads    rows output_rows correct_rows fetches
    ##   baseline     norm       1 1000000       1e+06        1e+06      17
    ##  candidate     norm       1 1000000       1e+06        1e+06      17
    ##   baseline     norm       4 1000000       1e+06        1e+06      17
    ##  candidate     norm       4 1000000       1e+06        1e+06      17
    ##   baseline    munge       1 1000000       1e+06        1e+06 1000000
    ##  candidate    munge       1 1000000       1e+06        1e+06      17
    ##   baseline    munge       4 1000000       1e+06        1e+06 1000000
    ##  candidate    munge       4 1000000       1e+06        1e+06      17
    ##  candidate     norm       8 1000000       1e+06        1e+06      17
    ##  candidate    munge       8 1000000       1e+06        1e+06      17
    ##  fetched_bases cache_hits handle_loads handle_destroys
    ##        1114112     999983            1               0
    ##        1114112     999983            1               0
    ##        1114112     999983            4               0
    ##        1114112     999983            4               0
    ##        1000000          0            1               0
    ##        1114112     999983            1               0
    ##        1000000          0            4               0
    ##        1114112     999983            4               0
    ##        1114112     999983            8               0
    ##        1114112     999983            8               0

## Scope

The nearest earlier reports are [normalization](benchmark_norm.md) and
[munging](benchmark_munge.md). They use different workloads, so their
timings are not an identical-input baseline for this change. The paired
baseline above rebuilds the pre-change code with the candidate toolchain
and dependencies. This measures clustered SNVs on a real plain FASTA,
not arbitrary indels, random genome-wide access, remote transport, BGZF
throughput, or BAM generation. No result here establishes a general “no
regression” or “fastest” claim.

The first post-transport-policy run is retained in
[benchmark_reference_cache.md](benchmark_reference_cache.md). It
recorded slower one-worker munging and faster four-worker munging. A
fixed-CPU-placement repeat, when present as
`benchmark_reference_cache_affinity.md`, is a separate experiment, not a
replacement for that result. Its process affinity is recorded above; all
subprocesses inherit it. Neither experiment establishes the cause of a
difference by itself.

The shared cache retains at most eight handles and eight 65,536-base
sequence windows **per thread**, plus terminators and metadata. HTSlib
index/transport memory and caller-owned fetch results are separate;
remote BGZF tuning can retain additional per-handle blocks. Tuning
remains caller-specific and is part of cache identity; this patch does
not add public I/O options. This is not a sealed-allocation claim.

Independent cache tests use seed 171, twelve mixed-case/IUPAC references
(alternating plain and BGZF), eight concurrent workers, and 16,000
random requests checked against retained source bytes. They also test
missing-index creation, exact versus clipped contig ends, explicit
FAI/GZI and tuning keys, cache eviction, caller-copy isolation, requests
at and above 65,536 bases, and TLS destruction. ASan/UBSan and
`make test-reference-cache-tsan CC=clang` pass. The GCC TSAN runtime
cannot initialize on this host (unexpected memory mapping); Clang TSAN
passes both the standard probe and a directly HTSlib-linked build.
HTSlib itself is not sanitizer-instrumented in these probes; they do not
prove the absence of every dependency race. SQL/R checks additionally
exercise mixed normalization/munging calls, one million rows, input
reversal, and 1/4/8 workers.
