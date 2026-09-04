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
    ##  candidate 62f9b8786852817b8db9747bd185666dac080d33
    ##                     extension_md5
    ##  bcf15693b036abeccc30941db47821b1
    ##  1895545a05cf4548e66d772efe2d15fd

    ## R: R version 4.6.0 (2026-04-24)
    ## DuckDB: 1.5.3

    ## Linux 6.8.0-78-generic x86_64 GNU/Linux

    ## Model name:                           13th Gen Intel(R) Core(TM) i5-13500 BIOS Model name:                      13th Gen Intel(R) Core(TM) i5-13500 To Be Filled By O.E.M. CPU @ 2.4GHz

    ## Rows per query: 1000000 ; measured repetitions: 3

    ## Ubuntu clang version 18.1.3 (1ubuntu1)

## Timings

Each function/thread/revision combination uses a fresh R process;
repetitions within that process reuse worker caches when DuckDB retains
workers. Table construction, FASTA staging, indexing, and R startup are
excluded from elapsed and CPU time. CPU is user + system process CPU
during the query. Peak RSS is the **whole-process high-water mark**,
including setup, not cache-only memory. The operating-system page cache
is warm; this is not a cold-storage benchmark.

    ##   revision operator threads median_elapsed_s median_cpu_s median_peak_rss_kib
    ##   baseline    munge       1            0.287        0.287              231264
    ##  candidate    munge       1            0.252        0.251              231732
    ##   baseline     norm       1            0.326        0.326              231392
    ##  candidate     norm       1            0.323        0.322              231504
    ##   baseline    munge       4            0.154        0.570              231496
    ##  candidate    munge       4            0.067        0.260              231340
    ##   baseline     norm       4            0.087        0.335              231956
    ##  candidate     norm       4            0.087        0.334              231880

    ##   revision operator threads run elapsed   cpu
    ##   baseline     norm       1   1   0.329 0.329
    ##   baseline     norm       1   2   0.323 0.323
    ##   baseline     norm       1   3   0.326 0.326
    ##  candidate     norm       1   1   0.323 0.322
    ##  candidate     norm       1   2   0.322 0.321
    ##  candidate     norm       1   3   0.326 0.325
    ##   baseline     norm       4   1   0.089 0.341
    ##   baseline     norm       4   2   0.087 0.326
    ##   baseline     norm       4   3   0.086 0.335
    ##  candidate     norm       4   1   0.090 0.341
    ##  candidate     norm       4   2   0.086 0.326
    ##  candidate     norm       4   3   0.087 0.334
    ##   baseline    munge       1   1   0.287 0.287
    ##   baseline    munge       1   2   0.284 0.283
    ##   baseline    munge       1   3   0.288 0.289
    ##  candidate    munge       1   1   0.255 0.254
    ##  candidate    munge       1   2   0.252 0.251
    ##  candidate    munge       1   3   0.251 0.251
    ##   baseline    munge       4   1   0.182 0.656
    ##   baseline    munge       4   2   0.154 0.570
    ##   baseline    munge       4   3   0.150 0.554
    ##  candidate    munge       4   1   0.077 0.282
    ##  candidate    munge       4   2   0.067 0.260
    ##  candidate    munge       4   3   0.067 0.259

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

The shared cache retains at most eight handles and eight 65,536-base
sequence windows **per thread**, plus terminators and metadata. HTSlib
index/transport memory and caller-owned fetch results are separate;
remote BGZF tuning can retain additional per-handle blocks. This is not
a sealed-allocation claim.

Independent cache tests use seed 171, twelve mixed-case/IUPAC references
(alternating plain and BGZF), eight concurrent workers, and 16,000
random requests checked against retained source bytes. They also test
missing-index creation, exact versus clipped contig ends, explicit
FAI/GZI keys, cache eviction, caller-copy isolation, requests at and
above 65,536 bases, and TLS destruction. ASan/UBSan and
`make test-reference-cache-tsan CC=clang` pass. The GCC TSAN runtime
cannot initialize on this host (unexpected memory mapping); Clang TSAN
passes both the standard probe and a directly HTSlib-linked build.
HTSlib itself is not sanitizer-instrumented in these probes; they do not
prove the absence of every dependency race. SQL/R checks additionally
exercise mixed normalization/munging calls, one million rows, input
reversal, and 1/4/8 workers.
