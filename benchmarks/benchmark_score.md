DuckHTS Score Benchmark
================

<!-- benchmark_score.md is generated from benchmark_score.Rmd. -->

# Benchmark

This benchmark validates and times DuckHTS `bcftools_score()` and the
test-only query in `test/sql/score_reference.sql` against the
`bcftools +score` plugin bundled by
[`RGenomicsETL/RBCFTools`](https://github.com/RGenomicsETL/RBCFTools).
It covers both a single TSV/PLINK summary and multi-PRS TSV scoring
where several summary files are scored in one genotype scan. Because
upstream `bcftools +score` writes scores with `%#.6g`, the comparison
formats DuckHTS native scores the same way before checking cell-level
equality. SQL sums float32 products in double precision. Its comparison
uses the forward error bound
`(gamma(n, 2^-24) + gamma(n, 2^-53)) * A / (1-gamma(n, 2^-53))`, where
`A` is SQL `sum(abs(term))`, `gamma(n, u) = n*u/(1-n*u)` and `n` is the
number of contributing markers. The denominator accounts for rounding in
`A`; using `n` rather than `n-1` leaves slack for evaluating the bound.
This assumes finite normal arithmetic and `n*2^-24 < 1`. Counts must
match exactly. All engines read the same BGZF VCF and TSV summaries,
with counts enabled. Parity is checked before any timed repeat. DuckDB
uses one thread, a 4 GB memory limit and a 3 GB spill limit; upstream
uses one thread.

# Run

`make bench-score`

Useful overrides:

- `SCORE_BENCH_ROWS`: synthetic VCF variants, default `500000`
- `SCORE_BENCH_SAMPLES`: genotype samples, default `100`
- `SCORE_BENCH_PRS`: TSV summary files for the multi-PRS case, default
  `4`
- `SCORE_BENCH_RUNS`: timed repeats, default `3`
- `BCFTOOLS_BIN`: optional override for the bcftools executable
- `SCORE_PLUGIN_DIR`: optional override for the bcftools plugin
  directory

``` r
# Counts compare as integers, independently of R's numeric printing conventions.
control <- data.frame(SAMPLE = "S", prs = 1, prs_CNT = 500000)
integer_control <- transform(control, prs_CNT = 500000L)
stopifnot(compare_outputs(control, integer_control)$outputs_match)
stopifnot(!compare_outputs(control, transform(control, prs_CNT = 500001L))$outputs_match)
stopifnot(!compare_outputs(control, transform(control, prs = 2))$outputs_match)
stopifnot(!compare_outputs(control, transform(control, SAMPLE = "other"))$outputs_match)
```

## Tools

    #>               tool
    #> 1         bcftools
    #> 2 bcftools_plugins
    #>                                                                path
    #> 1     /usr/local/lib/R/site-library/RBCFTools/bcftools/bin/bcftools
    #> 2 /usr/local/lib/R/site-library/RBCFTools/bcftools/libexec/bcftools

## Settings

    #>   variants samples multi_prs_files runs threads memory_limit spill_limit
    #> 1   500000     100               4    3       1          4GB         3GB
    #>                            source_revision
    #> 1 cf399608566f190dcee801800686c23ca50c3eb3

## Synthetic Single-PRS Case

    #>               engine runs median_sec min_sec max_sec variants samples
    #> 1     duckhts_native    3      0.975   0.971   0.979   500000     100
    #> 2        duckhts_sql    3     19.376  19.372  19.386   500000     100
    #> 3 bcftools_RBCFTools    3      1.160   1.158   1.173   500000     100
    #>   genotype_cells output_samples score_columns score_cells
    #> 1       50000000            100             1         100
    #> 2       50000000            100             1         100
    #> 3       50000000            100             1         100

    #>   duck_samples bcf_samples score_columns count_columns compared_score_cells
    #> 1          100         100             1             1                  100
    #>   compared_count_cells cell_matches cell_mismatches only_duck_samples
    #> 1                  100          200               0                 0
    #>   only_bcftools_samples max_abs_diff_after_bcftools_text
    #> 1                     0                                0
    #>   max_rel_diff_after_bcftools_text outputs_match
    #> 1                                0          TRUE
    #>   score_cells count_cells max_abs_difference max_bound outputs_match
    #> 1         100         100        0.008542642  3997.236          TRUE

## Synthetic Multi-PRS TSV Case

    #>               engine runs median_sec min_sec max_sec variants samples
    #> 1     duckhts_native    3      1.532   1.530   1.533   500000     100
    #> 2        duckhts_sql    3     51.003  50.897  51.040   500000     100
    #> 3 bcftools_RBCFTools    3      1.634   1.633   1.635   500000     100
    #>   genotype_cells output_samples score_columns score_cells
    #> 1       50000000            100             4         400
    #> 2       50000000            100             4         400
    #> 3       50000000            100             4         400

    #>   duck_samples bcf_samples score_columns count_columns compared_score_cells
    #> 1          100         100             4             4                  400
    #>   compared_count_cells cell_matches cell_mismatches only_duck_samples
    #> 1                  400          800               0                 0
    #>   only_bcftools_samples max_abs_diff_after_bcftools_text
    #> 1                     0                                0
    #>   max_rel_diff_after_bcftools_text outputs_match
    #> 1                                0          TRUE
    #>   score_cells count_cells max_abs_difference max_bound outputs_match
    #> 1         400         400         0.01521208  4005.554          TRUE

## Scaling

The growth dimensions are variants (V), samples (S) and PRS files (P).
The base case is V = 250,000, S = 100, P = 1. Each ladder doubles one
dimension twice, holding the others at the base, and a joint case
doubles all three. Every process is a fresh DuckDB CLI (`DUCKDB`,
default `duckdb`) with a 4 GB memory limit and a 3 GB spill limit. Three
processes run per engine, shape and thread count (`SCORE_SCALING_REPS`).
Seconds are the CLI `.timer` wall time of the measured statement. Peak
RSS comes from GNU `time` and covers the whole process. DuckDB’s peak
buffer memory and peak spill come from its profiler.

- **Output consumed:** each query returns the row count, the sum of
  every count column and whether every score is finite. The SQL
  reference result must equal the native result for the same shape
  exactly. Score values are checked against native within the float32
  bound in the sections above.
- **Work model:** `read_geno` streams V records and unnests V × S calls.
  The calls probe a hash join whose build side is the V × P
  first-per-key markers. The aggregate has S × P groups. The live state
  is the materialized summary rows (V × P), the window that keeps the
  first marker per key (V × P) and the hash build (V × P). None of it
  grows with S. Native keeps one marker hash per summary (V × P) and
  decodes one record at a time.
- **Memory ceiling, set before measuring:** 256 MiB of fixed runtime
  overhead plus three times the live state. The live state is taken as
  768 bytes per marker per PRS for SQL: a summary row of about 200
  bytes, a window row of about 250 bytes and a hash-build row of about
  275 bytes, each with its path string. For native it is 128 bytes: a
  16-byte marker, the allele string, and the growth slack of the khash
  table and marker array. A process above its ceiling fails its cell.
- **Exponents:** `log2(t(2x) / t(x))` per doubling, and the per-doubling
  mean for the joint case. SQL work grows with the V × S × P matched
  calls, so each ladder should stay near 1. Native runs take less than 5
  s at every shape, so they get a memory verdict only.

### Cells

A cell passes when all its processes finish, return the native result
for the shape, report finite scores and stay under the ceiling.

| ladder   | scale | variants | samples | prs | engine         | threads | status | median_sec | range_sec   | peak_rss_mib | ceiling_mib | peak_buffer_mib | peak_spill_mib |
|:---------|:------|---------:|--------:|----:|:---------------|--------:|:-------|-----------:|:------------|-------------:|------------:|----------------:|---------------:|
| base     | 1x    |   250000 |     100 |   1 | duckhts_native |       1 | PASS   |       0.48 | 0.48-0.49   |           58 |         348 |               0 |              0 |
| base     | 1x    |   250000 |     100 |   1 | duckhts_native |       4 | PASS   |       0.50 | 0.48-0.50   |           58 |         348 |               0 |              0 |
| base     | 1x    |   250000 |     100 |   1 | duckhts_sql    |       1 | PASS   |       7.86 | 7.84-7.96   |          238 |         805 |             400 |              0 |
| base     | 1x    |   250000 |     100 |   1 | duckhts_sql    |       4 | PASS   |       7.80 | 7.79-7.84   |          293 |         805 |             370 |              0 |
| variants | 2x    |   500000 |     100 |   1 | duckhts_native |       1 | PASS   |       0.98 | 0.97-0.99   |           76 |         439 |               0 |              0 |
| variants | 2x    |   500000 |     100 |   1 | duckhts_native |       4 | PASS   |       1.01 | 0.97-1.01   |           76 |         439 |               0 |              0 |
| variants | 2x    |   500000 |     100 |   1 | duckhts_sql    |       1 | PASS   |      15.61 | 15.58-15.64 |          406 |        1355 |             716 |              0 |
| variants | 2x    |   500000 |     100 |   1 | duckhts_sql    |       4 | PASS   |      15.46 | 15.44-15.55 |          451 |        1355 |             704 |              0 |
| variants | 4x    |  1000000 |     100 |   1 | duckhts_native |       1 | PASS   |       1.97 | 1.96-1.97   |          111 |         622 |               0 |              0 |
| variants | 4x    |  1000000 |     100 |   1 | duckhts_native |       4 | PASS   |       2.01 | 2.00-2.01   |          111 |         622 |               0 |              0 |
| variants | 4x    |  1000000 |     100 |   1 | duckhts_sql    |       1 | PASS   |      31.31 | 31.26-31.52 |          739 |        2453 |            1358 |              0 |
| variants | 4x    |  1000000 |     100 |   1 | duckhts_sql    |       4 | PASS   |      30.93 | 30.83-30.95 |          794 |        2453 |            1455 |              0 |
| samples  | 2x    |   250000 |     200 |   1 | duckhts_native |       1 | PASS   |       0.80 | 0.80-0.81   |           58 |         348 |               0 |              0 |
| samples  | 2x    |   250000 |     200 |   1 | duckhts_native |       4 | PASS   |       0.83 | 0.80-0.83   |           58 |         348 |               0 |              0 |
| samples  | 2x    |   250000 |     200 |   1 | duckhts_sql    |       1 | PASS   |      15.02 | 14.87-15.04 |          255 |         805 |             427 |              0 |
| samples  | 2x    |   250000 |     200 |   1 | duckhts_sql    |       4 | PASS   |      14.54 | 14.51-14.57 |          323 |         805 |             385 |              0 |
| samples  | 4x    |   250000 |     400 |   1 | duckhts_native |       1 | PASS   |       1.46 | 1.46-1.48   |           58 |         348 |               0 |              0 |
| samples  | 4x    |   250000 |     400 |   1 | duckhts_native |       4 | PASS   |       1.49 | 1.48-1.51   |           59 |         348 |               0 |              0 |
| samples  | 4x    |   250000 |     400 |   1 | duckhts_sql    |       1 | PASS   |      29.11 | 29.02-29.25 |          296 |         805 |             481 |              0 |
| samples  | 4x    |   250000 |     400 |   1 | duckhts_sql    |       4 | PASS   |      28.17 | 28.03-28.17 |          351 |         805 |             441 |              0 |
| prs      | 2x    |   250000 |     100 |   2 | duckhts_native |       1 | PASS   |       0.58 | 0.58-0.59   |           76 |         439 |               0 |              0 |
| prs      | 2x    |   250000 |     100 |   2 | duckhts_native |       4 | PASS   |       0.59 | 0.59-0.59   |           76 |         439 |               0 |              0 |
| prs      | 2x    |   250000 |     100 |   2 | duckhts_sql    |       1 | PASS   |      12.14 | 12.10-12.17 |          408 |        1355 |             779 |              0 |
| prs      | 2x    |   250000 |     100 |   2 | duckhts_sql    |       4 | PASS   |      10.78 | 10.78-10.87 |          454 |        1355 |             722 |              0 |
| prs      | 4x    |   250000 |     100 |   4 | duckhts_native |       1 | PASS   |       0.77 | 0.77-0.78   |          111 |         622 |               0 |              0 |
| prs      | 4x    |   250000 |     100 |   4 | duckhts_native |       4 | PASS   |       0.78 | 0.77-0.80   |          111 |         622 |               0 |              0 |
| prs      | 4x    |   250000 |     100 |   4 | duckhts_sql    |       1 | PASS   |      21.37 | 21.30-21.48 |          758 |        2453 |            1537 |              0 |
| prs      | 4x    |   250000 |     100 |   4 | duckhts_sql    |       4 | PASS   |      18.60 | 18.54-18.68 |          774 |        2453 |            1405 |              0 |
| joint    | 2x2x2 |   500000 |     200 |   2 | duckhts_native |       1 | PASS   |       1.88 | 1.87-1.89   |          111 |         622 |               0 |              0 |
| joint    | 2x2x2 |   500000 |     200 |   2 | duckhts_native |       4 | PASS   |       1.87 | 1.86-1.93   |          111 |         622 |               0 |              0 |
| joint    | 2x2x2 |   500000 |     200 |   2 | duckhts_sql    |       1 | PASS   |      45.22 | 44.70-45.45 |          762 |        2453 |            1410 |              0 |
| joint    | 2x2x2 |   500000 |     200 |   2 | duckhts_sql    |       4 | PASS   |      41.79 | 41.77-41.81 |          783 |        2453 |            1385 |              0 |

### Exponents

SQL reference median seconds along each ladder, starting from the base
case.

| threads | ladder   | seconds              | exponents  |
|--------:|:---------|:---------------------|:-----------|
|       1 | variants | 7.86 / 15.61 / 31.31 | 0.99, 1    |
|       1 | samples  | 7.86 / 15.02 / 29.11 | 0.93, 0.95 |
|       1 | prs      | 7.86 / 12.14 / 21.37 | 0.63, 0.82 |
|       1 | joint    | 7.86 / 45.22         | 0.84       |
|       4 | variants | 7.8 / 15.46 / 30.93  | 0.99, 1    |
|       4 | samples  | 7.8 / 14.54 / 28.17  | 0.9, 0.95  |
|       4 | prs      | 7.8 / 10.78 / 18.6   | 0.47, 0.79 |
|       4 | joint    | 7.8 / 41.79          | 0.81       |
