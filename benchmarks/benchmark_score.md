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
    #> 1 2801319c4187037e65a11cc1a375223a71f61472

## Synthetic Single-PRS Case

    #>               engine runs median_sec min_sec max_sec variants samples
    #> 1     duckhts_native    3      0.969   0.964   0.982   500000     100
    #> 2        duckhts_sql    3     19.719  19.641  19.721   500000     100
    #> 3 bcftools_RBCFTools    3      1.166   1.151   1.166   500000     100
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
    #> 1     duckhts_native    3      1.523   1.520   1.561   500000     100
    #> 2        duckhts_sql    3     52.189  52.069  52.490   500000     100
    #> 3 bcftools_RBCFTools    3      1.629   1.620   1.630   500000     100
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
| base     | 1x    |   250000 |     100 |   1 | duckhts_native |       1 | PASS   |       0.48 | 0.48-0.50   |           58 |         348 |               0 |              0 |
| base     | 1x    |   250000 |     100 |   1 | duckhts_native |       4 | PASS   |       0.50 | 0.49-0.50   |           59 |         348 |               0 |              0 |
| base     | 1x    |   250000 |     100 |   1 | duckhts_sql    |       1 | PASS   |       8.20 | 8.15-8.46   |          217 |         805 |             339 |              0 |
| base     | 1x    |   250000 |     100 |   1 | duckhts_sql    |       4 | PASS   |       7.99 | 7.99-8.02   |          284 |         805 |             364 |              0 |
| variants | 2x    |   500000 |     100 |   1 | duckhts_native |       1 | PASS   |       0.98 | 0.98-1.04   |           76 |         439 |               0 |              0 |
| variants | 2x    |   500000 |     100 |   1 | duckhts_native |       4 | PASS   |       1.00 | 0.98-1.01   |           76 |         439 |               0 |              0 |
| variants | 2x    |   500000 |     100 |   1 | duckhts_sql    |       1 | PASS   |      16.02 | 15.99-16.28 |          373 |        1355 |             654 |              0 |
| variants | 2x    |   500000 |     100 |   1 | duckhts_sql    |       4 | PASS   |      15.87 | 15.85-15.98 |          466 |        1355 |             714 |              0 |
| variants | 4x    |  1000000 |     100 |   1 | duckhts_native |       1 | PASS   |       1.96 | 1.96-1.97   |          111 |         622 |               0 |              0 |
| variants | 4x    |  1000000 |     100 |   1 | duckhts_native |       4 | PASS   |       2.00 | 1.98-2.01   |          111 |         622 |               0 |              0 |
| variants | 4x    |  1000000 |     100 |   1 | duckhts_sql    |       1 | PASS   |      32.07 | 32.05-32.23 |          687 |        2453 |            1296 |              0 |
| variants | 4x    |  1000000 |     100 |   1 | duckhts_sql    |       4 | PASS   |      31.77 | 31.67-32.50 |          739 |        2453 |            1386 |              0 |
| samples  | 2x    |   250000 |     200 |   1 | duckhts_native |       1 | PASS   |       0.81 | 0.81-0.81   |           58 |         348 |               0 |              0 |
| samples  | 2x    |   250000 |     200 |   1 | duckhts_native |       4 | PASS   |       0.83 | 0.82-0.83   |           59 |         348 |               0 |              0 |
| samples  | 2x    |   250000 |     200 |   1 | duckhts_sql    |       1 | PASS   |      15.14 | 15.11-15.22 |          235 |         805 |             368 |              0 |
| samples  | 2x    |   250000 |     200 |   1 | duckhts_sql    |       4 | PASS   |      15.09 | 14.69-15.25 |          299 |         805 |             383 |              0 |
| samples  | 4x    |   250000 |     400 |   1 | duckhts_native |       1 | PASS   |       1.47 | 1.46-1.56   |           58 |         348 |               0 |              0 |
| samples  | 4x    |   250000 |     400 |   1 | duckhts_native |       4 | PASS   |       1.50 | 1.48-1.54   |           59 |         348 |               0 |              0 |
| samples  | 4x    |   250000 |     400 |   1 | duckhts_sql    |       1 | PASS   |      29.39 | 29.28-29.53 |          272 |         805 |             422 |              0 |
| samples  | 4x    |   250000 |     400 |   1 | duckhts_sql    |       4 | PASS   |      29.34 | 28.90-30.20 |          337 |         805 |             435 |              0 |
| prs      | 2x    |   250000 |     100 |   2 | duckhts_native |       1 | PASS   |       0.58 | 0.58-0.59   |           76 |         439 |               0 |              0 |
| prs      | 2x    |   250000 |     100 |   2 | duckhts_native |       4 | PASS   |       0.59 | 0.58-0.60   |           76 |         439 |               0 |              0 |
| prs      | 2x    |   250000 |     100 |   2 | duckhts_sql    |       1 | PASS   |      12.52 | 12.47-12.66 |          370 |        1355 |             656 |              0 |
| prs      | 2x    |   250000 |     100 |   2 | duckhts_sql    |       4 | PASS   |      11.11 | 10.90-11.23 |          441 |        1355 |             707 |              0 |
| prs      | 4x    |   250000 |     100 |   4 | duckhts_native |       1 | PASS   |       0.78 | 0.77-0.79   |          111 |         622 |               0 |              0 |
| prs      | 4x    |   250000 |     100 |   4 | duckhts_native |       4 | PASS   |       0.77 | 0.77-0.79   |          111 |         622 |               0 |              0 |
| prs      | 4x    |   250000 |     100 |   4 | duckhts_sql    |       1 | PASS   |      22.22 | 22.12-22.25 |          685 |        2453 |            1290 |              0 |
| prs      | 4x    |   250000 |     100 |   4 | duckhts_sql    |       4 | PASS   |      18.78 | 18.73-19.23 |          733 |        2453 |            1393 |              0 |
| joint    | 2x2x2 |   500000 |     200 |   2 | duckhts_native |       1 | PASS   |       1.88 | 1.88-1.91   |          111 |         622 |               0 |              0 |
| joint    | 2x2x2 |   500000 |     200 |   2 | duckhts_native |       4 | PASS   |       1.92 | 1.89-1.92   |          111 |         622 |               0 |              0 |
| joint    | 2x2x2 |   500000 |     200 |   2 | duckhts_sql    |       1 | PASS   |      45.78 | 45.33-45.99 |          696 |        2453 |            1287 |              0 |
| joint    | 2x2x2 |   500000 |     200 |   2 | duckhts_sql    |       4 | PASS   |      42.85 | 42.74-43.56 |          748 |        2453 |            1387 |              0 |

### Exponents

SQL reference median seconds along each ladder, starting from the base
case.

| threads | ladder   | seconds              | exponents  |
|--------:|:---------|:---------------------|:-----------|
|       1 | variants | 8.2 / 16.02 / 32.07  | 0.97, 1    |
|       1 | samples  | 8.2 / 15.14 / 29.39  | 0.88, 0.96 |
|       1 | prs      | 8.2 / 12.52 / 22.22  | 0.61, 0.83 |
|       1 | joint    | 8.2 / 45.78          | 0.83       |
|       4 | variants | 7.99 / 15.87 / 31.77 | 0.99, 1    |
|       4 | samples  | 7.99 / 15.09 / 29.34 | 0.92, 0.96 |
|       4 | prs      | 7.99 / 11.11 / 18.78 | 0.48, 0.76 |
|       4 | joint    | 7.99 / 42.85         | 0.81       |
