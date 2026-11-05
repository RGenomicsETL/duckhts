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
    #> 1 1c9acb8d3f9491347156b26bd15783d85fd2e93a

## Synthetic Single-PRS Case

    #>               engine runs median_sec min_sec max_sec variants samples
    #> 1     duckhts_native    3      0.980   0.980   0.982   500000     100
    #> 2        duckhts_sql    3     19.403  19.394  19.424   500000     100
    #> 3 bcftools_RBCFTools    3      1.165   1.159   1.166   500000     100
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
    #> 1     duckhts_native    3      1.536   1.534   1.540   500000     100
    #> 2        duckhts_sql    3     51.083  51.031  51.361   500000     100
    #> 3 bcftools_RBCFTools    3      1.628   1.624   1.638   500000     100
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
| base     | 1x    |   250000 |     100 |   1 | duckhts_native |       4 | PASS   |       0.49 | 0.49-0.49   |           59 |         348 |               0 |              0 |
| base     | 1x    |   250000 |     100 |   1 | duckhts_sql    |       1 | PASS   |       8.03 | 8.03-8.09   |          237 |         805 |             402 |              0 |
| base     | 1x    |   250000 |     100 |   1 | duckhts_sql    |       4 | PASS   |       7.75 | 7.75-7.76   |          280 |         805 |             371 |              0 |
| variants | 2x    |   500000 |     100 |   1 | duckhts_native |       1 | PASS   |       0.98 | 0.97-0.99   |           76 |         439 |               0 |              0 |
| variants | 2x    |   500000 |     100 |   1 | duckhts_native |       4 | PASS   |       1.00 | 0.99-1.00   |           76 |         439 |               0 |              0 |
| variants | 2x    |   500000 |     100 |   1 | duckhts_sql    |       1 | PASS   |      15.80 | 15.70-15.80 |          410 |        1355 |             717 |              0 |
| variants | 2x    |   500000 |     100 |   1 | duckhts_sql    |       4 | PASS   |      15.46 | 15.38-15.52 |          456 |        1355 |             703 |              0 |
| variants | 4x    |  1000000 |     100 |   1 | duckhts_native |       1 | PASS   |       1.97 | 1.97-1.98   |          111 |         622 |               0 |              0 |
| variants | 4x    |  1000000 |     100 |   1 | duckhts_native |       4 | PASS   |       2.00 | 2.00-2.01   |          111 |         622 |               0 |              0 |
| variants | 4x    |  1000000 |     100 |   1 | duckhts_sql    |       1 | PASS   |      31.68 | 31.51-31.72 |          739 |        2453 |            1360 |              0 |
| variants | 4x    |  1000000 |     100 |   1 | duckhts_sql    |       4 | PASS   |      30.74 | 30.69-30.95 |          802 |        2453 |            1459 |              0 |
| samples  | 2x    |   250000 |     200 |   1 | duckhts_native |       1 | PASS   |       0.81 | 0.81-0.81   |           59 |         348 |               0 |              0 |
| samples  | 2x    |   250000 |     200 |   1 | duckhts_native |       4 | PASS   |       0.83 | 0.82-0.84   |           59 |         348 |               0 |              0 |
| samples  | 2x    |   250000 |     200 |   1 | duckhts_sql    |       1 | PASS   |      14.96 | 14.94-15.06 |          255 |         805 |             427 |              0 |
| samples  | 2x    |   250000 |     200 |   1 | duckhts_sql    |       4 | PASS   |      14.43 | 14.41-14.51 |          289 |         805 |             385 |              0 |
| samples  | 4x    |   250000 |     400 |   1 | duckhts_native |       1 | PASS   |       1.46 | 1.45-1.47   |           58 |         348 |               0 |              0 |
| samples  | 4x    |   250000 |     400 |   1 | duckhts_native |       4 | PASS   |       1.49 | 1.48-1.50   |           59 |         348 |               0 |              0 |
| samples  | 4x    |   250000 |     400 |   1 | duckhts_sql    |       1 | PASS   |      29.26 | 29.24-29.35 |          291 |         805 |             481 |              0 |
| samples  | 4x    |   250000 |     400 |   1 | duckhts_sql    |       4 | PASS   |      28.11 | 28.03-28.13 |          316 |         805 |             438 |              0 |
| prs      | 2x    |   250000 |     100 |   2 | duckhts_native |       1 | PASS   |       0.58 | 0.58-0.59   |           76 |         439 |               0 |              0 |
| prs      | 2x    |   250000 |     100 |   2 | duckhts_native |       4 | PASS   |       0.59 | 0.59-0.60   |           76 |         439 |               0 |              0 |
| prs      | 2x    |   250000 |     100 |   2 | duckhts_sql    |       1 | PASS   |      12.16 | 12.14-12.19 |          408 |        1355 |             780 |              0 |
| prs      | 2x    |   250000 |     100 |   2 | duckhts_sql    |       4 | PASS   |      10.82 | 10.78-10.98 |          472 |        1355 |             706 |              0 |
| prs      | 4x    |   250000 |     100 |   4 | duckhts_native |       1 | PASS   |       0.78 | 0.77-0.79   |          111 |         622 |               0 |              0 |
| prs      | 4x    |   250000 |     100 |   4 | duckhts_native |       4 | PASS   |       0.79 | 0.77-0.80   |          111 |         622 |               0 |              0 |
| prs      | 4x    |   250000 |     100 |   4 | duckhts_sql    |       1 | PASS   |      21.38 | 21.36-21.44 |          753 |        2453 |            1539 |              0 |
| prs      | 4x    |   250000 |     100 |   4 | duckhts_sql    |       4 | PASS   |      18.55 | 18.45-18.70 |          779 |        2453 |            1413 |              0 |
| joint    | 2x2x2 |   500000 |     200 |   2 | duckhts_native |       1 | PASS   |       1.87 | 1.86-1.89   |          111 |         622 |               0 |              0 |
| joint    | 2x2x2 |   500000 |     200 |   2 | duckhts_native |       4 | PASS   |       1.87 | 1.87-1.89   |          111 |         622 |               0 |              0 |
| joint    | 2x2x2 |   500000 |     200 |   2 | duckhts_sql    |       1 | PASS   |      44.79 | 44.72-45.28 |          764 |        2453 |            1413 |              0 |
| joint    | 2x2x2 |   500000 |     200 |   2 | duckhts_sql    |       4 | PASS   |      41.94 | 41.56-41.96 |          792 |        2453 |            1390 |              0 |

### Exponents

SQL reference median seconds along each ladder, starting from the base
case.

| threads | ladder   | seconds              | exponents  |
|--------:|:---------|:---------------------|:-----------|
|       1 | variants | 8.03 / 15.8 / 31.68  | 0.98, 1    |
|       1 | samples  | 8.03 / 14.96 / 29.26 | 0.9, 0.97  |
|       1 | prs      | 8.03 / 12.16 / 21.38 | 0.6, 0.81  |
|       1 | joint    | 8.03 / 44.79         | 0.83       |
|       4 | variants | 7.75 / 15.46 / 30.74 | 1, 0.99    |
|       4 | samples  | 7.75 / 14.43 / 28.11 | 0.9, 0.96  |
|       4 | prs      | 7.75 / 10.82 / 18.55 | 0.48, 0.78 |
|       4 | joint    | 7.75 / 41.94         | 0.81       |
