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
    #> 1 0537c5e02a5ba03962bced7884bab9e2a70596a0

## Synthetic Single-PRS Case

    #>               engine runs median_sec min_sec max_sec variants samples
    #> 1     duckhts_native    3      0.991   0.991   0.992   500000     100
    #> 2        duckhts_sql    3     37.581  37.368  37.865   500000     100
    #> 3 bcftools_RBCFTools    3      1.164   1.162   1.173   500000     100
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
    #> 1     duckhts_native    3      1.514   1.507   1.532   500000     100
    #> 2        duckhts_sql    3     69.739  69.619  70.032   500000     100
    #> 3 bcftools_RBCFTools    3      1.632   1.631   1.644   500000     100
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
