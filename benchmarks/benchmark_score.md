DuckHTS Score Benchmark
================

<!-- benchmark_score.md is generated from benchmark_score.Rmd. -->

# Benchmark

This benchmark:

- times DuckHTS `bcftools_score()` against `bcftools +score`
- compares per-sample score outputs with numeric tolerance

# Run

`make bench-score`

Override stress size with `SCORE_BENCH_ROWS` (number of VCF variants)
and repeats with `SCORE_BENCH_RUNS`.

## Synthetic Stress Case

    #>     engine runs median_sec min_sec max_sec output_rows
    #> 1  duckhts    3      0.070   0.069   0.070          10
    #> 2 bcftools    3      0.093   0.090   0.096          10

    #>   duck_samples bcf_samples score_columns cell_matches cell_mismatches
    #> 1           10          10             1           10               0
    #>   only_duck_samples max_abs_diff max_rel_diff outputs_match
    #> 1                 0  0.001803492 1.767224e-05          TRUE
