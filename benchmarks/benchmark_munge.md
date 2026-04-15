DuckHTS Munge Benchmark
================

<!-- benchmark_munge.md is generated from benchmark_munge.Rmd. -->

# Benchmark

This benchmark:

- times DuckHTS munge against `bcftools +munge`
- compares normalized output groups (`chrom`, `pos`, `ref`, `alt`,
  `filter`)

# Run

`make bench-munge`

Override stress size with `MUNGE_BENCH_ROWS` and repeats with
`MUNGE_BENCH_RUNS`.

## Synthetic Stress Case

    #>     engine threads runs median_sec min_sec max_sec output_rows
    #> 1  duckhts       1    3      0.972   0.384   0.996       2e+05
    #> 2  duckhts       4    3      0.308   0.298   0.330       2e+05
    #> 3 bcftools      NA    3      0.181   0.175   0.192       2e+05

    #>   duck_rows bcftools_rows only_duck_groups only_bcftools_groups
    #> 1     2e+05         2e+05                0                    0
    #>   mismatched_count_groups outputs_match
    #> 1                       0          TRUE

    #>   speedup_4_vs_1
    #> 1       3.155844
