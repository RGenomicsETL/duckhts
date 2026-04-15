DuckHTS Liftover Benchmark
================

<!-- benchmark_liftover.md is generated from benchmark_liftover.Rmd. -->

# Benchmark

This benchmark:

- times DuckHTS liftover against `bcftools +liftover`
- compares normalized lifted outputs

# Run

`make bench-lift`

Override stress size with `LIFTOVER_BENCH_ROWS` and repeats with
`LIFTOVER_BENCH_RUNS`.

## Synthetic Stress Case

    #>     engine runs median_sec min_sec max_sec output_rows
    #> 1  duckhts    3      0.014   0.013   0.021       1e+05
    #> 2 bcftools    3      0.457   0.446   0.459       1e+05

    #>   duck_rows bcftools_rows only_duck_groups only_bcftools_groups
    #> 1     1e+05         1e+05                0                    0
    #>   mismatched_count_groups outputs_match
    #> 1                       0          TRUE
