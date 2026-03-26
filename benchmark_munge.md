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

    #>     engine runs median_sec min_sec max_sec output_rows
    #> 1  duckhts    3      1.225   1.212   1.272       1e+06
    #> 2 bcftools    3      0.826   0.824   0.888       1e+06

    #>   duck_rows bcftools_rows only_duck_groups only_bcftools_groups
    #> 1     1e+06         1e+06                2                    2
    #>   mismatched_count_groups outputs_match
    #> 1                       0         FALSE
