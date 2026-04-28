DuckHTS Score Benchmark
================

<!-- benchmark_score.md is generated from benchmark_score.Rmd. -->

# Benchmark

This benchmark validates and times DuckHTS `bcftools_score()` against
the `bcftools +score` plugin bundled by
[`RGenomicsETL/RBCFTools`](https://github.com/RGenomicsETL/RBCFTools).
It covers both a single TSV/PLINK summary and multi-PRS TSV scoring
where several summary files are scored in one genotype scan.

# Run

`make bench-score`

Useful overrides:

- `SCORE_BENCH_ROWS`: synthetic VCF variants, default `100000`
- `SCORE_BENCH_SAMPLES`: genotype samples, default `10`
- `SCORE_BENCH_PRS`: TSV summary files for the multi-PRS case, default
  `4`
- `SCORE_BENCH_RUNS`: timed repeats, default `3`
- `BCFTOOLS_BIN`: optional override for the bcftools executable
- `SCORE_PLUGIN_DIR`: optional override for the bcftools plugin
  directory

## Tools

    #>               tool
    #> 1         bcftools
    #> 2 bcftools_plugins
    #>                                                                path
    #> 1     /usr/local/lib/R/site-library/RBCFTools/bcftools/bin/bcftools
    #> 2 /usr/local/lib/R/site-library/RBCFTools/bcftools/libexec/bcftools

## Settings

    #>   variants samples multi_prs_files runs
    #> 1    50000      10               4    3

## Synthetic Single-PRS Case

    #>               engine runs median_sec min_sec max_sec output_rows
    #> 1            duckhts    3      0.034   0.034   0.034          10
    #> 2 bcftools_RBCFTools    3      0.045   0.043   0.046          10

    #>   duck_samples bcf_samples score_columns cell_matches cell_mismatches
    #> 1           10          10             1           10               0
    #>   only_duck_samples only_bcftools_samples max_abs_diff max_rel_diff
    #> 1                 0                     0 0.0007000857 5.150188e-06
    #>   outputs_match
    #> 1          TRUE

## Synthetic Multi-PRS TSV Case

    #>               engine runs median_sec min_sec max_sec output_rows
    #> 1            duckhts    3      0.072    0.07   0.075          10
    #> 2 bcftools_RBCFTools    3      0.082    0.08   0.086          10

    #>   duck_samples bcf_samples score_columns cell_matches cell_mismatches
    #> 1           10          10             4           40               0
    #>   only_duck_samples only_bcftools_samples max_abs_diff max_rel_diff
    #> 1                 0                     0 0.0008962599 1.621158e-05
    #>   outputs_match
    #> 1          TRUE
