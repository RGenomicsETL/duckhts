DuckHTS Score Benchmark
================

<!-- benchmark_score.md is generated from benchmark_score.Rmd. -->

# Benchmark

This benchmark validates and times DuckHTS `bcftools_score()` against
the `bcftools +score` plugin bundled by
[`RGenomicsETL/RBCFTools`](https://github.com/RGenomicsETL/RBCFTools).
It covers both a single TSV/PLINK summary and multi-PRS TSV scoring
where several summary files are scored in one genotype scan. Because
upstream `bcftools +score` writes scores with `%#.6g`, the comparison
formats DuckHTS scores the same way before checking cell-level equality.

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

## Tools

    #>               tool
    #> 1         bcftools
    #> 2 bcftools_plugins
    #>                                                                path
    #> 1     /usr/local/lib/R/site-library/RBCFTools/bcftools/bin/bcftools
    #> 2 /usr/local/lib/R/site-library/RBCFTools/bcftools/libexec/bcftools

## Settings

    #>   variants samples multi_prs_files runs
    #> 1   500000     100               4    3

## Synthetic Single-PRS Case

    #>               engine runs median_sec min_sec max_sec variants samples
    #> 1            duckhts    3      0.967   0.961   0.968   500000     100
    #> 2 bcftools_RBCFTools    3      1.132   1.132   1.132   500000     100
    #>   genotype_cells output_samples score_columns score_cells
    #> 1       50000000            100             1         100
    #> 2       50000000            100             1         100

    #>   duck_samples bcf_samples score_columns compared_score_cells cell_matches
    #> 1          100         100             1                  100          100
    #>   cell_mismatches only_duck_samples only_bcftools_samples
    #> 1               0                 0                     0
    #>   max_abs_diff_after_bcftools_text max_rel_diff_after_bcftools_text
    #> 1                                0                                0
    #>   outputs_match
    #> 1          TRUE

## Synthetic Multi-PRS TSV Case

    #>               engine runs median_sec min_sec max_sec variants samples
    #> 1            duckhts    3      1.493   1.492   1.511   500000     100
    #> 2 bcftools_RBCFTools    3      1.588   1.587   1.599   500000     100
    #>   genotype_cells output_samples score_columns score_cells
    #> 1       50000000            100             4         400
    #> 2       50000000            100             4         400

    #>   duck_samples bcf_samples score_columns compared_score_cells cell_matches
    #> 1          100         100             4                  400          400
    #>   cell_mismatches only_duck_samples only_bcftools_samples
    #> 1               0                 0                     0
    #>   max_abs_diff_after_bcftools_text max_rel_diff_after_bcftools_text
    #> 1                                0                                0
    #>   outputs_match
    #> 1          TRUE
