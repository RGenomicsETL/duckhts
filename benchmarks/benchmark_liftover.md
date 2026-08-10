DuckHTS Liftover Benchmark
================

<!-- benchmark_liftover.md is generated from benchmark_liftover.Rmd. -->

# Benchmark

This benchmark:

- times DuckHTS liftover against `bcftools +liftover`
- compares normalized lifted outputs

# Run

`make bench-lift`

Useful overrides:

- `LIFTOVER_BENCH_ROWS`: synthetic VCF records, default `1000000`
- `LIFTOVER_BENCH_RUNS`: timed repeats, default `3`
- `BCFTOOLS_BIN`: optional override for the bcftools executable
- `BCFTOOLS_PLUGIN_DIR`: optional override for the bcftools plugin
  directory
- `LIFTOVER_REAL_VCF`, `LIFTOVER_REAL_PROVENANCE`,
  `LIFTOVER_REAL_REGION`, `LIFTOVER_REAL_CHAIN`,
  `LIFTOVER_REAL_SRC_FASTA`, `LIFTOVER_REAL_DST_FASTA`: cached
  real-callset inputs; use `make stage-giab-v4.2.1` and
  `make stage-liftover-references` to populate their defaults

For a shell-first real-callset comparison, see
`bash scripts/liftover_conformance.sh ...`.

For curated GIAB benchmark slices, see
`scripts/conformance_case_table.tsv` and run them in batch with
`bash scripts/liftover_conformance_batch.sh`.

## Tools

    #>               tool
    #> 1         bcftools
    #> 2 bcftools_plugins
    #>                                                                path
    #> 1     /usr/local/lib/R/site-library/RBCFTools/bcftools/bin/bcftools
    #> 2 /usr/local/lib/R/site-library/RBCFTools/bcftools/libexec/bcftools

## Settings

    #>   synthetic_rows runs real_region
    #> 1            100    1        <NA>

## Synthetic Stress Case

    #>               engine runs median_sec min_sec max_sec output_rows
    #> 1            duckhts    1      0.003   0.003   0.003         100
    #> 2 bcftools_RBCFTools    1      0.009   0.009   0.009         100

    #>   duck_rows bcftools_rows only_duck_groups only_bcftools_groups
    #> 1       100           100                0                    0
    #>   mismatched_count_groups outputs_match
    #> 1                       0          TRUE

## Real-Callset Case

    #>                                   status
    #> 1 skipped; set all LIFTOVER_REAL_* paths
