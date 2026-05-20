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
- `LIFTOVER_REAL_VCF`, `LIFTOVER_REAL_REGION`, `LIFTOVER_REAL_CHAIN`,
  `LIFTOVER_REAL_SRC_FASTA`, `LIFTOVER_REAL_DST_FASTA`: optional
  real-callset conformance case (for example GIAB HG001 GRCh37 chr20)

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
    #> 1        1000000    3        <NA>

## Synthetic Stress Case

    #>               engine runs median_sec min_sec max_sec output_rows
    #> 1            duckhts    3      0.047   0.046   0.054       1e+06
    #> 2 bcftools_RBCFTools    3    109.419 107.901 112.497       1e+06

    #>   duck_rows bcftools_rows only_duck_groups only_bcftools_groups
    #> 1     1e+06         1e+06                0                    0
    #>   mismatched_count_groups outputs_match
    #> 1                       0          TRUE

## Optional Real-Callset Case

    #>                                                                                                                           note
    #> 1 Set LIFTOVER_REAL_VCF, LIFTOVER_REAL_CHAIN, LIFTOVER_REAL_SRC_FASTA, and LIFTOVER_REAL_DST_FASTA to run a real-callset case.
