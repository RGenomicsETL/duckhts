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
  `LIFTOVER_REAL_SRC_FASTA`, `LIFTOVER_REAL_DST_FASTA`: override the
  rendered real-callset case; by default the benchmark renders the local
  full-file GIAB HG001 GRCh37 case from `/root/giab_norm/` against
  `/root/GRCh37/GRCh37_to_GRCh38.chain.gz`

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
    #> 1            duckhts    3      0.048   0.047   0.049       1e+06
    #> 2 bcftools_RBCFTools    3    101.570 101.306 104.393       1e+06

    #>   duck_rows bcftools_rows only_duck_groups only_bcftools_groups
    #> 1     1e+06         1e+06                0                    0
    #>   mismatched_count_groups outputs_match
    #> 1                       0          TRUE

## Real-Callset Case

    #>                                                   input_vcf     region
    #> 1 /root/giab_norm/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz full_input
    #>                                    chain                        src_fasta
    #> 1 /root/GRCh37/GRCh37_to_GRCh38.chain.gz /root/GRCh37/human_g1k_v37.fasta
    #>                                                      dst_fasta
    #> 1 /root/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

    #>               engine runs median_sec min_sec max_sec output_rows
    #> 1            duckhts    3      5.639   5.615   5.704     3886470
    #> 2 bcftools_RBCFTools    3     12.104  12.087  12.212     3886470

    #>   duck_rows bcftools_rows only_duck_groups only_bcftools_groups
    #> 1   3886470       3886470                0                    0
    #>   mismatched_count_groups outputs_match
    #> 1                       0          TRUE
