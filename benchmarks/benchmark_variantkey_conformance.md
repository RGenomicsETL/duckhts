DuckHTS VariantKey Conformance
================

<!-- benchmark_variantkey_conformance.md is generated from benchmark_variantkey_conformance.Rmd. -->

# Conformance

This report compares DuckHTS VariantKey output against installed
`bcftools` `%VKX` formatting on real whole-genome GIAB VCF files.

DuckHTS uses the vendored official VariantKey C API and follows the
bcftools `%VKX` / `+add-variantkey` convention of accepting 1-based VCF
`POS` at the SQL boundary while encoding the upstream 0-based field
internally.

VariantKey / RegionKey are based on Nicola Asuni’s work:
<https://doi.org/10.1101/473744>.

## Run

``` sh
make bench-variantkey
```

Useful overrides:

- `BCFTOOLS_BIN`: optional bcftools executable override
- `VARIANTKEY_REAL_CASES`: optional comma-separated local VCF/BCF paths
- `VARIANTKEY_OUT_DIR`: optional directory for emitted TSVs and
  summaries
- `DUCKHTS_EXTENSION`: optional extension path override

If `VARIANTKEY_REAL_CASES` is unset, this report looks for the locally
downloaded whole-WGS GIAB files already used in earlier DuckHTS
conformance work under `/root/giab_norm/`.

``` r
results <- do.call(
  rbind,
  lapply(seq_len(nrow(cases)), function(i) run_case(cases$label[i], cases$path[i]))
)
results
#>                    label
#> 1 giab_hg001_grch37_full
#> 2 giab_hg002_grch37_full
#> 3 giab_hg006_grch37_full
#>                                                        path bcftools_rows
#> 1 /root/giab_norm/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz       3891440
#> 2 /root/giab_norm/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz       4033796
#> 3 /root/giab_norm/HG006_GRCh37_1_22_v4.2.1_benchmark.vcf.gz       3878664
#>   duckhts_rows mismatch_groups elapsed_seconds
#> 1      3891440               0           3.494
#> 2      4033796               0           3.863
#> 3      3878664               0           3.324
```

All listed cases should have `mismatch_groups == 0`.

## Notes

- This framework compares against bcftools `%VKX`, which is the same
  VariantKey hexadecimal output exposed by `+add-variantkey`.
- DuckHTS currently compares `ALT[1]` from `read_bcf(...)`, mirroring
  the biallelic `%VKX` field that bcftools derives from the first ALT
  allele on the current record.
- Large, ambiguous, and symbolic alleles are expected to fall back to
  the official hashed nonreversible VariantKey mode. That still
  preserves bcftools `%VKX` parity, but it does not encode `END`,
  `SVLEN`, mate breakend coordinates, or other SV metadata.
