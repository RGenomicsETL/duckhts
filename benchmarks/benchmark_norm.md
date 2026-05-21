DuckHTS Normalization Benchmark
================

<!-- benchmark_norm.md is generated from benchmark_norm.Rmd. -->

# Benchmark

This benchmark:

- times DuckHTS `duckhts_bcftools_norm(...)` against installed
  `bcftools norm`
- compares normalized outputs after grouping identical normalized
  alleles
- covers both site-preserving normalization and
  `split_multiallelic := TRUE`

# Run

`make bench-norm`

Useful overrides:

- `NORM_BENCH_ROWS`: synthetic VCF records, default `100000`
- `NORM_BENCH_RUNS`: timed repeats, default `3`
- `BCFTOOLS_BIN`: optional override for the bcftools executable
- `NORM_REAL_VCF`, `NORM_REAL_REGION`, `NORM_REAL_FASTA`: optional
  real-callset conformance case
- `NORM_REAL_SPLIT`: `1` to run the optional real-callset case in split
  mode too
- `NORM_REAL_LABEL`: optional label for the real-callset section

For a shell-first comparison, see
`bash scripts/norm_conformance.sh ...`.

# Synthetic fixture

``` r
fixture <- write_norm_fixture(file.path(tempdir(), "duckhts_norm_bench_fixture"), bench_rows)
site_case <- run_case(fixture$input_path, fixture$fasta_path, split = FALSE, n = bench_runs)
split_case <- run_case(fixture$input_path, fixture$fasta_path, split = TRUE, n = bench_runs)

synthetic_bench <- rbind(site_case$bench, split_case$bench)
knitr::kable(synthetic_bench, digits = 3)
```

| engine        | mode               | runs | median_sec | min_sec | max_sec | output_rows |
|:--------------|:-------------------|-----:|-----------:|--------:|--------:|------------:|
| duckhts       | site_preserving    |    3 |      0.162 |   0.154 |   0.163 |      100000 |
| bcftools_norm | site_preserving    |    3 |     66.333 |  66.036 |  66.672 |      100000 |
| duckhts       | split_multiallelic |    3 |      0.202 |   0.202 |   0.203 |      120000 |
| bcftools_norm | split_multiallelic |    3 |     99.937 |  99.063 | 101.207 |      120000 |

Synthetic fixture mix:

- normalized insertions/deletions
- a multiallelic `TT,TTT` row
- a ref-only `ALT='.'` row that exercises split-mode row preservation

## Synthetic conformance

### Site-preserving mode

``` r
knitr::kable(site_case$compare_summary)
```

| status |   n |
|:-------|----:|
| match  |   5 |

``` r
knitr::kable(site_case$duck_status)
```

| norm_status |     n |
|:------------|------:|
| Normalized  | 80000 |
| RefOnly     | 20000 |

``` r
knitr::kable(site_case$mismatches)
```

| chrom | pos_normed | end_pos_normed | ref_normed | alt_normed | duck_n | bcf_n | status |
|:------|-----------:|---------------:|:-----------|:-----------|-------:|------:|:-------|

### Split mode

``` r
knitr::kable(split_case$compare_summary)
```

| status |   n |
|:-------|----:|
| match  |   4 |

``` r
knitr::kable(split_case$duck_status)
```

| norm_status |     n |
|:------------|------:|
| Normalized  | 1e+05 |
| RefOnly     | 2e+04 |

``` r
knitr::kable(split_case$mismatches)
```

| chrom | pos_normed | end_pos_normed | ref_normed | alt_normed | duck_n | bcf_n | status |
|:------|-----------:|---------------:|:-----------|:-----------|-------:|------:|:-------|

# Optional real-callset case

``` r
if (!nzchar(real_input_path) || !nzchar(real_fasta)) {
  cat("Real-callset benchmark skipped. Set `NORM_REAL_VCF` and `NORM_REAL_FASTA` to enable it.\n")
} else {
  real_fixture <- prepare_real_fixture(real_input_path, real_fasta, real_region)
  real_site <- run_case(real_fixture$input_path, real_fixture$fasta_path, split = FALSE, n = bench_runs)
  cat("## ", real_label, " (site-preserving)\n\n", sep = "")
  print(knitr::kable(real_site$bench, digits = 3))
  cat("\n")
  print(knitr::kable(real_site$compare_summary))
  cat("\n")
  print(knitr::kable(real_site$duck_status))
  cat("\n")
  print(knitr::kable(real_site$mismatches))
  cat("\n")

  if (real_split) {
    real_split_case <- run_case(real_fixture$input_path, real_fixture$fasta_path, split = TRUE, n = bench_runs)
    cat("## ", real_label, " (split_multiallelic)\n\n", sep = "")
    print(knitr::kable(real_split_case$bench, digits = 3))
    cat("\n")
    print(knitr::kable(real_split_case$compare_summary))
    cat("\n")
    print(knitr::kable(real_split_case$duck_status))
    cat("\n")
    print(knitr::kable(real_split_case$mismatches))
    cat("\n")
  }
}
```

Real-callset benchmark skipped. Set `NORM_REAL_VCF` and
`NORM_REAL_FASTA` to enable it.
