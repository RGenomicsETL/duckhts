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
- uses `bcftools norm -Ou` plus `bcftools query` so output compression
  does not dominate the upstream timing

# Run

`make bench-norm`

Useful overrides:

- `NORM_BENCH_ROWS`: synthetic VCF records, default `100000`
- `NORM_BENCH_RUNS`: timed repeats, default `3`
- `BCFTOOLS_BIN`: optional override for the bcftools executable
- `NORM_DUCK_THREADS`: DuckDB execution threads for DuckHTS, default `1`
- `NORM_DUCK_DECOMPRESSION_THREADS`: `read_bcf(...)` htslib worker
  threads, default `0`
- `NORM_BCFTOOLS_THREADS`: optional `bcftools norm --threads`, default
  `0`
- `NORM_REAL_VCF`, `NORM_REAL_REGION`, `NORM_REAL_FASTA`: optional
  real-callset conformance case
- `NORM_REAL_SPLIT`: `1` to run the optional real-callset case in split
  mode too
- `NORM_REAL_LABEL`: optional label for the real-callset section

For a shell-first comparison, see
`bash scripts/norm_conformance.sh ...`.

# Settings

``` r
knitr::kable(benchmark_settings)
```

| setting                        | value |
|:-------------------------------|:------|
| duckdb_threads                 | 1     |
| read_bcf_decompression_threads | 0     |
| bcftools_threads               | 0     |
| bcftools_output_mode           | -Ou   |

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
| duckhts       | site_preserving    |    3 |      0.151 |   0.150 |   0.156 |      100000 |
| bcftools_norm | site_preserving    |    3 |     70.449 |  67.853 |  74.957 |      100000 |
| duckhts       | split_multiallelic |    3 |      0.207 |   0.207 |   0.207 |      120000 |
| bcftools_norm | split_multiallelic |    3 |    101.604 | 101.530 | 102.502 |      120000 |

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

## giab_hg001_grch37_full (site-preserving)

| engine        | mode            | runs | median_sec | min_sec | max_sec | output_rows |
|:--------------|:----------------|-----:|-----------:|--------:|--------:|------------:|
| duckhts       | site_preserving |    3 |     15.405 |  15.256 |  16.009 |     3891440 |
| bcftools_norm | site_preserving |    3 |     34.662 |  34.575 |  39.775 |     3891440 |

| status |       n |
|:-------|--------:|
| match  | 3891439 |

| norm_status      |       n |
|:-----------------|--------:|
| Normalized       |      43 |
| SpanningDeletion |     345 |
| Unchanged        | 3891052 |

| chrom | pos_normed | end_pos_normed | ref_normed | alt_normed | duck_n | bcf_n | status |
|:------|-----------:|---------------:|:-----------|:-----------|-------:|------:|:-------|

## giab_hg001_grch37_full (split_multiallelic)

| engine        | mode               | runs | median_sec | min_sec | max_sec | output_rows |
|:--------------|:-------------------|-----:|-----------:|--------:|--------:|------------:|
| duckhts       | split_multiallelic |    3 |     23.607 |  23.153 |  24.313 |     3933714 |
| bcftools_norm | split_multiallelic |    3 |     32.456 |  32.431 |  36.844 |     3933714 |

| status |       n |
|:-------|--------:|
| match  | 3933713 |

| norm_status      |       n |
|:-----------------|--------:|
| Normalized       |   26837 |
| SpanningDeletion |     345 |
| Unchanged        | 3906532 |

| chrom | pos_normed | end_pos_normed | ref_normed | alt_normed | duck_n | bcf_n | status |
|:------|-----------:|---------------:|:-----------|:-----------|-------:|------:|:-------|
