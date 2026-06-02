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
- `NORM_GVCF_VCF`, `NORM_GVCF_FASTA`, `NORM_GVCF_LABEL`: optional
  locally staged gVCF slice benchmark; the current local default is a 10
  Mb 1000 Genomes DRAGEN HG00096 chr22 slice when present
- `NORM_GVCF_RUNS`: timed repeats for the optional gVCF slice, default
  `1`
- `NORM_GVCF_SPLIT`: `1` to also run split-multiallelic mode for the
  optional gVCF slice
- `NORM_GIAB_DIR`, `NORM_GIAB_FASTA`, `NORM_GIAB_SAMPLES`: full GIAB
  benchmark inputs; defaults target locally staged HG001/HG002/HG006
  GRCh37 files
- `NORM_GIAB_RUNS`: timed repeats for each full GIAB sample, default `1`
- `NORM_GIAB_SPLIT`: `1` to also run split-multiallelic mode for full
  GIAB inputs
- `NORM_REQUIRE_GIAB`: `1` to fail if any configured full GIAB input is
  missing

For a shell-first comparison, see
`bash scripts/norm_conformance.sh ...`.

# Settings

``` r
knitr::kable(benchmark_settings)
```

| setting                        | value                                                                            |
|:-------------------------------|:---------------------------------------------------------------------------------|
| duckdb_threads                 | 1                                                                                |
| read_bcf_decompression_threads | 0                                                                                |
| bcftools_threads               | 0                                                                                |
| bcftools_output_mode           | -Ou                                                                              |
| gvcf_input_path                | .tmp/1000g_dragen_hg00096/HG00096.hard-filtered.chr22_20000000_30000000.g.vcf.gz |
| gvcf_fasta                     | /root/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna                     |
| gvcf_runs                      | 1                                                                                |
| giab_dir                       | /root/giab_norm                                                                  |
| giab_fasta                     | /root/GRCh37/human_g1k_v37.fasta                                                 |
| giab_samples                   | HG001,HG002,HG006                                                                |
| giab_runs                      | 1                                                                                |

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
| duckhts       | site_preserving    |    3 |      0.193 |   0.190 |   0.200 |      100000 |
| bcftools_norm | site_preserving    |    3 |      0.604 |   0.595 |   0.728 |      100000 |
| duckhts       | split_multiallelic |    3 |      0.374 |   0.373 |   0.379 |      120000 |
| bcftools_norm | split_multiallelic |    3 |      0.776 |   0.755 |   0.792 |      120000 |

Synthetic fixture mix:

- normalized insertions/deletions
- a multiallelic `TT,TTT` row
- a ref-only `ALT='.'` row that exercises split-mode row preservation
- one variant per repeated reference motif, so the synthetic benchmark
  avoids same-position duplicate-path behavior and should be read as a
  compact kernel/conformance workload rather than a real cohort speed
  claim

## Synthetic conformance

### Site-preserving mode

``` r
knitr::kable(site_case$compare_summary)
```

| status |     n |
|:-------|------:|
| match  | 1e+05 |

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

| status |      n |
|:-------|-------:|
| match  | 120000 |

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

# Optional public gVCF slice benchmark

This section is non-network-dependent during rendering. Stage a local
slice first with:

``` sh
make stage-norm-1000g-dragen-gvcf
```

``` r
if (!file.exists(gvcf_input_path) || !file.exists(gvcf_fasta)) {
  cat("gVCF slice benchmark skipped. Set `NORM_GVCF_VCF` and `NORM_GVCF_FASTA`, or stage the default 1000G DRAGEN slice locally.\n")
} else {
  gvcf_case <- run_case(gvcf_input_path, gvcf_fasta, split = FALSE, n = gvcf_runs)
  cat("## ", gvcf_label, " (site-preserving)\n\n", sep = "")
  print(knitr::kable(gvcf_case$bench, digits = 3))
  cat("\n")
  print(knitr::kable(gvcf_case$compare_summary))
  cat("\n")
  print(knitr::kable(gvcf_case$duck_status))
  cat("\n")
  print(knitr::kable(gvcf_case$mismatches))
  cat("\n")

  if (gvcf_split) {
    gvcf_split_case <- run_case(gvcf_input_path, gvcf_fasta, split = TRUE, n = gvcf_runs)
    cat("## ", gvcf_label, " (split_multiallelic)\n\n", sep = "")
    print(knitr::kable(gvcf_split_case$bench, digits = 3))
    cat("\n")
    print(knitr::kable(gvcf_split_case$compare_summary))
    cat("\n")
    print(knitr::kable(gvcf_split_case$duck_status))
    cat("\n")
    print(knitr::kable(gvcf_split_case$mismatches))
    cat("\n")
  }
}
```

## 1000G_DRAGEN_HG00096_chr22_20_30Mb_gVCF (site-preserving)

| engine        | mode            | runs | median_sec | min_sec | max_sec | output_rows |
|:--------------|:----------------|-----:|-----------:|--------:|--------:|------------:|
| duckhts       | site_preserving |    1 |      1.458 |   1.458 |   1.458 |      677444 |
| bcftools_norm | site_preserving |    1 |      6.043 |   6.043 |   6.043 |      677444 |

| status |      n |
|:-------|-------:|
| match  | 677444 |

| norm_status        |      n |
|:-------------------|-------:|
| GVCFReferenceBlock | 659313 |
| Unchanged          |  18131 |

| chrom | pos_normed | end_pos_normed | ref_normed | alt_normed | duck_n | bcf_n | status |
|:------|-----------:|---------------:|:-----------|:-----------|-------:|------:|:-------|

# Full GIAB sample benchmarks

This section benchmarks full locally staged GIAB benchmark VCFs without
slicing. DuckHTS counts normalized rows and status classes through
`duckhts_bcftools_norm(...)`; `bcftools norm` is timed with
`-Ou > /dev/null` so output compression does not dominate.

``` r
giab_inputs <- if (length(giab_samples)) {
  data.frame(
    sample = giab_samples,
    path = file.path(giab_dir, paste0(giab_samples, "_GRCh37_1_22_v4.2.1_benchmark.vcf.gz")),
    stringsAsFactors = FALSE
  )
} else {
  data.frame(sample = character(), path = character(), stringsAsFactors = FALSE)
}
giab_inputs$exists <- file.exists(giab_inputs$path)
giab_ready <- length(giab_samples) > 0 && file.exists(giab_fasta) && all(giab_inputs$exists)

if (!giab_ready) {
  missing <- c(giab_inputs$path[!giab_inputs$exists], if (!file.exists(giab_fasta)) giab_fasta else character())
  msg <- if (!length(giab_samples)) {
    "Full GIAB benchmark skipped; set `NORM_GIAB_SAMPLES` to one or more staged sample IDs to enable it."
  } else {
    paste("Full GIAB benchmark skipped; missing:", paste(missing, collapse = ", "))
  }
  if (require_giab) stop(msg, call. = FALSE)
  cat(msg, "\n")
} else {
  giab_bench <- list()
  giab_status <- list()
  for (i in seq_len(nrow(giab_inputs))) {
    label <- giab_inputs$sample[[i]]
    cat("## ", label, " full GIAB (site-preserving)\n\n", sep = "")
    case <- run_giab_bench_case(label, giab_inputs$path[[i]], giab_fasta, split = FALSE, n = giab_runs)
    giab_bench[[length(giab_bench) + 1L]] <- case$bench
    giab_status[[label]] <- case$status
    print(knitr::kable(case$bench, digits = 3))
    cat("\n")
    print(knitr::kable(case$status))
    cat("\n")

    if (giab_split) {
      cat("## ", label, " full GIAB (split_multiallelic)\n\n", sep = "")
      split_case <- run_giab_bench_case(paste0(label, "_split"), giab_inputs$path[[i]], giab_fasta, split = TRUE, n = giab_runs)
      giab_bench[[length(giab_bench) + 1L]] <- split_case$bench
      print(knitr::kable(split_case$bench, digits = 3))
      cat("\n")
      print(knitr::kable(split_case$status))
      cat("\n")
    }
  }

  cat("## Full GIAB summary\n\n")
  print(knitr::kable(do.call(rbind, giab_bench), digits = 3))
}
```

## HG001 full GIAB (site-preserving)

| case  | engine        | mode            | runs | median_sec | min_sec | max_sec | output_rows | normalized_rows |
|:------|:--------------|:----------------|-----:|-----------:|--------:|--------:|------------:|----------------:|
| HG001 | duckhts       | site_preserving |    1 |      8.408 |   8.408 |   8.408 |     3891440 |              43 |
| HG001 | bcftools_norm | site_preserving |    1 |      4.960 |   4.960 |   4.960 |          NA |              NA |

| norm_status |       n |
|:------------|--------:|
| Normalized  |      43 |
| Unchanged   | 3891397 |

## HG002 full GIAB (site-preserving)

| case  | engine        | mode            | runs | median_sec | min_sec | max_sec | output_rows | normalized_rows |
|:------|:--------------|:----------------|-----:|-----------:|--------:|--------:|------------:|----------------:|
| HG002 | duckhts       | site_preserving |    1 |      9.366 |   9.366 |   9.366 |     4033796 |              18 |
| HG002 | bcftools_norm | site_preserving |    1 |      5.556 |   5.556 |   5.556 |          NA |              NA |

| norm_status |       n |
|:------------|--------:|
| Normalized  |      18 |
| Unchanged   | 4033778 |

## HG006 full GIAB (site-preserving)

| case  | engine        | mode            | runs | median_sec | min_sec | max_sec | output_rows | normalized_rows |
|:------|:--------------|:----------------|-----:|-----------:|--------:|--------:|------------:|----------------:|
| HG006 | duckhts       | site_preserving |    1 |      8.998 |   8.998 |   8.998 |     3878664 |              18 |
| HG006 | bcftools_norm | site_preserving |    1 |      4.718 |   4.718 |   4.718 |          NA |              NA |

| norm_status |       n |
|:------------|--------:|
| Normalized  |      18 |
| Unchanged   | 3878646 |

## Full GIAB summary

| case  | engine        | mode            | runs | median_sec | min_sec | max_sec | output_rows | normalized_rows |
|:------|:--------------|:----------------|-----:|-----------:|--------:|--------:|------------:|----------------:|
| HG001 | duckhts       | site_preserving |    1 |      8.408 |   8.408 |   8.408 |     3891440 |              43 |
| HG001 | bcftools_norm | site_preserving |    1 |      4.960 |   4.960 |   4.960 |          NA |              NA |
| HG002 | duckhts       | site_preserving |    1 |      9.366 |   9.366 |   9.366 |     4033796 |              18 |
| HG002 | bcftools_norm | site_preserving |    1 |      5.556 |   5.556 |   5.556 |          NA |              NA |
| HG006 | duckhts       | site_preserving |    1 |      8.998 |   8.998 |   8.998 |     3878664 |              18 |
| HG006 | bcftools_norm | site_preserving |    1 |      4.718 |   4.718 |   4.718 |          NA |              NA |
