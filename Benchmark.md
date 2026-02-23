DuckHTS Benchmark
================

<!-- Benchmark.md is generated from Benchmark.Rmd. -->

# Goal

Benchmark `read_bcf()` performance on realistic files (ClinVar and
larger VEP-annotated VCF/BCF), with focus on:

- Full scan throughput
- Projection pushdown impact
- INFO and FORMAT parsing cost
- Tidy vs wide FORMAT output
- Region query performance

# Setup

``` r
library(DBI)
library(duckdb)

drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
con <- dbConnect(drv, dbdir = ":memory:")

# Load extension from local build output (adjust if needed)
ext_path <- normalizePath("build/release/duckhts.duckdb_extension", mustWork = TRUE)
# DuckDB may block unsigned local extensions by default.
try(DBI::dbExecute(con, "SET allow_unsigned_extensions=true;"), silent = TRUE)
DBI::dbExecute(con, sprintf("LOAD '%s';", gsub("\\\\", "/", ext_path)))
#> [1] 0

# Optional: control parallelism for reproducibility
DBI::dbExecute(con, "PRAGMA threads=4;")
#> [1] 0

sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Berlin
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices datasets  utils     methods   base     
#> 
#> other attached packages:
#> [1] duckdb_1.4.3 DBI_1.2.3   
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39   fastmap_1.2.0   xfun_0.56       glue_1.8.0     
#>  [5] bspm_0.5.7      knitr_1.51      htmltools_0.5.9 rmarkdown_2.30 
#>  [9] lifecycle_1.0.5 cli_3.6.5       vctrs_0.7.1     compiler_4.5.2 
#> [13] tools_4.5.2     evaluate_1.0.5  pillar_1.11.1   yaml_2.3.12    
#> [17] otel_0.2.0      rlang_1.1.7
```

# File Paths

``` r
clinvar_vcf <- "clinvar.vcf.gz"
vep_vcf <- Sys.getenv("VEP_VCF", unset = "")

stopifnot(file.exists(clinvar_vcf))
has_vep <- nzchar(vep_vcf) && file.exists(vep_vcf)
```

# Helpers

``` r
build_read_bcf <- function(path, ..., tidy = FALSE) {
  named <- list(...)
  args <- c(sprintf("'%s'", gsub("'", "''", path)))
  if (length(named) > 0) {
    kv <- vapply(names(named), function(k) {
      v <- named[[k]]
      if (is.character(v)) sprintf("%s := '%s'", k, gsub("'", "''", v))
      else if (isTRUE(v) || identical(v, FALSE)) sprintf("%s := %s", k, tolower(as.character(v)))
      else sprintf("%s := %s", k, as.character(v))
    }, FUN.VALUE = character(1))
    args <- c(args, kv)
  }
  if (isTRUE(tidy)) args <- c(args, "tidy_format := true")
  sprintf("read_bcf(%s)", paste(args, collapse = ", "))
}

run_bench <- function(con, name, sql, iterations = 5, warmup = 1) {
  cat("\\n---\\n", name, "\\n", sep = "")
  cat(sql, "\\n")

  # Warmup
  for (i in seq_len(warmup)) DBI::dbGetQuery(con, sql)

  times <- numeric(iterations)
  rows <- integer(iterations)
  for (i in seq_len(iterations)) {
    gc()
    t0 <- proc.time()[["elapsed"]]
    out <- DBI::dbGetQuery(con, sql)
    t1 <- proc.time()[["elapsed"]]
    times[i] <- t1 - t0
    rows[i] <- if (is.data.frame(out)) nrow(out) else NA_integer_
  }

  data.frame(
    case = name,
    iterations = iterations,
    min_sec = min(times),
    median_sec = stats::median(times),
    mean_sec = mean(times),
    max_sec = max(times),
    rows = stats::median(rows)
  )
}

get_bcf_columns <- function(con, path, tidy = FALSE) {
  src <- build_read_bcf(path, tidy = tidy)
  q <- sprintf("DESCRIBE SELECT * FROM %s", src)
  DBI::dbGetQuery(con, q)$column_name
}
```

# ClinVar Benchmarks

``` r
clinvar_src <- build_read_bcf(clinvar_vcf)
clinvar_cols <- get_bcf_columns(con, clinvar_vcf)

cases <- list(
  list(
    name = "clinvar_count_all",
    sql = sprintf("SELECT COUNT(*) AS n FROM %s", clinvar_src)
  ),
  list(
    name = "clinvar_core_projection",
    sql = sprintf(
      "SELECT CHROM, POS, REF, ALT FROM %s WHERE POS > 0 LIMIT 200000",
      clinvar_src
    )
  )
)

# Add an INFO-heavy case if available
info_cols <- grep("^INFO_", clinvar_cols, value = TRUE)
if (length(info_cols) > 0) {
  pick <- paste(head(info_cols, 6), collapse = ", ")
  cases[[length(cases) + 1]] <- list(
    name = "clinvar_info_projection",
    sql = sprintf("SELECT %s FROM %s LIMIT 200000", pick, clinvar_src)
  )
}

# Region query case (update region to contig naming style in your file)
cases[[length(cases) + 1]] <- list(
  name = "clinvar_region_count",
  sql = sprintf("SELECT COUNT(*) AS n FROM %s", build_read_bcf(clinvar_vcf, region = "chr1:1-5000000"))
)

clinvar_results <- do.call(
  rbind,
  lapply(cases, function(x) run_bench(con, x$name, x$sql, iterations = 5, warmup = 1))
)
#> \n---\nclinvar_count_all\nSELECT COUNT(*) AS n FROM read_bcf('clinvar.vcf.gz') \n\n---\nclinvar_core_projection\nSELECT CHROM, POS, REF, ALT FROM read_bcf('clinvar.vcf.gz') WHERE POS > 0 LIMIT 200000 \n\n---\nclinvar_info_projection\nSELECT INFO_AF_ESP, INFO_AF_EXAC, INFO_AF_TGP, INFO_ALLELEID, INFO_CLNDN, INFO_CLNDNINCL FROM read_bcf('clinvar.vcf.gz') LIMIT 200000 \n\n---\nclinvar_region_count\nSELECT COUNT(*) AS n FROM read_bcf('clinvar.vcf.gz', region := 'chr1:1-5000000') \n
clinvar_results
#>                      case iterations min_sec median_sec mean_sec max_sec   rows
#> 1       clinvar_count_all          5   1.250      1.254   1.2582   1.274      1
#> 2 clinvar_core_projection          5   0.256      0.257   0.2586   0.265 200000
#> 3 clinvar_info_projection          5   0.295      0.297   0.2970   0.301 200000
#> 4    clinvar_region_count          5   0.020      0.020   0.0206   0.022      1
```

# VEP Benchmarks

``` r
if (has_vep) {
  vep_src <- build_read_bcf(vep_vcf)
  vep_tidy_src <- build_read_bcf(vep_vcf, tidy = TRUE)
  vep_cols <- get_bcf_columns(con, vep_vcf)
  vep_tidy_cols <- get_bcf_columns(con, vep_vcf, tidy = TRUE)

  vep_cases <- list(
    list(
      name = "vep_count_all",
      sql = sprintf("SELECT COUNT(*) AS n FROM %s", vep_src)
    )
  )

  vep_ann_cols <- grep("^VEP_", vep_cols, value = TRUE)
  if (length(vep_ann_cols) > 0) {
    pick_vep <- paste(head(vep_ann_cols, 8), collapse = ", ")
    vep_cases[[length(vep_cases) + 1]] <- list(
      name = "vep_annotation_projection",
      sql = sprintf("SELECT %s FROM %s LIMIT 200000", pick_vep, vep_src)
    )
  }

  fmt_cols_wide <- grep("^FORMAT_", vep_cols, value = TRUE)
  fmt_cols_tidy <- grep("^FORMAT_", vep_tidy_cols, value = TRUE)
  if (length(fmt_cols_wide) > 0 && length(fmt_cols_tidy) > 0) {
    pick_wide <- paste(head(fmt_cols_wide, 8), collapse = ", ")
    pick_tidy <- paste(head(fmt_cols_tidy, 8), collapse = ", ")

    vep_cases[[length(vep_cases) + 1]] <- list(
      name = "vep_format_wide_projection",
      sql = sprintf("SELECT %s FROM %s LIMIT 100000", pick_wide, vep_src)
    )
    vep_cases[[length(vep_cases) + 1]] <- list(
      name = "vep_format_tidy_projection",
      sql = sprintf("SELECT SAMPLE_ID, %s FROM %s LIMIT 100000", pick_tidy, vep_tidy_src)
    )
  }

  vep_results <- do.call(
    rbind,
    lapply(vep_cases, function(x) run_bench(con, x$name, x$sql, iterations = 5, warmup = 1))
  )
  vep_results
} else {
  vep_results <- data.frame()
  message("Skipping VEP benchmarks. Set env var VEP_VCF to a local VEP VCF/BCF file.")
}
#> Skipping VEP benchmarks. Set env var VEP_VCF to a local VEP VCF/BCF file.
```

# Optional: Export Results

``` r
out <- rbind(clinvar_results, vep_results)
write.csv(out, "benchmark_results.csv", row.names = FALSE)
out
#>                      case iterations min_sec median_sec mean_sec max_sec   rows
#> 1       clinvar_count_all          5   1.250      1.254   1.2582   1.274      1
#> 2 clinvar_core_projection          5   0.256      0.257   0.2586   0.265 200000
#> 3 clinvar_info_projection          5   0.295      0.297   0.2970   0.301 200000
#> 4    clinvar_region_count          5   0.020      0.020   0.0206   0.022      1
```

# Notes For Reproducibility

- Record CPU model, core count, RAM, and storage type
  (NVMe/SATA/network).
- Keep `PRAGMA threads` fixed across runs.
- Run each case multiple times and report median.
- Prefer local files for stable throughput comparisons.

# Larger VEP VCF Candidates

Good options for larger VEP-annotated files:

- Publicly released VCF/BCF from projects that include `INFO/CSQ` (or
  `ANN`/`BCSQ`) in headers.
- Any cohort VCF you already have can be VEP-annotated with `vep`
  (offline cache mode) to generate a controlled benchmark input.
- Keep one “wide FORMAT” dataset (many samples) and one “deep
  annotation” dataset (many VEP subfields) to separate bottlenecks.

Quick check that a file has VEP-like annotation fields:

``` r
if (has_vep) {
  q <- sprintf(
    "SELECT DISTINCT id FROM read_hts_header('%s') WHERE record_type = 'INFO' AND id IN ('CSQ','ANN','BCSQ')",
    gsub("'", "''", vep_vcf)
  )
  DBI::dbGetQuery(con, q)
}
```
