DuckHTS Benchmark
================

<!-- Benchmark.md is generated from Benchmark.Rmd. -->

# Goal

Benchmark DuckHTS with a small set of “napkin numbers” that answer:

- How fast can we scan a file end-to-end?
- How much slower are richer projections than a pure count?
- How much extra cost comes from conversion to Parquet?
- How close are we to a simple external baseline?
- How do DuckHTS and Exon compare on the same workload shape?

The benchmark structure is:

1.  Baseline: file size and simple throughput calculations
2.  DuckHTS scan/query benchmarks
3.  DuckHTS conversion benchmarks
4.  Optional external comparison points (for example Exon and bcftools)

# Setup

``` r
library(DBI)
library(duckdb)
library(tools)

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
#> [1] tools     stats     graphics  grDevices datasets  utils     methods  
#> [8] base     
#> 
#> other attached packages:
#> [1] duckdb_1.4.3 DBI_1.2.3   
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39   fastmap_1.2.0   xfun_0.56       glue_1.8.0     
#>  [5] bspm_0.5.7      knitr_1.51      htmltools_0.5.9 rmarkdown_2.30 
#>  [9] lifecycle_1.0.5 cli_3.6.5       vctrs_0.7.1     compiler_4.5.2 
#> [13] evaluate_1.0.5  pillar_1.11.1   yaml_2.3.12     otel_0.2.0     
#> [17] rlang_1.1.7
```

# File Paths

``` r
clinvar_vcf <- "clinvar.vcf.gz"
vep_vcf <- Sys.getenv("VEP_VCF", unset = "")
bam_path <- "HG00106.chrom11.ILLUMINA.bwa.GBR.exome.20130415.bam"
bam_index_path <- paste0(bam_path, ".bai")
bench_dir <- normalizePath(tempdir(), mustWork = TRUE)

stopifnot(file.exists(clinvar_vcf))
has_vep <- nzchar(vep_vcf) && file.exists(vep_vcf)
has_bam <- file.exists(bam_path) && file.exists(bam_index_path)

clinvar_bytes <- file.info(clinvar_vcf)$size
vep_bytes <- if (has_vep) file.info(vep_vcf)$size else NA_real_
bam_bytes <- if (has_bam) file.info(bam_path)$size else NA_real_
bcftools_bin <- Sys.which("bcftools")
has_bcftools <- nzchar(bcftools_bin)
samtools_bin <- Sys.which("samtools")
has_samtools <- nzchar(samtools_bin)
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

build_read_bam <- function(path, ...) {
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
  sprintf("read_bam(%s)", paste(args, collapse = ", "))
}

run_bench <- function(con, name, sql, iterations = 5, warmup = 1, bytes_read = NA_real_, rows_hint = NA_real_) {
  message(paste0("\n---\n", name, "\n", sql, "\n"))

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
    rows = stats::median(rows),
    bytes_read = bytes_read,
    compressed_mb = bytes_read / (1024 * 1024),
    compressed_mb_per_sec = if (is.finite(bytes_read)) (bytes_read / (1024 * 1024)) / stats::median(times) else NA_real_,
    rows_hint = rows_hint,
    rows_per_sec = if (is.finite(rows_hint)) rows_hint / stats::median(times) else NA_real_
  )
}

get_bcf_columns <- function(con, path, tidy = FALSE) {
  src <- build_read_bcf(path, tidy = tidy)
  q <- sprintf("DESCRIBE SELECT * FROM %s", src)
  DBI::dbGetQuery(con, q)$column_name
}

get_variant_count <- function(con, path, ...) {
  src <- build_read_bcf(path, ...)
  DBI::dbGetQuery(con, sprintf("SELECT COUNT(*) AS n FROM %s", src))$n[[1]]
}

get_header_ids <- function(con, path, record_type) {
  q <- sprintf(
    "SELECT id FROM read_hts_header('%s') WHERE record_type = '%s'",
    gsub("'", "''", path),
    gsub("'", "''", record_type)
  )
  DBI::dbGetQuery(con, q)$id
}

run_bcftools_bench <- function(name, args, iterations = 5, warmup = 1, bytes_read = NA_real_, rows_hint = NA_real_) {
  if (!has_bcftools) {
    return(data.frame())
  }

  shell_quote <- function(x) {
    paste0("'", gsub("'", "'\"'\"'", x, fixed = TRUE), "'")
  }

  cmd_display <- paste(shQuote(bcftools_bin), paste(vapply(args, shell_quote, character(1)), collapse = " "))
  message(paste0("\n---\n", name, "\n", cmd_display, "\n"))

  run_once <- function() {
    out_file <- tempfile("bcftools-out-")
    err_file <- tempfile("bcftools-err-")
    on.exit(unlink(c(out_file, err_file), force = TRUE), add = TRUE)

    status <- system(
      sprintf("%s > %s 2> %s", cmd_display, shQuote(out_file), shQuote(err_file)),
      ignore.stdout = TRUE,
      ignore.stderr = TRUE
    )
    if (!identical(status, 0L)) {
      err_lines <- readLines(err_file, warn = FALSE)
      stop(
        sprintf(
          "bcftools command failed for %s (exit %s)%s",
          name,
          status,
          if (length(err_lines)) paste0(": ", paste(err_lines, collapse = " ")) else ""
        ),
        call. = FALSE
      )
    }
  }

  for (i in seq_len(warmup)) {
    run_once()
  }

  times <- numeric(iterations)
  for (i in seq_len(iterations)) {
    gc()
    t0 <- proc.time()[["elapsed"]]
    run_once()
    t1 <- proc.time()[["elapsed"]]
    times[i] <- t1 - t0
  }

  data.frame(
    engine = "bcftools",
    case = name,
    iterations = iterations,
    min_sec = min(times),
    median_sec = stats::median(times),
    mean_sec = mean(times),
    max_sec = max(times),
    rows = rows_hint,
    bytes_read = bytes_read,
    compressed_mb = bytes_read / (1024 * 1024),
    compressed_mb_per_sec = if (is.finite(bytes_read)) (bytes_read / (1024 * 1024)) / stats::median(times) else NA_real_,
    rows_hint = rows_hint,
    rows_per_sec = if (is.finite(rows_hint)) rows_hint / stats::median(times) else NA_real_
  )
}

run_samtools_bench <- function(name, args, iterations = 5, warmup = 1, bytes_read = NA_real_, rows_hint = NA_real_) {
  if (!has_samtools) {
    return(data.frame())
  }

  cmd_display <- paste(shQuote(samtools_bin), paste(args, collapse = " "))
  message(paste0("\n---\n", name, "\n", cmd_display, "\n"))

  run_once <- function() {
    out_file <- tempfile("samtools-out-")
    err_file <- tempfile("samtools-err-")
    on.exit(unlink(c(out_file, err_file), force = TRUE), add = TRUE)

    status <- system2(
      samtools_bin,
      args = args,
      stdout = out_file,
      stderr = err_file
    )
    if (!identical(status, 0L)) {
      err_lines <- readLines(err_file, warn = FALSE)
      stop(
        sprintf(
          "samtools command failed for %s (exit %s)%s",
          name,
          status,
          if (length(err_lines)) paste0(": ", paste(err_lines, collapse = " ")) else ""
        ),
        call. = FALSE
      )
    }
  }

  for (i in seq_len(warmup)) {
    run_once()
  }

  times <- numeric(iterations)
  for (i in seq_len(iterations)) {
    gc()
    t0 <- proc.time()[["elapsed"]]
    run_once()
    t1 <- proc.time()[["elapsed"]]
    times[i] <- t1 - t0
  }

  data.frame(
    engine = "samtools",
    case = name,
    iterations = iterations,
    min_sec = min(times),
    median_sec = stats::median(times),
    mean_sec = mean(times),
    max_sec = max(times),
    rows = rows_hint,
    bytes_read = bytes_read,
    compressed_mb = bytes_read / (1024 * 1024),
    compressed_mb_per_sec = if (is.finite(bytes_read)) (bytes_read / (1024 * 1024)) / stats::median(times) else NA_real_,
    rows_hint = rows_hint,
    rows_per_sec = if (is.finite(rows_hint)) rows_hint / stats::median(times) else NA_real_
  )
}

make_case_key <- function(case_name) {
  key <- sub("^clinvar_", "", case_name)
  key <- sub("^vep_", "", key)
  key <- sub("^bcftools_", "", key)
  key
}

add_case_keys <- function(df) {
  if (!is.data.frame(df) || nrow(df) == 0) return(df)
  df$case_key <- vapply(df$case, make_case_key, character(1))
  df
}

run_copy_bench <- function(con, name, select_sql, out_path, iterations = 3, warmup = 1, bytes_read = NA_real_, rows_hint = NA_real_) {
  copy_sql <- sprintf(
    "COPY (%s) TO '%s' (FORMAT parquet, COMPRESSION zstd)",
    select_sql,
    gsub("\\\\", "/", out_path)
  )

  for (i in seq_len(warmup)) {
    if (file.exists(out_path)) file.remove(out_path)
    DBI::dbExecute(con, copy_sql)
  }

  times <- numeric(iterations)
  out_sizes <- numeric(iterations)
  for (i in seq_len(iterations)) {
    if (file.exists(out_path)) file.remove(out_path)
    gc()
    t0 <- proc.time()[["elapsed"]]
    DBI::dbExecute(con, copy_sql)
    t1 <- proc.time()[["elapsed"]]
    times[i] <- t1 - t0
    out_sizes[i] <- if (file.exists(out_path)) file.info(out_path)$size else NA_real_
  }

  data.frame(
    case = name,
    iterations = iterations,
    min_sec = min(times),
    median_sec = stats::median(times),
    mean_sec = mean(times),
    max_sec = max(times),
    rows = rows_hint,
    bytes_read = bytes_read,
    compressed_mb = bytes_read / (1024 * 1024),
    compressed_mb_per_sec = if (is.finite(bytes_read)) (bytes_read / (1024 * 1024)) / stats::median(times) else NA_real_,
    rows_hint = rows_hint,
    rows_per_sec = if (is.finite(rows_hint)) rows_hint / stats::median(times) else NA_real_,
    parquet_bytes = stats::median(out_sizes),
    parquet_mb = stats::median(out_sizes) / (1024 * 1024),
    parquet_mb_per_sec = (stats::median(out_sizes) / (1024 * 1024)) / stats::median(times)
  )
}

label_workload <- function(case_name) {
  if (grepl("copy_.*parquet$", case_name)) return("conversion")
  if (grepl("region_", case_name)) return("region")
  if (grepl("core_projection|info_projection|format_|annotation_projection", case_name)) return("projection")
  if (grepl("count_all$", case_name)) return("full_scan")
  "other"
}

add_relative_metrics <- function(df, baseline_case) {
  ref <- df[df$case == baseline_case, , drop = FALSE]
  if (nrow(ref) != 1) return(df)

  ref_scan_mb_s <- ref$compressed_mb_per_sec[[1]]
  ref_rows_s <- ref$rows_per_sec[[1]]

  df$vs_scan_baseline_mb <- df$compressed_mb_per_sec / ref_scan_mb_s
  df$vs_scan_baseline_rows <- df$rows_per_sec / ref_rows_s
  df
}

bind_result_frames <- function(...) {
  frames <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, list(...))
  if (length(frames) == 0) return(data.frame())

  all_cols <- unique(unlist(lapply(frames, names), use.names = FALSE))
  aligned <- lapply(frames, function(df) {
    missing <- setdiff(all_cols, names(df))
    for (nm in missing) df[[nm]] <- NA
    df[all_cols]
  })
  do.call(rbind, aligned)
}

finalize_results <- function(df) {
  if (!is.data.frame(df) || nrow(df) == 0) return(df)
  df$workload <- vapply(df$case, label_workload, character(1))
  for (i in seq_len(nrow(df))) {
    if (df$workload[[i]] != "full_scan") {
      df$compressed_mb[[i]] <- NA_real_
      df$compressed_mb_per_sec[[i]] <- NA_real_
    }
  }
  df
}
```

# Baseline Dataset Facts

``` r
data.frame(
  dataset = c("clinvar", if (has_vep) "vep" else NULL, if (has_bam) "bam" else NULL),
  path = c(clinvar_vcf, if (has_vep) vep_vcf else NULL, if (has_bam) bam_path else NULL),
  bytes = c(clinvar_bytes, if (has_vep) vep_bytes else NULL, if (has_bam) bam_bytes else NULL),
  compressed_mb = c(clinvar_bytes, if (has_vep) vep_bytes else NULL, if (has_bam) bam_bytes else NULL) / (1024 * 1024)
)
#>   dataset                                                path     bytes
#> 1 clinvar                                      clinvar.vcf.gz 189031261
#> 2     bam HG00106.chrom11.ILLUMINA.bwa.GBR.exome.20130415.bam 330528619
#>   compressed_mb
#> 1      180.2743
#> 2      315.2167
```

# DuckHTS Scan Benchmarks

``` r
clinvar_src <- build_read_bcf(clinvar_vcf)
clinvar_cols <- get_bcf_columns(con, clinvar_vcf)
clinvar_variants <- get_variant_count(con, clinvar_vcf)
clinvar_contigs <- get_header_ids(con, clinvar_vcf, "CONTIG")
clinvar_region <- NA_character_
if ("chr1" %in% clinvar_contigs) clinvar_region <- "chr1:1-5000000"
if (is.na(clinvar_region) && "1" %in% clinvar_contigs) clinvar_region <- "1:1-5000000"

cases <- list(
  list(
    name = "clinvar_count_all",
    sql = sprintf("SELECT COUNT(*) AS n FROM %s", clinvar_src),
    rows_hint = clinvar_variants
  ),
  list(
    name = "clinvar_core_projection",
    sql = sprintf(
      "SELECT CHROM, POS, REF, ALT FROM %s WHERE POS > 0 LIMIT 200000",
      clinvar_src
    ),
    rows_hint = 200000
  )
)

# Add an INFO-heavy case if available
info_cols <- grep("^INFO_", clinvar_cols, value = TRUE)
if (length(info_cols) > 0) {
  pick <- paste(head(info_cols, 6), collapse = ", ")
  cases[[length(cases) + 1]] <- list(
    name = "clinvar_info_projection",
    sql = sprintf("SELECT %s FROM %s LIMIT 200000", pick, clinvar_src),
    rows_hint = 200000
  )
}

if (!is.na(clinvar_region)) {
  cases[[length(cases) + 1]] <- list(
    name = "clinvar_region_count",
    sql = sprintf("SELECT COUNT(*) AS n FROM %s", build_read_bcf(clinvar_vcf, region = clinvar_region)),
    rows_hint = NA_real_
  )
} else {
  message("Skipping ClinVar region benchmark: no chr1/1 contig found in header.")
}
#> Skipping ClinVar region benchmark: no chr1/1 contig found in header.

clinvar_results <- do.call(
  rbind,
  lapply(cases, function(x) run_bench(
    con,
    x$name,
    x$sql,
    iterations = 5,
    warmup = 1,
    bytes_read = clinvar_bytes,
    rows_hint = x$rows_hint
  ))
)
#> 
#> ---
#> clinvar_count_all
#> SELECT COUNT(*) AS n FROM read_bcf('clinvar.vcf.gz')
#> 
#> ---
#> clinvar_core_projection
#> SELECT CHROM, POS, REF, ALT FROM read_bcf('clinvar.vcf.gz') WHERE POS > 0 LIMIT 200000
#> 
#> ---
#> clinvar_info_projection
#> SELECT INFO_AF_ESP, INFO_AF_EXAC, INFO_AF_TGP, INFO_ALLELEID, INFO_CLNDN, INFO_CLNDNINCL FROM read_bcf('clinvar.vcf.gz') LIMIT 200000
clinvar_results$engine <- "duckhts"
clinvar_results <- add_relative_metrics(clinvar_results, "clinvar_count_all")
clinvar_results <- add_case_keys(clinvar_results)
clinvar_results <- finalize_results(clinvar_results)
clinvar_results
#>                      case iterations min_sec median_sec mean_sec max_sec   rows
#> 1       clinvar_count_all          5   1.119      1.127   1.1246   1.128      1
#> 2 clinvar_core_projection          5   0.243      0.244   0.2452   0.252 200000
#> 3 clinvar_info_projection          5   0.299      0.306   0.3060   0.314 200000
#>   bytes_read compressed_mb compressed_mb_per_sec rows_hint rows_per_sec  engine
#> 1  189031261      180.2743              159.9594   4352930    3862404.6 duckhts
#> 2  189031261            NA                    NA    200000     819672.1 duckhts
#> 3  189031261            NA                    NA    200000     653594.8 duckhts
#>   vs_scan_baseline_mb vs_scan_baseline_rows        case_key   workload
#> 1            1.000000             1.0000000       count_all  full_scan
#> 2            4.618852             0.2122181 core_projection projection
#> 3            3.683007             0.1692197 info_projection projection
```

# VEP Scan Benchmarks

``` r
if (has_vep) {
  vep_src <- build_read_bcf(vep_vcf)
  vep_tidy_src <- build_read_bcf(vep_vcf, tidy = TRUE)
  vep_cols <- get_bcf_columns(con, vep_vcf)
  vep_tidy_cols <- get_bcf_columns(con, vep_vcf, tidy = TRUE)
  vep_variants <- get_variant_count(con, vep_vcf)

  vep_cases <- list(
    list(
      name = "vep_count_all",
      sql = sprintf("SELECT COUNT(*) AS n FROM %s", vep_src),
      rows_hint = vep_variants
    )
  )

  vep_ann_cols <- grep("^VEP_", vep_cols, value = TRUE)
  if (length(vep_ann_cols) > 0) {
    pick_vep <- paste(head(vep_ann_cols, 8), collapse = ", ")
    vep_cases[[length(vep_cases) + 1]] <- list(
      name = "vep_annotation_projection",
      sql = sprintf("SELECT %s FROM %s LIMIT 200000", pick_vep, vep_src),
      rows_hint = 200000
    )
  }

  fmt_cols_wide <- grep("^FORMAT_", vep_cols, value = TRUE)
  fmt_cols_tidy <- grep("^FORMAT_", vep_tidy_cols, value = TRUE)
  if (length(fmt_cols_wide) > 0 && length(fmt_cols_tidy) > 0) {
    pick_wide <- paste(head(fmt_cols_wide, 8), collapse = ", ")
    pick_tidy <- paste(head(fmt_cols_tidy, 8), collapse = ", ")

    vep_cases[[length(vep_cases) + 1]] <- list(
      name = "vep_format_wide_projection",
      sql = sprintf("SELECT %s FROM %s LIMIT 100000", pick_wide, vep_src),
      rows_hint = 100000
    )
    vep_cases[[length(vep_cases) + 1]] <- list(
      name = "vep_format_tidy_projection",
      sql = sprintf("SELECT SAMPLE_ID, %s FROM %s LIMIT 100000", pick_tidy, vep_tidy_src),
      rows_hint = 100000
    )
  }

  vep_results <- do.call(
    rbind,
    lapply(vep_cases, function(x) run_bench(
      con,
      x$name,
      x$sql,
      iterations = 5,
      warmup = 1,
      bytes_read = vep_bytes,
      rows_hint = x$rows_hint
    ))
  )
  vep_results$engine <- "duckhts"
  vep_results <- add_relative_metrics(vep_results, "vep_count_all")
  vep_results <- add_case_keys(vep_results)
  vep_results <- finalize_results(vep_results)
  vep_results
} else {
  vep_results <- data.frame()
  message("Skipping VEP benchmarks. Set env var VEP_VCF to a local VEP VCF/BCF file.")
}
#> Skipping VEP benchmarks. Set env var VEP_VCF to a local VEP VCF/BCF file.
```

# DuckHTS Conversion Benchmarks

``` r
parquet_cases <- list(
  list(
    name = "clinvar_copy_core_parquet",
    sql = sprintf("SELECT CHROM, POS, REF, ALT FROM %s", clinvar_src),
    out_path = file.path(bench_dir, "clinvar_core.parquet"),
    bytes_read = clinvar_bytes,
    rows_hint = clinvar_variants
  )
)

if (length(info_cols) > 0) {
  pick <- paste(head(info_cols, 6), collapse = ", ")
  parquet_cases[[length(parquet_cases) + 1]] <- list(
    name = "clinvar_copy_info_parquet",
    sql = sprintf("SELECT CHROM, POS, %s FROM %s", pick, clinvar_src),
    out_path = file.path(bench_dir, "clinvar_info.parquet"),
    bytes_read = clinvar_bytes,
    rows_hint = clinvar_variants
  )
}

parquet_results <- do.call(
  rbind,
  lapply(parquet_cases, function(x) run_copy_bench(
    con,
    x$name,
    x$sql,
    x$out_path,
    iterations = 3,
    warmup = 1,
    bytes_read = x$bytes_read,
    rows_hint = x$rows_hint
  ))
)
parquet_results$engine <- "duckhts"
parquet_results <- add_case_keys(parquet_results)
parquet_results <- finalize_results(parquet_results)
parquet_results
#>                        case iterations min_sec median_sec mean_sec max_sec
#> 1 clinvar_copy_core_parquet          3   4.852      4.881 4.886333   4.926
#> 2 clinvar_copy_info_parquet          3   6.352      6.367 6.368000   6.385
#>      rows bytes_read compressed_mb compressed_mb_per_sec rows_hint rows_per_sec
#> 1 4352930  189031261            NA                    NA   4352930     891811.1
#> 2 4352930  189031261            NA                    NA   4352930     683670.5
#>   parquet_bytes parquet_mb parquet_mb_per_sec  engine          case_key
#> 1      11846170   11.29739           2.314564 duckhts copy_core_parquet
#> 2      29680804   28.30582           4.445708 duckhts copy_info_parquet
#>     workload
#> 1 conversion
#> 2 conversion
```

# DuckHTS BAM Benchmarks

``` r
bam_results <- data.frame()
bam_total_rows <- NA_real_
bam_region <- NA_character_
bam_region_rows <- NA_real_

if (has_bam) {
  bam_src <- build_read_bam(bam_path)
  bam_total_rows <- DBI::dbGetQuery(con, sprintf("SELECT COUNT(*) AS n FROM %s", bam_src))$n[[1]]

  if (has_samtools) {
    idxstats <- system2(samtools_bin, args = c("idxstats", bam_path), stdout = TRUE, stderr = FALSE)
    idx_fields <- strsplit(idxstats, "\t", fixed = TRUE)
    idx_rows <- Filter(function(x) length(x) >= 3 && x[[1]] != "*" && suppressWarnings(as.numeric(x[[3]])) > 0, idx_fields)
    if (length(idx_rows) > 0) {
      contig <- idx_rows[[1]][[1]]
      contig_len <- suppressWarnings(as.numeric(idx_rows[[1]][[2]]))
      region_end <- if (is.finite(contig_len)) min(contig_len, 5000000) else 5000000
      bam_region <- sprintf("%s:1-%d", contig, as.integer(region_end))
      bam_region_rows <- DBI::dbGetQuery(
        con,
        sprintf("SELECT COUNT(*) AS n FROM %s", build_read_bam(bam_path, region = bam_region))
      )$n[[1]]
    }
  }

  bam_cases <- list()
  if (is.finite(bam_total_rows) && bam_total_rows > 0) {
    bam_cases[[length(bam_cases) + 1]] <- list(
      name = "bam_count_all",
      sql = sprintf("SELECT COUNT(*) AS n FROM %s", bam_src),
      rows_hint = bam_total_rows
    )
    bam_cases[[length(bam_cases) + 1]] <- list(
      name = "bam_core_projection",
      sql = sprintf("SELECT QNAME, RNAME, POS, MAPQ, CIGAR FROM %s LIMIT 200000", bam_src),
      rows_hint = min(bam_total_rows, 200000)
    )
  } else {
    message("Skipping BAM full-scan benchmarks: read_bam() returned zero rows for the whole file.")
  }

  if (!is.na(bam_region)) {
    bam_cases[[length(bam_cases) + 1]] <- list(
      name = "bam_region_count",
      sql = sprintf("SELECT COUNT(*) AS n FROM %s", build_read_bam(bam_path, region = bam_region)),
      rows_hint = bam_region_rows
    )
    bam_cases[[length(bam_cases) + 1]] <- list(
      name = "bam_region_core_projection",
      sql = sprintf(
        "SELECT QNAME, RNAME, POS, MAPQ, CIGAR FROM %s LIMIT 200000",
        build_read_bam(bam_path, region = bam_region)
      ),
      rows_hint = min(bam_region_rows, 200000)
    )
  } else {
    message("Skipping BAM region benchmarks: no indexed contig with reads found via samtools idxstats.")
  }

  if (length(bam_cases) > 0) {
    bam_results <- do.call(
      rbind,
      lapply(bam_cases, function(x) run_bench(
        con,
        x$name,
        x$sql,
        iterations = 5,
        warmup = 1,
        bytes_read = bam_bytes,
        rows_hint = x$rows_hint
      ))
    )
    bam_results$engine <- "duckhts"
    if ("bam_count_all" %in% bam_results$case) {
      bam_results <- add_relative_metrics(bam_results, "bam_count_all")
    }
    bam_results <- add_case_keys(bam_results)
    bam_results <- finalize_results(bam_results)
    bam_results
  } else {
    bam_results
  }
} else {
  message("Skipping BAM benchmarks: BAM/BAM index not found.")
}
#> Skipping BAM full-scan benchmarks: read_bam() returned zero rows for the whole file.
#> 
#> ---
#> bam_region_count
#> SELECT COUNT(*) AS n FROM read_bam('HG00106.chrom11.ILLUMINA.bwa.GBR.exome.20130415.bam', region := '11:1-5000000')
#> 
#> ---
#> bam_region_core_projection
#> SELECT QNAME, RNAME, POS, MAPQ, CIGAR FROM read_bam('HG00106.chrom11.ILLUMINA.bwa.GBR.exome.20130415.bam', region := '11:1-5000000') LIMIT 200000
#>                         case iterations min_sec median_sec mean_sec max_sec
#> 1           bam_region_count          5   0.046      0.046   0.0472   0.051
#> 2 bam_region_core_projection          5   0.061      0.062   0.0626   0.066
#>     rows bytes_read compressed_mb compressed_mb_per_sec rows_hint rows_per_sec
#> 1      1  330528619            NA                    NA    240068      5218870
#> 2 200000  330528619            NA                    NA    200000      3225806
#>    engine                   case_key workload
#> 1 duckhts           bam_region_count   region
#> 2 duckhts bam_region_core_projection   region
```

# Napkin Summary

``` r
scan_summary <- rbind(clinvar_results, vep_results)

summary_cols <- c(
  "case",
  "median_sec",
  "compressed_mb_per_sec",
  "rows_per_sec",
  "vs_scan_baseline_mb",
  "vs_scan_baseline_rows"
)

scan_summary[, intersect(summary_cols, names(scan_summary))]
#>                      case median_sec compressed_mb_per_sec rows_per_sec
#> 1       clinvar_count_all      1.127              159.9594    3862404.6
#> 2 clinvar_core_projection      0.244                    NA     819672.1
#> 3 clinvar_info_projection      0.306                    NA     653594.8
#>   vs_scan_baseline_mb vs_scan_baseline_rows
#> 1            1.000000             1.0000000
#> 2            4.618852             0.2122181
#> 3            3.683007             0.1692197
```

# bcftools Baselines

``` r
bcftools_results <- data.frame()

if (has_bcftools) {
  bcftools_cases <- list(
    list(
      name = "bcftools_count_all",
      args = c("view", "-Ou", "-o", "/dev/null", clinvar_vcf),
      bytes_read = clinvar_bytes,
      rows_hint = clinvar_variants
    ),
    list(
      name = "bcftools_core_projection",
      args = c("query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n", clinvar_vcf),
      bytes_read = clinvar_bytes,
      rows_hint = clinvar_variants
    )
  )

  if (length(info_cols) > 0) {
    info_tags <- sub("^INFO_", "", head(info_cols, 6))
    info_format <- paste0("%INFO/", info_tags, collapse = "\t")
    bcftools_cases[[length(bcftools_cases) + 1]] <- list(
      name = "bcftools_info_projection",
      args = c(
        "query",
        "-f",
        paste0(info_format, "\n"),
        clinvar_vcf
      ),
      bytes_read = clinvar_bytes,
      rows_hint = clinvar_variants
    )
  }

  if (!is.na(clinvar_region)) {
    clinvar_region_variants <- DBI::dbGetQuery(
      con,
      sprintf("SELECT COUNT(*) AS n FROM %s", build_read_bcf(clinvar_vcf, region = clinvar_region))
    )$n[[1]]

    bcftools_cases[[length(bcftools_cases) + 1]] <- list(
      name = "bcftools_region_count",
      args = c("view", "-Ou", "-r", clinvar_region, "-o", "/dev/null", clinvar_vcf),
      bytes_read = NA_real_,
      rows_hint = clinvar_region_variants
    )
    bcftools_cases[[length(bcftools_cases) + 1]] <- list(
      name = "bcftools_region_core_projection",
      args = c("query", "-r", clinvar_region, "-f", "%CHROM\t%POS\t%REF\t%ALT\n", clinvar_vcf),
      bytes_read = clinvar_bytes,
      rows_hint = clinvar_region_variants
    )
  }

  bcftools_results <- do.call(
    rbind,
    lapply(bcftools_cases, function(x) run_bcftools_bench(
      x$name,
      x$args,
      iterations = 5,
      warmup = 1,
      bytes_read = x$bytes_read,
      rows_hint = x$rows_hint
    ))
  )
  bcftools_results <- add_case_keys(bcftools_results)
  bcftools_results <- finalize_results(bcftools_results)
  bcftools_results
} else {
  message("Skipping bcftools baselines: bcftools not found on PATH.")
}
#> 
#> ---
#> bcftools_count_all
#> '/usr/local/bin/bcftools' 'view' '-Ou' '-o' '/dev/null' 'clinvar.vcf.gz'
#> 
#> ---
#> bcftools_core_projection
#> '/usr/local/bin/bcftools' 'query' '-f' '%CHROM   %POS    %REF    %ALT
#> ' 'clinvar.vcf.gz'
#> 
#> ---
#> bcftools_info_projection
#> '/usr/local/bin/bcftools' 'query' '-f' '%INFO/AF_ESP %INFO/AF_EXAC   %INFO/AF_TGP    %INFO/ALLELEID  %INFO/CLNDN %INFO/CLNDNINCL
#> ' 'clinvar.vcf.gz'
#>     engine                     case iterations min_sec median_sec mean_sec
#> 1 bcftools       bcftools_count_all          5   4.816      4.871   4.8742
#> 2 bcftools bcftools_core_projection          5   2.812      2.828   2.8282
#> 3 bcftools bcftools_info_projection          5   5.366      5.429   5.4250
#>   max_sec    rows bytes_read compressed_mb compressed_mb_per_sec rows_hint
#> 1   4.950 4352930  189031261      180.2743               37.0097   4352930
#> 2   2.848 4352930  189031261            NA                    NA   4352930
#> 3   5.455 4352930  189031261            NA                    NA   4352930
#>   rows_per_sec        case_key   workload
#> 1     893642.0       count_all  full_scan
#> 2    1539225.6 core_projection projection
#> 3     801792.2 info_projection projection
```

# samtools Baselines

``` r
samtools_results <- data.frame()

if (has_bam && has_samtools) {
  samtools_cases <- list()
  if (is.finite(bam_total_rows) && bam_total_rows > 0) {
    samtools_cases[[length(samtools_cases) + 1]] <- list(
      name = "samtools_bam_count_all",
      args = c("view", "-u", "-o", "/dev/null", bam_path),
      bytes_read = bam_bytes,
      rows_hint = bam_total_rows
    )
    samtools_cases[[length(samtools_cases) + 1]] <- list(
      name = "samtools_bam_core_projection",
      args = c("view", bam_path),
      bytes_read = NA_real_,
      rows_hint = bam_total_rows
    )
  }

  if (!is.na(bam_region)) {
    samtools_cases[[length(samtools_cases) + 1]] <- list(
      name = "samtools_bam_region_count",
      args = c("view", "-u", "-o", "/dev/null", bam_path, bam_region),
      bytes_read = NA_real_,
      rows_hint = bam_region_rows
    )
    samtools_cases[[length(samtools_cases) + 1]] <- list(
      name = "samtools_bam_region_core_projection",
      args = c("view", bam_path, bam_region),
      bytes_read = bam_bytes,
      rows_hint = bam_region_rows
    )
  }

  if (length(samtools_cases) > 0) {
    samtools_results <- do.call(
      rbind,
      lapply(samtools_cases, function(x) run_samtools_bench(
        x$name,
        x$args,
        iterations = 5,
        warmup = 1,
        bytes_read = x$bytes_read,
        rows_hint = x$rows_hint
      ))
    )
    samtools_results$case_key <- sub("^samtools_", "", samtools_results$case)
    samtools_results <- finalize_results(samtools_results)
    samtools_results
  } else {
    samtools_results
  }
} else {
  message("Skipping samtools baselines: samtools or BAM inputs not available.")
}
#> 
#> ---
#> samtools_bam_region_count
#> '/usr/local/bin/samtools' view -u -o /dev/null HG00106.chrom11.ILLUMINA.bwa.GBR.exome.20130415.bam 11:1-5000000
#> 
#> ---
#> samtools_bam_region_core_projection
#> '/usr/local/bin/samtools' view HG00106.chrom11.ILLUMINA.bwa.GBR.exome.20130415.bam 11:1-5000000
#>     engine                                case iterations min_sec median_sec
#> 1 samtools           samtools_bam_region_count          5   0.220       0.23
#> 2 samtools samtools_bam_region_core_projection          5   0.309       0.31
#>   mean_sec max_sec   rows bytes_read compressed_mb compressed_mb_per_sec
#> 1    0.227   0.233 240068         NA            NA                    NA
#> 2    0.310   0.311 240068  330528619            NA                    NA
#>   rows_hint rows_per_sec                   case_key workload
#> 1    240068    1043773.9           bam_region_count   region
#> 2    240068     774412.9 bam_region_core_projection   region
```

# Matched Comparison Table

``` r
comparison_results <- bind_result_frames(
  clinvar_results,
  vep_results,
  parquet_results,
  bam_results,
  bcftools_results,
  samtools_results
)

comparison_cols <- c(
  "engine",
  "case",
  "case_key",
  "workload",
  "median_sec",
  "compressed_mb_per_sec",
  "rows_per_sec",
  "parquet_mb_per_sec"
)

comparison_results[, intersect(comparison_cols, names(comparison_results))]
#>      engine                                case                   case_key
#> 1   duckhts                   clinvar_count_all                  count_all
#> 2   duckhts             clinvar_core_projection            core_projection
#> 3   duckhts             clinvar_info_projection            info_projection
#> 4   duckhts           clinvar_copy_core_parquet          copy_core_parquet
#> 5   duckhts           clinvar_copy_info_parquet          copy_info_parquet
#> 6   duckhts                    bam_region_count           bam_region_count
#> 7   duckhts          bam_region_core_projection bam_region_core_projection
#> 8  bcftools                  bcftools_count_all                  count_all
#> 9  bcftools            bcftools_core_projection            core_projection
#> 10 bcftools            bcftools_info_projection            info_projection
#> 11 samtools           samtools_bam_region_count           bam_region_count
#> 12 samtools samtools_bam_region_core_projection bam_region_core_projection
#>      workload median_sec compressed_mb_per_sec rows_per_sec parquet_mb_per_sec
#> 1   full_scan      1.127              159.9594    3862404.6                 NA
#> 2  projection      0.244                    NA     819672.1                 NA
#> 3  projection      0.306                    NA     653594.8                 NA
#> 4  conversion      4.881                    NA     891811.1           2.314564
#> 5  conversion      6.367                    NA     683670.5           4.445708
#> 6      region      0.046                    NA    5218869.6                 NA
#> 7      region      0.062                    NA    3225806.5                 NA
#> 8   full_scan      4.871               37.0097     893642.0                 NA
#> 9  projection      2.828                    NA    1539225.6                 NA
#> 10 projection      5.429                    NA     801792.2                 NA
#> 11     region      0.230                    NA    1043773.9                 NA
#> 12     region      0.310                    NA     774412.9                 NA
```

# Optional External Comparison

Use this section to place side-by-side numbers from Exon, bcftools, or
other tools on the same file and workload shape.

Recommended comparison cases:

- Full scan count
- Narrow projection
- Indexed region query
- Conversion to Parquet

Suggested output columns:

- `engine`
- `case`
- `median_sec`
- `compressed_mb_per_sec`
- `rows_per_sec`
- `vs_duckhts`

Keep the file path, region, projection width, and thread count fixed
across tools.

# Optional: Export Results

``` r
out <- comparison_results
write.csv(out, "benchmark_results.csv", row.names = FALSE)
out
#>                                   case iterations min_sec median_sec mean_sec
#> 1                    clinvar_count_all          5   1.119      1.127 1.124600
#> 2              clinvar_core_projection          5   0.243      0.244 0.245200
#> 3              clinvar_info_projection          5   0.299      0.306 0.306000
#> 4            clinvar_copy_core_parquet          3   4.852      4.881 4.886333
#> 5            clinvar_copy_info_parquet          3   6.352      6.367 6.368000
#> 6                     bam_region_count          5   0.046      0.046 0.047200
#> 7           bam_region_core_projection          5   0.061      0.062 0.062600
#> 8                   bcftools_count_all          5   4.816      4.871 4.874200
#> 9             bcftools_core_projection          5   2.812      2.828 2.828200
#> 10            bcftools_info_projection          5   5.366      5.429 5.425000
#> 11           samtools_bam_region_count          5   0.220      0.230 0.227000
#> 12 samtools_bam_region_core_projection          5   0.309      0.310 0.310000
#>    max_sec    rows bytes_read compressed_mb compressed_mb_per_sec rows_hint
#> 1    1.128       1  189031261      180.2743              159.9594   4352930
#> 2    0.252  200000  189031261            NA                    NA    200000
#> 3    0.314  200000  189031261            NA                    NA    200000
#> 4    4.926 4352930  189031261            NA                    NA   4352930
#> 5    6.385 4352930  189031261            NA                    NA   4352930
#> 6    0.051       1  330528619            NA                    NA    240068
#> 7    0.066  200000  330528619            NA                    NA    200000
#> 8    4.950 4352930  189031261      180.2743               37.0097   4352930
#> 9    2.848 4352930  189031261            NA                    NA   4352930
#> 10   5.455 4352930  189031261            NA                    NA   4352930
#> 11   0.233  240068         NA            NA                    NA    240068
#> 12   0.311  240068  330528619            NA                    NA    240068
#>    rows_per_sec   engine vs_scan_baseline_mb vs_scan_baseline_rows
#> 1     3862404.6  duckhts            1.000000             1.0000000
#> 2      819672.1  duckhts            4.618852             0.2122181
#> 3      653594.8  duckhts            3.683007             0.1692197
#> 4      891811.1  duckhts                  NA                    NA
#> 5      683670.5  duckhts                  NA                    NA
#> 6     5218869.6  duckhts                  NA                    NA
#> 7     3225806.5  duckhts                  NA                    NA
#> 8      893642.0 bcftools                  NA                    NA
#> 9     1539225.6 bcftools                  NA                    NA
#> 10     801792.2 bcftools                  NA                    NA
#> 11    1043773.9 samtools                  NA                    NA
#> 12     774412.9 samtools                  NA                    NA
#>                      case_key   workload parquet_bytes parquet_mb
#> 1                   count_all  full_scan            NA         NA
#> 2             core_projection projection            NA         NA
#> 3             info_projection projection            NA         NA
#> 4           copy_core_parquet conversion      11846170   11.29739
#> 5           copy_info_parquet conversion      29680804   28.30582
#> 6            bam_region_count     region            NA         NA
#> 7  bam_region_core_projection     region            NA         NA
#> 8                   count_all  full_scan            NA         NA
#> 9             core_projection projection            NA         NA
#> 10            info_projection projection            NA         NA
#> 11           bam_region_count     region            NA         NA
#> 12 bam_region_core_projection     region            NA         NA
#>    parquet_mb_per_sec
#> 1                  NA
#> 2                  NA
#> 3                  NA
#> 4            2.314564
#> 5            4.445708
#> 6                  NA
#> 7                  NA
#> 8                  NA
#> 9                  NA
#> 10                 NA
#> 11                 NA
#> 12                 NA
```

# Notes For Reproducibility

- Record CPU model, core count, RAM, and storage type
  (NVMe/SATA/network).
- Keep `PRAGMA threads` fixed across runs.
- Run each case multiple times and report median.
- Prefer local files for stable throughput comparisons.
- Distinguish compressed input MB/s from Parquet output MB/s.
- Compare only identical workload shapes across engines.
- Treat `COUNT(*)` as the scan baseline, not as a parser microbenchmark.

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
    "SELECT DISTINCT id FROM read_hts_header('%s') WHERE record_type = 'INFO' AND lower(id) IN ('csq','ann','bcsq','vep')",
    gsub("'", "''", vep_vcf)
  )
  DBI::dbGetQuery(con, q)
}
```
