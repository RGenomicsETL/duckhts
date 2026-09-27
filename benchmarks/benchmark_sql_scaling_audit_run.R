# Run from the repository root. Each measured query gets a fresh R/DuckDB process.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0L && identical(args[[1L]], "--child")) {
  library(DBI)
  case <- args[[2L]]
  stopifnot(length(args) == if (case == "norm") 7L else 9L)
  threads <- as.integer(args[[3L]])
  size <- as.integer(args[[4L]])
  extension <- args[[5L]]
  paths <- if (case == "norm") {
    stats::setNames(args[6:7], c("vcf", "fasta"))
  } else {
    stats::setNames(args[6:9], c("panel", "frequency", "evidence", "pairs"))
  }
  con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true"),
                                  shared_home = FALSE))
  sql_string <- function(value) as.character(dbQuoteString(con, value))
  dbExecute(con, paste("LOAD", sql_string(extension)))
  dbExecute(con, sprintf("SET threads=%d", threads))
  if (case == "norm") {
    dbExecute(con, paste0("CREATE TEMP TABLE bench_input AS SELECT CHROM, POS, REF, ALT ",
                          "FROM read_bcf(", sql_string(paths[["vcf"]]),
                          ", scan_mode := 'sequential') LIMIT ", size))
  } else {
    for (name in names(paths)) {
      site_limit <- if (case %in% c("panel_hash", "frequency_hash") &&
                        name %in% c("panel", "frequency")) {
        sprintf(" WHERE site_index < %d", size * 17000L / 32L)
      } else {
        ""
      }
      dbExecute(con, sprintf("CREATE TEMP VIEW %s AS SELECT * FROM read_parquet(%s)%s",
                             name, sql_string(paths[[name]]), site_limit))
    }
  }
  queries <- c(
    panel_hash = "SELECT duckhts_somalier_panel_sha256('panel') AS digest",
    frequency_hash = paste0("SELECT duckhts_somalier_frequency_sha256(",
                            "'frequency', 'panel') AS digest"),
    sketches = paste0("CREATE TEMP TABLE result AS SELECT * FROM ",
                      "duckhts_somalier_prepare_sketches(",
                      "'evidence', 'panel', 7, 0.3, 0.01, max_sites := 17000)"),
    charr = paste0("CREATE TEMP TABLE result AS SELECT * FROM ",
                   "duckhts_somalier_charr(",
                   "'evidence', 'panel', 'frequency', max_sites := 17000)"),
    matched = paste0("CREATE TEMP TABLE result AS SELECT * FROM ",
                     "duckhts_somalier_matched_contamination(",
                     "'evidence', 'panel', 'frequency', 'pairs', max_sites := 17000)"))
  if (case == "norm") {
    queries <- c(norm = paste0(
      "CREATE TEMP TABLE result AS SELECT * FROM duckhts_bcftools_norm(",
      "'bench_input', ", sql_string(paths[["fasta"]]),
      ", chrom_col := 'CHROM', pos_col := 'POS', ",
      "ref_col := 'REF', alt_col := 'ALT')"))
  }
  stopifnot(case %in% names(queries))
  profile_file <- tempfile(fileext = ".json")
  dbExecute(con, "PRAGMA enable_profiling='json'")
  dbExecute(con, paste("PRAGMA profiling_output=", sql_string(profile_file)))
  started <- proc.time()[["elapsed"]]
  if (case %in% c("panel_hash", "frequency_hash")) {
    result <- dbGetQuery(con, queries[[case]])
    result_rows <- nrow(result)
  } else {
    dbExecute(con, queries[[case]])
    result_rows <- NA_integer_
  }
  seconds <- proc.time()[["elapsed"]] - started
  dbExecute(con, "PRAGMA disable_profiling")
  profile <- jsonlite::fromJSON(profile_file)
  if (is.na(result_rows)) {
    result_rows <- dbGetQuery(con, "SELECT count(*) AS n FROM result")$n[[1L]]
  }
  status <- readLines("/proc/self/status", warn = FALSE)
  rss <- as.numeric(gsub("[^0-9]", "", status[startsWith(status, "VmHWM:")]))
  write.table(data.frame(case, size, threads, seconds,
                         peak_rss_kib = rss,
                         peak_buffer_mib = profile$system_peak_buffer_memory / 1048576,
                         result_rows), stdout(), sep = "\t", row.names = FALSE,
              col.names = FALSE, quote = FALSE)
  dbDisconnect(con, shutdown = TRUE)
  unlink(profile_file)
} else {
  stopifnot(length(args) %in% c(1L, 2L), file.exists(args[[1L]]))
  source("r/duckhtsbench/R/registry.R")
  Sys.setenv(DUCKHTSBENCH_REGISTRY = normalizePath(
    "r/duckhtsbench/inst/benchmark_registry.tsv", mustWork = TRUE))
  extension <- normalizePath(args[[1L]], mustWork = TRUE)
  norm_run <- length(args) == 2L && identical(args[[2L]], "--norm-run")
  fixed_pairs <- length(args) == 2L && identical(args[[2L]], "--fixed-pairs")
  stopifnot(length(args) == 1L || norm_run || fixed_pairs)
  if (norm_run) {
    cases <- "norm"
    sizes <- c(10000L, 20000L, 40000L)
    source("r/duckhtsbench/R/stage.R")
    vcf <- duckhts_bench_artifact_path("geno_giab_phased_source")
    fasta <- duckhts_bench_artifact_path("liftover_grch38_fasta")
    stopifnot(file.exists(vcf), file.exists(fasta))
    duckhts_bench_validate_identity("geno_giab_phased_source", vcf)
    duckhts_bench_validate_identity("liftover_grch38_fasta", fasta)
    paths <- c(vcf = vcf, fasta = fasta)
    output <- "benchmarks/sql_scaling_audit_norm.tsv"
  } else {
    source("r/duckhtsbench/R/stage.R")
    source("r/duckhtsbench/R/somalier.R")
    cases <- if (fixed_pairs) "matched" else
      c("panel_hash", "frequency_hash", "sketches", "charr", "matched")
    sizes <- c(8L, 16L, 32L)
    output <- if (fixed_pairs) "benchmarks/sql_scaling_audit_fixed_pairs.tsv" else
      "benchmarks/sql_scaling_audit_measurements.tsv"
  }
  measurements <- list()
  for (size in sizes) {
    if (!norm_run) {
      paths <- duckhts_bench_stage_somalier(
        samples = size,
        selected_samples = if (fixed_pairs) 3L else size,
        selected_pairs = if (fixed_pairs) 2L else size - 1L)
    }
    for (threads in c(1L, 4L)) {
      for (case in cases) {
        command <- c("benchmarks/benchmark_sql_scaling_audit_run.R", "--child",
                     case, threads, size, extension, unname(paths))
        stderr_file <- tempfile("sql-scaling-stderr-")
        result <- suppressWarnings(system2("Rscript", shQuote(command),
                                           stdout = TRUE, stderr = stderr_file))
        if (length(result) != 1L || !is.null(attr(result, "status"))) {
          stop(paste("Failed:", case, size, threads,
                     paste(readLines(stderr_file, warn = FALSE), collapse = "\n")))
        }
        measurements[[length(measurements) + 1L]] <- result
        message(result)
        unlink(stderr_file)
      }
    }
  }
  header <- "case\tsize\tthreads\tseconds\tpeak_rss_kib\tpeak_buffer_mib\tresult_rows"
  writeLines(c(header, unlist(measurements)), output)
  message("Wrote ", output)
}
