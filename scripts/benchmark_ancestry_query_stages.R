#!/usr/bin/env Rscript
# Profile each physical phase of the wide ancestry query.
main <- function(args) {
  if (length(args) < 2L || length(args) > 4L) {
    stop("Usage: benchmark_ancestry_query_stages.R input.parquet reference.parquet [threads] [samples]")
  }
  threads <- if (length(args) >= 3L) as.integer(args[[3L]]) else 4L
  samples <- if (length(args) == 4L) as.integer(args[[4L]]) else 1L
  stopifnot(threads %in% c(1L, 4L), samples %in% c(1L, 11L, 22L, 44L))
  one_sample <- samples == 1L
  library(DBI)
  library(Rduckhts)
  library(duckhtsbench)
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  dbExecute(con, paste0("SET threads=", threads))
  dbExecute(con, "SET memory_limit='1GB'")
  dbExecute(con, "SET preserve_insertion_order=false")
  quote_str <- function(x) as.character(dbQuoteString(con, x))
  spill <- tempfile("ancestry-stage-spill-")
  aligned <- tempfile("ancestry-stage-aligned-", fileext = ".parquet")
  profile_path <- tempfile("ancestry-stage-profile-", fileext = ".json")
  on.exit(unlink(c(spill, aligned, profile_path), recursive = TRUE), add = TRUE)
  dbExecute(con, paste0("SET temp_directory=", quote_str(spill)))
  input <- paste0("read_parquet(", quote_str(normalizePath(args[[1L]])), ")")
  if (!one_sample) {
    input <- paste0("(SELECT * FROM ", input, " WHERE sample_id IN ",
                    "(SELECT DISTINCT sample_id FROM ", input,
                    " ORDER BY sample_id LIMIT ", samples, "))")
  }
  reference <- paste0("read_parquet(", quote_str(normalizePath(args[[2L]])), ")")
  classification <- Rduckhts:::.ancestry_classification_query(
    c(input, reference), "frequency")
  rows <- list()
  record <- function(name, query) {
    dbExecute(con, "PRAGMA enable_profiling='json'")
    dbExecute(con, paste0("SET profiling_output=", quote_str(profile_path)))
    start <- proc.time()[["elapsed"]]
    dbExecute(con, query)
    seconds <- proc.time()[["elapsed"]] - start
    dbExecute(con, "PRAGMA disable_profiling")
    metrics <- jsonlite::fromJSON(profile_path)
    profile_dir <- Sys.getenv("DUCKHTS_ANCESTRY_PROFILE_DIR", "")
    if (nzchar(profile_dir)) {
      file.copy(profile_path, file.path(profile_dir, paste0(name, ".json")),
                overwrite = TRUE)
    }
    rss <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
    rows[[name]] <<- data.frame(
      phase = name, seconds = seconds,
      peak_buffer_mib = metrics$system_peak_buffer_memory / 1048576,
      peak_temp_mib = metrics$system_peak_temp_dir_size / 1048576,
      process_high_water_mib = as.numeric(strsplit(trimws(rss), "[[:space:]]+")[[1L]][[2L]]) / 1024)
  }
  record("audit", paste0("CREATE TEMP TABLE audit AS ",
                         Rduckhts:::.ancestry_indexed_audit_query(classification)))
  stopifnot(dbGetQuery(con, "SELECT count(*) AS n FROM audit")$n == samples)
  record("alignment", paste0("COPY (",
    Rduckhts:::.ancestry_aligned_query(classification), ") TO ", quote_str(aligned),
    " (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 32768)"))
  aligned_source <- paste0("read_parquet(", quote_str(aligned), ")")
  fields <- dbGetQuery(con, paste0("DESCRIBE SELECT * FROM ", reference))$column_name
  groups <- sort(setdiff(fields, c("chromosome", "position", "allele_a", "allele_b",
                                    paste0("PC", seq_len(16L)))))
  correction_path <- duckhts_bench_stage_repository_fixtures(
    normalizePath("."), "ancestry-projection")[["ancestry_correction"]]
  coefficients <- utils::read.delim(correction_path)$coefficient
  indexed_source <- if (one_sample) aligned_source else paste0(
    "(SELECT a.*, i.sample_index FROM ", aligned_source,
    " a JOIN audit i USING (sample_id))")
  record("moments", paste0("CREATE TEMP TABLE moments AS ",
    Rduckhts:::.ancestry_moments_query(con, reference, groups, coefficients,
                                      indexed_source, one_sample)))
  record("solver", paste0("CREATE TEMP TABLE solved AS SELECT m.sample_id, ",
    Rduckhts:::.ancestry_solver_expression(groups, TRUE), " AS q FROM moments m"))
  proportions <- if (one_sample) dbGetQuery(con, "SELECT q FROM solved")$q[[1L]] else NULL
  record("prediction", paste0("CREATE TEMP TABLE prediction AS ",
    Rduckhts:::.ancestry_prediction_query(con, reference, groups, aligned_source,
                                          "solved", one_sample, proportions)))
  record("output", paste0("CREATE TEMP TABLE output AS ",
    Rduckhts:::.ancestry_wide_query(con, groups, "audit", "moments", "solved",
                                   "prediction", TRUE, 0.4)))
  matched <- dbGetQuery(con, paste0("SELECT count(*) AS n, sum(length(sample_id) + ",
    "length(allele_a) + length(allele_b) + 4 + 8 + 8) AS bytes FROM ",
    aligned_source))
  stopifnot(dbGetQuery(con, "SELECT count(*) AS n FROM output")$n == 21L * samples)
  build_bytes <- matched$bytes + if (one_sample) 0 else 4 * matched$n
  cat(sprintf(paste0("matched_rows=%s aligned_payload_mib=%.2f ",
                     "build_payload_mib=%.2f aligned_parquet_mib=%.2f\n"),
              matched$n, matched$bytes / 1048576, build_bytes / 1048576,
              file.info(aligned)$size / 1048576))
  utils::write.table(do.call(rbind, rows), stdout(), sep = "\t",
                     row.names = FALSE, quote = FALSE)
}
main(commandArgs(trailingOnly = TRUE))
