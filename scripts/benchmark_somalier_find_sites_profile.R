#!/usr/bin/env Rscript

# Selection-only whole-chromosome timing with a declared, already-staged input.
library(DBI)
library(Rduckhts)
library(duckhtsbench)

args <- commandArgs(TRUE)
if (length(args) != 3L || !args[[1L]] %in% c("chr22", "chrX") ||
    !args[[2L]] %in% c("0", "2") || !args[[3L]] %in% c("input", "empty")) {
  stop("usage: benchmark_somalier_find_sites_profile.R chr22|chrX 0|2 input|empty")
}
chromosome <- args[[1L]]
decompression_threads <- as.integer(args[[2L]])
empty <- identical(args[[3L]], "empty")
input <- duckhts_bench_artifact_path(paste0(
  "somalier_find_sites_", tolower(chromosome), "_source"
))
if (!file.exists(input)) stop("stage the declared VCF before timing: ", input)

con <- rduckhts_connect()
on.exit(dbDisconnect(con, shutdown = TRUE))
invisible(dbExecute(con, "SET threads = 1"))
invisible(dbExecute(con, sprintf(
  "CREATE TEMP VIEW population AS SELECT * FROM read_bcf(%s, samples := '', scan_mode := 'sequential', decompression_threads := %d)%s",
  as.character(dbQuoteString(con, input)), decompression_threads,
  if (empty) " WHERE false" else ""
)))
profiled <- Sys.getenv("DUCKHTS_FIND_SITES_PROFILE", "1") == "1"
if (profiled) {
  profile_path <- tempfile("somalier-profile-", fileext = ".json")
  on.exit(unlink(profile_path), add = TRUE)
  invisible(dbExecute(con, "SET enable_profiling = 'json'"))
  invisible(dbExecute(con, sprintf(
    "SET profiling_output = %s", as.character(dbQuoteString(con, profile_path))
  )))
}
source_path <- Sys.getenv("DUCKHTS_FIND_SITES_SOURCE")
query_function <- Rduckhts:::.somalier_find_sites_query
if (nzchar(source_path)) {
  source_environment <- new.env(parent = asNamespace("Rduckhts"))
  sys.source(source_path, envir = source_environment)
  query_function <- source_environment$.somalier_find_sites_query
}
query <- query_function(
  con, "population", "GRCh38", 0.15, 6000, "AF", "AN", 10000, 0.48,
  list(include = NULL, exclude = NULL, gnotate = NULL),
  "somalier_v0.3.4", 65535, 10001, 5001
)
start <- proc.time()[["elapsed"]]
selected <- dbGetQuery(con, query)
output_path <- Sys.getenv("DUCKHTS_FIND_SITES_OUTPUT")
if (nzchar(output_path)) {
  utils::write.table(
    selected[, c("region", "position", "source_ref", "source_alt")],
    output_path, sep = "\t", row.names = FALSE, quote = FALSE
  )
}
seconds <- proc.time()[["elapsed"]] - start
peak_buffer <- NA_real_
if (profiled) {
  peak_buffer <- jsonlite::fromJSON(profile_path)$system_peak_buffer_memory
  invisible(dbExecute(con, "SET enable_profiling = 'no_output'"))
}
usage <- dbGetQuery(con, "SELECT sum(memory_usage_bytes) AS bytes FROM duckdb_memory()")
cat("chromosome:", chromosome, "decompression_workers:",
    decompression_threads, "input:", if (empty) "empty" else "VCF",
    "selected:", nrow(selected), "selection_seconds:", seconds,
    "query_peak_buffer_bytes:", peak_buffer,
    "post_query_buffer_bytes:", usage$bytes, "\n")
