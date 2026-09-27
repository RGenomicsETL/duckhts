#!/usr/bin/env Rscript

# Selection-only whole-chromosome timing with a declared, already-staged input.
library(DBI)
library(Rduckhts)
library(duckhtsbench)

args <- commandArgs(TRUE)
if (length(args) != 2L || !args[[1L]] %in% c("chr22", "chrX") ||
    !args[[2L]] %in% c("0", "2")) {
  stop("usage: benchmark_somalier_find_sites_profile.R chr22|chrX 0|2")
}
chromosome <- args[[1L]]
decompression_threads <- as.integer(args[[2L]])
input <- duckhts_bench_artifact_path(paste0(
  "somalier_find_sites_", tolower(chromosome), "_source"
))
if (!file.exists(input)) stop("stage the declared VCF before timing: ", input)

con <- rduckhts_connect()
on.exit(dbDisconnect(con, shutdown = TRUE))
invisible(dbExecute(con, "SET threads = 1"))
invisible(dbExecute(con, sprintf(
  "CREATE TEMP VIEW population AS SELECT * FROM read_bcf(%s, samples := '', scan_mode := 'sequential', decompression_threads := %d)",
  as.character(dbQuoteString(con, input)), decompression_threads
)))
query <- Rduckhts:::.somalier_find_sites_query(
  con, "population", "GRCh38", 0.15, 6000, "AF", "AN", 10000, 0.48,
  list(include = NULL, exclude = NULL, gnotate = NULL),
  "somalier_v0.3.4", 65535, 10001, 5001
)
start <- proc.time()[["elapsed"]]
selected <- dbGetQuery(con, query)
cat("chromosome:", chromosome, "decompression_workers:",
    decompression_threads, "selected:", nrow(selected),
    "selection_seconds:", proc.time()[["elapsed"]] - start, "\n")
