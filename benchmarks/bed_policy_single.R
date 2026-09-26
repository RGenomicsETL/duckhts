#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: Rscript benchmarks/bed_policy_single.R EXTENSION LABEL RUN")
}
source("r/duckhtsbench/R/registry.R")
Sys.setenv(DUCKHTSBENCH_REGISTRY = "r/duckhtsbench/inst/benchmark_registry.tsv")
file <- file.path(
  duckhts_bench_artifact_path("duckbedqc_118fc21"), "data", "GRCh38_exons.bed"
)
if (!file.exists(file) || !file.exists(args[[1L]])) {
  stop("stage corpus and build extension first")
}
library(DBI)
library(duckdb)
con <- dbConnect(
  duckdb(config = list(allow_unsigned_extensions = "true"), shared_home = FALSE),
  dbdir = ":memory:"
)
on.exit(dbDisconnect(con, shutdown = TRUE))
invisible(dbExecute(con, sprintf(
  "LOAD %s", as.character(dbQuoteString(con, normalizePath(args[[1L]])))
)))
invisible(dbExecute(con, "SET threads = 1"))
query <- sprintf(
  "SELECT sum(\"start\") + sum(\"end\") + count(*) AS result FROM read_bed(%s, scan_mode := 'sequential')",
  as.character(dbQuoteString(con, file))
)
start <- proc.time()[["elapsed"]]
result <- dbGetQuery(con, query)
seconds <- proc.time()[["elapsed"]] - start
if (nrow(result) != 1L || is.na(result$result[[1L]])) {
  stop("BED scan did not yield an aggregate")
}
write.table(
  data.frame(build = args[[2L]], run = as.integer(args[[3L]]),
             seconds = seconds, result = as.character(result$result[[1L]])),
  stdout(), sep = "\t", row.names = FALSE, quote = FALSE
)
