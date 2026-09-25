#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: Rscript benchmarks/bed_policy_run.R EXTENSION LABEL REVISION OUTPUT")
}
source("r/duckhtsbench/R/registry.R")
Sys.setenv(DUCKHTSBENCH_REGISTRY = "r/duckhtsbench/inst/benchmark_registry.tsv")
file <- file.path(
  duckhts_bench_artifact_path("duckbedqc_118fc21"), "data", "GRCh38_exons.bed"
)
if (!file.exists(file) || !file.exists(args[[1L]])) {
  stop("stage corpus and build extension first")
}
input_lines <- length(readLines(file, warn = FALSE))
library(DBI)
library(duckdb)
policies <- if (args[[2L]] == "baseline") "error" else c("error", "skip", "report")
results <- list()
for (threads in c("1", "default")) {
  con <- dbConnect(
    duckdb(config = list(allow_unsigned_extensions = "true"), shared_home = FALSE),
    dbdir = ":memory:"
  )
  dbExecute(con, sprintf(
    "LOAD %s", as.character(dbQuoteString(con, normalizePath(args[[1L]])))
  ))
  if (threads == "1") dbExecute(con, "SET threads = 1")
  effective_threads <- dbGetQuery(
    con, "SELECT current_setting('threads') AS threads"
  )$threads[[1L]]
  for (policy in policies) {
    policy_arg <- if (policy == "error") "" else sprintf(
      ", error_policy := %s", as.character(dbQuoteString(con, policy))
    )
    query <- sprintf(
      "SELECT count(*) AS n, sum(start) AS start_sum, sum(\"end\") AS end_sum FROM read_bed(%s, scan_mode := 'sequential'%s)",
      as.character(dbQuoteString(con, file)), policy_arg
    )
    expected <- dbGetQuery(con, query)
    if (expected$n[[1L]] != input_lines) {
      stop("BED scan did not return every physical row")
    }
    for (run in seq_len(11L)) {
      start <- proc.time()[["elapsed"]]
      observed <- dbGetQuery(con, query)
      seconds <- proc.time()[["elapsed"]] - start
      if (!identical(observed, expected)) stop("inconsistent aggregation")
      results[[length(results) + 1L]] <- data.frame(
        revision = args[[3L]], build = args[[2L]], threads = threads,
        effective_threads = effective_threads, policy = policy, run = run,
        seconds = seconds, input_lines = input_lines,
        output_rows = observed$n[[1L]],
        start_sum = format(observed$start_sum[[1L]], scientific = FALSE, digits = 16),
        end_sum = format(observed$end_sum[[1L]], scientific = FALSE, digits = 16)
      )
    }
  }
  dbDisconnect(con, shutdown = TRUE)
}
write.table(do.call(rbind, results), args[[4L]], sep = "\t", row.names = FALSE, quote = FALSE)
