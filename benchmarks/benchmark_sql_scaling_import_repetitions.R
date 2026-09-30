# Run from the repository root after stage_sql_scaling_audit.R.
# Source staging and keyed differentials are separate from the timed query.
source("r/duckhtsbench/R/registry.R")
source("r/duckhtsbench/R/stage.R")
Sys.setenv(DUCKHTSBENCH_REGISTRY = normalizePath(
  "r/duckhtsbench/inst/benchmark_registry.tsv", mustWork = TRUE))

extension <- normalizePath("build/release/duckhts.duckdb_extension", mustWork = TRUE)
output <- "benchmarks/sql_scaling_import_repetitions.tsv"
columns <- c("case", "multiplier", "threads", "input_records", "seconds",
             "peak_rss_kib", "peak_buffer_mib", "peak_spill_mib", "result_rows")

# Budget per process: 256 MiB DuckDB/R runtime + 3 times 64 bytes of
# decoded source and 128 bytes of materialized result per physical site.
# No query spill is allowed. Both caps are set before the measurements.
runtime_mib <- 256
live_bytes_per_site <- 64 + 128
spill_limit_mib <- 0
measurements <- vector("list", 18L)
index <- 0L
for (multiplier in c(1L, 2L, 4L)) {
  id <- sprintf("sql_scaling_phase3_%dx", multiplier)
  input <- duckhts_bench_artifact_path(id)
  duckhts_bench_validate_identity(id)
  stopifnot(file.exists(paste0(input, ".tbi")))
  for (threads in c(1L, 4L)) {
    for (repetition in seq_len(3L)) {
      command <- c("benchmarks/benchmark_sql_scaling_real_run.R", "--child",
                   "import_sites", threads, multiplier, extension, input)
      stderr <- tempfile("sql-import-stderr-")
      result <- suppressWarnings(system2("Rscript", shQuote(command),
                                         stdout = TRUE, stderr = stderr))
      if (length(result) != 1L || !is.null(attr(result, "status"))) {
        stop(paste("Import measurement failed:",
                   paste(readLines(stderr, warn = FALSE), collapse = "\n")))
      }
      unlink(stderr)
      row <- read.delim(text = result, header = FALSE, col.names = columns)
      stopifnot(row$result_rows == row$input_records)
      row$repetition <- repetition
      row$rss_limit_mib <- runtime_mib +
        3 * live_bytes_per_site * row$input_records / 1048576
      row$spill_limit_mib <- spill_limit_mib
      row$rss_verdict <- if (row$peak_rss_kib / 1024 <= row$rss_limit_mib) {
        "within_budget"
      } else {
        "over_budget"
      }
      row$spill_verdict <- if (row$peak_spill_mib <= spill_limit_mib) {
        "within_budget"
      } else {
        "over_budget"
      }
      index <- index + 1L
      measurements[[index]] <- row
      message(sprintf("%d x, %d threads, repetition %d: %.0f / %.0f MiB RSS",
                      multiplier, threads, repetition,
                      row$peak_rss_kib / 1024, row$rss_limit_mib))
    }
  }
}
write.table(do.call(rbind, measurements), output, sep = "\t", row.names = FALSE,
            quote = FALSE)
