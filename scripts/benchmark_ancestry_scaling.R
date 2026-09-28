#!/usr/bin/env Rscript
# Repeated, fresh-process measurements of the public ancestry wrapper.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: benchmark_ancestry_scaling.R epilepsy.parquet genotypes.parquet")
}
Sys.setenv(DUCKHTSBENCH_REGISTRY = normalizePath(
  "r/duckhtsbench/inst/benchmark_registry.tsv"))
reference <- duckhtsbench::duckhts_bench_stage_ancestry_parquet()
inputs <- c(epilepsy = args[[1L]], genotypes = args[[2L]], joint = args[[2L]])
stopifnot(all(file.exists(inputs)))
rows <- list()
for (workload in names(inputs)) {
  for (threads in c(1L, 4L)) {
    for (scale in c(1L, 2L, 4L)) {
      for (repetition in seq_len(3L)) {
        output <- tempfile("ancestry-scaling-")
        log <- tempfile("ancestry-scaling-log-")
        status <- system2("Rscript", c("--vanilla", "scripts/benchmark_ancestry_memory.R",
          workload, inputs[[workload]], reference, threads, scale),
          stdout = output, stderr = log)
        if (status != 0L) {
          stop("Scaling failure at ", workload, "/", threads, "/", scale,
               " repetition ", repetition, ": ", paste(readLines(log), collapse = "\n"))
        }
        lines <- readLines(output)
        header <- grep("^workload\\tthreads\\tscale\\t", lines)
        stopifnot(length(header) == 1L, length(lines) == header + 1L)
        row <- utils::read.delim(text = lines[header:length(lines)])
        stopifnot(nrow(row) == 1L, row$input_rows == row$used_variants,
                  row$peak_temp_mib == 0)
        row$repetition <- repetition
        decoded_full_mib <- if (workload == "epilepsy") 95.65 else 27.69
        row$decoded_required_column_mib <- if (workload == "epilepsy") {
          decoded_full_mib * scale / 4
        } else if (workload == "joint") {
          decoded_full_mib * scale^2 / 16
        } else {
          decoded_full_mib * scale / 4
        }
        row$budget_mib <- min(1024, 512 + 3 * decoded_full_mib)
        if (row$peak_rss_mib > row$budget_mib) {
          stop("Peak RSS exceeds declared ", workload, " budget at ",
               threads, " threads, ", scale, "x: ", row$peak_rss_mib,
               " > ", row$budget_mib, " MiB")
        }
        rows[[length(rows) + 1L]] <- row
        unlink(c(output, log))
        message(workload, " ", threads, "t ", scale, "x rep ", repetition,
                ": ", round(row$seconds, 3), "s / ", round(row$peak_rss_mib), " MiB")
      }
    }
  }
}
measurements <- do.call(rbind, rows)
utils::write.table(measurements, "benchmarks/ancestry_memory_repetitions.tsv",
                   sep = "\t", row.names = FALSE, quote = FALSE)
