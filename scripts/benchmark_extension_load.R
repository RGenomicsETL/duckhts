#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && identical(args[[1L]], "--once")) {
  path <- normalizePath(args[[2L]], mustWork = TRUE)
  suppressWarnings(suppressPackageStartupMessages(library(DBI)))
  suppressWarnings(suppressPackageStartupMessages(library(duckdb)))
  driver <- duckdb::duckdb(dbdir = ":memory:", shared_home = FALSE,
                           config = list(allow_unsigned_extensions = "true",
                                         autoinstall_known_extensions = "false",
                                         threads = "1"))
  connection <- DBI::dbConnect(driver)
  on.exit(DBI::dbDisconnect(connection, shutdown = TRUE), add = TRUE)
  sql <- paste0("LOAD '", gsub("'", "''", path, fixed = TRUE), "'")
  elapsed <- unname(system.time(DBI::dbExecute(connection, sql))[["elapsed"]])
  writeLines(sprintf("%.6f", elapsed))
} else if (length(args) == 4L) {
  baseline <- normalizePath(args[[1L]], mustWork = TRUE)
  candidate <- normalizePath(args[[2L]], mustWork = TRUE)
  repetitions <- suppressWarnings(as.integer(args[[3L]]))
  if (is.na(repetitions) || repetitions < 1L || repetitions > 100L) {
    stop("repetitions must be an integer from 1 through 100", call. = FALSE)
  }
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  output <- data.frame(build = character(), repetition = integer(), seconds = numeric())
  for (repetition in seq_len(repetitions)) {
    order <- if (repetition %% 2L == 1L) c("baseline", "candidate") else c("candidate", "baseline")
    for (build in order) {
      path <- if (build == "baseline") baseline else candidate
      value <- suppressWarnings(system2("Rscript", shQuote(c(script, "--once", path)),
                                        stdout = TRUE, stderr = TRUE))
      if (!is.null(attr(value, "status")) || length(value) != 1L) {
        stop("extension LOAD failed: ", paste(value, collapse = "\n"), call. = FALSE)
      }
      seconds <- suppressWarnings(as.numeric(value))
      if (!is.finite(seconds)) stop("extension LOAD did not return elapsed seconds", call. = FALSE)
      output <- rbind(output, data.frame(build = build, repetition = repetition,
                                         seconds = seconds))
    }
  }
  utils::write.csv(output, args[[4L]], row.names = FALSE)
  print(stats::aggregate(seconds ~ build, output, stats::median), row.names = FALSE)
} else {
  stop("usage: benchmark_extension_load.R baseline_extension candidate_extension repetitions output.csv",
       call. = FALSE)
}
