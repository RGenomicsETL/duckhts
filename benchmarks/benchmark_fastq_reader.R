#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
})

sql_string <- function(x) {
  sprintf("'%s'", gsub("'", "''", x, fixed = TRUE))
}

set_affinity <- function(cpu) {
  output <- system2(
    "taskset", c("-pc", cpu, Sys.getpid()), stdout = TRUE, stderr = TRUE
  )
  status <- attr(output, "status")
  stopifnot(is.null(status) || identical(status, 0L))
}

sha256 <- function(path) {
  output <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  stopifnot(is.null(attr(output, "status")))
  strsplit(output[[1]], "[[:space:]]+")[[1]][[1]]
}

peak_rss_kib <- function() {
  status <- readLines("/proc/self/status", warn = FALSE)
  line <- grep("^VmHWM:", status, value = TRUE)
  if (length(line) != 1L) return(NA_real_)
  as.double(strsplit(trimws(line), "[[:space:]]+")[[1]][[2]])
}

warm_file <- function(path) {
  input <- file(path, open = "rb")
  on.exit(close(input), add = TRUE)
  while (length(readBin(input, what = "raw", n = 16L * 1024L * 1024L)) > 0L) {
    # Populate the page cache outside the timed interval.
  }
  invisible(path)
}

read_exome_manifest <- function(manifest_path, input_dir) {
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  required <- c(
    "file_id", "sample", "library", "lane", "read_end",
    "compressed_bytes", "sha256", "url"
  )
  stopifnot(all(required %in% names(manifest)))
  manifest$path <- file.path(input_dir, basename(manifest$url))
  stopifnot(all(file.exists(manifest$path)))
  observed_size <- as.double(file.info(manifest$path)$size)
  stopifnot(identical(observed_size, as.double(manifest$compressed_bytes)))
  observed_sha <- vapply(manifest$path, sha256, character(1))
  stopifnot(identical(observed_sha, manifest$sha256))
  manifest
}

write_fixture <- function(path, reads = 2000000L, read_length = 150L) {
  if (file.exists(path)) {
    expected_lines <- as.double(reads) * 4
    wc_output <- system2("wc", c("-l", path), stdout = TRUE)
    fields <- strsplit(trimws(wc_output[[1]]), "[[:space:]]+")[[1]]
    if (as.double(fields[[1]]) == expected_lines) return(invisible(path))
    unlink(path)
  }

  sequence_a <- substr(
    paste(rep("ACGT", ceiling(read_length / 4)), collapse = ""),
    1L, read_length
  )
  sequence_b <- substr(
    paste(rep("TGCA", ceiling(read_length / 4)), collapse = ""),
    1L, read_length
  )
  quality_a <- paste(rep("I", read_length), collapse = "")
  quality_b <- paste(rep("G", read_length), collapse = "")
  output <- file(path, open = "wt")
  on.exit(close(output), add = TRUE)

  chunk_size <- 10000L
  for (first in seq.int(1L, reads, by = chunk_size)) {
    last <- min(reads, first + chunk_size - 1L)
    index <- seq.int(first, last)
    sequence <- ifelse(index %% 2L == 0L, sequence_a, sequence_b)
    quality <- ifelse(index %% 2L == 0L, quality_a, quality_b)
    records <- as.vector(rbind(
      sprintf("@bench:%09d/1", index), sequence, "+", quality
    ))
    writeLines(records, output, sep = "\n", useBytes = TRUE)
  }
  invisible(path)
}

open_duckhts <- function(extension) {
  con <- dbConnect(duckdb(config = list(allow_unsigned_extensions = "true")))
  dbExecute(con, sprintf("LOAD %s", sql_string(normalizePath(extension))))
  dbExecute(con, "PRAGMA threads=1")
  con
}

measure_queries <- function(build, revision, extension, fixture, queries,
                            repetitions, cpu, input_reads, input_bases) {
  set_affinity(cpu)
  con <- open_duckhts(extension)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  input_bytes <- as.double(file.info(fixture)$size)
  extension_digest <- sha256(extension)

  rows <- list()
  for (workload in names(queries)) {
    warm_file(fixture)
    expected <- as.double(dbGetQuery(con, queries[[workload]])$answer[[1]])
    for (rep in seq_len(repetitions)) {
      start <- proc.time()[["elapsed"]]
      answer <- as.double(dbGetQuery(con, queries[[workload]])$answer[[1]])
      elapsed <- proc.time()[["elapsed"]] - start
      stopifnot(answer == expected)
      rows[[length(rows) + 1L]] <- data.frame(
        engine = "DuckHTS",
        build = build,
        source_revision = revision,
        artifact_sha256 = extension_digest,
        workload = workload,
        repetition = rep,
        elapsed_seconds = elapsed,
        max_rss_kib = peak_rss_kib(),
        input_reads = input_reads,
        input_bases = input_bases,
        input_bytes = input_bytes,
        result_checksum = answer,
        configured_threads = 1L,
        cpu_affinity = cpu,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

run_timed <- function(program, args, time_path, stderr_path,
                      environment = character()) {
  format <- "elapsed_seconds=%e max_rss_kib=%M percent_cpu=%P"
  time_args <- c(
    "-f", shQuote(format), "-o", shQuote(time_path), "--",
    shQuote(normalizePath(program)), vapply(args, shQuote, character(1))
  )
  status <- system2(
    "/usr/bin/time", time_args, stdout = FALSE, stderr = stderr_path,
    env = environment
  )
  stopifnot(identical(status, 0L))
  fields <- strsplit(readLines(time_path, warn = FALSE)[[1]], " ")[[1]]
  values <- strsplit(fields, "=", fixed = TRUE)
  parsed <- setNames(vapply(values, `[[`, character(1), 2L),
                     vapply(values, `[[`, character(1), 1L))
  c(
    elapsed_seconds = as.double(parsed[["elapsed_seconds"]]),
    max_rss_kib = as.double(parsed[["max_rss_kib"]]),
    percent_cpu = as.double(sub("%$", "", parsed[["percent_cpu"]]))
  )
}

fastp_args <- function(input, json_path, html_path) {
  c(
    "-i", normalizePath(input), "-w", "1",
    "--disable_adapter_trimming", "--disable_quality_filtering",
    "--disable_length_filtering", "--disable_trim_poly_g",
    "--dont_eval_duplication", "-j", json_path, "-h", html_path
  )
}

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  stopifnot(length(file_arg) == 1L)
  normalizePath(sub("^--file=", "", file_arg[[1]]))
}

run_r_worker <- function(arguments) {
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(script_path()), arguments)
  )
  stopifnot(identical(status, 0L))
}

mode_synthetic_worker <- function(args) {
  stopifnot(length(args) == 10L)
  reads <- as.integer(args[[8]])
  read_length <- as.integer(args[[9]])
  fixture <- normalizePath(args[[5]])
  fixture_sql <- sql_string(fixture)
  queries <- c(
    count = sprintf(
      "SELECT count(*)::HUGEINT AS answer FROM read_fastq(%s, scan_mode := 'sequential')",
      fixture_sql
    ),
    name = sprintf(
      "SELECT sum(length(NAME))::HUGEINT AS answer FROM read_fastq(%s)",
      fixture_sql
    ),
    strings = sprintf(
      paste0(
        "SELECT sum(length(NAME) + length(SEQUENCE) + length(QUALITY))::HUGEINT AS answer ",
        "FROM read_fastq(%s)"
      ), fixture_sql
    ),
    packed = sprintf(
      paste0(
        "SELECT sum(length(NAME) + length(SEQUENCE) + length(QUALITY))::HUGEINT AS answer ",
        "FROM read_fastq(%s, sequence_encoding := 'nt16', ",
        "quality_representation := 'phred')"
      ), fixture_sql
    )
  )
  result <- measure_queries(
    args[[2]], args[[3]], args[[4]], fixture, queries,
    as.integer(args[[6]]), args[[7]], reads,
    as.double(reads) * read_length
  )
  write.csv(result, args[[10]], row.names = FALSE, quote = TRUE)
}

mode_synthetic <- function(args) {
  stopifnot(length(args) == 10L)
  baseline_extension <- normalizePath(args[[2]])
  current_extension <- normalizePath(args[[3]])
  fastp_binary <- normalizePath(args[[4]])
  fastp_library_path <- normalizePath(args[[5]])
  output_csv <- args[[6]]
  fixture_dir <- args[[7]]
  repetitions <- as.integer(args[[8]])
  cpu <- args[[9]]
  current_revision <- args[[10]]
  reads <- 2000000L
  read_length <- 150L

  dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)
  fixture <- file.path(fixture_dir, "fastq-reader-2m-150bp.fq")
  write_fixture(fixture, reads, read_length)

  previous_csv <- tempfile("fastq-previous-", fileext = ".csv")
  direct_csv <- tempfile("fastq-direct-", fileext = ".csv")
  run_r_worker(c(
    "synthetic-worker", "previous", "247e5c8", shQuote(baseline_extension),
    shQuote(fixture), repetitions, cpu, reads, read_length,
    shQuote(previous_csv)
  ))
  run_r_worker(c(
    "synthetic-worker", "direct", current_revision,
    shQuote(current_extension), shQuote(fixture), repetitions, cpu, reads,
    read_length, shQuote(direct_csv)
  ))

  set_affinity(cpu)
  report_dir <- file.path(fixture_dir, "fastp-synthetic")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  fastp_rows <- list()
  for (rep in seq_len(repetitions)) {
    warm_file(fixture)
    stem <- file.path(report_dir, sprintf("run-%02d", rep))
    json_path <- paste0(stem, ".json")
    timing <- run_timed(
      fastp_binary, fastp_args(fixture, json_path, paste0(stem, ".html")),
      paste0(stem, ".time"), paste0(stem, ".log"),
      sprintf("LD_LIBRARY_PATH=%s", fastp_library_path)
    )
    report <- jsonlite::fromJSON(json_path)
    stopifnot(report$summary$before_filtering$total_reads == reads)
    fastp_rows[[rep]] <- data.frame(
      engine = "fastp", build = "current",
      source_revision = "d517536b021bca0916cf33cb456f4e4b8aa24456",
      artifact_sha256 = sha256(fastp_binary), workload = "qc_no_output",
      repetition = rep, elapsed_seconds = timing[["elapsed_seconds"]],
      max_rss_kib = timing[["max_rss_kib"]], input_reads = reads,
      input_bases = as.double(reads) * read_length,
      input_bytes = as.double(file.info(fixture)$size),
      result_checksum = reads, configured_threads = 1L,
      cpu_affinity = cpu, stringsAsFactors = FALSE
    )
  }
  result <- rbind(
    read.csv(previous_csv, stringsAsFactors = FALSE),
    read.csv(direct_csv, stringsAsFactors = FALSE),
    do.call(rbind, fastp_rows)
  )
  expected <- c(
    count = reads,
    name = as.double(reads) * 15,
    strings = as.double(reads) * (15 + 2 * read_length),
    packed = as.double(reads) * (15 + 2 * read_length)
  )
  duckhts_rows <- result$engine == "DuckHTS"
  stopifnot(all(
    result$result_checksum[duckhts_rows] ==
      expected[result$workload[duckhts_rows]]
  ))
  write.csv(result, output_csv, row.names = FALSE, quote = TRUE)
}

mode_gzip_worker <- function(args) {
  stopifnot(length(args) == 9L)
  fixture <- normalizePath(args[[5]])
  fixture_sql <- sql_string(fixture)
  queries <- c(
    gzip_strings = sprintf(
      paste0(
        "SELECT sum(length(NAME) + length(SEQUENCE) + length(QUALITY))::HUGEINT AS answer ",
        "FROM read_fastq(%s)"
      ), fixture_sql
    ),
    gzip_packed = sprintf(
      paste0(
        "SELECT sum(length(NAME) + length(SEQUENCE) + length(QUALITY))::HUGEINT AS answer ",
        "FROM read_fastq(%s, sequence_encoding := 'nt16', ",
        "quality_representation := 'phred')"
      ), fixture_sql
    )
  )
  result <- measure_queries(
    args[[2]], args[[3]], args[[4]], fixture, queries,
    as.integer(args[[6]]), args[[7]], NA_real_, NA_real_
  )
  write.csv(result, args[[9]], row.names = FALSE, quote = TRUE)
}

mode_gzip <- function(args) {
  stopifnot(length(args) == 9L)
  manifest <- read_exome_manifest(args[[4]], args[[5]])
  fixture <- manifest$path[[1]]
  previous_csv <- tempfile("fastq-gzip-previous-", fileext = ".csv")
  direct_csv <- tempfile("fastq-gzip-direct-", fileext = ".csv")
  run_r_worker(c(
    "gzip-worker", "previous", "247e5c8", shQuote(normalizePath(args[[2]])),
    shQuote(fixture), args[[7]], args[[8]], "unused", shQuote(previous_csv)
  ))
  run_r_worker(c(
    "gzip-worker", "direct", args[[9]], shQuote(normalizePath(args[[3]])),
    shQuote(fixture), args[[7]], args[[8]], "unused", shQuote(direct_csv)
  ))
  result <- rbind(
    read.csv(previous_csv, stringsAsFactors = FALSE),
    read.csv(direct_csv, stringsAsFactors = FALSE)
  )
  result$file_id <- manifest$file_id[[1]]
  write.csv(result, args[[6]], row.names = FALSE, quote = TRUE)
}

mode_fastp_exome <- function(args) {
  stopifnot(length(args) == 11L)
  manifest <- read_exome_manifest(args[[2]], args[[3]])
  fastp_binary <- normalizePath(args[[4]])
  fastp_library_path <- normalizePath(args[[5]])
  report_dir <- args[[6]]
  run_csv <- args[[7]]
  metric_csv <- args[[8]]
  cycle_csv <- args[[9]]
  repetitions <- as.integer(args[[10]])
  cpu <- args[[11]]
  set_affinity(cpu)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

  run_rows <- list()
  metric_rows <- list()
  cycle_rows <- list()
  for (i in seq_len(nrow(manifest))) {
    item <- manifest[i, ]
    report <- NULL
    for (rep in seq_len(repetitions)) {
      warm_file(item$path)
      stem <- file.path(report_dir, sprintf("%s-run-%02d", item$file_id, rep))
      json_path <- paste0(stem, ".json")
      timing <- run_timed(
        fastp_binary,
        fastp_args(item$path, json_path, paste0(stem, ".html")),
        paste0(stem, ".time"), paste0(stem, ".log"),
        sprintf("LD_LIBRARY_PATH=%s", fastp_library_path)
      )
      report <- jsonlite::fromJSON(json_path)
      before <- report$read1_before_filtering
      run_rows[[length(run_rows) + 1L]] <- data.frame(
        engine = "fastp", build = "current",
        source_revision = "d517536b021bca0916cf33cb456f4e4b8aa24456",
        artifact_sha256 = sha256(fastp_binary), workload = "fastp_qc",
        file_id = item$file_id, repetition = rep,
        elapsed_seconds = timing[["elapsed_seconds"]],
        max_rss_kib = timing[["max_rss_kib"]],
        input_reads = before$total_reads, input_bases = before$total_bases,
        input_bytes = item$compressed_bytes,
        result_checksum = before$total_reads + before$total_bases,
        configured_threads = 1L, cpu_affinity = cpu,
        stringsAsFactors = FALSE
      )
    }

    summary <- report$summary$before_filtering
    before <- report$read1_before_filtering
    metric_rows[[i]] <- data.frame(
      file_id = item$file_id, read_end = item$read_end,
      total_reads = before$total_reads, total_bases = before$total_bases,
      q20_bases = before$q20_bases, q30_bases = before$q30_bases,
      q40_bases = before$q40_bases, total_cycles = before$total_cycles,
      mean_length = summary$read1_mean_length, gc_rate = summary$gc_content,
      stringsAsFactors = FALSE
    )
    cycles <- seq_len(before$total_cycles)
    cycle_rows[[i]] <- data.frame(
      file_id = item$file_id, read_end = item$read_end, cycle = cycles,
      quality_mean = before$quality_curves$mean,
      quality_a = before$quality_curves$A,
      quality_t = before$quality_curves$T,
      quality_c = before$quality_curves$C,
      quality_g = before$quality_curves$G,
      quality_n = before$quality_curves$N,
      content_a = before$content_curves$A,
      content_t = before$content_curves$T,
      content_c = before$content_curves$C,
      content_g = before$content_curves$G,
      content_n = before$content_curves$N,
      content_gc = before$content_curves$GC,
      stringsAsFactors = FALSE
    )
  }
  write.csv(do.call(rbind, run_rows), run_csv, row.names = FALSE, quote = TRUE)
  write.csv(do.call(rbind, metric_rows), metric_csv, row.names = FALSE, quote = TRUE)
  write.csv(do.call(rbind, cycle_rows), cycle_csv, row.names = FALSE, quote = TRUE)
}

mode_duckhts_exome <- function(args) {
  stopifnot(length(args) == 9L)
  manifest <- read_exome_manifest(args[[2]], args[[3]])
  extension <- normalizePath(args[[4]])
  run_csv <- args[[5]]
  cycle_csv <- args[[6]]
  cpu <- args[[7]]
  revision <- args[[8]]
  build <- args[[9]]
  set_affinity(cpu)
  invisible(lapply(manifest$path, warm_file))

  scans <- vapply(seq_len(nrow(manifest)), function(i) {
    sprintf(
      paste0(
        "SELECT %s::VARCHAR AS file_id, %d::UTINYINT AS read_end, ",
        "SEQUENCE, QUALITY FROM read_fastq(%s, ",
        "sequence_encoding := 'nt16', quality_representation := 'phred')"
      ),
      sql_string(manifest$file_id[[i]]), manifest$read_end[[i]],
      sql_string(normalizePath(manifest$path[[i]]))
    )
  }, character(1))
  query <- sprintf(paste0(
    "WITH reads AS (%s), bases AS (",
    "SELECT file_id, read_end, ",
    "generate_subscripts(SEQUENCE, 1)::USMALLINT AS cycle, ",
    "unnest(SEQUENCE)::UTINYINT AS base, ",
    "unnest(QUALITY)::UTINYINT AS quality FROM reads) ",
    "SELECT file_id, read_end, cycle, count(*)::UBIGINT AS total_bases, ",
    "sum(quality)::HUGEINT AS quality_sum, ",
    "sum((quality >= 20)::UTINYINT)::HUGEINT AS q20_bases, ",
    "sum((quality >= 30)::UTINYINT)::HUGEINT AS q30_bases, ",
    "sum((quality >= 40)::UTINYINT)::HUGEINT AS q40_bases, ",
    "sum((base = 1)::UTINYINT)::HUGEINT AS a_bases, ",
    "sum((base = 2)::UTINYINT)::HUGEINT AS c_bases, ",
    "sum((base = 4)::UTINYINT)::HUGEINT AS g_bases, ",
    "sum((base = 8)::UTINYINT)::HUGEINT AS t_bases, ",
    "sum((base = 15)::UTINYINT)::HUGEINT AS n_bases, ",
    "sum((base NOT IN (1, 2, 4, 8, 15))::UTINYINT)::HUGEINT AS other_bases, ",
    "sum(if(base = 1, quality, 0))::HUGEINT AS a_quality_sum, ",
    "sum(if(base = 2, quality, 0))::HUGEINT AS c_quality_sum, ",
    "sum(if(base = 4, quality, 0))::HUGEINT AS g_quality_sum, ",
    "sum(if(base = 8, quality, 0))::HUGEINT AS t_quality_sum, ",
    "sum(if(base = 15, quality, 0))::HUGEINT AS n_quality_sum ",
    "FROM bases GROUP BY file_id, read_end, cycle ",
    "ORDER BY file_id, cycle"
  ), paste(scans, collapse = " UNION ALL "))

  con <- open_duckhts(extension)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  start <- proc.time()[["elapsed"]]
  cycles <- dbGetQuery(con, query)
  elapsed <- proc.time()[["elapsed"]] - start
  stopifnot(nrow(cycles) > 0L, all(cycles$other_bases == 0))
  write.csv(cycles, cycle_csv, row.names = FALSE, quote = TRUE)

  total_reads <- sum(cycles$total_bases[cycles$cycle == 1])
  total_bases <- sum(cycles$total_bases)
  run <- data.frame(
    engine = "DuckHTS", build = build, source_revision = revision,
    artifact_sha256 = sha256(extension), workload = "pure_sql_fastp_cycles",
    file_id = "all-eight-files", repetition = 1L,
    elapsed_seconds = elapsed, max_rss_kib = peak_rss_kib(),
    input_reads = total_reads, input_bases = total_bases,
    input_bytes = sum(manifest$compressed_bytes),
    result_checksum = total_reads + total_bases + sum(cycles$q30_bases),
    configured_threads = 1L, cpu_affinity = cpu,
    stringsAsFactors = FALSE
  )
  write.csv(run, run_csv, row.names = FALSE, quote = TRUE)
}

usage <- paste(
  "modes:",
  "synthetic BASELINE CURRENT FASTP FASTP_LIB OUT_CSV FIXTURE_DIR REPS CPU REV;",
  "gzip BASELINE CURRENT MANIFEST INPUT_DIR OUT_CSV REPS CPU REV;",
  "fastp-exome MANIFEST INPUT_DIR FASTP FASTP_LIB REPORT_DIR RUN_CSV",
  "METRIC_CSV CYCLE_CSV REPS CPU;",
  "duckhts-exome MANIFEST INPUT_DIR EXT RUN_CSV CYCLE_CSV CPU REV BUILD"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) stop(usage)
switch(
  args[[1]],
  "synthetic-worker" = mode_synthetic_worker(args),
  "synthetic" = mode_synthetic(args),
  "gzip-worker" = mode_gzip_worker(args),
  "gzip" = mode_gzip(args),
  "fastp-exome" = mode_fastp_exome(args),
  "duckhts-exome" = mode_duckhts_exome(args),
  stop(usage)
)
