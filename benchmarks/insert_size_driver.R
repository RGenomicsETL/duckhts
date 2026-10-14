# Direct-tool insert-size benchmark. Smoke executes all arms without clocks.
insert_size_read_tsv <- function(path) {
  utils::read.delim(path, colClasses = "character", check.names = FALSE)
}

insert_size_run_tool <- function(command, args, scratch, smoke) {
  timing <- file.path(scratch, "time.tsv")
  log <- file.path(scratch, "tool.log")
  status <- if (smoke) system2(command, shQuote(args), stdout = log, stderr = log) else
    system2("/usr/bin/time", c("-f", shQuote("%e\t%M"), "-o", shQuote(timing),
      shQuote(command), shQuote(args)), stdout = log, stderr = log)
  if (status != 0L) {
    stop(paste(c("benchmark tool failed", readLines(log)), collapse = "\n"), call. = FALSE)
  }
  if (smoke) c(NA_real_, NA_real_) else scan(timing, quiet = TRUE)
}

insert_size_normalize <- function(frame, columns, keys) {
  frame <- frame[, columns, drop = FALSE]
  floats <- intersect(columns, c("mean_insert_size", "standard_deviation",
                                 "median_insert_size", "median_absolute_deviation"))
  for (column in floats) frame[[column]] <- sprintf("%.2f", as.numeric(frame[[column]]))
  frame[] <- lapply(frame, as.character)
  frame <- frame[do.call(order, frame[keys]), , drop = FALSE]
  rownames(frame) <- NULL
  frame
}

insert_size_riker <- function(config, prefix, raw = FALSE) {
  args <- c("isize", "-i", config$path, "-r", config$reference, "-o", prefix,
            "--threads", config$threads, "--min-frac", if (raw) 1.1 else 0.05)
  status <- system2(config$riker, shQuote(args), stdout = FALSE, stderr = FALSE)
  if (status != 0L) stop("riker failed", call. = FALSE)
  read_output <- function(suffix) insert_size_read_tsv(paste0(prefix, suffix))
  list(metrics = if (raw) NULL else read_output(".isize-metrics.txt"),
       histogram = read_output(".isize-histogram.txt"))
}

insert_size_tool_run <- function(config, duckdb, smoke) {
  prefixes <- file.path(config$scratch, paste0("pass-", seq_len(config$scale)))
  if (config$arm == "riker") {
    timings <- vapply(prefixes, function(prefix) {
      insert_size_run_tool(config$riker,
        c("isize", "-i", config$path, "-r", config$reference, "-o", prefix,
          "--threads", config$threads, "--min-frac", 0.05), config$scratch, smoke)
    }, numeric(2))
    measured <- c(sum(timings[1L, ]), max(timings[2L, ]))
    passes <- lapply(prefixes, function(prefix) list(
      metrics = insert_size_read_tsv(paste0(prefix, ".isize-metrics.txt")),
      histogram = insert_size_read_tsv(paste0(prefix, ".isize-histogram.txt"))))
    return(list(passes = passes, measured = measured, work_seconds = NA_real_))
  }
  quote <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
  copy <- function(query, path) sprintf(
    "COPY (%s) TO %s (FORMAT CSV, HEADER true, DELIMITER '\t');", query, quote(path))
  sql <- c("SET autoinstall_known_extensions = false;",
    "SET autoload_known_extensions = false;",
    paste0("LOAD ", quote(config$extension), ";"),
    sprintf("SET threads = %d;", config$threads),
    "SET memory_limit = '2GiB';", "SET max_temp_directory_size = '0B';",
    paste0("SET temp_directory = ", quote(file.path(config$scratch, "spill")), ";"))
  if (!smoke) sql <- c(sql,
    "SET VARIABLE work_start_ns = (SELECT epoch_ns(current_timestamp));")
  for (prefix in prefixes) {
    query <- if (config$arm == "duckhts_bam_insert_size_counts") sprintf(
      "SELECT * FROM duckhts_bam_insert_size_counts(%s, reference := %s, threads := %d)",
      quote(config$path), quote(config$reference), config$threads) else sprintf(
      "SELECT * FROM duckhts_bam_insert_size_counts_sql(%s, reference := %s)",
      quote(config$path), quote(config$reference))
    sql <- c(sql, paste0("CREATE OR REPLACE TEMP TABLE counts AS ", query, ";"),
      sprintf(paste0("CREATE OR REPLACE TEMP TABLE metrics AS ",
        "SELECT * FROM duckhts_insert_size_metrics('counts', sample_name := %s);"),
        quote(config$sample)),
      copy("SELECT * FROM counts", paste0(prefix, ".counts.tsv")),
      copy("SELECT * FROM metrics", paste0(prefix, ".metrics.tsv")),
      copy("SELECT insert_size,
        coalesce(sum(read_pairs) FILTER (WHERE pair_orientation='FR'),0) AS fr_count,
        coalesce(sum(read_pairs) FILTER (WHERE pair_orientation='RF'),0) AS rf_count,
        coalesce(sum(read_pairs) FILTER (WHERE pair_orientation='TANDEM'),0) AS tandem_count
        FROM counts WHERE insert_size <=
          (SELECT coalesce(max(histogram_cutoff), 4294967295) FROM metrics)
        GROUP BY insert_size", paste0(prefix, ".histogram.tsv")))
  }
  work_path <- file.path(config$scratch, "work.tsv")
  if (!smoke) sql <- c(sql, copy(
    "SELECT (epoch_ns(current_timestamp) - getvariable('work_start_ns')) / 1e9 AS seconds",
    work_path))
  script <- file.path(config$scratch, "run.sql")
  writeLines(sql, script)
  measured <- insert_size_run_tool(duckdb,
    c("-unsigned", "-no-init", "-batch", "-bail", "-f", script), config$scratch, smoke)
  if (length(list.files(file.path(config$scratch, "spill"), recursive = TRUE)) > 0L) {
    stop("undeclared spill", call. = FALSE)
  }
  passes <- lapply(prefixes, function(prefix) list(
    counts = insert_size_read_tsv(paste0(prefix, ".counts.tsv")),
    metrics = insert_size_read_tsv(paste0(prefix, ".metrics.tsv")),
    histogram = insert_size_read_tsv(paste0(prefix, ".histogram.tsv"))))
  work_seconds <- if (smoke) NA_real_ else as.numeric(insert_size_read_tsv(work_path)$seconds)
  if (!smoke) stopifnot(length(work_seconds) == 1L, is.finite(work_seconds), work_seconds >= 0)
  list(passes = passes, measured = measured, work_seconds = work_seconds)
}

run_insert_size_benchmark <- function(repo, extension, riker, smoke = TRUE, duckdb = "duckdb") {
  repo <- normalizePath(repo, mustWork = TRUE)
  extension <- normalizePath(extension, mustWork = TRUE)
  riker <- normalizePath(riker, mustWork = TRUE)
  duckdb <- normalizePath(Sys.which(duckdb), mustWork = TRUE)
  duckdb_version <- system2(duckdb, "--version", stdout = TRUE)
  build_settings <- readLines(file.path(repo, "Makefile"))
  api_target <- sub(".*=", "", grep("^TARGET_DUCKDB_VERSION=", build_settings, value = TRUE))
  header_version <- sub(".*=", "", grep("^DUCKDB_HEADER_VERSION=", build_settings, value = TRUE))
  stopifnot(any(grepl("0.4.1", system2(riker, "--version", stdout = TRUE))))
  load_sql <- paste0("LOAD '", gsub("'", "''", extension, fixed = TRUE), "';")
  status <- system2(duckdb, c("-unsigned", "-no-init", "-bail", "-c", shQuote(load_sql)),
                    stdout = FALSE)
  if (status != 0L) stop("DuckDB CLI cannot load the extension", call. = FALSE)
  pkgload::load_all(file.path(repo, "r/duckhtsbench"), quiet = TRUE)
  staged <- duckhtsbench::duckhts_bench_stage_repository_fixtures(repo, "insert-size-smoke")
  reference <- if (smoke) staged[["insert_size_smoke_reference_fa"]] else
    duckhtsbench::duckhts_bench_artifact_path("liftover_grch38_fasta")
  paths <- if (smoke) c(BAM = staged[["insert_size_smoke_edges_bam"]],
                        CRAM = staged[["insert_size_smoke_edges_cram"]]) else c(
    BAM = duckhtsbench::duckhts_bench_artifact_path("insert_size_hg00403_chr20_bam"),
    CRAM = duckhtsbench::duckhts_bench_artifact_path("roh_mixture_hg00403_chr20_cram"))
  stopifnot(all(file.exists(c(paths, reference))))
  scratch <- tempfile("insert-size-benchmark-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE))
  metrics_columns <- c("sample", "pair_orientation", "read_pairs", "mean_insert_size",
    "standard_deviation", "median_insert_size", "median_absolute_deviation",
    "mode_insert_size", "min_insert_size", "max_insert_size")
  histogram_columns <- c("insert_size", "fr_count", "rf_count", "tandem_count")
  rows <- list()
  for (format in names(paths)) {
    config <- list(path = paths[[format]], reference = reference, extension = extension,
      riker = riker, threads = 1L, scratch = scratch)
    expected <- insert_size_riker(config, file.path(scratch, "expected"))
    raw <- insert_size_riker(config, file.path(scratch, "raw"), raw = TRUE)$histogram
    config$sample <- unique(expected$metrics$sample)
    stopifnot(length(config$sample) == 1L)
    pairs <- sum(as.numeric(as.matrix(raw[c("fr_count", "rf_count", "tandem_count")])))
    bins <- sum(as.numeric(as.matrix(raw[c("fr_count", "rf_count", "tandem_count")])) > 0)
    # Identity denominators come from the committed SAM fixture or registry.
    registry <- duckhtsbench::duckhts_bench_registry()
    identity_id <- if (smoke) paste0("insert_size_smoke_edges_", tolower(format)) else
      if (format == "BAM") "insert_size_hg00403_chr20_bam" else "roh_mixture_hg00403_chr20_cram"
    identity <- getFromNamespace("duckhts_bench_identity_fields", "duckhtsbench")(
      registry$supplier_identity[registry$id == identity_id])
    records <- as.numeric(identity[["alignments"]])
    stopifnot(length(records) == 1L, is.finite(records), records > 0)
    for (threads in c(1L, 2L, 4L, 8L)) {
      for (scale in if (smoke) 1L else c(1L, 2L, 4L)) {
        for (arm in c("riker", "duckhts_bam_insert_size_counts",
                      "duckhts_bam_insert_size_counts_sql")) {
          for (replicate in seq_len(if (smoke) 1L else 3L)) {
            config$threads <- threads
            config$scale <- scale
            config$arm <- arm
            config$scratch <- tempfile("run-", tmpdir = scratch)
            dir.create(config$scratch)
            run <- insert_size_tool_run(config, duckdb, smoke)
            for (result in run$passes) {
              stopifnot(identical(
                insert_size_normalize(result$metrics, metrics_columns, "pair_orientation"),
                insert_size_normalize(expected$metrics, metrics_columns, "pair_orientation")),
                identical(insert_size_normalize(result$histogram, histogram_columns, "insert_size"),
                  insert_size_normalize(expected$histogram, histogram_columns, "insert_size")))
              if (arm != "riker") {
                actual_raw <- data.frame(insert_size = raw$insert_size,
                  fr_count = "0", rf_count = "0", tandem_count = "0")
                labels <- c(FR = "fr_count", RF = "rf_count", TANDEM = "tandem_count")
                for (i in seq_len(nrow(result$counts))) {
                  index <- match(as.character(result$counts$insert_size[[i]]), actual_raw$insert_size)
                  stopifnot(!is.na(index))
                  actual_raw[[labels[[result$counts$pair_orientation[[i]]]]]][[index]] <-
                    as.character(result$counts$read_pairs[[i]])
                }
                stopifnot(identical(insert_size_normalize(actual_raw, histogram_columns, "insert_size"),
                  insert_size_normalize(raw, histogram_columns, "insert_size")))
              }
            }
            measured <- run$measured
            # Hash keys/values/flags: 16.25 bytes per allocated slot; load <= .77.
            # Budget includes worst-case worker overlap and 8 MiB decoder state
            # per handle, plus a merged histogram and result rows (24 bytes/bin).
            capacity <- 2^ceiling(log2(max(4, bins / 0.77)))
            live <- threads * (16.25 * capacity + 8 * 1024^2) +
                    16.25 * capacity + 24 * bins
            rows[[length(rows) + 1L]] <- data.frame(format, arm, threads, scale, replicate,
              records = records * scale, read_pairs = pairs * scale, bins,
              metric_rows = nrow(result$metrics), histogram_rows = nrow(result$histogram),
              seconds = measured[[1L]], work_seconds = run$work_seconds,
              peak_rss_mib = measured[[2L]] / 1024,
              rss_budget_mib = 256 + 3 * live / 1024^2, parity = TRUE)
            unlink(config$scratch, recursive = TRUE)
          }
        }
      }
    }
  }
  list(revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
       dirty = system2("git", c("diff", "--quiet", "HEAD")) != 0L,
       driver_sha256 = sub(" .*", "", system2("sha256sum",
         shQuote(file.path(repo, "benchmarks/insert_size_driver.R")), stdout = TRUE)),
       riker_version = system2(riker, "--version", stdout = TRUE),
       riker_sha256 = sub(" .*", "", system2("sha256sum", shQuote(riker), stdout = TRUE)),
       duckdb_cli = duckdb, duckdb_version = duckdb_version,
       api_target = api_target, header_version = header_version,
       runtime = paste("R controller:", R.version.string, "; timed DuckDB CLI:", duckdb_version),
       smoke = smoke, rows = do.call(rbind, rows))
}
