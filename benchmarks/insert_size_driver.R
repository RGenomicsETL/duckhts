# Fresh-process insert-size benchmark. Smoke executes all arms without clocks.
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
  read_output <- function(suffix) utils::read.delim(paste0(prefix, suffix),
    colClasses = "character", check.names = FALSE)
  list(metrics = if (raw) NULL else read_output(".isize-metrics.txt"),
       histogram = read_output(".isize-histogram.txt"))
}

insert_size_worker <- function(config) {
  result <- NULL
  if (config$arm == "riker") {
    for (pass in seq_len(config$scale)) {
      result <- insert_size_riker(config, file.path(config$scratch, paste0("riker-", pass)))
    }
    return(result)
  }
  con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true",
    autoinstall_known_extensions = "false", autoload_known_extensions = "false"), shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
  quote <- function(x) as.character(DBI::dbQuoteString(con, x))
  DBI::dbExecute(con, paste("LOAD", quote(config$extension)))
  DBI::dbExecute(con, sprintf("SET threads = %d", config$threads))
  DBI::dbExecute(con, "SET memory_limit = '2GiB'")
  DBI::dbExecute(con, "SET max_temp_directory_size = '0B'")
  DBI::dbExecute(con, paste("SET temp_directory =", quote(file.path(config$scratch, "spill"))))
  DBI::dbExecute(con, paste("SET VARIABLE insert_path =", quote(config$path)))
  DBI::dbExecute(con, paste("SET VARIABLE insert_reference =", quote(config$reference)))
  DBI::dbExecute(con, "SET VARIABLE insert_duplicates = false")
  for (pass in seq_len(config$scale)) {
    query <- if (config$arm == "native") sprintf(
      "SELECT * FROM duckhts_bam_insert_size_counts(%s, reference := %s, threads := %d)",
      quote(config$path), quote(config$reference), config$threads) else config$reference_sql
    DBI::dbExecute(con, paste("CREATE OR REPLACE TEMP TABLE counts AS", sub(";[[:space:]]*$", "", query)))
    metrics <- DBI::dbGetQuery(con, sprintf(
      "SELECT * FROM duckhts_insert_size_metrics('counts', sample_name := %s)", quote(config$sample)))
    cutoff <- if (nrow(metrics)) max(metrics$histogram_cutoff) else 4294967295
    histogram <- DBI::dbGetQuery(con, sprintf(
      "SELECT insert_size, coalesce(sum(read_pairs) FILTER (WHERE pair_orientation='FR'),0) AS fr_count,
       coalesce(sum(read_pairs) FILTER (WHERE pair_orientation='RF'),0) AS rf_count,
       coalesce(sum(read_pairs) FILTER (WHERE pair_orientation='TANDEM'),0) AS tandem_count
       FROM counts WHERE insert_size <= %.0f GROUP BY insert_size", cutoff))
    result <- list(counts = DBI::dbGetQuery(con, "SELECT * FROM counts"),
                   metrics = metrics, histogram = histogram)
  }
  if (length(list.files(file.path(config$scratch, "spill"), recursive = TRUE))) {
    stop("undeclared spill", call. = FALSE)
  }
  result
}

run_insert_size_benchmark <- function(repo, extension, riker, smoke = TRUE) {
  repo <- normalizePath(repo, mustWork = TRUE)
  extension <- normalizePath(extension, mustWork = TRUE)
  riker <- normalizePath(riker, mustWork = TRUE)
  stopifnot(any(grepl("0.4.1", system2(riker, "--version", stdout = TRUE))))
  pkgload::load_all(file.path(repo, "r/duckhtsbench"), quiet = TRUE)
  staged <- duckhtsbench::duckhts_bench_stage_repository_fixtures(repo, "insert-size-smoke")
  reference <- if (smoke) staged[["insert_size_smoke_reference_fa"]] else
    duckhtsbench::duckhts_bench_artifact_path("liftover_grch38_fasta")
  paths <- if (smoke) c(BAM = staged[["insert_size_smoke_edges_bam"]],
                        CRAM = staged[["insert_size_smoke_edges_cram"]]) else c(
    BAM = duckhtsbench::duckhts_bench_artifact_path("insert_size_hg00403_chr20_bam"),
    CRAM = duckhtsbench::duckhts_bench_artifact_path("roh_mixture_hg00403_chr20_cram"))
  stopifnot(all(file.exists(c(paths, reference))))
  reference_sql <- paste(readLines(staged[["insert_size_smoke_insert_size_reference_sql"]]), collapse = "\n")
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
      riker = riker, threads = 1L, scratch = scratch, reference_sql = reference_sql)
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
        for (arm in c("riker", "native", "projected_sql")) {
          for (replicate in seq_len(if (smoke) 1L else 3L)) {
            config$threads <- threads
            config$scale <- scale
            config$arm <- arm
            config_path <- file.path(scratch, "config.rds")
            result_path <- file.path(scratch, "result.rds")
            saveRDS(list(config = config, result_path = result_path), config_path)
            worker_args <- c(file.path(repo, "benchmarks/insert_size_driver.R"), "--worker", config_path)
            timing <- file.path(scratch, "time.tsv")
            status <- if (smoke) system2(file.path(R.home("bin"), "Rscript"), shQuote(worker_args)) else
              system2("/usr/bin/time", c("-f", shQuote("%e\t%M"), "-o", shQuote(timing),
                shQuote(file.path(R.home("bin"), "Rscript")), shQuote(worker_args)))
            if (status != 0L) stop("benchmark worker failed", call. = FALSE)
            result <- readRDS(result_path)
            stopifnot(identical(insert_size_normalize(result$metrics, metrics_columns, "pair_orientation"),
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
            measured <- if (smoke) c(NA_real_, NA_real_) else scan(timing, quiet = TRUE)
            # Hash keys/values/flags: 16.25 bytes per allocated slot; load <= .77.
            # Budget includes worst-case worker overlap and 8 MiB decoder state
            # per handle, plus a merged histogram and result rows (24 bytes/bin).
            capacity <- 2^ceiling(log2(max(4, bins / 0.77)))
            live <- threads * (16.25 * capacity + 8 * 1024^2) +
                    16.25 * capacity + 24 * bins
            rows[[length(rows) + 1L]] <- data.frame(format, arm, threads, scale, replicate,
              records = records * scale, read_pairs = pairs * scale, bins,
              metric_rows = nrow(result$metrics), histogram_rows = nrow(result$histogram),
              seconds = measured[[1L]], peak_rss_mib = measured[[2L]] / 1024,
              rss_budget_mib = 256 + 3 * live / 1024^2, parity = TRUE)
          }
        }
      }
    }
  }
  list(revision = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
       riker_version = system2(riker, "--version", stdout = TRUE),
       riker_sha256 = sub(" .*", "", system2("sha256sum", shQuote(riker), stdout = TRUE)),
       runtime = paste(R.version.string, "DuckDB", as.character(utils::packageVersion("duckdb"))),
       smoke = smoke, rows = do.call(rbind, rows))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L && args[[1L]] == "--worker") {
  request <- readRDS(args[[2L]])
  saveRDS(insert_size_worker(request$config), request$result_path)
}
