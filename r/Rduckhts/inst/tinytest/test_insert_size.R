library(tinytest)
library(DBI)
library(Rduckhts)

local({
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  fixtures <- system.file("extdata", "insert_size", package = "Rduckhts")
  reference <- file.path(fixtures, "reference.fa.gz")
  oracle_sql <- paste(readLines(file.path(fixtures, "insert_size_reference.sql")), collapse = "\n")
  cases <- read.delim(file.path(fixtures, "cases.tsv"))
  oracle_metrics <- read.delim(file.path(fixtures, "metrics.tsv"), colClasses = "character")
  oracle_histograms <- read.delim(file.path(fixtures, "histograms.tsv"), colClasses = "character")
  float_columns <- c("mean_insert_size", "standard_deviation", "median_insert_size",
                     "median_absolute_deviation")
  metric_columns <- c("sample", "pair_orientation", "read_pairs", float_columns,
                      "mode_insert_size", "min_insert_size", "max_insert_size")
  columns <- c("pair_orientation", "insert_size", "read_pairs")
  key <- c("pair_orientation", "insert_size")
  normalize <- function(frame, columns, key) {
    frame <- frame[, columns, drop = FALSE]
    frame[] <- lapply(frame, as.character)
    frame <- frame[do.call(order, frame[key]), , drop = FALSE]
    rownames(frame) <- NULL
    frame
  }
  # SQL conformance owns the larger fixture matrix. Package tests exercise
  # both wrappers, bundled BAM/CRAM access and DBI summary/histogram workflows.
  for (fixture in c("edges", "empty")) {
    for (format in c("bam", "cram")) {
      path <- file.path(fixtures, paste(fixture, format, sep = "."))
      for (duplicates in c(FALSE, TRUE)) {
        dbExecute(con, paste("SET VARIABLE insert_path =", dbQuoteString(con, path)))
        dbExecute(con, paste("SET VARIABLE insert_reference =", dbQuoteString(con, reference)))
        dbExecute(con, paste("SET VARIABLE insert_duplicates =", if (duplicates) "TRUE" else "FALSE"))
        dbExecute(con, "SET threads = 1")
        baseline <- dbGetQuery(con, oracle_sql)
        raw_case <- if (duplicates) "raw_duplicates" else "raw"
        raw <- oracle_histograms[oracle_histograms$fixture == fixture &
          oracle_histograms$format == format & oracle_histograms$case_name == raw_case, ]
        raw_counts <- do.call(rbind, lapply(c(FR = "fr_count", RF = "rf_count", TANDEM = "tandem_count"),
          function(column) data.frame(insert_size = raw$insert_size, read_pairs = raw[[column]])))
        raw_counts$pair_orientation <- rep(c("FR", "RF", "TANDEM"), each = nrow(raw))
        raw_counts <- raw_counts[as.numeric(raw_counts$read_pairs) > 0, , drop = FALSE]
        for (arm in c("sql", "native")) {
          reader <- if (arm == "sql") rduckhts_bam_insert_size_counts_sql else
            rduckhts_bam_insert_size_counts
          counts <- reader(con, path, reference, include_duplicates = duplicates)
          expect_identical(normalize(counts, columns, key), normalize(baseline, columns, key),
                           info = paste(fixture, format, duplicates, arm, "projected SQL"))
          expect_identical(normalize(counts, columns, key), normalize(raw_counts, columns, key),
                           info = paste(fixture, format, duplicates, arm, "riker exact integers"))
          modes <- if (arm == "sql") c("auto", "sequential") else c("auto", "indexed", "sequential")
          for (mode in modes) {
            for (threads in c(1, 2, 4, 8)) {
              previous <- dbGetQuery(con, "SELECT current_setting('threads') AS n")$n
              parallel <- reader(con, path, reference, include_duplicates = duplicates,
                                 threads = threads, scan_mode = mode)
              expect_identical(normalize(parallel, columns, key), normalize(counts, columns, key),
                               info = paste(fixture, format, duplicates, arm, mode, threads))
              expect_identical(dbGetQuery(con, "SELECT current_setting('threads') AS n")$n, previous)
            }
          }
          dbWriteTable(con, "counts", counts, overwrite = TRUE)
          for (i in which(cases$duplicates == duplicates)) {
            case <- cases[i, ]
            metrics <- dbGetQuery(con, sprintf(
              "SELECT * FROM duckhts_insert_size_metrics('counts', deviations := %.17g, min_frac := %.17g, sample_name := %s)",
              case$deviations, case$min_frac, dbQuoteString(con, "edges")))
            expected <- oracle_metrics[oracle_metrics$fixture == fixture &
              oracle_metrics$format == format & oracle_metrics$case_name == case$name, ]
            metrics[float_columns] <- lapply(metrics[float_columns], function(x) sprintf("%.2f", x))
            expect_identical(normalize(metrics, metric_columns, "pair_orientation"),
                             normalize(expected, metric_columns, "pair_orientation"),
                             info = paste(fixture, format, arm, case$name, "riker exact 2dp"))
            cutoff <- if (nrow(metrics)) max(as.numeric(metrics$histogram_cutoff)) else Inf
            histogram <- dbGetQuery(con, sprintf(
              "SELECT insert_size, coalesce(sum(read_pairs) FILTER (WHERE pair_orientation='FR'),0) AS fr_count,
               coalesce(sum(read_pairs) FILTER (WHERE pair_orientation='RF'),0) AS rf_count,
               coalesce(sum(read_pairs) FILTER (WHERE pair_orientation='TANDEM'),0) AS tandem_count
               FROM counts WHERE insert_size <= %s GROUP BY insert_size ORDER BY insert_size",
              if (is.finite(cutoff)) as.character(cutoff) else "4294967295"))
            expected_hist <- oracle_histograms[oracle_histograms$fixture == fixture &
              oracle_histograms$format == format & oracle_histograms$case_name == case$name, ]
            expect_identical(normalize(histogram, names(histogram), "insert_size"),
                             normalize(expected_hist, names(histogram), "insert_size"),
                             info = paste(fixture, format, arm, case$name, "riker histogram"))
          }
        }
      }
    }
  }
  path <- file.path(fixtures, "edges.bam")
  for (reader in list(rduckhts_bam_insert_size_counts, rduckhts_bam_insert_size_counts_sql)) {
    expect_error(reader(con, path, threads = 0), "threads")
    expect_error(reader(con, path, threads = 9), "threads")
    expect_error(reader(con, path, include_duplicates = NA), "include_duplicates")
    expect_error(reader(con, path, scan_mode = "other"))
    expect_error(reader(con, "missing.bam"))
    expect_identical(normalize(reader(con, file.path(fixtures, "unindexed.bam")), columns, key),
                     normalize(reader(con, path), columns, key))
  }
  expect_error(rduckhts_bam_insert_size_counts(con, file.path(fixtures, "unindexed.bam"),
                                                    scan_mode = "indexed"), "index")
  expect_error(dbGetQuery(con, "SELECT * FROM duckhts_insert_size_metrics('counts', deviations := -1)"), "nonnegative")
  dbExecute(con, paste(
    "CREATE TEMP TABLE cutoff_counts AS SELECT 'FR' AS pair_orientation, x AS insert_size,",
    "1::BIGINT AS read_pairs FROM (VALUES (1::BIGINT), (3::BIGINT), (5::BIGINT)) t(x)"
  ))
  for (deviations in c(1e20, 1e308)) {
    metrics <- dbGetQuery(con, sprintf(
      "SELECT * FROM duckhts_insert_size_metrics('cutoff_counts', deviations := %.17g)", deviations))
    expect_equal(as.numeric(metrics$read_pairs), 3)
    expect_equal(metrics$mean_insert_size, 3)
    expect_equal(as.numeric(metrics$histogram_cutoff), 4294967295)
  }
})
