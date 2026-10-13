library(tinytest)
library(DBI)
library(Rduckhts)

local({
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  fixtures <- system.file("extdata", "insert_size", package = "Rduckhts")
  reference <- file.path(fixtures, "reference.fa")
  oracle_sql <- paste(readLines(file.path(fixtures, "insert_size_reference.sql")), collapse = "\n")
  cases <- read.delim(file.path(fixtures, "cases.tsv"))
  float_columns <- c("mean_insert_size", "standard_deviation", "median_insert_size",
                     "median_absolute_deviation")
  metric_columns <- c("sample", "pair_orientation", "read_pairs", float_columns,
                      "mode_insert_size", "min_insert_size", "max_insert_size")
  normalize <- function(frame, columns, key) {
    frame <- frame[, columns, drop = FALSE]
    frame[] <- lapply(frame, as.character)
    frame <- frame[do.call(order, frame[key]), , drop = FALSE]
    rownames(frame) <- NULL
    frame
  }
  read_oracle <- function(path, columns) {
    if (file.info(path)$size == 0) {
      return(as.data.frame(setNames(rep(list(character()), length(columns)), columns)))
    }
    read.delim(path, colClasses = "character", check.names = FALSE)
  }
  for (fixture in c("edges", "empty", "singleton", "task_overlap", "long_cigar", "even", "mad_zero", "minority")) {
    for (format in c("bam", "cram")) {
      path <- file.path(fixtures, paste(fixture, format, sep = "."))
      for (duplicates in c(FALSE, TRUE)) {
        counts <- rduckhts_bam_insert_size_counts(con, path, reference,
                                                 include_duplicates = duplicates)
        dbWriteTable(con, "counts", counts, overwrite = TRUE)
        dbExecute(con, paste("SET VARIABLE insert_path =", dbQuoteString(con, path)))
        dbExecute(con, paste("SET VARIABLE insert_reference =", dbQuoteString(con, reference)))
        dbExecute(con, paste("SET VARIABLE insert_duplicates =", if (duplicates) "TRUE" else "FALSE"))
        dbExecute(con, "SET threads = 1")
        baseline <- dbGetQuery(con, oracle_sql)
        columns <- c("pair_orientation", "insert_size", "read_pairs")
        key <- c("pair_orientation", "insert_size")
        expect_equal(normalize(counts, columns, key), normalize(baseline, columns, key),
                     info = paste(fixture, format, duplicates, "projected SQL"))
        for (mode in c("sequential", "indexed")) {
          for (threads in c(1, 2, 4, 8)) {
            parallel <- rduckhts_bam_insert_size_counts(con, path, reference,
              include_duplicates = duplicates, threads = threads, scan_mode = mode)
            expect_equal(normalize(parallel, columns, key), normalize(counts, columns, key),
                         info = paste(fixture, format, mode, threads, duplicates))
          }
        }
        raw_case <- if (duplicates) "raw_duplicates" else "raw"
        raw <- read_oracle(paste0(path, ".", raw_case, ".isize-histogram.txt"),
                          c("sample", "insert_size", "fr_count", "rf_count", "tandem_count"))
        raw_counts <- do.call(rbind, lapply(c(FR = "fr_count", RF = "rf_count", TANDEM = "tandem_count"),
          function(column) data.frame(insert_size = raw$insert_size, read_pairs = raw[[column]])))
        raw_counts$pair_orientation <- rep(c("FR", "RF", "TANDEM"), each = nrow(raw))
        raw_counts <- raw_counts[as.numeric(raw_counts$read_pairs) > 0, , drop = FALSE]
        expect_equal(normalize(counts, columns, key), normalize(raw_counts, columns, key),
                     info = paste(fixture, format, duplicates, "riker raw integer counts"))
        for (i in which(cases$duplicates == duplicates)) {
          case <- cases[i, ]
          metrics <- dbGetQuery(con, sprintf(
            "SELECT * FROM duckhts_insert_size_metrics('counts', deviations := %.17g, min_frac := %.17g, sample_name := %s)",
            case$deviations, case$min_frac, dbQuoteString(con, "edges")))
          expected <- read_oracle(paste0(path, ".", case$name, ".isize-metrics.txt"), metric_columns)
          metrics[float_columns] <- lapply(metrics[float_columns], function(x) sprintf("%.2f", x))
          expect_equal(normalize(metrics, metric_columns, "pair_orientation"),
                       normalize(expected, metric_columns, "pair_orientation"),
                       info = paste(fixture, format, case$name, "riker metrics, exact 2dp"))
          cutoff <- if (nrow(metrics)) max(as.numeric(metrics$histogram_cutoff)) else Inf
          histogram <- dbGetQuery(con, sprintf(
            "SELECT insert_size, sum(read_pairs) FILTER (WHERE pair_orientation='FR') AS fr_count,
             sum(read_pairs) FILTER (WHERE pair_orientation='RF') AS rf_count,
             sum(read_pairs) FILTER (WHERE pair_orientation='TANDEM') AS tandem_count
             FROM counts WHERE insert_size <= %s GROUP BY insert_size ORDER BY insert_size",
            if (is.finite(cutoff)) as.character(cutoff) else "4294967295"))
          histogram[is.na(histogram)] <- 0
          expected_hist <- read_oracle(paste0(path, ".", case$name, ".isize-histogram.txt"),
                                      c("sample", names(histogram)))
          expect_equal(normalize(histogram, names(histogram), "insert_size"),
                       normalize(expected_hist, names(histogram), "insert_size"),
                       info = paste(fixture, format, case$name, "riker histogram"))
        }
      }
    }
  }
  path <- file.path(fixtures, "edges.bam")
  expect_error(rduckhts_bam_insert_size_counts(con, path, threads = 0), "threads")
  expect_error(rduckhts_bam_insert_size_counts(con, path, threads = 9), "threads")
  expect_error(rduckhts_bam_insert_size_counts(con, path, include_duplicates = NA), "include_duplicates")
  expect_error(rduckhts_bam_insert_size_counts(con, path, scan_mode = "other"))
  expect_error(rduckhts_bam_insert_size_counts(con, "missing.bam"))
  unindexed <- file.path(fixtures, "unindexed.bam")
  expect_equal(rduckhts_bam_insert_size_counts(con, unindexed),
               rduckhts_bam_insert_size_counts(con, path))
  expect_error(rduckhts_bam_insert_size_counts(con, unindexed, scan_mode = "indexed"), "index")
  expect_error(dbGetQuery(con, "SELECT * FROM duckhts_insert_size_metrics('counts', deviations := -1)"), "nonnegative")
})
