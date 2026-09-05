#!/usr/bin/env Rscript
# POSIX fault-injection complement to installed-package tinytest contracts.
# Usage: Rscript test/scripts/reader_alloc_test.R /path/to/test-only/probe.so
library(DBI)
library(Rduckhts)

test_installed_reader_allocations <- function(probe_path) {
  shim <- dyn.load(normalizePath(probe_path, mustWork = TRUE))
  on.exit(dyn.unload(shim[["path"]]), add = TRUE, after = FALSE)
  con <- rduckhts_connect(config = list(threads = "1"))
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE, after = FALSE)
  extension <- system.file("duckhts_extension", "build", "duckhts.duckdb_extension", package = "Rduckhts")
  status <- .C("reader_alloc_r_open", extension, status = 0L, PACKAGE = shim[["name"]])$status
  stopifnot(status == 0L)
  on.exit(.C("reader_alloc_close", PACKAGE = shim[["name"]]), add = TRUE, after = FALSE)
  arm <- function(nth) {
    stopifnot(.C("reader_alloc_r_arm", as.integer(nth), status = 0L,
                PACKAGE = shim[["name"]])$status == 0L)
  }
  stats <- function() .C("reader_alloc_r_stats", count = 0L, remaining = 0L,
                         failed = 0L, PACKAGE = shim[["name"]])
  disarm <- function() invisible(.C("reader_alloc_disarm", PACKAGE = shim[["name"]]))
  fixtures <- c(read_bcf = "bcf_cache_lifecycle.vcf", read_bam = "bam_read_groups.sam",
                read_fasta = "region_names.fa", read_fastq = "r1.fq")
  failures <- 0L
  for (reader in names(fixtures)) {
    path <- system.file("extdata", fixtures[[reader]], package = "Rduckhts")
    stopifnot(nzchar(path))
    options <- if (reader == "read_bam") ", decompression_threads := 0" else ""
    sql <- sprintf("SELECT * FROM %s(%s, scan_mode := 'sequential'%s)",
                   reader, dbQuoteString(con, path), options)
    arm(0L)
    expected <- dbGetQuery(con, sql)
    dbGetQuery(con, "SELECT 4242")
    gc() # DBI prepared-result external pointers are finalized by R.
    control <- stats()
    stopifnot(control$count > 0L, control$remaining == 0L)
    disarm()
    for (nth in seq_len(control$count)) {
      cat(reader, "fail", nth, "of", control$count, "\n")
      arm(nth)
      error <- tryCatch({dbGetQuery(con, sql); NULL}, error = identity)
      stopifnot(inherits(error, "error"), grepl("out of memory", conditionMessage(error), fixed = TRUE))
      stopifnot(dbGetQuery(con, "SELECT 4242 AS n")$n == 4242L)
      gc()
      state <- stats()
      stopifnot(state$failed == 1L, state$remaining == 0L)
      disarm()
      stopifnot(identical(dbGetQuery(con, sql), expected))
      dbGetQuery(con, "SELECT 4242")
      failures <- failures + 1L
    }
  }
  cat("Installed R reader allocation failures:", failures,
      "errors, zero tracked leaks, exact DBI recovery: OK\n")
}

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
test_installed_reader_allocations(args[[1L]])
