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
                read_bcf_indexed = "bcf_scan_contigs.partial.vcf.gz",
                read_bam_materialize = "bam_materialize.sam",
                read_fasta = "region_names.fa", read_fastq = "r1.fq")
  failures <- 0L
  for (reader in names(fixtures)) {
    path <- system.file("extdata", fixtures[[reader]], package = "Rduckhts")
    stopifnot(nzchar(path))
    function_name <- switch(reader, read_bam_materialize = "read_bam",
                            read_bcf_indexed = "read_bcf", reader)
    options <- if (function_name == "read_bam") ", decompression_threads := 0" else ""
    if (reader == "read_bam_materialize")
      options <- paste0(options, ", standard_tags := true, auxiliary_tags := true")
    mode <- if (reader == "read_bcf_indexed") "auto" else "sequential"
    if (reader == "read_bcf_indexed")
      options <- paste0(options, ", index_path := ", dbQuoteString(con, paste0(path, ".index.tbi")))
    sql <- sprintf("SELECT * FROM %s(%s, scan_mode := '%s'%s) ORDER BY ALL",
                   function_name, dbQuoteString(con, path), mode, options)
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
  path <- dbQuoteString(con, system.file("extdata", "bam_materialize.sam", package = "Rduckhts"))
  projections <- c("CIGAR", "SEQ", "QUAL", "ML, FZ, CG", "AUXILIARY_TAGS", "*")
  options <- c("cigar_representation := 'binary'", "sequence_encoding := 'nt16'",
               "quality_representation := 'phred'", "standard_tags := true", "auxiliary_tags := true",
               paste("standard_tags := true, auxiliary_tags := true, cigar_representation := 'binary',",
                     "sequence_encoding := 'nt16', quality_representation := 'phred'"))
  list_arm <- function(kind, nth) stopifnot(.C("reader_list_r_arm", as.integer(kind),
    as.integer(nth), status = 0L, PACKAGE = shim[["name"]])$status == 0L)
  list_stats <- function() .C("reader_list_r_stats", count = 0L, failed = 0L,
                             unsafe_access = 0L, PACKAGE = shim[["name"]])
  list_disarm <- function() invisible(.C("reader_list_disarm", PACKAGE = shim[["name"]]))
  on.exit(list_disarm(), add = TRUE, after = FALSE)
  failures <- 0L
  for (i in seq_along(projections)) {
    sql <- sprintf("SELECT %s FROM read_bam(%s, %s, decompression_threads := 0)",
                   projections[[i]], path, options[[i]])
    for (kind in 1:2) {
      list_arm(kind, 0L)
      expected <- dbGetQuery(con, sql)
      count <- list_stats()$count
      stopifnot(count > 0L)
      list_disarm()
      for (nth in seq_len(count)) {
        list_arm(kind, nth)
        error <- tryCatch({dbGetQuery(con, sql); NULL}, error = identity)
        stopifnot(inherits(error, "error"), grepl("failed to grow output list", conditionMessage(error), fixed = TRUE))
        state <- list_stats()
        stopifnot(state$failed == 1L, state$unsafe_access == 0L)
        list_disarm()
        stopifnot(dbGetQuery(con, "SELECT 4242 AS n")$n == 4242L)
        stopifnot(identical(dbGetQuery(con, sql), expected))
        failures <- failures + 1L
      }
    }
  }
  cat("Installed R BAM list failures:", failures, "errors, no post-failure data access, exact DBI recovery: OK\n")
}

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
test_installed_reader_allocations(args[[1L]])
