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
                read_geno = "geno_calls.bcf", read_bcf_samples = "geno_calls.vcf.gz",
                read_geno_format = "geno_format.bcf",
                read_bam_materialize = "bam_materialize.sam",
                read_fasta = "region_names.fa", read_fastq = "r1.fq")
  failures <- 0L
  for (reader in names(fixtures)) {
    path <- system.file("extdata", fixtures[[reader]], package = "Rduckhts")
    stopifnot(nzchar(path))
    function_name <- switch(reader, read_bam_materialize = "read_bam",
                            read_bcf_indexed = "read_bcf", read_geno_format = "read_geno", reader)
    options <- if (function_name == "read_bam") ", decompression_threads := 0" else ""
    if (reader == "read_bam_materialize")
      options <- paste0(options, ", standard_tags := true, auxiliary_tags := true")
    if (reader == "read_geno_format")
      options <- ", format_fields := ['AD', 'DP', 'GQ', 'VI', 'VF', 'ST']"
    mode <- if (reader == "read_bcf_indexed") "auto" else "sequential"
    if (reader == "read_bcf_indexed")
      options <- paste0(options, ", index_path := ", dbQuoteString(con, paste0(path, ".index.tbi")))
    scan_option <- if (function_name == "read_bcf_samples") "" else sprintf(", scan_mode := '%s'", mode)
    sql <- sprintf("SELECT * FROM %s(%s%s%s) ORDER BY ALL",
                   function_name, dbQuoteString(con, path), scan_option, options)
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
  queries <- vapply(seq_along(projections), function(i)
    sprintf("SELECT %s FROM read_bam(%s, %s, decompression_threads := 0)",
            projections[[i]], path, options[[i]]), character(1))
  for (fixture in c("geno_calls.bcf", "geno_format.bcf", "bcf_cache_lifecycle.vcf", "mapping_number_families.vcf")) {
    path <- dbQuoteString(con, system.file("extdata", fixture, package = "Rduckhts"))
    queries <- c(queries, sprintf("SELECT * FROM read_bcf(%s, scan_mode := 'sequential')", path),
                 sprintf("SELECT * FROM read_bcf(%s, tidy_format := true, scan_mode := 'sequential')", path))
    if (fixture == "geno_calls.bcf") for (projection in c("ALT", "calls", "*")) {
      queries <- c(queries, sprintf("SELECT %s FROM read_geno(%s)", projection, path))
    }
    if (fixture == "geno_format.bcf") {
      queries <- c(queries, sprintf(paste("SELECT calls FROM read_geno(%s,",
        "format_fields := ['AD', 'DP', 'GQ', 'VI', 'VF', 'ST'])"), path))
    }
  }
  error_patterns <- rep("failed to grow output list", length(queries))
  tx <- paste("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,",
    "100::UBIGINT transcript_start,111::UBIGINT transcript_end,1::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,",
    "111::UBIGINT cds_end,'ATGGCTGCTTAA'::BLOB cds_sequence,1::UTINYINT codon_table,",
    "''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence")
  ex <- paste("SELECT 0::UINTEGER transcript_index,100::UBIGINT exon_start,111::UBIGINT exon_end,",
    "1::UBIGINT exon_cdna_start,12::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase")
  stopifnot(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('list_probe',",
    "'SELECT 0::UINTEGER seq_region',", dbQuoteString(con, tx), ",", dbQuoteString(con, ex), ")"))$loaded)
  for (raw in c(FALSE, TRUE)) {
    calls <- paste("SELECT 1 AS event_index,0 AS seq_region,104 AS position,'C' AS reference,",
      "0 transcript_index,s sample_index,", if (raw) paste(
        "['A'] alternates,CASE WHEN s=0 THEN '1|0' WHEN s=1 THEN '0|1' ELSE '.|1' END gt") else paste(
        "'A' alternate,1 alt_index,CASE WHEN s=0 THEN [1,0] WHEN s=1 THEN [0,1] ELSE [NULL,1] END alleles,",
        "[true,true] phase_before,NULL::BIGINT phase_set"), "FROM range(3) samples(s)")
    queries <- c(queries, paste0("SELECT * FROM duckvep_haplotypes(", dbQuoteString(con, calls),
      ",'list_probe'", if (raw) ",input_mode:='source_records',phase_policy:='vep116_compat'", ") ORDER BY ALL"))
    error_patterns <- c(error_patterns,
      paste0("duckvep_haplotypes: (cannot reset (output|block event) list|output list allocation failed)",
        if (!raw) "|duckvep_phase_call: could not reserve output allele slots"))
  }
  failures <- 0L
  for (i in seq_along(queries)) {
    sql <- queries[[i]]
    for (kind in 1:2) {
      list_arm(kind, 0L)
      expected <- dbGetQuery(con, sql)
      count <- list_stats()$count
      stopifnot(count > 0L)
      list_disarm()
      for (nth in seq_len(count)) {
        list_arm(kind, nth)
        error <- tryCatch({dbGetQuery(con, sql); NULL}, error = identity)
        stopifnot(inherits(error, "error"))
        if (!grepl(error_patterns[[i]], conditionMessage(error)))
          stop("list query ", i, ", kind ", kind, ", failure ", nth, ": ", conditionMessage(error))
        state <- list_stats()
        stopifnot(state$failed == 1L, state$unsafe_access == 0L)
        list_disarm()
        stopifnot(dbGetQuery(con, "SELECT 4242 AS n")$n == 4242L)
        stopifnot(identical(dbGetQuery(con, sql), expected))
        failures <- failures + 1L
      }
    }
  }
  cat("Installed R BAM/BCF/genotype/haplotype list failures:", failures,
      "errors, no post-failure data access, exact DBI recovery: OK\n")
}

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
test_installed_reader_allocations(args[[1L]])
