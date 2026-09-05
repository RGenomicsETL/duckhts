# Source this driver from benchmark Rmds. Each backend gets a fresh R/DuckDB
# process: process-wide SIMD selection must not leak between measurements.
benchmark_simd_bam_gc <- function(extension, bam, max_reads = 0L, iterations = 5L,
                                  threads = 1L,
                                  backends = c("scalar", "auto", "avx2", "avx512"),
                                  modes = c("bam_scan", "materialized_seq")) {
  extension <- normalizePath(extension, mustWork = TRUE)
  bam <- normalizePath(bam, mustWork = TRUE)
  for (value in list(max_reads, iterations, threads)) {
    stopifnot(length(value) == 1L, is.finite(value), value == floor(value))
  }
  stopifnot(max_reads >= 0, iterations > 0, threads > 0,
            length(backends) > 0L, !anyNA(backends), all(nzchar(backends)),
            length(modes) > 0L, !anyDuplicated(modes),
            all(modes %in% c("bam_scan", "bam_scan_offset", "materialized_seq",
                             "bam_materialize_text", "bam_materialize_numeric")))

  results <- lapply(backends, function(backend) {
    callr::r(function(extension, bam, max_reads, iterations, threads, backend, modes) {
      con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")),
                           bigint = "integer64")
      on.exit(DBI::dbDisconnect(con, shutdown = TRUE))
      DBI::dbExecute(con, sprintf("SET threads = %.0f", threads))
      query <- function(sql) DBI::dbGetQuery(con, sql)
      quote <- function(text) as.character(DBI::dbQuoteString(con, text))
      DBI::dbExecute(con, paste("LOAD", quote(extension)))
      available <- backend == "auto" || isTRUE(query(paste0(
        "SELECT duckhts_simd_backend_available(", quote(backend), ") AS available"
      ))$available)
      common <- list(backend_request = backend, available = available,
                     skipped = !available, bam_path = bam, max_reads = max_reads,
                     threads = threads, iterations = iterations)
      if (!available) {
        return(lapply(modes, function(mode) c(common, list(
          benchmark = mode, skip_reason = "backend request is not available in this process"
        ))))
      }
      selected <- query(paste0("SELECT backend FROM duckhts_simd_set_backend(", quote(backend), ")"))$backend
      diagnostics <- query("SELECT duckhts_simd_backend() AS selected, duckhts_simd_requested_backend() AS requested")
      stopifnot(identical(selected, diagnostics$selected))
      kernel <- query("SELECT selected_backend, selected_capability, scalar_fallback
                       FROM duckhts_simd_kernel_info() WHERE kernel = 'seq_base_counts'")
      stopifnot(nrow(kernel) == 1L)
      common <- c(common, list(requested_backend = diagnostics$requested,
                               selected_backend = selected,
                               kernel_backend = kernel$selected_backend,
                               kernel_capability = kernel$selected_capability,
                               kernel_scalar_fallback = kernel$scalar_fallback))
      limit <- if (max_reads > 0) sprintf(" LIMIT %.0f", max_reads) else ""
      source <- paste0("SELECT SEQ FROM read_bam(", quote(bam), ") WHERE SEQ IS NOT NULL")
      aggregate_sql <- function(source) paste0(
        "SELECT count(*) AS reads, coalesce(sum(length(SEQ)), 0) AS total_bases, ",
        "coalesce(sum(coalesce(seq_gc_content(SEQ), 0.0)), 0.0) AS gc_fraction_sum FROM ", source
      )
      same_aggregate <- function(actual, expected) {
        stopifnot(identical(actual$reads, expected$reads),
                  identical(actual$total_bases, expected$total_bases),
                  abs(actual$gc_fraction_sum - expected$gc_fraction_sum) < 1e-9)
      }

      lapply(modes, function(mode) {
        load_sec <- NULL
        materialize <- startsWith(mode, "bam_materialize_")
        if (materialize) {
          options <- ", standard_tags := true, auxiliary_tags := true"
          if (mode == "bam_materialize_numeric") {
            options <- paste0(options, ", cigar_representation := 'binary', ",
                              "sequence_encoding := 'nt16', quality_representation := 'phred'")
          }
          sql <- paste0("CREATE OR REPLACE TEMP TABLE bam_materialized AS SELECT * FROM read_bam(",
                        quote(bam), options, ")", limit)
        } else {
          selected_source <- if (mode == "bam_scan_offset") {
            paste0(source, " AND FILE_OFFSET IS NOT NULL")
          } else source
          selected_source <- paste0("(", selected_source, limit, ")")
          if (mode == "materialized_seq") {
            start <- Sys.time()
            DBI::dbExecute(con, paste0("CREATE TEMP TABLE bam_gc_seqs AS ", source, limit))
            load_sec <- as.numeric(difftime(Sys.time(), start, units = "secs"))
            selected_source <- "bam_gc_seqs"
          }
          sql <- aggregate_sql(selected_source)
        }
        timings <- numeric(iterations)
        expected <- expected_digest <- NULL
        for (repetition in 0:iterations) {
          start <- Sys.time()
          if (materialize) DBI::dbExecute(con, sql) else observed <- query(sql)
          elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
          row_digest <- NULL
          if (materialize) {
            # Sort and hash every typed row, including duplicates, outside the
            # timer. Arrival order is not a reader contract, even with one worker.
            observed <- query(aggregate_sql("bam_materialized"))
            row_digest <- digest::digest(list(
              schema = query("DESCRIBE bam_materialized"),
              rows = query("SELECT * FROM bam_materialized ORDER BY ALL")
            ), algo = "sha256", serializeVersion = 3)
          }
          if (repetition == 0L) {
            expected <- observed
            expected_digest <- row_digest
          } else {
            same_aggregate(observed, expected)
            stopifnot(identical(row_digest, expected_digest))
            timings[repetition] <- elapsed
          }
        }
        median_sec <- median(timings)
        reads <- as.numeric(expected$reads)
        bases <- as.numeric(expected$total_bases)
        c(common, list(benchmark = mode, reads = reads, total_bases = bases,
                       gc_fraction_sum = expected$gc_fraction_sum, load_sec = load_sec,
                       timings_sec = timings, median_sec = median_sec, min_sec = min(timings),
                       reads_per_sec_median = reads / median_sec,
                       mbases_per_sec_median = bases / median_sec / 1e6,
                       row_multiset_sha256 = expected_digest))
      })
    }, args = list(extension = extension, bam = bam, max_reads = max_reads,
                   iterations = iterations, threads = threads, backend = backend, modes = modes))
  })
  results <- do.call(c, results)
  for (mode in modes) {
    measured <- Filter(function(x) x$benchmark == mode && !x$skipped, results)
    if (!length(measured)) next
    for (actual in measured) {
      stopifnot(identical(actual$reads, measured[[1L]]$reads),
                identical(actual$total_bases, measured[[1L]]$total_bases),
                abs(actual$gc_fraction_sum - measured[[1L]]$gc_fraction_sum) < 1e-9,
                identical(actual$row_multiset_sha256, measured[[1L]]$row_multiset_sha256))
    }
  }
  results
}
