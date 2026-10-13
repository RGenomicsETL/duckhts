#' BAM/CRAM Insert-Size Counts
#'
#' Reduces eligible records into sparse orientation and absolute-TLEN counts
#' using Fulcrum Genomics riker 0.4.1 semantics (MIT). Only non-first segments
#' count; QC-fail records are excluded. No SEQ, QUAL or auxiliary tags are read.
#' For weighted metrics, materialize the counts and call
#' `duckhts_insert_size_metrics('counts_table')` in SQL.
#'
#' @param con DuckDB connection with DuckHTS loaded.
#' @param path BAM or CRAM path.
#' @param reference Reference FASTA for CRAM, or NULL.
#' @param include_duplicates Include duplicate-flagged records.
#' @param threads Indexed histogram workers, from 1 to 8.
#' @param decompression_threads htslib decompression workers per file handle.
#' @param scan_mode One of auto, indexed, sequential. Sequential uses one worker.
#' @param index_path Optional BAM/CRAM index path.
#' @param reference_index_path Optional reference FASTA index path.
#' @return Data frame with pair_orientation, insert_size and read_pairs.
#'   Rows are unordered; absent bins are omitted.
#' @export
rduckhts_bam_insert_size_counts <- function(
  con, path, reference = NULL, include_duplicates = FALSE, threads = 1,
  decompression_threads = 0, scan_mode = c("auto", "indexed", "sequential"),
  index_path = NULL, reference_index_path = NULL
) {
  if (!is.logical(include_duplicates) || length(include_duplicates) != 1L ||
      is.na(include_duplicates)) {
    stop("include_duplicates must be TRUE or FALSE", call. = FALSE)
  }
  scan_mode <- match.arg(scan_mode)
  params <- list(
    include_duplicates = if (include_duplicates) "TRUE" else "FALSE",
    threads = sql_whole_number(threads, "threads", 8, low = 1),
    decompression_threads = sql_whole_number(
      decompression_threads, "decompression_threads", .Machine$integer.max),
    scan_mode = sql_quote_string(con, scan_mode)
  )
  paths <- list(reference = reference, index_path = index_path,
                reference_index_path = reference_index_path)
  for (name in names(paths)) {
    if (!is.null(paths[[name]])) params[[name]] <- sql_quote_string(con, paths[[name]])
  }
  DBI::dbGetQuery(con, sprintf(
    "SELECT * FROM duckhts_bam_insert_size_counts(%s%s)",
    sql_quote_string(con, path), build_param_str(params)))
}
