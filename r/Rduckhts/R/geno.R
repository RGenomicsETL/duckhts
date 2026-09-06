#' Read Record-Major Genotypes
#'
#' Read typed GT/PS calls without repeating variant text per sample. Each
#' record has a zero-based scan-local `record_index` and a list of calls with
#' original-header sample indices, nullable allele indices, per-slot phase bits,
#' and nullable scalar phase sets. An absent GT has NULL allele/phase lists;
#' a missing allele still occupies a slot. Phase bits follow HTSlib decoding,
#' including its leading-slot convention for VCF versions before 4.4.
#'
#' Full scans preserve the input stream when assigning ordinals. Indexed regions
#' use HTSlib's union order and start a new ordinal at zero. Use SQL `ORDER BY
#' record_index` when order matters; the ordinal is not a persistent file locator.
#' Empty sample selection and sparse calls preserve every selected variant row.
#'
#' @inheritParams rduckhts_bcf
#' @param table_name Optional table to create; `NULL` returns a data frame.
#' @param non_reference_only Omit calls without any called alternate allele.
#'   This does not infer phase or remove variant records.
#' @return A data frame when `table_name` is `NULL`, otherwise invisible `TRUE`.
#' @seealso [rduckhts_bcf_samples()]
#' @examples
#' con <- rduckhts_connect()
#' path <- system.file("extdata", "geno_calls.bcf", package = "Rduckhts")
#' rduckhts_geno(con, "calls", path, non_reference_only = TRUE)
#' DBI::dbGetQuery(con, "SELECT record_index, len(calls) AS n FROM calls ORDER BY record_index")
#' DBI::dbDisconnect(con, shutdown = TRUE)
#' @export
rduckhts_geno <- function(con, table_name = NULL, path, region = NULL,
                          index_path = NULL, samples = NULL,
                          non_reference_only = FALSE, scan_mode = "auto",
                          decompression_threads = 0, decode_error_policy = "null",
                          overwrite = FALSE) {
  params <- list()
  if (!is.null(region)) params$region <- sql_quote_string(con, region)
  if (!is.null(index_path)) params$index_path <- sql_quote_string(con, index_path)
  if (!is.null(samples)) params$samples <- sql_quote_string(con, samples)
  if (!is.logical(non_reference_only) || length(non_reference_only) != 1L || is.na(non_reference_only)) {
    stop("non_reference_only must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE", call. = FALSE)
  }
  params$non_reference_only <- if (non_reference_only) "true" else "false"
  params$scan_mode <- sql_quote_string(con, .validate_scan_mode_param(scan_mode))
  params$decompression_threads <- .validate_nonnegative_integer_param(
    decompression_threads, "decompression_threads")
  params$decode_error_policy <- sql_quote_string(con, decode_error_policy)
  query <- paste0("SELECT * FROM read_geno(", sql_quote_string(con, path), build_param_str(params), ")")
  if (is.null(table_name)) return(DBI::dbGetQuery(con, query))
  prefix <- if (overwrite) "CREATE OR REPLACE TABLE " else "CREATE TABLE "
  DBI::dbExecute(con, paste0(prefix, sql_quote_identifier(con, table_name), " AS ", query))
  invisible(TRUE)
}

#' Read the VCF/BCF Sample Catalog
#'
#' Return one row per selected sample with its original-header zero-based
#' `sample_index` and `sample_name`. Join this relation to `read_geno()` calls
#' from the same unchanged file; names are not duplicated into every call.
#'
#' @inheritParams rduckhts_bcf
#' @return A data frame with `sample_index` and `sample_name` columns.
#' @export
rduckhts_bcf_samples <- function(con, path, samples = NULL) {
  params <- if (is.null(samples)) list() else list(samples = sql_quote_string(con, samples))
  DBI::dbGetQuery(con, paste0("SELECT * FROM read_bcf_samples(", sql_quote_string(con, path),
                             build_param_str(params), ") ORDER BY sample_index"))
}
