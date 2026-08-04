# ---------------------------------------------------------------------------
# Load: connect to DuckDB and load the extension
# ---------------------------------------------------------------------------

#' Load the duckhts extension into a DuckDB connection
#'
#' @description
#' Compatibility wrapper around \code{rduckhts_connect()} and
#' \code{rduckhts_load()}. Prefer the \code{rduckhts_*} functions in new code.
#'
#' @param con An existing DuckDB connection, or \code{NULL} to create a
#'   package-owned connection with \code{rduckhts_connect()}.
#' @param extension_path Explicit path to the \code{.duckdb_extension} file.
#'   If \code{NULL}, uses the default location in the installed package.
#' @return The DuckDB connection (invisibly).
#' @export
duckhts_load <- function(con = NULL, extension_path = NULL) {
  if (is.null(con)) {
    con <- rduckhts_connect(extension_path = extension_path)
  } else {
    rduckhts_load(con, extension_path = extension_path)
  }
  invisible(con)
}
