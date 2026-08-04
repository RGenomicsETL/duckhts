#' List DuckHTS Extension Functions
#'
#' Returns the package-bundled function catalog generated from the top-level
#' \code{functions.yaml} manifest in the duckhts repository.
#'
#' @param category Optional function category filter.
#' @param kind Optional function kind filter such as \code{"scalar"},
#'   \code{"table"}, or \code{"table_macro"}.
#'
#' @return A data frame describing the extension functions, including the
#'   DuckDB function name, kind, category, signature, return type, optional
#'   R helper wrapper, short description, and example SQL.
#'
#' @examples
#' catalog <- rduckhts_functions()
#' subset(catalog, category == "Sequence UDFs", select = c("name", "description"))
#' subset(rduckhts_functions(kind = "table"), select = c("name", "r_wrapper"))
#'
#' @export
rduckhts_functions <- function(category = NULL, kind = NULL) {
  catalog_path <- system.file(
    "function_catalog",
    "functions.tsv",
    package = "Rduckhts",
    mustWork = TRUE
  )
  catalog <- utils::read.delim(
    catalog_path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!is.null(category)) {
    catalog <- catalog[catalog$category %in% category, , drop = FALSE]
  }
  if (!is.null(kind)) {
    catalog <- catalog[catalog$kind %in% kind, , drop = FALSE]
  }
  rownames(catalog) <- NULL
  catalog
}

.simd_backend_arg <- function(backend) {
  if (!is.character(backend) || length(backend) != 1L || is.na(backend)) {
    stop("backend must be a single non-missing character string", call. = FALSE)
  }
  backend
}
