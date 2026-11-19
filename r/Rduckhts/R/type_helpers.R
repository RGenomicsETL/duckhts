#' DuckDB to R Type Mappings
#'
#' Returns named character vectors that map DuckDB column types to the R
#' storage type that \code{duckdb} returns, and R type names to the DuckDB
#' types used for tabix column types (see \code{\link{normalize_tabix_types}}).
#'
#' \code{duckdb_to_r} gives the \code{typeof()} of the returned column, with
#' one exception: \code{"MAP"} maps to \code{"data.frame"}, the class of each
#' cell. The class is not part of the mapping:
#' \itemize{
#'   \item DATE, TIMESTAMP and TIME columns have storage type double and class
#'     \code{Date}, \code{POSIXct} and \code{difftime} (seconds).
#'   \item BIGINT and UBIGINT columns are double under the default
#'     \code{bigint = "numeric"} connection setting.
#'   \item A LIST column such as \code{VARCHAR[]} is a list with one vector per
#'     row. A MAP column is a list with one data frame per row, with columns
#'     \code{key} and \code{value}.
#' }
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{duckdb_to_r}{Named character vector mapping DuckDB types to R types}
#'   \item{r_to_duckdb}{Named character vector mapping R types to DuckDB types}
#' }
#'
#' @examples
#' mappings <- duckdb_type_mappings()
#' mappings$duckdb_to_r["BIGINT"]
#' mappings$r_to_duckdb["integer"]
#'
#' @export
duckdb_type_mappings <- function() {
  duckdb_to_r <- c(
    "BIGINT" = "double",
    "DOUBLE" = "double",
    "VARCHAR" = "character",
    "BOOLEAN" = "logical",
    "DATE" = "double", # class Date
    "TIMESTAMP" = "double", # class POSIXct
    "TIME" = "double", # class difftime, seconds
    "BLOB" = "list", # one raw vector per row
    "INTEGER" = "integer",
    "FLOAT" = "double",
    "SMALLINT" = "integer",
    "USMALLINT" = "integer",
    "UBIGINT" = "double",
    "VARCHAR[]" = "list",
    "BIGINT[]" = "list",
    "INTEGER[]" = "list",
    "DOUBLE[]" = "list",
    "FLOAT[]" = "list",
    "BOOLEAN[]" = "list",
    "MAP" = "data.frame" # class of each cell; the column is a list
  )

  r_to_duckdb <- c(
    "integer" = "INTEGER",
    "int" = "INTEGER",
    "int32" = "INTEGER",
    "int64" = "BIGINT",
    "numeric" = "DOUBLE",
    "double" = "DOUBLE",
    "float" = "DOUBLE",
    "character" = "VARCHAR",
    "string" = "VARCHAR",
    "chr" = "VARCHAR",
    "logical" = "BOOLEAN",
    "bool" = "BOOLEAN",
    "boolean" = "BOOLEAN",
    "Date" = "DATE",
    "POSIXct" = "TIMESTAMP",
    "POSIXt" = "TIMESTAMP",
    "hms" = "TIME",
    "raw" = "BLOB"
  )

  list(
    duckdb_to_r = duckdb_to_r,
    r_to_duckdb = r_to_duckdb
  )
}

#' Detect Complex Types in DuckDB Table
#'
#' Lists the columns of a table whose DuckDB type contains \code{[} or
#' \code{MAP}. \code{duckdb} returns such columns as R lists.
#'
#' @param con A DuckDB connection
#' @param table_name Name of the table to analyze
#'
#' @return A data frame with columns \code{column_name}, \code{column_type},
#'   \code{r_type} (\code{"vector"} or \code{"data.frame"}, the class of each
#'   cell) and \code{description}; a zero-column data frame when the table has
#'   no such columns.
#'
#' @examples
#' library(DBI)
#' library(duckdb)
#'
#' con <- rduckhts_connect()
#' bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
#' rduckhts_bcf(con, "variants", bcf_path, overwrite = TRUE)
#' complex_cols <- detect_complex_types(con, "variants")
#' print(complex_cols)
#' dbDisconnect(con, shutdown = TRUE)
#'
#' @export
detect_complex_types <- function(con, table_name) {
  # Get table schema
  schema <- DBI::dbGetQuery(
    con,
    sprintf("DESCRIBE %s", sql_quote_identifier(con, table_name))
  )

  # Find complex types (containing [ or MAP)
  complex_mask <- grepl("\\[|MAP", schema$column_type)
  complex_cols <- schema[complex_mask, ]

  if (nrow(complex_cols) == 0) {
    return(data.frame())
  }

  # Add R type description
  complex_cols$r_type <- ifelse(
    grepl("\\[", complex_cols$column_type),
    "vector",
    "data.frame"
  )
  complex_cols$description <- ifelse(
    grepl("\\[", complex_cols$column_type),
    "ARRAY type - will be R vector",
    "MAP type - will be R data frame"
  )

  complex_cols[, c("column_name", "column_type", "r_type", "description")]
}

#' Extract Array Elements Safely
#'
#' Extracts one element from each row of a DuckDB LIST or ARRAY column,
#' which \code{duckdb} returns as an R list.
#'
#' @param array_col A list column from DuckDB array data
#' @param index Numeric index (1-based). If NULL, returns full list
#' @param default Default value if index is out of bounds
#'
#' @return The array element at the specified index, or full array if index is NULL
#'
#' @examples
#' library(DBI)
#' library(duckdb)
#'
#' con <- rduckhts_connect()
#' bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
#' rduckhts_bcf(con, "variants", bcf_path, overwrite = TRUE)
#' data <- dbGetQuery(con, "SELECT ALT FROM variants LIMIT 5")
#' first_alt <- extract_array_element(data$ALT, 1)
#' all_alts <- extract_array_element(data$ALT)
#' dbDisconnect(con, shutdown = TRUE)
#'
#' @export
extract_array_element <- function(array_col, index = NULL, default = NA) {
  if (is.null(index)) {
    return(array_col)
  }

  # Safe extraction with bounds checking
  sapply(
    array_col,
    function(x) {
      if (is.null(x) || length(x) < index) {
        return(default)
      }
      return(x[[index]])
    },
    USE.NAMES = FALSE
  )
}

#' Extract MAP Keys and Values
#'
#' Extracts keys, values, or the value of one key from each row of a DuckDB
#' MAP column. \code{duckdb} returns a MAP column as a list of data frames with
#' columns \code{key} and \code{value}.
#'
#' @param map_col A MAP column: a list of data frames with columns \code{key}
#'   and \code{value}
#' @param operation What to extract: "keys", "values", or a specific key name
#' @param default Default value if key is not found (only used when operation is a key name)
#'
#' @return Extracted data based on the operation
#'
#' @examples
#' library(DBI)
#' library(duckdb)
#'
#' con <- rduckhts_connect()
#' gff_path <- system.file("extdata", "gff_file.gff.gz", package = "Rduckhts")
#' rduckhts_gff(con, "annotations", gff_path, attributes_map = TRUE, overwrite = TRUE)
#' data <- dbGetQuery(con, "SELECT attributes FROM annotations LIMIT 5")
#' keys <- extract_map_data(data$attributes, "keys")
#' name_values <- extract_map_data(data$attributes, "Name")
#' dbDisconnect(con, shutdown = TRUE)
#'
#' @export
extract_map_data <- function(map_col, operation = "keys", default = NA) {
  is_empty_map <- function(x) {
    if (is.null(x)) {
      return(TRUE)
    }
    if (length(x) == 0) {
      return(TRUE)
    }
    if (length(x) == 1 && is.na(x)) {
      return(TRUE)
    }
    if (is.data.frame(x) && nrow(x) == 0) {
      return(TRUE)
    }
    FALSE
  }

  extract_kv <- function(x, field) {
    if (is_empty_map(x)) {
      return(character(0))
    }
    if (is.data.frame(x) && field %in% names(x)) {
      return(x[[field]])
    }
    if (is.list(x) && !is.data.frame(x) && field %in% names(x)) {
      return(x[[field]])
    }
    if (is.character(x) || is.numeric(x)) {
      return(character(0))
    }
    character(0)
  }

  if (operation == "keys") {
    return(lapply(map_col, function(x) extract_kv(x, "key")))
  }

  if (operation == "values") {
    return(lapply(map_col, function(x) extract_kv(x, "value")))
  }

  # Search for specific key
  return(lapply(map_col, function(x) {
    if (is_empty_map(x)) {
      return(default)
    }
    if (is.data.frame(x) && "key" %in% names(x) && "value" %in% names(x)) {
      key_pos <- which(x$key == operation)
      if (length(key_pos) == 0) {
        return(default)
      }
      return(x$value[key_pos[1]])
    }
    if (
      is.list(x) &&
        !is.data.frame(x) &&
        "key" %in% names(x) &&
        "value" %in% names(x)
    ) {
      key_pos <- which(x[["key"]] == operation)
      if (length(key_pos) == 0) {
        return(default)
      }
      return(x[["value"]][key_pos[1]])
    }
    return(default)
  }))
}
