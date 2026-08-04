#' DuckDB to R Type Mappings
#'
#' Returns a named list mapping between DuckDB and R data types.
#' This is useful for understanding type conversions when reading
#' HTS files or when specifying column types in tabix functions.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{duckdb_to_r}{Named character vector mapping DuckDB types to R types}
#'   \item{r_to_duckdb}{Named character vector mapping R types to DuckDB types}
#' }
#'
#' @description
#' The mapping covers the most common data types used in HTS file processing:
#' \itemize{
#'   \item BIGINT <-> double (not integer due to 64-bit overflow protection)
#'   \item DOUBLE <-> numeric/double
#'   \item VARCHAR <-> character/string
#'   \item BOOLEAN <-> logical
#'   \item ARRAY types (e.g., VARCHAR[], BIGINT[]) <-> list
#'   \item MAP types (e.g., MAP(VARCHAR, VARCHAR)) <-> data.frame
#' }
#'
#' Important notes:
#' \itemize{
#'   \item 64-bit integers (BIGINT, UBIGINT) become double to prevent overflow
#'   \item DATE/TIME values return as Unix epoch numbers (double)
#'   \item MAP types become data frames with 'key' and 'value' columns
#'   \item ARRAY types become vectors (which are lists in R terminology)
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
    "BIGINT" = "double", # Returns double due to 64-bit integer overflow protection
    "DOUBLE" = "double", # Returns double, not numeric
    "VARCHAR" = "character",
    "BOOLEAN" = "logical",
    "DATE" = "double", # Returns double (Unix epoch days)
    "TIMESTAMP" = "double", # Returns double (Unix epoch seconds)
    "TIME" = "double", # Returns double (nanoseconds since midnight)
    "BLOB" = "list", # Returns list, not raw
    # Additional integer types
    "INTEGER" = "integer",
    "FLOAT" = "double",
    "SMALLINT" = "integer",
    "USMALLINT" = "integer",
    "UBIGINT" = "double", # Returns double due to 64-bit overflow protection
    # Array types - DuckDB arrays become R lists
    "VARCHAR[]" = "list",
    "BIGINT[]" = "list",
    "INTEGER[]" = "list",
    "DOUBLE[]" = "list",
    "FLOAT[]" = "list",
    "BOOLEAN[]" = "list",
    # MAP types - DuckDB maps become R data frames
    "MAP" = "data.frame"
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
#' Identifies columns in a DuckDB table that contain complex types
#' (ARRAY or MAP) that will be returned as R lists.
#'
#' @param con A DuckDB connection
#' @param table_name Name of the table to analyze
#'
#' @return A data frame with columns that have complex types, showing
#'   column_name, column_type, and a description of R type.
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
  schema <- DBI::dbGetQuery(con, sprintf("DESCRIBE %s", table_name))

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
#' Helper function to safely extract elements from DuckDB arrays
#' (returned as R lists) with proper error handling.
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
#' Helper function to work with DuckDB MAP data (returned as data frames).
#' Can extract keys, values, or search for specific key-value pairs.
#'
#' @param map_col A data frame column from DuckDB MAP data
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
