#' DuckDB HTS File Reader Extension for R
#'
#' The Rduckhts package provides an interface to the DuckDB HTS (High Throughput Sequencing)
#' file reader extension from within R. It enables reading common bioinformatics file formats
#' such as VCF/BCF, SAM/BAM/CRAM, FASTA, FASTQ, GFF, GTF, and tabix-indexed files
#' directly from R using SQL queries via DuckDB.
#'
#' @docType package
#' @name Rduckhts
#' @author DuckHTS Contributors
#' @references \url{https://github.com/RGenomicsETL/duckhts}
#' @keywords internal
#' @importFrom DBI dbExecute dbExistsTable dbRemoveTable
#' @importFrom duckdb duckdb
"_PACKAGE"

build_param_str <- function(params) {
  if (length(params) == 0) {
    return("")
  }
  param_str <- paste(paste0(names(params), " := ", params), collapse = ", ")
  if (nchar(param_str) > 0) {
    return(paste0(", ", param_str))
  }
  ""
}

sql_quote_string <- function(x) {
  sprintf("'%s'", gsub("'", "''", x, fixed = TRUE))
}

sql_map_literal <- function(x) {
  if (is.null(x) || length(x) == 0) {
    stop("column_map must be non-empty", call. = FALSE)
  }
  nm <- names(x)
  if (is.null(nm) || any(!nzchar(nm))) {
    stop("column_map must be a named character vector with canonical names as names", call. = FALSE)
  }
  keys <- paste(vapply(nm, sql_quote_string, character(1)), collapse = ", ")
  vals <- paste(vapply(as.character(x), sql_quote_string, character(1)), collapse = ", ")
  sprintf("map([%s], [%s])", keys, vals)
}

read_munge_column_map_file <- function(path) {
  tbl <- utils::read.delim(
    path,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  if (ncol(tbl) < 2) {
    stop("column_map_file must be a two-column TSV with source and canonical names", call. = FALSE)
  }
  out <- structure(tbl[[1]], names = toupper(tbl[[2]]))
  out[nzchar(names(out))]
}

resolve_munge_column_map <- function(raw_map, available_columns) {
  if (is.null(raw_map) || length(raw_map) == 0 || is.null(available_columns) || length(available_columns) == 0) {
    return(character())
  }
  available_columns <- as.character(available_columns)
  available_upper <- toupper(available_columns)
  canonical_names <- unique(names(raw_map))
  selected <- character(length(canonical_names))
  names(selected) <- canonical_names
  for (i in seq_along(canonical_names)) {
    canonical <- canonical_names[[i]]
    candidates <- as.character(raw_map[names(raw_map) == canonical])
    for (candidate in candidates) {
      idx <- match(toupper(candidate), available_upper)
      if (!is.na(idx)) {
        selected[[i]] <- available_columns[[idx]]
        break
      }
    }
  }
  selected[nzchar(selected)]
}

read_munge_preset_map <- function(con, preset) {
  preset_sql <- sql_quote_string(preset)
  out <- DBI::dbGetQuery(con, sprintf("SELECT duckdb_munge_preset_map(%s) AS m", preset_sql))
  if (nrow(out) != 1 || is.null(out$m[[1]])) {
    stop("duckdb_munge: unknown preset", call. = FALSE)
  }
  m <- out$m[[1]]
  if (!is.data.frame(m) || !all(c("key", "value") %in% names(m))) {
    stop("duckdb_munge: failed to read preset map", call. = FALSE)
  }
  structure(as.character(m$value), names = toupper(as.character(m$key)))
}

#' Setup HTSlib Environment
#'
#' Sets the `HTS_PATH` environment variable to point to the bundled htslib
#' plugins directory. This enables remote file access via libcurl plugins
#' (e.g., s3://, gs://, http://) when plugins are available.
#'
#' @param plugins_dir Optional path to the htslib plugins directory. When NULL,
#'   uses the bundled plugins directory if available.
#'
#' @return Invisibly returns the previous value of `HTS_PATH` (or `NA` if unset).
#'
#' @details
#' Call this before querying remote URLs to allow htslib to locate its plugins.
#'
#' @examples
#' \dontrun{
#' setup_hts_env()
#'
#' plugins_path <- tempfile("hts_plugins_")
#' dir.create(plugins_path)
#' setup_hts_env(plugins_dir = plugins_path)
#' unlink(plugins_path, recursive = TRUE)
#' }
#' @export
setup_hts_env <- function(plugins_dir = NULL) {
  old_value <- Sys.getenv("HTS_PATH", unset = NA)
  if (is.null(plugins_dir)) {
    plugins_dir <- duckhts_htslib_plugins_dir()
  }
  if (nzchar(plugins_dir) && dir.exists(plugins_dir)) {
    Sys.setenv(HTS_PATH = plugins_dir)
  }
  invisible(old_value)
}

#' @keywords internal
duckhts_htslib_plugins_dir <- function() {
  ext_dir <- duckhts_extension_dir()
  candidates <- c(
    file.path(ext_dir, "htslib", "libexec", "htslib"),
    file.path(ext_dir, "htslib", "plugins"),
    file.path(ext_dir, "htslib")
  )
  for (p in candidates) {
    if (dir.exists(p)) {
      return(normalizePath(p))
    }
  }
  ""
}

#' Load DuckHTS Extension
#'
#' Loads the DuckHTS extension into a DuckDB connection. This must be called
#' before using any of the HTS reader functions.
#'
#' @details
#' The DuckDB connection must be created with
#' \code{allow_unsigned_extensions = "true"}.
#'
#' @param con A DuckDB connection object
#' @param extension_path Optional path to the duckhts extension file. If NULL,
#'   will try to use the bundled extension.
#'
#' @return TRUE if the extension was loaded successfully
#'
#' @examples
#' library(DBI)
#' library(duckdb)
#'
#' con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
#' rduckhts_load(con)
#' dbDisconnect(con, shutdown = TRUE)
#'
#' @export
rduckhts_load <- function(con, extension_path = NULL) {
  if (is.null(extension_path)) {
    # Try to use the bundled extension build directory
    ext_dir <- duckhts_extension_dir()
    if (is.null(ext_dir)) {
      stop(
        "duckhts_extension directory not found in installed package.",
        call. = FALSE
      )
    }
    extension_path <- file.path(ext_dir, "build", "duckhts.duckdb_extension")
  }

  if (!file.exists(extension_path)) {
    stop(
      "DuckHTS extension not found at: ",
      extension_path,
      "\nThis suggests the package was not built correctly during installation.",
      "\nTry recompiling the package with R CMD INSTALL --preclean Rduckhts"
    )
  }

  # Ensure unsigned extensions are allowed (must be set at connection creation)
  setting <- DBI::dbGetQuery(
    con,
    "SELECT value FROM duckdb_settings() WHERE name = 'allow_unsigned_extensions'"
  )
  if (nrow(setting) == 1 && tolower(setting$value[1]) != "true") {
    stop(
      "DuckDB connection must allow unsigned extensions. Recreate the ",
      "connection with duckdb::duckdb(config = list(allow_unsigned_extensions = \"true\")).",
      call. = FALSE
    )
  }
  DBI::dbExecute(con, "SET enable_progress_bar = false")

  # Load the extension
  result <- DBI::dbExecute(con, sprintf("LOAD '%s'", extension_path))
  return(result == 0)
}

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
#' con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
#' rduckhts_load(con)
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
#' con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
#' rduckhts_load(con)
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
#' con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
#' rduckhts_load(con)
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
    if (is.list(x) && !is.data.frame(x) && "key" %in% names(x) && "value" %in% names(x)) {
      key_pos <- which(x[["key"]] == operation)
      if (length(key_pos) == 0) {
        return(default)
      }
      return(x[["value"]][key_pos[1]])
    }
    return(default)
  }))
}


#' @keywords internal
duckhts_extension_dir <- function() {
  ext_path <- system.file(
    "duckhts_extension",
    package = "Rduckhts",
    mustWork = FALSE
  )
  if (nzchar(ext_path) && dir.exists(ext_path)) {
    return(ext_path)
  }
  return(NULL)
}

#' Create VCF/BCF Table
#'
#' Creates a DuckDB table from a VCF or BCF file using the DuckHTS extension.
#' This follows the RBCFTools pattern of creating a table that can be queried.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param table_name Name for the created table
#' @param path Path to the VCF/BCF file
#' @param region Optional genomic region (e.g., "chr1:1000-2000")
#' @param index_path Optional explicit path to index file (.csi/.tbi)
#' @param tidy_format Logical. If TRUE, FORMAT columns are returned in tidy format
#' @param additional_csq_column_types Optional bcftools-style `PATTERN TYPE`
#'   overrides for CSQ/ANN/BCSQ subfield typing, separated by newlines or `;`
#' @param overwrite Logical. If TRUE, overwrites existing table
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' library(DBI)
#' library(duckdb)
#'
#' con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
#' rduckhts_load(con)
#' bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
#' rduckhts_bcf(con, "variants", bcf_path, overwrite = TRUE)
#' dbGetQuery(con, "SELECT * FROM variants LIMIT 2")
#' dbDisconnect(con, shutdown = TRUE)
#'
#' @export
rduckhts_bcf <- function(
  con,
  table_name,
  path,
  region = NULL,
  index_path = NULL,
  tidy_format = FALSE,
  additional_csq_column_types = NULL,
  overwrite = FALSE
) {
  if (!missing(table_name) && !is.null(table_name)) {
    if (DBI::dbExistsTable(con, table_name) && !overwrite) {
      stop(
        "Table '",
        table_name,
        "' already exists. Use overwrite = TRUE to replace it."
      )
    }
    if (DBI::dbExistsTable(con, table_name)) {
      DBI::dbRemoveTable(con, table_name)
    }
  }

  # Build the CREATE TABLE query
  params <- list()
  if (!is.null(region)) {
    params$region <- sprintf("'%s'", region)
  }
  if (!is.null(index_path)) {
    params$index_path <- sprintf("'%s'", index_path)
  }
  if (tidy_format) {
    params$tidy_format <- "true"
  }
  if (!is.null(additional_csq_column_types)) {
    params$additional_csq_column_types <- sprintf("'%s'", additional_csq_column_types)
  }

  param_str <- build_param_str(params)

  if (!is.null(table_name)) {
    create_query <- sprintf(
      "CREATE TABLE %s AS SELECT * FROM read_bcf('%s'%s)",
      table_name,
      path,
      param_str
    )
  } else {
    create_query <- sprintf(
      "CREATE VIEW bcf_data AS SELECT * FROM read_bcf('%s'%s)",
      path,
      param_str
    )
  }

  DBI::dbExecute(con, create_query)
  invisible(TRUE)
}

#' Create SAM/BAM/CRAM Table
#'
#' Creates a DuckDB table from SAM, BAM, or CRAM files using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param table_name Name for the created table
#' @param path Path to the SAM/BAM/CRAM file
#' @param region Optional genomic region (e.g., "chr1:1000-2000")
#' @param index_path Optional explicit path to index file (.bai/.csi/.crai)
#' @param reference Optional reference file path for CRAM files
#' @param standard_tags Logical. If TRUE, include typed standard SAMtags columns. Default FALSE.
#' @param auxiliary_tags Logical. If TRUE, include AUXILIARY_TAGS map of non-standard tags. Default FALSE.
#' @param sequence_encoding Character. Sequence encoding for the SEQ column:
#'   \code{"string"} (default) returns decoded bases as \code{VARCHAR};
#'   \code{"nt16"} returns raw htslib nt16 4-bit codes as \code{UTINYINT[]}.
#' @param quality_representation Character. Quality representation for the QUAL column:
#'   \code{"string"} (default) returns canonical Phred+33 text;
#'   \code{"phred"} returns raw Phred values as \code{UTINYINT[]}.
#' @param overwrite Logical. If TRUE, overwrites existing table
#'
#' @return Invisible TRUE on success
#'
#' @examples
#' library(DBI)
#' library(duckdb)
#'
#' con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
#' rduckhts_load(con)
#' bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
#' rduckhts_bam(con, "reads", bam_path, overwrite = TRUE)
#' dbGetQuery(con, "SELECT COUNT(*) FROM reads WHERE FLAG & 4 = 0")
#' dbDisconnect(con, shutdown = TRUE)
#'
#' @export
rduckhts_bam <- function(
  con,
  table_name,
  path,
  region = NULL,
  index_path = NULL,
  reference = NULL,
  standard_tags = FALSE,
  auxiliary_tags = FALSE,
  sequence_encoding = NULL,
  quality_representation = NULL,
  overwrite = FALSE
) {
  if (!missing(table_name) && !is.null(table_name)) {
    if (DBI::dbExistsTable(con, table_name) && !overwrite) {
      stop(
        "Table '",
        table_name,
        "' already exists. Use overwrite = TRUE to replace it."
      )
    }
    if (DBI::dbExistsTable(con, table_name)) {
      DBI::dbRemoveTable(con, table_name)
    }
  }

  params <- list()
  if (!is.null(region)) {
    params$region <- sprintf("'%s'", region)
  }
  if (!is.null(index_path)) {
    params$index_path <- sprintf("'%s'", index_path)
  }
  if (!is.null(reference)) {
    params$reference <- sprintf("'%s'", reference)
  }
  if (!is.null(standard_tags)) {
    params$standard_tags <- if (isTRUE(standard_tags)) "true" else "false"
  }
  if (!is.null(auxiliary_tags)) {
    params$auxiliary_tags <- if (isTRUE(auxiliary_tags)) "true" else "false"
  }
  if (!is.null(sequence_encoding)) {
    params$sequence_encoding <- sprintf("'%s'", sequence_encoding)
  }
  if (!is.null(quality_representation)) {
    params$quality_representation <- sprintf("'%s'", quality_representation)
  }

  param_str <- build_param_str(params)

  if (!is.null(table_name)) {
    create_query <- sprintf(
      "CREATE TABLE %s AS SELECT * FROM read_bam('%s'%s)",
      table_name,
      path,
      param_str
    )
  } else {
    create_query <- sprintf(
      "CREATE VIEW bam_data AS SELECT * FROM read_bam('%s'%s)",
      path,
      param_str
    )
  }

  DBI::dbExecute(con, create_query)
  invisible(TRUE)
}

#' Normalize R Data Types to DuckDB Types for Tabix
#'
#' Normalizes R data type names to their corresponding DuckDB types for use with
#' tabix readers. This function handles common R type name variations and maps them
#' to appropriate DuckDB column types.
#'
#' @param types A character vector of R data type names to be normalized.
#'
#' @return A character vector of normalized DuckDB type names suitable for tabix columns.
#'
#' @details
#' The function performs the following normalizations:
#' \itemize{
#' \item Integer types (integer, int, int32, int64) -> BIGINT
#' \item Numeric types (numeric, double, float) -> DOUBLE
#' \item Character types (character, string, chr) -> VARCHAR
#' \item Logical types (logical, bool, boolean) -> BOOLEAN
#' \item Other types -> Converted to uppercase as-is
#' }
#' If an empty vector is provided, it returns the empty vector unchanged.
#'
#' @examples
#' normalize_tabix_types(c("integer", "character", "numeric"))
#' normalize_tabix_types(c("int", "string", "float"))
#'
#' @seealso
#' \code{\link{rduckhts_tabix}} for using normalized types with tabix readers,
#' \code{\link{duckdb_type_mappings}} for the complete type mapping table.
#'
#' @export
normalize_tabix_types <- function(types) {
  if (length(types) == 0) {
    return(types)
  }
  mappings <- duckdb_type_mappings()$r_to_duckdb
  cleaned <- trimws(types)
  lowered <- tolower(cleaned)
  mapped <- character(length(cleaned))
  for (i in seq_along(cleaned)) {
    mapped[i] <- switch(
      lowered[i],
      "integer" = "BIGINT",
      "int" = "BIGINT",
      "int32" = "BIGINT",
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
      toupper(cleaned[i])
    )
  }
  mapped
}

#' Create FASTA Table
#'
#' Creates a DuckDB table from FASTA files using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param table_name Name for the created table
#' @param path Path to the FASTA file
#' @param region Optional genomic region (e.g., "chr1:1000-2000" or "chr1:1-10,chr2:5-20")
#' @param index_path Optional explicit path to FASTA index file (.fai)
#' @param sequence_encoding Character. Sequence encoding for the SEQUENCE column:
#'   \code{"string"} (default) returns decoded bases as \code{VARCHAR};
#'   \code{"nt16"} returns raw htslib nt16 4-bit codes as \code{UTINYINT[]}.
#' @param overwrite Logical. If TRUE, overwrites existing table
#'
#' @return Invisible TRUE on success
#'
#' @export
rduckhts_fasta <- function(
  con,
  table_name,
  path,
  region = NULL,
  index_path = NULL,
  sequence_encoding = NULL,
  overwrite = FALSE
) {
  if (!missing(table_name) && !is.null(table_name)) {
    if (DBI::dbExistsTable(con, table_name) && !overwrite) {
      stop(
        "Table '",
        table_name,
        "' already exists. Use overwrite = TRUE to replace it."
      )
    }
    if (DBI::dbExistsTable(con, table_name)) {
      DBI::dbRemoveTable(con, table_name)
    }
  }

  params <- list()
  if (!is.null(region)) {
    params$region <- sprintf("'%s'", region)
  }
  if (!is.null(index_path)) {
    params$index_path <- sprintf("'%s'", index_path)
  }
  if (!is.null(sequence_encoding)) {
    params$sequence_encoding <- sprintf("'%s'", sequence_encoding)
  }
  param_str <- build_param_str(params)

  if (!is.null(table_name)) {
    create_query <- sprintf(
      "CREATE TABLE %s AS SELECT * FROM read_fasta('%s'%s)",
      table_name,
      path,
      param_str
    )
  } else {
    create_query <- sprintf(
      "CREATE VIEW fasta_data AS SELECT * FROM read_fasta('%s'%s)",
      path,
      param_str
    )
  }

  DBI::dbExecute(con, create_query)
  invisible(TRUE)
}

#' Create BED Table
#'
#' Creates a DuckDB table from a BED file using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param table_name Name for the created table
#' @param path Path to the BED file
#' @param region Optional genomic region for tabix-backed BED queries
#' @param index_path Optional explicit path to a BED tabix index
#' @param overwrite Logical. If TRUE, overwrites an existing table
#'
#' @return Invisible TRUE on success
#'
#' @export
rduckhts_bed <- function(
  con,
  table_name,
  path,
  region = NULL,
  index_path = NULL,
  overwrite = FALSE
) {
  if (!missing(table_name) && !is.null(table_name)) {
    if (DBI::dbExistsTable(con, table_name) && !overwrite) {
      stop(
        "Table '",
        table_name,
        "' already exists. Use overwrite = TRUE to replace it."
      )
    }
    if (DBI::dbExistsTable(con, table_name)) {
      DBI::dbRemoveTable(con, table_name)
    }
  }

  params <- list()
  if (!is.null(region)) {
    params$region <- sprintf("'%s'", region)
  }
  if (!is.null(index_path)) {
    params$index_path <- sprintf("'%s'", index_path)
  }
  param_str <- build_param_str(params)

  if (!is.null(table_name)) {
    create_query <- sprintf(
      "CREATE TABLE %s AS SELECT * FROM read_bed('%s'%s)",
      table_name,
      path,
      param_str
    )
  } else {
    create_query <- sprintf(
      "CREATE VIEW bed_data AS SELECT * FROM read_bed('%s'%s)",
      path,
      param_str
    )
  }

  DBI::dbExecute(con, create_query)
  invisible(TRUE)
}

#' Compute FASTA Interval Nucleotide Composition
#'
#' Computes bedtools nuc-style nucleotide composition over either a BED file or
#' generated fixed-width bins.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param path Path to the FASTA file
#' @param bed_path Optional BED path. Supply exactly one of `bed_path` or `bin_width`.
#' @param bin_width Optional fixed bin width in base pairs
#' @param region Optional FASTA region filter
#' @param index_path Optional explicit FASTA index path
#' @param bed_index_path Optional explicit BED tabix index path
#' @param include_seq Include the fetched interval sequence
#'
#' @return A data frame with interval composition statistics
#'
#' @export
rduckhts_fasta_nuc <- function(
  con,
  path,
  bed_path = NULL,
  bin_width = NULL,
  region = NULL,
  index_path = NULL,
  bed_index_path = NULL,
  include_seq = FALSE
) {
  params <- list()
  if (!is.null(bed_path)) params$bed_path <- sprintf("'%s'", bed_path)
  if (!is.null(bin_width)) params$bin_width <- bin_width
  if (!is.null(region)) params$region <- sprintf("'%s'", region)
  if (!is.null(index_path)) params$index_path <- sprintf("'%s'", index_path)
  if (!is.null(bed_index_path)) params$bed_index_path <- sprintf("'%s'", bed_index_path)
  if (include_seq) params$include_seq <- "true"
  param_str <- build_param_str(params)
  query <- sprintf("SELECT * FROM fasta_nuc('%s'%s)", path, param_str)
  DBI::dbGetQuery(con, query)
}

#' Build FASTA Index
#'
#' Builds a FASTA index (.fai) using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param path Path to the FASTA file
#' @param index_path Optional explicit output path for FASTA index file (.fai)
#'
#' @return A data frame with columns `success` and `index_path`
#'
#' @export
rduckhts_fasta_index <- function(con, path, index_path = NULL) {
  params <- list()
  if (!is.null(index_path)) {
    params$index_path <- sprintf("'%s'", index_path)
  }
  param_str <- build_param_str(params)
  query <- sprintf(
    "SELECT * FROM fasta_index('%s'%s)",
    path,
    param_str
  )
  DBI::dbGetQuery(con, query)
}

#' BGZF Compress a File
#'
#' Compresses a plain file to BGZF using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param path Path to the input file
#' @param output_path Optional explicit output path
#' @param threads BGZF worker thread count
#' @param level Compression level, or -1 for the htslib default
#' @param keep Keep the original input file after compression
#' @param overwrite Overwrite an existing output file
#'
#' @return A data frame describing the created BGZF file
#'
#' @export
rduckhts_bgzip <- function(
  con,
  path,
  output_path = NULL,
  threads = 4,
  level = -1,
  keep = TRUE,
  overwrite = FALSE
) {
  params <- list(threads = threads, level = level)
  if (!is.null(output_path)) params$output_path <- sprintf("'%s'", output_path)
  if (keep) params$keep <- "true"
  if (overwrite) params$overwrite <- "true"
  param_str <- build_param_str(params)
  query <- sprintf("SELECT * FROM bgzip('%s'%s)", path, param_str)
  DBI::dbGetQuery(con, query)
}

#' BGZF Decompress a File
#'
#' Decompresses a BGZF file using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param path Path to the BGZF-compressed input file
#' @param output_path Optional explicit output path
#' @param threads BGZF worker thread count
#' @param keep Keep the compressed input file after decompression
#' @param overwrite Overwrite an existing output file
#'
#' @return A data frame describing the created output file
#'
#' @export
rduckhts_bgunzip <- function(
  con,
  path,
  output_path = NULL,
  threads = 4,
  keep = TRUE,
  overwrite = FALSE
) {
  params <- list(threads = threads)
  if (!is.null(output_path)) params$output_path <- sprintf("'%s'", output_path)
  if (keep) params$keep <- "true"
  if (overwrite) params$overwrite <- "true"
  param_str <- build_param_str(params)
  query <- sprintf("SELECT * FROM bgunzip('%s'%s)", path, param_str)
  DBI::dbGetQuery(con, query)
}

#' Build BAM or CRAM Index
#'
#' Builds a BAM or CRAM index using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param path Path to the input BAM or CRAM file
#' @param index_path Optional explicit output path for the created index
#' @param min_shift Index format selector used by htslib
#' @param threads htslib indexing thread count
#'
#' @return A data frame with `success`, `index_path`, and `index_format`
#'
#' @export
rduckhts_bam_index <- function(
  con,
  path,
  index_path = NULL,
  min_shift = 0,
  threads = 4
) {
  params <- list(min_shift = min_shift, threads = threads)
  if (!is.null(index_path)) params$index_path <- sprintf("'%s'", index_path)
  param_str <- build_param_str(params)
  query <- sprintf("SELECT * FROM bam_index('%s'%s)", path, param_str)
  DBI::dbGetQuery(con, query)
}

#' Build VCF or BCF Index
#'
#' Builds a TBI or CSI index for a VCF/BCF file using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param path Path to the input VCF/BCF file
#' @param index_path Optional explicit output path for the created index
#' @param min_shift Optional explicit min_shift passed to htslib
#' @param threads htslib indexing thread count
#'
#' @return A data frame with `success`, `index_path`, and `index_format`
#'
#' @export
rduckhts_bcf_index <- function(
  con,
  path,
  index_path = NULL,
  min_shift = NULL,
  threads = 4
) {
  params <- list(threads = threads)
  if (!is.null(index_path)) params$index_path <- sprintf("'%s'", index_path)
  if (!is.null(min_shift)) params$min_shift <- min_shift
  param_str <- build_param_str(params)
  query <- sprintf("SELECT * FROM bcf_index('%s'%s)", path, param_str)
  DBI::dbGetQuery(con, query)
}

#' Build Tabix Index
#'
#' Builds a tabix index for a BGZF-compressed text file using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param path Path to the BGZF-compressed input file
#' @param preset Optional preset such as `"vcf"`, `"bed"`, `"gff"`, or `"sam"`
#' @param index_path Optional explicit output path for the created index
#' @param min_shift Index format selector used by htslib
#' @param threads htslib indexing thread count
#' @param seq_col,start_col,end_col Optional explicit tabix coordinate columns
#' @param comment_char Optional tabix comment/header prefix
#' @param skip_lines Optional fixed number of header lines to skip
#'
#' @return A data frame with `success`, `index_path`, and `index_format`
#'
#' @export
rduckhts_tabix_index <- function(
  con,
  path,
  preset = "vcf",
  index_path = NULL,
  min_shift = 0,
  threads = 4,
  seq_col = NULL,
  start_col = NULL,
  end_col = NULL,
  comment_char = NULL,
  skip_lines = NULL
) {
  params <- list(
    preset = sprintf("'%s'", preset),
    min_shift = min_shift,
    threads = threads
  )
  if (!is.null(index_path)) params$index_path <- sprintf("'%s'", index_path)
  if (!is.null(seq_col)) params$seq_col <- seq_col
  if (!is.null(start_col)) params$start_col <- start_col
  if (!is.null(end_col)) params$end_col <- end_col
  if (!is.null(comment_char)) params$comment_char <- sprintf("'%s'", comment_char)
  if (!is.null(skip_lines)) params$skip_lines <- skip_lines
  param_str <- build_param_str(params)
  query <- sprintf("SELECT * FROM tabix_index('%s'%s)", path, param_str)
  DBI::dbGetQuery(con, query)
}

#' Create FASTQ Table
#'
#' Creates a DuckDB table from FASTQ files using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param table_name Name for the created table
#' @param path Path to the FASTQ file
#' @param mate_path Optional path to mate file for paired reads
#' @param interleaved Logical indicating if file is interleaved paired reads
#' @param sequence_encoding Character. Sequence encoding for the SEQUENCE column:
#'   \code{"string"} (default) returns decoded bases as \code{VARCHAR};
#'   \code{"nt16"} returns raw htslib nt16 4-bit codes as \code{UTINYINT[]}.
#' @param quality_representation Character. Quality representation for the QUALITY column:
#'   \code{"string"} (default) returns canonical Phred+33 text;
#'   \code{"phred"} returns raw Phred values as \code{UTINYINT[]}.
#' @param input_quality_encoding Character. Input FASTQ quality encoding:
#'   \code{"phred33"} (default FASTQ convention), \code{"auto"}, \code{"phred64"},
#'   or \code{"solexa64"}.
#' @param overwrite Logical. If TRUE, overwrites existing table
#'
#' @return Invisible TRUE on success
#'
#' @export
rduckhts_fastq <- function(
  con,
  table_name,
  path,
  mate_path = NULL,
  interleaved = FALSE,
  sequence_encoding = NULL,
  quality_representation = NULL,
  input_quality_encoding = NULL,
  overwrite = FALSE
) {
  if (!missing(table_name) && !is.null(table_name)) {
    if (DBI::dbExistsTable(con, table_name) && !overwrite) {
      stop(
        "Table '",
        table_name,
        "' already exists. Use overwrite = TRUE to replace it."
      )
    }
    if (DBI::dbExistsTable(con, table_name)) {
      DBI::dbRemoveTable(con, table_name)
    }
  }

  params <- list()
  if (!is.null(mate_path)) {
    params$mate_path <- sprintf("'%s'", mate_path)
  }
  if (interleaved) {
    params$interleaved <- "true"
  }
  if (!is.null(sequence_encoding)) {
    params$sequence_encoding <- sprintf("'%s'", sequence_encoding)
  }
  if (!is.null(quality_representation)) {
    params$quality_representation <- sprintf("'%s'", quality_representation)
  }
  if (!is.null(input_quality_encoding)) {
    params$input_quality_encoding <- sprintf("'%s'", input_quality_encoding)
  }

  param_str <- build_param_str(params)

  if (!is.null(table_name)) {
    create_query <- sprintf(
      "CREATE TABLE %s AS SELECT * FROM read_fastq('%s'%s)",
      table_name,
      path,
      param_str
    )
  } else {
    create_query <- sprintf(
      "CREATE VIEW fastq_data AS SELECT * FROM read_fastq('%s'%s)",
      path,
      param_str
    )
  }

  DBI::dbExecute(con, create_query)
  invisible(TRUE)
}

#' Detect FASTQ Quality Encoding
#'
#' Inspects a FASTQ file's observed quality ASCII range and reports compatible
#' legacy encodings with a heuristic guessed encoding.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param path Path to the FASTQ file
#' @param max_records Maximum number of records to inspect
#'
#' @return A data frame with the detected quality encoding summary
#'
#' @export
rduckhts_detect_quality_encoding <- function(con, path, max_records = 10000) {
  params <- list(max_records = max_records)
  param_str <- build_param_str(params)
  query <- sprintf(
    "SELECT * FROM detect_quality_encoding('%s'%s)",
    path,
    param_str
  )
  DBI::dbGetQuery(con, query)
}

#' Create GFF3 Table
#'
#' Creates a DuckDB table from GFF3 files using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param table_name Name for the created table
#' @param path Path to the GFF3 file
#' @param region Optional genomic region (e.g., "chr1:1000-2000")
#' @param index_path Optional explicit path to index file (.tbi/.csi)
#' @param header Logical. If TRUE, use first non-meta line as column names
#' @param header_names Character vector to override column names
#' @param auto_detect Logical. If TRUE, infer basic numeric column types
#' @param column_types Character vector of column types (e.g. "BIGINT", "VARCHAR")
#' @param attributes_map Logical. If TRUE, returns attributes as a MAP column
#' @param overwrite Logical. If TRUE, overwrites existing table
#'
#' @return Invisible TRUE on success
#'
#' @export
rduckhts_gff <- function(
  con,
  table_name,
  path,
  region = NULL,
  index_path = NULL,
  header = NULL,
  header_names = NULL,
  auto_detect = NULL,
  column_types = NULL,
  attributes_map = FALSE,
  overwrite = FALSE
) {
  if (!missing(table_name) && !is.null(table_name)) {
    if (DBI::dbExistsTable(con, table_name) && !overwrite) {
      stop(
        "Table '",
        table_name,
        "' already exists. Use overwrite = TRUE to replace it."
      )
    }
    if (DBI::dbExistsTable(con, table_name)) {
      DBI::dbRemoveTable(con, table_name)
    }
  }

  params <- list()
  if (!is.null(region)) {
    params$region <- sprintf("'%s'", region)
  }
  if (!is.null(index_path)) {
    params$index_path <- sprintf("'%s'", index_path)
  }
  if (!is.null(header)) {
    params$header <- if (isTRUE(header)) "true" else "false"
  }
  if (!is.null(auto_detect)) {
    params$auto_detect <- if (isTRUE(auto_detect)) "true" else "false"
  }
  if (!is.null(header_names)) {
    if (!is.character(header_names)) {
      stop("header_names must be a character vector")
    }
    params$header_names <- sprintf(
      "[%s]",
      paste(sprintf("'%s'", header_names), collapse = ", ")
    )
  }
  if (!is.null(column_types)) {
    if (!is.character(column_types)) {
      stop("column_types must be a character vector")
    }
    normalized_types <- normalize_tabix_types(column_types)
    params$column_types <- sprintf(
      "[%s]",
      paste(sprintf("'%s'", normalized_types), collapse = ", ")
    )
  }
  if (attributes_map) {
    params$attributes_map <- "true"
  }

  param_str <- build_param_str(params)

  if (!is.null(table_name)) {
    create_query <- sprintf(
      "CREATE TABLE %s AS SELECT * FROM read_gff('%s'%s)",
      table_name,
      path,
      param_str
    )
  } else {
    create_query <- sprintf(
      "CREATE VIEW gff_data AS SELECT * FROM read_gff('%s'%s)",
      path,
      param_str
    )
  }

  DBI::dbExecute(con, create_query)
  invisible(TRUE)
}

#' Create GTF Table
#'
#' Creates a DuckDB table from GTF files using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param table_name Name for the created table
#' @param path Path to the GTF file
#' @param region Optional genomic region (e.g., "chr1:1000-2000")
#' @param index_path Optional explicit path to index file (.tbi/.csi)
#' @param header Logical. If TRUE, use first non-meta line as column names
#' @param header_names Character vector to override column names
#' @param auto_detect Logical. If TRUE, infer basic numeric column types
#' @param column_types Character vector of column types (e.g. "BIGINT", "VARCHAR")
#' @param attributes_map Logical. If TRUE, returns attributes as a MAP column
#' @param overwrite Logical. If TRUE, overwrites existing table
#'
#' @return Invisible TRUE on success
#'
#' @export
rduckhts_gtf <- function(
  con,
  table_name,
  path,
  region = NULL,
  index_path = NULL,
  header = NULL,
  header_names = NULL,
  auto_detect = NULL,
  column_types = NULL,
  attributes_map = FALSE,
  overwrite = FALSE
) {
  if (!missing(table_name) && !is.null(table_name)) {
    if (DBI::dbExistsTable(con, table_name) && !overwrite) {
      stop(
        "Table '",
        table_name,
        "' already exists. Use overwrite = TRUE to replace it."
      )
    }
    if (DBI::dbExistsTable(con, table_name)) {
      DBI::dbRemoveTable(con, table_name)
    }
  }

  params <- list()
  if (!is.null(region)) {
    params$region <- sprintf("'%s'", region)
  }
  if (!is.null(index_path)) {
    params$index_path <- sprintf("'%s'", index_path)
  }
  if (!is.null(header)) {
    params$header <- if (isTRUE(header)) "true" else "false"
  }
  if (!is.null(auto_detect)) {
    params$auto_detect <- if (isTRUE(auto_detect)) "true" else "false"
  }
  if (!is.null(header_names)) {
    if (!is.character(header_names)) {
      stop("header_names must be a character vector")
    }
    params$header_names <- sprintf(
      "[%s]",
      paste(sprintf("'%s'", header_names), collapse = ", ")
    )
  }
  if (!is.null(column_types)) {
    if (!is.character(column_types)) {
      stop("column_types must be a character vector")
    }
    normalized_types <- normalize_tabix_types(column_types)
    params$column_types <- sprintf(
      "[%s]",
      paste(sprintf("'%s'", normalized_types), collapse = ", ")
    )
  }
  if (attributes_map) {
    params$attributes_map <- "true"
  }

  param_str <- build_param_str(params)

  if (!is.null(table_name)) {
    create_query <- sprintf(
      "CREATE TABLE %s AS SELECT * FROM read_gtf('%s'%s)",
      table_name,
      path,
      param_str
    )
  } else {
    create_query <- sprintf(
      "CREATE VIEW gtf_data AS SELECT * FROM read_gtf('%s'%s)",
      path,
      param_str
    )
  }

  DBI::dbExecute(con, create_query)
  invisible(TRUE)
}

#' Create Tabix-Indexed File Table
#'
#' Creates a DuckDB table from any tabix-indexed file using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param table_name Name for the created table
#' @param path Path to the tabix-indexed file
#' @param region Optional genomic region (e.g., "chr1:1000-2000")
#' @param index_path Optional explicit path to index file (.tbi/.csi)
#' @param header Logical. If TRUE, use first non-meta line as column names
#' @param header_names Character vector to override column names
#' @param auto_detect Logical. If TRUE, infer basic numeric column types
#' @param column_types Character vector of column types (e.g. "BIGINT", "VARCHAR")
#' @param overwrite Logical. If TRUE, overwrites existing table
#'
#' @return Invisible TRUE on success
#'
#' @export
rduckhts_tabix <- function(
  con,
  table_name,
  path,
  region = NULL,
  index_path = NULL,
  header = NULL,
  header_names = NULL,
  auto_detect = NULL,
  column_types = NULL,
  overwrite = FALSE
) {
  if (!missing(table_name) && !is.null(table_name)) {
    if (DBI::dbExistsTable(con, table_name) && !overwrite) {
      stop(
        "Table '",
        table_name,
        "' already exists. Use overwrite = TRUE to replace it."
      )
    }
    if (DBI::dbExistsTable(con, table_name)) {
      DBI::dbRemoveTable(con, table_name)
    }
  }

  params <- list()
  if (!is.null(region)) {
    params$region <- sprintf("'%s'", region)
  }
  if (!is.null(index_path)) {
    params$index_path <- sprintf("'%s'", index_path)
  }
  if (!is.null(header)) {
    params$header <- if (isTRUE(header)) "true" else "false"
  }
  if (!is.null(auto_detect)) {
    params$auto_detect <- if (isTRUE(auto_detect)) "true" else "false"
  }
  if (!is.null(header_names)) {
    if (!is.character(header_names)) {
      stop("header_names must be a character vector")
    }
    params$header_names <- sprintf(
      "[%s]",
      paste(sprintf("'%s'", header_names), collapse = ", ")
    )
  }
  if (!is.null(column_types)) {
    if (!is.character(column_types)) {
      stop("column_types must be a character vector")
    }
    normalized_types <- normalize_tabix_types(column_types)
    params$column_types <- sprintf(
      "[%s]",
      paste(sprintf("'%s'", normalized_types), collapse = ", ")
    )
  }
  param_str <- build_param_str(params)

  if (!is.null(table_name)) {
    create_query <- sprintf(
      "CREATE TABLE %s AS SELECT * FROM read_tabix('%s'%s)",
      table_name,
      path,
      param_str
    )
  } else {
    create_query <- sprintf(
      "CREATE VIEW tabix_data AS SELECT * FROM read_tabix('%s'%s)",
      path,
      param_str
    )
  }

  DBI::dbExecute(con, create_query)
  invisible(TRUE)
}

#' Read HTS Header Metadata
#'
#' Reads file header records from HTS-supported formats using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param path Path to input HTS file
#' @param format Optional format hint (e.g., "auto", "vcf", "bcf", "bam", "cram", "tabix")
#' @param mode Header output mode: "parsed" (default), "raw", or "both"
#'
#' @return A data frame with parsed header metadata.
#'
#' @export
rduckhts_hts_header <- function(con, path, format = NULL, mode = NULL) {
  params <- list()
  if (!is.null(format)) {
    params$format <- sprintf("'%s'", format)
  }
  if (!is.null(mode)) {
    params$mode <- sprintf("'%s'", mode)
  }
  param_str <- build_param_str(params)
  query <- sprintf("SELECT * FROM read_hts_header('%s'%s)", path, param_str)
  DBI::dbGetQuery(con, query)
}

#' Read HTS Index Metadata
#'
#' Reads index metadata from HTS-supported index files via DuckHTS.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param path Path to input HTS file
#' @param format Optional format hint (e.g., "auto", "vcf", "bcf", "bam", "cram", "tabix")
#' @param index_path Optional explicit path to index file
#'
#' @return A data frame with index metadata.
#'
#' @export
rduckhts_hts_index <- function(con, path, format = NULL, index_path = NULL) {
  params <- list()
  if (!is.null(format)) {
    params$format <- sprintf("'%s'", format)
  }
  if (!is.null(index_path)) {
    params$index_path <- sprintf("'%s'", index_path)
  }
  param_str <- build_param_str(params)
  query <- sprintf("SELECT * FROM read_hts_index('%s'%s)", path, param_str)
  DBI::dbGetQuery(con, query)
}

#' Read HTS Index Spans
#'
#' Returns index span-oriented metadata for planning range workloads.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param path Path to input HTS file
#' @param format Optional format hint
#' @param index_path Optional explicit path to index file
#'
#' @return A data frame with span-oriented index metadata.
#'
#' @export
rduckhts_hts_index_spans <- function(con, path, format = NULL, index_path = NULL) {
  params <- list()
  if (!is.null(format)) {
    params$format <- sprintf("'%s'", format)
  }
  if (!is.null(index_path)) {
    params$index_path <- sprintf("'%s'", index_path)
  }
  param_str <- build_param_str(params)
  query <- sprintf("SELECT * FROM read_hts_index_spans('%s'%s)", path, param_str)
  DBI::dbGetQuery(con, query)
}

#' Read Raw HTS Index Blob
#'
#' Returns raw index metadata blob data for a file index.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param path Path to input HTS file
#' @param format Optional format hint
#' @param index_path Optional explicit path to index file
#'
#' @return A data frame with raw index blob metadata.
#'
#' @export
rduckhts_hts_index_raw <- function(con, path, format = NULL, index_path = NULL) {
  params <- list()
  if (!is.null(format)) {
    params$format <- sprintf("'%s'", format)
  }
  if (!is.null(index_path)) {
    params$index_path <- sprintf("'%s'", index_path)
  }
  param_str <- build_param_str(params)
  query <- sprintf("SELECT * FROM read_hts_index_raw('%s'%s)", path, param_str)
  DBI::dbGetQuery(con, query)
}

#' Lift Over Variant Coordinates Against a Query
#'
#' Applies the DuckHTS `duckdb_liftover(...)` table macro to rows from a SQL query or
#' table expression with chromosome and position columns, plus optional reference
#' and alternate alleles.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param query SQL query or table expression to lift over
#' @param chain_path Path to a UCSC chain file
#' @param dst_fasta_ref Path to the destination FASTA reference
#' @param chrom_col Source chromosome column name
#' @param pos_col Source 1-based position column name
#' @param ref_col Optional reference allele column name
#' @param alt_col Optional alternate allele column name
#' @param src_fasta_ref Optional source FASTA reference
#' @param max_snp_gap Maximum chain block merge gap
#' @param max_indel_inc Maximum indel anchor expansion
 #' @param lift_mt If FALSE (default), mitochondrial variants with matching
#'   source/destination contig lengths are passed through with only contig
#'   rename. If TRUE, MT variants are lifted through the chain like any
#'   other contig.
#' @param end_pos_col Optional column name containing INFO/END positions
#'   (1-based) to lift alongside the primary position. When provided, the
#'   output includes a `dest_end` column with the lifted end position.
#' @param no_left_align If FALSE (default), lifted indels are left-aligned
#'   against the destination reference. Set TRUE to skip left-alignment,
#'   mirroring \code{--no-left-align} in \code{bcftools +liftover}.
#'
#' @return A data frame with source columns, lifted coordinates/alleles, and warnings.
#'
#' @export
rduckhts_liftover <- function(
  con,
  query,
  chain_path,
  dst_fasta_ref,
  chrom_col = "chrom",
  pos_col = "pos",
  ref_col = NULL,
  alt_col = NULL,
  src_fasta_ref = NULL,
  max_snp_gap = 1,
  max_indel_inc = 250,
  lift_mt = FALSE,
  end_pos_col = NULL,
  no_left_align = FALSE
) {
  if (!is.numeric(max_snp_gap) || length(max_snp_gap) != 1 || is.na(max_snp_gap) || max_snp_gap < 0) {
    stop("max_snp_gap must be >= 0", call. = FALSE)
  }
  if (!is.numeric(max_indel_inc) || length(max_indel_inc) != 1 || is.na(max_indel_inc) || max_indel_inc < 0) {
    stop("max_indel_inc must be >= 0", call. = FALSE)
  }

  table_expr <- query
  if (grepl("^\\s*select\\b", table_expr, ignore.case = TRUE)) {
    table_expr <- sprintf("(%s) AS duckhts_src", table_expr)
  }
  table_sql <- gsub("'", "''", table_expr, fixed = TRUE)
  params <- list(
    sprintf("'%s'", table_sql),
    sprintf("'%s'", chrom_col),
    sprintf("'%s'", pos_col)
  )
  if (!is.null(ref_col)) params <- c(params, sprintf("ref_col := '%s'", ref_col))
  if (!is.null(alt_col)) params <- c(params, sprintf("alt_col := '%s'", alt_col))
  params <- c(
    params,
    sprintf("chain_path := '%s'", chain_path),
    sprintf("dst_fasta_ref := '%s'", dst_fasta_ref)
  )
  if (!is.null(src_fasta_ref)) params <- c(params, sprintf("src_fasta_ref := '%s'", src_fasta_ref))
  params <- c(
    params,
    sprintf("max_snp_gap := %d", as.integer(max_snp_gap)),
    sprintf("max_indel_inc := %d", as.integer(max_indel_inc)),
    sprintf("lift_mt := %s", tolower(as.character(lift_mt)))
  )
  if (!is.null(end_pos_col)) params <- c(params, sprintf("end_pos_col := '%s'", end_pos_col))
  params <- c(params, sprintf("no_left_align := %s", tolower(as.character(no_left_align))))
  sql <- sprintf("SELECT * FROM duckdb_liftover(%s)", paste(params, collapse = ", "))
  DBI::dbGetQuery(con, sql)
}

#' Munge Summary Statistics Rows
#'
#' Applies the DuckHTS `duckdb_munge(...)` table macro to rows from a SQL query or
#' table expression, using either an upstream-style preset, a named column map,
#' or a two-column mapping file. When no mapping mode is provided, the bundled
#' `colheaders.tsv` alias file is used by default.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param query SQL query or table expression to normalize
#' @param fasta_ref Path to the reference FASTA. When NULL (default), operates
#'   in fai-only mode: alleles pass through as-is without reference matching
#'   or allele swapping, matching upstream `--fai`-only behavior.
#' @param preset Optional preset such as `"PLINK"`, `"PLINK2"`, `"REGENIE"`,
#'   `"SAIGE"`, `"BOLT"`, `"METAL"`, `"PGS"`, or `"SSF"`
#' @param column_map Optional named character vector mapping canonical munge names
#'   such as `"CHR"`, `"BP"`, `"A1"`, `"A2"` to source column names
#' @param column_map_file Optional path to a two-column TSV mapping file in the
#'   upstream `source<TAB>canonical` format
#' @param iffy_tag FILTER tag for ambiguous reference resolution
#' @param mismatch_tag FILTER tag for reference mismatches
#' @param ns,nc,ne Optional global overrides for sample counts
#'
#' @return A data frame with normalized GWAS-VCF-style variant/effect columns.
#'
#' @export
rduckhts_munge <- function(
  con,
  query,
  fasta_ref = NULL,
  preset = NULL,
  column_map = NULL,
  column_map_file = NULL,
  iffy_tag = "IFFY",
  mismatch_tag = "REF_MISMATCH",
  ns = NULL,
  nc = NULL,
  ne = NULL
) {
  table_expr <- query
  if (grepl("^\\s*select\\b", table_expr, ignore.case = TRUE)) {
    table_expr <- sprintf("(%s) AS duckhts_src", table_expr)
  }
  provided_modes <- sum(!vapply(list(preset, column_map, column_map_file), is.null, logical(1)))
  if (provided_modes > 1) {
    stop("duckdb_munge: specify only one of preset, column_map, or column_map_file", call. = FALSE)
  }

  available_columns <- names(DBI::dbGetQuery(con, sprintf("SELECT * FROM %s LIMIT 0", table_expr)))

  if (!is.null(preset)) {
    default_map_path <- system.file("extdata", "colheaders.tsv", package = "Rduckhts", mustWork = TRUE)
    alias_map <- resolve_munge_column_map(read_munge_column_map_file(default_map_path), available_columns)
    preset_map <- resolve_munge_column_map(read_munge_preset_map(con, preset), available_columns)
    column_map <- c(alias_map, preset_map)
    column_map <- column_map[!duplicated(names(column_map), fromLast = TRUE)]
    preset <- NULL
  }
  if (is.null(preset) && is.null(column_map) && is.null(column_map_file)) {
    column_map_file <- system.file("extdata", "colheaders.tsv", package = "Rduckhts", mustWork = TRUE)
  }
  if (!is.null(column_map_file)) {
    file_map <- resolve_munge_column_map(read_munge_column_map_file(column_map_file), available_columns)
    if (!is.null(column_map)) {
      column_map <- c(file_map, column_map)
      column_map <- column_map[!duplicated(names(column_map), fromLast = TRUE)]
    } else {
      column_map <- file_map
    }
    column_map_file <- NULL
  }
  params <- list(sql_quote_string(table_expr))
  if (!is.null(preset)) params <- c(params, sprintf("preset := '%s'", preset))
  if (!is.null(column_map)) params <- c(params, sprintf("column_map := %s", sql_map_literal(column_map)))
  if (!is.null(column_map_file)) params <- c(params, sprintf("column_map_file := '%s'", column_map_file))
  if (!is.null(fasta_ref)) {
    params <- c(params, sprintf("fasta_ref := '%s'", fasta_ref))
  }
  params <- c(
    params,
    sprintf("iffy_tag := '%s'", iffy_tag),
    sprintf("mismatch_tag := '%s'", mismatch_tag)
  )
  if (!is.null(ns)) params <- c(params, sprintf("ns := %s", format(ns, scientific = FALSE, trim = TRUE)))
  if (!is.null(nc)) params <- c(params, sprintf("nc := %s", format(nc, scientific = FALSE, trim = TRUE)))
  if (!is.null(ne)) params <- c(params, sprintf("ne := %s", format(ne, scientific = FALSE, trim = TRUE)))
  metal_keys <- c("INFO", "HET_I2", "HET_P", "HET_LP", "DIRE")
  has_metal <- !is.null(column_map) && any(metal_keys %in% names(column_map))
  macro_name <- if (has_metal) "duckdb_munge_metal" else "duckdb_munge"
  sql <- sprintf("SELECT * FROM %s(%s)", macro_name, paste(params, collapse = ", "))
  DBI::dbGetQuery(con, sql)
}

#' Compute Polygenic Scores
#'
#' Calls the DuckHTS `bcftools_score(...)` table function to compute sample-level
#' polygenic scores from one genotype VCF/BCF file and one summary-statistics file.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param bcf_path Path to genotype VCF/BCF file
#' @param summary_path Path to summary-statistics file
#' @param use Optional dosage source (`"GT"`, `"DS"`, `"HDS"`, `"AP"`, `"GP"`, `"AS"`)
#' @param columns Optional summary preset (`"PLINK"`, `"PLINK2"`, `"REGENIE"`, `"SAIGE"`,
#'   `"BOLT"`, `"METAL"`, `"PGS"`, `"SSF"`, `"GWAS-SSF"`)
#' @param columns_file Optional two-column summary header mapping file
#' @param q_score_thr Optional comma-separated p-value thresholds (e.g. `"1e-8,1e-6,1e-4"`)
#' @param use_variant_id Logical; if TRUE, match variants by ID instead of CHR+BP
#' @param counts Logical; if TRUE, include per-threshold matched-variant counts
#' @param samples Optional comma-separated list of sample names to subset (e.g. `"SAMP1,SAMP2"`)
#' @param force_samples Logical; if TRUE, ignore missing samples instead of erroring
#' @param regions Optional comma-separated region list (e.g. `"1:1000-2000,2:50-90"`)
#' @param regions_file Optional path to a regions file
#' @param regions_overlap Overlap mode for regions (`0`, `1`, or `2`). Default 1 (trim to region).
#' @param targets Optional comma-separated targets list
#' @param targets_file Optional path to a targets file
#' @param targets_overlap Overlap mode for targets (`0`, `1`, or `2`). Default 0 (record must start in region).
#' @param apply_filters Optional comma-separated FILTER names to keep (e.g. `"PASS,."`)
#' @param include Optional site expression (currently unsupported)
#' @param exclude Optional site expression (currently unsupported)
#'
#' @return A data frame with one row per sample and score/count columns.
#'
#' @export
rduckhts_score <- function(
  con,
  bcf_path,
  summary_path,
  use = NULL,
  columns = "PLINK",
  columns_file = NULL,
  q_score_thr = NULL,
  use_variant_id = FALSE,
  counts = FALSE,
  samples = NULL,
  force_samples = FALSE,
  regions = NULL,
  regions_file = NULL,
  regions_overlap = 1,
  targets = NULL,
  targets_file = NULL,
  targets_overlap = 0,
  apply_filters = NULL,
  include = NULL,
  exclude = NULL
) {
  if (!is.null(use)) {
    use <- toupper(use)
    if (!(use %in% c("GT", "DS", "HDS", "AP", "GP", "AS"))) {
      stop("use must be one of GT, DS, HDS, AP, GP, AS", call. = FALSE)
    }
  }
  if (!is.null(columns)) {
    columns <- toupper(columns)
    if (!(columns %in% c("PLINK", "PLINK2", "REGENIE", "SAIGE", "BOLT", "METAL", "PGS", "SSF", "GWAS-SSF"))) {
      stop("columns must be one of PLINK, PLINK2, REGENIE, SAIGE, BOLT, METAL, PGS, SSF, GWAS-SSF", call. = FALSE)
    }
  }
  params <- list(
    sql_quote_string(bcf_path),
    sql_quote_string(summary_path)
  )
  if (!is.null(use)) params <- c(params, sprintf("use := '%s'", use))
  if (!is.null(columns)) params <- c(params, sprintf("columns := '%s'", columns))
  if (!is.null(columns_file)) params <- c(params, sprintf("columns_file := '%s'", columns_file))
  if (!is.null(q_score_thr)) params <- c(params, sprintf("q_score_thr := '%s'", q_score_thr))
  if (!is.null(samples)) params <- c(params, sprintf("samples := '%s'", samples))
  if (!is.null(regions)) params <- c(params, sprintf("regions := '%s'", regions))
  if (!is.null(regions_file)) params <- c(params, sprintf("regions_file := '%s'", regions_file))
  if (!is.null(regions_overlap)) params <- c(params, sprintf("regions_overlap := %d", as.integer(regions_overlap)))
  if (!is.null(targets)) params <- c(params, sprintf("targets := '%s'", targets))
  if (!is.null(targets_file)) params <- c(params, sprintf("targets_file := '%s'", targets_file))
  if (!is.null(targets_overlap)) params <- c(params, sprintf("targets_overlap := %d", as.integer(targets_overlap)))
  if (!is.null(apply_filters)) params <- c(params, sprintf("apply_filters := '%s'", apply_filters))
  if (!is.null(include)) params <- c(params, sprintf("include := '%s'", include))
  if (!is.null(exclude)) params <- c(params, sprintf("exclude := '%s'", exclude))
  params <- c(
    params,
    sprintf("use_variant_id := %s", if (isTRUE(use_variant_id)) "true" else "false"),
    sprintf("counts := %s", if (isTRUE(counts)) "true" else "false"),
    sprintf("force_samples := %s", if (isTRUE(force_samples)) "true" else "false")
  )
  sql <- sprintf("SELECT * FROM bcftools_score(%s)", paste(params, collapse = ", "))
  DBI::dbGetQuery(con, sql)
}
