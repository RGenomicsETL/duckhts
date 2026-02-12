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
#' # Use bundled plugins if present
#' setup_hts_env()
#'
#' # Or set an explicit plugins directory
#' # plugins_path <- "/path/to/htslib/plugins"
#' # setup_hts_env(plugins_dir = plugins_path)
#'
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
#' @param tidy_format Logical. If TRUE, FORMAT columns are returned in tidy format
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
  tidy_format = FALSE,
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
  if (tidy_format) {
    params$tidy_format <- "true"
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
#' @param reference Optional reference file path for CRAM files
#' @param standard_tags Logical. If TRUE, include typed standard SAMtags columns
#' @param auxiliary_tags Logical. If TRUE, include AUXILIARY_TAGS map of non-standard tags
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
  reference = NULL,
  standard_tags = NULL,
  auxiliary_tags = NULL,
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
  if (!is.null(reference)) {
    params$reference <- sprintf("'%s'", reference)
  }
  if (!is.null(standard_tags)) {
    params$standard_tags <- if (isTRUE(standard_tags)) "true" else "false"
  }
  if (!is.null(auxiliary_tags)) {
    params$auxiliary_tags <- if (isTRUE(auxiliary_tags)) "true" else "false"
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
#' # Normalize mixed type names
#' normalize_tabix_types(c("integer", "character", "numeric"))
#' # Returns: c("BIGINT", "VARCHAR", "DOUBLE")
#'
#' # Handle variations
#' normalize_tabix_types(c("int", "string", "float"))
#' # Returns: c("BIGINT", "VARCHAR", "DOUBLE")
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
#' @param overwrite Logical. If TRUE, overwrites existing table
#'
#' @return Invisible TRUE on success
#'
#' @export
rduckhts_fasta <- function(con, table_name, path, overwrite = FALSE) {
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

  if (!is.null(table_name)) {
    create_query <- sprintf(
      "CREATE TABLE %s AS SELECT * FROM read_fasta('%s')",
      table_name,
      path
    )
  } else {
    create_query <- sprintf(
      "CREATE VIEW fasta_data AS SELECT * FROM read_fasta('%s')",
      path
    )
  }

  DBI::dbExecute(con, create_query)
  invisible(TRUE)
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

#' Create GFF3 Table
#'
#' Creates a DuckDB table from GFF3 files using the DuckHTS extension.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param table_name Name for the created table
#' @param path Path to the GFF3 file
#' @param region Optional genomic region (e.g., "chr1:1000-2000")
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
