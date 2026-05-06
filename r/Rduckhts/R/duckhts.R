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

.validate_nonnegative_integer_param <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      value < 0 || value > .Machine$integer.max || value != floor(value)) {
    stop(
      sprintf(
        "%s must be a single whole number between 0 and %d",
        name, .Machine$integer.max
      ),
      call. = FALSE
    )
  }
  as.integer(value)
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
#' @param decompression_threads Integer. Number of htslib decompression worker
#'   threads per file handle. Default \code{2}. Use \code{0} to disable worker
#'   threads.
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
  decompression_threads = 2,
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
  if (!is.null(decompression_threads)) {
    params$decompression_threads <- sprintf(
      "%d",
      .validate_nonnegative_integer_param(
        decompression_threads,
        "decompression_threads"
      )
    )
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

#' Native Fixed-Width BAM/CRAM Bin Counts
#'
#' Count read starts into fixed-width genomic bins with optional duplicate
#' handling and optional per-bin GC and MAPQ summary statistics.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param path Path to the input BAM or CRAM file
#' @param bin_width Positive fixed bin width in bases
#' @param chrom Optional chromosome name filter
#' @param include_unmapped Logical. If `TRUE`, append one synthetic row for
#'   unmapped/no-coordinate records with `chrom = "*"`, and `start`, `end`, and
#'   `bin_id` set to `NA`.
#' @param reference Optional reference FASTA path for CRAM input when required,
#'   and for reference-GC output when `stats` includes `"gc"`
#' @param index_path Optional explicit BAM/CRAM index path
#' @param mapq Minimum mapping quality threshold applied after duplicate logic
#' @param require_flags Required SAM flag mask
#' @param exclude_flags Excluded SAM flag mask
#' @param rmdup Duplicate handling mode: `"none"`, `"flag"`, or `"streaming"`
#' @param stats Optional comma-separated subset of `"gc"` and `"mq"`
#'
#' @return A data frame with one row per fixed-width bin across the selected
#'   contig span, including zero-count bins, plus total, forward, reverse, and
#'   optional GC/MAPQ summary columns
#'
#' @export
rduckhts_bam_bin_counts <- function(
  con,
  path,
  bin_width,
  chrom = NULL,
  include_unmapped = FALSE,
  reference = NULL,
  index_path = NULL,
  mapq = 0,
  require_flags = 0,
  exclude_flags = 0,
  rmdup = "none",
  stats = NULL
) {
  params <- list(
    mapq = mapq,
    require_flags = require_flags,
    exclude_flags = exclude_flags,
    rmdup = sql_quote_string(rmdup)
  )
  if (!is.null(chrom)) params$chrom <- sql_quote_string(chrom)
  if (isTRUE(include_unmapped)) params$include_unmapped <- "true"
  if (!is.null(reference)) params$reference <- sql_quote_string(reference)
  if (!is.null(index_path)) params$index_path <- sql_quote_string(index_path)
  if (!is.null(stats)) params$stats <- sql_quote_string(stats)
  param_str <- build_param_str(params)
  query <- sprintf(
    "SELECT * FROM bam_bin_counts(%s, %s%s)",
    sql_quote_string(path),
    as.character(bin_width),
    param_str
  )
  DBI::dbGetQuery(con, query)
}

#' Native BAM/CRAM BED Regional Coverage Summary
#'
#' Computes samtools coverage-like regional summaries for BAM or CRAM input over
#' a BED target set, with DuckHTS-specific pre/post-filter and strand-aware
#' post-filter outputs.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param path Path to the input BAM or CRAM file
#' @param bed_path Path to the input BED file
#' @param reference Optional reference FASTA path for CRAM input when required
#' @param index_path Optional explicit BAM/CRAM index path
#' @param bed_index_path Optional explicit BED index path (reserved for future use)
#' @param mapq Minimum mapping quality threshold for post-filter summaries
#' @param min_baseq Minimum base quality threshold for post-filter base-level summaries
#' @param min_read_len Minimum read length threshold for post-filter summaries
#' @param require_flags Required SAM flag mask
#' @param exclude_flags Excluded SAM flag mask. Defaults to samtools coverage's
#'   `UNMAP|SECONDARY|QCFAIL|DUP` mask.
#' @param min_depth Minimum depth threshold for covered-base and mean-depth summaries
#' @param max_depth Maximum per-position depth cap. Set `0` to remove the cap.
#' @param decompression_threads Integer. Number of htslib decompression worker
#'   threads to use for BAM/CRAM input. `0` disables htslib worker threads.
#' @param fragment_mode Logical. Reserved for future fragment-level semantics.
#' @param strand_outputs Logical. Emit forward/reverse post-filter summary columns.
#' @param processing_threads Reserved for future parallel interval processing.
#'
#' @return A data frame with one row per BED interval and pre/post regional summaries
#'
#' @export
rduckhts_bam_bed_coverage <- function(
  con,
  path,
  bed_path,
  reference = NULL,
  index_path = NULL,
  bed_index_path = NULL,
  mapq = 0,
  min_baseq = 0,
  min_read_len = 0,
  require_flags = 0,
  exclude_flags = 1796,
  min_depth = 1,
  max_depth = 1000000,
  decompression_threads = 0,
  fragment_mode = FALSE,
  strand_outputs = TRUE,
  processing_threads = 0
) {
  params <- list(
    mapq = mapq,
    min_baseq = min_baseq,
    min_read_len = min_read_len,
    require_flags = require_flags,
    exclude_flags = exclude_flags,
    min_depth = min_depth,
    max_depth = max_depth,
    decompression_threads = .validate_nonnegative_integer_param(
      decompression_threads,
      "decompression_threads"
    ),
    fragment_mode = if (isTRUE(fragment_mode)) "true" else "false",
    strand_outputs = if (isTRUE(strand_outputs)) "true" else "false",
    processing_threads = processing_threads
  )
  if (!is.null(reference)) params$reference <- sql_quote_string(reference)
  if (!is.null(index_path)) params$index_path <- sql_quote_string(index_path)
  if (!is.null(bed_index_path)) params$bed_index_path <- sql_quote_string(bed_index_path)
  param_str <- build_param_str(params)
  query <- sprintf(
    "SELECT * FROM duckhts_bam_bed_coverage(%s, %s%s)",
    sql_quote_string(path),
    sql_quote_string(bed_path),
    param_str
  )
  DBI::dbGetQuery(con, query)
}

#' Native mosdepth-Compatible Coverage Outputs
#'
#' Writes native mosdepth-compatible coverage outputs for indexed BAM or CRAM input.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param prefix Output prefix for the mosdepth-style files
#' @param path Path to the input BAM or CRAM file
#' @param chrom Optional chromosome name filter
#' @param by Optional fixed-width window size as a string or a BED file path
#' @param fasta Optional reference FASTA path for CRAM input when required
#' @param read_groups Optional comma-separated read-group IDs, matching mosdepth's `-R`
#' @param no_per_base Skip writing `\{prefix\}.per-base.bed.gz`
#' @param threads Number of BAM decompression threads
#' @param processing_threads Number of parallel contig processing threads (0 = sequential)
#' @param flag Excluded SAM flag mask, matching mosdepth's `-F`
#' @param include_flag Required SAM flag mask, matching mosdepth's `-i`
#' @param fast_mode Logical. If `TRUE`, use mosdepth fast mode. Defaults to
#'   `FALSE`, matching upstream mosdepth.
#' @param fragment_mode Logical. If `TRUE`, count full fragment insert spans for
#'   proper pairs, matching mosdepth's `-a`. Cannot be combined with
#'   `fast_mode = TRUE`.
#' @param use_median Logical. If `TRUE`, write `by` region values as medians
#'   instead of means, matching mosdepth's `-m`.
#' @param mapq Minimum mapping quality threshold
#' @param min_frag_len Minimum absolute template length to keep, matching
#'   mosdepth's `-l`
#' @param max_frag_len Maximum absolute template length to keep, matching
#'   mosdepth's `-u`
#' @param precision_digits Number of decimal places to write in the text outputs
#' @param quantize Optional mosdepth-style quantize specification such as `":1:4:"`
#' @param thresholds Optional comma-separated coverage thresholds for `by`, matching mosdepth's `-T`
#' @param index_path Optional explicit BAM index path
#' @param overwrite Overwrite existing output files
#'
#' @return A data frame describing the written output paths
#'
#' @export
rduckhts_mosdepth <- function(
  con,
  prefix,
  path,
  chrom = NULL,
  by = NULL,
  fasta = NULL,
  read_groups = NULL,
  no_per_base = FALSE,
  threads = 2,
  processing_threads = 2,
  flag = 1796,
  include_flag = 0,
  fast_mode = FALSE,
  fragment_mode = FALSE,
  use_median = FALSE,
  mapq = 0,
  min_frag_len = -1,
  max_frag_len = -1,
  precision_digits = 2,
  quantize = NULL,
  thresholds = NULL,
  index_path = NULL,
  overwrite = FALSE
) {
  params <- list(
    no_per_base = if (isTRUE(no_per_base)) "true" else "false",
    threads = threads,
    processing_threads = processing_threads,
    flag = flag,
    include_flag = include_flag,
    fast_mode = if (isTRUE(fast_mode)) "true" else "false",
    fragment_mode = if (isTRUE(fragment_mode)) "true" else "false",
    use_median = if (isTRUE(use_median)) "true" else "false",
    mapq = mapq,
    min_frag_len = min_frag_len,
    max_frag_len = max_frag_len,
    precision_digits = precision_digits
  )
  if (!is.null(chrom)) params$chrom <- sql_quote_string(chrom)
  if (!is.null(by)) params$by <- sql_quote_string(by)
  if (!is.null(fasta)) params$fasta <- sql_quote_string(fasta)
  if (!is.null(read_groups)) params$read_groups <- sql_quote_string(read_groups)
  if (!is.null(quantize)) params$quantize <- sql_quote_string(quantize)
  if (!is.null(thresholds)) params$thresholds <- sql_quote_string(thresholds)
  if (!is.null(index_path)) params$index_path <- sql_quote_string(index_path)
  if (isTRUE(overwrite)) params$overwrite <- "true"
  param_str <- build_param_str(params)
  query <- sprintf(
    "SELECT * FROM duckhts_mosdepth(%s, %s%s)",
    sql_quote_string(prefix),
    sql_quote_string(path),
    param_str
  )
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

#' samtools idxstats-Compatible Alignment Summary
#'
#' Writes samtools idxstats-compatible alignment summary output for BAM, CRAM,
#' or SAM input.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param path Path to the input alignment file
#' @param output Optional output path for the written idxstats text file
#' @param index_path Optional explicit BAM/CRAM index path
#' @param threads htslib decompression thread count for scan fallback
#' @param overwrite Overwrite an existing output file
#'
#' @return A data frame with `success`, `path`, `output_path`,
#'   `used_index_fast_path`, and `error_message`
#'
#' @export
rduckhts_samtools_idxstats <- function(
  con,
  path,
  output = NULL,
  index_path = NULL,
  threads = 0,
  overwrite = FALSE
) {
  params <- list(threads = threads)
  if (!is.null(output)) params$output <- sql_quote_string(output)
  if (!is.null(index_path)) params$index_path <- sql_quote_string(index_path)
  if (isTRUE(overwrite)) params$overwrite <- "true"
  param_str <- build_param_str(params)
  query <- sprintf(
    "SELECT * FROM duckhts_samtools_idxstats(%s%s)",
    sql_quote_string(path),
    param_str
  )
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
#' @param attributes_map Logical. If TRUE, returns raw attributes as a scalar MAP column
#' @param attributes_list Logical. If TRUE, returns attributes as MAP(VARCHAR, VARCHAR[])
#' @param attributes_pairs Logical. If TRUE, returns attributes as a LIST of key/value/index structs
#' @param strict Logical. If TRUE, enforce GFF3 structural validation while scanning
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
  attributes_list = FALSE,
  attributes_pairs = FALSE,
  strict = FALSE,
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
  if (attributes_list) {
    params$attributes_list <- "true"
  }
  if (attributes_pairs) {
    params$attributes_pairs <- "true"
  }
  if (strict) {
    params$strict <- "true"
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
#' @param attributes_map Logical. If TRUE, returns raw attributes as a scalar MAP column
#' @param attributes_list Logical. If TRUE, returns attributes as MAP(VARCHAR, VARCHAR[])
#' @param attributes_pairs Logical. If TRUE, returns attributes as a LIST of key/value/index structs
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
  attributes_list = FALSE,
  attributes_pairs = FALSE,
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
  if (attributes_list) {
    params$attributes_list <- "true"
  }
  if (attributes_pairs) {
    params$attributes_pairs <- "true"
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
#' polygenic scores from one genotype VCF/BCF file and one or more summary-statistics files.
#'
#' @param con A DuckDB connection with DuckHTS loaded
#' @param bcf_path Path to genotype VCF/BCF file
#' @param summary_path Path(s) to summary-statistics file(s). A character vector
#'   computes multiple TSV/SSF PRS columns in one genotype scan. Use `NULL`
#'   with `summaries_list_file` to read paths from a file.
#' @param use Optional dosage source (`"GT"`, `"DS"`, `"HDS"`, `"AP"`, `"GP"`, `"AS"`)
#' @param columns Optional summary preset (`"PLINK"`, `"PLINK2"`, `"REGENIE"`, `"SAIGE"`,
#'   `"BOLT"`, `"METAL"`, `"PGS"`, `"SSF"`, `"GWAS-SSF"`)
#' @param columns_file Optional two-column summary header mapping file
#' @param q_score_thr Optional comma-separated p-value thresholds (e.g. `"1e-8,1e-6,1e-4"`)
#' @param summaries_list_file Optional path to a file (one summary path per line)
#'   or directory of summary files, matching upstream `bcftools +score --summaries`.
#' @param log_path Optional path for a matching/audit log with loaded, matched,
#'   allele-mismatch, and duplicate-marker counts per PRS.
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
  summary_path = NULL,
  use = NULL,
  columns = "PLINK",
  columns_file = NULL,
  q_score_thr = NULL,
  summaries_list_file = NULL,
  log_path = NULL,
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
  if (is.null(summary_path) && is.null(summaries_list_file)) {
    stop("summary_path or summaries_list_file must be supplied", call. = FALSE)
  }
  if (!is.null(summary_path) && !is.character(summary_path)) {
    stop("summary_path must be a character vector or NULL", call. = FALSE)
  }
  if (!is.null(summary_path) && !length(summary_path)) {
    stop("summary_path must not be empty", call. = FALSE)
  }
  if (!is.null(summary_path) && anyNA(summary_path)) {
    stop("summary_path must not contain NA", call. = FALSE)
  }
  summary_sql <- if (is.null(summary_path)) {
    "NULL"
  } else if (length(summary_path) == 1L) {
    sql_quote_string(summary_path)
  } else {
    sprintf("[%s]", paste(sql_quote_string(summary_path), collapse = ", "))
  }
  params <- list(
    sql_quote_string(bcf_path),
    summary_sql
  )
  if (!is.null(use)) params <- c(params, sprintf("use := '%s'", use))
  if (!is.null(columns)) params <- c(params, sprintf("columns := '%s'", columns))
  if (!is.null(columns_file)) params <- c(params, sprintf("columns_file := %s", sql_quote_string(columns_file)))
  if (!is.null(q_score_thr)) params <- c(params, sprintf("q_score_thr := '%s'", q_score_thr))
  if (!is.null(summaries_list_file)) params <- c(params, sprintf("summaries_list_file := %s", sql_quote_string(summaries_list_file)))
  if (!is.null(log_path)) params <- c(params, sprintf("log_path := %s", sql_quote_string(log_path)))
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

# --------------------------------------------------------------------------
# Multi-file reading helpers (internal)
# --------------------------------------------------------------------------

.format_hts_param <- function(name, value) {
  if (identical(name, "decompression_threads")) {
    value <- .validate_nonnegative_integer_param(value, name)
    return(sprintf("%s := %d", name, value))
  }

  if (is.logical(value)) {
    return(sprintf("%s := %s", name, if (isTRUE(value)) "true" else "false"))
  }
  if (is.numeric(value)) {
    return(sprintf("%s := %s", name, value))
  }
  if (is.character(value) && length(value) > 1) {
    # LIST literal
    items <- paste(sql_quote_string(value), collapse = ", ")
    return(sprintf("%s := [%s]", name, items))
  }
  if (is.character(value) && length(value) == 1) {
    return(sprintf("%s := %s", name, sql_quote_string(value)))
  }
  stop(sprintf("Unsupported parameter type for '%s': %s", name, class(value)[1]),
       call. = FALSE)
}

.build_hts_arm <- function(reader, file, params) {
  # params is a named list of already-validated non-NULL values
  param_parts <- character(0)
  if (length(params) > 0) {
    param_parts <- vapply(names(params), function(nm) {
      .format_hts_param(nm, params[[nm]])
    }, character(1))
  }
  param_str <- if (length(param_parts) > 0) {
    paste0(", ", paste(param_parts, collapse = ", "))
  } else {
    ""
  }
  quoted_file <- sql_quote_string(file)
  sprintf("SELECT %s AS filename, t.* FROM %s(%s%s) t",
          quoted_file, reader, quoted_file, param_str)
}

.expand_hts_files <- function(con, files) {
  # Use DuckDB glob() to expand each pattern (works with local and S3 paths)
  all_files <- character(0)
  for (pattern in files) {
    sql <- sprintf("SELECT file FROM glob(%s) g(file)", sql_quote_string(pattern))
    res <- DBI::dbGetQuery(con, sql)
    if (nrow(res) == 0) {
      warning(sprintf("Pattern '%s' matched no files", pattern), call. = FALSE)
    } else {
      all_files <- c(all_files, res$file)
    }
  }
  if (length(all_files) == 0) {
    stop("No files matched any of the supplied patterns", call. = FALSE)
  }
  unique(all_files)
}

.hts_multi_read <- function(con, table_name, reader, files, uniform_params,
                            .params, overwrite) {
  # Table guard — same pattern as rduckhts_bam, rduckhts_bcf, etc.
  if (!missing(table_name) && !is.null(table_name)) {
    if (DBI::dbExistsTable(con, table_name) && !overwrite) {
      stop(
        "Table '", table_name,
        "' already exists. Use overwrite = TRUE to replace it."
      )
    }
    if (DBI::dbExistsTable(con, table_name)) {
      DBI::dbRemoveTable(con, table_name)
    }
  }

  # Validate .params
  if (!is.null(.params)) {
    if (!is.data.frame(.params)) {
      stop(".params must be a data.frame or NULL", call. = FALSE)
    }
    if (!"file" %in% names(.params)) {
      stop(".params must contain a 'file' column", call. = FALSE)
    }
  }

  if (!is.null(.params)) {
    # Per-file mode: each row of .params specifies a file and optional overrides
    expanded <- character(0)
    row_map <- list()
    for (i in seq_len(nrow(.params))) {
      pat <- .params$file[i]
      row_files <- .expand_hts_files(con, pat)
      for (f in row_files) {
        expanded <- c(expanded, f)
        row_map[[f]] <- i
      }
    }
    if (length(expanded) == 0) {
      stop("No files matched any patterns in .params$file", call. = FALSE)
    }
    override_cols <- setdiff(names(.params), "file")
    arms <- vapply(expanded, function(f) {
      row_idx <- row_map[[f]]
      merged <- uniform_params
      for (col in override_cols) {
        val <- .params[[col]][row_idx]
        if (!is.na(val) && !is.null(val)) {
          merged[[col]] <- val
        }
      }
      .build_hts_arm(reader, f, merged)
    }, character(1), USE.NAMES = FALSE)
  } else {
    # Uniform mode: expand all globs and apply same params
    expanded <- .expand_hts_files(con, files)
    arms <- vapply(expanded, function(f) {
      .build_hts_arm(reader, f, uniform_params)
    }, character(1), USE.NAMES = FALSE)
  }

  union_sql <- paste(arms, collapse = " UNION ALL BY NAME ")
  create_sql <- sprintf("CREATE TABLE %s AS %s", table_name, union_sql)
  DBI::dbExecute(con, create_sql)
  invisible(TRUE)
}

# --------------------------------------------------------------------------
# Multi-file reading wrappers (exported)
# --------------------------------------------------------------------------

#' Read multiple BAM/SAM files into a DuckDB table
#'
#' Read and combine multiple BAM/SAM files via UNION ALL BY NAME,
#' materialising the result as a DuckDB table.
#' Each row includes a \code{filename} column identifying its source file.
#'
#' @param con A DBI connection to DuckDB with the duckhts extension loaded.
#' @param table_name Name of the DuckDB table to create.
#' @param files Character vector of file paths or glob patterns.
#' @param region Optional region string (e.g. \code{"chr1:1-1000"}).
#' @param index_path Optional index file path.
#' @param reference Optional reference FASTA path (for CRAM).
#' @param standard_tags Logical; include standard SAM tag columns.
#' @param auxiliary_tags Logical; include auxiliary tag map column.
#' @param sequence_encoding Optional sequence encoding (e.g. \code{"twoBit"}).
#' @param quality_representation Optional quality representation.
#' @param decompression_threads Integer. Number of htslib decompression worker
#'   threads per file handle. Default \code{2}. Use \code{0} to disable worker
#'   threads.
#' @param .params Optional data.frame with per-file parameter overrides.
#'   Must contain a \code{file} column; other columns override uniform parameters.
#'   \code{NA} values use the uniform default.
#' @param overwrite Logical; if \code{TRUE}, replace an existing table.
#' @return Invisible \code{TRUE} on success.
#' @export
rduckhts_bam_multi <- function(con, table_name, files, region = NULL,
                               index_path = NULL, reference = NULL,
                               standard_tags = FALSE, auxiliary_tags = FALSE,
                               sequence_encoding = NULL,
                               quality_representation = NULL,
                               decompression_threads = 2,
                               .params = NULL, overwrite = FALSE) {
  params <- list()
  if (!is.null(region)) params$region <- region
  if (!is.null(index_path)) params$index_path <- index_path
  if (!is.null(reference)) params$reference <- reference
  params$standard_tags <- standard_tags
  params$auxiliary_tags <- auxiliary_tags
  if (!is.null(sequence_encoding)) params$sequence_encoding <- sequence_encoding
  if (!is.null(quality_representation)) params$quality_representation <- quality_representation
  if (!is.null(decompression_threads)) params$decompression_threads <- decompression_threads
  .hts_multi_read(con, table_name, "read_bam", files, params, .params, overwrite)
}

#' Read multiple VCF/BCF files into a DuckDB table
#'
#' Read and combine multiple VCF/BCF files via UNION ALL BY NAME,
#' materialising the result as a DuckDB table.
#' Each row includes a \code{filename} column identifying its source file.
#'
#' @param con A DBI connection to DuckDB with the duckhts extension loaded.
#' @param table_name Name of the DuckDB table to create.
#' @param files Character vector of file paths or glob patterns.
#' @param region Optional region string.
#' @param index_path Optional index file path.
#' @param tidy_format Logical; use tidy FORMAT column output.
#' @param additional_csq_column_types Optional CSQ type override string.
#' @param .params Optional data.frame with per-file parameter overrides.
#' @param overwrite Logical; if \code{TRUE}, replace an existing table.
#' @return Invisible \code{TRUE} on success.
#' @export
rduckhts_bcf_multi <- function(con, table_name, files, region = NULL,
                               index_path = NULL, tidy_format = FALSE,
                               additional_csq_column_types = NULL,
                               .params = NULL, overwrite = FALSE) {
  params <- list()
  if (!is.null(region)) params$region <- region
  if (!is.null(index_path)) params$index_path <- index_path
  if (isTRUE(tidy_format)) params$tidy_format <- TRUE
  if (!is.null(additional_csq_column_types)) {
    params$additional_csq_column_types <- additional_csq_column_types
  }
  .hts_multi_read(con, table_name, "read_bcf", files, params, .params, overwrite)
}

#' Read multiple FASTQ files into a DuckDB table
#'
#' Read and combine multiple FASTQ files via UNION ALL BY NAME,
#' materialising the result as a DuckDB table.
#' Each row includes a \code{filename} column identifying its source file.
#'
#' @param con A DBI connection to DuckDB with the duckhts extension loaded.
#' @param table_name Name of the DuckDB table to create.
#' @param files Character vector of file paths or glob patterns.
#' @param mate_path Optional mate file path (for paired-end).
#' @param interleaved Logical; TRUE if file contains interleaved paired reads.
#' @param sequence_encoding Optional sequence encoding.
#' @param quality_representation Optional quality representation.
#' @param input_quality_encoding Optional input quality encoding override.
#' @param .params Optional data.frame with per-file parameter overrides.
#' @param overwrite Logical; if \code{TRUE}, replace an existing table.
#' @return Invisible \code{TRUE} on success.
#' @export
rduckhts_fastq_multi <- function(con, table_name, files, mate_path = NULL,
                                 interleaved = FALSE, sequence_encoding = NULL,
                                 quality_representation = NULL,
                                 input_quality_encoding = NULL,
                                 .params = NULL, overwrite = FALSE) {
  params <- list()
  if (!is.null(mate_path)) params$mate_path <- mate_path
  if (isTRUE(interleaved)) params$interleaved <- TRUE
  if (!is.null(sequence_encoding)) params$sequence_encoding <- sequence_encoding
  if (!is.null(quality_representation)) params$quality_representation <- quality_representation
  if (!is.null(input_quality_encoding)) params$input_quality_encoding <- input_quality_encoding
  .hts_multi_read(con, table_name, "read_fastq", files, params, .params, overwrite)
}

#' Read multiple FASTA files into a DuckDB table
#'
#' Read and combine multiple FASTA files via UNION ALL BY NAME,
#' materialising the result as a DuckDB table.
#' Each row includes a \code{filename} column identifying its source file.
#'
#' @param con A DBI connection to DuckDB with the duckhts extension loaded.
#' @param table_name Name of the DuckDB table to create.
#' @param files Character vector of file paths or glob patterns.
#' @param region Optional region string.
#' @param index_path Optional index file path.
#' @param sequence_encoding Optional sequence encoding.
#' @param .params Optional data.frame with per-file parameter overrides.
#' @param overwrite Logical; if \code{TRUE}, replace an existing table.
#' @return Invisible \code{TRUE} on success.
#' @export
rduckhts_fasta_multi <- function(con, table_name, files, region = NULL,
                                 index_path = NULL, sequence_encoding = NULL,
                                 .params = NULL, overwrite = FALSE) {
  params <- list()
  if (!is.null(region)) params$region <- region
  if (!is.null(index_path)) params$index_path <- index_path
  if (!is.null(sequence_encoding)) params$sequence_encoding <- sequence_encoding
  .hts_multi_read(con, table_name, "read_fasta", files, params, .params, overwrite)
}

#' Read multiple BED files into a DuckDB table
#'
#' Read and combine multiple BED files via UNION ALL BY NAME,
#' materialising the result as a DuckDB table.
#' Each row includes a \code{filename} column identifying its source file.
#'
#' @param con A DBI connection to DuckDB with the duckhts extension loaded.
#' @param table_name Name of the DuckDB table to create.
#' @param files Character vector of file paths or glob patterns.
#' @param region Optional region string.
#' @param index_path Optional index file path.
#' @param .params Optional data.frame with per-file parameter overrides.
#' @param overwrite Logical; if \code{TRUE}, replace an existing table.
#' @return Invisible \code{TRUE} on success.
#' @export
rduckhts_bed_multi <- function(con, table_name, files, region = NULL,
                               index_path = NULL, .params = NULL,
                               overwrite = FALSE) {
  params <- list()
  if (!is.null(region)) params$region <- region
  if (!is.null(index_path)) params$index_path <- index_path
  .hts_multi_read(con, table_name, "read_bed", files, params, .params, overwrite)
}

#' Read multiple tabix-indexed files into a DuckDB table
#'
#' Read and combine multiple tabix-indexed files via UNION ALL BY NAME,
#' materialising the result as a DuckDB table.
#' Each row includes a \code{filename} column identifying its source file.
#'
#' @param con A DBI connection to DuckDB with the duckhts extension loaded.
#' @param table_name Name of the DuckDB table to create.
#' @param files Character vector of file paths or glob patterns.
#' @param region Optional region string.
#' @param index_path Optional index file path.
#' @param header Logical or NULL; whether the file has a header line.
#' @param header_names Character vector of column names.
#' @param auto_detect Logical or NULL; enable type auto-detection.
#' @param column_types Character vector of column type names.
#' @param .params Optional data.frame with per-file parameter overrides.
#' @param overwrite Logical; if \code{TRUE}, replace an existing table.
#' @return Invisible \code{TRUE} on success.
#' @export
rduckhts_tabix_multi <- function(con, table_name, files, region = NULL,
                                 index_path = NULL, header = NULL,
                                 header_names = NULL, auto_detect = NULL,
                                 column_types = NULL, .params = NULL,
                                 overwrite = FALSE) {
  params <- list()
  if (!is.null(region)) params$region <- region
  if (!is.null(index_path)) params$index_path <- index_path
  if (!is.null(header)) params$header <- isTRUE(header)
  if (!is.null(auto_detect)) params$auto_detect <- isTRUE(auto_detect)
  if (!is.null(header_names)) {
    if (!is.character(header_names)) stop("header_names must be a character vector", call. = FALSE)
    params$header_names <- header_names
  }
  if (!is.null(column_types)) {
    if (!is.character(column_types)) stop("column_types must be a character vector", call. = FALSE)
    params$column_types <- normalize_tabix_types(column_types)
  }
  .hts_multi_read(con, table_name, "read_tabix", files, params, .params, overwrite)
}

#' Read multiple GFF files into a DuckDB table
#'
#' Read and combine multiple GFF3 files via UNION ALL BY NAME,
#' materialising the result as a DuckDB table.
#' Each row includes a \code{filename} column identifying its source file.
#'
#' @param con A DBI connection to DuckDB with the duckhts extension loaded.
#' @param table_name Name of the DuckDB table to create.
#' @param files Character vector of file paths or glob patterns.
#' @param region Optional region string.
#' @param index_path Optional index file path.
#' @param header Logical or NULL; whether the file has a header line.
#' @param header_names Character vector of column names.
#' @param auto_detect Logical or NULL; enable type auto-detection.
#' @param column_types Character vector of column type names.
#' @param attributes_map Logical; return raw attributes as a scalar MAP.
#' @param attributes_list Logical; return attributes as MAP(VARCHAR, VARCHAR[]).
#' @param attributes_pairs Logical; return attributes as a LIST of key/value/index structs.
#' @param strict Logical; enforce GFF3 structural validation while scanning.
#' @param .params Optional data.frame with per-file parameter overrides.
#' @param overwrite Logical; if \code{TRUE}, replace an existing table.
#' @return Invisible \code{TRUE} on success.
#' @export
rduckhts_gff_multi <- function(con, table_name, files, region = NULL,
                               index_path = NULL, header = NULL,
                               header_names = NULL, auto_detect = NULL,
                               column_types = NULL, attributes_map = FALSE,
                               attributes_list = FALSE, attributes_pairs = FALSE,
                               strict = FALSE,
                               .params = NULL, overwrite = FALSE) {
  params <- list()
  if (!is.null(region)) params$region <- region
  if (!is.null(index_path)) params$index_path <- index_path
  if (!is.null(header)) params$header <- isTRUE(header)
  if (!is.null(auto_detect)) params$auto_detect <- isTRUE(auto_detect)
  if (!is.null(header_names)) {
    if (!is.character(header_names)) stop("header_names must be a character vector", call. = FALSE)
    params$header_names <- header_names
  }
  if (!is.null(column_types)) {
    if (!is.character(column_types)) stop("column_types must be a character vector", call. = FALSE)
    params$column_types <- normalize_tabix_types(column_types)
  }
  if (isTRUE(attributes_map)) params$attributes_map <- TRUE
  if (isTRUE(attributes_list)) params$attributes_list <- TRUE
  if (isTRUE(attributes_pairs)) params$attributes_pairs <- TRUE
  if (isTRUE(strict)) params$strict <- TRUE
  .hts_multi_read(con, table_name, "read_gff", files, params, .params, overwrite)
}

#' Read multiple GTF files into a DuckDB table
#'
#' Read and combine multiple GTF files via UNION ALL BY NAME,
#' materialising the result as a DuckDB table.
#' Each row includes a \code{filename} column identifying its source file.
#'
#' @param con A DBI connection to DuckDB with the duckhts extension loaded.
#' @param table_name Name of the DuckDB table to create.
#' @param files Character vector of file paths or glob patterns.
#' @param region Optional region string.
#' @param index_path Optional index file path.
#' @param header Logical or NULL; whether the file has a header line.
#' @param header_names Character vector of column names.
#' @param auto_detect Logical or NULL; enable type auto-detection.
#' @param column_types Character vector of column type names.
#' @param attributes_map Logical; return raw attributes as a scalar MAP.
#' @param attributes_list Logical; return attributes as MAP(VARCHAR, VARCHAR[]).
#' @param attributes_pairs Logical; return attributes as a LIST of key/value/index structs.
#' @param .params Optional data.frame with per-file parameter overrides.
#' @param overwrite Logical; if \code{TRUE}, replace an existing table.
#' @return Invisible \code{TRUE} on success.
#' @export
rduckhts_gtf_multi <- function(con, table_name, files, region = NULL,
                               index_path = NULL, header = NULL,
                               header_names = NULL, auto_detect = NULL,
                               column_types = NULL, attributes_map = FALSE,
                               attributes_list = FALSE, attributes_pairs = FALSE,
                               .params = NULL, overwrite = FALSE) {
  params <- list()
  if (!is.null(region)) params$region <- region
  if (!is.null(index_path)) params$index_path <- index_path
  if (!is.null(header)) params$header <- isTRUE(header)
  if (!is.null(auto_detect)) params$auto_detect <- isTRUE(auto_detect)
  if (!is.null(header_names)) {
    if (!is.character(header_names)) stop("header_names must be a character vector", call. = FALSE)
    params$header_names <- header_names
  }
  if (!is.null(column_types)) {
    if (!is.character(column_types)) stop("column_types must be a character vector", call. = FALSE)
    params$column_types <- normalize_tabix_types(column_types)
  }
  if (isTRUE(attributes_map)) params$attributes_map <- TRUE
  if (isTRUE(attributes_list)) params$attributes_list <- TRUE
  if (isTRUE(attributes_pairs)) params$attributes_pairs <- TRUE
  .hts_multi_read(con, table_name, "read_gtf", files, params, .params, overwrite)
}
