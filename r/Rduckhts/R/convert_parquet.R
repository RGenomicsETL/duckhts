# Parquet converter wrappers. Generic SQL literal/path helpers live in duckhts.R.

run_convert_parquet_sql <- function(con, sql_function, path, output, params, overwrite) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("path must be a single non-empty character string", call. = FALSE)
  }
  if (!is.character(output) || length(output) != 1L || is.na(output) || !nzchar(output)) {
    stop("output must be a single non-empty character string", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE", call. = FALSE)
  }

  output_exists <- duckdb_path_exists(con, output)
  if (isTRUE(output_exists) && !overwrite) {
    stop("Output already exists: ", output, ". Use overwrite = TRUE to replace it.", call. = FALSE)
  }
  if (is.na(output_exists) && !overwrite) {
    stop(
      "Could not check output path with DuckDB: ", output,
      ". Use overwrite = TRUE to let DuckDB COPY handle the destination.",
      call. = FALSE
    )
  }
  params$overwrite <- if (overwrite) "true" else "false"

  query <- sprintf(
    "SELECT %s(%s, %s%s) AS copy_sql",
    sql_function,
    sql_quote_string(path),
    sql_quote_string(output),
    build_param_str(params)
  )
  copy_sql <- DBI::dbGetQuery(con, query)$copy_sql
  if (length(copy_sql) != 1L || is.na(copy_sql) || !nzchar(copy_sql)) {
    stop("failed to generate DuckHTS Parquet COPY statement", call. = FALSE)
  }
  DBI::dbExecute(con, copy_sql)
  invisible(output)
}

convert_parquet_common_params <- function(
  columns,
  where,
  compression,
  row_group_size,
  partition_by,
  include_metadata,
  header_text,
  metadata,
  metadata_json_file,
  write_format_version
) {
  compression <- match.arg(
    tolower(compression),
    c("zstd", "uncompressed", "none", "snappy", "gzip", "brotli", "lz4", "lz4_raw")
  )
  if (identical(compression, "none")) compression <- "uncompressed"

  row_group_size <- .validate_nonnegative_integer_param(row_group_size, "row_group_size")
  if (row_group_size == 0L) {
    stop("row_group_size must be a positive whole number", call. = FALSE)
  }
  if (!is.logical(include_metadata) || length(include_metadata) != 1L || is.na(include_metadata)) {
    stop("include_metadata must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.null(header_text)) {
    if (!is.character(header_text) || length(header_text) == 0L || anyNA(header_text)) {
      stop("header_text must be a character string or character vector without NA values", call. = FALSE)
    }
    header_text <- paste(header_text, collapse = "\n")
  }
  if (!is.null(metadata_json_file)) {
    if (!is.character(metadata_json_file) || length(metadata_json_file) != 1L || is.na(metadata_json_file) || !nzchar(metadata_json_file)) {
      stop("metadata_json_file must be a single non-empty character string", call. = FALSE)
    }
  }
  if (!is.null(where)) {
    if (!is.character(where) || length(where) != 1L || is.na(where) || !nzchar(where)) {
      stop("where must be a single non-empty character string", call. = FALSE)
    }
  }
  if (!is.character(write_format_version) || length(write_format_version) != 1L || is.na(write_format_version) || !nzchar(write_format_version)) {
    stop("write_format_version must be a single non-empty character string", call. = FALSE)
  }

  if (!is.null(metadata)) {
    if (is.null(names(metadata)) || any(!nzchar(names(metadata)))) {
      stop("metadata must be a named list or named character vector", call. = FALSE)
    }
    metadata <- as.list(metadata)
    metadata <- vapply(names(metadata), function(key) {
      value <- metadata[[key]]
      if (is.null(value) || length(value) == 0L) return("")
      if (anyNA(value)) stop("metadata value for '", key, "' must not be NA", call. = FALSE)
      paste(as.character(value), collapse = "\n")
    }, character(1), USE.NAMES = TRUE)
  }

  list(
    columns = sql_varchar_list_literal(columns, "columns"),
    where_sql = if (is.null(where)) "NULL" else sql_quote_string(where),
    compression = sql_quote_string(compression),
    row_group_size = as.character(row_group_size),
    partition_by = sql_varchar_list_literal(partition_by, "partition_by"),
    include_metadata = if (include_metadata) "true" else "false",
    header_text = if (is.null(header_text)) "NULL" else sql_quote_string(header_text),
    metadata = sql_map_literal(metadata, name = "metadata", allow_empty = TRUE),
    metadata_json_file = if (is.null(metadata_json_file)) "NULL" else sql_quote_string(metadata_json_file),
    write_format_version = sql_quote_string(write_format_version)
  )
}

#' Convert VCF/BCF reader output to Parquet with DuckHTS metadata
#'
#' `rduckhts_bcf_convert_parquet()` is a thin DBI wrapper around the extension
#' macro `duckhts_bcf_convert_parquet_sql(...)`. The extension macro builds the
#' `COPY ... TO ... (FORMAT PARQUET, KV_METADATA ...)` statement, including
#' DuckHTS write-format metadata, optional user metadata maps, selected columns,
#' optional SQL filter text, `tidy_format`, and the VCF header under the
#' `vcf_header` key. The wrapper executes the generated statement.
#'
#' @param con A DuckDB connection with DuckHTS loaded.
#' @param path Path or URI to the input VCF/BCF file.
#' @param output Path to the output Parquet file or partitioned directory.
#' @param columns Optional character vector of columns to include. Defaults to all columns.
#' @param region Optional genomic region string for indexed inputs.
#' @param index_path Optional explicit index path.
#' @param tidy_format Logical; request `read_bcf(..., tidy_format := true)`.
#' @param additional_csq_column_types Optional CSQ type override string.
#' @param decompression_threads Integer htslib decompression worker threads.
#' @param where Optional SQL predicate applied to the reader output before conversion.
#' @param compression Parquet compression, default `"zstd"`.
#' @param row_group_size Parquet row group size.
#' @param partition_by Optional character vector of output partition columns.
#' @param include_metadata Logical; include DuckHTS Parquet KV metadata.
#' @param header_text Optional corrected header text to store instead of the source header.
#' @param metadata Optional named list/vector of extra metadata. This is the
#'   primary CRAN/offline-safe path for arbitrary metadata; values with the same
#'   names as DuckHTS defaults override the default values.
#' @param metadata_json_file Optional path to a JSON file containing a top-level
#'   object of extra metadata. This requires DuckDB's `json` extension to be
#'   available when the conversion SQL is generated; otherwise DuckDB will report
#'   its normal missing-extension error. Use `metadata` for offline-safe metadata.
#' @param write_format_version DuckHTS Parquet write-format version string.
#' @param overwrite Logical; replace an existing output path. The wrapper checks
#'   existence through DuckDB `glob(...)` where possible and passes the same flag
#'   to DuckDB `COPY` for partitioned-output overwrite handling.
#'
#' @return Invisibly returns `output`.
#' @export
rduckhts_bcf_convert_parquet <- function(
  con,
  path,
  output,
  columns = NULL,
  region = NULL,
  index_path = NULL,
  tidy_format = FALSE,
  additional_csq_column_types = NULL,
  decompression_threads = 0,
  where = NULL,
  compression = "zstd",
  row_group_size = 100000L,
  partition_by = NULL,
  include_metadata = TRUE,
  header_text = NULL,
  metadata = NULL,
  metadata_json_file = NULL,
  write_format_version = "1",
  overwrite = FALSE
) {
  if (!is.logical(tidy_format) || length(tidy_format) != 1L || is.na(tidy_format)) {
    stop("tidy_format must be TRUE or FALSE", call. = FALSE)
  }
  params <- convert_parquet_common_params(
    columns, where, compression, row_group_size, partition_by, include_metadata,
    header_text, metadata, metadata_json_file, write_format_version
  )
  if (!is.null(region)) params$region <- sql_quote_string(region)
  if (!is.null(index_path)) params$index_path <- sql_quote_string(index_path)
  params$tidy_format <- if (tidy_format) "true" else "false"
  if (!is.null(additional_csq_column_types)) {
    params$additional_csq_column_types <- sql_quote_string(additional_csq_column_types)
  }
  params$decompression_threads <- sprintf(
    "%d",
    .validate_nonnegative_integer_param(decompression_threads, "decompression_threads")
  )
  run_convert_parquet_sql(con, "duckhts_bcf_convert_parquet_sql", path, output, params, overwrite)
}

#' Convert SAM/BAM/CRAM reader output to Parquet with DuckHTS metadata
#'
#' Thin DBI wrapper around extension macro `duckhts_bam_convert_parquet_sql(...)`.
#'
#' @inheritParams rduckhts_bcf_convert_parquet
#' @param path Path or URI to the input SAM/BAM/CRAM file.
#' @param reference Optional reference FASTA path for CRAM.
#' @param standard_tags Logical; include typed standard SAM tag columns.
#' @param auxiliary_tags Logical; include non-standard auxiliary tags as a map.
#' @param sequence_encoding Optional sequence encoding (`"string"` or `"nt16"`).
#' @param quality_representation Optional quality representation (`"string"` or `"phred"`).
#' @param cigar_representation Optional CIGAR representation (`"string"` or `"binary"`).
#' @param decompression_threads Integer htslib decompression worker threads.
#'
#' @return Invisibly returns `output`.
#' @export
rduckhts_bam_convert_parquet <- function(
  con,
  path,
  output,
  columns = NULL,
  region = NULL,
  index_path = NULL,
  reference = NULL,
  standard_tags = FALSE,
  auxiliary_tags = FALSE,
  sequence_encoding = NULL,
  quality_representation = NULL,
  cigar_representation = NULL,
  decompression_threads = 2,
  where = NULL,
  compression = "zstd",
  row_group_size = 100000L,
  partition_by = NULL,
  include_metadata = TRUE,
  header_text = NULL,
  metadata = NULL,
  metadata_json_file = NULL,
  write_format_version = "1",
  overwrite = FALSE
) {
  if (!is.logical(standard_tags) || length(standard_tags) != 1L || is.na(standard_tags)) {
    stop("standard_tags must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(auxiliary_tags) || length(auxiliary_tags) != 1L || is.na(auxiliary_tags)) {
    stop("auxiliary_tags must be TRUE or FALSE", call. = FALSE)
  }
  params <- convert_parquet_common_params(
    columns, where, compression, row_group_size, partition_by, include_metadata,
    header_text, metadata, metadata_json_file, write_format_version
  )
  if (!is.null(region)) params$region <- sql_quote_string(region)
  if (!is.null(index_path)) params$index_path <- sql_quote_string(index_path)
  if (!is.null(reference)) params$reference <- sql_quote_string(reference)
  params$standard_tags <- if (standard_tags) "true" else "false"
  params$auxiliary_tags <- if (auxiliary_tags) "true" else "false"
  if (!is.null(sequence_encoding)) params$sequence_encoding <- sql_quote_string(sequence_encoding)
  if (!is.null(quality_representation)) params$quality_representation <- sql_quote_string(quality_representation)
  if (!is.null(cigar_representation)) params$cigar_representation <- sql_quote_string(cigar_representation)
  params$decompression_threads <- sprintf(
    "%d",
    .validate_nonnegative_integer_param(decompression_threads, "decompression_threads")
  )
  run_convert_parquet_sql(con, "duckhts_bam_convert_parquet_sql", path, output, params, overwrite)
}

#' Convert GFF3 reader output to Parquet with DuckHTS metadata
#'
#' Thin DBI wrapper around extension macro `duckhts_gff_convert_parquet_sql(...)`.
#'
#' @inheritParams rduckhts_bcf_convert_parquet
#' @param path Path or URI to the input GFF3 file.
#' @param header Logical; pass `header := true/false` to `read_gff(...)`.
#' @param header_names Optional character vector of column names.
#' @param auto_detect Logical; request type auto-detection.
#' @param column_types Optional character vector of DuckDB column types.
#' @param attributes_map Logical; expose GFF attributes as `MAP(VARCHAR, VARCHAR)`.
#' @param attributes_list Logical; expose attributes as `MAP(VARCHAR, VARCHAR[])`.
#' @param attributes_pairs Logical; expose attributes as key/value/index structs.
#' @param strict Logical; enable strict GFF validation while scanning.
#'
#' @return Invisibly returns `output`.
#' @export
rduckhts_gff_convert_parquet <- function(
  con,
  path,
  output,
  columns = NULL,
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
  where = NULL,
  compression = "zstd",
  row_group_size = 100000L,
  partition_by = NULL,
  include_metadata = TRUE,
  header_text = NULL,
  metadata = NULL,
  metadata_json_file = NULL,
  write_format_version = "1",
  overwrite = FALSE
) {
  if (!is.logical(attributes_map) || length(attributes_map) != 1L || is.na(attributes_map)) {
    stop("attributes_map must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(attributes_list) || length(attributes_list) != 1L || is.na(attributes_list)) {
    stop("attributes_list must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(attributes_pairs) || length(attributes_pairs) != 1L || is.na(attributes_pairs)) {
    stop("attributes_pairs must be TRUE or FALSE", call. = FALSE)
  }
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("strict must be TRUE or FALSE", call. = FALSE)
  }
  params <- convert_parquet_common_params(
    columns, where, compression, row_group_size, partition_by, include_metadata,
    header_text, metadata, metadata_json_file, write_format_version
  )
  if (!is.null(region)) params$region <- sql_quote_string(region)
  if (!is.null(index_path)) params$index_path <- sql_quote_string(index_path)
  if (!is.null(header)) {
    if (!is.logical(header) || length(header) != 1L || is.na(header)) {
      stop("header must be TRUE or FALSE", call. = FALSE)
    }
    params$header <- if (header) "true" else "false"
  }
  params$header_names <- sql_varchar_list_literal(header_names, "header_names")
  if (!is.null(auto_detect)) {
    if (!is.logical(auto_detect) || length(auto_detect) != 1L || is.na(auto_detect)) {
      stop("auto_detect must be TRUE or FALSE", call. = FALSE)
    }
    params$auto_detect <- if (auto_detect) "true" else "false"
  }
  if (!is.null(column_types)) column_types <- normalize_tabix_types(column_types)
  params$column_types <- sql_varchar_list_literal(column_types, "column_types")
  params$attributes_map <- if (attributes_map) "true" else "false"
  params$attributes_list <- if (attributes_list) "true" else "false"
  params$attributes_pairs <- if (attributes_pairs) "true" else "false"
  params$strict <- if (strict) "true" else "false"
  run_convert_parquet_sql(con, "duckhts_gff_convert_parquet_sql", path, output, params, overwrite)
}

#' Convert generic tabix reader output to Parquet with DuckHTS metadata
#'
#' Thin DBI wrapper around extension macro `duckhts_tabix_convert_parquet_sql(...)`.
#'
#' @inheritParams rduckhts_gff_convert_parquet
#' @param path Path or URI to the input tabix-indexed text file.
#' @param header Logical; pass `header := true/false` to `read_tabix(...)`.
#'
#' @return Invisibly returns `output`.
#' @export
rduckhts_tabix_convert_parquet <- function(
  con,
  path,
  output,
  columns = NULL,
  region = NULL,
  index_path = NULL,
  header = NULL,
  header_names = NULL,
  auto_detect = NULL,
  column_types = NULL,
  where = NULL,
  compression = "zstd",
  row_group_size = 100000L,
  partition_by = NULL,
  include_metadata = TRUE,
  header_text = NULL,
  metadata = NULL,
  metadata_json_file = NULL,
  write_format_version = "1",
  overwrite = FALSE
) {
  params <- convert_parquet_common_params(
    columns, where, compression, row_group_size, partition_by, include_metadata,
    header_text, metadata, metadata_json_file, write_format_version
  )
  if (!is.null(region)) params$region <- sql_quote_string(region)
  if (!is.null(index_path)) params$index_path <- sql_quote_string(index_path)
  if (!is.null(header)) {
    if (!is.logical(header) || length(header) != 1L || is.na(header)) {
      stop("header must be TRUE or FALSE", call. = FALSE)
    }
    params$header <- if (header) "true" else "false"
  }
  params$header_names <- sql_varchar_list_literal(header_names, "header_names")
  if (!is.null(auto_detect)) {
    if (!is.logical(auto_detect) || length(auto_detect) != 1L || is.na(auto_detect)) {
      stop("auto_detect must be TRUE or FALSE", call. = FALSE)
    }
    params$auto_detect <- if (auto_detect) "true" else "false"
  }
  if (!is.null(column_types)) column_types <- normalize_tabix_types(column_types)
  params$column_types <- sql_varchar_list_literal(column_types, "column_types")
  run_convert_parquet_sql(con, "duckhts_tabix_convert_parquet_sql", path, output, params, overwrite)
}
