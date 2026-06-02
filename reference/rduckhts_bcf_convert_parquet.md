# Convert VCF/BCF reader output to Parquet with DuckHTS metadata

\`rduckhts_bcf_convert_parquet()\` is a thin DBI wrapper around the
extension macro \`duckhts_bcf_convert_parquet_sql(...)\`. The extension
macro builds the \`COPY ... TO ... (FORMAT PARQUET, KV_METADATA ...)\`
statement, including DuckHTS write-format metadata, optional user
metadata maps, selected columns, optional SQL filter text,
\`tidy_format\`, and the VCF header under the \`vcf_header\` key. The
wrapper executes the generated statement.

## Usage

``` r
rduckhts_bcf_convert_parquet(
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
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded.

- path:

  Path or URI to the input VCF/BCF file.

- output:

  Path to the output Parquet file or partitioned directory.

- columns:

  Optional character vector of columns to include. Defaults to all
  columns.

- region:

  Optional genomic region string for indexed inputs.

- index_path:

  Optional explicit index path.

- tidy_format:

  Logical; request \`read_bcf(..., tidy_format := true)\`.

- additional_csq_column_types:

  Optional CSQ type override string.

- decompression_threads:

  Integer htslib decompression worker threads.

- where:

  Optional SQL predicate applied to the reader output before conversion.

- compression:

  Parquet compression, default \`"zstd"\`.

- row_group_size:

  Parquet row group size.

- partition_by:

  Optional character vector of output partition columns.

- include_metadata:

  Logical; include DuckHTS Parquet KV metadata.

- header_text:

  Optional corrected header text to store instead of the source header.

- metadata:

  Optional named list/vector of extra metadata. This is the primary
  CRAN/offline-safe path for arbitrary metadata; values with the same
  names as DuckHTS defaults override the default values.

- metadata_json_file:

  Optional path to a JSON file containing a top-level object of extra
  metadata. This requires DuckDB's \`json\` extension to be available
  when the conversion SQL is generated; otherwise DuckDB will report its
  normal missing-extension error. Use \`metadata\` for offline-safe
  metadata.

- write_format_version:

  DuckHTS Parquet write-format version string.

- overwrite:

  Logical; replace an existing output path. The wrapper checks existence
  through DuckDB \`glob(...)\` where possible and passes the same flag
  to DuckDB \`COPY\` for partitioned-output overwrite handling.

## Value

Invisibly returns \`output\`.
