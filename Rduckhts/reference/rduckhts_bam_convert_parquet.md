# Convert SAM/BAM/CRAM reader output to Parquet with DuckHTS metadata

Thin DBI wrapper around extension macro
\`duckhts_bam_convert_parquet_sql(...)\`.

## Usage

``` r
rduckhts_bam_convert_parquet(
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
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded.

- path:

  Path or URI to the input SAM/BAM/CRAM file.

- output:

  Path to the output Parquet file or partitioned directory.

- columns:

  Optional character vector of columns to include. Defaults to all
  columns.

- region:

  Optional genomic region string for indexed inputs.

- index_path:

  Optional explicit index path.

- reference:

  Optional reference FASTA path for CRAM.

- standard_tags:

  Logical; include typed standard SAM tag columns.

- auxiliary_tags:

  Logical; include non-standard auxiliary tags as a map.

- sequence_encoding:

  Optional sequence encoding (\`"string"\` or \`"nt16"\`).

- quality_representation:

  Optional quality representation (\`"string"\` or \`"phred"\`).

- cigar_representation:

  Optional CIGAR representation (\`"string"\` or \`"binary"\`).

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
