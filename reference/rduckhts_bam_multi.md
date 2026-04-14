# Read multiple BAM/SAM files into a DuckDB table

Read and combine multiple BAM/SAM files via UNION ALL BY NAME,
materialising the result as a DuckDB table. Each row includes a
`filename` column identifying its source file.

## Usage

``` r
rduckhts_bam_multi(
  con,
  table_name,
  files,
  region = NULL,
  index_path = NULL,
  reference = NULL,
  standard_tags = FALSE,
  auxiliary_tags = FALSE,
  sequence_encoding = NULL,
  quality_representation = NULL,
  decompression_threads = 2,
  .params = NULL,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DBI connection to DuckDB with the duckhts extension loaded.

- table_name:

  Name of the DuckDB table to create.

- files:

  Character vector of file paths or glob patterns.

- region:

  Optional region string (e.g. `"chr1:1-1000"`).

- index_path:

  Optional index file path.

- reference:

  Optional reference FASTA path (for CRAM).

- standard_tags:

  Logical; include standard SAM tag columns.

- auxiliary_tags:

  Logical; include auxiliary tag map column.

- sequence_encoding:

  Optional sequence encoding (e.g. `"twoBit"`).

- quality_representation:

  Optional quality representation.

- decompression_threads:

  Integer. Number of htslib decompression worker threads per file
  handle. Default `2`. Use `0` to disable worker threads.

- .params:

  Optional data.frame with per-file parameter overrides. Must contain a
  `file` column; other columns override uniform parameters. `NA` values
  use the uniform default.

- overwrite:

  Logical; if `TRUE`, replace an existing table.

## Value

Invisible `TRUE` on success.
