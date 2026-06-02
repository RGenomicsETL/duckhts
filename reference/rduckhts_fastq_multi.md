# Read multiple FASTQ files into a DuckDB table

Read and combine multiple FASTQ files via UNION ALL BY NAME,
materialising the result as a DuckDB table. Each row includes a
`filename` column identifying its source file.

## Usage

``` r
rduckhts_fastq_multi(
  con,
  table_name,
  files,
  mate_path = NULL,
  interleaved = FALSE,
  sequence_encoding = NULL,
  quality_representation = NULL,
  input_quality_encoding = NULL,
  scan_mode = NULL,
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

- mate_path:

  Optional mate file path (for paired-end).

- interleaved:

  Logical; TRUE if file contains interleaved paired reads.

- sequence_encoding:

  Optional sequence encoding.

- quality_representation:

  Optional quality representation.

- input_quality_encoding:

  Optional input quality encoding override.

- scan_mode:

  Optional scan mode (`"auto"` or `"sequential"`).

- .params:

  Optional data.frame with per-file parameter overrides.

- overwrite:

  Logical; if `TRUE`, replace an existing table.

## Value

Invisible `TRUE` on success.
