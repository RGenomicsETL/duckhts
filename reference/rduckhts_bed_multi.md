# Read multiple BED files into a DuckDB table

Read and combine multiple BED files via UNION ALL BY NAME, materialising
the result as a DuckDB table. Each row includes a `filename` column
identifying its source file.

## Usage

``` r
rduckhts_bed_multi(
  con,
  table_name,
  files,
  region = NULL,
  index_path = NULL,
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

- region:

  Optional region string.

- index_path:

  Optional index file path.

- scan_mode:

  Optional scan mode (`"auto"` or `"sequential"`).

- .params:

  Optional data.frame with per-file parameter overrides.

- overwrite:

  Logical; if `TRUE`, replace an existing table.

## Value

Invisible `TRUE` on success.
