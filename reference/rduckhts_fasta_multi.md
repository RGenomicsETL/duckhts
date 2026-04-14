# Read multiple FASTA files into a DuckDB table

Read and combine multiple FASTA files via UNION ALL BY NAME,
materialising the result as a DuckDB table. Each row includes a
`filename` column identifying its source file.

## Usage

``` r
rduckhts_fasta_multi(
  con,
  table_name,
  files,
  region = NULL,
  index_path = NULL,
  sequence_encoding = NULL,
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

- sequence_encoding:

  Optional sequence encoding.

- .params:

  Optional data.frame with per-file parameter overrides.

- overwrite:

  Logical; if `TRUE`, replace an existing table.

## Value

Invisible `TRUE` on success.
