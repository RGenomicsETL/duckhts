# Read multiple GFF files into a DuckDB table

Read and combine multiple GFF3 files via UNION ALL BY NAME,
materialising the result as a DuckDB table. Each row includes a
`filename` column identifying its source file.

## Usage

``` r
rduckhts_gff_multi(
  con,
  table_name,
  files,
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

- header:

  Logical or NULL; whether the file has a header line.

- header_names:

  Character vector of column names.

- auto_detect:

  Logical or NULL; enable type auto-detection.

- column_types:

  Character vector of column type names.

- attributes_map:

  Logical; return raw attributes as a scalar MAP.

- attributes_list:

  Logical; return attributes as MAP(VARCHAR, VARCHAR\[\]).

- attributes_pairs:

  Logical; return attributes as a LIST of key/value/index structs.

- strict:

  Logical; enforce GFF3 structural validation while scanning.

- .params:

  Optional data.frame with per-file parameter overrides.

- overwrite:

  Logical; if `TRUE`, replace an existing table.

## Value

Invisible `TRUE` on success.
