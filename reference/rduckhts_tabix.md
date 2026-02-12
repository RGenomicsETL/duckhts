# Create Tabix-Indexed File Table

Creates a DuckDB table from any tabix-indexed file using the DuckHTS
extension.

## Usage

``` r
rduckhts_tabix(
  con,
  table_name,
  path,
  region = NULL,
  header = NULL,
  header_names = NULL,
  auto_detect = NULL,
  column_types = NULL,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- table_name:

  Name for the created table

- path:

  Path to the tabix-indexed file

- region:

  Optional genomic region (e.g., "chr1:1000-2000")

- header:

  Logical. If TRUE, use first non-meta line as column names

- header_names:

  Character vector to override column names

- auto_detect:

  Logical. If TRUE, infer basic numeric column types

- column_types:

  Character vector of column types (e.g. "BIGINT", "VARCHAR")

- overwrite:

  Logical. If TRUE, overwrites existing table

## Value

Invisible TRUE on success
