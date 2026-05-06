# Create GTF Table

Creates a DuckDB table from GTF files using the DuckHTS extension.

## Usage

``` r
rduckhts_gtf(
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
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- table_name:

  Name for the created table

- path:

  Path to the GTF file

- region:

  Optional genomic region (e.g., "chr1:1000-2000")

- index_path:

  Optional explicit path to index file (.tbi/.csi)

- header:

  Logical. If TRUE, use first non-meta line as column names

- header_names:

  Character vector to override column names

- auto_detect:

  Logical. If TRUE, infer basic numeric column types

- column_types:

  Character vector of column types (e.g. "BIGINT", "VARCHAR")

- attributes_map:

  Logical. If TRUE, returns raw attributes as a scalar MAP column

- attributes_list:

  Logical. If TRUE, returns attributes as MAP(VARCHAR, VARCHAR\[\])

- attributes_pairs:

  Logical. If TRUE, returns attributes as a LIST of key/value/index
  structs

- overwrite:

  Logical. If TRUE, overwrites existing table

## Value

Invisible TRUE on success
