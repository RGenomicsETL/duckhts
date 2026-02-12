# Create GTF Table

Creates a DuckDB table from GTF files using the DuckHTS extension.

## Usage

``` r
rduckhts_gtf(
  con,
  table_name,
  path,
  region = NULL,
  attributes_map = FALSE,
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

- attributes_map:

  Logical. If TRUE, returns attributes as a MAP column

- overwrite:

  Logical. If TRUE, overwrites existing table

## Value

Invisible TRUE on success
