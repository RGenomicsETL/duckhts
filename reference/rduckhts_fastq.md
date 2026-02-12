# Create FASTQ Table

Creates a DuckDB table from FASTQ files using the DuckHTS extension.

## Usage

``` r
rduckhts_fastq(
  con,
  table_name,
  path,
  mate_path = NULL,
  interleaved = FALSE,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- table_name:

  Name for the created table

- path:

  Path to the FASTQ file

- mate_path:

  Optional path to mate file for paired reads

- interleaved:

  Logical indicating if file is interleaved paired reads

- overwrite:

  Logical. If TRUE, overwrites existing table

## Value

Invisible TRUE on success
