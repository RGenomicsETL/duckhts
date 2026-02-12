# Create FASTA Table

Creates a DuckDB table from FASTA files using the DuckHTS extension.

## Usage

``` r
rduckhts_fasta(con, table_name, path, overwrite = FALSE)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- table_name:

  Name for the created table

- path:

  Path to the FASTA file

- overwrite:

  Logical. If TRUE, overwrites existing table

## Value

Invisible TRUE on success
