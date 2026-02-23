# Create FASTA Table

Creates a DuckDB table from FASTA files using the DuckHTS extension.

## Usage

``` r
rduckhts_fasta(
  con,
  table_name,
  path,
  region = NULL,
  index_path = NULL,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- table_name:

  Name for the created table

- path:

  Path to the FASTA file

- region:

  Optional genomic region (e.g., "chr1:1000-2000" or
  "chr1:1-10,chr2:5-20")

- index_path:

  Optional explicit path to FASTA index file (.fai)

- overwrite:

  Logical. If TRUE, overwrites existing table

## Value

Invisible TRUE on success
