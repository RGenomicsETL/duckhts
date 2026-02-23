# Build FASTA Index

Builds a FASTA index (.fai) using the DuckHTS extension.

## Usage

``` r
rduckhts_fasta_index(con, path, index_path = NULL)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to the FASTA file

- index_path:

  Optional explicit output path for FASTA index file (.fai)

## Value

A data frame with columns \`success\` and \`index_path\`
