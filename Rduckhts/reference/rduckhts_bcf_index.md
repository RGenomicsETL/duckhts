# Build VCF or BCF Index

Builds a TBI or CSI index for a VCF/BCF file using the DuckHTS
extension.

## Usage

``` r
rduckhts_bcf_index(con, path, index_path = NULL, min_shift = NULL, threads = 4)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to the input VCF/BCF file

- index_path:

  Optional explicit output path for the created index

- min_shift:

  Optional explicit min_shift passed to htslib

- threads:

  htslib indexing thread count

## Value

A data frame with \`success\`, \`index_path\`, and \`index_format\`
