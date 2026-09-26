# Build BAM or CRAM Index

Builds a BAM or CRAM index using the DuckHTS extension.

## Usage

``` r
rduckhts_bam_index(con, path, index_path = NULL, min_shift = 0, threads = 4)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to the input BAM or CRAM file

- index_path:

  Optional explicit output path for the created index

- min_shift:

  Index format selector used by htslib

- threads:

  htslib indexing thread count

## Value

A data frame with \`success\`, \`index_path\`, and \`index_format\`
