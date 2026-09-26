# BGZF Decompress a File

Decompresses a BGZF file using the DuckHTS extension.

## Usage

``` r
rduckhts_bgunzip(
  con,
  path,
  output_path = NULL,
  threads = 4,
  keep = TRUE,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to the BGZF-compressed input file

- output_path:

  Optional explicit output path

- threads:

  BGZF worker thread count

- keep:

  Keep the compressed input file after decompression

- overwrite:

  Overwrite an existing output file

## Value

A data frame describing the created output file
