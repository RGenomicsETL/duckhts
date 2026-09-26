# BGZF Compress a File

Compresses a plain file to BGZF using the DuckHTS extension.

## Usage

``` r
rduckhts_bgzip(
  con,
  path,
  output_path = NULL,
  threads = 4,
  level = -1,
  keep = TRUE,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to the input file

- output_path:

  Optional explicit output path

- threads:

  BGZF worker thread count

- level:

  Compression level, or -1 for the htslib default

- keep:

  Keep the original input file after compression

- overwrite:

  Overwrite an existing output file

## Value

A data frame describing the created BGZF file
