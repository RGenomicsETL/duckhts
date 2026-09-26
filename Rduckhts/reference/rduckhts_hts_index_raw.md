# Read Raw HTS Index Blob

Returns raw index metadata blob data for a file index.

## Usage

``` r
rduckhts_hts_index_raw(con, path, format = NULL, index_path = NULL)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to input HTS file

- format:

  Optional format hint

- index_path:

  Optional explicit path to index file

## Value

A data frame with raw index blob metadata.
