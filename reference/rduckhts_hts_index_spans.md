# Read HTS Index Spans

Returns index span-oriented metadata for planning range workloads.

## Usage

``` r
rduckhts_hts_index_spans(con, path, format = NULL, index_path = NULL)
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

A data frame with span-oriented index metadata.
