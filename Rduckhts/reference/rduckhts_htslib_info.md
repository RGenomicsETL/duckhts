# Inspect the Loaded htslib Build

Returns the version and build features reported by the htslib library
that is actually loaded with DuckHTS.

## Usage

``` r
rduckhts_htslib_info(con)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded.

## Value

A one-row data frame with \`version\`, \`feature_bits\`, and
\`feature_string\` columns.

## Examples

``` r
if (FALSE) { # \dontrun{
con <- rduckhts_connect()
rduckhts_htslib_info(con)
DBI::dbDisconnect(con, shutdown = TRUE)
} # }
```
