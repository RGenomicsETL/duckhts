# Load DuckHTS Extension

Loads the DuckHTS extension into a DuckDB connection. This must be
called before using any of the HTS reader functions.

## Usage

``` r
rduckhts_load(con, extension_path = NULL)
```

## Arguments

- con:

  A DuckDB connection object

- extension_path:

  Optional path to the duckhts extension file. If NULL, will try to use
  the bundled extension.

## Value

TRUE if the extension was loaded successfully

## Details

The DuckDB connection must be created with
`allow_unsigned_extensions = "true"`.

## Examples

``` r
library(DBI)
library(duckdb)

con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
rduckhts_load(con)
#> [1] TRUE
dbDisconnect(con, shutdown = TRUE)
```
