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
#> Error in duckdb_result(connection = conn, stmt_lst = stmt_lst, arrow = arrow): Invalid Error: IO Error: Extension "/tmp/Rtmph0wBxf/temp_libpath1b84b3ec7f5/Rduckhts/duckhts_extension/build/duckhts.duckdb_extension" could not be loaded: libhts.so.3: cannot open shared object file: No such file or directory
#> ℹ Context: rapi_execute
#> ℹ Error type: INVALID
dbDisconnect(con, shutdown = TRUE)
```
