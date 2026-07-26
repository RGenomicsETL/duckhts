# Load DuckHTS Extension

Loads the DuckHTS extension into an existing DuckDB connection. This
must be called before using HTS reader functions on a connection not
created by
[`rduckhts_connect()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_connect.md).

## Usage

``` r
rduckhts_load(con, extension_path = NULL)
```

## Arguments

- con:

  A DuckDB connection object.

- extension_path:

  Optional path to the DuckHTS extension file. If `NULL`, uses the
  extension bundled with Rduckhts.

## Value

`TRUE` if the extension was loaded successfully.

## Details

The connection must permit unsigned extension loading. With current
versions of the duckdb R package, its driver must also permit extension
loading. Prefer
[`rduckhts_connect()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_connect.md)
when Rduckhts owns the connection.

## Examples

``` r
con <- rduckhts_connect()
DBI::dbDisconnect(con, shutdown = TRUE)
```
