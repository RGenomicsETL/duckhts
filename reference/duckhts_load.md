# Load the duckhts extension into a DuckDB connection

Compatibility wrapper around
[`rduckhts_connect()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_connect.md)
and
[`rduckhts_load()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_load.md).
Prefer the `rduckhts_*` functions in new code.

## Usage

``` r
duckhts_load(con = NULL, extension_path = NULL)
```

## Arguments

- con:

  An existing DuckDB connection, or `NULL` to create a package-owned
  connection with
  [`rduckhts_connect()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_connect.md).

- extension_path:

  Explicit path to the `.duckdb_extension` file. If `NULL`, uses the
  default location in the installed package.

## Value

The DuckDB connection (invisibly).
