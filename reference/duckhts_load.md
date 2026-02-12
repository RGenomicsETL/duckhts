# Load the duckhts extension into a DuckDB connection

Load the duckhts extension into a DuckDB connection

## Usage

``` r
duckhts_load(con = NULL, extension_path = NULL)
```

## Arguments

- con:

  An existing DuckDB connection, or `NULL` to create one.

- extension_path:

  Explicit path to the `.duckdb_extension` file. If `NULL`, uses the
  default location in the installed package.

## Value

The DuckDB connection (invisibly).
