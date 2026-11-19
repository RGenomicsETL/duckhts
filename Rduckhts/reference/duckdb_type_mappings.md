# DuckDB to R Type Mappings

Returns named character vectors that map DuckDB column types to the R
storage type that `duckdb` returns, and R type names to the DuckDB types
used for tabix column types (see
[`normalize_tabix_types`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/normalize_tabix_types.md)).

## Usage

``` r
duckdb_type_mappings()
```

## Value

A named list with two elements:

- duckdb_to_r:

  Named character vector mapping DuckDB types to R types

- r_to_duckdb:

  Named character vector mapping R types to DuckDB types

## Details

`duckdb_to_r` gives the [`typeof()`](https://rdrr.io/r/base/typeof.html)
of the returned column, with one exception: `"MAP"` maps to
`"data.frame"`, the class of each cell. The class is not part of the
mapping:

- DATE, TIMESTAMP and TIME columns have storage type double and class
  `Date`, `POSIXct` and `difftime` (seconds).

- BIGINT and UBIGINT columns are double under the default
  `bigint = "numeric"` connection setting.

- A LIST column such as `VARCHAR[]` is a list with one vector per row. A
  MAP column is a list with one data frame per row, with columns `key`
  and `value`.

## Examples

``` r
mappings <- duckdb_type_mappings()
mappings$duckdb_to_r["BIGINT"]
#>   BIGINT 
#> "double" 
mappings$r_to_duckdb["integer"]
#>   integer 
#> "INTEGER" 
```
