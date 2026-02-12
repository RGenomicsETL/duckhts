# DuckDB to R Type Mappings

The mapping covers the most common data types used in HTS file
processing:

- BIGINT ↔ double (not integer due to 64-bit overflow protection)

- DOUBLE ↔ numeric/double

- VARCHAR ↔ character/string

- BOOLEAN ↔ logical

- ARRAY types (e.g., VARCHAR\[\], BIGINT\[\]) ↔ list

- MAP types (e.g., MAP(VARCHAR, VARCHAR)) ↔ data.frame

Important notes:

- 64-bit integers (BIGINT, UBIGINT) become double to prevent overflow

- DATE/TIME values return as Unix epoch numbers (double)

- MAP types become data frames with 'key' and 'value' columns

- ARRAY types become vectors (which are lists in R terminology)

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

Returns a named list mapping between DuckDB and R data types. This is
useful for understanding type conversions when reading HTS files or when
specifying column types in tabix functions.

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
