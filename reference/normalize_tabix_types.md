# Normalize R Data Types to DuckDB Types for Tabix

Normalizes R data type names to their corresponding DuckDB types for use
with tabix readers. This function handles common R type name variations
and maps them to appropriate DuckDB column types.

## Usage

``` r
normalize_tabix_types(types)
```

## Arguments

- types:

  A character vector of R data type names to be normalized.

## Value

A character vector of normalized DuckDB type names suitable for tabix
columns.

## Details

The function performs the following normalizations:

- Integer types (integer, int, int32, int64) -\> BIGINT

- Numeric types (numeric, double, float) -\> DOUBLE

- Character types (character, string, chr) -\> VARCHAR

- Logical types (logical, bool, boolean) -\> BOOLEAN

- Other types -\> Converted to uppercase as-is

If an empty vector is provided, it returns the empty vector unchanged.

## See also

[`rduckhts_tabix`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_tabix.md)
for using normalized types with tabix readers,
[`duckdb_type_mappings`](https://rgenomicsetl.github.io/duckhts/reference/duckdb_type_mappings.md)
for the complete type mapping table.

## Examples

``` r
# Normalize mixed type names
normalize_tabix_types(c("integer", "character", "numeric"))
#> [1] "BIGINT"  "VARCHAR" "DOUBLE" 
# Returns: c("BIGINT", "VARCHAR", "DOUBLE")

# Handle variations
normalize_tabix_types(c("int", "string", "float"))
#> [1] "BIGINT"  "VARCHAR" "DOUBLE" 
# Returns: c("BIGINT", "VARCHAR", "DOUBLE")
```
