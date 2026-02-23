# Detect Complex Types in DuckDB Table

Identifies columns in a DuckDB table that contain complex types (ARRAY
or MAP) that will be returned as R lists.

## Usage

``` r
detect_complex_types(con, table_name)
```

## Arguments

- con:

  A DuckDB connection

- table_name:

  Name of the table to analyze

## Value

A data frame with columns that have complex types, showing column_name,
column_type, and a description of R type.

## Examples

``` r
library(DBI)
library(duckdb)

con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
rduckhts_load(con)
#> [1] TRUE
bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
rduckhts_bcf(con, "variants", bcf_path, overwrite = TRUE)
complex_cols <- detect_complex_types(con, "variants")
print(complex_cols)
#>    column_name column_type r_type                   description
#> 5          ALT   VARCHAR[] vector ARRAY type - will be R vector
#> 7       FILTER   VARCHAR[] vector ARRAY type - will be R vector
#> 10     INFO_AC   INTEGER[] vector ARRAY type - will be R vector
#> 14 FORMAT_TT_A   INTEGER[] vector ARRAY type - will be R vector
#> 18 FORMAT_GL_A     FLOAT[] vector ARRAY type - will be R vector
#> 19 FORMAT_TT_B   INTEGER[] vector ARRAY type - will be R vector
#> 23 FORMAT_GL_B     FLOAT[] vector ARRAY type - will be R vector
dbDisconnect(con, shutdown = TRUE)
```
