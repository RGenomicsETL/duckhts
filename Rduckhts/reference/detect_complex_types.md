# Detect Complex Types in DuckDB Table

Lists the columns of a table whose DuckDB type contains `[` or `MAP`.
`duckdb` returns such columns as R lists.

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

A data frame with columns `column_name`, `column_type`, `r_type`
(`"vector"` or `"data.frame"`, the class of each cell) and
`description`; a zero-column data frame when the table has no such
columns.

## Examples

``` r
library(DBI)
library(duckdb)

con <- rduckhts_connect()
bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
rduckhts_bcf(con, "variants", bcf_path, overwrite = TRUE)
complex_cols <- detect_complex_types(con, "variants")
print(complex_cols)
#>    column_name column_type r_type                   description
#> 5          ALT   VARCHAR[] vector ARRAY type - will be R vector
#> 7       FILTER   VARCHAR[] vector ARRAY type - will be R vector
#> 9     INFO_DP4   INTEGER[] vector ARRAY type - will be R vector
#> 10     INFO_AC   INTEGER[] vector ARRAY type - will be R vector
#> 14 FORMAT_TT_A   INTEGER[] vector ARRAY type - will be R vector
#> 18 FORMAT_GL_A     FLOAT[] vector ARRAY type - will be R vector
#> 19 FORMAT_TT_B   INTEGER[] vector ARRAY type - will be R vector
#> 23 FORMAT_GL_B     FLOAT[] vector ARRAY type - will be R vector
dbDisconnect(con, shutdown = TRUE)
```
