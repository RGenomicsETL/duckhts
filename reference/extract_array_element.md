# Extract Array Elements Safely

Helper function to safely extract elements from DuckDB arrays (returned
as R lists) with proper error handling.

## Usage

``` r
extract_array_element(array_col, index = NULL, default = NA)
```

## Arguments

- array_col:

  A list column from DuckDB array data

- index:

  Numeric index (1-based). If NULL, returns full list

- default:

  Default value if index is out of bounds

## Value

The array element at the specified index, or full array if index is NULL

## Examples

``` r
library(DBI)
library(duckdb)

con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
rduckhts_load(con)
#> [1] TRUE
bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
rduckhts_bcf(con, "variants", bcf_path, overwrite = TRUE)
data <- dbGetQuery(con, "SELECT ALT FROM variants LIMIT 5")
first_alt <- extract_array_element(data$ALT, 1)
all_alts <- extract_array_element(data$ALT)
dbDisconnect(con, shutdown = TRUE)
```
