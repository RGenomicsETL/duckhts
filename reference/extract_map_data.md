# Extract MAP Keys and Values

Helper function to work with DuckDB MAP data (returned as data frames).
Can extract keys, values, or search for specific key-value pairs.

## Usage

``` r
extract_map_data(map_col, operation = "keys", default = NA)
```

## Arguments

- map_col:

  A data frame column from DuckDB MAP data

- operation:

  What to extract: "keys", "values", or a specific key name

- default:

  Default value if key is not found (only used when operation is a key
  name)

## Value

Extracted data based on the operation

## Examples

``` r
library(DBI)
library(duckdb)

con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
rduckhts_load(con)
#> [1] TRUE
gff_path <- system.file("extdata", "gff_file.gff.gz", package = "Rduckhts")
rduckhts_gff(con, "annotations", gff_path, attributes_map = TRUE, overwrite = TRUE)
data <- dbGetQuery(con, "SELECT attributes FROM annotations LIMIT 5")
keys <- extract_map_data(data$attributes, "keys")
name_values <- extract_map_data(data$attributes, "Name")
dbDisconnect(con, shutdown = TRUE)
```
