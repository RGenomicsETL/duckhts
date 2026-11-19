# Extract MAP Keys and Values

Extracts keys, values, or the value of one key from each row of a DuckDB
MAP column. `duckdb` returns a MAP column as a list of data frames with
columns `key` and `value`.

## Usage

``` r
extract_map_data(map_col, operation = "keys", default = NA)
```

## Arguments

- map_col:

  A MAP column: a list of data frames with columns `key` and `value`

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

con <- rduckhts_connect()
gff_path <- system.file("extdata", "gff_file.gff.gz", package = "Rduckhts")
rduckhts_gff(con, "annotations", gff_path, attributes_map = TRUE, overwrite = TRUE)
data <- dbGetQuery(con, "SELECT attributes FROM annotations LIMIT 5")
keys <- extract_map_data(data$attributes, "keys")
name_values <- extract_map_data(data$attributes, "Name")
dbDisconnect(con, shutdown = TRUE)
```
