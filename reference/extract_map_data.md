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
#> Error in duckdb_result(connection = conn, stmt_lst = stmt_lst, arrow = arrow): Invalid Error: IO Error: Extension "/tmp/Rtmph0wBxf/temp_libpath1b84b3ec7f5/Rduckhts/duckhts_extension/build/duckhts.duckdb_extension" could not be loaded: libhts.so.3: cannot open shared object file: No such file or directory
#> ℹ Context: rapi_execute
#> ℹ Error type: INVALID
gff_path <- system.file("extdata", "gff_file.gff.gz", package = "Rduckhts")
rduckhts_gff(con, "annotations", gff_path, attributes_map = TRUE, overwrite = TRUE)
#> Error in dbSendQuery(conn, statement, ...): Catalog Error: Table Function with name read_gff does not exist!
#> Did you mean "read_csv"?
#> 
#> LINE 1: CREATE TABLE annotations AS SELECT * FROM read_gff('/tmp/Rtmph0wBxf/temp_libpath1b84b3ec7f5/Rduckhts...
#>                                                   ^
#> ℹ Context: rapi_prepare
#> ℹ Error type: CATALOG
data <- dbGetQuery(con, "SELECT attributes FROM annotations LIMIT 5")
#> Error in dbSendQuery(conn, statement, ...): Catalog Error: Table with name annotations does not exist!
#> Did you mean "pg_settings"?
#> 
#> LINE 1: SELECT attributes FROM annotations LIMIT 5
#>                                ^
#> ℹ Context: rapi_prepare
#> ℹ Error type: CATALOG
keys <- extract_map_data(data$attributes, "keys")
#> Error in data$attributes: object of type 'closure' is not subsettable
name_values <- extract_map_data(data$attributes, "Name")
#> Error in data$attributes: object of type 'closure' is not subsettable
dbDisconnect(con, shutdown = TRUE)
```
