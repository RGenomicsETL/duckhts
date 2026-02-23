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
#> Error in duckdb_result(connection = conn, stmt_lst = stmt_lst, arrow = arrow): Invalid Error: IO Error: Extension "/tmp/RtmpWan7sJ/temp_libpath1b8825a39529/Rduckhts/duckhts_extension/build/duckhts.duckdb_extension" could not be loaded: libhts.so.3: cannot open shared object file: No such file or directory
#> ℹ Context: rapi_execute
#> ℹ Error type: INVALID
bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
rduckhts_bcf(con, "variants", bcf_path, overwrite = TRUE)
#> Error in dbSendQuery(conn, statement, ...): Catalog Error: Table Function with name read_bcf does not exist!
#> Did you mean "read_csv"?
#> 
#> LINE 1: CREATE TABLE variants AS SELECT * FROM read_bcf('/tmp/RtmpWan7sJ/temp_libpath1b8825a39529/Rduckhts...
#>                                                ^
#> ℹ Context: rapi_prepare
#> ℹ Error type: CATALOG
data <- dbGetQuery(con, "SELECT ALT FROM variants LIMIT 5")
#> Error in dbSendQuery(conn, statement, ...): Catalog Error: Table with name variants does not exist!
#> Did you mean "pg_views"?
#> 
#> LINE 1: SELECT ALT FROM variants LIMIT 5
#>                         ^
#> ℹ Context: rapi_prepare
#> ℹ Error type: CATALOG
first_alt <- extract_array_element(data$ALT, 1)
#> Error in data$ALT: object of type 'closure' is not subsettable
all_alts <- extract_array_element(data$ALT)
#> Error in data$ALT: object of type 'closure' is not subsettable
dbDisconnect(con, shutdown = TRUE)
```
