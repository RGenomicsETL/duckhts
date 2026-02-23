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
#> Error in duckdb_result(connection = conn, stmt_lst = stmt_lst, arrow = arrow): Invalid Error: IO Error: Extension "/tmp/Rtmph0wBxf/temp_libpath1b84b3ec7f5/Rduckhts/duckhts_extension/build/duckhts.duckdb_extension" could not be loaded: libhts.so.3: cannot open shared object file: No such file or directory
#> ℹ Context: rapi_execute
#> ℹ Error type: INVALID
bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
rduckhts_bcf(con, "variants", bcf_path, overwrite = TRUE)
#> Error in dbSendQuery(conn, statement, ...): Catalog Error: Table Function with name read_bcf does not exist!
#> Did you mean "read_csv"?
#> 
#> LINE 1: CREATE TABLE variants AS SELECT * FROM read_bcf('/tmp/Rtmph0wBxf/temp_libpath1b84b3ec7f5/Rduckhts...
#>                                                ^
#> ℹ Context: rapi_prepare
#> ℹ Error type: CATALOG
complex_cols <- detect_complex_types(con, "variants")
#> Error in dbSendQuery(conn, statement, ...): Catalog Error: Table with name variants does not exist!
#> Did you mean "pg_views"?
#> 
#> LINE 1: DESCRIBE variants
#>                  ^
#> ℹ Context: rapi_prepare
#> ℹ Error type: CATALOG
print(complex_cols)
#> Error: object 'complex_cols' not found
dbDisconnect(con, shutdown = TRUE)
```
