# Create VCF/BCF Table

Creates a DuckDB table from a VCF or BCF file using the DuckHTS
extension. This follows the RBCFTools pattern of creating a table that
can be queried.

## Usage

``` r
rduckhts_bcf(
  con,
  table_name,
  path,
  region = NULL,
  index_path = NULL,
  tidy_format = FALSE,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- table_name:

  Name for the created table

- path:

  Path to the VCF/BCF file

- region:

  Optional genomic region (e.g., "chr1:1000-2000")

- index_path:

  Optional explicit path to index file (.csi/.tbi)

- tidy_format:

  Logical. If TRUE, FORMAT columns are returned in tidy format

- overwrite:

  Logical. If TRUE, overwrites existing table

## Value

Invisible TRUE on success

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
dbGetQuery(con, "SELECT * FROM variants LIMIT 2")
#> Error in dbSendQuery(conn, statement, ...): Catalog Error: Table with name variants does not exist!
#> Did you mean "pg_views"?
#> 
#> LINE 1: SELECT * FROM variants LIMIT 2
#>                       ^
#> ℹ Context: rapi_prepare
#> ℹ Error type: CATALOG
dbDisconnect(con, shutdown = TRUE)
```
