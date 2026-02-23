# Create SAM/BAM/CRAM Table

Creates a DuckDB table from SAM, BAM, or CRAM files using the DuckHTS
extension.

## Usage

``` r
rduckhts_bam(
  con,
  table_name,
  path,
  region = NULL,
  index_path = NULL,
  reference = NULL,
  standard_tags = NULL,
  auxiliary_tags = NULL,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- table_name:

  Name for the created table

- path:

  Path to the SAM/BAM/CRAM file

- region:

  Optional genomic region (e.g., "chr1:1000-2000")

- index_path:

  Optional explicit path to index file (.bai/.csi/.crai)

- reference:

  Optional reference file path for CRAM files

- standard_tags:

  Logical. If TRUE, include typed standard SAMtags columns

- auxiliary_tags:

  Logical. If TRUE, include AUXILIARY_TAGS map of non-standard tags

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
#> Error in duckdb_result(connection = conn, stmt_lst = stmt_lst, arrow = arrow): Invalid Error: IO Error: Extension "/tmp/RtmpWan7sJ/temp_libpath1b8825a39529/Rduckhts/duckhts_extension/build/duckhts.duckdb_extension" could not be loaded: libhts.so.3: cannot open shared object file: No such file or directory
#> ℹ Context: rapi_execute
#> ℹ Error type: INVALID
bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
rduckhts_bam(con, "reads", bam_path, overwrite = TRUE)
#> Error in dbSendQuery(conn, statement, ...): Catalog Error: Table Function with name read_bam does not exist!
#> Did you mean "read_blob"?
#> 
#> LINE 1: CREATE TABLE reads AS SELECT * FROM read_bam('/tmp/RtmpWan7sJ/temp_libpath1b8825a39529/Rduckhts...
#>                                             ^
#> ℹ Context: rapi_prepare
#> ℹ Error type: CATALOG
dbGetQuery(con, "SELECT COUNT(*) FROM reads WHERE FLAG & 4 = 0")
#> Error in dbSendQuery(conn, statement, ...): Catalog Error: Table with name reads does not exist!
#> Did you mean "pg_prepared_statements"?
#> 
#> LINE 1: SELECT COUNT(*) FROM reads WHERE FLAG & 4 = 0
#>                              ^
#> ℹ Context: rapi_prepare
#> ℹ Error type: CATALOG
dbDisconnect(con, shutdown = TRUE)
```
