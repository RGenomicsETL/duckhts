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
#> [1] TRUE
bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
rduckhts_bam(con, "reads", bam_path, overwrite = TRUE)
dbGetQuery(con, "SELECT COUNT(*) FROM reads WHERE FLAG & 4 = 0")
#>   count_star()
#> 1          112
dbDisconnect(con, shutdown = TRUE)
```
