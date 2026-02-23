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
#> [1] TRUE
bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
rduckhts_bcf(con, "variants", bcf_path, overwrite = TRUE)
dbGetQuery(con, "SELECT * FROM variants LIMIT 2")
#>   CHROM     POS   ID REF ALT QUAL FILTER INFO_TEST INFO_DP4 INFO_AC INFO_AN
#> 1     1 3000150 <NA>   C   T 59.2   PASS        NA       NA       2       4
#> 2     1 3000151 <NA>   C   T 59.2   PASS        NA       NA       2       4
#>   INFO_INDEL INFO_STR FORMAT_TT_A FORMAT_GT_A FORMAT_GQ_A FORMAT_DP_A
#> 1      FALSE     <NA>        NULL         0/1         245          NA
#> 2      FALSE     <NA>        NULL         0/1         245          32
#>   FORMAT_GL_A FORMAT_TT_B FORMAT_GT_B FORMAT_GQ_B FORMAT_DP_B FORMAT_GL_B
#> 1        NULL        NULL         0/1         245          NA        NULL
#> 2        NULL        NULL         0/1         245          32        NULL
dbDisconnect(con, shutdown = TRUE)
```
