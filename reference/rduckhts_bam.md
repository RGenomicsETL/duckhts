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
  standard_tags = FALSE,
  auxiliary_tags = FALSE,
  sequence_encoding = NULL,
  quality_representation = NULL,
  cigar_representation = NULL,
  decompression_threads = 2,
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

  Logical. If TRUE, include typed standard SAMtags columns. Default
  FALSE.

- auxiliary_tags:

  Logical. If TRUE, include AUXILIARY_TAGS map of non-standard tags.
  Default FALSE.

- sequence_encoding:

  Character. Sequence encoding for the SEQ column: `"string"` (default)
  returns decoded bases as `VARCHAR`; `"nt16"` returns raw htslib nt16
  4-bit codes as `UTINYINT[]`.

- quality_representation:

  Character. Quality representation for the QUAL column: `"string"`
  (default) returns canonical Phred+33 text; `"phred"` returns raw Phred
  values as `UTINYINT[]`.

- cigar_representation:

  Character. CIGAR representation for the CIGAR column: `"string"`
  (default) returns SAM text such as `"36M"`; `"binary"` returns packed
  BAM operations as `UINTEGER[]` where each element is
  `(len << 4) | op`.

- decompression_threads:

  Integer. Number of htslib decompression worker threads per file
  handle. Default `2`. Use `0` to disable worker threads.

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
