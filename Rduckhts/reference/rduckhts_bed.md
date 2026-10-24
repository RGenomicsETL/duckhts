# Create BED Table

Creates a DuckDB table from a BED file using the DuckHTS extension.

## Usage

``` r
rduckhts_bed(
  con,
  table_name,
  path,
  region = NULL,
  index_path = NULL,
  scan_mode = NULL,
  overwrite = FALSE,
  error_policy = NULL
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- table_name:

  Name for the created table

- path:

  Path to the BED file

- region:

  Optional genomic region for tabix-backed BED queries

- index_path:

  Optional explicit path to a BED tabix index

- scan_mode:

  Optional scan mode. Use `"auto"` (default extension behavior) or
  `"sequential"` to force full-file streaming/counting instead of
  index-backed count paths. Sequential mode is incompatible with
  `region`.

- overwrite:

  Logical. If TRUE, overwrites an existing table

- error_policy:

  Optional short-line policy: `"error"` (default), `"skip"`, or
  `"report"`. Report adds diagnostic columns and requires a full-file
  scan rather than a region query.

## Value

Invisible TRUE on success

## Examples

``` r
con <- rduckhts_connect()

# Intervals are 0-based and half-open, like BED: keep BED starts as they are.
DBI::dbExecute(con, 'CREATE TABLE tiles AS SELECT * FROM (VALUES
  (\'CHROMOSOME_I\', 0, 10), (\'CHROMOSOME_I\', 10, 20)) t(chrom, start, "end")')
#> [1] 2
DBI::dbGetQuery(con, "SELECT * FROM duckhts_cgranges_from_table(
  'tiles_idx', 'tiles', 'chrom', 'start', 'end')")
#>   indexed
#> 1    TRUE

# end is a reserved SQL word: quote the column as "end" in the select list and
# pass that column, not the string 'end', to the scalar functions.
rduckhts_bed(con, "targets", system.file("extdata", "targets.bed", package = "Rduckhts"))
DBI::dbGetQuery(con, 'SELECT chrom, start, "end",
  duckhts_cgranges_count_overlaps(\'tiles_idx\', chrom, start, "end") AS tiles
  FROM targets')
#>            chrom start end tiles
#> 1   CHROMOSOME_I     0  10     1
#> 2   CHROMOSOME_I    10  20     1
#> 3  CHROMOSOME_II     0   8     0
#> 4 CHROMOSOME_III     0   6     0

DBI::dbDisconnect(con, shutdown = TRUE)
```
