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
