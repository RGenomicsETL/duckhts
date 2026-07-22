# Create a BigWig Signal Table

Materializes stored zero-based, half-open BigWig intervals through the
DuckHTS extension. Region filters use htslib's one-based inclusive
syntax; a character vector is combined into one multi-region request and
overlapping requests emit each stored interval once.

## Usage

``` r
rduckhts_bigwig(
  con,
  table_name,
  path,
  region = NULL,
  blocks_per_iteration = 64L,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded.

- table_name:

  Name for the created table, or `NULL` to create the `bigwig_data`
  view.

- path:

  Local path or URL to a BigWig file.

- region:

  Optional character vector of genomic regions such as
  `c("chr1:1000-2000", "chr2:1-500")`. A single comma-separated string
  is also accepted.

- blocks_per_iteration:

  Positive integer number of indexed BigWig data blocks decoded per
  iterator batch.

- overwrite:

  Logical. If `TRUE`, replace an existing table.

## Value

Invisibly returns `TRUE`.
