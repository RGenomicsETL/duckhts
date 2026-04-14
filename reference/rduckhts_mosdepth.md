# Native mosdepth-Compatible Fast-Mode Coverage Outputs

Writes an initial native mosdepth-compatible fast-mode output set for an
indexed BAM file.

## Usage

``` r
rduckhts_mosdepth(
  con,
  prefix,
  path,
  chrom = NULL,
  by = NULL,
  no_per_base = FALSE,
  threads = 0,
  flag = 1796,
  include_flag = 0,
  fast_mode = TRUE,
  mapq = 0,
  index_path = NULL,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- prefix:

  Output prefix for the mosdepth-style files

- path:

  Path to the input BAM file

- chrom:

  Optional chromosome name filter

- by:

  Optional fixed-width window size as a string or a BED file path

- no_per_base:

  Skip writing \`{prefix}.per-base.bed.gz\`

- threads:

  Number of BAM decompression threads

- flag:

  Excluded SAM flag mask, matching mosdepth's \`-F\`

- include_flag:

  Required SAM flag mask, matching mosdepth's \`-i\`

- fast_mode:

  Must currently remain \`TRUE\`

- mapq:

  Minimum mapping quality threshold

- index_path:

  Optional explicit BAM index path

- overwrite:

  Overwrite existing output files

## Value

A data frame describing the written output paths
