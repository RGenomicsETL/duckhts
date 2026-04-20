# samtools idxstats-Compatible Alignment Summary

Writes samtools idxstats-compatible alignment summary output for BAM,
CRAM, or SAM input.

## Usage

``` r
rduckhts_samtools_idxstats(
  con,
  path,
  output = NULL,
  index_path = NULL,
  threads = 0,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to the input alignment file

- output:

  Optional output path for the written idxstats text file

- index_path:

  Optional explicit BAM/CRAM index path

- threads:

  htslib decompression thread count for scan fallback

- overwrite:

  Overwrite an existing output file

## Value

A data frame with \`success\`, \`path\`, \`output_path\`,
\`used_index_fast_path\`, and \`error_message\`
