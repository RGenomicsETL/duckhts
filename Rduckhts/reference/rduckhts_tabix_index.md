# Build Tabix Index

Builds a tabix index for a BGZF-compressed text file using the DuckHTS
extension.

## Usage

``` r
rduckhts_tabix_index(
  con,
  path,
  preset = "vcf",
  index_path = NULL,
  min_shift = 0,
  threads = 4,
  seq_col = NULL,
  start_col = NULL,
  end_col = NULL,
  comment_char = NULL,
  skip_lines = NULL
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to the BGZF-compressed input file

- preset:

  Optional preset such as \`"vcf"\`, \`"bed"\`, \`"gff"\`, or \`"sam"\`

- index_path:

  Optional explicit output path for the created index

- min_shift:

  Index format selector used by htslib

- threads:

  htslib indexing thread count

- seq_col, start_col, end_col:

  Optional explicit tabix coordinate columns

- comment_char:

  Optional tabix comment/header prefix

- skip_lines:

  Optional fixed number of header lines to skip

## Value

A data frame with \`success\`, \`index_path\`, and \`index_format\`
