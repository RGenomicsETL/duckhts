# Native mosdepth-Compatible Coverage Outputs

Writes native mosdepth-compatible coverage outputs for indexed BAM or
CRAM input.

## Usage

``` r
rduckhts_mosdepth(
  con,
  prefix,
  path,
  chrom = NULL,
  by = NULL,
  fasta = NULL,
  read_groups = NULL,
  no_per_base = FALSE,
  threads = 2,
  processing_threads = 2,
  flag = 1796,
  include_flag = 0,
  fast_mode = FALSE,
  fragment_mode = FALSE,
  use_median = FALSE,
  mapq = 0,
  min_frag_len = -1,
  max_frag_len = -1,
  precision_digits = 2,
  quantize = NULL,
  thresholds = NULL,
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

  Path to the input BAM or CRAM file

- chrom:

  Optional chromosome name filter

- by:

  Optional fixed-width window size as a string or a BED file path

- fasta:

  Optional reference FASTA path for CRAM input when required

- read_groups:

  Optional comma-separated read-group IDs, matching mosdepth's \`-R\`

- no_per_base:

  Skip writing \`{prefix}.per-base.bed.gz\`

- threads:

  Number of BAM decompression threads

- processing_threads:

  Number of parallel contig processing threads (0 = sequential)

- flag:

  Excluded SAM flag mask, matching mosdepth's \`-F\`

- include_flag:

  Required SAM flag mask, matching mosdepth's \`-i\`

- fast_mode:

  Logical. If \`TRUE\`, use mosdepth fast mode. Defaults to \`FALSE\`,
  matching upstream mosdepth.

- fragment_mode:

  Logical. If \`TRUE\`, count full fragment insert spans for proper
  pairs, matching mosdepth's \`-a\`. Cannot be combined with \`fast_mode
  = TRUE\`.

- use_median:

  Logical. If \`TRUE\`, write \`by\` region values as medians instead of
  means, matching mosdepth's \`-m\`.

- mapq:

  Minimum mapping quality threshold

- min_frag_len:

  Minimum absolute template length to keep, matching mosdepth's \`-l\`

- max_frag_len:

  Maximum absolute template length to keep, matching mosdepth's \`-u\`

- precision_digits:

  Number of decimal places to write in the text outputs

- quantize:

  Optional mosdepth-style quantize specification such as \`":1:4:"\`

- thresholds:

  Optional comma-separated coverage thresholds for \`by\`, matching
  mosdepth's \`-T\`

- index_path:

  Optional explicit BAM index path

- overwrite:

  Overwrite existing output files

## Value

A data frame describing the written output paths
