# Native Fixed-Width BAM/CRAM Bin Counts

Count read starts into fixed-width genomic bins with optional duplicate
handling and optional per-bin GC and MAPQ summary statistics.

## Usage

``` r
rduckhts_bam_bin_counts(
  con,
  path,
  bin_width,
  chrom = NULL,
  reference = NULL,
  index_path = NULL,
  mapq = 0,
  require_flags = 0,
  exclude_flags = 0,
  rmdup = "none",
  stats = NULL
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to the input BAM or CRAM file

- bin_width:

  Positive fixed bin width in bases

- chrom:

  Optional chromosome name filter

- reference:

  Optional reference FASTA path for CRAM input when required, and for
  reference-GC output when \`stats\` includes \`"gc"\`

- index_path:

  Optional explicit BAM/CRAM index path

- mapq:

  Minimum mapping quality threshold applied after duplicate logic

- require_flags:

  Required SAM flag mask

- exclude_flags:

  Excluded SAM flag mask

- rmdup:

  Duplicate handling mode: \`"none"\`, \`"flag"\`, or \`"streaming"\`

- stats:

  Optional comma-separated subset of \`"gc"\` and \`"mq"\`

## Value

A data frame with one row per fixed-width bin across the selected contig
span, including zero-count bins, plus total, forward, reverse, and
optional GC/MAPQ summary columns
