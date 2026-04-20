# Native BAM/CRAM BED Regional Coverage Summary

Computes samtools coverage-like regional summaries for BAM or CRAM input
over a BED target set, with DuckHTS-specific pre/post-filter and
strand-aware post-filter outputs.

## Usage

``` r
rduckhts_bam_bed_coverage(
  con,
  path,
  bed_path,
  reference = NULL,
  index_path = NULL,
  bed_index_path = NULL,
  mapq = 0,
  min_baseq = 0,
  min_read_len = 0,
  require_flags = 0,
  exclude_flags = 1796,
  min_depth = 1,
  max_depth = 1e+06,
  decompression_threads = 0,
  fragment_mode = FALSE,
  strand_outputs = TRUE,
  processing_threads = 0
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to the input BAM or CRAM file

- bed_path:

  Path to the input BED file

- reference:

  Optional reference FASTA path for CRAM input when required

- index_path:

  Optional explicit BAM/CRAM index path

- bed_index_path:

  Optional explicit BED index path (reserved for future use)

- mapq:

  Minimum mapping quality threshold for post-filter summaries

- min_baseq:

  Minimum base quality threshold for post-filter base-level summaries

- min_read_len:

  Minimum read length threshold for post-filter summaries

- require_flags:

  Required SAM flag mask

- exclude_flags:

  Excluded SAM flag mask. Defaults to samtools coverage's
  \`UNMAP\|SECONDARY\|QCFAIL\|DUP\` mask.

- min_depth:

  Minimum depth threshold for covered-base and mean-depth summaries

- max_depth:

  Maximum per-position depth cap. Set \`0\` to remove the cap.

- decompression_threads:

  Integer. Number of htslib decompression worker threads to use for
  BAM/CRAM input. \`0\` disables htslib worker threads.

- fragment_mode:

  Logical. Reserved for future fragment-level semantics.

- strand_outputs:

  Logical. Emit forward/reverse post-filter summary columns.

- processing_threads:

  Reserved for future parallel interval processing.

## Value

A data frame with one row per BED interval and pre/post regional
summaries
