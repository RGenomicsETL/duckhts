# Count Aligned Bases Against the Reference

Counts the aligned bases of a SAM, BAM or CRAM file against the
reference by mate, cycle, base quality and substitution. With a mask of
the known variants of the sample or its population, the mismatches that
remain are mostly errors of the read, the library or the alignment. A
variant that the mask does not hold is still counted, so their rate is
an upper bound on the error rate.

## Usage

``` r
rduckhts_bam_mismatch_counts(
  con,
  path,
  reference,
  region = NULL,
  mask = NULL,
  index_path = NULL,
  reference_index_path = NULL,
  mask_index_path = NULL,
  min_mapq = 20,
  require_flags = 0,
  exclude_flags = 3844,
  indel_flank = 5
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to the SAM, BAM or CRAM file

- reference:

  Path to the indexed reference FASTA

- region:

  Optional comma-separated contigs or intervals. It needs an alignment
  index and selects alignments by overlap

- mask:

  Optional path to an indexed BCF or a bgzip-compressed VCF with a tabix
  index. Its positions are left out

- index_path, reference_index_path, mask_index_path:

  Optional explicit index paths

- min_mapq:

  Minimum mapping quality, 0 to 255

- require_flags, exclude_flags:

  SAM flag masks. The default \`exclude_flags\` leaves out unmapped,
  secondary, failed, duplicate and supplementary alignments

- indel_flank:

  Read bases to leave out on each side of an insertion, a deletion, a
  reference skip or a soft clip, and reference bases on each side of a
  mask record with an indel or a symbolic allele; 0 to 1000

## Value

A data frame with \`mate\` (1 or 2 for a paired read, 0 otherwise),
\`cycle\`, \`base_quality\` (\`NA\` when the read stores no qualities),
\`reference_base\` and \`read_base\` on the strand of the read, and
\`bases\`. Rows with equal bases are the matches.

## Details

Contig names of the alignment header, the reference and the mask are
compared byte for byte. A contig that the reference or the mask does not
know is an error.

## Examples

``` r
library(DBI)
library(duckdb)

con <- rduckhts_connect()
extdata <- function(name) system.file("extdata", name, package = "Rduckhts")
counts <- rduckhts_bam_mismatch_counts(
  con, extdata("bam_mismatch.bam"), extdata("bam_mismatch.fa"),
  mask = extdata("bam_mismatch.mask.bcf")
)
sum(counts$bases[counts$reference_base != counts$read_base]) / sum(counts$bases)
#> [1] 0.04040404
dbDisconnect(con, shutdown = TRUE)
```
