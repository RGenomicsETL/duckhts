# Count Agreement Between Overlapping Mates

Compares read 1 and read 2 of one fragment at each reference position
where both mates have an aligned base, and counts the pairs of base
qualities and bases. A real allele sits on the molecule and both mates
show it; a sequencing error usually shows in one mate only. So the
disagreement rate is about twice the instrument error rate, with no
reference, mask or contamination model. Errors made before sequencing,
such as PCR or damage errors, show in both mates and are not seen.

## Usage

``` r
rduckhts_bam_mate_overlap_counts(
  con,
  path,
  region = NULL,
  index_path = NULL,
  reference = NULL,
  reference_index_path = NULL,
  min_mapq = 20,
  require_flags = 0,
  exclude_flags = 3844,
  max_pending = 1e+06
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to the paired-end SAM, BAM or CRAM file

- region:

  Optional comma-separated contigs or intervals. It needs an alignment
  index and selects alignments by overlap; a pair counts only when both
  mates are selected

- index_path:

  Optional explicit alignment index path

- reference:

  Optional indexed reference FASTA. Required for CRAM; it is otherwise
  used only to resolve read bases stored as \`=\`

- reference_index_path:

  Optional explicit FASTA index path

- min_mapq:

  Minimum mapping quality, 0 to 255

- require_flags, exclude_flags:

  SAM flag masks that both mates must pass. The default
  \`exclude_flags\` leaves out unmapped, secondary, failed, duplicate
  and supplementary alignments

- max_pending:

  Most mates held at once while they wait for their mate, 1 to 100000000

## Value

A data frame with \`read1_quality\` and \`read2_quality\` (\`NA\` when
the read stores no qualities), \`read1_base\` and \`read2_base\` on the
forward strand of the reference, and \`positions\`. Read 1 is the mate
with flag 0x40. Rows with different bases are the disagreements.

## Details

Alignments must be in coordinate order. The first mate of a pair is held
until its mate arrives; more than \`max_pending\` mates held at once is
an error.

## Examples

``` r
library(DBI)
library(duckdb)

con <- rduckhts_connect()
bam <- system.file("extdata", "empty-tids.bam", package = "Rduckhts")
counts <- rduckhts_bam_mate_overlap_counts(con, bam)
sum(counts$positions[counts$read1_base != counts$read2_base]) / sum(counts$positions)
#> [1] 0.004428698
dbDisconnect(con, shutdown = TRUE)
```
