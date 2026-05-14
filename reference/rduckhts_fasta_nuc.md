# Compute FASTA Interval Nucleotide Composition

Computes bedtools nuc-style nucleotide composition over either a BED
file or generated fixed-width bins.

## Usage

``` r
rduckhts_fasta_nuc(
  con,
  path,
  bed_path = NULL,
  bin_width = NULL,
  region = NULL,
  index_path = NULL,
  gzi_path = NULL,
  bed_index_path = NULL,
  include_seq = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to the FASTA file

- bed_path:

  Optional BED path. Supply exactly one of \`bed_path\` or
  \`bin_width\`.

- bin_width:

  Optional fixed bin width in base pairs

- region:

  Optional FASTA region filter

- index_path:

  Optional explicit FASTA index path

- gzi_path:

  Optional explicit BGZF FASTA block index path (.gzi) for bgzipped
  FASTA inputs when the sidecar is not colocated with the FASTA.

- bed_index_path:

  Optional explicit BED tabix index path

- include_seq:

  Include the fetched interval sequence

## Value

A data frame with interval composition statistics
