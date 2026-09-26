# Create BAM Pileup Table

Creates a DuckDB table from a region-scoped BAM pileup using the DuckHTS
extension. This is a compact base/quality pileup view backed by htslib's
pileup engine; it is not samtools mpileup text parity.

## Usage

``` r
rduckhts_pileup(
  con,
  table_name,
  path,
  region,
  index_path = NULL,
  min_mapq = 0,
  flag_mask = 1796,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- table_name:

  Name for the created table

- path:

  Path to the BAM file

- region:

  Required genomic region (e.g., "chr1:1000-2000")

- index_path:

  Optional explicit path to the BAM index (.bai/.csi)

- min_mapq:

  Minimum mapping quality to include in the pileup

- flag_mask:

  Bitmask of SAM flags to exclude before pileup construction. The
  default `1796` matches samtools depth-style filtering of unmapped,
  secondary, QC-fail, and duplicate reads.

- overwrite:

  Logical. If TRUE, overwrites existing table

## Value

Invisible TRUE on success
