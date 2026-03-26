# Lift Over Variant Coordinates Against a Query

Applies the DuckHTS \`duckdb_liftover(...)\` table macro to rows from a
SQL query or table expression with chromosome and position columns, plus
optional reference and alternate alleles.

## Usage

``` r
rduckhts_liftover(
  con,
  query,
  chain_path,
  dst_fasta_ref,
  chrom_col = "chrom",
  pos_col = "pos",
  ref_col = NULL,
  alt_col = NULL,
  src_fasta_ref = NULL,
  max_snp_gap = 1,
  max_indel_inc = 250
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- query:

  SQL query or table expression to lift over

- chain_path:

  Path to a UCSC chain file

- dst_fasta_ref:

  Path to the destination FASTA reference

- chrom_col:

  Source chromosome column name

- pos_col:

  Source 1-based position column name

- ref_col:

  Optional reference allele column name

- alt_col:

  Optional alternate allele column name

- src_fasta_ref:

  Optional source FASTA reference

- max_snp_gap:

  Maximum chain block merge gap

- max_indel_inc:

  Maximum indel anchor expansion

## Value

A data frame with source columns, lifted coordinates/alleles, and
warnings.
