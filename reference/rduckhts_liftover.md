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
  max_indel_inc = 250,
  lift_mt = FALSE,
  end_pos_col = NULL,
  no_left_align = FALSE
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

- lift_mt:

  If FALSE (default), mitochondrial variants with matching
  source/destination contig lengths are passed through with only contig
  rename. If TRUE, MT variants are lifted through the chain like any
  other contig.

- end_pos_col:

  Optional column name containing INFO/END positions (1-based) to lift
  alongside the primary position. When provided, the output includes a
  \`dest_end\` column with the lifted end position.

- no_left_align:

  If FALSE (default), lifted indels are left-aligned against the
  destination reference. Set TRUE to skip left-alignment, mirroring
  `--no-left-align` in `bcftools +liftover`.

## Value

A data frame with source columns, lifted coordinates/alleles, and
warnings.
