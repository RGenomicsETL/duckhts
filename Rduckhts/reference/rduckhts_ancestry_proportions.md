# Estimate projected ancestry proportions from frequency or dosage relations

Input columns are sample_id, chromosome, position, allele_a, allele_b
and frequency (or dosage), describing allele_b. Positions are positive
one-based whole numbers. Reference data may be a keyed wide relation
with consecutive PC1..PCn columns and one frequency column per group, or
a long relation with group_id and frequency plus a keyed pc/loading
relation. Corrections contain pc and coefficient. Groups (1..30), PCs
(1..64), frequencies and loadings must be complete and finite at each
contributing reference locus; group frequencies must lie in \[0, 1\].
Duplicate loci matching the input are rejected. Long references require
complete keyed rows during pivoting. Inputs and references must share an
assembly, uppercase biallelic SNV alleles and reference orientation.
Chromosome names are compared with \`duckhts_contig_key()\`: one leading
\`chr\` is removed in any letter case, \`M\` and \`MT\` become \`MT\`,
and \`X\` and \`Y\` are uppercased. Any other name must match byte for
byte, so \`01\` does not match \`1\`. A reference, wide or long, whose
chromosome column has an integer type compares that key as an integer,
so \`chr01\` joins \`1\`. Input duplicates at a sample/locus,
palindromic alleles and missing frequencies are dropped; other alleles
match directly, reversed, strand complemented or both. Reversed alleles
use 1-frequency. Missing genotypes are not imputed. Audit counts include
dropped physical input rows. The solver and correlation gates use full
precision; returned proportions are rounded to seven decimals. Aligned
scratch Parquet in \`tempdir()\` is removed on return. Correlation
failures return NULL proportions and a status.

## Usage

``` r
rduckhts_ancestry_proportions(
  con,
  input_table,
  reference_table,
  loadings_table,
  correction_table = NULL,
  input_kind = c("frequency", "dosage"),
  sum_to_one = TRUE,
  min_cor = 0.4,
  table_name = NULL,
  overwrite = FALSE,
  group_ids = NULL
)
```

## Arguments

- con:

  DuckDB connection with DuckHTS loaded.

- input_table, reference_table:

  Caller-owned relation names.

- loadings_table:

  Long-format PC loading relation; for wide references, pass the
  correction relation here with \`correction_table = NULL\`.

- correction_table:

  PC correction relation for long references.

- input_kind:

  \`frequency\` or diploid \`dosage\` (divided by two).

- sum_to_one:

  Require coefficients to sum to one; otherwise at most one.

- min_cor:

  Minimum predicted-frequency correlation, default 0.4 as in bigsnpr
  1.12.21. A stricter gate may suit a specific panel.

- table_name:

  Optional destination; NULL returns a data frame.

- overwrite:

  Replace an existing destination.

- group_ids:

  Optional named character vector for wide references, mapping frequency
  column names to output group IDs. Use distinct frequency column
  aliases when a group ID is also a PC column name or differs only by
  case.

## Value

A row per sample and reference group, with status and matching audit.
Data-frame results carry \`aligned_bytes\`, the size of temporary
aligned Parquet written during the call.
