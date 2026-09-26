# Normalize Variant Alleles with bcftools-style Semantics

Applies the DuckHTS \`duckhts_bcftools_norm(...)\` table macro to rows
from a SQL query or table expression. Input rows must expose chromosome,
1-based position, reference allele, and alternate allele columns.
Alternate alleles may be supplied either as a comma-delimited
\`VARCHAR\` or as a \`VARCHAR\[\]\` list, matching the common DuckDB
representations used by plain tables and \`read_bcf(...)\`.

## Usage

``` r
rduckhts_bcftools_norm(
  con,
  query,
  fasta_ref,
  chrom_col = "chrom",
  pos_col = "pos",
  ref_col = "ref",
  alt_col = "alt",
  split_multiallelic = FALSE,
  end_pos_col = NULL,
  svlen_col = NULL,
  fasta_index_path = NULL,
  gzi_path = NULL
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- query:

  SQL query or table expression to normalize

- fasta_ref:

  Path to the reference FASTA

- chrom_col:

  Source chromosome column name

- pos_col:

  Source 1-based position column name

- ref_col:

  Source reference allele column name

- alt_col:

  Source alternate allele column name (\`VARCHAR\` or \`VARCHAR\[\]\`)

- split_multiallelic:

  If \`TRUE\`, split multiallelic sites before normalization so
  \`alt_normed\` is emitted as \`VARCHAR\` plus \`alt_index\`. If
  \`FALSE\` (default), keep sites intact and emit \`alt_normed\` as
  \`VARCHAR\[\]\`.

- end_pos_col:

  Optional source column name containing an END-like 1-based end
  coordinate for symbolic deletions.

- svlen_col:

  Optional source column name containing an SVLEN-like signed length for
  symbolic duplications.

- fasta_index_path:

  Optional explicit \`.fai\` sidecar path.

- gzi_path:

  Optional explicit \`.gzi\` sidecar path for bgzipped FASTA.

## Value

A data frame with the original columns plus \`pos_normed\`,
\`end_pos_normed\`, \`ref_normed\`, \`alt_normed\`, \`normed\`, and
\`norm_status\`. In split mode the result additionally includes
\`alt_index\`.
