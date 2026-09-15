# Prepare Somalier-Derived Sample Sketches

Build packed, panel-verified relatedness sketches from measured
A/B/other count evidence. The panel must contain \`assembly\`,
zero-based \`site_index\`, \`region\`, one-based \`position\`, and
uppercase single-base \`allele_a\` and \`allele_b\`. Evidence must
contain the same site identity columns plus \`sample_id\` and nullable
count columns \`a\`, \`b\`, and \`other\`. All three counts are NULL for
unavailable evidence; three measured zeros are not unavailable. The
native SQL preparation checks every evidence site's geometry and A/B
orientation against the ordered panel before computing its digest. Panel
alleles must be distinct uppercase single-base A/C/G/T with lexical A \<
B; the exact X/Y aliases excluded by Somalier v0.3.4 are rejected. Other
contig aliases cannot be classified biologically from the region string.
The three-state calculation assumes diploid sites; count evidence alone
does not prove sample ploidy.

## Usage

``` r
rduckhts_somalier_sketches(
  con,
  evidence_table = NULL,
  evidence_parquet = NULL,
  panel_table = NULL,
  panel_parquet = NULL,
  table_name = NULL,
  sample_ids = NULL,
  min_depth = 7,
  min_het_balance = 0.3,
  hom_balance_cutoff = 0.01,
  max_sites = 1e+06,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded.

- evidence_table:

  Name of an ordinary evidence table or view.

- evidence_parquet:

  Path to an evidence Parquet file, instead of \`evidence_table\`.

- panel_table:

  Name of the required ordered panel table or view. Panel assembly and
  region values are each limited to 1,024 bytes.

- panel_parquet:

  Path to the required ordered panel Parquet file, instead of
  \`panel_table\`.

- table_name:

  Optional output table. \`NULL\` returns a data frame.

- sample_ids:

  Optional nonempty vector of distinct sample IDs to retain.

- min_depth:

  Minimum A+B count depth for a relatedness genotype call.

- min_het_balance:

  Lower B/(A+B) balance accepted as heterozygous.

- hom_balance_cutoff:

  B/(A+B) balance below which a site is homozygous A; its upper
  symmetric limit determines homozygous B.

- max_sites:

  Positive per-sample panel capacity, at most 100,000,000.

- overwrite:

  Whether an existing output table may be replaced.

## Value

A data frame if \`table_name\` is \`NULL\`; otherwise invisible
\`TRUE\`.

## Details

Supply each source as either a table/view name or an ordinary Parquet
path. Parquet inputs are exposed through query-scoped temporary views;
no private sketch format or user-supplied panel digest is involved. The
result has one \`sketch\` struct per selected sample. Its packed words
can be persisted with DuckDB's usual Parquet \`COPY\` statement.
