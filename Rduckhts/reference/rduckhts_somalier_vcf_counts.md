# Extract Panel-Aligned Counts from VCF or BCF

Produce the complete sample-by-panel count relation consumed by the
Somalier-derived relatedness and contamination functions. \`FORMAT/AD\`
must declare \`Number=R,Type=Integer\`; A and B slots are matched by
exact REF/ALT identity, and \`other\` sums only the remaining
declared-allele slots. Missing sites and unavailable AD remain rows with
three NULL counts, distinct from measured zero depth. The panel can be
any typed table/view or ordinary Parquet file with the canonical six
panel identity columns.

## Usage

``` r
rduckhts_somalier_vcf_counts(
  con,
  path,
  panel_table = NULL,
  panel_parquet = NULL,
  samples = NULL,
  filter_policy = c("pass_or_unapplied", "include_all", "error"),
  table_name = NULL,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded.

- path:

  One VCF/BCF path or URI.

- panel_table:

  Name of the ordered panel table or view.

- panel_parquet:

  Ordinary panel Parquet path, instead of \`panel_table\`.

- samples:

  Optional HTSlib sample selector: comma-separated inclusion, leading
  \`^\` exclusion, \`"-"\` for all, or \`""\` for none.

- filter_policy:

  Record FILTER policy: \`"pass_or_unapplied"\` makes named failures
  unavailable, \`"include_all"\` uses their AD, and \`"error"\` rejects
  a selected panel record with a named failure.

- table_name:

  Optional output table. \`NULL\` returns a data frame.

- overwrite:

  Whether an existing output table may be replaced.

## Value

A data frame if \`table_name\` is \`NULL\`; otherwise invisible
\`TRUE\`.
