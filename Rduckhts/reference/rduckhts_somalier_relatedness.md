# Compare Somalier-Derived Sample Sketches

Compute the named relatedness and concordance statistics for either
every distinct pair in one sketch relation or the ordered pairs in
\`pairs_table\`. The sketch relation must contain one non-NULL
\`sketch\` struct per sample. Selected pairs are an ordinary relation
with \`sample_a\` and \`sample_b\` columns. Missing or duplicate sample
and pair identities error instead of silently dropping or multiplying
requested comparisons. The native kernel checks assembly, ordered-panel
digest, classification settings, mask shape, and mask contents for each
comparison. No SQL row-order guarantee is implied.

## Usage

``` r
rduckhts_somalier_relatedness(
  con,
  sketches_table = NULL,
  sketches_parquet = NULL,
  pairs_table = NULL,
  table_name = NULL,
  max_sites = 1e+06,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded.

- sketches_table:

  Name of a prepared sketch table or view.

- sketches_parquet:

  Path to ordinary Parquet-persisted sketches instead of
  \`sketches_table\`.

- pairs_table:

  Optional name of an ordered-pair table or view; \`NULL\` requests all
  distinct unordered sample pairs.

- table_name:

  Optional output table. \`NULL\` returns a data frame, which should be
  used only for a result small enough to fit in R memory.

- max_sites:

  Positive per-pair panel capacity, at most 100,000,000.

- overwrite:

  Whether an existing output table may be replaced.

## Value

A data frame if \`table_name\` is \`NULL\`; otherwise invisible
\`TRUE\`.
