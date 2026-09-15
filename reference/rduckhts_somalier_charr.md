# Estimate Per-Sample Contamination with CHARR

Apply the Somalier-derived CHARR estimator to measured A/B/other count
evidence aligned to an ordered panel and population-B allele
frequencies. \`frequency_table\` contains the panel's six identity
columns plus \`population_b_af\`; it must cover every panel site exactly
once. Evidence and frequency identities, coordinates, and A/B
orientation are checked against the panel before their digests are
derived. A result with no usable homozygous-like evidence has status
\`no_evidence\` and a NULL estimate.

## Usage

``` r
rduckhts_somalier_charr(
  con,
  evidence_table = NULL,
  evidence_parquet = NULL,
  panel_table = NULL,
  panel_parquet = NULL,
  frequency_table = NULL,
  frequency_parquet = NULL,
  table_name = NULL,
  sample_ids = NULL,
  min_depth = 15,
  max_depth = 1e+06,
  hom_minor_rate = 0.12,
  hom_tail_alpha = 0.002,
  max_threshold_work = 1.6e+07,
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

- frequency_table:

  Name of the required population-frequency table or view.

- frequency_parquet:

  Path to the required population-frequency Parquet file, instead of
  \`frequency_table\`.

- table_name:

  Optional output table. \`NULL\` returns a data frame.

- sample_ids:

  Optional nonempty vector of distinct sample IDs to retain.

- min_depth:

  Minimum measured A+B depth for CHARR homozygous-like eligibility.

- max_depth:

  Maximum supported A+B depth for the binomial eligibility test, at most
  1,000,000.

- hom_minor_rate:

  Expected minor-read rate used to recognize homozygous-like anchors.

- hom_tail_alpha:

  Binomial upper-tail threshold for homozygous-like eligibility.

- max_threshold_work:

  Positive cumulative limit on exact binomial certification steps, at
  most 100,000,000. Distinct observed depths are certified once per call
  and shared across samples.

- max_sites:

  Positive per-sample panel capacity, at most 100,000,000.

- overwrite:

  Whether an existing output table may be replaced.

## Value

A data frame if \`table_name\` is \`NULL\`; otherwise invisible
\`TRUE\`.

## Details

Each of the evidence, panel, and frequency inputs is supplied as exactly
one named table/view or ordinary Parquet path. Optional sample selection
is exact: every requested ID must occur. Results contain the panel and
frequency digests, usable-site denominators, numerical status, and all
filter settings.
