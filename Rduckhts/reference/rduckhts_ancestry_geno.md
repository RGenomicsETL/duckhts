# Estimate ancestry from a materialized \`read_geno\` relation

Converts each diploid, biallelic call to ALT dosage / 2. A haploid,
partially missing, or multiallelic call has NULL dosage and is counted
in the missing-site audit. Original sample indices identify the calls
unless a \`read_bcf_samples\` relation is provided for original-header
sample names. That relation must give every selected sample index one
distinct non-null name.

## Usage

``` r
rduckhts_ancestry_geno(
  con,
  geno_table,
  reference_table,
  loadings_table,
  correction_table,
  samples_table = NULL,
  sum_to_one = TRUE,
  min_cor = 0.4,
  non_reference_only
)
```

## Arguments

- con:

  Connection with DuckHTS loaded.

- geno_table:

  Materialized \`read_geno\` rows with CHROM, POS, REF, ALT and calls
  columns. Each site must contain a call for every selected sample,
  including homozygous-reference calls.

- reference_table, loadings_table, correction_table:

  Reference products.

- samples_table:

  Optional \`read_bcf_samples\` relation with sample_index and
  sample_name.

- sum_to_one, min_cor:

  Quality/constraint arguments for the solver.

- non_reference_only:

  Required declaration of the \`read_geno\` setting used to create
  \`geno_table\`. Must be \`FALSE\`: sparse call lists omit zero dosages
  and cannot be identified from the materialized relation's schema.

## Value

A row per sample and group with the matching audit.
