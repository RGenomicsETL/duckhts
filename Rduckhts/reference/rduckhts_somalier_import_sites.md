# Import an Already Selected Somalier Sites VCF or BCF

Convert an existing Somalier-compatible sites file into the canonical
typed panel and population-frequency relation used by DuckHTS
extraction, relatedness, and contamination functions. REF and ALT are
oriented into lexical A/B order and alternate-allele frequency is
flipped with the alleles, so \`population_b_af\` always describes
\`allele_b\`. Exact Somalier v0.3.4 X/Y aliases are excluded, matching
its autosomal frequency importer.

## Usage

``` r
rduckhts_somalier_import_sites(
  con,
  path,
  assembly,
  max_sites = 1e+06,
  table_name = NULL,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded.

- path:

  One already selected sites VCF/BCF path or URI. \`INFO/AF\` must be
  declared \`Number=A,Type=Float\` and every retained record must be one
  canonical biallelic SNV with one finite AF value.

- assembly:

  Nonempty assembly identifier attached to every site.

- max_sites:

  Positive input-site limit, at most 100,000,000.

- table_name:

  Optional output table. \`NULL\` returns a data frame.

- overwrite:

  Whether an existing output table may be replaced.

## Value

A data frame if \`table_name\` is \`NULL\`; otherwise invisible
\`TRUE\`.

## Details

This function does not select sites from a population VCF. Somalier's
\`find-sites\` algorithm has separate AF/AN, QC, interval-exclusion, and
spacing semantics and is not implied by importing its output.
