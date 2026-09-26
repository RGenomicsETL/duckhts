# Read the VCF/BCF Sample Catalog

Return one row per selected sample with its original-header zero-based
\`sample_index\` and \`sample_name\`. Join this relation to
\`read_geno()\` calls from the same unchanged file; names are not
duplicated into every call.

## Usage

``` r
rduckhts_bcf_samples(con, path, samples = NULL)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- path:

  Path to the VCF/BCF file

- samples:

  Optional HTSlib sample selector: \`NULL\` or \`"-"\` keeps all, \`""\`
  keeps none, comma-separated names include samples, and a leading
  \`"^"\` excludes them. Unknown names error; selected samples retain
  header order.

## Value

A data frame with \`sample_index\` and \`sample_name\` columns.
