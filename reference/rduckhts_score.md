# Compute Polygenic Scores

Calls the DuckHTS \`bcftools_score(...)\` table function to compute
sample-level polygenic scores from one genotype VCF/BCF file and one or
more summary-statistics files.

## Usage

``` r
rduckhts_score(
  con,
  bcf_path,
  summary_path = NULL,
  use = NULL,
  columns = "PLINK",
  columns_file = NULL,
  q_score_thr = NULL,
  summaries_list_file = NULL,
  log_path = NULL,
  use_variant_id = FALSE,
  counts = FALSE,
  samples = NULL,
  force_samples = FALSE,
  regions = NULL,
  regions_file = NULL,
  regions_overlap = 1,
  targets = NULL,
  targets_file = NULL,
  targets_overlap = 0,
  apply_filters = NULL,
  include = NULL,
  exclude = NULL
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- bcf_path:

  Path to genotype VCF/BCF file

- summary_path:

  Path(s) to summary-statistics file(s). A character vector computes
  multiple TSV/SSF PRS columns in one genotype scan. Use \`NULL\` with
  \`summaries_list_file\` to read paths from a file.

- use:

  Optional dosage source (\`"GT"\`, \`"DS"\`, \`"HDS"\`, \`"AP"\`,
  \`"GP"\`, \`"AS"\`)

- columns:

  Optional summary preset (\`"PLINK"\`, \`"PLINK2"\`, \`"REGENIE"\`,
  \`"SAIGE"\`, \`"BOLT"\`, \`"METAL"\`, \`"PGS"\`, \`"SSF"\`,
  \`"GWAS-SSF"\`)

- columns_file:

  Optional two-column summary header mapping file

- q_score_thr:

  Optional comma-separated p-value thresholds (e.g.
  \`"1e-8,1e-6,1e-4"\`)

- summaries_list_file:

  Optional path to a file (one summary path per line) or directory of
  summary files, matching upstream \`bcftools +score –summaries\`.

- log_path:

  Optional path for a matching/audit log with loaded, matched,
  allele-mismatch, and duplicate-marker counts per PRS.

- use_variant_id:

  Logical; if TRUE, match variants by ID instead of CHR+BP

- counts:

  Logical; if TRUE, include per-threshold matched-variant counts

- samples:

  Optional comma-separated list of sample names to subset (e.g.
  \`"SAMP1,SAMP2"\`)

- force_samples:

  Logical; if TRUE, ignore missing samples instead of erroring

- regions:

  Optional comma-separated region list (e.g. \`"1:1000-2000,2:50-90"\`)

- regions_file:

  Optional path to a regions file

- regions_overlap:

  Overlap mode for regions (\`0\`, \`1\`, or \`2\`). Default 1 (trim to
  region).

- targets:

  Optional comma-separated targets list

- targets_file:

  Optional path to a targets file

- targets_overlap:

  Overlap mode for targets (\`0\`, \`1\`, or \`2\`). Default 0 (record
  must start in region).

- apply_filters:

  Optional comma-separated FILTER names to keep (e.g. \`"PASS,."\`)

- include:

  Optional site expression (currently unsupported)

- exclude:

  Optional site expression (currently unsupported)

## Value

A data frame with one row per sample and score/count columns.
