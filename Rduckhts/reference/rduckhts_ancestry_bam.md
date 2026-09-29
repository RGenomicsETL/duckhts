# Project ancestry directly from an indexed BAM or CRAM

The panel must be a committed, versioned canonical panel intersected
with the reference products. It has Somalier panel columns assembly,
site_index, region, position, allele_a, and allele_b. Each count is
extracted by the DuckHTS panel reader on the caller's connection. The
count relation and frequency view are discarded on return; no private
connection is created.

## Usage

``` r
rduckhts_ancestry_bam(
  con,
  source_path,
  sample_id,
  reference_path,
  panel_table,
  reference_table,
  loadings_table,
  correction_table = NULL,
  frequency_method = c("allele_fraction", "called_genotype"),
  min_depth = 7,
  min_cor = 0.4,
  ...
)
```

## Arguments

- con:

  Connection with DuckHTS loaded.

- source_path:

  Indexed BAM/CRAM path.

- sample_id:

  Sample identifier in the output.

- reference_path:

  FASTA path for reference checks and CRAM decoding.

- panel_table:

  Committed panel relation with canonical site identity.

- reference_table:

  Keyed wide or long reference relation.

- loadings_table:

  PC loadings for long references; pass the correction relation here for
  wide references.

- correction_table:

  PC correction relation for long references; NULL for wide.

- frequency_method:

  \`"allele_fraction"\` for B/(A+B) or \`"called_genotype"\` for the
  Somalier balance-rule genotype divided by two.

- min_depth:

  Minimum measured A+B depth to retain a site.

- min_cor:

  Minimum predicted-frequency correlation, default 0.4.

- ...:

  Additional arguments to \`rduckhts_somalier_bam_counts\`, such as
  index_path, map/base quality thresholds and worker_count.

## Value

Proportions and matching audit with frequency_method and min_depth.
