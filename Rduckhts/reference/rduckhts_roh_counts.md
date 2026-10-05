# Find Runs of Homozygosity from Read Counts

Detect runs of homozygosity from allele read counts, for example counts
at panel sites from a BAM or CRAM, with the two-state model and
transitions of \[rduckhts_roh()\]. Genotype likelihoods come from a
binomial read model with a per-read sequencing error and an optional
contamination fraction instead of FORMAT/PL or GT. This emission is a
DuckHTS extension; \`bcftools roh\` has no read-count mode.

## Usage

``` r
rduckhts_roh_counts(
  con,
  counts_table,
  seq_error = 0.001,
  contamination = 0,
  genetic_map = NULL,
  hw_to_az = 6.7e-08,
  az_to_hw = 5e-09,
  rec_rate = NULL,
  table_name = NULL,
  overwrite = FALSE,
  max_sites = 2e+07,
  max_site_bytes = 2^32
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded.

- counts_table:

  Name of a table or view of read counts (see Details).

- seq_error:

  Per-read probability of showing the other allele, in \`(0, 0.5)\`.

- contamination:

  Fraction of reads from a contaminating individual, in \`\[0, 1)\`.

- genetic_map:

  Optional name of a genetic-map table or view.

- hw_to_az, az_to_hw:

  Per-base-pair transition probabilities Hardy-Weinberg to autozygous
  (\`-a\`) and back (\`-H\`), each in \`\[0, 1\]\`.

- rec_rate:

  Optional constant recombination rate per base pair (\`-M\`).

- table_name:

  Optional output table. \`NULL\` returns a data frame.

- overwrite:

  Whether an existing output table may be replaced.

- max_sites:

  Most sites one sample and chromosome may hold, from 1 to 100,000,000.
  A larger group is an error.

- max_site_bytes:

  Bound on the bytes the site buffers hold at once (16 bytes per site,
  24 for read counts). A decode does not grow a buffer when the bytes
  held by all ROH decodes in the process would pass its own
  \`max_site_bytes\`; decodes that run at the same time with different
  values each apply their own. Exceeding it is an error: decode fewer
  samples per call or raise it.

## Value

A data frame (or invisible \`TRUE\` when \`table_name\` is given) with
the columns of \[rduckhts_roh()\]: one row per run, unordered.

## Details

\`counts_table\` names a table or view with columns \`sample_id\`,
\`chrom\`, \`pos\` (one-based), \`ref_count\`, \`alt_count\` and \`af\`,
one row per sample and site. \`alt_count\` counts reads that show the
allele whose population frequency is \`af\`; \`ref_count\` counts reads
that show the other allele. The frequency can come from any source
joined in SQL beforehand: an INFO tag, a frequency table or
ancestry-weighted frequencies. Somalier site counts map as \`ref_count =
a\`, \`alt_count = b\`, with \`af\` the frequency of \`allele_b\`.

Each read shows the counted allele with probability \`(1 -
contamination) \* q + contamination \* c\`, where \`q\` is
\`seq_error\`, \`1/2\` or \`1 - seq_error\` for zero, one or two copies,
and \`c = af \* (1 - seq_error) + (1 - af) \* seq_error\` is the chance
that a read from a contaminating individual of the same population shows
it. Without the contamination term, contaminant reads at homozygous
sites look like heterozygous evidence, which shortens runs and lowers
FROH.

Sites with no reads, a NULL count, or an \`af\` that is NULL, NaN or 0
are skipped; a repeated position keeps one row. Negative counts are an
error. The model assumes biallelic sites, independent reads (no
duplicate or strand modelling) and one error rate for all reads. A
heterozygote is assumed to show each allele in half its reads;
reference-biased mapping makes the other allele rarer, more so where
flanking heterozygosity is high, so such heterozygotes can look
homozygous and lengthen runs in divergent regions. The contamination
fraction is supplied by the caller, for example from a separate
estimate; it is not estimated here. A contaminant from a different
population than \`af\` describes is approximated by \`af\`.

\`seq_error\` is an error of one read, and the reads of a site multiply.
A site with balanced reads excludes both homozygous genotypes whatever
\`seq_error\` is. The \`gt_error\` of \[rduckhts_roh()\] is an error of
the site, so a run can continue through one heterozygous call. The
read-count decode therefore corresponds to the genotype decode with no
tolerated genotype error, and reports fewer bases in runs than
\`rduckhts_roh(gt_error = 30)\` on the same sample.

## Examples

``` r
con <- rduckhts_connect()
DBI::dbExecute(con, paste(
  "CREATE TEMP TABLE counts AS SELECT 'S1' AS sample_id, '1' AS chrom,",
  "i * 25000 AS pos, 0.5 AS af,",
  "CASE WHEN i BETWEEN 40 AND 120 OR i % 2 = 0 THEN 30 ELSE 15 END AS ref_count,",
  "CASE WHEN i BETWEEN 40 AND 120 OR i % 2 = 0 THEN 0 ELSE 15 END AS alt_count",
  "FROM range(1, 161) AS t(i)"
))
#> [1] 160
rduckhts_roh_counts(con, "counts")
#>   sample chrom start   end  length n_markers  quality
#> 1     S1     1 1e+06 3e+06 2000001        81 47.67115
DBI::dbDisconnect(con, shutdown = TRUE)
```
