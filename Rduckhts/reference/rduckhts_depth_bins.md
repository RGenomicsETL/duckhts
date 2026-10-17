# Correct Bin Read Counts for GC Content

Gives each bin of each sample an expected read count from the median
read rate of the sample's autosomal training bins in the same GC
stratum, so that \`ratio = reads / expected\` is about 1 for two copies,
0.5 for one and 1.5 for three. The output feeds
\[rduckhts_depth_region_fraction()\] and
\[rduckhts_depth_sex_fraction()\]. Normalization uses eligible autosomal
bins with positive reads; empty bins retain zero ratios when their bin
and stratum rules pass.

## Usage

``` r
rduckhts_depth_bins(
  con,
  bins_table,
  min_mappability = 0.9,
  gc_step = 0.01,
  min_stratum_bins = 10,
  outlier_quantile = 0.01,
  max_ratio = 4,
  table_name = NULL,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded.

- bins_table:

  Name of a table or view with \`sample_id\`, \`chrom\`, \`start\`,
  \`end\` (zero-based, half-open), \`reads\` and \`gc\` (the G+C share
  of the A, C, G and T reference bases, from 0 to 1), and optionally
  \`mappability\` (0 to 1).

- min_mappability:

  Bins with a lower mappability are not kept.

- gc_step:

  Width of a GC stratum, in \\1/2147483647, 1\\, so the stratum number
  fits the integer output column.

- min_stratum_bins:

  Training bins a stratum needs, a whole number.

- outlier_quantile:

  Training bins outside this quantile of the read rate, and its
  complement, set no stratum median; in \\0, 0.5).

- max_ratio:

  Bins with a higher ratio are not kept; above 1.

- table_name:

  Optional table to create instead of returning a data frame.

- overwrite:

  Whether \`table_name\` may replace an existing table.

## Value

A data frame (or invisible \`TRUE\` when \`table_name\` is given) with
one row per input row: \`sample_id\`, \`chrom\`, \`start\`, \`end\`,
\`reads\`, \`gc\`, \`mappability\`, \`gc_stratum\`, \`stratum_bins\`,
\`expected\`, \`ratio\`, \`log2_ratio\`, \`reason\`, \`kept\` and
\`method\`.

## Details

The GC stratum of a bin is \`floor(gc / gc_step)\`. Training bins are
autosomal bins that pass the bin rules, hold at least one read and have
a read rate inside the \`outlier_quantile\` and \`1 - outlier_quantile\`
quantiles of the sample. A bin is kept unless it lacks a count or a GC
value, has a mappability below \`min_mappability\` (when the relation
has a \`mappability\` column), falls in a stratum with fewer than
\`min_stratum_bins\` training bins, or has a ratio above \`max_ratio\`;
\`reason\` names the rule. Only \`max_ratio\` looks at a bin's count, so
a bin with no reads stays kept. Autosomal means a contig whose
\`duckhts_contig_key()\` is not X, Y or MT. Restrict contigs and remove
problem bins in SQL first.

## See also

\[rduckhts_bam_bin_counts()\], \[rduckhts_fasta_nuc()\]

## Examples

``` r
con <- rduckhts_connect()
DBI::dbExecute(con, paste(
  "CREATE TEMP TABLE bins AS SELECT 'S1' AS sample_id, CAST(c AS VARCHAR) AS chrom,",
  "b * 1000000 AS start, (b + 1) * 1000000 AS \"end\",",
  "0.35 + 0.25 * ((b * 37 + c * 11) % 100) / 100.0 AS gc,",
  "CAST(600 + 4 * ((b * 13 + c * 7) % 11) AS BIGINT) AS reads",
  "FROM range(1, 23) t(c), range(0, 100) r(b)"
))
#> [1] 2200
head(rduckhts_depth_bins(con, "bins")[, c("chrom", "start", "reads", "ratio", "kept")])
#>   chrom   start reads     ratio kept
#> 1    21 4.8e+07   604 0.9741935 TRUE
#> 2    21 4.9e+07   612 0.9870968 TRUE
#> 3    21 5.0e+07   620 1.0000000 TRUE
#> 4    21 5.1e+07   628 1.0129032 TRUE
#> 5    21 5.2e+07   636 1.0225080 TRUE
#> 6    21 5.3e+07   600 0.9677419 TRUE
DBI::dbDisconnect(con, shutdown = TRUE)
```
