# Estimate the Share of Cells with a One-Copy Change of a Region

For each sample and named region, compares the depth of the region's
kept bins with the depth of the sample's other kept autosomal bins:
\`ratio = (C_r / E_r) / (C_b / E_b)\` with \`C\` and \`E\` the sums of
reads and expected counts, and \`shift = 2 \* (ratio - 1)\`, the share
of cells with a one-copy gain (positive) or loss (negative) relative to
two copies when the baseline is unchanged. A one-copy change in every
cell gives 1 or -1.

## Usage

``` r
rduckhts_depth_region_fraction(
  con,
  bins_table,
  regions_table,
  reference_table = NULL,
  table_name = NULL,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded.

- bins_table:

  Name of a table or view of corrected bins with \`sample_id\`,
  \`chrom\`, \`start\`, \`end\`, \`reads\`, \`expected\` and \`kept\`,
  such as the output of \[rduckhts_depth_bins()\].

- regions_table:

  Name of a table or view with \`region\`, \`chrom\`, \`start\` and
  \`end\` (zero-based, half-open); rows with one name form one region,
  and a bin belongs to it when its midpoint lies in one of its
  intervals.

- reference_table:

  Optional name of a table or view of corrected bins of reference
  samples, in the layout of \`bins_table\`.

- table_name:

  Optional table to create instead of returning a data frame.

- overwrite:

  Whether \`table_name\` may replace an existing table.

## Value

A data frame (or invisible \`TRUE\` when \`table_name\` is given) with
one row per sample and region: \`sample_id\`, \`region\`, \`bins\`,
\`baseline_bins\`, \`ratio\`, \`shift\`, \`shift_se\`, \`shift_low\`,
\`shift_high\`, \`dispersion\`, \`reference_shift\`,
\`reference_samples\`, \`status\` (\`ok\`, \`no_bins\`, \`no_baseline\`
or \`few_reference_samples\`) and \`method\`.

## Details

Without \`reference_table\`, the interval treats bin counts as
independent and quasi-Poisson, with the dispersion of the baseline bins,
and is \`shift +- 1.959964 shift_se\`. With \`reference_table\`,
corrected bins of reference samples, the reported \`shift\` is the
sample's shift minus the mean shift of the references (leaving out the
tested sample if it is a reference), and the interval uses Student's t
with \`reference_samples - 1\` degrees of freedom. Other changes in the
sample move the baseline.

## Examples

``` r
con <- rduckhts_connect()
DBI::dbExecute(con, paste(
  "CREATE TEMP TABLE bins AS SELECT 'S1' AS sample_id, CAST(c AS VARCHAR) AS chrom,",
  "b * 1000000 AS start, (b + 1) * 1000000 AS \"end\",",
  "0.35 + 0.25 * ((b * 37 + c * 11) % 100) / 100.0 AS gc,",
  "CAST((600 + 4 * ((b * 13 + c * 7) % 11)) *",
  "CASE WHEN c = 8 THEN 1.1 ELSE 1 END AS BIGINT) AS reads",
  "FROM range(1, 23) t(c), range(0, 100) r(b)"
))
#> [1] 2200
rduckhts_depth_bins(con, "bins", table_name = "corrected")
DBI::dbExecute(con, paste(
  "CREATE TEMP TABLE chromosomes AS SELECT CAST(c AS VARCHAR) AS region,",
  "CAST(c AS VARCHAR) AS chrom, 0 AS start, 100000000 AS \"end\" FROM range(1, 23) t(c)"
))
#> [1] 22
fractions <- rduckhts_depth_region_fraction(con, "corrected", "chromosomes")
fractions[fractions$region == "8", c("region", "shift", "shift_low", "shift_high")]
#>   region     shift shift_low shift_high
#> 1      8 0.1993923 0.1824568  0.2163279
DBI::dbDisconnect(con, shutdown = TRUE)
```
