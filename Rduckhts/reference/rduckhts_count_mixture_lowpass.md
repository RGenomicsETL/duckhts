# Estimate a Conditional Second-Genome Share from Low-Pass Counts

Uses first-read and discordant-pair moments with fixed external
calibration. This estimator is separate from
[`rduckhts_count_error_fit()`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/rduckhts_count_error_fit.md).
It is suitable for low-pass, cell-free DNA or low-input DNA counts only
under the supplied population frequencies and calibration; the share is
not an unconditional biological fraction. Matched-frequency unrelated
and IBD1 hypotheses can have identical moments and are reported as
unresolved.

## Usage

``` r
rduckhts_count_mixture_lowpass(
  con,
  counts_table,
  seq_error,
  artefact_weight,
  homozygosity_excess,
  block_bases = 1e+07,
  max_depth = 8,
  min_blocks = 30,
  max_blocks = 512,
  max_candidates = 16,
  max_state_bytes = 2^28,
  table_name = NULL,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded.

- counts_table:

  Table or view containing `sample_id`, `candidate`, `chrom`, one-based
  `pos`, `ref_count`, `alt_count`, `other_count`, `sample_af` and
  `candidate_af`. Counts and `sample_af` must agree across candidates;
  site keys must be unique within each sample and candidate. Frequencies
  refer to the counted alternate allele. Sites missing any candidate
  frequency are unavailable to every candidate of that sample.

- seq_error:

  Required fixed sequencing error in `[0, 0.5)`.

- artefact_weight:

  Required fixed uniform-count artefact weight in `[0, 1)`.

- homozygosity_excess:

  Required fixed excess of homozygosity in `[0, 1]`.

- block_bases:

  Bases per chromosome-local block, at least 1.

- max_depth:

  Maximum ref plus alt depth, 1 to 65535. Deeper sites are excluded;
  `other_count` contributes only to QC totals.

- min_blocks:

  Informative blocks required for 95 percent limits, at least 2,
  default 30. Matched frequencies require pair-bearing blocks; otherwise
  read-bearing blocks count. Lower values support small-fixture checks,
  not well-calibrated inference. Below this floor the share is reported
  but limits are NA and status is `insufficient_blocks`.

- max_blocks:

  Maximum physical blocks per sample, 1 to 4096.

- max_candidates:

  Maximum candidates per sample, 1 to 64.

- max_state_bytes:

  Maximum bytes of concurrently owned native moment buffers and fit
  copies. DuckDB owns SQL sorting, materialization and spill.

- table_name:

  Optional output table; otherwise return a data frame.

- overwrite:

  Whether to replace an existing output table.

## Value

One row per sample, candidate and relationship model. IBD1 rows require
exact `candidate_af == sample_af` at every input site, including
zero-read and excluded sites. Columns include the bounded share (0 to
0.5), block-jackknife standard error and score-inverted limits, singles
share, pair residual, relationship diagnostic, QC read/pair totals,
calibration, score rank and status. The discordance ratio divides
discordant pairs by
`SUM(choose(depth, 2) * sample_af * (1 - sample_af))`. When `table_name`
is supplied, returns invisible `TRUE`.
