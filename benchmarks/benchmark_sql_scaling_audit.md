Real-input SQL scaling audit: measured paths and open coverage
================

<!-- benchmark_sql_scaling_audit.md is rendered from this file. -->

## Measurement contract

Source: `a0e7949dae4e7958d87cb9954f3cfa8fce59943a` (candidate);
`0b9f32a8f572c4b01f2346f00086602ece473c67` (legacy VCF-count baseline
only). DuckDB R 1.5.5; Linux; these previously recorded broad-run cells
use one fresh DuckDB process and one observation per case/size/thread.
The capped VCF-count matrix below uses three fresh processes per cell.
`seconds` measures the query, including source reading and
temporary-table materialization, but not source VCF subset derivation,
panel preparation, connection, or COPY used for keyed differential. RSS
is process `VmHWM` (includes preparation); DuckDB buffer is the query
JSON profile’s `system_peak_buffer_memory`. Numbers are MiB. These
measurements are not cross-machine performance comparisons.

Source identity is checked against
`r/duckhtsbench/inst/benchmark_registry.tsv`.
`stage_sql_scaling_audit.R` derives real chromosome prefixes of the GIAB
HG002 NIST v4.2.1 GRCh38 phased VCF (319,349 / 646,810 / 1,236,049
records), and unique biallelic SNV windows of the 1000 Genomes phase-3
GRCh37 chr22 VCF (124,658 / 274,567 / 562,061 records). VCF counts use
2,000 panel sites selected from each corresponding public VCF; import
sites return one row per staged site. GRCh37 epilepsy METAL summary
statistics supply genuine physical prefixes of 200,000 / 400,000 /
800,000 rows via a streaming `LIMIT`, with the registered GRCh37 FASTA
supplied to both munge macros. Liftover uses 250,000 / 500,000 /
1,000,000 real phase-3 chr22 variants and the registered GRCh37→GRCh38
chain and FASTAs. Exponent is `log(t4/t1)/log(input4/input1)` on the
displayed actual row counts. Source and output row counts are separate
denominators.

Reproduce from repository root with `bcftools` and the registered inputs
already in the cache:

``` sh
Rscript benchmarks/stage_sql_scaling_audit.R
make release -j2
Rscript benchmarks/benchmark_sql_scaling_real_run.R build/release/duckhts.duckdb_extension
repo=$PWD
mkdir -p /tmp/duckhts-scaling-r-lib
(cd r/Rduckhts && Rscript bootstrap.R "$repo" && R_LIBS_USER=/tmp/duckhts-scaling-r-lib THREADS=4 make test)
R_LIBS_USER=/tmp/duckhts-scaling-r-lib Rscript benchmarks/benchmark_sql_scaling_real_run.R build/release/duckhts.duckdb_extension --wrappers
Rscript -e 'rmarkdown::render("benchmarks/benchmark_sql_scaling_audit.Rmd", knit_root_dir = getwd(), quiet = TRUE)'
```

Pass `--download` to the staging script to fetch missing direct-download
inputs through the registry. The `0b9f32a8` VCF-count baseline is in
`sql_scaling_audit_vcf_before.tsv`; `--wrappers` writes
`sql_scaling_audit_wrappers.tsv`. The synthetic diagnostic inputs in
`sql_scaling_audit_measurements.tsv`,
`sql_scaling_audit_fixed_pairs.tsv`, and `sql_scaling_audit_norm.tsv`
fall below the timing floor and are excluded from scaling conclusions.

## Earlier above-floor probes (one observation per cell)

Columns give 1× / 2× / 4× observations. The single-thread 1× workloads
take 6–8 seconds. Both SQL macros and their R wrappers use
reference-backed allele validation and emit all selected rows. The
wrappers fetch the full result into R; their peak RSS includes the
returned data frame, while the SQL timings materialize a DuckDB
temporary table. The displayed exponents are exploratory: one
observation per cell does not meet the three-repetition timing evidence
requirement. These are distinct output-materialization workloads, not a
wrapper-overhead comparison.

| entry_point   | threads | input_records               | output_rows                 | seconds                 | peak_rss_mib    | peak_buffer_mib | exponent |
|:--------------|--------:|:----------------------------|:----------------------------|:------------------------|:----------------|:----------------|:---------|
| munge         |       1 | 200,000 / 400,000 / 800,000 | 200,000 / 400,000 / 800,000 | 6.665 / 12.745 / 25.161 | 353 / 421 / 563 | 297 / 457 / 779 | 0.96     |
| munge         |       4 | 200,000 / 400,000 / 800,000 | 200,000 / 400,000 / 800,000 | 6.836 / 12.765 / 25.378 | 338 / 413 / 675 | 318 / 468 / 799 | 0.95     |
| munge_metal   |       1 | 200,000 / 400,000 / 800,000 | 200,000 / 400,000 / 800,000 | 6.800 / 12.897 / 25.576 | 361 / 421 / 564 | 298 / 458 / 780 | 0.96     |
| munge_metal   |       4 | 200,000 / 400,000 / 800,000 | 200,000 / 400,000 / 800,000 | 6.944 / 13.123 / 26.162 | 344 / 409 / 707 | 319 / 469 / 800 | 0.96     |
| r_munge       |       1 | 200,000 / 400,000 / 800,000 | 200,000 / 400,000 / 800,000 | 7.360 / 14.039 / 26.606 | 390 / 454 / 645 | 264 / 278 / 557 | 0.93     |
| r_munge       |       4 | 200,000 / 400,000 / 800,000 | 200,000 / 400,000 / 800,000 | 7.229 / 14.472 / 26.123 | 359 / 462 / 653 | 278 / 278 / 556 | 0.93     |
| r_munge_metal |       1 | 200,000 / 400,000 / 800,000 | 200,000 / 400,000 / 800,000 | 7.303 / 13.761 / 24.817 | 402 / 488 / 710 | 281 / 281 / 562 | 0.88     |
| r_munge_metal |       4 | 200,000 / 400,000 / 800,000 | 200,000 / 400,000 / 800,000 | 7.212 / 13.300 / 26.388 | 372 / 489 / 714 | 279 / 279 / 559 | 0.94     |

## Real-input probes below the timing floor

These probes identify allocation behavior but their time exponents are
**not scaling verdicts**: all one-thread 1× queries finish in less than
two seconds. The liftover denominator is capped by the single staged
chr22 VCF; the count and site-import prefixes do not reach a 5-second 1×
while retaining a 4× prefix of the same public source.

| entry_point  | threads | input_records                 | output_rows                   | seconds               | peak_rss_mib    | peak_buffer_mib | exponent |
|:-------------|--------:|:------------------------------|:------------------------------|:----------------------|:----------------|:----------------|:---------|
| vcf_counts   |       1 | 319,349 / 646,810 / 1,236,049 | 2,000 / 2,000 / 2,000         | 0.592 / 1.166 / 2.198 | 210 / 289 / 430 | 178 / 333 / 612 | 0.97     |
| vcf_counts   |       4 | 319,349 / 646,810 / 1,236,049 | 2,000 / 2,000 / 2,000         | 0.585 / 1.158 / 2.197 | 220 / 286 / 427 | 186 / 341 / 620 | 0.98     |
| import_sites |       1 | 124,658 / 274,567 / 562,061   | 124,658 / 274,567 / 562,061   | 0.760 / 1.669 / 3.466 | 254 / 399 / 677 | 222 / 469 / 942 | 1.01     |
| import_sites |       4 | 124,658 / 274,567 / 562,061   | 124,658 / 274,567 / 562,061   | 0.736 / 1.553 / 3.156 | 255 / 407 / 693 | 228 / 506 / 986 | 0.97     |
| liftover     |       1 | 250,000 / 500,000 / 1,000,000 | 250,000 / 500,000 / 1,000,000 | 1.774 / 3.442 / 6.919 | 247 / 348 / 550 | 176 / 331 / 647 | 0.98     |
| liftover     |       4 | 250,000 / 500,000 / 1,000,000 | 250,000 / 500,000 / 1,000,000 | 1.431 / 2.815 / 5.622 | 281 / 399 / 616 | 192 / 378 / 690 | 0.99     |

## VCF counts: a materialized scan worth removing

At 1,236,049 physical GIAB VCF records and 2,000 output panel rows, the
legacy baseline JSON EXPLAIN ANALYZE shows `READ_GENO` (1,236,049 rows)
under a materialized CTE feeding `HASH_JOIN`. Removing only
`AS MATERIALIZED` from `__dht_source_records` lets the scan feed the
join directly. At one thread, in that earlier uncapped one-observation
run, before → after: **2.327 → 2.198 s**, **749 → 430 MiB peak RSS**,
**1,367 → 612 MiB DuckDB buffer**. Across 319,349 / 646,810 / 1,236,049
records, before → after buffers (MiB) are **373 / 728 / 1,367 → 178 /
333 / 612**. Whole-result keyed differential on
`(sample_id, site_index)` for the 4× input (full results from the
0b9f32a8 and candidate builds, exported separately to Parquet): 2,000
baseline and 2,000 candidate rows, zero duplicate keys; `EXCEPT ALL` in
both directions returns zero rows.
`check_sql_scaling_vcf_differential.R` reproduces the assertion from the
two Parquet files; the checked-in receipt is
`sql_scaling_vcf_keyed_differential.tsv`. This earlier allocation
comparison is not a bounded-memory claim; the capped repeated matrix
below applies the declared 256 MiB budgets.

## duckhts_somalier_vcf_counts: capped before/after matrix

This matrix uses the frozen pre-change extension and the extension built
from the recorded branch revision. Every cell runs in three fresh serial
processes at 1 and 4 DuckDB threads, and the query result is fetched
with every output column consumed. The file-axis cells hold a 2,000-site
panel fixed across the registered 1× / 2× / 4× GIAB inputs; the
panel-axis cells hold the 4× file fixed and use nested 2,000 / 4,000 /
8,000-site panels. The joint-growth cell uses the 2× file and 4,000-site
panel. The duplicate-error cell reads one million copies of a single
registered chr1 biallelic SNP with its GT and AD calls and uses the
matching one-site panel; the staging test exercises this deterministic
transformation from a local VCF record without network access.

Before these measurements, the benchmark source declared **process RSS
≤256 MiB**, **query-wide DuckDB buffer ≤256 MiB**, **spill exactly 0
MiB**, and for each file doubling at fixed panel **both median RSS and
median query buffer increase by no more than 16 MiB**. Each child sets
`memory_limit='256MB'`, `max_temp_directory_size='0B'`, and a private
temporary directory. The raw run file records the per-process settings,
all input and panel checksums, output row counts, timings, memory
observations, versions, build identities, and the frozen/candidate macro
hashes. Any failed or unmeasured metric remains a budget
failure/incomplete result. The fastest 1× run is below the 5-second
timing floor, so no timing exponent or timing-scaling verdict is made.

### Frozen pre-change implementation

| cell                   | threads | input_records | panel_sites | output_rows | execution     | seconds_median_min_max   | rss_mib_median_min_max   | buffer_mib_median_min_max | spill_mib_median_min_max | budget_verdict                      |
|:-----------------------|--------:|:--------------|:------------|:------------|:--------------|:-------------------------|:-------------------------|:--------------------------|:-------------------------|:------------------------------------|
| panel_1x_fixed_file_4x |       1 | 1236049       | 2000        | 0           | out_of_memory | 1.356 \[ 1.348, 1.370 \] | 358.6 \[ 358.3, 359.4 \] | NA                        | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |
| panel_1x_fixed_file_4x |       4 | 1236049       | 2000        | 0           | out_of_memory | 1.303 \[ 1.296, 1.315 \] | 357.0 \[ 350.7, 357.8 \] | NA                        | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |
| panel_2x_fixed_file_4x |       1 | 1236049       | 4000        | 0           | out_of_memory | 1.352 \[ 1.349, 1.374 \] | 357.9 \[ 357.6, 358.4 \] | NA                        | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |
| panel_2x_fixed_file_4x |       4 | 1236049       | 4000        | 0           | out_of_memory | 1.268 \[ 1.260, 1.270 \] | 354.9 \[ 353.7, 355.7 \] | NA                        | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |
| panel_4x_fixed_file_4x |       1 | 1236049       | 8000        | 0           | out_of_memory | 1.344 \[ 1.339, 1.345 \] | 357.0 \[ 356.7, 357.3 \] | NA                        | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |
| panel_4x_fixed_file_4x |       4 | 1236049       | 8000        | 0           | out_of_memory | 1.229 \[ 1.220, 1.233 \] | 353.2 \[ 347.1, 354.2 \] | NA                        | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |
| file_1x_fixed_panel_1x |       1 | 319349        | 2000        | 2000        | success       | 0.637 \[ 0.630, 0.645 \] | 318.2 \[ 318.0, 318.4 \] | 390.5 \[ 390.5, 391.0 \]  | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |
| file_1x_fixed_panel_1x |       4 | 319349        | 2000        | 2000        | success       | 0.580 \[ 0.579, 0.583 \] | 320.7 \[ 317.5, 322.4 \] | 413.5 \[ 412.1, 427.6 \]  | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |
| file_2x_fixed_panel_1x |       1 | 646810        | 2000        | 0           | out_of_memory | 1.168 \[ 1.166, 1.171 \] | 353.6 \[ 353.2, 353.7 \] | NA                        | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |
| file_2x_fixed_panel_1x |       4 | 646810        | 2000        | 0           | out_of_memory | 1.127 \[ 1.117, 1.144 \] | 331.9 \[ 330.8, 333.7 \] | NA                        | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |
| file_4x_fixed_panel_1x |       1 | 1236049       | 2000        | 0           | out_of_memory | 1.360 \[ 1.359, 1.361 \] | 358.2 \[ 358.1, 359.0 \] | NA                        | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |
| file_4x_fixed_panel_1x |       4 | 1236049       | 2000        | 0           | out_of_memory | 1.292 \[ 1.272, 1.299 \] | 356.1 \[ 352.9, 357.9 \] | NA                        | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |
| joint_file_2x_panel_2x |       1 | 646810        | 4000        | 0           | out_of_memory | 1.159 \[ 1.153, 1.160 \] | 353.8 \[ 353.3, 354.2 \] | NA                        | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |
| joint_file_2x_panel_2x |       4 | 646810        | 4000        | 0           | out_of_memory | 1.147 \[ 1.145, 1.148 \] | 336.5 \[ 334.5, 336.7 \] | NA                        | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |
| duplicate_error_1m     |       1 | 1000000       | 1           | 0           | out_of_memory | 0.719 \[ 0.714, 0.732 \] | 374.0 \[ 373.8, 374.3 \] | NA                        | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |
| duplicate_error_1m     |       4 | 1000000       | 1           | 0           | out_of_memory | 0.702 \[ 0.695, 0.711 \] | 373.4 \[ 367.4, 374.4 \] | NA                        | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |

### Branch implementation

| cell                   | threads | input_records | panel_sites | output_rows | execution       | seconds_median_min_max   | rss_mib_median_min_max   | buffer_mib_median_min_max | spill_mib_median_min_max | budget_verdict                      |
|:-----------------------|--------:|:--------------|:------------|:------------|:----------------|:-------------------------|:-------------------------|:--------------------------|:-------------------------|:------------------------------------|
| panel_1x_fixed_file_4x |       1 | 1236049       | 2000        | 2000        | success         | 1.962 \[ 1.962, 1.976 \] | 178.5 \[ 176.2, 179.3 \] | 59.7 \[ 59.6, 59.8 \]     | 0.0 \[ 0.0, 0.0 \]       | PASS                                |
| panel_1x_fixed_file_4x |       4 | 1236049       | 2000        | 2000        | success         | 1.965 \[ 1.948, 1.965 \] | 182.8 \[ 182.1, 184.3 \] | 121.6 \[ 117.4, 124.9 \]  | 0.0 \[ 0.0, 0.0 \]       | PASS                                |
| panel_2x_fixed_file_4x |       1 | 1236049       | 4000        | 4000        | success         | 1.972 \[ 1.968, 1.972 \] | 180.9 \[ 179.4, 181.3 \] | 67.9 \[ 67.9, 67.9 \]     | 0.0 \[ 0.0, 0.0 \]       | PASS                                |
| panel_2x_fixed_file_4x |       4 | 1236049       | 4000        | 4000        | success         | 1.973 \[ 1.951, 2.004 \] | 191.6 \[ 191.4, 191.9 \] | 129.6 \[ 129.1, 129.8 \]  | 0.0 \[ 0.0, 0.0 \]       | PASS                                |
| panel_4x_fixed_file_4x |       1 | 1236049       | 8000        | 8000        | success         | 2.032 \[ 2.024, 2.038 \] | 186.1 \[ 185.5, 186.6 \] | 97.2 \[ 97.2, 97.2 \]     | 0.0 \[ 0.0, 0.0 \]       | PASS                                |
| panel_4x_fixed_file_4x |       4 | 1236049       | 8000        | 8000        | success         | 2.030 \[ 2.016, 2.035 \] | 205.0 \[ 204.3, 207.2 \] | 215.5 \[ 213.0, 222.8 \]  | 0.0 \[ 0.0, 0.0 \]       | PASS                                |
| file_1x_fixed_panel_1x |       1 | 319349        | 2000        | 2000        | success         | 0.540 \[ 0.529, 0.545 \] | 175.4 \[ 175.4, 176.0 \] | 59.7 \[ 59.6, 59.7 \]     | 0.0 \[ 0.0, 0.0 \]       | PASS                                |
| file_1x_fixed_panel_1x |       4 | 319349        | 2000        | 2000        | success         | 0.541 \[ 0.534, 0.543 \] | 180.1 \[ 179.0, 180.6 \] | 121.3 \[ 117.9, 122.7 \]  | 0.0 \[ 0.0, 0.0 \]       | PASS                                |
| file_2x_fixed_panel_1x |       1 | 646810        | 2000        | 2000        | success         | 1.055 \[ 1.049, 1.062 \] | 176.4 \[ 175.7, 176.8 \] | 59.7 \[ 59.6, 59.7 \]     | 0.0 \[ 0.0, 0.0 \]       | PASS                                |
| file_2x_fixed_panel_1x |       4 | 646810        | 2000        | 2000        | success         | 1.044 \[ 1.040, 1.067 \] | 179.6 \[ 178.5, 179.6 \] | 122.8 \[ 118.6, 123.9 \]  | 0.0 \[ 0.0, 0.0 \]       | PASS                                |
| file_4x_fixed_panel_1x |       1 | 1236049       | 2000        | 2000        | success         | 1.951 \[ 1.942, 1.956 \] | 178.4 \[ 177.5, 178.4 \] | 59.6 \[ 59.6, 59.7 \]     | 0.0 \[ 0.0, 0.0 \]       | PASS                                |
| file_4x_fixed_panel_1x |       4 | 1236049       | 2000        | 2000        | success         | 1.961 \[ 1.953, 1.967 \] | 183.1 \[ 183.0, 184.5 \] | 123.4 \[ 122.7, 123.8 \]  | 0.0 \[ 0.0, 0.0 \]       | PASS                                |
| joint_file_2x_panel_2x |       1 | 646810        | 4000        | 4000        | success         | 1.071 \[ 1.064, 1.081 \] | 178.5 \[ 178.5, 179.2 \] | 68.0 \[ 67.9, 68.0 \]     | 0.0 \[ 0.0, 0.0 \]       | PASS                                |
| joint_file_2x_panel_2x |       4 | 646810        | 4000        | 4000        | success         | 1.058 \[ 1.057, 1.060 \] | 189.7 \[ 189.0, 190.4 \] | 130.2 \[ 130.1, 130.2 \]  | 0.0 \[ 0.0, 0.0 \]       | PASS                                |
| duplicate_error_1m     |       1 | 1000000       | 1           | 0           | duplicate_error | 0.740 \[ 0.738, 0.742 \] | 161.3 \[ 161.0, 161.4 \] | NA                        | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |
| duplicate_error_1m     |       4 | 1000000       | 1           | 0           | duplicate_error | 0.759 \[ 0.754, 0.762 \] | 162.6 \[ 162.1, 163.8 \] | NA                        | 0.0 \[ 0.0, 0.0 \]       | FAIL/INCOMPLETE (0/3 within budget) |

### File-doubling memory gates

The `1×→2×` and `2×→4×` verdicts require both median memory deltas to
remain at or below 16 MiB. Failures and unavailable profiles remain
visible.

| implementation | threads | doubling | median_rss_delta_mib | median_buffer_delta_mib | verdict                      |
|:---------------|--------:|:---------|---------------------:|------------------------:|:-----------------------------|
| prechange      |       1 | 1x to 2x |           35.4414062 |                      NA | FAIL/INCOMPLETE (unmeasured) |
| prechange      |       1 | 2x to 4x |            4.5898438 |                      NA | FAIL/INCOMPLETE (unmeasured) |
| prechange      |       4 | 1x to 2x |           11.1289062 |                      NA | FAIL/INCOMPLETE (unmeasured) |
| prechange      |       4 | 2x to 4x |           24.2812500 |                      NA | FAIL/INCOMPLETE (unmeasured) |
| branch         |       1 | 1x to 2x |            0.9960938 |               0.0312500 | PASS                         |
| branch         |       1 | 2x to 4x |            1.9609375 |              -0.0546875 | PASS                         |
| branch         |       4 | 1x to 2x |           -0.5507812 |               1.5195312 | PASS                         |
| branch         |       4 | 2x to 4x |            3.4960938 |               0.5703125 | PASS                         |

### Per-cell output and provenance checks

| cell                   | threads | output_identity               |
|:-----------------------|--------:|:------------------------------|
| panel_1x_fixed_file_4x |       1 | comparison unavailable        |
| panel_1x_fixed_file_4x |       4 | comparison unavailable        |
| panel_2x_fixed_file_4x |       1 | comparison unavailable        |
| panel_2x_fixed_file_4x |       4 | comparison unavailable        |
| panel_4x_fixed_file_4x |       1 | comparison unavailable        |
| panel_4x_fixed_file_4x |       4 | comparison unavailable        |
| file_1x_fixed_panel_1x |       1 | equal serialized full results |
| file_1x_fixed_panel_1x |       4 | equal serialized full results |
| file_2x_fixed_panel_1x |       1 | comparison unavailable        |
| file_2x_fixed_panel_1x |       4 | comparison unavailable        |
| file_4x_fixed_panel_1x |       1 | comparison unavailable        |
| file_4x_fixed_panel_1x |       4 | comparison unavailable        |
| joint_file_2x_panel_2x |       1 | comparison unavailable        |
| joint_file_2x_panel_2x |       4 | comparison unavailable        |
| duplicate_error_1m     |       1 | comparison unavailable        |
| duplicate_error_1m     |       4 | comparison unavailable        |

The raw per-run records are `sql_scaling_vcf_counts_runs.tsv`. The
following tables identify the exact builds and staged inputs used by
those runs; the baseline source commit is
2c5b69b6d23fcd3cb7928465c1d886d23ac8a26b and its empty dirty-diff hash
identifies the clean commit tree.

| implementation | source_git_sha                           | source_dirty_diff_sha256                                         | extension_sha256                                                 | macro_text_sha256                                                | duckdb_version | duckdb_source_id | R_version                    | DBI_version | duckdb_R_version | Rduckhts_installed_version | Rduckhts_source_version | duckhtsbench_source_version | htslib_version |
|:---------------|:-----------------------------------------|:-----------------------------------------------------------------|:-----------------------------------------------------------------|:-----------------------------------------------------------------|:---------------|:-----------------|:-----------------------------|:------------|:-----------------|:---------------------------|:------------------------|:----------------------------|---------------:|
| prechange      | 2c5b69b6d23fcd3cb7928465c1d886d23ac8a26b | e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855 | f9f90d05bb5ae7b0265047652ca1e0184875aeb3535f89bb1bf08def578dde56 | cc5f8aeae2a6db3210a547cebc3fcf092da6d9cdeae82185c6fb3ecce5fbe557 | v1.5.5         | d8cdaa33fda      | R version 4.6.0 (2026-04-24) | 1.3.0       | 1.5.5            | 0.1.3.0.0.2                | 1.5.2.9007-0.1.5        | 0.0.0.9001                  |           1.24 |
| branch         | ceb2341ff9e746379955b5ea6462f05e9bd6fa8a | e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855 | 8a11aea8e552b3953dd1b10410b4670e0f72e6999fcbcd878cc04b26b10f0087 | ee2361335884564210aa91a7e8260720cc2d72333c56d8ddb153392ae9256cc8 | v1.5.5         | d8cdaa33fda      | R version 4.6.0 (2026-04-24) | 1.3.0       | 1.5.5            | 0.1.3.0.0.2                | 1.5.2.9007-0.1.5        | 0.0.0.9001                  |           1.24 |

| input_id                                 | input_sha256                                                     | input_bytes | input_records | panel_source_id     | panel_source_sha256                                              | panel_source_bytes | panel_sites | panel_sha256                                                     |
|:-----------------------------------------|:-----------------------------------------------------------------|------------:|--------------:|:--------------------|:-----------------------------------------------------------------|-------------------:|------------:|:-----------------------------------------------------------------|
| sql_scaling_giab_4x                      | 932617388d0dd702066a5bcd55fec3e443946bf2d93c1796ed82312a55a1560f |    47804032 |       1236049 | sql_scaling_giab_1x | a3f740ed18f64ceb129365c5ada886391cdc85fd3450f9f8f72941282d77f71b |           12430977 |        2000 | 3745830e874f19e927203bacfd2d09ac834b8c2d05055997e834750038ccace9 |
| sql_scaling_giab_4x                      | 932617388d0dd702066a5bcd55fec3e443946bf2d93c1796ed82312a55a1560f |    47804032 |       1236049 | sql_scaling_giab_1x | a3f740ed18f64ceb129365c5ada886391cdc85fd3450f9f8f72941282d77f71b |           12430977 |        4000 | 303abd2b4af8a7c46c145e1c02c030857fecb913554b234528d9f0ea727bbd7c |
| sql_scaling_giab_4x                      | 932617388d0dd702066a5bcd55fec3e443946bf2d93c1796ed82312a55a1560f |    47804032 |       1236049 | sql_scaling_giab_1x | a3f740ed18f64ceb129365c5ada886391cdc85fd3450f9f8f72941282d77f71b |           12430977 |        8000 | 66f69b2054964f223bb73eb125e717cb53db87722ab165aede74c3c6d336a5e0 |
| sql_scaling_giab_1x                      | a3f740ed18f64ceb129365c5ada886391cdc85fd3450f9f8f72941282d77f71b |    12430977 |        319349 | sql_scaling_giab_1x | a3f740ed18f64ceb129365c5ada886391cdc85fd3450f9f8f72941282d77f71b |           12430977 |        2000 | 3745830e874f19e927203bacfd2d09ac834b8c2d05055997e834750038ccace9 |
| sql_scaling_giab_2x                      | 7e525c659ec98de8994c7aaa55f5a662d7f6457625901578b77ef98f0b99e9f4 |    25107378 |        646810 | sql_scaling_giab_1x | a3f740ed18f64ceb129365c5ada886391cdc85fd3450f9f8f72941282d77f71b |           12430977 |        2000 | 3745830e874f19e927203bacfd2d09ac834b8c2d05055997e834750038ccace9 |
| sql_scaling_giab_2x                      | 7e525c659ec98de8994c7aaa55f5a662d7f6457625901578b77ef98f0b99e9f4 |    25107378 |        646810 | sql_scaling_giab_1x | a3f740ed18f64ceb129365c5ada886391cdc85fd3450f9f8f72941282d77f71b |           12430977 |        4000 | 303abd2b4af8a7c46c145e1c02c030857fecb913554b234528d9f0ea727bbd7c |
| sql_scaling_vcf_counts_duplicate_million | 1bdd840714bf8c18a8fd8a9f356530ffde82f45510bce0635c0c697a47721921 |    42006757 |       1000000 | sql_scaling_giab_1x | a3f740ed18f64ceb129365c5ada886391cdc85fd3450f9f8f72941282d77f71b |           12430977 |           1 | 4ffeb81ef45b97e60d208c4e4f4618b6b24ec274383e02c0bef46826460dab9c |

## Site import: real-input materialization

The keyed differential and the paired baseline → candidate figures below
come from candidate `a2e81c6f`, based on
`509190790a3ce9d51ec735b2313d0c46b0969a70`. The repetitions table was
measured again on the commit containing this report, after merging the
Somalier X/Y panel change; chr22 input has no X/Y records, so the output
rows are unchanged. The public phase-3 GRCh37 chr22 VCF supplies 124,658
/ 274,567 / 562,061 biallelic SNVs; the importer emits exactly one row
per site. Offline registry validation and VCF derivation
(`stage_sql_scaling_audit.R`) are outside query timing. The full
562,061-row keyed differential on `(assembly, site_index)` is in
`sql_scaling_import_keyed_differential.tsv`; zero duplicate keys,
missing rows, or column disagreements. Its baseline and candidate
Parquet exports were produced in separate fresh processes. The baseline
physical plan materializes `__dht_source` (562,061 records),
`__dht_autosomal`, and `__dht_oriented` before sorting; the candidate
retains only the shared source, allowing validation and orientation to
consume it without two further wide copies. `READ_BCF` runs once in
either plan.

Before measurement, the per-process budget was set at **256 MiB
runtime + 3 × (64 bytes decoded source + 128 bytes materialized output)
per physical input row**; query spill limit **0 MiB**. The three
fresh-process repetitions per size at each of 1 and 4 threads are in
`sql_scaling_import_repetitions.tsv`. Each measurement materializes
every output row in a temporary table. Peak RSS includes connection and
extension overhead; buffer/spill are taken from the measured query’s
JSON profile, before any row-count probe. Staging and oracle work are
reported separately. The 1× one-thread time is \<5 s, so **no timing
exponent or timing scaling verdict** is asserted. Memory stays under the
predetermined RSS budget with zero spill at both thread counts.

| input_rows | output_rows | threads | seconds               | rss_mib         | rss_limit_mib | buffer_mib      | spill_mib |
|-----------:|------------:|--------:|:----------------------|:----------------|:--------------|:----------------|:----------|
|     124658 |      124658 |       1 | 0.784 / 0.775 / 0.755 | 219 / 218 / 218 | 324           | 145 / 145 / 145 | 0 / 0 / 0 |
|     124658 |      124658 |       4 | 0.749 / 0.742 / 0.741 | 221 / 223 / 221 | 324           | 150 / 150 / 150 | 0 / 0 / 0 |
|     274567 |      274567 |       1 | 1.703 / 1.675 / 1.681 | 322 / 322 / 322 | 407           | 303 / 303 / 303 | 0 / 0 / 0 |
|     274567 |      274567 |       4 | 1.564 / 1.573 / 1.575 | 329 / 330 / 329 | 407           | 332 / 335 / 329 | 0 / 0 / 0 |
|     562061 |      562061 |       1 | 3.485 / 3.416 / 3.475 | 519 / 519 / 519 | 565           | 606 / 606 / 606 | 0 / 0 / 0 |
|     562061 |      562061 |       4 | 3.199 / 3.190 / 3.224 | 531 / 532 / 534 | 565           | 643 / 642 / 646 | 0 / 0 / 0 |

At one thread, baseline → candidate buffer for 1× / 2× / 4× is **222 /
469 / 942 → 149 / 309 / 614 MiB**. At 562,061 rows, the paired
keyed-export runs show **687 → 527 MiB RSS**; both builds produced all
562,061 rows. The repeated candidate 4-thread 4× RSS is **531–534 MiB**
against the **565 MiB** budget, with **0 MiB** spill. Reproduce query
repetitions with
`Rscript benchmarks/benchmark_sql_scaling_import_repetitions.R`; compare
separately exported results with
`Rscript benchmarks/check_sql_scaling_import_differential.R before.parquet after.parquet`.
Timing excludes input staging, keyed exports, and result-count probes.

## Open findings, ranked by measured exposure

1.  **VCF-count join (addressed on this branch):** the pre-change join
    built its hash table from every decoded VCF record (612 MiB of
    DuckDB buffer for 1.24 million GIAB records and 2,000 output rows).
    The branch keeps call payloads only for panel coordinates. In the
    capped matrix above, the branch completes 42 of 42 non-error runs
    within the 256 MiB budgets, and the pre-change build runs out of
    memory in 36 of 42. Timings stay below the 5-second floor, so this
    is a memory result, not a timing verdict.
2.  **Site import:** the output grows one-to-one with input; the
    candidate reduces its 4× one-thread buffer from 942 to 614 MiB, but
    the 1× query is still below the timing floor. A larger public source
    is required before a time-scaling verdict.
3.  **Munge/Metal reference-backed modes:** the single-observation 200k
    / 400k / 800k real epilepsy probe approaches linear time, but does
    not meet the three-repetition timing gate. The 800k-row SQL output
    occupies about 779 MiB DuckDB buffer; COPY-to-disk output and any
    bounded streaming requirement are separate questions. The R wrapper
    800k-row result reaches 645 MiB RSS (`r_munge`) and 710 MiB
    (`r_munge_metal`), consistent with returning a large R data frame.
4.  **Liftover:** the entire staged chr22 GRCh37 source holds only ~1.1
    million variants; the current 250k baseline is below the timing
    floor. A larger public source, allele/mapped-output verification,
    and EXPLAIN ANALYZE at that scale are required.
5.  **Not measured with real multi-sample evidence:** BAM/CRAM
    extraction, sketch preparation/verification, Charr, matched
    contamination and R all-pairs relatedness. The phased 1000 Genomes
    chr22 VCFs have GT but no FORMAT/AD. Public indexed 30× CRAM
    regional access was probed on HG00188 chr22:16,000,000–16,020,000
    (2,680 alignments; 1.75 s; 37 MiB RSS) using the registered GRCh38
    reference. This establishes feasibility, not a checksum-validated
    16/32/64-sample staged cohort or a Somalier performance measurement.
    All-pairs relatedness needs a separately capped quadratic
    output-size and memory contract.
6.  **Not measured:** executed BAM/VCF/GFF `*_convert_parquet_sql`
    COPYs. These require physical multi-size public inputs, output
    denominator checks and independent disk-space limits.

The first-pass synthetic and 10k–40k-row rows have no extrapolatable
scaling conclusion. This report does not claim completion of the
unmeasured entry points.
