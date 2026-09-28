Real-input SQL scaling audit: measured paths and open coverage
================

<!-- benchmark_sql_scaling_audit.md is rendered from this file. -->

## Measurement contract

Source: `a0e7949dae4e7958d87cb9954f3cfa8fce59943a` (candidate);
`0b9f32a8f572c4b01f2346f00086602ece473c67` (VCF-count baseline only).
DuckDB R 1.5.5; Linux; one fresh DuckDB process per case/size/thread
cell, one observation per cell. `seconds` measures the query, including
source reading and temporary-table materialization, but not source VCF
subset derivation, panel preparation, connection, or COPY used for keyed
differential. RSS is process `VmHWM` (includes preparation); DuckDB
buffer is the query JSON profile’s `system_peak_buffer_memory`. Numbers
are MiB. These are exploratory measurements, not confidence intervals or
a cross-machine performance comparison.

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
baseline JSON EXPLAIN ANALYZE shows `READ_GENO` (1,236,049 rows) under a
materialized CTE feeding `HASH_JOIN`. Removing only `AS MATERIALIZED`
from `__dht_source_records` lets the scan feed the join directly. At one
thread, before → after: **2.327 → 2.198 s**, **749 → 430 MiB peak RSS**,
**1,367 → 612 MiB DuckDB buffer**. Across 319,349 / 646,810 / 1,236,049
records, before → after buffers (MiB) are **373 / 728 / 1,367 → 178 /
333 / 612**. Whole-result keyed differential on
`(sample_id, site_index)` for the 4× input (full results from the
0b9f32a8 and candidate builds, exported separately to Parquet): 2,000
baseline and 2,000 candidate rows, zero duplicate keys; `EXCEPT ALL` in
both directions returns zero rows.
`check_sql_scaling_vcf_differential.R` reproduces the assertion from the
two Parquet files; the checked-in receipt is
`sql_scaling_vcf_keyed_differential.tsv`. This is an allocation
improvement, not a claim of bounded memory: the remaining join still
scales with input, and the 1× query is below the timing floor.

## Site import: real-input materialization

The candidate build is the commit containing this report, based on
`509190790a3ce9d51ec735b2313d0c46b0969a70`. The public phase-3 GRCh37
chr22 VCF supplies 124,658 / 274,567 / 562,061 biallelic SNVs; the
importer emits exactly one row per site. Offline registry validation and
VCF derivation (`stage_sql_scaling_audit.R`) are outside query timing.
The full 562,061-row keyed differential on `(assembly, site_index)` is
in `sql_scaling_import_keyed_differential.tsv`; zero duplicate keys,
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
|     124658 |      124658 |       1 | 0.740 / 0.738 / 0.744 | 219 / 220 / 220 | 324           | 149 / 149 / 149 | 0 / 0 / 0 |
|     124658 |      124658 |       4 | 0.726 / 0.738 / 0.734 | 222 / 224 / 223 | 324           | 155 / 155 / 155 | 0 / 0 / 0 |
|     274567 |      274567 |       1 | 1.614 / 1.621 / 1.630 | 323 / 323 / 323 | 407           | 309 / 309 / 309 | 0 / 0 / 0 |
|     274567 |      274567 |       4 | 1.544 / 1.509 / 1.550 | 332 / 332 / 332 | 407           | 342 / 346 / 342 | 0 / 0 / 0 |
|     562061 |      562061 |       1 | 3.421 / 3.387 / 3.371 | 521 / 520 / 521 | 565           | 614 / 614 / 614 | 0 / 0 / 0 |
|     562061 |      562061 |       4 | 3.120 / 3.097 / 3.146 | 537 / 539 / 538 | 565           | 659 / 660 / 660 | 0 / 0 / 0 |

At one thread, baseline → candidate buffer for 1× / 2× / 4× is **222 /
469 / 942 → 149 / 309 / 614 MiB**. At 562,061 rows, the paired
keyed-export runs show **687 → 527 MiB RSS**; both builds produced all
562,061 rows. The repeated candidate 4-thread 4× RSS is **537–539 MiB**
against the **565 MiB** budget, with **0 MiB** spill. Reproduce query
repetitions with
`Rscript benchmarks/benchmark_sql_scaling_import_repetitions.R`; compare
separately exported results with
`Rscript benchmarks/check_sql_scaling_import_differential.R before.parquet after.parquet`.
Timing excludes input staging, keyed exports, and result-count probes.

## Open findings, ranked by measured exposure

1.  **VCF-count join:** 612 MiB of DuckDB buffer for 1.24 million GIAB
    records and 2,000 output rows. `READ_GENO` is estimated at one row;
    the left join therefore builds a hash table from the full VCF
    source. An exploratory panel `SEMI JOIN` was planned as
    `RIGHT_SEMI`, still building the source (658 MiB buffer in a
    single-process diagnostic), so it was not retained. A genuinely
    selective indexed input or keyed stream join needs a separate design
    and keyed differential; this diagnostic is not a timing verdict.
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
6.  **Not measured:** executed BAM/VCF/GFF `*_convert_parquet_sql` COPYs
    and DuckVEP. Both require physical multi-size public inputs, output
    denominator checks and independent disk-space limits.

The first-pass synthetic and 10k–40k-row rows have no extrapolatable
scaling conclusion. This report does not claim completion of the
unmeasured entry points.
