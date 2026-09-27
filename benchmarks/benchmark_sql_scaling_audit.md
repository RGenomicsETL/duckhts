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

## Qualified 1× workloads (5–10 seconds, one thread)

Columns give 1× / 2× / 4× observations. The single-thread 1× workloads
take 6–8 seconds. Both SQL macros and their R wrappers use
reference-backed allele validation and emit all selected rows. The
wrappers fetch the full result into R; their peak RSS includes the
returned data frame, while the SQL timings materialize a DuckDB
temporary table. Near-linear time includes CSV decompression and FASTA
access. These are distinct output-materialization workloads, not a
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

## Open findings, ranked by measured exposure

1.  **VCF-count join:** 612 MiB of DuckDB buffer for 1.24 million GIAB
    records and 2,000 output rows after the materialization change. A
    selective indexed input or keyed stream join merits a separate
    design and larger real-input differential; an optimizer hint alone
    does not make memory bounded.
2.  **Site import:** buffer rises 222 → 942 MiB as unique phase-3 chr22
    SNPs rise 125k → 562k; output also rises one-to-one, so this is not
    by itself evidence of a leak. A larger panel/source is required
    before a time-scaling verdict.
3.  **Munge/Metal reference-backed modes:** 200k / 400k / 800k real
    epilepsy records scale close to linearly in time, but the 800k-row
    SQL output occupies about 779 MiB DuckDB buffer; COPY-to-disk output
    and any bounded streaming requirement are separate questions. The R
    wrapper 800k-row result reaches 645 MiB RSS (`r_munge`) and 710 MiB
    (`r_munge_metal`), consistent with returning a large R data frame.
4.  **Liftover:** the entire staged chr22 GRCh37 source holds only ~1.1
    million variants; the current 250k baseline is below the timing
    floor. A larger public source, allele/mapped-output verification,
    and EXPLAIN ANALYZE at that scale are required.
5.  **Not measured with real multi-sample evidence:** BAM/CRAM
    extraction, `verify_sketches`, matched contamination and the other
    SQL-building R wrappers. The staged phased 1000 Genomes chr22 VCFs
    have GT but no FORMAT/AD, so they cannot supply the observed allele
    depths required for `vcf_counts`; the available GIAB FORMAT/AD VCF
    has one sample. The existing three indexed chr22 CRAM extracts are
    not a many-sample cohort. A new real evidence cohort and a checked
    sampling denominator are needed.
6.  **Not measured:** executed BAM/VCF/GFF `*_convert_parquet_sql` COPYs
    and DuckVEP. Both require physical multi-size public inputs, output
    denominator checks and independent disk-space limits.

The first-pass synthetic and 10k–40k-row rows have no extrapolatable
scaling conclusion. This report does not claim completion of the
unmeasured entry points.
