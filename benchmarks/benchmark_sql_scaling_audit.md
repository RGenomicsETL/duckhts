SQL entry-point scaling audit: measured subset
================

<!-- benchmark_sql_scaling_audit.md is rendered from this file. -->

## Workload and measurement contract

Source revision: `ecf3d4a6a48dc81872ffb65e6a2a6e90346384e6` (no
extension or SQL changes). DuckDB R 1.5.5; one or four DuckDB threads;
one cold DuckDB connection per timed query. Linux `VmHWM` is process
high-water resident memory in KiB, converted to MiB; DuckDB JSON
profiling supplies `system_peak_buffer_memory`, converted to MiB. RSS
includes process startup and input-table preparation; DuckDB buffer
memory is the profiled query’s peak. Wall time covers only the query
execution, including output materialization to a temporary table; input
staging and reading into the GIAB input table are excluded. Each cell is
one run, not a confidence interval. This is an exploratory scaling
check, not a performance regression claim.

The Somalier input is the `somalier-synthetic` registry workload: 17,000
autosomal sites and `samples × 17,000` count-evidence records, with one
directional pair per additional sample. For the hashes, the dense
panel/frequency prefix is 4,250, 8,500 or 17,000 sites; other Somalier
calls use the whole panel. The public normalization input is the
registry’s GIAB NIST v4.2.1 HG002 GRCh38 phased VCF
(`geno_giab_phased_source`) with the registered GRCh38 FASTA
(`liftover_grch38_fasta`), both checksum-validated. Prefixes of 10,000,
20,000 and 40,000 physical VCF records are materialized before timing.
Output denominators are one digest for each hash, `size` sketches or
charr rows, `size - 1` matched pairs, and `size` normalized rows. The
normalization probe uses the default site-preserving mode;
split-multiallelic and reference-block cases are not covered.

Reproduce from the repository root after `make release`:

``` sh
DUCKHTS_CACHE_DIR=/tmp/duckhts-sql-scaling-cache Rscript benchmarks/benchmark_sql_scaling_audit_run.R build/release/duckhts.duckdb_extension
DUCKHTS_CACHE_DIR=/tmp/duckhts-sql-scaling-cache Rscript benchmarks/benchmark_sql_scaling_audit_run.R build/release/duckhts.duckdb_extension --fixed-pairs
Rscript benchmarks/benchmark_sql_scaling_audit_run.R build/release/duckhts.duckdb_extension --norm-run
Rscript -e 'rmarkdown::render("benchmarks/benchmark_sql_scaling_audit.Rmd", knit_root_dir = getwd(), quiet = TRUE)'
```

The default cache must contain the registered public VCF and FASTA for
the `--norm-run` command; see the `genotype-phase-set` and
`liftover-reference-bundle` registry workloads for staging. Raw
per-process observations are in `sql_scaling_audit_measurements.tsv`,
`sql_scaling_audit_fixed_pairs.tsv`, and `sql_scaling_audit_norm.tsv`.

## Measured scaling

The exponent is `log(time at 4× / time at 1×) / log(4)`; startup and
planner cost dominate the fastest queries. Time, RSS and buffer columns
give 1× / 2× / 4× values. A verdict applies only to the displayed input
range and workload.

| entry_point                            | threads | input_1x_2x_4x                       | time_s                | peak_rss_mib    | peak_buffer_mib | exponent | verdict                              |
|:---------------------------------------|--------:|:-------------------------------------|:----------------------|:----------------|:----------------|:---------|:-------------------------------------|
| duckhts_somalier_panel_sha256          |       1 | 4,250 / 8,500 / 17,000 sites         | 0.013 / 0.017 / 0.024 | 136 / 136 / 139 | 13 / 21 / 23    | 0.44     | No superlinear signal at these sizes |
| duckhts_somalier_panel_sha256          |       4 | 4,250 / 8,500 / 17,000 sites         | 0.013 / 0.013 / 0.017 | 139 / 146 / 148 | 14 / 24 / 36    | 0.19     | No superlinear signal at these sizes |
| duckhts_somalier_frequency_sha256      |       1 | 4,250 / 8,500 / 17,000 sites         | 0.022 / 0.031 / 0.047 | 143 / 143 / 146 | 21 / 35 / 46    | 0.55     | No superlinear signal at these sizes |
| duckhts_somalier_frequency_sha256      |       4 | 4,250 / 8,500 / 17,000 sites         | 0.018 / 0.019 / 0.026 | 154 / 159 / 165 | 33 / 70 / 80    | 0.27     | No superlinear signal at these sizes |
| duckhts_somalier_prepare_sketches      |       1 | 8 / 16 / 32 samples (17k sites each) | 0.049 / 0.065 / 0.099 | 146 / 146 / 147 | 40 / 42 / 42    | 0.51     | No superlinear signal at these sizes |
| duckhts_somalier_prepare_sketches      |       4 | 8 / 16 / 32 samples (17k sites each) | 0.036 / 0.049 / 0.050 | 160 / 161 / 163 | 73 / 81 / 87    | 0.24     | No superlinear signal at these sizes |
| duckhts_somalier_charr                 |       1 | 8 / 16 / 32 samples (17k sites each) | 0.129 / 0.164 / 0.246 | 177 / 173 / 173 | 119 / 99 / 99   | 0.47     | No superlinear signal at these sizes |
| duckhts_somalier_charr                 |       4 | 8 / 16 / 32 samples (17k sites each) | 0.077 / 0.074 / 0.106 | 200 / 198 / 202 | 197 / 248 / 234 | 0.23     | No superlinear signal at these sizes |
| duckhts_somalier_matched_contamination |       1 | 8 / 16 / 32 samples (17k sites each) | 0.906 / 1.833 / 3.703 | 193 / 215 / 226 | 144 / 165 / 192 | 1.02     | No superlinear signal at these sizes |
| duckhts_somalier_matched_contamination |       4 | 8 / 16 / 32 samples (17k sites each) | 0.870 / 1.668 / 3.461 | 220 / 249 / 304 | 281 / 399 / 485 | 1.00     | No superlinear signal at these sizes |
| duckhts_bcftools_norm                  |       1 | 10k / 20k / 40k VCF records          | 0.115 / 0.120 / 0.136 | 135 / 135 / 136 | 10 / 19 / 21    | 0.12     | No superlinear signal at these sizes |
| duckhts_bcftools_norm                  |       4 | 10k / 20k / 40k VCF records          | 0.091 / 0.124 / 0.145 | 135 / 138 / 141 | 10 / 19 / 22    | 0.34     | No superlinear signal at these sizes |

## Candidate-set control

For `duckhts_somalier_matched_contamination`, a second staged run grows
the evidence from 8 to 32 samples (136k / 272k / 544k input records)
while holding the candidate set at three selected samples and two output
pairs. Peak RSS (MiB) at one thread is 187 / 185 / 185; at four threads
it is 206 / 208 / 216. DuckDB peak buffer memory (MiB) at one thread is
127 / 127 / 121; at four threads it is 269 / 248 / 270. All six runs
produce two pairs. This probe does not show memory tracking unselected
evidence. Matched contamination’s full-candidate run is close to linear
in candidate pairs at these sizes; its ordered `list()` aggregates hold
per-pair site vectors, so larger candidate or panel sweeps remain
necessary.

## Open audit work (ranked by likely exposure, not confirmed defect)

1.  **Somalier VCF/BAM extraction, `duckhts_somalier_vcf_counts`,
    `duckhts_somalier_import_sites` and the R extraction wrappers:** no
    public multi-size input was exercised here. The VCF-count SQL has
    materialized validation and a per-record alternate/AD path; real
    indexed VCF and BAM workloads need 1×/2×/4× records, thread sweeps
    and operator profiles.
2.  **`duckhts_somalier_verify_sketches`, R sketch wrappers, and matched
    contamination on a real multi-sample cohort:** synthetic count
    evidence cannot resolve transport/decompression or
    high-panel-cardinality costs. Hold selected pairs fixed and then
    grow candidates on a publicly staged cohort.
3.  **`duckdb_munge`, `duckdb_munge_metal`, `duckdb_liftover` and their
    R builders:** their public input and reference workloads need
    independent staging, output-count controls, two thread settings and
    complete plan profiling. No scaling verdict is assigned.
4.  **BAM/BCF/GFF/tabix `*_convert_parquet_sql` execution:** these
    macros produce COPY statements; scaling must include executing the
    generated statements, including partition/cardinality and disk
    output denominators. SQL-string generation alone would not test the
    writer.
5.  **DuckVEP annotate, transcript projection and Ensembl macros:**
    lower priority, not measured; high-cardinality transcript cohorts
    and repeated exon/reference joins need a separate audit.

No query in the measured subset exceeded a 1.2 time-scaling exponent.
Accordingly there is no flagged plan to diagnose, semantic rewrite, or
before/after differential in this report. This is not an audit verdict
on the unmeasured entry points. A source search of `src/` and
`r/Rduckhts/R/` found no `row_number() OVER ()` without ordering; the
`row_number()` call in `duckhts_parquet_metadata_args` orders by
metadata priority. No order-dependent SQL was changed.
