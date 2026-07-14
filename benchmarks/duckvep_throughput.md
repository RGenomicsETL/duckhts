DuckVEP sorted-stream throughput
================

<!-- duckvep_throughput.md is generated from duckvep_throughput.Rmd. -->

This benchmark drives either the rich `duckvep_annotate` result or the
numeric `duckvep_annotate_compact` result through DuckDB’s stable
extension API with nondecreasing sequence-region and position keys. It
records candidate sweep, consequence classification, list
materialization, `unnest`, and final aggregation together. Model loading
and input staging are outside the timed pass.

The checked-in one-transcript workload is a repeatable adapter floor.
The production-density workload uses a pre-staged DuckDB database
containing `bench_regions`, `bench_transcripts`, `bench_exons`, and
`bench_variants`. Its recorded Ensembl-116 GRCh38 model contains 644,292
chromosome transcripts and 5,078,384 exon memberships; the input is a
deterministic 1-in-40 hash sample of 100,957 sorted GIAB sites. The
external staging database is not committed. The production query
restates `ORDER BY` in the timed input CTE, so its result includes the
cost of enforcing the kernel’s nondecreasing-input contract even though
the staging table was already written in that order. The `79cfaaf2` rows
predate that review correction and relied on the single-thread physical
scan preserving insertion order; they remain only as historical
measurements.

## Recorded runs

| run_date   | revision | workload                            | output_mode | threads | variants   | transcripts | exons     | annotated_rows | passes | min_seconds | median_seconds | max_seconds | variants_per_second | ns_per_variant | cpu                                 |
|:-----------|:---------|:------------------------------------|:------------|--------:|:-----------|:------------|:----------|:---------------|-------:|------------:|---------------:|------------:|--------------------:|---------------:|:------------------------------------|
| 2026-07-11 | 8cc22218 | fixture_one_transcript_sorted       | rich        |       1 | 10,000,000 | 1           | 2         | 10,000,000     |      3 |       3.225 |          3.334 |       3.371 |             2999400 |          333.4 | 13th Gen Intel(R) Core(TM) i5-13500 |
| 2026-07-13 | 87f03a2a | fixture_one_transcript_sorted       | rich        |       1 | 10,000,000 | 1           | 2         | 10,000,000     |      3 |       2.812 |          2.824 |       2.825 |             3541076 |          282.4 | 13th Gen Intel(R) Core(TM) i5-13500 |
| 2026-07-14 | 79cfaaf2 | ensembl116_grch38_giab_sites_hash40 | compact     |       1 | 100,957    | 644,292     | 5,078,384 | 1,179,465      |      9 |       0.107 |          0.108 |       0.110 |              934787 |         1069.8 | 13th Gen Intel(R) Core(TM) i5-13500 |
| 2026-07-14 | 79cfaaf2 | ensembl116_grch38_giab_sites_hash40 | rich        |       1 | 100,957    | 644,292     | 5,078,384 | 1,179,465      |      9 |       0.225 |          0.228 |       0.232 |              442794 |         2258.4 | 13th Gen Intel(R) Core(TM) i5-13500 |
| 2026-07-14 | bcca6f6c | ensembl116_grch38_giab_sites_hash40 | compact     |       1 | 100,957    | 644,292     | 5,078,384 | 1,179,465      |      9 |       0.106 |          0.107 |       0.109 |              943523 |         1059.9 | 13th Gen Intel(R) Core(TM) i5-13500 |
| 2026-07-14 | bcca6f6c | ensembl116_grch38_giab_sites_hash40 | rich        |       1 | 100,957    | 644,292     | 5,078,384 | 1,179,465      |      9 |       0.221 |          0.225 |       0.231 |              448698 |         2228.7 | 13th Gen Intel(R) Core(TM) i5-13500 |

Each pass consumes every staged input and checks output cardinality plus
either the rendered consequence-byte total or the numeric
consequence-mask sum. Rows with different workload, output mode, thread
count, host, or variant count are separate measurements and should not
be compared as if only the engine changed.

The production form is run with:

``` sh
Rscript benchmarks/duckvep_throughput.R \
  --database /path/to/staged.duckdb \
  --variants 100957 --passes 9 --warmup 10000 --threads 1 \
  --output compact
```
