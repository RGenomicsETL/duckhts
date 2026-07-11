DuckVEP sorted-stream throughput
================

<!-- duckvep_throughput.md is generated from duckvep_throughput.Rmd. -->

This benchmark drives `duckvep_annotate` through DuckDB’s stable
extension API with nondecreasing sequence-region and position keys. It
records vector adapter, sweep reuse, coding classification, list
materialization, `unnest`, and final aggregation together.

The current checked-in workload has one transcript and two exons. It is
a repeatable hot-transcript floor, not a whole-genome throughput claim.
Production qualification also requires an imported Ensembl model and a
sorted gnomAD-scale variant stream spanning realistic transcript
density.

## Recorded runs

| run_date   | revision | workload                      | threads | variants   | passes | min_seconds | median_seconds | max_seconds | variants_per_second | ns_per_variant | cpu                                 |
|:-----------|:---------|:------------------------------|--------:|:-----------|-------:|------------:|---------------:|------------:|--------------------:|---------------:|:------------------------------------|
| 2026-07-11 | 8cc22218 | fixture_one_transcript_sorted |       1 | 10,000,000 |      3 |       3.225 |          3.334 |       3.371 |             2999400 |          333.4 | 13th Gen Intel(R) Core(TM) i5-13500 |

Each pass consumes every generated input and checks output cardinality
and a consequence-string checksum. Rows with different workload, thread
count, host, or variant count are separate measurements and should not
be compared as if only the engine changed.
