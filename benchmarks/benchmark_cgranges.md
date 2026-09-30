DuckHTS cgranges Streaming-Provider Benchmark
================

<!-- benchmark_cgranges.md is generated from benchmark_cgranges.Rmd. -->

# Benchmark

This benchmark measures the cgranges path that matters for provider
streaming:

- build one session-scoped cgranges target index with
  `duckhts_cgranges_from_table(...)`
- stream query intervals from `read_bed(...)`
- filter query rows with the vectorized scalar predicate
  `duckhts_cgranges_has_overlap(...)`
- annotate query rows with vectorized scalar counts via
  `duckhts_cgranges_count_overlaps(...)`
- expand streaming provider rows to one row per hit with
  `duckhts_cgranges_overlaps_list(...)` plus `UNNEST(...)`
- compare overlap-existence against `bedtk flt`
- compare overlap-existence, overlap counts, and one-row-per-hit
  expansion against `bedtools intersect -u`, `bedtools intersect -c`,
  and `bedtools intersect -wa -wb` when `bedtools` is installed

The old benchmark shape generated thousands of `UNION ALL` calls to
`duckhts_cgranges_overlaps(...)`, or issued one SQL statement per probe.
That measured bind/SQL dispatch overhead rather than the desired
streaming provider path. The current scripts avoid that pattern.

# Run

Synthetic, deterministic default used for the rendered table:

``` sh
Rscript r/duckhtsbench/scripts/cgranges_benchmark.R \
  --extension build/release/duckhts.duckdb_extension \
  --bedtk .sync/bedtk/bedtk \
  --bedtools bedtools \
  --subjects 50000 \
  --queries 5000 \
  --passes 3
```

Real DuckBedQC BED files can be benchmarked after staging their pinned
source revision:

``` sh
make stage-duckbedqc-data
python3 scripts/cgranges_benchmark_real.py \
  --extension build/release/duckhts.duckdb_extension \
  --passes 1
```

For a shell/CLI-only smoke path that still avoids generated bulk-query
SQL:

``` sh
scripts/cgranges_benchmark_cli.sh
```

# Configuration

| parameter | value                        |
|:----------|:-----------------------------|
| dataset   | synthetic deterministic BED4 |
| subjects  | 50000                        |
| queries   | 5000                         |
| passes    | 3                            |
| bedtk     | available                    |
| bedtools  | available                    |

# Results

| tool     | variant         | subject_intervals | query_intervals | passes | build_index_sec | query_total_sec | query_pass_1_sec | total_elapsed_sec | peak_rss_mb | matched_query_intervals | total_hits | time_per_query_ms |
|:---------|:----------------|------------------:|----------------:|-------:|----------------:|----------------:|-----------------:|------------------:|------------:|------------------------:|-----------:|------------------:|
| duckhts  | scalar_filter   |             50000 |            5000 |      3 |           0.019 |           0.006 |            0.002 |             0.025 |          NA |                    2261 |         NA |            0.0004 |
| duckhts  | scalar_count    |             50000 |            5000 |      3 |           0.019 |           0.009 |            0.002 |             0.028 |          NA |                    2261 |       2942 |            0.0006 |
| duckhts  | scalar_expand   |             50000 |            5000 |      3 |           0.019 |           0.023 |            0.008 |             0.042 |          NA |                    2261 |       2942 |            0.0015 |
| bedtk    | flt             |             50000 |            5000 |      3 |           0.000 |           0.050 |            0.015 |             0.050 |          NA |                    2261 |         NA |            0.0033 |
| bedtools | intersect_u     |             50000 |            5000 |      3 |           0.000 |           0.083 |            0.026 |             0.083 |          NA |                    2261 |         NA |            0.0055 |
| bedtools | intersect_c     |             50000 |            5000 |      3 |           0.000 |           0.085 |            0.028 |             0.085 |          NA |                    2261 |       2942 |            0.0057 |
| bedtools | intersect_wa_wb |             50000 |            5000 |      3 |           0.000 |           0.114 |            0.038 |             0.114 |          NA |                    2261 |       2942 |            0.0076 |

# Semantic checks

| check                                                                                                        | result |
|:-------------------------------------------------------------------------------------------------------------|:-------|
| matched query intervals agree for DuckHTS scalar_filter, bedtk flt, and bedtools -u when present             | TRUE   |
| overlap counts agree for DuckHTS scalar_count, scalar_expand, bedtools -c, and bedtools -wa -wb when present | TRUE   |

# Notes

- `duckhts:scalar_filter` is the row-preserving provider-streaming
  predicate path.
- `duckhts:scalar_count` computes one cgranges overlap count per
  streamed query row and is comparable to `bedtools intersect -c`.
- `duckhts:scalar_expand` uses `duckhts_cgranges_overlaps_list(...)`
  plus `UNNEST(...)` to emit one row per hit while preserving streamed
  provider rows; it is comparable to `bedtools intersect -wa -wb`.
- `bedtk flt` and `bedtools intersect -u` are overlap-existence
  comparators; they do not report total hit counts.
- DuckHTS reports explicit target-index build time because the SQL API
  exposes index construction as a session-scoped operation. The external
  tools build any internal structures inside their command runtime.
- The rendered table is intentionally a modest synthetic run. Use
  `scripts/cgranges_benchmark_real.py` for the larger DuckBedQC
  WGS/exome BED workload, and optional `--limit-subjects` /
  `--limit-queries` for quick smoke runs.
