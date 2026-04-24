DuckHTS cgranges Streaming-Provider Benchmark
================

<!-- benchmark_cgranges.md is generated from benchmark_cgranges.Rmd. -->

# Benchmark

This benchmark measures the cgranges path that matters for provider
streaming:

- build one session-scoped cgranges target index with
  `duckhts_cgranges_from_query(...)` / `duckhts_cgranges_index(...)`
- stream query intervals from `read_bed(...)`
- filter query rows with the vectorized scalar predicate
  `duckhts_cgranges_has_overlap(...)`
- annotate query rows with vectorized scalar counts via
  `duckhts_cgranges_count_overlaps(...)`
- compare overlap-existence against `bedtk flt`
- compare overlap-existence and overlap counts against
  `bedtools intersect -u` and `bedtools intersect -c` when `bedtools` is
  installed

The old benchmark shape generated thousands of `UNION ALL` calls to
`duckhts_cgranges_overlaps(...)`, or issued one SQL statement per probe.
That measured bind/SQL dispatch overhead rather than the desired
streaming provider path. The current scripts avoid that pattern.

# Run

Synthetic, deterministic default used for the rendered table:

``` sh
python3 scripts/cgranges_benchmark.py \
  --extension build/release/duckhts.duckdb_extension \
  --bedtk .sync/bedtk/bedtk \
  --bedtools bedtools \
  --subjects 50000 \
  --queries 5000 \
  --passes 3 \
  --out-dir .tmp/cgranges_benchmark
```

Real DuckBedQC BED files can be benchmarked with:

``` sh
python3 scripts/cgranges_benchmark_real.py \
  --extension build/release/duckhts.duckdb_extension \
  --passes 1 \
  --out-dir .tmp/cgranges_benchmark_real
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

| tool     | variant       | subject_intervals | query_intervals | passes | build_index_sec | query_total_sec | query_pass_1_sec | total_elapsed_sec | peak_rss_mb | matched_query_intervals | total_hits | time_per_query_ms |
|:---------|:--------------|------------------:|----------------:|-------:|----------------:|----------------:|-----------------:|------------------:|------------:|------------------------:|-----------:|------------------:|
| duckhts  | scalar_filter |             50000 |            5000 |      3 |           0.011 |           0.004 |            0.002 |             0.102 |        68.4 |                    2239 |         NA |            0.0003 |
| duckhts  | scalar_count  |             50000 |            5000 |      3 |           0.013 |           0.005 |            0.002 |             0.105 |        68.0 |                    2239 |       3000 |            0.0003 |
| bedtk    | flt           |             50000 |            5000 |      3 |           0.000 |           0.022 |            0.007 |             0.054 |        15.0 |                    2239 |         NA |            0.0015 |
| bedtools | intersect_u   |             50000 |            5000 |      3 |           0.000 |           0.057 |            0.019 |             0.094 |        21.9 |                    2239 |         NA |            0.0038 |
| bedtools | intersect_c   |             50000 |            5000 |      3 |           0.000 |           0.063 |            0.020 |             0.095 |        21.7 |                    2239 |       3000 |            0.0042 |

# Semantic checks

| check                                                                                            | result |
|:-------------------------------------------------------------------------------------------------|:-------|
| matched query intervals agree for DuckHTS scalar_filter, bedtk flt, and bedtools -u when present | TRUE   |
| overlap counts agree for DuckHTS scalar_count and bedtools -c when present                       | TRUE   |

# Notes

- `duckhts:scalar_filter` is the row-preserving provider-streaming
  predicate path.
- `duckhts:scalar_count` computes one cgranges overlap count per
  streamed query row and is comparable to `bedtools intersect -c`.
- `bedtk flt` and `bedtools intersect -u` are overlap-existence
  comparators; they do not report total hit counts.
- DuckHTS reports explicit target-index build time because the SQL API
  exposes index construction as a session-scoped operation. The external
  tools build any internal structures inside their command runtime.
- The rendered table is intentionally a modest synthetic run. Use
  `scripts/cgranges_benchmark_real.py` for the larger DuckBedQC
  WGS/exome BED workload, and optional `--limit-subjects` /
  `--limit-queries` for quick smoke runs.
