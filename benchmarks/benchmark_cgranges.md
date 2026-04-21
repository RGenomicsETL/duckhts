DuckHTS cgranges Benchmark
================

<!-- benchmark_cgranges.md is generated from benchmark_cgranges.Rmd. -->

# Benchmark

This benchmark measures:

- `duckhts_cgranges_from_query(...)` + `duckhts_cgranges_index(...)`
- repeated bound overlap probes through `duckhts_cgranges_overlaps(...)`
- upstream `bedtk flt` as the overlap-existence comparator
- wall time and peak RSS for both engines

It uses a deterministic synthetic interval set so the benchmark is
stable and rerunnable.

# Run

``` sh
python3 scripts/cgranges_benchmark.py \
  --extension build/release/duckhts.duckdb_extension \
  --bedtk .sync/bedtk/bedtk \
  --subjects 50000 \
  --queries 5000 \
  --passes 3 \
  --out-dir .tmp/cgranges_benchmark
```

# Configuration

| parameter | value |
|:----------|------:|
| subjects  | 50000 |
| queries   |  5000 |
| passes    |     3 |

# Results

| tool             | subjects | queries | passes | build_index_sec | query_total_sec | query_pass_1_sec | total_elapsed_sec | peak_rss_mb | matched_query_intervals | total_hits | time_per_query_ms |
|:-----------------|---------:|--------:|-------:|----------------:|----------------:|-----------------:|------------------:|------------:|------------------------:|-----------:|------------------:|
| duckhts_cgranges |    50000 |    5000 |      3 |            0.03 |           1.733 |            0.660 |             1.873 |       125.4 |                    2239 |       3000 |            0.1156 |
| bedtk            |    50000 |    5000 |      3 |            0.00 |           0.030 |            0.007 |             0.060 |        14.5 |                    2239 |         NA |            0.0020 |

# Notes

- `bedtk` reports overlap-existence for the whole query BED in one pass.
- DuckHTS builds an immutable registry index once, then reuses it across
  probe passes.
- `matched_query_intervals` should agree across both engines.
- `build_index_sec` is only meaningful for DuckHTS because the public
  API exposes an explicit build/finalize step.
