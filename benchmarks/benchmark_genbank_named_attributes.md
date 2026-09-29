GenBank named-attribute projections
================

`read_genbank(attributes := [...])` adds one VARCHAR column per
qualifier key. The C API table function learns which columns a query
projects, never which MAP keys, so a query that needs one qualifier used
to build the whole `attributes_map` for every feature. This report
compares that MAP projection with named columns.

The input is the registered NCBI RefSeq record for *Escherichia coli*
K-12 MG1655 (`genbank_ecoli_k12_gbff`, staged by
`duckhtsbench::duckhts_bench_stage_genbank()`), 11,450,954 bytes
yielding 9,306 feature rows. Each measurement is one query in a fresh R
process with one DuckDB thread: the extension is loaded, the query runs
once, and the process reports elapsed seconds and its peak resident set
(`VmHWM`). Three fresh processes per workload run in an interleaved
order; the table reports medians, and the individual runs are in
`benchmark_genbank_named_attributes_runs.csv`. Every query returns one
aggregate row, so transport is not measured. Named columns carry the
same values as the map lookups (`n` and `metric` agree between paired
workloads, and the SQL and R tests compare every column with the
lookup).

Source revision: 14b5f13a59153aee59dd3dff78d060ccb10a3b40. Extension
SHA-256:
73bdb3d6384498437b791dedda387cc3abd5c17cb8f782ac57441f6823de9465. Input
SHA-256:
8a50dc9bc68b8d1b2222f022c086c6a177004692b7a0318a4e2225aacd190673. DuckDB
runtime: 1.5.5. Host: 13th Gen Intel(R) Core(TM) i5-13500; `uptime`
before rendering: 20:22:28 up 357 days, 5:16, 22 users, load average:
3.89, 3.89, 2.97.

| workload             | input_rows | output_rows | metric      | median_seconds | median_peak_rss_mib | input_mib_per_second | runs |
|:---------------------|-----------:|------------:|:------------|---------------:|--------------------:|---------------------:|-----:|
| plain_columns        |       9306 |           1 | 42860429501 |          0.017 |               124.3 |              642.381 |    3 |
| attributes_string    |       9306 |           1 | 1708355     |          0.022 |               123.5 |              496.385 |    3 |
| map_one              |       9306 |           1 | 4556        |          0.025 |               124.3 |              436.819 |    3 |
| named_one            |       9306 |           1 | 4556        |          0.019 |               123.9 |              574.762 |    3 |
| map_three            |       9306 |           1 | 23044       |          0.024 |               124.8 |              455.020 |    3 |
| named_three          |       9306 |           1 | 23044       |          0.020 |               124.1 |              546.024 |    3 |
| named_none_projected |       9306 |           1 | 42860429501 |          0.018 |               123.8 |              606.693 |    3 |

| keys | map_seconds | named_seconds | map_peak_rss_mib | named_peak_rss_mib | speedup |
|-----:|------------:|--------------:|-----------------:|-------------------:|--------:|
|    1 |       0.025 |         0.019 |            124.3 |              123.9 |   1.316 |
|    3 |       0.024 |         0.020 |            124.8 |              124.1 |   1.200 |

## Reading the result

This is a single 1x input. Its runtime is below the 5 s floor that
`STYLE.md` (Scale) sets for timing claims, so the report makes no claim
about how time scales with input size and has no 2x/4x series. It
reports the measured seconds, peak resident memory and input throughput
for this file only. Memory is the process peak, dominated by DuckDB and
R baseline plus the parser’s one-record buffers;
`benchmark_genbank_memory.md` is the report for the memory budget across
input sizes.

`named_none_projected` requests three keys at bind but projects none, so
the scan computes no qualifier values. `map_one` builds the entire MAP
for each feature to read a single key, which is the cost named columns
avoid; `named_one` and `named_three` compute values only for matching
keys. The paired workloads must return identical metrics, and the render
stops if they differ. Workloads with different aggregation expressions
are not isolated parsing costs. No earlier recorded GenBank
named-attribute workload exists; the nearest report for this change is
`benchmark_gff_named_attributes.md`, which measures a different format
on a different input and is not comparable.
