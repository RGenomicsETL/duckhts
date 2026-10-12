GENCODE GFF3 named-attribute projections
================

The source is the public GENCODE mouse vM25 BASIC GFF3, staged by
`scripts/stage_tabix_split.R --offline` from the `tabix_split_gff3`
registry artifact (24,025,067 compressed bytes; 1,299,172 annotation
rows). Each query streams the whole file in sequential mode. Measured
outputs are scalar aggregates of the stated columns, not `count(*)`
shortcuts. Five measured runs follow one warm-up per workload and thread
setting. Times include gzip decode, row scan, aggregation and R/DBI
query execution. Default threads are the runtime’s own default; the GFF
reader uses one worker in either setting. The host is a 13th Gen Intel
Core i5-13500 (20 logical CPUs) with two core types, so the render
records the CPU affinity of its process in the runs file and compares
with the previous render only on the same CPUs. Choose them with
`taskset -c`. `GFF_NAMED_SIZE=reduced` runs a driver check with one
measured run per cell, and writes its runs to a temporary file instead
of `benchmark_gff_named_attributes_runs.csv`; it is not a measurement.

``` sh
taskset -c 8-11 Rscript -e 'rmarkdown::render("benchmarks/benchmark_gff_named_attributes.Rmd")'
```

Source revision: 68dfa0332ad9ea34eea556eddea231cd83d842f6. DuckDB
runtime: 1.5.5. Host load before the run: 1.36, 1.36, 0.95. CPU
affinity: 8-11. Extension SHA-256:
66598566776aa44faffb68a8e0ba62961346c11db160962d115c30291d33586e. Input
SHA-256:
e8ed48bef6a44fdf0db7c10a551d4398aa341318d00fbd9efd69530593106846.

| threads      | workload             | input_rows | output_rows | metric          | median_seconds | runs |
|:-------------|:---------------------|-----------:|------------:|:----------------|---------------:|-----:|
| 1            | plain_columns        |    1299172 |           1 | 193775078663005 |          0.479 |    5 |
| 1            | named_none_projected |    1299172 |           1 | 193775078663005 |          0.480 |    5 |
| 1            | attributes_string    |    1299172 |           1 | 577618955       |          0.692 |    5 |
| 1            | map_parent           |    1299172 |           1 | 1243771         |          2.588 |    5 |
| 1            | regexp_three         |    1299172 |           1 | 76627121        |          2.373 |    5 |
| 1            | named_three          |    1299172 |           1 | 76627121        |          0.842 |    5 |
| 1            | named_one_projected  |    1299172 |           1 | 25136049        |          0.732 |    5 |
| default (20) | plain_columns        |    1299172 |           1 | 193775078663005 |          0.482 |    5 |
| default (20) | named_none_projected |    1299172 |           1 | 193775078663005 |          0.479 |    5 |
| default (20) | attributes_string    |    1299172 |           1 | 577618955       |          0.688 |    5 |
| default (20) | map_parent           |    1299172 |           1 | 1243771         |          2.582 |    5 |
| default (20) | regexp_three         |    1299172 |           1 | 76627121        |          2.376 |    5 |
| default (20) | named_three          |    1299172 |           1 | 76627121        |          0.845 |    5 |
| default (20) | named_one_projected  |    1299172 |           1 | 25136049        |          0.725 |    5 |

The one-key projection requests three keys at bind but scans only
`Parent`. `named_none_projected` requests three keys at bind but reads
only the physical coordinate columns; its scan skips attribute parsing
entirely. Comparisons across workloads reflect different aggregation and
expression costs, not isolated parsing cycles.
`make wasm-playwright-test` built the extension from this checkout for
duckdb-wasm 1.31.0; all six Chromium tests passed, including
`test/wasm/gff-attributes.spec.ts`. Browser wasm timings were not
measured.

## Change from the previous render

The previous render, 60b36911, ran the same queries on the same input
with the same driver and host. The previous render recorded no CPU
affinity, so the CPUs it ran on are unknown and the ratios may include a
change of core type. Host load differs between renders. Median seconds
per cell:

| setting | workload             | seconds_before | seconds_after | time_ratio |
|:--------|:---------------------|---------------:|--------------:|-----------:|
| 1       | plain_columns        |          0.479 |         0.479 |      1.000 |
| 1       | named_none_projected |          0.479 |         0.480 |      1.002 |
| 1       | attributes_string    |          0.692 |         0.692 |      1.000 |
| 1       | map_parent           |          2.242 |         2.588 |      1.154 |
| 1       | regexp_three         |          2.358 |         2.373 |      1.006 |
| 1       | named_three          |          0.850 |         0.842 |      0.991 |
| 1       | named_one_projected  |          0.731 |         0.732 |      1.001 |
| default | plain_columns        |          0.485 |         0.482 |      0.994 |
| default | named_none_projected |          0.480 |         0.479 |      0.998 |
| default | attributes_string    |          0.688 |         0.688 |      1.000 |
| default | map_parent           |          2.262 |         2.582 |      1.141 |
| default | regexp_three         |          2.392 |         2.376 |      0.993 |
| default | named_three          |          0.853 |         0.845 |      0.991 |
| default | named_one_projected  |          0.732 |         0.725 |      0.990 |
