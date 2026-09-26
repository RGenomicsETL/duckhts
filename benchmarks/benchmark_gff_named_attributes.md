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
Core i5-13500 (20 logical CPUs); `uptime` before the run reported load
averages of 1.51, 1.09, and 0.79.

Source revision: 58e6e2b1ea404e246df84e55ee589e14c3048aeb. DuckDB
runtime: 1.5.5. Previous medians: source revision
c4728030d2770b591f0875206aac0af9be9d2c26; extension SHA-256
1d868ba87527209147d86dfdbf796f00ee29ecd0cbc5d81951e8ee21a8331421.
Extension SHA-256:
4e0997a92ef2cfc6dd7a6777dba4b6a9b6edce80e32b3be9c523bedfcdcb2e82. Input
SHA-256:
e8ed48bef6a44fdf0db7c10a551d4398aa341318d00fbd9efd69530593106846.

| threads      | workload             | input_rows | output_rows | metric          | median_seconds | runs | previous_median_seconds |
|:-------------|:---------------------|-----------:|------------:|:----------------|---------------:|-----:|------------------------:|
| 1            | plain_columns        |    1299172 |           1 | 193775078663005 |          0.479 |    5 |                   0.540 |
| 1            | named_none_projected |    1299172 |           1 | 193775078663005 |          0.479 |    5 |                   0.542 |
| 1            | attributes_string    |    1299172 |           1 | 577618955       |          0.692 |    5 |                   0.690 |
| 1            | map_parent           |    1299172 |           1 | 1243771         |          2.242 |    5 |                   2.326 |
| 1            | regexp_three         |    1299172 |           1 | 76627121        |          2.358 |    5 |                   2.371 |
| 1            | named_three          |    1299172 |           1 | 76627121        |          0.850 |    5 |                   0.900 |
| 1            | named_one_projected  |    1299172 |           1 | 25136049        |          0.731 |    5 |                   0.788 |
| default (20) | plain_columns        |    1299172 |           1 | 193775078663005 |          0.485 |    5 |                   0.542 |
| default (20) | named_none_projected |    1299172 |           1 | 193775078663005 |          0.480 |    5 |                   0.543 |
| default (20) | attributes_string    |    1299172 |           1 | 577618955       |          0.688 |    5 |                   0.691 |
| default (20) | map_parent           |    1299172 |           1 | 1243771         |          2.262 |    5 |                   2.327 |
| default (20) | regexp_three         |    1299172 |           1 | 76627121        |          2.392 |    5 |                   2.376 |
| default (20) | named_three          |    1299172 |           1 | 76627121        |          0.853 |    5 |                   0.905 |
| default (20) | named_one_projected  |    1299172 |           1 | 25136049        |          0.732 |    5 |                   0.784 |

The one-key projection requests three keys at bind but scans only
`Parent`. `named_none_projected` requests three keys at bind but reads
only the physical coordinate columns; its scan skips attribute parsing
entirely. Comparisons across workloads reflect different aggregation and
expression costs, not isolated parsing cycles.
`make wasm-playwright-test` built the extension from this checkout for
duckdb-wasm 1.31.0; all six Chromium tests passed, including
`test/wasm/gff-attributes.spec.ts`. Browser wasm timings were not
measured.
