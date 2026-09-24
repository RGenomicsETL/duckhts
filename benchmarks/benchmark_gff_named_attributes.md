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
reader uses one worker in either setting.

Source revision: c4728030d2770b591f0875206aac0af9be9d2c26. DuckDB
runtime: 1.5.5. Extension SHA-256:
1d868ba87527209147d86dfdbf796f00ee29ecd0cbc5d81951e8ee21a8331421. Input
SHA-256:
e8ed48bef6a44fdf0db7c10a551d4398aa341318d00fbd9efd69530593106846.

| threads      | workload             | input_rows | output_rows | metric          | median_seconds | runs |
|:-------------|:---------------------|-----------:|------------:|:----------------|---------------:|-----:|
| 1            | plain_columns        |    1299172 |           1 | 193775078663005 |          0.540 |    5 |
| 1            | named_none_projected |    1299172 |           1 | 193775078663005 |          0.542 |    5 |
| 1            | attributes_string    |    1299172 |           1 | 577618955       |          0.690 |    5 |
| 1            | map_parent           |    1299172 |           1 | 1243771         |          2.326 |    5 |
| 1            | regexp_three         |    1299172 |           1 | 76627121        |          2.371 |    5 |
| 1            | named_three          |    1299172 |           1 | 76627121        |          0.900 |    5 |
| 1            | named_one_projected  |    1299172 |           1 | 25136049        |          0.788 |    5 |
| default (20) | plain_columns        |    1299172 |           1 | 193775078663005 |          0.542 |    5 |
| default (20) | named_none_projected |    1299172 |           1 | 193775078663005 |          0.543 |    5 |
| default (20) | attributes_string    |    1299172 |           1 | 577618955       |          0.691 |    5 |
| default (20) | map_parent           |    1299172 |           1 | 1243771         |          2.327 |    5 |
| default (20) | regexp_three         |    1299172 |           1 | 76627121        |          2.376 |    5 |
| default (20) | named_three          |    1299172 |           1 | 76627121        |          0.905 |    5 |
| default (20) | named_one_projected  |    1299172 |           1 | 25136049        |          0.784 |    5 |

The one-key projection requests three keys at bind but scans only
`Parent`. `named_none_projected` requests three keys at bind but reads
only the physical coordinate columns; its scan skips attribute parsing
entirely. Comparisons across workloads reflect different aggregation and
expression costs, not isolated parsing cycles. The local duckdb-wasm
browser test could not load this extension: `duckdb_get_list_size` has
an incompatible wasm import type with the staged browser runtime.
Browser wasm timings are unavailable.
