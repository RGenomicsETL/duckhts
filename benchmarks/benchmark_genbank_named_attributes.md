GenBank named-qualifier projections
================

This report measures `read_genbank(..., attributes := [...])` against the two
ways a query could obtain the same qualifiers before it existed: the
`attributes_map` MAP and the raw `attributes` string. The input is NCBI RefSeq
release 237’s plasmid division (workload `genbank-plasmid`): its first
`plasmid.N.genomic.gbff.gz` part alone, its first two parts joined, and its
first four parts joined, each part pinned by NCBI’s published MD5 and each
joined file by its own SHA-256, byte size and record count. Joining consecutive
release parts grows the input while the records themselves stay what the
release supplies. The reader’s work grows with feature rows and input bytes,
not LOCUS records: feature density differs between parts, so the scaling table
compares time with the feature-row and byte ratios. A record whose only
feature is `source` emits no row, so row counts are not a multiple of record
counts. Stage all four parts once; the render itself makes no network access:

``` sh
Rscript -e 'duckhtsbench::duckhts_bench_stage_genbank_plasmid()'
```

NCBI serves only the current RefSeq release’s sequence files and archives the
release catalog, not the files, so the parts can be staged only while release
237 is current. Once NCBI moves on, staging stops before downloading and names
the archived `release237.files.installed` catalog that still lists the
registered MD5s; re-pinning to the current release and re-rendering restores
the workload.

Every workload is a DuckDB table macro whose parameters and output columns
carry a `p_` prefix, so no name can resolve to a reader column. The input path
enters as a SQL variable read through `getvariable()`, so the macro body never
contains it. Each workload returns one row: a row count and a scalar metric
over the stated columns.

Every observation is one fresh DuckDB CLI process under `/usr/bin/time`: it
loads the extension, runs one workload with one thread, writes its row to a CSV
file, and exits. The JSON profiler supplies the query latency and the peak
buffer-manager and temporary-directory bytes; `/usr/bin/time` supplies the
process’s peak RSS. `max_temp_directory_size` is zero, so any spill fails the
render. Before timing, one fresh process per input asserts row-level parity
(every named column equals `attributes_map[key]` on every row). Each workload’s
metric then has an independent expectation computed through the other
representation in a separate scan: a MAP-only scan gives the coordinate,
`named_one_projected` and `named_three` metrics, and a named-only scan gives
the `map_gene` and `map_three` metrics. `attributes_string` has no second
representation and is checked for its row count. One unrecorded warm-up per
workload and input must match its expectation, every timed observation must
reproduce the warm-up answer, and `named_three` and `map_three` must agree on
the same three-key metric. Five timed repetitions follow. Each
repetition visits the inputs in turn and runs the workloads in a rotated
order, reversed on even repetitions, so no workload always runs first and
every pair of workloads runs in both orders. The one-part run exceeds the
five-second floor `STYLE.md` sets for a timing verdict.

The memory budget is the reader’s contract: peak memory proportional to the
largest record, not to records consumed. It is declared before any timed run:
the fixed overhead (median peak RSS of three empty runs that load the extension
and scan nothing) plus three times the reader’s live state, taken as four times
the largest record’s bytes plus one 1 MiB output chunk. The ceiling applies to
every observation of every workload at every input size. The one-thread condition is pinned to the highest logical CPU in the
rendering process’s affinity mask, which each CLI process inherits, so wrap the
render in `taskset -c` to choose it. `read_genbank` declares one worker, so a
multi-thread condition would time the same scan and is not recorded.
`DUCKDB_CLI` selects the CLI (default: `duckdb` on `PATH`) and
`DUCKHTS_EXTENSION` the extension build.

``` sh
taskset -c 8-11 Rscript -e 'rmarkdown::render("benchmarks/benchmark_genbank_named_attributes.Rmd")'
```

Source revision: c4a73210ab8461465ffe306e5e5cac0541c04971; measured `src` tree: b0474026510f692da693fe57fd8d1c56d79cd70b (compare with `git rev-parse HEAD:src` at the reviewed head, since committing this render moves the commit hash but not the tree). DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI).
Extension SHA-256: 5fd33730e8af91c365ac3dcf58091e3faf3b2d49c7849ef283a8db8450bf4bfe.
Host load before the run: 0.62, 0.32, 0.22. One thread, affinity 19, 5 timed fresh-process repetitions per cell; `seconds` is the profiler’s query latency, and every observation is in [`benchmark_genbank_named_attributes_runs.csv`](benchmark_genbank_named_attributes_runs.csv).

Inputs, with the registered byte sizes and the feature rows each emits:

| parts | records |     rows |      bytes |
|------:|--------:|---------:|-----------:|
|     1 |   23713 |  2563736 | 1693174640 |
|     2 |   43922 |  5263269 | 3392019192 |
|     4 |   71677 | 10670134 | 6585648536 |

| parts | records | workload             | input_rows | output_rows | metric        | min_seconds | median_seconds | max_seconds | repetitions |
|------:|--------:|:---------------------|-----------:|------------:|:--------------|------------:|---------------:|------------:|------------:|
|     1 |   23713 | plain_columns        |    2563736 |           1 | 1301781980348 |       6.085 |          6.458 |       6.647 |           5 |
|     1 |   23713 | named_none_projected |    2563736 |           1 | 1301781980348 |       6.132 |          6.480 |       6.843 |           5 |
|     1 |   23713 | attributes_string    |    2563736 |           1 | 712833807     |       8.254 |          8.310 |       8.571 |           5 |
|     1 |   23713 | map_gene             |    2563736 |           1 | 422264        |       8.997 |          9.082 |       9.375 |           5 |
|     1 |   23713 | map_three            |    2563736 |           1 | 77807365      |       9.419 |          9.456 |       9.635 |           5 |
|     1 |   23713 | named_three          |    2563736 |           1 | 77807365      |       7.062 |          7.139 |       8.143 |           5 |
|     1 |   23713 | named_one_projected  |    2563736 |           1 | 1714086       |       6.477 |          6.687 |       7.909 |           5 |
|     2 |   43922 | plain_columns        |    5263269 |           1 | 2542095936640 |      12.271 |         12.661 |      12.994 |           5 |
|     2 |   43922 | named_none_projected |    5263269 |           1 | 2542095936640 |      12.435 |         13.056 |      13.404 |           5 |
|     2 |   43922 | attributes_string    |    5263269 |           1 | 1440737469    |      16.549 |         16.628 |      17.210 |           5 |
|     2 |   43922 | map_gene             |    5263269 |           1 | 849114        |      18.017 |         18.145 |      18.234 |           5 |
|     2 |   43922 | map_three            |    5263269 |           1 | 159110522     |      18.954 |         19.479 |      19.728 |           5 |
|     2 |   43922 | named_three          |    5263269 |           1 | 159110522     |      14.151 |         14.254 |      14.634 |           5 |
|     2 |   43922 | named_one_projected  |    5263269 |           1 | 3461839       |      13.240 |         13.690 |      14.734 |           5 |
|     4 |   71677 | plain_columns        |   10670134 |           1 | 4287043938404 |      24.352 |         25.715 |      27.183 |           5 |
|     4 |   71677 | named_none_projected |   10670134 |           1 | 4287043938404 |      24.065 |         24.349 |      26.493 |           5 |
|     4 |   71677 | attributes_string    |   10670134 |           1 | 2840780134    |      32.355 |         32.545 |      35.403 |           5 |
|     4 |   71677 | map_gene             |   10670134 |           1 | 1829615       |      35.443 |         35.684 |      36.326 |           5 |
|     4 |   71677 | map_three            |   10670134 |           1 | 320600172     |      37.059 |         37.267 |      37.912 |           5 |
|     4 |   71677 | named_three          |   10670134 |           1 | 320600172     |      27.852 |         27.974 |      30.172 |           5 |
|     4 |   71677 | named_one_projected  |   10670134 |           1 | 7548571       |      25.961 |         27.358 |      27.700 |           5 |

Row-level parity of the three named columns against `attributes_map`, asserted before timing:

| parts | records |     rows | mismatches |
|------:|--------:|---------:|-----------:|
|     1 |   23713 |  2563736 |          0 |
|     2 |   43922 |  5263269 |          0 |
|     4 |   71677 | 10670134 |          0 |

Growth with input size: the ratio of the two- and four-part medians to the one-part median, beside the ratios of feature rows and input bytes (release parts differ in size and feature density, so linear growth means the time ratio tracks those ratios):

| workload      | seconds_1_part | time_ratio_2_parts | rows_ratio_2_parts | bytes_ratio_2_parts | time_ratio_4_parts | rows_ratio_4_parts | bytes_ratio_4_parts |
|:--------------|---------------:|-------------------:|-------------------:|--------------------:|-------------------:|-------------------:|--------------------:|
| named_three   |          7.139 |              1.997 |              2.053 |               2.003 |              3.918 |              4.162 |                3.89 |
| map_three     |          9.456 |              2.060 |              2.053 |               2.003 |              3.941 |              4.162 |                3.89 |
| plain_columns |          6.458 |              1.961 |              2.053 |               2.003 |              3.982 |              4.162 |                3.89 |

Memory over the same fresh-process observations, spill forbidden. The peak-RSS ceiling, declared before measuring, is 124.1 MiB: 39 MiB fixed overhead (median of three empty runs) plus three times a live state of 28.4 MiB (four times the largest record, 7,176,865 bytes, plus one 1 MiB output chunk). The render stops if any observation exceeds it or spills:

| parts | workload             | observations | min_max_rss_mib | median_max_rss_mib | max_max_rss_mib | max_duckdb_peak_buffer_mib | max_duckdb_peak_temp_bytes |
|------:|:---------------------|-------------:|----------------:|-------------------:|----------------:|---------------------------:|---------------------------:|
|     1 | plain_columns        |            5 |            46.0 |               46.2 |            46.3 |                       0.25 |                          0 |
|     1 | named_none_projected |            5 |            46.0 |               46.2 |            46.2 |                       0.25 |                          0 |
|     1 | attributes_string    |            5 |            48.4 |               48.4 |            49.6 |                       2.56 |                          0 |
|     1 | map_gene             |            5 |            48.1 |               48.1 |            48.1 |                       1.66 |                          0 |
|     1 | map_three            |            5 |            48.6 |               48.7 |            48.7 |                       1.62 |                          0 |
|     1 | named_three          |            5 |            46.9 |               47.1 |            47.1 |                       0.49 |                          0 |
|     1 | named_one_projected  |            5 |            46.4 |               46.4 |            46.6 |                       0.25 |                          0 |
|     2 | plain_columns        |            5 |            47.7 |               47.9 |            47.9 |                       0.25 |                          0 |
|     2 | named_none_projected |            5 |            47.9 |               47.9 |            47.9 |                       0.25 |                          0 |
|     2 | attributes_string    |            5 |            49.8 |               49.8 |            49.9 |                       2.56 |                          0 |
|     2 | map_gene             |            5 |            50.0 |               50.2 |            50.3 |                       3.34 |                          0 |
|     2 | map_three            |            5 |            50.3 |               50.9 |            51.1 |                       3.28 |                          0 |
|     2 | named_three          |            5 |            48.9 |               48.9 |            48.9 |                       0.49 |                          0 |
|     2 | named_one_projected  |            5 |            48.0 |               48.1 |            48.1 |                       0.25 |                          0 |
|     4 | plain_columns        |            5 |            49.0 |               49.0 |            49.1 |                       0.25 |                          0 |
|     4 | named_none_projected |            5 |            49.0 |               49.1 |            49.1 |                       0.25 |                          0 |
|     4 | attributes_string    |            5 |            51.0 |               51.0 |            51.2 |                       2.56 |                          0 |
|     4 | map_gene             |            5 |            51.1 |               51.2 |            51.3 |                       3.34 |                          0 |
|     4 | map_three            |            5 |            52.3 |               52.3 |            52.5 |                       3.28 |                          0 |
|     4 | named_three          |            5 |            49.7 |               49.9 |            49.9 |                       0.49 |                          0 |
|     4 | named_one_projected  |            5 |            48.9 |               49.1 |            49.3 |                       0.25 |                          0 |

`named_none_projected` requests three keys at bind but reads only the
coordinate columns; its scan resolves no qualifiers. `named_one_projected`
requests three keys and scans only `gene`. `map_three` and `named_three` are
the same metric through the two representations. Comparisons across workloads
reflect different aggregation and expression costs, not isolated lookup cycles.
The nearest recorded GenBank measurements are in
[`benchmark_genbank_reader.md`](benchmark_genbank_reader.md), which times the
plain-column and MAP scans of a single E. coli record without named columns
and is not directly comparable because its input and workloads differ.

## Change from the previous render

The previous render, b6f799e8, measured the same workloads on the
same inputs with the same driver, host and CPU (19), resolving the
named qualifiers once per output row. The cache saves work only on the rows after a feature’s first
location segment, so the workloads that read named columns are the ones it can
change. Host load differs between renders. Median seconds of five repetitions:

| parts | workload             | seconds_before | seconds_after | time_ratio |
|------:|:---------------------|---------------:|--------------:|-----------:|
|     1 | plain_columns        |          6.509 |         6.458 |      0.992 |
|     1 | named_none_projected |          6.004 |         6.480 |      1.079 |
|     1 | attributes_string    |          8.236 |         8.310 |      1.009 |
|     1 | map_gene             |          8.912 |         9.082 |      1.019 |
|     1 | map_three            |          9.359 |         9.456 |      1.010 |
|     1 | named_three          |          7.128 |         7.139 |      1.001 |
|     1 | named_one_projected  |          6.678 |         6.687 |      1.001 |
|     2 | plain_columns        |         12.421 |        12.661 |      1.019 |
|     2 | named_none_projected |         12.111 |        13.056 |      1.078 |
|     2 | attributes_string    |         16.227 |        16.628 |      1.025 |
|     2 | map_gene             |         18.002 |        18.145 |      1.008 |
|     2 | map_three            |         18.756 |        19.479 |      1.039 |
|     2 | named_three          |         14.094 |        14.254 |      1.011 |
|     2 | named_one_projected  |         13.082 |        13.690 |      1.046 |
|     4 | plain_columns        |         23.671 |        25.715 |      1.086 |
|     4 | named_none_projected |         23.878 |        24.349 |      1.020 |
|     4 | attributes_string    |         31.522 |        32.545 |      1.032 |
|     4 | map_gene             |         34.779 |        35.684 |      1.026 |
|     4 | map_three            |         36.601 |        37.267 |      1.018 |
|     4 | named_three          |         27.319 |        27.974 |      1.024 |
|     4 | named_one_projected  |         25.233 |        27.358 |      1.084 |
