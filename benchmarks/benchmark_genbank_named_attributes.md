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
taskset -c 16-19 Rscript -e 'rmarkdown::render("benchmarks/benchmark_genbank_named_attributes.Rmd")'
```

Source revision: 3d7fb282cc9d9f79c187ddf573c1eb5b3cee1a5c; measured `src` tree: 9ada6d903f188bc52198290590ce9b0599c765ca (compare with `git rev-parse HEAD:src` at the reviewed head, since committing this render moves the commit hash but not the tree). DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI).
Extension SHA-256: 0682f5dded2dda4a5a3a720e615769d250029ad1a02449cb5d192ae959e61ce5.
Host load before the run: 1.29, 1.59, 1.49. One thread, affinity 19, 5 timed fresh-process repetitions per cell; `seconds` is the profiler’s query latency, and every observation is in [`benchmark_genbank_named_attributes_runs.csv`](benchmark_genbank_named_attributes_runs.csv).

Inputs, with the registered byte sizes and the feature rows each emits:

| parts | records |     rows |      bytes |
|------:|--------:|---------:|-----------:|
|     1 |   23713 |  2563736 | 1693174640 |
|     2 |   43922 |  5263269 | 3392019192 |
|     4 |   71677 | 10670134 | 6585648536 |

| parts | records | workload             | input_rows | output_rows | metric        | min_seconds | median_seconds | max_seconds | repetitions |
|------:|--------:|:---------------------|-----------:|------------:|:--------------|------------:|---------------:|------------:|------------:|
|     1 |   23713 | plain_columns        |    2563736 |           1 | 1301781980348 |       6.124 |          6.195 |       6.536 |           5 |
|     1 |   23713 | named_none_projected |    2563736 |           1 | 1301781980348 |       6.053 |          6.364 |       6.693 |           5 |
|     1 |   23713 | attributes_string    |    2563736 |           1 | 712833807     |       8.244 |          8.289 |       8.591 |           5 |
|     1 |   23713 | map_gene             |    2563736 |           1 | 422264        |       9.173 |          9.205 |       9.516 |           5 |
|     1 |   23713 | map_three            |    2563736 |           1 | 77807365      |       9.751 |          9.776 |       9.935 |           5 |
|     1 |   23713 | named_three          |    2563736 |           1 | 77807365      |       7.062 |          7.090 |       7.148 |           5 |
|     1 |   23713 | named_one_projected  |    2563736 |           1 | 1714086       |       6.433 |          6.457 |       6.556 |           5 |
|     2 |   43922 | plain_columns        |    5263269 |           1 | 2542095936640 |      12.201 |         12.551 |      12.860 |           5 |
|     2 |   43922 | named_none_projected |    5263269 |           1 | 2542095936640 |      12.124 |         12.289 |      13.036 |           5 |
|     2 |   43922 | attributes_string    |    5263269 |           1 | 1440737469    |      16.593 |         16.644 |      17.094 |           5 |
|     2 |   43922 | map_gene             |    5263269 |           1 | 849114        |      18.510 |         18.595 |      19.448 |           5 |
|     2 |   43922 | map_three            |    5263269 |           1 | 159110522     |      19.180 |         19.706 |      20.150 |           5 |
|     2 |   43922 | named_three          |    5263269 |           1 | 159110522     |      14.265 |         14.296 |      14.529 |           5 |
|     2 |   43922 | named_one_projected  |    5263269 |           1 | 3461839       |      12.855 |         13.241 |      14.254 |           5 |
|     4 |   71677 | plain_columns        |   10670134 |           1 | 4287043938404 |      23.485 |         23.719 |      25.329 |           5 |
|     4 |   71677 | named_none_projected |   10670134 |           1 | 4287043938404 |      23.721 |         24.302 |      25.474 |           5 |
|     4 |   71677 | attributes_string    |   10670134 |           1 | 2840780134    |      31.989 |         32.339 |      33.055 |           5 |
|     4 |   71677 | map_gene             |   10670134 |           1 | 1829615       |      36.139 |         36.266 |      38.272 |           5 |
|     4 |   71677 | map_three            |   10670134 |           1 | 320600172     |      37.995 |         38.158 |      38.573 |           5 |
|     4 |   71677 | named_three          |   10670134 |           1 | 320600172     |      27.753 |         27.877 |      28.690 |           5 |
|     4 |   71677 | named_one_projected  |   10670134 |           1 | 7548571       |      25.301 |         25.618 |      25.653 |           5 |

Row-level parity of the three named columns against `attributes_map`, asserted before timing:

| parts | records |     rows | mismatches |
|------:|--------:|---------:|-----------:|
|     1 |   23713 |  2563736 |          0 |
|     2 |   43922 |  5263269 |          0 |
|     4 |   71677 | 10670134 |          0 |

Growth with input size: the ratio of the two- and four-part medians to the one-part median, beside the ratios of feature rows and input bytes (release parts differ in size and feature density, so linear growth means the time ratio tracks those ratios):

| workload      | seconds_1_part | time_ratio_2_parts | rows_ratio_2_parts | bytes_ratio_2_parts | time_ratio_4_parts | rows_ratio_4_parts | bytes_ratio_4_parts |
|:--------------|---------------:|-------------------:|-------------------:|--------------------:|-------------------:|-------------------:|--------------------:|
| named_three   |          7.090 |              2.016 |              2.053 |               2.003 |              3.932 |              4.162 |                3.89 |
| map_three     |          9.776 |              2.016 |              2.053 |               2.003 |              3.903 |              4.162 |                3.89 |
| plain_columns |          6.195 |              2.026 |              2.053 |               2.003 |              3.829 |              4.162 |                3.89 |

Memory over the same fresh-process observations, spill forbidden. The peak-RSS ceiling, declared before measuring, is 124.7 MiB: 39.6 MiB fixed overhead (median of three empty runs) plus three times a live state of 28.4 MiB (four times the largest record, 7,176,865 bytes, plus one 1 MiB output chunk). The render stops if any observation exceeds it or spills:

| parts | workload             | observations | min_max_rss_mib | median_max_rss_mib | max_max_rss_mib | max_duckdb_peak_buffer_mib | max_duckdb_peak_temp_bytes |
|------:|:---------------------|-------------:|----------------:|-------------------:|----------------:|---------------------------:|---------------------------:|
|     1 | plain_columns        |            5 |            46.9 |               46.9 |            47.1 |                       0.25 |                          0 |
|     1 | named_none_projected |            5 |            46.9 |               47.1 |            47.1 |                       0.25 |                          0 |
|     1 | attributes_string    |            5 |            47.7 |               47.9 |            47.9 |                       2.56 |                          0 |
|     1 | map_gene             |            5 |            48.4 |               48.6 |            48.7 |                       1.66 |                          0 |
|     1 | map_three            |            5 |            49.2 |               49.2 |            49.2 |                       1.62 |                          0 |
|     1 | named_three          |            5 |            46.8 |               46.8 |            47.0 |                       0.49 |                          0 |
|     1 | named_one_projected  |            5 |            46.5 |               46.5 |            46.5 |                       0.25 |                          0 |
|     2 | plain_columns        |            5 |            48.1 |               48.2 |            48.3 |                       0.25 |                          0 |
|     2 | named_none_projected |            5 |            48.1 |               48.1 |            48.2 |                       0.25 |                          0 |
|     2 | attributes_string    |            5 |            49.9 |               50.1 |            50.7 |                       2.56 |                          0 |
|     2 | map_gene             |            5 |            50.6 |               50.7 |            50.8 |                       3.34 |                          0 |
|     2 | map_three            |            5 |            51.0 |               51.2 |            51.8 |                       3.28 |                          0 |
|     2 | named_three          |            5 |            48.6 |               48.7 |            48.7 |                       0.49 |                          0 |
|     2 | named_one_projected  |            5 |            48.2 |               48.2 |            48.2 |                       0.25 |                          0 |
|     4 | plain_columns        |            5 |            49.1 |               49.2 |            49.3 |                       0.25 |                          0 |
|     4 | named_none_projected |            5 |            49.3 |               49.3 |            49.3 |                       0.25 |                          0 |
|     4 | attributes_string    |            5 |            51.0 |               51.1 |            51.2 |                       2.56 |                          0 |
|     4 | map_gene             |            5 |            51.4 |               52.0 |            52.0 |                       3.34 |                          0 |
|     4 | map_three            |            5 |            52.2 |               52.2 |            52.8 |                       3.28 |                          0 |
|     4 | named_three          |            5 |            49.7 |               49.7 |            49.8 |                       0.49 |                          0 |
|     4 | named_one_projected  |            5 |            49.5 |               49.6 |            49.6 |                       0.25 |                          0 |

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

The previous render, dcda99c7, measured the same workloads on the
same inputs with the same driver, host and CPU (19), before the
attributes MAP was written through the shared row writer. The workloads that
build `attributes_map` (`map_gene`, `map_three`) are the ones the writer can
change; the others are the control. Host load differs between renders. Median seconds of five repetitions:

| parts | workload             | seconds_before | seconds_after | time_ratio |
|------:|:---------------------|---------------:|--------------:|-----------:|
|     1 | plain_columns        |          6.458 |         6.195 |      0.959 |
|     1 | named_none_projected |          6.480 |         6.364 |      0.982 |
|     1 | attributes_string    |          8.310 |         8.289 |      0.997 |
|     1 | map_gene             |          9.082 |         9.205 |      1.014 |
|     1 | map_three            |          9.456 |         9.776 |      1.034 |
|     1 | named_three          |          7.139 |         7.090 |      0.993 |
|     1 | named_one_projected  |          6.687 |         6.457 |      0.966 |
|     2 | plain_columns        |         12.661 |        12.551 |      0.991 |
|     2 | named_none_projected |         13.056 |        12.289 |      0.941 |
|     2 | attributes_string    |         16.628 |        16.644 |      1.001 |
|     2 | map_gene             |         18.145 |        18.595 |      1.025 |
|     2 | map_three            |         19.479 |        19.706 |      1.012 |
|     2 | named_three          |         14.254 |        14.296 |      1.003 |
|     2 | named_one_projected  |         13.690 |        13.241 |      0.967 |
|     4 | plain_columns        |         25.715 |        23.719 |      0.922 |
|     4 | named_none_projected |         24.349 |        24.302 |      0.998 |
|     4 | attributes_string    |         32.545 |        32.339 |      0.994 |
|     4 | map_gene             |         35.684 |        36.266 |      1.016 |
|     4 | map_three            |         37.267 |        38.158 |      1.024 |
|     4 | named_three          |         27.974 |        27.877 |      0.997 |
|     4 | named_one_projected  |         27.358 |        25.618 |      0.936 |
