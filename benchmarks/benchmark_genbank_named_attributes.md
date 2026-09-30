GenBank named-qualifier projections
================

This report measures `read_genbank(..., attributes := [...])` against the two
ways a query could obtain the same qualifiers before it existed: the
`attributes_map` MAP and the raw `attributes` string. The input is NCBI RefSeq
release 237’s plasmid division (workload `genbank-plasmid`): its first
`plasmid.N.genomic.gbff.gz` part alone, its first two parts joined, and its
first four parts joined, each part pinned by NCBI’s published MD5 and each
joined file by its own SHA-256 and record count. Joining consecutive release
parts grows the record count while every other dimension stays what the
release supplies; the `records` column in every table is the LOCUS count of
the input. A record whose only feature is `source` emits no row, so row counts
are below record counts. Stage all four parts without network access at render
time:

``` sh
Rscript -e 'duckhtsbench::duckhts_bench_stage_genbank_plasmid()'
```

Every workload is a DuckDB table macro whose parameters and output columns
carry a `p_` prefix, so no name can resolve to a reader column, and the input
path is bound at execution rather than pasted into SQL. Each workload returns
one row: a row count and a scalar metric over the stated columns. Before any
timing, the render asserts row-level parity: on every input, every named column
equals `attributes_map[key]` on every row, or the render stops. The
`named_three` and `map_three` workloads then compute the same metric over the
same three keys through the two representations, and the render stops if they
disagree, so the timed queries are known to return the same answer. Five timed
repetitions follow one verified warm-up run per workload and input; the
one-part run exceeds the five-second floor `STYLE.md` sets for a timing
verdict. Memory is measured the way `benchmark_genbank_memory.Rmd` does: one
fresh DuckDB CLI process per observation under `/usr/bin/time`, with the JSON
profiler’s peak buffer-manager and temporary-directory bytes, for the three-key
MAP and named-column workloads at every scale, three observations each, with
`max_temp_directory_size` set to zero so any spill is an error. The budget is
the reader’s contract, peak memory proportional to the largest record, not to
records consumed: the ceiling is three times the one-part fresh-process peak
RSS, applied to every observation. The one-thread condition is pinned to the
highest logical CPU in the rendering process’s affinity mask, so wrap the
render in `taskset -c` to choose it; the table records the affinity.
`read_genbank` declares one worker, so a multi-thread condition would time the
same scan and is not recorded.

``` sh
taskset -c 8-11 Rscript -e 'rmarkdown::render("benchmarks/benchmark_genbank_named_attributes.Rmd")'
```

Source revision: 2d01deb2f314e9dba2d30a48b253efa932afae47; measured `src` tree: a6db7399672b87534728a9c354744007fba1c838 (compare with `git rev-parse HEAD:src` at the reviewed head, since committing this render moves the commit hash but not the tree). DuckDB runtime: v1.5.6.
Extension SHA-256: d308cf4c1c615e1172589696f520a6cec6460d76a4574dc940956143ccfc3a93.
Host load before the run: 1.14, 0.76, 0.51. One thread, affinity 11, 5 timed repetitions per cell.

| parts | records | workload             | affinity | input_rows | output_rows | metric        | min_seconds | median_seconds | max_seconds | repetitions |
|------:|--------:|:---------------------|:---------|-----------:|------------:|:--------------|------------:|---------------:|------------:|------------:|
|     1 |   23713 | plain_columns        | 11       |    2563736 |           1 | 1301781980348 |       6.299 |          6.381 |       6.392 |           5 |
|     1 |   23713 | named_none_projected | 11       |    2563736 |           1 | 1301781980348 |       6.203 |          6.235 |       6.407 |           5 |
|     1 |   23713 | attributes_string    | 11       |    2563736 |           1 | 712833807     |       8.270 |          8.325 |       9.180 |           5 |
|     1 |   23713 | map_gene             | 11       |    2563736 |           1 | 422264        |       8.896 |          8.959 |       9.099 |           5 |
|     1 |   23713 | map_three            | 11       |    2563736 |           1 | 77807365      |       9.930 |         10.041 |      10.614 |           5 |
|     1 |   23713 | named_three          | 11       |    2563736 |           1 | 77807365      |       7.549 |          7.630 |       7.674 |           5 |
|     1 |   23713 | named_one_projected  | 11       |    2563736 |           1 | 1714086       |       6.882 |          6.992 |       7.173 |           5 |
|     2 |   43922 | plain_columns        | 11       |    5263269 |           1 | 2542095936640 |      12.966 |         13.370 |      13.653 |           5 |
|     2 |   43922 | named_none_projected | 11       |    5263269 |           1 | 2542095936640 |      12.319 |         12.391 |      12.692 |           5 |
|     2 |   43922 | attributes_string    | 11       |    5263269 |           1 | 1440737469    |      16.110 |         16.121 |      16.145 |           5 |
|     2 |   43922 | map_gene             | 11       |    5263269 |           1 | 849114        |      17.878 |         17.961 |      17.971 |           5 |
|     2 |   43922 | map_three            | 11       |    5263269 |           1 | 159110522     |      18.609 |         18.731 |      19.332 |           5 |
|     2 |   43922 | named_three          | 11       |    5263269 |           1 | 159110522     |      15.050 |         15.109 |      15.191 |           5 |
|     2 |   43922 | named_one_projected  | 11       |    5263269 |           1 | 3461839       |      13.871 |         13.939 |      14.587 |           5 |
|     4 |   71677 | plain_columns        | 11       |   10670134 |           1 | 4287043938404 |      24.090 |         24.826 |      25.864 |           5 |
|     4 |   71677 | named_none_projected | 11       |   10670134 |           1 | 4287043938404 |      24.028 |         24.114 |      24.328 |           5 |
|     4 |   71677 | attributes_string    | 11       |   10670134 |           1 | 2840780134    |      31.227 |         31.331 |      31.369 |           5 |
|     4 |   71677 | map_gene             | 11       |   10670134 |           1 | 1829615       |      34.821 |         34.904 |      35.142 |           5 |
|     4 |   71677 | map_three            | 11       |   10670134 |           1 | 320600172     |      36.227 |         36.375 |      36.506 |           5 |
|     4 |   71677 | named_three          | 11       |   10670134 |           1 | 320600172     |      28.977 |         29.260 |      30.212 |           5 |
|     4 |   71677 | named_one_projected  | 11       |   10670134 |           1 | 7548571       |      25.936 |         26.331 |      26.753 |           5 |

Row-level parity of the three named columns against `attributes_map`, asserted before timing:

| parts | records |     rows | mismatches |
|------:|--------:|---------:|-----------:|
|     1 |   23713 |  2563736 |          0 |
|     2 |   43922 |  5263269 |          0 |
|     4 |   71677 | 10670134 |          0 |

Growth with record count: the ratio of the two- and four-part medians to the one-part median, beside the ratio of their record counts (release parts are not equal in size, so linear growth means the time ratio tracks the record ratio):

| workload      | seconds_1_part | ratio_2_parts | ratio_4_parts | records_ratio_2_parts | records_ratio_4_parts |
|:--------------|---------------:|--------------:|--------------:|----------------------:|----------------------:|
| named_three   |          7.630 |         1.980 |         3.835 |                 1.852 |                 3.023 |
| map_three     |         10.041 |         1.865 |         3.623 |                 1.852 |                 3.023 |
| plain_columns |          6.381 |         2.095 |         3.891 |                 1.852 |                 3.023 |

Memory, one fresh DuckDB CLI process per observation, spill forbidden. The peak-RSS ceiling is 174.3 MiB (three times the one-part named-column median); the render stops if any observation exceeds it or spills:

| parts | records | workload    | input_rows | median_max_rss_mib | max_max_rss_mib | duckdb_peak_buffer_mib | duckdb_peak_temp_bytes | observations |
|------:|--------:|:------------|-----------:|-------------------:|----------------:|-----------------------:|-----------------------:|-------------:|
|     1 |   23713 | map_three   |    2563736 |               58.6 |            59.3 |                   1.62 |                      0 |            3 |
|     1 |   23713 | named_three |    2563736 |               58.1 |            58.1 |                   0.49 |                      0 |            3 |
|     2 |   43922 | map_three   |    5263269 |               55.6 |            59.2 |                   3.28 |                      0 |            3 |
|     2 |   43922 | named_three |    5263269 |               57.5 |            57.5 |                   0.49 |                      0 |            3 |
|     4 |   71677 | map_three   |   10670134 |               60.3 |            60.4 |                   3.28 |                      0 |            3 |
|     4 |   71677 | named_three |   10670134 |               58.0 |            58.0 |                   0.49 |                      0 |            3 |

`named_none_projected` requests three keys at bind but reads only the
coordinate columns; its scan resolves no qualifiers. `named_one_projected`
requests three keys and scans only `gene`. `map_three` and `named_three` are
the same metric through the two representations. Comparisons across workloads
reflect different aggregation and expression costs, not isolated lookup cycles.
The nearest recorded GenBank measurements are in
[`benchmark_genbank_reader.md`](benchmark_genbank_reader.md), which times the
plain-column and MAP scans of a single E. coli record without named columns
and is not directly comparable because its input and workloads differ.
