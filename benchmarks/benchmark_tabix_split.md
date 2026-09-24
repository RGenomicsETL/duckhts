Tabix GFF/GTF single-split reader benchmark
================

This report compares identical release builds of origin/develop and the single-split
reader on the same host. Each query returns one aggregate row but hashes every
requested column so DuckDB must materialize it. The two builds must agree on
row count, XOR fingerprint and hash sum for each workload. Warm-up, load, and staging are
outside the five timed repetitions. The input is GENCODE mouse vM25 BASIC GFF3
and GTF; the BED-like file is the first 100,000 GFF3 data records rearranged
into 16 columns (chromosome, zero-based start, end, feature, score, strand,
source, phase, attributes, then seven repeated source fields). The `gff_plain`
projection selects `seqname`, `start` and `end`; `gff_all_nine` and
`gtf_all_nine` select all nine physical columns. The BED-like query selects
all 16 fields. This derived file tests generic column projection, not BED
semantics or the BED reader.

## Reproduction

The inputs are registered under workload `tabix-single-split` in
`r/duckhtsbench/inst/benchmark_registry.tsv`. Stage with network access only
in the explicit staging command; `--offline` validates both source checksums
and reconstructs the derived input from an already-cached GFF3:

``` sh
Rscript scripts/stage_tabix_split.R
Rscript scripts/stage_tabix_split.R --offline
```

Build origin/develop in a separate worktree with `make configure && make release -j4`.
Set `DUCKHTS_TABIX_BASELINE_EXTENSION` to its release binary and
`DUCKHTS_TABIX_BASELINE_REVISION` to its exact commit; the candidate defaults
to `build/release/duckhts.duckdb_extension`. Render with
`Rscript -e 'rmarkdown::render("benchmarks/benchmark_tabix_split.Rmd")'`.
The `default` condition leaves DuckDB’s `threads` setting untouched. Both
conditions run on the same process affinity (reported below); the reader is
not a parallel table scan, so the thread conditions are not a scaling claim.

    #> Warning: package 'duckdb' was built under R version 4.6.1

| build   | revision                                 | binary                                                                |
|:--------|:-----------------------------------------|:----------------------------------------------------------------------|
| develop | 81edf9bff85ab2daabeccd1daaf13e4bfa39ee47 | /tmp/duckhts-tabsplit-baseline/build/release/duckhts.duckdb_extension |
| branch  | 5cc6914da50aa6fe765415f165d5818de4639e91 | /root/duckhts-tabsplit/build/release/duckhts.duckdb_extension         |

Release binaries (built from the named revisions).

    #> Host: Ubuntu-2404-noble-amd64-base ; affinity: pid 1011537's current affinity list: 0-3 ; candidate src tree: 5b7798c9b48470677a6229fab444565d44f2d28a

| workload     | build   | threads | effective_threads | input_bytes | input_rows | result_rows | median_seconds |
|:-------------|:--------|:--------|------------------:|------------:|:-----------|------------:|---------------:|
| gff_plain    | develop | 1       |                 1 |    24025067 | 1299172    |           1 |          0.513 |
| gff_plain    | branch  | 1       |                 1 |    24025067 | 1299172    |           1 |          0.586 |
| gff_all_nine | develop | 1       |                 1 |    24025067 | 1299172    |           1 |          0.956 |
| gff_all_nine | branch  | 1       |                 1 |    24025067 | 1299172    |           1 |          0.797 |
| gtf_all_nine | develop | 1       |                 1 |    19833013 | 1302163    |           1 |          0.910 |
| gtf_all_nine | branch  | 1       |                 1 |    19833013 | 1302163    |           1 |          0.752 |
| tabix_bed16  | develop | 1       |                 1 |    97197357 | 100000     |           1 |          0.181 |
| tabix_bed16  | branch  | 1       |                 1 |    97197357 | 100000     |           1 |          0.070 |
| gff_plain    | develop | default |                20 |    24025067 | 1299172    |           1 |          0.511 |
| gff_plain    | branch  | default |                20 |    24025067 | 1299172    |           1 |          0.578 |
| gff_all_nine | develop | default |                20 |    24025067 | 1299172    |           1 |          0.949 |
| gff_all_nine | branch  | default |                20 |    24025067 | 1299172    |           1 |          0.796 |
| gtf_all_nine | develop | default |                20 |    19833013 | 1302163    |           1 |          0.907 |
| gtf_all_nine | branch  | default |                20 |    19833013 | 1302163    |           1 |          0.747 |
| tabix_bed16  | develop | default |                20 |    97197357 | 100000     |           1 |          0.183 |
| tabix_bed16  | branch  | default |                20 |    97197357 | 100000     |           1 |          0.070 |

Five-run median wall time (seconds); output fingerprints and hash sums match across builds.

In this run the narrow GFF projection is slower while all-column GFF/GTF
and wide tabix improve. Scanning even the unprojected final GFF attribute
field costs time in the narrow projection; this workload does not establish
a universal speedup.

The raw input checksum, byte count and derivation are pinned in the registry
and validated before measurement. `input_rows` is the count returned by each
scan; `result_rows` counts aggregate rows materialized by the query. Times include scanning, decompression
of GENCODE inputs, aggregation and DuckDB execution but exclude opening the
connection and loading the extension. Fingerprint and hash-sum equality check output values,
not ordering; no output rows are filtered by the benchmark.
