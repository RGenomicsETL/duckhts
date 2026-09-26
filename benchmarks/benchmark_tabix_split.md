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
projection selects `seqname`, `start` and `end`; `gff_five` selects `seqname`,
`start`, `end`, `feature` and `strand`. `gff_all_nine` and `gtf_all_nine`
select all nine physical columns. The BED-like query selects all 16 fields.
This derived file tests generic column projection, not BED semantics or the
BED reader.

## Reproduction

The inputs are registered under workload `tabix-single-split` in
`r/duckhtsbench/inst/benchmark_registry.tsv`. Stage with network access only
in the explicit staging command; `--offline` validates both source checksums
and reconstructs the derived input from an already-cached GFF3:

``` sh
Rscript scripts/stage_tabix_split.R
Rscript scripts/stage_tabix_split.R --offline
```

Build origin/develop from a source archive in a separate directory with
`make configure && make release -j4` (copy the pinned extension-ci-tools
submodule into that directory before building).
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
| branch  | 7f02880d99f48c7873d2c85f394870a343e51e4e | /root/duckhts-tabsplit/build/release/duckhts.duckdb_extension         |

Release binaries (built from the named revisions).

    #> Host: Ubuntu-2404-noble-amd64-base ; affinity: pid 1539597's current affinity list: 0-3 ; candidate src tree: fc4e3615ee77d4b491eec438666b8f15c1c35ed6

| workload     | build   | threads | effective_threads | input_bytes | input_rows | result_rows | median_seconds |
|:-------------|:--------|:--------|------------------:|------------:|:-----------|------------:|---------------:|
| gff_plain    | develop | 1       |                 1 |    24025067 | 1299172    |           1 |          0.539 |
| gff_plain    | branch  | 1       |                 1 |    24025067 | 1299172    |           1 |          0.521 |
| gff_five     | develop | 1       |                 1 |    24025067 | 1299172    |           1 |          0.595 |
| gff_five     | branch  | 1       |                 1 |    24025067 | 1299172    |           1 |          0.604 |
| gff_all_nine | develop | 1       |                 1 |    24025067 | 1299172    |           1 |          0.998 |
| gff_all_nine | branch  | 1       |                 1 |    24025067 | 1299172    |           1 |          0.821 |
| gtf_all_nine | develop | 1       |                 1 |    19833013 | 1302163    |           1 |          0.912 |
| gtf_all_nine | branch  | 1       |                 1 |    19833013 | 1302163    |           1 |          0.742 |
| tabix_bed16  | develop | 1       |                 1 |    97197357 | 100000     |           1 |          0.182 |
| tabix_bed16  | branch  | 1       |                 1 |    97197357 | 100000     |           1 |          0.071 |
| gff_plain    | develop | default |                20 |    24025067 | 1299172    |           1 |          0.519 |
| gff_plain    | branch  | default |                20 |    24025067 | 1299172    |           1 |          0.514 |
| gff_five     | develop | default |                20 |    24025067 | 1299172    |           1 |          0.610 |
| gff_five     | branch  | default |                20 |    24025067 | 1299172    |           1 |          0.602 |
| gff_all_nine | develop | default |                20 |    24025067 | 1299172    |           1 |          0.975 |
| gff_all_nine | branch  | default |                20 |    24025067 | 1299172    |           1 |          0.828 |
| gtf_all_nine | develop | default |                20 |    19833013 | 1302163    |           1 |          1.030 |
| gtf_all_nine | branch  | default |                20 |    19833013 | 1302163    |           1 |          0.759 |
| tabix_bed16  | develop | default |                20 |    97197357 | 100000     |           1 |          0.188 |
| tabix_bed16  | branch  | default |                20 |    97197357 | 100000     |           1 |          0.072 |

Five-run median wall time (seconds); output fingerprints and hash sums match across builds.

## Narrow GFF profile

Two additional five-run passes at affinity 0–3 measured `gff_plain` at
0.711/0.772 s and 0.573/0.643 s (develop/branch) with one thread; the second
pass also measured 0.522/0.538 s with default threads. The narrow projection
is slower in both one-thread passes. A separate paired `perf stat`/`perf record`
run pinned to CPU 0 used the same query and input, with one warm-up and 12
repeated scans per build. Develop used 7.705 s CPU task-clock and 92.99 billion
core instructions; the branch used 7.471 s and 89.47 billion. Sampled core
cycles were about 59% in zlib in each build. The branch spent 3.5% in libc
`strcspn` (not present in the develop profile); sampled reader share was
17.3% versus 14.3% and libc share was 11.8% versus 15.3%. Develop’s narrow reader
uses an inline byte scan per projected field; the branch pays for five
`strcspn` calls and field-span storage, with little redundant scanning to
eliminate. Total core cycles were 33.02 billion and 33.55 billion,
respectively. The binaries were stripped, so samples inside the reader cannot
be assigned to individual C lines. CPU-0 task-clock does not reproduce the
0–3 wall-time regression, and the first pass had competing jobs on the
pinned CPUs. The remaining difference cannot be apportioned between parser
cost, scheduling and frequency variation; no narrow-projection speedup is
claimed.

The raw input checksum, byte count and derivation are pinned in the registry
and validated before measurement. `input_rows` is the count returned by each
scan; `result_rows` counts aggregate rows materialized by the query. Times
include scanning, decompression of GENCODE inputs, aggregation and DuckDB
execution but exclude opening the connection and loading the extension.
Fingerprint and hash-sum equality check output values, not ordering; no
output rows are filtered by the benchmark.
