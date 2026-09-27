DuckHTS vs GFFBase GFF3 Conformance Benchmark
================

<!-- benchmark_gffbase_conformance.md is generated from benchmark_gffbase_conformance.Rmd. -->

# Benchmark

This benchmark compares DuckHTS `read_gff(...)` with
[GFFBase](https://github.com/Kuanhao-Chao/gffbase) across parser,
conformance, scan and feature-database workloads:

1.  **GFFBase audit:** record whether GFFBase is using its Rust parser
    or Python fallback, check Rust/Python parser parity when both are
    available, and verify basic GFF3/GTF dialect detection. We do not
    treat marketing claims as proof.
2.  **GFF3 conformance:** compare DuckHTS default permissive
    `read_gff(...)` and `read_gff(..., strict := true)` against observed
    GFFBase strict-parser behavior. The cases are adapted from GFFBase’s
    NCBI/GFF3 compliance tests at commit
    `78714cf30a9d799eab544e00a79a4da9754987ca`.
3.  **Direct scan/parser throughput:** DuckHTS scans the same GFF3 files
    through `read_gff(...)`; GFFBase parses them with
    `parse_gff(strict = TRUE)`.
4.  **Feature-database core:** SQL over `read_gff(...)` materializes
    features, parent edges and descendants for comparison with
    `gffbase.create_db`.

DuckHTS supports GFF3 through `read_gff(...)`, GTF through
`read_gtf(...)`, and plain tabix-style GFF/GTF-like rows through
`read_tabix(...)`. The new `strict := true` option is specifically GFF3
structural validation for `read_gff(...)`; it is not a full
feature-database hierarchy validator.

# Run

Stage GFFBase into the DuckHTS cache, then render this report:

``` sh
make stage-gffbase
make bench-gffbase
```

Override defaults with `GFFBASE_BENCH_ROWS`, `GFFBASE_BENCH_PASSES`,
`GFFBASE_BENCH_FORCE=1`, `GFFBASE_BENCH_INCLUDE_CREATE_DB=1`, or
`GFFBASE_HUMAN_GFF=/path/to/gencode.gff3.gz[,/path/to/refseq.gff3.gz]`
for real human-scale files.

# Configuration

| parameter               | value                                  |
|:------------------------|:---------------------------------------|
| DuckHTS git rev         | 27e47c86d052+dirty                     |
| DuckHTS extension       | build/release/duckhts.duckdb_extension |
| GFFBase version         | 0.1.0                                  |
| GFFBase native parser   | TRUE                                   |
| GFFBase upstream commit | 78714cf30a9d                           |
| DuckDB Python           | 1.5.2                                  |
| server hostname         | Ubuntu-2404-noble-amd64-base           |
| server OS               | Ubuntu 24.04.3 LTS                     |
| kernel / machine        | 6.8.0-78-generic / x86_64              |
| CPU model               | 13th Gen Intel(R) Core(TM) i5-13500    |
| CPU logical / affinity  | 20 / 20                                |
| memory                  | 62.58 GiB                              |
| DuckDB threads          | 4                                      |
| synthetic rows          | 200,000                                |
| timed passes            | 3                                      |

# GFFBase audit

| check                     | status | detail                                                                                                                                                                                                                                                                       | value |
|:--------------------------|:-------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------|
| import                    | ok     | /root/.cache/duckhts/benchmarks/gffbase/python-site/gffbase/**init**.py                                                                                                                                                                                                      | 0.1.0 |
| native_available          | ok     | auto engine resolves to rust                                                                                                                                                                                                                                                 | True  |
| rust_python_parser_parity | ok     | 25/25 strict cases matched by status,row_count,error_kind                                                                                                                                                                                                                    | 0     |
| detect_dialect_gff3       | ok     | {“field separator”: “;”, “fmt”: “gff3”, “keyval separator”: “=”, “leading semicolon”: false, “multival separator”: “,”, “order”: \[“ID”, “Name”\], “quoted GFF2 values”: false, “repeated keys”: false, “semicolon in quotes”: false, “trailing semicolon”: false}           | gff3  |
| detect_dialect_gtf        | ok     | {“field separator”: “;”, “fmt”: “gtf”, “keyval separator”: ” “,”leading semicolon”: false, “multival separator”: “,”, “order”: \[“gene_id”, “transcript_id”\], “quoted GFF2 values”: true, “repeated keys”: false, “semicolon in quotes”: false, “trailing semicolon”: true} | gtf   |

# GFF3 conformance summary

| metric                                          | value |
|:------------------------------------------------|------:|
| cases                                           |    25 |
| GFFBase matched expected strict behavior        |    25 |
| DuckHTS default matched GFFBase strict behavior |     7 |
| DuckHTS default gaps vs GFFBase strict behavior |    18 |
| DuckHTS strict matched GFFBase strict behavior  |    25 |
| DuckHTS strict gaps vs GFFBase strict behavior  |     0 |

## Remaining DuckHTS strict gaps

| note                                                                                                          |
|:--------------------------------------------------------------------------------------------------------------|
| No remaining gaps in the local strict-conformance cases when DuckHTS uses attributes_list / attributes_pairs. |

# Direct scan/parser throughput

| dataset               | tool    | variant                            | rows    | passes | median_sec | rows_per_sec | mb_per_sec | vs_gffbase_parse |
|:----------------------|:--------|:-----------------------------------|:--------|-------:|-----------:|-------------:|-----------:|-----------------:|
| synthetic_200000.gff3 | DuckHTS | read_gff COUNT(\*)                 | 200,000 |      3 |     0.0108 |   18568013.9 |     2577.6 |            51.85 |
| synthetic_200000.gff3 | DuckHTS | read_gff strict COUNT(\*)          | 200,000 |      3 |     0.0372 |    5379576.4 |      746.8 |            15.02 |
| synthetic_200000.gff3 | DuckHTS | read_gff projected filter/sum      | 200,000 |      3 |     0.0260 |    7706230.9 |     1069.8 |            21.52 |
| synthetic_200000.gff3 | DuckHTS | read_gff attributes_map            | 200,000 |      3 |     0.0999 |    2002228.9 |      277.9 |             5.59 |
| synthetic_200000.gff3 | DuckHTS | read_gff attributes_list           | 200,000 |      3 |     0.1726 |    1158461.4 |      160.8 |             3.24 |
| synthetic_200000.gff3 | DuckHTS | read_gff attributes_pairs          | 200,000 |      3 |     0.1213 |    1649459.0 |      229.0 |             4.61 |
| synthetic_200000.gff3 | DuckHTS | read_gff attributes_map+list+pairs | 200,000 |      3 |     0.3728 |     536488.2 |       74.5 |             1.50 |
| synthetic_200000.gff3 | GFFBase | parse_gff strict (rust)            | 200,000 |      3 |     0.5585 |     358102.1 |       49.7 |             1.00 |

# Interpretation

- GFFBase should be read as a comparison implementation, not an
  unquestioned oracle. The audit table exposes whether `auto` used Rust
  or Python and whether the two parsers agree on the local strict cases.
- DuckHTS default reading remains permissive for existing ingestion
  workflows. `strict := true` now rejects the structural GFF3 failures
  covered here.
- `attributes_list := true` and `attributes_pairs := true` close the
  local strict-conformance attribute gaps by representing
  repeated/comma-split GFF3 values and URL-decoded values. The older
  `attributes_map := true` remains a backward-compatible scalar
  convenience column and intentionally cannot express duplicate keys or
  multi-valued attributes losslessly.
- Requesting `attributes_map`, `attributes_list`, and `attributes_pairs`
  together materializes three nested representations independently. That
  is useful as a stress test and compatibility escape hatch, but callers
  should normally choose the one representation their query needs.
- The timing table compares direct scans/parsing, not GFFBase’s full
  persistent feature database API. Set `GFFBASE_HUMAN_GFF` to real
  GENCODE/RefSeq/MANE files for human-scale rows, and enable
  `GFFBASE_BENCH_INCLUDE_CREATE_DB=1` when the workload of interest is
  database materialization plus hierarchy indexing.

# Exercise the thesis: a feature database from `read_gff`

DuckHTS observes GFF3 records; SQL builds the persistent,
position-sorted feature table (including ordered, repeated attributes),
parent edges and transitive closure. This section uses the registered
MANE and GENCODE inputs and GFFBase v0.2.1; the conformance and parser
sections above use v0.1.0.

``` sql
CREATE TABLE features AS
WITH input AS (
  SELECT row_number() OVER () AS rid, ID, Parent, seqname, source, feature,
         start, "end", score, strand, frame, attributes_pairs
  FROM read_gff('{input}', attributes_pairs := true, attributes := ['ID', 'Parent'])
), numbered AS (
  SELECT *, row_number() OVER (PARTITION BY ID ORDER BY rid) - 1 AS occurrence
  FROM input
)
SELECT rid,
       CASE WHEN ID IS NULL THEN 'row:' || rid
            WHEN occurrence = 0 THEN ID
            ELSE ID || '_' || occurrence END AS id,
       seqname AS seqid, source, feature AS featuretype, start, "end", score,
       strand, frame, Parent AS parents, attributes_pairs AS attributes
FROM numbered
ORDER BY seqid, start, "end";

CREATE TABLE relations AS
SELECT DISTINCT f.id AS child, trim(p) AS parent, 1 AS level
FROM features f, unnest(string_split(f.parents, ',')) AS u(p)
WHERE f.parents IS NOT NULL;

CREATE TABLE closure AS
WITH RECURSIVE c(ancestor, descendant, level) AS (
  SELECT parent, child, level FROM relations
  UNION ALL
  SELECT c.ancestor, r.child, c.level + 1
  FROM c JOIN relations r ON r.parent = c.descendant
)
SELECT DISTINCT ancestor, descendant, level FROM c;

CREATE INDEX closure_ancestor ON closure(ancestor);
CHECKPOINT;
```

`gffbase.create_db` uses `force=True`, `merge_strategy="create_unique"`,
`force_gff=True`, `pragmas={"threads": threads}` and `validation=None`.
Validation modes are not run. SQL assigns suffixes to repeated physical
IDs in source order to match `create_unique`; both sides retain every
input row.

## Parity before timing

Seeded genes and windows are identical across engines; sorted descendant
rows compare anchor, descendant ID, chromosome, source, type,
coordinates, score, strand, frame, source order and depth (missing
numeric score is rendered as `'.'` on the SQL side). Region hits compare
sorted feature IDs for every window, not just the total count.

| dataset                | features | anchors | windows | descendants | region_hits | parity | passes |
|:-----------------------|---------:|--------:|--------:|------------:|------------:|:-------|-------:|
| mane_v15_ensembl_gff3  |   524834 |    5000 |    1000 |      129527 |       31233 | PASS   |      3 |
| gencode_v49_basic_gff3 |  5866158 |    5000 |    1000 |      375777 |      101445 | PASS   |      1 |

## Materialization and queries

Wall time and peak resident memory include process startup, build, and
both queries; build time covers database creation only. Each pass starts
with a new process and database. File sizes follow a closed checkpoint.
The GENCODE GFFBase process runs in a `systemd-run` scope with memory
cap 45G. GENCODE has a single pass, so it supplies no run-to-run
variance estimate.

| dataset                | engine  | passes | build_s | wall_s | peak_rss_GiB | database_GiB | descendants_s | windows_s |
|:-----------------------|:--------|-------:|--------:|-------:|-------------:|-------------:|--------------:|----------:|
| mane_v15_ensembl_gff3  | DuckHTS |      3 |    3.95 |   5.42 |         1.20 |          0.1 |          0.22 |      0.76 |
| mane_v15_ensembl_gff3  | gffbase |      3 |   20.83 |  23.14 |         1.13 |          0.6 |          0.44 |      1.45 |
| gencode_v49_basic_gff3 | DuckHTS |      1 |   40.80 |  43.53 |         9.42 |          1.0 |          0.46 |      1.00 |
| gencode_v49_basic_gff3 | gffbase |      1 |  314.89 | 319.44 |         9.02 |          6.9 |          0.69 |      2.65 |

| dataset                | input_sha256                                                     | gffbase_memory_cap |
|:-----------------------|:-----------------------------------------------------------------|:-------------------|
| mane_v15_ensembl_gff3  | 69089bbc84d1d3c3ce31c2ed3f85b6c3169fb8836d092a082623c59a43fd22ef | none               |
| gencode_v49_basic_gff3 | 639c16217b6341abe802041118178cd65b176c2cfd1b0c07415deb6b269c6138 | 45G                |

| parameter            | value                                                            |
|:---------------------|:-----------------------------------------------------------------|
| revision             | 27e47c86d052364c235e0a54f134f54de1b7929d                         |
| extension_sha256     | 847efdefbde57dcb00b56dd0df960558160163059bd48f9ddb773cab24dcaf97 |
| host                 | Ubuntu-2404-noble-amd64-base                                     |
| os                   | Linux-6.8.0-78-generic-x86_64-with-glibc2.39                     |
| cpu                  | 13th Gen Intel(R) Core(TM) i5-13500                              |
| threads              | 4                                                                |
| duckdb_version       | 1.5.2                                                            |
| gffbase_version      | 0.2.1                                                            |
| duckhts_memory_limit | 8GB                                                              |
| validation           | not run                                                          |

This is the storage and query core, not GFFBase’s validation modes, GTF
synthesis, gffutils API or R-tree. The recorded DuckHTS batched
descendant query is faster: its one SQL join handles all anchors;
GFFBase also batches through its dedicated `children_batched` API. This
result does not describe per-gene SQL lookups, where GFFBase’s indexed
batched API can avoid repeated query overhead. DuckDB’s default memory
limit is 80% of RAM; these DuckHTS runs set a limit explicitly. GENCODE
DuckHTS peak RSS exceeds both that limit and GFFBase’s measured RSS
despite a smaller database. The DuckDB limit does not cap all process
allocations.

To restage inputs and the pinned wheel, repeat the measurements and
render:

``` sh
make bench-gffbase-featuredb
```

`GFFBASE_FEATUREDB_PYTHON` selects a Python with DuckDB and PyArrow
installed; `DUCKHTS_CACHE_DIR` selects the registry cache. Set
`GFFBASE_FEATUREDB_MANE_PASSES` and `GFFBASE_FEATUREDB_GENCODE_PASSES`
to increase repetitions. Python dependencies are installed from the
verified wheel without network access during the measurement; staging
may download inputs.
