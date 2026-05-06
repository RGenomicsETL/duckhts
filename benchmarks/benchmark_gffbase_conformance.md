DuckHTS vs GFFBase GFF3 Conformance Benchmark
================

<!-- benchmark_gffbase_conformance.md is generated from benchmark_gffbase_conformance.Rmd. -->

# Benchmark

This benchmark compares DuckHTS `read_gff(...)` with
[GFFBase](https://github.com/Kuanhao-Chao/gffbase) around two related
but different workload shapes:

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

DuckHTS supports GFF3 through `read_gff(...)`, GTF through
`read_gtf(...)`, and plain tabix-style GFF/GTF-like rows through
`read_tabix(...)`. The new `strict := true` option is specifically GFF3
structural validation for `read_gff(...)`; it is not a full
feature-database hierarchy validator.

# Run

Install GFFBase into a local target directory, then render this report:

``` sh
python3 -m pip install --target .tmp/gffbase_site gffbase==0.1.0
PYTHONPATH=.tmp/gffbase_site python3 scripts/gffbase_conformance_benchmark.py \
  --extension build/release/duckhts.duckdb_extension \
  --rows 200000 \
  --passes 3 \
  --out-dir .tmp/gffbase_conformance
make bench-gffbase
```

Override defaults with `GFFBASE_BENCH_ROWS`, `GFFBASE_BENCH_PASSES`,
`GFFBASE_BENCH_FORCE=1`, `GFFBASE_BENCH_INCLUDE_CREATE_DB=1`, or
`GFFBASE_HUMAN_GFF=/path/to/gencode.gff3.gz[,/path/to/refseq.gff3.gz]`
for real human-scale files.

# Configuration

| parameter               | value                                                |
|:------------------------|:-----------------------------------------------------|
| DuckHTS git rev         | 78532b859c12                                         |
| DuckHTS extension       | /root/duckhts/build/release/duckhts.duckdb_extension |
| GFFBase version         | 0.1.0                                                |
| GFFBase native parser   | TRUE                                                 |
| GFFBase upstream commit | 78714cf30a9d                                         |
| DuckDB Python           | 1.5.2                                                |
| server hostname         | Ubuntu-2404-noble-amd64-base                         |
| server OS               | Ubuntu 24.04.3 LTS                                   |
| kernel / machine        | 6.8.0-78-generic / x86_64                            |
| CPU model               | 13th Gen Intel(R) Core(TM) i5-13500                  |
| CPU logical / affinity  | 20 / 20                                              |
| memory                  | 62.58 GiB                                            |
| DuckDB threads          | 4                                                    |
| synthetic rows          | 200,000                                              |
| timed passes            | 3                                                    |

# GFFBase audit

| check                     | status | detail                                                                                                                                                                                                                                                                       | value |
|:--------------------------|:-------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------|
| import                    | ok     | /root/duckhts/.tmp/gffbase_site/gffbase/**init**.py                                                                                                                                                                                                                          | 0.1.0 |
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

| dataset                              | tool    | variant                            | rows      | passes | median_sec | rows_per_sec | mb_per_sec | vs_gffbase_parse |
|:-------------------------------------|:--------|:-----------------------------------|:----------|-------:|-----------:|-------------:|-----------:|-----------------:|
| synthetic_200000.gff3                | DuckHTS | read_gff COUNT(\*)                 | 200,000   |      3 |     0.0117 |   17090326.8 |     2372.4 |            46.44 |
| synthetic_200000.gff3                | DuckHTS | read_gff strict COUNT(\*)          | 200,000   |      3 |     0.0637 |    3139960.5 |      435.9 |             8.53 |
| synthetic_200000.gff3                | DuckHTS | read_gff projected filter/sum      | 200,000   |      3 |     0.0329 |    6075425.0 |      843.4 |            16.51 |
| synthetic_200000.gff3                | DuckHTS | read_gff attributes_map            | 200,000   |      3 |     0.1367 |    1462720.5 |      203.1 |             3.98 |
| synthetic_200000.gff3                | DuckHTS | read_gff attributes_list           | 200,000   |      3 |     0.1851 |    1080222.0 |      150.0 |             2.94 |
| synthetic_200000.gff3                | DuckHTS | read_gff attributes_pairs          | 200,000   |      3 |     0.1207 |    1657507.5 |      230.1 |             4.50 |
| synthetic_200000.gff3                | DuckHTS | read_gff attributes_map+list+pairs | 200,000   |      3 |     0.3930 |     508842.8 |       70.6 |             1.38 |
| synthetic_200000.gff3                | GFFBase | parse_gff strict (rust)            | 200,000   |      3 |     0.5435 |     367972.6 |       51.1 |             1.00 |
| gencode.v49.basic.annotation.gff3.gz | DuckHTS | read_gff COUNT(\*)                 | 5,866,158 |      3 |     2.0129 |    2914263.7 |       40.7 |            19.98 |
| gencode.v49.basic.annotation.gff3.gz | DuckHTS | read_gff strict COUNT(\*)          | 5,866,158 |      3 |     4.4549 |    1316797.3 |       18.4 |             9.03 |
| gencode.v49.basic.annotation.gff3.gz | DuckHTS | read_gff projected filter/sum      | 5,866,158 |      3 |     2.9410 |    1994602.0 |       27.9 |            13.68 |
| gencode.v49.basic.annotation.gff3.gz | DuckHTS | read_gff attributes_map            | 5,866,158 |      3 |    10.0948 |     581105.6 |        8.1 |             3.98 |
| gencode.v49.basic.annotation.gff3.gz | DuckHTS | read_gff attributes_list           | 5,866,158 |      3 |    16.9759 |     345558.7 |        4.8 |             2.37 |
| gencode.v49.basic.annotation.gff3.gz | DuckHTS | read_gff attributes_pairs          | 5,866,158 |      3 |    14.2517 |     411611.4 |        5.8 |             2.82 |
| gencode.v49.basic.annotation.gff3.gz | DuckHTS | read_gff attributes_map+list+pairs | 5,866,158 |      3 |    37.5871 |     156068.4 |        2.2 |             1.07 |
| gencode.v49.basic.annotation.gff3.gz | GFFBase | parse_gff strict (rust)            | 5,866,158 |      3 |    40.2213 |     145846.9 |        2.0 |             1.00 |

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
