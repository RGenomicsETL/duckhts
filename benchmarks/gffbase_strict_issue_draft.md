Draft upstream GFFBase strict-mode issue
================

<!-- gffbase_strict_issue_draft.md is generated from gffbase_strict_issue_draft.Rmd. -->

# Copy-ready issue body

## Summary

`parse_gff(..., strict=True)` accepts several GFF3 rows that look
structurally malformed: rows with a 10th tab-separated field, rows with
a partially malformed attribute segment such as `ID=ok;broken`, and rows
with an unknown start but concrete `end=0`.

## Environment

|                             |                                                                                                                                                                                         |
|-----------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **OS**                      | Ubuntu 24.04.3 LTS                                                                                                                                                                      |
| **Python version**          | 3.13.12                                                                                                                                                                                 |
| **GFFBase version**         | 0.1.0                                                                                                                                                                                   |
| **Install method**          | `pip install --target .tmp/gffbase_site gffbase==0.1.0` for this reproduction; also cross-checked against the local source mirror at commit `78714cf30a9d799eab544e00a79a4da9754987ca`. |
| **DuckDB version**          | 1.5.2                                                                                                                                                                                   |
| **PyArrow version**         | 24.0.0                                                                                                                                                                                  |
| **Native extension built?** | True                                                                                                                                                                                    |

## Minimal reproducible example

``` python
import tempfile
from pathlib import Path
import gffbase

cases = {
    "extra_field": b"chr1\tsrc\texon\t1\t10\t.\t+\t.\tID=ok\textra\n",
    "partial_malformed_attribute": b"chr1\tsrc\texon\t1\t10\t.\t+\t.\tID=ok;broken\n",
    "nonpositive_end_with_unknown_start": b"chr1\tsrc\tregion\t.\t0\t.\t.\t.\tID=x\n",
}

for name, text in cases.items():
    path = Path(tempfile.mkdtemp()) / f"{name}.gff3"
    path.write_bytes(text)
    try:
        rows = list(gffbase.parse_gff(str(path), strict=True, engine="rust"))
        first = rows[0] if rows else None
        print(name, "ACCEPTED", "rows=", len(rows))
        if first is not None:
            print("  start/end:", first.start, first.end)
            print("  attributes_pairs:", first.attributes_pairs)
            print("  extra:", first.extra)
    except Exception as exc:
        print(name, "REJECTED", type(exc).__name__, str(exc))
```

Minimal GFF3 snippets:

``` text
##gff-version 3
chr1    src    exon      1    10    .    +    .    ID=ok    extra
chr1    src    exon      1    10    .    +    .    ID=ok;broken
chr1    src    region    .    0     .    .    .    ID=x
```

## What you expected

I expected `strict=True` either to reject these rows or for the
documentation to explicitly describe this as intentionally permissive /
gffutils-compatible parsing.

My reading is:

- GFF3 feature rows are structurally nine tab-separated columns; a 10th
  field is not part of the core GFF3 row shape.
- Non-empty GFF3 attribute segments are `tag=value` pairs separated by
  semicolons; `ID=ok;broken` has one valid pair plus one malformed
  segment.
- Concrete coordinates are 1-based positive integers. If a parser allows
  unknown `start='.'`, a concrete `end=0` still seems nonpositive.
- The current `parse_gff()` docstring says `strict=True` raises
  `GFFFormatError` on the first malformed line.

## What actually happened

On `gffbase==0.1.0` with the native Rust parser available, all three
rows were accepted:

``` text
extra_field ACCEPTED rows= 1
  start/end: 1 10
  attributes_pairs: [('ID', 'ok', 0)]
  extra: ['extra']
partial_malformed_attribute ACCEPTED rows= 1
  start/end: 1 10
  attributes_pairs: [('ID', 'ok', 0), ('broken', '', 0)]
  extra: []
nonpositive_end_with_unknown_start ACCEPTED rows= 1
  start/end: None 0
  attributes_pairs: [('ID', 'x', 0)]
  extra: []
```

## Additional context

We found this while implementing strict GFF3 validation in DuckHTS.
GFFBase’s parser and NCBI/GFF3-oriented tests were useful inspiration
for DuckHTS’ own strict mode, so this report is meant as
cross-validation feedback rather than a criticism. DuckHTS now
intentionally rejects the three cases above under
`read_gff(..., strict := true)`.

For context, we also ran direct parser/scan benchmarks against GFFBase.
These numbers compare DuckHTS direct DuckDB scans with
`gffbase.parse_gff(strict=True, engine="rust")`; they do **not** compare
GFFBase’s heavier `create_db()` annotation-database ingest path.

| dataset                              | tool    | variant                            | rows      | passes | median_sec | vs_gffbase_parse |
|:-------------------------------------|:--------|:-----------------------------------|:----------|-------:|-----------:|-----------------:|
| synthetic_200000.gff3                | DuckHTS | read_gff strict COUNT(\*)          | 200,000   |      3 |     0.0672 |             8.25 |
| synthetic_200000.gff3                | DuckHTS | read_gff attributes_list           | 200,000   |      3 |     0.1790 |             3.09 |
| synthetic_200000.gff3                | DuckHTS | read_gff attributes_pairs          | 200,000   |      3 |     0.1246 |             4.45 |
| synthetic_200000.gff3                | DuckHTS | read_gff attributes_map+list+pairs | 200,000   |      3 |     0.3943 |             1.41 |
| synthetic_200000.gff3                | GFFBase | parse_gff strict (rust)            | 200,000   |      3 |     0.5540 |             1.00 |
| gencode.v49.basic.annotation.gff3.gz | DuckHTS | read_gff strict COUNT(\*)          | 5,866,158 |      3 |     5.1975 |             7.86 |
| gencode.v49.basic.annotation.gff3.gz | DuckHTS | read_gff attributes_list           | 5,866,158 |      3 |    17.2256 |             2.37 |
| gencode.v49.basic.annotation.gff3.gz | DuckHTS | read_gff attributes_pairs          | 5,866,158 |      3 |    13.7290 |             2.97 |
| gencode.v49.basic.annotation.gff3.gz | DuckHTS | read_gff attributes_map+list+pairs | 5,866,158 |      3 |    38.2677 |             1.07 |
| gencode.v49.basic.annotation.gff3.gz | GFFBase | parse_gff strict (rust)            | 5,866,158 |      3 |    40.8327 |             1.00 |

Benchmark machine:

| parameter               | value                               |
|:------------------------|:------------------------------------|
| OS                      | Ubuntu 24.04.3 LTS                  |
| kernel / machine        | 6.8.0-78-generic / x86_64           |
| CPU                     | 13th Gen Intel(R) Core(TM) i5-13500 |
| logical / affinity CPUs | 20 / 20                             |
| memory                  | 62.58 GiB                           |
| DuckDB threads          | 4                                   |

One performance observation: GFFBase’s Arrow staging path is a good
Python-to-DuckDB ingestion strategy, especially for `create_db()`. For
maximum DuckDB throughput, however, a future Rust-native path could
avoid the
`Rust parser -> Python objects/lists -> PyArrow Table -> DuckDB INSERT`
boundary. If `duckdb-rs` exposes the needed
appender/table-function/vector APIs for the supported DuckDB versions,
writing parsed records directly into DuckDB vectors from Rust may be a
simpler high-performance path. If not, Arrow remains the pragmatic
Python boundary.
