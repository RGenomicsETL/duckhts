GenBank reader benchmark
================

This report measures `read_genbank(...)` and `genbank_to_fasta(...)` on
a real annotation: the NCBI RefSeq record for *Escherichia coli* K-12
MG1655 (`NC_000913.3`), a 4,641,652 bp circular genome carrying 9,306
features across 8 feature kinds. Fetching the record is staging work and
is outside every timed interval.

There is no earlier checked-in DuckHTS GenBank benchmark. The nearest
reader report is
[`benchmark_gffbase_conformance.md`](benchmark_gffbase_conformance.md),
which measures `read_gff(...)`. It is the closest analogue because
`read_genbank` emits `read_gff`’s exact column shape, but the two are
not numerically comparable: that workload reads tabix-style tabular
rows, while this parses a GenBank flat file line by line through hFILE.

`read_genbank` declares no parallel scan and parses the whole record in
`init` before the first chunk, so the one- and four-thread conditions
exist to show whether threading adds overhead rather than to show it
scales. Memory therefore tracks the annotation size rather than the
projection; that limit is visible in the input denominators below but is
not itself measured here.

## Reproduction

Stage the record from NCBI. `gbwithparts` returns the full FEATURES
table and the ORIGIN sequence, which is what both functions read.

``` sh
mkdir -p /tmp/ncbi-genbank-bench
curl -fsSL -o /tmp/ncbi-genbank-bench/NC_000913.3.gb \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=gbwithparts&retmode=text"
```

Then render this report from a tracked-clean revision:

``` sh
Rscript -e 'rmarkdown::render("benchmarks/benchmark_genbank_reader.Rmd")'
```

The benchmark pins the one-thread condition to logical CPU 12 and the
four-thread condition to logical CPUs 12–15 on the recorded host. Input
staging, extension loading, connection setup, warm-up, result
verification, and removal of the previous FASTA output are excluded from
timing.

    #> duckdb is storing downloaded extensions and secrets under ~/.duckdb:
    #> ℹ /home/ryan.ward/.duckdb
    #> This persists across sessions and is shared with the DuckDB CLI and other clients.
    #> ℹ Run duckdb(shared_home = FALSE) to use a temporary directory instead.
    #> ℹ See ?duckdb_storage for details and alternatives.

| Property               | Recorded value                       |
|:-----------------------|:-------------------------------------|
| source revision        | f1aad54e4655                         |
| run date               | 2026-09-18                           |
| input source           | NCBI RefSeq efetch (gbwithparts)     |
| input accession        | NC_000913.3                          |
| GenBank bytes          | 11,450,917                           |
| genome bases           | 4,641,652                            |
| emitted features       | 9,306                                |
| distinct feature kinds | 8                                    |
| feature fingerprint    | 231BDB64ADB2A55                      |
| DuckDB version         | v1.5.5                               |
| htslib version         | 1.24                                 |
| CPU                    | 13th Gen Intel(R) Core(TM) i9-13900K |

## Results

The aggregate workloads parse the FEATURES table and return one row to
R, so the timing covers parsing and decoding rather than transport. The
FASTA workload reads ORIGIN instead and writes every base as a real
file. Every timed pass is verified: the aggregates by exact row count
and an order-independent full-row fingerprint, with the summed span
checked to a relative tolerance of `1e-12`, and the FASTA by reading it
back and asserting both the record count and 4,641,652 bases.

| Workload | Threads | CPU affinity | Runs | Output rows | FASTA bytes | Minimum seconds | Median seconds | Maximum seconds | Median rows/s |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| features + aggregate | 1 | 12 | 9 | 9,306 | – | 0.069 | 0.072 | 0.073 | 129,250 |
| features + aggregate | 4 | 12-15 | 9 | 9,306 | – | 0.071 | 0.072 | 0.078 | 129,250 |
| features + attribute MAP | 1 | 12 | 9 | 9,306 | – | 0.074 | 0.076 | 0.078 | 122,447 |
| features + attribute MAP | 4 | 12-15 | 9 | 9,306 | – | 0.075 | 0.076 | 0.079 | 122,447 |
| ORIGIN to FASTA | 1 | 12 | 5 | 1 | 4,708,035 | 0.013 | 0.014 | 0.015 | 331,546,571 |
| ORIGIN to FASTA | 4 | 12-15 | 5 | 1 | 4,708,035 | 0.013 | 0.013 | 0.014 | 357,050,154 |

The features and ORIGIN rows count different things: the aggregate
workloads emit one row per feature, while the FASTA rate is expressed in
bases written, since the converter produces a single record for this
genome. Adding `attributes_map := TRUE` builds a `MAP<VARCHAR,VARCHAR>`
per feature on top of the raw attribute string and is the more expensive
projection, which is why both are recorded.

The measurement does not include the NCBI fetch, extension loading, or
any downstream join against the emitted intervals. It also does not
bound memory: the reader materializes every feature before the first
chunk, so a much larger annotation is a different question than this
report answers.
