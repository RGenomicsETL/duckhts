BCF appender contig-parallel benchmark
================

This report measures `read_bcf_appender(...)` after making
`region_threads` partition work by primary contig. Each contig keeps its
complete interval set in one htslib 1.24 multi-region iterator, so
repeated and overlapping intervals produce the same target schema and
row multiset at every thread count.

The previous parallel implementation scanned each requested interval
independently, emitted duplicate rows, and added `duckhts_region_idx`.
Its row denominator and write schema therefore differ from this
workload. There is no valid before/after throughput comparison. The
nearest checked-in measurement is
[`benchmark_multi_region_readers.md`](benchmark_multi_region_readers.md),
which measures native iterator decoding without appender
materialization.

## Reproduction

Build and install the current `Rduckhts` tarball, then render with the
process allowed on all CPUs named below. The script pins the serial run
to one CPU and the four-job run to four CPUs.

``` bash
DUCKHTS_BENCH_SERIAL_CPU=16 \
DUCKHTS_BENCH_PARALLEL_CPUS=16-19 \
taskset -c 16-19 Rscript -e "rmarkdown::render('benchmarks/benchmark_bcf_appender_contigs.Rmd', output_format = 'github_document')"
```

The fixture contains four contigs with consecutive SNVs. Each contig
receives eight 50,000-base intervals at a 25,000-base step. Fixture
creation, BCF conversion, indexing, connection setup, and warm-up are
outside the timed interval. The timed `seconds` value covers
target-table replacement, contig-job planning, indexed reads, narrow-row
materialization through DuckDB appenders, and commit.

| Property                             | Recorded value                                    |
|:-------------------------------------|:--------------------------------------------------|
| source revision                      | 408fbcb5bc7c                                      |
| run date                             | 2026-07-18                                        |
| DuckDB version                       | v1.5.3                                            |
| htslib version                       | 1.24                                              |
| CPU                                  | 13th Gen Intel(R) Core(TM) i5-13500               |
| serial CPU                           | 16                                                |
| parallel CPUs                        | 16-19                                             |
| DuckDB threads                       | 1                                                 |
| htslib decompression threads per job | 0                                                 |
| input BCF records                    | 1,000,000                                         |
| primary contigs                      | 4                                                 |
| requested intervals                  | 32                                                |
| interval width                       | 50,000                                            |
| interval step                        | 25,000                                            |
| raw interval-record hits             | 1,600,000                                         |
| unique output records                | 900,000                                           |
| overlap work factor                  | 1.778                                             |
| uncompressed VCF bytes               | 23,555,760                                        |
| BCF bytes                            | 2,314,452                                         |
| bcftools version                     | bcftools 1.23.1-70-g6dbd8fef; Using htslib 1.22.1 |

## Results

The denominator is the 900,000 unique target rows after htslib merges
the overlapping intervals. Both runs materialize `FILE_OFFSET`, validate
the exact row count, position checksum, and distinct record offsets, and
finish with a bidirectional `EXCEPT ALL` comparison of zero rows.

| Scheduling                         | region_threads | Recorded affinity | Runs | Unique rows | Median seconds | Median unique rows/s |
|:-----------------------------------|---------------:|:------------------|-----:|------------:|---------------:|---------------------:|
| one complete multi-region iterator |              1 | 16                |   15 |     900,000 |          0.300 |            2,999,997 |
| one iterator per primary contig    |              4 | 16-19             |   15 |     900,000 |          0.176 |            5,111,454 |

This is a warm-page-cache local-BCF materialization measurement. It does
not measure cold storage, remote object-store latency, rich INFO/FORMAT
decoding, or full `read_bcf(...)` output. The parallel result measures
four independent contig iterators and four DuckDB appenders; a request
confined to one contig correctly remains a one-iterator workload.
