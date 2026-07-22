BigWig reader benchmark
================

This report measures `read_bigwig(...)` on a dense, real conservation
track. The input is a local 10 Mb slice of the UCSC GRCh38 phyloP
100-way BigWig. Network transfer and slice construction are staging work
and are outside every timed interval.

There is no earlier checked-in DuckHTS BigWig benchmark. The indexed
multi-region BCF/VCF report in
[`benchmark_multi_region_readers.md`](benchmark_multi_region_readers.md)
is the nearest reader benchmark, but it exercises htslib variant
iterators and is not numerically comparable.

## Reproduction

Stage the real track with the UCSC utilities. Their interval flags are
zero-based and half-open, matching the stored BigWig coordinates.

``` bash
mkdir -p /tmp/ucsc-bigwig-bench
curl -fsSL https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph \
  -o /tmp/ucsc-bigwig-bench/bigWigToBedGraph
curl -fsSL https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig \
  -o /tmp/ucsc-bigwig-bench/bedGraphToBigWig
chmod +x /tmp/ucsc-bigwig-bench/bigWigToBedGraph \
  /tmp/ucsc-bigwig-bench/bedGraphToBigWig

/tmp/ucsc-bigwig-bench/bigWigToBedGraph \
  -chrom=chr22 -start=20000000 -end=30000000 \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw \
  /tmp/ucsc-bigwig-bench/chr22_20m_30m.bedGraph
printf 'chr22\t50818468\n' > /tmp/ucsc-bigwig-bench/hg38.chr22.sizes
/tmp/ucsc-bigwig-bench/bedGraphToBigWig \
  /tmp/ucsc-bigwig-bench/chr22_20m_30m.bedGraph \
  /tmp/ucsc-bigwig-bench/hg38.chr22.sizes \
  /tmp/ucsc-bigwig-bench/hg38_chr22_20m_30m_phyloP100way.bw
```

Build and install `Rduckhts` from the source revision under test, then
render:

``` bash
DUCKHTS_BIGWIG_BENCHMARK_INPUT=/tmp/ucsc-bigwig-bench/hg38_chr22_20m_30m_phyloP100way.bw \
  Rscript -e "rmarkdown::render('benchmarks/benchmark_bigwig_reader.Rmd', output_format='github_document')"
```

The benchmark pins the one-thread condition to CPU 12 and the
four-thread condition to physical CPUs 12–15 on the recorded host. Input
construction, extension loading, connection setup, warm-up, result
verification, and removal of the previous Parquet output are excluded
from timing.

| Property                   | Recorded value                      |
|:---------------------------|:------------------------------------|
| source revision            | 669d2c4d9b91                        |
| run date                   | 2026-07-22                          |
| input source               | UCSC hg38.phyloP100way.bw           |
| input interval             | chr22:\[20,000,000, 30,000,000)     |
| local BigWig bytes         | 73,950,347                          |
| full-scan stored intervals | 9,314,560                           |
| full-scan value sum        | 79960.8974669                       |
| full-scan fingerprint      | A2B271AC831ABD4A                    |
| indexed ranges             | 20                                  |
| indexed output intervals   | 4,656,131                           |
| DuckDB version             | v1.5.3                              |
| htslib version             | 1.24                                |
| CPU                        | 13th Gen Intel(R) Core(TM) i5-13500 |

## Results

The aggregate workloads decode the reader output but return one row to
R. The Parquet workload decodes and writes every interval as a real
ZSTD-compressed four-column relation. Every timed pass is followed by
exact row-count and order-independent full-row fingerprint validation;
floating-point sums are also checked to a relative tolerance of `1e-12`.

| Workload                      | Threads | CPU affinity | Runs | Output intervals | Parquet bytes | Minimum seconds | Median seconds | Maximum seconds | Median intervals/s |
|:------------------------------|--------:|-------------:|-----:|-----------------:|--------------:|----------------:|---------------:|----------------:|-------------------:|
| full scan + aggregate         |       1 |           12 |    9 |        9,314,560 |             – |           0.989 |          0.994 |           0.998 |          9,370,785 |
| full scan + aggregate         |       4 |        12-15 |    9 |        9,314,560 |             – |           0.987 |          0.993 |           1.005 |          9,380,222 |
| 20 indexed ranges + aggregate |       1 |           12 |    9 |        4,656,131 |             – |           0.497 |          0.498 |           0.503 |          9,349,661 |
| 20 indexed ranges + aggregate |       4 |        12-15 |    9 |        4,656,131 |             – |           0.130 |          0.130 |           0.131 |         35,816,392 |
| full scan + ZSTD Parquet      |       1 |           12 |    5 |        9,314,560 |    69,313,917 |           1.747 |          1.868 |           1.869 |          4,986,381 |
| full scan + ZSTD Parquet      |       4 |        12-15 |    5 |        9,314,560 |    69,313,917 |           1.716 |          1.716 |           1.728 |          5,428,065 |

The full scan contains one nonempty contig, so it deliberately tests
whether threading overhead remains bounded when contig-level parallel
work is absent. The separated indexed ranges expose independent tasks.
The measurement does not include remote latency, BigWig construction,
extension loading, or a downstream conservation reduction over variant
spans.
