Indexed multi-region reader benchmark
================

This report measures the htslib 1.24 native multi-region iterator used
by `read_bcf(...)`. The same generated records are read from indexed
BCF, which uses `bcf_itr_regarray(...)`, and indexed VCF text, which
uses `tbx_itr_regarray(...)`. The requested intervals overlap by design;
each query must emit every matching VCF record once.

There is no earlier checked-in multi-region reader benchmark with the
same workload. The full-file `read_bcf(...)` measurements in
[`benchmark_variantkey_conformance.md`](benchmark_variantkey_conformance.md)
are the nearest reader baseline, but they are not directly comparable
because they do not create indexed region iterators.

## Reproduction

Build and install the current `Rduckhts` tarball, then render on one
allowed CPU. For example:

``` bash
taskset -c 2 Rscript -e "rmarkdown::render('benchmarks/benchmark_multi_region_readers.Rmd', output_format = 'github_document')"
```

The fixture is generated locally and deterministically. Fixture
creation, BGZF compression, indexing, connection setup, and warm-up are
outside the timed interval.

| Property                     | Recorded value                                                                                                                                                                                                                                         |
|:-----------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| source revision              | 2abd25645223                                                                                                                                                                                                                                           |
| run date                     | 2026-07-18                                                                                                                                                                                                                                             |
| DuckDB version               | v1.5.3                                                                                                                                                                                                                                                 |
| htslib version               | 1.24                                                                                                                                                                                                                                                   |
| CPU                          | 13th Gen Intel(R) Core(TM) i5-13500                                                                                                                                                                                                                    |
| process affinity             | pid 2811511’s current affinity list: 2                                                                                                                                                                                                                 |
| DuckDB threads               | 1                                                                                                                                                                                                                                                      |
| htslib decompression threads | 0                                                                                                                                                                                                                                                      |
| input records per format     | 2,000,000                                                                                                                                                                                                                                              |
| requested intervals          | 64                                                                                                                                                                                                                                                     |
| interval width               | 50,000                                                                                                                                                                                                                                                 |
| interval step                | 25,000                                                                                                                                                                                                                                                 |
| raw interval-record hits     | 3,200,000                                                                                                                                                                                                                                              |
| unique output records        | 1,625,000                                                                                                                                                                                                                                              |
| overlap work factor          | 1.969                                                                                                                                                                                                                                                  |
| uncompressed VCF bytes       | 48,888,987                                                                                                                                                                                                                                             |
| VCF.gz bytes                 | 3,391,626                                                                                                                                                                                                                                              |
| BCF bytes                    | 4,626,191                                                                                                                                                                                                                                              |
| bgzip version                | bgzip (htslib) 1.19 Copyright (C) 2023 Genome Research Ltd.                                                                                                                                                                                            |
| tabix version                | tabix (htslib) 1.19 Copyright (C) 2023 Genome Research Ltd.                                                                                                                                                                                            |
| bcftools version             | bcftools 1.23.1-70-g6dbd8fef Using htslib 1.22.1 Copyright (C) 2025 Genome Research Ltd. License Expat: The MIT/Expat license This is free software: you are free to change and redistribute it. There is NO WARRANTY, to the extent permitted by law. |

## Results

The denominator is unique VCF records returned after htslib merges the
64 overlapping intervals. `raw interval-record hits` is what independent
per-region scans would encounter before duplicate removal. Both formats
are validated against the exact expected row count and position checksum
on every timed run.

| Format and native iterator      | Runs | Unique rows | Minimum seconds | Median seconds | Maximum seconds | Median unique rows/s |
|:--------------------------------|-----:|------------:|----------------:|---------------:|----------------:|---------------------:|
| BCF + CSI (bcf_itr_regarray)    |    9 |   1,625,000 |           0.084 |          0.086 |           0.086 |           18,895,349 |
| VCF.gz + TBI (tbx_itr_regarray) |    9 |   1,625,000 |           0.282 |          0.284 |           0.289 |            5,721,831 |

This is a warm-page-cache, one-core iterator/decode measurement. It does
not measure cold storage, remote object-store latency, rich INFO/FORMAT
decoding, or DuckDB materialization of every returned column. It
specifically protects the indexed overlapping-region path changed by the
htslib 1.24 integration.
