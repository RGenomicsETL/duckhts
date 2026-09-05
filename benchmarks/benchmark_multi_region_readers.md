Indexed multi-region reader benchmark
================

This compares the canonical `read_bcf` indexed region path before and
after shared region-list validation and retirement of the appender
experiment. HTSlib 1.24 still owns interval union and record visitation:
BCF uses `bcf_itr_regarray`, and indexed VCF text uses
`tbx_itr_regarray`.

The same-run baseline and candidate read identical cached files, with
equal schema, record and aggregate-output denominators. The nearest
identical logical workload was recorded at `2abd25645223` on July 18,
2026: 0.086 s BCF and 0.284 s VCF.gz medians, one DuckDB thread and no
HTSlib workers. That historical run used different compression tools; it
is context, not an environment-matched regression comparison. The old
rendered report remains in git history.

## Reproduction and identity

Build the two named revisions separately. Set
`DUCKHTS_REGION_BASELINE_EXTENSION` and
`DUCKHTS_REGION_BASELINE_REVISION` to the pre-fix binary and commit.
`DUCKHTS_EXTENSION` selects the candidate (default: release build).
Render on one allowed CPU, for example:

``` bash
taskset -c 2 Rscript -e "rmarkdown::render('benchmarks/benchmark_multi_region_readers.Rmd')"
```

The network-free `duckhts_bench_stage_multi_region()` stage resolves the
synthetic VCF and its derived BCF/BGZF/index artifacts through the
benchmark registry. It regenerates the exact consecutive-SNV workload
before timing; creation, compression, indexing, connection setup and
warm-up are excluded. Each build uses a separate in-memory database
connection, with baseline first.

| property                  | value                                                                                                                                                                                                                                                  |
|:--------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| baseline revision         | 9b4b232127d5043716a1ec34cb66b96f2c7b2c03                                                                                                                                                                                                               |
| candidate revision        | d6278e154aef6a877039b2d5fb9ba80acc512bb2                                                                                                                                                                                                               |
| candidate src tree        | b78accf31f330d42712944a3244fdb3422d21b31                                                                                                                                                                                                               |
| baseline binary MD5       | 725b64e3c894a1ef045245d53c4a98a6                                                                                                                                                                                                                       |
| candidate binary MD5      | 5704c4f173be9553e6e0d94fe1bc720d                                                                                                                                                                                                                       |
| run date                  | 2026-09-05                                                                                                                                                                                                                                             |
| DuckDB R package          | 1.5.3                                                                                                                                                                                                                                                  |
| HTSlib                    | 1.24                                                                                                                                                                                                                                                   |
| CPU                       | 13th Gen Intel(R) Core(TM) i5-13500                                                                                                                                                                                                                    |
| CPU affinity              | pid 468465’s current affinity list: 2                                                                                                                                                                                                                  |
| DuckDB threads            | 1                                                                                                                                                                                                                                                      |
| HTSlib workers per handle | 0                                                                                                                                                                                                                                                      |
| input records per format  | 2000000                                                                                                                                                                                                                                                |
| requested intervals       | 64                                                                                                                                                                                                                                                     |
| interval width            | 50000                                                                                                                                                                                                                                                  |
| interval step             | 25000                                                                                                                                                                                                                                                  |
| raw interval-record hits  | 3200000                                                                                                                                                                                                                                                |
| unique reader rows        | 1625000                                                                                                                                                                                                                                                |
| SQL aggregate rows        | 1                                                                                                                                                                                                                                                      |
| bgzip version             | bgzip (htslib) 1.19 Copyright (C) 2023 Genome Research Ltd.                                                                                                                                                                                            |
| tabix version             | tabix (htslib) 1.19 Copyright (C) 2023 Genome Research Ltd.                                                                                                                                                                                            |
| bcftools version          | bcftools 1.23.1-70-g6dbd8fef Using htslib 1.22.1 Copyright (C) 2025 Genome Research Ltd. License Expat: The MIT/Expat license This is free software: you are free to change and redistribute it. There is NO WARRANTY, to the extent permitted by law. |

| artifact |    bytes | md5                              |
|:---------|---------:|:---------------------------------|
| vcf      | 48888987 | 4ba8f8a4b430d5a4a583b5b26fdad307 |
| vcf_gz   |  3391626 | 84e23cd87a0ab5a05fbe730ad38a2c7f |
| tbi      |     1371 | 1e62024e109cbe839370f85379b419b8 |
| bcf      |  4626191 | 19e549c948bd7ea285d4610eda32d3f4 |
| csi      |     1222 | 4fa8c74bc1a3a62178bd3a286847afe9 |

## Results

The reader-row denominator is the 1,625,000 stored records returned
after HTSlib unions 64 overlapping intervals, from 2,000,000 input
records per format. Independent interval scans would encounter 3,200,000
record hits. Each SQL query returns one aggregate row. All nine timed
runs validate both count and position checksum; separate SQL/R tests
cover long REF/END spans, missing targets and physical duplicate
records.

| build     | format       | runs | reader_rows | sql_output_rows | min_seconds | median_seconds | max_seconds |
|:----------|:-------------|-----:|------------:|----------------:|------------:|---------------:|------------:|
| baseline  | BCF + CSI    |    9 |     1625000 |               1 |       0.078 |          0.080 |       0.082 |
| baseline  | VCF.gz + TBI |    9 |     1625000 |               1 |       0.303 |          0.307 |       0.310 |
| candidate | BCF + CSI    |    9 |     1625000 |               1 |       0.078 |          0.078 |       0.082 |
| candidate | VCF.gz + TBI |    9 |     1625000 |               1 |       0.301 |          0.304 |       0.311 |

| format       | median_seconds_baseline | median_seconds_candidate | percent_change |
|:-------------|------------------------:|-------------------------:|---------------:|
| BCF + CSI    |                   0.080 |                    0.078 |         -2.500 |
| VCF.gz + TBI |                   0.307 |                    0.304 |         -0.977 |

This is a warm-page-cache, one-core iterator/decode measurement,
including region parsing and initialization. It does not measure remote
seeks, cold storage, rich INFO/FORMAT decoding, every-column
materialization, BAM/tabix/FASTA reader performance, or staged
pipelining. The removed appender cannot be benchmarked at this revision;
its historical narrow-row timings in
[benchmark_bcf_appender_contigs.md](benchmark_bcf_appender_contigs.md)
are not canonical-reader throughput measurements.
