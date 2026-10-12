Indexed BCF/VCF scan initialization: small-fixture cost
================

This paired measurement exercises automatic full-file scan planning,
complete and projected materialization, metadata counts and overlapping
indexed regions on committed BCF/VCF fixtures staged through
duckhtsbench’s registry. The registry pins encoded file bytes; staging
is a network-free verified copy, not fixture reconstruction. The
original two-record inputs have one occupied contig and 256 empty header
contigs. A second five-record corpus has two occupied contigs, two
samples, equal-position distinct alleles, and a physical duplicate; its
VCF header declares all, some, or none of the contigs. It uses
non-colocated indexes. Timings measure setup and small-result
materialization, not cohort throughput.

Both builds plan VCF work using the Tabix dictionary. The
one-occupied-contig VCF therefore selects a single full-file stream
instead of 257 indexed contig tasks; the five-record VCFs use two contig
tasks irrespective of header completeness. BCF retains its header
dictionary. The requested DuckDB thread limit is not a promise about the
number of active workers. Both builds retain parsed indexes at bind. The
candidate shares its file/header/iterator and decode-validation code
with a host-neutral scanner. Data and indexes remain unchanged during
timing.

The nearest larger full-column report is
[`benchmark_reader_allocations.md`](benchmark_reader_allocations.md):
4,048,342 GIAB HG002 records in wide/tidy form, one sequential scan
handle and zero HTSlib decompression workers. It does not exercise
automatic parallel index reloads and is not an identical-workload
baseline. The previous identical two/five-record full-materialization
workload is retained at
[`4929162:benchmark_bcf_scan_init.md`](https://github.com/RGenomicsETL/duckhts/blob/4929162aa103fe703c848a9b87bfc062cd53aa15/benchmarks/benchmark_bcf_scan_init.md).
Count-only, projected and explicit-region cases below are new paired
workloads; they do not have a previous identical recorded baseline. This
report measures the pre-fix binary alongside the candidate on the same
host.

## Reproduction and identity

Install duckhtsbench from this checkout. Set
`BCF_SCAN_BASELINE_EXTENSION` and `BCF_SCAN_BASELINE_REVISION` to a
separately built baseline. `DUCKHTS_EXTENSION` selects the candidate.
From the repository root:

``` sh
taskset -c 0,2,4,6 Rscript -e 'rmarkdown::render("benchmarks/benchmark_bcf_scan_init.Rmd")'
```

| property             | value                                                                                |
|:---------------------|:-------------------------------------------------------------------------------------|
| baseline revision    | a5ab5147                                                                             |
| candidate revision   | 3624158b5141082af6545e19a1be9e8f03f06a81                                             |
| candidate src tree   | 9ada6d903f188bc52198290590ce9b0599c765ca                                             |
| baseline binary MD5  | 3ee8f7c34478550c542334e13e927bfb                                                     |
| candidate binary MD5 | a8cbff3d6b339c1865237221e7f42502                                                     |
| R                    | R version 4.6.0 (2026-04-24)                                                         |
| DuckDB R package     | 1.5.5                                                                                |
| HTSlib               | 1.24                                                                                 |
| CPU affinity         | pid 320001’s current affinity list: 0,2,4,6                                          |
| repetitions          | 9 timed calls after warm-up; fresh R process per build/input/thread/workload setting |

| file                                      | bytes | md5                              |
|:------------------------------------------|------:|:---------------------------------|
| parallel_empty_contigs.bcf                |  1351 | e7576013a2a4cfc843130eff617ea527 |
| parallel_empty_contigs.vcf.gz             |   372 | 8583a37cb25fe73d7dc8925652804942 |
| bcf_scan_contigs.full.vcf.gz              |   424 | abd4cc644200e4812df5ccbf185e3a4f |
| bcf_scan_contigs.partial.vcf.gz           |   421 | 4dcaa2047c179520b1259d7671c8a0e9 |
| bcf_scan_contigs.none.vcf.gz              |   404 | 01c9b6bbc25b859aebf67f74d9e265c3 |
| parallel_empty_contigs.bcf.csi            |   117 | 1cb71253d787009c2867430c349f2769 |
| parallel_empty_contigs.vcf.gz.tbi         |   106 | 9f978fe03cf1a2d047683b8c2e578eae |
| bcf_scan_contigs.full.vcf.gz.index.tbi    |   126 | fbb48c2c7cd524e891c35b7bc11f6401 |
| bcf_scan_contigs.partial.vcf.gz.index.tbi |   126 | 953769f59bf8df49d7ea3ee4de9c4bb7 |
| bcf_scan_contigs.none.vcf.gz.index.tbi    |   126 | 698f40d722cfee51932a386cb0b2f937 |
| bcf_scan_contigs.bcf                      |   470 | 8e006ec9582009fc56a1616e5dc67010 |

| build     | workload  | format     | header      | duckdb_threads | htslib_workers_per_handle | header_contigs | input_records | selected_records | materialized_rows | materialized_columns | exact_rows | median_ms | min_ms | max_ms |
|:----------|:----------|:-----------|:------------|---------------:|--------------------------:|---------------:|--------------:|-----------------:|------------------:|---------------------:|:-----------|----------:|-------:|-------:|
| baseline  | full      | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     1.160 |  1.121 |  1.413 |
| candidate | full      | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     1.207 |  1.176 |  1.509 |
| baseline  | count     | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     1.141 |  1.041 |  1.283 |
| candidate | count     | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     1.143 |  1.005 |  1.297 |
| baseline  | projected | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     1.185 |  1.110 |  1.332 |
| candidate | projected | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     1.168 |  1.058 |  1.321 |
| baseline  | regions   | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     1.180 |  1.147 |  1.382 |
| candidate | regions   | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     1.165 |  1.110 |  1.405 |
| baseline  | full      | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     1.256 |  1.148 |  1.519 |
| candidate | full      | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     1.177 |  1.158 |  1.287 |
| baseline  | count     | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     1.184 |  1.079 |  1.290 |
| candidate | count     | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     1.170 |  1.061 |  1.300 |
| baseline  | projected | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     1.213 |  1.108 |  1.382 |
| candidate | projected | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     1.263 |  1.087 |  1.516 |
| baseline  | regions   | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     1.205 |  1.155 |  1.539 |
| candidate | regions   | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     1.225 |  1.117 |  1.408 |
| baseline  | full      | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     0.928 |  0.892 |  0.975 |
| candidate | full      | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     0.911 |  0.875 |  0.955 |
| baseline  | count     | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     0.832 |  0.761 |  0.890 |
| candidate | count     | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     0.817 |  0.762 |  0.842 |
| baseline  | projected | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     0.840 |  0.804 |  0.896 |
| candidate | projected | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     0.820 |  0.813 |  0.848 |
| baseline  | regions   | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     0.926 |  0.904 |  0.946 |
| candidate | regions   | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     0.932 |  0.900 |  0.981 |
| baseline  | full      | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     0.909 |  0.874 |  0.977 |
| candidate | full      | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     0.910 |  0.872 |  1.013 |
| baseline  | count     | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     0.791 |  0.772 |  0.860 |
| candidate | count     | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     0.844 |  0.816 |  1.013 |
| baseline  | projected | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     0.875 |  0.840 |  0.896 |
| candidate | projected | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     0.838 |  0.800 |  0.857 |
| baseline  | regions   | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     0.935 |  0.908 |  0.950 |
| candidate | regions   | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     0.886 |  0.839 |  0.906 |
| baseline  | full      | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 5 |                   10 | TRUE       |     0.951 |  0.921 |  0.997 |
| candidate | full      | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 5 |                   10 | TRUE       |     0.918 |  0.894 |  0.955 |
| baseline  | count     | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 1 |                    1 | TRUE       |     0.832 |  0.786 |  0.859 |
| candidate | count     | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 1 |                    1 | TRUE       |     0.798 |  0.773 |  0.839 |
| baseline  | projected | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 5 |                    4 | TRUE       |     0.843 |  0.824 |  0.880 |
| candidate | projected | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 5 |                    4 | TRUE       |     0.858 |  0.830 |  0.915 |
| baseline  | regions   | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                1 |                 1 |                   10 | TRUE       |     0.929 |  0.919 |  0.982 |
| candidate | regions   | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                1 |                 1 |                   10 | TRUE       |     0.881 |  0.854 |  0.971 |
| baseline  | full      | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 5 |                   10 | TRUE       |     0.953 |  0.911 |  1.006 |
| candidate | full      | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 5 |                   10 | TRUE       |     0.965 |  0.924 |  1.100 |
| baseline  | count     | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 1 |                    1 | TRUE       |     0.789 |  0.763 |  0.849 |
| candidate | count     | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 1 |                    1 | TRUE       |     0.762 |  0.738 |  0.874 |
| baseline  | projected | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 5 |                    4 | TRUE       |     0.830 |  0.803 |  0.869 |
| candidate | projected | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 5 |                    4 | TRUE       |     0.842 |  0.810 |  0.913 |
| baseline  | regions   | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                1 |                 1 |                   10 | TRUE       |     0.946 |  0.924 |  1.015 |
| candidate | regions   | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                1 |                 1 |                   10 | TRUE       |     0.965 |  0.937 |  1.006 |
| baseline  | full      | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 5 |                   10 | TRUE       |     0.950 |  0.891 |  1.044 |
| candidate | full      | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 5 |                   10 | TRUE       |     0.924 |  0.888 |  0.981 |
| baseline  | count     | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 1 |                    1 | TRUE       |     0.855 |  0.846 |  0.877 |
| candidate | count     | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 1 |                    1 | TRUE       |     0.820 |  0.790 |  0.842 |
| baseline  | projected | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 5 |                    4 | TRUE       |     0.834 |  0.804 |  0.883 |
| candidate | projected | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 5 |                    4 | TRUE       |     0.879 |  0.834 |  0.897 |
| baseline  | regions   | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                1 |                 1 |                   10 | TRUE       |     0.923 |  0.876 |  1.004 |
| candidate | regions   | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                1 |                 1 |                   10 | TRUE       |     0.881 |  0.853 |  0.994 |
| baseline  | full      | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 5 |                   10 | TRUE       |     0.910 |  0.885 |  1.019 |
| candidate | full      | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 5 |                   10 | TRUE       |     0.926 |  0.902 |  1.440 |
| baseline  | count     | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 1 |                    1 | TRUE       |     0.823 |  0.740 |  0.839 |
| candidate | count     | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 1 |                    1 | TRUE       |     0.820 |  0.785 |  0.893 |
| baseline  | projected | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 5 |                    4 | TRUE       |     0.853 |  0.836 |  0.912 |
| candidate | projected | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 5 |                    4 | TRUE       |     0.857 |  0.826 |  0.931 |
| baseline  | regions   | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                1 |                 1 |                   10 | TRUE       |     0.942 |  0.901 |  0.991 |
| candidate | regions   | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                1 |                 1 |                   10 | TRUE       |     0.905 |  0.849 |  0.958 |
| baseline  | full      | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 5 |                   10 | TRUE       |     0.954 |  0.927 |  1.017 |
| candidate | full      | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 5 |                   10 | TRUE       |     0.957 |  0.923 |  0.980 |
| baseline  | count     | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 1 |                    1 | TRUE       |     0.881 |  0.854 |  0.928 |
| candidate | count     | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 1 |                    1 | TRUE       |     0.873 |  0.860 |  0.880 |
| baseline  | projected | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 5 |                    4 | TRUE       |     0.860 |  0.841 |  0.883 |
| candidate | projected | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 5 |                    4 | TRUE       |     0.872 |  0.865 |  0.924 |
| baseline  | regions   | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                1 |                 1 |                   10 | TRUE       |     0.908 |  0.885 |  1.012 |
| candidate | regions   | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                1 |                 1 |                   10 | TRUE       |     0.933 |  0.888 |  1.003 |
| baseline  | full      | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 5 |                   10 | TRUE       |     1.009 |  0.977 |  1.098 |
| candidate | full      | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 5 |                   10 | TRUE       |     0.992 |  0.970 |  1.177 |
| baseline  | count     | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 1 |                    1 | TRUE       |     0.829 |  0.786 |  0.880 |
| candidate | count     | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 1 |                    1 | TRUE       |     0.828 |  0.788 |  0.884 |
| baseline  | projected | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 5 |                    4 | TRUE       |     0.860 |  0.830 |  0.900 |
| candidate | projected | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 5 |                    4 | TRUE       |     0.879 |  0.839 |  0.931 |
| baseline  | regions   | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                1 |                 1 |                   10 | TRUE       |     0.921 |  0.907 |  0.965 |
| candidate | regions   | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                1 |                 1 |                   10 | TRUE       |     0.954 |  0.922 |  1.016 |

Each timed `CREATE OR REPLACE TEMP TABLE AS SELECT` materializes the
requested result. Full/projected cases return two or five rows, regions
return one selected record, and count returns one aggregate row
representing two or five records. Count may read only index statistics;
it is not a physical decode throughput test. Elapsed `Sys.time()`
includes the DBI call; connection creation, warm-up, and complete
typed-row verification are outside the timer. Every candidate result and
every valid baseline agrees exactly with the complete BCF sequential
rows, including GT, ALT and physical duplicate multiplicity, without
asserting scan arrival order. The BCF oracle also agrees across builds.

There is no public sample-selection option on canonical `read_bcf`;
retired experimental selectors are not reinstated to produce a timing.
Native selected sample tests are correctness evidence, not a
selected-cohort benchmark.

This does not measure the latency of an index failure, remote I/O,
large-cohort throughput, RSS, or maximum worker concurrency. Thread
settings are DuckDB’s requested limits. Small timings can vary with
scheduling; no general speedup or no-regression claim follows from them.
