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
| candidate revision   | 68dfa0332ad9ea34eea556eddea231cd83d842f6                                             |
| candidate src tree   | 2de704e472ceeda276caa0eeef7093221378f3c4                                             |
| baseline binary MD5  | 3ee8f7c34478550c542334e13e927bfb                                                     |
| candidate binary MD5 | fe0858278cd4db241ba3dac8f6323d6e                                                     |
| R                    | R version 4.6.0 (2026-04-24)                                                         |
| DuckDB R package     | 1.5.5                                                                                |
| HTSlib               | 1.24                                                                                 |
| CPU affinity         | pid 548625’s current affinity list: 0,2,4,6                                          |
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
| baseline  | full      | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     1.170 |  1.136 |  1.500 |
| candidate | full      | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     1.189 |  1.157 |  1.429 |
| baseline  | count     | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     1.203 |  1.066 |  1.276 |
| candidate | count     | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     1.139 |  1.050 |  1.256 |
| baseline  | projected | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     1.193 |  1.113 |  1.314 |
| candidate | projected | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     1.234 |  1.084 |  1.871 |
| baseline  | regions   | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     1.196 |  1.151 |  1.524 |
| candidate | regions   | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     1.204 |  1.131 |  1.304 |
| baseline  | full      | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     1.206 |  1.139 |  1.618 |
| candidate | full      | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     1.220 |  1.172 |  1.463 |
| baseline  | count     | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     1.205 |  1.066 |  1.415 |
| candidate | count     | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     1.226 |  1.047 |  1.277 |
| baseline  | projected | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     1.207 |  1.110 |  1.328 |
| candidate | projected | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     1.154 |  1.078 |  1.453 |
| baseline  | regions   | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     1.280 |  1.185 |  1.422 |
| candidate | regions   | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     1.177 |  1.152 |  1.511 |
| baseline  | full      | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     0.900 |  0.880 |  0.949 |
| candidate | full      | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     0.908 |  0.895 |  0.928 |
| baseline  | count     | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     0.834 |  0.768 |  0.898 |
| candidate | count     | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     0.823 |  0.759 |  0.837 |
| baseline  | projected | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     0.849 |  0.805 |  1.048 |
| candidate | projected | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     0.842 |  0.796 |  0.885 |
| baseline  | regions   | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     0.914 |  0.906 |  1.001 |
| candidate | regions   | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     0.932 |  0.886 |  1.018 |
| baseline  | full      | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     0.901 |  0.884 |  0.948 |
| candidate | full      | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     0.921 |  0.882 |  0.955 |
| baseline  | count     | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     0.816 |  0.746 |  0.855 |
| candidate | count     | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     0.809 |  0.762 |  0.873 |
| baseline  | projected | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     0.853 |  0.811 |  0.890 |
| candidate | projected | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     0.836 |  0.813 |  0.874 |
| baseline  | regions   | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     0.922 |  0.896 |  0.962 |
| candidate | regions   | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     0.916 |  0.875 |  1.034 |
| baseline  | full      | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 5 |                   10 | TRUE       |     0.943 |  0.910 |  1.000 |
| candidate | full      | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 5 |                   10 | TRUE       |     0.934 |  0.901 |  0.975 |
| baseline  | count     | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 1 |                    1 | TRUE       |     0.827 |  0.810 |  1.369 |
| candidate | count     | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 1 |                    1 | TRUE       |     0.837 |  0.824 |  0.863 |
| baseline  | projected | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 5 |                    4 | TRUE       |     0.834 |  0.799 |  1.036 |
| candidate | projected | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 5 |                    4 | TRUE       |     0.843 |  0.821 |  0.882 |
| baseline  | regions   | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                1 |                 1 |                   10 | TRUE       |     0.943 |  0.916 |  1.019 |
| candidate | regions   | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                1 |                 1 |                   10 | TRUE       |     0.901 |  0.879 |  0.959 |
| baseline  | full      | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 5 |                   10 | TRUE       |     0.971 |  0.916 |  1.026 |
| candidate | full      | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 5 |                   10 | TRUE       |     1.033 |  0.894 |  1.376 |
| baseline  | count     | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 1 |                    1 | TRUE       |     0.801 |  0.740 |  0.851 |
| candidate | count     | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 1 |                    1 | TRUE       |     0.798 |  0.772 |  0.890 |
| baseline  | projected | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 5 |                    4 | TRUE       |     0.887 |  0.808 |  0.940 |
| candidate | projected | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 5 |                    4 | TRUE       |     0.839 |  0.805 |  0.911 |
| baseline  | regions   | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                1 |                 1 |                   10 | TRUE       |     0.915 |  0.886 |  0.994 |
| candidate | regions   | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                1 |                 1 |                   10 | TRUE       |     0.895 |  0.848 |  0.965 |
| baseline  | full      | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 5 |                   10 | TRUE       |     0.981 |  0.964 |  1.239 |
| candidate | full      | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 5 |                   10 | TRUE       |     0.962 |  0.930 |  1.030 |
| baseline  | count     | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 1 |                    1 | TRUE       |     0.849 |  0.836 |  0.890 |
| candidate | count     | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 1 |                    1 | TRUE       |     0.846 |  0.817 |  0.879 |
| baseline  | projected | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 5 |                    4 | TRUE       |     0.844 |  0.816 |  0.918 |
| candidate | projected | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 5 |                    4 | TRUE       |     0.901 |  0.878 |  0.930 |
| baseline  | regions   | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                1 |                 1 |                   10 | TRUE       |     0.938 |  0.899 |  0.996 |
| candidate | regions   | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                1 |                 1 |                   10 | TRUE       |     0.905 |  0.872 |  0.980 |
| baseline  | full      | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 5 |                   10 | TRUE       |     0.957 |  0.891 |  1.038 |
| candidate | full      | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 5 |                   10 | TRUE       |     0.963 |  0.939 |  1.011 |
| baseline  | count     | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 1 |                    1 | TRUE       |     0.831 |  0.786 |  0.921 |
| candidate | count     | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 1 |                    1 | TRUE       |     0.818 |  0.774 |  0.885 |
| baseline  | projected | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 5 |                    4 | TRUE       |     0.848 |  0.823 |  0.908 |
| candidate | projected | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 5 |                    4 | TRUE       |     0.872 |  0.854 |  0.915 |
| baseline  | regions   | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                1 |                 1 |                   10 | TRUE       |     0.943 |  0.901 |  1.026 |
| candidate | regions   | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                1 |                 1 |                   10 | TRUE       |     0.961 |  0.900 |  1.017 |
| baseline  | full      | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 5 |                   10 | TRUE       |     0.986 |  0.949 |  1.017 |
| candidate | full      | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 5 |                   10 | TRUE       |     1.009 |  0.994 |  1.250 |
| baseline  | count     | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 1 |                    1 | TRUE       |     0.878 |  0.861 |  0.896 |
| candidate | count     | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 1 |                    1 | TRUE       |     0.881 |  0.864 |  0.891 |
| baseline  | projected | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 5 |                    4 | TRUE       |     0.907 |  0.856 |  0.938 |
| candidate | projected | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 5 |                    4 | TRUE       |     0.868 |  0.836 |  0.905 |
| baseline  | regions   | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                1 |                 1 |                   10 | TRUE       |     1.037 |  1.009 |  1.149 |
| candidate | regions   | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                1 |                 1 |                   10 | TRUE       |     0.942 |  0.933 |  1.002 |
| baseline  | full      | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 5 |                   10 | TRUE       |     1.002 |  0.982 |  1.117 |
| candidate | full      | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 5 |                   10 | TRUE       |     1.000 |  0.970 |  1.049 |
| baseline  | count     | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 1 |                    1 | TRUE       |     0.837 |  0.806 |  0.910 |
| candidate | count     | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 1 |                    1 | TRUE       |     0.840 |  0.749 |  0.936 |
| baseline  | projected | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 5 |                    4 | TRUE       |     0.908 |  0.863 |  0.949 |
| candidate | projected | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 5 |                    4 | TRUE       |     0.884 |  0.839 |  0.955 |
| baseline  | regions   | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                1 |                 1 |                   10 | TRUE       |     0.975 |  0.920 |  1.270 |
| candidate | regions   | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                1 |                 1 |                   10 | TRUE       |     0.952 |  0.935 |  1.073 |

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
