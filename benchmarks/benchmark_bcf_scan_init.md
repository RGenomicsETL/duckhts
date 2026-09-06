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
| baseline revision    | 4929162aa103fe703c848a9b87bfc062cd53aa15                                             |
| candidate revision   | e87abf20bd603cb0cce8dfcc3bf45a6048edc914                                             |
| candidate src tree   | 05b2af83394afd3e6a83b4ce28ab387e1d6c0bfd                                             |
| baseline binary MD5  | 37996465549069014fa3ae36bdb7589c                                                     |
| candidate binary MD5 | a8d0f2b125cb5dbd29568e1364c59f2f                                                     |
| R                    | R version 4.6.0 (2026-04-24)                                                         |
| DuckDB R package     | 1.5.3                                                                                |
| HTSlib               | 1.24                                                                                 |
| CPU affinity         | pid 1275682’s current affinity list: 0,2,4,6                                         |
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
| baseline  | full      | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     1.170 |  1.139 |  1.422 |
| candidate | full      | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     1.201 |  1.177 |  1.234 |
| baseline  | count     | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     1.067 |  1.033 |  1.247 |
| candidate | count     | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     1.162 |  1.054 |  1.249 |
| baseline  | projected | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     1.116 |  1.095 |  1.342 |
| candidate | projected | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     1.206 |  1.062 |  1.320 |
| baseline  | regions   | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     1.169 |  1.133 |  1.341 |
| candidate | regions   | BCF/CSI    | 257 contigs |              1 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     1.158 |  1.117 |  1.272 |
| baseline  | full      | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     1.187 |  1.150 |  1.438 |
| candidate | full      | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     1.268 |  1.200 |  1.512 |
| baseline  | count     | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     1.140 |  1.038 |  1.271 |
| candidate | count     | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     1.271 |  1.058 |  1.306 |
| baseline  | projected | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     1.173 |  1.065 |  1.377 |
| candidate | projected | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     1.211 |  1.072 |  1.455 |
| baseline  | regions   | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     1.176 |  1.109 |  1.377 |
| candidate | regions   | BCF/CSI    | 257 contigs |              4 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     1.266 |  1.128 |  1.288 |
| baseline  | full      | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     0.914 |  0.906 |  0.960 |
| candidate | full      | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     0.917 |  0.872 |  0.957 |
| baseline  | count     | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     0.826 |  0.786 |  0.901 |
| candidate | count     | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     0.821 |  0.776 |  0.885 |
| baseline  | projected | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     0.862 |  0.800 |  0.907 |
| candidate | projected | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     0.807 |  0.777 |  0.863 |
| baseline  | regions   | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     0.924 |  0.891 |  0.977 |
| candidate | regions   | VCF.gz/TBI | 257 contigs |              1 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     0.900 |  0.866 |  0.918 |
| baseline  | full      | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     0.889 |  0.866 |  0.924 |
| candidate | full      | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    8 | TRUE       |     0.903 |  0.866 |  0.940 |
| baseline  | count     | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     0.801 |  0.712 |  0.868 |
| candidate | count     | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 1 |                    1 | TRUE       |     0.776 |  0.721 |  0.869 |
| baseline  | projected | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     0.807 |  0.761 |  0.863 |
| candidate | projected | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                2 |                 2 |                    4 | TRUE       |     0.816 |  0.777 |  0.887 |
| baseline  | regions   | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     0.932 |  0.895 |  1.044 |
| candidate | regions   | VCF.gz/TBI | 257 contigs |              4 |                         0 |            257 |             2 |                1 |                 1 |                    8 | TRUE       |     0.927 |  0.901 |  0.976 |
| baseline  | full      | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 5 |                   10 | TRUE       |     0.936 |  0.906 |  0.974 |
| candidate | full      | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 5 |                   10 | TRUE       |     0.921 |  0.883 |  1.050 |
| baseline  | count     | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 1 |                    1 | TRUE       |     0.816 |  0.783 |  0.843 |
| candidate | count     | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 1 |                    1 | TRUE       |     0.820 |  0.785 |  0.842 |
| baseline  | projected | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 5 |                    4 | TRUE       |     0.818 |  0.779 |  0.854 |
| candidate | projected | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                5 |                 5 |                    4 | TRUE       |     0.866 |  0.836 |  0.919 |
| baseline  | regions   | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                1 |                 1 |                   10 | TRUE       |     0.931 |  0.895 |  0.968 |
| candidate | regions   | VCF.gz/TBI | full        |              1 |                         0 |              3 |             5 |                1 |                 1 |                   10 | TRUE       |     0.919 |  0.882 |  1.023 |
| baseline  | full      | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 5 |                   10 | TRUE       |     0.934 |  0.882 |  1.000 |
| candidate | full      | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 5 |                   10 | TRUE       |     0.887 |  0.859 |  0.967 |
| baseline  | count     | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 1 |                    1 | TRUE       |     0.767 |  0.725 |  0.855 |
| candidate | count     | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 1 |                    1 | TRUE       |     0.812 |  0.783 |  0.860 |
| baseline  | projected | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 5 |                    4 | TRUE       |     0.805 |  0.784 |  0.868 |
| candidate | projected | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                5 |                 5 |                    4 | TRUE       |     0.841 |  0.797 |  0.874 |
| baseline  | regions   | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                1 |                 1 |                   10 | TRUE       |     0.960 |  0.918 |  1.040 |
| candidate | regions   | VCF.gz/TBI | full        |              4 |                         0 |              3 |             5 |                1 |                 1 |                   10 | TRUE       |     0.941 |  0.885 |  1.000 |
| baseline  | full      | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 5 |                   10 | TRUE       |     0.923 |  0.904 |  0.979 |
| candidate | full      | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 5 |                   10 | TRUE       |     0.955 |  0.894 |  1.031 |
| baseline  | count     | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 1 |                    1 | TRUE       |     0.821 |  0.795 |  0.901 |
| candidate | count     | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 1 |                    1 | TRUE       |     0.854 |  0.816 |  0.878 |
| baseline  | projected | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 5 |                    4 | TRUE       |     0.883 |  0.842 |  0.928 |
| candidate | projected | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                5 |                 5 |                    4 | TRUE       |     0.837 |  0.821 |  0.890 |
| baseline  | regions   | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                1 |                 1 |                   10 | TRUE       |     0.909 |  0.879 |  0.956 |
| candidate | regions   | VCF.gz/TBI | partial     |              1 |                         0 |              2 |             5 |                1 |                 1 |                   10 | TRUE       |     0.930 |  0.908 |  0.950 |
| baseline  | full      | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 5 |                   10 | TRUE       |     0.978 |  0.948 |  1.070 |
| candidate | full      | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 5 |                   10 | TRUE       |     0.954 |  0.932 |  1.017 |
| baseline  | count     | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 1 |                    1 | TRUE       |     0.827 |  0.794 |  0.918 |
| candidate | count     | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 1 |                    1 | TRUE       |     0.822 |  0.784 |  0.869 |
| baseline  | projected | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 5 |                    4 | TRUE       |     0.919 |  0.847 |  1.070 |
| candidate | projected | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                5 |                 5 |                    4 | TRUE       |     0.884 |  0.847 |  0.929 |
| baseline  | regions   | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                1 |                 1 |                   10 | TRUE       |     0.904 |  0.865 |  0.986 |
| candidate | regions   | VCF.gz/TBI | partial     |              4 |                         0 |              2 |             5 |                1 |                 1 |                   10 | TRUE       |     0.936 |  0.899 |  0.991 |
| baseline  | full      | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 5 |                   10 | TRUE       |     1.006 |  0.985 |  1.047 |
| candidate | full      | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 5 |                   10 | TRUE       |     0.988 |  0.946 |  1.044 |
| baseline  | count     | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 1 |                    1 | TRUE       |     0.881 |  0.873 |  0.895 |
| candidate | count     | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 1 |                    1 | TRUE       |     0.852 |  0.839 |  0.887 |
| baseline  | projected | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 5 |                    4 | TRUE       |     0.884 |  0.874 |  0.930 |
| candidate | projected | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                5 |                 5 |                    4 | TRUE       |     0.887 |  0.847 |  0.921 |
| baseline  | regions   | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                1 |                 1 |                   10 | TRUE       |     0.971 |  0.911 |  0.982 |
| candidate | regions   | VCF.gz/TBI | none        |              1 |                         0 |              0 |             5 |                1 |                 1 |                   10 | TRUE       |     0.938 |  0.896 |  1.016 |
| baseline  | full      | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 5 |                   10 | TRUE       |     0.992 |  0.968 |  1.070 |
| candidate | full      | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 5 |                   10 | TRUE       |     0.968 |  0.951 |  1.034 |
| baseline  | count     | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 1 |                    1 | TRUE       |     0.864 |  0.781 |  0.904 |
| candidate | count     | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 1 |                    1 | TRUE       |     0.836 |  0.819 |  0.904 |
| baseline  | projected | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 5 |                    4 | TRUE       |     0.897 |  0.872 |  0.926 |
| candidate | projected | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                5 |                 5 |                    4 | TRUE       |     0.881 |  0.865 |  0.933 |
| baseline  | regions   | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                1 |                 1 |                   10 | TRUE       |     0.945 |  0.935 |  1.008 |
| candidate | regions   | VCF.gz/TBI | none        |              4 |                         0 |              0 |             5 |                1 |                 1 |                   10 | TRUE       |     0.919 |  0.903 |  0.975 |

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
