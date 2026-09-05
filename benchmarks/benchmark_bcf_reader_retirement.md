DuckHTS BCF record-cache ownership
================

<!-- Reports are generated from benchmark_bcf_record_cache.Rmd. -->

This paired workload materializes **every declared column** of the
registered GIAB HG002 GRCh38 NIST v4.2.1 benchmark VCF with `read_bcf`,
in wide and tidy form. The file has one sample, so input records and
output rows have the same denominator. This exercises typed INFO, GT,
and FORMAT decoding across many output chunks, not multi-sample
expansion or CSQ.

Stage `variantkey_giab_hg002_v421` explicitly with the `duckhtsbench`
`duckhts_bench_fetch()` helper before running; rendering never downloads
data. Then build both revisions with the same toolchain and run from the
repo root:

``` sh
BCF_CACHE_BASELINE_EXTENSION=/path/to/baseline/duckhts.duckdb_extension \
BCF_CACHE_BASELINE_REV=<commit> \
Rscript -e 'rmarkdown::render("benchmarks/benchmark_bcf_record_cache.Rmd", knit_root_dir=getwd())'
```

The earlier `benchmark_bcf_record_cache.md` records the experiment at
its named source revisions. Preserve that report when measuring
retirement: render this Rmd with
`output_file = "benchmark_bcf_reader_retirement.md"`. Set
`BCF_CACHE_READERS=read_bcf,read_bcf_v2` only to reproduce the
historical comparison with two old binaries that both expose the
experimental reader.

## Revisions, input, and environment

    ##   revision                                   commit
    ##   baseline 419cba3255ef2f4c7f49031673ac2b39f4fc94f4
    ##  candidate 8b8458a70b156694160bdb81a21fb8527e2a6582
    ##                                  src_tree                    extension_md5
    ##  55a6ff4cfaf380dd488b39262cb3e752aa27d1a4 20700c4efda8dbb99cdc51dd358956f1
    ##  e62b8014a51d8564bc599048c4831d93ee3385be 92aef5253e9920fe0b59df156412991d

    ## Artifact: variantkey_giab_hg002_v421; registered bytes: 156252944

    ## Observed input MD5: dc750b3807d4af1f7ffec852e9c2f771

    ## Input records: 4048342 ; sample: HG002; output rows per scan: 4048342

    ## DuckDB threads: 1; sequential scan handles: 1; HTSlib decompression workers: 0

    ## Measured repetitions per reader/mode/revision: 3

    ## R: R version 4.6.0 (2026-04-24) ; DuckDB: 1.5.3

    ## Linux 6.8.0-78-generic x86_64 GNU/Linux

    ## Ubuntu clang version 18.1.3 (1ubuntu1)

    ## bcftools 1.23.1-70-g6dbd8fef
    ## Using htslib 1.22.1

    ## 
    ##  pid 87044's current affinity list: 0,2,4,6

    ## Model name:                           13th Gen Intel(R) Core(TM) i5-13500 BIOS Model name:                      13th Gen Intel(R) Core(TM) i5-13500 To Be Filled By O.E.M. CPU @ 2.4GHz

## Full typed materialization

Elapsed and user + system CPU times include binding, initialization,
VCF/BGZF decoding, typed vectors, and `CREATE TABLE AS` materialization.
Every repetition uses a fresh R/DuckDB process. Startup, staging,
independent counting, checksum validation, Parquet snapshot writing, and
exact comparison are outside timing. The OS page cache is warm. RSS is
the whole-process high-water mark immediately after materialization, not
an allocation count or cache-only measurement.

    ##  mode   reader  revision elapsed    cpu peak_rss_kib input_records output_rows
    ##  tidy read_bcf  baseline  15.060 15.060      3792628       4048342     4048342
    ##  wide read_bcf  baseline  14.888 14.887      3748436       4048342     4048342
    ##  tidy read_bcf candidate  14.821 14.819      3792600       4048342     4048342
    ##  wide read_bcf candidate  14.611 14.610      3748168       4048342     4048342
    ##  rows_per_second
    ##           268814
    ##           271920
    ##           273149
    ##           277075

    ##  mode   reader  revision run    rows columns elapsed    cpu
    ##  wide read_bcf  baseline   1 4048342      29  14.872 14.872
    ##  wide read_bcf candidate   1 4048342      29  14.597 14.597
    ##  wide read_bcf candidate   2 4048342      29  14.730 14.730
    ##  wide read_bcf  baseline   2 4048342      29  14.907 14.908
    ##  wide read_bcf  baseline   3 4048342      29  14.888 14.887
    ##  wide read_bcf candidate   3 4048342      29  14.611 14.610
    ##  tidy read_bcf  baseline   1 4048342      30  15.060 15.060
    ##  tidy read_bcf candidate   1 4048342      30  14.821 14.819
    ##  tidy read_bcf candidate   2 4048342      30  14.821 14.821
    ##  tidy read_bcf  baseline   2 4048342      30  15.009 15.008
    ##  tidy read_bcf  baseline   3 4048342      30  15.070 15.071
    ##  tidy read_bcf candidate   3 4048342      30  14.791 14.790

## Output agreement and scope

    ##  mode              reference                compared different_rows
    ##  wide wide_read_bcf_baseline wide_read_bcf_candidate              0
    ##  tidy tidy_read_bcf_baseline tidy_read_bcf_candidate              0

    ##  mode    rows columns                   checksum
    ##  wide 4048342      29 37346610631881826082171150
    ##  tidy 4048342      30 37336419453416803666352226

Each first-run snapshot is compared against the pre-change `read_bcf`
output with `EXCEPT ALL` in both directions. Every repetition also
checks the independent record count and a sum of hashes covering every
materialized field. Hash agreement alone is not an exact oracle; the
separate first-run comparisons check complete typed row multisets. The
registry pins the supplier byte count; the input MD5 above is observed
at measurement time, not a supplier-authenticated checksum.

The earlier [appender contig report](benchmark_bcf_appender_contigs.md)
exercises a narrow scalar appender with different data and scheduling.
It is not an identical-workload baseline; the paired pre-change build
here supplies that baseline. This report does not measure indexed or
remote scans, appender pipelining, many-sample throughput, CSQ parsing
throughput, or bounded/static allocation. It establishes neither
universal performance parity nor a fastest-reader claim.

The independent SQL/R regression fixtures separately exercise two-sample
CSQ replication, all 2,053 sample rows of a record crossing an output
chunk, missing annotations and FORMAT values, growing list/string
buffers, and changing ploidy. This workload has no CSQ: its timings do
not establish annotation correctness. The paired revisions above identify
the actual comparison; experimental API retirement and canonical
record-cache correctness are separate changes. No incorrect annotation
output is accepted as a throughput oracle.
