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
its named source revisions. Preserve historical reports when measuring
subsequent changes: render this Rmd with a change-specific output file
(for allocation-error handling,
`output_file = "benchmark_reader_allocations.md"`). The paired
pre-change binary provides the identical-workload baseline. This
workload exercises BCF construction and materialization, not BAM
SAMPLE_ID caching or FASTA/FASTQ throughput. Set
`BCF_CACHE_READERS=read_bcf,read_bcf_v2` only to reproduce the
historical comparison with two old binaries that both expose the
experimental reader.

## Revisions, input, and environment

    ##   revision                                   commit
    ##   baseline ac0e2528850ca224cf8af3dbac415a25ce0f3123
    ##  candidate 6609e9bd93418f9bbd67b24d00c1f796d4df3370
    ##                                  src_tree                    extension_md5
    ##  59b510f543bde31e0667602ce214a5559a914db5 3b197043d39f11b500fff05a2e24d1db
    ##  74bd979a13e3272dd3f475437fd7235341f8cefa 0bd702b6cbddb0e393f675eeeb2c9aa5

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

    ## pid 1937594's current affinity list: 2

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
    ##  tidy read_bcf  baseline  15.227 15.226      3793312       4048342     4048342
    ##  wide read_bcf  baseline  15.169 15.169      3749088       4048342     4048342
    ##  tidy read_bcf candidate  15.321 15.320      3792720       4048342     4048342
    ##  wide read_bcf candidate  15.634 15.633      3748852       4048342     4048342
    ##  rows_per_second
    ##           265866
    ##           266883
    ##           264235
    ##           258945

    ##  mode   reader  revision run    rows columns elapsed    cpu
    ##  wide read_bcf  baseline   1 4048342      29  15.111 15.110
    ##  wide read_bcf candidate   1 4048342      29  15.832 15.754
    ##  wide read_bcf candidate   2 4048342      29  15.634 15.633
    ##  wide read_bcf  baseline   2 4048342      29  15.694 15.629
    ##  wide read_bcf  baseline   3 4048342      29  15.169 15.169
    ##  wide read_bcf candidate   3 4048342      29  15.330 15.328
    ##  tidy read_bcf  baseline   1 4048342      30  15.227 15.226
    ##  tidy read_bcf candidate   1 4048342      30  15.321 15.320
    ##  tidy read_bcf candidate   2 4048342      30  15.309 15.309
    ##  tidy read_bcf  baseline   2 4048342      30  15.280 15.279
    ##  tidy read_bcf  baseline   3 4048342      30  15.205 15.205
    ##  tidy read_bcf candidate   3 4048342      30  15.389 15.388

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
not establish annotation correctness. The paired revisions above
identify the actual comparison; experimental API retirement and
canonical record-cache correctness are separate changes. No incorrect
annotation output is accepted as a throughput oracle.
