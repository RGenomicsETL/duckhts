DuckHTS BCF record-cache ownership
================

<!-- benchmark_bcf_record_cache.md is generated from this Rmd. -->

This paired workload materializes **every declared column** of the
registered GIAB HG002 GRCh38 NIST v4.2.1 benchmark VCF, with both
`read_bcf` and `read_bcf_v2`, in wide and tidy form. The file has one
sample, so input records and output rows have the same denominator. This
exercises typed INFO, GT, and FORMAT decoding across many output chunks,
not multi-sample expansion or CSQ.

Stage `variantkey_giab_hg002_v421` explicitly with the `duckhtsbench`
`duckhts_bench_fetch()` helper before running; rendering never downloads
data. Then build both revisions with the same toolchain and run from the
repo root:

``` sh
BCF_CACHE_BASELINE_EXTENSION=/path/to/baseline/duckhts.duckdb_extension \
BCF_CACHE_BASELINE_REV=<commit> \
Rscript -e 'rmarkdown::render("benchmarks/benchmark_bcf_record_cache.Rmd", knit_root_dir=getwd())'
```

## Revisions, input, and environment

    ##   revision                                   commit
    ##   baseline 19b4c04405c98d18a69e4bf98ad5ff7f7ea30503
    ##  candidate d32de982510ae2a79c79983888bc7bd0aaba16be
    ##                                  src_tree                    extension_md5
    ##  ad13a3f1612ac69e3cd47189d2425039572c7fe6 e83dbba4a843ebd84125a71cee449b94
    ##  55a6ff4cfaf380dd488b39262cb3e752aa27d1a4 20700c4efda8dbb99cdc51dd358956f1

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
    ##  pid 13250's current affinity list: 0,2,4,6

    ## Model name:                           13th Gen Intel(R) Core(TM) i5-13500 BIOS Model name:                      13th Gen Intel(R) Core(TM) i5-13500 To Be Filled By O.E.M. CPU @ 2.4GHz

## Full typed materialization

Elapsed and user + system CPU times include binding, initialization,
VCF/BGZF decoding, typed vectors, and `CREATE TABLE AS` materialization.
Every repetition uses a fresh R/DuckDB process. Startup, staging,
independent counting, checksum validation, Parquet snapshot writing, and
exact comparison are outside timing. The OS page cache is warm. RSS is
the whole-process high-water mark immediately after materialization, not
an allocation count or cache-only measurement.

    ##  mode      reader  revision elapsed    cpu peak_rss_kib input_records
    ##  tidy    read_bcf  baseline  15.211 15.212      3792468       4048342
    ##  wide    read_bcf  baseline  15.053 15.052      3747852       4048342
    ##  tidy read_bcf_v2  baseline  15.252 15.252      3792424       4048342
    ##  wide read_bcf_v2  baseline  15.062 15.063      3747668       4048342
    ##  tidy    read_bcf candidate  15.195 15.194      3792272       4048342
    ##  wide    read_bcf candidate  14.851 14.850      3747524       4048342
    ##  tidy read_bcf_v2 candidate  15.073 15.073      3793468       4048342
    ##  wide read_bcf_v2 candidate  14.961 14.962      3747816       4048342
    ##  output_rows rows_per_second
    ##      4048342          266146
    ##      4048342          268939
    ##      4048342          265430
    ##      4048342          268779
    ##      4048342          266426
    ##      4048342          272597
    ##      4048342          268582
    ##      4048342          270593

    ##  mode      reader  revision run    rows columns elapsed    cpu
    ##  wide    read_bcf  baseline   1 4048342      29  15.039 15.038
    ##  wide    read_bcf candidate   1 4048342      29  14.851 14.850
    ##  wide    read_bcf candidate   2 4048342      29  14.838 14.838
    ##  wide    read_bcf  baseline   2 4048342      29  15.067 15.066
    ##  wide    read_bcf  baseline   3 4048342      29  15.053 15.052
    ##  wide    read_bcf candidate   3 4048342      29  14.970 14.969
    ##  wide read_bcf_v2  baseline   1 4048342      29  14.968 14.968
    ##  wide read_bcf_v2 candidate   1 4048342      29  14.961 14.962
    ##  wide read_bcf_v2 candidate   2 4048342      29  14.827 14.826
    ##  wide read_bcf_v2  baseline   2 4048342      29  15.063 15.063
    ##  wide read_bcf_v2  baseline   3 4048342      29  15.062 15.063
    ##  wide read_bcf_v2 candidate   3 4048342      29  15.044 15.042
    ##  tidy    read_bcf  baseline   1 4048342      30  15.293 15.291
    ##  tidy    read_bcf candidate   1 4048342      30  15.196 15.195
    ##  tidy    read_bcf candidate   2 4048342      30  15.195 15.194
    ##  tidy    read_bcf  baseline   2 4048342      30  15.211 15.212
    ##  tidy    read_bcf  baseline   3 4048342      30  15.204 15.203
    ##  tidy    read_bcf candidate   3 4048342      30  15.168 15.166
    ##  tidy read_bcf_v2  baseline   1 4048342      30  15.165 15.164
    ##  tidy read_bcf_v2 candidate   1 4048342      30  15.123 15.124
    ##  tidy read_bcf_v2 candidate   2 4048342      30  15.073 15.073
    ##  tidy read_bcf_v2  baseline   2 4048342      30  15.252 15.252
    ##  tidy read_bcf_v2  baseline   3 4048342      30  15.256 15.256
    ##  tidy read_bcf_v2 candidate   3 4048342      30  15.055 15.055

## Output agreement and scope

    ##  mode              reference                   compared different_rows
    ##  wide wide_read_bcf_baseline    wide_read_bcf_candidate              0
    ##  wide wide_read_bcf_baseline  wide_read_bcf_v2_baseline              0
    ##  wide wide_read_bcf_baseline wide_read_bcf_v2_candidate              0
    ##  tidy tidy_read_bcf_baseline    tidy_read_bcf_candidate              0
    ##  tidy tidy_read_bcf_baseline  tidy_read_bcf_v2_baseline              0
    ##  tidy tidy_read_bcf_baseline tidy_read_bcf_v2_candidate              0

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
buffers, and changing ploidy. The pre-change legacy reader returns
incorrect CSQ rows; that incorrect output is not used as a valid
throughput baseline here. Reader names, scan defaults, sample-selection
policies, and output ordering guarantees remain unchanged.
