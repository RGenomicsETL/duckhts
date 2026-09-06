Record-major genotypes and carrier expansion on an HPRC cohort
================

The registered HPRC v2.0 GRCh38 cohort region `chr22:20000000-21000000`
retains every source sample and allele. This is not the
consequence-only, split-allele corpus. Stage explicitly with
`duckhtsbench::duckhts_bench_stage_genotypes()` before rendering;
rendering is network-free. BCF and VCF.gz contain the same complete
records, verified by the staging test and full HTSlib text comparison.

This compares `read_geno()` with the canonical
`read_bcf(..., tidy_format:=true)` GT projection. The experiment named
`read_bcf_v2` in the original request was deleted; it is not
reintroduced as a benchmark dependency. Full cohorts, the first eight
header samples, sparse non-reference calls, and downstream typed carrier
expansion have separate denominators. The input has GT but no PS; these
timings do not measure PS decoding, which SQL/R/native fixtures test.

## Source and workload identity

    ## Source revision: 48ac22cbd06a1dd19ed425489d9b1e883209ada6 
    ## Source tree: 74bd979a13e3272dd3f475437fd7235341f8cefa

    ## Extension MD5: 0bd702b6cbddb0e393f675eeeb2c9aa5

    ## R: R version 4.6.0 (2026-04-24) ; DuckDB: 1.5.3

    ## Linux 6.8.0-78-generic x86_64 GNU/Linux

    ## pid 1908847's current affinity list: 2

    ## Model name:                           13th Gen Intel(R) Core(TM) i5-13500 BIOS Model name:                      13th Gen Intel(R) Core(TM) i5-13500 To Be Filled By O.E.M. CPU @ 2.4GHz

    ## DuckDB threads: 1; scan handles: 1; HTSlib decompression workers: 0; repetitions: 3

| format | artifact        |   bytes | md5                              |
|:-------|:----------------|--------:|:---------------------------------|
| VCF    | geno_hprc_vcfgz | 6223471 | 8d7132125abf081ffe4a235de3bc0873 |
| BCF    | geno_hprc_bcf   | 6179639 | 32e9a194b62ffbdf49c141ca8721a2b7 |

| workload | input_records | selected_samples | selected_calls | selected_slots | non_reference_calls | sparse_slots | alt_slots |
|:---------|--------------:|-----------------:|---------------:|---------------:|--------------------:|-------------:|----------:|
| full     |         14365 |              232 |        3332680 |        6665360 |              437639 |       875278 |    614611 |
| selected |         14365 |                8 |         114920 |         229840 |               14078 |        28156 |     19388 |
| sparse   |         14365 |              232 |        3332680 |        6665360 |              437639 |       875278 |    614611 |
| carriers |         14365 |              232 |        3332680 |        6665360 |              437639 |       875278 |    614611 |

## Materialization and carrier expansion

Timing includes binding, decoding, selected output vectors and
`CREATE TABLE AS`. Carrier timing also includes allele-slot expansion
and ALT selection. Every measurement uses a fresh R/DuckDB process;
source data is warm in the OS cache. Startup, independent counting,
checksums, checkpoint, snapshots and exact output comparison are outside
timing. RSS is the process high-water mark immediately after
materialization. Materialized database bytes are measured after
checkpoint and include the tiny sample catalog and database overhead;
they are physical compressed storage, not vector-memory bytes. The fixed
DuckDB memory limit is 8 GiB; these are not sealed/static-allocation
measurements.

| format | reader    | workload | elapsed | records_per_second | calls_per_second | allele_slots_per_second | input_bytes_per_second | peak_rss_kib | materialized_database_bytes |
|:-------|:----------|:---------|--------:|-------------------:|-----------------:|------------------------:|-----------------------:|-------------:|----------------------------:|
| BCF    | read_bcf  | carriers |   6.806 |           2110.638 |         489667.9 |                979335.9 |               907969.3 |      8471464 |                    65548288 |
| VCF    | read_bcf  | carriers |   7.018 |           2046.879 |         474876.0 |                949752.1 |               886787.0 |      8504656 |                    65548288 |
| BCF    | read_geno | carriers |   0.901 |          15943.396 |        3698867.9 |               7397735.8 |              6858644.8 |       380488 |                    65548288 |
| VCF    | read_geno | carriers |   1.161 |          12372.954 |        2870525.4 |               5741050.8 |              5360440.1 |       401352 |                    65548288 |
| BCF    | read_bcf  | full     |  24.130 |            595.317 |         138113.6 |                276227.1 |               256097.8 |      7556348 |                  5269106688 |
| VCF    | read_bcf  | full     |  22.839 |            628.968 |         145920.6 |                291841.1 |               272493.1 |      7593976 |                  5269368832 |
| BCF    | read_geno | full     |   0.938 |          15314.499 |        3552963.8 |               7105927.5 |              6588101.3 |       433280 |                    30420992 |
| VCF    | read_geno | full     |   0.709 |          20260.931 |        4700536.0 |               9401071.9 |              8777815.2 |       434200 |                    30420992 |
| BCF    | read_bcf  | selected |   0.953 |          15073.452 |         120587.6 |                241175.2 |              6484406.1 |       617888 |                   191901696 |
| VCF    | read_bcf  | selected |   0.749 |          19178.905 |         153431.2 |                306862.5 |              8309040.1 |       617624 |                   191901696 |
| BCF    | read_geno | selected |   0.197 |          72918.782 |         583350.3 |               1166700.5 |             31368725.9 |       247312 |                    24915968 |
| VCF    | read_geno | selected |   0.233 |          61652.361 |         493218.9 |                986437.8 |             26710176.0 |       268628 |                    24915968 |
| BCF    | read_bcf  | sparse   |  24.691 |            581.791 |         134975.5 |                269951.0 |               250279.0 |     11156896 |                  5438451712 |
| VCF    | read_bcf  | sparse   |  24.651 |            582.735 |         135194.5 |                270389.0 |               252463.2 |     11175080 |                  5438451712 |
| BCF    | read_geno | sparse   |   0.309 |          46488.673 |       10785372.2 |              21570744.3 |             19998831.7 |       260216 |                    25702400 |
| VCF    | read_geno | sparse   |   0.361 |          39792.244 |        9231800.6 |              18463601.1 |             17239531.9 |       280912 |                    25702400 |

| format | reader    | workload | repetition | output_rows | output_calls | output_slots | elapsed |    cpu |
|:-------|:----------|:---------|-----------:|------------:|-------------:|-------------:|--------:|-------:|
| VCF    | read_bcf  | full     |          1 |     3332680 |      3332680 |      6665360 |  20.026 |  9.496 |
| VCF    | read_geno | full     |          1 |       14365 |      3332680 |      6665360 |   0.709 |  0.693 |
| VCF    | read_geno | full     |          2 |       14365 |      3332680 |      6665360 |   0.693 |  0.674 |
| VCF    | read_bcf  | full     |          2 |     3332680 |      3332680 |      6665360 |  22.839 |  9.411 |
| VCF    | read_bcf  | full     |          3 |     3332680 |      3332680 |      6665360 |  23.268 |  9.593 |
| VCF    | read_geno | full     |          3 |       14365 |      3332680 |      6665360 |   1.087 |  0.718 |
| VCF    | read_bcf  | selected |          1 |      114920 |       114920 |       229840 |   0.749 |  0.431 |
| VCF    | read_geno | selected |          1 |       14365 |       114920 |       229840 |   0.214 |  0.140 |
| VCF    | read_geno | selected |          2 |       14365 |       114920 |       229840 |   0.233 |  0.142 |
| VCF    | read_bcf  | selected |          2 |      114920 |       114920 |       229840 |   0.831 |  0.431 |
| VCF    | read_bcf  | selected |          3 |      114920 |       114920 |       229840 |   0.601 |  0.423 |
| VCF    | read_geno | selected |          3 |       14365 |       114920 |       229840 |   0.250 |  0.144 |
| VCF    | read_bcf  | sparse   |          1 |      437639 |       437639 |       875278 |  22.657 | 10.532 |
| VCF    | read_geno | sparse   |          1 |       14365 |       437639 |       875278 |   0.361 |  0.244 |
| VCF    | read_geno | sparse   |          2 |       14365 |       437639 |       875278 |   0.244 |  0.237 |
| VCF    | read_bcf  | sparse   |          2 |      437639 |       437639 |       875278 |  26.722 |  9.918 |
| VCF    | read_bcf  | sparse   |          3 |      437639 |       437639 |       875278 |  24.651 |  9.998 |
| VCF    | read_geno | sparse   |          3 |       14365 |       437639 |       875278 |   0.367 |  0.238 |
| VCF    | read_bcf  | carriers |          1 |      614611 |           NA |       614611 |   7.018 |  6.791 |
| VCF    | read_geno | carriers |          1 |      614611 |           NA |       614611 |   1.087 |  0.951 |
| VCF    | read_geno | carriers |          2 |      614611 |           NA |       614611 |   1.161 |  0.942 |
| VCF    | read_bcf  | carriers |          2 |      614611 |           NA |       614611 |   7.024 |  6.775 |
| VCF    | read_bcf  | carriers |          3 |      614611 |           NA |       614611 |   6.997 |  6.837 |
| VCF    | read_geno | carriers |          3 |      614611 |           NA |       614611 |   1.165 |  0.953 |
| BCF    | read_bcf  | full     |          1 |     3332680 |      3332680 |      6665360 |  20.431 |  9.138 |
| BCF    | read_geno | full     |          1 |       14365 |      3332680 |      6665360 |   1.198 |  0.644 |
| BCF    | read_geno | full     |          2 |       14365 |      3332680 |      6665360 |   0.634 |  0.622 |
| BCF    | read_bcf  | full     |          2 |     3332680 |      3332680 |      6665360 |  25.406 |  8.989 |
| BCF    | read_bcf  | full     |          3 |     3332680 |      3332680 |      6665360 |  24.130 |  9.036 |
| BCF    | read_geno | full     |          3 |       14365 |      3332680 |      6665360 |   0.938 |  0.658 |
| BCF    | read_bcf  | selected |          1 |      114920 |       114920 |       229840 |   0.901 |  0.418 |
| BCF    | read_geno | selected |          1 |       14365 |       114920 |       229840 |   0.197 |  0.117 |
| BCF    | read_geno | selected |          2 |       14365 |       114920 |       229840 |   0.223 |  0.115 |
| BCF    | read_bcf  | selected |          2 |      114920 |       114920 |       229840 |   0.953 |  0.405 |
| BCF    | read_bcf  | selected |          3 |      114920 |       114920 |       229840 |   1.003 |  0.409 |
| BCF    | read_geno | selected |          3 |       14365 |       114920 |       229840 |   0.192 |  0.116 |
| BCF    | read_bcf  | sparse   |          1 |      437639 |       437639 |       875278 |  27.979 | 10.148 |
| BCF    | read_geno | sparse   |          1 |       14365 |       437639 |       875278 |   0.309 |  0.189 |
| BCF    | read_geno | sparse   |          2 |       14365 |       437639 |       875278 |   0.314 |  0.187 |
| BCF    | read_bcf  | sparse   |          2 |      437639 |       437639 |       875278 |  20.561 |  9.483 |
| BCF    | read_bcf  | sparse   |          3 |      437639 |       437639 |       875278 |  24.691 |  9.594 |
| BCF    | read_geno | sparse   |          3 |       14365 |       437639 |       875278 |   0.275 |  0.185 |
| BCF    | read_bcf  | carriers |          1 |      614611 |           NA |       614611 |   6.886 |  6.698 |
| BCF    | read_geno | carriers |          1 |      614611 |           NA |       614611 |   1.070 |  0.909 |
| BCF    | read_geno | carriers |          2 |      614611 |           NA |       614611 |   0.901 |  0.895 |
| BCF    | read_bcf  | carriers |          2 |      614611 |           NA |       614611 |   6.806 |  6.766 |
| BCF    | read_bcf  | carriers |          3 |      614611 |           NA |       614611 |   6.670 |  6.647 |
| BCF    | read_geno | carriers |          3 |      614611 |           NA |       614611 |   0.878 |  0.870 |

## Complete output comparison

First-run results are normalized outside the timer to typed variant
identity, original-header sample index, every allele slot, every phase
flag and phase set. For the VCF 4.2 source, the comparison reconstructs
the leading phase flag from HTSlib’s string representation; VCF 4.4
explicit-prefix semantics are tested in the committed SQL/R witnesses,
not inferred from this dataset. Every field and duplicate contributes to
`EXCEPT ALL` in both directions. Carrier comparisons retain variant
identity, sample, allele slot/index, selected ALT, phase and PS. Each
repetition separately checks the complete materialized-row checksum.

| format | workload | different_typed_rows |
|:-------|:---------|---------------------:|
| VCF    | full     |                    0 |
| VCF    | selected |                    0 |
| VCF    | sparse   |                    0 |
| VCF    | carriers |                    0 |
| BCF    | full     |                    0 |
| BCF    | selected |                    0 |
| BCF    | sparse   |                    0 |
| BCF    | carriers |                    0 |

The nearest earlier full-materialization report is [the shared BCF
scanner report](benchmark_bcf_shared_scan.md): a single-sample GIAB
workload, not an identical cohort baseline. This report compares the two
current interfaces, not pre-change and post-change revisions. The
separately rendered paired BCF regression workload measures changes to
the existing reader. These results do not establish universal
fastest-reader performance, remote-I/O throughput, multi-worker scaling,
PS throughput, or phased annotation correctness.
