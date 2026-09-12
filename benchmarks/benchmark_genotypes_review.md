Record-major genotypes and carrier expansion on an HPRC cohort
================

The registered HPRC v2.0 GRCh38 cohort region `chr22:20000000-21000000`
retains every source sample and allele. This is not the
consequence-only, split-allele corpus. Stage explicitly with
`duckhtsbench::duckhts_bench_stage_genotypes()` before rendering;
rendering is network-free. BCF and VCF.gz contain the same complete
records, verified by the staging test and full HTSlib text comparison.

This compares `read_geno()` with the canonical
`read_bcf(..., tidy_format:=true)` GT projection. Full cohorts, the
first eight header samples, sparse non-reference calls, and downstream
typed carrier expansion have separate denominators. The input has GT but
no PS; these timings do not measure PS decoding, which SQL/R/native
fixtures test.

`benchmark_genotypes_format_shared.md` and
`benchmark_genotypes_scalar_counts.md` are revision-specific snapshots
of this Rmd. Preserve their measured revisions; render another snapshot
from the repository root with a new report name:

``` sh
BENCHMARK_RMD=benchmark_genotypes.Rmd \
BENCHMARK_REPORT=benchmark_genotypes_review.md taskset -c 2 make bench-snapshot
```

## Source and workload identity

    ## Source revision: be6c4383b7c6fc8c4e7fd5b7f1bfc2cd3181d8b1
    ## Source tree: 4a29fc0f39f8bf15fe7bb9c37f4e047c6554b2f0

    ## Extension MD5: e9288246a7e03a2159aa6bd761370552

    ## R: R version 4.6.0 (2026-04-24) ; DuckDB: 1.5.3

    ## Linux 6.8.0-78-generic x86_64 GNU/Linux

    ## pid 1538596's current affinity list: 2

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
| BCF    | read_bcf  | carriers |   6.814 |           2108.160 |         489093.0 |                978186.1 |               906903.3 |      8472244 |                    65810432 |
| VCF    | read_bcf  | carriers |   6.997 |           2053.023 |         476301.3 |                952602.5 |               889448.5 |      8505484 |                    65810432 |
| BCF    | read_geno | carriers |   1.124 |          12780.249 |        2965017.8 |               5930035.6 |              5497899.5 |       380400 |                    65810432 |
| VCF    | read_geno | carriers |   1.169 |          12288.281 |        2850881.1 |               5701762.2 |              5323756.2 |       401164 |                    66072576 |
| BCF    | read_bcf  | full     |  23.469 |            612.084 |         142003.5 |                284007.0 |               263310.7 |      7556860 |                  5269893120 |
| VCF    | read_bcf  | full     |  24.760 |            580.170 |         134599.4 |                269198.7 |               251351.8 |      7594636 |                  5269106688 |
| BCF    | read_geno | full     |   1.119 |          12837.355 |        2978266.3 |               5956532.6 |              5522465.6 |       433872 |                    30420992 |
| VCF    | read_geno | full     |   0.706 |          20347.025 |        4720509.9 |               9441019.8 |              8815114.7 |       434776 |                    30420992 |
| BCF    | read_bcf  | selected |   0.980 |          14658.163 |         117265.3 |                234530.6 |              6305754.1 |       618472 |                   192163840 |
| VCF    | read_bcf  | selected |   0.969 |          14824.561 |         118596.5 |                237193.0 |              6422570.7 |       618332 |                   192163840 |
| BCF    | read_geno | selected |   0.221 |          65000.000 |         520000.0 |               1040000.0 |             27962167.4 |       248212 |                    24915968 |
| VCF    | read_geno | selected |   0.228 |          63004.386 |         504035.1 |               1008070.2 |             27295925.4 |       269304 |                    24915968 |
| BCF    | read_bcf  | sparse   |  24.739 |            580.662 |         134713.6 |                269427.2 |               249793.4 |     11157248 |                  5438451712 |
| VCF    | read_bcf  | sparse   |  25.007 |            574.439 |         133269.9 |                266539.8 |               248869.2 |     11175412 |                  5438451712 |
| BCF    | read_geno | sparse   |   0.281 |          51120.996 |       11860071.2 |              23720142.3 |             21991597.9 |       260884 |                    25702400 |
| VCF    | read_geno | sparse   |   0.362 |          39682.320 |        9206298.3 |              18412596.7 |             17191908.8 |       281336 |                    25702400 |

| format | reader    | workload | repetition | output_rows | output_calls | output_slots | elapsed |    cpu |
|:-------|:----------|:---------|-----------:|------------:|-------------:|-------------:|--------:|-------:|
| VCF    | read_bcf  | full     |          1 |     3332680 |      3332680 |      6665360 |  24.642 |  9.936 |
| VCF    | read_geno | full     |          1 |       14365 |      3332680 |      6665360 |   0.706 |  0.687 |
| VCF    | read_geno | full     |          2 |       14365 |      3332680 |      6665360 |   0.696 |  0.675 |
| VCF    | read_bcf  | full     |          2 |     3332680 |      3332680 |      6665360 |  24.760 |  9.702 |
| VCF    | read_bcf  | full     |          3 |     3332680 |      3332680 |      6665360 |  26.462 |  9.956 |
| VCF    | read_geno | full     |          3 |       14365 |      3332680 |      6665360 |   1.184 |  0.730 |
| VCF    | read_bcf  | selected |          1 |      114920 |       114920 |       229840 |   0.994 |  0.468 |
| VCF    | read_geno | selected |          1 |       14365 |       114920 |       229840 |   0.270 |  0.149 |
| VCF    | read_geno | selected |          2 |       14365 |       114920 |       229840 |   0.228 |  0.143 |
| VCF    | read_bcf  | selected |          2 |      114920 |       114920 |       229840 |   0.966 |  0.434 |
| VCF    | read_bcf  | selected |          3 |      114920 |       114920 |       229840 |   0.969 |  0.458 |
| VCF    | read_geno | selected |          3 |       14365 |       114920 |       229840 |   0.151 |  0.139 |
| VCF    | read_bcf  | sparse   |          1 |      437639 |       437639 |       875278 |  25.207 |  9.730 |
| VCF    | read_geno | sparse   |          1 |       14365 |       437639 |       875278 |   0.362 |  0.251 |
| VCF    | read_geno | sparse   |          2 |       14365 |       437639 |       875278 |   0.296 |  0.241 |
| VCF    | read_bcf  | sparse   |          2 |      437639 |       437639 |       875278 |  25.007 |  9.693 |
| VCF    | read_bcf  | sparse   |          3 |      437639 |       437639 |       875278 |  23.185 |  9.752 |
| VCF    | read_geno | sparse   |          3 |       14365 |       437639 |       875278 |   0.391 |  0.250 |
| VCF    | read_bcf  | carriers |          1 |      614611 |           NA |       614611 |   6.824 |  6.762 |
| VCF    | read_geno | carriers |          1 |      614611 |           NA |       614611 |   1.169 |  0.963 |
| VCF    | read_geno | carriers |          2 |      614611 |           NA |       614611 |   1.102 |  0.957 |
| VCF    | read_bcf  | carriers |          2 |      614611 |           NA |       614611 |   7.043 |  6.825 |
| VCF    | read_bcf  | carriers |          3 |      614611 |           NA |       614611 |   6.997 |  6.834 |
| VCF    | read_geno | carriers |          3 |      614611 |           NA |       614611 |   1.217 |  0.969 |
| BCF    | read_bcf  | full     |          1 |     3332680 |      3332680 |      6665360 |  20.354 |  8.964 |
| BCF    | read_geno | full     |          1 |       14365 |      3332680 |      6665360 |   1.119 |  0.662 |
| BCF    | read_geno | full     |          2 |       14365 |      3332680 |      6665360 |   0.639 |  0.624 |
| BCF    | read_bcf  | full     |          2 |     3332680 |      3332680 |      6665360 |  23.469 |  9.192 |
| BCF    | read_bcf  | full     |          3 |     3332680 |      3332680 |      6665360 |  23.751 |  9.172 |
| BCF    | read_geno | full     |          3 |       14365 |      3332680 |      6665360 |   1.121 |  0.664 |
| BCF    | read_bcf  | selected |          1 |      114920 |       114920 |       229840 |   1.052 |  0.408 |
| BCF    | read_geno | selected |          1 |       14365 |       114920 |       229840 |   0.229 |  0.116 |
| BCF    | read_geno | selected |          2 |       14365 |       114920 |       229840 |   0.199 |  0.116 |
| BCF    | read_bcf  | selected |          2 |      114920 |       114920 |       229840 |   0.897 |  0.408 |
| BCF    | read_bcf  | selected |          3 |      114920 |       114920 |       229840 |   0.980 |  0.420 |
| BCF    | read_geno | selected |          3 |       14365 |       114920 |       229840 |   0.221 |  0.123 |
| BCF    | read_bcf  | sparse   |          1 |      437639 |       437639 |       875278 |  23.720 |  9.585 |
| BCF    | read_geno | sparse   |          1 |       14365 |       437639 |       875278 |   0.281 |  0.189 |
| BCF    | read_geno | sparse   |          2 |       14365 |       437639 |       875278 |   0.190 |  0.183 |
| BCF    | read_bcf  | sparse   |          2 |      437639 |       437639 |       875278 |  25.645 | 10.090 |
| BCF    | read_bcf  | sparse   |          3 |      437639 |       437639 |       875278 |  24.739 |  9.852 |
| BCF    | read_geno | sparse   |          3 |       14365 |       437639 |       875278 |   0.298 |  0.188 |
| BCF    | read_bcf  | carriers |          1 |      614611 |           NA |       614611 |   6.947 |  6.693 |
| BCF    | read_geno | carriers |          1 |      614611 |           NA |       614611 |   1.149 |  0.904 |
| BCF    | read_geno | carriers |          2 |      614611 |           NA |       614611 |   1.124 |  0.898 |
| BCF    | read_bcf  | carriers |          2 |      614611 |           NA |       614611 |   6.814 |  6.734 |
| BCF    | read_bcf  | carriers |          3 |      614611 |           NA |       614611 |   6.760 |  6.740 |
| BCF    | read_geno | carriers |          3 |      614611 |           NA |       614611 |   0.930 |  0.920 |

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
