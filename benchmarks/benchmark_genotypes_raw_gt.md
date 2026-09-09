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

    ## Source revision: f99a34c3bbc9549aea2773e6c8c5f2a600a6eb0f
    ## Source tree: a0b18d333f0e1a6dd446ab042d8067db60e3ec4f

    ## Extension MD5: 674bcbc743a651b8fdcb6e199536141b

    ## R: R version 4.6.0 (2026-04-24) ; DuckDB: 1.5.3

    ## Linux 6.8.0-78-generic x86_64 GNU/Linux

    ## pid 1968656's current affinity list: 2

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
| BCF    | read_bcf  | carriers |   6.935 |           2071.377 |         480559.5 |                961119.0 |               891079.9 |      8472356 |                    65810432 |
| VCF    | read_bcf  | carriers |   6.931 |           2072.573 |         480836.8 |                961673.6 |               897918.2 |      8505328 |                    65810432 |
| BCF    | read_geno | carriers |   1.140 |          12600.877 |        2923403.5 |               5846807.0 |              5420736.0 |       381044 |                    65810432 |
| VCF    | read_geno | carriers |   1.034 |          13892.650 |        3223094.8 |               6446189.6 |              6018830.8 |       401584 |                    65810432 |
| BCF    | read_bcf  | full     |  22.893 |            627.484 |         145576.4 |                291152.8 |               269935.7 |      7557136 |                  5269368832 |
| VCF    | read_bcf  | full     |  20.557 |            698.789 |         162119.0 |                324238.0 |               302742.2 |      7594900 |                  5269368832 |
| BCF    | read_geno | full     |   1.125 |          12768.889 |        2962382.2 |               5924764.4 |              5493012.4 |       433760 |                    30420992 |
| VCF    | read_geno | full     |   1.079 |          13313.253 |        3088674.7 |               6177349.4 |              5767813.7 |       434808 |                    30420992 |
| BCF    | read_bcf  | selected |   0.868 |          16549.539 |         132396.3 |                264792.6 |              7119399.8 |       618400 |                   192163840 |
| VCF    | read_bcf  | selected |   1.053 |          13641.975 |         109135.8 |                218271.6 |              5910228.9 |       618560 |                   192163840 |
| BCF    | read_geno | selected |   0.201 |          71467.662 |         571741.3 |               1143482.6 |             30744472.6 |       248044 |                    24915968 |
| VCF    | read_geno | selected |   0.238 |          60357.143 |         482857.1 |                965714.3 |             26149037.8 |       269440 |                    24915968 |
| BCF    | read_bcf  | sparse   |  25.626 |            560.563 |         130050.7 |                260101.5 |               241147.2 |     11158160 |                  5438451712 |
| VCF    | read_bcf  | sparse   |  23.672 |            606.835 |         140785.7 |                281571.5 |               262904.3 |     11175796 |                  5438451712 |
| BCF    | read_geno | sparse   |   0.303 |          47409.241 |       10998943.9 |              21997887.8 |             20394848.2 |       260868 |                    25702400 |
| VCF    | read_geno | sparse   |   0.258 |          55678.295 |       12917364.3 |              25834728.7 |             24121980.6 |       281368 |                    25702400 |

| format | reader    | workload | repetition | output_rows | output_calls | output_slots | elapsed |    cpu |
|:-------|:----------|:---------|-----------:|------------:|-------------:|-------------:|--------:|-------:|
| VCF    | read_bcf  | full     |          1 |     3332680 |      3332680 |      6665360 |  18.730 |  9.287 |
| VCF    | read_geno | full     |          1 |       14365 |      3332680 |      6665360 |   1.079 |  0.745 |
| VCF    | read_geno | full     |          2 |       14365 |      3332680 |      6665360 |   0.700 |  0.681 |
| VCF    | read_bcf  | full     |          2 |     3332680 |      3332680 |      6665360 |  29.411 |  9.628 |
| VCF    | read_bcf  | full     |          3 |     3332680 |      3332680 |      6665360 |  20.557 |  9.596 |
| VCF    | read_geno | full     |          3 |       14365 |      3332680 |      6665360 |   1.115 |  0.775 |
| VCF    | read_bcf  | selected |          1 |      114920 |       114920 |       229840 |   1.053 |  0.472 |
| VCF    | read_geno | selected |          1 |       14365 |       114920 |       229840 |   0.238 |  0.142 |
| VCF    | read_geno | selected |          2 |       14365 |       114920 |       229840 |   0.217 |  0.140 |
| VCF    | read_bcf  | selected |          2 |      114920 |       114920 |       229840 |   0.914 |  0.427 |
| VCF    | read_bcf  | selected |          3 |      114920 |       114920 |       229840 |   1.117 |  0.449 |
| VCF    | read_geno | selected |          3 |       14365 |       114920 |       229840 |   0.255 |  0.145 |
| VCF    | read_bcf  | sparse   |          1 |      437639 |       437639 |       875278 |  23.434 | 10.760 |
| VCF    | read_geno | sparse   |          1 |       14365 |       437639 |       875278 |   0.258 |  0.250 |
| VCF    | read_geno | sparse   |          2 |       14365 |       437639 |       875278 |   0.247 |  0.241 |
| VCF    | read_bcf  | sparse   |          2 |      437639 |       437639 |       875278 |  24.779 | 10.303 |
| VCF    | read_bcf  | sparse   |          3 |      437639 |       437639 |       875278 |  23.672 |  9.870 |
| VCF    | read_geno | sparse   |          3 |       14365 |       437639 |       875278 |   0.381 |  0.256 |
| VCF    | read_bcf  | carriers |          1 |      614611 |           NA |       614611 |   6.931 |  6.842 |
| VCF    | read_geno | carriers |          1 |      614611 |           NA |       614611 |   1.034 |  0.969 |
| VCF    | read_geno | carriers |          2 |      614611 |           NA |       614611 |   1.226 |  0.975 |
| VCF    | read_bcf  | carriers |          2 |      614611 |           NA |       614611 |   7.049 |  6.830 |
| VCF    | read_bcf  | carriers |          3 |      614611 |           NA |       614611 |   6.730 |  6.708 |
| VCF    | read_geno | carriers |          3 |      614611 |           NA |       614611 |   0.973 |  0.962 |
| BCF    | read_bcf  | full     |          1 |     3332680 |      3332680 |      6665360 |  22.893 |  8.953 |
| BCF    | read_geno | full     |          1 |       14365 |      3332680 |      6665360 |   1.125 |  0.664 |
| BCF    | read_geno | full     |          2 |       14365 |      3332680 |      6665360 |   0.690 |  0.670 |
| BCF    | read_bcf  | full     |          2 |     3332680 |      3332680 |      6665360 |  21.615 |  9.142 |
| BCF    | read_bcf  | full     |          3 |     3332680 |      3332680 |      6665360 |  26.195 |  9.047 |
| BCF    | read_geno | full     |          3 |       14365 |      3332680 |      6665360 |   1.210 |  0.675 |
| BCF    | read_bcf  | selected |          1 |      114920 |       114920 |       229840 |   1.037 |  0.412 |
| BCF    | read_geno | selected |          1 |       14365 |       114920 |       229840 |   0.217 |  0.117 |
| BCF    | read_geno | selected |          2 |       14365 |       114920 |       229840 |   0.190 |  0.115 |
| BCF    | read_bcf  | selected |          2 |      114920 |       114920 |       229840 |   0.857 |  0.414 |
| BCF    | read_bcf  | selected |          3 |      114920 |       114920 |       229840 |   0.868 |  0.413 |
| BCF    | read_geno | selected |          3 |       14365 |       114920 |       229840 |   0.201 |  0.120 |
| BCF    | read_bcf  | sparse   |          1 |      437639 |       437639 |       875278 |  23.593 |  9.883 |
| BCF    | read_geno | sparse   |          1 |       14365 |       437639 |       875278 |   0.303 |  0.189 |
| BCF    | read_geno | sparse   |          2 |       14365 |       437639 |       875278 |   0.194 |  0.185 |
| BCF    | read_bcf  | sparse   |          2 |      437639 |       437639 |       875278 |  25.626 | 10.020 |
| BCF    | read_bcf  | sparse   |          3 |      437639 |       437639 |       875278 |  27.978 | 10.524 |
| BCF    | read_geno | sparse   |          3 |       14365 |       437639 |       875278 |   0.345 |  0.195 |
| BCF    | read_bcf  | carriers |          1 |      614611 |           NA |       614611 |   7.134 |  6.846 |
| BCF    | read_geno | carriers |          1 |      614611 |           NA |       614611 |   1.140 |  0.911 |
| BCF    | read_geno | carriers |          2 |      614611 |           NA |       614611 |   1.077 |  0.914 |
| BCF    | read_bcf  | carriers |          2 |      614611 |           NA |       614611 |   6.892 |  6.797 |
| BCF    | read_bcf  | carriers |          3 |      614611 |           NA |       614611 |   6.935 |  6.833 |
| BCF    | read_geno | carriers |          3 |      614611 |           NA |       614611 |   1.149 |  0.926 |

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
