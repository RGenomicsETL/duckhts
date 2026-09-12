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

## Source and workload identity

    ## Source revision: 615df08f86d314befb00eca08e9f8d489dbcd314 
    ## Source tree: e63b99672649882467c07c4f6af0fee3f10ba867

    ## Extension MD5: 567e6bca77b5016a8085f9e495b7c3ac

    ## R: R version 4.6.0 (2026-04-24) ; DuckDB: 1.5.3

    ## Linux 6.8.0-78-generic x86_64 GNU/Linux

    ## pid 1250369's current affinity list: 2

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
| BCF    | read_bcf  | carriers |   6.911 |           2078.570 |         482228.3 |                964456.7 |               894174.4 |      8472084 |                    65810432 |
| VCF    | read_bcf  | carriers |   6.716 |           2138.922 |         496229.9 |                992459.8 |               926663.3 |      8505152 |                    65810432 |
| BCF    | read_geno | carriers |   1.094 |          13130.713 |        3046325.4 |               6092650.8 |              5648664.5 |       381012 |                    65810432 |
| VCF    | read_geno | carriers |   1.012 |          14194.664 |        3293162.1 |               6586324.1 |              6149674.9 |       401620 |                    65810432 |
| BCF    | read_bcf  | full     |  22.852 |            628.610 |         145837.6 |                291675.1 |               270420.1 |      7556880 |                  5269106688 |
| VCF    | read_bcf  | full     |  22.894 |            627.457 |         145570.0 |                291140.0 |               271838.5 |      7594908 |                  5269106688 |
| BCF    | read_geno | full     |   1.002 |          14336.327 |        3326027.9 |               6652055.9 |              6167304.4 |       433768 |                    30420992 |
| VCF    | read_geno | full     |   0.730 |          19678.082 |        4565315.1 |               9130630.1 |              8525302.7 |       434676 |                    30420992 |
| BCF    | read_bcf  | selected |   0.919 |          15631.121 |         125049.0 |                250097.9 |              6724307.9 |       618624 |                   192163840 |
| VCF    | read_bcf  | selected |   0.939 |          15298.190 |         122385.5 |                244771.0 |              6627764.6 |       618348 |                   192163840 |
| BCF    | read_geno | selected |   0.205 |          70073.171 |         560585.4 |               1121170.7 |             30144580.5 |       248088 |                    24915968 |
| VCF    | read_geno | selected |   0.223 |          64417.040 |         515336.3 |               1030672.6 |             27907941.7 |       269272 |                    24915968 |
| BCF    | read_bcf  | sparse   |  27.149 |            529.117 |         122755.2 |                245510.3 |               227619.4 |     11157964 |                  5438451712 |
| VCF    | read_bcf  | sparse   |  25.340 |            566.890 |         131518.5 |                263037.1 |               245598.7 |     11175776 |                  5438451712 |
| BCF    | read_geno | sparse   |   0.296 |          48530.405 |       11259054.1 |              22518108.1 |             20877158.8 |       260876 |                    25702400 |
| VCF    | read_geno | sparse   |   0.371 |          38719.677 |        8982965.0 |              17965929.9 |             16774854.4 |       281376 |                    25702400 |

| format | reader    | workload | repetition | output_rows | output_calls | output_slots | elapsed |    cpu |
|:-------|:----------|:---------|-----------:|------------:|-------------:|-------------:|--------:|-------:|
| VCF    | read_bcf  | full     |          1 |     3332680 |      3332680 |      6665360 |  22.894 | 10.063 |
| VCF    | read_geno | full     |          1 |       14365 |      3332680 |      6665360 |   0.730 |  0.699 |
| VCF    | read_geno | full     |          2 |       14365 |      3332680 |      6665360 |   0.700 |  0.682 |
| VCF    | read_bcf  | full     |          2 |     3332680 |      3332680 |      6665360 |  22.874 | 10.078 |
| VCF    | read_bcf  | full     |          3 |     3332680 |      3332680 |      6665360 |  26.302 | 10.506 |
| VCF    | read_geno | full     |          3 |       14365 |      3332680 |      6665360 |   1.054 |  0.720 |
| VCF    | read_bcf  | selected |          1 |      114920 |       114920 |       229840 |   0.970 |  0.447 |
| VCF    | read_geno | selected |          1 |       14365 |       114920 |       229840 |   0.232 |  0.140 |
| VCF    | read_geno | selected |          2 |       14365 |       114920 |       229840 |   0.223 |  0.139 |
| VCF    | read_bcf  | selected |          2 |      114920 |       114920 |       229840 |   0.939 |  0.419 |
| VCF    | read_bcf  | selected |          3 |      114920 |       114920 |       229840 |   0.914 |  0.445 |
| VCF    | read_geno | selected |          3 |       14365 |       114920 |       229840 |   0.219 |  0.140 |
| VCF    | read_bcf  | sparse   |          1 |      437639 |       437639 |       875278 |  21.312 |  9.777 |
| VCF    | read_geno | sparse   |          1 |       14365 |       437639 |       875278 |   0.408 |  0.257 |
| VCF    | read_geno | sparse   |          2 |       14365 |       437639 |       875278 |   0.354 |  0.243 |
| VCF    | read_bcf  | sparse   |          2 |      437639 |       437639 |       875278 |  35.150 |  9.835 |
| VCF    | read_bcf  | sparse   |          3 |      437639 |       437639 |       875278 |  25.340 |  9.615 |
| VCF    | read_geno | sparse   |          3 |       14365 |       437639 |       875278 |   0.371 |  0.242 |
| VCF    | read_bcf  | carriers |          1 |      614611 |           NA |       614611 |   7.090 |  6.816 |
| VCF    | read_geno | carriers |          1 |      614611 |           NA |       614611 |   1.129 |  0.943 |
| VCF    | read_geno | carriers |          2 |      614611 |           NA |       614611 |   1.012 |  0.934 |
| VCF    | read_bcf  | carriers |          2 |      614611 |           NA |       614611 |   6.716 |  6.669 |
| VCF    | read_bcf  | carriers |          3 |      614611 |           NA |       614611 |   6.715 |  6.694 |
| VCF    | read_geno | carriers |          3 |      614611 |           NA |       614611 |   0.966 |  0.950 |
| BCF    | read_bcf  | full     |          1 |     3332680 |      3332680 |      6665360 |  22.852 |  9.076 |
| BCF    | read_geno | full     |          1 |       14365 |      3332680 |      6665360 |   1.096 |  0.674 |
| BCF    | read_geno | full     |          2 |       14365 |      3332680 |      6665360 |   0.636 |  0.616 |
| BCF    | read_bcf  | full     |          2 |     3332680 |      3332680 |      6665360 |  22.489 |  9.079 |
| BCF    | read_bcf  | full     |          3 |     3332680 |      3332680 |      6665360 |  28.435 |  9.220 |
| BCF    | read_geno | full     |          3 |       14365 |      3332680 |      6665360 |   1.002 |  0.656 |
| BCF    | read_bcf  | selected |          1 |      114920 |       114920 |       229840 |   0.949 |  0.437 |
| BCF    | read_geno | selected |          1 |       14365 |       114920 |       229840 |   0.208 |  0.116 |
| BCF    | read_geno | selected |          2 |       14365 |       114920 |       229840 |   0.202 |  0.112 |
| BCF    | read_bcf  | selected |          2 |      114920 |       114920 |       229840 |   0.494 |  0.375 |
| BCF    | read_bcf  | selected |          3 |      114920 |       114920 |       229840 |   0.919 |  0.405 |
| BCF    | read_geno | selected |          3 |       14365 |       114920 |       229840 |   0.205 |  0.114 |
| BCF    | read_bcf  | sparse   |          1 |      437639 |       437639 |       875278 |  27.967 |  9.797 |
| BCF    | read_geno | sparse   |          1 |       14365 |       437639 |       875278 |   0.197 |  0.184 |
| BCF    | read_geno | sparse   |          2 |       14365 |       437639 |       875278 |   0.307 |  0.176 |
| BCF    | read_bcf  | sparse   |          2 |      437639 |       437639 |       875278 |  25.824 |  9.925 |
| BCF    | read_bcf  | sparse   |          3 |      437639 |       437639 |       875278 |  27.149 |  9.679 |
| BCF    | read_geno | sparse   |          3 |       14365 |       437639 |       875278 |   0.296 |  0.192 |
| BCF    | read_bcf  | carriers |          1 |      614611 |           NA |       614611 |   6.942 |  6.708 |
| BCF    | read_geno | carriers |          1 |      614611 |           NA |       614611 |   1.094 |  0.896 |
| BCF    | read_geno | carriers |          2 |      614611 |           NA |       614611 |   1.118 |  0.895 |
| BCF    | read_bcf  | carriers |          2 |      614611 |           NA |       614611 |   6.911 |  6.730 |
| BCF    | read_bcf  | carriers |          3 |      614611 |           NA |       614611 |   6.736 |  6.721 |
| BCF    | read_geno | carriers |          3 |      614611 |           NA |       614611 |   0.910 |  0.901 |

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
