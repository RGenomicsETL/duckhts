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

    ## Source revision: 76bd0f0c90d8d52baed0278b7e7e9c21a81330a4 
    ## Source tree: f1078b0250b804f7a403f977181ea4f29b265e95

    ## Extension MD5: 35f62d3595a03a112e1903d6b0d7eda5

    ## R: R version 4.6.0 (2026-04-24) ; DuckDB: 1.5.3

    ## Linux 6.8.0-78-generic x86_64 GNU/Linux

    ## pid 1169607's current affinity list: 2

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
| BCF    | read_bcf  | carriers |   6.984 |           2056.844 |         477187.9 |                954375.7 |               884828.0 |      8472260 |                    65810432 |
| VCF    | read_bcf  | carriers |   7.038 |           2041.063 |         473526.6 |                947053.1 |               884267.0 |      8505344 |                    65810432 |
| BCF    | read_geno | carriers |   1.145 |          12545.852 |        2910637.6 |               5821275.1 |              5397064.6 |       380944 |                    65810432 |
| VCF    | read_geno | carriers |   1.196 |          12010.870 |        2786521.7 |               5573043.5 |              5203571.1 |       401164 |                    65810432 |
| BCF    | read_bcf  | full     |  25.593 |            561.286 |         130218.4 |                260436.8 |               241458.2 |      7556560 |                  5269368832 |
| VCF    | read_bcf  | full     |  26.365 |            544.851 |         126405.5 |                252810.9 |               236050.5 |      7594716 |                  5269630976 |
| BCF    | read_geno | full     |   0.809 |          17756.489 |        4119505.6 |               8239011.1 |              7638614.3 |       433908 |                    30420992 |
| VCF    | read_geno | full     |   1.163 |          12351.677 |        2865589.0 |               5731178.0 |              5351221.8 |       434768 |                    30420992 |
| BCF    | read_bcf  | selected |   0.984 |          14598.577 |         116788.6 |                233577.2 |              6280120.9 |       618384 |                   192163840 |
| VCF    | read_bcf  | selected |   1.054 |          13629.032 |         109032.3 |                218064.5 |              5904621.4 |       618488 |                   192163840 |
| BCF    | read_geno | selected |   0.207 |          69396.135 |         555169.1 |               1110338.2 |             29853328.5 |       248044 |                    24915968 |
| VCF    | read_geno | selected |   0.231 |          62186.147 |         497489.2 |                994978.4 |             26941432.9 |       269312 |                    24915968 |
| BCF    | read_bcf  | sparse   |  27.798 |            516.764 |         119889.2 |                239778.4 |               222305.2 |     11157792 |                  5438451712 |
| VCF    | read_bcf  | sparse   |  27.157 |            528.961 |         122719.0 |                245438.0 |               229166.4 |     11175524 |                  5438451712 |
| BCF    | read_geno | sparse   |   0.294 |          48860.544 |       11335646.3 |              22671292.5 |             21019180.3 |       260844 |                    25702400 |
| VCF    | read_geno | sparse   |   0.412 |          34866.505 |        8089029.1 |              16178058.3 |             15105512.1 |       281240 |                    25702400 |

| format | reader    | workload | repetition | output_rows | output_calls | output_slots | elapsed |    cpu |
|:-------|:----------|:---------|-----------:|------------:|-------------:|-------------:|--------:|-------:|
| VCF    | read_bcf  | full     |          1 |     3332680 |      3332680 |      6665360 |  26.365 | 10.258 |
| VCF    | read_geno | full     |          1 |       14365 |      3332680 |      6665360 |   1.163 |  0.717 |
| VCF    | read_geno | full     |          2 |       14365 |      3332680 |      6665360 |   0.707 |  0.684 |
| VCF    | read_bcf  | full     |          2 |     3332680 |      3332680 |      6665360 |  25.761 |  9.637 |
| VCF    | read_bcf  | full     |          3 |     3332680 |      3332680 |      6665360 |  27.840 |  9.780 |
| VCF    | read_geno | full     |          3 |       14365 |      3332680 |      6665360 |   1.286 |  0.824 |
| VCF    | read_bcf  | selected |          1 |      114920 |       114920 |       229840 |   1.054 |  0.478 |
| VCF    | read_geno | selected |          1 |       14365 |       114920 |       229840 |   0.268 |  0.145 |
| VCF    | read_geno | selected |          2 |       14365 |       114920 |       229840 |   0.215 |  0.144 |
| VCF    | read_bcf  | selected |          2 |      114920 |       114920 |       229840 |   1.177 |  0.483 |
| VCF    | read_bcf  | selected |          3 |      114920 |       114920 |       229840 |   0.834 |  0.456 |
| VCF    | read_geno | selected |          3 |       14365 |       114920 |       229840 |   0.231 |  0.144 |
| VCF    | read_bcf  | sparse   |          1 |      437639 |       437639 |       875278 |  27.009 | 10.732 |
| VCF    | read_geno | sparse   |          1 |       14365 |       437639 |       875278 |   0.376 |  0.251 |
| VCF    | read_geno | sparse   |          2 |       14365 |       437639 |       875278 |   0.412 |  0.246 |
| VCF    | read_bcf  | sparse   |          2 |      437639 |       437639 |       875278 |  27.157 | 10.248 |
| VCF    | read_bcf  | sparse   |          3 |      437639 |       437639 |       875278 |  28.434 | 10.074 |
| VCF    | read_geno | sparse   |          3 |       14365 |       437639 |       875278 |   0.444 |  0.265 |
| VCF    | read_bcf  | carriers |          1 |      614611 |           NA |       614611 |   7.074 |  6.848 |
| VCF    | read_geno | carriers |          1 |      614611 |           NA |       614611 |   1.130 |  0.957 |
| VCF    | read_geno | carriers |          2 |      614611 |           NA |       614611 |   1.224 |  0.960 |
| VCF    | read_bcf  | carriers |          2 |      614611 |           NA |       614611 |   6.971 |  6.760 |
| VCF    | read_bcf  | carriers |          3 |      614611 |           NA |       614611 |   7.038 |  6.770 |
| VCF    | read_geno | carriers |          3 |      614611 |           NA |       614611 |   1.196 |  0.953 |
| BCF    | read_bcf  | full     |          1 |     3332680 |      3332680 |      6665360 |  25.593 |  9.324 |
| BCF    | read_geno | full     |          1 |       14365 |      3332680 |      6665360 |   1.080 |  0.666 |
| BCF    | read_geno | full     |          2 |       14365 |      3332680 |      6665360 |   0.809 |  0.796 |
| BCF    | read_bcf  | full     |          2 |     3332680 |      3332680 |      6665360 |  27.119 | 11.726 |
| BCF    | read_bcf  | full     |          3 |     3332680 |      3332680 |      6665360 |  25.139 |  9.237 |
| BCF    | read_geno | full     |          3 |       14365 |      3332680 |      6665360 |   0.768 |  0.660 |
| BCF    | read_bcf  | selected |          1 |      114920 |       114920 |       229840 |   0.984 |  0.403 |
| BCF    | read_geno | selected |          1 |       14365 |       114920 |       229840 |   0.250 |  0.118 |
| BCF    | read_geno | selected |          2 |       14365 |       114920 |       229840 |   0.188 |  0.115 |
| BCF    | read_bcf  | selected |          2 |      114920 |       114920 |       229840 |   0.959 |  0.436 |
| BCF    | read_bcf  | selected |          3 |      114920 |       114920 |       229840 |   1.002 |  0.410 |
| BCF    | read_geno | selected |          3 |       14365 |       114920 |       229840 |   0.207 |  0.123 |
| BCF    | read_bcf  | sparse   |          1 |      437639 |       437639 |       875278 |  27.798 |  9.484 |
| BCF    | read_geno | sparse   |          1 |       14365 |       437639 |       875278 |   0.333 |  0.198 |
| BCF    | read_geno | sparse   |          2 |       14365 |       437639 |       875278 |   0.292 |  0.187 |
| BCF    | read_bcf  | sparse   |          2 |      437639 |       437639 |       875278 |  22.732 |  9.914 |
| BCF    | read_bcf  | sparse   |          3 |      437639 |       437639 |       875278 |  33.258 |  9.899 |
| BCF    | read_geno | sparse   |          3 |       14365 |       437639 |       875278 |   0.294 |  0.189 |
| BCF    | read_bcf  | carriers |          1 |      614611 |           NA |       614611 |   6.984 |  6.758 |
| BCF    | read_geno | carriers |          1 |      614611 |           NA |       614611 |   1.145 |  0.920 |
| BCF    | read_geno | carriers |          2 |      614611 |           NA |       614611 |   1.190 |  0.950 |
| BCF    | read_bcf  | carriers |          2 |      614611 |           NA |       614611 |   6.996 |  6.782 |
| BCF    | read_bcf  | carriers |          3 |      614611 |           NA |       614611 |   6.940 |  6.712 |
| BCF    | read_geno | carriers |          3 |      614611 |           NA |       614611 |   1.115 |  0.896 |

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
