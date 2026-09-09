Selected FORMAT materialization in record-major genotype calls
================

This measures `read_geno` on every record of the registered GIAB HG002
GRCh38 NIST v4.2.1 benchmark VCF.gz. The single sample has declared GT,
PS, AD, DP and GQ fields. GT/PS, +AD and +AD/DP/GQ have separate
measurements; each selection also runs with `calls` unprojected. No
record, allele or sample is filtered.

Stage the input explicitly with
`duckhtsbench::duckhts_bench_fetch("variantkey_giab_hg002_v421")`.
Rendering is network-free and resolves the input through the artifact
registry. The output is a single-sample workload, not a large-cohort or
remote-I/O result.

## Source and denominators

    ## Source revision: c4c848553a5c7abc64560b7e7a145b4a7886f67f
    ## Source tree: 8a0e318e6e27be8c5e96c0544689629af4ee4ac6

    ## Extension SHA-256: 1905733b5a8b96274e1232c1f2f20749f8412f739eb2219c4cfbdf247f5b2b5b  /root/duckhts/build/release/duckhts.duckdb_extension

    ## Input artifact: variantkey_giab_hg002_v421 ; bytes: 156252944 ; observed MD5: dc750b3807d4af1f7ffec852e9c2f771

    ## bcftools: bcftools 1.23.1-70-g6dbd8fef
    ## R: R version 4.6.0 (2026-04-24) ; DuckDB: 1.5.3

    ## Linux 6.8.0-78-generic x86_64 GNU/Linux

    ## pid 1312381's current affinity list: 2

    ## Model name:                           13th Gen Intel(R) Core(TM) i5-13500 BIOS Model name:                      13th Gen Intel(R) Core(TM) i5-13500 To Be Filled By O.E.M. CPU @ 2.4GHz

    ## DuckDB threads: 1; scan handles: 1; HTSlib decompression workers: 0; repetitions: 3

| denominator |   count |
|:------------|--------:|
| records     | 4048342 |
| calls       | 4048342 |
| gt_slots    | 8096684 |
| ps_values   |       0 |
| ad_slots    | 8144465 |
| ad_values   | 8144465 |
| dp_values   | 4022526 |
| gq_values   | 4022526 |

Input denominators come from complete `bcftools view -H` records in
bounded R batches. Every physical duplicate counts. AD slots include
missing list elements; AD values count only nonmissing elements. DP/GQ
values and PS values count nonmissing scalars. An absent AD list has no
child slots. An emitted scalar field occupies one slot per call even
when NULL. No PS throughput claim is made when its observed count is
zero.

## Materialization cost

Each timer includes bind, sequential read/decompression, decoding,
vector writing and `CREATE TABLE AS` in a fresh R/DuckDB process with an
8 GiB memory limit. Calls-unprojected cases still materialize all site
columns and scan-local ordinals; they are not index-only counts.
Downstream counts, checksums, checkpoints, snapshots and comparisons are
outside timing. Peak RSS is the process high-water mark immediately
after CTAS, not first-party workspace size. Database bytes are measured
after checkpoint. Case order reverses on alternate repetitions; input
data is warm in the OS cache.

| selection      | calls_projected | elapsed |   cpu | peak_rss_kib | database_bytes | records |   calls | gt_slots | extra_slots | ad_values | dp_values | gq_values |
|:---------------|:----------------|--------:|------:|-------------:|---------------:|--------:|--------:|---------:|------------:|----------:|----------:|----------:|
| GT_PS          | FALSE           |   5.768 | 5.748 |       206028 |       27275264 | 4048342 |       0 |        0 |           0 |         0 |         0 |         0 |
| GT_PS_AD       | FALSE           |   5.751 | 5.734 |       209968 |       27275264 | 4048342 |       0 |        0 |           0 |         0 |         0 |         0 |
| GT_PS_AD_DP_GQ | FALSE           |   5.771 | 5.752 |       209940 |       27537408 | 4048342 |       0 |        0 |           0 |         0 |         0 |         0 |
| GT_PS          | TRUE            |   6.913 | 6.897 |       252712 |       29372416 | 4048342 | 4048342 |  8096684 |           0 |         0 |         0 |         0 |
| GT_PS_AD       | TRUE            |   7.422 | 7.396 |       277820 |       40644608 | 4048342 | 4048342 |  8096684 |     8144465 |   8144465 |         0 |         0 |
| GT_PS_AD_DP_GQ | TRUE            |   7.980 | 7.959 |       297016 |       53227520 | 4048342 | 4048342 |  8096684 |    16241149 |   8144465 |   4022526 |   4022526 |

|     | repetition | selection      | calls_projected | elapsed |   cpu | peak_rss_kib |
|:----|-----------:|:---------------|:----------------|--------:|------:|-------------:|
| 1   |          1 | GT_PS          | TRUE            |   6.913 | 6.897 |       252712 |
| 2   |          1 | GT_PS_AD       | TRUE            |   7.422 | 7.396 |       276320 |
| 3   |          1 | GT_PS_AD_DP_GQ | TRUE            |   7.973 | 7.949 |       296972 |
| 4   |          1 | GT_PS          | FALSE           |   5.831 | 5.809 |       204716 |
| 5   |          1 | GT_PS_AD       | FALSE           |   5.751 | 5.732 |       209716 |
| 6   |          1 | GT_PS_AD_DP_GQ | FALSE           |   5.826 | 5.807 |       210624 |
| 61  |          2 | GT_PS_AD_DP_GQ | FALSE           |   5.734 | 5.715 |       209940 |
| 51  |          2 | GT_PS_AD       | FALSE           |   5.751 | 5.734 |       213688 |
| 41  |          2 | GT_PS          | FALSE           |   5.763 | 5.744 |       206992 |
| 31  |          2 | GT_PS_AD_DP_GQ | TRUE            |   8.005 | 7.988 |       297016 |
| 21  |          2 | GT_PS_AD       | TRUE            |   7.385 | 7.364 |       277820 |
| 11  |          2 | GT_PS          | TRUE            |   6.939 | 6.916 |       254408 |
| 12  |          3 | GT_PS          | TRUE            |   6.890 | 6.871 |       251888 |
| 22  |          3 | GT_PS_AD       | TRUE            |   7.439 | 7.415 |       278392 |
| 32  |          3 | GT_PS_AD_DP_GQ | TRUE            |   7.980 | 7.959 |       298124 |
| 42  |          3 | GT_PS          | FALSE           |   5.768 | 5.748 |       206028 |
| 52  |          3 | GT_PS_AD       | FALSE           |   5.764 | 5.749 |       209968 |
| 62  |          3 | GT_PS_AD_DP_GQ | FALSE           |   5.771 | 5.752 |       207844 |

## Complete comparisons and scope

All six layouts preserve site fields and scan-local record ordinals.
Projected layouts preserve the complete default GT/PS call state.
Selected AD/DP/GQ values are compared with projected `read_bcf` across
every record using `EXCEPT ALL` in both directions. Per-pass
complete-row fingerprints must also agree within each layout. Four
corruption controls test that missing-slot compaction, dropped
duplicates, changed depth and NULL-to-zero substitution are rejected.

| comparison                 | differences |
|:---------------------------|------------:|
| sites GT_PS TRUE           |           0 |
| GT/PS GT_PS TRUE           |           0 |
| sites GT_PS_AD TRUE        |           0 |
| GT/PS GT_PS_AD TRUE        |           0 |
| sites GT_PS_AD_DP_GQ TRUE  |           0 |
| GT/PS GT_PS_AD_DP_GQ TRUE  |           0 |
| sites GT_PS FALSE          |           0 |
| sites GT_PS_AD FALSE       |           0 |
| sites GT_PS_AD_DP_GQ FALSE |           0 |
| read_bcf GT_PS_AD          |           0 |
| read_bcf GT_PS_AD_DP_GQ    |           0 |

    ## Rejected corruption controls: 4

This is the selected-FORMAT baseline, not a pre-change/post-change
speedup claim. The nearest prior [GIAB full-reader
report](benchmark_bcf_shared_scan.md) reads all declared columns with
`read_bcf`; it does not measure this nested schema or projection set.
The [matched HPRC report](benchmark_genotypes_scalar_counts.md)
separately exercises the existing GT-only interfaces against the
[shared-decoder baseline](benchmark_genotypes_format_shared.md) and
[earlier GT-only baseline](benchmark_genotypes.md). These workloads do
not establish universal fastest-reader performance, multi-sample FORMAT
scaling, remote-I/O throughput or sealed allocation.

## Matched GT-only cohort comparisons

The following comparison requires identical registered input hashes and
all four input-denominator rows in the rendered HPRC reports. Both use
one scan thread, no HTSlib decompression workers and three fresh-process
passes on CPU 2. Values are medians from those reports; the percentage
is descriptive, not a statistical significance test or a general
no-regression claim. The reports retain complete within-revision typed
comparisons, not cross-revision output snapshots.

    ## benchmark_genotypes_format_shared.md :     ## Source revision: 76bd0f0c90d8d52baed0278b7e7e9c21a81330a4
    ## benchmark_genotypes.md :     ## Source revision: 48ac22cbd06a1dd19ed425489d9b1e883209ada6
    ## benchmark_genotypes_scalar_counts.md :     ## Source revision: 615df08f86d314befb00eca08e9f8d489dbcd314

    ## Comparison baseline: benchmark_genotypes_format_shared.md
    ##
    ##
    ## |format |reader    |workload | elapsed_baseline| peak_rss_kib_baseline| elapsed_current| peak_rss_kib_current| elapsed_change_percent|
    ## |:------|:---------|:--------|----------------:|---------------------:|---------------:|--------------------:|----------------------:|
    ## |BCF    |read_bcf  |carriers |            6.984|               8472260|           6.911|              8472084|                 -1.045|
    ## |BCF    |read_bcf  |full     |           25.593|               7556560|          22.852|              7556880|                -10.710|
    ## |BCF    |read_bcf  |selected |            0.984|                618384|           0.919|               618624|                 -6.606|
    ## |BCF    |read_bcf  |sparse   |           27.798|              11157792|          27.149|             11157964|                 -2.335|
    ## |BCF    |read_geno |carriers |            1.145|                380944|           1.094|               381012|                 -4.454|
    ## |BCF    |read_geno |full     |            0.809|                433908|           1.002|               433768|                 23.857|
    ## |BCF    |read_geno |selected |            0.207|                248044|           0.205|               248088|                 -0.966|
    ## |BCF    |read_geno |sparse   |            0.294|                260844|           0.296|               260876|                  0.680|
    ## |VCF    |read_bcf  |carriers |            7.038|               8505344|           6.716|              8505152|                 -4.575|
    ## |VCF    |read_bcf  |full     |           26.365|               7594716|          22.894|              7594908|                -13.165|
    ## |VCF    |read_bcf  |selected |            1.054|                618488|           0.939|               618348|                -10.911|
    ## |VCF    |read_bcf  |sparse   |           27.157|              11175524|          25.340|             11175776|                 -6.691|
    ## |VCF    |read_geno |carriers |            1.196|                401164|           1.012|               401620|                -15.385|
    ## |VCF    |read_geno |full     |            1.163|                434768|           0.730|               434676|                -37.231|
    ## |VCF    |read_geno |selected |            0.231|                269312|           0.223|               269272|                 -3.463|
    ## |VCF    |read_geno |sparse   |            0.412|                281240|           0.371|               281376|                 -9.951|
    ## Comparison baseline: benchmark_genotypes.md
    ##
    ##
    ## |format |reader    |workload | elapsed_baseline| peak_rss_kib_baseline| elapsed_current| peak_rss_kib_current| elapsed_change_percent|
    ## |:------|:---------|:--------|----------------:|---------------------:|---------------:|--------------------:|----------------------:|
    ## |BCF    |read_bcf  |carriers |            6.806|               8471464|           6.911|              8472084|                  1.543|
    ## |BCF    |read_bcf  |full     |           24.130|               7556348|          22.852|              7556880|                 -5.296|
    ## |BCF    |read_bcf  |selected |            0.953|                617888|           0.919|               618624|                 -3.568|
    ## |BCF    |read_bcf  |sparse   |           24.691|              11156896|          27.149|             11157964|                  9.955|
    ## |BCF    |read_geno |carriers |            0.901|                380488|           1.094|               381012|                 21.421|
    ## |BCF    |read_geno |full     |            0.938|                433280|           1.002|               433768|                  6.823|
    ## |BCF    |read_geno |selected |            0.197|                247312|           0.205|               248088|                  4.061|
    ## |BCF    |read_geno |sparse   |            0.309|                260216|           0.296|               260876|                 -4.207|
    ## |VCF    |read_bcf  |carriers |            7.018|               8504656|           6.716|              8505152|                 -4.303|
    ## |VCF    |read_bcf  |full     |           22.839|               7593976|          22.894|              7594908|                  0.241|
    ## |VCF    |read_bcf  |selected |            0.749|                617624|           0.939|               618348|                 25.367|
    ## |VCF    |read_bcf  |sparse   |           24.651|              11175080|          25.340|             11175776|                  2.795|
    ## |VCF    |read_geno |carriers |            1.161|                401352|           1.012|               401620|                -12.834|
    ## |VCF    |read_geno |full     |            0.709|                434200|           0.730|               434676|                  2.962|
    ## |VCF    |read_geno |selected |            0.233|                268628|           0.223|               269272|                 -4.292|
    ## |VCF    |read_geno |sparse   |            0.361|                280912|           0.371|               281376|                  2.770|
