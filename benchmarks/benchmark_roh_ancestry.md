Ancestry-tuned allele frequencies for ROH
================

## Evaluation protocol

This report evaluates whether ancestry-tuned allele frequencies reduce
unsupported long runs of homozygosity (ROH) in admixed 1000 Genomes trio
children without changing non-admixed controls. The acceptance criterion
is declared before calculating or comparing any arm:

- In admixed children, ancestry-tuned AF must yield fewer unsupported
  ROH of at least 1 Mb than both pooled and single-population AF.
- Its fraction of frequency-free truth-run length covered must be no
  more than 5% relatively below either comparator.
- In controls, both unsupported-ROH count and truth coverage must remain
  within 5% relative of the single-population arm. For a zero
  comparator, equality is required.

No amendment to the truth threshold has been made. A truth segment is a
merged run of windows (100 kb each) with at most one heterozygous call,
supported when at least 90% of its length meets that window criterion.
Truth calls use all biallelic SNVs in each full source VCF, not just
ancestry-reference sites. The window distribution and separation check
must be reviewed before any ROH-arm comparison; if inadequate, the
threshold may only be amended here before computing arm results. The
primary length threshold is 1 Mb, with 2 Mb and 5 Mb also reported. FROH
is total called ROH length divided by callable autosomal length.

## Data identity and processing status

The phased high-coverage VCF source is
`https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr{N}.filtered.shapeit2-duohmm-phased.vcf.gz`.
The pedigree source is
`https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt`.
The ancestry reference is registry artifact
`ancestry_reference_grch38_parquet` (GRCh38, 5,012,826 sites).

| Input                     | Identity / SHA-256                                                                                                                                                     | Status                                                                                                  |
|:--------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------|
| chr20 phased VCF          | Source SHA-256 `c2624c726c6cb27e288fe5908621125740379a3d523299dd07c13abf97645e4b`; remote ETag `2fbb7ed2-5b2d2762b9fc0`; Last-Modified `Thu, 29 Oct 2020 17:17:59 GMT` | Local source checksum verified; 3,202-sample source used to make the 377-child BCF                      |
| chr20 children BCF        | 849,143 records; SHA-256 `a53f83cebd7de6fb4661c14ecd93e5346022ff453392095a085e0b7cda116541`                                                                            | All 377 eligible children present; biallelic SNVs with INFO/AF, INFO/AF_AFR/AMR/EAS/EUR retained        |
| GRCh38 ancestry reference | SHA-256 `4982810aa69b9c71924e9a02dcd454ccedccc34bd8bcf71506ce9752589847fa`                                                                                             | Chr20 wide rows converted to long form; 108,757 non-palindromic, matched dosage loci per child          |
| Pedigree                  | SHA-256 `4e164b717433fc9bfc16e29b20e5b4659b4216185b3275322d421fe3d323c132`                                                                                             | Downloaded; expected child counts matched in all ten populations; all 377 found in the chr20 VCF header |

The planned trio counts are ACB 20, ASW 13, CLM 35, MXL 32, PEL 35 and
PUR 35 (admixed), and YRI 56, ESN 43, CEU 57 and CHS 51 (controls).
Single-population AF uses AFR for ACB, ASW, YRI and ESN; AMR for CLM,
MXL, PEL and PUR; EUR for CEU; and EAS for CHS.

Chr20 has the complete truth, q and arm evaluation. All 22 children BCFs
are now staged, but truth windows, genome-wide q and arm comparisons
remain pending for the other chromosomes. The population child counts
match the declared expectations: ACB 20, ASW 13, CLM 35, MXL 32, PEL 35,
PUR 35, YRI 56, ESN 43, CEU 57 and CHS 51. The chr20 source VCF was
copied from the supplied local file after its checksum was verified; its
remote metadata is recorded above. bcftools reported
`1.23.1-70-g6dbd8fef`; the source header had 3,202 samples and the
children BCF has 377. The DuckHTS extension was rebuilt from revision
`86f53052` for the arm measurements. The remote chr1 input was streamed
through bcftools; its ETag, Last-Modified value, record count and output
SHA-256 are recorded in [its
receipt](results/roh-ancestry/chr1/input_receipt.tsv).

## Autosome BCF staging

Every children BCF has 377 samples. The per-chromosome source URL,
Last-Modified value, ETag, record count and output SHA-256 are listed in
[the staging
summary](results/roh-ancestry/autosome_staging_summary.csv); raw source
VCFs were streamed and not cached.

| chromosome | source_last_modified          | source_etag              | source_samples | child_samples | child_records | output_sha256                                                    |
|-----------:|:------------------------------|:-------------------------|---------------:|--------------:|--------------:|:-----------------------------------------------------------------|
|          1 | Thu, 29 Oct 2020 17:17:41 GMT | “a63897e0-5b2d27518f740” |           3202 |           377 |       2930166 | 96a9c74a3e3051f74ac1b46b4acf716983ba16713c1e18c000ac0b92eed1fe43 |
|          2 | Thu, 29 Oct 2020 17:19:45 GMT | “ae80101c-5b2d27c7d0e40” |           3202 |           377 |       3139480 | 1ce29c20d828e99420f5f8488c928cc487d309e3aea41f7482dfc75221c46485 |
|          3 | Thu, 29 Oct 2020 17:20:33 GMT | “917b0cb8-5b2d27f597a40” |           3202 |           377 |       2598706 | c9c6c3da847123e437ebb99fd5b10e53c45019ababd4f938ef26e394bb834cc3 |
|          4 | Thu, 29 Oct 2020 17:21:15 GMT | “906bc41b-5b2d281da58c0” |           3202 |           377 |       2575009 | 93127e17f01d424696804b47e996544f524d3113a704e6dcc47f102d3418b067 |
|          5 | Thu, 29 Oct 2020 17:28:14 GMT | “83dd2c4e-5b2d29ad3c780” |           3202 |           377 |       2387322 | 359132848a57db194347242dc78ea9ca0360ea9cd68ba0af61373886c8be9b4d |
|          6 | Thu, 29 Oct 2020 17:28:50 GMT | “81a4595f-5b2d29cf91880” |           3202 |           377 |       2265577 | 32b17e13e0239131c980cbd879ea686d11ea857e69c6f45b59815c73917152db |
|          7 | Thu, 29 Oct 2020 17:29:25 GMT | “7972f934-5b2d29f0f2740” |           3202 |           377 |       2149504 | 7ba8639570bd38b44bee9935fcd4eba3f9f2948b25e3c3646fe1e10389893550 |
|          8 | Thu, 29 Oct 2020 17:51:19 GMT | “705b4bda-5b2d2ed6133c0” |           3202 |           377 |       2036799 | 74c65b9f3dc65d110be11abc0927813260c5410c2609f14a2f731d2a46a935da |
|          9 | Thu, 29 Oct 2020 17:51:45 GMT | “5bba6f2f-5b2d2eeedee40” |           3202 |           377 |       1651537 | bc7ce3156c2561c43c2126e1622d0637c92355222596d4ed25f428e9eb9e501d |
|         10 | Thu, 29 Oct 2020 16:36:26 GMT | “676779b4-5b2d1e1937680” |           3202 |           377 |       1823834 | 05938cee249d0c17acc96d1de0ec0a03d7e84a6a6a4d22a4a0d7ad164082c83a |
|         11 | Thu, 29 Oct 2020 16:36:57 GMT | “64264546-5b2d1e36c7c40” |           3202 |           377 |       1798385 | 63f1a70b3dcfa65cf06ae485600478fa695be0706a87fd071cecdc590f3711a9 |
|         12 | Thu, 29 Oct 2020 16:38:48 GMT | “6231f666-5b2d1ea0a3600” |           3202 |           377 |       1729266 | f071c7dae1c6c76bb29fe06a2bd91d36d101e74522fd45d9dca123b8479075e4 |
|         13 | Thu, 29 Oct 2020 17:06:58 GMT | “4a472041-5b2d24ec59080” |           3202 |           377 |       1303393 | 19c980ecbd9ef1e2dfb1fcae6aa442ccd74502b48f42a0d30b3705f9e22b3afc |
|         14 | Thu, 29 Oct 2020 17:09:12 GMT | “43049970-5b2d256c23e00” |           3202 |           377 |       1188783 | 167df8da0f3b365fc84fb8328ddb4eecf82c430375c824de15d31d0be1f7fcf5 |
|         15 | Thu, 29 Oct 2020 17:09:32 GMT | “3d7a622e-5b2d257f36b00” |           3202 |           377 |       1091118 | 3f14c31a2b5e21dff79b80301f1fb970e0d90570b55d1caaa8968e96dda838a9 |
|         16 | Thu, 29 Oct 2020 17:12:05 GMT | “430843cd-5b2d261120340” |           3202 |           377 |       1212714 | 9f365903619abcdd476bfa0766812ba02e878191f36cbe5a6e2b642dfbdc2940 |
|         17 | Thu, 29 Oct 2020 17:15:50 GMT | “3bdfac68-5b2d26e7b3d80” |           3202 |           377 |       1041119 | 18aa70829e5ee0dbf9086bfffb7265dd9abd1fc9edee5e7c51520b3d3160946b |
|         18 | Thu, 29 Oct 2020 17:16:08 GMT | “39c4fcb9-5b2d26f8de600” |           3202 |           377 |       1028543 | e681322da2b28d9c6d13a308ea1fb99fbe5774c6ef5a3103f86db73e9e915859 |
|         19 | Thu, 29 Oct 2020 17:17:00 GMT | “31a21292-5b2d272a75b00” |           3202 |           377 |        838108 | 892ae4ab2d21b19fb2564edf8b752d1e13b7a55bcf4b31697c722131df5da456 |
|         20 | Thu, 29 Oct 2020 17:17:59 GMT | 2fbb7ed2-5b2d2762b9fc0   |           3202 |           377 |        849143 | a53f83cebd7de6fb4661c14ecd93e5346022ff453392095a085e0b7cda116541 |
|         21 | Thu, 29 Oct 2020 17:18:31 GMT | “1d92bc4c-5b2d27813e7c0” |           3202 |           377 |        526031 | 576c8232cbddf3a6dff5f3daa7a402f5b5a91cfb125446f84bcd04db1858dfe1 |
|         22 | Thu, 29 Oct 2020 17:19:03 GMT | “1efd81b1-5b2d279fc2fc0” |           3202 |           377 |        540617 | 8f6129e3587be10a8ca11dd2904d5e5a8931c2deeef37548a002353261f0be7b |

## Truth-window separation

The complete chr20 truth-window distribution contains 243,165
child-windows (645 windows per child, including the terminal partial
window). Of these, 6,572 have zero heterozygotes and 2,840 have one;
9,412 (3.87%) meet the declared `<=1` criterion, while 233,753 (96.13%)
have at least two. The binned distribution is
[here](results/roh-ancestry/chr20/truth_window_bins.csv), and the exact
per-count distribution is
[here](results/roh-ancestry/chr20/truth_window_distribution.csv). This
separates the low-heterozygosity tail from the remaining windows
sufficiently to retain the preregistered cutoff; no amendment was made.
The retained truth is 102 runs of at least 1 Mb across 92 children; the
run counts and truth-base denominators by population are
[here](results/roh-ancestry/chr20/truth_runs_by_population.csv).

The binned distribution is summarized below; the exact per-count table
is [here](results/roh-ancestry/chr20/truth_window_distribution.csv).

| bin   | child_windows |
|:------|--------------:|
| 0     |          6572 |
| 1     |          2840 |
| 2-4   |          8374 |
| 5-9   |         11801 |
| 10-19 |         13120 |
| 20-49 |         31570 |
| 50+   |        168888 |

The pooled and single-population comparisons began only after this
distribution was reviewed.

## Ancestry proportions

Proportions were estimated from chr20 dosage input only, using the
GRCh38 wide reference and the bigsnpr 1.12.21 correction. Each child
contributed 108,757 matched, non-palindromic loci. Every child has
`status=ok`; the `cor_pred` range and top three groups are reported per
child in [this table](results/roh-ancestry/chr20/q_by_child.csv). This
is a chr20-only estimate, not the final genome-wide q.

The three ROH arms used the same 108,757 sites for every child. The AF
and ancestry reference site keys have zero entries in either direction
of their symmetric difference; the exact comparison is
[here](results/roh-ancestry/chr20/site_set_equivalence.csv). The
`INFO_AF*` values are array-valued Number=A fields; because the input
BCF is biallelic, the AF relation takes its sole value and clamps it to
`[1e-3, 0.999]`. The ancestry arm uses the same clamp. Per-child site
counts are [here](results/roh-ancestry/chr20/site_counts_by_child.csv).

| sample_id | status | cor_pred | used_variants | top1                |     q1 | top2                |     q2 | top3                |     q3 | Population | Superpopulation |
|:----------|:-------|---------:|--------------:|:--------------------|-------:|:--------------------|-------:|:--------------------|-------:|:-----------|:----------------|
| HG01881   | ok     |   0.6665 |        108757 | Africa (West)       | 0.6719 | Africa (South)      | 0.1531 | Europe (North East) | 0.0901 | ACB        | AFR             |
| HG01884   | ok     |   0.7092 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ACB        | AFR             |
| HG01887   | ok     |   0.6719 |        108757 | Africa (West)       | 0.8261 | Europe (North East) | 0.1220 | United Kingdom      | 0.0520 | ACB        | AFR             |
| HG01888   | ok     |   0.6750 |        108757 | Africa (West)       | 0.7409 | Africa (North)      | 0.1261 | Africa (South)      | 0.0614 | ACB        | AFR             |
| HG01891   | ok     |   0.7147 |        108757 | Africa (West)       | 0.9711 | Middle East         | 0.0194 | Philippines         | 0.0094 | ACB        | AFR             |
| HG01895   | ok     |   0.6764 |        108757 | Africa (West)       | 0.8771 | Europe (South West) | 0.1229 | Africa (East)       | 0.0000 | ACB        | AFR             |
| HG01897   | ok     |   0.7084 |        108757 | Africa (West)       | 0.9908 | Asia (East)         | 0.0092 | Africa (East)       | 0.0000 | ACB        | AFR             |
| HG01916   | ok     |   0.6545 |        108757 | Africa (West)       | 0.4087 | Africa (South)      | 0.2464 | United Kingdom      | 0.1977 | ACB        | AFR             |
| HG01959   | ok     |   0.6916 |        108757 | Africa (West)       | 0.9086 | Middle East         | 0.0690 | Japan               | 0.0127 | ACB        | AFR             |
| HG01960   | ok     |   0.6326 |        108757 | Africa (West)       | 0.5565 | Europe (South West) | 0.2835 | United Kingdom      | 0.0975 | ACB        | AFR             |
| HG01987   | ok     |   0.7036 |        108757 | Africa (West)       | 0.9440 | Africa (South)      | 0.0560 | Africa (East)       | 0.0000 | ACB        | AFR             |
| HG02011   | ok     |   0.6559 |        108757 | Africa (West)       | 0.5504 | Africa (South)      | 0.1853 | United Kingdom      | 0.0937 | ACB        | AFR             |
| HG02055   | ok     |   0.6817 |        108757 | Africa (West)       | 0.8710 | Ashkenazi           | 0.0951 | Africa (South)      | 0.0182 | ACB        | AFR             |
| HG02145   | ok     |   0.7101 |        108757 | Africa (West)       | 0.9808 | Africa (South)      | 0.0192 | Africa (East)       | 0.0000 | ACB        | AFR             |
| HG02257   | ok     |   0.6438 |        108757 | Africa (West)       | 0.6058 | Ireland             | 0.1926 | Pakistan            | 0.1159 | ACB        | AFR             |
| HG02258   | ok     |   0.6940 |        108757 | Africa (West)       | 0.8066 | Africa (South)      | 0.1226 | Ireland             | 0.0708 | ACB        | AFR             |
| HG02280   | ok     |   0.6631 |        108757 | Africa (West)       | 0.7630 | Ireland             | 0.1631 | Africa (South)      | 0.0415 | ACB        | AFR             |
| HG02316   | ok     |   0.6920 |        108757 | Africa (West)       | 0.9007 | Ireland             | 0.0796 | United Kingdom      | 0.0181 | ACB        | AFR             |
| HG02321   | ok     |   0.6927 |        108757 | Africa (West)       | 0.7796 | Africa (East)       | 0.1386 | Ireland             | 0.0563 | ACB        | AFR             |
| HG02451   | ok     |   0.7110 |        108757 | Africa (West)       | 0.9746 | Finland             | 0.0231 | Asia (East)         | 0.0022 | ACB        | AFR             |
| NA19702   | ok     |   0.6893 |        108757 | Africa (South)      | 0.4740 | Africa (West)       | 0.3646 | Africa (North)      | 0.1614 | ASW        | AFR             |
| NA19705   | ok     |   0.7154 |        108757 | Africa (West)       | 0.9913 | Africa (South)      | 0.0087 | Africa (East)       | 0.0000 | ASW        | AFR             |
| NA19828   | ok     |   0.6405 |        108757 | Africa (West)       | 0.3319 | Africa (South)      | 0.2274 | Sri Lanka           | 0.1274 | ASW        | AFR             |
| NA19836   | ok     |   0.6720 |        108757 | Africa (West)       | 0.8487 | Finland             | 0.0893 | Europe (North East) | 0.0216 | ASW        | AFR             |
| NA19902   | ok     |   0.6896 |        108757 | Africa (West)       | 0.6270 | Africa (South)      | 0.2822 | Ireland             | 0.0434 | ASW        | AFR             |
| NA19918   | ok     |   0.7045 |        108757 | Africa (West)       | 0.7602 | Africa (South)      | 0.2398 | Africa (East)       | 0.0000 | ASW        | AFR             |
| NA19919   | ok     |   0.6699 |        108757 | Africa (West)       | 0.6600 | Africa (South)      | 0.1938 | Middle East         | 0.0751 | ASW        | AFR             |
| NA19924   | ok     |   0.6403 |        108757 | Africa (West)       | 0.4296 | Scandinavia         | 0.2159 | Africa (South)      | 0.1798 | ASW        | AFR             |
| NA19983   | ok     |   0.6748 |        108757 | Africa (West)       | 0.4777 | Africa (South)      | 0.2844 | Ireland             | 0.2089 | ASW        | AFR             |
| NA20128   | ok     |   0.6754 |        108757 | Africa (West)       | 0.8412 | Ashkenazi           | 0.1443 | Africa (North)      | 0.0106 | ASW        | AFR             |
| NA20129   | ok     |   0.6527 |        108757 | Africa (West)       | 0.5033 | Europe (South West) | 0.4121 | Ashkenazi           | 0.0595 | ASW        | AFR             |
| NA20279   | ok     |   0.6638 |        108757 | Africa (West)       | 0.6908 | Ireland             | 0.1279 | Finland             | 0.0974 | ASW        | AFR             |
| NA20358   | ok     |   0.6755 |        108757 | Africa (West)       | 0.4711 | Africa (South)      | 0.4217 | South America       | 0.0590 | ASW        | AFR             |
| NA06991   | ok     |   0.6646 |        108757 | Ireland             | 0.4924 | Finland             | 0.2503 | Italy               | 0.1275 | CEU        | EUR             |
| NA06995   | ok     |   0.6566 |        108757 | United Kingdom      | 0.3996 | Ireland             | 0.2656 | Europe (South West) | 0.1655 | CEU        | EUR             |
| NA06997   | ok     |   0.6894 |        108757 | Ireland             | 0.5416 | United Kingdom      | 0.2030 | Europe (South West) | 0.1699 | CEU        | EUR             |
| NA07014   | ok     |   0.6427 |        108757 | Scandinavia         | 0.4501 | Ireland             | 0.4224 | Europe (North East) | 0.0497 | CEU        | EUR             |
| NA07019   | ok     |   0.6695 |        108757 | United Kingdom      | 0.5328 | Middle East         | 0.2835 | Ireland             | 0.0998 | CEU        | EUR             |
| NA07029   | ok     |   0.6585 |        108757 | United Kingdom      | 0.4458 | Scandinavia         | 0.3009 | Africa (North)      | 0.1670 | CEU        | EUR             |
| NA07048   | ok     |   0.6809 |        108757 | Scandinavia         | 0.9083 | Europe (North East) | 0.0531 | Philippines         | 0.0386 | CEU        | EUR             |
| NA07348   | ok     |   0.6574 |        108757 | Scandinavia         | 0.8821 | Finland             | 0.0949 | Philippines         | 0.0213 | CEU        | EUR             |
| NA07349   | ok     |   0.6622 |        108757 | Europe (North East) | 0.3470 | Europe (South West) | 0.2226 | Finland             | 0.1682 | CEU        | EUR             |
| NA10830   | ok     |   0.6871 |        108757 | Ireland             | 0.4292 | Finland             | 0.2568 | Europe (North East) | 0.1845 | CEU        | EUR             |
| NA10831   | ok     |   0.6541 |        108757 | Europe (South West) | 0.3909 | Ireland             | 0.3399 | United Kingdom      | 0.2478 | CEU        | EUR             |
| NA10835   | ok     |   0.6810 |        108757 | Scandinavia         | 0.3663 | Europe (North East) | 0.3385 | Ireland             | 0.2953 | CEU        | EUR             |
| NA10836   | ok     |   0.6649 |        108757 | United Kingdom      | 0.6476 | Ireland             | 0.2759 | Europe (South West) | 0.0674 | CEU        | EUR             |
| NA10837   | ok     |   0.6510 |        108757 | Scandinavia         | 0.4560 | Europe (South West) | 0.3208 | United Kingdom      | 0.1007 | CEU        | EUR             |
| NA10838   | ok     |   0.6500 |        108757 | Scandinavia         | 0.4143 | Europe (North East) | 0.1857 | Ireland             | 0.1099 | CEU        | EUR             |
| NA10839   | ok     |   0.6788 |        108757 | Ireland             | 0.7326 | United Kingdom      | 0.1291 | Africa (North)      | 0.0779 | CEU        | EUR             |
| NA10840   | ok     |   0.6783 |        108757 | Europe (North East) | 0.6352 | Europe (South West) | 0.3648 | Africa (East)       | 0.0000 | CEU        | EUR             |
| NA10842   | ok     |   0.6471 |        108757 | Ireland             | 0.4624 | Italy               | 0.2509 | Europe (South East) | 0.1142 | CEU        | EUR             |
| NA10843   | ok     |   0.6567 |        108757 | Scandinavia         | 0.4544 | Europe (North East) | 0.2375 | Pakistan            | 0.1700 | CEU        | EUR             |
| NA10845   | ok     |   0.6558 |        108757 | Ireland             | 0.5233 | Europe (North East) | 0.1869 | United Kingdom      | 0.0983 | CEU        | EUR             |
| NA10846   | ok     |   0.6577 |        108757 | United Kingdom      | 0.5444 | Europe (South West) | 0.2208 | Finland             | 0.1626 | CEU        | EUR             |
| NA10847   | ok     |   0.6337 |        108757 | United Kingdom      | 0.5961 | Italy               | 0.2927 | Scandinavia         | 0.0607 | CEU        | EUR             |
| NA10851   | ok     |   0.6861 |        108757 | Scandinavia         | 0.8879 | Ireland             | 0.1121 | Africa (East)       | 0.0000 | CEU        | EUR             |
| NA10852   | ok     |   0.6498 |        108757 | Scandinavia         | 0.5878 | United Kingdom      | 0.2658 | Asia (East)         | 0.0628 | CEU        | EUR             |
| NA10854   | ok     |   0.6895 |        108757 | United Kingdom      | 0.4021 | Finland             | 0.3123 | Europe (South West) | 0.1581 | CEU        | EUR             |
| NA10855   | ok     |   0.6610 |        108757 | United Kingdom      | 0.7579 | Ashkenazi           | 0.1381 | Africa (South)      | 0.0473 | CEU        | EUR             |
| NA10856   | ok     |   0.6413 |        108757 | Ireland             | 0.4521 | Scandinavia         | 0.2867 | Europe (South West) | 0.1201 | CEU        | EUR             |
| NA10857   | ok     |   0.6846 |        108757 | Europe (South West) | 0.5401 | Scandinavia         | 0.2765 | Finland             | 0.0961 | CEU        | EUR             |
| NA10859   | ok     |   0.6508 |        108757 | Ireland             | 0.8004 | Ashkenazi           | 0.1507 | South America       | 0.0281 | CEU        | EUR             |
| NA10860   | ok     |   0.6754 |        108757 | Scandinavia         | 0.6334 | Europe (North East) | 0.2497 | Finland             | 0.0926 | CEU        | EUR             |
| NA10861   | ok     |   0.6662 |        108757 | Ireland             | 0.5553 | United Kingdom      | 0.2555 | Scandinavia         | 0.1488 | CEU        | EUR             |
| NA10863   | ok     |   0.6711 |        108757 | Ireland             | 0.9723 | Japan               | 0.0216 | Africa (South)      | 0.0061 | CEU        | EUR             |
| NA10864   | ok     |   0.6892 |        108757 | Ireland             | 0.4370 | Finland             | 0.2481 | Italy               | 0.2136 | CEU        | EUR             |
| NA10865   | ok     |   0.6843 |        108757 | United Kingdom      | 0.5033 | Europe (South West) | 0.2615 | Finland             | 0.1115 | CEU        | EUR             |
| NA12329   | ok     |   0.6910 |        108757 | Ireland             | 0.7380 | United Kingdom      | 0.2019 | Middle East         | 0.0369 | CEU        | EUR             |
| NA12335   | ok     |   0.6743 |        108757 | Europe (South West) | 0.5385 | Ireland             | 0.2638 | Finland             | 0.1896 | CEU        | EUR             |
| NA12336   | ok     |   0.6655 |        108757 | Europe (South West) | 0.4362 | United Kingdom      | 0.4162 | Europe (North East) | 0.0911 | CEU        | EUR             |
| NA12344   | ok     |   0.6484 |        108757 | Europe (South West) | 0.5153 | Italy               | 0.4317 | Sri Lanka           | 0.0425 | CEU        | EUR             |
| NA12376   | ok     |   0.6591 |        108757 | Ireland             | 0.6537 | Europe (South West) | 0.3143 | Philippines         | 0.0320 | CEU        | EUR             |
| NA12386   | ok     |   0.6397 |        108757 | Italy               | 0.5016 | Europe (South West) | 0.3337 | United Kingdom      | 0.0852 | CEU        | EUR             |
| NA12485   | ok     |   0.6739 |        108757 | Ireland             | 0.5270 | Europe (South East) | 0.3377 | Europe (South West) | 0.1353 | CEU        | EUR             |
| NA12707   | ok     |   0.6752 |        108757 | Europe (South West) | 0.6394 | Ireland             | 0.3024 | Asia (East)         | 0.0405 | CEU        | EUR             |
| NA12739   | ok     |   0.6558 |        108757 | United Kingdom      | 0.8967 | Ireland             | 0.0681 | Africa (West)       | 0.0213 | CEU        | EUR             |
| NA12740   | ok     |   0.6652 |        108757 | Ireland             | 0.6885 | Pakistan            | 0.1587 | Finland             | 0.0983 | CEU        | EUR             |
| NA12752   | ok     |   0.6595 |        108757 | Ireland             | 0.6146 | Finland             | 0.2700 | Middle East         | 0.1025 | CEU        | EUR             |
| NA12753   | ok     |   0.6610 |        108757 | Europe (South West) | 0.9117 | South America       | 0.0511 | Ireland             | 0.0250 | CEU        | EUR             |
| NA12766   | ok     |   0.6586 |        108757 | Europe (South West) | 0.4242 | Scandinavia         | 0.3695 | South America       | 0.0787 | CEU        | EUR             |
| NA12767   | ok     |   0.6734 |        108757 | United Kingdom      | 0.3157 | Italy               | 0.3121 | Europe (South West) | 0.1717 | CEU        | EUR             |
| NA12801   | ok     |   0.6821 |        108757 | Europe (South East) | 0.4875 | Europe (South West) | 0.4079 | Ireland             | 0.0788 | CEU        | EUR             |
| NA12802   | ok     |   0.6737 |        108757 | Scandinavia         | 0.5939 | Europe (North East) | 0.2075 | Ireland             | 0.1883 | CEU        | EUR             |
| NA12817   | ok     |   0.6813 |        108757 | Scandinavia         | 0.8159 | Ashkenazi           | 0.1601 | Asia (East)         | 0.0240 | CEU        | EUR             |
| NA12818   | ok     |   0.6630 |        108757 | Ireland             | 0.4911 | Europe (North East) | 0.1901 | Scandinavia         | 0.1812 | CEU        | EUR             |
| NA12832   | ok     |   0.6562 |        108757 | Ireland             | 0.6803 | Europe (South West) | 0.2656 | Finland             | 0.0541 | CEU        | EUR             |
| NA12864   | ok     |   0.6514 |        108757 | Europe (South West) | 0.7864 | Ireland             | 0.1920 | Finland             | 0.0216 | CEU        | EUR             |
| NA12865   | ok     |   0.6668 |        108757 | Europe (South East) | 0.6190 | Ireland             | 0.2435 | Europe (North East) | 0.1113 | CEU        | EUR             |
| NA12877   | ok     |   0.6971 |        108757 | Scandinavia         | 0.4812 | United Kingdom      | 0.4442 | Europe (South West) | 0.0746 | CEU        | EUR             |
| NA12878   | ok     |   0.6664 |        108757 | Ireland             | 0.8512 | Finland             | 0.0941 | Sri Lanka           | 0.0441 | CEU        | EUR             |
| HG00405   | ok     |   0.6756 |        108757 | Asia (East)         | 0.6358 | Japan               | 0.3390 | Philippines         | 0.0244 | CHS        | EAS             |
| HG00408   | ok     |   0.7036 |        108757 | Asia (East)         | 0.8757 | Japan               | 0.0957 | Africa (East)       | 0.0286 | CHS        | EAS             |
| HG00420   | ok     |   0.6858 |        108757 | Asia (East)         | 0.8753 | Philippines         | 0.0777 | South America       | 0.0470 | CHS        | EAS             |
| HG00423   | ok     |   0.7242 |        108757 | Asia (East)         | 0.6303 | Philippines         | 0.2004 | Japan               | 0.1693 | CHS        | EAS             |
| HG00429   | ok     |   0.7113 |        108757 | Asia (East)         | 0.8821 | Japan               | 0.0986 | Africa (West)       | 0.0124 | CHS        | EAS             |
| HG00438   | ok     |   0.7195 |        108757 | Asia (East)         | 0.9551 | Japan               | 0.0449 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00444   | ok     |   0.7025 |        108757 | Asia (East)         | 0.7298 | Japan               | 0.1868 | Sri Lanka           | 0.0811 | CHS        | EAS             |
| HG00447   | ok     |   0.7121 |        108757 | Asia (East)         | 0.8974 | Japan               | 0.0916 | Africa (East)       | 0.0110 | CHS        | EAS             |
| HG00450   | ok     |   0.7160 |        108757 | Asia (East)         | 0.9194 | Japan               | 0.0264 | Ashkenazi           | 0.0259 | CHS        | EAS             |
| HG00453   | ok     |   0.6927 |        108757 | Asia (East)         | 0.7767 | Japan               | 0.2233 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00459   | ok     |   0.7010 |        108757 | Asia (East)         | 0.8908 | Japan               | 0.1092 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00465   | ok     |   0.7019 |        108757 | Asia (East)         | 0.7419 | Japan               | 0.2144 | South America       | 0.0437 | CHS        | EAS             |
| HG00474   | ok     |   0.7164 |        108757 | Asia (East)         | 0.8789 | Japan               | 0.1211 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00477   | ok     |   0.7055 |        108757 | Asia (East)         | 0.8559 | Japan               | 0.1224 | Finland             | 0.0217 | CHS        | EAS             |
| HG00480   | ok     |   0.7077 |        108757 | Asia (East)         | 0.7751 | Japan               | 0.2249 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00502   | ok     |   0.7221 |        108757 | Asia (East)         | 0.7903 | Japan               | 0.2097 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00514   | ok     |   0.7028 |        108757 | Asia (East)         | 0.6607 | Japan               | 0.2783 | Europe (North East) | 0.0610 | CHS        | EAS             |
| HG00526   | ok     |   0.7174 |        108757 | Asia (East)         | 0.9656 | Japan               | 0.0344 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00532   | ok     |   0.7289 |        108757 | Asia (East)         | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | CHS        | EAS             |
| HG00535   | ok     |   0.6858 |        108757 | Asia (East)         | 0.8653 | Japan               | 0.0861 | Finland             | 0.0486 | CHS        | EAS             |
| HG00538   | ok     |   0.6863 |        108757 | Asia (East)         | 0.7714 | Japan               | 0.1724 | Sri Lanka           | 0.0382 | CHS        | EAS             |
| HG00544   | ok     |   0.6762 |        108757 | Asia (East)         | 0.8196 | Japan               | 0.1396 | Europe (North East) | 0.0210 | CHS        | EAS             |
| HG00558   | ok     |   0.7111 |        108757 | Asia (East)         | 0.8476 | Japan               | 0.1524 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00561   | ok     |   0.7012 |        108757 | Asia (East)         | 0.7754 | Japan               | 0.1939 | Africa (East)       | 0.0307 | CHS        | EAS             |
| HG00567   | ok     |   0.7048 |        108757 | Asia (East)         | 0.7622 | Japan               | 0.2378 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00579   | ok     |   0.6792 |        108757 | Asia (East)         | 0.8178 | Japan               | 0.0977 | Sri Lanka           | 0.0845 | CHS        | EAS             |
| HG00582   | ok     |   0.6994 |        108757 | Asia (East)         | 0.7985 | Japan               | 0.1399 | Sri Lanka           | 0.0333 | CHS        | EAS             |
| HG00585   | ok     |   0.6919 |        108757 | Asia (East)         | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | CHS        | EAS             |
| HG00591   | ok     |   0.7023 |        108757 | Asia (East)         | 0.6015 | Japan               | 0.3016 | Ashkenazi           | 0.0330 | CHS        | EAS             |
| HG00594   | ok     |   0.7116 |        108757 | Asia (East)         | 0.7219 | Japan               | 0.2516 | Bangladesh          | 0.0265 | CHS        | EAS             |
| HG00597   | ok     |   0.7148 |        108757 | Asia (East)         | 0.8793 | Japan               | 0.1207 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00609   | ok     |   0.7337 |        108757 | Asia (East)         | 0.9353 | South America       | 0.0410 | Japan               | 0.0238 | CHS        | EAS             |
| HG00612   | ok     |   0.7389 |        108757 | Asia (East)         | 0.8389 | Japan               | 0.1611 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00615   | ok     |   0.6961 |        108757 | Asia (East)         | 0.9741 | Ashkenazi           | 0.0154 | Japan               | 0.0071 | CHS        | EAS             |
| HG00621   | ok     |   0.7064 |        108757 | Asia (East)         | 0.7556 | Japan               | 0.2332 | Italy               | 0.0112 | CHS        | EAS             |
| HG00627   | ok     |   0.7169 |        108757 | Asia (East)         | 0.6780 | Philippines         | 0.2718 | Sri Lanka           | 0.0501 | CHS        | EAS             |
| HG00630   | ok     |   0.7089 |        108757 | Asia (East)         | 0.9284 | Japan               | 0.0716 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00636   | ok     |   0.7150 |        108757 | Asia (East)         | 0.7946 | Japan               | 0.2054 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00652   | ok     |   0.7125 |        108757 | Asia (East)         | 0.9852 | Sri Lanka           | 0.0148 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00655   | ok     |   0.7241 |        108757 | Asia (East)         | 0.8205 | Japan               | 0.1311 | Philippines         | 0.0484 | CHS        | EAS             |
| HG00658   | ok     |   0.7136 |        108757 | Asia (East)         | 0.9402 | Sri Lanka           | 0.0470 | Finland             | 0.0098 | CHS        | EAS             |
| HG00664   | ok     |   0.7099 |        108757 | Asia (East)         | 0.9898 | South America       | 0.0102 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00673   | ok     |   0.7139 |        108757 | Asia (East)         | 0.9591 | Japan               | 0.0409 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00685   | ok     |   0.7010 |        108757 | Asia (East)         | 0.9429 | Sri Lanka           | 0.0571 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00691   | ok     |   0.6958 |        108757 | Asia (East)         | 0.9839 | Ashkenazi           | 0.0083 | Africa (South)      | 0.0078 | CHS        | EAS             |
| HG00694   | ok     |   0.7299 |        108757 | Asia (East)         | 0.8852 | Philippines         | 0.1083 | Ireland             | 0.0065 | CHS        | EAS             |
| HG00700   | ok     |   0.7242 |        108757 | Asia (East)         | 0.8991 | Japan               | 0.1009 | Africa (East)       | 0.0000 | CHS        | EAS             |
| HG00702   | ok     |   0.7198 |        108757 | Asia (East)         | 0.9026 | Japan               | 0.0380 | Philippines         | 0.0308 | CHS        | EAS             |
| HG00703   | ok     |   0.7097 |        108757 | Asia (East)         | 0.9426 | Philippines         | 0.0363 | Pakistan            | 0.0211 | CHS        | EAS             |
| HG00706   | ok     |   0.6786 |        108757 | Asia (East)         | 0.8826 | Japan               | 0.1108 | Africa (West)       | 0.0066 | CHS        | EAS             |
| HG00709   | ok     |   0.6922 |        108757 | Asia (East)         | 0.9398 | Japan               | 0.0515 | South America       | 0.0087 | CHS        | EAS             |
| HG01114   | ok     |   0.6576 |        108757 | Europe (South West) | 0.5189 | Europe (South East) | 0.3021 | Sri Lanka           | 0.0803 | CLM        | AMR             |
| HG01126   | ok     |   0.6301 |        108757 | South America       | 0.9506 | Philippines         | 0.0472 | Africa (South)      | 0.0022 | CLM        | AMR             |
| HG01135   | ok     |   0.6518 |        108757 | Europe (South West) | 0.7012 | South America       | 0.1953 | Europe (North East) | 0.0900 | CLM        | AMR             |
| HG01138   | ok     |   0.6641 |        108757 | Italy               | 0.7266 | South America       | 0.1517 | Europe (South East) | 0.0670 | CLM        | AMR             |
| HG01141   | ok     |   0.6375 |        108757 | Europe (South West) | 0.4728 | Africa (West)       | 0.2168 | South America       | 0.2161 | CLM        | AMR             |
| HG01150   | ok     |   0.6358 |        108757 | Europe (South West) | 0.7103 | South America       | 0.1320 | Asia (East)         | 0.0555 | CLM        | AMR             |
| HG01252   | ok     |   0.6882 |        108757 | South America       | 0.9857 | Europe (South West) | 0.0143 | Africa (East)       | 0.0000 | CLM        | AMR             |
| HG01255   | ok     |   0.6246 |        108757 | South America       | 0.6377 | Europe (South West) | 0.2146 | Italy               | 0.0560 | CLM        | AMR             |
| HG01258   | ok     |   0.6664 |        108757 | Europe (South West) | 0.4549 | Europe (North East) | 0.3272 | South America       | 0.1473 | CLM        | AMR             |
| HG01261   | ok     |   0.6222 |        108757 | South America       | 0.4875 | Africa (West)       | 0.1572 | Africa (North)      | 0.1214 | CLM        | AMR             |
| HG01273   | ok     |   0.6395 |        108757 | South America       | 0.4652 | Italy               | 0.2117 | Ashkenazi           | 0.1079 | CLM        | AMR             |
| HG01276   | ok     |   0.6402 |        108757 | South America       | 0.7541 | Scandinavia         | 0.1036 | Ashkenazi           | 0.0837 | CLM        | AMR             |
| HG01279   | ok     |   0.6580 |        108757 | South America       | 0.6508 | Europe (South West) | 0.2226 | Ashkenazi           | 0.0816 | CLM        | AMR             |
| HG01343   | ok     |   0.6957 |        108757 | South America       | 0.9711 | Japan               | 0.0289 | Africa (East)       | 0.0000 | CLM        | AMR             |
| HG01346   | ok     |   0.6714 |        108757 | South America       | 0.4208 | Europe (South West) | 0.2998 | Finland             | 0.1075 | CLM        | AMR             |
| HG01349   | ok     |   0.6687 |        108757 | South America       | 0.9117 | Japan               | 0.0728 | Finland             | 0.0155 | CLM        | AMR             |
| HG01352   | ok     |   0.7021 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | CLM        | AMR             |
| HG01355   | ok     |   0.6398 |        108757 | South America       | 0.5695 | Europe (South West) | 0.3391 | Italy               | 0.0817 | CLM        | AMR             |
| HG01358   | ok     |   0.6498 |        108757 | South America       | 0.6398 | Europe (South West) | 0.2470 | Africa (West)       | 0.0670 | CLM        | AMR             |
| HG01361   | ok     |   0.6552 |        108757 | Europe (South West) | 0.3183 | Africa (West)       | 0.1923 | United Kingdom      | 0.1723 | CLM        | AMR             |
| HG01367   | ok     |   0.6449 |        108757 | South America       | 0.7270 | Africa (West)       | 0.1775 | Africa (East)       | 0.0471 | CLM        | AMR             |
| HG01376   | ok     |   0.6676 |        108757 | Europe (North East) | 0.3398 | Europe (South West) | 0.2950 | South America       | 0.2857 | CLM        | AMR             |
| HG01379   | ok     |   0.6692 |        108757 | South America       | 0.4766 | Europe (South West) | 0.3615 | Ireland             | 0.1378 | CLM        | AMR             |
| HG01385   | ok     |   0.6510 |        108757 | South America       | 0.5708 | Europe (South West) | 0.2377 | Finland             | 0.1077 | CLM        | AMR             |
| HG01391   | ok     |   0.6539 |        108757 | South America       | 0.3682 | Italy               | 0.3585 | Africa (West)       | 0.1454 | CLM        | AMR             |
| HG01433   | ok     |   0.6579 |        108757 | South America       | 0.5059 | Europe (South West) | 0.4874 | Philippines         | 0.0067 | CLM        | AMR             |
| HG01439   | ok     |   0.6703 |        108757 | Ireland             | 0.4568 | Europe (South West) | 0.2105 | South America       | 0.1339 | CLM        | AMR             |
| HG01457   | ok     |   0.6586 |        108757 | South America       | 0.3421 | Italy               | 0.3215 | Europe (South West) | 0.3165 | CLM        | AMR             |
| HG01463   | ok     |   0.6460 |        108757 | South America       | 0.5370 | Scandinavia         | 0.2565 | Africa (South)      | 0.0965 | CLM        | AMR             |
| HG01466   | ok     |   0.6434 |        108757 | Europe (South West) | 0.4869 | Italy               | 0.2047 | South America       | 0.1667 | CLM        | AMR             |
| HG01490   | ok     |   0.6406 |        108757 | South America       | 0.3055 | Europe (South West) | 0.2700 | Scandinavia         | 0.2411 | CLM        | AMR             |
| HG01493   | ok     |   0.6514 |        108757 | Europe (South West) | 0.6034 | South America       | 0.3754 | Finland             | 0.0213 | CLM        | AMR             |
| HG01496   | ok     |   0.6749 |        108757 | South America       | 0.6873 | Middle East         | 0.1054 | Europe (North East) | 0.0977 | CLM        | AMR             |
| HG01499   | ok     |   0.6487 |        108757 | South America       | 0.5003 | Europe (South West) | 0.2692 | Scandinavia         | 0.1968 | CLM        | AMR             |
| HG01552   | ok     |   0.6324 |        108757 | Europe (South West) | 0.4521 | South America       | 0.3471 | Africa (South)      | 0.1984 | CLM        | AMR             |
| HG02924   | ok     |   0.7182 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ESN        | AFR             |
| HG02945   | ok     |   0.7223 |        108757 | Africa (West)       | 0.9687 | Africa (South)      | 0.0280 | South America       | 0.0033 | ESN        | AFR             |
| HG02948   | ok     |   0.7151 |        108757 | Africa (West)       | 0.8464 | Africa (South)      | 0.1536 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG02954   | ok     |   0.7063 |        108757 | Africa (West)       | 0.8480 | Africa (South)      | 0.1520 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG02966   | ok     |   0.7194 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ESN        | AFR             |
| HG02972   | ok     |   0.7097 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ESN        | AFR             |
| HG02975   | ok     |   0.7097 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ESN        | AFR             |
| HG02978   | ok     |   0.6994 |        108757 | Africa (West)       | 0.9839 | Philippines         | 0.0146 | Japan               | 0.0016 | ESN        | AFR             |
| HG02980   | ok     |   0.7031 |        108757 | Africa (West)       | 0.7221 | Africa (South)      | 0.2779 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG03101   | ok     |   0.7023 |        108757 | Africa (West)       | 0.9839 | Africa (South)      | 0.0161 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG03110   | ok     |   0.7092 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ESN        | AFR             |
| HG03113   | ok     |   0.7110 |        108757 | Africa (West)       | 0.8737 | Africa (South)      | 0.1263 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG03116   | ok     |   0.7045 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ESN        | AFR             |
| HG03119   | ok     |   0.6993 |        108757 | Africa (West)       | 0.9590 | Africa (South)      | 0.0373 | Philippines         | 0.0021 | ESN        | AFR             |
| HG03122   | ok     |   0.7169 |        108757 | Africa (West)       | 0.9982 | Philippines         | 0.0018 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG03125   | ok     |   0.7167 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ESN        | AFR             |
| HG03128   | ok     |   0.7093 |        108757 | Africa (West)       | 0.8702 | Africa (South)      | 0.1298 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG03131   | ok     |   0.7014 |        108757 | Africa (West)       | 0.8759 | Africa (South)      | 0.1139 | Africa (East)       | 0.0101 | ESN        | AFR             |
| HG03134   | ok     |   0.7125 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ESN        | AFR             |
| HG03137   | ok     |   0.7192 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ESN        | AFR             |
| HG03161   | ok     |   0.6991 |        108757 | Africa (West)       | 0.9659 | Europe (South West) | 0.0234 | Ireland             | 0.0106 | ESN        | AFR             |
| HG03164   | ok     |   0.7073 |        108757 | Africa (West)       | 0.7492 | Africa (South)      | 0.2444 | Philippines         | 0.0064 | ESN        | AFR             |
| HG03170   | ok     |   0.7005 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ESN        | AFR             |
| HG03191   | ok     |   0.7107 |        108757 | Africa (West)       | 0.6306 | Africa (South)      | 0.3694 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG03197   | ok     |   0.7108 |        108757 | Africa (West)       | 0.8390 | Africa (South)      | 0.1610 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG03200   | ok     |   0.7019 |        108757 | Africa (West)       | 0.8454 | Africa (South)      | 0.1045 | Philippines         | 0.0280 | ESN        | AFR             |
| HG03269   | ok     |   0.7053 |        108757 | Africa (West)       | 0.6693 | Africa (South)      | 0.3110 | Philippines         | 0.0197 | ESN        | AFR             |
| HG03272   | ok     |   0.7163 |        108757 | Africa (West)       | 0.9152 | Africa (South)      | 0.0848 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG03296   | ok     |   0.7021 |        108757 | Africa (West)       | 0.9201 | Africa (South)      | 0.0799 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG03299   | ok     |   0.7183 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ESN        | AFR             |
| HG03302   | ok     |   0.7146 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ESN        | AFR             |
| HG03305   | ok     |   0.7132 |        108757 | Africa (West)       | 0.9890 | Africa (South)      | 0.0110 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG03308   | ok     |   0.7060 |        108757 | Africa (West)       | 0.9993 | Pakistan            | 0.0007 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG03310   | ok     |   0.6875 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ESN        | AFR             |
| HG03314   | ok     |   0.6985 |        108757 | Africa (West)       | 0.9992 | Japan               | 0.0008 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG03344   | ok     |   0.7128 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ESN        | AFR             |
| HG03350   | ok     |   0.7153 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ESN        | AFR             |
| HG03368   | ok     |   0.7167 |        108757 | Africa (West)       | 0.8155 | Africa (South)      | 0.1845 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG03371   | ok     |   0.7228 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | ESN        | AFR             |
| HG03374   | ok     |   0.7166 |        108757 | Africa (West)       | 0.9576 | Africa (South)      | 0.0424 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG03516   | ok     |   0.6994 |        108757 | Africa (West)       | 0.9698 | Africa (South)      | 0.0302 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG03519   | ok     |   0.6950 |        108757 | Africa (West)       | 0.9913 | Africa (South)      | 0.0087 | Africa (East)       | 0.0000 | ESN        | AFR             |
| HG03522   | ok     |   0.7013 |        108757 | Africa (West)       | 0.9149 | Africa (South)      | 0.0692 | Philippines         | 0.0159 | ESN        | AFR             |
| NA19650   | ok     |   0.6629 |        108757 | United Kingdom      | 0.2864 | Ashkenazi           | 0.2564 | Ireland             | 0.1667 | MXL        | AMR             |
| NA19653   | ok     |   0.6559 |        108757 | South America       | 0.4704 | Italy               | 0.3498 | Middle East         | 0.1125 | MXL        | AMR             |
| NA19656   | ok     |   0.6476 |        108757 | South America       | 0.6772 | Africa (West)       | 0.1092 | Italy               | 0.0811 | MXL        | AMR             |
| NA19659   | ok     |   0.6630 |        108757 | South America       | 0.4043 | Europe (South West) | 0.2284 | Europe (South East) | 0.1853 | MXL        | AMR             |
| NA19662   | ok     |   0.6348 |        108757 | Italy               | 0.3652 | South America       | 0.3444 | Ashkenazi           | 0.1173 | MXL        | AMR             |
| NA19665   | ok     |   0.6550 |        108757 | South America       | 0.7919 | Africa (West)       | 0.1545 | Asia (East)         | 0.0408 | MXL        | AMR             |
| NA19671   | ok     |   0.6474 |        108757 | South America       | 0.6523 | Italy               | 0.1549 | Finland             | 0.0754 | MXL        | AMR             |
| NA19675   | ok     |   0.6458 |        108757 | South America       | 0.4463 | Middle East         | 0.1666 | Africa (West)       | 0.1334 | MXL        | AMR             |
| NA19677   | ok     |   0.6394 |        108757 | South America       | 0.8994 | Africa (West)       | 0.1006 | Africa (East)       | 0.0000 | MXL        | AMR             |
| NA19680   | ok     |   0.6205 |        108757 | United Kingdom      | 0.3389 | South America       | 0.2805 | Africa (West)       | 0.1731 | MXL        | AMR             |
| NA19683   | ok     |   0.6809 |        108757 | South America       | 0.9565 | Finland             | 0.0363 | Philippines         | 0.0071 | MXL        | AMR             |
| NA19685   | ok     |   0.6359 |        108757 | Europe (South West) | 0.3632 | South America       | 0.2773 | Ashkenazi           | 0.0995 | MXL        | AMR             |
| NA19686   | ok     |   0.6413 |        108757 | South America       | 0.5399 | Europe (South West) | 0.3143 | Sri Lanka           | 0.0611 | MXL        | AMR             |
| NA19718   | ok     |   0.6704 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | MXL        | AMR             |
| NA19721   | ok     |   0.6648 |        108757 | South America       | 0.8791 | Japan               | 0.0928 | Africa (West)       | 0.0281 | MXL        | AMR             |
| NA19724   | ok     |   0.6311 |        108757 | South America       | 0.6179 | Africa (North)      | 0.2077 | Europe (South West) | 0.0657 | MXL        | AMR             |
| NA19727   | ok     |   0.6671 |        108757 | South America       | 0.9340 | Europe (North East) | 0.0410 | Ireland             | 0.0251 | MXL        | AMR             |
| NA19730   | ok     |   0.7095 |        108757 | South America       | 0.8229 | Japan               | 0.1771 | Africa (East)       | 0.0000 | MXL        | AMR             |
| NA19733   | ok     |   0.7047 |        108757 | South America       | 0.9228 | Japan               | 0.0772 | Africa (East)       | 0.0000 | MXL        | AMR             |
| NA19748   | ok     |   0.6626 |        108757 | South America       | 0.6922 | Europe (South East) | 0.1662 | Ireland             | 0.0753 | MXL        | AMR             |
| NA19751   | ok     |   0.6381 |        108757 | South America       | 0.4060 | Europe (South West) | 0.2443 | Finland             | 0.1435 | MXL        | AMR             |
| NA19757   | ok     |   0.6648 |        108757 | South America       | 0.8854 | Japan               | 0.1146 | Africa (East)       | 0.0000 | MXL        | AMR             |
| NA19760   | ok     |   0.6606 |        108757 | South America       | 0.9634 | Japan               | 0.0243 | Asia (East)         | 0.0123 | MXL        | AMR             |
| NA19763   | ok     |   0.6717 |        108757 | South America       | 0.9816 | Asia (East)         | 0.0184 | Africa (East)       | 0.0000 | MXL        | AMR             |
| NA19772   | ok     |   0.6542 |        108757 | South America       | 0.3148 | Ireland             | 0.2459 | Italy               | 0.1755 | MXL        | AMR             |
| NA19775   | ok     |   0.6723 |        108757 | South America       | 0.7403 | Italy               | 0.2597 | Africa (East)       | 0.0000 | MXL        | AMR             |
| NA19778   | ok     |   0.6785 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | MXL        | AMR             |
| NA19781   | ok     |   0.6765 |        108757 | South America       | 0.9417 | Japan               | 0.0583 | Africa (East)       | 0.0000 | MXL        | AMR             |
| NA19784   | ok     |   0.7008 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | MXL        | AMR             |
| NA19787   | ok     |   0.6710 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | MXL        | AMR             |
| NA19790   | ok     |   0.7076 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | MXL        | AMR             |
| NA19796   | ok     |   0.6616 |        108757 | South America       | 0.6678 | Italy               | 0.1922 | Finland             | 0.1137 | MXL        | AMR             |
| HG01567   | ok     |   0.6723 |        108757 | South America       | 0.8698 | Ireland             | 0.0832 | Middle East         | 0.0470 | PEL        | AMR             |
| HG01573   | ok     |   0.6962 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01579   | ok     |   0.6747 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01898   | ok     |   0.6866 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01919   | ok     |   0.7124 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01922   | ok     |   0.7023 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01925   | ok     |   0.7035 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01928   | ok     |   0.7115 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01934   | ok     |   0.7002 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01937   | ok     |   0.7087 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01940   | ok     |   0.7186 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01943   | ok     |   0.7078 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01946   | ok     |   0.7070 |        108757 | South America       | 0.7794 | Japan               | 0.1833 | Asia (East)         | 0.0374 | PEL        | AMR             |
| HG01949   | ok     |   0.6676 |        108757 | South America       | 0.9179 | Finland             | 0.0767 | Europe (North East) | 0.0055 | PEL        | AMR             |
| HG01952   | ok     |   0.7048 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01955   | ok     |   0.7085 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01969   | ok     |   0.7100 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01972   | ok     |   0.6782 |        108757 | South America       | 0.9204 | Japan               | 0.0796 | Africa (East)       | 0.0000 | PEL        | AMR             |
| HG01975   | ok     |   0.6892 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01978   | ok     |   0.7014 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01981   | ok     |   0.6922 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01984   | ok     |   0.6840 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01993   | ok     |   0.6888 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG01998   | ok     |   0.6698 |        108757 | South America       | 0.9191 | Philippines         | 0.0644 | Ashkenazi           | 0.0092 | PEL        | AMR             |
| HG02004   | ok     |   0.7137 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG02091   | ok     |   0.7040 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG02106   | ok     |   0.6970 |        108757 | South America       | 0.9808 | Japan               | 0.0192 | Africa (East)       | 0.0000 | PEL        | AMR             |
| HG02148   | ok     |   0.6977 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG02261   | ok     |   0.6980 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG02273   | ok     |   0.7154 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG02279   | ok     |   0.7013 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG02287   | ok     |   0.7007 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG02293   | ok     |   0.7057 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG02300   | ok     |   0.6945 |        108757 | South America       | 0.9249 | Japan               | 0.0751 | Africa (East)       | 0.0000 | PEL        | AMR             |
| HG02303   | ok     |   0.7070 |        108757 | South America       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | PEL        | AMR             |
| HG00552   | ok     |   0.6585 |        108757 | Europe (South West) | 0.3551 | Ireland             | 0.2098 | Africa (North)      | 0.1261 | PUR        | AMR             |
| HG00555   | ok     |   0.6235 |        108757 | South America       | 0.4795 | Africa (West)       | 0.2162 | Europe (South West) | 0.1903 | PUR        | AMR             |
| HG00639   | ok     |   0.6623 |        108757 | Europe (South West) | 0.4297 | Italy               | 0.2288 | Ashkenazi           | 0.1053 | PUR        | AMR             |
| HG00642   | ok     |   0.6317 |        108757 | Europe (South West) | 0.3514 | Europe (North East) | 0.2891 | South America       | 0.1379 | PUR        | AMR             |
| HG00733   | ok     |   0.6688 |        108757 | Italy               | 0.4782 | South America       | 0.1924 | Europe (South West) | 0.1477 | PUR        | AMR             |
| HG00735   | ok     |   0.6560 |        108757 | Europe (South West) | 0.5078 | South America       | 0.4473 | Japan               | 0.0449 | PUR        | AMR             |
| HG00738   | ok     |   0.6297 |        108757 | Italy               | 0.4219 | Ireland             | 0.3159 | Africa (West)       | 0.1995 | PUR        | AMR             |
| HG00741   | ok     |   0.6404 |        108757 | Europe (South West) | 0.3887 | Africa (South)      | 0.3089 | South America       | 0.1159 | PUR        | AMR             |
| HG01050   | ok     |   0.6406 |        108757 | South America       | 0.4780 | Africa (North)      | 0.1692 | Middle East         | 0.1660 | PUR        | AMR             |
| HG01053   | ok     |   0.6460 |        108757 | South America       | 0.3709 | Africa (West)       | 0.2940 | Ashkenazi           | 0.1937 | PUR        | AMR             |
| HG01056   | ok     |   0.6476 |        108757 | Italy               | 0.4856 | South America       | 0.2644 | Africa (North)      | 0.1307 | PUR        | AMR             |
| HG01062   | ok     |   0.6492 |        108757 | South America       | 0.4362 | Europe (South West) | 0.3105 | Finland             | 0.2035 | PUR        | AMR             |
| HG01068   | ok     |   0.6322 |        108757 | Europe (South West) | 0.4928 | Europe (South East) | 0.1627 | Africa (North)      | 0.1001 | PUR        | AMR             |
| HG01071   | ok     |   0.6495 |        108757 | Italy               | 0.4537 | Africa (North)      | 0.3090 | Finland             | 0.1715 | PUR        | AMR             |
| HG01074   | ok     |   0.6587 |        108757 | South America       | 0.3859 | Italy               | 0.2634 | Europe (South West) | 0.2369 | PUR        | AMR             |
| HG01081   | ok     |   0.6357 |        108757 | Italy               | 0.3869 | Europe (South West) | 0.3261 | Africa (East)       | 0.1956 | PUR        | AMR             |
| HG01084   | ok     |   0.6375 |        108757 | South America       | 0.3011 | Ireland             | 0.2431 | Italy               | 0.2306 | PUR        | AMR             |
| HG01087   | ok     |   0.6465 |        108757 | South America       | 0.4387 | Scandinavia         | 0.2565 | Ashkenazi           | 0.0881 | PUR        | AMR             |
| HG01096   | ok     |   0.6446 |        108757 | Europe (South West) | 0.3881 | Ashkenazi           | 0.1920 | Finland             | 0.1733 | PUR        | AMR             |
| HG01099   | ok     |   0.6570 |        108757 | Europe (South West) | 0.6086 | Italy               | 0.1458 | Ireland             | 0.1258 | PUR        | AMR             |
| HG01100   | ok     |   0.6509 |        108757 | Africa (North)      | 0.2710 | Italy               | 0.1925 | Ireland             | 0.1522 | PUR        | AMR             |
| HG01103   | ok     |   0.6283 |        108757 | Europe (South West) | 0.6334 | Africa (South)      | 0.1100 | Africa (East)       | 0.0890 | PUR        | AMR             |
| HG01106   | ok     |   0.6439 |        108757 | Ireland             | 0.2272 | South America       | 0.1935 | Finland             | 0.1761 | PUR        | AMR             |
| HG01109   | ok     |   0.6479 |        108757 | Africa (West)       | 0.4509 | South America       | 0.3240 | Africa (South)      | 0.1624 | PUR        | AMR             |
| HG01169   | ok     |   0.6146 |        108757 | South America       | 0.4682 | Italy               | 0.2394 | Ashkenazi           | 0.1490 | PUR        | AMR             |
| HG01172   | ok     |   0.6307 |        108757 | Europe (South East) | 0.5496 | Europe (South West) | 0.2355 | South America       | 0.0724 | PUR        | AMR             |
| HG01175   | ok     |   0.6279 |        108757 | Italy               | 0.6059 | South America       | 0.2613 | Africa (East)       | 0.0655 | PUR        | AMR             |
| HG01178   | ok     |   0.6361 |        108757 | South America       | 0.3097 | Middle East         | 0.1665 | Italy               | 0.1645 | PUR        | AMR             |
| HG01184   | ok     |   0.6505 |        108757 | Europe (South West) | 0.6819 | South America       | 0.2070 | Africa (West)       | 0.0805 | PUR        | AMR             |
| HG01189   | ok     |   0.6347 |        108757 | South America       | 0.5956 | Africa (West)       | 0.1644 | Europe (North East) | 0.1264 | PUR        | AMR             |
| HG01192   | ok     |   0.6281 |        108757 | Italy               | 0.8386 | South America       | 0.1566 | Philippines         | 0.0047 | PUR        | AMR             |
| HG01199   | ok     |   0.6424 |        108757 | Ireland             | 0.3275 | Africa (North)      | 0.1839 | Europe (South West) | 0.1721 | PUR        | AMR             |
| HG01206   | ok     |   0.6565 |        108757 | Europe (South West) | 0.3573 | South America       | 0.2806 | Ireland             | 0.1623 | PUR        | AMR             |
| HG01243   | ok     |   0.6620 |        108757 | Africa (West)       | 0.4106 | Africa (South)      | 0.2391 | South America       | 0.1719 | PUR        | AMR             |
| HG01249   | ok     |   0.6283 |        108757 | Italy               | 0.4559 | South America       | 0.2677 | Scandinavia         | 0.2033 | PUR        | AMR             |
| NA18484   | ok     |   0.7114 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA18485   | ok     |   0.7183 |        108757 | Africa (West)       | 0.9697 | Africa (South)      | 0.0303 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA18497   | ok     |   0.7179 |        108757 | Africa (West)       | 0.9968 | Japan               | 0.0032 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA18500   | ok     |   0.7107 |        108757 | Africa (West)       | 0.9622 | Africa (South)      | 0.0378 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA18503   | ok     |   0.7102 |        108757 | Africa (West)       | 0.9963 | Africa (East)       | 0.0037 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA18506   | ok     |   0.6996 |        108757 | Africa (West)       | 0.9992 | Africa (East)       | 0.0008 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA18515   | ok     |   0.7245 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA18518   | ok     |   0.7146 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA18521   | ok     |   0.7204 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA18854   | ok     |   0.7188 |        108757 | Africa (West)       | 0.9997 | Asia (East)         | 0.0003 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA18857   | ok     |   0.6846 |        108757 | Africa (West)       | 0.9243 | Africa (East)       | 0.0757 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA18860   | ok     |   0.7104 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA18863   | ok     |   0.7124 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA18869   | ok     |   0.7108 |        108757 | Africa (West)       | 0.9502 | Africa (East)       | 0.0498 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA18872   | ok     |   0.7150 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA18875   | ok     |   0.6920 |        108757 | Africa (West)       | 0.9351 | Africa (South)      | 0.0383 | Asia (East)         | 0.0115 | YRI        | AFR             |
| NA18906   | ok     |   0.6982 |        108757 | Africa (West)       | 0.7249 | Africa (South)      | 0.2573 | Philippines         | 0.0078 | YRI        | AFR             |
| NA18911   | ok     |   0.7043 |        108757 | Africa (West)       | 0.9534 | Asia (East)         | 0.0265 | Africa (East)       | 0.0201 | YRI        | AFR             |
| NA18914   | ok     |   0.7061 |        108757 | Africa (West)       | 0.9945 | Sri Lanka           | 0.0039 | Asia (East)         | 0.0016 | YRI        | AFR             |
| NA18925   | ok     |   0.7015 |        108757 | Africa (West)       | 0.9968 | South America       | 0.0032 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA18930   | ok     |   0.7219 |        108757 | Africa (West)       | 0.9524 | Africa (South)      | 0.0476 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA18935   | ok     |   0.7147 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA19094   | ok     |   0.6992 |        108757 | Africa (West)       | 0.9829 | Sri Lanka           | 0.0171 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA19097   | ok     |   0.7164 |        108757 | Africa (West)       | 0.9920 | Japan               | 0.0060 | Africa (South)      | 0.0020 | YRI        | AFR             |
| NA19100   | ok     |   0.7099 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA19103   | ok     |   0.7155 |        108757 | Africa (West)       | 0.9855 | Finland             | 0.0145 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA19109   | ok     |   0.7052 |        108757 | Africa (West)       | 0.8242 | Africa (South)      | 0.1758 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA19115   | ok     |   0.7105 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA19120   | ok     |   0.7077 |        108757 | Africa (West)       | 0.9717 | Finland             | 0.0273 | Africa (South)      | 0.0010 | YRI        | AFR             |
| NA19123   | ok     |   0.6971 |        108757 | Africa (West)       | 0.8961 | Africa (South)      | 0.0991 | Finland             | 0.0036 | YRI        | AFR             |
| NA19129   | ok     |   0.7085 |        108757 | Africa (West)       | 0.9645 | Africa (South)      | 0.0355 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA19132   | ok     |   0.7206 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA19139   | ok     |   0.7031 |        108757 | Africa (West)       | 0.9106 | Africa (South)      | 0.0851 | Philippines         | 0.0026 | YRI        | AFR             |
| NA19142   | ok     |   0.7226 |        108757 | Africa (West)       | 0.8647 | Africa (South)      | 0.1353 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA19145   | ok     |   0.7059 |        108757 | Africa (West)       | 0.8398 | Africa (South)      | 0.1493 | Sri Lanka           | 0.0110 | YRI        | AFR             |
| NA19148   | ok     |   0.6920 |        108757 | Africa (West)       | 0.7918 | Africa (South)      | 0.1653 | Pakistan            | 0.0385 | YRI        | AFR             |
| NA19151   | ok     |   0.7060 |        108757 | Africa (West)       | 0.9669 | Africa (South)      | 0.0253 | Japan               | 0.0078 | YRI        | AFR             |
| NA19154   | ok     |   0.7105 |        108757 | Africa (West)       | 0.9967 | South America       | 0.0033 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA19161   | ok     |   0.6904 |        108757 | Africa (West)       | 0.9878 | Africa (North)      | 0.0122 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA19173   | ok     |   0.7039 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA19174   | ok     |   0.7043 |        108757 | Africa (West)       | 0.9879 | Africa (South)      | 0.0121 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA19177   | ok     |   0.7110 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA19186   | ok     |   0.7134 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA19191   | ok     |   0.7123 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA19199   | ok     |   0.7078 |        108757 | Africa (West)       | 0.9500 | Africa (South)      | 0.0500 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA19202   | ok     |   0.7208 |        108757 | Africa (West)       | 0.9911 | South America       | 0.0075 | Ashkenazi           | 0.0015 | YRI        | AFR             |
| NA19205   | ok     |   0.7123 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA19208   | ok     |   0.7009 |        108757 | Africa (West)       | 0.8988 | Africa (South)      | 0.1012 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA19211   | ok     |   0.7129 |        108757 | Africa (West)       | 0.8420 | Africa (South)      | 0.1580 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA19215   | ok     |   0.7070 |        108757 | Africa (West)       | 0.9978 | Japan               | 0.0022 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA19221   | ok     |   0.7209 |        108757 | Africa (West)       | 0.9800 | Ireland             | 0.0200 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA19224   | ok     |   0.7047 |        108757 | Africa (West)       | 0.8925 | Africa (South)      | 0.1075 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA19237   | ok     |   0.7119 |        108757 | Africa (West)       | 0.9980 | Ireland             | 0.0020 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA19240   | ok     |   0.7144 |        108757 | Africa (West)       | 0.9840 | Sri Lanka           | 0.0160 | Africa (East)       | 0.0000 | YRI        | AFR             |
| NA19249   | ok     |   0.7136 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |
| NA19258   | ok     |   0.7103 |        108757 | Africa (West)       | 1.0000 | Africa (East)       | 0.0000 | Africa (North)      | 0.0000 | YRI        | AFR             |

## ROH results

All count and length summaries are per-child quantities and are
aggregated by population only after retaining each child’s observations.
The chr20-only FROH column is reported as a chromosome-length fraction;
it is not the declared genome-wide FROH denominator.

| Population | Arm               | Unsupported ROH count / child (\>=1 Mb; \>=2 Mb; \>=5 Mb) | Unsupported ROH length / child (\>=1 Mb; \>=2 Mb; \>=5 Mb) | Supported count / child | Supported length / child |     FROH | Truth-run length covered |
|:-----------|:------------------|:----------------------------------------------------------|:-----------------------------------------------------------|------------------------:|-------------------------:|---------:|-------------------------:|
| ACB        | Pooled            | 0.350; 0.100; 0.000                                       | 670900; 356078; 0                                          |                   0.000 |                        0 | 0.151113 |                    21.52 |
| ACB        | Single population | 0.300; 0.050; 0.000                                       | 468436; 153856; 0                                          |                   0.000 |                        0 | 0.133340 |                    21.52 |
| ACB        | Ancestry tuned    | 0.300; 0.050; 0.000                                       | 468436; 153856; 0                                          |                   0.000 |                        0 | 0.132315 |                    21.52 |
| ASW        | Pooled            | 0.308; 0.231; 0.000                                       | 1031233; 950841; 0                                         |                   0.000 |                        0 | 0.161486 |                     0.00 |
| ASW        | Single population | 0.231; 0.154; 0.000                                       | 641463; 561376; 0                                          |                   0.000 |                        0 | 0.143329 |                     0.00 |
| ASW        | Ancestry tuned    | 0.308; 0.231; 0.000                                       | 956017; 875930; 0                                          |                   0.000 |                        0 | 0.147672 |                     0.00 |
| CLM        | Pooled            | 0.914; 0.200; 0.029                                       | 1447021; 565503; 192631                                    |                   0.029 |                   101071 | 0.259322 |                    30.95 |
| CLM        | Single population | 0.943; 0.200; 0.029                                       | 1474631; 565486; 192631                                    |                   0.029 |                   101071 | 0.244852 |                    30.95 |
| CLM        | Ancestry tuned    | 0.886; 0.200; 0.029                                       | 1413754; 565486; 192631                                    |                   0.029 |                   101071 | 0.244794 |                    30.95 |
| MXL        | Pooled            | 0.875; 0.125; 0.031                                       | 1325084; 421403; 211752                                    |                   0.031 |                    51818 | 0.281746 |                    28.30 |
| MXL        | Single population | 0.750; 0.125; 0.031                                       | 1185224; 420997; 211473                                    |                   0.031 |                    51818 | 0.267161 |                    28.30 |
| MXL        | Ancestry tuned    | 0.750; 0.125; 0.031                                       | 1188555; 420997; 211473                                    |                   0.031 |                    51818 | 0.267304 |                    28.30 |
| PEL        | Pooled            | 1.886; 0.257; 0.000                                       | 2765371; 637938; 0                                         |                   0.000 |                        0 | 0.361562 |                     7.81 |
| PEL        | Single population | 1.800; 0.257; 0.000                                       | 2662697; 637723; 0                                         |                   0.000 |                        0 | 0.341866 |                     7.81 |
| PEL        | Ancestry tuned    | 1.829; 0.257; 0.000                                       | 2704393; 637290; 0                                         |                   0.000 |                        0 | 0.339094 |                     7.81 |
| PUR        | Pooled            | 0.514; 0.057; 0.000                                       | 741818; 175457; 0                                          |                   0.057 |                   346568 | 0.229328 |                    59.65 |
| PUR        | Single population | 0.429; 0.057; 0.000                                       | 645491; 175714; 0                                          |                   0.057 |                   343779 | 0.216111 |                    59.65 |
| PUR        | Ancestry tuned    | 0.514; 0.057; 0.000                                       | 748088; 175714; 0                                          |                   0.057 |                   343779 | 0.217421 |                    59.65 |
| YRI        | Pooled            | 0.429; 0.179; 0.018                                       | 1048185; 703952; 145594                                    |                   0.000 |                        0 | 0.166384 |                    46.76 |
| YRI        | Single population | 0.393; 0.143; 0.018                                       | 897098; 561710; 146193                                     |                   0.000 |                        0 | 0.140664 |                    46.76 |
| YRI        | Ancestry tuned    | 0.393; 0.143; 0.018                                       | 897079; 561710; 146193                                     |                   0.000 |                        0 | 0.140369 |                    46.76 |
| ESN        | Pooled            | 0.512; 0.186; 0.093                                       | 1269975; 888192; 591068                                    |                   0.000 |                        0 | 0.170931 |                    33.85 |
| ESN        | Single population | 0.419; 0.163; 0.070                                       | 1071104; 764335; 467527                                    |                   0.000 |                        0 | 0.143510 |                    15.38 |
| ESN        | Ancestry tuned    | 0.419; 0.163; 0.070                                       | 1071104; 764335; 467527                                    |                   0.000 |                        0 | 0.142799 |                    15.38 |
| CEU        | Pooled            | 0.649; 0.070; 0.018                                       | 993521; 293245; 105846                                     |                   0.000 |                        0 | 0.286503 |                    12.17 |
| CEU        | Single population | 0.614; 0.070; 0.018                                       | 953874; 293245; 105846                                     |                   0.000 |                        0 | 0.260306 |                    12.17 |
| CEU        | Ancestry tuned    | 0.632; 0.070; 0.018                                       | 966819; 293245; 105846                                     |                   0.000 |                        0 | 0.261963 |                    12.17 |
| CHS        | Pooled            | 0.843; 0.078; 0.000                                       | 1194930; 170283; 0                                         |                   0.000 |                        0 | 0.329422 |                     0.00 |
| CHS        | Single population | 0.824; 0.078; 0.000                                       | 1150715; 170482; 0                                         |                   0.000 |                        0 | 0.292096 |                     0.00 |
| CHS        | Ancestry tuned    | 0.824; 0.078; 0.000                                       | 1150742; 170510; 0                                         |                   0.000 |                        0 | 0.294230 |                     0.00 |

## Criterion

The chr20 criterion **did not hold**. The comparison for each population
is shown below.

| Population | comparison        | ancestry_unsupported_count | comparator_unsupported_count | ancestry_truth_coverage | comparator_truth_coverage | count_pass | truth_coverage_pass | criterion_pass |
|:-----------|:------------------|---------------------------:|-----------------------------:|------------------------:|--------------------------:|:-----------|:--------------------|:---------------|
| ACB        | pooled            |                     0.3000 |                       0.3500 |                  0.2152 |                    0.2152 | TRUE       | TRUE                | TRUE           |
| ACB        | single_population |                     0.3000 |                       0.3000 |                  0.2152 |                    0.2152 | FALSE      | TRUE                | FALSE          |
| ASW        | pooled            |                     0.3077 |                       0.3077 |                  0.0000 |                    0.0000 | FALSE      | TRUE                | FALSE          |
| ASW        | single_population |                     0.3077 |                       0.2308 |                  0.0000 |                    0.0000 | FALSE      | TRUE                | FALSE          |
| CEU        | single_population |                     0.6316 |                       0.6140 |                  0.1217 |                    0.1217 | TRUE       | TRUE                | TRUE           |
| CHS        | single_population |                     0.8235 |                       0.8235 |                  0.0000 |                    0.0000 | TRUE       | TRUE                | TRUE           |
| CLM        | pooled            |                     0.8857 |                       0.9143 |                  0.3095 |                    0.3095 | TRUE       | TRUE                | TRUE           |
| CLM        | single_population |                     0.8857 |                       0.9429 |                  0.3095 |                    0.3095 | TRUE       | TRUE                | TRUE           |
| ESN        | single_population |                     0.4186 |                       0.4186 |                  0.1538 |                    0.1538 | TRUE       | TRUE                | TRUE           |
| MXL        | pooled            |                     0.7500 |                       0.8750 |                  0.2830 |                    0.2830 | TRUE       | TRUE                | TRUE           |
| MXL        | single_population |                     0.7500 |                       0.7500 |                  0.2830 |                    0.2830 | FALSE      | TRUE                | FALSE          |
| PEL        | pooled            |                     1.8286 |                       1.8857 |                  0.0781 |                    0.0781 | TRUE       | TRUE                | TRUE           |
| PEL        | single_population |                     1.8286 |                       1.8000 |                  0.0781 |                    0.0781 | FALSE      | TRUE                | FALSE          |
| PUR        | pooled            |                     0.5143 |                       0.5143 |                  0.5965 |                    0.5965 | FALSE      | TRUE                | FALSE          |
| PUR        | single_population |                     0.5143 |                       0.4286 |                  0.5965 |                    0.5965 | FALSE      | TRUE                | FALSE          |
| YRI        | single_population |                     0.3929 |                       0.3929 |                  0.4676 |                    0.4676 | TRUE       | TRUE                | TRUE           |

## Runtime and memory

Each value is from a fresh R process using four DuckDB threads. The
ancestry arm used 16-child batches after the 64-child batch exceeded the
12 GB DuckDB memory limit; other arms used 64-child batches.
Per-replicate runtime and peak resident memory are reported below.

| chromosome | arm        | repetition | elapsed_seconds | peak_rss_gib | threads | batch_size |
|-----------:|:-----------|-----------:|----------------:|-------------:|--------:|-----------:|
|         20 | pooled     |          1 |          126.59 |        12.44 |       4 |         64 |
|         20 | pooled     |          2 |          126.74 |        12.47 |       4 |         64 |
|         20 | pooled     |          3 |          125.21 |        12.43 |       4 |         64 |
|         20 | single_AFR |          1 |           44.90 |        12.53 |       4 |         64 |
|         20 | single_AFR |          2 |           44.88 |        12.86 |       4 |         64 |
|         20 | single_AFR |          3 |           44.55 |        12.18 |       4 |         64 |
|         20 | single_AMR |          1 |           46.85 |        12.27 |       4 |         64 |
|         20 | single_AMR |          2 |           46.69 |        12.35 |       4 |         64 |
|         20 | single_AMR |          3 |           46.84 |        12.40 |       4 |         64 |
|         20 | single_EUR |          1 |           19.37 |        10.75 |       4 |         64 |
|         20 | single_EUR |          2 |           19.25 |        10.77 |       4 |         64 |
|         20 | single_EUR |          3 |           19.06 |        10.75 |       4 |         64 |
|         20 | single_EAS |          1 |           17.08 |        10.14 |       4 |         64 |
|         20 | single_EAS |          2 |           17.37 |        10.16 |       4 |         64 |
|         20 | single_EAS |          3 |           17.38 |        10.12 |       4 |         64 |
|         20 | ancestry   |          1 |          208.83 |         4.58 |       4 |         16 |
|         20 | ancestry   |          2 |          208.96 |         4.58 |       4 |         16 |
|         20 | ancestry   |          3 |          208.90 |         4.54 |       4 |         16 |
