Base mismatch count scaling
================

This report gives the 1×/2×/4× scaling evidence `STYLE.md` asks for, for
`duckhts_bam_mismatch_counts()`. The function reads alignments once and counts
their aligned bases against the reference by mate, cycle, base quality and
substitution. Its state is meant to be fixed, so the memory verdict is the main
one. No earlier report measures this function; this is its first baseline.

## Workload

The workload is real and public: chr1, the largest contig, of the
1000 Genomes 30× CRAM of HG00403, staged as a CRAM slice (registry artifact
`bam_mismatch_hg00403_chr1_cram`), with every chr1 record of the phased
1000 Genomes panel as the mask (`bam_mismatch_chr1_mask_bcf`) and the
registered GRCh38 FASTA as the reference (`liftover_grch38_fasta`).
`benchmarks/bam_mismatch_stage.R` stages both derived files and checks them
against their registered identities. This report checks them again before it
measures, and it checks the reference against its registered size and MD5.
Those checks read all three files in full, so every measured run reads from a
warm page cache.

There are two growth dimensions: the aligned bases and the mask records.

- The `mask` series doubles both together. Its 1× and 2× regions end at the
  25% and 50% quantiles of the alignment start positions of the slice, and its
  4× region is the whole contig. The mask records under a region grow with it.
- The `no mask` series reads the same regions without a mask, so the aligned
  bases grow alone.

The mask records are not doubled alone. The difference between the two series
is the whole cost of the mask, and the report gives it at each size.

The cuts and the input denominators of each region are computed at render
time. That is untimed staging, as is the one-time staging of the two files.

Every observation is one fresh DuckDB CLI process under `/usr/bin/time`. The
JSON profiler gives the query latency and DuckDB’s peak buffer and temporary
bytes, and `/usr/bin/time` gives the peak RSS. The host is shared, so the
render command pins the processes to eight cores and each measured process
runs at scheduling priority −15 (`nice`), which makes other load on those
cores yield to it; the host load at the start is recorded, and the spread of
the repetitions shows what interference remained. `max_temp_directory_size` is
zero, so a spill fails the render; DuckDB’s memory limit is left at its
default, which the report records. Each cell has three repetitions at one
thread and at four threads, run in an order that alternates sizes. The
function has one worker, so the thread setting should not change it, and the
four-thread cells show whether that holds for time and memory. Each query
reduces the rows of the function to their number, the counted bases and the
mismatches, so every row is produced and consumed. The R wrapper only builds
this SQL call; it is not timed here.

## Budget, declared before measuring

The overhead is the median peak RSS of three empty runs that start the CLI,
load the extension and scan nothing.

The ceiling is overhead + 3 × the live state. It is the same at every size,
and it is a gate: the render stops if any observation exceeds it or spills.
The live state is:

- the counters: 3 mates × 1,001 cycles × 95 base qualities × 16 substitutions
  × 8 bytes = 36,516,480 bytes;
- one reference window of 1 MiB and its mask of 1 MiB;
- htslib’s state for the three files: the CRAM container in work with its
  decoded alignments, the reference of that container’s span, the indexes and
  the file buffers. 16 MiB is allowed for it. The function sets no decoder
  threads, and for sorted input without them htslib’s `cram_get_ref()` loads
  the reference of one slice’s span, not the contig.

That is 52.8 MiB, so the ceiling is overhead + 158.5 MiB.

The render also stops if a time exponent is above 1.25 at either doubling, or
if the 1× run on one thread takes less than 5 s, which is the least that
`STYLE.md` accepts for a timing verdict.

``` sh
taskset -c 8-15 Rscript -e 'rmarkdown::render("benchmarks/benchmark_bam_mismatch_counts.Rmd")'
```

Source revision: 37b438ca55a09d335c6d2fb218ca053730930ec4; `src` tree 1f7c4e3dbbd535d8b1ea6f59cedf0993ca1eb4e4. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: f8e0276f8d94671a91700e8142c281cdb551f9cc9efbb158dd9f05cc375d4675. Host load at start: 6.17, 6.90, 7.42. DuckDB memory limit: 50.0 GiB (the default); temporary directory limit: 0.

Staged inputs: the slice has 62,660,978 alignments (position sum 7728419844834032, 9,399,146,700 stored bases) in 6,270 CRAM slices, of which the widest starts at 125,183,655 and spans 18,001,341 reference bases; the mask has 5,769,087 records.

Overhead: 38.6 MiB. Ceiling: 197 MiB.

## Inputs of each region

`alignments` and `stored_bases` are what the region selects before the
function’s own filters; `decoded_base_mib` is their bases and qualities as
htslib hands them over, at 1.5 bytes for each stored base. `mask_records` are
the mask records that the region selects.

| scale | region           | alignments |  stored_bases | decoded_base_mib | mask_records |
|------:|:-----------------|-----------:|--------------:|-----------------:|-------------:|
|     1 | chr1:1-61075319  | 15,665,245 | 2,349,786,750 |          3,361.4 |    1,567,146 |
|     2 | chr1:1-123524417 | 31,330,489 | 4,699,573,350 |          6,722.8 |    3,074,594 |
|     4 | chr1:1-248956422 | 62,660,978 | 9,399,146,700 |         13,445.6 |    5,769,087 |

## Measurements

`bases`, `mismatches` and `output_rows` are the result of the query. Times are
the query latency of three repetitions.

| series  | scale | threads |         bases | mismatches | output_rows | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | max_buffer_mib | budget_mib |
|:--------|------:|--------:|--------------:|-----------:|------------:|---------------:|------------:|------------:|---------------:|------------:|---------------:|-----------:|
| mask    |     1 |       1 | 1,920,445,937 |  4,034,678 |      17,756 |         16.552 |      15.451 |      20.093 |           67.6 |        67.6 |            0.3 |        197 |
| mask    |     2 |       1 | 3,733,164,470 |  7,164,101 |      17,810 |         32.380 |      30.514 |      38.069 |           67.4 |        67.5 |            0.3 |        197 |
| mask    |     4 |       1 | 7,053,095,954 | 14,515,783 |      17,839 |         65.996 |      59.774 |      80.851 |           83.5 |        83.5 |            0.3 |        197 |
| mask    |     1 |       4 | 1,920,445,937 |  4,034,678 |      17,756 |         15.969 |      15.414 |      31.107 |           68.8 |        68.8 |            0.3 |        197 |
| mask    |     2 |       4 | 3,733,164,470 |  7,164,101 |      17,810 |         35.274 |      31.428 |      42.297 |           67.8 |        68.5 |            0.3 |        197 |
| mask    |     4 |       4 | 7,053,095,954 | 14,515,783 |      17,839 |         68.034 |      61.035 |      77.057 |           86.3 |        86.6 |            0.3 |        197 |
| no mask |     1 |       1 | 2,015,065,094 |  6,059,797 |      17,807 |         15.296 |      14.131 |      21.515 |           65.7 |        66.1 |            0.3 |        197 |
| no mask |     2 |       1 | 3,908,455,277 | 11,109,173 |      17,832 |         31.173 |      27.690 |      36.799 |           65.8 |        67.0 |            0.3 |        197 |
| no mask |     4 |       1 | 7,378,304,736 | 21,938,306 |      17,839 |         60.987 |      54.096 |      72.025 |           80.0 |        81.3 |            0.3 |        197 |
| no mask |     1 |       4 | 2,015,065,094 |  6,059,797 |      17,807 |         14.494 |      14.372 |      19.483 |           66.8 |        67.4 |            0.3 |        197 |
| no mask |     2 |       4 | 3,908,455,277 | 11,109,173 |      17,832 |         28.354 |      28.297 |      41.230 |           67.3 |        67.4 |            0.3 |        197 |
| no mask |     4 |       4 | 7,378,304,736 | 21,938,306 |      17,839 |         55.035 |      54.060 |      70.933 |           83.1 |        85.7 |            0.3 |        197 |

Growth of the median time against the bases counted (an exponent of 1 is linear), and peak-RSS ratios to 1×:

| series  | threads | one_x_seconds | bases_ratio_2x | bases_ratio_4x | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x |
|:--------|--------:|--------------:|---------------:|---------------:|-----------------:|-----------------:|-------------:|-------------:|
| mask    |       1 |        16.552 |          1.944 |          1.889 |            1.010 |            1.119 |        0.997 |        1.235 |
| mask    |       4 |        15.969 |          1.944 |          1.889 |            1.192 |            1.032 |        0.985 |        1.254 |
| no mask |       1 |        15.296 |          1.940 |          1.888 |            1.075 |            1.056 |        1.002 |        1.218 |
| no mask |       4 |        14.494 |          1.940 |          1.888 |            1.013 |            1.044 |        1.007 |        1.244 |

The cost of the mask on one thread, as the difference of the two series:

| scale | mask_records | seconds_with_mask | seconds_without_mask | mask_seconds | bases_left_out | rss_mib_with_mask | rss_mib_without_mask |
|------:|-------------:|------------------:|---------------------:|-------------:|---------------:|------------------:|---------------------:|
|     1 |    1,567,146 |            16.552 |               15.296 |        1.255 |     94,619,157 |              67.6 |                 65.7 |
|     2 |    3,074,594 |            32.380 |               31.173 |        1.207 |    175,290,807 |              67.4 |                 65.8 |
|     4 |    5,769,087 |            65.996 |               60.987 |        5.009 |    325,208,782 |              83.5 |                 80.0 |

## Findings

Every observation is within the ceiling and none spills. The 1× runs on one
thread take at least 5 s and no time exponent is above 1.25. The render
enforces all of this. It also stops unless the number of rows, the counted
bases and the mismatches of a region are the same in every repetition and at
both thread settings. This report does not check the counts themselves; the
SQL tests do, against counts derived by hand and against the same counts made
independently in SQL.

With the mask on one thread, the whole contig takes
66 s for 7,053,095,954
counted bases, which is 107 million counted bases
a second, at 83.5 MiB peak RSS. The 1× region takes
16.6 s at 67.6 MiB.
Peak RSS at 4× is 1.24
times the 1× value, and the largest peak RSS of any run is
86.6 MiB against the ceiling of 197 MiB.

The difference between the sizes in peak RSS is htslib’s, not the function’s.
The decoder holds the reference of the CRAM slice in hand, and the widest
slice of this file spans 18,001,341 bases
(17.2 MiB of reference) from position
125,183,655, over the sparse region of the
contig. Only the whole-contig runs read that slice. That one span is above the
16 MiB the budget allows for all of htslib’s state, so the allowance was too
small for this file; the ceiling, which is what gates the render, held with a
wide margin. The step depends on how the file’s slices are laid out, not on
the number of alignments.

The operating scale this report supports is what it measures: one contig of
one 30× short-read sample, read from CRAM, with a mask of
5,769,087 records. It does not measure
a whole genome in one call, BAM input, unsorted input, long reads, or a mask
whose records are doubled alone.

## What the counts say

The whole contig with the mask, by base quality, from one untimed run:

| base_quality |         bases | mismatches | mismatch_rate | share_of_bases | share_of_mismatches |
|-------------:|--------------:|-----------:|--------------:|---------------:|--------------------:|
|            3 |        72,739 |     30,288 |      0.416393 |         0.0000 |              0.0021 |
|            4 |       932,553 |    355,352 |      0.381053 |         0.0001 |              0.0245 |
|            5 |       774,048 |    276,661 |      0.357421 |         0.0001 |              0.0191 |
|            6 |     3,648,377 |    652,581 |      0.178869 |         0.0005 |              0.0450 |
|           10 |    95,321,180 |  9,174,850 |      0.096252 |         0.0135 |              0.6321 |
|           20 |    48,665,990 |    959,063 |      0.019707 |         0.0069 |              0.0661 |
|           30 | 6,903,681,067 |  3,066,988 |      0.000444 |         0.9788 |              0.2113 |

The mask leaves out every position of a panel record, and the sample is one of
the panel’s samples. So these mismatches are errors of the read, the library
or the alignment, together with the variants of the sample that the panel
does not hold.
