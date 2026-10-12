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

There are two growth dimensions, the aligned bases and the mask records, and
three series.

- The `mask` series doubles both together. Its 1× and 2× regions end at the
  25% and 50% quantiles of the alignment start positions of the slice, and its
  4× region is the whole contig. The mask records under a region grow with it.
- The `no mask` series reads the same regions without a mask, so the aligned
  bases grow alone. The difference between the two series is the whole cost
  of the mask at each size.
- The `mask records` series doubles the mask records alone. Its alignments
  are fixed: one two-base region every 2 MiB of the contig. The function reads
  the mask one reference window at a time, and a window starts at the first
  kept alignment that the last window does not hold and is 1 MiB long. A
  region 2 MiB after the last one is therefore outside the last window, and
  each region opens one window of its own. The render derives these windows
  from the first kept alignment of each region, checks that no alignment of a
  region leaves its window or reaches the next one, and counts the mask
  records that the windows hold; that count, not the size of the mask file,
  is the denominator of the series. The windows cover about half of the
  contig. The masks are the same panel records thinned to the positions
  divisible by 4 (`bam_mismatch_chr1_mask_quarter_bcf`), by 2
  (`bam_mismatch_chr1_mask_half_bcf`) and not at all, so 1×, 2× and 4×
  records; a run without a mask gives the cost of the alignments alone.

The cuts and the input denominators of each region are computed at render
time. That is untimed staging, as is the one-time staging of the files.

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
- one reference window of 1 MiB and its mask of 1 MiB. The window moves along
  an alignment that reaches past it, so no alignment needs more;
- htslib’s state for the three files: the CRAM container in work with its
  decoded alignments, the reference of that container’s span, the indexes and
  the file buffers. 16 MiB is allowed for it. The function sets no decoder
  threads, and for sorted input without them htslib’s `cram_get_ref()` loads
  the reference of one slice’s span, not the contig.

That is 52.8 MiB, so the ceiling is overhead + 158.5 MiB.

The render also stops if a time exponent of the `mask` or `no mask` series is
above 1.25 at either doubling, or if their 1× run on one thread takes less
than 5 s, which is the least that `STYLE.md` accepts for a timing verdict.
The `mask records` series reads a few million records in a few seconds, so
its 1× run is expected under 5 s: it gets the memory verdict in full and its
times are reported without a verdict, as `STYLE.md` provides.

``` sh
taskset -c 8-15 Rscript -e 'rmarkdown::render("benchmarks/benchmark_bam_mismatch_counts.Rmd")'
```

Source revision: d1a678ce2830e60e48ecb226993207cc9f6116cd; `src` tree b0474026510f692da693fe57fd8d1c56d79cd70b. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: 5fd33730e8af91c365ac3dcf58091e3faf3b2d49c7849ef283a8db8450bf4bfe. Host load at start: 0.75, 0.27, 0.19. DuckDB memory limit: 50.0 GiB (the default); temporary directory limit: 0.

Staged inputs: the slice has 62,660,978 alignments (position sum 7728419844834032, 9,399,146,700 stored bases) in 6,270 CRAM slices, of which the widest starts at 125,183,655 and spans 18,001,341 reference bases; the mask has 5,769,087 records, its half 2,884,829 and its quarter 1,441,811.

Overhead: 38.7 MiB. Ceiling: 197.2 MiB.

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
the query latency of three repetitions. The two alignment series:

| series  | scale | threads |         bases | mismatches | output_rows | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | max_buffer_mib | budget_mib |
|:--------|------:|--------:|--------------:|-----------:|------------:|---------------:|------------:|------------:|---------------:|------------:|---------------:|-----------:|
| mask    |     1 |       1 | 1,920,445,937 |  4,034,678 |      17,756 |         16.055 |      15.869 |      17.060 |           68.0 |        68.2 |            0.3 |      197.2 |
| mask    |     2 |       1 | 3,733,164,470 |  7,164,101 |      17,810 |         31.507 |      31.119 |      33.833 |           67.9 |        68.1 |            0.3 |      197.2 |
| mask    |     4 |       1 | 7,053,095,954 | 14,515,783 |      17,839 |         61.220 |      61.075 |      65.124 |           83.0 |        83.2 |            0.3 |      197.2 |
| mask    |     1 |       4 | 1,920,445,937 |  4,034,678 |      17,756 |         16.362 |      16.096 |      16.971 |           68.9 |        69.4 |            0.2 |      197.2 |
| mask    |     2 |       4 | 3,733,164,470 |  7,164,101 |      17,810 |         34.083 |      31.603 |      34.770 |           69.1 |        69.2 |            0.3 |      197.2 |
| mask    |     4 |       4 | 7,053,095,954 | 14,515,783 |      17,839 |         65.418 |      61.490 |      65.682 |           83.5 |        86.4 |            0.3 |      197.2 |
| no mask |     1 |       1 | 2,015,065,094 |  6,059,797 |      17,807 |         14.602 |      14.534 |      15.047 |           66.0 |        66.2 |            0.3 |      197.2 |
| no mask |     2 |       1 | 3,908,455,277 | 11,109,173 |      17,832 |         28.836 |      28.437 |      31.227 |           66.2 |        66.2 |            0.3 |      197.2 |
| no mask |     4 |       1 | 7,378,304,736 | 21,938,306 |      17,839 |         56.816 |      56.287 |      58.224 |           81.8 |        81.9 |            0.3 |      197.2 |
| no mask |     1 |       4 | 2,015,065,094 |  6,059,797 |      17,807 |         15.224 |      14.861 |      16.410 |           67.4 |        67.6 |            0.3 |      197.2 |
| no mask |     2 |       4 | 3,908,455,277 | 11,109,173 |      17,832 |         29.728 |      29.404 |      31.938 |           66.3 |        66.4 |            0.3 |      197.2 |
| no mask |     4 |       4 | 7,378,304,736 | 21,938,306 |      17,839 |         57.704 |      56.657 |      62.228 |           82.4 |        83.9 |            0.3 |      197.2 |

Growth of the median time against the bases counted (an exponent of 1 is linear), and peak-RSS ratios to 1×:

| series  | threads | one_x_seconds | bases_ratio_2x | bases_ratio_4x | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x |
|:--------|--------:|--------------:|---------------:|---------------:|-----------------:|-----------------:|-------------:|-------------:|
| mask    |       1 |        16.055 |          1.944 |          1.889 |            1.014 |            1.044 |        0.999 |        1.221 |
| mask    |       4 |        16.362 |          1.944 |          1.889 |            1.104 |            1.025 |        1.003 |        1.212 |
| no mask |       1 |        14.602 |          1.940 |          1.888 |            1.027 |            1.067 |        1.003 |        1.239 |
| no mask |       4 |        15.224 |          1.940 |          1.888 |            1.010 |            1.044 |        0.984 |        1.223 |

The cost of the mask on one thread, as the difference of the two series:

| scale | mask_records | seconds_with_mask | seconds_without_mask | mask_seconds | bases_left_out | rss_mib_with_mask | rss_mib_without_mask |
|------:|-------------:|------------------:|---------------------:|-------------:|---------------:|------------------:|---------------------:|
|     1 |    1,567,146 |            16.055 |               14.602 |        1.453 |     94,619,157 |              68.0 |                 66.0 |
|     2 |    3,074,594 |            31.507 |               28.836 |        2.672 |    175,290,807 |              67.9 |                 66.2 |
|     4 |    5,769,087 |            61.220 |               56.816 |        4.404 |    325,208,782 |              83.0 |                 81.8 |

The mask-records series, at fixed alignments. Its 109 windows
hold 114,294,784 reference
bases, and `mask_records` are the records of each mask that those windows
hold, out of 5,769,087 in the full
mask file. Scale 0 is the same regions without a mask:

| scale | threads | mask_records |   bases | mismatches | output_rows | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | budget_mib |
|------:|--------:|-------------:|--------:|-----------:|------------:|---------------:|------------:|------------:|---------------:|------------:|-----------:|
|     0 |       1 |            0 | 566,202 |      1,784 |       4,907 |          0.569 |       0.546 |       0.580 |           76.1 |        76.1 |      197.2 |
|     1 |       1 |      716,832 | 560,511 |      1,646 |       4,781 |          1.100 |       1.072 |       1.103 |           77.5 |        77.7 |      197.2 |
|     2 |       1 |    1,434,126 | 554,837 |      1,566 |       4,704 |          1.609 |       1.559 |       1.638 |           77.4 |        77.7 |      197.2 |
|     4 |       1 |    2,866,198 | 546,018 |      1,292 |       4,468 |          2.615 |       2.567 |       2.674 |           78.7 |        78.8 |      197.2 |
|     0 |       4 |            0 | 566,202 |      1,784 |       4,907 |          0.584 |       0.577 |       0.603 |           76.3 |        76.9 |      197.2 |
|     1 |       4 |      716,832 | 560,511 |      1,646 |       4,781 |          1.099 |       1.098 |       1.100 |           78.0 |        78.9 |      197.2 |
|     2 |       4 |    1,434,126 | 554,837 |      1,566 |       4,704 |          1.608 |       1.602 |       1.672 |           77.7 |        79.1 |      197.2 |
|     4 |       4 |    2,866,198 | 546,018 |      1,292 |       4,468 |          2.630 |       2.596 |       2.648 |           77.7 |        77.8 |      197.2 |

Its growth against the mask records, for the whole run and for the mask’s own cost (the run without a mask taken off), without a timing verdict:

| threads | no_mask_seconds | one_x_seconds | records_ratio_2x | records_ratio_4x | time_exponent_2x | time_exponent_4x | mask_cost_exponent_2x | mask_cost_exponent_4x | rss_ratio_4x |
|--------:|----------------:|--------------:|-----------------:|-----------------:|-----------------:|-----------------:|----------------------:|----------------------:|-------------:|
|       1 |           0.569 |         1.100 |            2.001 |            1.999 |            0.548 |            0.701 |                 0.969 |                 0.977 |        1.015 |
|       4 |           0.584 |         1.099 |            2.001 |            1.999 |            0.549 |            0.711 |                 0.993 |                 1.000 |        0.996 |

## Findings

Every observation is within the ceiling and none spills. The 1× runs of the
two alignment series on one thread take at least 5 s and none of their time
exponents is above 1.25. The render enforces all of this. It also stops
unless the number of rows, the counted bases and the mismatches of a region
are the same in every repetition and at both thread settings. This report
does not check the counts themselves; the SQL tests do, against counts
derived by hand and against the same counts made independently in SQL.

Mask records alone: reading the 2,866,198
records that the 109 windows hold in the full mask costs
2.05 s
on one thread (1.4 million records a
second) on top of the 0.57 s of the fixed
alignments, and peak RSS is 76.1 MiB without a mask
and 78.7 MiB with the full mask: the mask window is
1 MiB whatever the mask holds. The whole-run time exponents at one thread are
0.55 and
0.7, and
the mask’s own cost grows with exponents
0.97 and
0.98.
These times carry no verdict, because the 1× run takes
1.1 s.

With the mask on one thread, the whole contig takes
61.2 s for 7,053,095,954
counted bases, which is 115 million counted bases
a second, at 83 MiB peak RSS. The 1× region takes
16.1 s at 68 MiB.
Peak RSS at 4× is 1.22
times the 1× value, and the largest peak RSS of any run is
86.4 MiB against the ceiling of 197.2 MiB.

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
a whole genome in one call, BAM input, unsorted input or long reads. It also
does not time an alignment that reaches past the reference window, which the
function serves by moving the window: this file has
0 such alignments, and
its widest alignment spans 216
reference bases. That path costs a window fetch on each side of the alignment
and no memory; the SQL tests check its counts on a reference skip of
1,100,000 bases.

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
does not hold. Their rate is an upper bound on the error rate, not the error
rate.

## Change from the previous render

The previous render, 66899c92, measured the same cells with the
same driver on the same host and cores. Host load differs between renders.
Median of the repetitions:

| series       | scale | threads | seconds_before | rss_mib_before | seconds_after | rss_mib_after | time_ratio | rss_ratio |
|:-------------|------:|--------:|---------------:|---------------:|--------------:|--------------:|-----------:|----------:|
| mask         |     1 |       1 |         15.829 |         68.148 |        16.055 |        67.953 |      1.014 |     0.997 |
| mask         |     1 |       4 |         15.980 |         69.637 |        16.362 |        68.949 |      1.024 |     0.990 |
| mask         |     2 |       1 |         31.053 |         67.957 |        31.507 |        67.938 |      1.015 |     1.000 |
| mask         |     2 |       4 |         30.954 |         67.953 |        34.083 |        69.070 |      1.101 |     1.016 |
| mask         |     4 |       1 |         60.498 |         81.867 |        61.220 |        83.000 |      1.012 |     1.014 |
| mask         |     4 |       4 |         62.225 |         84.758 |        65.418 |        83.457 |      1.051 |     0.985 |
| mask records |     0 |       1 |          0.548 |         76.066 |         0.569 |        76.090 |      1.038 |     1.000 |
| mask records |     0 |       4 |          0.562 |         76.125 |         0.584 |        76.266 |      1.040 |     1.002 |
| mask records |     1 |       1 |          1.069 |         77.477 |         1.100 |        77.496 |      1.029 |     1.000 |
| mask records |     1 |       4 |          1.093 |         78.156 |         1.099 |        77.965 |      1.005 |     0.998 |
| mask records |     2 |       1 |          1.567 |         77.844 |         1.609 |        77.414 |      1.027 |     0.994 |
| mask records |     2 |       4 |          1.606 |         77.957 |         1.608 |        77.676 |      1.001 |     0.996 |
| mask records |     4 |       1 |          2.570 |         77.855 |         2.615 |        78.723 |      1.018 |     1.011 |
| mask records |     4 |       4 |          2.635 |         77.902 |         2.630 |        77.656 |      0.998 |     0.997 |
| no mask      |     1 |       1 |         14.426 |         66.398 |        14.602 |        66.047 |      1.012 |     0.995 |
| no mask      |     1 |       4 |         14.795 |         67.730 |        15.224 |        67.355 |      1.029 |     0.994 |
| no mask      |     2 |       1 |         28.250 |         67.074 |        28.836 |        66.156 |      1.021 |     0.986 |
| no mask      |     2 |       4 |         28.464 |         66.875 |        29.728 |        66.262 |      1.044 |     0.991 |
| no mask      |     4 |       1 |         55.900 |         81.914 |        56.816 |        81.828 |      1.016 |     0.999 |
| no mask      |     4 |       4 |         57.026 |         82.363 |        57.704 |        82.441 |      1.012 |     1.001 |
