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
- one reference window of 1 MiB and its mask of 1 MiB. An alignment that
  spans more reference than that grows both to its span, up to the 64 MiB
  limit; the render checks from the CIGARs that no alignment of this file
  spans more than 1 MiB, so the window stays at 1 MiB;
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

Source revision: 79a839cc3c1b91e9ce5222b7a69210f5ae650b07; `src` tree 9d1d5e9d39cd3d3bef47844938d031636fd16573. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: 1615882a792790b46e1f87cd2c80890efc748a315ca64e341064aab1eca51bb1. Host load at start: 1.18, 0.97, 0.78. DuckDB memory limit: 50.0 GiB (the default); temporary directory limit: 0.

Staged inputs: the slice has 62,660,978 alignments (position sum 7728419844834032, 9,399,146,700 stored bases) in 6,270 CRAM slices, of which the widest starts at 125,183,655 and spans 18,001,341 reference bases; the mask has 5,769,087 records, its half 2,884,829 and its quarter 1,441,811.

Overhead: 39 MiB. Ceiling: 197.5 MiB.

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
| mask    |     1 |       1 | 1,920,445,937 |  4,034,678 |      17,756 |         15.894 |      15.802 |      15.934 |           67.9 |        68.0 |            0.3 |      197.5 |
| mask    |     2 |       1 | 3,733,164,470 |  7,164,101 |      17,810 |         31.074 |      30.959 |      31.327 |           68.1 |        68.2 |            0.3 |      197.5 |
| mask    |     4 |       1 | 7,053,095,954 | 14,515,783 |      17,839 |         61.313 |      61.195 |      61.338 |           83.0 |        83.1 |            0.3 |      197.5 |
| mask    |     1 |       4 | 1,920,445,937 |  4,034,678 |      17,756 |         16.051 |      15.885 |      16.320 |           69.1 |        69.6 |            0.3 |      197.5 |
| mask    |     2 |       4 | 3,733,164,470 |  7,164,101 |      17,810 |         31.642 |      31.631 |      32.072 |           69.3 |        69.6 |            0.3 |      197.5 |
| mask    |     4 |       4 | 7,053,095,954 | 14,515,783 |      17,839 |         61.594 |      61.415 |      63.511 |           83.8 |        84.2 |            0.3 |      197.5 |
| no mask |     1 |       1 | 2,015,065,094 |  6,059,797 |      17,807 |         14.266 |      14.220 |      14.437 |           66.3 |        66.5 |            0.3 |      197.5 |
| no mask |     2 |       1 | 3,908,455,277 | 11,109,173 |      17,832 |         28.650 |      28.367 |      28.753 |           67.3 |        67.4 |            0.3 |      197.5 |
| no mask |     4 |       1 | 7,378,304,736 | 21,938,306 |      17,839 |         55.982 |      55.894 |      56.029 |           82.0 |        82.1 |            0.2 |      197.5 |
| no mask |     1 |       4 | 2,015,065,094 |  6,059,797 |      17,807 |         14.589 |      14.566 |      14.619 |           67.4 |        67.8 |            0.3 |      197.5 |
| no mask |     2 |       4 | 3,908,455,277 | 11,109,173 |      17,832 |         28.688 |      28.632 |      28.994 |           67.4 |        68.7 |            0.3 |      197.5 |
| no mask |     4 |       4 | 7,378,304,736 | 21,938,306 |      17,839 |         56.639 |      56.271 |      57.310 |           82.5 |        85.3 |            0.3 |      197.5 |

Growth of the median time against the bases counted (an exponent of 1 is linear), and peak-RSS ratios to 1×:

| series  | threads | one_x_seconds | bases_ratio_2x | bases_ratio_4x | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x |
|:--------|--------:|--------------:|---------------:|---------------:|-----------------:|-----------------:|-------------:|-------------:|
| mask    |       1 |        15.894 |          1.944 |          1.889 |            1.009 |            1.068 |        1.003 |        1.222 |
| mask    |       4 |        16.051 |          1.944 |          1.889 |            1.021 |            1.047 |        1.003 |        1.213 |
| no mask |       1 |        14.266 |          1.940 |          1.888 |            1.053 |            1.054 |        1.015 |        1.237 |
| no mask |       4 |        14.589 |          1.940 |          1.888 |            1.021 |            1.071 |        1.000 |        1.224 |

The cost of the mask on one thread, as the difference of the two series:

| scale | mask_records | seconds_with_mask | seconds_without_mask | mask_seconds | bases_left_out | rss_mib_with_mask | rss_mib_without_mask |
|------:|-------------:|------------------:|---------------------:|-------------:|---------------:|------------------:|---------------------:|
|     1 |    1,567,146 |            15.894 |               14.266 |        1.628 |     94,619,157 |              67.9 |                 66.3 |
|     2 |    3,074,594 |            31.074 |               28.650 |        2.424 |    175,290,807 |              68.1 |                 67.3 |
|     4 |    5,769,087 |            61.313 |               55.982 |        5.331 |    325,208,782 |              83.0 |                 82.0 |

The mask-records series, at fixed alignments. Its 109 windows
hold 114,294,784 reference
bases, and `mask_records` are the records of each mask that those windows
hold, out of 5,769,087 in the full
mask file. Scale 0 is the same regions without a mask:

| scale | threads | mask_records |   bases | mismatches | output_rows | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | budget_mib |
|------:|--------:|-------------:|--------:|-----------:|------------:|---------------:|------------:|------------:|---------------:|------------:|-----------:|
|     0 |       1 |            0 | 566,202 |      1,784 |       4,907 |          0.543 |       0.542 |       0.561 |           76.1 |        76.9 |      197.5 |
|     1 |       1 |      716,832 | 560,511 |      1,646 |       4,781 |          1.085 |       1.075 |       1.107 |           77.7 |        77.7 |      197.5 |
|     2 |       1 |    1,434,126 | 554,837 |      1,566 |       4,704 |          1.593 |       1.575 |       1.621 |           77.8 |        77.8 |      197.5 |
|     4 |       1 |    2,866,198 | 546,018 |      1,292 |       4,468 |          2.592 |       2.573 |       2.633 |           77.9 |        78.3 |      197.5 |
|     0 |       4 |            0 | 566,202 |      1,784 |       4,907 |          0.575 |       0.573 |       0.580 |           76.2 |        76.2 |      197.5 |
|     1 |       4 |      716,832 | 560,511 |      1,646 |       4,781 |          1.105 |       1.093 |       1.116 |           77.8 |        78.0 |      197.5 |
|     2 |       4 |    1,434,126 | 554,837 |      1,566 |       4,704 |          1.602 |       1.596 |       1.629 |           77.9 |        78.2 |      197.5 |
|     4 |       4 |    2,866,198 | 546,018 |      1,292 |       4,468 |          2.636 |       2.618 |       2.644 |           77.8 |        78.9 |      197.5 |

Its growth against the mask records, for the whole run and for the mask’s own cost (the run without a mask taken off), without a timing verdict:

| threads | no_mask_seconds | one_x_seconds | records_ratio_2x | records_ratio_4x | time_exponent_2x | time_exponent_4x | mask_cost_exponent_2x | mask_cost_exponent_4x | rss_ratio_4x |
|--------:|----------------:|--------------:|-----------------:|-----------------:|-----------------:|-----------------:|----------------------:|----------------------:|-------------:|
|       1 |           0.543 |         1.085 |            2.001 |            1.999 |            0.554 |            0.703 |                 0.954 |                 0.965 |        1.003 |
|       4 |           0.575 |         1.105 |            2.001 |            1.999 |            0.536 |            0.719 |                 0.953 |                 1.006 |        1.000 |

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
second) on top of the 0.54 s of the fixed
alignments, and peak RSS is 76.1 MiB without a mask
and 77.9 MiB with the full mask: the mask window is
1 MiB whatever the mask holds. The whole-run time exponents at one thread are
0.55 and
0.7, and
the mask’s own cost grows with exponents
0.95 and
0.97.
These times carry no verdict, because the 1× run takes
1.08 s.

With the mask on one thread, the whole contig takes
61.3 s for 7,053,095,954
counted bases, which is 115 million counted bases
a second, at 83 MiB peak RSS. The 1× region takes
15.9 s at 67.9 MiB.
Peak RSS at 4× is 1.22
times the 1× value, and the largest peak RSS of any run is
85.3 MiB against the ceiling of 197.5 MiB.

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
a whole genome in one call, BAM input, unsorted input or long reads.

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
