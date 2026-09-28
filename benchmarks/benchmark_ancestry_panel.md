Ancestry panel building at scale
================

`rduckhts_ancestry_panel()` reads a long reference (loci × groups
frequency rows and loci × PCs loading rows), keeps the eligible loci,
and caps the panel at `max_sites`. This report measures three builders
on the same inputs:

- **develop** keeps the first `max_sites` eligible loci in genomic order
  and checks eligibility with distinct aggregates (one allele pair,
  every group and PC once).
- **distinct** allocates the cap per contig (one site each, the rest in
  proportion to each contig’s eligible sites minus one, spread along the
  contig) with the same distinct-aggregate eligibility checks.
- **final** keeps the per-contig allocation and checks eligibility with
  constant per-locus state: `min = max` for each allele, and a row count
  equal to the group (or PC) count whose ordinal bits cover the full
  mask, so each group and PC occurs exactly once.

## Method

`benchmark_ancestry_panel_run.R` stages bigsnpr’s reference through the
benchmark registry (`ancestry_reference_parquet`: 5,816,590 loci, 21
groups, 16 PCs) and keeps 1, 2 or 4 quarters of its loci by hash, the
one growth dimension (groups and PCs are fixed by the reference). Each
run expands them into the long relations the builder reads, builds a
committed 17,000-site panel through the public function, and records the
panel hash. `spaced` uses the default 5 kb spacing; `candidates` offers
every locus as a candidate, the largest input the allocation windows
see. Every run is a fresh process with `memory_limit = '16GB'` and its
own empty temporary directory, at one and four threads, three
repetitions interleaved across builders. The reference Parquet file is
read from a warm page cache. DuckDB peak buffer memory and spill are the
maxima of the per-statement profiles over every statement the builder
runs; peak RSS is the process high-water mark. Other jobs shared the
machine; builders were interleaved so they saw the same conditions.

The memory budget is fixed overhead (111 MiB: packages, connection and
extension load, measured by runs that build no panel) plus three times
the eligibility live state, 170 bytes per locus of group key and
constant aggregate state; spill must be zero. Decoded input is the
required long columns as DuckDB vectors.

## Results

| implementation | threads | mode       | quarters |    loci | long_reference_rows | decoded_input_mib | completed | median_s |  min_s |  max_s | peak_rss_mib | peak_buffer_mib | spill_mib | budget_mib | contigs |
|:---------------|--------:|:-----------|---------:|--------:|--------------------:|------------------:|:----------|---------:|-------:|-------:|-------------:|----------------:|----------:|-----------:|--------:|
| develop        |       1 | candidates |        1 | 1452923 |            30511383 |              3630 | 3/3       |    24.28 |  23.28 |  32.16 |         2813 |            4752 |     0.000 |        818 |       1 |
| distinct       |       1 | candidates |        1 | 1452923 |            30511383 |              3630 | 3/3       |    24.52 |  23.89 |  25.46 |         2816 |            4922 |     0.000 |        818 |      22 |
| final          |       1 | candidates |        1 | 1452923 |            30511383 |              3630 | 3/3       |    13.45 |  13.42 |  16.21 |          596 |             705 |     0.000 |        818 |      22 |
| develop        |       1 | candidates |        2 | 2906994 |            61046874 |              7263 | 3/3       |    48.41 |  46.17 |  61.87 |         5384 |            9542 |     0.000 |       1525 |       1 |
| distinct       |       1 | candidates |        2 | 2906994 |            61046874 |              7263 | 3/3       |    47.98 |  47.14 |  48.99 |         5387 |            9554 |     0.000 |       1525 |      22 |
| final          |       1 | candidates |        2 | 2906994 |            61046874 |              7263 | 3/3       |    22.33 |  22.02 |  25.07 |          958 |            1340 |     0.000 |       1525 |      22 |
| develop        |       1 | candidates |        4 | 5816590 |           122148390 |             14533 | 3/3       |   114.97 | 114.91 | 120.33 |        16223 |           29818 |  2555.312 |       2940 |       1 |
| distinct       |       1 | candidates |        4 | 5816590 |           122148390 |             14533 | 3/3       |   119.74 | 117.43 | 123.58 |        16251 |           29817 |  2522.562 |       2940 |      22 |
| final          |       1 | candidates |        4 | 5816590 |           122148390 |             14533 | 3/3       |    41.11 |  40.05 |  41.77 |         1770 |            2605 |     0.000 |       2940 |      22 |
| develop        |       1 | spaced     |        1 | 1452923 |            30511383 |              3630 | 3/3       |    24.19 |  23.37 |  25.72 |         2801 |            4752 |     0.000 |        818 |       1 |
| distinct       |       1 | spaced     |        1 | 1452923 |            30511383 |              3630 | 3/3       |    23.43 |  23.32 |  24.09 |         2804 |            4922 |     0.000 |        818 |      22 |
| final          |       1 | spaced     |        1 | 1452923 |            30511383 |              3630 | 3/3       |    12.64 |  12.15 |  12.94 |          585 |             596 |     0.000 |        818 |      22 |
| develop        |       1 | spaced     |        2 | 2906994 |            61046874 |              7263 | 3/3       |    49.01 |  47.39 |  90.64 |         5372 |            9527 |     0.000 |       1525 |       1 |
| distinct       |       1 | spaced     |        2 | 2906994 |            61046874 |              7263 | 3/3       |    48.91 |  47.13 |  72.64 |         5373 |            9553 |     0.000 |       1525 |      22 |
| final          |       1 | spaced     |        2 | 2906994 |            61046874 |              7263 | 3/3       |    29.01 |  20.98 |  29.79 |          948 |            1003 |     0.000 |       1525 |      22 |
| develop        |       1 | spaced     |        4 | 5816590 |           122148390 |             14533 | 3/3       |   117.83 | 113.43 | 151.03 |        15967 |           29815 |  2486.938 |       2940 |       1 |
| distinct       |       1 | spaced     |        4 | 5816590 |           122148390 |             14533 | 3/3       |   116.95 | 115.46 | 119.85 |        16275 |           29814 |  2551.844 |       2940 |      22 |
| final          |       1 | spaced     |        4 | 5816590 |           122148390 |             14533 | 3/3       |    37.81 |  36.88 |  41.28 |         1760 |            1942 |     0.000 |       2940 |      22 |
| develop        |       4 | candidates |        1 | 1452923 |            30511383 |              3630 | 3/3       |     7.91 |   7.70 |   8.49 |         4523 |            8111 |     0.000 |        818 |       1 |
| distinct       |       4 | candidates |        1 | 1452923 |            30511383 |              3630 | 3/3       |     7.96 |   7.88 |  10.42 |         4522 |            8120 |     0.000 |        818 |      22 |
| final          |       4 | candidates |        1 | 1452923 |            30511383 |              3630 | 3/3       |     3.80 |   3.79 |   4.13 |          659 |             857 |     0.000 |        818 |      22 |
| develop        |       4 | candidates |        2 | 2906994 |            61046874 |              7263 | 3/3       |    13.72 |  13.62 |  15.25 |         8743 |           16052 |     0.000 |       1525 |       1 |
| distinct       |       4 | candidates |        2 | 2906994 |            61046874 |              7263 | 3/3       |    15.70 |  15.11 |  16.02 |         8816 |           15974 |     0.000 |       1525 |      22 |
| final          |       4 | candidates |        2 | 2906994 |            61046874 |              7263 | 3/3       |     6.59 |   6.34 |   6.74 |         1030 |            1487 |     0.000 |       1525 |      22 |
| develop        |       4 | candidates |        4 | 5816590 |           122148390 |             14533 | 0/3       |       NA |     NA |     NA |           NA |              NA |        NA |       2940 |      NA |
| distinct       |       4 | candidates |        4 | 5816590 |           122148390 |             14533 | 0/3       |       NA |     NA |     NA |           NA |              NA |        NA |       2940 |      NA |
| final          |       4 | candidates |        4 | 5816590 |           122148390 |             14533 | 3/3       |    11.29 |  11.29 |  13.24 |         1739 |            2740 |     0.000 |       2940 |      22 |
| develop        |       4 | spaced     |        1 | 1452923 |            30511383 |              3630 | 3/3       |     7.88 |   7.59 |   9.29 |         4513 |            8135 |     0.000 |        818 |       1 |
| distinct       |       4 | spaced     |        1 | 1452923 |            30511383 |              3630 | 3/3       |     7.74 |   7.63 |  11.41 |         4505 |            8110 |     0.000 |        818 |      22 |
| final          |       4 | spaced     |        1 | 1452923 |            30511383 |              3630 | 3/3       |     3.64 |   3.63 |   3.79 |          644 |             712 |     0.000 |        818 |      22 |
| develop        |       4 | spaced     |        2 | 2906994 |            61046874 |              7263 | 3/3       |    14.24 |  14.00 |  14.68 |         8730 |           16055 |     0.000 |       1525 |       1 |
| distinct       |       4 | spaced     |        2 | 2906994 |            61046874 |              7263 | 3/3       |    14.98 |  13.88 |  15.24 |         8737 |           15961 |     0.000 |       1525 |      22 |
| final          |       4 | spaced     |        2 | 2906994 |            61046874 |              7263 | 3/3       |     6.21 |   5.86 |   6.32 |         1010 |            1121 |     0.000 |       1525 |      22 |
| develop        |       4 | spaced     |        4 | 5816590 |           122148390 |             14533 | 0/3       |       NA |     NA |     NA |           NA |              NA |        NA |       2940 |      NA |
| distinct       |       4 | spaced     |        4 | 5816590 |           122148390 |             14533 | 0/3       |       NA |     NA |     NA |           NA |              NA |        NA |       2940 |      NA |
| final          |       4 | spaced     |        4 | 5816590 |           122148390 |             14533 | 3/3       |    10.50 |  10.24 |  11.59 |         1730 |            2103 |     0.000 |       2940 |      22 |

These runs stopped with DuckDB’s out-of-memory error at the 16GB limit:
develop at 4 threads, candidates, 4 quarters; distinct at 4 threads,
candidates, 4 quarters; develop at 4 threads, spaced, 4 quarters;
distinct at 4 threads, spaced, 4 quarters.

Scaling exponents of the final builder over 1, 2 and 4 quarters:

| threads | mode       | time_exponent | rss_exponent |
|--------:|:-----------|--------------:|-------------:|
|       1 | candidates |          0.81 |         0.79 |
|       4 | candidates |          0.79 |         0.70 |
|       1 | spaced     |          0.79 |         0.79 |
|       4 | spaced     |          0.76 |         0.71 |

## Reading

The final builder produces the same panel hash as the distinct builder
wherever the distinct builder completed, so the constant-state checks
select exactly the same loci. At the full reference with one thread it
takes 37.8 s and 1,760 MiB peak RSS, against 117.0 s and 16,275 MiB with
distinct aggregates, which hold a set of alleles, groups or PCs per
locus across the whole reference: they exceed the budget from the first
quarter, spill up to 2,555 MiB at full scale with one thread, and run
out of memory at four threads. Every final run is within its budget with
no spill; with four threads the full reference takes 10.5 s. Time and
memory grow no faster than the reference.

The per-contig allocation itself is cheap: the develop and distinct
builders differ only in allocation and take similar time and memory.
Every develop panel comes from chromosome 1; every final panel spans all
22 autosomes.
