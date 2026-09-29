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

Each builder is an installed Rduckhts library verified by
`benchmark_ancestry_panel_identity.R` against the commit it was built
from: every tracked code file of `r/Rduckhts` at that commit is
byte-identical in the tree the library was built from, which holds
nothing else but build outputs, and every function in the installed
namespace equals the function that commit’s R sources define, and the
installed extension binary equals one rebuilt from that commit, apart
from R’s temporary install path embedded in it and the GNU build-id
derived from that path. The hashes of the installed package code and
extension then tie each measured run to its library. Every measured
process records the hashes of the build it loaded, and the report
requires them to match. The develop library was built from the
DuckVEP-removal branch’s merge of `develop`, whose panel builder is
byte-identical to `develop` at `d8475811`; the tag
`bench/ancestry-panel-baseline-develop` keeps that commit reachable. The
final builder’s source equals this tree’s.

| implementation | commit       | builder_blob | verified_code_files | verified_functions | extension_rebuilt_equal | package_code_sha256 | extension_sha256 |
|:---------------|:-------------|:-------------|:--------------------|:-------------------|:------------------------|:--------------------|:-----------------|
| develop        | af8dd4ac705c | 69ade4e3e9a2 | 1444                | 141                | TRUE                    | f0c752ff0a72        | c7876673b03f     |
| distinct       | 36c212ba6112 | 6e72fea378d0 | 1543                | 143                | TRUE                    | 36bd0e3dada4        | 4e97043b5372     |
| final          | 5090db2465dc | 9fdc2acf13d7 | 1543                | 143                | TRUE                    | 31427abdbb52        | 88d173d8be15     |

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
| develop        |       1 | candidates |        1 | 1452923 |            30511383 |              3630 | 3/3       |    23.63 |  22.81 |  25.00 |         2814 |            4752 |     0.000 |        818 |       1 |
| distinct       |       1 | candidates |        1 | 1452923 |            30511383 |              3630 | 3/3       |    24.02 |  23.60 |  32.61 |         2816 |            4922 |     0.000 |        818 |      22 |
| final          |       1 | candidates |        1 | 1452923 |            30511383 |              3630 | 3/3       |    14.33 |  12.79 |  15.00 |          596 |             705 |     0.000 |        818 |      22 |
| develop        |       1 | candidates |        2 | 2906994 |            61046874 |              7263 | 3/3       |    46.83 |  46.01 |  48.25 |         5384 |            9527 |     0.000 |       1525 |       1 |
| distinct       |       1 | candidates |        2 | 2906994 |            61046874 |              7263 | 3/3       |    46.88 |  46.67 |  52.68 |         5386 |            9554 |     0.000 |       1525 |      22 |
| final          |       1 | candidates |        2 | 2906994 |            61046874 |              7263 | 3/3       |    23.71 |  21.99 |  26.63 |          960 |            1340 |     0.000 |       1525 |      22 |
| develop        |       1 | candidates |        4 | 5816590 |           122148390 |             14533 | 3/3       |   120.02 | 113.16 | 134.22 |        16251 |           29818 |  2511.906 |       2940 |       1 |
| distinct       |       1 | candidates |        4 | 5816590 |           122148390 |             14533 | 3/3       |   125.99 | 116.96 | 126.69 |        16263 |           29817 |  2552.125 |       2940 |      22 |
| final          |       1 | candidates |        4 | 5816590 |           122148390 |             14533 | 3/3       |    42.80 |  40.02 |  45.11 |         1770 |            2605 |     0.000 |       2940 |      22 |
| develop        |       1 | spaced     |        1 | 1452923 |            30511383 |              3630 | 3/3       |    22.94 |  22.43 |  23.88 |         2802 |            4752 |     0.000 |        818 |       1 |
| distinct       |       1 | spaced     |        1 | 1452923 |            30511383 |              3630 | 3/3       |    23.09 |  22.80 |  23.94 |         2803 |            4767 |     0.000 |        818 |      22 |
| final          |       1 | spaced     |        1 | 1452923 |            30511383 |              3630 | 3/3       |    12.42 |  12.00 |  12.75 |          585 |             596 |     0.000 |        818 |      22 |
| develop        |       1 | spaced     |        2 | 2906994 |            61046874 |              7263 | 3/3       |    44.99 |  44.15 |  47.82 |         5372 |            9527 |     0.000 |       1525 |       1 |
| distinct       |       1 | spaced     |        2 | 2906994 |            61046874 |              7263 | 3/3       |    45.40 |  44.49 |  51.77 |         5374 |            9554 |     0.000 |       1525 |      22 |
| final          |       1 | spaced     |        2 | 2906994 |            61046874 |              7263 | 3/3       |    20.40 |  20.23 |  21.82 |          949 |            1003 |     0.000 |       1525 |      22 |
| develop        |       1 | spaced     |        4 | 5816590 |           122148390 |             14533 | 3/3       |   114.04 | 111.36 | 118.16 |        16241 |           29815 |  2481.688 |       2940 |       1 |
| distinct       |       1 | spaced     |        4 | 5816590 |           122148390 |             14533 | 3/3       |   114.58 | 113.23 | 121.23 |        16277 |           29814 |  2553.281 |       2940 |      22 |
| final          |       1 | spaced     |        4 | 5816590 |           122148390 |             14533 | 3/3       |    38.29 |  37.14 |  50.10 |         1759 |            1942 |     0.000 |       2940 |      22 |
| develop        |       4 | candidates |        1 | 1452923 |            30511383 |              3630 | 3/3       |     8.14 |   7.67 |   8.90 |         4516 |            8120 |     0.000 |        818 |       1 |
| distinct       |       4 | candidates |        1 | 1452923 |            30511383 |              3630 | 3/3       |     8.52 |   7.91 |  10.21 |         4518 |            8111 |     0.000 |        818 |      22 |
| final          |       4 | candidates |        1 | 1452923 |            30511383 |              3630 | 3/3       |     3.88 |   3.83 |   5.00 |          654 |             860 |     0.000 |        818 |      22 |
| develop        |       4 | candidates |        2 | 2906994 |            61046874 |              7263 | 3/3       |    14.55 |  13.67 |  17.94 |         8731 |           16055 |     0.000 |       1525 |       1 |
| distinct       |       4 | candidates |        2 | 2906994 |            61046874 |              7263 | 3/3       |    15.17 |  14.38 |  17.61 |         8808 |           15959 |     0.000 |       1525 |      22 |
| final          |       4 | candidates |        2 | 2906994 |            61046874 |              7263 | 3/3       |     6.67 |   6.30 |   8.32 |         1019 |            1491 |     0.000 |       1525 |      22 |
| develop        |       4 | candidates |        4 | 5816590 |           122148390 |             14533 | 0/3       |       NA |     NA |     NA |           NA |              NA |        NA |       2940 |      NA |
| distinct       |       4 | candidates |        4 | 5816590 |           122148390 |             14533 | 0/3       |       NA |     NA |     NA |           NA |              NA |        NA |       2940 |      NA |
| final          |       4 | candidates |        4 | 5816590 |           122148390 |             14533 | 3/3       |    12.02 |  11.32 |  14.28 |         1743 |            2741 |     0.000 |       2940 |      22 |
| develop        |       4 | spaced     |        1 | 1452923 |            30511383 |              3630 | 3/3       |     8.87 |   7.66 |  10.25 |         4518 |            8169 |     0.000 |        818 |       1 |
| distinct       |       4 | spaced     |        1 | 1452923 |            30511383 |              3630 | 3/3       |     9.11 |   7.69 |  11.99 |         4508 |            8270 |     0.000 |        818 |      22 |
| final          |       4 | spaced     |        1 | 1452923 |            30511383 |              3630 | 3/3       |     5.09 |   4.42 |   5.12 |          646 |             705 |     0.000 |        818 |      22 |
| develop        |       4 | spaced     |        2 | 2906994 |            61046874 |              7263 | 3/3       |    16.13 |  14.06 |  18.65 |         8748 |           16048 |     0.000 |       1525 |       1 |
| distinct       |       4 | spaced     |        2 | 2906994 |            61046874 |              7263 | 3/3       |    17.66 |  17.47 |  17.86 |         8772 |           16106 |     0.000 |       1525 |      22 |
| final          |       4 | spaced     |        2 | 2906994 |            61046874 |              7263 | 3/3       |     7.60 |   6.33 |   8.20 |         1007 |            1123 |     0.000 |       1525 |      22 |
| develop        |       4 | spaced     |        4 | 5816590 |           122148390 |             14533 | 0/3       |       NA |     NA |     NA |           NA |              NA |        NA |       2940 |      NA |
| distinct       |       4 | spaced     |        4 | 5816590 |           122148390 |             14533 | 0/3       |       NA |     NA |     NA |           NA |              NA |        NA |       2940 |      NA |
| final          |       4 | spaced     |        4 | 5816590 |           122148390 |             14533 | 3/3       |    12.49 |  10.35 |  13.53 |         1730 |            2068 |     0.000 |       2940 |      22 |

These runs stopped with DuckDB’s out-of-memory error at the 16GB limit:
develop at 4 threads, candidates, 4 quarters; distinct at 4 threads,
candidates, 4 quarters; develop at 4 threads, spaced, 4 quarters;
distinct at 4 threads, spaced, 4 quarters.

Scaling exponents of the final builder over 1, 2 and 4 quarters:

| threads | mode       | time_exponent | rss_exponent |
|--------:|:-----------|--------------:|-------------:|
|       1 | candidates |          0.79 |         0.79 |
|       4 | candidates |          0.82 |         0.71 |
|       1 | spaced     |          0.81 |         0.79 |
|       4 | spaced     |          0.65 |         0.71 |

## Reading

The final builder produces the same panel hash as the distinct builder
wherever the distinct builder completed, so the constant-state checks
select exactly the same loci. At the full reference with one thread it
takes 38.3 s and 1,759 MiB peak RSS, against 114.6 s and 16,277 MiB with
distinct aggregates, which hold a set of alleles, groups or PCs per
locus across the whole reference: they exceed the budget from the first
quarter, spill up to 2,553 MiB at full scale with one thread, and run
out of memory at four threads. Every final run is within its budget with
no spill; with four threads the full reference takes 12.5 s. Time and
memory grow no faster than the reference.

The per-contig allocation itself is cheap: the develop and distinct
builders differ only in allocation and take similar time and memory.
Every develop panel comes from chromosome 1; every final panel spans all
22 autosomes.
