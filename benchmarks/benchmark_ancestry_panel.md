Ancestry panel building: per-contig allocation at scale
================

`rduckhts_ancestry_panel()` caps a panel at `max_sites`. The previous
builder kept the first `max_sites` candidates in genomic order, so a
17,000-site panel over bigsnpr’s reference came entirely from
chromosome 1. The current builder gives each contig one site, shares the
rest in proportion to each contig’s eligible sites minus that one, and
spreads each contig’s sites evenly along it. That adds per-contig row
numbers and counts over every candidate, where the old builder kept a
bounded top-N. This report measures what that costs.

## Method

`benchmark_ancestry_panel_run.R` stages bigsnpr’s reference and loadings
(`ancestry_reference_parquet`, 5,816,590 loci, 21 groups, 16 PCs) and
keeps 1, 2 or 4 quarters of its loci by hash. Each run expands them into
the long frequency and loading relations the builder reads, builds a
committed 17,000-site panel, and counts it. Every run is a fresh
single-thread process; peak RSS comes from GNU time. There are three
repetitions of each implementation, mode and scale, in interleaved
order. `spaced` uses the default 5 kb window spacing; `candidates`
offers every locus as a candidate, the largest input the allocation
windows can see. `old` is the builder at the commit before this change
and `new` is the per-contig builder.

## Results

Medians of three runs:

| implementation | mode       | quarters | sites | contigs | seconds | peak_rss_gib |
|:---------------|:-----------|---------:|------:|--------:|--------:|-------------:|
| old            | candidates |        1 | 17000 |       1 |   23.09 |         2.75 |
| new            | candidates |        1 | 17000 |      22 |   23.60 |         2.75 |
| old            | candidates |        2 | 17000 |       1 |   41.44 |         5.26 |
| new            | candidates |        2 | 17000 |      22 |   42.75 |         5.26 |
| old            | candidates |        4 | 17000 |       1 |   77.11 |        10.35 |
| new            | candidates |        4 | 17000 |      22 |   81.97 |        10.35 |
| old            | spaced     |        1 | 17000 |       1 |   22.73 |         2.73 |
| new            | spaced     |        1 | 17000 |      22 |   22.97 |         2.74 |
| old            | spaced     |        2 | 17000 |       1 |   40.45 |         5.24 |
| new            | spaced     |        2 | 17000 |      22 |   40.81 |         5.25 |
| old            | spaced     |        4 | 17000 |       1 |   76.82 |        10.34 |
| new            | spaced     |        4 | 17000 |      22 |   77.58 |        10.34 |

The new builder against the old one at each scale:

| mode       | quarters | time_ratio | rss_delta_mib |
|:-----------|---------:|-----------:|--------------:|
| candidates |        1 |      1.022 |           2.0 |
| candidates |        2 |      1.032 |           1.9 |
| candidates |        4 |      1.063 |           3.0 |
| spaced     |        1 |      1.011 |           2.0 |
| spaced     |        2 |      1.009 |           1.7 |
| spaced     |        4 |      1.010 |           2.3 |

Scaling exponents of the new builder over 1, 2 and 4 quarters:

| mode       | time_exponent | memory_exponent |
|:-----------|--------------:|----------------:|
| candidates |          0.90 |            0.96 |
| spaced     |          0.88 |            0.96 |

## Reading

The allocation adds at most 6.3% to the median time and 3 MiB to peak
memory at any scale, and every new panel spans all 22 autosomes. Time
grows linearly with the reference in both builders.

Peak memory also grows linearly in both builders, to about 10.4 GiB at
full scale. That growth predates this change: the builder groups the
whole long reference (loci × groups) to find eligible sites before any
capping happens. Reducing it belongs to the builder’s site-eligibility
stage and is tracked in
[\#293](https://github.com/RGenomicsETL/duckhts/issues/293).
