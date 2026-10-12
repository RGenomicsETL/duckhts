Mate overlap count scaling
================

This report gives the 1×/2×/4× scaling evidence `STYLE.md` asks for, for
`duckhts_bam_mate_overlap_counts()`. The function reads coordinate-sorted
alignments once, on one worker. It holds the first mate of a pair whose mate
starts inside it until the mate arrives, and it counts the positions that both
mates align to in fixed counters. So its work grows with the alignments it
reads, and its memory grows with the mates it holds at once. No earlier report
measures this function; this is its first baseline.

## Workloads

- **Real.** chr1 of the 1000 Genomes 30× CRAM of HG00403, staged as a CRAM
  slice (registry artifact `bam_mismatch_hg00403_chr1_cram`, staged by
  `benchmarks/bam_mismatch_stage.R`), read with the registered GRCh38 FASTA
  (`liftover_grch38_fasta`) as its reference. The render checks the slice
  against its registered identity and the FASTA against its registered size
  and MD5 before it measures. Those checks read both files in full, so every
  measured run reads from a warm page cache. Three regions double the
  alignments: chr1 up to the 25% and the 50% quantiles of the alignment start
  positions, and the whole contig. The mates held at once are what the data
  give. The render derives them, and they are not fixed.
- **Synthetic.** Generated at render time, untimed: DuckDB writes SAM text, and
  the samtools that the RBCFTools package bundles sorts it and indexes it as
  BAM. A contig holds fragment families 1,000 bases apart. A family is P pairs
  of 150-base reads whose read 1s all start at one position. Read 2 of pair j
  starts d = 1 + (97 j mod 149) bases after its read 1. So it shares 150 − d
  positions (1 to 149) with its read 1, and every read 1 is held until its
  read 2 arrives. All P read 1s of a family are held at once: the held-mate
  count of a contig is P, and the contig has 2 × families × P alignments. One
  pair in 16 has one substituted read-2 base, a known disagreement. One pair
  in 16 has a read 2 with MAPQ 10, which the default `min_mapq` of 20 leaves
  out; its read 1 is held and then dropped when the scan passes its mate’s
  start. Seven contigs give three series that share the 1× contig
  (families × P):
  - alignments at a fixed held-mate count: 4 × 100,000, 8 × 100,000 and
    16 × 100,000;
  - held mates at fixed alignments: 4 × 100,000, 2 × 200,000 and
    1 × 400,000;
  - both: 4 × 100,000, 4 × 200,000 and 4 × 400,000.

  A family stands for a deep amplicon, or a deep pileup of short fragments:
  the stress shape for held mates. The amplicon sequence and the base
  qualities are fixed, so this BAM compresses far better than real data.

## Held mates, derived from the input

The function holds a kept alignment when its mate starts inside it and no
mate of the same name is held. It releases the held mate when its mate
arrives, or when the scan passes the mate’s start without it. A record with
the name of a held mate that does not pair with it is left out, and the held
mate stays held. The render derives the number of mates held at once from
each input, with SQL over `read_bam()` and without the function:

- an alignment is kept by the function’s default filters: MAPQ at least 20, no
  flag of 3844, paired, mate mapped, and mate on the same contig;
- the first kept alignment of a name, by position, is held when its mate
  starts inside it (`PNEXT` from `POS` up to its reference end);
- the held alignment is released at its mate’s position when the name has
  two kept alignments that pair: one is read 1 and the other read 2, and each
  names the other’s start. Otherwise it is dropped before the first alignment
  past the position it names.

A drop comes before every hold at its position. A hold and a release by an
arriving mate at one position come in file order, which `read_bam()` does not
give for CRAM. So the derived count is bounded. The lower bound is the most
mates held after all the alignments of a position. The upper bound adds back
the mates that arriving mates release at that position, as if every hold
there came first. The render also runs the function untimed with
`max_pending` at the upper bound, which must succeed, and at one below the
lower bound, which must fail with its limit error. The budget uses the upper
bound. The derivation assumes that no name has more than two kept alignments
and that every kept alignment has a CIGAR that matches its sequence; the
render checks both.

## Budget, declared before measuring

The overhead is the median peak RSS of three empty runs that start the CLI,
load the extension and scan nothing.

The ceiling of each input is overhead + 3 × its live state. It is a gate on
the peak RSS: the render stops if any observation exceeds its ceiling or
spills. The live state is:

- the counters: 95 × 95 × 16 cells of 8 bytes, 1,155,200 bytes;
- the held mates, from the held-mate count H of the input (its upper bound)
  and D, the bytes of CIGAR, sequence and qualities of its longest kept
  alignment:
  - the entries, 320 bytes each and allocated 1,024 at a time, each with
    16 bytes in the heap and the free list: 1,024 × ⌈H / 1,024⌉ × 336 bytes;
  - the name table, 2 × C + 1 slots of 24 bytes, where C is the smallest
    1,024 × 2<sup>k</sup> that is at least H. While the table doubles, the
    half-size table it replaces is live too, so 24 × (3 C + 2) bytes, or
    24 × 2,049 bytes when C is 1,024;
  - the CIGAR, sequence and qualities of each held mate, H × (D + 24) bytes;
    the 24 bytes cover the allocator’s header and rounding. An entry is
    reused and keeps its largest buffer, so D is the longest kept alignment;
- htslib’s state for the file. For the CRAM, 32 MiB is allowed for the
  container in hand (its compressed and uncompressed blocks and its decoded
  records); this part is an allowance, not derived from the input. The
  reference of the slice in hand and of the slice before it is added to it.
  Without decoder threads, htslib asks for the reference of one slice’s span
  (`third_party/htslib/cram/cram_decode.c` lines 2452–2456).
  `cram_get_ref()` loads only that span, unless the span is at least half the
  contig (`third_party/htslib/cram/cram_io.c` lines 3248–3256). It frees the
  previous span after it loads the next one (lines 3340–3347). The render
  takes the widest pair of consecutive slices in each region from the CRAM
  index. For the BAM, 4 MiB is allowed for its blocks, the record and the
  index.

The render also stops if a time exponent of the real series is above 1.25 at
either doubling, or if its 1× run on one thread takes less than 5 s, the least
that `STYLE.md` accepts for a timing verdict. The synthetic series read their
1× contig in well under 5 s. So they get the memory verdict in full, and their
times are reported without a verdict, as `STYLE.md` provides.

## Protocol

Every observation is one fresh DuckDB CLI process under `/usr/bin/time`. The
JSON profiler gives the query latency and DuckDB’s peak buffer and temporary
bytes, and `/usr/bin/time` gives the peak RSS. The host is shared, so the
render is pinned to eight cores (8 to 15) and runs at scheduling priority −15
(`nice`), which every measured process inherits. The host load at the start
is recorded, and the spread of the repetitions shows what interference
remained. `max_temp_directory_size` is zero, so a spill fails the render;
DuckDB’s memory limit is left at its default, which the report records. Each
cell has three repetitions at one thread and at four threads, in an order
that alternates sizes. The function has one worker, so the thread setting
should not change it, and the four-thread cells show whether that holds. Each
query reduces the rows of the function to their number, the shared positions
and the disagreements, so every row is produced and consumed. The R wrapper
only builds this SQL call; it is not timed here.

The render also checks the counts. On each synthetic contig, every row of the
output must equal the counts that follow from the design, which R computes
from the amplicon, the qualities and the pair rules alone. On the real data,
every row of the output in the 400 kb around the position with the most held
mates must equal the same count made in SQL from `read_bam()` and
`cigar_aligned_blocks()`, as the SQL tests make it. That SQL count keeps only
alignments that can share a position with their mate: each mate must start
before the other ends.

``` sh
nice -n -15 taskset -c 8-15 Rscript -e 'rmarkdown::render("benchmarks/benchmark_bam_mate_overlap_counts.Rmd")'
```

The render needs RBCFTools 1.24-1.1.0 or later for its samtools. The
parameter `reduced` (`params = list(reduced = TRUE)`, with an `output_dir`
outside the repository) runs every step on less data: a 2 Mb window of the
CRAM around the most held mates (chr1:124,000,001–126,000,000) and synthetic
families of a hundredth of the pairs. It keeps every gate except the timing
verdict, writes its runs to a temporary file and says at its top that it is
not the report. It checks the code end to end in a few minutes.

    #> [1] "/tmp/claude-0/-root-duckhts/324b1460-ce1f-42c1-b0e3-6e2d29c7b777/scratchpad/rlib-rbcftools/RBCFTools/samtools/bin/samtools"

Source revision: 57c27de3446a8482d03eae179ce490294f06f4d4; `src` tree fda943b6dcbb090d387fc4a934091e6d17589e62. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). Extension SHA-256: ea551be30ca37c6e5a672b7cf7c2c48396cc762c3611c2c4a005b6a5ebc5dbdc. samtools 1.24 from RBCFTools 1.24.1.1.1.9000. Host load at start: 0.94, 1.46, 2.20. DuckDB memory limit: 50.0 GiB (the default); temporary directory limit: 0.

Staged input: the slice has 62,660,978 alignments (position sum 7728419844834032, 9,399,146,700 stored bases) in 6,270 CRAM slices. Its widest slice starts at 125,183,655 and spans 18,001,341 reference bases.

Overhead: 38.8 MiB.

## Inputs and budgets

Real regions. `alignments` are what the region selects, before the function’s
filters; `kept_alignments` pass them. `pairs_held` are the pairs whose first
mate is held, and `unpaired_names` the names with two kept alignments that do
not pair. `held_lower` and `held_upper` bound the mates held at once, and
`limit_agrees` is TRUE when the function succeeds with `max_pending` at the
upper bound and fails one below the lower bound. `peak_pos` is where the most
mates are held. `longest_bytes` is D, `reference_span` the bases of the widest
pair of consecutive slices, and the last columns the budget in MiB: the held
mates’ share, the live state, overhead + live state (`model_mib`) and the
ceiling.

| scale | region           | alignments | kept_alignments | pairs_held | unpaired_names | held_lower | held_upper | limit_agrees |    peak_pos | longest_bytes | reference_span | held_mib | live_mib | model_mib | ceiling_mib |
|------:|:-----------------|-----------:|----------------:|-----------:|---------------:|-----------:|-----------:|:-------------|------------:|--------------:|---------------:|---------:|---------:|----------:|------------:|
|     1 | chr1:1-61075319  | 15,665,245 |      13,456,359 |    391,357 |              0 |        115 |        119 | TRUE         |   2,652,264 |           293 |        175,413 |      0.4 |     33.7 |      72.5 |       139.8 |
|     2 | chr1:1-123524417 | 31,330,489 |      26,076,226 |    781,384 |              0 |        355 |        356 | TRUE         | 123,262,329 |           293 |        350,566 |      0.5 |     33.9 |      72.7 |       140.6 |
|     4 | chr1:1-248956422 | 62,660,978 |      49,228,489 |  1,509,814 |              0 |      4,189 |      4,240 | TRUE         | 125,180,069 |           293 |     18,001,706 |      3.5 |     53.8 |      92.6 |       200.1 |

Synthetic contigs, with the same columns. The derivation finds P mates held at
once on every contig, as the design says.

| cell | families | pairs | alignments | kept_alignments | pairs_held | unpaired_names | held_lower | held_upper | limit_agrees | longest_bytes | held_mib | live_mib | model_mib | ceiling_mib |
|:-----|---------:|------:|-----------:|----------------:|-----------:|---------------:|-----------:|-----------:|:-------------|--------------:|---------:|---------:|----------:|------------:|
| c11  |        4 | 1e+05 |    800,000 |         775,000 |    400,000 |              0 |      1e+05 |      1e+05 | TRUE         |           229 |     65.3 |     70.4 |     109.2 |       250.0 |
| c21  |        8 | 1e+05 |  1,600,000 |       1,550,000 |    800,000 |              0 |      1e+05 |      1e+05 | TRUE         |           229 |     65.3 |     70.4 |     109.2 |       250.0 |
| c41  |       16 | 1e+05 |  3,200,000 |       3,100,000 |  1,600,000 |              0 |      1e+05 |      1e+05 | TRUE         |           229 |     65.3 |     70.4 |     109.2 |       250.0 |
| c12  |        2 | 2e+05 |    800,000 |         775,000 |    400,000 |              0 |      2e+05 |      2e+05 | TRUE         |           229 |    130.6 |    135.7 |     174.5 |       445.8 |
| c14  |        1 | 4e+05 |    800,000 |         775,000 |    400,000 |              0 |      4e+05 |      4e+05 | TRUE         |           229 |    260.8 |    265.9 |     304.7 |       836.5 |
| c22  |        4 | 2e+05 |  1,600,000 |       1,550,000 |    800,000 |              0 |      2e+05 |      2e+05 | TRUE         |           229 |    130.6 |    135.7 |     174.5 |       445.8 |
| c44  |        4 | 4e+05 |  3,200,000 |       3,100,000 |  1,600,000 |              0 |      4e+05 |      4e+05 | TRUE         |           229 |    260.8 |    265.9 |     304.7 |       836.5 |

## Count checks

Each synthetic contig against the counts that follow from its design; the
measured runs reproduce the same totals:

| cell | expected_rows | observed_rows | expected_positions | expected_disagreements | rows_equal |
|:-----|--------------:|--------------:|-------------------:|-----------------------:|:-----------|
| c11  |            32 |            32 |         28,125,328 |                 25,000 | TRUE       |
| c21  |            32 |            32 |         56,250,656 |                 50,000 | TRUE       |
| c41  |            32 |            32 |        112,501,312 |                100,000 | TRUE       |
| c12  |            32 |            32 |         28,125,190 |                 25,000 | TRUE       |
| c14  |            32 |            32 |         28,125,063 |                 25,000 | TRUE       |
| c22  |            32 |            32 |         56,250,380 |                 50,000 | TRUE       |
| c44  |            32 |            32 |        112,500,252 |                100,000 | TRUE       |

The real data in chr1:124980069-125380069, the 400,001 bases
around the position with the most held mates: the function gives
279 rows and
the SQL count 279 rows of
921,511 shared positions, with
0 rows only in the function’s output and
0 only in the SQL count.

## Measurements

The real series. `positions`, `disagreements` and `output_rows` are the result
of the query. Times are the query latency of three repetitions.

| scale | threads | alignments | held_upper | output_rows |  positions | disagreements | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | max_buffer_mib | model_mib | ceiling_mib |
|------:|--------:|-----------:|-----------:|------------:|-----------:|--------------:|---------------:|------------:|------------:|---------------:|------------:|---------------:|----------:|------------:|
|     1 |       1 | 15,665,245 |        119 |         383 | 15,317,803 |        96,673 |          7.928 |       7.844 |       7.970 |           62.4 |        62.5 |            0.3 |      72.5 |       139.8 |
|     2 |       1 | 31,330,489 |        356 |         396 | 30,569,572 |       167,128 |         15.738 |      15.712 |      15.833 |           62.4 |        62.5 |            0.3 |      72.7 |       140.6 |
|     4 |       1 | 62,660,978 |      4,240 |         409 | 59,427,497 |       328,228 |         32.222 |      32.096 |      32.280 |           78.1 |        80.1 |            0.2 |      92.6 |       200.1 |
|     1 |       4 | 15,665,245 |        119 |         383 | 15,317,803 |        96,673 |          8.138 |       8.096 |       8.186 |           63.7 |        64.2 |            0.3 |      72.5 |       139.8 |
|     2 |       4 | 31,330,489 |        356 |         396 | 30,569,572 |       167,128 |         16.226 |      16.139 |      16.318 |           63.9 |        64.6 |            0.3 |      72.7 |       140.6 |
|     4 |       4 | 62,660,978 |      4,240 |         409 | 59,427,497 |       328,228 |         33.267 |      33.140 |      33.284 |           83.1 |        83.2 |            0.3 |      92.6 |       200.1 |

Growth of the median time against the alignments (an exponent of 1 is
linear), and peak-RSS ratios to 1×:

| threads | one_x_seconds | alignments_ratio_2x | alignments_ratio_4x | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x |
|--------:|--------------:|--------------------:|--------------------:|-----------------:|-----------------:|-------------:|-------------:|
|       1 |         7.928 |                   2 |                   2 |            0.989 |            1.034 |        1.000 |        1.252 |
|       4 |         8.138 |                   2 |                   2 |            0.995 |            1.036 |        1.003 |        1.305 |

The synthetic series:

| series     | scale | cell | threads | alignments | held_upper | output_rows |   positions | disagreements | median_seconds | min_seconds | max_seconds | median_rss_mib | max_rss_mib | max_buffer_mib | model_mib | ceiling_mib |
|:-----------|------:|:-----|--------:|-----------:|-----------:|------------:|------------:|--------------:|---------------:|------------:|------------:|---------------:|------------:|---------------:|----------:|------------:|
| alignments |     1 | c11  |       1 |    800,000 |      1e+05 |          32 |  28,125,328 |        25,000 |          0.342 |       0.341 |       0.343 |          100.9 |       101.1 |            0.3 |     109.2 |       250.0 |
| alignments |     2 | c21  |       1 |  1,600,000 |      1e+05 |          32 |  56,250,656 |        50,000 |          0.663 |       0.652 |       0.665 |          101.0 |       101.1 |            0.3 |     109.2 |       250.0 |
| alignments |     4 | c41  |       1 |  3,200,000 |      1e+05 |          32 | 112,501,312 |       100,000 |          1.302 |       1.299 |       1.304 |          100.8 |       101.1 |            0.3 |     109.2 |       250.0 |
| alignments |     1 | c11  |       4 |    800,000 |      1e+05 |          32 |  28,125,328 |        25,000 |          0.350 |       0.346 |       0.364 |          100.4 |       101.0 |            0.3 |     109.2 |       250.0 |
| alignments |     2 | c21  |       4 |  1,600,000 |      1e+05 |          32 |  56,250,656 |        50,000 |          0.675 |       0.658 |       0.683 |          100.7 |       101.0 |            0.3 |     109.2 |       250.0 |
| alignments |     4 | c41  |       4 |  3,200,000 |      1e+05 |          32 | 112,501,312 |       100,000 |          1.324 |       1.298 |       1.332 |          100.5 |       101.2 |            0.3 |     109.2 |       250.0 |
| held mates |     1 | c11  |       1 |    800,000 |      1e+05 |          32 |  28,125,328 |        25,000 |          0.342 |       0.341 |       0.343 |          100.9 |       101.1 |            0.3 |     109.2 |       250.0 |
| held mates |     2 | c12  |       1 |    800,000 |      2e+05 |          32 |  28,125,190 |        25,000 |          0.404 |       0.401 |       0.410 |          163.2 |       163.3 |            0.3 |     174.5 |       445.8 |
| held mates |     4 | c14  |       1 |    800,000 |      4e+05 |          32 |  28,125,063 |        25,000 |          0.495 |       0.493 |       0.497 |          292.2 |       292.4 |            0.3 |     304.7 |       836.5 |
| held mates |     1 | c11  |       4 |    800,000 |      1e+05 |          32 |  28,125,328 |        25,000 |          0.350 |       0.346 |       0.364 |          100.4 |       101.0 |            0.3 |     109.2 |       250.0 |
| held mates |     2 | c12  |       4 |    800,000 |      2e+05 |          32 |  28,125,190 |        25,000 |          0.425 |       0.418 |       0.432 |          164.9 |       165.0 |            0.3 |     174.5 |       445.8 |
| held mates |     4 | c14  |       4 |    800,000 |      4e+05 |          32 |  28,125,063 |        25,000 |          0.553 |       0.536 |       0.570 |          290.4 |       290.8 |            0.3 |     304.7 |       836.5 |
| both       |     1 | c11  |       1 |    800,000 |      1e+05 |          32 |  28,125,328 |        25,000 |          0.342 |       0.341 |       0.343 |          100.9 |       101.1 |            0.3 |     109.2 |       250.0 |
| both       |     2 | c22  |       1 |  1,600,000 |      2e+05 |          32 |  56,250,380 |        50,000 |          0.766 |       0.765 |       0.768 |          163.2 |       163.4 |            0.3 |     174.5 |       445.8 |
| both       |     4 | c44  |       1 |  3,200,000 |      4e+05 |          32 | 112,500,252 |       100,000 |          1.684 |       1.678 |       1.704 |          292.4 |       292.5 |            0.3 |     304.7 |       836.5 |
| both       |     1 | c11  |       4 |    800,000 |      1e+05 |          32 |  28,125,328 |        25,000 |          0.350 |       0.346 |       0.364 |          100.4 |       101.0 |            0.3 |     109.2 |       250.0 |
| both       |     2 | c22  |       4 |  1,600,000 |      2e+05 |          32 |  56,250,380 |        50,000 |          0.786 |       0.784 |       0.797 |          164.9 |       165.0 |            0.3 |     174.5 |       445.8 |
| both       |     4 | c44  |       4 |  3,200,000 |      4e+05 |          32 | 112,500,252 |       100,000 |          1.733 |       1.731 |       1.757 |          290.3 |       290.3 |            0.3 |     304.7 |       836.5 |

Their time ratios between doublings, reported without a verdict, and their
median peak RSS:

| series     | threads | one_x_seconds | time_ratio_2x | time_ratio_4x | rss_mib_1x | rss_mib_2x | rss_mib_4x |
|:-----------|--------:|--------------:|--------------:|--------------:|-----------:|-----------:|-----------:|
| alignments |       1 |         0.342 |         1.940 |         1.964 |      100.9 |      101.0 |      100.8 |
| alignments |       4 |         0.350 |         1.928 |         1.961 |      100.4 |      100.7 |      100.5 |
| held mates |       1 |         0.342 |         1.182 |         1.227 |      100.9 |      163.2 |      292.2 |
| held mates |       4 |         0.350 |         1.215 |         1.300 |      100.4 |      164.9 |      290.4 |
| both       |       1 |         0.342 |         2.243 |         2.197 |      100.9 |      163.2 |      292.4 |
| both       |       4 |         0.350 |         2.244 |         2.206 |      100.4 |      164.9 |      290.3 |

## Findings

Every observation is within its ceiling and none spills.
The 1× run of the real series on one thread takes at least 5 s, and none of its time exponents is above 1.25.
The render enforces all of this. It also stops unless a cell gives the same
rows, positions and disagreements in every repetition and at both thread
settings, unless the synthetic counts equal the design’s, and unless the real
counts around the most held mates equal the SQL count.

Real data: the 4× region, chr1:1-248956422, takes
32.2 s on one thread for
62,660,978 alignments, which is
1.94 million
alignments a second, and compares 59,427,497 shared
positions at 78.1 MiB peak RSS. The 1× region takes
7.9 s at
62.4 MiB. The time exponents at one thread are
0.99 and
1.03.

The held-mate count is not fixed in this series: it is at most
119, 356 and
4,240 in the three regions. The most are held at
125,180,069, past the end of
the 2× region, where the 400,001 bases of the count
check hold 1,500,739 alignments, a depth of about
563×. At 4× the model charges
3.5 MiB for the held mates, and
17.2 MiB for the reference of the
widest pair of consecutive slices. The widest slice spans
18,001,341 reference bases from
125,183,655; of the three regions, only the 4× region reads it.
On one thread, the median peak RSS rises by
15.7 MiB
from 2× to 4×.

Synthetic data: at a fixed held-mate count of
100,000, four times the alignments take the median
peak RSS on one thread from
100.9 to
100.8 MiB.
At fixed alignments, four times the held mates take it from
100.9 to
292.2 MiB.
Each further held mate costs 669 bytes of peak
RSS, and the model charges 683 bytes.
No observation, real or synthetic, has a peak RSS above overhead + live state.
The synthetic runs take 0.34 to
1.73 s, so their times carry no
verdict.

The operating scale this report supports is what it measures: one contig of
one 30× short-read sample from CRAM, and up to
400,000 mates held at once from BAM, under the default
`max_pending` of 1,000,000. It does not measure a whole genome in one call,
long reads, or more held mates than that. The synthetic times are too short to
show whether the time per alignment grows with the held mates.

## What the counts say

The 4× region by the lower of the two base qualities of a shared position,
from one untimed run:

| lower_quality |  positions | disagreements | disagreement_rate | share_of_positions |
|--------------:|-----------:|--------------:|------------------:|-------------------:|
|             3 |      2,379 |           760 |          0.319462 |             0.0000 |
|             4 |     26,191 |         7,950 |          0.303539 |             0.0004 |
|             5 |     25,133 |         7,510 |          0.298810 |             0.0004 |
|             6 |     77,367 |        17,865 |          0.230912 |             0.0013 |
|            10 |  2,071,871 |       232,150 |          0.112048 |             0.0349 |
|            20 |    942,913 |        20,908 |          0.022174 |             0.0159 |
|            30 | 56,281,643 |        41,085 |          0.000730 |             0.9471 |

A real allele sits on the molecule and both mates show it; a sequencing error
usually shows in one mate. So the disagreement rate is about twice the
instrument error rate. Errors made before sequencing, in the library, by PCR
or by DNA damage, show in both mates and are not seen. Only fragments shorter
than the two reads overlap, so the count leans toward short fragments and read
ends.
