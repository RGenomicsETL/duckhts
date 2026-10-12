read_bam CRAM projection
================

This report measures commit `2c36deb1` (“Ask htslib only for the CRAM fields
that read_bam projects”). `read_bam()` on a CRAM now sets htslib’s
`CRAM_OPT_REQUIRED_FIELDS` from the projected columns
(`bam_cram_required_fields()` in `src/bam_reader.c`). It always asks for QNAME
to TLEN. It adds SEQ and QUAL only when SEQ, QUAL, MD, NM or AUXILIARY_TAGS is
projected. It adds AUX and RGAUX only when a tag column, READ_GROUP_ID or
SAMPLE_ID is projected. Before the change, htslib decoded every field of every
record. The report compares the build of this revision with a baseline build
that lacks the change, on the same input, in the same render. It gives the
before/after ratio of time and peak RSS for each projection and size, and the
1×/2×/4× scaling evidence that `STYLE.md` asks for. No earlier report measures
`read_bam()` projections on a CRAM; this is their first baseline.

## Workload

The input is chr1 of the 1000 Genomes 30× CRAM of HG00403, staged as a CRAM
slice (registry artifact `bam_mismatch_hg00403_chr1_cram`, staged by
`benchmarks/bam_mismatch_stage.R`, 150-base reads). It is read with the
registered GRCh38 FASTA (`liftover_grch38_fasta`) as its reference. The render
checks the FASTA against its registered size and MD5. It checks the CRAM
against its registered content: the alignments, the sum of their positions
and the sum of their bases. Those checks read both files in full, so every
measured run reads from a warm page cache.

Three regions double the alignments read: chr1 up to the 25% and the 50%
quantiles of the alignment start positions, and the whole contig. The cuts
are the cuts of `benchmark_convert_parquet_scaling.Rmd`: `quantile_disc(POS, f)`
over every alignment, and the region `chr1:1-<cut>`.

Four queries select four projections. Each reduces its rows to one row, so
every row is decoded and consumed:

- **fixed**: `count(*)`, `sum(POS)` and `sum(MAPQ)`. The branch build asks for
  the fixed fields only.
- **sequence**: `count(*)`, `sum(len(SEQ))` and `sum(len(QUAL))`. The branch
  build adds SEQ and QUAL.
- **read_group**: `count(*)`, `count(DISTINCT READ_GROUP_ID)` and
  `sum(hash(READ_GROUP_ID))`. The branch build adds AUX and RGAUX.
- **all**: `count(*)` and `sum(hash(COLUMNS(*)))`, one hash sum per column of
  `read_bam()` with its defaults. Both builds ask for every field, so this
  projection is the control: its ratios should be near 1.

Each query runs with `decompression_threads := 2`, the default, and with
`decompression_threads := 0`. DuckDB runs one thread (`SET threads = 1`).

## Builds

- **branch**: the extension built from the source revision of the render,
  `build/release/duckhts.duckdb_extension` or `DUCKHTS_EXTENSION`.
- **baseline**: the extension in `CRAM_PROJECTION_BASE_EXTENSION`, built from
  the revision in `CRAM_PROJECTION_BASE_REVISION`. That revision must be an
  ancestor of the branch revision, and its `src/bam_reader.c` must lack
  `bam_cram_required_fields`. DuckDB finds an extension’s entry point from its
  file name, so the render copies the baseline file to a temporary directory
  as `duckhts.duckdb_extension`. The report lists every commit between the two
  revisions that changes `src/bam_reader.c`. Only the CRAM change touches the
  columns these queries project.

## Protocol

Every observation is one fresh DuckDB CLI process under `/usr/bin/time`, which
gives the peak RSS. The JSON profiler gives the latency of the query, its CPU
time, and DuckDB’s peak buffer and temporary bytes. `max_temp_directory_size`
is zero, so a spill fails the render. DuckDB’s memory limit is left at its
default, which the report records. Each cell (build, projection, decompression
threads, size) has three repetitions. Each repetition runs both builds for
each projection and thread setting. The build order and the size order
alternate between repetitions. The render command pins the processes to eight
cores and runs them at scheduling priority −15. The host load at the start is
recorded.

The render stops unless:

- every cell returns the same result in both builds, at both thread settings
  and in every repetition;
- every projection returns the same alignments at a size, and the whole
  contig returns the registered alignments, position sum (fixed) and base sum
  (sequence);
- no observation spills;
- every observation of the branch build is within its ceiling below;
- every time exponent of a branch series is at most 1.25 at either doubling,
  for each series whose 1× run takes at least 5 s on one thread. That is the
  least `STYLE.md` accepts for a timing verdict. A series under 5 s gets the
  memory verdict alone.

The baseline build gets no ceiling and no timing gate. It is the before
state, and its numbers are reported.

A branch series is **flat in memory** when its median peak RSS at 4× is at
most 10% above its median peak RSS at 1×. This is a finding, not a gate.

## Budget, declared before measuring

The overhead is the median peak RSS of three empty runs that start the CLI,
load the branch build and run `SELECT 1`. The ceiling of a branch observation
is overhead + 3 × the live state of its cell. The live state follows htslib
1.24 (`third_party/htslib`), with T the decompression threads.

**Per-record bytes.** For projection p, a record takes m<sub>p</sub> bytes:
36 for the BAM core, 4 × ⌈(name length + 1) / 4⌉ for the name, 4 for each
CIGAR operation, plus 2 × bases when p decodes SEQ and QUAL, plus
Σ (3 + max(value length + 1, 4)) over the tags when p decodes tags. This
bounds the BAM record and the slice’s decoded blocks. A base takes half a byte
in BAM and one byte in the slice’s sequence block, and a quality takes one
byte in both. A tag takes 2 bytes of key and 1 of type. A string value takes
its length and a NUL. An integer, float or character value takes at most 4
bytes. An array tag would break the bound, so the render stops if any tag
value has the array form. The tags are those of a scan that decodes every
field, so MD and NM, which htslib makes from the reference, are counted for
every projection that decodes tags. D<sub>p</sub> is the largest
m<sub>p</sub> on the contig. B<sub>p</sub> is the sum of its N largest
m<sub>p</sub>, where N = 10,000 records per slice (`SEQS_PER_SLICE`,
`cram/cram_structs.h:87`). The render checks that the CRAM index has at least
alignments / N slices. A<sub>p</sub> is D<sub>p</sub> rounded up to a power of
two, plus 32 bytes: the buffer that `sam_realloc_bam_data()` allocates
(`sam.c:438-444`).

**Slices in flight.** S = 2T + 1 slices with threads, 1 without. With T
threads, htslib’s decode queue has a size of 2T (`cram/cram_io.c:5581`). A
worker starts a job only while the decoded slices and the slices being decoded
number fewer than that (`thread_pool.c:537`, `579` and `681`). The slice in
hand adds one. Each slice in flight holds:

- its compressed container, at most C bytes, the largest container in the
  CRAM index;
- N decoded records of 152 bytes (`sizeof(cram_record)` on x86-64,
  `cram/cram_structs.h:473-513`);
- its decoded name, CIGAR, sequence, quality and tag blocks
  (`cram/cram_decode.c:3012-3016`), at most B<sub>p</sub>.

Each of max(T, 1) decoding workers also holds the slice’s uncompressed blocks
until the slice is decoded (`cram/cram_decode.c:2999-3006`), at most
B<sub>p</sub> more. The base-difference features of a slice, 16 bytes each,
are not counted; the factor of 3 must cover them.

**Record buffers.** With threads, each slice in flight converts its records
into a recycled list of N BAM records (`cram/cram_decode.c:2303-2351`), and the
slice in hand returns its list for reuse (`cram/cram_decode.c:3596-3600`). A
record’s buffer only grows, so each slot keeps the largest record it has held.
The bound is S × N × A<sub>p</sub>. It is reached only when every slot has
held a record as large as the largest record; how close the data come grows
with the slices read. Without threads, htslib converts one record at a time
into `read_bam()`’s own record: A<sub>p</sub>.

**Reference.** Without SEQ and QUAL, htslib loads no reference for a slice
that covers one reference (`cram/cram_decode.c:2452-2458`). With SEQ and QUAL
and threads, htslib shares one copy of the reference between its threads
(`cram/cram_io.c:5581-5582`) and loads the whole contig (`cram/cram_io.c:3253`):
248,956,422 bytes for chr1. Without threads, it loads the span of one slice
(`cram/cram_decode.c:2452-2456`), unless the span is at least half the contig
(`cram/cram_io.c:3248-3256`). It frees the previous span after it loads the
next one (`cram/cram_io.c:3340-3347`). The render takes the widest pair of
consecutive slices that each region reads from the CRAM index.

**DuckDB and the rest.** One output chunk of 2,048 rows, 2,048 × (16 × 14 +
D<sub>p</sub>) bytes, and an allowance of 8 MiB for the header, the CRAM index
and the file buffers. The allowance is not derived from the input.

The live state of a cell is the sum: S × (C + 152 N + B<sub>p</sub>) +
max(T, 1) × B<sub>p</sub> + the record buffers + the reference + the chunk +
8 MiB. The baseline build decodes every field, so the live state of the
**all** projection is also the model of every baseline cell; the report shows
it, but does not enforce it.

## Render

``` sh
CRAM_PROJECTION_BASE_EXTENSION=/path/to/baseline/duckhts.duckdb_extension \
CRAM_PROJECTION_BASE_REVISION=<baseline revision> \
nice -n -15 taskset -c 8-15 Rscript -e 'rmarkdown::render("benchmarks/benchmark_read_bam_cram_projection.Rmd")'
```

`CRAM_PROJECTION_SIZE=reduced` runs the same code on small regions, to check
the report end to end: chr1 up to the 0.2%, 0.4% and 0.8% quantiles, with one
repetition. The identity checks and the quantile pass still read the whole
contig. That run writes its observations to a temporary directory, so pass
`output_dir` to keep its report out of `benchmarks/` too; its numbers are not
evidence.

``` sh
CRAM_PROJECTION_SIZE=reduced CRAM_PROJECTION_BASE_EXTENSION=... CRAM_PROJECTION_BASE_REVISION=... \
taskset -c 0-3 Rscript -e 'rmarkdown::render("benchmarks/benchmark_read_bam_cram_projection.Rmd", output_dir = "/tmp")'
```

Source revision: 188934a282ca20a4a7883d93a4b853efb8952f4b; `src` tree 9ada6d903f188bc52198290590ce9b0599c765ca. Baseline revision: a5ab5147abb7abe398583e26ded5183ae7ecd82f. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). htslib 1.24. Host load at start: 1.12, 1.20, 1.31. CPU affinity: pid 258774’s current affinity list: 8-15. DuckDB memory limit: 50.0 GiB (the default); temporary directory limit: 0. Every observation is in [`benchmark_read_bam_cram_projection_runs.csv`](benchmark_read_bam_cram_projection_runs.csv).

| build    | sha256                                                           | empty_run_mib |
|:---------|:-----------------------------------------------------------------|--------------:|
| branch   | 0682f5dded2dda4a5a3a720e615769d250029ad1a02449cb5d192ae959e61ce5 |          38.8 |
| baseline | 274dee1c390b25474a251d89508badacc15414a9d4d9695700f47dbd96a2a279 |          38.5 |

Commits between the baseline and the branch revision that change `src/bam_reader.c`:

    2c36deb1 Ask htslib only for the CRAM fields that read_bam projects
    b22f4385 Write attribute MAP rows through one shared helper

## Input

The CRAM slice has 62,660,978 alignments (position sum 7728419844834032, 9,399,146,700 bases) in 6,270 slices. Its largest container is 820,458 bytes. No tag value has the array form.

| scale | region           | alignments | widest_slice_pair |
|------:|:-----------------|-----------:|------------------:|
|     1 | chr1:1-61075319  | 15,665,245 |           175,413 |
|     2 | chr1:1-123524417 | 31,330,489 |           350,566 |
|     4 | chr1             | 62,660,978 |        18,001,706 |

Per-record bytes of each projection, from the pass that decodes every field: `largest_bytes` is D<sub>p</sub>, `top_bytes` B<sub>p</sub>, the sum of the 10,000 largest, and `record_buffer_bytes` A<sub>p</sub>.

| projection | largest_bytes |  top_bytes | record_buffer_bytes |
|:-----------|--------------:|-----------:|--------------------:|
| fixed      |           152 |  1,098,672 |                 288 |
| sequence   |           452 |  4,098,672 |                 544 |
| read_group |         7,070 | 48,756,078 |               8,224 |
| all        |         7,370 | 51,756,078 |               8,224 |

The live state and the ceiling of each branch cell, in MiB:

| projection | threads | scale | slices_in_flight | slices_mib | record_buffers_mib | reference_mib | live_state_mib | budget_mib |
|:-----------|--------:|------:|-----------------:|-----------:|-------------------:|--------------:|---------------:|-----------:|
| fixed      |       2 |     1 |                5 |       18.5 |               13.7 |           0.0 |           41.0 |      161.7 |
| fixed      |       2 |     2 |                5 |       18.5 |               13.7 |           0.0 |           41.0 |      161.7 |
| fixed      |       2 |     4 |                5 |       18.5 |               13.7 |           0.0 |           41.0 |      161.7 |
| fixed      |       0 |     1 |                1 |        4.3 |                0.0 |           0.0 |           13.1 |       78.0 |
| fixed      |       0 |     2 |                1 |        4.3 |                0.0 |           0.0 |           13.1 |       78.0 |
| fixed      |       0 |     4 |                1 |        4.3 |                0.0 |           0.0 |           13.1 |       78.0 |
| sequence   |       2 |     1 |                5 |       38.5 |               25.9 |         237.4 |          311.2 |      972.4 |
| sequence   |       2 |     2 |                5 |       38.5 |               25.9 |         237.4 |          311.2 |      972.4 |
| sequence   |       2 |     4 |                5 |       38.5 |               25.9 |         237.4 |          311.2 |      972.4 |
| sequence   |       0 |     1 |                1 |       10.0 |                0.0 |           0.2 |           19.5 |       97.4 |
| sequence   |       0 |     2 |                1 |       10.0 |                0.0 |           0.3 |           19.7 |       97.9 |
| sequence   |       0 |     4 |                1 |       10.0 |                0.0 |          17.2 |           36.5 |      148.4 |
| read_group |       2 |     1 |                5 |      336.6 |              392.2 |           0.0 |          751.0 |     2291.9 |
| read_group |       2 |     2 |                5 |      336.6 |              392.2 |           0.0 |          751.0 |     2291.9 |
| read_group |       2 |     4 |                5 |      336.6 |              392.2 |           0.0 |          751.0 |     2291.9 |
| read_group |       0 |     1 |                1 |       95.2 |                0.0 |           0.0 |          117.5 |      391.3 |
| read_group |       0 |     2 |                1 |       95.2 |                0.0 |           0.0 |          117.5 |      391.3 |
| read_group |       0 |     4 |                1 |       95.2 |                0.0 |           0.0 |          117.5 |      391.3 |
| all        |       2 |     1 |                5 |      356.7 |              392.2 |         237.4 |         1009.1 |     3066.0 |
| all        |       2 |     2 |                5 |      356.7 |              392.2 |         237.4 |         1009.1 |     3066.0 |
| all        |       2 |     4 |                5 |      356.7 |              392.2 |         237.4 |         1009.1 |     3066.0 |
| all        |       0 |     1 |                1 |      100.9 |                0.0 |           0.2 |          124.0 |      410.7 |
| all        |       0 |     2 |                1 |      100.9 |                0.0 |           0.3 |          124.1 |      411.2 |
| all        |       0 |     4 |                1 |      100.9 |                0.0 |          17.2 |          141.0 |      461.7 |

## Measurements

`seconds` are the latency of the query over the repetitions. `cpu_per_second` is DuckDB’s CPU time per second of latency; htslib’s decompression threads are not in it. `budget_mib` is the ceiling of a branch cell.

| build    | projection | threads | scale | alignments | median_seconds | min_seconds | max_seconds | cpu_per_second | median_rss_mib | max_rss_mib | max_buffer_mib | budget_mib |
|:---------|:-----------|--------:|------:|-----------:|---------------:|------------:|------------:|---------------:|---------------:|------------:|---------------:|-----------:|
| baseline | fixed      |       2 |     1 | 15,665,245 |          4.524 |       4.476 |       4.662 |           0.99 |          415.6 |       422.2 |            0.1 |         NA |
| baseline | fixed      |       2 |     2 | 31,330,489 |          8.975 |       8.921 |       9.155 |           0.99 |          477.7 |       487.6 |            0.1 |         NA |
| baseline | fixed      |       2 |     4 | 62,660,978 |         18.013 |      17.703 |      18.318 |           1.00 |          542.5 |       552.3 |            0.1 |         NA |
| branch   | fixed      |       2 |     1 | 15,665,245 |          1.865 |       1.860 |       2.086 |           0.99 |           75.1 |        75.8 |            0.2 |      161.7 |
| branch   | fixed      |       2 |     2 | 31,330,489 |          3.749 |       3.624 |       4.048 |           0.99 |           77.1 |        80.8 |            0.1 |      161.7 |
| branch   | fixed      |       2 |     4 | 62,660,978 |          7.720 |       7.487 |       9.307 |           1.00 |           74.3 |        75.8 |            0.2 |      161.7 |
| baseline | fixed      |       0 |     1 | 15,665,245 |          7.558 |       7.532 |       7.583 |           1.00 |           60.0 |        60.9 |            0.1 |         NA |
| baseline | fixed      |       0 |     2 | 31,330,489 |         15.078 |      15.045 |      15.111 |           1.00 |           61.4 |        61.7 |            0.1 |         NA |
| baseline | fixed      |       0 |     4 | 62,660,978 |         30.553 |      30.539 |      30.555 |           1.00 |           74.6 |        75.0 |            0.1 |         NA |
| branch   | fixed      |       0 |     1 | 15,665,245 |          4.750 |       4.718 |       4.774 |           1.00 |           48.2 |        48.3 |            0.1 |       78.0 |
| branch   | fixed      |       0 |     2 | 31,330,489 |          9.418 |       9.390 |       9.459 |           1.00 |           48.2 |        48.4 |            0.1 |       78.0 |
| branch   | fixed      |       0 |     4 | 62,660,978 |         19.336 |      19.220 |      19.367 |           1.00 |           48.3 |        48.7 |            0.1 |       78.0 |
| baseline | sequence   |       2 |     1 | 15,665,245 |          5.676 |       5.671 |       5.724 |           0.99 |          496.5 |       496.8 |            3.2 |         NA |
| baseline | sequence   |       2 |     2 | 31,330,489 |         11.598 |      11.361 |      12.064 |           1.00 |          566.2 |       571.6 |            3.2 |         NA |
| baseline | sequence   |       2 |     4 | 62,660,978 |         23.608 |      23.395 |      24.101 |           1.00 |          645.4 |       651.5 |            3.3 |         NA |
| branch   | sequence   |       2 |     1 | 15,665,245 |          3.965 |       3.954 |       4.022 |           0.99 |          358.3 |       359.1 |            3.2 |      972.4 |
| branch   | sequence   |       2 |     2 | 31,330,489 |          7.820 |       7.804 |       7.867 |           1.00 |          356.5 |       358.0 |            3.3 |      972.4 |
| branch   | sequence   |       2 |     4 | 62,660,978 |         15.976 |      15.673 |      16.626 |           1.00 |          361.2 |       363.2 |            3.2 |      972.4 |
| baseline | sequence   |       0 |     1 | 15,665,245 |         10.639 |      10.636 |      10.679 |           1.00 |           62.1 |        62.2 |            3.2 |         NA |
| baseline | sequence   |       0 |     2 | 31,330,489 |         21.190 |      21.157 |      21.209 |           1.00 |           63.2 |        65.2 |            3.2 |         NA |
| baseline | sequence   |       0 |     4 | 62,660,978 |         42.718 |      42.685 |      43.069 |           1.00 |           75.5 |        76.3 |            3.3 |         NA |
| branch   | sequence   |       0 |     1 | 15,665,245 |         10.038 |       9.710 |      10.487 |           1.00 |           53.4 |        54.6 |            3.2 |       97.4 |
| branch   | sequence   |       0 |     2 | 31,330,489 |         20.161 |      19.992 |      20.231 |           1.00 |           54.0 |        54.0 |            3.2 |       97.9 |
| branch   | sequence   |       0 |     4 | 62,660,978 |         37.786 |      37.544 |      38.763 |           1.00 |           70.0 |        70.2 |            3.2 |      148.4 |
| baseline | read_group |       2 |     1 | 15,665,245 |          4.692 |       4.687 |       4.725 |           0.99 |          446.0 |       448.4 |            5.6 |         NA |
| baseline | read_group |       2 |     2 | 31,330,489 |          9.427 |       9.294 |       9.680 |           0.99 |          509.5 |       521.3 |            5.5 |         NA |
| baseline | read_group |       2 |     4 | 62,660,978 |         18.726 |      18.687 |      18.908 |           1.00 |          574.2 |       576.0 |            5.6 |         NA |
| branch   | read_group |       2 |     1 | 15,665,245 |          3.498 |       3.486 |       3.549 |           0.99 |          219.1 |       225.6 |            5.6 |    2,291.9 |
| branch   | read_group |       2 |     2 | 31,330,489 |          7.391 |       7.297 |       7.482 |           0.99 |          274.5 |       304.8 |            5.6 |    2,291.9 |
| branch   | read_group |       2 |     4 | 62,660,978 |         15.200 |      15.051 |      15.712 |           1.00 |          332.0 |       338.0 |            5.6 |    2,291.9 |
| baseline | read_group |       0 |     1 | 15,665,245 |          8.375 |       8.347 |       8.384 |           1.00 |           64.4 |        64.8 |            5.5 |         NA |
| baseline | read_group |       0 |     2 | 31,330,489 |         16.606 |      16.556 |      16.654 |           1.00 |           64.1 |        64.2 |            5.5 |         NA |
| baseline | read_group |       0 |     4 | 62,660,978 |         33.897 |      33.876 |      33.927 |           1.00 |           79.3 |        79.6 |            5.6 |         NA |
| branch   | read_group |       0 |     1 | 15,665,245 |          6.421 |       6.414 |       6.436 |           1.00 |           60.8 |        62.5 |            3.1 |      391.3 |
| branch   | read_group |       0 |     2 | 31,330,489 |         12.748 |      12.734 |      12.782 |           1.00 |           62.2 |        62.5 |            5.5 |      391.3 |
| branch   | read_group |       0 |     4 | 62,660,978 |         25.791 |      25.618 |      25.944 |           1.00 |           63.1 |        63.6 |            5.5 |      391.3 |
| baseline | all        |       2 |     1 | 15,665,245 |          9.685 |       9.641 |       9.746 |           1.00 |          514.5 |       519.1 |            6.2 |         NA |
| baseline | all        |       2 |     2 | 31,330,489 |         19.698 |      19.569 |      19.812 |           1.00 |          580.2 |       584.6 |            6.1 |         NA |
| baseline | all        |       2 |     4 | 62,660,978 |         40.120 |      39.971 |      40.220 |           1.00 |          649.6 |       656.0 |            6.1 |         NA |
| branch   | all        |       2 |     1 | 15,665,245 |          9.626 |       9.584 |       9.635 |           1.00 |          511.0 |       512.6 |            6.1 |    3,066.0 |
| branch   | all        |       2 |     2 | 31,330,489 |         19.346 |      19.345 |      19.482 |           1.00 |          581.6 |       581.9 |            6.1 |    3,066.0 |
| branch   | all        |       2 |     4 | 62,660,978 |         39.528 |      39.267 |      39.641 |           1.00 |          660.4 |       668.0 |            6.2 |    3,066.0 |
| baseline | all        |       0 |     1 | 15,665,245 |         15.533 |      15.525 |      15.578 |           1.00 |           66.1 |        66.1 |            6.1 |         NA |
| baseline | all        |       0 |     2 | 31,330,489 |         31.178 |      31.174 |      31.235 |           1.00 |           65.2 |        65.6 |            6.1 |         NA |
| baseline | all        |       0 |     4 | 62,660,978 |         62.462 |      62.403 |      62.699 |           1.00 |           80.0 |        81.9 |            6.2 |         NA |
| branch   | all        |       0 |     1 | 15,665,245 |         15.543 |      15.537 |      15.554 |           1.00 |           63.4 |        63.5 |            6.1 |      410.7 |
| branch   | all        |       0 |     2 | 31,330,489 |         30.955 |      30.875 |      30.978 |           1.00 |           63.8 |        66.4 |            6.2 |      411.2 |
| branch   | all        |       0 |     4 | 62,660,978 |         62.651 |      62.389 |      62.680 |           1.00 |           81.0 |        81.2 |            6.1 |      461.7 |

Before/after: `time_ratio` is the baseline median time over the branch median time, `rss_ratio` the baseline median peak RSS over the branch’s, and `rss_saved_mib` their difference.

| projection | threads | scale | median_seconds_before | median_rss_mib_before | median_seconds_after | median_rss_mib_after | time_ratio | rss_ratio | rss_saved_mib |
|:-----------|--------:|------:|----------------------:|----------------------:|---------------------:|---------------------:|-----------:|----------:|--------------:|
| fixed      |       2 |     1 |                  4.52 |                 415.6 |                 1.86 |                 75.1 |       2.43 |      5.53 |         340.5 |
| fixed      |       2 |     2 |                  8.97 |                 477.7 |                 3.75 |                 77.1 |       2.39 |      6.20 |         400.6 |
| fixed      |       2 |     4 |                 18.01 |                 542.5 |                 7.72 |                 74.3 |       2.33 |      7.30 |         468.2 |
| fixed      |       0 |     1 |                  7.56 |                  60.0 |                 4.75 |                 48.2 |       1.59 |      1.24 |          11.8 |
| fixed      |       0 |     2 |                 15.08 |                  61.4 |                 9.42 |                 48.2 |       1.60 |      1.27 |          13.2 |
| fixed      |       0 |     4 |                 30.55 |                  74.6 |                19.34 |                 48.3 |       1.58 |      1.54 |          26.3 |
| sequence   |       2 |     1 |                  5.68 |                 496.5 |                 3.96 |                358.3 |       1.43 |      1.39 |         138.2 |
| sequence   |       2 |     2 |                 11.60 |                 566.2 |                 7.82 |                356.5 |       1.48 |      1.59 |         209.7 |
| sequence   |       2 |     4 |                 23.61 |                 645.4 |                15.98 |                361.2 |       1.48 |      1.79 |         284.2 |
| sequence   |       0 |     1 |                 10.64 |                  62.1 |                10.04 |                 53.4 |       1.06 |      1.16 |           8.7 |
| sequence   |       0 |     2 |                 21.19 |                  63.2 |                20.16 |                 54.0 |       1.05 |      1.17 |           9.2 |
| sequence   |       0 |     4 |                 42.72 |                  75.5 |                37.79 |                 70.0 |       1.13 |      1.08 |           5.5 |
| read_group |       2 |     1 |                  4.69 |                 446.0 |                 3.50 |                219.1 |       1.34 |      2.04 |         226.9 |
| read_group |       2 |     2 |                  9.43 |                 509.5 |                 7.39 |                274.5 |       1.28 |      1.86 |         235.0 |
| read_group |       2 |     4 |                 18.73 |                 574.2 |                15.20 |                332.0 |       1.23 |      1.73 |         242.2 |
| read_group |       0 |     1 |                  8.38 |                  64.4 |                 6.42 |                 60.8 |       1.30 |      1.06 |           3.6 |
| read_group |       0 |     2 |                 16.61 |                  64.1 |                12.75 |                 62.2 |       1.30 |      1.03 |           1.9 |
| read_group |       0 |     4 |                 33.90 |                  79.3 |                25.79 |                 63.1 |       1.31 |      1.26 |          16.2 |
| all        |       2 |     1 |                  9.68 |                 514.5 |                 9.63 |                511.0 |       1.01 |      1.01 |           3.5 |
| all        |       2 |     2 |                 19.70 |                 580.2 |                19.35 |                581.6 |       1.02 |      1.00 |          -1.4 |
| all        |       2 |     4 |                 40.12 |                 649.6 |                39.53 |                660.4 |       1.01 |      0.98 |         -10.8 |
| all        |       0 |     1 |                 15.53 |                  66.1 |                15.54 |                 63.4 |       1.00 |      1.04 |           2.7 |
| all        |       0 |     2 |                 31.18 |                  65.2 |                30.95 |                 63.8 |       1.01 |      1.02 |           1.4 |
| all        |       0 |     4 |                 62.46 |                  80.0 |                62.65 |                 81.0 |       1.00 |      0.99 |          -1.0 |

Growth of the median time against the alignments read (an exponent of 1 is linear), the size ratio at each doubling, and peak-RSS ratios to 1×. `timed` marks a series whose 1× run takes at least 5 s; only the branch series that are timed have a gate.

| build    | projection | threads | one_x_seconds | ratio_2x | ratio_4x | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | rss_growth_mib | flat  | timed |
|:---------|:-----------|--------:|--------------:|---------:|---------:|-----------------:|-----------------:|-------------:|-------------:|---------------:|:------|:------|
| branch   | fixed      |       2 |         1.865 |        2 |        2 |            1.007 |            1.042 |        1.027 |        0.989 |           -0.8 | TRUE  | FALSE |
| branch   | fixed      |       0 |         4.750 |        2 |        2 |            0.987 |            1.038 |        1.000 |        1.002 |            0.1 | TRUE  | FALSE |
| branch   | sequence   |       2 |         3.965 |        2 |        2 |            0.980 |            1.031 |        0.995 |        1.008 |            2.9 | TRUE  | FALSE |
| branch   | sequence   |       0 |        10.038 |        2 |        2 |            1.006 |            0.906 |        1.011 |        1.311 |           16.6 | FALSE | TRUE  |
| branch   | read_group |       2 |         3.498 |        2 |        2 |            1.080 |            1.040 |        1.253 |        1.515 |          112.9 | FALSE | FALSE |
| branch   | read_group |       0 |         6.421 |        2 |        2 |            0.989 |            1.017 |        1.023 |        1.038 |            2.3 | TRUE  | TRUE  |
| branch   | all        |       2 |         9.626 |        2 |        2 |            1.007 |            1.031 |        1.138 |        1.292 |          149.4 | FALSE | TRUE  |
| branch   | all        |       0 |        15.543 |        2 |        2 |            0.994 |            1.017 |        1.006 |        1.278 |           17.6 | FALSE | TRUE  |
| baseline | fixed      |       2 |         4.524 |        2 |        2 |            0.988 |            1.005 |        1.149 |        1.305 |          126.9 | FALSE | FALSE |
| baseline | fixed      |       0 |         7.558 |        2 |        2 |            0.996 |            1.019 |        1.023 |        1.243 |           14.6 | FALSE | TRUE  |
| baseline | sequence   |       2 |         5.676 |        2 |        2 |            1.031 |            1.025 |        1.140 |        1.300 |          148.9 | FALSE | TRUE  |
| baseline | sequence   |       0 |        10.639 |        2 |        2 |            0.994 |            1.011 |        1.018 |        1.216 |           13.4 | FALSE | TRUE  |
| baseline | read_group |       2 |         4.692 |        2 |        2 |            1.007 |            0.990 |        1.142 |        1.287 |          128.2 | FALSE | FALSE |
| baseline | read_group |       0 |         8.375 |        2 |        2 |            0.987 |            1.029 |        0.995 |        1.231 |           14.9 | FALSE | TRUE  |
| baseline | all        |       2 |         9.685 |        2 |        2 |            1.024 |            1.026 |        1.128 |        1.263 |          135.1 | FALSE | TRUE  |
| baseline | all        |       0 |        15.533 |        2 |        2 |            1.005 |            1.002 |        0.986 |        1.210 |           13.9 | FALSE | TRUE  |

## Findings

Every cell returns the same result in both builds, at both thread settings and in every repetition. No observation spills, and every branch observation is within its ceiling. The render enforces all of this.

**fixed.** Over the three sizes, the baseline takes 2.33 to 2.43 times the branch’s time with two decompression threads and 1.58 to 1.60 times without them. Its peak RSS is 5.53 to 7.30 times the branch’s with two threads (baseline minus branch: +468.2 MiB at 4×) and 1.24 to 1.54 times without them (+26.3 MiB at 4×).

**sequence.** Over the three sizes, the baseline takes 1.43 to 1.48 times the branch’s time with two decompression threads and 1.05 to 1.13 times without them. Its peak RSS is 1.39 to 1.79 times the branch’s with two threads (baseline minus branch: +284.2 MiB at 4×) and 1.08 to 1.17 times without them (+5.5 MiB at 4×).

**read_group.** Over the three sizes, the baseline takes 1.23 to 1.34 times the branch’s time with two decompression threads and 1.30 to 1.31 times without them. Its peak RSS is 1.73 to 2.04 times the branch’s with two threads (baseline minus branch: +242.2 MiB at 4×) and 1.03 to 1.26 times without them (+16.2 MiB at 4×).

**all.** Over the three sizes, the baseline takes 1.01 to 1.02 times the branch’s time with two decompression threads and 1.00 to 1.01 times without them. Its peak RSS is 0.98 to 1.01 times the branch’s with two threads (baseline minus branch: -10.8 MiB at 4×) and 0.99 to 1.04 times without them (-1.0 MiB at 4×).

Flat in memory (at most 10% more peak RSS at 4× than at 1×), branch build: `fixed` with 2 threads (-0.8 MiB from 1× to 4×, a ratio of 0.99); `fixed` with 0 threads (0.1 MiB from 1× to 4×, a ratio of 1.00); `sequence` with 2 threads (2.9 MiB from 1× to 4×, a ratio of 1.01); `read_group` with 0 threads (2.3 MiB from 1× to 4×, a ratio of 1.04).

Not flat, branch build: `sequence` with 0 threads (16.6 MiB from 1× to 4×, a ratio of 1.31); `read_group` with 2 threads (112.9 MiB from 1× to 4×, a ratio of 1.52); `all` with 2 threads (149.4 MiB from 1× to 4×, a ratio of 1.29); `all` with 0 threads (17.6 MiB from 1× to 4×, a ratio of 1.28).

The peak RSS of the branch build above the overhead is between 0.19 (`read_group`, 0 threads, 1×) and 1.04 (`sequence`, 2 threads, 4×) times the declared live state, against the factor of 3 that the ceiling allows. The record-buffer and slice terms use the largest records of the contig, so the model is an upper bound, not a prediction.

The branch series with a timing verdict (`sequence` with 0 threads, `read_group` with 0 threads, `all` with 2 threads, `all` with 0 threads) have time exponents of at most 1.03. `fixed` with 2 threads, `fixed` with 0 threads, `sequence` with 2 threads, `read_group` with 2 threads take less than 5 s at 1× and get no timing verdict.

The `all` projection asks htslib for every field in both builds. At 4× with two threads, its baseline-to-branch time ratio is 1.01 and its peak-RSS ratio 0.98; these show how far the two builds differ when the projection does not change what htslib decodes.
