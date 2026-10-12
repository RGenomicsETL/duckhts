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

Source revision: 68dfa0332ad9ea34eea556eddea231cd83d842f6; `src` tree 2de704e472ceeda276caa0eeef7093221378f3c4. Baseline revision: a5ab5147abb7abe398583e26ded5183ae7ecd82f. DuckDB runtime: v1.5.1 (Variegata) 7dbb2e646f (CLI). htslib 1.24. Host load at start: 1.15, 1.28, 0.98. CPU affinity: pid 510410’s current affinity list: 8-15. DuckDB memory limit: 50.0 GiB (the default); temporary directory limit: 0. Every observation is in [`benchmark_read_bam_cram_projection_runs.csv`](benchmark_read_bam_cram_projection_runs.csv).

| build    | sha256                                                           | empty_run_mib |
|:---------|:-----------------------------------------------------------------|--------------:|
| branch   | 66598566776aa44faffb68a8e0ba62961346c11db160962d115c30291d33586e |          38.7 |
| baseline | 274dee1c390b25474a251d89508badacc15414a9d4d9695700f47dbd96a2a279 |          38.2 |

Commits between the baseline and the branch revision that change `src/bam_reader.c`:

    68dfa033 Decode CRAM tags and sequence when read_bam projects CIGAR
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
| fixed      |       2 |     1 |                5 |       18.5 |               13.7 |           0.0 |           41.0 |      161.6 |
| fixed      |       2 |     2 |                5 |       18.5 |               13.7 |           0.0 |           41.0 |      161.6 |
| fixed      |       2 |     4 |                5 |       18.5 |               13.7 |           0.0 |           41.0 |      161.6 |
| fixed      |       0 |     1 |                1 |        4.3 |                0.0 |           0.0 |           13.1 |       77.9 |
| fixed      |       0 |     2 |                1 |        4.3 |                0.0 |           0.0 |           13.1 |       77.9 |
| fixed      |       0 |     4 |                1 |        4.3 |                0.0 |           0.0 |           13.1 |       77.9 |
| sequence   |       2 |     1 |                5 |       38.5 |               25.9 |         237.4 |          311.2 |      972.4 |
| sequence   |       2 |     2 |                5 |       38.5 |               25.9 |         237.4 |          311.2 |      972.4 |
| sequence   |       2 |     4 |                5 |       38.5 |               25.9 |         237.4 |          311.2 |      972.4 |
| sequence   |       0 |     1 |                1 |       10.0 |                0.0 |           0.2 |           19.5 |       97.4 |
| sequence   |       0 |     2 |                1 |       10.0 |                0.0 |           0.3 |           19.7 |       97.9 |
| sequence   |       0 |     4 |                1 |       10.0 |                0.0 |          17.2 |           36.5 |      148.4 |
| read_group |       2 |     1 |                5 |      336.6 |              392.2 |           0.0 |          751.0 |     2291.9 |
| read_group |       2 |     2 |                5 |      336.6 |              392.2 |           0.0 |          751.0 |     2291.9 |
| read_group |       2 |     4 |                5 |      336.6 |              392.2 |           0.0 |          751.0 |     2291.9 |
| read_group |       0 |     1 |                1 |       95.2 |                0.0 |           0.0 |          117.5 |      391.2 |
| read_group |       0 |     2 |                1 |       95.2 |                0.0 |           0.0 |          117.5 |      391.2 |
| read_group |       0 |     4 |                1 |       95.2 |                0.0 |           0.0 |          117.5 |      391.2 |
| all        |       2 |     1 |                5 |      356.7 |              392.2 |         237.4 |         1009.1 |     3066.0 |
| all        |       2 |     2 |                5 |      356.7 |              392.2 |         237.4 |         1009.1 |     3066.0 |
| all        |       2 |     4 |                5 |      356.7 |              392.2 |         237.4 |         1009.1 |     3066.0 |
| all        |       0 |     1 |                1 |      100.9 |                0.0 |           0.2 |          124.0 |      410.6 |
| all        |       0 |     2 |                1 |      100.9 |                0.0 |           0.3 |          124.1 |      411.1 |
| all        |       0 |     4 |                1 |      100.9 |                0.0 |          17.2 |          141.0 |      461.6 |

## Measurements

`seconds` are the latency of the query over the repetitions. `cpu_per_second` is DuckDB’s CPU time per second of latency; htslib’s decompression threads are not in it. `budget_mib` is the ceiling of a branch cell.

| build    | projection | threads | scale | alignments | median_seconds | min_seconds | max_seconds | cpu_per_second | median_rss_mib | max_rss_mib | max_buffer_mib | budget_mib |
|:---------|:-----------|--------:|------:|-----------:|---------------:|------------:|------------:|---------------:|---------------:|------------:|---------------:|-----------:|
| baseline | fixed      |       2 |     1 | 15,665,245 |          4.602 |       4.571 |       4.639 |           0.99 |          420.6 |       439.0 |            0.1 |         NA |
| baseline | fixed      |       2 |     2 | 31,330,489 |          8.946 |       8.899 |       8.982 |           0.99 |          476.5 |       483.9 |            0.1 |         NA |
| baseline | fixed      |       2 |     4 | 62,660,978 |         18.171 |      18.030 |      18.205 |           1.00 |          531.6 |       543.5 |            0.1 |         NA |
| branch   | fixed      |       2 |     1 | 15,665,245 |          1.909 |       1.889 |       1.927 |           0.99 |           76.2 |        78.7 |            0.1 |      161.6 |
| branch   | fixed      |       2 |     2 | 31,330,489 |          3.750 |       3.694 |       4.064 |           0.99 |           78.0 |        80.9 |            0.2 |      161.6 |
| branch   | fixed      |       2 |     4 | 62,660,978 |          7.574 |       7.559 |       7.998 |           1.00 |           77.2 |        77.5 |            0.1 |      161.6 |
| baseline | fixed      |       0 |     1 | 15,665,245 |          7.575 |       7.491 |       8.682 |           1.00 |           60.8 |        61.3 |            0.1 |         NA |
| baseline | fixed      |       0 |     2 | 31,330,489 |         15.087 |      15.003 |      15.097 |           1.00 |           60.8 |        61.6 |            0.1 |         NA |
| baseline | fixed      |       0 |     4 | 62,660,978 |         30.813 |      30.580 |      31.504 |           1.00 |           74.6 |        74.9 |            0.2 |         NA |
| branch   | fixed      |       0 |     1 | 15,665,245 |          4.748 |       4.745 |       4.758 |           1.00 |           48.3 |        48.4 |            0.1 |       77.9 |
| branch   | fixed      |       0 |     2 | 31,330,489 |          9.541 |       9.415 |       9.584 |           1.00 |           48.4 |        48.6 |            0.1 |       77.9 |
| branch   | fixed      |       0 |     4 | 62,660,978 |         19.371 |      19.300 |      19.439 |           1.00 |           48.5 |        48.6 |            0.1 |       77.9 |
| baseline | sequence   |       2 |     1 | 15,665,245 |          5.715 |       5.677 |       5.983 |           0.99 |          498.6 |       499.4 |            3.3 |         NA |
| baseline | sequence   |       2 |     2 | 31,330,489 |         11.608 |      11.358 |      11.629 |           1.00 |          563.8 |       564.7 |            3.3 |         NA |
| baseline | sequence   |       2 |     4 | 62,660,978 |         23.482 |      23.366 |      24.670 |           1.00 |          644.4 |       647.7 |            3.3 |         NA |
| branch   | sequence   |       2 |     1 | 15,665,245 |          3.997 |       3.997 |       4.078 |           0.99 |          359.3 |       359.6 |            3.2 |      972.4 |
| branch   | sequence   |       2 |     2 | 31,330,489 |          7.833 |       7.746 |       7.892 |           1.00 |          359.2 |       360.8 |            3.2 |      972.4 |
| branch   | sequence   |       2 |     4 | 62,660,978 |         15.679 |      15.501 |      15.857 |           1.00 |          357.4 |       357.5 |            3.2 |      972.4 |
| baseline | sequence   |       0 |     1 | 15,665,245 |         10.659 |      10.643 |      10.696 |           1.00 |           60.3 |        62.1 |            3.2 |         NA |
| baseline | sequence   |       0 |     2 | 31,330,489 |         21.195 |      21.157 |      21.248 |           1.00 |           63.1 |        63.4 |            3.2 |         NA |
| baseline | sequence   |       0 |     4 | 62,660,978 |         42.768 |      42.725 |      43.019 |           1.00 |           77.2 |        77.8 |            3.2 |         NA |
| branch   | sequence   |       0 |     1 | 15,665,245 |          9.923 |       9.880 |       9.946 |           1.00 |           54.0 |        54.1 |            3.2 |       97.4 |
| branch   | sequence   |       0 |     2 | 31,330,489 |         19.962 |      19.821 |      20.061 |           1.00 |           52.9 |        53.1 |            3.2 |       97.9 |
| branch   | sequence   |       0 |     4 | 62,660,978 |         37.858 |      37.812 |      38.018 |           1.00 |           70.1 |        70.3 |            3.2 |      148.4 |
| baseline | read_group |       2 |     1 | 15,665,245 |          4.685 |       4.675 |       4.724 |           0.99 |          450.4 |       453.2 |            5.6 |         NA |
| baseline | read_group |       2 |     2 | 31,330,489 |          9.279 |       9.267 |       9.776 |           0.99 |          508.9 |       516.9 |            5.6 |         NA |
| baseline | read_group |       2 |     4 | 62,660,978 |         18.985 |      18.809 |      19.210 |           1.00 |          573.7 |       581.8 |            5.6 |         NA |
| branch   | read_group |       2 |     1 | 15,665,245 |          3.559 |       3.486 |       3.639 |           0.99 |          209.9 |       210.5 |            5.6 |    2,291.9 |
| branch   | read_group |       2 |     2 | 31,330,489 |          7.444 |       7.336 |       7.459 |           0.99 |          273.0 |       297.8 |            5.6 |    2,291.9 |
| branch   | read_group |       2 |     4 | 62,660,978 |         15.411 |      15.394 |      15.549 |           1.00 |          322.6 |       332.9 |            5.6 |    2,291.9 |
| baseline | read_group |       0 |     1 | 15,665,245 |          8.391 |       8.352 |       8.426 |           1.00 |           64.7 |        66.0 |            5.5 |         NA |
| baseline | read_group |       0 |     2 | 31,330,489 |         16.669 |      16.653 |      16.792 |           1.00 |           64.2 |        64.2 |            5.5 |         NA |
| baseline | read_group |       0 |     4 | 62,660,978 |         33.909 |      33.893 |      34.162 |           1.00 |           79.4 |        79.6 |            5.5 |         NA |
| branch   | read_group |       0 |     1 | 15,665,245 |          6.448 |       6.389 |       6.487 |           1.00 |           62.4 |        62.6 |            5.5 |      391.2 |
| branch   | read_group |       0 |     2 | 31,330,489 |         12.748 |      12.738 |      12.937 |           1.00 |           62.3 |        62.7 |            5.5 |      391.2 |
| branch   | read_group |       0 |     4 | 62,660,978 |         25.995 |      25.923 |      26.182 |           1.00 |           62.5 |        65.1 |            5.6 |      391.2 |
| baseline | all        |       2 |     1 | 15,665,245 |          9.728 |       9.692 |      10.282 |           1.00 |          504.5 |       509.3 |            6.2 |         NA |
| baseline | all        |       2 |     2 | 31,330,489 |         19.643 |      19.635 |      19.699 |           1.00 |          581.4 |       587.9 |            6.1 |         NA |
| baseline | all        |       2 |     4 | 62,660,978 |         39.977 |      39.717 |      40.091 |           1.00 |          638.0 |       652.7 |            6.1 |         NA |
| branch   | all        |       2 |     1 | 15,665,245 |          9.591 |       9.561 |       9.596 |           1.00 |          515.6 |       523.0 |            6.1 |    3,066.0 |
| branch   | all        |       2 |     2 | 31,330,489 |         19.429 |      19.413 |      19.656 |           1.00 |          587.1 |       598.3 |            6.2 |    3,066.0 |
| branch   | all        |       2 |     4 | 62,660,978 |         39.488 |      39.231 |      39.608 |           1.00 |          653.3 |       654.8 |            6.1 |    3,066.0 |
| baseline | all        |       0 |     1 | 15,665,245 |         15.610 |      15.581 |      15.757 |           1.00 |           65.8 |        66.0 |            6.1 |         NA |
| baseline | all        |       0 |     2 | 31,330,489 |         31.223 |      31.207 |      31.259 |           1.00 |           65.6 |        65.9 |            6.2 |         NA |
| baseline | all        |       0 |     4 | 62,660,978 |         62.665 |      62.386 |      63.279 |           1.00 |           80.0 |        81.8 |            6.1 |         NA |
| branch   | all        |       0 |     1 | 15,665,245 |         15.767 |      15.748 |      15.912 |           1.00 |           66.0 |        66.7 |            6.1 |      410.6 |
| branch   | all        |       0 |     2 | 31,330,489 |         31.188 |      31.074 |      31.268 |           1.00 |           65.6 |        66.6 |            6.2 |      411.1 |
| branch   | all        |       0 |     4 | 62,660,978 |         62.595 |      62.588 |      62.882 |           1.00 |           80.6 |        80.7 |            6.1 |      461.6 |

Before/after: `time_ratio` is the baseline median time over the branch median time, `rss_ratio` the baseline median peak RSS over the branch’s, and `rss_saved_mib` their difference.

| projection | threads | scale | median_seconds_before | median_rss_mib_before | median_seconds_after | median_rss_mib_after | time_ratio | rss_ratio | rss_saved_mib |
|:-----------|--------:|------:|----------------------:|----------------------:|---------------------:|---------------------:|-----------:|----------:|--------------:|
| fixed      |       2 |     1 |                  4.60 |                 420.6 |                 1.91 |                 76.2 |       2.41 |      5.52 |         344.4 |
| fixed      |       2 |     2 |                  8.95 |                 476.5 |                 3.75 |                 78.0 |       2.39 |      6.11 |         398.5 |
| fixed      |       2 |     4 |                 18.17 |                 531.6 |                 7.57 |                 77.2 |       2.40 |      6.89 |         454.4 |
| fixed      |       0 |     1 |                  7.58 |                  60.8 |                 4.75 |                 48.3 |       1.60 |      1.26 |          12.5 |
| fixed      |       0 |     2 |                 15.09 |                  60.8 |                 9.54 |                 48.4 |       1.58 |      1.26 |          12.4 |
| fixed      |       0 |     4 |                 30.81 |                  74.6 |                19.37 |                 48.5 |       1.59 |      1.54 |          26.1 |
| sequence   |       2 |     1 |                  5.72 |                 498.6 |                 4.00 |                359.3 |       1.43 |      1.39 |         139.3 |
| sequence   |       2 |     2 |                 11.61 |                 563.8 |                 7.83 |                359.2 |       1.48 |      1.57 |         204.6 |
| sequence   |       2 |     4 |                 23.48 |                 644.4 |                15.68 |                357.4 |       1.50 |      1.80 |         287.0 |
| sequence   |       0 |     1 |                 10.66 |                  60.3 |                 9.92 |                 54.0 |       1.07 |      1.12 |           6.3 |
| sequence   |       0 |     2 |                 21.19 |                  63.1 |                19.96 |                 52.9 |       1.06 |      1.19 |          10.2 |
| sequence   |       0 |     4 |                 42.77 |                  77.2 |                37.86 |                 70.1 |       1.13 |      1.10 |           7.1 |
| read_group |       2 |     1 |                  4.68 |                 450.4 |                 3.56 |                209.9 |       1.32 |      2.15 |         240.5 |
| read_group |       2 |     2 |                  9.28 |                 508.9 |                 7.44 |                273.0 |       1.25 |      1.86 |         235.9 |
| read_group |       2 |     4 |                 18.99 |                 573.7 |                15.41 |                322.6 |       1.23 |      1.78 |         251.1 |
| read_group |       0 |     1 |                  8.39 |                  64.7 |                 6.45 |                 62.4 |       1.30 |      1.04 |           2.3 |
| read_group |       0 |     2 |                 16.67 |                  64.2 |                12.75 |                 62.3 |       1.31 |      1.03 |           1.9 |
| read_group |       0 |     4 |                 33.91 |                  79.4 |                25.99 |                 62.5 |       1.30 |      1.27 |          16.9 |
| all        |       2 |     1 |                  9.73 |                 504.5 |                 9.59 |                515.6 |       1.01 |      0.98 |         -11.1 |
| all        |       2 |     2 |                 19.64 |                 581.4 |                19.43 |                587.1 |       1.01 |      0.99 |          -5.7 |
| all        |       2 |     4 |                 39.98 |                 638.0 |                39.49 |                653.3 |       1.01 |      0.98 |         -15.3 |
| all        |       0 |     1 |                 15.61 |                  65.8 |                15.77 |                 66.0 |       0.99 |      1.00 |          -0.2 |
| all        |       0 |     2 |                 31.22 |                  65.6 |                31.19 |                 65.6 |       1.00 |      1.00 |           0.0 |
| all        |       0 |     4 |                 62.66 |                  80.0 |                62.59 |                 80.6 |       1.00 |      0.99 |          -0.6 |

Growth of the median time against the alignments read (an exponent of 1 is linear), the size ratio at each doubling, and peak-RSS ratios to 1×. `timed` marks a series whose 1× run takes at least 5 s; only the branch series that are timed have a gate.

| build    | projection | threads | one_x_seconds | ratio_2x | ratio_4x | time_exponent_2x | time_exponent_4x | rss_ratio_2x | rss_ratio_4x | rss_growth_mib | flat  | timed |
|:---------|:-----------|--------:|--------------:|---------:|---------:|-----------------:|-----------------:|-------------:|-------------:|---------------:|:------|:------|
| branch   | fixed      |       2 |         1.909 |        2 |        2 |            0.974 |            1.014 |        1.024 |        1.013 |            1.0 | TRUE  | FALSE |
| branch   | fixed      |       0 |         4.748 |        2 |        2 |            1.007 |            1.022 |        1.002 |        1.004 |            0.2 | TRUE  | FALSE |
| branch   | sequence   |       2 |         3.997 |        2 |        2 |            0.971 |            1.001 |        1.000 |        0.995 |           -1.9 | TRUE  | FALSE |
| branch   | sequence   |       0 |         9.923 |        2 |        2 |            1.008 |            0.923 |        0.980 |        1.298 |           16.1 | FALSE | TRUE  |
| branch   | read_group |       2 |         3.559 |        2 |        2 |            1.065 |            1.050 |        1.301 |        1.537 |          112.7 | FALSE | FALSE |
| branch   | read_group |       0 |         6.448 |        2 |        2 |            0.983 |            1.028 |        0.998 |        1.002 |            0.1 | TRUE  | TRUE  |
| branch   | all        |       2 |         9.591 |        2 |        2 |            1.018 |            1.023 |        1.139 |        1.267 |          137.7 | FALSE | TRUE  |
| branch   | all        |       0 |        15.767 |        2 |        2 |            0.984 |            1.005 |        0.994 |        1.221 |           14.6 | FALSE | TRUE  |
| baseline | fixed      |       2 |         4.602 |        2 |        2 |            0.959 |            1.022 |        1.133 |        1.264 |          111.0 | FALSE | FALSE |
| baseline | fixed      |       0 |         7.575 |        2 |        2 |            0.994 |            1.030 |        1.000 |        1.227 |           13.8 | FALSE | TRUE  |
| baseline | sequence   |       2 |         5.715 |        2 |        2 |            1.022 |            1.016 |        1.131 |        1.292 |          145.8 | FALSE | TRUE  |
| baseline | sequence   |       0 |        10.659 |        2 |        2 |            0.992 |            1.013 |        1.046 |        1.280 |           16.9 | FALSE | TRUE  |
| baseline | read_group |       2 |         4.685 |        2 |        2 |            0.986 |            1.033 |        1.130 |        1.274 |          123.3 | FALSE | FALSE |
| baseline | read_group |       0 |         8.391 |        2 |        2 |            0.990 |            1.024 |        0.992 |        1.227 |           14.7 | FALSE | TRUE  |
| baseline | all        |       2 |         9.728 |        2 |        2 |            1.014 |            1.025 |        1.152 |        1.265 |          133.5 | FALSE | TRUE  |
| baseline | all        |       0 |        15.610 |        2 |        2 |            1.000 |            1.005 |        0.997 |        1.216 |           14.2 | FALSE | TRUE  |

## Findings

Every cell returns the same result in both builds, at both thread settings and in every repetition. No observation spills, and every branch observation is within its ceiling. The render enforces all of this.

**fixed.** Over the three sizes, the baseline takes 2.39 to 2.41 times the branch’s time with two decompression threads and 1.58 to 1.60 times without them. Its peak RSS is 5.52 to 6.89 times the branch’s with two threads (baseline minus branch: +454.4 MiB at 4×) and 1.26 to 1.54 times without them (+26.1 MiB at 4×).

**sequence.** Over the three sizes, the baseline takes 1.43 to 1.50 times the branch’s time with two decompression threads and 1.06 to 1.13 times without them. Its peak RSS is 1.39 to 1.80 times the branch’s with two threads (baseline minus branch: +287.0 MiB at 4×) and 1.10 to 1.19 times without them (+7.1 MiB at 4×).

**read_group.** Over the three sizes, the baseline takes 1.23 to 1.32 times the branch’s time with two decompression threads and 1.30 to 1.31 times without them. Its peak RSS is 1.78 to 2.15 times the branch’s with two threads (baseline minus branch: +251.1 MiB at 4×) and 1.03 to 1.27 times without them (+16.9 MiB at 4×).

**all.** Over the three sizes, the baseline takes 1.01 to 1.01 times the branch’s time with two decompression threads and 0.99 to 1.00 times without them. Its peak RSS is 0.98 to 0.99 times the branch’s with two threads (baseline minus branch: -15.3 MiB at 4×) and 0.99 to 1.00 times without them (-0.6 MiB at 4×).

Flat in memory (at most 10% more peak RSS at 4× than at 1×), branch build: `fixed` with 2 threads (1.0 MiB from 1× to 4×, a ratio of 1.01); `fixed` with 0 threads (0.2 MiB from 1× to 4×, a ratio of 1.00); `sequence` with 2 threads (-1.9 MiB from 1× to 4×, a ratio of 0.99); `read_group` with 0 threads (0.1 MiB from 1× to 4×, a ratio of 1.00).

Not flat, branch build: `sequence` with 0 threads (16.1 MiB from 1× to 4×, a ratio of 1.30); `read_group` with 2 threads (112.7 MiB from 1× to 4×, a ratio of 1.54); `all` with 2 threads (137.7 MiB from 1× to 4×, a ratio of 1.27); `all` with 0 threads (14.6 MiB from 1× to 4×, a ratio of 1.22).

The peak RSS of the branch build above the overhead is between 0.20 (`read_group`, 0 threads, 2×) and 1.03 (`sequence`, 2 threads, 1×) times the declared live state, against the factor of 3 that the ceiling allows. The record-buffer and slice terms use the largest records of the contig, so the model is an upper bound, not a prediction.

The branch series with a timing verdict (`sequence` with 0 threads, `read_group` with 0 threads, `all` with 2 threads, `all` with 0 threads) have time exponents of at most 1.03. `fixed` with 2 threads, `fixed` with 0 threads, `sequence` with 2 threads, `read_group` with 2 threads take less than 5 s at 1× and get no timing verdict.

The `all` projection asks htslib for every field in both builds. At 4× with two threads, its baseline-to-branch time ratio is 1.01 and its peak-RSS ratio 0.98; these show how far the two builds differ when the projection does not change what htslib decodes.
