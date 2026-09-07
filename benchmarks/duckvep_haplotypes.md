Phased replay: native stream and public SQL
================

<!-- duckvep_haplotypes.md is generated from duckvep_haplotypes.Rmd. -->

Source: b6f0016ffd287336e8c1777740f67457016d7cbe. This is the first
recorded baseline for the current **literal phased replay**
implementation, not combined SO/HGVS. All runs use one thread pinned to
CPU 2 on an Intel Core i5-13500, DuckDB 1.5.3, and a source-bound clean
extension build. The native bridge uses the same kernel sources compiled
with `-O3 -DNDEBUG`; receipts retain compiler/version, binary and
fixture hashes, worker jobs, GNU-time logs and complete result
fingerprints.

## Workload and verification

The [registry](../r/duckhtsbench/inst/benchmark_registry.tsv) stages two
committed, checksum-verified fixtures without network access: a 180-base
CDS and four ALT events. Two substitutions share a codon; an
insertion/deletion pair displaces and restores the frame. Four diploid
sample templates repeat across the cohort, producing exactly three
occupied paths per transcript. Strict phase uses PS=10. Transcript
groups share genomic intervals, with non-overlapping groups 1,000 bases
apart; `overlap` is the number of simultaneously active transcripts.

This deliberately favourable prefix-sharing workload is not a real
population, a long-transcript benchmark, a reader benchmark, or evidence
for arbitrary phase and structural configurations. Native input is a
factorized, already-ordered feed; SQL input is its complete
event/candidate/sample relation stored in hash order. Both retain every
source identity and called carrier.

Before timing, native validation checks every transcript, complete
CDS/protein against independently rebuilt/Biostrings-translated
sequences, contributors and sample/lane/PS/ploidy membership. SQL
validation compares the full carrier multiset with decoded input, checks
sequences and edit provenance, and rejects an injected duplicate. Native
validation also rejects a corrupted expected sequence. Every measured
pass must retain all common native/SQL denominators and SQL’s full,
canonicalized nested-row fingerprint; failed runs are not promoted.

| transcripts | samples | overlap | input_records | projected_events | input_calls | output_leaves | output_carriers | prefixes_created | translated_bases |
|------------:|--------:|--------:|--------------:|-----------------:|------------:|--------------:|----------------:|-----------------:|-----------------:|
|        1024 |       4 |       1 |          4096 |             4096 |       16384 |          3072 |            7168 |             6144 |           552960 |
|        1024 |      64 |       1 |          4096 |             4096 |      262144 |          3072 |          114688 |             6144 |           552960 |
|        1024 |      64 |      16 |           256 |             4096 |      262144 |          3072 |          114688 |             6144 |           552960 |
|        1024 |      64 |      64 |            64 |             4096 |      262144 |          3072 |          114688 |             6144 |           552960 |
|        1024 |     256 |       1 |          4096 |             4096 |     1048576 |          3072 |          458752 |             6144 |           552960 |
|       10240 |      64 |      16 |          2560 |            40960 |     2621440 |         30720 |         1146880 |            61440 |          5529600 |

`input_records` counts distinct source ALT events, not candidate/sample
rows. `projected_events` counts event/transcript pairs.
`output_carriers` counts memberships of sample/lane keys in completed
leaves; `output_leaves` counts unique occupied prefixes, not unique
peptide strings. Prefix creation and active pool peaks below are
observed in the native stream, not inferred SQL counters.

## Timing and process memory

| transcripts | samples | overlap | mode   | min_s | median_s | max_s | max_process_rss_mib |
|------------:|--------:|--------:|:-------|------:|---------:|------:|--------------------:|
|        1024 |       4 |       1 | native | 0.003 |    0.003 | 0.003 |              74.477 |
|        1024 |       4 |       1 | sql    | 0.083 |    0.084 | 0.084 |             210.406 |
|        1024 |      64 |       1 | native | 0.013 |    0.013 | 0.013 |              74.477 |
|        1024 |      64 |       1 | sql    | 0.842 |    0.850 | 0.856 |             392.992 |
|        1024 |      64 |      16 | native | 0.017 |    0.017 | 0.017 |              74.480 |
|        1024 |      64 |      16 | sql    | 0.875 |    0.878 | 0.892 |             393.516 |
|        1024 |      64 |      64 | native | 0.018 |    0.018 | 0.018 |              74.480 |
|        1024 |      64 |      64 | sql    | 0.850 |    0.853 | 0.857 |             393.062 |
|        1024 |     256 |       1 | native | 0.046 |    0.047 | 0.048 |              74.480 |
|        1024 |     256 |       1 | sql    | 3.438 |    3.446 | 3.488 |            1030.109 |
|       10240 |      64 |      16 | native | 0.167 |    0.167 | 0.168 |              74.480 |
|       10240 |      64 |      16 | sql    | 9.105 |    9.162 | 9.300 |            2213.879 |

Each recorded pass uses a fresh process and a full warm-up. Native
timing starts after workspace initialization and includes ordered-feed
construction, GT routing, projection, replay/translation, block
construction and a count/length sink. It does **not** include alignment,
DuckDB, or nested output materialization.

SQL timing covers
`CREATE OR REPLACE TABLE ... AS SELECT * FROM duckvep_haplotypes(...)`:
preparation, phase-domain discovery, mandatory sort, workspace
initialization, replay, both sequence alignments and **all** nested
output columns. Input/model staging is outside timing. The warm-up
output table is replaced by the measured table. JSON sizing and
fingerprint aggregation occur afterward, not in the timer; no
file-output throughput is claimed.

GNU `/usr/bin/time -v` measures whole-process maximum RSS, including the
R host, model/input tables, warm-up, measured materialization and
post-run aggregates. The heavier SQL correctness audit runs separately
and is excluded from recorded RSS. RSS is not an isolated operator
allocation measurement. These two execution scopes are intentionally
different: their times are not a speedup comparison.

| transcripts | samples | overlap | output_leaves | cds_bytes | protein_bytes | json_bytes |
|------------:|--------:|--------:|--------------:|----------:|--------------:|-----------:|
|        1024 |       4 |       1 |          3072 |    552960 |        184320 |    4622074 |
|        1024 |      64 |       1 |          3072 |    552960 |        184320 |   11665146 |
|        1024 |      64 |      16 |          3072 |    552960 |        184320 |   11637374 |
|        1024 |      64 |      64 |          3072 |    552960 |        184320 |   11621374 |
|        1024 |     256 |       1 |          3072 |    552960 |        184320 |   34508538 |
|       10240 |      64 |      16 |         30720 |   5529600 |       1843200 |  116648334 |

CDS/protein bytes count the sequence views returned once per occupied
leaf. JSON bytes are the measured UTF-8 size of complete canonicalized
SQL rows, including carriers, contributors, blocks and differences; they
are not bytes written to disk.

## Sparse native state

| transcripts | samples | overlap | peak_transcripts | peak_carriers | peak_prefixes | peak_events | peak_projections | peak_allele_bytes | workspace_bytes | model_bytes |
|------------:|--------:|--------:|-----------------:|--------------:|--------------:|------------:|-----------------:|------------------:|----------------:|------------:|
|        1024 |       4 |       1 |                1 |             7 |             6 |           4 |                4 |                10 |            4484 |       34136 |
|        1024 |      64 |       1 |                1 |           112 |             6 |           4 |                4 |                10 |           14084 |       34136 |
|        1024 |      64 |      16 |               16 |          1792 |            96 |           4 |               64 |                10 |          189104 |       34136 |
|        1024 |      64 |      64 |               64 |          7168 |           384 |           4 |              256 |                10 |          749168 |       34136 |
|        1024 |     256 |       1 |                1 |           448 |             6 |           4 |                4 |                10 |           44804 |       34136 |
|       10240 |      64 |      16 |               16 |          1792 |            96 |           4 |               64 |                10 |          189104 |      338264 |

Native workspace bytes count actual preallocated buffer capacities and
state, separately from immutable model allocations. The fixed fixture is
borrowed; validation additionally owns `transcripts + 2 * samples` bytes
outside the timed kernel workspace. Peak pool occupancy is measured by
the production stream. Increasing total transcripts at fixed cohort size
and overlap does not increase that workspace; immutable model storage
still grows with transcript count.

Both paths use the same active-pool limits. SQL additionally allows
65,536 alignment cells and 64 differing runs per axis, with a 64 MiB
native workspace ceiling and a 4 GB DuckDB memory limit. SQL input,
sorting and output memory belong to DuckDB and may grow with the
complete relation. The process RSS table must not be presented as a
constant-total-memory guarantee.

There is no identical predecessor workload for a regression claim. Three
repeated passes, one machine and this deliberately shared synthetic
cohort do not establish production throughput or a statistical
performance change. Combined SO/HGVS and typed structural composition
will require renewed measurements when implemented.

## Reproduction

Run [the R driver](duckvep_haplotypes.R) with a clean source-bound
extension receipt from the shared evidence builder; it stages inputs
through `duckhtsbench`. All six configurations above use
`--passes 3 --cpu 2`. To prepare the receipt from a clean checkout (no
network access in the build):

Pass `--extension-receipt /tmp/haplotype-extension.tsv` to the R driver.
For a network-free smoke test without a clean-build certificate:

Diagnostic results are retained but are not accepted into the recorded
baseline.
