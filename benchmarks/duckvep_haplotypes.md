Phased replay: native stream and public SQL
================

<!-- duckvep_haplotypes.md is generated from duckvep_haplotypes.Rmd. -->

Current source: 8f9987e3826018cfa73c155eebac7a956b6dc024; same-input
predecessor: d1c591b76f8a9a07036736ac0666a004eb58e0eb. Native
measurements cover literal phased replay. SQL materializes all current
fields, including local coding-block SO; whole-haplotype SO/HGVS is
unfinished. All runs use one thread pinned to CPU 2 on an Intel Core
i5-13500, DuckDB 1.5.3, and a source-bound clean extension build. The
native bridge uses the same kernel sources compiled with `-O3 -DNDEBUG`;
receipts retain compiler/version, binary and fixture hashes, worker
jobs, GNU-time logs and complete result fingerprints.

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
canonicalized nested-row fingerprint; failed runs are not promoted. Full
fingerprints must also agree across revisions sharing an output
contract. An additional projection of every pre-existing field must
agree across the additive coding-block schema change; it does not
replace validation of the new full output.

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
|        1024 |       4 |       1 | native | 0.003 |    0.003 | 0.003 |              74.637 |
|        1024 |       4 |       1 | sql    | 0.051 |    0.052 | 0.053 |             216.789 |
|        1024 |      64 |       1 | native | 0.013 |    0.013 | 0.013 |              74.633 |
|        1024 |      64 |       1 | sql    | 0.343 |    0.345 | 0.345 |             375.930 |
|        1024 |      64 |      16 | native | 0.017 |    0.018 | 0.018 |              74.633 |
|        1024 |      64 |      16 | sql    | 0.352 |    0.355 | 0.356 |             376.910 |
|        1024 |      64 |      64 | native | 0.018 |    0.018 | 0.024 |              74.480 |
|        1024 |      64 |      64 | sql    | 0.347 |    0.348 | 0.351 |             377.922 |
|        1024 |     256 |       1 | native | 0.047 |    0.047 | 0.047 |              74.480 |
|        1024 |     256 |       1 | sql    | 1.319 |    1.319 | 1.329 |             957.906 |
|       10240 |      64 |      16 | native | 0.167 |    0.168 | 0.171 |              74.633 |
|       10240 |      64 |      16 | sql    | 4.146 |    4.194 | 4.340 |            2135.238 |

The current SQL output adds each block’s local SO mask, coding status
and position relative to the first stop. It opens a shared coding
context on completed replay and evaluates physical blocks; it does not
rebuild each leaf again. Consequently this is a same-input cost
comparison for additional work and wider output, not an identical-output
speedup or a claim of no regression. The native count sink still omits
local SO evaluation. Prior sequence, provenance and difference fields
retain their complete projected fingerprints; new full fingerprints and
byte counts are recorded separately.

| transcripts | samples | overlap | median_s_before | median_s_after | median_change_percent | max_process_rss_mib_before | max_process_rss_mib_after |
|------------:|--------:|--------:|----------------:|---------------:|----------------------:|---------------------------:|--------------------------:|
|        1024 |       4 |       1 |           0.050 |          0.052 |                 4.000 |                    213.789 |                   216.789 |
|        1024 |      64 |       1 |           0.335 |          0.345 |                 2.985 |                    375.457 |                   375.930 |
|        1024 |      64 |      16 |           0.353 |          0.355 |                 0.567 |                    376.910 |                   376.910 |
|        1024 |      64 |      64 |           0.343 |          0.348 |                 1.458 |                    377.098 |                   377.922 |
|        1024 |     256 |       1 |           1.316 |          1.319 |                 0.228 |                    957.355 |                   957.906 |
|       10240 |      64 |      16 |           4.184 |          4.194 |                 0.239 |                   2130.625 |                  2135.238 |

Each recorded pass uses a fresh process and a full warm-up. Native
timing starts after workspace initialization and includes ordered-feed
construction, GT routing, projection, replay/translation, block
construction and a count/length sink. It does **not** include alignment,
DuckDB, or nested output materialization.

SQL timing covers
`CREATE OR REPLACE TABLE ... AS SELECT * FROM duckvep_haplotypes(...)`:
preparation, phase-domain discovery, mandatory sort, workspace
initialization, replay, both sequence alignments, local coding-block SO
and **all** nested output columns. Input/model staging is outside
timing. The warm-up output table is replaced by the measured table. JSON
sizing and fingerprint aggregation occur afterward, not in the timer; no
file-output throughput is claimed.

GNU `/usr/bin/time -v` measures whole-process maximum RSS, including the
R host, model/input tables, warm-up, measured materialization and
post-run aggregates. The heavier SQL correctness audit runs separately
and is excluded from recorded RSS. RSS is not an isolated operator
allocation measurement. These two execution scopes are intentionally
different: their times are not a speedup comparison.

| transcripts | samples | overlap | output_leaves | cds_bytes | protein_bytes | json_bytes |
|------------:|--------:|--------:|--------------:|----------:|--------------:|-----------:|
|        1024 |       4 |       1 |          3072 |    552960 |        184320 |    4933370 |
|        1024 |      64 |       1 |          3072 |    552960 |        184320 |   11976442 |
|        1024 |      64 |      16 |          3072 |    552960 |        184320 |   11948670 |
|        1024 |      64 |      64 |          3072 |    552960 |        184320 |   11932670 |
|        1024 |     256 |       1 |          3072 |    552960 |        184320 |   34819834 |
|       10240 |      64 |      16 |         30720 |   5529600 |       1843200 |  119761294 |

CDS/protein bytes count the sequence views returned once per occupied
leaf. JSON bytes are the measured UTF-8 size of complete canonicalized
SQL rows, including carriers, contributors, blocks and differences; they
are not bytes written to disk.

## Sparse native state

| transcripts | samples | overlap | peak_transcripts | peak_carriers | peak_prefixes | peak_events | peak_projections | peak_allele_bytes | workspace_bytes | model_bytes |
|------------:|--------:|--------:|-----------------:|--------------:|--------------:|------------:|-----------------:|------------------:|----------------:|------------:|
|        1024 |       4 |       1 |                1 |             7 |             6 |           4 |                4 |                10 |            4516 |       34136 |
|        1024 |      64 |       1 |                1 |           112 |             6 |           4 |                4 |                10 |           14116 |       34136 |
|        1024 |      64 |      16 |               16 |          1792 |            96 |           4 |               64 |                10 |          189136 |       34136 |
|        1024 |      64 |      64 |               64 |          7168 |           384 |           4 |              256 |                10 |          749200 |       34136 |
|        1024 |     256 |       1 |                1 |           448 |             6 |           4 |                4 |                10 |           44836 |       34136 |
|       10240 |      64 |      16 |               16 |          1792 |            96 |           4 |               64 |                10 |          189136 |      338264 |

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

This matched before/after comparison retains three repeated passes per
workload and revision. One machine and this deliberately shared
synthetic cohort do not establish production throughput or statistical
significance. Whole-haplotype SO/HGVS and typed structural composition
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
