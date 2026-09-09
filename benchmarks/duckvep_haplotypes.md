Phased replay: native stream and public SQL
================

<!-- duckvep_haplotypes.md is generated from duckvep_haplotypes.Rmd. -->

Current source: 16afccd0fd795df318bbd8b3e66f89b8321da022; same-input
baseline: 4119d55c43fe0649ffe8325135bbdf26c1f37e94. Native measurements
cover literal phased replay. SQL materializes all current fields,
including local coding-block SO. HGVS generation is disabled;
whole-haplotype SO/HGVS is unfinished and its computation is not timed
here. Both paths consume standalone ALT events and decoded calls; native
raw-record parsing and conditional replay are not timed by this
workload. All runs use one thread pinned to CPU 2 on an Intel Core
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
contract. Local coding-block and literal-replay projections provide
separate comparisons with their recorded output contracts; full current
output is always checked. The HGVS-status schema adds `hgvsp` (NULL) and
`hgvsp_status` (`not_requested`) in this workload. Only the separate
prior-schema projections exclude these columns.

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
|        1024 |       4 |       1 | native | 0.003 |    0.003 | 0.004 |              74.477 |
|        1024 |       4 |       1 | sql    | 0.053 |    0.053 | 0.054 |             215.789 |
|        1024 |      64 |       1 | native | 0.014 |    0.015 | 0.024 |              74.633 |
|        1024 |      64 |       1 | sql    | 0.361 |    0.365 | 0.403 |             376.344 |
|        1024 |      64 |      16 | native | 0.018 |    0.018 | 0.018 |              74.480 |
|        1024 |      64 |      16 | sql    | 0.361 |    0.366 | 0.368 |             377.344 |
|        1024 |      64 |      64 | native | 0.020 |    0.021 | 0.021 |              74.480 |
|        1024 |      64 |      64 | sql    | 0.367 |    0.368 | 0.368 |             375.516 |
|        1024 |     256 |       1 | native | 0.051 |    0.051 | 0.060 |              74.480 |
|        1024 |     256 |       1 | sql    | 1.396 |    1.402 | 1.426 |             957.293 |
|       10240 |      64 |      16 | native | 0.184 |    0.185 | 0.188 |              74.477 |
|       10240 |      64 |      16 | sql    | 4.389 |    4.445 | 4.539 |            2134.574 |

Both compared revisions return each block’s local SO mask, coding status
and position relative to the first stop. A shared coding context
evaluates physical blocks on completed replay. Full-output fingerprints
and input/output denominators match across the compared configurations.
The native count sink omits local SO evaluation. HGVS-enabled
performance requires a separate workload.

| transcripts | samples | overlap | median_s_before | median_s_after | median_change_percent | max_process_rss_mib_before | max_process_rss_mib_after |
|------------:|--------:|--------:|----------------:|---------------:|----------------------:|---------------------------:|--------------------------:|
|        1024 |       4 |       1 |           0.051 |          0.053 |                 3.922 |                    216.574 |                   215.789 |
|        1024 |      64 |       1 |           0.330 |          0.365 |                10.606 |                    376.836 |                   376.344 |
|        1024 |      64 |      16 |           0.340 |          0.366 |                 7.647 |                    377.082 |                   377.344 |
|        1024 |      64 |      64 |           0.343 |          0.368 |                 7.289 |                    376.406 |                   375.516 |
|        1024 |     256 |       1 |           1.296 |          1.402 |                 8.179 |                    957.453 |                   957.293 |
|       10240 |      64 |      16 |           4.069 |          4.445 |                 9.241 |                   2133.770 |                  2134.574 |

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
|        1024 |       4 |       1 |                1 |             7 |             6 |           4 |                4 |                10 |            4810 |       34136 |
|        1024 |      64 |       1 |                1 |           112 |             6 |           4 |                4 |                10 |           14410 |       34136 |
|        1024 |      64 |      16 |               16 |          1792 |            96 |           4 |               64 |                10 |          189430 |       34136 |
|        1024 |      64 |      64 |               64 |          7168 |           384 |           4 |              256 |                10 |          749494 |       34136 |
|        1024 |     256 |       1 |                1 |           448 |             6 |           4 |                4 |                10 |           45130 |       34136 |
|       10240 |      64 |      16 |               16 |          1792 |            96 |           4 |               64 |                10 |          189430 |      338264 |

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

This matched comparison retains three repeated passes per workload and
revision. Sanitizer and conformance jobs did not overlap either
benchmark campaign. Both ran at different times on a shared machine; CPU
pinning does not isolate the cause of timing differences between
campaigns. One machine and this deliberately shared synthetic cohort do
not establish production throughput, statistical significance or a
general absence of regression. Whole-haplotype SO/HGVS and typed
structural composition require their own measurements.

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

### Raw source-record SQL lane

`--modes sql_records` measures the public
`input_mode := 'source_records'` path under `vep116_compat`; the default
`--modes native,sql` selects the decoded lanes reported above. All modes
use the same registered reference and four biallelic events. Source rows
carry the ordered ALT list and literal diploid GT text. The timed query
includes source validation, replay-order planning, raw-GT parsing,
allele expansion, sorting, replay and full nested output. File decoding
and model/input staging are outside timing. This is not a standalone raw
native kernel benchmark or a rare-configuration conformance campaign.

The source lane returns four paths per transcript and two carriers per
sample/transcript, including the occupied reference path. Decoded lanes
return three paths per transcript and 7/4 carriers per
sample/transcript. Complete CDS/protein strings, carrier multisets,
contributor IDs, ALT ordinals and block/edit counts are checked against
the fixture. Eleven corruption controls cover changed CDS, protein,
missing CDS, carrier identity, missing ploidy, contributor identity,
missing ALT ordinal, invented HGVS, changed HGVS status, invalid source
GT and duplicate output. They run in the separate correctness worker.

`input_physical_records`, `input_record_sample_calls` and
`input_candidate_sample_rows` are separately counted from the staged
source relation. REF, ALT and undefined-slot descriptors are execution
work, not extra physical records. For this biallelic workload the raw
lane reserves 13 active event slots, `12 * overlap + 1` projection slots
and 192 allele bytes; other limits match the decoded SQL lane.
Whole-process RSS includes both the decoded fixture table and the
derived source table. Output fingerprints include the raw lane’s
complete schema and reference paths; they are not substituted for
decoded-output fingerprints.

No source-bound raw-lane measurement is recorded in the tables above. A
clean-build run with three passes is required before a timing
comparison; different output denominators preclude treating raw versus
decoded time as an implementation speedup.
