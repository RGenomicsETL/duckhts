Phased replay: native stream and public SQL
================

<!-- duckvep_haplotypes.md is generated from duckvep_haplotypes.Rmd. -->

Compound-replay source: ac423de1e6ec0f449b2f5ff2367763e009115b0a;
same-input baseline: 16afccd0fd795df318bbd8b3e66f89b8321da022. Native
measurements cover literal phased replay. SQL materializes all current
fields, including local coding-block SO. HGVS generation is disabled;
whole-haplotype SO/HGVS is unfinished and its computation is not timed
here. Native and decoded SQL consume standalone ALT events and decoded
calls. The additional source-record SQL lane measures literal GT parsing
and replay; standalone native raw-record execution and conditional
source omissions are not timed. All runs use one thread pinned to CPU 2
on an Intel Core i5-13500, DuckDB 1.5.3, and a source-bound clean
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
decoded pass must retain all common native/SQL denominators and SQL’s
full, canonicalized nested-row fingerprint; failed runs are not
promoted. Full fingerprints must also agree across revisions sharing an
output contract. Local coding-block and literal-replay projections
provide separate comparisons with their recorded output contracts; full
current output is always checked. The HGVS-status schema adds `hgvsp`
(NULL) and `hgvsp_status` (`not_requested`) in this workload. Only the
separate prior-schema projections exclude these columns.

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

| transcripts | samples | overlap | mode        | min_s | median_s | max_s | max_process_rss_mib |
|------------:|--------:|--------:|:------------|------:|---------:|------:|--------------------:|
|        1024 |       4 |       1 | native      | 0.003 |    0.003 | 0.003 |              74.324 |
|        1024 |       4 |       1 | sql         | 0.049 |    0.050 | 0.050 |             219.137 |
|        1024 |      64 |       1 | native      | 0.013 |    0.013 | 0.014 |              74.320 |
|        1024 |      64 |       1 | sql         | 0.329 |    0.334 | 0.334 |             378.934 |
|        1024 |      64 |      16 | native      | 0.016 |    0.017 | 0.017 |              74.473 |
|        1024 |      64 |      16 | sql         | 0.336 |    0.337 | 0.337 |             379.340 |
|        1024 |      64 |      16 | sql_records | 0.711 |    0.711 | 0.723 |             610.504 |
|        1024 |      64 |      64 | native      | 0.018 |    0.019 | 0.019 |              74.320 |
|        1024 |      64 |      64 | sql         | 0.338 |    0.341 | 0.342 |             379.176 |
|        1024 |     256 |       1 | native      | 0.045 |    0.047 | 0.047 |              74.465 |
|        1024 |     256 |       1 | sql         | 1.321 |    1.328 | 1.331 |             960.242 |
|       10240 |      64 |      16 | native      | 0.164 |    0.166 | 0.168 |              74.473 |
|       10240 |      64 |      16 | sql         | 4.184 |    4.188 | 4.221 |            2133.344 |

Both compared revisions return each block’s local SO mask, coding status
and position relative to the first stop. A shared coding context
evaluates physical blocks on completed replay. The current output adds
HGVS columns. All prior local-coding-block fields and input/output
denominators match, but complete rows have different schemas and byte
counts. This is a same-input comparison, not identical-output work. The
native count sink omits local SO evaluation. HGVS-enabled performance
requires a separate workload.

| transcripts | samples | overlap | median_s_before | median_s_after | median_change_percent | max_process_rss_mib_before | max_process_rss_mib_after |
|------------:|--------:|--------:|----------------:|---------------:|----------------------:|---------------------------:|--------------------------:|
|        1024 |       4 |       1 |           0.053 |          0.050 |                -5.660 |                    215.789 |                   219.137 |
|        1024 |      64 |       1 |           0.365 |          0.334 |                -8.493 |                    376.344 |                   378.934 |
|        1024 |      64 |      16 |           0.366 |          0.337 |                -7.923 |                    377.344 |                   379.340 |
|        1024 |      64 |      64 |           0.368 |          0.341 |                -7.337 |                    375.516 |                   379.176 |
|        1024 |     256 |       1 |           1.402 |          1.328 |                -5.278 |                    957.293 |                   960.242 |
|       10240 |      64 |      16 |           4.445 |          4.188 |                -5.782 |                   2134.574 |                  2133.344 |

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
different: their times are not a speedup comparison. Post-run
fingerprint aggregates depend on the output schema, so a schema change
also prevents interpreting process RSS as an identical-work memory
comparison.

| transcripts | samples | overlap | output_leaves | cds_bytes | protein_bytes | json_bytes |
|------------:|--------:|--------:|--------------:|----------:|--------------:|-----------:|
|        1024 |       4 |       1 |          3072 |    552960 |        184320 |    5068538 |
|        1024 |      64 |       1 |          3072 |    552960 |        184320 |   12111610 |
|        1024 |      64 |      16 |          3072 |    552960 |        184320 |   12083838 |
|        1024 |      64 |      64 |          3072 |    552960 |        184320 |   12067838 |
|        1024 |     256 |       1 |          3072 |    552960 |        184320 |   34955002 |
|       10240 |      64 |      16 |         30720 |   5529600 |       1843200 |  121112974 |

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

``` r
source("scripts/duckvep_evidence.R")
root <- normalizePath(".")
revision <- duckvep_evidence_revision(root)
extension <- duckvep_evidence_build_extension(root,
  "build/release/duckhts.duckdb_extension", revision)
duckvep_evidence_write_extension_receipt("/tmp/haplotype-extension.tsv", revision, extension)
```

Pass `--extension-receipt /tmp/haplotype-extension.tsv` to the R driver.
For a network-free smoke test without a clean-build certificate:

``` bash
Rscript benchmarks/duckvep_haplotypes.R \
  --transcripts 16 --samples 4 --overlap 4 --passes 1 --diagnostic
```

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

``` bash
Rscript benchmarks/duckvep_haplotypes.R \
  --transcripts 16 --samples 4 --overlap 4 --passes 1 \
  --modes native,sql,sql_records --diagnostic
```

The source-bound raw lane uses three fresh-process passes at the current
source, with one thread pinned to CPU 2. Its timing and process RSS
appear in the timing table above; its input and output denominators are:

| transcripts | samples | overlap | input_physical_records | input_record_sample_calls | input_candidate_sample_rows | output_leaves | output_carriers |
|------------:|--------:|--------:|-----------------------:|--------------------------:|----------------------------:|--------------:|----------------:|
|        1024 |      64 |      16 |                    256 |                     16384 |                      262144 |          4096 |          131072 |

| output_leaves | cds_bytes | protein_bytes | json_bytes |
|--------------:|----------:|--------------:|-----------:|
|          4096 |    737280 |        245760 |   14052904 |

There is no earlier source-bound raw-lane measurement for a
same-contract comparison. Different output denominators preclude
treating raw versus decoded time as an implementation speedup.

## Singleton HGVS materialization

`--modes sql_singletons,sql_singletons_hgvs` compares the same decoded
input with protein HGVS disabled and enabled. Each sample carries one of
the four registered ALT events heterozygously. Four paths per transcript
retain all non-reference carriers; reference paths remain implicit. The
two substitutions and two frame-changing indels exercise the
independent-event HGVS consumer, not compound HGVS. This is a different
carrier distribution from the compound replay workload above; their
timings are not interchangeable.

The driver repeats the registered reference CDS at the model’s genomic
coordinates, pads the inter-transcript sequence with A, and indexes that
per-run FASTA. Both modes load exactly the same model and reference. A
transcript has no UTR: explicitly empty pre-/post-CDS sequences
distinguish a complete zero-length flank from an unavailable flank.
Reference/model/input preparation is outside timing. The timed query
includes sorting, bounded query initialization, sequence replay,
differences, local coding facts and every output column. HGVS-enabled
execution also opens its query-local faidx handle, retrieves genomic
context and computes protein HGVS.

Each mode has a separate correctness process, followed by three
fresh-process measurement passes with a full warm-up. Direct sequence
replacement and Biostrings translation check every expected CDS/protein.
Multiset comparisons check complete carrier identities, phase sets,
ploidies and source-event membership. The enabled path also compares
every event/transcript HGVSp with independent annotation, allowing only
its documented prediction-parenthesis difference. That internal
consistency check is not an independent VEP oracle or compound-HGVS
certificate. The four fixed event shapes do not measure broad HGVS
conformance.

Nine shared output corruptions and duplicate-output insertion must fail
the correctness process. HGVS-enabled checks additionally reject missing
HGVS and malformed internal parentheses. Whole-output byte counts and
fingerprints must match every pass of the same mode. All non-HGVS fields
must match between modes. HGVS-on/off complete outputs intentionally
differ and are never reported as identical work. The worker receipts
retain the derived FASTA/index, actual versus independent HGVS,
source/binary/input hashes, jobs and process logs. Process RSS includes
setup, warm-up and post-query aggregates, not just native workspace.

Measured source: `e6c59820dd5089fc205eadb1bbfe6d899c5138fc`. One thread,
CPU 2, Intel i5-13500, DuckDB 1.5.3.

| transcripts | samples | overlap | mode                | median_s | min_s | max_s | peak_process_rss_mib |
|------------:|--------:|--------:|:--------------------|---------:|------:|------:|---------------------:|
|        1024 |      64 |      16 | sql_singletons      |    0.343 | 0.339 | 0.351 |             409.8438 |
|       10240 |      64 |      16 | sql_singletons      |    4.132 | 4.098 | 4.224 |            2192.3711 |
|        1024 |      64 |      16 | sql_singletons_hgvs |    0.345 | 0.344 | 0.347 |             410.0391 |
|       10240 |      64 |      16 | sql_singletons_hgvs |    4.176 | 4.162 | 4.199 |            2193.1875 |

| transcripts | samples | overlap | mode                | input_physical_records | input_record_sample_calls | input_candidate_sample_rows | output_leaves | output_carriers | cds_bytes | protein_bytes | json_bytes |
|------------:|--------:|--------:|:--------------------|-----------------------:|--------------------------:|----------------------------:|--------------:|----------------:|----------:|--------------:|-----------:|
|        1024 |      64 |      16 | sql_singletons      |                    256 |                     16384 |                      262144 |          4096 |           65536 |    737280 |        244736 |    8907368 |
|        1024 |      64 |      16 | sql_singletons_hgvs |                    256 |                     16384 |                      262144 |          4096 |           65536 |    737280 |        244736 |    8915560 |
|       10240 |      64 |      16 | sql_singletons      |                   2560 |                    163840 |                     2621440 |         40960 |          655360 |   7372800 |       2447360 |   89236552 |
|       10240 |      64 |      16 | sql_singletons_hgvs |                   2560 |                    163840 |                     2621440 |         40960 |          655360 |   7372800 |       2447360 |   89318472 |

This is the first source-bound measurement of this HGVS-enabled
workload, not a before/after implementation comparison. Compound HGVS,
raw-input omissions, multi-exon models and heterogeneous genomic
contexts still need their own matched measurements. The compound replay
and raw-record baselines above retain their original source revisions
and denominators.

``` bash
Rscript benchmarks/duckvep_haplotypes.R \
  --transcripts 16 --samples 4 --overlap 4 --passes 1 \
  --modes sql_singletons,sql_singletons_hgvs --diagnostic
```
