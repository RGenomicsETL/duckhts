Phased replay: native stream and public SQL
================

<!-- duckvep_haplotypes.md is generated from duckvep_haplotypes.Rmd. -->

Compound-replay source: 677d684810a78ba2749e5958954e7f745eb7d278;
same-input baseline: c69e0e0cf124c65cc8e23c25ca80691bd446a2cd. Native
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
current output is always checked. Contracts ending in
`_nominal_length_diff` include the exact signed sum of source
replacement-length changes. The verifier calculates that sum from
fixture alleles, independently of replayed CDS and block lengths, and
rejects a +3 corruption that preserves frame category. Historical
contracts do not certify this numeric field.

`hgvs_so_*` fingerprints retain every local-SO/HGVS-schema field and
exclude only `nominal_length_diff`. `local_so_*` additionally excludes
`hgvsp` and `hgvsp_status`; `replay_*` also projects each block to its
literal replay fields. `non_hgvs_*` retains nominal length for the
current HGVS-on/off comparison. The recorded contract determines which
fields were measured; historical full fingerprints are used only for
projections identical to their recorded schema.

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
|        1024 |       4 |       1 | native      | 0.003 |    0.003 | 0.003 |              73.699 |
|        1024 |       4 |       1 | sql         | 0.051 |    0.051 | 0.052 |             217.879 |
|        1024 |      64 |       1 | native      | 0.014 |    0.014 | 0.014 |              73.699 |
|        1024 |      64 |       1 | sql         | 0.335 |    0.337 | 0.337 |             377.520 |
|        1024 |      64 |      16 | native      | 0.017 |    0.017 | 0.017 |              73.852 |
|        1024 |      64 |      16 | sql         | 0.346 |    0.352 | 0.360 |             378.555 |
|        1024 |      64 |      16 | sql_records | 0.725 |    0.733 | 0.738 |             612.230 |
|        1024 |      64 |      64 | native      | 0.019 |    0.019 | 0.019 |              73.699 |
|        1024 |      64 |      64 | sql         | 0.343 |    0.345 | 0.350 |             378.305 |
|        1024 |     256 |       1 | native      | 0.045 |    0.046 | 0.048 |              73.695 |
|        1024 |     256 |       1 | sql         | 1.299 |    1.303 | 1.303 |             958.102 |
|       10240 |      64 |      16 | native      | 0.170 |    0.172 | 0.173 |              73.852 |
|       10240 |      64 |      16 | sql         | 4.100 |    4.102 | 4.106 |            2143.492 |

Both compared revisions return each block’s local SO mask, coding status
and position relative to the first stop. A shared coding context
evaluates physical blocks on completed replay. Full-output fingerprints
and input/output denominators match across the compared configurations.
The native count sink omits local SO evaluation. HGVS-enabled
performance requires a separate workload.

| transcripts | samples | overlap | median_s_before | median_s_after | median_change_percent | max_process_rss_mib_before | max_process_rss_mib_after |
|------------:|--------:|--------:|----------------:|---------------:|----------------------:|---------------------------:|--------------------------:|
|        1024 |       4 |       1 |           0.055 |          0.051 |                -7.273 |                    219.172 |                   217.879 |
|        1024 |      64 |       1 |           0.364 |          0.337 |                -7.418 |                    378.012 |                   377.520 |
|        1024 |      64 |      16 |           0.372 |          0.352 |                -5.376 |                    378.488 |                   378.555 |
|        1024 |      64 |      64 |           0.379 |          0.345 |                -8.971 |                    377.684 |                   378.305 |
|        1024 |     256 |       1 |           1.467 |          1.303 |               -11.179 |                    957.824 |                   958.102 |
|       10240 |      64 |      16 |           4.841 |          4.102 |               -15.265 |                   2135.391 |                  2143.492 |

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
|        1024 |       4 |       1 |          3072 |    552960 |        184320 |    5142266 |
|        1024 |      64 |       1 |          3072 |    552960 |        184320 |   12185338 |
|        1024 |      64 |      16 |          3072 |    552960 |        184320 |   12157566 |
|        1024 |      64 |      64 |          3072 |    552960 |        184320 |   12141566 |
|        1024 |     256 |       1 |          3072 |    552960 |        184320 |   35028730 |
|       10240 |      64 |      16 |         30720 |   5529600 |       1843200 |  121850254 |

CDS/protein bytes count the sequence views returned once per occupied
leaf. JSON bytes are the measured UTF-8 size of complete canonicalized
SQL rows, including carriers, contributors, blocks and differences; they
are not bytes written to disk.

## Sparse native state

| transcripts | samples | overlap | peak_transcripts | peak_carriers | peak_prefixes | peak_events | peak_projections | peak_allele_bytes | workspace_bytes | model_bytes |
|------------:|--------:|--------:|-----------------:|--------------:|--------------:|------------:|-----------------:|------------------:|----------------:|------------:|
|        1024 |       4 |       1 |                1 |             7 |             6 |           4 |                4 |                10 |            4944 |       34136 |
|        1024 |      64 |       1 |                1 |           112 |             6 |           4 |                4 |                10 |           14544 |       34136 |
|        1024 |      64 |      16 |               16 |          1792 |            96 |           4 |               64 |                10 |          189564 |       34136 |
|        1024 |      64 |      64 |               64 |          7168 |           384 |           4 |              256 |                10 |          749628 |       34136 |
|        1024 |     256 |       1 |                1 |           448 |             6 |           4 |                4 |                10 |           45264 |       34136 |
|       10240 |      64 |      16 |               16 |          1792 |            96 |           4 |               64 |                10 |          189564 |      338264 |

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

This same-input comparison retains three repeated passes per workload
and revision. The `c69e0e0` baseline overlapped reference and model
conformance jobs pinned to other CPUs. The current campaign did not
deliberately overlap conformance work. Both ran at different times on a
shared machine. CPU pinning does not isolate memory, thermal or
scheduling effects, so timing differences do not establish the cost of
the code change. One machine and this deliberately shared synthetic
cohort do not establish production throughput, statistical significance
or a general absence of regression. Whole-haplotype SO/HGVS and typed
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
the fixture. Corruption controls cover changed CDS, protein, missing
CDS, carrier identity, missing ploidy, contributor identity, missing ALT
ordinal, invented HGVS, changed HGVS status, invalid source GT and
duplicate output. They run in the separate correctness worker. Contracts
with `nominal_length_diff` also reject a signed length change of +3
while all sequence flags remain untouched.

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
|          4096 |    737280 |        245760 |   14151208 |

The raw-lane comparison uses the same source revisions as the decoded
comparison. Input hashes, complete output fingerprints, all listed
denominators, thread count and output contract match between revisions.
Different output denominators preclude treating raw versus decoded time
as an implementation speedup.

| transcripts | samples | overlap | median_s_before | median_s_after | median_change_percent | max_process_rss_mib_before | max_process_rss_mib_after |
|------------:|--------:|--------:|----------------:|---------------:|----------------------:|---------------------------:|--------------------------:|
|        1024 |      64 |      16 |           0.897 |          0.733 |               -18.283 |                    612.504 |                    612.23 |

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

Shared output corruptions and duplicate-output insertion must fail the
correctness process. Nominal-length contracts also reject a same-frame
+3 numeric corruption. HGVS-enabled checks additionally reject missing
HGVS and malformed internal parentheses. Whole-output byte counts and
fingerprints must match every pass of the same mode. All non-HGVS fields
must match between modes. HGVS-on/off complete outputs intentionally
differ and are never reported as identical work. The worker receipts
retain the derived FASTA/index, actual versus independent HGVS,
source/binary/input hashes, jobs and process logs. Process RSS includes
setup, warm-up and post-query aggregates, not just native workspace.

Measured source: `677d684810a78ba2749e5958954e7f745eb7d278`. One thread,
CPU 2, Intel i5-13500, DuckDB 1.5.3.

| transcripts | samples | overlap | mode                | median_s | min_s | max_s | peak_process_rss_mib |
|------------:|--------:|--------:|:--------------------|---------:|------:|------:|---------------------:|
|        1024 |      64 |      16 | sql_singletons      |    0.340 | 0.339 | 0.344 |             410.9961 |
|       10240 |      64 |      16 | sql_singletons      |    4.131 | 4.127 | 4.254 |            2190.8164 |
|        1024 |      64 |      16 | sql_singletons_hgvs |    0.347 | 0.345 | 0.347 |             408.7773 |
|       10240 |      64 |      16 | sql_singletons_hgvs |    4.207 | 4.206 | 4.210 |            2191.1172 |

| transcripts | samples | overlap | mode                | input_physical_records | input_record_sample_calls | input_candidate_sample_rows | output_leaves | output_carriers | cds_bytes | protein_bytes | json_bytes |
|------------:|--------:|--------:|:--------------------|-----------------------:|--------------------------:|----------------------------:|--------------:|----------------:|----------:|--------------:|-----------:|
|        1024 |      64 |      16 | sql_singletons      |                    256 |                     16384 |                      262144 |          4096 |           65536 |    737280 |        244736 |    9006696 |
|        1024 |      64 |      16 | sql_singletons_hgvs |                    256 |                     16384 |                      262144 |          4096 |           65536 |    737280 |        244736 |    9014888 |
|       10240 |      64 |      16 | sql_singletons      |                   2560 |                    163840 |                     2621440 |         40960 |          655360 |   7372800 |       2447360 |   90229832 |
|       10240 |      64 |      16 | sql_singletons_hgvs |                   2560 |                    163840 |                     2621440 |         40960 |          655360 |   7372800 |       2447360 |   90311752 |

Nearest same-input source: `c69e0e0cf124c65cc8e23c25ca80691bd446a2cd`.
Full-output fingerprints, input identities and every listed denominator
match within each mode.

| transcripts | samples | overlap | mode                | median_s_before | min_s_before | max_s_before | peak_process_rss_mib_before | median_s_after | min_s_after | max_s_after | peak_process_rss_mib_after | median_change_percent |
|------------:|--------:|--------:|:--------------------|----------------:|-------------:|-------------:|----------------------------:|---------------:|------------:|------------:|---------------------------:|----------------------:|
|        1024 |      64 |      16 | sql_singletons      |           0.397 |        0.396 |        0.406 |                    411.0273 |          0.340 |       0.339 |       0.344 |                   410.9961 |              -14.3577 |
|        1024 |      64 |      16 | sql_singletons_hgvs |           0.399 |        0.398 |        0.412 |                    409.7812 |          0.347 |       0.345 |       0.347 |                   408.7773 |              -13.0326 |
|       10240 |      64 |      16 | sql_singletons      |           5.229 |        4.837 |        6.093 |                   2192.9805 |          4.131 |       4.127 |       4.254 |                  2190.8164 |              -20.9983 |
|       10240 |      64 |      16 | sql_singletons_hgvs |           5.168 |        5.089 |        5.451 |                   2191.1953 |          4.207 |       4.206 |       4.210 |                  2191.1172 |              -18.5952 |

Singleton HGVS borrows the decoded stream’s successful physical CDS
projection. VEP’s uploaded-feature interpretation and shifted-HGVS
placement remain distinct; raw-source HGVS prepares and projects its own
minimized allele. Leaf and contributor storage is included in the query
workspace limit. These whole-process measurements do not isolate that
storage cost or the projection cost. Three passes on one machine do not
establish statistical significance or a general absence of regression.
Compound HGVS, raw-input omissions, multi-exon models and heterogeneous
genomic contexts still need their own matched measurements. The compound
replay and raw-record baselines above retain their original source
revisions and denominators.

``` bash
Rscript benchmarks/duckvep_haplotypes.R \
  --transcripts 16 --samples 4 --overlap 4 --passes 1 \
  --modes sql_singletons,sql_singletons_hgvs --diagnostic
```
