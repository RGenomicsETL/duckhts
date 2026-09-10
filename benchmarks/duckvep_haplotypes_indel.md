Indel checkpoint: replay benchmark with unequal outputs
================

<!-- duckvep_haplotypes_indel.md is generated from duckvep_haplotypes_indel.Rmd. -->

Source: 20efcf2af33b38c5be7596e2db7b243a3be53c47. Nearest same-input
baseline: e3ec6d231cc7e5c769685761a4739c8fcffc30ff. All 51 timed
observations are retained in [the checkpoint
data](data/duckvep_haplotypes_indel.csv), including every full-output
fingerprint. **The cross-revision full-output equality check failed.**
None of these rows was appended to the existing benchmark history, and
its equality guard remains unchanged. Rendering this diagnostic report
is not a passing benchmark-promotion verdict.

## Workload and execution

[The benchmark driver](duckvep_haplotypes.R) uses the same registered
180-base canonical-DNA CDS, four original ALT events and diploid sample
templates as [the baseline report](duckvep_haplotypes.md). It checks
complete carrier multisets, sequences, edit provenance, corrupted inputs
and full output fingerprints within each mode. The singleton HGVS check
compares independent and phased DuckHTS paths; it is not an independent
VEP oracle.

Both campaigns ran on the Intel i5-13500 with DuckDB 1.5.3, one thread
pinned to CPU 2, three fresh-process timed passes and separate
correctness workers. The current extension was rebuilt from a clean
tracked checkout with HTSlib distclean; its SHA-256 is
`a5dceecf83f96b484339fec50372953465dc0543199508bb368df34a9c647f13`.
Retained receipts identify source files, native bridge, extension,
registered inputs, derived FASTA/index and every worker artifact. They
are local source-bound records, not signed CI attestations.

Native timing covers initialized ordered-feed replay and a count sink.
SQL timing covers full table materialization, including preparation,
sort, replay, sequence differences, local coding facts and all nested
output. Model/input staging and post-query fingerprints are outside
timing. Whole-process maximum RSS includes the R/DuckDB host, setup,
warm-up, output and post-query checks.

## All measured configurations

| transcripts | samples | overlap | mode                | baseline_median_s | current_median_s | current_max_rss_mib | full_output_equal | replay_equal |
|:------------|:--------|:--------|:--------------------|------------------:|-----------------:|--------------------:|:------------------|:-------------|
| 1024        | 256     | 1       | native              |            0.0475 |           0.0478 |             73.8516 | NA                | NA           |
| 1024        | 256     | 1       | sql                 |            1.3070 |           1.3110 |            957.9570 | TRUE              | TRUE         |
| 1024        | 4       | 1       | native              |            0.0032 |           0.0032 |             73.8516 | NA                | NA           |
| 1024        | 4       | 1       | sql                 |            0.0500 |           0.0500 |            218.0352 | TRUE              | TRUE         |
| 1024        | 64      | 1       | native              |            0.0138 |           0.0135 |             73.6953 | NA                | NA           |
| 1024        | 64      | 1       | sql                 |            0.3370 |           0.3380 |            376.6367 | TRUE              | TRUE         |
| 1024        | 64      | 16      | native              |            0.0177 |           0.0177 |             73.6992 | NA                | NA           |
| 1024        | 64      | 16      | sql                 |            0.3430 |           0.3480 |            378.4023 | TRUE              | TRUE         |
| 1024        | 64      | 16      | sql_records         |            0.7220 |           0.7300 |            611.9102 | TRUE              | TRUE         |
| 1024        | 64      | 16      | sql_singletons      |            0.3390 |           0.3470 |            408.9727 | FALSE             | TRUE         |
| 1024        | 64      | 16      | sql_singletons_hgvs |            0.3480 |           0.3490 |            409.5078 | FALSE             | TRUE         |
| 1024        | 64      | 64      | native              |            0.0189 |           0.0186 |             73.8555 | NA                | NA           |
| 1024        | 64      | 64      | sql                 |            0.3440 |           0.3450 |            377.8203 | TRUE              | TRUE         |
| 10240       | 64      | 16      | native              |            0.1713 |           0.1705 |             73.6992 | NA                | NA           |
| 10240       | 64      | 16      | sql                 |            4.0950 |           4.1620 |           2137.1797 | TRUE              | TRUE         |
| 10240       | 64      | 16      | sql_singletons      |            4.1070 |           4.1820 |           2190.9141 | FALSE             | TRUE         |
| 10240       | 64      | 16      | sql_singletons_hgvs |            4.1790 |           4.2270 |           2191.0586 | FALSE             | TRUE         |

`NA` marks the native count sink, which has no serialized SQL-row
fingerprint. Counts, input identities, output contracts, sequence bytes
and JSON byte lengths match between revisions. The unequal singleton
fingerprints mean those timings are **same-input, different-output
measurements**, not identical-output work. Matching byte lengths do not
establish matching values.

| transcripts | samples | overlap | mode                | input_records | projected_events | input_calls | input_physical_records | input_record_sample_calls | input_candidate_sample_rows | output_leaves | output_carriers | cds_bytes | protein_bytes | json_bytes |
|:------------|:--------|:--------|:--------------------|:--------------|:-----------------|:------------|:-----------------------|:--------------------------|:----------------------------|:--------------|:----------------|:----------|:--------------|:-----------|
| 1024        | 256     | 1       | native              | 4096          | 4096             | 1048576     | NA                     | NA                        | NA                          | 3072          | 458752          | 552960    | 184320        | NA         |
| 1024        | 256     | 1       | sql                 | NA            | NA               | NA          | NA                     | NA                        | NA                          | 3072          | 458752          | 552960    | 184320        | 35028730   |
| 1024        | 4       | 1       | native              | 4096          | 4096             | 16384       | NA                     | NA                        | NA                          | 3072          | 7168            | 552960    | 184320        | NA         |
| 1024        | 4       | 1       | sql                 | NA            | NA               | NA          | NA                     | NA                        | NA                          | 3072          | 7168            | 552960    | 184320        | 5142266    |
| 1024        | 64      | 1       | native              | 4096          | 4096             | 262144      | NA                     | NA                        | NA                          | 3072          | 114688          | 552960    | 184320        | NA         |
| 1024        | 64      | 1       | sql                 | NA            | NA               | NA          | NA                     | NA                        | NA                          | 3072          | 114688          | 552960    | 184320        | 12185338   |
| 1024        | 64      | 16      | native              | 256           | 4096             | 262144      | NA                     | NA                        | NA                          | 3072          | 114688          | 552960    | 184320        | NA         |
| 1024        | 64      | 16      | sql                 | NA            | NA               | NA          | NA                     | NA                        | NA                          | 3072          | 114688          | 552960    | 184320        | 12157566   |
| 1024        | 64      | 16      | sql_records         | NA            | NA               | NA          | 256                    | 16384                     | 262144                      | 4096          | 131072          | 737280    | 245760        | 14151208   |
| 1024        | 64      | 16      | sql_singletons      | NA            | NA               | NA          | 256                    | 16384                     | 262144                      | 4096          | 65536           | 737280    | 244736        | 9006696    |
| 1024        | 64      | 16      | sql_singletons_hgvs | NA            | NA               | NA          | 256                    | 16384                     | 262144                      | 4096          | 65536           | 737280    | 244736        | 9014888    |
| 1024        | 64      | 64      | native              | 64            | 4096             | 262144      | NA                     | NA                        | NA                          | 3072          | 114688          | 552960    | 184320        | NA         |
| 1024        | 64      | 64      | sql                 | NA            | NA               | NA          | NA                     | NA                        | NA                          | 3072          | 114688          | 552960    | 184320        | 12141566   |
| 10240       | 64      | 16      | native              | 2560          | 40960            | 2621440     | NA                     | NA                        | NA                          | 30720         | 1146880         | 5529600   | 1843200       | NA         |
| 10240       | 64      | 16      | sql                 | NA            | NA               | NA          | NA                     | NA                        | NA                          | 30720         | 1146880         | 5529600   | 1843200       | 121850254  |
| 10240       | 64      | 16      | sql_singletons      | NA            | NA               | NA          | 2560                   | 163840                    | 2621440                     | 40960         | 655360          | 7372800   | 2447360       | 90229832   |
| 10240       | 64      | 16      | sql_singletons_hgvs | NA            | NA               | NA          | 2560                   | 163840                    | 2621440                     | 40960         | 655360          | 7372800   | 2447360       | 90311752   |

## Retained output differences

All three measured passes agree within each mode. Across revisions, 12
of 51 timed rows have unequal fingerprints: the two singleton modes at
both transcript counts, with three passes each. Literal replay
fingerprints agree. The complete observations below retain the unequal
hash values; no field was normalized away to obtain a match.

| transcripts | samples | overlap | mode                | pass | field             | baseline_value           | current_value            |
|:------------|:--------|:--------|:--------------------|:-----|:------------------|:-------------------------|:-------------------------|
| 1024        | 64      | 16      | sql_singletons      | 1    | xor_hash          | 126084453999582039       | 700836657811329227       |
| 1024        | 64      | 16      | sql_singletons_hgvs | 1    | xor_hash          | 2000475952795925397      | 6726181131062349297      |
| 10240       | 64      | 16      | sql_singletons      | 1    | xor_hash          | 8822573126786218751      | 13918678999014856763     |
| 10240       | 64      | 16      | sql_singletons_hgvs | 1    | xor_hash          | 2105622761814596718      | 11764073491056183881     |
| 1024        | 64      | 16      | sql_singletons      | 1    | sum_hash          | 37582730671888802355371  | 37573251854262826388211  |
| 1024        | 64      | 16      | sql_singletons_hgvs | 1    | sum_hash          | 37824893397053031140415  | 37715148776440476629199  |
| 10240       | 64      | 16      | sql_singletons      | 1    | sum_hash          | 378298692363800329214139 | 377624590454045657951587 |
| 10240       | 64      | 16      | sql_singletons_hgvs | 1    | sum_hash          | 378148534262711101761592 | 378079038706658146159825 |
| 1024        | 64      | 16      | sql_singletons      | 1    | local_so_xor_hash | 654628605784468392       | 2335123865651189926      |
| 1024        | 64      | 16      | sql_singletons_hgvs | 1    | local_so_xor_hash | 654628605784468392       | 2335123865651189926      |
| 10240       | 64      | 16      | sql_singletons      | 1    | local_so_xor_hash | 14224601500847809748     | 5814153375303733475      |
| 10240       | 64      | 16      | sql_singletons_hgvs | 1    | local_so_xor_hash | 14224601500847809748     | 5814153375303733475      |
| 1024        | 64      | 16      | sql_singletons      | 1    | local_so_sum_hash | 37446723966248878265056  | 37469739201441226011254  |
| 1024        | 64      | 16      | sql_singletons_hgvs | 1    | local_so_sum_hash | 37446723966248878265056  | 37469739201441226011254  |
| 10240       | 64      | 16      | sql_singletons      | 1    | local_so_sum_hash | 377356556873241402705208 | 376754957532195857701877 |
| 10240       | 64      | 16      | sql_singletons_hgvs | 1    | local_so_sum_hash | 377356556873241402705208 | 376754957532195857701877 |
| 1024        | 64      | 16      | sql_singletons      | 1    | hgvs_so_xor_hash  | 14345163728266764224     | 4389172685361975006      |
| 1024        | 64      | 16      | sql_singletons_hgvs | 1    | hgvs_so_xor_hash  | 4662359414007789427      | 2986948762521831407      |
| 10240       | 64      | 16      | sql_singletons      | 1    | hgvs_so_xor_hash  | 17221334368169699363     | 18088277829535065931     |
| 10240       | 64      | 16      | sql_singletons_hgvs | 1    | hgvs_so_xor_hash  | 453159649544162334       | 17894173913969220412     |
| 1024        | 64      | 16      | sql_singletons      | 1    | hgvs_so_sum_hash  | 37210505766262944422094  | 37319589750075617092114  |
| 1024        | 64      | 16      | sql_singletons_hgvs | 1    | hgvs_so_sum_hash  | 37310727840131450206907  | 37380943966290453737923  |
| 10240       | 64      | 16      | sql_singletons      | 1    | hgvs_so_sum_hash  | 375900478728341038687991 | 376659840003214785142417 |
| 10240       | 64      | 16      | sql_singletons_hgvs | 1    | hgvs_so_sum_hash  | 379046938936219164021386 | 379285744703700591171964 |
| 1024        | 64      | 16      | sql_singletons      | 1    | non_hgvs_xor_hash | 6052207445973056023      | 18003099339576597538     |
| 1024        | 64      | 16      | sql_singletons_hgvs | 1    | non_hgvs_xor_hash | 6052207445973056023      | 18003099339576597538     |
| 10240       | 64      | 16      | sql_singletons      | 1    | non_hgvs_xor_hash | 9551608422294852071      | 2779620574718580229      |
| 10240       | 64      | 16      | sql_singletons_hgvs | 1    | non_hgvs_xor_hash | 9551608422294852071      | 2779620574718580229      |
| 1024        | 64      | 16      | sql_singletons      | 1    | non_hgvs_sum_hash | 37933509573515682488189  | 38315864159327482580440  |
| 1024        | 64      | 16      | sql_singletons_hgvs | 1    | non_hgvs_sum_hash | 37933509573515682488189  | 38315864159327482580440  |
| 10240       | 64      | 16      | sql_singletons      | 1    | non_hgvs_sum_hash | 377195884730473591625225 | 378704497758867019612841 |
| 10240       | 64      | 16      | sql_singletons_hgvs | 1    | non_hgvs_sum_hash | 377195884730473591625225 | 378704497758867019612841 |

A [retained four-event, one-transcript
reproduction](data/duckvep_haplotypes_indel_witness.jsonl.gz) isolates
the deletion `GC>G` at CDS position 67: its block
`local_consequence_mask` changes from 8,388,608 (`frameshift`) to
8,396,800 (`frameshift` plus raw `missense`). Every other field in those
four complete rows is unchanged. The physical deletion changes local
`GCT` to partial `GT`, giving peptide operands `A` and `X`. The
documented local-block contract evaluates raw coding predicates without
an independent source event’s uploaded-feature class gates. See
[ERRATA.md](../ERRATA.md) for the pinned predicate evidence and its
limits. The compact witness preserves both complete native outputs and
all original files from the exact pinned-VEP observation, including
source, observation scripts, commands, receipts and failed-attempt
diagnostics. Its 38 files are stored as byte-exact base64 payloads with
per-file hashes. Rendering verifies the bundle identity, the complete
four-row difference and the independent raw-predicate observation. It is
unsigned diagnostic evidence.

The [complete-output
audit](data/duckvep_haplotypes_indel_full_audit.jsonl.gz) retains all
eight baseline/current Parquet relations for the two singleton modes at
both transcript counts. Each untimed capture matches all 27 non-timing
result fields from its original correctness worker and all three timed
workers, including every full-output fingerprint; Parquet round trips
preserve those metrics and types. Rendering verifies the pinned bundle
and receipt identities, binds the retained metrics to both ledgers, and
checks the exact workload configuration. The [original-input
capsule](data/duckvep_haplotypes_indel_inputs.jsonl.gz) preserves 32
original worker jobs, four benchmark receipts, their four FASTA/index
pairs, and both registered fixtures. Each original receipt is bound to
its capture and ledger; every job’s actual bytes must match that
receipt, and the correctness job must also match its capture’s input
hash. Decoded jobs must contain the exact registered reference and event
table. Their complete semantic fields agree across revisions and passes,
with each pass’s correctness flag checked separately; only named binary,
output and reference-file paths are excluded from that equality.
Extension identities are checked against their own revision, not equated
across revisions. The driver remains independently pinned, and retained
FASTA/index bytes must match both receipt layers. Original paths are
identifiers only: rendering does not reopen historical files. Input
controls reject missing, malformed or altered job hashes and
self-consistently rehashed jobs with changed events or configuration.
Rendering then reruns the complete typed, exact-key comparison.

| transcripts | mode                | joined_rows | missing_keys | changed_rows | changed_other_events | unexplained_rows |
|------------:|:--------------------|------------:|-------------:|-------------:|---------------------:|-----------------:|
|        1024 | sql_singletons      |        4096 |            0 |         1024 |                    0 |                0 |
|        1024 | sql_singletons_hgvs |        4096 |            0 |         1024 |                    0 |                0 |
|       10240 | sql_singletons      |       40960 |            0 |        10240 |                    0 |                0 |
|       10240 | sql_singletons_hgvs |       40960 |            0 |        10240 |                    0 |                0 |

Across 90112 keyed row pairs, 22528 differ only in the declared deletion
mask above. There are no missing keys, changes to other events or
unexplained field changes. The audit checks the full expected
transcript/event inventory and deletion geometry, not only aggregate
hashes. Its seven corruption controls reject NULL blocks,
count-preserving key shifts, altered integer types, changed protein
sequence, omitted or duplicate rows and an incorrect deletion mask;
rendering reruns them. The bundle retains the observation/comparison
scripts, receipts, complete corrupted inputs and original logs. These
remain unsigned local observations. The outputs are still unequal:
explaining their difference does not turn these measurements into
identical-output work or change the historical promotion guard.

Three passes on one shared machine do not establish statistical
significance, a causal performance effect or a general absence of
regression. This fixture does not measure N-rich annotation,
whole-haplotype consequences or compound HGVS.

## Reproduction

Use the clean extension-receipt procedure in the baseline report, then
run `benchmarks/duckvep_haplotypes.R` with `--passes 3 --cpu 2` and
these six workloads:

| transcripts | samples | overlap | modes                                                     |
|:------------|:--------|:--------|:----------------------------------------------------------|
| 1024        | 256     | 1       | native,sql                                                |
| 1024        | 4       | 1       | native,sql                                                |
| 1024        | 64      | 1       | native,sql                                                |
| 1024        | 64      | 16      | native,sql,sql_records,sql_singletons,sql_singletons_hgvs |
| 1024        | 64      | 64      | native,sql                                                |
| 10240       | 64      | 16      | native,sql,sql_singletons,sql_singletons_hgvs             |

Pass `--extension-receipt PATH` and the listed `--transcripts`,
`--samples`, `--overlap` and `--modes` to each invocation. Preserve
every result and receipt; do not append a campaign that fails the
historical full-output check.

## Canonical workload at the translation checkpoints

Sources `52ce78513b3f5f7c2a46fe82dc44dd458333092b` and
`de3d008004f615c198a422a42a773bd7c9de0b70` each use the nearest recorded
identical workload as their timing baseline:
`20efcf2af33b38c5be7596e2db7b243a3be53c47` for `52ce785`, and `52ce785`
for `de3d008`. The workload has 1,024 transcripts, 64 samples, 16
overlapping transcripts, 256 physical events and 262,144
candidate/sample input rows. All three runs use one thread on CPU 2,
DuckDB 1.5.3 and three fresh-process timed passes. No DuckHTS test or
conformance jobs ran concurrently with the `52ce785` timers. The
`de3d008` campaign ran from 00:59:19 to 00:59:59 +0200 on September 11,
2026; a previously launched small read-only classification/oracle
diagnostic completed at 01:00:01, so concurrent diagnostic activity
cannot be excluded. This measured run is retained without replacement.
Other shared-host activity was not controlled.

| revision | baseline_revision | mode                | input_candidate_sample_rows | output_leaves | output_carriers | baseline_median_s | checkpoint_median_s | checkpoint_max_rss_mib |
|:---------|:------------------|:--------------------|----------------------------:|:--------------|:----------------|------------------:|--------------------:|-----------------------:|
| 52ce785  | 20efcf2           | native              |                      262144 | 3072          | 114688          |          0.017662 |            0.017174 |               73.69531 |
| 52ce785  | 20efcf2           | sql                 |                      262144 | 3072          | 114688          |          0.348000 |            0.346000 |              378.93750 |
| 52ce785  | 20efcf2           | sql_records         |                      262144 | 4096          | 131072          |          0.730000 |            0.731000 |              611.99219 |
| 52ce785  | 20efcf2           | sql_singletons      |                      262144 | 4096          | 65536           |          0.347000 |            0.343000 |              411.19531 |
| 52ce785  | 20efcf2           | sql_singletons_hgvs |                      262144 | 4096          | 65536           |          0.349000 |            0.353000 |              408.82812 |
| de3d008  | 52ce785           | native              |                      262144 | 3072          | 114688          |          0.017174 |            0.017235 |               73.69922 |
| de3d008  | 52ce785           | sql                 |                      262144 | 3072          | 114688          |          0.346000 |            0.344000 |              379.13672 |
| de3d008  | 52ce785           | sql_records         |                      262144 | 4096          | 131072          |          0.731000 |            0.734000 |              612.45703 |
| de3d008  | 52ce785           | sql_singletons      |                      262144 | 4096          | 65536           |          0.343000 |            0.343000 |              409.75391 |
| de3d008  | 52ce785           | sql_singletons_hgvs |                      262144 | 4096          | 65536           |          0.353000 |            0.348000 |              409.91797 |

Each checkpoint’s 15 matched rows have identical non-timing result
fields, including every SQL full-output fingerprint, output denominator
and native workspace/count metric. The [30 checkpoint
observations](data/duckvep_haplotypes_indel_translation.csv) retain
individual passes. The [264-file evidence
capsule](data/duckvep_haplotypes_indel_translation.jsonl.gz) retains all
three runs’ original worker jobs, results, logs, process-time reports,
benchmark receipts and FASTA/index pairs, all three clean-build
receipts, the registered fixtures and the unchanged driver. Rendering
verifies every retained file against its receipt before comparing actual
jobs and results. It requires complete job semantics to agree across
revisions, excluding only the explicitly checked correctness flag and
named binary/output/reference-file paths. All 60 result objects are
checked for within-mode repeatability; all timed values and RSS
measurements are bound to the corresponding ledger rows. Thirty
job-mutation controls reject changed sample counts or event alleles.
Original paths remain identifiers only; rendering does not reopen
historical files.

Compiled extension and native-bridge payloads are not included in this
capsule; their identities are retained and checked through the original
receipts. The `de3d008` extension SHA-256 is
`cb34f65f738f5d26ce5b849bbcf9628e95db180e62a199a153aacbdac44f25b0`.
These are unsigned local source-bound measurements. The historical
unequal-output comparison and its complete-output controls remain
separate and unchanged. The checkpoint comparisons establish fingerprint
equality, not complete typed-row equality.

This canonical table-1 workload measures shared replay and
singleton-HGVS paths. It does not measure the newly accepted N-bearing
REF or nonstandard-codon-table branches. Three passes on a shared host
do not establish statistical significance, causality or a general
absence of regression. Reproduce with the clean-build receipt procedure
above and this command:

``` bash
Rscript benchmarks/duckvep_haplotypes.R \
  --transcripts 1024 --samples 64 --overlap 16 --passes 3 --cpu 2 \
  --modes native,sql,sql_records,sql_singletons,sql_singletons_hgvs \
  --extension-receipt /tmp/duckhts-indel-de3d008-extension.tsv
```
