DuckVEP sorted-stream throughput
================

<!-- duckvep_throughput.md is generated from duckvep_throughput.Rmd. -->

This benchmark drives either the rich `duckvep_annotate` result or the
numeric `duckvep_annotate_compact` result through DuckDB’s stable
extension API with nondecreasing sequence-region and position keys. It
records candidate sweep, consequence classification, list
materialization, `unnest`, and final aggregation together. Model loading
and input staging are outside the timed pass.

The checked-in one-transcript workload is a repeatable adapter floor.
The final production-density workload uses the receipt-hashed Ensembl
116 GRCh38 model: 194 sequence regions, 644,427 transcripts, 5,068,416
exon memberships, and 369,631 sequence-backed coding transcripts. Its
ordinary topology workload is a deterministic 1-in-40 hash sample of
100,957 sorted GIAB sites. The coding workloads deliberately repeat a
frozen set to stabilize one-core kernel timing: 200,000 SNVs repeat
1,000 distinct alleles; the 36,000 non-SNV workload repeats 1,000 MNVs
and 8,000 length-changing alleles four times; the 200,000 mixed workload
is 10% SNV, 10% MNV, 40% lengthening, and 40% shortening. Repetition
favours warm model and sequence caches, so these rows measure the
resident execution lanes, not file ingestion or an end-to-end genome
run.

Earlier chromosome-only and complete-primary staging rows remain in the
ledger as historical measurements. They have different transcript
inventories and checksums and must not be treated as revisions of the
final receipt. The production query restates `ORDER BY` in the timed
input CTE, so its result includes the cost of enforcing the kernel’s
nondecreasing-input contract even though the staging table was already
written in that order. The `79cfaaf2` rows predate that review
correction and relied on the single-thread physical scan preserving
insertion order; they remain only as historical measurements.

## Recorded runs

| run_date   | revision | workload                                    | output_mode | threads | variants   | transcripts | exons     | annotated_rows | passes | min_seconds | median_seconds | max_seconds | variants_per_second | ns_per_variant | cpu                                 | output_rows_per_second |
|:-----------|:---------|:--------------------------------------------|:------------|--------:|:-----------|:------------|:----------|:---------------|-------:|------------:|---------------:|------------:|--------------------:|---------------:|:------------------------------------|:-----------------------|
| 2026-07-11 | 8cc22218 | fixture_one_transcript_sorted               | rich        |       1 | 10,000,000 | 1           | 2         | 10,000,000     |      3 |       3.225 |          3.334 |       3.371 |             2999400 |          333.4 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,999,400              |
| 2026-07-13 | 87f03a2a | fixture_one_transcript_sorted               | rich        |       1 | 10,000,000 | 1           | 2         | 10,000,000     |      3 |       2.812 |          2.824 |       2.825 |             3541076 |          282.4 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,541,076              |
| 2026-07-14 | 751d3c74 | ensembl116_grch38_giab_sites_hash40         | compact     |       1 | 100,957    | 644,292     | 5,078,384 | 1,179,465      |      9 |       0.104 |          0.107 |       0.114 |              943523 |         1059.9 | 13th Gen Intel(R) Core(TM) i5-13500 | 11,023,037             |
| 2026-07-14 | 751d3c74 | ensembl116_grch38_giab_sites_hash40         | rich        |       1 | 100,957    | 644,292     | 5,078,384 | 1,179,465      |      9 |       0.222 |          0.224 |       0.231 |              450701 |         2218.8 | 13th Gen Intel(R) Core(TM) i5-13500 | 5,265,469              |
| 2026-07-14 | 79cfaaf2 | ensembl116_grch38_giab_sites_hash40         | compact     |       1 | 100,957    | 644,292     | 5,078,384 | 1,179,465      |      9 |       0.107 |          0.108 |       0.110 |              934787 |         1069.8 | 13th Gen Intel(R) Core(TM) i5-13500 | 10,920,972             |
| 2026-07-14 | 79cfaaf2 | ensembl116_grch38_giab_sites_hash40         | rich        |       1 | 100,957    | 644,292     | 5,078,384 | 1,179,465      |      9 |       0.225 |          0.228 |       0.232 |              442794 |         2258.4 | 13th Gen Intel(R) Core(TM) i5-13500 | 5,173,092              |
| 2026-07-14 | 81b9e2a7 | ensembl116_grch38_primary_giab_sites_hash40 | compact     |       1 | 100,957    | 646,577     | 5,087,789 | 1,179,465      |      9 |       0.112 |          0.115 |       0.125 |              877887 |         1139.1 | 13th Gen Intel(R) Core(TM) i5-13500 | 10,256,217             |
| 2026-07-14 | 81b9e2a7 | ensembl116_grch38_primary_giab_sites_hash40 | rich        |       1 | 100,957    | 646,577     | 5,087,789 | 1,179,465      |      9 |       0.233 |          0.236 |       0.242 |              427784 |         2337.6 | 13th Gen Intel(R) Core(TM) i5-13500 | 4,997,733              |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_coding_mixed_200k   | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,757,180      |      3 |      58.833 |         59.010 |      59.212 |                3389 |       295050.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 97,563                 |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_coding_nonsnv_36k   | compact     |       1 | 36,000     | 644,427     | 5,068,416 | 1,047,616      |      5 |      11.837 |         11.898 |      11.959 |                3026 |       330500.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 88,050                 |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_coding_nonsnv_36k   | rich        |       1 | 36,000     | 644,427     | 5,068,416 | 1,047,616      |      3 |      11.903 |         11.952 |      11.978 |                3012 |       332000.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 87,652                 |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_coding_snv_200k     | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,191,000      |      9 |       0.690 |          0.693 |       0.719 |              288600 |         3465.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 7,490,620              |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_coding_snv_200k     | rich        |       1 | 200,000    | 644,427     | 5,068,416 | 5,191,000      |      5 |       1.213 |          1.238 |       1.262 |              161551 |         6190.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 4,193,053              |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_giab_sites_hash40   | compact     |       1 | 100,957    | 644,427     | 5,068,416 | 1,174,245      |      9 |       0.116 |          0.117 |       0.119 |              862880 |         1158.9 | 13th Gen Intel(R) Core(TM) i5-13500 | 10,036,282             |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_giab_sites_hash40   | rich        |       1 | 100,957    | 644,427     | 5,068,416 | 1,174,245      |      9 |       0.239 |          0.244 |       0.263 |              413758 |         2416.9 | 13th Gen Intel(R) Core(TM) i5-13500 | 4,812,480              |
| 2026-07-14 | bcca6f6c | ensembl116_grch38_giab_sites_hash40         | compact     |       1 | 100,957    | 644,292     | 5,078,384 | 1,179,465      |      9 |       0.106 |          0.107 |       0.109 |              943523 |         1059.9 | 13th Gen Intel(R) Core(TM) i5-13500 | 11,023,037             |
| 2026-07-14 | bcca6f6c | ensembl116_grch38_giab_sites_hash40         | rich        |       1 | 100,957    | 644,292     | 5,078,384 | 1,179,465      |      9 |       0.221 |          0.225 |       0.231 |              448698 |         2228.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 5,242,067              |
| 2026-07-15 | 1903f7fb | ensembl116_grch38_final_coding_mixed_200k   | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,756,720      |     21 |       1.489 |          1.514 |       1.618 |              132100 |         7570.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,802,325              |
| 2026-07-15 | 1903f7fb | ensembl116_grch38_final_coding_nonsnv_36k   | compact     |       1 | 36,000     | 644,427     | 5,068,416 | 1,047,524      |     31 |       0.323 |          0.329 |       0.362 |              109422 |         9138.9 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,183,964              |
| 2026-07-15 | 1903f7fb | ensembl116_grch38_final_coding_snv_200k     | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,191,000      |     31 |       0.610 |          0.622 |       0.657 |              321543 |         3110.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 8,345,659              |
| 2026-07-15 | 57b5b566 | ensembl116_grch38_final_coding_mixed_200k   | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,756,720      |     15 |       1.467 |          1.480 |       1.489 |              135135 |         7400.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,889,676              |
| 2026-07-15 | 57b5b566 | ensembl116_grch38_final_coding_nonsnv_36k   | compact     |       1 | 36,000     | 644,427     | 5,068,416 | 1,047,524      |     21 |       0.310 |          0.322 |       0.351 |              111801 |         8944.4 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,253,180              |
| 2026-07-15 | 57b5b566 | ensembl116_grch38_final_coding_snv_200k     | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,191,000      |     15 |       0.609 |          0.615 |       0.638 |              325203 |         3075.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 8,440,650              |
| 2026-07-15 | 7dd90ce8 | ensembl116_grch38_final_coding_mixed_200k   | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,756,720      |      9 |       1.550 |          1.572 |       1.583 |              127226 |         7860.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,662,036              |
| 2026-07-15 | 7dd90ce8 | ensembl116_grch38_final_coding_nonsnv_36k   | compact     |       1 | 36,000     | 644,427     | 5,068,416 | 1,047,524      |     11 |       0.341 |          0.344 |       0.357 |              104651 |         9555.6 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,045,128              |
| 2026-07-15 | 7dd90ce8 | ensembl116_grch38_final_coding_snv_200k     | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,191,000      |      9 |       0.609 |          0.613 |       0.629 |              326264 |         3065.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 8,468,189              |
| 2026-07-15 | 7dd90ce8 | ensembl116_grch38_final_giab_sites_hash40   | compact     |       1 | 100,957    | 644,427     | 5,068,416 | 1,174,245      |      9 |       0.118 |          0.119 |       0.137 |              848378 |         1178.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 9,867,605              |
| 2026-07-15 | a346a7aa | ensembl116_grch38_final_coding_mixed_200k   | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,756,720      |     21 |       1.417 |          1.443 |       1.511 |              138600 |         7215.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,989,411              |
| 2026-07-15 | a346a7aa | ensembl116_grch38_final_coding_nonsnv_36k   | compact     |       1 | 36,000     | 644,427     | 5,068,416 | 1,047,524      |     31 |       0.307 |          0.315 |       0.376 |              114286 |         8750.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,325,473              |
| 2026-07-15 | a346a7aa | ensembl116_grch38_final_coding_snv_200k     | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,191,000      |     31 |       0.599 |          0.617 |       0.791 |              324149 |         3085.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 8,413,290              |
| 2026-07-15 | c361346f | ensembl116_grch38_final_coding_mixed_200k   | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,756,720      |     21 |       1.445 |          1.468 |       1.498 |              136240 |         7340.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,921,471              |
| 2026-07-15 | c361346f | ensembl116_grch38_final_coding_nonsnv_36k   | compact     |       1 | 36,000     | 644,427     | 5,068,416 | 1,047,524      |     31 |       0.320 |          0.323 |       0.341 |              111455 |         8972.2 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,243,108              |
| 2026-07-15 | c361346f | ensembl116_grch38_final_coding_snv_200k     | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,191,000      |     31 |       0.584 |          0.602 |       0.779 |              332226 |         3010.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 8,622,924              |
| 2026-07-15 | ecb7eaae | ensembl116_grch38_final_coding_mixed_200k   | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,756,720      |     15 |       1.467 |          1.495 |       1.512 |              133779 |         7475.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,850,649              |
| 2026-07-15 | ecb7eaae | ensembl116_grch38_final_coding_nonsnv_36k   | compact     |       1 | 36,000     | 644,427     | 5,068,416 | 1,047,524      |     21 |       0.311 |          0.324 |       0.356 |              111111 |         9000.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,233,099              |
| 2026-07-15 | ecb7eaae | ensembl116_grch38_final_coding_snv_200k     | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,191,000      |     15 |       0.605 |          0.610 |       0.624 |              327869 |         3050.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 8,509,836              |
| 2026-07-15 | f649d724 | ensembl116_grch38_final_coding_mixed_200k   | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,756,720      |      9 |       1.487 |          1.504 |       1.626 |              132979 |         7520.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,827,606              |
| 2026-07-15 | f649d724 | ensembl116_grch38_final_coding_nonsnv_36k   | compact     |       1 | 36,000     | 644,427     | 5,068,416 | 1,047,524      |     11 |       0.324 |          0.327 |       0.342 |              110092 |         9083.3 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,203,437              |
| 2026-07-15 | f649d724 | ensembl116_grch38_final_coding_snv_200k     | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,191,000      |     11 |       0.611 |          0.613 |       0.634 |              326264 |         3065.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 8,468,189              |
| 2026-07-15 | f649d724 | ensembl116_grch38_final_giab_sites_hash40   | compact     |       1 | 100,957    | 644,427     | 5,068,416 | 1,174,245      |     15 |       0.118 |          0.120 |       0.132 |              841308 |         1188.6 | 13th Gen Intel(R) Core(TM) i5-13500 | 9,785,375              |
| 2026-07-16 | 9c97cb07 | ensembl116_grch38_final_coding_mixed_200k   | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,756,720      |     15 |       1.451 |          1.481 |       1.638 |              135044 |         7405.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,887,049              |
| 2026-07-16 | 9c97cb07 | ensembl116_grch38_final_coding_nonsnv_36k   | compact     |       1 | 36,000     | 644,427     | 5,068,416 | 1,047,524      |     15 |       0.317 |          0.322 |       0.339 |              111801 |         8944.4 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,253,180              |
| 2026-07-16 | 9c97cb07 | ensembl116_grch38_final_coding_snv_200k     | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,191,000      |     15 |       0.593 |          0.603 |       0.709 |              331675 |         3015.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 8,608,624              |
| 2026-07-16 | 9c97cb07 | ensembl116_grch38_final_giab_sites_hash40   | compact     |       1 | 100,957    | 644,427     | 5,068,416 | 1,174,245      |     15 |       0.117 |          0.118 |       0.131 |              855568 |         1168.8 | 13th Gen Intel(R) Core(TM) i5-13500 | 9,951,229              |

Each pass consumes every staged input and checks output cardinality plus
either the rendered consequence-byte total or the numeric
consequence-mask sum. Rows with different workload, output mode, thread
count, host, or variant count are separate measurements and should not
be compared as if only the engine changed.

## Validated CDS projection cache

The `7dd90ce8` extension source is byte-identical to the immediate
pre-cache parent; the intervening commits changed only documentation and
benchmark data. The `f649d724` rows measure the cache with the same
final model, input tables, single pinned CPU, ordered query, compact
output, and warmup size. Positive change means less elapsed time. The
non-SNV and mixed lanes avoid repeated binary searches for both CDS
endpoints. The sorted SNV and ordinary GIAB lanes show the control cost
because they rarely or never consume that cache.

| workload       | variants | output rows | seconds before | variants/s before | seconds after | variants/s after | elapsed change (%) |
|:---------------|:---------|:------------|---------------:|------------------:|--------------:|-----------------:|-------------------:|
| mixed coding   | 200,000  | 5,756,720   |          1.572 |            127226 |         1.504 |           132979 |                4.3 |
| non-SNV coding | 36,000   | 1,047,524   |          0.344 |            104651 |         0.327 |           110092 |                4.9 |
| coding SNV     | 200,000  | 5,191,000   |          0.613 |            326264 |         0.613 |           326264 |                0.0 |
| GIAB sites     | 100,957  | 1,174,245   |          0.119 |            848378 |         0.120 |           841308 |               -0.8 |

## Reused intron coordinates for mismatch islands

The `ecb7eaae` rows are the merged projection-cache source immediately
before the splice-classifier change. The `57b5b566` rows reuse the
intron coordinates already resolved for each exon pair, avoid calling
the shared predicate when neither the intron nor a splice or
short-intron window can contribute, and keep that predicate inlined in
the sorted point classifier. Both revisions use the same final model,
input tables, pinned CPU, ordered query, compact output, and
10,000-variant warmup. Positive change means less elapsed time.

| workload           | variants | output rows | seconds before | variants/s before | seconds after | variants/s after | elapsed change (%) |
|:-------------------|:---------|:------------|---------------:|------------------:|--------------:|-----------------:|-------------------:|
| non-SNV coding     | 36,000   | 1,047,524   |          0.324 |            111111 |         0.322 |           111801 |                0.6 |
| mixed coding       | 200,000  | 5,756,720   |          1.495 |            133779 |         1.480 |           135135 |                1.0 |
| coding SNV control | 200,000  | 5,191,000   |          0.610 |            327869 |         0.615 |           325203 |               -0.8 |

The coding-SNV row is a control: that route does not use mismatch-island
classification. Its 5 ms movement is reported as measurement noise
rather than attributed to this change.

## Reused coding projection fact for NMD

The `1903f7fb` rows are the exact merged source before this change. The
`a346a7aa` rows record the first projection-cache optimization: retain
the `CDS end <= 101` comparison already established by the sequence
classifier and consume it when VEP 116 NMD prediction needs the same
fact. That revision is a measured historical intermediate. The later
executable-plugin audit proved that the fact is interchangeable only
when the minimized edit and VEP feature cover exactly the same genomic
bases. Both revisions below were rebuilt and measured adjacently on the
same final model, input tables, pinned CPU, ordered query, compact
output, and 10,000-variant warmup. Positive change means less elapsed
time.

| workload           | variants | output rows | seconds before | variants/s before | seconds after | variants/s after | elapsed change (%) |
|:-------------------|:---------|:------------|---------------:|------------------:|--------------:|-----------------:|-------------------:|
| mixed coding       | 200,000  | 5,756,720   |          1.514 |            132100 |         1.443 |           138600 |                4.7 |
| non-SNV coding     | 36,000   | 1,047,524   |          0.329 |            109422 |         0.315 |           114286 |                4.3 |
| coding SNV control | 200,000  | 5,191,000   |          0.622 |            321543 |         0.617 |           324149 |                0.8 |

The coding-SNV route does not carry the reused projection fact. Its 5 ms
movement is reported as a control rather than attributed to the NMD
change.

## Cost of exact VEP feature projection

The `c361346f` correction retains the cached fact only for identical
genomic spans. Wider equal-length features and insertions instead
project the same parent `TranscriptVariation` state read by `NMD.pm`.
These rows compare that exact implementation with the earlier over-broad
cache revision. Output row counts and consequence-mask checksums are
identical for every workload. Negative elapsed change means the exact
implementation took longer; the SNV lane is a control because it does
not use this NMD path.

| workload           | variants | output rows | seconds before | variants/s before | seconds exact | variants/s exact | elapsed change (%) |
|:-------------------|:---------|:------------|---------------:|------------------:|--------------:|-----------------:|-------------------:|
| mixed coding       | 200,000  | 5,756,720   |          1.443 |            138600 |         1.468 |           136240 |               -1.7 |
| non-SNV coding     | 36,000   | 1,047,524   |          0.315 |            114286 |         0.323 |           111455 |               -2.5 |
| coding SNV control | 200,000  | 5,191,000   |          0.617 |            324149 |         0.602 |           332226 |                2.4 |

## Final GRCh38 model resources

| measure                                | value         |
|:---------------------------------------|:--------------|
| revision                               | b204dd49      |
| model artifact                         | 1.608 GiB     |
| model load                             | 3.047 s       |
| process RSS before load                | 0.114 GiB     |
| process RSS after load                 | 4.125 GiB     |
| peak process RSS                       | 5.041 GiB     |
| process RSS after drop                 | 2.834 GiB     |
| DuckDB base-table cache                | 1.594 GiB     |
| planned scratch per worker             | 3.208 MiB     |
| observed RSS change on first workspace | 0.156 MiB     |
| largest transcript run                 | 62,116        |
| largest prepared CDS                   | 107,976 bases |

The artifact size is a compact DuckDB containing only `model_regions`,
`model_transcripts`, and `model_receipt`; benchmark inputs are excluded.
RSS is for one R/DuckDB process, not a C-heap estimate. The load peak
includes both the DuckDB scan/cache and construction of the second
immutable C representation. Dropping the model releases the C
representation, but DuckDB retains its table cache in the process. The
planned per-worker number is the exact sum of the C buffer capacities
(active and candidate transcript indices, one exon cursor per
transcript, and worst-case edit/CDS/peptide scratch); allocator metadata
is not included. The much smaller first-workspace RSS change reflects
lazy or reused pages, not a smaller capacity.

The standalone model ABI 0.12 cache adds three `uint32_t` arrays and one
`uint8_t` array per transcript. Its exact element payload for this model
is 8,377,551 bytes (7.989 MiB), plus four allocator records. This
deterministic allocation is not retroactively included in the older
`b204dd49` RSS measurement above.

The production form is run with:

``` sh
Rscript benchmarks/duckvep_throughput.R \
  --database /path/to/staged.duckdb \
  --variants-table bench_variants \
  --workload-name ensembl116_grch38_final_giab_sites_hash40 \
  --variants 100957 --passes 9 --warmup 10000 --threads 1 \
  --output compact
```
