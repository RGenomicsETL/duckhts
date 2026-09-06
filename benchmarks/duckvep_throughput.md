DuckVEP sorted-stream throughput
================

<!-- duckvep_throughput.md is generated from duckvep_throughput.Rmd. -->

This benchmark drives the internal rich `_duckvep_annotate_small_rich`,
numeric `_duckvep_annotate_small_compact`, independent-event
`_duckvep_annotate_small_hgvs`, or fused rich-plus-HGVS result through
DuckDB’s stable extension API with nondecreasing sequence-region and
position keys. It records candidate sweep, consequence classification,
list materialization, `unnest`, and final aggregation together. Model
loading and input staging are outside the timed pass.

The checked-in one-transcript workload is a repeatable adapter floor.
The final production-density workload uses the receipt-hashed Ensembl
116 GRCh38 model: 194 sequence regions, 644,427 transcripts, 5,068,416
exon memberships, and 369,631 sequence-backed coding transcripts. The
regulation-enabled form adds 380,818 RegulatoryFeature and 1,002,762
MotifFeature intervals in an independent resident SoA. Its ordinary
topology workload is a deterministic 1-in-40 hash sample of 100,957
sorted GIAB sites. The coding workloads deliberately repeat a frozen set
to stabilize one-core kernel timing: 200,000 SNVs repeat 1,000 distinct
alleles; the 36,000 non-SNV workload repeats 1,000 MNVs and 8,000
length-changing alleles four times; the 200,000 mixed workload is 10%
SNV, 10% MNV, 40% lengthening, and 40% shortening. Repetition favours
warm model and sequence caches, so these rows measure the resident
execution lanes, not file ingestion or an end-to-end genome run.

The paired-breakend workload repeats the 1,004-event multichromosome
seed-31 VEP differential to 200,000 sorted endpoint pairs. Both
endpoints participate in candidate discovery and each transcript
receives the union of the two endpoint consequence sets. Its much larger
transcript fan-out makes output rows/s as important as paired events/s.

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

| run_date   | revision | workload                                                         | output_mode | threads | input_partitions | transcript_distance | variants   | transcripts | exons     | regulation_features | annotated_rows | passes | min_seconds | median_seconds | max_seconds | variants_per_second | ns_per_variant | cpu                                 | cpu_affinity | output_rows_per_second |
|:-----------|:---------|:-----------------------------------------------------------------|:------------|--------:|-----------------:|--------------------:|:-----------|:------------|:----------|:--------------------|:---------------|-------:|------------:|---------------:|------------:|--------------------:|---------------:|:------------------------------------|:-------------|:-----------------------|
| 2026-07-11 | 8cc22218 | fixture_one_transcript_sorted                                    | rich        |       1 |                1 |                5000 | 10,000,000 | 1           | 2         | 0                   | 10,000,000     |      3 |       3.225 |          3.334 |       3.371 |             2999400 |          333.4 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 2,999,400              |
| 2026-07-13 | 87f03a2a | fixture_one_transcript_sorted                                    | rich        |       1 |                1 |                5000 | 10,000,000 | 1           | 2         | 0                   | 10,000,000     |      3 |       2.812 |          2.824 |       2.825 |             3541076 |          282.4 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,541,076              |
| 2026-07-14 | 751d3c74 | ensembl116_grch38_giab_sites_hash40                              | compact     |       1 |                1 |                5000 | 100,957    | 644,292     | 5,078,384 | 0                   | 1,179,465      |      9 |       0.104 |          0.107 |       0.114 |              943523 |         1059.9 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 11,023,037             |
| 2026-07-14 | 751d3c74 | ensembl116_grch38_giab_sites_hash40                              | rich        |       1 |                1 |                5000 | 100,957    | 644,292     | 5,078,384 | 0                   | 1,179,465      |      9 |       0.222 |          0.224 |       0.231 |              450701 |         2218.8 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 5,265,469              |
| 2026-07-14 | 79cfaaf2 | ensembl116_grch38_giab_sites_hash40                              | compact     |       1 |                1 |                5000 | 100,957    | 644,292     | 5,078,384 | 0                   | 1,179,465      |      9 |       0.107 |          0.108 |       0.110 |              934787 |         1069.8 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 10,920,972             |
| 2026-07-14 | 79cfaaf2 | ensembl116_grch38_giab_sites_hash40                              | rich        |       1 |                1 |                5000 | 100,957    | 644,292     | 5,078,384 | 0                   | 1,179,465      |      9 |       0.225 |          0.228 |       0.232 |              442794 |         2258.4 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 5,173,092              |
| 2026-07-14 | 81b9e2a7 | ensembl116_grch38_primary_giab_sites_hash40                      | compact     |       1 |                1 |                5000 | 100,957    | 646,577     | 5,087,789 | 0                   | 1,179,465      |      9 |       0.112 |          0.115 |       0.125 |              877887 |         1139.1 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 10,256,217             |
| 2026-07-14 | 81b9e2a7 | ensembl116_grch38_primary_giab_sites_hash40                      | rich        |       1 |                1 |                5000 | 100,957    | 646,577     | 5,087,789 | 0                   | 1,179,465      |      9 |       0.233 |          0.236 |       0.242 |              427784 |         2337.6 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 4,997,733              |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_coding_mixed_200k                        | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,757,180      |      3 |      58.833 |         59.010 |      59.212 |                3389 |       295050.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 97,563                 |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_coding_nonsnv_36k                        | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,616      |      5 |      11.837 |         11.898 |      11.959 |                3026 |       330500.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 88,050                 |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_coding_nonsnv_36k                        | rich        |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,616      |      3 |      11.903 |         11.952 |      11.978 |                3012 |       332000.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 87,652                 |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_coding_snv_200k                          | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |      9 |       0.690 |          0.693 |       0.719 |              288600 |         3465.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 7,490,620              |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_coding_snv_200k                          | rich        |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |      5 |       1.213 |          1.238 |       1.262 |              161551 |         6190.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 4,193,053              |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_giab_sites_hash40                        | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 0                   | 1,174,245      |      9 |       0.116 |          0.117 |       0.119 |              862880 |         1158.9 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 10,036,282             |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_giab_sites_hash40                        | rich        |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 0                   | 1,174,245      |      9 |       0.239 |          0.244 |       0.263 |              413758 |         2416.9 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 4,812,480              |
| 2026-07-14 | bcca6f6c | ensembl116_grch38_giab_sites_hash40                              | compact     |       1 |                1 |                5000 | 100,957    | 644,292     | 5,078,384 | 0                   | 1,179,465      |      9 |       0.106 |          0.107 |       0.109 |              943523 |         1059.9 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 11,023,037             |
| 2026-07-14 | bcca6f6c | ensembl116_grch38_giab_sites_hash40                              | rich        |       1 |                1 |                5000 | 100,957    | 644,292     | 5,078,384 | 0                   | 1,179,465      |      9 |       0.221 |          0.225 |       0.231 |              448698 |         2228.7 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 5,242,067              |
| 2026-07-15 | 1903f7fb | ensembl116_grch38_final_coding_mixed_200k                        | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |     21 |       1.489 |          1.514 |       1.618 |              132100 |         7570.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,802,325              |
| 2026-07-15 | 1903f7fb | ensembl116_grch38_final_coding_nonsnv_36k                        | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,524      |     31 |       0.323 |          0.329 |       0.362 |              109422 |         9138.9 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,183,964              |
| 2026-07-15 | 1903f7fb | ensembl116_grch38_final_coding_snv_200k                          | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |     31 |       0.610 |          0.622 |       0.657 |              321543 |         3110.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 8,345,659              |
| 2026-07-15 | 57b5b566 | ensembl116_grch38_final_coding_mixed_200k                        | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |     15 |       1.467 |          1.480 |       1.489 |              135135 |         7400.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,889,676              |
| 2026-07-15 | 57b5b566 | ensembl116_grch38_final_coding_nonsnv_36k                        | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,524      |     21 |       0.310 |          0.322 |       0.351 |              111801 |         8944.4 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,253,180              |
| 2026-07-15 | 57b5b566 | ensembl116_grch38_final_coding_snv_200k                          | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |     15 |       0.609 |          0.615 |       0.638 |              325203 |         3075.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 8,440,650              |
| 2026-07-15 | 7dd90ce8 | ensembl116_grch38_final_coding_mixed_200k                        | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |      9 |       1.550 |          1.572 |       1.583 |              127226 |         7860.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,662,036              |
| 2026-07-15 | 7dd90ce8 | ensembl116_grch38_final_coding_nonsnv_36k                        | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,524      |     11 |       0.341 |          0.344 |       0.357 |              104651 |         9555.6 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,045,128              |
| 2026-07-15 | 7dd90ce8 | ensembl116_grch38_final_coding_snv_200k                          | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |      9 |       0.609 |          0.613 |       0.629 |              326264 |         3065.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 8,468,189              |
| 2026-07-15 | 7dd90ce8 | ensembl116_grch38_final_giab_sites_hash40                        | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 0                   | 1,174,245      |      9 |       0.118 |          0.119 |       0.137 |              848378 |         1178.7 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 9,867,605              |
| 2026-07-15 | a346a7aa | ensembl116_grch38_final_coding_mixed_200k                        | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |     21 |       1.417 |          1.443 |       1.511 |              138600 |         7215.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,989,411              |
| 2026-07-15 | a346a7aa | ensembl116_grch38_final_coding_nonsnv_36k                        | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,524      |     31 |       0.307 |          0.315 |       0.376 |              114286 |         8750.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,325,473              |
| 2026-07-15 | a346a7aa | ensembl116_grch38_final_coding_snv_200k                          | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |     31 |       0.599 |          0.617 |       0.791 |              324149 |         3085.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 8,413,290              |
| 2026-07-15 | c361346f | ensembl116_grch38_final_coding_mixed_200k                        | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |     21 |       1.445 |          1.468 |       1.498 |              136240 |         7340.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,921,471              |
| 2026-07-15 | c361346f | ensembl116_grch38_final_coding_nonsnv_36k                        | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,524      |     31 |       0.320 |          0.323 |       0.341 |              111455 |         8972.2 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,243,108              |
| 2026-07-15 | c361346f | ensembl116_grch38_final_coding_snv_200k                          | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |     31 |       0.584 |          0.602 |       0.779 |              332226 |         3010.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 8,622,924              |
| 2026-07-15 | ecb7eaae | ensembl116_grch38_final_coding_mixed_200k                        | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |     15 |       1.467 |          1.495 |       1.512 |              133779 |         7475.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,850,649              |
| 2026-07-15 | ecb7eaae | ensembl116_grch38_final_coding_nonsnv_36k                        | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,524      |     21 |       0.311 |          0.324 |       0.356 |              111111 |         9000.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,233,099              |
| 2026-07-15 | ecb7eaae | ensembl116_grch38_final_coding_snv_200k                          | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |     15 |       0.605 |          0.610 |       0.624 |              327869 |         3050.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 8,509,836              |
| 2026-07-15 | f649d724 | ensembl116_grch38_final_coding_mixed_200k                        | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |      9 |       1.487 |          1.504 |       1.626 |              132979 |         7520.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,827,606              |
| 2026-07-15 | f649d724 | ensembl116_grch38_final_coding_nonsnv_36k                        | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,524      |     11 |       0.324 |          0.327 |       0.342 |              110092 |         9083.3 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,203,437              |
| 2026-07-15 | f649d724 | ensembl116_grch38_final_coding_snv_200k                          | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |     11 |       0.611 |          0.613 |       0.634 |              326264 |         3065.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 8,468,189              |
| 2026-07-15 | f649d724 | ensembl116_grch38_final_giab_sites_hash40                        | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 0                   | 1,174,245      |     15 |       0.118 |          0.120 |       0.132 |              841308 |         1188.6 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 9,785,375              |
| 2026-07-16 | 220c6ddc | ensembl116_grch38_final_giab_sites_hash40                        | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 0                   | 1,174,245      |     31 |       0.118 |          0.120 |       0.144 |              841308 |         1188.6 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 9,785,375              |
| 2026-07-16 | 220c6ddc | ensembl116_grch38_final_giab_sites_hash40_regulation             | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 1,383,580           | 1,179,329      |     31 |       0.125 |          0.126 |       0.162 |              801246 |         1248.1 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 9,359,754              |
| 2026-07-16 | 360619ed | ensembl116_grch38_final_breakend_200k                            | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 18,766,240     |     15 |       2.676 |          2.710 |       2.727 |               73801 |        13550.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 6,924,812              |
| 2026-07-16 | 96b4cd45 | ensembl116_grch38_final_giab_sites_hash40                        | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 0                   | 1,174,245      |     31 |       0.131 |          0.134 |       0.147 |              753410 |         1327.3 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 8,763,022              |
| 2026-07-16 | 96b4cd45 | ensembl116_grch38_final_giab_sites_hash40_regulation             | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 1,383,580           | 1,179,329      |     31 |       0.138 |          0.140 |       0.162 |              721121 |         1386.7 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 8,423,779              |
| 2026-07-16 | 9c97cb07 | ensembl116_grch38_final_coding_mixed_200k                        | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |     15 |       1.451 |          1.481 |       1.638 |              135044 |         7405.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,887,049              |
| 2026-07-16 | 9c97cb07 | ensembl116_grch38_final_coding_nonsnv_36k                        | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,524      |     15 |       0.317 |          0.322 |       0.339 |              111801 |         8944.4 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 3,253,180              |
| 2026-07-16 | 9c97cb07 | ensembl116_grch38_final_coding_snv_200k                          | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |     15 |       0.593 |          0.603 |       0.709 |              331675 |         3015.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 8,608,624              |
| 2026-07-16 | 9c97cb07 | ensembl116_grch38_final_giab_sites_hash40                        | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 0                   | 1,174,245      |     15 |       0.117 |          0.118 |       0.131 |              855568 |         1168.8 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 9,951,229              |
| 2026-07-16 | f3796ae9 | ensembl116_grch38_final_giab_sites_hash40                        | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 0                   | 1,174,245      |     31 |       0.130 |          0.132 |       0.141 |              764826 |         1307.5 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 8,895,795              |
| 2026-07-16 | f3796ae9 | ensembl116_grch38_final_giab_sites_hash40_regulation             | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 1,383,580           | 1,179,329      |     31 |       0.138 |          0.139 |       0.148 |              726309 |         1376.8 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 8,484,381              |
| 2026-07-19 | 2548d675 | ensembl116_grch38_final_breakend_200k                            | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 18,766,240     |     15 |       2.676 |          2.681 |       2.693 |               74599 |        13405.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 6,999,717              |
| 2026-07-19 | 2548d675 | ensembl116_grch38_final_breakend_regulation_200k                 | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 1,383,580           | 18,877,120     |     15 |       2.737 |          2.761 |       2.769 |               72438 |        13805.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 6,837,059              |
| 2026-07-19 | e7c3623d | ensembl116_grch38_final_breakend_200k                            | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 18,766,240     |     15 |       2.589 |          2.642 |       2.654 |               75700 |        13210.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 7,103,043              |
| 2026-07-19 | e7c3623d | ensembl116_grch38_final_breakend_regulation_200k                 | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 1,383,580           | 18,877,120     |     15 |       2.712 |          2.725 |       2.745 |               73394 |        13625.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 6,927,383              |
| 2026-07-19 | f97101e1 | ensembl116_grch38_final_breakend_200k                            | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 18,766,240     |     15 |       3.056 |          3.105 |       4.435 |               64412 |        15525.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 6,043,878              |
| 2026-07-19 | f97101e1 | ensembl116_grch38_final_breakend_regulation_200k                 | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 1,383,580           | 18,877,120     |     15 |       3.307 |          3.348 |       3.370 |               59737 |        16740.0 | 13th Gen Intel(R) Core(TM) i5-13500 |              | 5,638,327              |
| 2026-07-20 | 0714235a | ensembl116_grch38_final_coding_mixed_200k                        | hgvs        |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |      9 |      32.023 |         32.110 |      32.150 |                6229 |       160550.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 179,281                |
| 2026-07-20 | 0714235a | ensembl116_grch38_final_coding_mixed_200k                        | rich        |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |      9 |       2.031 |          2.038 |       2.047 |               98135 |        10190.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,824,691              |
| 2026-07-20 | 25328813 | ensembl116_grch38_final_coding_mixed_200k                        | rich        |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |      9 |       2.024 |          2.033 |       2.046 |               98377 |        10165.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,831,638              |
| 2026-07-20 | 25328813 | ensembl116_grch38_final_coding_mixed_200k                        | hgvs        |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |      9 |      32.047 |         32.106 |      32.204 |                6229 |       160530.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 179,304                |
| 2026-07-20 | 7dae50cd | ensembl116_grch38_final_coding_mixed_200k                        | rich        |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |      9 |       2.045 |          2.067 |       2.081 |               96759 |        10335.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,785,060              |
| 2026-07-20 | 7dae50cd | ensembl116_grch38_final_coding_mixed_200k                        | hgvs        |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |      9 |       4.565 |          4.579 |       4.610 |               43678 |        22895.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 1,257,200              |
| 2026-07-20 | e25c1513 | ensembl116_grch38_clinvar_annotation_dense_v2                    | compact     |       1 |                1 |                   0 | 517,097    | 644,427     | 5,068,416 | 1,383,580           | 18,163,993     |      5 |       1.828 |          1.834 |       1.845 |              281950 |         3546.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 9,904,031              |
| 2026-07-20 | e25c1513 | ensembl116_grch38_clinvar_annotation_dense_v2                    | compact     |       4 |                4 |                   0 | 517,097    | 644,427     | 5,068,416 | 1,383,580           | 18,163,993     |      5 |       0.535 |          0.540 |       0.545 |              957587 |         1044.3 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 33,637,024             |
| 2026-07-20 | e25c1513 | ensembl116_grch38_clinvar_annotation_dense_v2                    | compact     |       1 |                1 |                5000 | 517,097    | 644,427     | 5,068,416 | 1,383,580           | 26,518,787     |      5 |       2.081 |          2.106 |       2.112 |              245535 |         4072.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 12,592,017             |
| 2026-07-20 | e25c1513 | ensembl116_grch38_clinvar_annotation_dense_v2                    | compact     |       4 |                4 |                5000 | 517,097    | 644,427     | 5,068,416 | 1,383,580           | 26,518,787     |      5 |       0.626 |          0.632 |       0.662 |              818191 |         1222.2 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 41,960,106             |
| 2026-07-20 | e25c1513 | ensembl116_grch38_clinvar_annotation_dense_v2                    | compact     |       1 |                1 |               10000 | 517,097    | 644,427     | 5,068,416 | 1,383,580           | 34,248,323     |      5 |       2.359 |          2.370 |       2.418 |              218184 |         4583.3 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 14,450,769             |
| 2026-07-20 | e25c1513 | ensembl116_grch38_clinvar_annotation_dense_v2                    | compact     |       4 |                4 |               10000 | 517,097    | 644,427     | 5,068,416 | 1,383,580           | 34,248,323     |      5 |       0.713 |          0.719 |       0.758 |              719189 |         1390.5 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 47,633,273             |
| 2026-07-20 | e25c1513 | ensembl116_grch38_clinvar_annotation_dense_v2                    | compact     |       1 |                1 |               50000 | 517,097    | 644,427     | 5,068,416 | 1,383,580           | 88,784,213     |      5 |       4.130 |          4.161 |       4.185 |              124272 |         8046.8 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 21,337,230             |
| 2026-07-20 | e25c1513 | ensembl116_grch38_clinvar_annotation_dense_v2                    | compact     |       4 |                4 |               50000 | 517,097    | 644,427     | 5,068,416 | 1,383,580           | 88,784,213     |      5 |       1.326 |          1.342 |       1.436 |              385318 |         2595.3 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 66,158,132             |
| 2026-07-21 | 03803bd8 | ensembl116_grch38_clinvar_20260706_full_literal_regulation       | compact     |       1 |                1 |                5000 | 4,438,467  | 644,427     | 5,068,416 | 1,383,580           | 126,930,027    |      5 |      11.412 |         11.441 |      11.470 |              387944 |         2577.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 11,094,312             |
| 2026-07-21 | 03803bd8 | ensembl116_grch38_clinvar_20260706_full_literal_regulation       | rich        |       1 |                1 |                5000 | 4,438,467  | 644,427     | 5,068,416 | 1,383,580           | 126,930,027    |      5 |      25.451 |         25.591 |      25.659 |              173439 |         5765.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 4,959,948              |
| 2026-07-21 | 03803bd8 | ensembl116_grch38_clinvar_20260706_full_literal_regulation       | hgvs        |       1 |                1 |                5000 | 4,438,467  | 644,427     | 5,068,416 | 1,383,580           | 126,930,027    |      5 |      47.729 |         47.825 |      47.909 |               92806 |        10775.1 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,654,052              |
| 2026-07-21 | 03803bd8 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_regulation      | compact     |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       3.651 |          3.653 |       3.676 |             1121164 |          891.9 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 13,094,950             |
| 2026-07-21 | 03803bd8 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_regulation      | rich        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       8.718 |          8.747 |       8.794 |              468230 |         2135.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 5,468,829              |
| 2026-07-21 | 03803bd8 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_regulation      | hgvs        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |      16.384 |         16.448 |      16.527 |              249004 |         4016.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,908,308              |
| 2026-07-21 | 056e2d56 | ensembl116_grch38_clinvar_20260706_full_literal_regulation       | compact     |       1 |                1 |                5000 | 4,438,467  | 644,427     | 5,068,416 | 1,383,580           | 126,930,027    |      5 |      11.562 |         11.617 |      11.692 |              382067 |         2617.3 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 10,926,231             |
| 2026-07-21 | 056e2d56 | ensembl116_grch38_clinvar_20260706_full_literal_regulation       | rich        |       1 |                1 |                5000 | 4,438,467  | 644,427     | 5,068,416 | 1,383,580           | 126,930,027    |      5 |      25.804 |         25.899 |      25.988 |              171376 |         5835.1 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 4,900,962              |
| 2026-07-21 | 056e2d56 | ensembl116_grch38_clinvar_20260706_full_literal_regulation       | hgvs        |       1 |                1 |                5000 | 4,438,467  | 644,427     | 5,068,416 | 1,383,580           | 126,930,027    |      5 |      48.761 |         48.869 |      48.906 |               90824 |        11010.3 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,597,353              |
| 2026-07-21 | 056e2d56 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_regulation      | compact     |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       3.639 |          3.668 |       3.679 |             1116579 |          895.6 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 13,041,399             |
| 2026-07-21 | 056e2d56 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_regulation      | rich        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       8.776 |          8.788 |       8.836 |              466046 |         2145.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 5,443,315              |
| 2026-07-21 | 056e2d56 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_regulation      | hgvs        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |      16.765 |         16.790 |      16.872 |              243932 |         4099.5 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,849,068              |
| 2026-07-21 | 310890b5 | ensembl116_grch38_clinvar_20260706_full_literal_regulation       | compact     |       1 |                1 |                5000 | 4,438,467  | 644,427     | 5,068,416 | 1,383,580           | 126,930,027    |      5 |      11.550 |         11.667 |      11.687 |              380429 |         2628.6 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 10,879,406             |
| 2026-07-21 | 310890b5 | ensembl116_grch38_clinvar_20260706_full_literal_regulation       | rich        |       1 |                1 |                5000 | 4,438,467  | 644,427     | 5,068,416 | 1,383,580           | 126,930,027    |      5 |      25.854 |         25.939 |      25.982 |              171112 |         5844.1 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 4,893,405              |
| 2026-07-21 | 310890b5 | ensembl116_grch38_clinvar_20260706_full_literal_regulation       | hgvs        |       1 |                1 |                5000 | 4,438,467  | 644,427     | 5,068,416 | 1,383,580           | 126,930,027    |      5 |      48.586 |         48.676 |      48.927 |               91184 |        10966.8 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,607,651              |
| 2026-07-21 | 310890b5 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_regulation      | compact     |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       3.686 |          3.703 |       3.721 |             1106025 |          904.1 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 12,918,134             |
| 2026-07-21 | 310890b5 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_regulation      | rich        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       8.854 |          8.885 |       8.944 |              460958 |         2169.4 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 5,383,889              |
| 2026-07-21 | 310890b5 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_regulation      | hgvs        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |      16.693 |         16.742 |      16.791 |              244631 |         4087.8 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,857,236              |
| 2026-07-21 | 33a15026 | ensembl116_grch38_clinvar_20260706_full_literal                  | compact     |       1 |                1 |                5000 | 4,438,467  | 644,427     | 5,068,416 | 0                   | 126,733,707    |      5 |      11.317 |         11.327 |      11.370 |              391848 |         2552.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 11,188,638             |
| 2026-07-21 | 33a15026 | ensembl116_grch38_clinvar_20260706_full_literal                  | rich        |       1 |                1 |                5000 | 4,438,467  | 644,427     | 5,068,416 | 0                   | 126,733,707    |      5 |      25.287 |         25.365 |      25.415 |              174984 |         5714.8 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 4,996,401              |
| 2026-07-21 | 33a15026 | ensembl116_grch38_clinvar_20260706_full_literal                  | hgvs        |       1 |                1 |                5000 | 4,438,467  | 644,427     | 5,068,416 | 0                   | 126,733,707    |      5 |      48.625 |         48.993 |      51.699 |               90594 |        11038.3 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,586,772              |
| 2026-07-21 | 33a15026 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal                 | compact     |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 0                   | 47,629,345     |      5 |       3.592 |          4.446 |       6.477 |              921190 |         1085.6 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 10,712,853             |
| 2026-07-21 | 33a15026 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal                 | rich        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 0                   | 47,629,345     |      5 |       9.147 |          9.157 |       9.244 |              447266 |         2235.8 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 5,201,414              |
| 2026-07-21 | 33a15026 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal                 | hgvs        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 0                   | 47,629,345     |      5 |      16.701 |         16.770 |      17.858 |              244222 |         4094.6 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,840,152              |
| 2026-07-21 | 6c640890 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | compact     |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       4.162 |          4.247 |       4.427 |              964354 |         1037.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 11,263,445             |
| 2026-07-21 | 6c640890 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | compact     |       4 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       1.271 |          1.283 |       1.334 |             3192214 |          313.3 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 37,284,373             |
| 2026-07-21 | 6c640890 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | rich        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       9.566 |         11.836 |      35.602 |              346030 |         2889.9 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 4,041,556              |
| 2026-07-21 | 6c640890 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | rich        |       4 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       2.693 |          2.706 |       2.759 |             1513530 |          660.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 17,677,698             |
| 2026-07-21 | 6c640890 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | hgvs        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |      17.117 |         17.185 |      17.346 |              238325 |         4196.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,783,582              |
| 2026-07-21 | 6c640890 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | hgvs        |       4 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       4.725 |          4.770 |       4.884 |              858619 |         1164.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 10,028,480             |
| 2026-07-21 | 6c640890 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | rich_hgvs   |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |      22.609 |         22.691 |      22.890 |              180495 |         5540.3 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,108,142              |
| 2026-07-21 | 6c640890 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | rich_hgvs   |       4 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       6.228 |          6.421 |       6.447 |              637846 |         1567.8 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 7,449,907              |
| 2026-07-21 | 8da057b0 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal                 | compact     |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 0                   | 47,629,345     |      5 |       3.556 |          3.573 |       3.585 |             1146267 |          872.4 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 13,330,351             |
| 2026-07-21 | 8da057b0 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal                 | hgvs        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 0                   | 47,629,345     |      5 |      16.486 |         16.505 |      16.565 |              248144 |         4029.9 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,885,752              |
| 2026-07-21 | 8da057b0 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | compact     |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 0                   | 47,629,345     |      5 |       4.072 |          4.092 |       4.268 |             1000882 |          999.1 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 11,639,625             |
| 2026-07-21 | 8da057b0 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | compact     |       4 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 0                   | 47,629,345     |      5 |       1.206 |          1.305 |       1.339 |             3138399 |          318.6 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 36,497,582             |
| 2026-07-21 | 8da057b0 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | hgvs        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 0                   | 47,629,345     |      5 |      17.070 |         17.151 |      17.374 |              238797 |         4187.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,777,059              |
| 2026-07-21 | 8da057b0 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | hgvs        |       4 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 0                   | 47,629,345     |      5 |       5.491 |          5.553 |       5.823 |              737549 |         1355.8 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 8,577,228              |
| 2026-07-21 | c6b532b4 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal                 | compact     |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 0                   | 47,629,345     |      5 |       3.513 |          3.549 |       4.983 |             1154018 |          866.5 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 13,420,497             |
| 2026-07-21 | c6b532b4 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal                 | hgvs        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 0                   | 47,629,345     |      5 |      16.247 |         16.291 |      16.320 |              251403 |         3977.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,923,660              |
| 2026-07-21 | c6b532b4 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | compact     |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 0                   | 47,629,345     |      5 |       4.014 |          4.046 |       4.244 |             1012262 |          987.9 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 11,771,959             |
| 2026-07-21 | c6b532b4 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | compact     |       4 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 0                   | 47,629,345     |      5 |       1.202 |          1.215 |       1.265 |             3370873 |          296.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 39,201,107             |
| 2026-07-21 | c6b532b4 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | hgvs        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 0                   | 47,629,345     |      5 |      16.868 |         16.978 |      17.336 |              241230 |         4145.4 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,805,357              |
| 2026-07-21 | c6b532b4 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | hgvs        |       4 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 0                   | 47,629,345     |      5 |       4.622 |          4.777 |       5.277 |              857360 |         1166.4 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 9,970,556              |
| 2026-07-22 | ca35fd7b | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | compact     |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       4.245 |          4.260 |       4.430 |              961411 |         1040.1 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 11,229,073             |
| 2026-07-22 | ca35fd7b | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | compact     |       4 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       1.273 |          1.291 |       1.329 |             3172433 |          315.2 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 37,053,332             |
| 2026-07-22 | ca35fd7b | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | rich        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       9.510 |          9.520 |       9.717 |              430211 |         2324.4 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 5,024,774              |
| 2026-07-22 | ca35fd7b | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | rich        |       4 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       2.669 |          2.681 |       2.736 |             1527643 |          654.6 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 17,842,540             |
| 2026-07-22 | ca35fd7b | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | hgvs        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |      17.198 |         17.326 |      17.498 |              236385 |         4230.4 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,760,929              |
| 2026-07-22 | ca35fd7b | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | hgvs        |       4 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       4.736 |          4.883 |       4.934 |              838749 |         1192.3 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 9,796,406              |
| 2026-07-22 | ca35fd7b | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | rich_hgvs   |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |      22.753 |         22.771 |      23.063 |              179861 |         5559.9 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,100,736              |
| 2026-07-22 | ca35fd7b | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | rich_hgvs   |       4 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       6.251 |          6.487 |       6.537 |              631357 |         1583.9 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 7,374,110              |
| 2026-07-23 | 0eb1b440 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | compact     |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       4.264 |          4.287 |       4.501 |              955356 |         1046.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 11,158,351             |
| 2026-07-23 | 0eb1b440 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | rich        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       9.665 |          9.682 |       9.864 |              423013 |         2364.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 4,940,699              |
| 2026-07-23 | 0eb1b440 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | hgvs        |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |      17.992 |         18.005 |      18.152 |              227471 |         4396.2 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,656,809              |
| 2026-07-23 | 0eb1b440 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | rich_hgvs   |       1 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |      23.516 |         23.640 |      23.789 |              173249 |         5772.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,023,513              |
| 2026-07-23 | 0eb1b440 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | compact     |       4 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       1.295 |          1.300 |       1.357 |             3150470 |          317.4 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 36,796,808             |
| 2026-07-23 | 0eb1b440 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | rich        |       4 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       2.715 |          2.740 |       2.818 |             1494749 |          669.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 17,458,340             |
| 2026-07-23 | 0eb1b440 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | hgvs        |       4 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       4.916 |          5.034 |       5.094 |              813590 |         1229.1 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 9,502,553              |
| 2026-07-23 | 0eb1b440 | ensembl116_grch38_giab_hg002_v4_2_1_full_literal_public_relation | rich_hgvs   |       4 |                1 |                5000 | 4,095,611  | 644,427     | 5,068,416 | 1,383,580           | 47,835,851     |      5 |       6.410 |          6.602 |       6.608 |              620359 |         1612.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 7,245,661              |
| 2026-09-06 | 41e84614 | fixture_one_transcript_sorted                                    | rich        |       1 |                1 |                5000 | 10,000,000 | 1           | 2         | 0                   | 10,000,000     |      3 |       3.320 |          3.321 |       3.346 |             3011141 |          332.1 | 13th Gen Intel(R) Core(TM) i5-13500 | 0-19         | 3,011,141              |
| 2026-09-06 | fc67aef3 | fixture_one_transcript_sorted                                    | rich        |       1 |                1 |                5000 | 10,000,000 | 1           | 2         | 0                   | 10,000,000     |      3 |       3.316 |          3.333 |       3.339 |             3000300 |          333.3 | 13th Gen Intel(R) Core(TM) i5-13500 | 0-19         | 3,000,300              |
| 2026-09-06 | e5f1f269 | fixture_one_transcript_sorted_indels                             | hgvs        |       1 |                1 |                5000 | 1,000,000  | 1           | 2         | 0                   | 1,000,000      |      5 |       1.239 |          1.242 |       1.254 |              805153 |         1242.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 805,153                |
| 2026-09-06 | 10320db6 | fixture_one_transcript_sorted_indels                             | hgvs        |       1 |                1 |                5000 | 1,000,000  | 1           | 2         | 0                   | 1,000,000      |      5 |       1.233 |          1.240 |       1.249 |              806452 |         1240.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 806,452                |
| 2026-09-06 | 60a15f21 | fixture_one_transcript_sorted_breakend                           | compact     |       1 |                1 |                5000 | 1,000,000  | 1           | 2         | 0                   | 1,000,000      |      5 |       0.266 |          0.267 |       0.272 |             3745318 |          267.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 3,745,318              |
| 2026-09-06 | 9ed772f7 | fixture_one_transcript_sorted_breakend                           | compact     |       1 |                1 |                5000 | 1,000,000  | 1           | 2         | 0                   | 1,000,000      |      5 |       0.266 |          0.270 |       0.272 |             3703704 |          270.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 3,703,704              |
| 2026-09-06 | 9ed772f7 | fixture_one_transcript_sorted_breakend_raw_alt                   | compact     |       1 |                1 |                5000 | 1,000,000  | 1           | 2         | 0                   | 1,000,000      |      5 |       0.394 |          0.396 |       0.402 |             2525253 |          396.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,525,253              |

Each pass consumes every staged input and checks output cardinality plus
either the rendered consequence-byte total or the numeric
consequence-mask sum. Rows with different workload, output mode, thread
count, host, or variant count are separate measurements and should not
be compared as if only the engine changed. `cpu_affinity` is read from
the benchmark process itself; a single integer records an externally
pinned run, while a range records a scheduler-visible CPU set.

## Borrowed-CDS projection bounds check

The `--fixture-indels` workload alternates `T>TAC` insertions and `TA>T`
deletions at position 124 of the checked-in one-transcript model and
FASTA. Unlike the single-SNV adapter control, these length-changing
alleles exercise CDS edit projection and HGVS construction. They remain
independent events, not one compound haplotype. Model loading and input
generation are untimed.

| revision | variants | annotated_rows | threads | passes | min_seconds | median_seconds | max_seconds | checksum_value |
|:---------|---------:|---------------:|--------:|-------:|------------:|---------------:|------------:|:---------------|
| e5f1f269 |    1e+06 |          1e+06 |       1 |      5 |       1.239 |          1.242 |       1.254 | 26000000       |
| 10320db6 |    1e+06 |          1e+06 |       1 |      5 |       1.233 |          1.240 |       1.249 | 26000000       |

Both clean in-tree builds use DuckDB v1.5.3 on 13th Gen Intel(R)
Core(TM) i5-13500, pinned to CPU 2, with 100,000 warm-up events and a
5,000-base transcript halo. The baseline commit adds only the benchmark
case; the second commit adds the borrowed-slice checks and genomic
carrier differential. Their median difference is -0.16%. This narrow
valid-model workload does not measure malformed-model rejection, file
ingestion, cohort state, sort-plus-phased execution, or peak memory.

## Annotation-dense transcript distance and ordered parallel partitions

| revision | transcript_distance | threads | ordered_partitions | input_alleles | output_rows | median_seconds | input_alleles_per_second | output_rows_per_second | speedup | peak_RSS_GiB | RSS_above_one_thread_MiB |
|:---------|:--------------------|--------:|-------------------:|:--------------|:------------|---------------:|:-------------------------|:-----------------------|--------:|-------------:|-------------------------:|
| e25c1513 | 0                   |       1 |                  1 | 517,097       | 18,163,993  |          1.834 | 281,950                  | 9,904,031              |    1.00 |        5.194 |                     0.00 |
| e25c1513 | 0                   |       4 |                  4 | 517,097       | 18,163,993  |          0.540 | 957,587                  | 33,637,024             |    3.40 |        5.200 |                     6.54 |
| e25c1513 | 5,000               |       1 |                  1 | 517,097       | 26,518,787  |          2.106 | 245,535                  | 12,592,017             |    1.00 |        5.194 |                     0.00 |
| e25c1513 | 5,000               |       4 |                  4 | 517,097       | 26,518,787  |          0.632 | 818,191                  | 41,960,106             |    3.33 |        5.201 |                     7.03 |
| e25c1513 | 10,000              |       1 |                  1 | 517,097       | 34,248,323  |          2.370 | 218,184                  | 14,450,769             |    1.00 |        5.193 |                     0.00 |
| e25c1513 | 10,000              |       4 |                  4 | 517,097       | 34,248,323  |          0.719 | 719,189                  | 47,633,273             |    3.30 |        5.203 |                     9.80 |
| e25c1513 | 50,000              |       1 |                  1 | 517,097       | 88,784,213  |          4.161 | 124,272                  | 21,337,230             |    1.00 |        5.193 |                     0.00 |
| e25c1513 | 50,000              |       4 |                  4 | 517,097       | 88,784,213  |          1.342 | 385,318                  | 66,158,132             |    3.10 |        5.201 |                     8.16 |

This is the release-116 GRCh38 model with all 644,427 admitted
transcripts, 5,068,416 exon memberships, and 1,383,580 regulation/motif
intervals. The 517,097 ClinVar alleles were selected from 318
annotation-dense genomic tiles, not diluted with transcript-sparse
whole-genome positions. The largest REF and ALT alleles are 9,974 and
9,971 bases. Each row is the median of five complete passes after a
100,000-allele warmup. Model loading and corpus staging are outside the
timed pass; compact list production, `unnest`, and final aggregation are
included.

One-thread runs are pinned to physical P core 2. Four-thread runs are
pinned to physical P cores 2, 4, 6, and 8 and explicitly split the
already ordered input into four disjoint coordinate ranges. Merely
setting DuckDB `threads = 4` does not split one scalar-function branch.
Each concurrent scalar callback borrows an exclusive mutable workspace
from the shared immutable model; a partition is not permanently bound to
one workspace. `input_variant_index` is the global input ordinal, while
`annotation_index` is the one-based position inside each variant’s
result list. Consumers that require a canonical materialized row order
use `ORDER BY input_variant_index, annotation_index` after the parallel
union.

The one- and four-partition runs have identical XOR and sum fingerprints
over every `(input_variant_index, annotation_index, public row)` at each
tested distance. These commutative fingerprints prove row-multiset and
per-input-list equality, not global emission order. The 50,000-base case
emits 88,784,213 rows and still scales from 124,272 to 385,318 input
alleles/s. Peak RSS is GNU `/usr/bin/time -v` maximum resident set size
for the complete one-process R/DuckDB invocation, including model load,
warmup, and all five passes. Across the matrix, the largest observed
four-thread premium over its one-thread pair is 9.80 MiB. Source
ownership establishes that callbacks share the immutable model; the RSS
result shows no model-sized process increase for these four-branch runs
but is not allocation attribution. Wider transcript distance increases
candidate pairs and output cardinality, not model ownership. These
observations apply to the tested 0, 5,000, 10,000, and 50,000-base
configurations; arbitrary partition counts still require measurement.

The benchmark source accepts every transcript distance from zero through
`UINT32_MAX`; no classifier or sweep buffer is sized for a 5,000-base
special case. VEP’s separate fixed 5,000-base paired-breakend
close-to-feature rule is a biological compatibility predicate applied
after transcript discovery and must not be confused with the
configurable transcript-distance search.

## Routine optimization backlog

The declared fixed-event consequence contract no longer depends on
another semantic rewrite. Future speed work should preserve the
conformance and full-row-fingerprint gates while attacking measured
costs in this order:

1.  write compact SoA fields directly into DuckDB result vectors,
    avoiding the compact-list-to-`unnest` materialization round trip;
2.  make the common coding-SNV route leaner while retaining one
    generated effect composition authority;
3.  reuse genomic-to-transcript, exon, splice-window, and coding
    projection facts that are currently rediscovered by adjacent
    predicates;
4.  separate HGVS fact production from string rendering and reuse one
    normalized shift result across HGVSc/HGVSn/HGVSp;
5.  expose a persistent sorted table-function sweep so one logical input
    stream does not reseed at every DuckDB vector edge; and
6.  make disjoint ordered partitions a first-class planner surface,
    keeping stable identifiers numeric and joining presentation strings
    late.

Exploratory local profiling suggested low branch-miss and cache-miss
rates in the candidate sweep and a small candidate-discovery share of
dense runtime. That profile is not retained as a checked comparable
benchmark, so it is an avenue-selection hint, not a quantitative claim.
Broad SIMD work is therefore not the first optimization target. Each
item above needs an adjacent exact-revision comparison on the same
model, corpus, transcript distance, output denominator, core affinity,
and fingerprint before it can be called an improvement.

## Public relation API cost and scaling

| revision | output    | threads | input_alleles | output_rows | median_seconds | alleles_per_second | output_rows_per_second | scaling_vs_one_core | elapsed_vs_compact |
|:---------|:----------|--------:|:--------------|:------------|---------------:|:-------------------|:-----------------------|--------------------:|-------------------:|
| 0eb1b440 | compact   |       1 | 4,095,611     | 47,835,851  |          4.287 | 955,356            | 11,158,351             |                1.00 |               1.00 |
| 0eb1b440 | compact   |       4 | 4,095,611     | 47,835,851  |          1.300 | 3,150,470          | 36,796,808             |                3.30 |               1.00 |
| 0eb1b440 | rich      |       1 | 4,095,611     | 47,835,851  |          9.682 | 423,013            | 4,940,699              |                1.00 |               2.26 |
| 0eb1b440 | rich      |       4 | 4,095,611     | 47,835,851  |          2.740 | 1,494,749          | 17,458,340             |                3.53 |               2.11 |
| 0eb1b440 | hgvs      |       1 | 4,095,611     | 47,835,851  |         18.005 | 227,471            | 2,656,809              |                1.00 |               4.20 |
| 0eb1b440 | hgvs      |       4 | 4,095,611     | 47,835,851  |          5.034 | 813,590            | 9,502,553              |                3.58 |               3.87 |
| 0eb1b440 | rich_hgvs |       1 | 4,095,611     | 47,835,851  |         23.640 | 173,249            | 2,023,513              |                1.00 |               5.51 |
| 0eb1b440 | rich_hgvs |       4 | 4,095,611     | 47,835,851  |          6.602 | 620,359            | 7,245,661              |                3.58 |               5.08 |

### Adjacent revision comparison

| output    | threads | previous_revision | current_revision | previous_median_seconds | current_median_seconds | elapsed_change_percent | previous_alleles_per_second | current_alleles_per_second | throughput_change_percent |
|:----------|--------:|:------------------|:-----------------|------------------------:|-----------------------:|-----------------------:|:----------------------------|:---------------------------|--------------------------:|
| compact   |       1 | ca35fd7b          | 0eb1b440         |                   4.260 |                  4.287 |                   0.63 | 961,411                     | 955,356                    |                     -0.63 |
| compact   |       4 | ca35fd7b          | 0eb1b440         |                   1.291 |                  1.300 |                   0.70 | 3,172,433                   | 3,150,470                  |                     -0.69 |
| rich      |       1 | ca35fd7b          | 0eb1b440         |                   9.520 |                  9.682 |                   1.70 | 430,211                     | 423,013                    |                     -1.67 |
| rich      |       4 | ca35fd7b          | 0eb1b440         |                   2.681 |                  2.740 |                   2.20 | 1,527,643                   | 1,494,749                  |                     -2.15 |
| hgvs      |       1 | ca35fd7b          | 0eb1b440         |                  17.326 |                 18.005 |                   3.92 | 236,385                     | 227,471                    |                     -3.77 |
| hgvs      |       4 | ca35fd7b          | 0eb1b440         |                   4.883 |                  5.034 |                   3.09 | 838,749                     | 813,590                    |                     -3.00 |
| rich_hgvs |       1 | ca35fd7b          | 0eb1b440         |                  22.771 |                 23.640 |                   3.82 | 179,861                     | 173,249                    |                     -3.68 |
| rich_hgvs |       4 | ca35fd7b          | 0eb1b440         |                   6.487 |                  6.602 |                   1.77 | 631,357                     | 620,359                    |                     -1.74 |

The adjacent table compares the nearest complete ancestor on the same
host, corpus, model semantics, output denominator, core affinity, and
pass count. Positive elapsed change means slower execution; negative
throughput change means fewer input alleles per second. Identical output
checksums are required before a pair is admitted to this table.

This exact-revision comparison uses the complete Ensembl 116 GRCh38
model: 644,427 transcripts, 5,068,416 exon memberships, and all
1,383,580 resident RegulatoryFeature and MotifFeature intervals. It
annotates every one of the 4,095,611 model-addressable literal alleles
staged from GIAB HG002 v4.2.1 and emits 47,835,851 rows at a 5,000-base
transcript distance after a 100,000-allele warmup and five complete
passes. One-thread runs are pinned to CPU 2; four-thread runs are
restricted to CPUs 2, 4, 6, and 8. Input VCF decoding, coordinate
sorting, canonical event-table staging, and model loading are outside
the timed pass. Each public row times event validation and family
dispatch through `duckvep_annotate(...)`, native consequence/HGVS work,
fixed-schema row expansion, and checksum aggregation. The current
one-thread/four-thread pair requires the same aggregate checksum for
each output contract. The current rows do not retain the optional
untimed full-row fingerprint. `elapsed_vs_compact` compares projections
at the same thread count; it is not a comparison against a private
kernel lane.

## Historical full-corpus transcript-only diagnostic

| revision | corpus             | output  | input_alleles | output_rows | rows_per_allele | median_seconds | alleles_per_second | output_rows_per_second | elapsed_vs_compact | projected_700M_allele_hours |
|:---------|:-------------------|:--------|:--------------|:------------|----------------:|---------------:|:-------------------|:-----------------------|-------------------:|----------------------------:|
| 33a15026 | ClinVar 2026-07-06 | compact | 4,438,467     | 126,733,707 |           28.55 |         11.327 | 391,848            | 11,188,638             |               1.00 |                        0.50 |
| 33a15026 | ClinVar 2026-07-06 | rich    | 4,438,467     | 126,733,707 |           28.55 |         25.365 | 174,984            | 4,996,401              |               2.24 |                        1.11 |
| 33a15026 | ClinVar 2026-07-06 | hgvs    | 4,438,467     | 126,733,707 |           28.55 |         48.993 | 90,594             | 2,586,772              |               4.33 |                        2.15 |
| 33a15026 | GIAB HG002 v4.2.1  | compact | 4,095,611     | 47,629,345  |           11.63 |          4.446 | 921,190            | 10,712,853             |               1.00 |                        0.21 |
| 33a15026 | GIAB HG002 v4.2.1  | rich    | 4,095,611     | 47,629,345  |           11.63 |          9.157 | 447,266            | 5,201,414              |               2.06 |                        0.43 |
| 33a15026 | GIAB HG002 v4.2.1  | hgvs    | 4,095,611     | 47,629,345  |           11.63 |         16.770 | 244,222            | 2,840,152              |               3.77 |                        0.80 |

These are full prepared literal-allele relations, not annotation-density
samples. They are retained historical diagnostics from before core
regulation/motif loading; they are not the production comparison. The
ClinVar 2026-07-06 source contains 4,439,617 records and 4,438,515
source ALTs; 4,438,467 literal, model-addressable alleles enter the
benchmark (4,134,694 SNVs and 303,773 non-SNVs). The longest retained
REF and ALT are 9,983 and 9,971 bases. The GIAB HG002 GRCh38 v4.2.1
benchmark source contains 4,048,342 records and 4,096,123 source ALTs;
4,095,611 literal, model-addressable alleles enter the benchmark
(3,463,000 SNVs and 632,611 non-SNVs). Its longest retained REF and ALT
are 9,704 and 3,686 bases.

Every row uses the complete Ensembl 116 GRCh38 transcript model, a
5,000-base transcript distance, one ordered input partition, one thread
pinned to CPU 2, a 100,000-allele warmup, and five complete timed
passes. Model loading, VCF staging, and the untimed corpus receipt are
excluded. Regulatory and motif features are not loaded in these rows.
The timed query includes candidate discovery, consequence calculation,
list materialization, `unnest`, and the mode-specific checksum
aggregation.

The 126,733,707 ClinVar rows are 28.55 transcript rows per input allele;
the 47,629,345 GIAB rows are 11.63 per allele. Consequently, the same
HGVS engine reports 90,594 input alleles/s on ClinVar and 244,222/s on
GIAB while sustaining 2.59 and 2.84 million expanded rows/s
respectively. At those measured rates, 700 million input alleles project
to 2.15 or 0.80 single-core hours; that range is an input-density
sensitivity, not a gnomAD runtime claim. The full ClinVar HGVS pass also
supplied the 1,406-base-ALT regression that found the former
scratch-buffer over-read; all five recorded HGVS passes complete after
the exact-size retry fix.

## Full-corpus core VEP model

| revision | corpus             | output  | resident_regulation_motif | input_alleles | output_rows | rows_per_allele | median_seconds | alleles_per_second | output_rows_per_second | elapsed_vs_compact | projected_700M_allele_hours |
|:---------|:-------------------|:--------|:--------------------------|:--------------|:------------|----------------:|---------------:|:-------------------|:-----------------------|-------------------:|----------------------------:|
| 03803bd8 | ClinVar 2026-07-06 | compact | 1,383,580                 | 4,438,467     | 126,930,027 |           28.60 |         11.441 | 387,944            | 11,094,312             |               1.00 |                        0.50 |
| 03803bd8 | ClinVar 2026-07-06 | rich    | 1,383,580                 | 4,438,467     | 126,930,027 |           28.60 |         25.591 | 173,439            | 4,959,948              |               2.24 |                        1.12 |
| 03803bd8 | ClinVar 2026-07-06 | hgvs    | 1,383,580                 | 4,438,467     | 126,930,027 |           28.60 |         47.825 | 92,806             | 2,654,052              |               4.18 |                        2.10 |
| 03803bd8 | GIAB HG002 v4.2.1  | compact | 1,383,580                 | 4,095,611     | 47,835,851  |           11.68 |          3.653 | 1,121,164          | 13,094,950             |               1.00 |                        0.17 |
| 03803bd8 | GIAB HG002 v4.2.1  | rich    | 1,383,580                 | 4,095,611     | 47,835,851  |           11.68 |          8.747 | 468,230            | 5,468,829              |               2.39 |                        0.42 |
| 03803bd8 | GIAB HG002 v4.2.1  | hgvs    | 1,383,580                 | 4,095,611     | 47,835,851  |           11.68 |         16.448 | 249,004            | 2,908,308              |               4.50 |                        0.78 |

| revision | corpus             | output  | extension_build               | extension_binary | physical_model | logical_model | corpus_source | staged_corpus | public_row_xor       | public_row_sum               |
|:---------|:-------------------|:--------|:------------------------------|:-----------------|:---------------|:--------------|:--------------|:--------------|:---------------------|:-----------------------------|
| 03803bd8 | ClinVar 2026-07-06 | compact | htslib_distclean_make_release | 53af1444f11b     | 4c2077c83958   | 296bc9063356  | 59a83b34d425  | fde5aa206844  | 5434592916959930391  | 1170715483948138239433590223 |
| 03803bd8 | ClinVar 2026-07-06 | rich    | htslib_distclean_make_release | 53af1444f11b     | 4c2077c83958   | 296bc9063356  | 59a83b34d425  | fde5aa206844  | 15731069621491065835 | 1170781327306548312741272715 |
| 03803bd8 | ClinVar 2026-07-06 | hgvs    | htslib_distclean_make_release | 53af1444f11b     | 4c2077c83958   | 296bc9063356  | 59a83b34d425  | fde5aa206844  | 3854825828222035738  | 1170739873865994974587183880 |
| 03803bd8 | GIAB HG002 v4.2.1  | compact | htslib_distclean_make_release | 53af1444f11b     | 4c2077c83958   | 296bc9063356  | adb4d4a50048  | 9e8eccb0719e  | 18071298266035115485 | 441264289735621502248147381  |
| 03803bd8 | GIAB HG002 v4.2.1  | rich    | htslib_distclean_make_release | 53af1444f11b     | 4c2077c83958   | 296bc9063356  | adb4d4a50048  | 9e8eccb0719e  | 2464071972670265587  | 441226953508179288041384557  |
| 03803bd8 | GIAB HG002 v4.2.1  | hgvs    | htslib_distclean_make_release | 53af1444f11b     | 4c2077c83958   | 296bc9063356  | adb4d4a50048  | 9e8eccb0719e  | 13799808553808525439 | 441236299918702761797429749  |

This is the complete resident core-VEP model used for independent small
variants: 644,427 transcripts plus 380,818 RegulatoryFeature and
1,002,762 MotifFeature intervals. It is a production core-consequence
workload, not an end-to-end product benchmark: VCF/Parquet ingestion and
external exact, positional, gene, interval, frequency, and
clinical-resource joins remain separate DuckDB stages.

The core model emits 126,930,027 rows for ClinVar and 47,835,851 for
GIAB, including 196,320 and 206,506 regulation/motif rows respectively.
The pinned single-core rates are 387,944 compact, 173,439 rich, and
92,806 cumulative-HGVS alleles/s for ClinVar; GIAB reaches 1,121,164,
468,230, and 249,004/s. At the measured HGVS rates, 700 million input
alleles project to 2.10 or 0.78 single-core hours before supplementary
SQL joins.

Every row in the production-core table is tied to exact source revision
03803bd8dff47f0d41584e7fb7a79d42c9183f04, a vendored-htslib `distclean`
followed by an in-tree `make release` rebuild and its extension-binary
digest, the physical and logical model digests, the Ensembl
source-manifest and logical-reference digests, the model-region ordinal
digest, the staged-corpus digest, and the original VCF digest. The HGVS
rows additionally record the physical FASTA and FAI digests. The receipt
table shows shortened artifact digests for readability and the complete
schema-version-2 public-row XOR/sum fingerprints; the full SHA-256
values are retained in `benchmarks/data/duckvep_throughput.csv`.

`_duckvep_annotate_small_hgvs(...)` preserves interval-feature rows from
the same candidate sweep and leaves their transcript/protein HGVS fields
NULL. The HGVS byte checksum is therefore unchanged while output
cardinality includes every regulatory/motif row; SQL and R regressions
pin that pass-through. The earlier transcript-only GIAB run ranged from
3.592 to 6.477 seconds across five passes, whereas this core-model run
ranged from 3.651 to 3.676 seconds. The apparently faster compact result
is not attributed to loading more intervals; use the stable absolute
core-model result rather than treating non-adjacent noisy runs as an
optimization comparison.

## Independent-event HGVS materialization

| revision | workload                                  | output_mode | threads | cpu_affinity | variants | annotated_rows | median_seconds | variants_per_second | checksum_kind          | output_rows_per_second | elapsed_vs_rich |
|:---------|:------------------------------------------|:------------|--------:|:-------------|:---------|:---------------|---------------:|--------------------:|:-----------------------|:-----------------------|----------------:|
| 7dae50cd | ensembl116_grch38_final_coding_mixed_200k | rich        |       1 | 2            | 200,000  | 5,756,720      |          2.067 |               96759 | consequence_text_bytes | 2,785,060              |            1.00 |
| 7dae50cd | ensembl116_grch38_final_coding_mixed_200k | hgvs        |       1 | 2            | 200,000  | 5,756,720      |          4.579 |               43678 | hgvs_text_status_bytes | 1,257,200              |            2.22 |

| baseline_revision | current_revision | input_variants | output_rows | baseline_seconds | current_seconds | speedup | current_input_variants_per_second | current_output_rows_per_second |
|:------------------|:-----------------|:---------------|:------------|-----------------:|----------------:|--------:|----------------------------------:|:-------------------------------|
| 25328813          | 7dae50cd         | 200,000        | 5,756,720   |           32.106 |           4.579 |    7.01 |                             43678 | 1,257,200                      |

[Ferro-hgvs
v0.9.0](https://github.com/fulcrumgenomics/ferro-hgvs/tree/278e2c11134e3b49067d0c334f650c7c29db9cbe)
reports 5.1 million patterns/s for parsing already rendered HGVS and a
77,000 patterns/s full-population normalization peak on an Apple M2 Max.
Those are useful independent reference points but not the operation
timed here: DuckVEP begins with genomic alleles, discovers every
overlapping transcript, computes consequences, shifts and renders
transcript/protein HGVS, and emits 28.78 transcript rows per input in
this workload. The table therefore reports both input alleles/s and
expanded transcript-HGVS rows/s; it does not present a direct
cross-machine speed ranking.

The HGVS row materializes transcript/protein strings and status fields
and checks their total byte count. Compare it only with the rich
consequence row at the same source revision, host, workload, variant
count, thread count, and CPU affinity. Model loading, FASTA validation,
and input staging remain outside the timed pass; lazy per-worker faidx
setup occurs during warmup.

## Core regulation and motif sweep

| lane                           | revision | input_variants | resident_interval_features | output_rows | median_seconds | variants_per_second | output_rows_per_second |
|:-------------------------------|:---------|:---------------|:---------------------------|:------------|---------------:|:--------------------|:-----------------------|
| transcripts                    | 220c6ddc | 100,957        | 0                          | 1,174,245   |          0.120 | 841,308             | 9,785,375              |
| transcripts + regulation/motif | 220c6ddc | 100,957        | 1,383,580                  | 1,179,329   |          0.126 | 801,246             | 9,359,754              |

This is one pinned-core comparison at one source revision. Both lanes
load the same transcript model and consume the same coordinate-sorted
GIAB sample. The second lane additionally loads every VEP-admitted
GRCh38 RegulatoryFeature and MotifFeature interval and emits their
compact consequence rows. The comparison therefore measures candidate
traversal and the extra output together; it is not an isolated
interval-lookup microbenchmark.

## Paired breakends

### Raw ALT preparation and endpoint-union correction

The 2026-09-06 fixture comparison uses source
`60a15f2121c5674fcae7bba6a8ff22d0298988f8` (benchmark-only baseline) and
`9ed772f71376b9ee47fcac03ba3e010fb2bfc973` (native ALT parser and
endpoint-union fix). Each run repeats **1,000,000 input BND events into
1,000,000 compact transcript rows**, using one DuckDB thread pinned to
logical CPU 2 of the same i5-13500, DuckDB 1.5.3, five timed passes, a
100,000-event warmup, and a 5,000-base transcript window. Both revisions
are clean, in-tree release rebuilds.

| revision | workload                                       | min_seconds | median_seconds | max_seconds |
|:---------|:-----------------------------------------------|------------:|---------------:|------------:|
| 60a15f21 | fixture_one_transcript_sorted_breakend         |       0.266 |          0.267 |       0.272 |
| 9ed772f7 | fixture_one_transcript_sorted_breakend         |       0.266 |          0.270 |       0.272 |
| 9ed772f7 | fixture_one_transcript_sorted_breakend_raw_alt |       0.394 |          0.396 |       0.402 |

The identical typed-endpoint control changes from **0.267 s to 0.270 s**
(+1.1% median; both observed ranges 0.266–0.272 s). Adding
`--fixture-breakend-alt` measures **0.396 s** (0.394–0.402 s), including
all four raw bracket forms, `duckvep_breakend_geometry()`, the explicit
mate-name join, ordering, annotation, unnesting and aggregation. It is
additional preparation work, not a before/after implementation
comparison. The untimed full compact-row fingerprints match for all
three runs; no consequence, status, coordinate or input-identity field
is dropped from those fingerprints.

The model is the checked-in one-transcript, two-exon fixture, with local
position 159 and mate position 170. This is a paired-adapter floor, not
a measurement of large raw records, high-fan-out models, fusion
reconstruction, or phased work. It does not measure the newly corrected
empty-endpoint predicate state; fixed native/SQL/R witnesses and
executable VEP differentials cover that state. No prior checked-in
workload measured raw ALT preparation. The production BND baseline below
(200,000 events / 18,766,240 transcript rows) has a different model and
fan-out and is not a same-workload comparator for this fixture.

| lane                                     | revision | paired_events | consequence_rows | rows_per_event | median_seconds | paired_events_per_second | consequence_rows_per_second |
|:-----------------------------------------|:---------|:--------------|:-----------------|---------------:|---------------:|:-------------------------|:----------------------------|
| transcript baseline                      | 360619ed | 200,000       | 18,766,240       |          93.83 |          2.710 | 73,801                   | 6,924,812                   |
| rejected full-result sort                | f97101e1 | 200,000       | 18,766,240       |          93.83 |          3.105 | 64,412                   | 6,043,878                   |
| pre-literal transcript A/B               | e7c3623d | 200,000       | 18,766,240       |          93.83 |          2.642 | 75,700                   | 7,103,043                   |
| pre-literal transcript + regulation A/B  | e7c3623d | 200,000       | 18,877,120       |          94.39 |          2.725 | 73,394                   | 6,927,383                   |
| post-literal transcript A/B              | 2548d675 | 200,000       | 18,766,240       |          93.83 |          2.681 | 74,599                   | 6,999,717                   |
| post-literal transcript + regulation A/B | 2548d675 | 200,000       | 18,877,120       |          94.39 |          2.761 | 72,438                   | 6,837,059                   |

This measurement used the complete receipt-hashed GRCh38 model, compact
output, one DuckDB thread, 15 timed passes, and a 10,000-event warmup.
Model loading and staging were outside the timed region. The historical
baseline and rejected-sort rows used `taskset -c 2`; the direct
pre-literal/post-literal A/B rows were rebuilt from their exact
revisions and run consecutively with `taskset -c 4`. Here “post-literal”
is revision 2548d675, not the report’s current source; source `e25c151`
has not rerun this identical BND workload. The site-rate headline alone
would be misleading: the transcript-only workload emits 93.83 rows per
event, and the checksum covers every emitted mask. The `f97101e1`
intermediate sorted all 18.8 million already variant-major transcript
rows and was rejected. The accepted comparison revision instead sorts
only the smaller feature stream and merges the two sorted streams
linearly.

The direct transcript-only A/B is 2.681 seconds versus 2.642 seconds
before the long-literal predicate fix (+1.5%). The corresponding
regulation A/B is 2.761 versus 2.725 seconds (+1.3%). Loading 1,383,580
regulation/motif intervals and emitting 110,880 additional rows costs
+3.0% relative to the post-literal transcript-only lane and preserves
6.84 million output rows/s.

An earlier run of the comparison revision pinned to CPU 2 reported a
3.386-second median while its fastest pass still completed in 2.712
seconds. Host inspection found three long-running unpinned CPU-bound
workers, and a clean CPU-4 rerun returned to the 2.6-second range. That
contaminated row is not retained in the ledger. The remaining 1–2%
direct A/B difference is reported without a no-regression claim.

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

## Last measured transcript-only GRCh38 model resources

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

This resource row predates the resident regulation/motif relation and
remains a transcript-only baseline; the throughput table above records
the feature count for regulation-enabled runs. The artifact size is a
compact DuckDB containing only `model_regions`, `model_transcripts`, and
`model_receipt`; benchmark inputs are excluded. RSS is for one R/DuckDB
process, not a C-heap estimate. The load peak includes both the DuckDB
scan/cache and construction of the second immutable C representation.
Dropping the model releases the C representation, but DuckDB retains its
table cache in the process. The planned per-worker number is the exact
sum of the C buffer capacities (active and candidate transcript indices,
one exon cursor per transcript, and worst-case edit/CDS/peptide
scratch); allocator metadata is not included. The much smaller
first-workspace RSS change reflects lazy or reused pages, not a smaller
capacity.

The historical standalone model ABI 0.12 cache adds three `uint32_t`
arrays and one `uint8_t` array per transcript. Its exact element payload
for this model is 8,377,551 bytes (7.989 MiB), plus four allocator
records. This deterministic allocation is not retroactively included in
the older `b204dd49` RSS measurement above.

The production form is run with an externally pinned process when a
stable CPU comparison is required:

``` sh
taskset -c 2 Rscript benchmarks/duckvep_throughput.R \
  --database /path/to/model.duckdb \
  --variants-database /path/to/staged-corpus.duckdb \
  --variants-table bench_variants \
  --corpus-source /path/to/source.vcf.gz \
  --workload-name ensembl116_grch38_full_corpus_regulation \
  --regulatory --variants 4095611 --passes 5 --warmup 100000 \
  --threads 1 --input-partitions 1 --transcript-distance 5000 \
  --output compact --fingerprint /tmp/full-public-row-fingerprint.csv
```

For cumulative HGVS, add `--reference-fasta` and change `--output` to
`hgvs`. The source revision must be the current checkout; the history
row records the extension binary, physical/logical model, model source,
reference, staged corpus, original source, and public-row fingerprint
receipts.

## Real whole-genome composition run

The resident-engine tables above isolate comparable annotation surfaces.
The following separate run verifies the complete relational composition
shown in the root README on real human data. It scans the public HG002
40x PCR-free DeepVariant GRCh38 WGS callset sequentially, retains
literal alleles on chromosomes 1–22, X, Y, and MT, loads the complete
Ensembl 116 GRCh38 transcript plus regulation/motif model, requests rich
consequence and HGVSc/HGVSp, and joins dated ClinVar, ClinvArbitration,
AlphaMissense, gnomAD v2.1.1 gene constraint, and an Ensembl regulatory
interval Parquet relation. The final relation is written as ZSTD Parquet
rather than reduced to a scalar.

| stage                                       | result rows | seconds | threads | peak RSS (GiB) |
|:--------------------------------------------|:------------|--------:|--------:|:---------------|
| load Ensembl 116 model                      | 644,427     |   2.398 |       4 |                |
| read VCF, canonicalize and sort             | 7,378,240   |  12.200 |       4 |                |
| ClinVar exact join                          | 50,749      |   0.259 |       4 |                |
| ClinvArbitration exact join                 | 45,232      |   0.176 |       4 |                |
| AlphaMissense exact join                    | 14,850      |   0.523 |       4 |                |
| regulatory RegionKey IEJoin + write         | 414,813     |   0.770 |       4 | 1.52           |
| regulatory cgranges build/query + write     | 414,813     |   1.144 |       4 | 0.79           |
| rich + HGVS + core regulation, four writers | 88,392,840  |  24.524 |       4 | 4.98           |
| rich + HGVS + all providers, four writers   | 88,392,840  |  28.841 |       4 | 5.31           |

These are warm-page-cache, unpinned integration measurements, not
replacements for the pinned resident-engine comparisons above. DuckDB
CLI `.timer on` records the statements and GNU `/usr/bin/time -v`
records process high-water RSS. The corrected production query sets a 4
GB DuckDB memory limit and writes one Parquet file per worker so the
complete consequence relation is never retained or globally sorted.

| method                                     | matched alleles | overlap pairs | seconds | peak RSS (MiB) |
|:-------------------------------------------|:----------------|:--------------|--------:|:---------------|
| chromosome hash join + range residuals     | 414,813         | 745,252       |  83.960 |                |
| two packed RegionKey inequalities (IEJoin) | 414,813         | 745,252       |   0.770 | 1557           |
| cgranges bulk index + query                | 414,813         | 745,252       |   1.144 | 812            |

The first query was mislabeled as IEJoin in the initial receipt.
`EXPLAIN` showed a chromosome `HASH_JOIN` with all four range predicates
applied to the same-chromosome candidate products. The receipt remains
above rather than being hidden. Removing the redundant string equality
lets the two packed RegionKey inequalities encode both chromosome and
half-open overlap and selects `IE_JOIN`. It returns the same 745,252
pairs and 414,813 matched alleles in 109.04x less wall time. The
cgranges path is slightly slower here but uses roughly half the process
memory and supports arbitrary literal contig names.

The bounded four-writer composition emits 88,392,840 rows in 28.841
seconds and writes 285,798,217 bytes. Peak RSS is 5.31 GiB; the
high-water mark includes transient construction of the complete
immutable Ensembl model, whose C allocations are outside DuckDB’s 4 GB
buffer-manager limit.

The earlier retained-intermediate/single-file run is also preserved in
the CSV: it took 284.890 seconds and peaked at 46.87 GiB. The bounded
output has the same row count and the same order-independent XOR and sum
fingerprints over every projected column. Controlled resident and
FastVEP comparisons remain the authorities above.
