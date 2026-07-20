DuckVEP sorted-stream throughput
================

<!-- duckvep_throughput.md is generated from duckvep_throughput.Rmd. -->

This benchmark drives the rich `duckvep_annotate`, numeric
`duckvep_annotate_compact`, or independent-event `duckvep_annotate_hgvs`
result through DuckDB’s stable extension API with nondecreasing
sequence-region and position keys. It records candidate sweep,
consequence classification, list materialization, `unnest`, and final
aggregation together. Model loading and input staging are outside the
timed pass.

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

| run_date   | revision | workload                                             | output_mode | threads | input_partitions | transcript_distance | variants   | transcripts | exons     | regulation_features | annotated_rows | passes | min_seconds | median_seconds | max_seconds | variants_per_second | ns_per_variant | cpu                                 | cpu_affinity | output_rows_per_second |
|:-----------|:---------|:-----------------------------------------------------|:------------|--------:|-----------------:|--------------------:|:-----------|:------------|:----------|:--------------------|:---------------|-------:|------------:|---------------:|------------:|--------------------:|---------------:|:------------------------------------|:-------------|:-----------------------|
| 2026-07-11 | 8cc22218 | fixture_one_transcript_sorted                        | rich        |       1 |                1 |                5000 | 10,000,000 | 1           | 2         | 0                   | 10,000,000     |      3 |       3.225 |          3.334 |       3.371 |             2999400 |          333.4 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 2,999,400              |
| 2026-07-13 | 87f03a2a | fixture_one_transcript_sorted                        | rich        |       1 |                1 |                5000 | 10,000,000 | 1           | 2         | 0                   | 10,000,000     |      3 |       2.812 |          2.824 |       2.825 |             3541076 |          282.4 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,541,076              |
| 2026-07-14 | 751d3c74 | ensembl116_grch38_giab_sites_hash40                  | compact     |       1 |                1 |                5000 | 100,957    | 644,292     | 5,078,384 | 0                   | 1,179,465      |      9 |       0.104 |          0.107 |       0.114 |              943523 |         1059.9 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 11,023,037             |
| 2026-07-14 | 751d3c74 | ensembl116_grch38_giab_sites_hash40                  | rich        |       1 |                1 |                5000 | 100,957    | 644,292     | 5,078,384 | 0                   | 1,179,465      |      9 |       0.222 |          0.224 |       0.231 |              450701 |         2218.8 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 5,265,469              |
| 2026-07-14 | 79cfaaf2 | ensembl116_grch38_giab_sites_hash40                  | compact     |       1 |                1 |                5000 | 100,957    | 644,292     | 5,078,384 | 0                   | 1,179,465      |      9 |       0.107 |          0.108 |       0.110 |              934787 |         1069.8 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 10,920,972             |
| 2026-07-14 | 79cfaaf2 | ensembl116_grch38_giab_sites_hash40                  | rich        |       1 |                1 |                5000 | 100,957    | 644,292     | 5,078,384 | 0                   | 1,179,465      |      9 |       0.225 |          0.228 |       0.232 |              442794 |         2258.4 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 5,173,092              |
| 2026-07-14 | 81b9e2a7 | ensembl116_grch38_primary_giab_sites_hash40          | compact     |       1 |                1 |                5000 | 100,957    | 646,577     | 5,087,789 | 0                   | 1,179,465      |      9 |       0.112 |          0.115 |       0.125 |              877887 |         1139.1 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 10,256,217             |
| 2026-07-14 | 81b9e2a7 | ensembl116_grch38_primary_giab_sites_hash40          | rich        |       1 |                1 |                5000 | 100,957    | 646,577     | 5,087,789 | 0                   | 1,179,465      |      9 |       0.233 |          0.236 |       0.242 |              427784 |         2337.6 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 4,997,733              |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_coding_mixed_200k            | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,757,180      |      3 |      58.833 |         59.010 |      59.212 |                3389 |       295050.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 97,563                 |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_coding_nonsnv_36k            | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,616      |      5 |      11.837 |         11.898 |      11.959 |                3026 |       330500.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 88,050                 |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_coding_nonsnv_36k            | rich        |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,616      |      3 |      11.903 |         11.952 |      11.978 |                3012 |       332000.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 87,652                 |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_coding_snv_200k              | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |      9 |       0.690 |          0.693 |       0.719 |              288600 |         3465.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 7,490,620              |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_coding_snv_200k              | rich        |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |      5 |       1.213 |          1.238 |       1.262 |              161551 |         6190.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 4,193,053              |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_giab_sites_hash40            | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 0                   | 1,174,245      |      9 |       0.116 |          0.117 |       0.119 |              862880 |         1158.9 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 10,036,282             |
| 2026-07-14 | b204dd49 | ensembl116_grch38_final_giab_sites_hash40            | rich        |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 0                   | 1,174,245      |      9 |       0.239 |          0.244 |       0.263 |              413758 |         2416.9 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 4,812,480              |
| 2026-07-14 | bcca6f6c | ensembl116_grch38_giab_sites_hash40                  | compact     |       1 |                1 |                5000 | 100,957    | 644,292     | 5,078,384 | 0                   | 1,179,465      |      9 |       0.106 |          0.107 |       0.109 |              943523 |         1059.9 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 11,023,037             |
| 2026-07-14 | bcca6f6c | ensembl116_grch38_giab_sites_hash40                  | rich        |       1 |                1 |                5000 | 100,957    | 644,292     | 5,078,384 | 0                   | 1,179,465      |      9 |       0.221 |          0.225 |       0.231 |              448698 |         2228.7 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 5,242,067              |
| 2026-07-15 | 1903f7fb | ensembl116_grch38_final_coding_mixed_200k            | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |     21 |       1.489 |          1.514 |       1.618 |              132100 |         7570.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,802,325              |
| 2026-07-15 | 1903f7fb | ensembl116_grch38_final_coding_nonsnv_36k            | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,524      |     31 |       0.323 |          0.329 |       0.362 |              109422 |         9138.9 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,183,964              |
| 2026-07-15 | 1903f7fb | ensembl116_grch38_final_coding_snv_200k              | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |     31 |       0.610 |          0.622 |       0.657 |              321543 |         3110.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 8,345,659              |
| 2026-07-15 | 57b5b566 | ensembl116_grch38_final_coding_mixed_200k            | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |     15 |       1.467 |          1.480 |       1.489 |              135135 |         7400.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,889,676              |
| 2026-07-15 | 57b5b566 | ensembl116_grch38_final_coding_nonsnv_36k            | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,524      |     21 |       0.310 |          0.322 |       0.351 |              111801 |         8944.4 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,253,180              |
| 2026-07-15 | 57b5b566 | ensembl116_grch38_final_coding_snv_200k              | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |     15 |       0.609 |          0.615 |       0.638 |              325203 |         3075.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 8,440,650              |
| 2026-07-15 | 7dd90ce8 | ensembl116_grch38_final_coding_mixed_200k            | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |      9 |       1.550 |          1.572 |       1.583 |              127226 |         7860.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,662,036              |
| 2026-07-15 | 7dd90ce8 | ensembl116_grch38_final_coding_nonsnv_36k            | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,524      |     11 |       0.341 |          0.344 |       0.357 |              104651 |         9555.6 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,045,128              |
| 2026-07-15 | 7dd90ce8 | ensembl116_grch38_final_coding_snv_200k              | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |      9 |       0.609 |          0.613 |       0.629 |              326264 |         3065.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 8,468,189              |
| 2026-07-15 | 7dd90ce8 | ensembl116_grch38_final_giab_sites_hash40            | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 0                   | 1,174,245      |      9 |       0.118 |          0.119 |       0.137 |              848378 |         1178.7 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 9,867,605              |
| 2026-07-15 | a346a7aa | ensembl116_grch38_final_coding_mixed_200k            | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |     21 |       1.417 |          1.443 |       1.511 |              138600 |         7215.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,989,411              |
| 2026-07-15 | a346a7aa | ensembl116_grch38_final_coding_nonsnv_36k            | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,524      |     31 |       0.307 |          0.315 |       0.376 |              114286 |         8750.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,325,473              |
| 2026-07-15 | a346a7aa | ensembl116_grch38_final_coding_snv_200k              | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |     31 |       0.599 |          0.617 |       0.791 |              324149 |         3085.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 8,413,290              |
| 2026-07-15 | c361346f | ensembl116_grch38_final_coding_mixed_200k            | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |     21 |       1.445 |          1.468 |       1.498 |              136240 |         7340.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,921,471              |
| 2026-07-15 | c361346f | ensembl116_grch38_final_coding_nonsnv_36k            | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,524      |     31 |       0.320 |          0.323 |       0.341 |              111455 |         8972.2 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,243,108              |
| 2026-07-15 | c361346f | ensembl116_grch38_final_coding_snv_200k              | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |     31 |       0.584 |          0.602 |       0.779 |              332226 |         3010.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 8,622,924              |
| 2026-07-15 | ecb7eaae | ensembl116_grch38_final_coding_mixed_200k            | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |     15 |       1.467 |          1.495 |       1.512 |              133779 |         7475.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,850,649              |
| 2026-07-15 | ecb7eaae | ensembl116_grch38_final_coding_nonsnv_36k            | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,524      |     21 |       0.311 |          0.324 |       0.356 |              111111 |         9000.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,233,099              |
| 2026-07-15 | ecb7eaae | ensembl116_grch38_final_coding_snv_200k              | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |     15 |       0.605 |          0.610 |       0.624 |              327869 |         3050.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 8,509,836              |
| 2026-07-15 | f649d724 | ensembl116_grch38_final_coding_mixed_200k            | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |      9 |       1.487 |          1.504 |       1.626 |              132979 |         7520.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,827,606              |
| 2026-07-15 | f649d724 | ensembl116_grch38_final_coding_nonsnv_36k            | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,524      |     11 |       0.324 |          0.327 |       0.342 |              110092 |         9083.3 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,203,437              |
| 2026-07-15 | f649d724 | ensembl116_grch38_final_coding_snv_200k              | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |     11 |       0.611 |          0.613 |       0.634 |              326264 |         3065.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 8,468,189              |
| 2026-07-15 | f649d724 | ensembl116_grch38_final_giab_sites_hash40            | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 0                   | 1,174,245      |     15 |       0.118 |          0.120 |       0.132 |              841308 |         1188.6 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 9,785,375              |
| 2026-07-16 | 220c6ddc | ensembl116_grch38_final_giab_sites_hash40            | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 0                   | 1,174,245      |     31 |       0.118 |          0.120 |       0.144 |              841308 |         1188.6 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 9,785,375              |
| 2026-07-16 | 220c6ddc | ensembl116_grch38_final_giab_sites_hash40_regulation | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 1,383,580           | 1,179,329      |     31 |       0.125 |          0.126 |       0.162 |              801246 |         1248.1 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 9,359,754              |
| 2026-07-16 | 360619ed | ensembl116_grch38_final_breakend_200k                | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 18,766,240     |     15 |       2.676 |          2.710 |       2.727 |               73801 |        13550.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 6,924,812              |
| 2026-07-16 | 96b4cd45 | ensembl116_grch38_final_giab_sites_hash40            | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 0                   | 1,174,245      |     31 |       0.131 |          0.134 |       0.147 |              753410 |         1327.3 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 8,763,022              |
| 2026-07-16 | 96b4cd45 | ensembl116_grch38_final_giab_sites_hash40_regulation | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 1,383,580           | 1,179,329      |     31 |       0.138 |          0.140 |       0.162 |              721121 |         1386.7 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 8,423,779              |
| 2026-07-16 | 9c97cb07 | ensembl116_grch38_final_coding_mixed_200k            | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |     15 |       1.451 |          1.481 |       1.638 |              135044 |         7405.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,887,049              |
| 2026-07-16 | 9c97cb07 | ensembl116_grch38_final_coding_nonsnv_36k            | compact     |       1 |                1 |                5000 | 36,000     | 644,427     | 5,068,416 | 0                   | 1,047,524      |     15 |       0.317 |          0.322 |       0.339 |              111801 |         8944.4 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 3,253,180              |
| 2026-07-16 | 9c97cb07 | ensembl116_grch38_final_coding_snv_200k              | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,191,000      |     15 |       0.593 |          0.603 |       0.709 |              331675 |         3015.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 8,608,624              |
| 2026-07-16 | 9c97cb07 | ensembl116_grch38_final_giab_sites_hash40            | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 0                   | 1,174,245      |     15 |       0.117 |          0.118 |       0.131 |              855568 |         1168.8 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 9,951,229              |
| 2026-07-16 | f3796ae9 | ensembl116_grch38_final_giab_sites_hash40            | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 0                   | 1,174,245      |     31 |       0.130 |          0.132 |       0.141 |              764826 |         1307.5 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 8,895,795              |
| 2026-07-16 | f3796ae9 | ensembl116_grch38_final_giab_sites_hash40_regulation | compact     |       1 |                1 |                5000 | 100,957    | 644,427     | 5,068,416 | 1,383,580           | 1,179,329      |     31 |       0.138 |          0.139 |       0.148 |              726309 |         1376.8 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 8,484,381              |
| 2026-07-19 | 2548d675 | ensembl116_grch38_final_breakend_200k                | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 18,766,240     |     15 |       2.676 |          2.681 |       2.693 |               74599 |        13405.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 6,999,717              |
| 2026-07-19 | 2548d675 | ensembl116_grch38_final_breakend_regulation_200k     | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 1,383,580           | 18,877,120     |     15 |       2.737 |          2.761 |       2.769 |               72438 |        13805.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 6,837,059              |
| 2026-07-19 | e7c3623d | ensembl116_grch38_final_breakend_200k                | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 18,766,240     |     15 |       2.589 |          2.642 |       2.654 |               75700 |        13210.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 7,103,043              |
| 2026-07-19 | e7c3623d | ensembl116_grch38_final_breakend_regulation_200k     | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 1,383,580           | 18,877,120     |     15 |       2.712 |          2.725 |       2.745 |               73394 |        13625.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 6,927,383              |
| 2026-07-19 | f97101e1 | ensembl116_grch38_final_breakend_200k                | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 18,766,240     |     15 |       3.056 |          3.105 |       4.435 |               64412 |        15525.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 6,043,878              |
| 2026-07-19 | f97101e1 | ensembl116_grch38_final_breakend_regulation_200k     | compact     |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 1,383,580           | 18,877,120     |     15 |       3.307 |          3.348 |       3.370 |               59737 |        16740.0 | 13th Gen Intel(R) Core(TM) i5-13500 | NA           | 5,638,327              |
| 2026-07-20 | 0714235a | ensembl116_grch38_final_coding_mixed_200k            | hgvs        |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |      9 |      32.023 |         32.110 |      32.150 |                6229 |       160550.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 179,281                |
| 2026-07-20 | 0714235a | ensembl116_grch38_final_coding_mixed_200k            | rich        |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |      9 |       2.031 |          2.038 |       2.047 |               98135 |        10190.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,824,691              |
| 2026-07-20 | 25328813 | ensembl116_grch38_final_coding_mixed_200k            | rich        |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |      9 |       2.024 |          2.033 |       2.046 |               98377 |        10165.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,831,638              |
| 2026-07-20 | 25328813 | ensembl116_grch38_final_coding_mixed_200k            | hgvs        |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |      9 |      32.047 |         32.106 |      32.204 |                6229 |       160530.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 179,304                |
| 2026-07-20 | 7dae50cd | ensembl116_grch38_final_coding_mixed_200k            | rich        |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |      9 |       2.045 |          2.067 |       2.081 |               96759 |        10335.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 2,785,060              |
| 2026-07-20 | 7dae50cd | ensembl116_grch38_final_coding_mixed_200k            | hgvs        |       1 |                1 |                5000 | 200,000    | 644,427     | 5,068,416 | 0                   | 5,756,720      |      9 |       4.565 |          4.579 |       4.610 |               43678 |        22895.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 1,257,200              |
| 2026-07-20 | e25c1513 | ensembl116_grch38_clinvar_annotation_dense_v2        | compact     |       1 |                1 |                   0 | 517,097    | 644,427     | 5,068,416 | 1,383,580           | 18,163,993     |      5 |       1.828 |          1.834 |       1.845 |              281950 |         3546.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 9,904,031              |
| 2026-07-20 | e25c1513 | ensembl116_grch38_clinvar_annotation_dense_v2        | compact     |       4 |                4 |                   0 | 517,097    | 644,427     | 5,068,416 | 1,383,580           | 18,163,993     |      5 |       0.535 |          0.540 |       0.545 |              957587 |         1044.3 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 33,637,024             |
| 2026-07-20 | e25c1513 | ensembl116_grch38_clinvar_annotation_dense_v2        | compact     |       1 |                1 |                5000 | 517,097    | 644,427     | 5,068,416 | 1,383,580           | 26,518,787     |      5 |       2.081 |          2.106 |       2.112 |              245535 |         4072.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 12,592,017             |
| 2026-07-20 | e25c1513 | ensembl116_grch38_clinvar_annotation_dense_v2        | compact     |       4 |                4 |                5000 | 517,097    | 644,427     | 5,068,416 | 1,383,580           | 26,518,787     |      5 |       0.626 |          0.632 |       0.662 |              818191 |         1222.2 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 41,960,106             |
| 2026-07-20 | e25c1513 | ensembl116_grch38_clinvar_annotation_dense_v2        | compact     |       1 |                1 |               10000 | 517,097    | 644,427     | 5,068,416 | 1,383,580           | 34,248,323     |      5 |       2.359 |          2.370 |       2.418 |              218184 |         4583.3 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 14,450,769             |
| 2026-07-20 | e25c1513 | ensembl116_grch38_clinvar_annotation_dense_v2        | compact     |       4 |                4 |               10000 | 517,097    | 644,427     | 5,068,416 | 1,383,580           | 34,248,323     |      5 |       0.713 |          0.719 |       0.758 |              719189 |         1390.5 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 47,633,273             |
| 2026-07-20 | e25c1513 | ensembl116_grch38_clinvar_annotation_dense_v2        | compact     |       1 |                1 |               50000 | 517,097    | 644,427     | 5,068,416 | 1,383,580           | 88,784,213     |      5 |       4.130 |          4.161 |       4.185 |              124272 |         8046.8 | 13th Gen Intel(R) Core(TM) i5-13500 | 2            | 21,337,230             |
| 2026-07-20 | e25c1513 | ensembl116_grch38_clinvar_annotation_dense_v2        | compact     |       4 |                4 |               50000 | 517,097    | 644,427     | 5,068,416 | 1,383,580           | 88,784,213     |      5 |       1.326 |          1.342 |       1.436 |              385318 |         2595.3 | 13th Gen Intel(R) Core(TM) i5-13500 | 2,4,6,8      | 66,158,132             |

Each pass consumes every staged input and checks output cardinality plus
either the rendered consequence-byte total or the numeric
consequence-mask sum. Rows with different workload, output mode, thread
count, host, or variant count are separate measurements and should not
be compared as if only the engine changed. `cpu_affinity` is read from
the benchmark process itself; a single integer records an externally
pinned run, while a range records a scheduler-visible CPU set.

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
  --database /path/to/staged.duckdb \
  --variants-table bench_variants \
  --workload-name ensembl116_grch38_final_giab_sites_hash40 \
  --variants 100957 --passes 9 --warmup 10000 --threads 1 \
  --output compact
```

Add `--regulatory` and use the
`ensembl116_grch38_final_giab_sites_hash40_regulation` workload name to
load `duckvep_regulation_features` and record the integrated
core-feature lane.
