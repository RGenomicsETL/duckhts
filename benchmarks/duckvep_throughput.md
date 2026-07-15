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
| 2026-07-15 | 7dd90ce8 | ensembl116_grch38_final_coding_mixed_200k   | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,756,720      |      3 |       1.560 |          1.561 |       1.566 |              128123 |         7805.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,687,841              |
| 2026-07-15 | 7dd90ce8 | ensembl116_grch38_final_coding_nonsnv_36k   | compact     |       1 | 36,000     | 644,427     | 5,068,416 | 1,047,524      |      5 |       0.339 |          0.340 |       0.347 |              105882 |         9444.4 | 13th Gen Intel(R) Core(TM) i5-13500 | 3,080,953              |
| 2026-07-15 | 7dd90ce8 | ensembl116_grch38_final_coding_snv_200k     | compact     |       1 | 200,000    | 644,427     | 5,068,416 | 5,191,000      |      9 |       0.609 |          0.613 |       0.629 |              326264 |         3065.0 | 13th Gen Intel(R) Core(TM) i5-13500 | 8,468,189              |
| 2026-07-15 | 7dd90ce8 | ensembl116_grch38_final_giab_sites_hash40   | compact     |       1 | 100,957    | 644,427     | 5,068,416 | 1,174,245      |      9 |       0.118 |          0.119 |       0.137 |              848378 |         1178.7 | 13th Gen Intel(R) Core(TM) i5-13500 | 9,867,605              |

Each pass consumes every staged input and checks output cardinality plus
either the rendered consequence-byte total or the numeric
consequence-mask sum. Rows with different workload, output mode, thread
count, host, or variant count are separate measurements and should not
be compared as if only the engine changed.

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

The production form is run with:

``` sh
Rscript benchmarks/duckvep_throughput.R \
  --database /path/to/staged.duckdb \
  --variants-table bench_variants \
  --workload-name ensembl116_grch38_final_giab_sites_hash40 \
  --variants 100957 --passes 9 --warmup 10000 --threads 1 \
  --output compact
```
