DuckVEP SQL transcript presentation
================

This measures the complete 23-column result of
`duckvep_transcript_projection(events, annotations, transcripts)` as a
DuckDB temporary table. It is the SQL reference, not a fused native
implementation. The literal alleles and consequence rows are prepared
before timing. The timer includes the macro’s validation, joins,
coordinate projection, codon/peptide calculation and complete output
materialization. One warm-up precedes the reported measured passes;
preparation and full-output verification are excluded from those
timings.

The raw FASTA, index and minimal GFF are pinned in `r/duckhtsbench`
under the `duckvep-projection` workload. `projection_fixtures.R` derives
a forward model, a reverse model, and a three-exon model with two-base
CDS phase padding from those inputs. The unchanged independent-event
generator adds 1,000 seed-173 random alleles to its targeted witnesses.
Only record IDs are relabelled so repeated source alleles retain
separate event keys. Each set is repeated eight times with distinct
keys. These tiny warm models are a reproducible adapter workload, not a
genome-scale annotation or file-ingestion benchmark.

Each model/worker-count pair runs in a fresh R process. Peak RSS is the
Linux process high-water mark, including model preparation and
verification; it is not isolated DuckDB workspace usage. Output bytes
are the sum of the UTF-8 lengths of the complete rows serialized as
JSON, not DuckDB physical storage size. The SHA-256 covers the schema
and every typed row sorted by all columns, including duplicates and
NULLs. It must match across passes and worker counts.

## Measured materialization

| source_revision | case              | threads | input_alleles | expanded_transcript_rows | output_rows | json_utf8_bytes | iterations | median_seconds | min_seconds | max_seconds | peak_rss_bytes |
|:----------------|:------------------|--------:|--------------:|-------------------------:|------------:|----------------:|-----------:|---------------:|------------:|------------:|---------------:|
| f0184155b0fa    | forward           |       1 |         10144 |                    10144 |       10144 |         4618058 |          3 |          0.045 |       0.045 |       0.045 |      253714432 |
| f0184155b0fa    | forward           |       4 |         10144 |                    10144 |       10144 |         4618058 |          3 |          0.030 |       0.030 |       0.032 |      285921280 |
| f0184155b0fa    | reverse           |       1 |         10112 |                    10112 |       10112 |         4603562 |          3 |          0.046 |       0.045 |       0.047 |      258449408 |
| f0184155b0fa    | reverse           |       4 |         10112 |                    10112 |       10112 |         4603562 |          3 |          0.032 |       0.031 |       0.032 |      288108544 |
| f0184155b0fa    | three_exon_phase2 |       1 |         11112 |                    11112 |       11112 |         5056410 |          3 |          0.049 |       0.049 |       0.051 |      243851264 |
| f0184155b0fa    | three_exon_phase2 |       4 |         11112 |                    11112 |       11112 |         5056410 |          3 |          0.033 |       0.033 |       0.033 |      284536832 |
| 20b2902c4672    | forward           |       1 |         10144 |                    10144 |       10144 |         4618058 |          3 |          0.046 |       0.045 |       0.048 |      256720896 |
| 20b2902c4672    | forward           |       4 |         10144 |                    10144 |       10144 |         4618058 |          3 |          0.032 |       0.031 |       0.032 |      285495296 |
| 20b2902c4672    | reverse           |       1 |         10112 |                    10112 |       10112 |         4603562 |          3 |          0.046 |       0.045 |       0.047 |      255000576 |
| 20b2902c4672    | reverse           |       4 |         10112 |                    10112 |       10112 |         4603562 |          3 |          0.032 |       0.031 |       0.032 |      273793024 |
| 20b2902c4672    | three_exon_phase2 |       1 |         11112 |                    11112 |       11112 |         5056410 |          3 |          0.050 |       0.050 |       0.051 |      245338112 |
| 20b2902c4672    | three_exon_phase2 |       4 |         11112 |                    11112 |       11112 |         5056410 |          3 |          0.033 |       0.032 |       0.033 |      283779072 |

## Identity and complete output

|     | source_revision                          | extension_sha256                                                 | build_binding                 | machine                                                    | cpu                                 | duckdb_version |
|:----|:-----------------------------------------|:-----------------------------------------------------------------|:------------------------------|:-----------------------------------------------------------|:------------------------------------|:---------------|
| 1   | f0184155b0faa2486937c8b774648ce897c88f2f | 420e6d1b14a2831384117d2bec74cba72732769a4349ff9866eeb409d526c7d9 | htslib_distclean_make_release | Linux 6.8.0-78-generic x86_64 Ubuntu-2404-noble-amd64-base | 13th Gen Intel(R) Core(TM) i5-13500 | 1.5.3          |
| 7   | 20b2902c4672fb0401616a7f09ff8fc0e703909e | 420e6d1b14a2831384117d2bec74cba72732769a4349ff9866eeb409d526c7d9 | htslib_distclean_make_release | Linux 6.8.0-78-generic x86_64 Ubuntu-2404-noble-amd64-base | 13th Gen Intel(R) Core(TM) i5-13500 | 1.5.3          |

|     | case              | input_vcf_sha256                                                 | model_gff_sha256                                                 | reference_sha256                                                 | output_sha256                                                    |
|:----|:------------------|:-----------------------------------------------------------------|:-----------------------------------------------------------------|:-----------------------------------------------------------------|:-----------------------------------------------------------------|
| 1   | forward           | 3f79dffd03a5d552e89d15a3eb393c0dd2428ef86cf2f02ce32cf8ddbd1dfe3e | 41defe13bfea82d43afc45dd1016c678fc199e634eab215309477ee0bc821685 | 01d1f025213063a747cb0c53cdcbea67262f32ab1dbdb56aa4131d72142e6d26 | aa92aa02feb723b06fea518580209008730d8e7c11992341a1030f787c980e52 |
| 3   | reverse           | da6bbb2e1cd5c354a62a831aa37bcbef417dbefcfe4f55260650be03dc72cb7f | 7a30f5acd27936582d7a1856e95baea8b92de01e2a6c93312c915a7d16fb8c02 | 01d1f025213063a747cb0c53cdcbea67262f32ab1dbdb56aa4131d72142e6d26 | 3f297fdf19db2876b614b01cc77077aa82394810fd38fc17eacc2a686c1f5a3a |
| 5   | three_exon_phase2 | d7d75e06f61766d825b288f4846e827db1cedd7193c9ce722e59d0bb27d756ae | 4e554456fd7f761f54ea25ac8def4bcfd92dc8810e69ea8058e2d7c2917c6dbc | 01d1f025213063a747cb0c53cdcbea67262f32ab1dbdb56aa4131d72142e6d26 | 4c983975b15fe4d200158a9c1d87f25ce874345bdaeba45431679338744e4ce5 |

The first `f0184155b0fa` runs establish this complete-contract baseline.
Later rows retain that measurement rather than replace it. The following
comparison pairs the first and last recorded runs only when input hash,
model, output hash, copy count, worker count, DuckDB version and machine
agree.

| case              | threads | baseline_revision | current_revision | baseline_seconds | current_seconds | current_over_baseline |
|:------------------|--------:|:------------------|:-----------------|-----------------:|----------------:|----------------------:|
| forward           |       1 | f0184155b0fa      | 20b2902c4672     |            0.045 |           0.046 |              1.022222 |
| forward           |       4 | f0184155b0fa      | 20b2902c4672     |            0.030 |           0.032 |              1.066667 |
| reverse           |       1 | f0184155b0fa      | 20b2902c4672     |            0.046 |           0.046 |              1.000000 |
| reverse           |       4 | f0184155b0fa      | 20b2902c4672     |            0.032 |           0.032 |              1.000000 |
| three_exon_phase2 |       1 | f0184155b0fa      | 20b2902c4672     |            0.049 |           0.050 |              1.020408 |
| three_exon_phase2 |       4 | f0184155b0fa      | 20b2902c4672     |            0.033 |           0.033 |              1.000000 |

The [sorted-stream throughput report](duckvep_throughput.md) instead
measures consequence/HGVS output and is not an identical comparator. No
fused implementation is measured here. A production-density model and
comparison with a future fused implementation remain unmeasured; small
timing differences on this warm adapter workload are not proof of a
speedup or absence of regression.

## Pinned VEP field differential

The separate differential checks all 21 presentation fields, plus the
complete event/transcript key set, against VEP 116 JSON. Six models add
a partial-CDS end, a noncoding first exon and a noncoding transcript to
the three timed models. Each seed runs the unchanged generator’s
targeted cases plus 1,000 random cases; the 128-base differing-allele
campaign supplements the default 10-base campaign to exercise longer
exon/intron spans. Every physical input record is retained, including
duplicate alleles. Missing or duplicate keys, changed fields, and
changes to/from NULL must all trip deliberate verifier controls. Actual
and expected tables, original and relabelled VCFs, GFFs and oracle JSON
have separate hashes in the checked-in CSV; detailed artifacts remain in
the run directory.

| source_revision                          | case                 |     seed | max_random_length | input_alleles | expected_pairs | actual_pairs | mismatches | detected_controls |
|:-----------------------------------------|:---------------------|---------:|------------------:|--------------:|---------------:|-------------:|-----------:|------------------:|
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | forward              |      173 |                10 |          1268 |           1268 |         1268 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | reverse              |      173 |                10 |          1264 |           1264 |         1264 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | three_exon_phase2    |      173 |                10 |          1389 |           1389 |         1389 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | partial_cds_end      |      173 |                10 |          1385 |           1385 |         1385 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | noncoding_first_exon |      173 |                10 |          1389 |           1389 |         1389 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | noncoding            |      173 |                10 |          1198 |           1198 |         1198 |          0 |                54 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | forward              | 20260906 |                10 |          1268 |           1268 |         1268 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | reverse              | 20260906 |                10 |          1264 |           1264 |         1264 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | three_exon_phase2    | 20260906 |                10 |          1389 |           1389 |         1389 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | partial_cds_end      | 20260906 |                10 |          1385 |           1385 |         1385 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | noncoding_first_exon | 20260906 |                10 |          1389 |           1389 |         1389 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | noncoding            | 20260906 |                10 |          1198 |           1198 |         1198 |          0 |                54 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | forward              |      173 |               128 |          1268 |           1268 |         1268 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | reverse              |      173 |               128 |          1264 |           1264 |         1264 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | three_exon_phase2    |      173 |               128 |          1389 |           1389 |         1389 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | partial_cds_end      |      173 |               128 |          1385 |           1385 |         1385 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | noncoding_first_exon |      173 |               128 |          1389 |           1389 |         1389 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | noncoding            |      173 |               128 |          1198 |           1198 |         1198 |          0 |                54 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | forward              | 20260906 |               128 |          1268 |           1268 |         1268 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | reverse              | 20260906 |               128 |          1264 |           1264 |         1264 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | three_exon_phase2    | 20260906 |               128 |          1389 |           1389 |         1389 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | partial_cds_end      | 20260906 |               128 |          1385 |           1385 |         1385 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | noncoding_first_exon | 20260906 |               128 |          1389 |           1389 |         1389 |          0 |                62 |
| 20b2902c4672fb0401616a7f09ff8fc0e703909e | noncoding            | 20260906 |               128 |          1198 |           1198 |         1198 |          0 |                54 |

The oracle sources are Ensembl VEP
`57ea5c52340acc1f156267f810ad162e26597082` and Variation
`2fb834b987ede3824e200197a838ce11e91aeb4b`, with the exact package set
in
`test/duckvep/upstream/receipts/vep116_2026-07-22.conda-explicit.txt`.
GFF input does not import the full core transcript-quality and
translation-edit attributes. Positive CDS quality flags and reference
peptide edits are therefore tested as SQL model-attribute properties,
not represented as positive GFF oracle coverage. Symbolic structural and
phased events are outside this function’s literal-allele contract;
existing consequence/HGVS/phased acceptance histories and generator
settings are unchanged.

Reproduce from a clean source revision with an extension build receipt:

``` r
source("scripts/benchmark_duckvep_projection.R")
benchmark_duckvep_projection(
  normalizePath("."), "build/release/duckhts.duckdb_extension",
  "benchmarks/data/duckvep_projection.csv",
  extension_receipt = Sys.getenv("DUCKHTS_EXTENSION_RECEIPT")
)
rmarkdown::render("benchmarks/benchmark_duckvep_projection.Rmd",
  knit_root_dir = normalizePath("."))
```

For each reported seed and maximum length, record conformance
separately:

``` bash
make test-duckvep-projection ARGS="--seed 173 --random-cases 1000 --max-random-length 10 --extension-receipt $DUCKHTS_EXTENSION_RECEIPT --summary-output benchmarks/data/duckvep_projection_conformance.csv"
```

`VEP_PREFIX` names the environment matching the checked-in package lock.
