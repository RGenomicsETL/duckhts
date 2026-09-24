Reader performance across initialization compatibility changes
================

## Scope and method

This is a matched check of five reader/format combinations on one host,
with **one DuckDB thread**, automatic SIMD dispatch, and no htslib
decompression workers for BAM/VCF/BCF. Each observation starts a fresh R
process and database. All reader columns are decoded and materialized by
`CREATE TABLE AS SELECT *`; checksums and comparison snapshots are
outside the timed interval. Extension load time is recorded separately.
There are 7 measured observations per build/workload, in alternating
build order; no cold-cache or peak-memory claim is made.

Inputs resolve through `duckhtsbench`: the ONT E. coli BAM, its FASTQ
reads and reference, and the registered two-million-record synthetic
VCF/BCF. Expected counts come from the retained ONT staging evidence,
reference/read registry metadata, and synthetic generator contract. They
are not inferred from candidate output. Every repetition must match its
expected record count, complete schema, and all-field checksum.
First-run outputs additionally undergo exact multiset comparison in both
directions, including duplicate rows.

This workload does not measure DuckVEP, CRAM, indexed interval
selection, or reader scaling across thread counts. The nearest broader
reader baselines are [BCF record
caching](benchmark_bcf_record_cache.md), [multi-region
readers](benchmark_multi_region_readers.md), and [FASTQ
reading](benchmark_fastq_reader.md); their inputs/projections differ, so
the matched baseline in this report is the performance comparator.

## Artifact identity

| build     | source_revision                          | extension_sha256                                                 | duckdb | R                            | threads | measured_repetitions | run_utc                 | cpu_affinity |
|:----------|:-----------------------------------------|:-----------------------------------------------------------------|:-------|:-----------------------------|--------:|---------------------:|:------------------------|-------------:|
| baseline  | 183ae5fc5f464e68c697b990406a26a5d8138b51 | 77f881c8de8f19d457ff8227025a0cb2f1db3d2359f848807f9bfe183019c07f | 1.5.5  | R version 4.6.0 (2026-04-24) |       1 |                    7 | 2026-09-24 04:49:13 UTC |            0 |
| candidate | e602cd66796c9f1c17dce8f27cbf3fce5f779363 | e122aee4aec5d55130bd40f74c379a4e83f1b0e72600a00a95a3d09ebd886689 | 1.5.5  | R version 4.6.0 (2026-04-24) |       1 |                    7 | 2026-09-24 04:49:13 UTC |            0 |

    ## Linux / 6.8.0-78-generic / x86_64
    ## 13th Gen Intel(R) Core(TM) i5-13500

## Inputs and output denominators

| workload | artifact                     | expected_rows |    bytes | sha256                                                           |
|:---------|:-----------------------------|--------------:|---------:|:-----------------------------------------------------------------|
| BAM      | ont_ecoli_k12_bam            |         27377 | 56460404 | b9bf633c99fd27503f3413c9c3208dc8d39a05c675ff15395972b981fa2225b2 |
| FASTQ.gz | ont_ecoli_k12_reads_fastq_gz |         25950 | 53549982 | ce3dace293ab2e5171e3c5fc95a7863715214649628568393b37a4f465a61f51 |
| FASTA    | ont_ecoli_k12_reference_fna  |             1 |  4699745 | 53bb6a51b6e92139ced1e38f74b7938781027c52200922ff03718c2237d23bb4 |
| VCF.gz   | multi_region_vcf_gz          |       2000000 |  3391626 | a052d2f7836ce92c00cb3d8033edadf8e5db2372d78a46c44c1d6c0410555ad7 |
| BCF      | multi_region_bcf             |       2000000 |  4626191 | dce85f5c236290141e47401a837946a81e39e567a20c3dbd8fbd1dc223bb3461 |

| workload | mismatches | rows  | columns | checksum                   |
|:---------|:-----------|:------|:--------|:---------------------------|
| BAM      | 0          | 27377 | 14      | 253218916404808632371051   |
| FASTQ.gz | 0          | 25950 | 4       | 239238163988283001417796   |
| FASTA    | 0          | 1     | 3       | 17575931472513828760       |
| VCF.gz   | 0          | 2e+06 | 7       | 18437751706168571091955499 |
| BCF      | 0          | 2e+06 | 7       | 18437751706168571091955499 |

## Matched reader timing

Median seconds and candidate/baseline ratios. The single-record FASTA
has millisecond-scale variation and limited timing resolution; its 14–15
ms range is not evidence of a general reader speed difference.

| workload | elapsed_baseline | elapsed_candidate | ratio |
|:---------|-----------------:|------------------:|------:|
| BAM      |            0.312 |             0.310 | 0.994 |
| BCF      |            0.564 |             0.555 | 0.984 |
| FASTA    |            0.014 |             0.015 | 1.071 |
| FASTQ.gz |            0.498 |             0.499 | 1.002 |
| VCF.gz   |            0.687 |             0.687 | 1.000 |

| workload | build     | elapsed.min | elapsed.median | elapsed.max |
|:---------|:----------|------------:|---------------:|------------:|
| BAM      | baseline  |       0.310 |          0.312 |       0.314 |
| BCF      | baseline  |       0.551 |          0.564 |       0.578 |
| FASTA    | baseline  |       0.014 |          0.014 |       0.015 |
| FASTQ.gz | baseline  |       0.494 |          0.498 |       0.499 |
| VCF.gz   | baseline  |       0.678 |          0.687 |       0.694 |
| BAM      | candidate |       0.306 |          0.310 |       0.312 |
| BCF      | candidate |       0.545 |          0.555 |       0.564 |
| FASTA    | candidate |       0.014 |          0.015 |       0.015 |
| FASTQ.gz | candidate |       0.494 |          0.499 |       0.500 |
| VCF.gz   | candidate |       0.682 |          0.687 |       0.694 |

| build     | load_elapsed |
|:----------|-------------:|
| baseline  |        0.138 |
| candidate |        0.125 |

Raw timings, source/binary identities, input hashes and equality proofs
are kept in `data/init_readers_{timings,metadata,inputs,proofs}.csv`.
The default render measures and validates the workload; `measure: false`
renders those recorded CSVs without executing readers or creating new
measurements.
