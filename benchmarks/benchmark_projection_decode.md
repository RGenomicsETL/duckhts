Decode only what the query reads: VCF parsing, CRAM fields, mosdepth medians
================

This report measures three changes that remove work the query does not need:

- **VCF parsing.** `read_bcf` sets htslib’s `max_unpack` from the projected columns, so a
  text VCF record is parsed only as far as the selected columns need. A site-only query
  parses no genotypes.
- **CRAM fields.** `read_bam`, `bam_bin_counts`, `duckhts_bam_bed_coverage` and
  `read_pileup` ask the CRAM decoder only for the fields they read
  (`CRAM_OPT_REQUIRED_FIELDS`). `read_bam` now also leaves out read names, MAPQ and
  mate fields unless they are selected.
- **mosdepth medians.** With `use_median := TRUE`, each window clears and scans only
  the depth counters it used instead of all 65,536.

## Inputs

All inputs resolve through the `duckhtsbench` registry; rendering does not download.

| Input | Registry id                                                      | Shape                                                                                              |
|-------|------------------------------------------------------------------|----------------------------------------------------------------------------------------------------|
| VCF   | `ancestry_1000g_chr22`                                           | 1000 Genomes phase 3 chr22, 1,103,547 records, 2,504 samples, bgzipped text VCF with a tabix index |
| CRAM  | `insert_size_hg00403_chr20_cram_2x` with `liftover_grch38_fasta` | 30x HG00403 chr20 repeated twice, 33.6 M reads, CRAM 3 against GRCh38                              |
| BAM   | `insert_size_hg00403_chr20_bam`                                  | 30x HG00403 chr20, one populated reference                                                         |

## Method

- Engines: **baseline**, the extension built from `main` before these changes
  (`PROJECTION_DECODE_BASELINE_EXTENSION`), and **candidate**, the current build
  (`DUCKHTS_EXTENSION`, default `build/release/duckhts.duckdb_extension`).
- Every timed query runs in a fresh DuckDB CLI process (`DUCKDB`, default `duckdb`)
  with `SET threads` and a 16 GiB memory limit. Seconds are the CLI `.timer` wall time of
  the measured statement; peak RSS and CPU seconds (user plus system) are GNU `time`
  for the whole process. Three processes per engine, workload and thread count.
- Each query consumes its output through an exact fingerprint (count plus summed
  hashes), so DuckDB cannot prune the work under test.
- **Correctness before timing:** for every workload, the baseline’s and the candidate’s
  fingerprints at one thread must be identical; for mosdepth, the MD5 of every output
  file must be identical. A timed run whose result differs from that fingerprint fails
  its cell. Limits set before measuring: 900 s per process and 20 GiB peak RSS.

## Correctness

| workload        | equal |
|:----------------|:------|
| vcf_site        | TRUE  |
| vcf_info        | TRUE  |
| vcf_tidy_region | TRUE  |
| cram_count      | TRUE  |
| cram_pos        | TRUE  |
| cram_qname      | TRUE  |
| cram_bin_counts | TRUE  |
| mosdepth_median | TRUE  |

## Results

Median of three fresh processes, with the range; seconds are the measured statement
only, and peak RSS is the largest of the three. A cell with any failed, timed-out,
over-ceiling or mismatching process is marked failed and has no median or speedup.

| workload        | threads | status.baseline | seconds.baseline | range.baseline | peak_rss_mib.baseline | cpu_seconds.baseline | status.candidate | seconds.candidate | range.candidate | peak_rss_mib.candidate | cpu_seconds.candidate | speedup |
|:----------------|--------:|:----------------|-----------------:|:---------------|----------------------:|---------------------:|:-----------------|------------------:|:----------------|-----------------------:|----------------------:|--------:|
| vcf_site        |       1 | PASS            |            28.32 | 28.23-28.33    |                    43 |                 28.4 | PASS             |              4.34 | 4.33-4.36       |                     42 |                   4.4 |    6.53 |
| vcf_site        |       4 | PASS            |            28.20 | 28.19-28.24    |                    42 |                 28.2 | PASS             |              4.33 | 4.32-4.33       |                     42 |                   4.4 |    6.51 |
| vcf_info        |       1 | PASS            |            28.47 | 28.41-28.50    |                    42 |                 28.5 | PASS             |              4.87 | 4.84-4.87       |                     42 |                   4.9 |    5.85 |
| vcf_info        |       4 | PASS            |            28.42 | 28.39-28.44    |                    42 |                 28.5 | PASS             |              4.86 | 4.86-4.87       |                     42 |                   4.9 |    5.85 |
| vcf_tidy_region |       1 | PASS            |             2.54 | 2.54-2.56      |                    44 |                  2.6 | PASS             |              2.45 | 2.42-2.46       |                     45 |                   2.5 |    1.04 |
| vcf_tidy_region |       4 | PASS            |             2.57 | 2.52-2.58      |                    44 |                  5.2 | PASS             |              2.44 | 2.43-2.48       |                     44 |                   2.5 |    1.05 |
| cram_count      |       1 | PASS            |             4.23 | 4.09-4.89      |                    72 |                 10.0 | PASS             |              3.51 | 3.33-3.65       |                     76 |                   8.6 |    1.21 |
| cram_count      |       4 | PASS            |             4.58 | 4.02-5.27      |                    84 |                 10.8 | PASS             |              3.29 | 3.28-3.60       |                     81 |                   8.1 |    1.39 |
| cram_pos        |       1 | PASS            |             4.45 | 4.44-4.48      |                    78 |                 10.1 | PASS             |              1.74 | 1.64-1.94       |                     81 |                   4.6 |    2.56 |
| cram_pos        |       4 | PASS            |             4.30 | 4.00-4.52      |                    82 |                  9.7 | PASS             |              1.77 | 1.55-2.01       |                     91 |                   4.7 |    2.43 |
| cram_qname      |       1 | PASS            |             4.15 | 4.09-4.84      |                    74 |                 10.2 | PASS             |              3.34 | 3.29-3.79       |                     75 |                   8.5 |    1.24 |
| cram_qname      |       4 | PASS            |             4.37 | 4.14-4.47      |                    87 |                 10.6 | PASS             |              3.26 | 3.23-4.03       |                     83 |                   8.3 |    1.34 |
| cram_bin_counts |       1 | PASS            |             9.47 | 9.42-9.50      |                   292 |                 20.4 | PASS             |              5.54 | 5.51-5.59       |                    303 |                  12.3 |    1.71 |
| cram_bin_counts |       4 | PASS            |             9.49 | 9.38-9.51      |                   289 |                 20.4 | PASS             |              5.58 | 5.54-5.60       |                    303 |                  12.3 |    1.70 |
| mosdepth_median |       1 | PASS            |             6.17 | 6.14-6.18      |                  1921 |                  9.2 | PASS             |              5.33 | 5.29-5.36       |                   1921 |                   8.4 |    1.16 |
| mosdepth_median |       4 | PASS            |             6.14 | 6.12-6.17      |                  1922 |                  9.2 | PASS             |              5.32 | 5.32-5.39       |                   1921 |                   8.3 |    1.15 |

| key                      | value                                                            |
|:-------------------------|:-----------------------------------------------------------------|
| source_revision          | c335b7628eb06e95c32df140af896cae69c2395e                         |
| candidate_sha256         | 51df5145980674bb42143ebb568816fce91eb3fbb02b85295289fbc352ca8e8e |
| baseline_sha256          | 2cae6a33478cdd279fa7fcc6a2db882e73c83109ad14315d71baea22402ad20e |
| baseline_source_revision | e9ecdea51d44e731303fa90b73538c4b54fca637                         |
| duckdb                   | v1.5.1 (Variegata) 7dbb2e646f                                    |
| nproc                    | 20                                                               |
| cpu                      | 13th Gen Intel(R) Core(TM) i5-13500                              |
| rendered_utc             | 2026-11-03 20:41:45 UTC                                          |

`read_bam` and `bam_bin_counts` read one populated reference here, so their CRAM cells
run one reader at both thread counts. `duckhts_mosdepth` runs with `threads := 1`
(htslib decompression threads); its processing does not follow `SET threads`.
