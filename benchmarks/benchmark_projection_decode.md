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

| Input | Registry id                                                                   | Shape                                                                                              |
|-------|-------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------|
| VCF   | `ancestry_1000g_chr22`                                                        | 1000 Genomes phase 3 chr22, 1,103,547 records, 2,504 samples, bgzipped text VCF with a tabix index |
| CRAM  | `insert_size_hg00403_chr20_cram_2x`, `_cram_4x`, with `liftover_grch38_fasta` | 30x HG00403 chr20 repeated twice and four times, CRAM 3 against GRCh38                             |
| BAM   | `insert_size_hg00403_chr20_bam`, `_bam_2x`, `_bam_4x`                         | 30x HG00403 chr20 and its 2x and 4x repeats, one populated reference                               |

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
- **Scale ladders.** Each changed path is measured at 1×, 2× and 4× of its growth
  dimension with everything else fixed: VCF site scans and CRAM scans by records, through
  region queries over a quarter, half and all of chr22 or chr20 (on the 4x CRAM); mosdepth
  by depth, through the 1x/2x/4x BAMs. `bam_bin_counts` scans whole files and the registry
  has only the 2x and 4x CRAMs, so it has two points. The full-file VCF site and INFO scans
  and a tidy FORMAT region are kept as single points. The time exponent per doubling,
  `log2(t(2s) / t(s))`, is reported for the candidate; above 1.25 it would need
  investigation for these linear scans.
- **Memory budget, set before measuring:** 512 MiB for VCF scans (one record and its
  decoded fields live at a time), 1 GiB for CRAM scans and bin counts (decoder slices,
  reference and bins), and 3 GiB for mosdepth (a whole-contig int32 array per processing
  worker, 2 × 258 MB for chr20, plus output buffers and the median histogram). A process
  above its budget fails its cell.
- **Correctness before timing:** for every workload, the baseline’s and the candidate’s
  fingerprints at one thread must be identical; for mosdepth, the MD5 of every output
  file must be identical. A timed run whose result differs from that fingerprint fails
  its cell. Each process is limited to 900 s.

## Correctness

| key                    | records                                                                                                                                                                                                                                                                                                          | equal |
|:-----------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------|
| vcf_site 1x            | 249909                                                                                                                                                                                                                                                                                                           | TRUE  |
| vcf_site 2x            | 510323                                                                                                                                                                                                                                                                                                           | TRUE  |
| vcf_site 4x            | 1103547                                                                                                                                                                                                                                                                                                          | TRUE  |
| vcf_site_file full     | 1103547                                                                                                                                                                                                                                                                                                          | TRUE  |
| vcf_info full          | 1103547                                                                                                                                                                                                                                                                                                          | TRUE  |
| vcf_tidy_region 0.5 Mb | 29559720                                                                                                                                                                                                                                                                                                         | TRUE  |
| cram_count 1x          | 15435240                                                                                                                                                                                                                                                                                                         | TRUE  |
| cram_count 2x          | 34591792                                                                                                                                                                                                                                                                                                         | TRUE  |
| cram_count 4x          | 67265720                                                                                                                                                                                                                                                                                                         | TRUE  |
| cram_pos 1x            | 15435240                                                                                                                                                                                                                                                                                                         | TRUE  |
| cram_pos 2x            | 34591792                                                                                                                                                                                                                                                                                                         | TRUE  |
| cram_pos 4x            | 67265720                                                                                                                                                                                                                                                                                                         | TRUE  |
| cram_qname 1x          | 15435240                                                                                                                                                                                                                                                                                                         | TRUE  |
| cram_qname 2x          | 34591792                                                                                                                                                                                                                                                                                                         | TRUE  |
| cram_qname 4x          | 67265720                                                                                                                                                                                                                                                                                                         | TRUE  |
| cram_bin_counts 2x     | 35263                                                                                                                                                                                                                                                                                                            | TRUE  |
| cram_bin_counts 4x     | 35263                                                                                                                                                                                                                                                                                                            | TRUE  |
| mosdepth_median 1x     | PREFIX.mosdepth.global.dist.txt=92927c72ba1ff9d31214b60fc3fd9b18;PREFIX.mosdepth.region.dist.txt=0a84a21fbd883f521456c216cf174f72;PREFIX.mosdepth.summary.txt=d38de1afe1c1e4699b435076bb69038c;PREFIX.regions.bed.gz=f496d9cae30d310e813ed3fa70d5dbb1;PREFIX.regions.bed.gz.csi=7ceb1f2d15ffe4c2fe66af6b7b1761f4 | TRUE  |
| mosdepth_median 2x     | PREFIX.mosdepth.global.dist.txt=4f46958c73a26d3bcc933d18de6829c4;PREFIX.mosdepth.region.dist.txt=54bdf403a1531b6efc8d9e9c3c182958;PREFIX.mosdepth.summary.txt=8d8b4710e4dc8d34fe529aedc74a369d;PREFIX.regions.bed.gz=f55f5ab051eed566deb466359002326b;PREFIX.regions.bed.gz.csi=aa5350ed216b65daa99dd12860574703 | TRUE  |
| mosdepth_median 4x     | PREFIX.mosdepth.global.dist.txt=6f1642fc99d32329934a861cc698e373;PREFIX.mosdepth.region.dist.txt=445979615115fb976fa3c0d12721cf07;PREFIX.mosdepth.summary.txt=08a065dd992fe07abc972c43789c10c7;PREFIX.regions.bed.gz=371b297633df68b9cac6a4b8e0b2472e;PREFIX.regions.bed.gz.csi=7fb81b9b3e0794dbde0c9626f047ca82 | TRUE  |

## Results

Median of three fresh processes, with the range; seconds are the measured statement
only, and peak RSS is the largest of the three. A cell with any failed, timed-out,
over-ceiling or mismatching process is marked failed and has no median or speedup.

| workload        | scale  | threads | status.baseline | seconds.baseline | range.baseline | peak_rss_mib.baseline | cpu_seconds.baseline | status.candidate | seconds.candidate | range.candidate | peak_rss_mib.candidate | cpu_seconds.candidate | speedup |
|:----------------|:-------|--------:|:----------------|-----------------:|:---------------|----------------------:|---------------------:|:-----------------|------------------:|:----------------|-----------------------:|----------------------:|--------:|
| vcf_site        | 1x     |       1 | PASS            |             7.99 | 7.99-8.00      |                    43 |                  8.1 | PASS             |              1.03 | 1.01-1.03       |                     43 |                   1.1 |    7.76 |
| vcf_site        | 1x     |       4 | PASS            |             7.90 | 7.89-8.06      |                    43 |                  8.0 | PASS             |              1.01 | 1.01-1.03       |                     43 |                   1.1 |    7.82 |
| vcf_site        | 2x     |       1 | PASS            |            16.16 | 16.15-16.18    |                    43 |                 16.2 | PASS             |              2.06 | 2.06-2.06       |                     43 |                   2.1 |    7.84 |
| vcf_site        | 2x     |       4 | PASS            |            16.14 | 16.13-16.16    |                    43 |                 16.2 | PASS             |              2.07 | 2.06-2.10       |                     44 |                   2.1 |    7.80 |
| vcf_site        | 4x     |       1 | PASS            |            34.95 | 34.93-34.96    |                    43 |                 35.0 | PASS             |              4.48 | 4.45-4.49       |                     42 |                   4.5 |    7.80 |
| vcf_site        | 4x     |       4 | PASS            |            34.87 | 34.87-34.92    |                    43 |                 35.0 | PASS             |              4.49 | 4.46-4.56       |                     44 |                   4.6 |    7.77 |
| vcf_site_file   | full   |       1 | PASS            |            28.24 | 28.22-28.27    |                    43 |                 28.3 | PASS             |              4.34 | 4.32-4.36       |                     43 |                   4.4 |    6.51 |
| vcf_site_file   | full   |       4 | PASS            |            28.18 | 28.18-28.21    |                    43 |                 28.3 | PASS             |              4.33 | 4.33-4.35       |                     43 |                   4.4 |    6.51 |
| vcf_info        | full   |       1 | PASS            |            28.43 | 28.42-28.50    |                    43 |                 28.5 | PASS             |              4.86 | 4.86-4.86       |                     43 |                   4.9 |    5.85 |
| vcf_info        | full   |       4 | PASS            |            28.92 | 28.44-28.96    |                    44 |                 57.9 | PASS             |              4.90 | 4.85-4.95       |                     44 |                   5.0 |    5.90 |
| vcf_tidy_region | 0.5 Mb |       1 | PASS            |             2.51 | 2.50-2.52      |                    45 |                  2.6 | PASS             |              2.41 | 2.37-2.46       |                     45 |                   2.5 |    1.04 |
| vcf_tidy_region | 0.5 Mb |       4 | PASS            |             2.57 | 2.56-2.57      |                    44 |                  5.2 | PASS             |              2.41 | 2.39-2.44       |                     44 |                   4.9 |    1.07 |
| cram_count      | 1x     |       1 | PASS            |             1.49 | 1.48-1.63      |                    73 |                  3.7 | PASS             |              1.46 | 1.44-1.47       |                     74 |                   3.7 |    1.02 |
| cram_count      | 1x     |       4 | PASS            |             1.72 | 1.50-1.77      |                    78 |                  6.2 | PASS             |              1.44 | 1.42-1.50       |                     72 |                   3.6 |    1.19 |
| cram_count      | 2x     |       1 | PASS            |             3.62 | 3.54-3.70      |                    76 |                  8.8 | PASS             |              3.51 | 3.47-3.57       |                     76 |                   8.6 |    1.03 |
| cram_count      | 2x     |       4 | PASS            |             4.01 | 3.66-4.43      |                    83 |                 14.3 | PASS             |              3.54 | 3.47-3.90       |                     83 |                   8.7 |    1.13 |
| cram_count      | 4x     |       1 | PASS            |             6.80 | 6.56-6.99      |                    77 |                 16.6 | PASS             |              6.64 | 6.51-6.99       |                     75 |                  16.4 |    1.02 |
| cram_count      | 4x     |       4 | PASS            |             7.93 | 6.62-8.28      |                    78 |                 28.2 | PASS             |              6.87 | 6.83-6.88       |                     75 |                  16.7 |    1.15 |
| cram_pos        | 1x     |       1 | PASS            |             1.49 | 1.45-1.52      |                    72 |                  3.5 | PASS             |              0.78 | 0.76-0.80       |                     79 |                   2.1 |    1.91 |
| cram_pos        | 1x     |       4 | PASS            |             1.67 | 1.44-1.78      |                    76 |                  4.1 | PASS             |              0.77 | 0.72-0.79       |                     78 |                   2.1 |    2.17 |
| cram_pos        | 2x     |       1 | PASS            |             3.60 | 3.47-3.67      |                    74 |                  8.4 | PASS             |              1.68 | 1.55-1.83       |                     84 |                   4.5 |    2.14 |
| cram_pos        | 2x     |       4 | PASS            |             3.93 | 3.68-4.03      |                    80 |                 13.5 | PASS             |              1.60 | 1.55-2.00       |                     86 |                   4.3 |    2.46 |
| cram_pos        | 4x     |       1 | PASS            |             6.58 | 6.55-6.65      |                    75 |                 15.3 | PASS             |              3.41 | 3.19-3.59       |                     84 |                   9.0 |    1.93 |
| cram_pos        | 4x     |       4 | PASS            |             6.86 | 6.82-6.87      |                    74 |                 15.8 | PASS             |              3.26 | 3.22-3.36       |                     85 |                   8.7 |    2.10 |
| cram_qname      | 1x     |       1 | PASS            |             1.47 | 1.46-1.54      |                    75 |                  3.8 | PASS             |              1.49 | 1.45-1.61       |                     74 |                   3.8 |    0.99 |
| cram_qname      | 1x     |       4 | PASS            |             1.87 | 1.45-1.89      |                    78 |                  6.8 | PASS             |              1.82 | 1.41-1.92       |                     77 |                   6.7 |    1.03 |
| cram_qname      | 2x     |       1 | PASS            |             3.82 | 3.49-4.16      |                    75 |                  9.5 | PASS             |              3.62 | 3.48-3.73       |                     75 |                   9.1 |    1.06 |
| cram_qname      | 2x     |       4 | PASS            |             3.84 | 3.52-4.24      |                    77 |                  9.6 | PASS             |              4.31 | 3.55-4.37       |                     80 |                  15.7 |    0.89 |
| cram_qname      | 4x     |       1 | PASS            |             7.24 | 6.65-7.51      |                    76 |                 18.1 | PASS             |              6.61 | 6.46-6.67       |                     80 |                  16.7 |    1.10 |
| cram_qname      | 4x     |       4 | PASS            |             6.91 | 6.74-9.20      |                    78 |                 17.3 | PASS             |              6.57 | 6.50-7.15       |                     78 |                  16.6 |    1.05 |
| cram_bin_counts | 2x     |       1 | PASS            |             9.41 | 9.38-9.51      |                   275 |                 20.1 | PASS             |              5.56 | 5.52-5.64       |                    299 |                  12.3 |    1.69 |
| cram_bin_counts | 2x     |       4 | PASS            |             9.45 | 9.35-9.50      |                   282 |                 20.2 | PASS             |              5.74 | 5.59-5.75       |                    307 |                  18.9 |    1.65 |
| cram_bin_counts | 4x     |       1 | PASS            |            15.78 | 15.67-16.29    |                   307 |                 34.3 | PASS             |             10.10 | 10.06-10.16     |                    347 |                  22.5 |    1.56 |
| cram_bin_counts | 4x     |       4 | PASS            |            15.68 | 15.66-15.96    |                   305 |                 34.2 | PASS             |             10.19 | 10.17-10.64     |                    344 |                  22.6 |    1.54 |
| mosdepth_median | 1x     |       1 | PASS            |             6.07 | 6.06-6.13      |                  1922 |                  9.1 | PASS             |              5.31 | 5.29-5.32       |                   1923 |                   8.2 |    1.14 |
| mosdepth_median | 1x     |       4 | PASS            |             6.12 | 6.05-6.15      |                  1922 |                  9.1 | PASS             |              5.30 | 5.29-5.36       |                   1922 |                   8.2 |    1.15 |
| mosdepth_median | 2x     |       1 | PASS            |             7.64 | 7.58-7.71      |                  1922 |                 12.2 | PASS             |              6.86 | 6.86-6.89       |                   1922 |                  11.4 |    1.11 |
| mosdepth_median | 2x     |       4 | PASS            |             7.68 | 7.58-7.73      |                  1922 |                 12.3 | PASS             |              6.83 | 6.81-6.97       |                   1922 |                  11.5 |    1.12 |
| mosdepth_median | 4x     |       1 | PASS            |            10.47 | 10.41-10.53    |                  1922 |                 18.0 | PASS             |              9.72 | 9.70-9.74       |                   1923 |                  17.3 |    1.08 |
| mosdepth_median | 4x     |       4 | PASS            |            10.48 | 10.35-10.52    |                  1922 |                 17.9 | PASS             |              9.61 | 9.60-9.63       |                   1922 |                  17.0 |    1.09 |

| key                      | value                                                            |
|:-------------------------|:-----------------------------------------------------------------|
| source_revision          | 508026e3599655ece9aea93279b458529443dd10                         |
| candidate_sha256         | 4e193c5415a3e610099bff40c1a01687564219f0db8f91791c2ed3d8fde743a6 |
| baseline_sha256          | 2cae6a33478cdd279fa7fcc6a2db882e73c83109ad14315d71baea22402ad20e |
| baseline_source_revision | e9ecdea51d44e731303fa90b73538c4b54fca637                         |
| duckdb                   | v1.5.1 (Variegata) 7dbb2e646f                                    |
| nproc                    | 6                                                                |
| cpu                      | 13th Gen Intel(R) Core(TM) i5-13500                              |
| rendered_utc             | 2026-11-03 22:57:59 UTC                                          |

## Scaling

Candidate median seconds per scale. Region ladders do not double the records exactly
(the quarter, half and whole chromosome are cut by base pairs), so the exponent is taken
against the records each scan returned, `log(t2 / t1) / log(n2 / n1)`, from the
correctness fingerprints (reads counted, for bin counts). mosdepth’s ladder doubles depth through the 1x/2x/4x BAMs, so
its exponent uses n2 / n1 = 2. A ladder with a failed cell has no exponent.

| workload        | threads | scales   | records                              | seconds            | peak_rss_mib       | exponents  |
|:----------------|--------:|:---------|:-------------------------------------|:-------------------|:-------------------|:-----------|
| vcf_site        |       1 | 1x/2x/4x | 249,909 / 510,323 / 1,103,547        | 1.03 / 2.06 / 4.48 | 43 / 43 / 42       | 0.97, 1.01 |
| vcf_site        |       4 | 1x/2x/4x | 249,909 / 510,323 / 1,103,547        | 1.01 / 2.07 / 4.49 | 43 / 44 / 44       | 1.01, 1    |
| cram_count      |       1 | 1x/2x/4x | 15,435,240 / 34,591,792 / 67,265,720 | 1.46 / 3.51 / 6.64 | 74 / 76 / 75       | 1.09, 0.96 |
| cram_count      |       4 | 1x/2x/4x | 15,435,240 / 34,591,792 / 67,265,720 | 1.44 / 3.54 / 6.87 | 72 / 83 / 75       | 1.11, 1    |
| cram_pos        |       1 | 1x/2x/4x | 15,435,240 / 34,591,792 / 67,265,720 | 0.78 / 1.68 / 3.41 | 79 / 84 / 84       | 0.95, 1.06 |
| cram_pos        |       4 | 1x/2x/4x | 15,435,240 / 34,591,792 / 67,265,720 | 0.77 / 1.6 / 3.26  | 78 / 86 / 85       | 0.91, 1.07 |
| cram_qname      |       1 | 1x/2x/4x | 15,435,240 / 34,591,792 / 67,265,720 | 1.49 / 3.62 / 6.61 | 74 / 75 / 80       | 1.1, 0.91  |
| cram_qname      |       4 | 1x/2x/4x | 15,435,240 / 34,591,792 / 67,265,720 | 1.82 / 4.31 / 6.57 | 77 / 80 / 78       | 1.07, 0.63 |
| cram_bin_counts |       1 | 2x/4x    | 33,632,860 / 67,265,720              | 5.56 / 10.1        | 299 / 347          | 0.86       |
| cram_bin_counts |       4 | 2x/4x    | 33,632,860 / 67,265,720              | 5.74 / 10.19       | 307 / 344          | 0.83       |
| mosdepth_median |       1 | 1x/2x/4x | depth 1x / 2x / 4x                   | 5.31 / 6.86 / 9.72 | 1923 / 1922 / 1923 | 0.37, 0.5  |
| mosdepth_median |       4 | 1x/2x/4x | depth 1x / 2x / 4x                   | 5.3 / 6.83 / 9.61  | 1922 / 1922 / 1922 | 0.37, 0.49 |

`read_bam` and `bam_bin_counts` read one populated reference here, so their CRAM cells
run one reader at both thread counts. `duckhts_mosdepth` runs with `threads := 1`
(htslib decompression threads); its processing does not follow `SET threads`.
