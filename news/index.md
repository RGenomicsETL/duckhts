# Changelog

## Rduckhts 1.1.5-0.0.1

- Fix bundled Windows `windows_amd64_rtools` CMake builds: the upstream
  extension `Makefile` now pins `CC`/`AR`/`RANLIB` from `R CMD config`,
  avoiding mixed non-Rtools compiler and Rtools library selection when
  vendored `htslib` is configured; the vendored `htslib` CMake path also
  returns to separate configure/build steps on MinGW for simpler
  diagnostics and behavior, and MinGW static-libcurl builds now define
  `CURL_STATICLIB` to match Rtools `libcurl.a`.
- Fix bundled `read_bcf(...)` /
  [`rduckhts_bcf()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bcf.md)
  mapping of fixed-count INFO/FORMAT arrays: exact-cardinality fields
  such as `Number=2` and `Number=4` now materialize as DuckDB array/list
  columns instead of silently dropping all but the first value.
- Fix bundled `read_bcf(...)` /
  [`rduckhts_bcf()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bcf.md)
  handling of string FORMAT lists such as DRAGEN `FORMAT/LAA`:
  `Number != 1` string FORMAT fields now materialize as `VARCHAR[]`
  instead of triggering DuckDB internal assertion failures.
- Fix bundled `duckdb_munge(...)` /
  [`rduckhts_munge()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_munge.md)
  multithreaded FASTA lookups: FASTA index handles are now thread-local
  and FASTA fetches are synchronized in `munge`, avoiding intermittent
  `fai_retrieve` failures and aborts when `fasta_ref` is used with
  `PRAGMA threads > 1`.
- Add
  [`rduckhts_score()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_score.md):
  polygenic risk score computation backed by the `bcftools +score`
  plugin, supporting GT/DS/HDS/AP/GP/AS dosage modes, all major GWAS
  summary presets (PLINK, PLINK2, REGENIE, SAIGE, BOLT, METAL, PGS,
  SSF/GWAS-SSF), GWAS-VCF multi-PRS scoring, p-value thresholding,
  sample subsetting, and region/filter controls.
- Add
  [`rduckhts_munge()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_munge.md):
  GWAS summary statistics normalization backed by `bcftools +munge`,
  with FASTA reference allele resolution, swap-aware effect/frequency
  transforms, and METAL meta-analysis column support.
- Add
  [`rduckhts_liftover()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_liftover.md):
  variant coordinate liftover backed by `bcftools +liftover` using UCSC
  chain files, with full indel normalization, INFO/END lifting, and MT
  passthrough.
- Add
  [`rduckhts_bed()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bed.md)
  for BED3–BED12 interval files and
  [`rduckhts_fasta_nuc()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_fasta_nuc.md)
  for nucleotide composition over BED intervals or fixed-width bins.
- Add compression and index helpers:
  [`rduckhts_bgzip()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bgzip.md),
  [`rduckhts_bgunzip()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bgunzip.md),
  [`rduckhts_bam_index()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bam_index.md),
  [`rduckhts_bcf_index()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bcf_index.md),
  and
  [`rduckhts_tabix_index()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_tabix_index.md).
- Add HTS metadata readers:
  [`rduckhts_hts_header()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_hts_header.md),
  [`rduckhts_hts_index()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_hts_index.md),
  [`rduckhts_hts_index_spans()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_hts_index_spans.md),
  and
  [`rduckhts_hts_index_raw()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_hts_index_raw.md).
- Add quality encoding controls to
  [`rduckhts_bam()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bam.md)
  and
  [`rduckhts_fastq()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_fastq.md)
  (`quality_representation`, `input_quality_encoding`) and
  [`rduckhts_detect_quality_encoding()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_detect_quality_encoding.md)
  for heuristic FASTQ encoding detection.
- Add `sequence_encoding := 'nt16'` parameter to
  [`rduckhts_bam()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bam.md),
  [`rduckhts_fasta()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_fasta.md),
  and
  [`rduckhts_fastq()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_fastq.md)
  for raw htslib nt16 sequence output as `UTINYINT[]`.
- Add SAM flag helpers `sam_flag_bits()` and `sam_flag_has()`, CIGAR
  utility functions, and `is_forward_aligned()`.
- Bundle duckhts 1.1.5 extension.

## Rduckhts 0.1.3-0.0.2.9000

## Rduckhts 0.1.3-0.0.2

CRAN release: 2026-02-24

- Conditionaly enable plugins in windows

- Updates the configure script to avoid check faillure on CRAN MacOS

- Update the extension version to 0.1.3

## Rduckhts 0.1.2-0.1.5

- Fixed inadvertant removal of libexec
- Updated the plugin to add header table functions

## Rduckhts 0.1.2-0.1.4

CRAN release: 2026-02-23

- CRAN Submission

## Rduckhts 0.1.2-0.0.9000

- Different fixes for CRAN submission
  - Updated DESCRIPTION Title/Description formatting and added HTSlib
    reference.
  - Removed default write paths in bootstrap/build helpers; now require
    explicit paths.
  - setup_hts_env now accepts an explicit plugins_dir parameter.
  - duckhts_build now accepts a make argument (GNU make required).
- modified configure to attemp to support wasm
- Update bootstrapped extension code to match `duckhts` 0.1.2.
- Add SAMtags + auxiliary tag support (standard_tags, auxiliary_tags).
- Add tabix header/typing options (header, header_names, auto_detect,
  column_types).

## Rduckhts 0.1.1-0.0.3

- make the build single threaded

## Rduckhts 0.1.1-0.0.3

- misspeling correction

## Rduckhts 0.1.1-0.0.2

- CRAN resubmission: apply DuckDB C API header patch to avoid
  strict-prototypes warnings.

## Rduckhts 0.1.1-0.0.1

- CRAN Submission

- Bump bundled duckhts extension version to 0.1.1.

- Initial development release.

- Bundles the DuckHTS DuckDB extension and htslib for HTS file readers.

- Adds table-creation helpers for VCF/BCF, BAM/CRAM, FASTA/FASTQ,
  GFF/GTF, and tabix.
