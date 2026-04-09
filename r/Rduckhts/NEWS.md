
# Rduckhts 1.1.6-0.0.2 (2026-04-09)

- Fix `test_bam_file_offset`: cast `COUNT(*)` results to `INTEGER` in SQL so the DuckDB driver returns R `integer` rather than `numeric` (BIGINT maps to double in the duckdb R driver), restoring `expect_identical` assertions.

# Rduckhts 1.1.6-0.0.1 (2026-04-09)

- Fix bundled `read_hts_index_spans(...)` / `rduckhts_hts_index_spans()`: the span view now returns real chunk rows from CSI/TBI/BAI indexes, including populated `bin`, `chunk_beg_vo`, `chunk_end_vo`, `chunk_bytes`, `seq_start`, and `seq_end` values instead of placeholder `NA`s; BCF-backed calls also avoid the previous noisy `tbx` probe warning on `.csi` indexes.
- Add `FILE_OFFSET` column to `rduckhts_bam()` / `read_bam(...)`: exposes the BGZF virtual file offset after each record. Zero runtime overhead (macro over already-open struct fields). Enables `ORDER BY FILE_OFFSET` in SQL `LAG()` / `LAST_VALUE()` window functions to reproduce exact BAM file order for streaming deduplication algorithms. Together with the `//` integer-division operator and `LAST_VALUE(... IGNORE NULLS)`, this permits exact replication of WisecondorX's larp/larp2 state machine in pure SQL, confirmed at 0 mismatches across 25,115 non-zero bins on a real NIPT BAM.

# Rduckhts 1.1.5-0.0.1 (2026-04-08)

- Fix bundled `bcftools_liftover(...)` / `rduckhts_liftover()` cache and realignment hardening: per-thread chain/FASTA contexts are now bounded instead of accumulating for the lifetime of worker threads, and scalar left-alignment no longer reuses stale traceback state after failed/empty alignments.
- Fix bundled `read_bam(...)` / `rduckhts_bam()` and `read_bcf(...)` / `rduckhts_bcf()` indexed parallel full scans when headers contain leading empty contigs: contig claiming now retries iteratively instead of recursively, and the BAM reader no longer returns an empty chunk after successfully handing off to the next contig.
- Keep the top-level extension `README.Rmd` examples aligned with direct extension usage: the extension README now renders its example queries through a custom DuckDB SQL knitr engine instead of `R`/`DBI`, and its liftover example uses bundled fixtures rather than temporary `R`-generated FASTA/chain files.
- Fix bundled Windows GNU CMake builds: the vendored `htslib` configure step now distinguishes `windows_amd64_mingw` from `windows_amd64_rtools`; the MinGW path keeps the smaller `configure.win`-style library set, while the Rtools path restores the fuller static `libcurl` dependency closure required by its `htslib` feature probes. `CURL_STATICLIB` remains on the built objects rather than on `./configure` test probes.
- Fix bundled Windows `windows_amd64_rtools` CMake builds: the upstream extension `Makefile` now pins `CC`/`AR`/`RANLIB` from `R CMD config`, avoiding mixed non-Rtools compiler and Rtools library selection when vendored `htslib` is configured; the vendored `htslib` CMake path also returns to separate configure/build steps on MinGW for simpler diagnostics and behavior, and MinGW static-libcurl builds now define `CURL_STATICLIB` to match Rtools `libcurl.a`.
- Fix bundled `read_bcf(...)` / `rduckhts_bcf()` mapping of fixed-count INFO/FORMAT arrays: exact-cardinality fields such as `Number=2` and `Number=4` now materialize as DuckDB array/list columns instead of silently dropping all but the first value.
- Fix bundled `read_bcf(...)` / `rduckhts_bcf()` handling of string FORMAT lists such as DRAGEN `FORMAT/LAA`: `Number != 1` string FORMAT fields now materialize as `VARCHAR[]` instead of triggering DuckDB internal assertion failures.
- Fix bundled `duckdb_munge(...)` / `rduckhts_munge()` multithreaded FASTA lookups: FASTA index handles are now thread-local and FASTA fetches are synchronized in `munge`, avoiding intermittent `fai_retrieve` failures and aborts when `fasta_ref` is used with `PRAGMA threads > 1`.
- Add `rduckhts_score()`: polygenic risk score computation backed by the `bcftools +score` plugin, supporting GT/DS/HDS/AP/GP/AS dosage modes, all major GWAS summary presets (PLINK, PLINK2, REGENIE, SAIGE, BOLT, METAL, PGS, SSF/GWAS-SSF), GWAS-VCF multi-PRS scoring, p-value thresholding, sample subsetting, and region/filter controls.
- Add `rduckhts_munge()`: GWAS summary statistics normalization backed by `bcftools +munge`, with FASTA reference allele resolution, swap-aware effect/frequency transforms, and METAL meta-analysis column support.
- Add `rduckhts_liftover()`: variant coordinate liftover backed by `bcftools +liftover` using UCSC chain files, with full indel normalization, INFO/END lifting, and MT passthrough.
- Add `rduckhts_bed()` for BED3–BED12 interval files and `rduckhts_fasta_nuc()` for nucleotide composition over BED intervals or fixed-width bins.
- Add compression and index helpers: `rduckhts_bgzip()`, `rduckhts_bgunzip()`, `rduckhts_bam_index()`, `rduckhts_bcf_index()`, and `rduckhts_tabix_index()`.
- Add HTS metadata readers: `rduckhts_hts_header()`, `rduckhts_hts_index()`, `rduckhts_hts_index_spans()`, and `rduckhts_hts_index_raw()`.
- Add quality encoding controls to `rduckhts_bam()` and `rduckhts_fastq()` (`quality_representation`, `input_quality_encoding`) and `rduckhts_detect_quality_encoding()` for heuristic FASTQ encoding detection.
- Add `sequence_encoding := 'nt16'` parameter to `rduckhts_bam()`, `rduckhts_fasta()`, and `rduckhts_fastq()` for raw htslib nt16 sequence output as `UTINYINT[]`.
- Add SAM flag helpers `sam_flag_bits()` and `sam_flag_has()`, CIGAR utility functions, and `is_forward_aligned()`.
- Bundle duckhts 1.1.5 extension.

# Rduckhts 0.1.3-0.0.2.9000

# Rduckhts 0.1.3-0.0.2

- Conditionaly enable plugins in windows

- Updates the configure script to avoid check faillure on CRAN MacOS 

- Update the extension version to 0.1.3

# Rduckhts  0.1.2-0.1.5

- Fixed inadvertant removal of libexec
- Updated the plugin to add header table functions

# Rduckhts 0.1.2-0.1.4

- CRAN Submission

# Rduckhts 0.1.2-0.0.9000

- Different fixes for CRAN submission
    - Updated DESCRIPTION Title/Description formatting and added HTSlib reference.
    - Removed default write paths in bootstrap/build helpers; now require explicit paths.
    - setup_hts_env now accepts an explicit plugins_dir parameter.
    - duckhts_build now accepts a make argument (GNU make required).

- modified configure to attemp to support wasm
- Update bootstrapped extension code to match `duckhts` 0.1.2.
- Add SAMtags + auxiliary tag support (standard_tags, auxiliary_tags).
- Add tabix header/typing options (header, header_names, auto_detect, column_types).


# Rduckhts 0.1.1-0.0.3

- make the build single threaded

# Rduckhts 0.1.1-0.0.3

- misspeling correction

# Rduckhts 0.1.1-0.0.2

- CRAN resubmission: apply DuckDB C API header patch to avoid strict-prototypes warnings.

# Rduckhts 0.1.1-0.0.1

- CRAN Submission

- Bump bundled duckhts extension version to 0.1.1.

- Initial development release.
- Bundles the DuckHTS DuckDB extension and htslib for HTS file readers.
- Adds table-creation helpers for VCF/BCF, BAM/CRAM, FASTA/FASTQ, GFF/GTF, and tabix.
