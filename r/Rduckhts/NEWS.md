
# Rduckhts 1.1.5-0.0.1

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
