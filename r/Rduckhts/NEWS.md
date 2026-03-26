
# Rduckhts 0.1.3-0.0.2.9001

- harden bundled liftover diagnostics: `rduckhts_liftover()` now surfaces explicit SQL errors for invalid `chrom`/`pos` rows through `duckdb_liftover(...)`, and the package tests cover the wrapper-level invalid-row error paths
- add README liftover examples and broaden `rduckhts_liftover()` tinytest coverage for unmapped rows plus chain/FASTA and parameter validation failures
- expose the bundled `liftover(...)` table macro for score-style variant rows via a new `rduckhts_liftover()` helper that runs the macro against an input SQL query/table expression and returns lifted coordinates, alleles, and warnings
- Keep the generated community extension metadata in sync with the bundled extension version by sourcing the emitted top-level `version` field from the repo-level `description.yml`.
- Bundle the `duckhts` `0.1.3.9001` extension update.
- Add `quality_representation` to `rduckhts_bam()` and `rduckhts_fastq()` so qualities can be returned as raw `UTINYINT[]` Phred values.
- Add `input_quality_encoding` to `rduckhts_fastq()`, defaulting to modern `phred33` while allowing explicit legacy FASTQ decoding; expose `rduckhts_detect_quality_encoding()` for heuristic FASTQ encoding inspection.
- Add `sequence_encoding` parameter to `rduckhts_bam()`, `rduckhts_fasta()`, and `rduckhts_fastq()` wrappers, forwarding `sequence_encoding := 'nt16'` to the extension readers for raw nt16-encoded sequence output as `UTINYINT[]`.
- **Breaking**: `seq_encode_4bit()` no longer returns NULL for unknown characters — they now map to N (code 15), matching htslib `seq_nt16_table` behavior; UDF and reader nt16 paths use shared code.
- Remove the unsupported `attributes_map` argument from `rduckhts_tabix()` so the wrapper matches the generic `read_tabix(...)` surface; attribute maps remain available via `rduckhts_gff()` and `rduckhts_gtf()`.
- Expose the bundled `read_bcf(...)` CSQ/ANN/BCSQ typing cleanup, including centralized builtin rules and the `additional_csq_column_types := ...` override parameter.
- Add `read_bed(...)` and `fasta_nuc(...)` to the bundled extension surface, plus `rduckhts_bed()` and `rduckhts_fasta_nuc()` wrappers.
- Add `rduckhts_bgzip()`, `rduckhts_bgunzip()`, `rduckhts_bam_index()`, `rduckhts_bcf_index()`, and `rduckhts_tabix_index()` wrappers for the new extension compression and indexing functions.
- Expose the newer bundled extension surface in the package catalog, including HTS metadata readers, additional sequence helpers, `sam_flag_bits()`/`sam_flag_has()`, the new `CIGAR utils` helpers, and the expanded SAM/tag and tabix reader capabilities.
- Rename the ambiguous SAM flag helper names in the bundled extension/catalog to the clearer forms `is_paired()`, `is_proper_pair()`, `is_next_segment_unmapped()`, and `is_next_segment_reverse_complemented()`.
- Add `is_forward_aligned()` for mapped-strand checks in pure SQL workflows such as strand-split bin counting.
- Bootstrap the new extension sources into the package build and update `configure`/`configure.win` so the bundled extension compiles them on Unix and Windows.
- Regenerate the package-bundled function catalog and roxygen documentation for the new wrappers.
- Add installed-package tinytest coverage for BED reading, FASTA nucleotide composition, BGZF round-trips, tabix indexing, and BAM/BCF index creation.

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
