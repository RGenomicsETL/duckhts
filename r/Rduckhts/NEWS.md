
# Rduckhts 1.1.4-0.0.1

- expand `rduckhts_score()` test coverage: add FORMAT/AS integer dosage fixture, missing DS value fixture, and seven GWAS summary preset fixtures (REGENIE, SAIGE, BOLT, METAL, PGS, SSF/GWAS-SSF); add R tinytest cases for TSV counts without threshold, GWAS-VCF with `q_score_thr`, and mutual-exclusion error assertions (`include`+`exclude`, `regions`+`regions_file`, `targets`+`targets_file`)

- Bundle duckhts 1.1.4 extension.

- add `no_left_align` parameter to `rduckhts_liftover()` (default `FALSE`): when `TRUE`, skips post-liftover left-alignment of indels, mirroring `bcftools +liftover --no-left-align`
- harden `rduckhts_score()` filter-expression failures for installed-package CI runs on Linux/Windows: invalid `include`/`exclude` expressions now surface as normal package errors from the bundled score engine instead of relying on a cross-frame bcftools shim jump path
- extend `rduckhts_score()` with bcftools-style filtering arguments (`regions`, `regions_file`, `regions_overlap`, `targets`, `targets_file`, `targets_overlap`, `apply_filters`, `include`, `exclude`) and package-level tests for region/target/FILTER behavior; `include`/`exclude` expressions are now evaluated against core VCF fields (`POS`, `QUAL`, `CHROM`, `ID`, `REF`, `ALT`, `FILTER`)
- improve `rduckhts_liftover()` multiallelic semantic parity with upstream `bcftools +liftover`: preserve all ALT alleles through reference-introduction/swaps, apply dynamic-allele indel normalization/left-alignment, treat `ALT='.'` as no alternate alleles, and add tinytest coverage for forward/reverse multiallelic liftover outputs
- add GWAS-VCF multi-PRS scoring: `rduckhts_score()` now supports GWAS-VCF summary files where each sample in the VCF defines a separate PRS with `FORMAT/ES` and `FORMAT/LP` fields
- add INFO/END liftover: `rduckhts_liftover()` gains an `end_pos_col` parameter to lift INFO/END positions alongside the primary coordinate, producing a `dest_end` output column
- add mitochondria passthrough: `rduckhts_liftover()` gains a `lift_mt` parameter (default `FALSE`) that renames MT/chrM contigs with matching source/destination sizes without chain lifting, matching upstream `bcftools +liftover`
- harden `rduckhts_liftover()` wrapper validation: negative `max_snp_gap`/`max_indel_inc` now fail in R before SQL execution to keep tinytest/check-r-package behavior stable across CI platforms
- add fai-only mode: `rduckhts_munge()` `fasta_ref` parameter is now optional (default `NULL`); when omitted, alleles pass through as-is without reference matching, matching upstream `bcftools +munge --fai`-only behavior

- port full indel liftover pipeline from upstream `bcftools +liftover` (Tier 3 — Algorithmic Completeness): allele extension, Needleman-Wunsch re-alignment with BWA-MEM affine gap scoring, 3-phase `find_reference` with interior matching and SNP rescue, left-alignment post-liftover via VT normalization, and SNP retry-as-indel dispatch — `rduckhts_liftover()` now handles indels near chain gaps, difficult SNPs at chain boundaries, and padded variants with correct coordinate adjustment
- **breaking**: some previously-rejected edge-of-chain SNPs now map successfully with `note = 'Padded'` due to the upstream-compatible retry-as-indel dispatch; R test expectations updated accordingly

- harden `rduckhts_score()` preset column parity: add 45 missing column entries across all 8 presets so non-scoring columns (A2, FRQ, SE, N, INFO, etc.) are recognized and silently consumed during header matching, matching upstream `score.h` definitions
- add gzipped summary file support to `rduckhts_score()`: `.gz`/`.bgz` compressed PGS Catalog and GWAS summary files can now be passed directly without decompression
- fix `rduckhts_score()` PRS name derivation to strip multiple known extensions (`{gz,txt,tsv,vcf,bcf}`) iteratively, matching upstream behavior
- add validation for explicitly requested `use` tag in `rduckhts_score()`: error if the VCF header lacks the specified FORMAT field instead of silently scoring zero
- fix `rduckhts_score()` CNT column naming to match upstream: when `q_score_thr` and `counts` are both active, count columns are now named `<prs>_CNT_p<thr>` (upstream pattern) instead of `<prs>_p<thr>_CNT` (**breaking** column name change for q_score_thr + counts queries)
- fix `rduckhts_score()` threshold boundary precision to match upstream: parse `q_score_thr` values with `strtof`→double promotion before `-log10`, reproducing the exact float→double comparison asymmetry in upstream `bcftools +score`; markers at exact P-value boundaries may now be excluded where they were previously included
- fix `rduckhts_score()` wrong-result bugs: remove incorrect METAL `Zscore → P` and SSF `standard_error → P` column mappings, add haploid GT support for chrX male samples, and fix memory leak on skipped markers
- fix `rduckhts_liftover()` wrong-result bugs: skip reverse complement for symbolic alleles and detect insertions as indels in the liftover path
- fix `rduckhts_munge()` wrong-result bugs: emit `filter = 'MissingContig'` instead of aborting when FASTA ref fetch fails for unknown contigs, and propagate NAN correctly on AC swap with missing NS
- harden `rduckhts_score()` upstream parity: port flexible chromosome name resolution (chr prefix strip/prepend, 23→X, 24→Y, 26/MT→chrM aliases), numerically stable `-log10(p)` parsing for very small p-values, and NA/`.` missing value handling in summary stats fields (BETA, OR, P, LP)
- add comprehensive `rduckhts_score()` test coverage: DS/HDS/AP/GP dosage modes with real FORMAT values, OR-to-beta conversion, PLINK2 preset with LOG10_P, custom columns_file mapping, allele mismatch (zero-match), missing genotype (`./. `) handling, NA in summary stats, chr-prefix flexible matching, auto-detection priority, and small p-value threshold precision
- strengthen `rduckhts_score()` parity coverage for variant matching by adding rsID-based fixture/tests (`use_variant_id = TRUE`), explicit CHR/BP mismatch behavior checks, and an error assertion when marker IDs are required but missing
- add METAL meta-analysis output support: `rduckhts_munge()` auto-detects METAL/heterogeneity keys (`INFO`, `HET_I2`, `HET_P`, `HET_LP`, `DIRE`) in the resolved column map and dispatches to `duckdb_munge_metal(...)`, returning `si`, `i2`, `cq`, and `ed` columns; non-METAL presets continue to use the base 16-column schema without these fields
- align bundled munge outputs with the extension-level schema by renaming `swapped` to `alleles_swapped` with explicit allele-orientation semantics
- include munge APIs (`duckdb_munge` and `bcftools_munge_row`) in the generated bundled function catalog metadata and add a dedicated `benchmark_munge.Rmd` workflow (`make bench-munge`) for DuckHTS vs `bcftools +munge` benchmarking
- harden bundled liftover diagnostics: `rduckhts_liftover()` now surfaces explicit SQL errors for invalid `chrom`/`pos` rows through `duckdb_liftover(...)`, direct `bcftools_liftover(...)` calls also error on invalid required inputs, and the package tests cover both wrapper-level and scalar invalid-row paths
- harden bundled liftover upstream parity: full IUPAC complement table, IUPAC→N reference sequence sanitization, `25→X`/`26→MT` contig aliases, Ensembl chain ID dedup, chain coverage validation, and `uint64_t` regidx payload — all exact ports from upstream `bcftools +liftover`
- align bundled liftover outputs with `bcftools +liftover` semantics by replacing the old `warning` field with `reject_reason` for rejected rows and `note` for emitted rows that carry extra annotations
- add README liftover examples and broaden `rduckhts_liftover()` tinytest coverage for unmapped rows plus chain/FASTA and parameter validation failures
- expose the bundled `liftover(...)` table macro for score-style variant rows via a new `rduckhts_liftover()` helper that runs the macro against an input SQL query/table expression and returns lifted coordinates, alleles, reject reasons, and note annotations
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
