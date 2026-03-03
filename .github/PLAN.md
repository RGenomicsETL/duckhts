# DuckHTS Implementation Plan

## ✅ Phase 0 — Repository Setup (COMPLETED)
- [x] Add vendoring layout under `third_party/`.
- [x] Add initial vendoring scripts under `scripts/`.
- [x] Add license capture workflow for upstream sources.

## ✅ Phase 1 — htslib Integration (COMPLETED)
- [x] Add minimal build of htslib into extension build.
- [x] Provide a small C wrapper layer for reading records.
- [x] Add basic `read_bcf` table function.
- [x] Implement header parsing and structured metadata capture.

## ✅ Phase 2 — SAM/BAM/CRAM Readers (COMPLETED)
- [x] Add `read_bam` table function (covers SAM/BAM/CRAM).
- [x] Implement indexed region filtering.
- [x] Add tag handling strategy for optional fields.
- [x] Preserve read group and reference metadata for round-trip conversion.

## ✅ Phase 3 — Additional Readers (COMPLETED)
- [x] Add `read_fasta` and `read_fastq` table functions.
- [x] Add `read_gff` and `read_gtf` table functions.
- [x] Add `read_tabix` generic tabix reader.
- [x] Implement VEP annotation parser.

## ✅ Phase 4 — Conformance & CI (COMPLETED)
- [x] Add conformance datasets (small, indexed).
- [x] Add tests in `test/sql/` for core read paths.
- [x] Add CI job to run conformance tests offline.

## ✅ Phase 5 — R Package Harness (COMPLETED)
- [x] Add an R package subdirectory for testing.
- [x] Provide `bootstrap.R` to vendor and build the extension from R.
- [x] Add smoke tests that load and query the extension from R.

## 🔄 Phase 6 — CRAN Preparation (IN PROGRESS)
- [x] Basic R package structure with DESCRIPTION, NAMESPACE
- [x] Bootstrap script for building extension
- [x] Vendored htslib in R package structure
- [x] Adapt R package to use CMAKE and configure/configure.win
- [x] Add cleanup and cleanup.win scripts
- [x] Simplify package bootstrapping to copy needed extension files
- [x] Ensure self-contained package (assume cmake, make available)
- [x] Update package versioning scheme (duckhtsVersion-x)
- [x] Remove vcpkg dependency for CRAN compatibility
- [x] Add comprehensive R package tests
- [x] Update DESCRIPTION with proper dependencies and SystemRequirements
- [x] CRAN submission completed at 0.1.1-0.0.4 (extension baseline)

## 📋 Phase 7 — Enhanced Testing & Documentation
- [ ] Add more comprehensive edge case tests
- [ ] Add performance benchmarks
- [ ] Update README.Rmd with CRAN-specific installation instructions
- [ ] Add DuckDB COPY-to-Parquet usage example
- [ ] Add vignettes for common use cases
- [ ] Test across multiple platforms (Linux, macOS, Windows MinGW)

## 🚀 Phase 8 — CRAN Submission Preparation
- [ ] Run R CMD check --as-cran on all platforms
- [ ] Fix any NOTEs, WARNINGs, or ERRORs
- [ ] Prepare submission materials
- [ ] Address CRAN policy compliance
- [ ] Submit to CRAN (if desired) or prepare for community-extensions

## ⏸️ Phase 9 — Write Path (FUTURE)
- [ ] Implement write path for VCF/BCF (if needed)
- [ ] Consider based on community feedback after CRAN release

## 🧬 Phase 10 — Genomic Coverage & Ranges APIs (PROPOSED)

### 10.1 Design Principles
- [ ] Keep base-level read semantics in BAM/CRAM-native table functions (htslib-backed), not post-load SQL only.
- [ ] Keep interval analytics composable in DuckDB SQL (joins/windows/aggregations on loaded tables).
- [ ] Add C genomic ranges bindings for interval algebra primitives, not as a replacement for pileup/depth engines.
- [ ] Keep responsibilities explicit across Heng Li components:
  - `cgranges` for interval indexing/overlap/contain queries (BED/GRanges operations).
  - `seqtk` parsing style (`kseq`/buffered streaming patterns) for FASTA/FASTQ streaming and low-allocation readers.
  - Do not conflate `seqtk` with `cgranges`; they solve different layers.

### 10.2 New Table Functions (Draft)
- [ ] `bam_depth(path, [region, region_file, region_mode, region_merge, region_combine, min_mapq, min_baseq, overlap, max_depth, excl_flags, incl_flags, index_path, reference])`
- [ ] `bam_coverage(path, [region, region_file, region_mode, region_merge, region_combine, min_mapq, min_baseq, overlap, max_depth, bins, excl_flags, incl_flags, index_path, reference])`
- [ ] `bam_pileup(path, [region, region_file, region_mode, region_merge, region_combine, min_mapq, min_baseq, overlap, include_deletions, include_insertions, max_depth, excl_flags, incl_flags, index_path, reference])`

### 10.3 Region Input Contract
- [ ] Support both `region := 'chr:start-end,...'` and `region_file := 'targets.bed[.gz|.bgz]'`.
- [ ] Add `region_mode := 'union' | 'per_interval'`:
  - `union`: merge overlapping target intervals before counting.
  - `per_interval`: keep intervals separate and return per-target metrics.
- [ ] Add `region_merge := true|false` for explicit overlap merge control.
- [ ] Add `region_combine := 'intersect' | 'union' | 'error'` when both `region` and `region_file` are supplied.
- [ ] Normalize BED coordinates from 0-based half-open to 1-based closed internally.
- [ ] Carry target labels (`name`/`region_id`) in `per_interval` mode.

### 10.4 Overlap Semantics (Paired-End)
- [ ] Make overlap policy explicit in all depth/coverage/pileup APIs:
  - `overlap := 'none'` (read-level counting; overlap can be double counted)
  - `overlap := 'fragment'` (count fragment once in mate-overlap span)
  - `overlap := 'samtools'` (mpileup-style overlap handling)
- [ ] Document defaults per function and samtools parity rationale.
- [ ] Avoid hidden behavior changes: overlap removal must be opt-in/out via a clear argument.

### 10.5 Coverage Output Shape (GRanges-Compatible)
- [ ] Return GRanges-like columns directly from DuckDB:
  - `seqnames`, `start`, `end`, `width`, `strand`
  - `region_id`, `region_name` (nullable; populated in `per_interval`)
  - `numreads`, `covbases`, `coverage_pct`, `mean_depth`, `mean_baseq`, `mean_mapq`
- [ ] Keep naming stable and R-friendly; add R wrapper option `as_granges = TRUE` for conversion to `GenomicRanges::GRanges`.

### 10.6 Base Reader Enhancements (Proposed)
- [ ] Add `region_file` support to base readers where indexed region pushdown exists:
  - `read_bam`, `read_bcf`, `read_fasta`, `read_gff`, `read_gtf`, `read_tabix`
- [ ] Preserve existing `region` behavior for backward compatibility.
- [ ] Document duplicate semantics clearly for overlapping region lists/files.

### 10.7 Performance & Memory Considerations
- [ ] Require/encourage index-backed scanning for region file workflows.
- [ ] Pre-sort and optionally coalesce many BED intervals by contig to reduce iterator churn.
- [ ] Maintain streaming behavior (bounded memory), including mate-pair overlap bookkeeping.
- [ ] Benchmark BAM-native depth/coverage vs equivalent SQL-on-materialized reads.
- [ ] For parser-heavy code paths, prefer seqtk-like buffered/chunked line parsing over byte-wise reads.

### 10.8 Validation & Tests
- [ ] Add synthetic tests for paired-end overlap cases:
  - overlapping mates, non-overlapping mates, discordant pairs, cross-contig mates.
- [ ] Add BED region tests:
  - union vs per-interval mode, overlapping targets, duplicate targets, mixed contig ordering.
- [ ] Add coordinate conversion tests (BED 0-based half-open edge cases).
- [ ] Add regression tests for flags, mapping/base quality thresholds, and max-depth clipping.
- [ ] Add parity checks versus selected samtools behaviors and document intentional deviations.

### 10.9 Documentation & Rollout
- [ ] Update `README.Rmd` with new function signatures and deterministic examples using local test data.
- [ ] Add R package examples returning GRanges-compatible outputs.
- [ ] Mark APIs as experimental initially; freeze signatures after one release cycle of feedback.

## 📝 Notes
- **Scope**: Focus on READERS ONLY - application code will handle format conversion (e.g., to parquet)
- **Target**: DuckDB 1.4+ C API compatibility maintained
- **Documentation**: Primary documentation in README.Rmd
- **Build**: Self-contained package, no vcpkg on CRAN
- **Testing**: Comprehensive conformance and edge case coverage required
- **Versioning**: Extension bumps use .9000 dev suffix; R package tracks release as 0.1.x-0.0.y (e.g., next 0.1.2-0.0.1 after extension news release)

## 🔍 Review Feedback (2026-02-10)
- bcf_reader: region lookup currently errors for both “contig not found” and “no overlapping records” (TODO in `src/bcf_reader.c`); consider distinguishing to avoid false failures on empty regions.
- seq_reader: paired FASTQ path assumes reads are in lockstep but doesn’t validate QNAME pairing; mismatched mates will silently pair (recommend add name check + test).
- seq_reader: interleaved mode toggles mate 1/2 regardless of QNAME suffix or odd record count; consider handling trailing unpaired read and/or validating suffixes.
- tabix_reader (generic): column count inferred from first non-# line only; files with variable columns or non-# meta-char may mis-bind schema (consider using tabix conf/meta-char when indexed, and add tests for varying columns).
- behavior consistency: bcf_reader errors on “region not found” while tabix_reader returns empty; decide on a consistent contract and document it in README/tests.
