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
- [ ] Adapt R package to use CMAKE and configure/configure.win
- [ ] Add cleanup and cleanup.win scripts
- [ ] Simplify package bootstrapping to copy needed extension files
- [ ] Ensure self-contained package (assume cmake, make available)
- [ ] Update package versioning scheme (duckhtsVersion-x)
- [ ] Remove vcpkg dependency for CRAN compatibility
- [ ] Add comprehensive R package tests
- [ ] Update DESCRIPTION with proper dependencies and SystemRequirements

## 📋 Phase 7 — Enhanced Testing & Documentation
- [ ] Add more comprehensive edge case tests
- [ ] Add performance benchmarks
- [ ] Update README.Rmd with CRAN-specific installation instructions
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
- [ ] Add `COPY` integration where feasible
- [ ] Consider based on community feedback after CRAN release

## 📝 Notes
- **Scope**: Focus on READERS ONLY - application code will handle format conversion (e.g., to parquet)
- **Target**: DuckDB 1.4+ C API compatibility maintained
- **Documentation**: Primary documentation in README.Rmd
- **Build**: Self-contained package, no vcpkg on CRAN
- **Testing**: Comprehensive conformance and edge case coverage required
