# DuckHTS Implementation Plan (Draft)

## Phase 0 — Repository Setup
- [ ] Add vendoring layout under `third_party/`.
- [ ] Add initial vendoring scripts under `scripts/`.
- [ ] Add license capture workflow for upstream sources.

## Phase 1 — htslib Integration
- [ ] Add minimal build of htslib into extension build.
- [ ] Provide a small C wrapper layer for reading records.
- [ ] Add basic `bcf_read` table function.
- [ ] Implement header parsing and structured metadata capture.

## Phase 2 — SAM/BAM/CRAM Readers
- [ ] Add `sam_read`, `bam_read`, `cram_read` table functions.
- [ ] Implement indexed region filtering.
- [ ] Add tag handling strategy for optional fields.
- [ ] Preserve read group and reference metadata for round-trip conversion.

## Phase 3 — Write Path (Selective)
- [ ] Implement initial write path for VCF/BCF.
- [ ] Add `COPY` integration where feasible.

## Phase 4 — Conformance & CI
- [ ] Add conformance datasets (small, indexed).
- [ ] Add tests in `test/sql/` for core read paths.
- [ ] Add CI job to run conformance tests offline.

## Phase 5 — R Package Harness
- [ ] Add an R package subdirectory for testing.
- [ ] Provide `bootstrap.R` to vendor and build the extension from R.
- [ ] Add smoke tests that load and query the extension from R.

## Phase 6 — Performance & Polish
- [ ] Projection pushdown for FORMAT/INFO/tag columns.
- [ ] Tidy format optimization for sample-wise reads.
- [ ] Documentation and examples in README.
