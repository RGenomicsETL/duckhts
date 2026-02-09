# DuckHTS Agent Guidelines

This document provides guidance for AI agents working on the DuckHTS DuckDB extension.

## Project Goal
Build a DuckDB 1.4+ extension that **reads** HTS file formats using htslib, with a focus on CRAN package preparation. Application code will handle format conversion (e.g., to parquet).

## Current Status
- ✅ All readers implemented: VCF/BCF, SAM/BAM/CRAM, FASTA/FASTQ, GTF/GFF, tabix
- ✅ Basic R package structure exists
- ✅ Comprehensive SQL tests implemented
- 🔄 CRAN preparation in progress

## Target DuckDB Version
- **Minimum**: DuckDB 1.4 (C/C++ extension API)
- Keep all public APIs compatible with DuckDB 1.4+.

## Documentation Conventions
- **Primary documentation**: README.Rmd (R Markdown format)
- New specs go in [./.github/SPEC.md](SPEC.md).
- The implementation plan is tracked in [./.github/PLAN.md](PLAN.md).
- Update description.yml when needed for extension metadata.
- If a changelog is introduced, list **newest entries first**.

## Agent Working Instructions
**CRITICAL**: Always read existing code and implementation before making changes:
1. Read the relevant source files to understand current implementation
2. Check existing tests to understand expected behavior
3. Review current build system before making changes
4. Preserve existing functionality and API compatibility

## Source Layout Expectations
- Extension sources live under [src/](src/).
- Public headers live under [src/include/](src/include/).
- Vendor sources live under a dedicated top-level folder (see SPEC.md).
- Test SQL lives under [test/sql/](test/sql/).
- R package structure in [r/duckhts/](r/duckhts/).

## Vendoring Rules
- Vendor scripts must be **reproducible** and pin exact versions/commits.
- All downloads must be checksum-verified.
- Scripts must support **offline rebuilds** once sources are fetched.
- Do not modify vendored code directly; patch with explicit patch files.
- Capture upstream licenses in a dedicated licenses folder and document them.

## Build & CI Rules
- Keep build scripts deterministic and non-interactive.
- Do not add network access to the extension build step.
- Network access is allowed **only** in explicit vendoring scripts.
- **CRAN requirement**: No vcpkg usage - use CMAKE and traditional autotools
- Assume cmake, make, and standard build tools available on R platforms

## R Package Specific Rules
- Package should be self-contained as much as possible
- Use CMAKE and configure/configure.win scripts for cross-platform builds
- Include cleanup and cleanup.win scripts for proper cleanup
- Simplify package bootstrapping to copy necessary extension files
- Version scheme: duckhtsVersion-x format
- All R package modifications must maintain CRAN compatibility

## Testing Guidelines
- Current tests cover: VCF/BCF, BAM/CRAM, FASTA/FASTQ, GTF/GFF, tabix readers
- Tests should cover:
  - schema correctness
  - record counts  
  - region queries (indexed files)
  - FORMAT/INFO parsing (VCF/BCF)
  - paired FASTQ handling
  - CRAM with reference files
  - attributes_map functionality for GTF/GFF
- Add more edge case tests for CRAN preparation
- Test R package loading and function calls

## Current Focus Areas
1. **CRAN Preparation**: Adapt R package build system, remove vcpkg dependency
2. **Enhanced Testing**: Add more comprehensive tests and edge cases
3. **Documentation**: Update README.Rmd and add R package vignettes
4. **Cross-platform**: Ensure builds work on Linux, macOS, Windows MinGW

## Style
- Keep changes minimal and focused.
- Preserve existing code style and APIs unless the task explicitly requires changes.
- Follow R package conventions for CRAN submission.
- Test changes with R CMD check when modifying R package components.
