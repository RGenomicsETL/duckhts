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
- Extension function documentation and community-extension descriptor metadata are sourced from [functions.yaml](functions.yaml); keep it updated whenever functions/macros/UDFs are added, removed, renamed, or re-described.
- Regenerate the package-bundled function catalog and the checked-in community extension descriptor from `functions.yaml` with `python3 scripts/render_function_catalog.py`.
- Keep the generated files in `r/Rduckhts/inst/function_catalog/` and [community-extensions/extensions/duckhts/description.yml](community-extensions/extensions/duckhts/description.yml) in sync with the manifest.
- New specs go in [./.github/SPEC.md](SPEC.md).
- The implementation plan is tracked in [./.github/PLAN.md](PLAN.md).
- Do not hand-edit the generated `community-extensions/extensions/duckhts/description.yml`; update `functions.yaml` and re-render it instead.
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
- Keep `r/Rduckhts/README.Rmd` wired to the generated function catalog instead of duplicating extension function lists by hand.
- When upstream extension sources under `src/` change, update the bundled R package copy by running `Rscript bootstrap.R ~/duckhts/` from `r/Rduckhts/`, then reinstall the package before running tests.
- Never run `R CMD INSTALL .` from `r/Rduckhts/`: it mutates `inst/duckhts_extension/htslib` in place and removes the vendored htslib source tree needed for subsequent installs. Build a tarball with `R CMD build .` and install the resulting tarball instead.
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

## Best Practices Learned
- Keep Windows builds on MinGW/Rtools; match DuckDB platform strings (e.g., windows_amd64_mingw) in extension metadata.
- Use $ORIGIN rpaths on Linux so libhts.so can be found at runtime.
- Prefer pkg-config (plus Homebrew fallbacks) for OpenSSL/libcurl detection on macOS.
- Avoid plugins on Windows; keep plugin support behind configure checks on Unix.
- Bundle all runtime data under inst/extdata for R and keep README examples runnable with local files.
- Keep README examples deterministic and short; use eval=FALSE only when external network access is required.

## Heng Li Component Scope
- Treat `cgranges` as the interval engine (overlap/contain/index) for BED/GRanges-style operations.
- Treat `seqtk` ideas (`kseq`/buffered stream parsing patterns) as parser implementation guidance for FASTA/FASTQ-like read paths.
- Do not conflate `cgranges` and `seqtk`; they address different layers.
- When adding interval APIs, keep base-level BAM/CRAM pileup/depth logic in htslib-native code paths and use `cgranges` for interval algebra.

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
