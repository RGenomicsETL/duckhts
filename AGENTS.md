# DuckHTS Agent Guidelines

This document provides guidance for AI agents working on the DuckHTS DuckDB extension.

## Project Goal
Build a DuckDB 1.4+ extension that reads and writes HTS file formats using htslib, with optional use of bcftools and samtools for compatibility and conformance testing.

## Target DuckDB Version
- **Minimum**: DuckDB 1.4 (C/C++ extension API)
- Keep all public APIs compatible with DuckDB 1.4+.

## Documentation Conventions
- New specs go in [SPEC.md](SPEC.md).
- The implementation plan is tracked in [PLAN.md](PLAN.md).
- If a changelog is introduced, list **newest entries first**.

## Source Layout Expectations
- Extension sources live under [src/](src/).
- Public headers live under [src/include/](src/include/).
- Vendor sources live under a dedicated top-level folder (see SPEC.md).
- Test SQL lives under [test/sql/](test/sql/).

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

## Testing Guidelines
- Provide minimal conformance tests for VCF/BCF and BAM/CRAM.
- Tests should cover:
  - schema correctness
  - record counts
  - region queries (indexed files)
  - FORMAT/INFO parsing (VCF/BCF)

## Style
- Keep changes minimal and focused.
- Preserve existing code style and APIs unless the task explicitly requires changes.
