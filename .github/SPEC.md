# DuckHTS Extension Specification (Draft)

**Status**: Draft

## Overview
DuckHTS is a DuckDB 1.4+ extension that provides high-performance reading and writing of HTS file formats using htslib, with optional integration of bcftools and samtools for tooling parity and conformance testing.

The longer-term direction is **DuckHtsLake**, enabling round-trip conversions and metadata preservation for HTS-backed lakehouse workflows.

## Goals
- Provide table functions to read VCF/BCF, SAM/BAM/CRAM, and common annotation formats (GTF/GFF) into DuckDB tables.
- Support indexed region queries (tabix/CSI for VCF/BCF, BAM/CRAM indexes).
- Support column projection pushdown for performance.
- Preserve and expose file header metadata (including read groups) for round-trip conversions and auditing.
- Provide deterministic vendoring scripts for htslib, bcftools, samtools.
- Include conformance testing assets and scripts.

## Non-Goals (Initial)
- Full feature parity with bcftools/samtools command-line tools.
- Distributed / cloud-backed storage integration beyond what htslib already supports.
- Writer support for all formats (initial focus is read-path).

## Target DuckDB API
- DuckDB C/C++ extension API for **DuckDB 1.4+**.
- Use DuckDB table functions for file readers.
- Avoid deprecated APIs; update when DuckDB 1.4+ requires changes.

## Extension Surface (Proposed)
### Table Functions
- `bcf_read(path, region := NULL, tidy_format := false, allow_no_index := false)`
- `bam_read(path, region := NULL, allow_no_index := false)`
- `fastq_read(path, allow_no_index := false)`
- `fasta_read(path, allow_no_index := false)`
- `gtf_read(path, region := NULL, allow_no_index := false)`
- `gff_read(path, region := NULL, allow_no_index := false)`
- `gtb_read(path, region := NULL, allow_no_index := false)`
- `gtb_read(path, region := NULL, allow_no_index := false)`
- `tabix_read(path, region := NULL, allow_no_index := false)`

> htslib will detect the input format where applicable; function names map to expected schemas.

### Options
- `region`: genomic region string (e.g., `chr1:1000-2000`)
- `tidy_format` (VCF/BCF): emit one row per sample with `SAMPLE_ID`
- `allow_no_index`: allow sequential scan if no index is present

## Output Schema (Draft)
### VCF/BCF
- Fixed fields: `CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, `INFO_*`
- FORMAT fields expanded to `FORMAT_*`
- When `tidy_format = true`:
  - add `SAMPLE_ID`
  - expand FORMAT values per sample

### SAM/BAM/CRAM
- Core fields: `QNAME`, `FLAG`, `RNAME`, `POS`, `MAPQ`, `CIGAR`, `RNEXT`, `PNEXT`, `TLEN`, `SEQ`, `QUAL`
- Optional tags expanded to `TAG_*` columns where feasible
- Header metadata preserved and exposed, including `@RG` read groups and `@SQ` references

### Header Metadata (All Formats)
- Preserve raw header lines and structured key/value metadata
- Provide a consistent representation that can be embedded into Parquet metadata for round-trip conversions

### GTF/GFF
- Core fields: `SEQNAME`, `SOURCE`, `FEATURE`, `START`, `END`, `SCORE`, `STRAND`, `FRAME`, `ATTRIBUTES`
- `ATTRIBUTES` parsed into `ATTR_*` columns where feasible

### FASTA/FASTQ
- FASTA: `ID`, `DESCRIPTION`, `SEQUENCE`
- FASTQ: `ID`, `DESCRIPTION`, `SEQUENCE`, `QUALITY`

## `bcftools_score` Exact-Port Contract

DuckHTS will add a scoring surface that is intentionally aligned to upstream `bcftools +score` behavior under `.sync/score/score.c`.

### Parity Scope
- Preserve upstream option semantics and defaults, including `--use` tag handling (`GT`, `DS`, `HDS`, `AP`, `GP`, `AS`) and autodetection order.
- Preserve summary statistics ingest semantics for both:
  - VCF summary files with `FORMAT/ES` (and `FORMAT/LP` when thresholds are requested), and
  - TSV summaries with preset or custom header mappings.
- Preserve variant matching behavior (`variant_id` matching vs `chrom+pos` matching).
- Preserve missing-value handling and per-sample accumulation behavior (`score += ES * dosage(effect_allele)` for non-missing samples only).
- Preserve q-threshold behavior by comparing against `LP` (`-log10(p)` space) exactly as upstream.
- Preserve optional count outputs equivalent to `--counts`.

### DuckDB Surface
- Implement `bcftools_score(...)` as a DuckDB table function (matrix output in table form).
- Keep arguments 1:1 with upstream plugin semantics where possible (SQL-safe names for options).
- Keep table output deterministic and directly exportable via DuckDB `COPY`.

### Output Contract
- One row per target sample.
- One score column per PRS source (and per q-threshold when requested).
- Optional `*_CNT` columns when count mode is enabled.

### Conformance Requirements
- Golden parity tests must compare DuckHTS output to native `bcftools +score` on the same fixtures.
- Include fixtures for edge cases: multiallelic variants, missing dosage values, missing/ambiguous effect alleles, duplicate IDs, and q-thresholded runs.
- Numeric checks should use strict tolerances appropriate for float accumulation while ensuring rank/order and aggregate parity.

## Vendoring Layout (Proposed)
```
third_party/
  htslib/
  bcftools/
  samtools/
  cgranges/
  patches/
  licenses/
```

## Vendoring Scripts (Proposed)
- `scripts/vendor_htslib.sh`
- `scripts/vendor_bcftools.sh`
- `scripts/vendor_samtools.sh`
- `scripts/vendor_conformance_data.sh`

Each script must:
1. Pin to a version or commit.
2. Download official release tarballs.
3. Verify checksums.
4. Unpack into `third_party/`.
5. Apply patches (if any).
6. Capture licenses into `third_party/licenses/`.

Notes:
- bcftools and samtools release tarballs bundle their own htslib copy; vendoring must strip or ignore those duplicates and use the top-level `third_party/htslib`.

## Conformance Testing
- Include minimal datasets and checks:
  - VCF/BCF: small indexed cohort file(s)
  - BAM/CRAM: small alignment with index
  - Region queries + record counts
- Provide a test driver script that can run in CI with no network access.
- Prefer staging datasets from vendored htslib/bcftools/samtools test trees.

## Build Strategy
- htslib is the only required build dependency for the extension.
- bcftools and samtools are vendored for testing and output comparison only.
- CMake/Makefile integration to locate vendored headers/libs.
- Avoid network access at build time; vendored code only.

## R Package (Testing Harness)
- Provide a subdirectory R package under `r/duckhts` for testing and developer workflows.
- Include a `bootstrap.R` script (modeled after Arrow's approach) to vendor the extension build infrastructure for use from R.
- The R package should be able to build and load the extension locally for tests.

## Licensing
- Capture upstream licenses from htslib/bcftools/samtools in `third_party/licenses/`.
- Document third-party dependencies in README.
