## Submission

Rduckhts 1.4.0-0.1.0 — minor release. New features since 1.3.0: `read_bcf_v2`,
`read_bcf_appender`, `hts_region_union_query`, Parquet converter builders
(`duckhts_{bcf,bam,gff,tabix}_convert_parquet_sql`), runtime DuckDB type-support
probes, `decode_error_policy` for corrupt BCF handling, `scan_mode` reader
controls, gVCF-aware normalization, nt16/binary-CIGAR overloads for the sequence
and CIGAR helper UDFs, and an AVX2 SIMD kernel behind the nt16 GC path. Full list
in `NEWS.md`. No user-facing API breakage.

## Test environments

- Ubuntu 24.04 (x86_64), R 4.6.0 — local `R CMD check --as-cran` (tarball)
- win-builder — R-devel and R-release (`devtools::check_win_devel()` / `check_win_release()`)
- macOS builder — R-release (`devtools::check_mac_release()`)
- Fedora R-devel (`rhub` container `ghcr.io/r-hub/containers/gcc16`, mirrors the
  CRAN `r-devel-linux-x86_64-fedora-gcc` flavour), via a Docker check

## R CMD check results

0 errors | 0 warnings | 1 note

The expected note is installed size:

- installed size is ~28 Mb; largest directories `duckhts_extension` (bundled
  DuckDB extension C sources + vendored `htslib` 1.23.1) and `extdata` (small HTS
  fixtures).

The size is intentional: the package bundles the `duckhts` DuckDB extension
sources and vendored `htslib` so installation and tests run without network
access. No compiled objects or large binaries are shipped in the sources.

## Reverse dependencies

No reverse dependencies on CRAN.
