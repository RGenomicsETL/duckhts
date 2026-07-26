## Submission

This submission updates the bundled 'duckhts' DuckDB extension in Rduckhts
(1.4.0 -> 1.5.1). The 0.1.1 R packaging revision also responds to the current
CRAN check results for Rduckhts 1.5.0. Since the previous CRAN release, it:

- adds `rduckhts_connect()` so package-owned connections explicitly permit the
  locally compiled bundled extension while disabling implicit installation and
  loading of unrelated known DuckDB extensions. This resolves the Fedora-Clang
  example and tinytest errors caused when `duckdb` 1.5.5 correctly disabled
  extension loading by default on its Linux libc++ build;
- adds the DuckVEP consequence/HGVS relation, Ensembl model builder,
  regulatory/motif, structural-variant, breakend, and NMD support, with
  VEP-116 compatibility fixes;
- adds BigWig reading through the package's htslib transport and region
  semantics;
- updates bundled htslib to 1.24 and exposes versioned downstream linking
  metadata; and
- extends BCF parsing, normalization, FASTQ projection/QC, cross-platform
  package tests, wasm tests, and differential validation; and
- makes DuckHTS-owned native diagnostics fatal in package builds, fixes the
  reported Windows/macOS compiler diagnostics, and adds a Fedora 44 Clang 22
  CRAN-like check with checksum-verified R and compiler artifacts.

Full details are in `NEWS.md`. There is no user-facing API breakage.

## Test environments

- Ubuntu 24.04.3 (x86_64), R 4.6.0 — local 1.5.1 tarball
  `R CMD check --as-cran`
- Fedora 44 (x86_64), Fedora-packaged R 4.6.1, Clang 22.1.8 — local and
  GitHub CRAN-like tarball checks
- Ubuntu R-release and R-devel, Windows R-release, and macOS R-release —
  GitHub Actions
- Emscripten/webR package build and package test — GitHub Actions

## R CMD check results

0 errors | 0 warnings | 1 note on the completed native environments above.
The note reports eight CRAN updates in the preceding six months. It is not a
package-code, compilation, documentation, or test note.

The libBigWig strict-prototype warning in the current CRAN macOS and Fedora-
Clang results is fixed in 1.5.1 and enforced by the strict package build. The
current CRAN Fedora-Clang example/test errors are covered by package examples
and tinytests that use `rduckhts_connect()` with `duckdb` 1.5.5's driver-level
extension guard enabled. The earlier Windows compiler diagnostics are also
fixed. The Fedora 44 Clang 22.1.8 CRAN-like check completed with no package
warnings.

Installed size is reported as INFO (31.5 Mb). The size is intentional: the
package bundles the 'duckhts' DuckDB extension C sources and vendored 'htslib'
1.24 so it builds and tests without network access; no compiled objects or
large binaries are shipped in the sources. Largest directories are
`duckhts_extension` (26.7 Mb) and `extdata` (3.7 Mb). The source tarball is
approximately 3.3 Mb and contains no compiled objects or libraries. If CRAN's
incoming check reports installed size as a NOTE, it is expected and by design.

## Reverse dependencies

No reverse dependencies on CRAN.
