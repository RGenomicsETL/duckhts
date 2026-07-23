## Submission

This submission updates the bundled 'duckhts' DuckDB extension in Rduckhts
(1.4.0 -> 1.5.0). Since the previous CRAN release, it:

- adds the DuckVEP consequence/HGVS relation, Ensembl model builder,
  regulatory/motif, structural-variant, breakend, and NMD support, with
  VEP-116 compatibility fixes;
- adds BigWig reading through the package's htslib transport and region
  semantics;
- updates bundled htslib to 1.24 and exposes versioned downstream linking
  metadata; and
- extends BCF parsing, normalization, FASTQ projection/QC, cross-platform
  package tests, wasm tests, and differential validation.

Full details are in `NEWS.md`. There is no user-facing API breakage.

## Test environments

- Ubuntu 24.04 (x86_64), R 4.6.0 — local tarball
  `R CMD check --as-cran`
- Fedora 44 (x86_64), R-devel 2026-07-22 r90289, GCC 16.1.1 — GitHub
  CRAN-reproduction job
- Ubuntu R-release and R-devel, Windows R-release, and macOS R-release —
  GitHub Actions
- Emscripten/webR package build and package test — GitHub Actions

## R CMD check results

0 errors | 0 warnings | 1 note on the completed native environments above.
The note reports seven CRAN updates in the preceding six months. It is not a
package-code, compilation, documentation, or test note.

The Fedora warning reported for CRAN's Rduckhts 1.4.0 checks came from const
qualifier diagnostics in bundled htslib code. The 1.5.0 Fedora GCC 16.1.1
reproduction completed with no warnings.

Installed size is reported as INFO (31.5 Mb). The size is intentional: the
package bundles the 'duckhts' DuckDB extension C sources and vendored 'htslib'
1.24 so it builds and tests without network access; no compiled objects or
large binaries are shipped in the sources. Largest directories are
`duckhts_extension` (26.7 Mb) and `extdata` (3.7 Mb). The source tarball is
approximately 3.3 Mb and contains no compiled objects or libraries. If CRAN's
incoming check reports installed size as a NOTE, it is expected and by design.

## Reverse dependencies

No reverse dependencies on CRAN.
