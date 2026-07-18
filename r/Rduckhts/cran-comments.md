## Submission

This submission updates the bundled 'duckhts' DuckDB extension in Rduckhts
(1.3.0 -> 1.4.0). Minor release adding new reader/converter functions,
sequence-UDF overloads, and some bug fixes. Full list in `NEWS.md`. No
user-facing API breakage.

## Test environments

- Ubuntu 24.04 (x86_64), R 4.6.0 — local `R CMD check --as-cran` (tarball): Status OK
- win-builder — R-release (R 4.6.1) Status OK; R-devel submitted
- Fedora R-devel (rhub `ghcr.io/r-hub/containers/gcc16`, mirrors CRAN
  `r-devel-linux-x86_64-fedora-gcc`), R-devel 2026-07-17 and GCC 16.1.1:
  0 errors, 0 warnings, 1 expected incoming-feasibility note for the `.9000`
  development version and recent submission count
- macOS builder — R-release submitted

## R CMD check results

0 errors | 0 warnings on the completed environments above. The local Fedora
check of the development version has the single incoming-feasibility note
described above; it is not a package-code, compilation, documentation, or test
note.

Installed size is reported as INFO (~28 Mb). The size is intentional: the
package bundles the 'duckhts' DuckDB extension C sources and vendored 'htslib'
1.24 so it builds and tests without network access; no compiled objects or
large binaries are shipped in the sources. Largest directories are
`duckhts_extension` and `extdata` (small HTS fixtures). If CRAN's incoming
check reports installed size as a NOTE, it is expected and by design.

## Reverse dependencies

No reverse dependencies on CRAN.
