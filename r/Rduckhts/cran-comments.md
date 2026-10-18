This submission updates Rduckhts to 2.0.0.0.1.7 and the bundled DuckHTS
extension from 1.5.2 to 2.0.0. The packaging revision also keeps webR
builds on the bundled htslib headers.

This is a major release. Consequence prediction has moved to the separate
DuckVEP extension and Rduckvep package, and the old bundled DuckVEP interfaces
are removed. It adds genotype, depth, ancestry and sequencing-QC workflows,
and connection-local macro loading for database files. It also replaces the
deprecated SQL lambda syntax, addressing the DuckDB 2.0 reverse-dependency
failure reported in https://github.com/RGenomicsETL/duckhts/issues/352. Full
details are in NEWS.md.

## Test environments

The exact source tarball passed `R CMD check --as-cran` on Ubuntu 24.04 with
R 4.6.0: 0 errors, 0 warnings and 1 note. The note reports
`-mno-omit-leaf-frame-pointer`, which is supplied by R rather than the package.

GitHub Actions checked Ubuntu R-release and R-devel, Fedora 44 with GCC and
Clang, Windows x86_64 and ARM64, macOS ARM64, Linux ARM64, sanitizers,
webR/duckdb-wasm and DuckDB 2.0.

There are no reverse dependencies on CRAN.
