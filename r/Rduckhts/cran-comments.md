## Submission

<!-- Draft for the 2.0.0 submission. The package is 1.5.2.9013.0.1.5, a 2.0.0 preview, until release. -->

Rduckhts 2.0.0.0.1.5 packages DuckHTS 2.0.0. It is a major release with
breaking changes; details are in `NEWS.md`.

- Requires `duckdb` 1.5.0 or newer (previously 1.4.0).
- Consequence prediction has moved to the separate DuckVEP extension and its
  R package, Rduckvep. `rduckhts_haplotypes()` and the DuckVEP functions are
  removed from this package; CSQ, ANN and BCSQ parsing in `rduckhts_bcf()` is
  unchanged.
- The bundled extension installs SQL macros into the database catalog only for
  a writable in-memory database. `rduckhts_connect()` and `rduckhts_load()`
  install connection-local macros for database files, including read-only
  files, without modifying them.
- `duckhts_cgranges_from_query()` and `duckhts_cgranges_overlaps_bulk()` are
  removed; `duckhts_cgranges_from_table()` is a table macro.
- The Somalier panel now carries X/Y sites, which changes panel, frequency and
  sketch digests; `panel_table` is exported to a scratch Parquet file for
  `rduckhts_somalier_bam_counts()` and `rduckhts_ancestry_bam()`.

New features include indexed interval reads (`regions_var`,
`rduckhts_geno_sites()`), `rduckhts_somalier_sex()`, ancestry proportion
estimation, Somalier panel selection, named attribute columns for GenBank,
GFF and GTF, and `error_policy` for BED. A closed database file is now
released, which fixes reopening a database file in the same process on Windows.

`inst/COPYRIGHT` now inventories every bundled third-party component, and
`Authors@R` lists the copyright holders named in the bundled sources.

## Test environments

- Local source-tarball installation (private library) and
  `R CMD check --as-cran --no-manual` with `_R_CHECK_FORCE_SUGGESTS_=false`:
  Ubuntu 24.04, x86_64, R 4.6.0 (2026-04-24), GCC 13.3.0, `duckdb` R 1.5.5.

- GitHub Actions, R-CMD-check run 36750177321 on the release-preparation head:
  - Ubuntu 24.04, R 4.6.1 and R-devel (4.7.0): passed.
  - macOS 26 (aarch64), R 4.6.1: passed.
  - Windows Server 2022, R 4.6.1 (ucrt): passed.
  - Fedora 44, clang 22, R 4.6.1, CRAN-like with warnings as errors: 2 NOTEs.
    One is the update count; the URL check also could not connect from the
    container. The other lists compiler flags that Fedora's R configuration
    supplies (`-march=x86-64`, `-mtls-dialect=gnu2`, `-D_FORTIFY_SOURCE=3`
    and others), not flags the package adds.
  - Fedora, GCC 16, R-devel (2026-09-29 r90598): 1 NOTE, the update count.
  - Linux ARM64 and Windows ARM64 package contract checks: passed.

## R CMD check results

0 errors, 0 warnings, 2 notes:

- CRAN incoming feasibility reports the maintainer and "Number of updates in
  past 6 months: 7". This is a coordinated major release with the DuckDB
  community-extension release.
- The local R installation supplies `-mno-omit-leaf-frame-pointer`, which the
  compilation-flags check reports as non-portable. Rduckhts does not add this
  flag; it preserves the installation's compiler settings.

The check also reports an installed size of 31.3 Mb (`duckhts_extension`
25.7 Mb, `extdata` 4.0 Mb) as information. The source tarball is 3.4 MB. Both
are large because the package compiles the DuckHTS extension together with the
bundled HTSlib (with htscodecs), libBigWig, cgranges, the bcftools filter
engine and VariantKey sources during installation, and ships test fixtures.
Each component and its licence is listed in `inst/COPYRIGHT`. If the CRAN
incoming check reports a size NOTE on other platforms, this is the reason.

## Reverse dependencies

There are no reverse dependencies on CRAN. This was checked with
`tools::package_dependencies("Rduckhts", reverse = TRUE, db =
available.packages(repos = "https://cloud.r-project.org"))`, which returned
`character(0)` (also with `which = "all"`).
