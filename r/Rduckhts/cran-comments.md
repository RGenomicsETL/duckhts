## Submission

This is a patch release following Rduckhts 1.5.0-0.1.0. It updates the bundled
'duckhts' DuckDB extension from 1.5.0 to 1.5.1 and fixes the current CRAN check
issues.

It fixes the reported native diagnostics, including uninitialized normalization
cleanup state and liftover alias pointers, unused coverage, interval, and VEP
code, the libBigWig `bwCleanup` prototype, and MinGW visibility and
allocation-size checks.

This resubmission no longer promotes compiler warnings supplied by R or the
installation environment to errors. Strict package-owned diagnostics remain an
explicit CI check, and vendored htslib headers are treated as system headers in
Unix-like package builds. This fixes the reported M1Mac installation failure.
It also derives Windows extension metadata from R's target rather than the MSYS
shell architecture and validates native Windows ARM64 package loading.

It adds `rduckhts_connect()` to permit loading the bundled extension with
`duckdb` 1.5.5, fixing the Fedora-Clang example and tinytest errors. There is no
user-facing API breakage. Full details are in `NEWS.md`.

## Test environments

The package was checked on Ubuntu 24.04 with R-devel and on Fedora 44 with R
4.6.1 and Clang 22.1.8. GitHub Actions also checked Ubuntu R-release and
R-devel, Windows x86_64 and ARM64 R-release, and macOS ARM64 R-release.

## R CMD check results

Checks completed with 0 errors, 0 warnings, and 1 note. The note reports CRAN
package updates from the preceding six months.

## Reverse dependencies

There are no reverse dependencies on CRAN.
