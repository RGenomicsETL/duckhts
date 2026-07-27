## Submission

This is the requested resubmission of Rduckhts 1.5.1-0.1.1, now versioned
1.5.1-0.1.3. It is a patch release following the accepted Rduckhts
1.5.0-0.1.0 and updates the bundled 'duckhts' DuckDB extension from 1.5.0 to
1.5.1.

The reported M1Mac installation failure was caused by the package's former
default of promoting compiler warnings supplied by R and the installation
environment to errors. Released package builds now preserve that external
warning policy, while strict package-owned diagnostics remain an explicit CI
check. Vendored htslib headers are treated as system headers in Unix-like
package builds. The macOS ARM64 package check enables the conversion-warning
flags that produced the report; package installation emits no compiler warning
or error diagnostics, and `R CMD check` finishes with `Status: OK`.

This resubmission also derives Windows extension metadata from R's target
rather than the architecture reported by the MSYS shell. A native Windows
ARM64 job builds a clean source tarball, installs it, loads the resulting
`windows_arm64_mingw` extension, and validates the bundled htslib receipt.

The release fixes the other reported native diagnostics, including
uninitialized normalization cleanup state and liftover alias pointers, unused
coverage, interval, and VEP code, the libBigWig `bwCleanup` prototype, and
MinGW visibility and allocation-size checks. It also adds
`rduckhts_connect()` to permit loading the bundled extension with `duckdb`
1.5.5, fixing the Fedora-Clang example and tinytest errors. There is no
user-facing API breakage. Full details are in `NEWS.md`.

## Test environments

The package was checked on Ubuntu 24.04 with R-devel and on Fedora 44 with R
4.6.1 and Clang 22.1.8. GitHub Actions also checked Ubuntu R-release and
R-devel, Windows x86_64 R-release, and macOS ARM64 R-release, and built,
installed, and validated a clean source tarball on Windows ARM64 R-release.

## R CMD check results

Checks completed with 0 errors, 0 warnings, and 1 note. The note reports CRAN
package updates from the preceding six months.

## Reverse dependencies

There are no reverse dependencies on CRAN.
