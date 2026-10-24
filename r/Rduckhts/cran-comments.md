This is a corrective release. It updates Rduckhts from 2.0.0.0.1.9 to
2.0.1.0.1.1 and the bundled DuckHTS extension from 2.0.0 to 2.0.1. It fixes
the failures CRAN reports for 2.0.0.0.1.9 with DuckDB 2.0.0: the tests on both
Fedora flavours (`test_bcf.R`, "Invalid BGZF header at offset 447") and the
no-long-double validation of `af_clamp = 0.5`.

The test failure was in the package's own test, not the reader. DuckDB 2.0
binds a prepared statement again at every `EXECUTE`, so a test that replaced an
index file after `PREPARE` and expected the original plan to survive read the
replaced index. The test now expects a fresh bind on DuckDB 2.0 and keeps the old
expectation on DuckDB 1.5. The same rebinding exposed three real problems that
this release also fixes: `bgzip()`/`bgunzip()` did their work at bind time, so a
second bind failed ("already exists", or the input was already removed with
`keep := FALSE`); the low-pass macro put a volatile expression in a window
`ORDER BY`; and an index whose first intervals had NULL labels wrote to
unallocated storage when a typed label followed. A Parquet wrapper also checked
a local output path with a DuckDB `glob()` query that failed on one Windows
R-devel build; it now uses the file system for local paths. Public SQL and R
signatures are unchanged. Full details are in NEWS.md.

This replaces 2.0.1.0.1.0, which the incoming check answered with one NOTE,
"Overall checktime 11 min > 10 min" on r-devel-windows-x86_64. Three long tests
(`test_count_error_fit.R`, `test_depth_state_model.R`, `test_ancestry.R`) now exit
early unless `NOT_CRAN=true`; that is the only change. Win-builder's R-devel check
now takes 260 seconds (it was 694) with 190 seconds of installation.

## Test environments

The published 2.0.0.0.1.9 source reproduces CRAN's failure on this host, in
three Fedora 44 images with R-devel r90760: Clang 23.1.4 with libc++, GCC 16, and
GCC 16 with `--disable-long-double`. The exact tarball submitted here
(`Rduckhts_2.0.1.0.1.1.tar.gz`, SHA-256 `ac618b648eafac081b0541c205137de49be942ab4ebbdcd09c634152c1f27431`)
passes `R CMD check --as-cran` in those three images and in the R-hub
no-long-double image: 0 errors, 0 warnings, and only NOTEs (update frequency,
and CRAN's own `-Wp,-D_FORTIFY_SOURCE=3` flag on Clang). Dependencies were
refreshed to the newest CRAN versions, including DuckDB 2.0.0. Deviations from
CRAN's configuration: gfortran instead of flang, no recommended packages, no LTO.

Win-builder checked the same tarball: R-devel (r90760), R-release (4.6.1) and
R-oldrelease (4.5.3) each report Status: 1 NOTE, "Days since last update: 4; Number
of updates in past 6 months: 7", which follows from correcting the previous
release quickly. Tests pass on all three, in 87 to 176 seconds.

The same package sources passed GitHub Actions on Windows (release, devel and
4.5.3, including a debugger replay of the full test suite), macOS and Linux, and
the extension's sanitizer and DuckDB 2.0 gates. The R macOS builder could not be
used: it reported `DBI` and `duckdb` as not available in its R 4.6 library before
checking the package.

There are no reverse dependencies on CRAN.
