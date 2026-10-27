This is a corrective release. It updates Rduckhts from 2.0.1.0.1.0 to
2.0.2.0.1.0 and the bundled DuckHTS extension from 2.0.1 to 2.0.2. It fixes the
test failure CRAN reports for 2.0.1.0.1.0 on r-devel-linux-x86_64-debian-gcc and
r-oldrel-windows-x86_64: an out-of-memory error at `test_somalier_extraction.R`.

That test sets a 128MB memory limit to show that a query's memory does not grow
with the number of records. It did not fix the thread count. DuckDB gives each
worker thread its own buffers, and with duckdb 2.0 this query needs about 8 MiB
plus 2.3 MiB per thread, so the limit failed on hosts with more than about 50
threads. The test now sets two threads before the limit (about 14 MiB).

The release also fixes functions that wrote files when a query was bound.
DuckDB 2.0 binds a prepared statement again at every `EXECUTE`, and DBI prepares
and then executes, so `rduckhts_mosdepth()` and
`rduckhts_samtools_idxstats(output = ...)` failed with "already exists" on a fresh
path unless `overwrite = TRUE`. These functions, and the BAM, BCF/VCF, tabix and
FASTA index builders, now write once per execution. Three long tests run only with
`NOT_CRAN=true`, which keeps the check under 10 minutes (the incoming check of
2.0.1.0.1.0 reported 11 minutes on r-devel-windows-x86_64). Public SQL and R
signatures are unchanged. Full details are in NEWS.md.

## Test environments

`R CMD check --as-cran` of the submitted tarball
(`Rduckhts_2.0.2.0.1.0.tar.gz`, SHA-256
`567e9aaeede141e4666c47d0782e61c3b766035684fdee803982579812434f54`) with duckdb 2.0.0
passes in four local images: the R-hub Fedora GCC 16 and no-long-double images, and
Fedora 44 with R-devel r90760 built with Clang 23.1.4/libc++ and with GCC 16
`--disable-long-double`. Each reports 0 errors, 0 warnings and only NOTEs: the update
frequency ("Days since last update: 3; Number of updates in past 6 months: 8"),
which follows from correcting the previous release, and CRAN's own
`-Wp,-D_FORTIFY_SOURCE=3` flag in the two images that build R themselves. Tests take
67 to 71 seconds there. Deviations from CRAN's configuration: gfortran instead of flang,
no recommended packages, no LTO.

The package tests also pass with `NOT_CRAN=true` (8,456 results, 0 failures) on
duckdb 1.5.5 and on duckdb with engine 2.0.0. The memory threshold was measured at 1
to 128 threads with the DuckDB 2.0.0 and 1.5.2 Python wheels.

Win-builder and the R macOS builder results for this tarball are recorded here
before submission.

There are no reverse dependencies on CRAN.
