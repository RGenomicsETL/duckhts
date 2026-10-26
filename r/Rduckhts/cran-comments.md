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
path unless `overwrite = TRUE`. These functions, and the BAM, tabix and FASTA
index builders, now write once per execution. Three long tests run only with
`NOT_CRAN=true`, which keeps the check under 10 minutes (the incoming check of
2.0.1.0.1.0 reported 11 minutes on r-devel-windows-x86_64). Public SQL and R
signatures are unchanged. Full details are in NEWS.md.

## Test environments

The package tests pass on this host with `NOT_CRAN=true` (8,312 results, 0
failures) on duckdb 1.5.5 and on duckdb with engine 2.0.0. The memory threshold was
measured at 1 to 128 threads with the DuckDB 2.0.0 and 1.5.2 Python wheels.

`R CMD check --as-cran` of the submitted tarball, in the local CRAN-flavour images
and on Win-builder, is still to be run and recorded here before submission.

There are no reverse dependencies on CRAN.
