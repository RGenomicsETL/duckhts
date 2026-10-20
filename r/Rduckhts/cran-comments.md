This submission updates Rduckhts to 2.0.0.0.1.9 and the bundled DuckHTS
extension from 1.5.2 to 2.0.0.

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
webR/duckdb-wasm and DuckDB 2.0. The R macOS builder (r-release, arm64)
checked the same tarball: Status OK, tests included.

Earlier 2.0.0 revisions intermittently ended the Windows test process with
status 0xC0000005 at `dbDisconnect(shutdown = TRUE)`, on Win-builder
R-oldrelease and on Windows Server 2022. A debugger capture showed a C++
`thread_local` destructor inside the `duckdb` package running on freed memory
when a DuckDB worker thread exits, because the winpthreads in Rtools runs
thread-exit key destructors in creation order (fixed in mingw-w64 on
2024-12-05, after Rtools45). The bundled extension triggered it by closing a
private DuckDB instance on a worker thread; it now does that on the calling
thread only, and the Windows Server 2022 check replays the full test suite
under a debugger three times to confirm. Win-builder R-oldrelease (R 4.5.3,
Windows Server 2022) then checked this revision: tests OK in 367 s, with one
NOTE for an example above 10 s of CPU time; that example is shortened in the
submitted tarball.

There are no reverse dependencies on CRAN.
