This release fixes the tests ERROR reported for 2.0.1.0.1.0 on
r-devel-linux-x86_64-debian-gcc and r-oldrel-windows-x86_64, and shortens the check time.

Memory: `test_somalier_extraction.R` ran a query under a 128MB memory limit without
fixing the thread count. DuckDB allocates buffers per worker thread, so on hosts with
many threads the query ran out of memory. The test now sets two threads before the
limit.

Check time: the incoming check of 2.0.1.0.1.0 reported an overall check time of 11
minutes on r-devel-windows-x86_64. Three long tests now run only when `NOT_CRAN=true`.
On Win-builder the whole check now takes 217 to 392 seconds.

## Test environments

- Win-builder R-devel (r90766), R-release (4.6.1) and R-oldrelease (4.5.3): 1 NOTE
  (days since last update).
- R macOS builder, R 4.6.1, arm64: OK.
- Local R-devel images: Fedora GCC 16, Fedora Clang 23 with libc++, and two
  `--disable-long-double` builds: NOTEs only.

There are no reverse dependencies on CRAN.
