This corrective release updates Rduckhts from 2.0.0.0.1.9 to
2.0.0.0.1.10. It fixes CRAN's no-long-double validation failure:
`af_clamp = 0.5` is rejected with a direct open-endpoint comparison
before SQL evaluation.

The bundled extension also supports DuckDB 2.0 low-pass window binding,
explicit NULL handling in cgranges mutators, and execution-scoped BGZF
compression/decompression. Tests distinguish BCF index binding from
prepared-statement rebinding. Public SQL and R signatures are unchanged.
Tinytest runs directly in-process.

## Completed checks

Package sources at revision 7111af111e98d36eb387bde0f037461693a55648:

- Ubuntu 24.04, R 4.6.0, source-tarball `R CMD check --as-cran --no-manual`:
  0 errors, 0 warnings, 2 NOTEs. These concern update frequency and
  `-mno-omit-leaf-frame-pointer`, which is present in `R CMD config CFLAGS`.
- R-devel r90728 configured with `--disable-long-double`, source-tarball
  `R CMD check --as-cran --no-manual`: 0 errors, 0 warnings, 1 update-frequency NOTE.
- Full in-process package suite: 8,309 results passed on each of R DuckDB
  1.5.5 and the r-universe 1.99.99.9001 package with the v2.0.0 engine.
- All 67 extension SQL conformance files passed on DuckDB 1.5.2 and 2.0.0.
- AddressSanitizer and UndefinedBehaviorSanitizer extension gates passed.

## Pending release checks

Fresh Win-builder results and current-head GitHub Actions/Codex review are
required before submission. Earlier Windows checks do not certify this tarball.
No CRAN submission has been made for this release candidate.

The CRAN noLD failure is recorded at
https://www.stats.ox.ac.uk/pub/bdr/noLD/Rduckhts.out.

There are no reverse dependencies on CRAN.
