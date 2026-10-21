This corrective release updates Rduckhts from 2.0.0.0.1.9 to
2.0.0.0.1.10. It fixes the `rduckhts_roh()` validation error reported by
CRAN's no-long-double check: `af_clamp = 0.5` is rejected with a direct
open-endpoint comparison before SQL is evaluated. The bundled DuckHTS
extension and its public SQL API are unchanged.

## Test environments

- Ubuntu 24.04, R 4.6.0, `R CMD check --as-cran`: 0 errors, 0 warnings.
- R-devel r90728 configured with `--disable-long-double`: 0 errors,
  0 warnings.
- GitHub Actions: R-release and R-devel on Linux, macOS and Windows;
  Fedora GCC 16 and Clang 22; Linux and Windows ARM64; sanitizers;
  webR/duckdb-wasm; DuckDB 2.0.
- Win-builder R-devel: source tarball submitted on 2026-10-21; results are
  delivered to the maintainer by email.

The CRAN noLD failure is recorded at
https://www.stats.ox.ac.uk/pub/bdr/noLD/Rduckhts.out.

There are no reverse dependencies on CRAN.
