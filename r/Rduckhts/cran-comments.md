## Test environments

- Ubuntu 24.04.3 LTS (x86_64), R 4.6.0

## R CMD check results

0 errors | 0 warnings | 0 notes

Checked with:

```sh
THREADS=2 make check
```

The check reported an installed-size INFO:

- installed size: 26.3 Mb
- largest directories: `duckhts_extension` (22.6 Mb) and `extdata` (3.1 Mb)

This size is expected. The package bundles the `duckhts` DuckDB extension sources, vendored `htslib` 1.23.1, and small HTS fixtures so installation and tests can run without network access.

## Reverse dependencies

No reverse dependencies on CRAN.
