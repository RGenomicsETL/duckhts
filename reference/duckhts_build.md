# Build the duckhts DuckDB extension

Compiles htslib and the duckhts extension from the sources bundled in
the installed R package. The built `.duckdb_extension` file is placed in
the extension directory.

## Usage

``` r
duckhts_build(build_dir = NULL, make = NULL, force = FALSE, verbose = TRUE)
```

## Arguments

- build_dir:

  Where to build. Required. Use a writable location such as
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) when the installed
  package directory is read-only.

- make:

  Optional GNU make command to use (e.g., "gmake" or "make"). When NULL,
  auto-detects gmake or make. If a non-GNU make is used, htslib's
  configure step will fail.

- force:

  Rebuild even if the extension file already exists.

- verbose:

  Print build output.

## Value

Path to the built `duckhts.duckdb_extension` file.
