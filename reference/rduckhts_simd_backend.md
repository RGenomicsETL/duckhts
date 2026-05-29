# DuckHTS SIMD backend diagnostics

Inspect or explicitly select the SIMD dispatch policy used by bundled
DuckHTS byte-oriented helper kernels such as `seq_gc_content(...)`.

## Usage

``` r
rduckhts_simd_backend(con)

rduckhts_simd_requested_backend(con)

rduckhts_simd_backend_compiled(con, backend)

rduckhts_simd_backend_cpu_supported(con, backend)

rduckhts_simd_backend_available(con, backend)

rduckhts_simd_info(con)

rduckhts_simd_kernel_info(con)

rduckhts_simd_set_backend(con, backend = "auto")
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded via
  [`rduckhts_load()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_load.md).

- backend:

  A single backend request. `"auto"` selects the best available
  implementation independently for each logical kernel at runtime for
  `rduckhts_simd_set_backend()`. Backend inventory predicates such as
  `rduckhts_simd_backend_available()` refer to concrete backends such as
  `"scalar"`, `"sse2"`, `"sse41"`, `"avx2"`, `"avx512"`, `"neon"`, and
  `"wasm_simd128"`.

## Value

`rduckhts_simd_backend()`, `rduckhts_simd_requested_backend()`, and
`rduckhts_simd_set_backend()` return a character scalar.
`rduckhts_simd_backend_compiled()`,
`rduckhts_simd_backend_cpu_supported()`, and
`rduckhts_simd_backend_available()` return logical scalars.
`rduckhts_simd_info()` returns the extension-owned backend inventory
table with one row per known backend; availability means compiled and
CPU/runtime supported; SQL-level backend changes are process-wide and
use the one-row `duckhts_simd_set_backend(...)` table function; the
`selectable` column reports whether the backend has a selectable
implementation path. Explicit selection still requires
`available = TRUE`. `rduckhts_simd_kernel_info()` returns one row per
logical SIMD kernel and is the authoritative diagnostic for mixed
per-kernel auto-dispatch.

## Examples

``` r
if (FALSE) { # \dontrun{
con <- DBI::dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
rduckhts_load(con)
rduckhts_simd_info(con)
rduckhts_simd_kernel_info(con)
rduckhts_simd_backend_available(con, "scalar")
rduckhts_simd_set_backend(con, "scalar")
rduckhts_simd_set_backend(con, "auto")
DBI::dbDisconnect(con, shutdown = TRUE)
} # }
```
