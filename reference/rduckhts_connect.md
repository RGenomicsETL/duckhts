# Create a DuckDB connection with bundled DuckHTS loaded

Creates a DuckDB connection that explicitly permits loading the
package-built DuckHTS extension, disables automatic installation and
loading of unrelated known DuckDB extensions, and loads the bundled
DuckHTS extension.

## Usage

``` r
rduckhts_connect(
  dbdir = ":memory:",
  read_only = FALSE,
  bigint = "numeric",
  config = list(),
  extension_path = NULL
)
```

## Arguments

- dbdir:

  Path to a DuckDB database, or `":memory:"` for an in-memory database.

- read_only:

  Logical; open a file-backed database read-only.

- bigint:

  How DuckDB 64-bit integers are returned; passed to
  [`duckdb::duckdb()`](https://r.duckdb.org/reference/duckdb.html).

- config:

  Named list of additional DuckDB configuration settings.

- extension_path:

  Optional path to a DuckHTS `.duckdb_extension` file. If `NULL`, uses
  the extension bundled with Rduckhts.

## Value

A DuckDB connection with DuckHTS loaded.

## Details

Current versions of the duckdb R package can disable extension loading
on Linux builds that use a C++ standard library other than `libstdc++`.
DuckHTS is compiled locally with the same package toolchain, so this
helper explicitly enables extension loading for this package-owned
connection. Package-owned connections use per-session DuckDB
extension/secret storage by default when the installed duckdb version
supports that setting.

DuckDB reuses a live database instance when the same file-backed `dbdir`
is opened again, ignoring new driver and configuration settings. This
helper rejects such reuse rather than silently weakening its connection
policy. Use
[`rduckhts_load()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_load.md)
on the existing connection if it already permits unsigned, driver-level
extension loading, or close the existing database instance before
calling this helper.

The settings `allow_unsigned_extensions`,
`autoinstall_known_extensions`, and `autoload_known_extensions` are
controlled by this helper and override entries with those names in
`config`. Other named DuckDB settings are passed through unchanged.

## Examples

``` r
con <- rduckhts_connect()
DBI::dbGetQuery(con, "SELECT duckhts_htslib_version() AS version")
#>   version
#> 1    1.24
DBI::dbDisconnect(con, shutdown = TRUE)
```
