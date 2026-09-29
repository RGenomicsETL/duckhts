# Install DuckHTS macros on a connection

Executes the ordered definitions exported by the loaded DuckHTS
extension as connection-local TEMP macros. Use this on DBI or pool
connections that did not come from
[`rduckhts_connect()`](https://rgenomicsetl.github.io/duckhts/Rduckhts/reference/rduckhts_connect.md).
The operation is idempotent and shadows same-named persistent macros
without deleting them. Inspect old persistent macros in
`duckdb_functions()` using the file's `database_name`; drop them
explicitly if no longer wanted.

## Usage

``` r
rduckhts_install_macros(con)
```

## Arguments

- con:

  A DBI connection with DuckHTS loaded.

## Value

Invisibly returns the definitions SHA-256 digest.

## Details

This function does not begin or commit a transaction. Inside a caller
transaction its DDL participates in that transaction: rollback removes
newly installed TEMP macros, and the caller can install them again after
rollback. Call it separately for each connection, including pool
connections.
