# DuckDB C API headers

`duckdb.h` and `duckdb_extension.h` are copied unmodified from the DuckDB tag named by
`TARGET_DUCKDB_VERSION` in the `Makefile`. Refresh them with:

```shell
make update_duckdb_headers
```

The header's default `DUCKDB_EXTENSION_API_VERSION_*` is the C API version the extension
requests when it loads, and `TARGET_DUCKDB_VERSION` is the version written to the metadata
footer. CMake stops the build when the two differ.

Take the headers from the `v1.5.6` tag and no later one. v1.5.6 is the last version of the
v1 C API, and the `v2.0.0` header places `duckdb_create_timestamp_tz_ns` and
`duckdb_get_timestamp_tz_ns` under its 1.5.6 version guard although DuckDB 1.5.6 does not
provide them. An extension built from that header would read past the function table that
DuckDB 1.5.6 supplies.
