# DuckDB C API deprecation scan

Status: historical scan. Re-run and update this note after bundled DuckDB header/runtime bumps or broad extension C API refactors.

Date: 2026-04-21

## Scope

Scan this repository for use of DuckDB C API methods that DuckDB 1.5 marks as deprecated or scheduled for removal.

I verified the deprecation texts against a fresh clone of `duckdb/duckdb-web`:

- repo: `https://github.com/duckdb/duckdb-web`
- commit: `f125322f46a5a10d0340d49a025665f5bb77ec80`
- primary doc checked: `docs/current/clients/c/api.md`
- supporting doc checked: `docs/current/clients/c/query.md`

## Upstream deprecation text checked

DuckDB's C API docs currently state:

> The reference contains several deprecation notices. These concern methods whose long-term availability is not guaranteed as they may be removed in the future.

Representative exact notices verified in `duckdb-web`:

- `duckdb_row_count`: "Deprecation notice. This method is scheduled for removal in a future release."
- `duckdb_column_data`: "This method has been deprecated. Prefer using `duckdb_result_get_chunk` instead."
- `duckdb_nullmask_data`: same as above
- `duckdb_value_*` result accessors: deprecated / scheduled for removal
- `duckdb_execute_prepared_streaming`: scheduled for removal
- `duckdb_pending_prepared_streaming`: scheduled for removal
- `duckdb_appender_error`: scheduled for removal, use `duckdb_appender_error_data`
- Arrow result/scan API family (`duckdb_query_arrow`, `duckdb_query_arrow_schema`, `duckdb_prepared_arrow_schema`, `duckdb_result_arrow_array`, `duckdb_query_arrow_array`, `duckdb_arrow_column_count`, `duckdb_arrow_row_count`, `duckdb_arrow_rows_changed`, `duckdb_query_arrow_error`, `duckdb_destroy_arrow`, `duckdb_destroy_arrow_stream`, `duckdb_execute_prepared_arrow`, `duckdb_arrow_scan`, `duckdb_arrow_array_scan`): scheduled for removal
- `duckdb_stream_fetch_chunk`: scheduled for removal

`docs/current/clients/c/query.md` also says:

> The `duckdb_value` functions are deprecated and are scheduled for removal in a future release.

## Methods included in the repo scan

The scan looked for these deprecated method names:

- `duckdb_row_count`
- `duckdb_column_data`
- `duckdb_nullmask_data`
- `duckdb_result_get_chunk`
- `duckdb_result_is_streaming`
- `duckdb_result_chunk_count`
- `duckdb_value_boolean`
- `duckdb_value_int8`
- `duckdb_value_int16`
- `duckdb_value_int32`
- `duckdb_value_int64`
- `duckdb_value_hugeint`
- `duckdb_value_uhugeint`
- `duckdb_value_decimal`
- `duckdb_value_uint8`
- `duckdb_value_uint16`
- `duckdb_value_uint32`
- `duckdb_value_uint64`
- `duckdb_value_float`
- `duckdb_value_double`
- `duckdb_value_date`
- `duckdb_value_time`
- `duckdb_value_timestamp`
- `duckdb_value_interval`
- `duckdb_value_varchar`
- `duckdb_value_string`
- `duckdb_value_varchar_internal`
- `duckdb_value_string_internal`
- `duckdb_value_blob`
- `duckdb_value_is_null`
- `duckdb_execute_prepared_streaming`
- `duckdb_pending_prepared_streaming`
- `duckdb_appender_error`
- `duckdb_query_arrow`
- `duckdb_query_arrow_schema`
- `duckdb_prepared_arrow_schema`
- `duckdb_result_arrow_array`
- `duckdb_query_arrow_array`
- `duckdb_arrow_column_count`
- `duckdb_arrow_row_count`
- `duckdb_arrow_rows_changed`
- `duckdb_query_arrow_error`
- `duckdb_destroy_arrow`
- `duckdb_destroy_arrow_stream`
- `duckdb_execute_prepared_arrow`
- `duckdb_arrow_scan`
- `duckdb_arrow_array_scan`
- `duckdb_stream_fetch_chunk`

## Scan result

### Result for DuckHTS source and tests

No uses were found in project source, R package code, SQL tests, or build files after excluding vendored/generated DuckDB API headers.

Command used:

```bash
cd /root/duckhts && \
rg -n '\b(duckdb_row_count|duckdb_column_data|duckdb_nullmask_data|duckdb_result_get_chunk|duckdb_result_is_streaming|duckdb_result_chunk_count|duckdb_value_boolean|duckdb_value_int8|duckdb_value_int16|duckdb_value_int32|duckdb_value_int64|duckdb_value_hugeint|duckdb_value_uhugeint|duckdb_value_decimal|duckdb_value_uint8|duckdb_value_uint16|duckdb_value_uint32|duckdb_value_uint64|duckdb_value_float|duckdb_value_double|duckdb_value_date|duckdb_value_time|duckdb_value_timestamp|duckdb_value_interval|duckdb_value_varchar|duckdb_value_string|duckdb_value_varchar_internal|duckdb_value_string_internal|duckdb_value_blob|duckdb_value_is_null|duckdb_execute_prepared_streaming|duckdb_pending_prepared_streaming|duckdb_appender_error|duckdb_query_arrow|duckdb_query_arrow_schema|duckdb_prepared_arrow_schema|duckdb_result_arrow_array|duckdb_query_arrow_array|duckdb_arrow_column_count|duckdb_arrow_row_count|duckdb_arrow_rows_changed|duckdb_query_arrow_error|duckdb_destroy_arrow|duckdb_destroy_arrow_stream|duckdb_execute_prepared_arrow|duckdb_arrow_scan|duckdb_arrow_array_scan|duckdb_stream_fetch_chunk)\b' \
  . \
  -g '!third_party/**' \
  -g '!.sync/**' \
  -g '!duckdb_capi/**' \
  -g '!r/Rduckhts/inst/duckhts_extension/duckdb_capi/**' \
  -g '!r/Rduckhts/inst/duckhts_extension/htslib/**' \
  -g '!r/Rduckhts/src/htslib/**'
```

Observed result: no matches.

### Vendored/generated header matches

Matches do exist in the repository's DuckDB C API header mirrors:

- `duckdb_capi/duckdb.h`
- `duckdb_capi/duckdb_extension.h`
- `r/Rduckhts/inst/duckhts_extension/duckdb_capi/duckdb.h`
- `r/Rduckhts/inst/duckhts_extension/duckdb_capi/duckdb_extension.h`

Those are API definition surfaces, not call sites in DuckHTS logic.

## Interpretation

At the time of this scan:

- DuckHTS does **not** appear to call DuckDB C API methods that DuckDB 1.5 documents as deprecated.
- The deprecated symbols are present only in bundled DuckDB header copies, which is expected.
- There is no immediate source migration work required in `src/`, `r/`, `test/`, or build glue based on this scan.

## Follow-up suggestions

1. Keep this scan handy when bumping bundled DuckDB headers/runtime.
2. Re-run the same search after any DuckDB vendor update.
3. If future code adds result materialization helpers, prefer current non-deprecated chunk/vector APIs such as:
   - `duckdb_fetch_chunk`
   - vector accessors (`duckdb_data_chunk_get_vector`, `duckdb_vector_get_data`, validity helpers)
   - `duckdb_appender_error_data` instead of `duckdb_appender_error`

## Bottom line

No deprecated DuckDB C API call sites were found in DuckHTS code itself.
