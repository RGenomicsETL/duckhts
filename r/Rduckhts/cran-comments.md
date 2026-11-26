This release bundles htslib 1.25 and follows bcftools 1.25 for filter expressions and
variant normalization. It adds an opt-in route that opens files through the file
system of the DuckDB connection (`SET duckhts_file_system = 'duckdb'`).

## Test environments

- Local Linux, R 4.6.0 with `--as-cran`: 2 NOTEs. One is the number of updates in
  the past 6 months; the other is `-mno-omit-leaf-frame-pointer`, which comes from the
  local R build flags, not from the package.

There are no reverse dependencies on CRAN.
