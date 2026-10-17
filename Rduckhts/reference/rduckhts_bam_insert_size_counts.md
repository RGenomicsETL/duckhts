# BAM/CRAM Insert-Size Counts

Reduces eligible records into sparse orientation and absolute-TLEN
counts using Fulcrum Genomics riker 0.4.1 semantics (MIT). Only
non-first segments count; QC-fail records are excluded. The unsuffixed
wrapper uses the native reducer; the \`\_sql\` wrapper uses the SQL
macro. The SQL arm projects binary CIGAR and flag/coordinate columns
only; CRAM decoding also requests sequence and tags internally for
long-CIGAR restoration. The native arm avoids these fields. For weighted
metrics, materialize the counts and call
\`duckhts_insert_size_metrics('counts_table')\` in SQL.

## Usage

``` r
rduckhts_bam_insert_size_counts(
  con,
  path,
  reference = NULL,
  include_duplicates = FALSE,
  threads = 1,
  decompression_threads = 0,
  scan_mode = c("auto", "indexed", "sequential"),
  index_path = NULL,
  reference_index_path = NULL
)

rduckhts_bam_insert_size_counts_sql(
  con,
  path,
  reference = NULL,
  include_duplicates = FALSE,
  threads = 1,
  decompression_threads = 0,
  scan_mode = c("auto", "sequential"),
  index_path = NULL
)
```

## Arguments

- con:

  DuckDB connection with DuckHTS loaded.

- path:

  BAM or CRAM path.

- reference:

  Reference FASTA for CRAM, or NULL.

- include_duplicates:

  Include duplicate-flagged records.

- threads:

  Workers, from 1 to 8. The SQL wrapper temporarily sets and restores
  the connection's DuckDB threads setting; the native wrapper uses
  indexed histogram workers.

- decompression_threads:

  htslib decompression workers per file handle.

- scan_mode:

  SQL: auto or sequential. Native also accepts indexed; its sequential
  mode uses one histogram worker.

- index_path:

  Optional BAM/CRAM index path.

- reference_index_path:

  Optional reference FASTA index path (native only).

## Value

Data frame with pair_orientation, insert_size and read_pairs. Rows are
unordered; absent bins are omitted.
