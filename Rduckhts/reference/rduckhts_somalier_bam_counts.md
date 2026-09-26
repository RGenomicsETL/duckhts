# Extract Panel-Aligned Counts from BAM or CRAM

Count observed A, B, and other query bases at every site in a typed
panel. The panel is prepared once, then each scan worker owns one
indexed multi-region scan and independent BAM/CRAM, index, reference,
pileup, and overlap state. \`worker_count\` bounds the number of
DuckDB-scheduled panel shards; the connection's thread setting bounds
how many can run concurrently. \`decompression_threads\` separately
controls htslib decompression workers per source handle. Valid uncovered
sites are measured zero depth; reference or alignment-header mismatches
remain rows with NULL counts and a named status. The panel can be a
committed table/view or an ordinary Parquet file. Caller-local temporary
relations and uncommitted changes are not visible during panel
preparation. One retained-connection preparation slot is shared by
concurrent calls; nested or concurrent preparation errors and callers
may retry.

## Usage

``` r
rduckhts_somalier_bam_counts(
  con,
  source_path,
  sample_id,
  reference_path,
  panel_table = NULL,
  panel_parquet = NULL,
  index_path = NULL,
  reference_index_path = NULL,
  min_mapq = 1,
  min_baseq = 0,
  require_flags = 0,
  exclude_flags = 1796,
  overlap_policy = c("hileup_v0.1.0", "none"),
  decompression_threads = 0,
  worker_count = 1,
  max_depth = 1e+05,
  max_overlap_qnames = 1e+05,
  max_sites = 1e+06,
  max_region_bytes = 67108864,
  remote_block_bytes = 1048576,
  remote_cache_bytes = 67108864,
  reference_cache_bytes = 67108864,
  table_name = NULL,
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded.

- source_path:

  One indexed BAM/CRAM path or URI.

- sample_id:

  Nonempty sample identity assigned to the extracted rows.

- reference_path:

  Reference FASTA matching the panel and alignments.

- panel_table:

  Name of a committed ordered panel table or view.

- panel_parquet:

  Ordinary panel Parquet path, instead of \`panel_table\`.

- index_path:

  Optional explicit BAM/CRAM index path.

- reference_index_path:

  Optional explicit FASTA index path.

- min_mapq:

  Minimum alignment mapping quality.

- min_baseq:

  Minimum observed-base quality. Missing qualities pass only when this
  is zero.

- require_flags:

  SAM flag bits that every retained alignment must have.

- exclude_flags:

  SAM flag bits that exclude an alignment.

- overlap_policy:

  Either \`"hileup_v0.1.0"\` encounter-order mate suppression or
  \`"none"\`.

- decompression_threads:

  Number of htslib decompression worker threads.

- worker_count:

  Number of independently schedulable panel shards, from 1 through 64.
  Each nonempty shard opens its own reader/reference state.

- max_depth:

  Maximum admitted pileup depth before an explicit error.

- max_overlap_qnames:

  Per-site, per-job capacity for overlap-suppression names.

- max_sites:

  Maximum panel cardinality.

- max_region_bytes:

  Per-job capacity for the indexed multi-region request.

- remote_block_bytes:

  Remote alignment block size per worker handle.

- remote_cache_bytes:

  Remote alignment cache size per worker handle.

- reference_cache_bytes:

  Remote reference cache size per worker handle.

- table_name:

  Optional output table. \`NULL\` returns a data frame.

- overwrite:

  Whether an existing output table may be replaced.

## Value

A data frame ordered by \`site_index\` if \`table_name\` is \`NULL\`;
otherwise invisible \`TRUE\`.
