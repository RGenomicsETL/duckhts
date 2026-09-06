# Read Record-Major Genotypes

Read typed GT/PS calls without repeating variant text per sample. Each
record has a zero-based scan-local \`record_index\` and a list of calls
with original-header sample indices, nullable allele indices, per-slot
phase bits, and nullable scalar phase sets. An absent GT has NULL
allele/phase lists; a missing allele still occupies a slot. Phase bits
follow HTSlib decoding, including its leading-slot convention for VCF
versions before 4.4.

## Usage

``` r
rduckhts_geno(
  con,
  table_name = NULL,
  path,
  region = NULL,
  index_path = NULL,
  samples = NULL,
  non_reference_only = FALSE,
  scan_mode = "auto",
  decompression_threads = 0,
  decode_error_policy = "null",
  overwrite = FALSE
)
```

## Arguments

- con:

  A DuckDB connection with DuckHTS loaded

- table_name:

  Optional table to create; \`NULL\` returns a data frame.

- path:

  Path to the VCF/BCF file

- region:

  Optional genomic region (e.g., "chr1:1000-2000")

- index_path:

  Optional explicit path to index file (.csi/.tbi)

- samples:

  Optional HTSlib sample selector: \`NULL\` or \`"-"\` keeps all, \`""\`
  keeps none, comma-separated names include samples, and a leading
  \`"^"\` excludes them. Unknown names error; selected samples retain
  header order.

- non_reference_only:

  Omit calls without any called alternate allele. This does not infer
  phase or remove variant records.

- scan_mode:

  Optional scan mode. Use `"auto"` (default extension behavior) or
  `"sequential"` to force full-file streaming instead of index-backed
  count/parallel scan paths. Sequential mode is incompatible with
  \`region\`.

- decompression_threads:

  Integer. Number of htslib decompression worker threads per file
  handle. Default \`0\`. Use \`0\` to keep BCF/VCF reads
  single-threaded.

- decode_error_policy:

  Character. Dirty/corrupt BCF decode policy: `"null"` returns NULL for
  header-vs-payload type clashes, `"warn"` emits a DuckHTS warning and
  returns NULL, and `"error"` raises a DuckDB/R error.

- overwrite:

  Logical. If TRUE, overwrites existing table

## Value

A data frame when \`table_name\` is \`NULL\`, otherwise invisible
\`TRUE\`.

## Details

Full scans preserve the input stream when assigning ordinals. Indexed
regions use HTSlib's union order and start a new ordinal at zero. Use
SQL \`ORDER BY record_index\` when order matters; the ordinal is not a
persistent file locator. Empty sample selection and sparse calls
preserve every selected variant row.

## See also

\[rduckhts_bcf_samples()\]

## Examples

``` r
con <- rduckhts_connect()
path <- system.file("extdata", "geno_calls.bcf", package = "Rduckhts")
rduckhts_geno(con, "calls", path, non_reference_only = TRUE)
DBI::dbGetQuery(con, "SELECT record_index, len(calls) AS n FROM calls ORDER BY record_index")
#>   record_index n
#> 1            0 2
#> 2            1 0
#> 3            2 0
#> 4            3 0
#> 5            4 0
#> 6            5 2
DBI::dbDisconnect(con, shutdown = TRUE)
```
