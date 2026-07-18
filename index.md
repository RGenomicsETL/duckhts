# Rduckhts: DuckDB HTS File Reader Extension for R

[![CRAN
Status](https://www.r-pkg.org/badges/version/Rduckhts)](https://cran.r-project.org/package=Rduckhts)[![R-universe
version](https://RGenomicsETL.r-universe.dev/Rduckhts/badges/version)](https://RGenomicsETL.r-universe.dev/Rduckhts)

`Rduckhts` provides an R interface to a [DuckDB](https://duckdb.org/)
`HTS` (High Throughput Sequencing) file reader extension. This enables
reading common bioinformatics file formats such as `VCF`/`BCF`,
`SAM`/`BAM`/`CRAM`, `FASTA`, `FASTQ`, `GFF`, `GTF`, and tabix-indexed
files directly from `R` using `SQL` queries via
[`duckhts`](https://github.com/RGenomicsETL/duckhts).

## How it works

Following [RBCFTools](https://github.com/RGenomicsETL/RBCFTools), tables
are created and returned instead of data frames. `VCF`/`BCF`,
`SAM`/`BAM`/`CRAM`, `FASTA`, `FASTQ`, `GFF`, `GTF`, and `tabix` formats
can be queried. We support region queries for indexed files, and we
target Linux, macOS, and RTools.
[`htslib`](https://github.com/samtools/htslib) 1.24 is bundled so build
dependencies stay minimal. The package build adapts the generic
extension infrastructure to a GNU make-based R package workflow, while
the standalone community extension uses the submitted
[`duckhts`](https://github.com/RGenomicsETL/duckhts) build path.

## Installation

The package can be installed from r-universe

``` r
# Install 'Rduckhts' in R:
install.packages('Rduckhts', repos = c('https://rgenomicsetl.r-universe.dev', 'https://cloud.r-project.org'))
# When on CRAN
install.packages("Rduckhts")
```

## System Requirements

Installation requires `htslib` dependencies such as zlib and libbz2, and
optionally liblzma, libcurl, and OpenSSL for full functionality. The
package requires GNU make. On Windows Rtools builds, `htslib` plugins
are not enabled.

## Browser wasm/webR networking and setup

In browser wasm/webR builds, remote `http`/`https` access does not use
htslib `libcurl`.

- The bundled extension enables a custom htslib `hFILE` backend in
  `src/wasm_http_hfile.c` for `http` and `https`.
- This backend uses synchronous worker-side XHR for range reads and
  tabix index access.
- htslib `libcurl`/`S3`/`GCS` paths are intentionally disabled for wasm
  in the package build.

What this means in practice:

- Same-origin URLs are the simplest setup and work well for local
  browser testing.
- Remote URLs work only when browser CORS policy allows them.
- CORS must allow both the primary file and index sidecars
  (`.tbi`/`.csi`), including range requests.
- htslib may probe a `.csi` sidecar before falling back to `.tbi`; a
  `.csi` `404` is not a regression if the `.tbi` path succeeds.
- `ALL_PROXY` and websocket proxy settings do not affect this wasm XHR
  backend.
- The simple local `python3 -m http.server` setup ignores HTTP `Range`;
  in that case the backend warns and falls back to full-object fetch +
  local byte slicing. This is acceptable for small smoke tests, but
  production servers should support `Range`.

### `Module.duckhtsWasmHttpConfig`

The browser HTTP backend reads optional request/auth settings from
`Module.duckhtsWasmHttpConfig`.

- In plain browser JavaScript, set it in the page/worker before running
  queries.
- In webR, consumers can set it from R with `webr::eval_js()`; they do
  not need to hand-edit the host HTML/JS as long as they can run that
  call before the relevant HTTP reads.

Plain JavaScript example:

``` js
Module.duckhtsWasmHttpConfig = {
  headers: {
    Authorization: "Bearer <short-lived-token>",
    "X-Request-Source": "webr-local"
  },
  allowHosts: ["ftp.ebi.ac.uk", ".ebi.ac.uk"],
  enforceHostAllowlist: true,
  withCredentials: false,
  allowInsecureAuth: false
};
```

Inside webR, set the same config through `webr::eval_js()` because the
code runs inside the webR worker:

``` r
webr::eval_js("
  Module.duckhtsWasmHttpConfig = {
    headers: {
      Authorization: 'Bearer <short-lived-token>',
      'X-Request-Source': 'webr-local'
    },
    allowHosts: ['ftp.ebi.ac.uk', '.ebi.ac.uk'],
    enforceHostAllowlist: true,
    withCredentials: false,
    allowInsecureAuth: false
  };
")
```

Configuration fields:

- `headers`: named request headers to attach to matching hosts.
- `allowHosts`: hostname allowlist for header injection. Entries can be
  exact hosts (`"ftp.ebi.ac.uk"`) or suffix matches with a leading dot
  (`".ebi.ac.uk"`).
- `enforceHostAllowlist`: when `true`, block requests to hosts outside
  `allowHosts` instead of merely omitting configured headers.
- `withCredentials`: when `true`, send cookies/credentials with XHR
  requests.
- `allowInsecureAuth`: when `true`, allow `Authorization` headers on
  non-HTTPS URLs. The default is `false`.

Security behavior of this config:

- Headers are only attached when URL hostnames match `allowHosts`.
- `Authorization` is blocked for non-HTTPS URLs unless
  `allowInsecureAuth: true` is set.
- Cookies/credentials are only sent when `withCredentials: true` is set.
- The config can be updated or cleared between queries if different
  hosts need different policies.

## Quick Start

The extension is loaded with
`rduckhts_load(con, extension_path = NULL)`. The wrappers break down
into:

- readers:
  [`rduckhts_bcf()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bcf.md),
  [`rduckhts_bam()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bam.md),
  [`rduckhts_fasta()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_fasta.md),
  [`rduckhts_fastq()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_fastq.md),
  [`rduckhts_gff()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_gff.md),
  [`rduckhts_gtf()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_gtf.md),
  [`rduckhts_tabix()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_tabix.md),
  [`rduckhts_bed()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bed.md)
- multi-file readers:
  [`rduckhts_bam_multi()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bam_multi.md),
  [`rduckhts_bcf_multi()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bcf_multi.md),
  [`rduckhts_fastq_multi()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_fastq_multi.md),
  [`rduckhts_fasta_multi()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_fasta_multi.md),
  [`rduckhts_gff_multi()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_gff_multi.md),
  [`rduckhts_gtf_multi()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_gtf_multi.md),
  [`rduckhts_tabix_multi()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_tabix_multi.md),
  [`rduckhts_bed_multi()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bed_multi.md)
- reference helpers:
  [`rduckhts_fasta_index()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_fasta_index.md),
  [`rduckhts_fasta_nuc()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_fasta_nuc.md)
- compression/indexing:
  [`rduckhts_bgzip()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bgzip.md),
  [`rduckhts_bgunzip()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bgunzip.md),
  [`rduckhts_bam_index()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bam_index.md),
  [`rduckhts_bcf_index()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bcf_index.md),
  [`rduckhts_tabix_index()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_tabix_index.md)
- metadata helpers:
  [`rduckhts_hts_header()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_hts_header.md),
  [`rduckhts_hts_index()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_hts_index.md),
  [`rduckhts_hts_index_spans()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_hts_index_spans.md),
  [`rduckhts_hts_index_raw()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_hts_index_raw.md)
- Parquet converters:
  [`rduckhts_bcf_convert_parquet()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bcf_convert_parquet.md),
  [`rduckhts_bam_convert_parquet()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bam_convert_parquet.md),
  [`rduckhts_gff_convert_parquet()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_gff_convert_parquet.md),
  [`rduckhts_tabix_convert_parquet()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_tabix_convert_parquet.md)
- SIMD diagnostics:
  [`rduckhts_simd_backend()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_simd_backend.md),
  [`rduckhts_simd_requested_backend()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_simd_backend.md),
  [`rduckhts_simd_backend_available()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_simd_backend.md),
  [`rduckhts_simd_set_backend()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_simd_backend.md)

Start with one reader, then materialize tables and compose the richer
helpers around them.

``` r
library(DBI)
library(duckdb)
library(Rduckhts)

# Set HTS_PATH before htslib's first file access so bundled remote handlers
# are discoverable when this build uses dynamic plugins.
setup_hts_env()

fasta_path <- system.file("extdata", "ce.fa", package = "Rduckhts")
fastq_r1 <- system.file("extdata", "r1.fq", package = "Rduckhts")
fastq_r2 <- system.file("extdata", "r2.fq", package = "Rduckhts")
con <- dbConnect(duckdb::duckdb(config = list(allow_unsigned_extensions = "true")))
rduckhts_load(con)
#> [1] TRUE

rduckhts_fasta(con, "sequences", fasta_path, overwrite = TRUE)
rduckhts_fastq(con, "reads", fastq_r1, mate_path = fastq_r2, overwrite = TRUE)

dbGetQuery(con, "SELECT COUNT(*) AS n FROM sequences")
#>   n
#> 1 7
dbGetQuery(con, "SELECT COUNT(*) AS n FROM reads")
#>    n
#> 1 10
```

## SIMD dispatch flow

The bundled extension uses explicit runtime SIMD dispatch for
byte-oriented helper kernels, starting with `seq_gc_content(...)`.
`scalar` is always available and is the portable baseline. Optional
platform backends such as `avx2` or `avx512` should be checked with
[`rduckhts_simd_backend_available()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_simd_backend.md)
before being requested. The `auto` policy resolves each logical kernel
independently from the current compiled-and-CPU-supported capability
mask; use
[`rduckhts_simd_kernel_info()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_simd_backend.md)
for the per-kernel result and `rduckhts_simd_set_backend(con, "auto")`
to return to runtime auto-detection.

``` r
rduckhts_simd_info(con)[, c("backend", "selectable", "compiled", "cpu_supported", "available", "selected")]
#>        backend selectable compiled cpu_supported available selected
#> 1       scalar       TRUE     TRUE          TRUE      TRUE    FALSE
#> 2         sse2      FALSE    FALSE          TRUE     FALSE    FALSE
#> 3        sse41      FALSE    FALSE          TRUE     FALSE    FALSE
#> 4         avx2       TRUE     TRUE          TRUE      TRUE     TRUE
#> 5       avx512       TRUE     TRUE         FALSE     FALSE    FALSE
#> 6         neon       TRUE    FALSE         FALSE     FALSE    FALSE
#> 7 wasm_simd128       TRUE    FALSE         FALSE     FALSE    FALSE

rduckhts_simd_kernel_info(con)[, c("kernel", "selected_backend", "scalar_fallback")]
#>            kernel selected_backend scalar_fallback
#> 1 seq_base_counts             avx2           FALSE
#> 2 bam_nt16_counts             avx2           FALSE
#> 3  nt16_gc_counts             avx2           FALSE

rduckhts_simd_set_backend(con, "scalar")
#> [1] "scalar"

DBI::dbGetQuery(
  con,
  paste(
    "SELECT duckhts_simd_requested_backend() AS requested_backend,",
    "duckhts_simd_backend() AS selected_backend,",
    "round(seq_gc_content('ACGTNNacgtnn'), 3) AS gc_content"
  )
)
#>   requested_backend selected_backend gc_content
#> 1            scalar           scalar        0.5

data.frame(
  requested_backend = rduckhts_simd_requested_backend(con),
  selected_backend = rduckhts_simd_backend(con)
)
#>   requested_backend selected_backend
#> 1            scalar           scalar

restored_backend <- rduckhts_simd_set_backend(con, "auto")
data.frame(
  requested_backend = rduckhts_simd_requested_backend(con),
  selected_backend_known = nzchar(restored_backend)
)
#>   requested_backend selected_backend_known
#> 1              auto                   TRUE
```

## Multi-file Reading

The `rduckhts_*_multi` family reads multiple files into a single DuckDB
table with a `filename` column, following the same
`(con, table_name, ...)` convention as the single-file wrappers:

``` r
fq_files <- c(
  system.file("extdata", "r1.fq", package = "Rduckhts"),
  system.file("extdata", "r2.fq", package = "Rduckhts")
)
rduckhts_fastq_multi(con, "fq_multi", fq_files, overwrite = TRUE)
dbGetQuery(con, "SELECT filename, count(*) AS n FROM fq_multi GROUP BY ALL")
#>                                               filename n
#> 1 /usr/local/lib/R/site-library/Rduckhts/extdata/r1.fq 5
#> 2 /usr/local/lib/R/site-library/Rduckhts/extdata/r2.fq 5
```

Per-file parameters are supported via a `.params` data.frame with a
`file` column and columns matching reader arguments. `NA` values fall
back to the uniform default:

``` r
bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
bam_idx  <- system.file("extdata", "range.bam.bai", package = "Rduckhts")

params <- data.frame(
  file       = bam_path,
  region     = "CHROMOSOME_I:1-1000",
  index_path = bam_idx
)
rduckhts_bam_multi(con, "bam_multi", bam_path, .params = params,
                   overwrite = TRUE)
dbGetQuery(con, "SELECT count(*) AS n FROM bam_multi")
#>   n
#> 1 2
```

## Function Catalog

Use
[`rduckhts_functions()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_functions.md)
inside R to inspect the generated extension catalog.

Show generated function catalog

## Extension Function Catalog

This section is generated from `functions.yaml`.

### Diagnostics

| Function                             | Kind         | Returns                | R helper                              | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|--------------------------------------|--------------|------------------------|---------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `duckhts_simd_backend`               | scalar       | VARCHAR                | `rduckhts_simd_backend`               | Return the current DuckHTS SIMD dispatch label. For explicit scalar or concrete backend requests this is the requested policy; for auto it is the single selected backend when all logical kernels resolve to the same backend, or mixed when per-kernel auto-dispatch resolves to multiple backends. Use duckhts_simd_kernel_info() for per-kernel details.                                                                                                                                                                                                |
| `duckhts_simd_requested_backend`     | scalar       | VARCHAR                | `rduckhts_simd_requested_backend`     | Return the current explicit SIMD backend request, usually auto unless `SELECT backend FROM duckhts_simd_set_backend('auto'\|'scalar'\|backend)` was called. The selected per-kernel backend may differ under auto-dispatch across x86, ARM, wasm, and scalar-only builds.                                                                                                                                                                                                                                                                                   |
| `duckhts_simd_backend_compiled`      | scalar       | BOOLEAN                | `rduckhts_simd_backend_compiled`      | Return whether a concrete DuckHTS SIMD backend was compiled into this build. This is independent of whether the current CPU/runtime supports executing that backend; for example avx512 can be compiled but not CPU-supported on the running host.                                                                                                                                                                                                                                                                                                          |
| `duckhts_simd_backend_cpu_supported` | scalar       | BOOLEAN                | `rduckhts_simd_backend_cpu_supported` | Return whether the current CPU/runtime supports a concrete DuckHTS SIMD backend, independent of whether DuckHTS compiled an implementation for it. Availability is the intersection of compiled and CPU-supported.                                                                                                                                                                                                                                                                                                                                          |
| `duckhts_simd_backend_available`     | scalar       | BOOLEAN                | `rduckhts_simd_backend_available`     | Return whether a concrete SIMD backend is usable in the current process. Availability means the backend is compiled into DuckHTS and supported by the current CPU/runtime. auto is a selection request rather than a concrete backend and is not reported as available here.                                                                                                                                                                                                                                                                                |
| `duckhts_simd_info`                  | table        | table                  | `rduckhts_simd_info`                  | Return one row per known concrete DuckHTS SIMD backend with extension-owned selectable, compiled, CPU-supported, available, selected, requested, and dispatch-mode diagnostics. Availability is the intersection of compiled and CPU/runtime-supported. selectable reports whether the backend has a selectable implementation path; explicit selection still requires available = TRUE. selected is TRUE when the current dispatch table uses that backend for at least one logical kernel. auto is a selection request and is not a concrete backend row. |
| `duckhts_simd_kernel_info`           | table        | table                  | `rduckhts_simd_kernel_info`           | Return one row per logical DuckHTS SIMD kernel showing the concrete backend selected by the current immutable dispatch table, the selected capability, the requested backend policy, whether scalar was used as a per-kernel fallback, and the dispatch mode. This is the authoritative diagnostic for mixed auto-dispatch when different kernels resolve to different backends.                                                                                                                                                                            |
| `duckhts_simd_set_backend`           | table        | table(backend VARCHAR) | `rduckhts_simd_set_backend`           | Explicitly select the DuckHTS SIMD dispatch policy for this process using a one-row table-function call and return the current dispatch label in a backend column. Use auto for per-kernel runtime dispatch or scalar for a portable baseline; unavailable platform-specific requests such as avx512 on non-AVX-512 CPUs raise an error instead of silently falling back.                                                                                                                                                                                   |
| `duckhts_duckdb_type_supported`      | scalar_macro | BOOLEAN                |                                       | Return whether the currently open DuckDB runtime advertises a logical type with the given name through duckdb_types(). This is a catalog-level runtime probe for feature gating SQL/macros across DuckDB versions.                                                                                                                                                                                                                                                                                                                                          |
| `duckhts_duckdb_supports_variant`    | scalar_macro | BOOLEAN                |                                       | Return whether the currently open DuckDB runtime advertises the VARIANT logical type. Use this to gate optional SQL that depends on DuckDB VARIANT support.                                                                                                                                                                                                                                                                                                                                                                                                 |
| `duckhts_duckdb_supports_geometry`   | scalar_macro | BOOLEAN                |                                       | Return whether the currently open DuckDB runtime advertises the GEOMETRY logical type. Use this to gate optional SQL that depends on DuckDB GEOMETRY support.                                                                                                                                                                                                                                                                                                                                                                                               |

### Variant Annotation

| Function                              | Kind        | Returns                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | R helper | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|---------------------------------------|-------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `duckvep_ensembl_regions`             | table_macro | table                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |          | Match a tiled reference sequence to one Ensembl core assembly and assign dense DuckVEP sequence-region ordinals. core_schema must contain the Ensembl seq_region and coord_system tables. reference_chunks_table must contain chrom, zero-based start, half-open end, and seq columns, normally materialized from fasta_nuc(…, include_seq := TRUE). Chunks must be contiguous from zero, their sequence lengths must match their intervals, and every FASTA contig must match exactly one same-length Ensembl region. Regions absent from the supplied FASTA are deliberately excluded.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| `duckvep_ensembl_transcripts`         | table_macro | table                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |          | Build the validated DuckVEP transcript relation from Ensembl core tables and matching tiled FASTA sequence. The VEP-116 core selection is applied before dense ordinals are assigned: transcripts must be current and have a non-empty stable ID, and artifact-biotype or readthrough_tra transcripts are excluded. The result contains the exact resident-loader columns, stable/source identifiers, biotypes, versioned RefSeq accessions from MANE Select and MANE Plus Clinical attributes when present, a nested ranked-exon projection, mature miRNA cDNA attributes projected into genomic exon segments, and supported Translation SeqEdits. MANE accessions remain cold DuckDB columns while only selection flags enter the resident C model. Prepared sequence contains the phase-adjusted CDS plus the complete transcript-oriented pre-CDS and post-CDS spliced sequence; the legacy post_cds_bases field remains as the first three post-CDS bases. Codon-table IDs come from the Ensembl seq_region_attrib relation and default to table 1 exactly as VEP 116 does; every BioPerl/VEP-supported NCBI table ID is accepted and invalid or conflicting source attributes reject the import. Single-residue initial_met, \_selenocysteine, amino_acid_sub, and \_stop_codon_rt edits are retained as a sparse reference-peptide overlay. Other Translation SeqEdit shapes and transcript-level RNA edits fail closed by withholding sequence with an explicit reason.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| `duckvep_ensembl_regulation_features` | table_macro | table                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |          | Build the VEP-116-compatible interval relation for Ensembl RegulatoryFeature and MotifFeature rows on the sequence regions selected for a DuckVEP model. funcgen_schema must contain regulatory_feature, feature_type, and motif_feature; regions_table is the output of duckvep_ensembl_regions(…). Matching VEP’s database and cache source, epigenetically_modified_region (EMAR) RegulatoryFeature rows are excluded before dense ordinals are assigned. The result maps source sequence-region IDs to model ordinals, validates one-based inclusive coordinates and source feature types, and preserves stable/source IDs, feature metadata, binding-matrix IDs, regulatory-build IDs, and motif scores as cold relation columns. feature_kind is the compact resident code: 1 regulatory region or 2 transcription-factor binding site. Pass the five hot columns regulation_feature_index, seq_region, feature_start, feature_end, and feature_kind to duckvep_model_load(interval_feature_query := …); the remaining columns stay relational for late projection.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| `duckvep_model_receipt`               | table_macro | table                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |          | Create one deterministic receipt for prepared DuckVEP region and transcript relations plus an optional prepared regulation-feature relation. The original eight positional arguments remain valid. When regulation_features_table is supplied by name, the receipt also validates its dense ordinals, region agreement, interval geometry, and feature kinds; counts regulatory regions and motif features; and includes their five resident columns in the model hash. It returns the declared provenance as source_name, source_version, assembly, source_manifest_sha256, reference_sha256, and transcript_filter; counts transcript, sequence, mature-miRNA, and peptide-edit content; and hashes every semantic resident-model field in stable ordinal order. source_manifest_sha256 should identify a canonical manifest containing every exact core and funcgen input used. The receipt contains no clock time, so rebuilding identical inputs yields the same model hash.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| `duckvep_model_load`                  | table       | table(loaded BOOLEAN)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |          | Load one immutable consequence model into the current DuckDB database instance under a caller-chosen name. Three required query strings read committed, non-temporary DuckDB relations: sorted UINTEGER sequence-region ordinals; dense transcript rows with genomic span, strand, gene ordinal, flags, optional CDS span/sequence/codon table; and transcript-ordered exon rows with genomic and cDNA spans plus phase. The optional mature_mirna_query reads transcript_index UINTEGER, mature_mirna_start UBIGINT, and mature_mirna_end UBIGINT ordered by transcript and start. The optional peptide_edit_query reads transcript_index UINTEGER, protein_position UINTEGER, and alternate_amino_acid VARCHAR ordered uniquely by transcript and protein position. The optional interval_feature_query reads regulation_feature_index UINTEGER, seq_region UINTEGER, feature_start UINTEGER, feature_end UINTEGER, and feature_kind UTINYINT ordered by region, start, and index; kind 1 is RegulatoryFeature and kind 2 is MotifFeature. Loading keeps those hot columns in a compact SoA and builds a separate cgranges seed index, while identifiers and funcgen metadata remain ordinary DuckDB columns. The transcript projection accepts 11 base columns, 12 columns with the legacy post_cds_bases BLOB of up to three bases, or 13 columns ending in complete transcript-oriented pre_cds_sequence and post_cds_sequence BLOBs. Complete flanks enable VEP start/stop predicates for length-changing edits that cross a CDS start or end. The default is a partial transcript model: a variant with no loaded transcript is unresolved, not intergenic. Setting transcript_coverage_complete := TRUE requires sequence_length UBIGINT and permits supported intergenic results after coordinate bounds checks. Loading validates and narrows the model once and returns one TRUE row. Several named models may coexist.                                                                                                            |
| `duckvep_model_drop`                  | scalar      | BOOLEAN                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |          | Remove a named resident DuckVEP consequence model and release its transcript and regulation-feature interval indexes, sequences, and cached worker state. Returns FALSE when the name is absent or the model is in use by an annotation vector.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `duckvep_annotate`                    | scalar      | LIST\<STRUCT(transcript_index UINTEGER, gene_index UINTEGER, consequence VARCHAR, impact VARCHAR, region VARCHAR, status VARCHAR, reason VARCHAR, cdna_position UBIGINT, cds_position UBIGINT, protein_position UBIGINT, reference_amino_acid VARCHAR, alternate_amino_acid VARCHAR, nmd_prediction VARCHAR, nmd_escape_intronless BOOLEAN, nmd_escape_early_cds BOOLEAN, nmd_escape_last_exon BOOLEAN, nmd_escape_penultimate_exon_end BOOLEAN, regulation_feature_index UINTEGER, overlap_object VARCHAR)\> |          | Annotate a biallelic A/C/G/T/N variant against a named resident consequence model. position is one-based and seq_region is the model’s compact ordinal. With five arguments the upstream and downstream transcript windows both default to 5000 bases, matching VEP. A sixth distance applies to both directions; a seventh argument makes the sixth upstream_distance and the seventh downstream_distance. Zero disables that direction. RegulatoryFeature and MotifFeature use exact overlap with no flank. Apply UNNEST in the SELECT list over sorted variant columns, then expand the returned struct in an outer query; this retains DuckDB vector execution. Each sorted run seeds the transcript and regulation indexes once through cgranges and advances both sweeps. Transcript rows set transcript_index and overlap_object = transcript; regulation rows set regulation_feature_index and overlap_object = regulatory_region or transcription_factor_binding_site. Join either ordinal to its cold prepared relation for stable IDs and metadata. A no-transcript-overlap input is supported intergenic only for a model loaded with transcript_coverage_complete := TRUE; partial models return no_feature_in_loaded_model as unresolved even when separate regulation rows are present. Coding rows report stable reasons including missing_sequence, missing_transcript_tail, missing_transcript_flank, ambiguous_sequence, reference_mismatch, non_contiguous_cds_edit, unsupported_compound_consequence, invalid_model_projection, and internal_capacity_error. Amino-acid fields are NULL when the compact row has no scalar one-letter amino-acid value, independently of protein_position. For stop-gained, frameshift, splice-donor, and splice-acceptor consequences on coding transcripts, nmd_prediction applies the pinned VEP Plugins release/116 NMD.pm policy. These fields are separate from the NMD_transcript_variant SO term, which describes an already curated nonsense_mediated_decay transcript biotype. |
| `duckvep_annotate_compact`            | scalar      | LIST\<STRUCT(transcript_index UINTEGER, gene_index UINTEGER, consequence_mask UBIGINT, region_mask UINTEGER, impact_code UTINYINT, status_code UTINYINT, reason_code UTINYINT, cdna_position UINTEGER, cds_position UINTEGER, protein_position UINTEGER, reference_amino_acid_code UTINYINT, alternate_amino_acid_code UTINYINT, nmd_prediction_code UTINYINT, nmd_escape_reasons UTINYINT, regulation_feature_index UINTEGER, overlap_object_code UTINYINT)\>                                                |          | Run the same sorted transcript and regulation-feature sweeps and the same optional transcript-distance windows as duckvep_annotate(…), but keep the result numeric so filtering, aggregation, and joins happen before text rendering. overlap_object_code is 0 transcript, 1 regulatory region, or 2 transcription-factor binding site; exactly one of transcript_index and regulation_feature_index is present on a biological feature row. consequence_mask uses the stable bit assignments in src/duckvep/kernel/include/duckvep_so.h. region_mask uses 1 upstream, 2 downstream, 4 intron, 8 exon, 16 CDS, 32 UTR, and 64 splice. impact_code is 0 modifier, 1 low, 2 moderate, or 3 high; status_code is 0 supported or 1 unresolved. reason_code is 0 none, 1 no feature in a partial model, 2 missing sequence, 3 ambiguous sequence, 4 reference mismatch, 5 non-contiguous CDS edit, 6 unsupported compound consequence, 7 invalid model projection, 8 internal capacity, 9 missing transcript tail, or 10 missing transcript flank. Amino-acid codes are ASCII bytes and remain NULL when no scalar amino acid exists. nmd_prediction_code is 0 not applicable, 1 triggering, 2 escaping, or 3 unresolved; nmd_escape_reasons is a bit mask with 1 intronless, 2 early CDS, 4 last exon, and 8 penultimate-exon end. Coordinate fields are NULL when unavailable. The numeric codes are stable SQL output contracts, not model ordinals.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| `duckvep_annotate_sv`                 | scalar      | LIST\<STRUCT(transcript_index UINTEGER, gene_index UINTEGER, consequence VARCHAR, impact VARCHAR, region VARCHAR, status VARCHAR, reason VARCHAR, cdna_position UBIGINT, cds_position UBIGINT, protein_position UBIGINT, reference_amino_acid VARCHAR, alternate_amino_acid VARCHAR, nmd_prediction VARCHAR, nmd_escape_intronless BOOLEAN, nmd_escape_early_cds BOOLEAN, nmd_escape_last_exon BOOLEAN, nmd_escape_penultimate_exon_end BOOLEAN, regulation_feature_index UINTEGER, overlap_object VARCHAR)\> |          | Annotate an exact single-locus structural event against a named resident consequence model. DEL, DUP, TDUP, INV, CNV, and UNKNOWN use a one-based inclusive start/end span. INS uses start = end = P for the interbase site after reference base P; symbolic VCF preparation therefore maps a span event to start = POS + 1 and end = INFO/END, but maps an insertion to P = POS after removing its left anchor. structural_type is DEL, DUP, TDUP, INV, INS, CNV, or UNKNOWN; copy_change is UNKNOWN, LOSS, NEUTRAL, or GAIN. Contradictory pairs such as DEL with GAIN fail instead of silently choosing one meaning. The optional distance arguments have the same 5000-base default, symmetric-seventh-argument form, separate upstream/downstream form, and zero-to-disable behavior as duckvep_annotate(…). Input must be sorted by model, seq_region, and start to retain the forward transcript and regulation-feature sweeps. BND is rejected because it requires both loci through duckvep_annotate_breakend(…). This surface does not accept imprecise confidence intervals or unexpanded STR repeat units/counts.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `duckvep_annotate_sv_compact`         | scalar      | LIST\<STRUCT(transcript_index UINTEGER, gene_index UINTEGER, consequence_mask UBIGINT, region_mask UINTEGER, impact_code UTINYINT, status_code UTINYINT, reason_code UTINYINT, cdna_position UINTEGER, cds_position UINTEGER, protein_position UINTEGER, reference_amino_acid_code UTINYINT, alternate_amino_acid_code UTINYINT, nmd_prediction_code UTINYINT, nmd_escape_reasons UTINYINT, regulation_feature_index UINTEGER, overlap_object_code UTINYINT)\>                                                |          | Run the same typed exact single-locus structural-event path, insertion-site convention, and distance controls as duckvep_annotate_sv(…), returning the same numeric consequence, region, impact, status, reason, coordinate, amino-acid, and NMD codes as duckvep_annotate_compact(…). This is the join-oriented form for filtering and aggregation before text rendering. It rejects BND, imprecise structural coordinates, and unexpanded STR repeat units/counts.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| `duckvep_annotate_breakend`           | scalar      | LIST\<STRUCT(transcript_index UINTEGER, gene_index UINTEGER, consequence VARCHAR, impact VARCHAR, region VARCHAR, status VARCHAR, reason VARCHAR, cdna_position UBIGINT, cds_position UBIGINT, protein_position UBIGINT, reference_amino_acid VARCHAR, alternate_amino_acid VARCHAR, nmd_prediction VARCHAR, nmd_escape_intronless BOOLEAN, nmd_escape_early_cds BOOLEAN, nmd_escape_last_exon BOOLEAN, nmd_escape_penultimate_exon_end BOOLEAN, regulation_feature_index UINTEGER, overlap_object VARCHAR)\> |          | Annotate one exact paired BND against the transcript portion of a named resident consequence model. Both positions are raw one-based VCF endpoint coordinates. To reproduce VEP 116, ordinary topology uses local_position + 1 while the mate position is retained verbatim; both loci discover candidates, and the function returns the union of VEP’s endpoint-overlap consequence sets once per transcript. region is NULL when only the mate contributes a transcript consequence. VEP’s fixed 5000-base endpoint-to-transcript admission cap applies in addition to the same optional directional distances and zero-to-disable behavior as duckvep_annotate(…). Input must be sorted by model, local_seq_region, and local_position. This path emits transcript consequences only: regulation_feature_index is NULL and overlap_object is transcript even when the loaded model also contains regulatory or motif intervals. Paired-BND regulatory and motif consequences are not part of this contract. Keep the raw ALT, bracket orientation, event identifier, and provenance in the surrounding relation for HGVS, fusion, and round-trip consumers; this function does not parse them and does not approximate the pair as a span.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `duckvep_annotate_breakend_compact`   | scalar      | LIST\<STRUCT(transcript_index UINTEGER, gene_index UINTEGER, consequence_mask UBIGINT, region_mask UINTEGER, impact_code UTINYINT, status_code UTINYINT, reason_code UTINYINT, cdna_position UINTEGER, cds_position UINTEGER, protein_position UINTEGER, reference_amino_acid_code UTINYINT, alternate_amino_acid_code UTINYINT, nmd_prediction_code UTINYINT, nmd_escape_reasons UTINYINT, regulation_feature_index UINTEGER, overlap_object_code UTINYINT)\>                                                |          | Run the same paired two-locus BND candidate search, VEP-116 coordinate rules, transcript-level consequence union, and distance controls as duckvep_annotate_breakend(…), but return the numeric consequence contract used by duckvep_annotate_compact(…). A zero region_mask means that the transcript was affected only through the mate and has no local topology. This path emits transcript consequences only: regulation_feature_index is NULL and overlap_object_code is 0 even when the loaded model also contains regulatory or motif intervals. Paired-BND regulatory and motif consequences are not part of this contract. Keep event identity, raw ALT, orientation, and provenance as ordinary columns beside this join-oriented result.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |

### Readers

| Function                 | Kind         | Returns | R helper                                                                                                                                                               | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
|--------------------------|--------------|---------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `read_bcf`               | table        | table   | `rduckhts_bcf`                                                                                                                                                         | Read VCF and BCF variant data with typed INFO, FORMAT, typed CSQ/ANN/BCSQ subfields, optional tidy sample output, optional bcftools-style CSQ type overrides, explicit scan_mode control (‘auto’ or full-file ‘sequential’), optional htslib decompression worker threads via decompression_threads (default 0 for single-threaded reads), and decode_error_policy (‘null’, ‘warn’, or ‘error’) for corrupt BCF header-vs-payload decode mismatches. Comma-separated indexed regions use htslib’s native multi-region iterator and emit a record once when requested regions overlap.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| `read_bcf_v2`            | table        | table   |                                                                                                                                                                        | Experimental clone of read_bcf(…) for non-destructive benchmarking of current-schema scan optimizations. It preserves the read_bcf wide/tidy output schema by default, streams full-file scans by default instead of auto-iterating indexed contigs, accepts the same scan_mode and decode_error_policy (‘null’, ‘warn’, or ‘error’) contract as read_bcf, and shares read_bcf’s native deduplicating multi-region iterator. It sets htslib VCF text max_unpack from DuckDB projection pushdown so unprojected INFO/FORMAT text can be skipped earlier, keeps INFO/FORMAT/GT decode buffers and parsed VEP/CSQ records in local scan state across output chunks, repeats VEP-derived values correctly across tidy sample rows, supports v2-only htslib sample pushdown via comma-separated samples/exclude_samples strings or samples_file/exclude_samples_file lists, can bind only requested INFO/FORMAT/VEP fields through include_info/include_format/include_vep and comma-separated info_fields/format_fields/vep_fields controls, rejects format_fields when no samples are selected, treats tidy_format := true with include_format := false as one row per variant because no SAMPLE_ID/FORMAT columns are bound, and counts VCF text rows without parsing when no columns are projected.                                                                                                                                                                                                                                                                   |
| `read_bcf_appender`      | table        | table   |                                                                                                                                                                        | Experimental side-effecting benchmark helper that reads BCF/VCF with htslib and appends a narrow scalar tidy table through DuckDB’s chunk appender API. The target table is created in the default schema with columns CHROM, POS, REF, ALT, SAMPLE_ID, and FORMAT_GT; set include_file_offset := TRUE to add FILE_OFFSET as a file-order token. Comma-separated regions use htslib’s native deduplicating multi-region iterators. region_threads controls at most 16 worker handles claiming primary-contig jobs iteratively: every contig keeps its complete requested interval set in one native iterator, so repeated, overlapping, or record-spanning intervals produce the same target schema and row multiset for every thread count. Header contig count does not create worker threads or open handles, and a request limited to one contig uses one iterator. Parallel row arrival is unordered; use ORDER BY FILE_OFFSET when file order matters. decompression_threads configures htslib workers per open region-worker handle, so total helper threads can scale as region workers times decompression workers. The operation runs inside an internal transaction and surfaces malformed-record errors as DuckDB errors; decode_error_policy controls corrupt BCF FORMAT/GT header-vs-payload mismatches with ‘null’, ‘warn’, or ‘error’. target_table must be an unqualified SQL identifier. This is intended to compare appender-based materialization against read_bcf(…) query pipelines, not as a full replacement for the typed read_bcf surface. |
| `read_bam`               | table        | table   | `rduckhts_bam`                                                                                                                                                         | Read SAM, BAM, and CRAM alignments with optional typed SAMtags and auxiliary tag maps. Use sequence_encoding := ‘nt16’ to return SEQ as UTINYINT\[\], quality_representation := ‘phred’ to return QUAL as UTINYINT\[\], and cigar_representation := ‘binary’ to return packed BAM CIGAR operations as UINTEGER\[\] instead of SAM text. scan_mode := ‘sequential’ forces a full-file streaming scan instead of index-backed count/parallel paths and is incompatible with region. decompression_threads controls per-file htslib worker threads and defaults to 2; use 0 to disable worker threads.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `read_fasta`             | table        | table   | `rduckhts_fasta`                                                                                                                                                       | Read FASTA records or indexed FASTA regions as sequence rows. Use sequence_encoding := ‘nt16’ to return SEQUENCE as UTINYINT\[\] (htslib nt16 4-bit codes) instead of VARCHAR. For bgzipped FASTA, gzi_path may point to an explicit .gzi sidecar when it is not colocated with the FASTA. scan_mode := ‘sequential’ forces full-file streaming/counting instead of index-backed count paths and is incompatible with region.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| `read_bed`               | table        | table   | `rduckhts_bed`                                                                                                                                                         | Read BED3-BED12 interval files with canonical typed columns and optional tabix-backed region filtering. scan_mode := ‘sequential’ forces full-file streaming/counting instead of index-backed count paths and is incompatible with region.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| `fasta_nuc`              | table        | table   | `rduckhts_fasta_nuc`                                                                                                                                                   | Compute bedtools nuc-style nucleotide composition for supplied BED intervals or generated fixed-width bins over a FASTA reference. For bgzipped FASTA, gzi_path may point to an explicit .gzi sidecar when it is not colocated with the FASTA.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `read_fastq`             | table        | table   | `rduckhts_fastq`                                                                                                                                                       | Read single-end, paired-end, or interleaved FASTQ files with optional legacy quality decoding. By default, FASTQ qualities are interpreted as modern Phred+33 input. Use sequence_encoding := ‘nt16’ to return SEQUENCE as UTINYINT\[\] and quality_representation := ‘phred’ to return QUALITY as UTINYINT\[\] instead of VARCHAR. input_quality_encoding accepts ‘phred33’, ‘auto’, ‘phred64’, or ‘solexa64’. scan_mode := ‘sequential’ forces raw streaming/counting instead of index-backed count paths.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| `read_gff`               | table        | table   | `rduckhts_gff`                                                                                                                                                         | Read GFF annotations with optional raw scalar and richer list/pair parsed attribute columns, strict GFF3 structural validation, and indexed region filtering. Comma-separated indexed regions use htslib’s native multi-region iterator and emit a row once when requested regions overlap. scan_mode := ‘sequential’ forces full-file streaming/counting instead of index-backed count paths and is incompatible with region.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `read_gtf`               | table        | table   | `rduckhts_gtf`                                                                                                                                                         | Read GTF annotations with optional raw scalar and richer list/pair parsed attribute columns and indexed region filtering. Comma-separated indexed regions use htslib’s native multi-region iterator and emit a row once when requested regions overlap. scan_mode := ‘sequential’ forces full-file streaming/counting instead of index-backed count paths and is incompatible with region.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| `read_tabix`             | table        | table   | `rduckhts_tabix`                                                                                                                                                       | Read generic tabix-indexed text data with optional header handling and type inference. Comma-separated indexed regions use htslib’s native multi-region iterator and emit a row once when requested regions overlap. scan_mode := ‘sequential’ forces full-file streaming/counting instead of index-backed count paths and is incompatible with region.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| `fasta_index`            | table        | table   | `rduckhts_fasta_index`                                                                                                                                                 | Build a FASTA index (.fai) and return a single row with columns success (BOOLEAN) and index_path (VARCHAR).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| `hts_union_query`        | scalar_macro | VARCHAR | `rduckhts_bam_multi, rduckhts_bcf_multi, rduckhts_fastq_multi, rduckhts_fasta_multi, rduckhts_bed_multi, rduckhts_tabix_multi, rduckhts_gff_multi, rduckhts_gtf_multi` | Generate a UNION ALL BY NAME query string that reads every file matching a glob pattern through the named reader function. The result includes a ‘filename’ column identifying the source file for each row. Assign to a variable with SET VARIABLE and execute via query(getvariable(…)). Optional params string is appended to each reader call. In R, use the typed rduckhts\_\*\_multi() helpers instead, which accept file vectors with optional per-file parameters and create DuckDB tables directly.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| `hts_region_union_query` | scalar_macro | VARCHAR |                                                                                                                                                                        | Generate a UNION ALL BY NAME query string that reads one HTS file through separate per-region table-function scans. The result adds filename, duckhts_region_shard_id, and duckhts_region_shard columns so benchmark queries can compare region-sharded execution and explicit ordered writes. Pass regions as a DuckDB list of region strings; params is appended to each reader call and should not include region. The generated UNION ALL does not deduplicate records; adjacent or overlapping shards can duplicate BAM alignments or BCF/VCF records that span shard boundaries, so add explicit shard-local filters or downstream deduplication when exact once-only semantics are required.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |

### Converters

| Function                            | Kind         | Returns | R helper                         | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
|-------------------------------------|--------------|---------|----------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `duckhts_bcf_convert_parquet_sql`   | scalar_macro | VARCHAR | `rduckhts_bcf_convert_parquet`   | Build a DuckDB COPY statement that converts read_bcf(…) output to Parquet with DuckHTS key-value metadata, preserved or corrected VCF header text, user metadata from metadata := map(…), selected columns, filters, and partition columns. Optional metadata_json_file is caller-managed and requires DuckDB’s json extension to be available when the builder is called; use metadata maps for CRAN/offline-safe workflows. Clients execute the returned SQL string; the R wrapper is a thin DBI helper that does this execution.       |
| `duckhts_bam_convert_parquet_sql`   | scalar_macro | VARCHAR | `rduckhts_bam_convert_parquet`   | Build a DuckDB COPY statement that converts read_bam(…) output to Parquet with DuckHTS key-value metadata, preserved or corrected SAM header text, user metadata from metadata := map(…), selected columns, filters, and partition columns. Optional metadata_json_file is caller-managed and requires DuckDB’s json extension to be available when the builder is called; use metadata maps for CRAN/offline-safe workflows. Clients execute the returned SQL string; the R wrapper is a thin DBI helper that does this execution.       |
| `duckhts_gff_convert_parquet_sql`   | scalar_macro | VARCHAR | `rduckhts_gff_convert_parquet`   | Build a DuckDB COPY statement that converts read_gff(…) output to Parquet with DuckHTS key-value metadata, preserved or corrected GFF/tabix header text, user metadata from metadata := map(…), selected columns, filters, and partition columns. Optional metadata_json_file is caller-managed and requires DuckDB’s json extension to be available when the builder is called; use metadata maps for CRAN/offline-safe workflows. Clients execute the returned SQL string; the R wrapper is a thin DBI helper that does this execution. |
| `duckhts_tabix_convert_parquet_sql` | scalar_macro | VARCHAR | `rduckhts_tabix_convert_parquet` | Build a DuckDB COPY statement that converts read_tabix(…) output to Parquet with DuckHTS key-value metadata, preserved or corrected tabix header text, user metadata from metadata := map(…), selected columns, filters, and partition columns. Optional metadata_json_file is caller-managed and requires DuckDB’s json extension to be available when the builder is called; use metadata maps for CRAN/offline-safe workflows. Clients execute the returned SQL string; the R wrapper is a thin DBI helper that does this execution.   |

### Coverage

| Function                   | Kind  | Returns | R helper                    | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
|----------------------------|-------|---------|-----------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `read_pileup`              | table | table   | `rduckhts_pileup`           | Construct a region-scoped BAM pileup with one row per covered position, emitting chrom, 1-based position, depth, observed bases, and Phred+33 qualities after SAM flag and MAPQ filtering. This is a compact htslib pileup view, not samtools mpileup text parity.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `bam_bin_counts`           | table | table   | `rduckhts_bam_bin_counts`   | Count BAM or CRAM read starts into fixed-width bins. Returns one row per bin across the selected contig span, including zero-count bins, with total, forward, and reverse counts; `rmdup := 'streaming'` applies the WisecondorX-style larp/larp2 consecutive-position deduplication, `rmdup := 'flag'` drops SAM duplicate-flagged reads, and `stats := 'gc'`, `'mq'`, or `'gc,mq'` adds per-bin pre/post-filter GC and MAPQ sufficient statistics, including reference GC when `reference` is provided.                                                                                                                                                                                                                                                                                                                                                                                                                |
| `duckhts_bam_bed_coverage` | table | table   | `rduckhts_bam_bed_coverage` | Compute samtools coverage-like regional summaries for BAM or CRAM input over a BED target set, returning one row per BED interval with DuckHTS-specific pre/post-filter read counts, covered bases, percentage covered, mean depth, mean baseQ, mean mapQ, and strand-specific post-filter summaries in read mode. Indexed BAM/CRAM input is required in the current implementation. decompression_threads controls htslib worker threads for BAM/CRAM decoding; use 0 to disable them.                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `duckhts_mosdepth`         | table | table   | `rduckhts_mosdepth`         | Write native mosdepth-compatible coverage outputs for indexed BAM or CRAM input. Produces mosdepth-style summary, global distribution, per-base BED.gz + CSI, optional window/BED region outputs, optional quantized BED.gz + CSI, and optional threshold counts for `by`; `fast_mode` defaults to FALSE to match upstream mosdepth, default mode performs CIGAR-aware coverage with mate-overlap correction, `fragment_mode` switches coverage to full-fragment insert spans for proper pairs, `use_median` switches `by` outputs from mean to median, `read_groups` filters by comma-separated RG IDs, `min_frag_len` and `max_frag_len` filter on absolute template length, `fasta` is required for CRAM when htslib needs a reference, `precision_digits` controls decimal places in the text outputs, and `processing_threads` enables parallel contig processing (0 = sequential, \>0 = number of worker threads). |

### Intervals

| Function                           | Kind   | Returns                                                                                                                                      | R helper | Description                                                                                                                                                                                                                                                                                                                                                                                                                  |
|------------------------------------|--------|----------------------------------------------------------------------------------------------------------------------------------------------|----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `duckhts_cgranges_create`          | scalar | BOOLEAN                                                                                                                                      |          | Create an empty session-scoped cgranges registry entry that can be populated with intervals and finalized for overlap queries.                                                                                                                                                                                                                                                                                               |
| `duckhts_cgranges_add`             | scalar | BOOLEAN                                                                                                                                      |          | Append an interval to a session-scoped cgranges registry entry before finalization. Labels may be BIGINT-like, DOUBLE, VARCHAR, or BOOLEAN.                                                                                                                                                                                                                                                                                  |
| `duckhts_cgranges_index`           | scalar | BOOLEAN                                                                                                                                      |          | Finalize a populated cgranges registry entry and build its immutable overlap index for subsequent queries.                                                                                                                                                                                                                                                                                                                   |
| `duckhts_cgranges_destroy`         | scalar | BOOLEAN                                                                                                                                      |          | Destroy a session-scoped cgranges registry entry and release its indexed interval storage when it is not in active use.                                                                                                                                                                                                                                                                                                      |
| `duckhts_cgranges_from_query`      | scalar | BOOLEAN                                                                                                                                      |          | Execute a SQL query on an extension-owned DuckDB connection, append its interval rows into a session-scoped cgranges registry entry, and leave the populated index ready for explicit finalization with duckhts_cgranges_index(…).                                                                                                                                                                                           |
| `duckhts_cgranges_from_table`      | scalar | BOOLEAN                                                                                                                                      |          | Reserved convenience constructor for bulk cgranges population from a table name. The current implementation is intentionally deferred and directs callers to duckhts_cgranges_from_query(…).                                                                                                                                                                                                                                 |
| `duckhts_cgranges_has_overlap`     | scalar | BOOLEAN                                                                                                                                      |          | Vectorized scalar predicate for streaming provider rows through a finalized session-scoped cgranges index. Returns TRUE when the query interval overlaps at least one indexed interval, or when mode = ‘contain’ and it fully contains at least one indexed interval; NULL inputs return NULL.                                                                                                                               |
| `duckhts_cgranges_count_overlaps`  | scalar | BIGINT                                                                                                                                       |          | Vectorized scalar overlap counter for streaming provider rows through a finalized session-scoped cgranges index. Returns the number of indexed intervals that overlap the query interval, or with mode = ‘contain’ the number fully contained by it; NULL inputs return NULL.                                                                                                                                                |
| `duckhts_cgranges_overlaps_list`   | scalar | STRUCT(interval_ordinal BIGINT, label VARCHAR, label_type VARCHAR, interval_chrom VARCHAR, interval_start INTEGER, interval_end INTEGER)\[\] |          | Vectorized scalar overlap expander for streaming provider rows through a finalized session-scoped cgranges index. Returns a LIST of hit STRUCTs that can be expanded with UNNEST, preserving provider columns while emitting one row per matching indexed interval. Because scalar return types are fixed, labels are returned as text with label_type describing the original cgranges label kind; NULL inputs return NULL. |
| `duckhts_cgranges_overlaps`        | table  | table                                                                                                                                        |          | Query a finalized session-scoped cgranges registry entry and return one row per overlapping or containing indexed interval, preserving the original label type and interval coordinates.                                                                                                                                                                                                                                     |
| `duckhts_cgranges_overlaps_bulk`   | table  | table                                                                                                                                        |          | Run a SQL query that yields overlap probes, stream those rows through a finalized session-scoped cgranges registry entry, and return one row per matching indexed interval. The probe query runs on the extension-owned helper connection, so it must reference regular tables/views rather than connection-local temp tables. When query_row_id_col is omitted, query_row_id defaults to the 1-based probe row ordinal.     |
| `regionkey`                        | scalar | UBIGINT                                                                                                                                      |          | Encode a genomic interval as an official RegionKey-compatible 64-bit unsigned integer. Start and end use 0-based half-open interval semantics, matching BED-style coordinates; strand accepts -1, 0, or 1.                                                                                                                                                                                                                   |
| `regionkey_hex`                    | scalar | VARCHAR                                                                                                                                      |          | Render a RegionKey as its lowercase 16-character hexadecimal string representation.                                                                                                                                                                                                                                                                                                                                          |
| `parse_regionkey_hex`              | scalar | UBIGINT                                                                                                                                      |          | Parse a 16-character hexadecimal RegionKey string back into its UBIGINT code. Invalid or non-hex strings return NULL.                                                                                                                                                                                                                                                                                                        |
| `encode_regionkey`                 | scalar | UBIGINT                                                                                                                                      |          | Encode the raw upstream RegionKey fields directly: chromosome code, 0-based start, 0-based end, and strand code (0 = unknown, 1 = +, 2 = -).                                                                                                                                                                                                                                                                                 |
| `extract_regionkey_chrom`          | scalar | UTINYINT                                                                                                                                     |          | Extract the raw upstream RegionKey chromosome code.                                                                                                                                                                                                                                                                                                                                                                          |
| `extract_regionkey_startpos`       | scalar | UINTEGER                                                                                                                                     |          | Extract the raw upstream RegionKey 0-based start position.                                                                                                                                                                                                                                                                                                                                                                   |
| `extract_regionkey_endpos`         | scalar | UINTEGER                                                                                                                                     |          | Extract the raw upstream RegionKey 0-based end position.                                                                                                                                                                                                                                                                                                                                                                     |
| `extract_regionkey_strand`         | scalar | UTINYINT                                                                                                                                     |          | Extract the raw upstream RegionKey strand code (0 = unknown, 1 = +, 2 = -).                                                                                                                                                                                                                                                                                                                                                  |
| `decode_regionkey`                 | scalar | STRUCT                                                                                                                                       |          | Decode a RegionKey into its raw upstream numeric fields: chrom_code, start, end, and strand_code.                                                                                                                                                                                                                                                                                                                            |
| `reverse_regionkey`                | scalar | STRUCT                                                                                                                                       |          | Decode a RegionKey into a STRUCT with chrom, chrom_code, start, end, strand, and strand_code.                                                                                                                                                                                                                                                                                                                                |
| `extend_regionkey`                 | scalar | UBIGINT                                                                                                                                      |          | Extend a RegionKey interval by a fixed number of bases on both sides, clamping to the official 28-bit RegionKey position range.                                                                                                                                                                                                                                                                                              |
| `are_overlapping_regions`          | scalar | BOOLEAN                                                                                                                                      |          | Return TRUE when two explicit 0-based half-open intervals overlap on the same canonical chromosome.                                                                                                                                                                                                                                                                                                                          |
| `are_overlapping_region_regionkey` | scalar | BOOLEAN                                                                                                                                      |          | Return TRUE when a 0-based half-open interval overlaps the supplied RegionKey interval.                                                                                                                                                                                                                                                                                                                                      |
| `are_overlapping_regionkeys`       | scalar | BOOLEAN                                                                                                                                      |          | Return TRUE when two RegionKeys overlap.                                                                                                                                                                                                                                                                                                                                                                                     |

### Metadata

| Function                    | Kind        | Returns | R helper                           | Description                                                                                                                                                                                                                                                                |
|-----------------------------|-------------|---------|------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `detect_quality_encoding`   | table       | table   | `rduckhts_detect_quality_encoding` | Inspect a FASTQ file’s observed quality ASCII range and report compatible legacy encodings with a heuristic guessed encoding.                                                                                                                                              |
| `duckhts_samtools_idxstats` | table       | table   | `rduckhts_samtools_idxstats`       | Write samtools idxstats-compatible TAB-delimited output for BAM, CRAM, or SAM input. Indexed BAM uses `hts_idx_get_stat(...)` for the fast path; CRAM, SAM, and unindexed BAM fall back to a full scan while preserving samtools-style contig rows plus the final `*` row. |
| `read_hts_header`           | table       | table   | `rduckhts_hts_header`              | Inspect HTS headers in parsed, raw, or combined form across supported formats. Raw VCF/BCF mode includes the final `#CHROM` sample header line so the returned text is suitable for Parquet metadata and future VCF/BCF regeneration.                                      |
| `read_hts_index`            | table       | table   | `rduckhts_hts_index`               | Inspect high-level HTS index metadata such as sequence names and mapped counts.                                                                                                                                                                                            |
| `read_hts_index_spans`      | table       | table   | `rduckhts_hts_index_spans`         | Expand index metadata into span and chunk rows suitable for low-level index inspection.                                                                                                                                                                                    |
| `read_hts_index_raw`        | table_macro | table   | `rduckhts_hts_index_raw`           | Return the raw on-disk HTS index blob together with basic identifying metadata.                                                                                                                                                                                            |

### Compression

| Function  | Kind  | Returns | R helper           | Description                                                                           |
|-----------|-------|---------|--------------------|---------------------------------------------------------------------------------------|
| `bgzip`   | table | table   | `rduckhts_bgzip`   | Compress a plain file to BGZF and return the created output path and byte counts.     |
| `bgunzip` | table | table   | `rduckhts_bgunzip` | Decompress a BGZF-compressed file and return the created output path and byte counts. |

### Indexing

| Function      | Kind  | Returns | R helper               | Description                                                                                        |
|---------------|-------|---------|------------------------|----------------------------------------------------------------------------------------------------|
| `bam_index`   | table | table   | `rduckhts_bam_index`   | Build a BAM or CRAM index and report the written index path and format.                            |
| `bcf_index`   | table | table   | `rduckhts_bcf_index`   | Build a TBI or CSI index for a VCF or BCF file and report the written index path and format.       |
| `tabix_index` | table | table   | `rduckhts_tabix_index` | Build a tabix index for a BGZF-compressed text file using a preset or explicit coordinate columns. |

### Variants

| Function                    | Kind        | Returns  | R helper                 | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
|-----------------------------|-------------|----------|--------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `variantkey`                | scalar      | UBIGINT  |                          | Encode a normalized biallelic variant as an official VariantKey-compatible 64-bit unsigned integer. This DuckHTS wrapper accepts 1-based VCF/DuckHTS POS to match bcftools `%VKX` / `+add-variantkey`, internally converts to the upstream 0-based field, and preserves the official hashed nonreversible mode for large, ambiguous, and symbolic REF/ALT strings. Only CHROM, POS, REF, and ALT are encoded; END, SVLEN, mate breakend coordinates, and other SV metadata are not.                                                                                                                                                                                                                                                                                                                                                                                          |
| `variantkey_hex`            | scalar      | VARCHAR  |                          | Render a VariantKey as its lowercase 16-character hexadecimal string representation.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| `parse_variantkey_hex`      | scalar      | UBIGINT  |                          | Parse a 16-character hexadecimal VariantKey string back into its UBIGINT code. Invalid or non-hex strings return NULL.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `encode_variantkey`         | scalar      | UBIGINT  |                          | Encode the raw upstream VariantKey fields directly: chromosome code, 0-based position, and 31-bit REF+ALT code.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| `extract_variantkey_chrom`  | scalar      | UTINYINT |                          | Extract the raw upstream VariantKey chromosome code.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| `extract_variantkey_pos`    | scalar      | UINTEGER |                          | Extract the raw upstream VariantKey 0-based position field.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `extract_variantkey_refalt` | scalar      | UINTEGER |                          | Extract the raw upstream 31-bit VariantKey REF+ALT code.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| `decode_variantkey`         | scalar      | STRUCT   |                          | Decode a VariantKey into its raw upstream numeric fields: chrom_code, pos0, and refalt_code.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `reverse_variantkey`        | scalar      | STRUCT   |                          | Decode a VariantKey into a STRUCT with chrom, chrom_code, 1-based pos, upstream 0-based pos0, ref, alt, refalt_code, and reversible. For hashed nonreversible keys, reversible is FALSE and ref/alt are returned as NULL because DuckHTS v1 does not ship the optional NRVK lookup sidecar.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `variantkey_range`          | scalar      | STRUCT   |                          | Return the inclusive minimum and maximum VariantKey bounds for a chromosome plus 1-based VCF position range, suitable for numeric range filtering on precomputed VariantKeys.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| `bcftools_liftover`         | scalar      | STRUCT   | `rduckhts_liftover`      | Row-oriented liftover kernel intended to mirror bcftools +liftover semantics as closely as possible while returning one STRUCT per input row with fields: src_chrom, src_pos, src_ref, src_alt, dest_chrom, dest_pos, dest_end, dest_ref, dest_alt, mapped, reverse_complemented, swap, reject_reason, and note. Set no_left_align := true to skip post-liftover left-alignment of lifted indels (mirrors –no-left-align in bcftools +liftover).                                                                                                                                                                                                                                                                                                                                                                                                                             |
| `duckdb_liftover`           | table_macro | table    | `rduckhts_liftover`      | DuckDB-specific wrapper over bcftools_liftover that takes either a table name or a derived-table expression plus column-name strings for chrom/pos/ref/alt and returns the lifted table. The no_left_align parameter mirrors –no-left-align in bcftools +liftover.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| `bcftools_norm_row`         | scalar      | STRUCT   |                          | Normalize one variant row with bcftools/vt-style left-alignment semantics against a FASTA reference. The alt argument may be either a comma-delimited VARCHAR or a VARCHAR\[\] list. The returned STRUCT contains pos_normed, end_pos_normed, ref_normed, alt_normed (always VARCHAR\[\]), normed (TRUE/FALSE/NULL), and norm_status. Symbolic ~~rows can use end_pos, and symbolic rows can use svlen. gVCF /\<\*\> reference-block alleles pass through with GVCFReferenceBlock, and mixed real-plus-gVCF-symbolic alleles normalize the real alleles while preserving symbolic alleles and caller-supplied reference-block END in site-preserving output.~~                                                                                                                                                                                                               |
| `duckhts_bcftools_norm`     | table_macro | table    | `rduckhts_bcftools_norm` | DuckDB table macro wrapper over bcftools_norm_row that normalizes variants from a table or derived-table expression while preserving the original columns. The input ALT column may be either VARCHAR or VARCHAR\[\]. The result appends pos_normed, end_pos_normed, ref_normed, alt_normed, normed, and norm_status; with split_multiallelic := TRUE, multiallelic sites are split before normalization and alt_normed becomes VARCHAR plus alt_index. This is a vt/vcfnorm-style row/table transform, not a full-record VCF/BCF rewrite: genotype, PL/GP/DS, and PS fields are preserved as caller columns unless a separate full-record writer/remapper layer is used.                                                                                                                                                                                                    |
| `bcftools_score`            | table       | table    | `rduckhts_score`         | Compute polygenic scores from one genotype BCF/VCF and one or more summary-statistics files with bcftools +score-compatible GT/DS/HDS/AP/GP/AS dosage semantics, sample subsetting, and region/target/FILTER-string controls. The second argument accepts a scalar path or a DuckDB LIST/array of paths; TSV/SSF summaries produce one PRS column per file in a single genotype scan, while GWAS-VCF summaries still produce one PRS column per FORMAT sample. Use summaries_list_file with a NULL second argument to read paths from a file or directory; list-file entries are interpreted as written, matching upstream `bcftools +score --summaries` behavior, while directory inputs scan supported regular summary files in lexicographic order and ignore index sidecars. Use log_path to write per-PRS loaded/matched/allele-mismatch/duplicate-marker audit counts. |
| `bcftools_munge_row`        | scalar      | STRUCT   |                          | Normalize one summary-statistics row into GWAS-VCF-style fields (chrom/pos/ref/alt/effect metrics), resolving REF/ALT orientation against a FASTA reference and applying swap-aware sign/frequency/count transforms. The output flag `alleles_swapped` means REF/ALT orientation was swapped to match the FASTA reference.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `duckdb_munge`              | table_macro | table    | `rduckhts_munge`         | DuckDB macro wrapper over bcftools_munge_row that maps source columns (via preset or explicit map) and returns normalized GWAS-VCF-style rows with lean outputs and explicit `alleles_swapped` semantics. Output columns: chrom, pos, id, ref, alt, alleles_swapped, filter, ns, ez, nc, es, se, lp, af, ac, ne (16 columns). For METAL meta-analysis output with SI/I2/CQ/ED columns, use duckdb_munge_metal.                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `duckdb_munge_metal`        | table_macro | table    | `rduckhts_munge`         | Extended munge macro with METAL meta-analysis output columns. Same as duckdb_munge but additionally emits: si (imputation info, from INFO input), i2 (Cochran’s I² heterogeneity, from HET_I2), cq (Cochran’s Q -log10 p, from HET_LP or -log10(HET_P)), and ed (effect direction string, from DIRE; +/- flipped on allele swap). The R wrapper rduckhts_munge() auto-dispatches to this macro when metal keys (INFO, HET_I2, HET_P, HET_LP, DIRE) are present in the resolved column map.                                                                                                                                                                                                                                                                                                                                                                                   |

### Sequence UDFs

| Function          | Kind   | Returns      | R helper | Description                                                                                                                                                                                                                                                                                                                                                                                                             |
|-------------------|--------|--------------|----------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `seq_revcomp`     | scalar | VARCHAR      |          | Compute the reverse complement of a DNA sequence using A, C, G, T, and N bases. Overloaded: accepts either a VARCHAR text sequence (returns VARCHAR) or a UTINYINT\[\] of htslib nt16 codes as produced by read_bam(sequence_encoding := ‘nt16’) (returns UTINYINT\[\]); the nt16 overload is bit-identical to the text path after decoding, so BAM pipelines can reverse-complement without leaving the nt16 encoding. |
| `seq_canonical`   | scalar | VARCHAR      |          | Return the lexicographically smaller of a sequence and its reverse complement. Overloaded: accepts either a VARCHAR text sequence (returns VARCHAR) or a UTINYINT\[\] of htslib nt16 codes as produced by read_bam(sequence_encoding := ‘nt16’) (returns UTINYINT\[\]); the nt16 overload compares by decoded base order and is bit-identical to the text path after decoding.                                          |
| `seq_hash_2bit`   | scalar | UBIGINT      |          | Encode a short DNA sequence as a 2-bit unsigned integer hash. Overloaded to also accept a UTINYINT\[\] of htslib nt16 codes (from read_bam(sequence_encoding := ‘nt16’)); non-ACGT codes yield NULL, bit-identical to the text path.                                                                                                                                                                                    |
| `seq_encode_4bit` | scalar | UTINYINT\[\] |          | Encode an IUPAC DNA sequence as a list of 4-bit base codes, preserving ambiguity symbols including N.                                                                                                                                                                                                                                                                                                                   |
| `seq_decode_4bit` | scalar | VARCHAR      |          | Decode a list of 4-bit IUPAC DNA base codes back into a sequence string.                                                                                                                                                                                                                                                                                                                                                |
| `seq_gc_content`  | scalar | DOUBLE       |          | Compute GC fraction for a DNA sequence as a value between 0 and 1. Overloaded: accepts either a VARCHAR text sequence or a UTINYINT\[\] of htslib nt16 codes as produced by read_bam(sequence_encoding := ‘nt16’); the nt16 overload classifies codes directly and is bit-identical to the text path, so BAM pipelines can compute GC without decoding sequences back to text.                                          |
| `seq_kmers`       | table  | table        |          | Expand a sequence into positional k-mers with optional canonicalization.                                                                                                                                                                                                                                                                                                                                                |

### SAM Flag UDFs

| Function                               | Kind   | Returns | R helper | Description                                                                                                                                                                        |
|----------------------------------------|--------|---------|----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `sam_flag_bits`                        | scalar | STRUCT  |          | Decode a SAM flag into a struct of boolean bit fields using explicit SAM-oriented names such as `is_paired`, `is_proper_pair`, `is_next_segment_unmapped`, and `is_supplementary`. |
| `sam_flag_has`                         | scalar | BOOLEAN |          | Test whether any bits from the provided SAM flag mask are set in a flag value.                                                                                                     |
| `is_forward_aligned`                   | scalar | BOOLEAN |          | Test whether a mapped segment is aligned to the forward strand. Returns `NULL` for unmapped segments because SAM flag `0x10` does not define genomic strand when `0x4` is set.     |
| `is_paired`                            | scalar | BOOLEAN |          | Test whether the SAM flag indicates that the template has multiple segments in sequencing (`0x1`).                                                                                 |
| `is_proper_pair`                       | scalar | BOOLEAN |          | Test whether the SAM flag indicates that each segment is properly aligned according to the aligner (`0x2`).                                                                        |
| `is_unmapped`                          | scalar | BOOLEAN |          | Test whether the read itself is unmapped according to the SAM flag.                                                                                                                |
| `is_next_segment_unmapped`             | scalar | BOOLEAN |          | Test whether the next segment in the template is flagged as unmapped (`0x8`).                                                                                                      |
| `is_reverse_complemented`              | scalar | BOOLEAN |          | Test whether `SEQ` is stored reverse complemented (`0x10`); for mapped reads this corresponds to reverse-strand alignment.                                                         |
| `is_next_segment_reverse_complemented` | scalar | BOOLEAN |          | Test whether `SEQ` of the next segment in the template is stored reverse complemented (`0x20`).                                                                                    |
| `is_first_segment`                     | scalar | BOOLEAN |          | Test whether the read is marked as the first segment in the template.                                                                                                              |
| `is_last_segment`                      | scalar | BOOLEAN |          | Test whether the read is marked as the last segment in the template.                                                                                                               |
| `is_secondary`                         | scalar | BOOLEAN |          | Test whether the alignment is marked as secondary.                                                                                                                                 |
| `is_qc_fail`                           | scalar | BOOLEAN |          | Test whether the read failed vendor or pipeline quality checks.                                                                                                                    |
| `is_duplicate`                         | scalar | BOOLEAN |          | Test whether the alignment is flagged as a duplicate.                                                                                                                              |
| `is_supplementary`                     | scalar | BOOLEAN |          | Test whether the alignment is marked as supplementary.                                                                                                                             |

### CIGAR Utils

| Function                     | Kind   | Returns | R helper | Description                                                                                                                                                                                                                                                                                   |
|------------------------------|--------|---------|----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `cigar_has_soft_clip`        | scalar | BOOLEAN |          | Test whether a CIGAR string contains any soft-clipped segment (`S`). Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                                                |
| `cigar_has_hard_clip`        | scalar | BOOLEAN |          | Test whether a CIGAR string contains any hard-clipped segment (`H`). Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                                                |
| `cigar_left_soft_clip`       | scalar | BIGINT  |          | Return the left-end soft-clipped length from a CIGAR string, or zero if the alignment does not start with `S`. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.      |
| `cigar_right_soft_clip`      | scalar | BIGINT  |          | Return the right-end soft-clipped length from a CIGAR string, or zero if the alignment does not end with `S`. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.       |
| `cigar_query_length`         | scalar | BIGINT  |          | Return the query-consuming length from a CIGAR string, counting `M`, `I`, `S`, `=`, and `X`. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                        |
| `cigar_aligned_query_length` | scalar | BIGINT  |          | Return the aligned query length from a CIGAR string, counting `M`, `=`, and `X` but excluding clips and insertions. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path. |
| `cigar_reference_length`     | scalar | BIGINT  |          | Return the reference-consuming length from a CIGAR string, counting `M`, `D`, `N`, `=`, and `X`. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                    |
| `cigar_has_op`               | scalar | BOOLEAN |          | Test whether a CIGAR string contains at least one instance of the requested operator. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                               |

## Common Workflows

### Region-aware variant and alignment queries

``` r
bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
bcf_index_path <- system.file("extdata", "vcf_file.bcf.csi", package = "Rduckhts")
bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
bam_index_path <- system.file("extdata", "range.bam.bai", package = "Rduckhts")

rduckhts_bcf(
  con, "variants_idx", bcf_path,
  region = "1:3000150-3000151",
  index_path = bcf_index_path,
  overwrite = TRUE
)
dbGetQuery(con, "SELECT CHROM, POS, REF, ALT FROM variants_idx")
#>   CHROM     POS REF ALT
#> 1     1 3000150   C   T
#> 2     1 3000151   C   T

rduckhts_bam(
  con, "bam_idx_reads", bam_path,
  region = "CHROMOSOME_I:1-1000",
  index_path = bam_index_path,
  overwrite = TRUE
)
dbGetQuery(con, "SELECT QNAME, FLAG, POS, MAPQ FROM bam_idx_reads")
#>                           QNAME FLAG POS MAPQ
#> 1 HS18_09653:4:1315:19857:61712  145 914   23
#> 2 HS18_09653:4:1308:11522:27107  161 934    0
```

### Variant normalization

[`rduckhts_bcftools_norm()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bcftools_norm.md)
wraps the bundled `duckhts_bcftools_norm(...)` table macro and keeps the
original input columns alongside normalized position/reference/ALT
outputs.

``` r
norm_fa <- system.file("extdata", "liftover_repeat_src.fa", package = "Rduckhts")
invisible(DBI::dbExecute(
  con,
  paste(
    "CREATE OR REPLACE TEMP TABLE readme_norm AS SELECT * FROM (VALUES",
    "('chrS', 2, 'T', 'TT,TTT'),",
    "('chrS', 2, 'T', '*,TT')",
    ") AS t(chrom, pos, ref, alt)"
  )
))

norm_out <- rduckhts_bcftools_norm(
  con,
  "readme_norm",
  norm_fa,
  split_multiallelic = TRUE
)
norm_out[order(norm_out$alt, norm_out$alt_index),
         c("chrom", "pos", "ref", "alt", "alt_index", "pos_normed", "ref_normed", "alt_normed", "norm_status")]
#>   chrom pos ref    alt alt_index pos_normed ref_normed alt_normed
#> 4  chrS   2   T   *,TT         1          2          T          *
#> 2  chrS   2   T   *,TT         2          1          G         GT
#> 3  chrS   2   T TT,TTT         1          1          G         GT
#> 1  chrS   2   T TT,TTT         2          1          G        GTT
#>        norm_status
#> 4 SpanningDeletion
#> 2       Normalized
#> 3       Normalized
#> 1       Normalized
```

### VariantKey + RegionKey

DuckHTS now bundles the official VariantKey / RegionKey C API. The SQL
helper `variantkey(...)` accepts 1-based VCF positions to match bcftools
`%VKX`, while `regionkey(...)` keeps 0-based half-open interval
semantics for span work. Large, ambiguous, and symbolic alleles still
encode through the official hashed nonreversible VariantKey mode, but
those keys do not encode `END`, `SVLEN`, mate breakend coordinates, or
other SV metadata; use RegionKey explicitly for interval/SV span
indexing. See Nicola Asuni (2018) <https://doi.org/10.1101/473744>.

``` r
dbGetQuery(
  con,
  paste(
    "SELECT variantkey_hex(variantkey('1', 324684, 'C', 'G')) AS vkx,",
    "reverse_variantkey(parse_variantkey_hex('08027a2588b00000')) AS reversed"
  )
)
#>                vkx reversed.chrom reversed.chrom_code reversed.pos
#> 1 08027a2588b00000              1                   1       324684
#>   reversed.pos0 reversed.ref reversed.alt reversed.refalt_code
#> 1        324683            C            G            145752064
#>   reversed.reversible
#> 1                TRUE

dbGetQuery(
  con,
  paste(
    "SELECT regionkey_hex(regionkey('X', 1007, 1807, 1)) AS rkx,",
    "are_overlapping_regionkeys(regionkey('X', 1007, 1807, 1), parse_regionkey_hex('b80001f78000387a')) AS overlaps"
  )
)
#>                rkx overlaps
#> 1 b80001f78000387a     TRUE
```

### Interval + reference helpers

``` r
bed_path <- system.file("extdata", "targets.bed", package = "Rduckhts")
fai_path <- tempfile("duckhts_readme_", fileext = ".fai")
rduckhts_fasta_index(con, fasta_path, index_path = fai_path)
#>   success                                        index_path
#> 1    TRUE <tempfile>

rduckhts_bed(con, "targets", bed_path, overwrite = TRUE)
dbGetQuery(con, "SELECT chrom, start, \"end\", name, block_count FROM targets")
#>            chrom start end    name block_count
#> 1   CHROMOSOME_I     0  10 target1           2
#> 2   CHROMOSOME_I    10  20 target2           1
#> 3  CHROMOSOME_II     0   8 target3          NA
#> 4 CHROMOSOME_III     0   6 target4           1

rduckhts_fasta_nuc(con, fasta_path, bed_path = bed_path, index_path = fai_path)
#>            chrom start end pct_at pct_gc num_a num_c num_g num_t num_n
#> 1   CHROMOSOME_I     0  10  0.400  0.600     2     4     2     2     0
#> 2   CHROMOSOME_I    10  20  0.500  0.500     4     3     2     1     0
#> 3  CHROMOSOME_II     0   8  0.375  0.625     2     4     1     1     0
#> 4 CHROMOSOME_III     0   6  0.500  0.500     2     2     1     1     0
#>   num_other seq_len
#> 1         0      10
#> 2         0      10
#> 3         0       8
#> 4         0       6
rduckhts_fasta_nuc(con, fasta_path, bin_width = 10, region = "CHROMOSOME_I:1-20", index_path = fai_path)
#>          chrom start end pct_at pct_gc num_a num_c num_g num_t num_n num_other
#> 1 CHROMOSOME_I     0  10    0.4    0.6     2     4     2     2     0         0
#> 2 CHROMOSOME_I    10  20    0.5    0.5     4     3     2     1     0         0
#>   seq_len
#> 1      10
#> 2      10
unlink(fai_path)
```

### cgranges registry entry points

The bundled extension also exposes SQL-first `duckhts_cgranges_*` entry
points. These are session-scoped interval indexes that you can populate
either row-wise or in bulk from a SQL query, then query through
`duckhts_cgranges_overlaps(...)`. For row-preserving filters or count
annotations over provider rows, use the vectorized scalar helpers
`duckhts_cgranges_has_overlap(...)` and
`duckhts_cgranges_count_overlaps(...)` directly in queries over
`read_bed(...)`, `read_bam(...)`, `read_bcf(...)`, or regular tables.
For streaming one-row-per-hit expansion while keeping provider columns,
use `duckhts_cgranges_overlaps_list(...)` with `UNNEST(...)`. The older
`duckhts_cgranges_overlaps_bulk(...)` table function still accepts a
probe query and emits matching indexed intervals in one table-function
call; that bulk query runs on the extension-owned helper connection, so
use a regular table or view rather than a temp table. There is no
dedicated R wrapper yet, so use them through `DBI`.

``` r
DBI::dbGetQuery(con, "SELECT duckhts_cgranges_create('readme_idx') AS ok")
#>     ok
#> 1 TRUE
DBI::dbGetQuery(con, "SELECT duckhts_cgranges_add('readme_idx', 'chr1', 10, 20, 'a') AS ok")
#>     ok
#> 1 TRUE
DBI::dbGetQuery(con, "SELECT duckhts_cgranges_add('readme_idx', 'chr1', 30, 40, 'b') AS ok")
#>     ok
#> 1 TRUE
DBI::dbGetQuery(con, "SELECT duckhts_cgranges_index('readme_idx') AS ok")
#>     ok
#> 1 TRUE
DBI::dbGetQuery(
  con,
  paste(
    "SELECT interval_ordinal, label, interval_chrom, interval_start, interval_end",
    "FROM duckhts_cgranges_overlaps('readme_idx', 'chr1', 35, 36, query_row_id := 7)"
  )
)
#>   interval_ordinal label interval_chrom interval_start interval_end
#> 1                1     b           chr1             30           40

DBI::dbGetQuery(
  con,
  paste(
    "SELECT duckhts_cgranges_from_query(",
    "  'readme_qry_idx',",
    "  'SELECT * FROM (VALUES (''chr2'', 100, 110, ''alpha''), (''chr2'', 150, 170, ''beta'')) AS t(chrom, start, \"end\", label)',",
    "  'chrom', 'start', 'end', 'label'",
    ") AS ok"
  )
)
#>     ok
#> 1 TRUE
DBI::dbGetQuery(con, "SELECT duckhts_cgranges_index('readme_qry_idx') AS ok")
#>     ok
#> 1 TRUE
DBI::dbGetQuery(
  con,
  paste(
    "SELECT interval_ordinal, label, interval_chrom, interval_start, interval_end",
    "FROM duckhts_cgranges_overlaps('readme_qry_idx', 'chr2', 140, 170, mode := 'contain')"
  )
)
#>   interval_ordinal label interval_chrom interval_start interval_end
#> 1                1  beta           chr2            150          170
DBI::dbExecute(
  con,
  paste(
    "CREATE TABLE readme_probes AS SELECT * FROM (VALUES",
    "(10, 'chr2', 100, 105),",
    "(20, 'chr2', 160, 161),",
    "(30, 'chr2', 500, 510)",
    ") AS t(probe_id, chrom, start, \"end\")"
  )
)
#> [1] 3
DBI::dbGetQuery(
  con,
  paste(
    "SELECT p.probe_id, hit.interval_ordinal, hit.label, hit.label_type,",
    "  hit.interval_chrom, hit.interval_start, hit.interval_end",
    "FROM readme_probes AS p",
    "CROSS JOIN UNNEST(",
    "  duckhts_cgranges_overlaps_list('readme_qry_idx', p.chrom, p.start, p.\"end\")",
    ") AS u(hit)",
    "ORDER BY p.probe_id, hit.interval_ordinal"
  )
)
#>   probe_id interval_ordinal label label_type interval_chrom interval_start
#> 1       10                0 alpha    VARCHAR           chr2            100
#> 2       20                1  beta    VARCHAR           chr2            150
#>   interval_end
#> 1          110
#> 2          170
DBI::dbGetQuery(
  con,
  paste(
    "SELECT query_row_id, interval_ordinal, label, interval_chrom, interval_start, interval_end",
    "FROM duckhts_cgranges_overlaps_bulk(",
    "  'readme_qry_idx',",
    "  'SELECT probe_id, chrom, start, \"end\" FROM readme_probes',",
    "  'chrom', 'start', 'end',",
    "  query_row_id_col := 'probe_id'",
    ")",
    "ORDER BY query_row_id, interval_ordinal"
  )
)
#>   query_row_id interval_ordinal label interval_chrom interval_start
#> 1           10                0 alpha           chr2            100
#> 2           20                1  beta           chr2            150
#>   interval_end
#> 1          110
#> 2          170
DBI::dbGetQuery(con, "SELECT duckhts_cgranges_destroy('readme_idx') AS ok")
#>     ok
#> 1 TRUE
DBI::dbGetQuery(con, "SELECT duckhts_cgranges_destroy('readme_qry_idx') AS ok")
#>     ok
#> 1 TRUE
```

### Fixed-bin native counting

[`rduckhts_bam_bin_counts()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_bam_bin_counts.md)
exposes the native fixed-width read-start counting kernel. It returns
one row per bin across the selected contig span, including zero-count
bins, which makes it suitable as a dense CNV binning primitive. It can
also add one-pass GC and MAPQ summaries on the same scan.

``` r
mixed_cram <- system.file("extdata", "fixture_mixed.cram", package = "Rduckhts")
fixture_ref <- system.file("extdata", "fixture_ref.fa", package = "Rduckhts")

bin_counts <- rduckhts_bam_bin_counts(
  con,
  mixed_cram,
  5000,
  reference = fixture_ref,
  rmdup = "streaming",
  stats = "gc,mq"
)

bin_counts[, c(
  "bin_id", "count_total", "count_fwd", "count_rev",
  "count_pre", "gc_perc_pre", "gc_perc_post", "mean_mapq_post"
)]
#>    bin_id count_total count_fwd count_rev count_pre gc_perc_pre gc_perc_post
#> 1       0           2         1         1         4         0.5            0
#> 2       1           2         1         1         2         0.0            0
#> 3       2           1         1         0         2         1.0            1
#> 4       3           0         0         0         0          NA           NA
#> 5       4           0         0         0         0          NA           NA
#> 6       5           0         0         0         0          NA           NA
#> 7       6           0         0         0         0          NA           NA
#> 8       7           0         0         0         0          NA           NA
#> 9       8           0         0         0         0          NA           NA
#> 10      9           0         0         0         0          NA           NA
#>    mean_mapq_post
#> 1              60
#> 2              60
#> 3              60
#> 4              NA
#> 5              NA
#> 6              NA
#> 7              NA
#> 8              NA
#> 9              NA
#> 10             NA
```

### Mosdepth-compatible coverage outputs

[`rduckhts_mosdepth()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_mosdepth.md)
writes mosdepth-style outputs to disk and returns the paths it created.
This example writes windowed fragment coverage from the bundled BAM
fixture and previews the generated regions BED.gz.

``` r
mos_prefix <- tempfile("duckhts_readme_mosdepth_")
mos_out <- rduckhts_mosdepth(
  con,
  prefix = mos_prefix,
  path = bam_path,
  chrom = "CHROMOSOME_II",
  by = "1000",
  no_per_base = TRUE,
  fragment_mode = TRUE,
  use_median = TRUE,
  overwrite = TRUE
)

mos_out[, c("summary_path", "regions_path")]
#>                                                                  summary_path
#> 1 <tempfile>
#>                                                            regions_path
#> 1 <tempfile>

utils::read.delim(
  gzfile(mos_out$regions_path[[1]]),
  header = FALSE,
  sep = "\t",
  nrows = 3,
  col.names = c("chrom", "start", "end", "depth")
)
#>           chrom start  end depth
#> 1 CHROMOSOME_II     0 1000     0
#> 2 CHROMOSOME_II  1000 2000     5
#> 3 CHROMOSOME_II  2000 3000     3

unlink(
  c(
    paste0(mos_prefix, ".mosdepth.summary.txt"),
    paste0(mos_prefix, ".mosdepth.global.dist.txt"),
    paste0(mos_prefix, ".mosdepth.region.dist.txt"),
    paste0(mos_prefix, ".regions.bed.gz"),
    paste0(mos_prefix, ".regions.bed.gz.csi")
  ),
  force = TRUE
)
```

### Liftover score-style rows

``` r
lift_src <- tempfile("duckhts_liftover_src_", fileext = ".fa")
lift_dst <- tempfile("duckhts_liftover_dst_", fileext = ".fa")
lift_chain <- tempfile("duckhts_liftover_", fileext = ".chain")

writeLines(c(
  ">chrF",
  "ACGTACGTAA",
  ">chrR",
  "AACCGGTTAA"
), lift_src)
writeLines(c(
  ">chrLiftF",
  "ACGTACGTAA",
  ">chrLiftR",
  "TTAACCGGTT"
), lift_dst)
writeLines(c(
  "chain 100 chrF 10 + 0 10 chrLiftF 10 + 0 10 1",
  "10",
  "",
  "chain 100 chrR 10 + 0 10 chrLiftR 10 - 0 10 2",
  "10"
), lift_chain)

rduckhts_fasta_index(con, lift_src, index_path = paste0(lift_src, ".fai"))
#>   success                                                 index_path
#> 1    TRUE <tempfile>
rduckhts_fasta_index(con, lift_dst, index_path = paste0(lift_dst, ".fai"))
#>   success                                                 index_path
#> 1    TRUE <tempfile>

lifted <- rduckhts_liftover(
  con,
  query = paste(
    "SELECT * FROM (VALUES",
    "('chrF', 2, 'C', 'T'),",
    "('chrR', 2, 'A', 'G'),",
    "('chrF', 11, 'A', 'T')",
    ") AS t(chrom, pos, ref, alt)"
  ),
  chain_path = lift_chain,
  dst_fasta_ref = lift_dst,
  ref_col = "ref",
  alt_col = "alt",
  src_fasta_ref = lift_src
)

lifted[, c(
  "src_chrom", "src_pos", "dest_chrom", "dest_pos",
  "dest_ref", "dest_alt", "mapped", "reverse_complemented",
  "reject_reason", "note"
)]
#>   src_chrom src_pos dest_chrom dest_pos dest_ref dest_alt mapped
#> 1      chrF       2   chrLiftF        2        C        T   TRUE
#> 2      chrR       2   chrLiftR        9        T        C   TRUE
#> 3      chrF      11       <NA>       NA     <NA>     <NA>  FALSE
#>   reverse_complemented     reject_reason note
#> 1                FALSE              <NA> <NA>
#> 2                 TRUE              <NA> <NA>
#> 3                FALSE SourceRefMismatch <NA>

unlink(c(lift_src, paste0(lift_src, ".fai"), lift_dst, paste0(lift_dst, ".fai"), lift_chain))
```

### Munge score-style rows

``` r
munge_fasta <- tempfile("duckhts_munge_", fileext = ".fa")
writeLines(c(
  ">chrF",
  "ACGTACGTAA"
), munge_fasta)
rduckhts_fasta_index(con, munge_fasta, index_path = paste0(munge_fasta, ".fai"))
#>   success                                          index_path
#> 1    TRUE <tempfile>

munge_out <- rduckhts_munge(
  con,
  query = paste(
    "SELECT * FROM (VALUES",
    "('rs1', 2, 'chrF', 'A', 'C', 0.01, 1.10, 0.20, 0.98, 0.10, 0.01, 1000),",
    "('rs2', 2, 'chrF', 'C', 'A', 0.02, 0.90, -0.20, 0.98, 0.90, 0.01, 1000)",
    ") AS t(SNP, BP, CHR, A1, A2, P, OR_VALUE, BETA, INFO, FRQ, SE, N)"
  ),
  fasta_ref = munge_fasta,
  column_map = c(
    SNP = "SNP", BP = "BP", CHR = "CHR", A1 = "A1", A2 = "A2",
    P = "P", OR = "OR_VALUE", BETA = "BETA", INFO = "INFO", FRQ = "FRQ", SE = "SE", N = "N"
  )
)

munge_out[, c("chrom", "pos", "id", "ref", "alt", "alleles_swapped", "filter", "af", "es", "ns")]
#>   chrom pos  id ref alt alleles_swapped filter  af  es   ns
#> 1  chrF   2 rs2   C   A            TRUE   <NA> 0.1 0.2 1000
#> 2  chrF   2 rs1   C   A           FALSE   <NA> 0.1 0.2 1000

unlink(c(munge_fasta, paste0(munge_fasta, ".fai")))
```

### Polygenic risk scoring

[`rduckhts_score()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_score.md)
computes per-sample polygenic risk scores (PRS) from a genotype VCF/BCF
and one or more GWAS summary statistics files, wrapping
`bcftools_score`.

``` r
vcf_path    <- system.file("extdata", "score_input.vcf",        package = "Rduckhts")
dosage_path <- system.file("extdata", "score_dosage.vcf",       package = "Rduckhts")
sumf_path   <- system.file("extdata", "score_summary.tsv",      package = "Rduckhts")
gwas_path   <- system.file("extdata", "score_gwas_summary.vcf", package = "Rduckhts")

# Hard-call (GT) PRS with PLINK-format summary statistics
# S1: 0×0.5 + 1×(−0.2) + 2×1.0 = 1.8
# S2: 1×0.5 + 2×(−0.2) + 0×1.0 = 0.1
gt_prs <- rduckhts_score(con, vcf_path, sumf_path, use = "GT", columns = "PLINK")
gt_prs[, c("SAMPLE", "score_summary")]
#>   SAMPLE score_summary
#> 1     S1    1.79999995
#> 2     S2    0.09999999

# Multiple TSV/SSF summary files are scored in one genotype scan
sumf_na_path <- system.file("extdata", "score_summary_na.tsv", package = "Rduckhts")
multi_prs <- rduckhts_score(con, vcf_path, c(sumf_path, sumf_na_path),
                            use = "GT", columns = "PLINK")
multi_prs[, c("SAMPLE", "score_summary", "score_summary_na")]
#>   SAMPLE score_summary score_summary_na
#> 1     S1    1.79999995              2.0
#> 2     S2    0.09999999              0.5

# Optional audit log records loaded/matched/allele-mismatch marker counts
sumf_mismatch_path <- system.file("extdata", "score_summary_mismatch.tsv", package = "Rduckhts")
score_log <- tempfile("duckhts_score_", fileext = ".log")
invisible(rduckhts_score(con, vcf_path, c(sumf_path, sumf_mismatch_path),
                         use = "GT", columns = "PLINK", log_path = score_log))
read.delim(score_log, comment.char = "#")[, c("summary_name", "loaded_markers",
                                                "matched_markers", "allele_mismatch_markers")]
#>             summary_name loaded_markers matched_markers allele_mismatch_markers
#> 1          score_summary              3               3                       0
#> 2 score_summary_mismatch              3               0                       3

# Dosage-based PRS (DS field) for imputed genotypes
# S1: 0.1×0.5 + 0.8×(−0.2) + 1.8×1.0 = 1.69
# S2: 1.0×0.5 + 1.9×(−0.2) + 0.2×1.0 = 0.32
ds_prs <- rduckhts_score(con, dosage_path, sumf_path, use = "DS", columns = "PLINK")
ds_prs[, c("SAMPLE", "score_summary")]
#>   SAMPLE score_summary
#> 1     S1          1.69
#> 2     S2          0.32

# GWAS-VCF multi-PRS: each FORMAT/ES sample column becomes a separate PRS track
gwas_prs <- rduckhts_score(con, vcf_path, gwas_path, use = "GT")
gwas_prs[, c("SAMPLE", "PRS_A", "PRS_B")]
#>   SAMPLE      PRS_A PRS_B
#> 1     S1 1.79999995   1.0
#> 2     S2 0.09999999   0.3
```

### Compression + tabix round-trips

``` r
bed_src <- system.file("extdata", "targets.bed", package = "Rduckhts")
bam_src <- system.file("extdata", "range.bam", package = "Rduckhts")
bcf_src <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")

tmp_bed <- tempfile("duckhts_targets_", fileext = ".bed")
tmp_bgz <- paste0(tmp_bed, ".gz")
tmp_tbi <- paste0(tmp_bgz, ".tbi")
tmp_roundtrip <- tempfile("duckhts_targets_roundtrip_", fileext = ".bed")
tmp_bai <- tempfile("duckhts_range_", fileext = ".bam.bai")
tmp_csi <- tempfile("duckhts_variants_", fileext = ".bcf.csi")
file.copy(bed_src, tmp_bed, overwrite = TRUE)
#> [1] TRUE

bgzip_meta <- rduckhts_bgzip(
  con, tmp_bed,
  output_path = tmp_bgz,
  threads = 1,
  keep = TRUE,
  overwrite = TRUE
)
bgzip_meta[, c("success", "output_path", "bytes_out")]
#>   success                                           output_path bytes_out
#> 1    TRUE <tempfile>       169

bgunzip_meta <- rduckhts_bgunzip(
  con, tmp_bgz,
  output_path = tmp_roundtrip,
  threads = 1,
  keep = TRUE,
  overwrite = TRUE
)
bgunzip_meta[, c("success", "output_path", "bytes_out")]
#>   success                                                  output_path
#> 1    TRUE <tempfile>
#>   bytes_out
#> 1       194

bam_index_meta <- rduckhts_bam_index(
  con, bam_src,
  index_path = tmp_bai,
  threads = 1
)
bam_index_meta
#>   success                                           index_path index_format
#> 1    TRUE <tempfile>          BAI

bcf_index_meta <- rduckhts_bcf_index(
  con, bcf_src,
  index_path = tmp_csi,
  threads = 1
)
bcf_index_meta
#>   success                                             index_path index_format
#> 1    TRUE <tempfile>          CSI

tabix_meta <- rduckhts_tabix_index(
  con, tmp_bgz,
  preset = "bed",
  index_path = tmp_tbi,
  threads = 1
)
tabix_meta
#>   success                                                index_path
#> 1    TRUE <tempfile>
#>   index_format
#> 1          TBI

rduckhts_bed(con, "targets_idx", tmp_bgz, region = "CHROMOSOME_I:1-20", index_path = tmp_tbi, overwrite = TRUE)
dbGetQuery(con, "SELECT * FROM targets_idx")
#>          chrom start end    name score strand thick_start thick_end item_rgb
#> 1 CHROMOSOME_I     0  10 target1   100      +           0        10  255,0,0
#> 2 CHROMOSOME_I    10  20 target2   200      -          10        20  0,0,255
#>   block_count block_sizes block_starts extra
#> 1           2         5,5          0,5  <NA>
#> 2           1          10            0  <NA>

unlink(c(tmp_bed, tmp_bgz, tmp_tbi, tmp_roundtrip, tmp_bai, tmp_csi))
```

## Sequence UDFs

The extension also exposes sequence utility UDFs directly in DuckDB SQL,
including 4-bit IUPAC DNA encode/decode helpers. These can be applied to
`SEQUENCE` columns from FASTA and FASTQ scans.

``` r
dbGetQuery(
  con,
  "SELECT
     NAME,
     seq_hash_2bit(substr(SEQUENCE, 1, 12)) AS hash_2bit_prefix,
     seq_encode_4bit(substr(SEQUENCE, 1, 16)) AS codes,
     seq_decode_4bit(seq_encode_4bit(substr(SEQUENCE, 1, 16))) AS roundtrip
   FROM sequences
   LIMIT 2"
)
#>            NAME hash_2bit_prefix                                          codes
#> 1  CHROMOSOME_I          9898352 4, 2, 2, 8, 1, 1, 4, 2, 2, 8, 1, 1, 4, 2, 2, 8
#> 2 CHROMOSOME_II          6038978 2, 2, 8, 1, 1, 4, 2, 2, 8, 1, 1, 4, 2, 2, 8, 1
#>          roundtrip
#> 1 GCCTAAGCCTAAGCCT
#> 2 CCTAAGCCTAAGCCTA

dbGetQuery(
  con,
  "SELECT
     NAME,
     MATE,
     seq_encode_4bit(substr(SEQUENCE, 1, 12)) AS codes,
     seq_decode_4bit(seq_encode_4bit(substr(SEQUENCE, 1, 12))) AS roundtrip
   FROM reads
   LIMIT 2"
)
#>                              NAME MATE                              codes
#> 1 HS25_09827:2:1201:1505:59795#49    1 2, 2, 4, 8, 8, 1, 4, 1, 4, 2, 1, 8
#> 2 HS25_09827:2:1201:1505:59795#49    2 1, 1, 4, 4, 1, 1, 1, 4, 1, 1, 4, 4
#>      roundtrip
#> 1 CCGTTAGAGCAT
#> 2 AAGGAAAGAAGG
```

### FASTA region queries

`read_fasta` supports indexed region queries via
`rduckhts_fasta(..., region = ...)`.

``` r
fai_path <- tempfile("duckhts_readme_", fileext = ".fai")
fai_info <- rduckhts_fasta_index(con, fasta_path, index_path = fai_path)
fai_info
#>   success                                        index_path
#> 1    TRUE <tempfile>

rduckhts_fasta(
  con, "fasta_region", fasta_path,
  region = "CHROMOSOME_I:1-25",
  overwrite = TRUE
)
dbGetQuery(con, "SELECT NAME, length(SEQUENCE) AS n FROM fasta_region")
#>           NAME  n
#> 1 CHROMOSOME_I 25
unlink(fai_path)
```

## Examples

### Region Queries

Region queries can use implicit sidecar indexes or an explicit
`index_path` for custom index names/locations.

``` r
bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
bcf_index_path <- system.file("extdata", "vcf_file.bcf.csi", package = "Rduckhts")
rduckhts_bcf(con, "variants", bcf_path, overwrite = TRUE)
variants <- dbGetQuery(con, "SELECT * FROM variants LIMIT 5")
variants
#>   CHROM     POS    ID  REF  ALT  QUAL FILTER INFO_TEST   INFO_DP4 INFO_AC
#> 1     1 3000150  <NA>    C    T  59.2   PASS        NA       NULL       2
#> 2     1 3000151  <NA>    C    T  59.2   PASS        NA       NULL       2
#> 3     1 3062915  id3D GTTT    G  12.9    q10        NA 1, 2, 3, 4       2
#> 4     1 3062915 idSNP    G T, C  12.6   test         5 1, 2, 3, 4    1, 1
#> 5     1 3106154  <NA> CAAA    C 342.0   PASS        NA       NULL       2
#>   INFO_AN INFO_INDEL INFO_STR FORMAT_TT_A FORMAT_GT_A FORMAT_GQ_A FORMAT_DP_A
#> 1       4      FALSE     <NA>        NULL         0/1         245          NA
#> 2       4      FALSE     <NA>        NULL         0/1         245          32
#> 3       4       TRUE     test        NULL         0/1         409          35
#> 4       3      FALSE     <NA>        0, 1         0/1         409          35
#> 5       4      FALSE     <NA>        NULL         0/1         245          32
#>                  FORMAT_GL_A FORMAT_TT_B FORMAT_GT_B FORMAT_GQ_B FORMAT_DP_B
#> 1                       NULL        NULL         0/1         245          NA
#> 2                       NULL        NULL         0/1         245          32
#> 3               -20, -5, -20        NULL         0/1         409          35
#> 4 -20, -5, -20, -20, -5, -20        0, 1           2         409          35
#> 5                       NULL        NULL         0/1         245          32
#>    FORMAT_GL_B
#> 1         NULL
#> 2         NULL
#> 3 -20, -5, -20
#> 4 -20, -5, -20
#> 5         NULL

rduckhts_bcf(
  con, "variants_idx", bcf_path,
  region = "1:3000150-3000151",
  index_path = bcf_index_path,
  overwrite = TRUE
)
dbGetQuery(con, "SELECT count(*) AS n FROM variants_idx")
#>   n
#> 1 2

# Convert to a round-trip-aware Parquet copy with VCF header metadata.
# Extra named metadata is merged into the Parquet key-value metadata.
parquet_path <- tempfile(fileext = ".parquet")
rduckhts_bcf_convert_parquet(
  con,
  bcf_path,
  parquet_path,
  columns = c("CHROM", "POS", "REF", "ALT"),
  metadata = list(project = "demo-cohort", batch = "1"),
  overwrite = TRUE
)

metadata_preview <- dbGetQuery(con, sprintf(
  "SELECT key::VARCHAR AS key, left(value::VARCHAR, 40) AS value_prefix
   FROM parquet_kv_metadata('%s')
   WHERE key::VARCHAR IN ('duckhts_write_format_version', 'duckhts_reader', 'project', 'batch', 'vcf_header')
   ORDER BY key::VARCHAR",
  parquet_path
))
metadata_preview
#>                            key                              value_prefix
#> 1                        batch                                         1
#> 2               duckhts_reader                                  read_bcf
#> 3 duckhts_write_format_version                                         1
#> 4                      project                               demo-cohort
#> 5                   vcf_header ##fileformat=VCFv4.1\\x5Cn##FILTER=<ID=PA

# The same converter helpers support partitioned output for DuckLake-style
# registration of premade Parquet files.
gff_path <- system.file("extdata", "gff_file.gff.gz", package = "Rduckhts")
gff_parquet_dir <- tempfile("duckhts_gff_parquet_")
rduckhts_gff_convert_parquet(
  con,
  gff_path,
  gff_parquet_dir,
  columns = c("seqname", "source", "feature", "start", "end"),
  partition_by = "feature",
  overwrite = TRUE
)
length(list.files(gff_parquet_dir, pattern = "\\.parquet$", recursive = TRUE))
#> [1] 5
unlink(c(parquet_path, gff_parquet_dir), recursive = TRUE)

# Span-oriented index view from the same file
index_spans_preview <- rduckhts_hts_index_spans(con, bcf_path, index_path = bcf_index_path)
head(index_spans_preview[, c("seqname", "tid", "index_type", "chunk_beg_vo", "chunk_end_vo")], 5)
#>   seqname tid index_type chunk_beg_vo chunk_end_vo
#> 1       1   0        CSI         1586         1713
#> 2       1   0        CSI         1713         1973
#> 3       1   0        CSI         1973         2109
#> 4       1   0        CSI         2109         2242
#> 5       1   0        CSI         2242         2372
```

### Remote VCF on S3

S3 files can be query when `htslib` is built with plugins enable. This
is not the case on RTools

``` r
# Example S3 URL (1000 Genomes cohort VCF)
s3_base <- "s3://1000genomes-dragen-v3.7.6/data/cohorts/"
s3_path <- "gvcf-genotyper-dragen-3.7.6/hg19/3202-samples-cohort/"
s3_vcf_file <- "3202_samples_cohort_gg_chr22.vcf.gz"
s3_vcf_uri <- paste0(s3_base, s3_path, s3_vcf_file)

rduckhts_bcf(con, "s3_variants", s3_vcf_uri, region = "chr22:16050000-16050500", overwrite = TRUE)
dbGetQuery(con, "SELECT CHROM, COUNT(*) AS n FROM s3_variants GROUP BY CHROM")
#>   CHROM  n
#> 1 chr22 11
```

### FASTQ files

Three modes for fastq files, single, paired and interleaved

``` r
r1 <- system.file("extdata", "r1.fq", package = "Rduckhts")
r2 <- system.file("extdata", "r2.fq", package = "Rduckhts")
interleaved <- system.file("extdata", "interleaved.fq", package = "Rduckhts")
rduckhts_fastq(con, "paired_reads", r1, mate_path = r2, overwrite = TRUE)
rduckhts_fastq(con, "interleaved_reads", interleaved, interleaved = TRUE, overwrite = TRUE)
pairs <- dbGetQuery(con, "SELECT * FROM paired_reads WHERE MATE = 1 LIMIT 5")
pairs
#>                              NAME DESCRIPTION
#> 1 HS25_09827:2:1201:1505:59795#49        <NA>
#> 2 HS25_09827:2:1201:1559:70726#49        <NA>
#> 3 HS25_09827:2:1201:1564:39627#49        <NA>
#> 4 HS25_09827:2:1201:1565:91731#49        <NA>
#> 5 HS25_09827:2:1201:1624:69925#49        <NA>
#>                                                                                               SEQUENCE
#> 1 CCGTTAGAGCATTTGTTGAAAATGCTTTCCTTGCTCCATGTGATGACTCTGGTGCCCTTGTCAAAAGCCAGCTGGGCCTATTCGTGTGGGTCTGTTTCTG
#> 2 TTGTTAAAATGACCATACCCAAAGTGATCTACAGACTCAATACAATTTCTATTGAAATACCAATCACACTCTTCACAGAACTAGAAAAACAGTTCTAAAA
#> 3 ACGCGGCAATCCAATGTGTGAGTTGAGAAGCGGTGAGGAGGGAATCCTAATTTTATGAGCAGGTCAGGACCGTGGGAGATACCTGACACCTGAGATGGTA
#> 4 GACATGCCATAACATTCATGTTTTATGTGTACAAGTCAATGAATTTTAGTATATTTACAGAGTTGTATGACTGTCTCCACAATCTAATTTTAGGTTTCCA
#> 5 GCCAGCCTCCTTCTCAATGGTCTTTTTAAACATTATATGAAAACCAGACATTTACATTTGATTTCTTTTTCAATACTATACAGTTCTAAGAGAAAAAACA
#>                                                                                                QUALITY
#> 1 CABCFGDEEFFEFHGHGGFFGDIGIJFIFHHGHEIFGHBCGHDIFBE9GIAICGGICFIBFGGHGDGGGHE?GIGDFGGHEGIEJG>;FG<GGHACEFGH
#> 2 CABEFGFFGFHGGGGJGGFFGKIHHJFIEHHHGIEGGEHJGHDHFGHIGICIJEFIFGIF8GGHKFHGGFEI6GGGFIGHGGIE>EFCFHGGGHEJEAJE
#> 3 BACCFGBFGFHGGJGHGGFEGHIGIJHFEH:HHEHGHHBGGH9IAGHGFHIFJFFAFGIFDIGHKEIG<C>F,CGD66?7EFI5EEG>EGGGGD5=HH6E
#> 4 CABFFGFFJFHEGEGJGGDG?FIGHHHBGHHHGIIGHGHGGHDGHFHIDFCIKEGIFHGGII9HFFGGGEEIGGEEHGGEEGDEHFH>FGGGGHAFAHGE
#> 5 CABEFGFGIFGGGJGHGGFH?FDHGHDHGHEHHJCGHHFHDHDHFGHIGHIFFHGHFGGGI9GHF@IGGH;FICGEFEIHGGIEEFC:DEGGGBDJHHFF
#>   MATE                         PAIR_ID
#> 1    1 HS25_09827:2:1201:1505:59795#49
#> 2    1 HS25_09827:2:1201:1559:70726#49
#> 3    1 HS25_09827:2:1201:1564:39627#49
#> 4    1 HS25_09827:2:1201:1565:91731#49
#> 5    1 HS25_09827:2:1201:1624:69925#49
```

### FASTQ quality decoding and per-position histograms

FASTQ quality handling has two separate knobs:

- `input_quality_encoding` controls how incoming FASTQ ASCII is decoded.
  The default is modern `phred33`; use `phred64`, `solexa64`, or `auto`
  only for legacy data.
- `quality_representation` controls how qualities are returned to
  DuckDB: canonical Phred+33 text (`"string"`) or numeric Phred arrays
  (`"phred"`).

The flow is:

1.  Decode FASTQ ASCII using `input_quality_encoding`.
2.  Normalize to numeric Phred qualities.
3.  Return either numeric arrays or canonical Phred+33 text.

`BAM`/`CRAM` reads skip the text-decoding step because qualities are
already stored numerically.

``` r
legacy_fastq <- system.file("extdata", "legacy_phred64.fq", package = "Rduckhts")

rduckhts_detect_quality_encoding(con, legacy_fastq)
#>   format observed_ascii_min observed_ascii_max records_sampled
#> 1  fastq                104                104               1
#>       compatible_encodings guessed_encoding is_ambiguous
#> 1 phred33,phred64,solexa64          phred64         TRUE

quality_hist <- dbGetQuery(
  con,
  sprintf(
    "WITH q AS (
       SELECT NAME, QUALITY
       FROM read_fastq('%s', quality_representation := 'phred')
     ),
     expanded AS (
       SELECT
         NAME,
         generate_subscripts(QUALITY, 1) AS pos,
         unnest(QUALITY) AS q
       FROM q
     )
     SELECT pos, q AS phred, count(*) AS n_reads
     FROM expanded
     GROUP BY pos, phred
     ORDER BY pos, phred
     LIMIT 12",
    fastq_r1
  )
)
quality_hist
#>    pos phred n_reads
#> 1    1    33       1
#> 2    1    34       4
#> 3    2    32       5
#> 4    3    33       4
#> 5    3    34       1
#> 6    4    34       2
#> 7    4    36       2
#> 8    4    37       1
#> 9    5    37       5
#> 10   6    38       5
#> 11   7    33       1
#> 12   7    35       1
```

### GFF/GTF annotation attributes

GFF3 files are read with
[`rduckhts_gff()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_gff.md)
/ SQL `read_gff(...)`; GTF files are read with
[`rduckhts_gtf()`](https://rgenomicsetl.github.io/duckhts/reference/rduckhts_gtf.md)
/ SQL `read_gtf(...)`. `strict = TRUE` enables GFF3 structural
validation. Attribute decoding can be scalar and raw for legacy
convenience (`attributes_map`), grouped and lossless for multi-values
(`attributes_list`, a DuckDB `MAP(VARCHAR, VARCHAR[])`), or exact
parser-style pairs (`attributes_pairs`, a DuckDB
`LIST<STRUCT(key, value, idx)>`).

The extension-level GFF3 implementation is benchmarked and audited
against [GFFBase](https://github.com/Kuanhao-Chao/gffbase) in the
DuckHTS repo:
<https://github.com/RGenomicsETL/duckhts/blob/develop/benchmarks/benchmark_gffbase_conformance.md>.

``` r
gff_path <- system.file("extdata", "gff_file.gff.gz", package = "Rduckhts")
rduckhts_gff(con, "genes", gff_path, attributes_map = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT seqname, start, \"end\" FROM genes WHERE feature = 'gene' LIMIT 5")
#>   seqname   start     end
#> 1       X 2934816 2964270

gff_attrs_path <- system.file("extdata", "gff_attrs.gff3", package = "Rduckhts")
rduckhts_gff(
  con,
  "gff_attrs",
  gff_attrs_path,
  strict = TRUE,
  attributes_list = TRUE,
  attributes_pairs = TRUE,
  overwrite = TRUE
)
dbGetQuery(con, paste(
  "SELECT seqname, feature,",
  "list_extract(map_extract_value(attributes_list, 'Dbxref'), 1) AS first_dbxref,",
  "list_count(attributes_pairs) AS n_attr_pairs FROM gff_attrs"
))
#>   seqname feature first_dbxref n_attr_pairs
#> 1    chr1    gene     GeneID:1            6

gtf_attrs_path <- system.file("extdata", "gtf_attrs.gtf", package = "Rduckhts")
rduckhts_gtf(con, "gtf_attrs", gtf_attrs_path, attributes_list = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT list_extract(map_extract_value(attributes_list, 'note'), 1) AS note FROM gtf_attrs")
#>          note
#> 1 weird; semi
```

### BAM/CRAM

When built with htslib codec, `CRAM` can be opened in addition to `BAM`
files. `index_path` can also be passed for region scans with
non-standard index names.

``` r
cram_path <- system.file("extdata", "range.cram", package = "Rduckhts")
ref_path <- system.file("extdata", "ce.fa", package = "Rduckhts")
bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
bam_index_path <- system.file("extdata", "range.bam.bai", package = "Rduckhts")

rduckhts_bam(con, "cram_reads", cram_path, reference = ref_path, overwrite = TRUE)
cram_reads <- dbGetQuery(con, "SELECT QNAME, FLAG, POS, MAPQ FROM cram_reads LIMIT 5")
cram_reads
#>                           QNAME FLAG  POS MAPQ
#> 1 HS18_09653:4:1315:19857:61712  145  914   23
#> 2 HS18_09653:4:1308:11522:27107  161  934    0
#> 3 HS18_09653:4:2314:14991:85680   83 1020   10
#> 4 HS18_09653:4:2108:14085:93656  147 1122   60
#> 5  HS18_09653:4:1303:4347:38100   83 1137   37

rduckhts_bam(
  con, "bam_idx_reads", bam_path,
  region = "CHROMOSOME_I:1-1000",
  index_path = bam_index_path,
  overwrite = TRUE
)
dbGetQuery(con, "SELECT count(*) AS n FROM bam_idx_reads")
#>   n
#> 1 2
```

### SAMtags + auxiliary tags

Standard SAMtags can be exposed as typed columns, and any remaining tags
are available via `AUXILIARY_TAGS`:

``` r
aux_path <- system.file("extdata", "aux_tags.sam.gz", package = "Rduckhts")
rduckhts_bam(con, "aux_reads", aux_path, standard_tags = TRUE, auxiliary_tags = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT RG, NM, map_extract(AUXILIARY_TAGS, 'XZ') AS XZ FROM aux_reads LIMIT 1")
#>   RG NM  XZ
#> 1 x1  2 foo
```

### Tabix headers + types

Use `header = TRUE` to use the first non-meta row as column names, and
`auto_detect = TRUE` / `column_types` to control column typing:

``` r
tabix_header <- system.file("extdata", "header_tabix.tsv.gz", package = "Rduckhts")
tabix_meta <- system.file("extdata", "meta_tabix.tsv.gz", package = "Rduckhts")

rduckhts_tabix(con, "header_tabix", tabix_header, header = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT chrom, pos FROM header_tabix LIMIT 2")
#>   chrom pos
#> 1  chr1   1
#> 2  chr1   2

rduckhts_tabix(con, "typed_tabix", tabix_meta, auto_detect = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT typeof(column1) AS column1_type FROM typed_tabix LIMIT 1")
#>   column1_type
#> 1       BIGINT

rduckhts_tabix(con, "typed_tabix_explicit", tabix_header,
               header = TRUE,
               column_types = c("VARCHAR", "BIGINT", "VARCHAR"),
               overwrite = TRUE)
dbGetQuery(con, "SELECT pos + 1 AS pos_plus_one FROM typed_tabix_explicit LIMIT 1")
#>   pos_plus_one
#> 1            2
```

### HTS header and index metadata

Use metadata helpers to inspect parsed headers, raw header lines, index
summaries, span-oriented index views, and raw index blobs.

``` r
header_meta <- rduckhts_hts_header(con, bcf_path)
head(header_meta[, c("record_type", "id", "number", "value_type")], 5)
#>   record_type   id number value_type
#> 1  fileformat <NA>   <NA>       <NA>
#> 2      FILTER PASS   <NA>       <NA>
#> 3        INFO TEST      1    Integer
#> 4      FORMAT   TT      A    Integer
#> 5        INFO  DP4      4    Integer

header_raw <- rduckhts_hts_header(con, bcf_path, mode = "raw")
head(header_raw[, c("idx", "raw")], 5)
#>   idx
#> 1   0
#> 2   1
#> 3   2
#> 4   3
#> 5   4
#>                                                                                                                                                          raw
#> 1                                                                                                                                       ##fileformat=VCFv4.1
#> 2                                                                                                        ##FILTER=<ID=PASS,Description="All filters passed">
#> 3                                                                                           ##INFO=<ID=TEST,Number=1,Type=Integer,Description="Testing Tag">
#> 4 ##FORMAT=<ID=TT,Number=A,Type=Integer,Description="Testing Tag, with commas and \\"escapes\\" and escaped escapes combined with \\\\\\"quotes\\\\\\\\\\"">
#> 5                       ##INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">

index_meta <- rduckhts_hts_index(con, bcf_path, index_path = bcf_index_path)
head(index_meta[, c("seqname", "mapped", "unmapped", "index_type")], 5)
#>   seqname mapped unmapped index_type
#> 1       1     11        0        CSI
#> 2       2      1        0        CSI
#> 3       3      1        0        CSI
#> 4       4      2        0        CSI

index_spans <- rduckhts_hts_index_spans(con, bcf_path, index_path = bcf_index_path)
head(index_spans[, c("seqname", "tid", "index_type", "chunk_beg_vo", "chunk_end_vo")], 5)
#>   seqname tid index_type chunk_beg_vo chunk_end_vo
#> 1       1   0        CSI         1586         1713
#> 2       1   0        CSI         1713         1973
#> 3       1   0        CSI         1973         2109
#> 4       1   0        CSI         2109         2242
#> 5       1   0        CSI         2242         2372

index_raw <- rduckhts_hts_index_raw(con, bcf_path, index_path = bcf_index_path)
head(index_raw, 1)
#> [1] index_type
#> [2] '/usr/local/lib/R/site-library/Rduckhts/extdata/vcf_file.bcf.csi'
#> [3] raw
#> <0 rows> (or 0-length row.names)
```

### Remote GTEx tabix example

GTEx eQTL matrices on EBI are tabix-indexed. In browser wasm/webR, this
still depends on CORS policy on both the data object and index object.

``` r
gtex_url <- "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Brain_Cerebellar_Hemisphere.tsv.gz" 
rduckhts_tabix(con, "gtex_eqtl", gtex_url, region = "1:11868-14409",
                  header = TRUE, auto_detect = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT * FROM gtex_eqtl LIMIT 5")
#>          variant r2    pvalue molecular_trait_object_id molecular_trait_id
#> 1 chr1_13550_G_A NA 0.0204520           ENSG00000188290    ENSG00000188290
#> 2 chr1_13550_G_A NA 0.0303633           ENSG00000230699    ENSG00000230699
#> 3 chr1_13550_G_A NA 0.1057900           ENSG00000177757    ENSG00000177757
#> 4 chr1_13550_G_A NA 0.1617190           ENSG00000241860    ENSG00000241860
#> 5 chr1_13550_G_A NA 0.1919580           ENSG00000198744    ENSG00000198744
#>         maf         gene_id median_tpm      beta       se  an ac chromosome
#> 1 0.0114286 ENSG00000188290  6.3960000  0.633986 0.270285 350  4          1
#> 2 0.0114286 ENSG00000230699  0.0674459 -0.980082 0.447861 350  4          1
#> 3 0.0114286 ENSG00000177757  1.2659000  0.631359 0.387738 350  4          1
#> 4 0.0114286 ENSG00000241860  0.1081970 -0.791695 0.562674 350  4          1
#> 5 0.0114286 ENSG00000198744 21.6284000 -0.592354 0.451705 350  4          1
#>   position ref alt type        rsid
#> 1    13550   G   A  SNP rs554008981
#> 2    13550   G   A  SNP rs554008981
#> 3    13550   G   A  SNP rs554008981
#> 4    13550   G   A  SNP rs554008981
#> 5    13550   G   A  SNP rs554008981
```

``` r
dbDisconnect(con, shutdown = TRUE)
```

## References

- DuckDB: <https://duckdb.org/>
- DuckDB Extension API: <https://duckdb.org/docs/extensions/overview>
- DuckDB extension template (C):
  <https://github.com/duckdb/extension-template-c>
- htslib: <https://github.com/samtools/htslib>
- RBCFTools: <https://github.com/RGenomicsETL/RBCFTools>

## License

GPL-3.
