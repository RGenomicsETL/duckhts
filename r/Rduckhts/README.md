
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Rduckhts: DuckDB HTS File Reader Extension for R

[![CRAN
Status](https://www.r-pkg.org/badges/version/Rduckhts)](https://cran.r-project.org/package=Rduckhts)[![R-universe
version](https://RGenomicsETL.r-universe.dev/Rduckhts/badges/version)](https://RGenomicsETL.r-universe.dev/Rduckhts)

`Rduckhts` provides an R interface to a [DuckDB](https://duckdb.org/)
`HTS` (High Throughput Sequencing) file reader extension. This enables
reading common bioinformatics file formats such as `VCF`/`BCF`,
`SAM`/`BAM`/`CRAM`, `FASTA`, `FASTQ`, BigWig, `GFF`, `GTF`, and
tabix-indexed files directly from `R` using `SQL` queries via
[`duckhts`](https://github.com/RGenomicsETL/duckhts).

## How it works

Following [RBCFTools](https://github.com/RGenomicsETL/RBCFTools), tables
are created and returned instead of data frames. `VCF`/`BCF`,
`SAM`/`BAM`/`CRAM`, `FASTA`, `FASTQ`, BigWig, `GFF`, `GTF`, and `tabix`
formats can be queried. We support region queries for indexed files, and
we target Linux, macOS, and RTools.
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
- Local test servers should return `Accept-Ranges: bytes` and
  `206 Partial   Content` for both the data file and its index sidecars.
  An R-native option is
  [`goserveR`](https://github.com/sounkou-bioinfo/goserveR). Without
  byte-range support the backend warns and falls back to fetching the
  complete object before slicing it locally.

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

Create a package-owned connection with `rduckhts_connect()`, which
explicitly loads the bundled extension without enabling automatic
installation or loading of unrelated DuckDB extensions. Use
`rduckhts_load(con, extension_path = NULL)` only for an existing
connection whose driver permits extension loading. The wrappers break
down into:

- readers: `rduckhts_bcf()`, `rduckhts_bam()`, `rduckhts_fasta()`,
  `rduckhts_fastq()`, `rduckhts_bigwig()`, `rduckhts_gff()`,
  `rduckhts_gtf()`, `rduckhts_tabix()`, `rduckhts_bed()`
- multi-file readers: `rduckhts_bam_multi()`, `rduckhts_bcf_multi()`,
  `rduckhts_fastq_multi()`, `rduckhts_fasta_multi()`,
  `rduckhts_gff_multi()`, `rduckhts_gtf_multi()`,
  `rduckhts_tabix_multi()`, `rduckhts_bed_multi()`
- reference helpers: `rduckhts_fasta_index()`, `rduckhts_fasta_nuc()`
- compression/indexing: `rduckhts_bgzip()`, `rduckhts_bgunzip()`,
  `rduckhts_bam_index()`, `rduckhts_bcf_index()`,
  `rduckhts_tabix_index()`
- metadata helpers: `rduckhts_hts_header()`, `rduckhts_hts_index()`,
  `rduckhts_hts_index_spans()`, `rduckhts_hts_index_raw()`
- htslib diagnostics/linking: `rduckhts_htslib_info()`,
  `rduckhts_htslib_version()`, `rduckhts_htslib_config()`
- Parquet converters: `rduckhts_bcf_convert_parquet()`,
  `rduckhts_bam_convert_parquet()`, `rduckhts_gff_convert_parquet()`,
  `rduckhts_tabix_convert_parquet()`
- SIMD diagnostics: `rduckhts_simd_backend()`,
  `rduckhts_simd_requested_backend()`,
  `rduckhts_simd_backend_available()`, `rduckhts_simd_set_backend()`

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
con <- rduckhts_connect()

rduckhts_fasta(con, "sequences", fasta_path, overwrite = TRUE)
rduckhts_fastq(con, "reads", fastq_r1, mate_path = fastq_r2, overwrite = TRUE)

dbGetQuery(con, "SELECT COUNT(*) AS n FROM sequences")
#>   n
#> 1 7
dbGetQuery(con, "SELECT COUNT(*) AS n FROM reads")
#>    n
#> 1 10
```

### Link a downstream package to the bundled htslib

`rduckhts_htslib_config()` returns the exact installed headers, shared
or static library, required flags, enabled features, and build identity.
With no `link` argument it selects the mode chosen when Rduckhts was
configured. Validation loads DuckHTS and rejects a source/header/runtime
version mismatch. A downstream `configure` script can emit its
`PKG_CPPFLAGS` and `PKG_LIBS` from this one receipt.

``` r
hts_config <- rduckhts_htslib_config()
hts_config[c("contract_version", "htslib_version", "runtime_version", "link")]
#> $contract_version
#> [1] 1
#>
#> $htslib_version
#> [1] "1.24"
#>
#> $runtime_version
#> [1] "1.24"
#>
#> $link
#> [1] "shared"
hts_config$features
#> $cram
#> [1] TRUE
#>
#> $zlib
#> [1] TRUE
#>
#> $bzip2
#> [1] TRUE
#>
#> $lzma
#> [1] TRUE
#>
#> $libdeflate
#> [1] TRUE
#>
#> $curl
#> [1] TRUE
#>
#> $openssl
#> [1] TRUE
#>
#> $plugins
#> [1] TRUE
#>
#> $s3
#> [1] TRUE
#>
#> $gcs
#> [1] TRUE
```

## SIMD dispatch flow

The bundled extension uses explicit runtime SIMD dispatch for
byte-oriented helper kernels, starting with `seq_gc_content(...)`.
`scalar` is always available and is the portable baseline. Optional
platform backends such as `avx2` or `avx512` should be checked with
`rduckhts_simd_backend_available()` before being requested. The `auto`
policy resolves each logical kernel independently from the current
compiled-and-CPU-supported capability mask; use
`rduckhts_simd_kernel_info()` for the per-kernel result and
`rduckhts_simd_set_backend(con, "auto")` to return to runtime
auto-detection.

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
#> 4        fastq_qc             avx2           FALSE

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
dbGetQuery(con, "SELECT filename, count(*) AS n FROM fq_multi GROUP BY ALL ORDER BY filename")
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

Use `rduckhts_functions()` inside R to inspect the generated extension
catalog.

<details>
<summary>
Show generated function catalog
</summary>

## Extension Function Catalog

This section is generated from `functions.yaml`.

### Diagnostics

| Function                                                                                                      | Kind         | R helper                              | Description                                                                                                                                                                                                                                                                                                                                                                      |
|---------------------------------------------------------------------------------------------------------------|--------------|---------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`duckhts_htslib_version`](inst/function_catalog/reference.md#duckhts_htslib_version)                         | scalar       | `rduckhts_htslib_version`             | Return the runtime version reported by the htslib library loaded with DuckHTS. Rduckhts uses this value to reject a downstream linking receipt whose source/header version does not match the loaded library.                                                                                                                                                                    |
| [`duckhts_htslib_features`](inst/function_catalog/reference.md#duckhts_htslib_features)                       | scalar       |                                       | Return the htslib runtime feature bitfield reported by hts_features(). Use duckhts_htslib_feature_string() for the corresponding build description.                                                                                                                                                                                                                              |
| [`duckhts_htslib_feature_string`](inst/function_catalog/reference.md#duckhts_htslib_feature_string)           | scalar       |                                       | Return htslib’s runtime build-feature description, including configured transports, compression libraries, compiler, and build flags. DuckHTS snapshots it once while loading the extension so parallel SQL calls read immutable text.                                                                                                                                           |
| [`duckhts_simd_backend`](inst/function_catalog/reference.md#duckhts_simd_backend)                             | scalar       | `rduckhts_simd_backend`               | Return the current DuckHTS SIMD dispatch label. For explicit scalar or concrete backend requests this is the requested policy; for auto it is the single selected backend when all logical kernels resolve to the same backend, or mixed when per-kernel auto-dispatch resolves to multiple backends. Use duckhts_simd_kernel_info() for per-kernel details.                     |
| [`duckhts_simd_requested_backend`](inst/function_catalog/reference.md#duckhts_simd_requested_backend)         | scalar       | `rduckhts_simd_requested_backend`     | Return the current explicit SIMD backend request, usually auto unless `SELECT backend FROM duckhts_simd_set_backend('auto'\|'scalar'\|backend)` was called. The selected per-kernel backend may differ under auto-dispatch across x86, ARM, wasm, and scalar-only builds.                                                                                                        |
| [`duckhts_simd_backend_compiled`](inst/function_catalog/reference.md#duckhts_simd_backend_compiled)           | scalar       | `rduckhts_simd_backend_compiled`      | Return whether a concrete DuckHTS SIMD backend was compiled into this build. This is independent of whether the current CPU/runtime supports executing that backend; for example avx512 can be compiled but not CPU-supported on the running host.                                                                                                                               |
| [`duckhts_simd_backend_cpu_supported`](inst/function_catalog/reference.md#duckhts_simd_backend_cpu_supported) | scalar       | `rduckhts_simd_backend_cpu_supported` | Return whether the current CPU/runtime supports a concrete DuckHTS SIMD backend, independent of whether DuckHTS compiled an implementation for it. Availability is the intersection of compiled and CPU-supported.                                                                                                                                                               |
| [`duckhts_simd_backend_available`](inst/function_catalog/reference.md#duckhts_simd_backend_available)         | scalar       | `rduckhts_simd_backend_available`     | Return whether a concrete SIMD backend is usable in the current process. Availability means the backend is compiled into DuckHTS and supported by the current CPU/runtime. auto is a selection request rather than a concrete backend and is not reported as available here.                                                                                                     |
| [`duckhts_simd_info`](inst/function_catalog/reference.md#duckhts_simd_info)                                   | table        | `rduckhts_simd_info`                  | Report compiled, runtime-supported and selected status for each concrete DuckHTS SIMD backend.                                                                                                                                                                                                                                                                                   |
| [`duckhts_simd_kernel_info`](inst/function_catalog/reference.md#duckhts_simd_kernel_info)                     | table        | `rduckhts_simd_kernel_info`           | Return one row per logical DuckHTS SIMD kernel showing the concrete backend selected by the current immutable dispatch table, the selected capability, the requested backend policy, whether scalar was used as a per-kernel fallback, and the dispatch mode. This is the authoritative diagnostic for mixed auto-dispatch when different kernels resolve to different backends. |
| [`duckhts_simd_set_backend`](inst/function_catalog/reference.md#duckhts_simd_set_backend)                     | table        | `rduckhts_simd_set_backend`           | Explicitly select the DuckHTS SIMD dispatch policy for this process using a one-row table-function call and return the current dispatch label in a backend column. Use auto for per-kernel runtime dispatch or scalar for a portable baseline; unavailable platform-specific requests such as avx512 on non-AVX-512 CPUs raise an error instead of silently falling back.        |
| [`duckhts_duckdb_type_supported`](inst/function_catalog/reference.md#duckhts_duckdb_type_supported)           | scalar_macro |                                       | Return whether the currently open DuckDB runtime advertises a logical type with the given name through duckdb_types(). This is a catalog-level runtime probe for feature gating SQL/macros across DuckDB versions.                                                                                                                                                               |
| [`duckhts_duckdb_supports_variant`](inst/function_catalog/reference.md#duckhts_duckdb_supports_variant)       | scalar_macro |                                       | Return whether the currently open DuckDB runtime advertises the VARIANT logical type. Use this to gate optional SQL that depends on DuckDB VARIANT support.                                                                                                                                                                                                                      |
| [`duckhts_duckdb_supports_geometry`](inst/function_catalog/reference.md#duckhts_duckdb_supports_geometry)     | scalar_macro |                                       | Return whether the currently open DuckDB runtime advertises the GEOMETRY logical type. Use this to gate optional SQL that depends on DuckDB GEOMETRY support.                                                                                                                                                                                                                    |

### Variant Annotation

| Function                                                                                                        | Kind         | R helper              | Description                                                                                                                                                                                                                                     |
|-----------------------------------------------------------------------------------------------------------------|--------------|-----------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`duckvep_ensembl_regions`](inst/function_catalog/reference.md#duckvep_ensembl_regions)                         | table_macro  |                       | Match tiled FASTA sequence to one Ensembl core assembly and assign dense model-local sequence-region ordinals.                                                                                                                                  |
| [`duckvep_ensembl_transcripts`](inst/function_catalog/reference.md#duckvep_ensembl_transcripts)                 | table_macro  |                       | Build validated VEP-116 Ensembl core transcript models from core tables and matching tiled FASTA sequence.                                                                                                                                      |
| [`duckvep_ensembl_regulation_features`](inst/function_catalog/reference.md#duckvep_ensembl_regulation_features) | table_macro  |                       | Prepare VEP-116 RegulatoryFeature and MotifFeature intervals for a DuckVEP model.                                                                                                                                                               |
| [`duckvep_model_receipt`](inst/function_catalog/reference.md#duckvep_model_receipt)                             | table_macro  |                       | Create a deterministic provenance receipt and semantic hash for prepared DuckVEP model relations.                                                                                                                                               |
| [`duckvep_model_load`](inst/function_catalog/reference.md#duckvep_model_load)                                   | table        |                       | Load a validated immutable consequence model under a name in the current DuckDB database; return one TRUE row.                                                                                                                                  |
| [`duckvep_model_drop`](inst/function_catalog/reference.md#duckvep_model_drop)                                   | scalar       |                       | Remove a named resident DuckVEP consequence model and release its transcript and regulation-feature interval indexes, sequences, and cached worker state. Returns FALSE when the name is absent or the model is in use by an annotation vector. |
| [`duckvep_allele_geometry`](inst/function_catalog/reference.md#duckvep_allele_geometry)                         | scalar       |                       | Separate uploaded, VEP-116 feature and minimized-edit geometry for one literal biallelic allele.                                                                                                                                                |
| [`duckvep_transcript_projection`](inst/function_catalog/reference.md#duckvep_transcript_projection)             | table_macro  |                       | Project independent literal alleles and existing DuckVEP annotations into typed, unshifted VEP-116 transcript display fields.                                                                                                                   |
| [`duckvep_repeat_sequence`](inst/function_catalog/reference.md#duckvep_repeat_sequence)                         | scalar_macro |                       | Expand a caller-asserted exact ordered repeat description into a bounded literal sequence.                                                                                                                                                      |
| [`duckvep_breakend_geometry`](inst/function_catalog/reference.md#duckvep_breakend_geometry)                     | scalar       |                       | Parse one raw VCF 4.5 breakend ALT into mate coordinates, orientation and retained replacement sequence.                                                                                                                                        |
| [`duckvep_haplotypes`](inst/function_catalog/reference.md#duckvep_haplotypes)                                   | table        | `rduckhts_haplotypes` | Replay literal phased CDS/protein paths with carriers, source contributors, coding blocks, aligned differences and optional protein HGVS.                                                                                                       |
| [`duckvep_phase_call`](inst/function_catalog/reference.md#duckvep_phase_call)                                   | scalar       |                       | Assign decoded GT/PS allele slots to haplotype lanes under strict or pinned VEP-116 phase policy.                                                                                                                                               |
| [`duckvep_annotate`](inst/function_catalog/reference.md#duckvep_annotate)                                       | table        |                       | Annotate independent literal alleles, exact typed structural events and paired breakends against a resident VEP-116-compatible model.                                                                                                           |
| [`duckvep_so_terms`](inst/function_catalog/reference.md#duckvep_so_terms)                                       | table        |                       | Return VEP-116 Sequence Ontology terms, consequence-mask bits, impact, severity rank and evaluator tier.                                                                                                                                        |

### Readers

| Function                                                                              | Kind         | R helper                                                                                                                                                               | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
|---------------------------------------------------------------------------------------|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`read_bcf`](inst/function_catalog/reference.md#read_bcf)                             | table        | `rduckhts_bcf`                                                                                                                                                         | Read VCF/BCF with header-typed INFO/FORMAT, typed CSQ/ANN/BCSQ annotations, sample selection and optional tidy sample rows.                                                                                                                                                                                                                                                                                                                                                                                  |
| [`read_geno`](inst/function_catalog/reference.md#read_geno)                           | table        | `rduckhts_geno`                                                                                                                                                        | Read one row per VCF/BCF record with typed arbitrary-ploidy GT/PS calls, selected FORMAT fields and optional original VCF genotype text.                                                                                                                                                                                                                                                                                                                                                                     |
| [`read_bcf_samples`](inst/function_catalog/reference.md#read_bcf_samples)             | table        | `rduckhts_bcf_samples`                                                                                                                                                 | Read the typed VCF/BCF sample catalog as sample_index UINTEGER and sample_name VARCHAR without reading records. Indices are zero-based positions in the original header, remain stable under selection and join read_geno calls from the same unchanged file. NULL or ‘-’ selects all; an empty string selects none; comma-separated names include samples; ‘^’ excludes them. Names are validated by HTSlib, selected rows retain header order, and unknown names error.                                    |
| [`read_bam`](inst/function_catalog/reference.md#read_bam)                             | table        | `rduckhts_bam`                                                                                                                                                         | Read SAM/BAM/CRAM alignments with optional typed SAM tags, auxiliary maps and packed sequence, quality or CIGAR output.                                                                                                                                                                                                                                                                                                                                                                                      |
| [`read_fasta`](inst/function_catalog/reference.md#read_fasta)                         | table        | `rduckhts_fasta`                                                                                                                                                       | Read full FASTA records or indexed regions with text or packed sequence output.                                                                                                                                                                                                                                                                                                                                                                                                                              |
| [`read_bed`](inst/function_catalog/reference.md#read_bed)                             | table        | `rduckhts_bed`                                                                                                                                                         | Read BED3-BED12 interval files with canonical typed columns and optional tabix-backed region filtering. scan_mode := ‘sequential’ forces full-file streaming/counting instead of index-backed count paths and is incompatible with region.                                                                                                                                                                                                                                                                   |
| [`fasta_nuc`](inst/function_catalog/reference.md#fasta_nuc)                           | table        | `rduckhts_fasta_nuc`                                                                                                                                                   | Compute bedtools nuc-style nucleotide composition for supplied BED intervals or generated fixed-width bins over a FASTA reference. A failed reference fetch fails the query with the file and zero-based half-open interval; requested intervals are not silently omitted. For bgzipped FASTA, gzi_path may point to an explicit .gzi sidecar when it is not colocated with the FASTA.                                                                                                                       |
| [`read_fastq`](inst/function_catalog/reference.md#read_fastq)                         | table        | `rduckhts_fastq`                                                                                                                                                       | Read single-end, paired-end, or interleaved FASTQ files with optional legacy quality decoding. By default, FASTQ qualities are interpreted as modern Phred+33 input. Use sequence_encoding := ‘nt16’ to return SEQUENCE as UTINYINT\[\] and quality_representation := ‘phred’ to return QUALITY as UTINYINT\[\] instead of VARCHAR. input_quality_encoding accepts ‘phred33’, ‘auto’, ‘phred64’, or ‘solexa64’. scan_mode := ‘sequential’ forces raw streaming/counting instead of index-backed count paths. |
| [`read_bigwig`](inst/function_catalog/reference.md#read_bigwig)                       | table        | `rduckhts_bigwig`                                                                                                                                                      | Read stored BigWig signal intervals as CHROM, START0, END0 and VALUE.                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| [`read_gff`](inst/function_catalog/reference.md#read_gff)                             | table        | `rduckhts_gff`                                                                                                                                                         | Read GFF annotations with optional raw scalar and parsed list/pair attributes, strict GFF3 validation and indexed region selection.                                                                                                                                                                                                                                                                                                                                                                          |
| [`read_gtf`](inst/function_catalog/reference.md#read_gtf)                             | table        | `rduckhts_gtf`                                                                                                                                                         | Read GTF annotations with optional raw scalar and parsed list/pair attributes and indexed region selection.                                                                                                                                                                                                                                                                                                                                                                                                  |
| [`read_tabix`](inst/function_catalog/reference.md#read_tabix)                         | table        | `rduckhts_tabix`                                                                                                                                                       | Read tabix-indexed text with optional header handling, inferred types and region selection.                                                                                                                                                                                                                                                                                                                                                                                                                  |
| [`fasta_index`](inst/function_catalog/reference.md#fasta_index)                       | table        | `rduckhts_fasta_index`                                                                                                                                                 | Build a FASTA index (.fai) and return a single row with columns success (BOOLEAN) and index_path (VARCHAR).                                                                                                                                                                                                                                                                                                                                                                                                  |
| [`hts_union_query`](inst/function_catalog/reference.md#hts_union_query)               | scalar_macro | `rduckhts_bam_multi, rduckhts_bcf_multi, rduckhts_fastq_multi, rduckhts_fasta_multi, rduckhts_bed_multi, rduckhts_tabix_multi, rduckhts_gff_multi, rduckhts_gtf_multi` | Generate a UNION ALL BY NAME query string that reads every file matching a glob pattern through the named reader function. The result includes a ‘filename’ column identifying the source file for each row. Assign to a variable with SET VARIABLE and execute via query(getvariable(…)). Optional params string is appended to each reader call. In R, use the typed rduckhts\_\*\_multi() helpers instead, which accept file vectors with optional per-file parameters and create DuckDB tables directly. |
| [`hts_region_union_query`](inst/function_catalog/reference.md#hts_region_union_query) | scalar_macro |                                                                                                                                                                        | Generate UNION ALL BY NAME SQL over separate per-region scans of one HTS file.                                                                                                                                                                                                                                                                                                                                                                                                                               |

### Converters

| Function                                                                                                    | Kind         | R helper                         | Description                                                                                                                    |
|-------------------------------------------------------------------------------------------------------------|--------------|----------------------------------|--------------------------------------------------------------------------------------------------------------------------------|
| [`duckhts_bcf_convert_parquet_sql`](inst/function_catalog/reference.md#duckhts_bcf_convert_parquet_sql)     | scalar_macro | `rduckhts_bcf_convert_parquet`   | Build COPY SQL for read_bcf() output with Parquet metadata, VCF header text and selected columns, filters or partitions.       |
| [`duckhts_bam_convert_parquet_sql`](inst/function_catalog/reference.md#duckhts_bam_convert_parquet_sql)     | scalar_macro | `rduckhts_bam_convert_parquet`   | Build COPY SQL for read_bam() output with Parquet metadata, SAM header text and selected columns, filters or partitions.       |
| [`duckhts_gff_convert_parquet_sql`](inst/function_catalog/reference.md#duckhts_gff_convert_parquet_sql)     | scalar_macro | `rduckhts_gff_convert_parquet`   | Build COPY SQL for read_gff() output with Parquet metadata, GFF/tabix header text and selected columns, filters or partitions. |
| [`duckhts_tabix_convert_parquet_sql`](inst/function_catalog/reference.md#duckhts_tabix_convert_parquet_sql) | scalar_macro | `rduckhts_tabix_convert_parquet` | Build COPY SQL for read_tabix() output with Parquet metadata, header text and selected columns, filters or partitions.         |

### Coverage

| Function                                                                                  | Kind  | R helper                    | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
|-------------------------------------------------------------------------------------------|-------|-----------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`read_pileup`](inst/function_catalog/reference.md#read_pileup)                           | table | `rduckhts_pileup`           | Construct a region-scoped BAM pileup with one row per covered position, emitting chrom, 1-based position, depth, observed bases, and Phred+33 qualities after SAM flag and MAPQ filtering. This is a compact htslib pileup view, not samtools mpileup text parity.                                                                                                                                                                                                                                        |
| [`bam_bin_counts`](inst/function_catalog/reference.md#bam_bin_counts)                     | table | `rduckhts_bam_bin_counts`   | Count BAM or CRAM read starts into fixed-width bins. Returns one row per bin across the selected contig span, including zero-count bins, with total, forward, and reverse counts; `rmdup := 'streaming'` applies the WisecondorX-style larp/larp2 consecutive-position deduplication, `rmdup := 'flag'` drops SAM duplicate-flagged reads, and `stats := 'gc'`, `'mq'`, or `'gc,mq'` adds per-bin pre/post-filter GC and MAPQ sufficient statistics, including reference GC when `reference` is provided. |
| [`duckhts_bam_bed_coverage`](inst/function_catalog/reference.md#duckhts_bam_bed_coverage) | table | `rduckhts_bam_bed_coverage` | Compute samtools coverage-like regional summaries for BAM or CRAM input over a BED target set, returning one row per BED interval with DuckHTS-specific pre/post-filter read counts, covered bases, percentage covered, mean depth, mean baseQ, mean mapQ, and strand-specific post-filter summaries in read mode. Indexed BAM/CRAM input is required in the current implementation. decompression_threads controls htslib worker threads for BAM/CRAM decoding; use 0 to disable them.                   |
| [`duckhts_mosdepth`](inst/function_catalog/reference.md#duckhts_mosdepth)                 | table | `rduckhts_mosdepth`         | Write mosdepth-compatible coverage files from indexed BAM/CRAM.                                                                                                                                                                                                                                                                                                                                                                                                                                           |

### Intervals

| Function                                                                                                  | Kind   | R helper | Description                                                                                                                                                                                                                                                                                                                                                                                                                  |
|-----------------------------------------------------------------------------------------------------------|--------|----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`duckhts_cgranges_create`](inst/function_catalog/reference.md#duckhts_cgranges_create)                   | scalar |          | Create an empty session-scoped cgranges registry entry that can be populated with intervals and finalized for overlap queries.                                                                                                                                                                                                                                                                                               |
| [`duckhts_cgranges_add`](inst/function_catalog/reference.md#duckhts_cgranges_add)                         | scalar |          | Append an interval to a session-scoped cgranges registry entry before finalization. Labels may be BIGINT-like, DOUBLE, VARCHAR, or BOOLEAN.                                                                                                                                                                                                                                                                                  |
| [`duckhts_cgranges_index`](inst/function_catalog/reference.md#duckhts_cgranges_index)                     | scalar |          | Finalize a populated cgranges registry entry and build its immutable overlap index for subsequent queries.                                                                                                                                                                                                                                                                                                                   |
| [`duckhts_cgranges_destroy`](inst/function_catalog/reference.md#duckhts_cgranges_destroy)                 | scalar |          | Destroy a session-scoped cgranges registry entry and release its indexed interval storage when it is not in active use.                                                                                                                                                                                                                                                                                                      |
| [`duckhts_cgranges_from_query`](inst/function_catalog/reference.md#duckhts_cgranges_from_query)           | scalar |          | Execute a SQL query on an extension-owned DuckDB connection, append its interval rows into a session-scoped cgranges registry entry, and leave the populated index ready for explicit finalization with duckhts_cgranges_index(…).                                                                                                                                                                                           |
| [`duckhts_cgranges_from_table`](inst/function_catalog/reference.md#duckhts_cgranges_from_table)           | scalar |          | Reserved convenience constructor for bulk cgranges population from a table name. The current implementation is intentionally deferred and directs callers to duckhts_cgranges_from_query(…).                                                                                                                                                                                                                                 |
| [`duckhts_cgranges_has_overlap`](inst/function_catalog/reference.md#duckhts_cgranges_has_overlap)         | scalar |          | Vectorized scalar predicate for streaming provider rows through a finalized session-scoped cgranges index. Returns TRUE when the query interval overlaps at least one indexed interval, or when mode = ‘contain’ and it fully contains at least one indexed interval; NULL inputs return NULL.                                                                                                                               |
| [`duckhts_cgranges_count_overlaps`](inst/function_catalog/reference.md#duckhts_cgranges_count_overlaps)   | scalar |          | Vectorized scalar overlap counter for streaming provider rows through a finalized session-scoped cgranges index. Returns the number of indexed intervals that overlap the query interval, or with mode = ‘contain’ the number fully contained by it; NULL inputs return NULL.                                                                                                                                                |
| [`duckhts_cgranges_overlaps_list`](inst/function_catalog/reference.md#duckhts_cgranges_overlaps_list)     | scalar |          | Vectorized scalar overlap expander for streaming provider rows through a finalized session-scoped cgranges index. Returns a LIST of hit STRUCTs that can be expanded with UNNEST, preserving provider columns while emitting one row per matching indexed interval. Because scalar return types are fixed, labels are returned as text with label_type describing the original cgranges label kind; NULL inputs return NULL. |
| [`duckhts_cgranges_overlaps`](inst/function_catalog/reference.md#duckhts_cgranges_overlaps)               | table  |          | Query a finalized session-scoped cgranges registry entry and return one row per overlapping or containing indexed interval, preserving the original label type and interval coordinates.                                                                                                                                                                                                                                     |
| [`duckhts_cgranges_overlaps_bulk`](inst/function_catalog/reference.md#duckhts_cgranges_overlaps_bulk)     | table  |          | Run a SQL query that yields overlap probes, stream those rows through a finalized session-scoped cgranges registry entry, and return one row per matching indexed interval. The probe query runs on the extension-owned helper connection, so it must reference regular tables/views rather than connection-local temp tables. When query_row_id_col is omitted, query_row_id defaults to the 1-based probe row ordinal.     |
| [`regionkey`](inst/function_catalog/reference.md#regionkey)                                               | scalar |          | Encode a genomic interval as an official RegionKey-compatible 64-bit unsigned integer. Start and end use 0-based half-open interval semantics, matching BED-style coordinates; strand accepts -1, 0, or 1.                                                                                                                                                                                                                   |
| [`regionkey_hex`](inst/function_catalog/reference.md#regionkey_hex)                                       | scalar |          | Render a RegionKey as its lowercase 16-character hexadecimal string representation.                                                                                                                                                                                                                                                                                                                                          |
| [`parse_regionkey_hex`](inst/function_catalog/reference.md#parse_regionkey_hex)                           | scalar |          | Parse a 16-character hexadecimal RegionKey string back into its UBIGINT code. Invalid or non-hex strings return NULL.                                                                                                                                                                                                                                                                                                        |
| [`encode_regionkey`](inst/function_catalog/reference.md#encode_regionkey)                                 | scalar |          | Encode the raw upstream RegionKey fields directly: chromosome code, 0-based start, 0-based end, and strand code (0 = unknown, 1 = +, 2 = -).                                                                                                                                                                                                                                                                                 |
| [`extract_regionkey_chrom`](inst/function_catalog/reference.md#extract_regionkey_chrom)                   | scalar |          | Extract the raw upstream RegionKey chromosome code.                                                                                                                                                                                                                                                                                                                                                                          |
| [`extract_regionkey_startpos`](inst/function_catalog/reference.md#extract_regionkey_startpos)             | scalar |          | Extract the raw upstream RegionKey 0-based start position.                                                                                                                                                                                                                                                                                                                                                                   |
| [`extract_regionkey_endpos`](inst/function_catalog/reference.md#extract_regionkey_endpos)                 | scalar |          | Extract the raw upstream RegionKey 0-based end position.                                                                                                                                                                                                                                                                                                                                                                     |
| [`extract_regionkey_strand`](inst/function_catalog/reference.md#extract_regionkey_strand)                 | scalar |          | Extract the raw upstream RegionKey strand code (0 = unknown, 1 = +, 2 = -).                                                                                                                                                                                                                                                                                                                                                  |
| [`decode_regionkey`](inst/function_catalog/reference.md#decode_regionkey)                                 | scalar |          | Decode a RegionKey into its raw upstream numeric fields: chrom_code, start, end, and strand_code.                                                                                                                                                                                                                                                                                                                            |
| [`reverse_regionkey`](inst/function_catalog/reference.md#reverse_regionkey)                               | scalar |          | Decode a RegionKey into a STRUCT with chrom, chrom_code, start, end, strand, and strand_code.                                                                                                                                                                                                                                                                                                                                |
| [`extend_regionkey`](inst/function_catalog/reference.md#extend_regionkey)                                 | scalar |          | Extend a RegionKey interval by a fixed number of bases on both sides, clamping to the official 28-bit RegionKey position range.                                                                                                                                                                                                                                                                                              |
| [`are_overlapping_regions`](inst/function_catalog/reference.md#are_overlapping_regions)                   | scalar |          | Return TRUE when two explicit 0-based half-open intervals overlap on the same canonical chromosome.                                                                                                                                                                                                                                                                                                                          |
| [`are_overlapping_region_regionkey`](inst/function_catalog/reference.md#are_overlapping_region_regionkey) | scalar |          | Return TRUE when a 0-based half-open interval overlaps the supplied RegionKey interval.                                                                                                                                                                                                                                                                                                                                      |
| [`are_overlapping_regionkeys`](inst/function_catalog/reference.md#are_overlapping_regionkeys)             | scalar |          | Return TRUE when two RegionKeys overlap.                                                                                                                                                                                                                                                                                                                                                                                     |

### Quality Control

| Function                                                                  | Kind      | R helper | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
|---------------------------------------------------------------------------|-----------|----------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`duckhts_fastq_qc`](inst/function_catalog/reference.md#duckhts_fastq_qc) | aggregate |          | Aggregate canonical sequence and Phred+33 quality strings directly into exact read/base/Q20/Q30/Q40, nucleotide, quality-sum, and per-cycle sufficient statistics. The nested cycles list supports mean-quality, nucleotide-content, GC, and read-length curves without expanding one SQL row per base. Rows with any NULL input are ignored. Per-cycle state defaults to at most 1,048,576 cycles; pass a constant max_cycles per aggregate group to choose a larger explicit limit, up to 16,777,216. |

### Sample Identity

| Function                                                                                                              | Kind         | R helper                                  | Description                                                                                      |
|-----------------------------------------------------------------------------------------------------------------------|--------------|-------------------------------------------|--------------------------------------------------------------------------------------------------|
| [`duckhts_somalier_panel_sha256`](inst/function_catalog/reference.md#duckhts_somalier_panel_sha256)                   | scalar_macro |                                           | Derive a stable SHA-256 identity for an ordered biallelic sample-fingerprinting panel.           |
| [`duckhts_somalier_frequency_sha256`](inst/function_catalog/reference.md#duckhts_somalier_frequency_sha256)           | scalar_macro |                                           | Derive a stable identity for panel-aligned population-B allele frequencies.                      |
| [`duckhts_somalier_classify`](inst/function_catalog/reference.md#duckhts_somalier_classify)                           | scalar       |                                           | Classify one measured A/B/other count tuple for Somalier-derived autosomal relatedness.          |
| [`duckhts_somalier_prepare_sketches`](inst/function_catalog/reference.md#duckhts_somalier_prepare_sketches)           | table_macro  | `rduckhts_somalier_sketches`              | Build one panel-verified packed relatedness sketch per sample from typed count evidence.         |
| [`duckhts_somalier_verify_sketches`](inst/function_catalog/reference.md#duckhts_somalier_verify_sketches)             | scalar_macro |                                           | Verify persisted relatedness sketches against their retained raw count evidence.                 |
| [`duckhts_somalier_relatedness`](inst/function_catalog/reference.md#duckhts_somalier_relatedness)                     | scalar       | `rduckhts_somalier_relatedness`           | Compute fused Somalier-derived relatedness and concordance statistics for two prepared sketches. |
| [`duckhts_somalier_verify_relatedness`](inst/function_catalog/reference.md#duckhts_somalier_verify_relatedness)       | scalar       |                                           | Verify a typed relatedness result against its two sealed sketches.                               |
| [`duckhts_somalier_charr`](inst/function_catalog/reference.md#duckhts_somalier_charr)                                 | table_macro  | `rduckhts_somalier_charr`                 | Estimate per-sample contamination with a bounded Somalier-derived CHARR reduction.               |
| [`duckhts_somalier_matched_contamination`](inst/function_catalog/reference.md#duckhts_somalier_matched_contamination) | table_macro  | `rduckhts_somalier_matched_contamination` | Estimate directional contamination for explicitly selected receiver/anchor sample pairs.         |

### Metadata

| Function                                                                                    | Kind        | R helper                           | Description                                                                                                                                                                                                                                                                |
|---------------------------------------------------------------------------------------------|-------------|------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`detect_quality_encoding`](inst/function_catalog/reference.md#detect_quality_encoding)     | table       | `rduckhts_detect_quality_encoding` | Inspect a FASTQ file’s observed quality ASCII range and report compatible legacy encodings with a heuristic guessed encoding.                                                                                                                                              |
| [`duckhts_samtools_idxstats`](inst/function_catalog/reference.md#duckhts_samtools_idxstats) | table       | `rduckhts_samtools_idxstats`       | Write samtools idxstats-compatible TAB-delimited output for BAM, CRAM, or SAM input. Indexed BAM uses `hts_idx_get_stat(...)` for the fast path; CRAM, SAM, and unindexed BAM fall back to a full scan while preserving samtools-style contig rows plus the final `*` row. |
| [`read_hts_header`](inst/function_catalog/reference.md#read_hts_header)                     | table       | `rduckhts_hts_header`              | Inspect HTS headers in parsed, raw, or combined form across supported formats. Raw VCF/BCF mode includes the final `#CHROM` sample header line so the returned text is suitable for Parquet metadata and future VCF/BCF regeneration.                                      |
| [`read_hts_index`](inst/function_catalog/reference.md#read_hts_index)                       | table       | `rduckhts_hts_index`               | Inspect high-level HTS index metadata such as sequence names and mapped counts.                                                                                                                                                                                            |
| [`read_hts_index_spans`](inst/function_catalog/reference.md#read_hts_index_spans)           | table       | `rduckhts_hts_index_spans`         | Expand index metadata into span and chunk rows suitable for low-level index inspection.                                                                                                                                                                                    |
| [`read_hts_index_raw`](inst/function_catalog/reference.md#read_hts_index_raw)               | table_macro | `rduckhts_hts_index_raw`           | Return the raw on-disk HTS index blob together with basic identifying metadata.                                                                                                                                                                                            |

### Compression

| Function                                                | Kind  | R helper           | Description                                                                           |
|---------------------------------------------------------|-------|--------------------|---------------------------------------------------------------------------------------|
| [`bgzip`](inst/function_catalog/reference.md#bgzip)     | table | `rduckhts_bgzip`   | Compress a plain file to BGZF and return the created output path and byte counts.     |
| [`bgunzip`](inst/function_catalog/reference.md#bgunzip) | table | `rduckhts_bgunzip` | Decompress a BGZF-compressed file and return the created output path and byte counts. |

### Indexing

| Function                                                        | Kind  | R helper               | Description                                                                                        |
|-----------------------------------------------------------------|-------|------------------------|----------------------------------------------------------------------------------------------------|
| [`bam_index`](inst/function_catalog/reference.md#bam_index)     | table | `rduckhts_bam_index`   | Build a BAM or CRAM index and report the written index path and format.                            |
| [`bcf_index`](inst/function_catalog/reference.md#bcf_index)     | table | `rduckhts_bcf_index`   | Build a TBI or CSI index for a VCF or BCF file and report the written index path and format.       |
| [`tabix_index`](inst/function_catalog/reference.md#tabix_index) | table | `rduckhts_tabix_index` | Build a tabix index for a BGZF-compressed text file using a preset or explicit coordinate columns. |

### Variants

| Function                                                                                    | Kind        | R helper                 | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
|---------------------------------------------------------------------------------------------|-------------|--------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`variantkey`](inst/function_catalog/reference.md#variantkey)                               | scalar      |                          | Encode a normalized biallelic variant as an official VariantKey-compatible 64-bit unsigned integer. This DuckHTS wrapper accepts 1-based VCF/DuckHTS POS to match bcftools `%VKX` / `+add-variantkey`, internally converts to the upstream 0-based field, and preserves the official hashed nonreversible mode for large, ambiguous, and symbolic REF/ALT strings. Only CHROM, POS, REF, and ALT are encoded; END, SVLEN, mate breakend coordinates, and other SV metadata are not.        |
| [`variantkey_hex`](inst/function_catalog/reference.md#variantkey_hex)                       | scalar      |                          | Render a VariantKey as its lowercase 16-character hexadecimal string representation.                                                                                                                                                                                                                                                                                                                                                                                                       |
| [`parse_variantkey_hex`](inst/function_catalog/reference.md#parse_variantkey_hex)           | scalar      |                          | Parse a 16-character hexadecimal VariantKey string back into its UBIGINT code. Invalid or non-hex strings return NULL.                                                                                                                                                                                                                                                                                                                                                                     |
| [`encode_variantkey`](inst/function_catalog/reference.md#encode_variantkey)                 | scalar      |                          | Encode the raw upstream VariantKey fields directly: chromosome code, 0-based position, and 31-bit REF+ALT code.                                                                                                                                                                                                                                                                                                                                                                            |
| [`extract_variantkey_chrom`](inst/function_catalog/reference.md#extract_variantkey_chrom)   | scalar      |                          | Extract the raw upstream VariantKey chromosome code.                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| [`extract_variantkey_pos`](inst/function_catalog/reference.md#extract_variantkey_pos)       | scalar      |                          | Extract the raw upstream VariantKey 0-based position field.                                                                                                                                                                                                                                                                                                                                                                                                                                |
| [`extract_variantkey_refalt`](inst/function_catalog/reference.md#extract_variantkey_refalt) | scalar      |                          | Extract the raw upstream 31-bit VariantKey REF+ALT code.                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| [`decode_variantkey`](inst/function_catalog/reference.md#decode_variantkey)                 | scalar      |                          | Decode a VariantKey into its raw upstream numeric fields: chrom_code, pos0, and refalt_code.                                                                                                                                                                                                                                                                                                                                                                                               |
| [`reverse_variantkey`](inst/function_catalog/reference.md#reverse_variantkey)               | scalar      |                          | Decode a VariantKey into a STRUCT with chrom, chrom_code, 1-based pos, upstream 0-based pos0, ref, alt, refalt_code, and reversible. For hashed nonreversible keys, reversible is FALSE and ref/alt are returned as NULL because DuckHTS v1 does not ship the optional NRVK lookup sidecar.                                                                                                                                                                                                |
| [`variantkey_range`](inst/function_catalog/reference.md#variantkey_range)                   | scalar      |                          | Return the inclusive minimum and maximum VariantKey bounds for a chromosome plus 1-based VCF position range, suitable for numeric range filtering on precomputed VariantKeys.                                                                                                                                                                                                                                                                                                              |
| [`duckhts_contig_key`](inst/function_catalog/reference.md#duckhts_contig_key)               | scalar      |                          | Return a conservative contig join key by removing one non-empty leading chr prefix case-insensitively and normalizing M/MT to MT. X and Y are uppercased; all other suffixes are preserved. This does not map numeric sex chromosomes, accessions, patches, or alternate loci.                                                                                                                                                                                                             |
| [`bcftools_liftover`](inst/function_catalog/reference.md#bcftools_liftover)                 | scalar      | `rduckhts_liftover`      | Row-oriented liftover kernel intended to mirror bcftools +liftover semantics as closely as possible while returning one STRUCT per input row with fields: src_chrom, src_pos, src_ref, src_alt, dest_chrom, dest_pos, dest_end, dest_ref, dest_alt, mapped, reverse_complemented, swap, reject_reason, and note. Set no_left_align := true to skip post-liftover left-alignment of lifted indels (mirrors –no-left-align in bcftools +liftover).                                           |
| [`duckdb_liftover`](inst/function_catalog/reference.md#duckdb_liftover)                     | table_macro | `rduckhts_liftover`      | DuckDB-specific wrapper over bcftools_liftover that takes either a table name or a derived-table expression plus column-name strings for chrom/pos/ref/alt and returns the lifted table. The no_left_align parameter mirrors –no-left-align in bcftools +liftover.                                                                                                                                                                                                                         |
| [`bcftools_norm_row`](inst/function_catalog/reference.md#bcftools_norm_row)                 | scalar      |                          | Normalize one variant against FASTA with bcftools/vt-style left alignment.                                                                                                                                                                                                                                                                                                                                                                                                                 |
| [`duckhts_bcftools_norm`](inst/function_catalog/reference.md#duckhts_bcftools_norm)         | table_macro | `rduckhts_bcftools_norm` | Normalize variants from a table or derived-table expression while preserving input columns.                                                                                                                                                                                                                                                                                                                                                                                                |
| [`bcftools_score`](inst/function_catalog/reference.md#bcftools_score)                       | table       | `rduckhts_score`         | Compute polygenic scores from genotype VCF/BCF and summary statistics using bcftools +score dosage semantics.                                                                                                                                                                                                                                                                                                                                                                              |
| [`bcftools_munge_row`](inst/function_catalog/reference.md#bcftools_munge_row)               | scalar      |                          | Normalize one summary-statistics row into GWAS-VCF-style fields (chrom/pos/ref/alt/effect metrics), resolving REF/ALT orientation against a FASTA reference and applying swap-aware sign/frequency/count transforms. The output flag `alleles_swapped` means REF/ALT orientation was swapped to match the FASTA reference.                                                                                                                                                                 |
| [`duckdb_munge`](inst/function_catalog/reference.md#duckdb_munge)                           | table_macro | `rduckhts_munge`         | DuckDB macro wrapper over bcftools_munge_row that maps source columns (via preset or explicit map) and returns normalized GWAS-VCF-style rows with lean outputs and explicit `alleles_swapped` semantics. Output columns: chrom, pos, id, ref, alt, alleles_swapped, filter, ns, ez, nc, es, se, lp, af, ac, ne (16 columns). For METAL meta-analysis output with SI/I2/CQ/ED columns, use duckdb_munge_metal.                                                                             |
| [`duckdb_munge_metal`](inst/function_catalog/reference.md#duckdb_munge_metal)               | table_macro | `rduckhts_munge`         | Extended munge macro with METAL meta-analysis output columns. Same as duckdb_munge but additionally emits: si (imputation info, from INFO input), i2 (Cochran’s I² heterogeneity, from HET_I2), cq (Cochran’s Q -log10 p, from HET_LP or -log10(HET_P)), and ed (effect direction string, from DIRE; +/- flipped on allele swap). The R wrapper rduckhts_munge() auto-dispatches to this macro when metal keys (INFO, HET_I2, HET_P, HET_LP, DIRE) are present in the resolved column map. |

### Sequence UDFs

| Function                                                                | Kind   | R helper | Description                                                                                                                                                                                                                                                                                                                                                                                                             |
|-------------------------------------------------------------------------|--------|----------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`seq_revcomp`](inst/function_catalog/reference.md#seq_revcomp)         | scalar |          | Compute the reverse complement of a DNA sequence using A, C, G, T, and N bases. Overloaded: accepts either a VARCHAR text sequence (returns VARCHAR) or a UTINYINT\[\] of htslib nt16 codes as produced by read_bam(sequence_encoding := ‘nt16’) (returns UTINYINT\[\]); the nt16 overload is bit-identical to the text path after decoding, so BAM pipelines can reverse-complement without leaving the nt16 encoding. |
| [`seq_canonical`](inst/function_catalog/reference.md#seq_canonical)     | scalar |          | Return the lexicographically smaller of a sequence and its reverse complement. Overloaded: accepts either a VARCHAR text sequence (returns VARCHAR) or a UTINYINT\[\] of htslib nt16 codes as produced by read_bam(sequence_encoding := ‘nt16’) (returns UTINYINT\[\]); the nt16 overload compares by decoded base order and is bit-identical to the text path after decoding.                                          |
| [`seq_hash_2bit`](inst/function_catalog/reference.md#seq_hash_2bit)     | scalar |          | Encode a short DNA sequence as a 2-bit unsigned integer hash. Overloaded to also accept a UTINYINT\[\] of htslib nt16 codes (from read_bam(sequence_encoding := ‘nt16’)); non-ACGT codes yield NULL, bit-identical to the text path.                                                                                                                                                                                    |
| [`seq_encode_4bit`](inst/function_catalog/reference.md#seq_encode_4bit) | scalar |          | Encode an IUPAC DNA sequence as a list of 4-bit base codes, preserving ambiguity symbols including N.                                                                                                                                                                                                                                                                                                                   |
| [`seq_decode_4bit`](inst/function_catalog/reference.md#seq_decode_4bit) | scalar |          | Decode a list of 4-bit IUPAC DNA base codes back into a sequence string.                                                                                                                                                                                                                                                                                                                                                |
| [`seq_gc_content`](inst/function_catalog/reference.md#seq_gc_content)   | scalar |          | Compute GC fraction for a DNA sequence as a value between 0 and 1. Overloaded: accepts either a VARCHAR text sequence or a UTINYINT\[\] of htslib nt16 codes as produced by read_bam(sequence_encoding := ‘nt16’); the nt16 overload classifies codes directly and is bit-identical to the text path, so BAM pipelines can compute GC without decoding sequences back to text.                                          |
| [`seq_kmers`](inst/function_catalog/reference.md#seq_kmers)             | table  |          | Expand a sequence into positional k-mers with optional canonicalization.                                                                                                                                                                                                                                                                                                                                                |

### SAM Flag UDFs

| Function                                                                                                          | Kind   | R helper | Description                                                                                                                                                                        |
|-------------------------------------------------------------------------------------------------------------------|--------|----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`sam_flag_bits`](inst/function_catalog/reference.md#sam_flag_bits)                                               | scalar |          | Decode a SAM flag into a struct of boolean bit fields using explicit SAM-oriented names such as `is_paired`, `is_proper_pair`, `is_next_segment_unmapped`, and `is_supplementary`. |
| [`sam_flag_has`](inst/function_catalog/reference.md#sam_flag_has)                                                 | scalar |          | Test whether any bits from the provided SAM flag mask are set in a flag value.                                                                                                     |
| [`is_forward_aligned`](inst/function_catalog/reference.md#is_forward_aligned)                                     | scalar |          | Test whether a mapped segment is aligned to the forward strand. Returns `NULL` for unmapped segments because SAM flag `0x10` does not define genomic strand when `0x4` is set.     |
| [`is_paired`](inst/function_catalog/reference.md#is_paired)                                                       | scalar |          | Test whether the SAM flag indicates that the template has multiple segments in sequencing (`0x1`).                                                                                 |
| [`is_proper_pair`](inst/function_catalog/reference.md#is_proper_pair)                                             | scalar |          | Test whether the SAM flag indicates that each segment is properly aligned according to the aligner (`0x2`).                                                                        |
| [`is_unmapped`](inst/function_catalog/reference.md#is_unmapped)                                                   | scalar |          | Test whether the read itself is unmapped according to the SAM flag.                                                                                                                |
| [`is_next_segment_unmapped`](inst/function_catalog/reference.md#is_next_segment_unmapped)                         | scalar |          | Test whether the next segment in the template is flagged as unmapped (`0x8`).                                                                                                      |
| [`is_reverse_complemented`](inst/function_catalog/reference.md#is_reverse_complemented)                           | scalar |          | Test whether `SEQ` is stored reverse complemented (`0x10`); for mapped reads this corresponds to reverse-strand alignment.                                                         |
| [`is_next_segment_reverse_complemented`](inst/function_catalog/reference.md#is_next_segment_reverse_complemented) | scalar |          | Test whether `SEQ` of the next segment in the template is stored reverse complemented (`0x20`).                                                                                    |
| [`is_first_segment`](inst/function_catalog/reference.md#is_first_segment)                                         | scalar |          | Test whether the read is marked as the first segment in the template.                                                                                                              |
| [`is_last_segment`](inst/function_catalog/reference.md#is_last_segment)                                           | scalar |          | Test whether the read is marked as the last segment in the template.                                                                                                               |
| [`is_secondary`](inst/function_catalog/reference.md#is_secondary)                                                 | scalar |          | Test whether the alignment is marked as secondary.                                                                                                                                 |
| [`is_qc_fail`](inst/function_catalog/reference.md#is_qc_fail)                                                     | scalar |          | Test whether the read failed vendor or pipeline quality checks.                                                                                                                    |
| [`is_duplicate`](inst/function_catalog/reference.md#is_duplicate)                                                 | scalar |          | Test whether the alignment is flagged as a duplicate.                                                                                                                              |
| [`is_supplementary`](inst/function_catalog/reference.md#is_supplementary)                                         | scalar |          | Test whether the alignment is marked as supplementary.                                                                                                                             |

### CIGAR Utils

| Function                                                                                      | Kind   | R helper | Description                                                                                                                                                                                                                                                                                   |
|-----------------------------------------------------------------------------------------------|--------|----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`cigar_has_soft_clip`](inst/function_catalog/reference.md#cigar_has_soft_clip)               | scalar |          | Test whether a CIGAR string contains any soft-clipped segment (`S`). Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                                                |
| [`cigar_has_hard_clip`](inst/function_catalog/reference.md#cigar_has_hard_clip)               | scalar |          | Test whether a CIGAR string contains any hard-clipped segment (`H`). Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                                                |
| [`cigar_left_soft_clip`](inst/function_catalog/reference.md#cigar_left_soft_clip)             | scalar |          | Return the left-end soft-clipped length from a CIGAR string, or zero if the alignment does not start with `S`. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.      |
| [`cigar_right_soft_clip`](inst/function_catalog/reference.md#cigar_right_soft_clip)           | scalar |          | Return the right-end soft-clipped length from a CIGAR string, or zero if the alignment does not end with `S`. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.       |
| [`cigar_query_length`](inst/function_catalog/reference.md#cigar_query_length)                 | scalar |          | Return the query-consuming length from a CIGAR string, counting `M`, `I`, `S`, `=`, and `X`. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                        |
| [`cigar_aligned_query_length`](inst/function_catalog/reference.md#cigar_aligned_query_length) | scalar |          | Return the aligned query length from a CIGAR string, counting `M`, `=`, and `X` but excluding clips and insertions. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path. |
| [`cigar_reference_length`](inst/function_catalog/reference.md#cigar_reference_length)         | scalar |          | Return the reference-consuming length from a CIGAR string, counting `M`, `D`, `N`, `=`, and `X`. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                    |
| [`cigar_has_op`](inst/function_catalog/reference.md#cigar_has_op)                             | scalar |          | Test whether a CIGAR string contains at least one instance of the requested operator. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                               |

</details>

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

### BigWig signals

The package bundles libBigWig’s upstream test file so the wrapper,
stored zero-based half-open coordinates, and htslib-style multi-region
semantics are executable offline. `rduckhts_bigwig()` materializes a
table; pass `table_name = NULL` to create the `bigwig_data` view
instead.

``` r
bigwig_path <- system.file(
  "extdata", "libbigwig_test.bw", package = "Rduckhts", mustWork = TRUE
)
rduckhts_bigwig(
  con, "bigwig_intervals", bigwig_path,
  region = c("1:1-150", "10:201-300"),
  overwrite = TRUE
)
dbGetQuery(
  con,
  paste(
    "SELECT CHROM, START0, END0, round(VALUE::DOUBLE, 1) AS VALUE",
    "FROM bigwig_intervals ORDER BY CHROM, START0"
  )
)
#>   CHROM START0 END0 VALUE
#> 1     1      0    1   0.1
#> 2     1      1    2   0.2
#> 3     1      2    3   0.3
#> 4     1    100  150   1.4
#> 5    10    200  300   2.0
```

### Variant normalization

`rduckhts_bcftools_norm()` wraps the bundled
`duckhts_bcftools_norm(...)` table macro and keeps the original input
columns alongside normalized position/reference/ALT outputs.

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

### DuckVEP consequence and HGVS annotation

#### Design and validation

DuckVEP annotates a narrow event relation against an immutable named
transcript model. Production models come from receipted Ensembl
relations; this compact non-coding transcript keeps the README
deterministic while exercising the same public DBI surface. The bundled
extension shares one read-only model across DuckDB workers and gives
each worker private mutable annotation state. Compact masks/codes are
the high-cardinality interface; rich consequence text and
independent-event HGVS are optional projections of the same candidate
sweep.

DuckVEP returns a relation rather than owning a closed annotation-cache
format. Any source DuckDB can scan–local or S3 Parquet, a DuckLake
snapshot, VCF/BCF, tabix text, or an attached database–can provide exact
allele, interval, gene, or disease annotations through ordinary joins.
The root [DuckHTS
README](https://github.com/RGenomicsETL/duckhts#duckvep-a-resident-ensembl-vep-consequence-engine)
shows the complete real HG002 WGS + Ensembl 116 +
ClinVar/ClinvArbitration + AlphaMissense + gnomAD + interval workflow
and its materialized timing evidence; the illustrated [*DuckVEP: the
fastest Ensembl VEP-compatible consequence predictor in the
West?*](https://github.com/RGenomicsETL/duckhts/blob/main/benchmarks/benchmark_duckvep_fastvep.md)
report presents the whole-genome speed, memory, and VEP-conformance
evidence. It also directly benchmarks the same dated ClinVar payload
through FastVEP fastSA and DuckDB’s typed collision-safe join, then
covers AlphaMissense, assembly-correct REVEL, BigWig conservation,
algorithm design, and the layered fuzz/differential/sanitizer
infrastructure. The small example below remains only the CRAN/offline
executable check.

The release gate combines pure-C property/statistical tests, SQL and R
end-to-end tests, source-labelled Ensembl fixtures, and fail-closed
differentials against executable VEP 116. Complete ClinVar/GIAB and
selected GRCh37 campaigns retain every compared allele/object pair. The
large-campaign `{targets}`/`blit`/micromamba workflow lives in the
DuckHTS source repository; it is not reimplemented in the R package.

#### Current scope

The current surface covers independent small variants, exact span SVs,
paired BNDs, regulation/motif consequences, NMD, and independent-event
HGVSc/HGVSn/HGVSp. Combined phased haplotype consequences/HGVS,
imprecise-SV confidence intervals, STR-specific rules, and arbitrary
producer-specific symbolic interpretation remain explicit follow-up
work. Supplementary population and clinical annotations remain typed
DuckDB/Parquet/DuckLake relations joined after consequence filtering.

The resident kernel accepts any transcript catalog that satisfies its
checked model contract. The bundled builder currently compiles the
Ensembl core transcript set, not VEP’s `--refseq` or `--merged` sets. It
retains GENCODE Basic and Primary flags for SQL filtering but does not
prefilter the model to either set. The builder withholds sequence for
source rows marked with transcript-level sequence corrections. A model
with inserted or deleted transcript bases relative to its genomic exons
also needs a richer coordinate map than the current resident interface.
Mitochondrial genetic codes and ordinary MT coordinates are supported;
coordinates that cross a circular sequence origin are not.

#### Minimal R example

``` r
duckvep_reference <- system.file(
  "extdata", "fixture_ref.fa", package = "Rduckhts", mustWork = TRUE
)

dbExecute(con, "
  CREATE OR REPLACE TABLE readme_duckvep_regions AS
  SELECT 1::UINTEGER AS seq_region, 50000::UBIGINT AS sequence_length,
         '11'::VARCHAR AS seq_region_name
")
#> [1] 1
dbExecute(con, "
  CREATE OR REPLACE TABLE readme_duckvep_transcripts AS
  SELECT 0::UINTEGER AS transcript_index, 1::UINTEGER AS seq_region,
         100::UBIGINT AS transcript_start, 150::UBIGINT AS transcript_end,
         1::TINYINT AS strand, 0::UINTEGER AS gene_index,
         0::UBIGINT AS transcript_flags, NULL::UBIGINT AS cds_start,
         NULL::UBIGINT AS cds_end, NULL::BLOB AS cds_sequence,
         NULL::UTINYINT AS codon_table, NULL::BLOB AS pre_cds_sequence,
         NULL::BLOB AS post_cds_sequence
")
#> [1] 1
dbExecute(con, "
  CREATE OR REPLACE TABLE readme_duckvep_exons AS
  SELECT 0::UINTEGER AS transcript_index, 100::UBIGINT AS exon_start,
         150::UBIGINT AS exon_end, 1::UBIGINT AS exon_cdna_start,
         51::UBIGINT AS exon_cdna_end, -1::TINYINT AS phase,
         -1::TINYINT AS end_phase
")
#> [1] 1

duckvep_model_queries <- c(
  "SELECT * FROM readme_duckvep_regions ORDER BY seq_region",
  "SELECT * FROM readme_duckvep_transcripts
   ORDER BY seq_region, transcript_start, transcript_index",
  "SELECT * FROM readme_duckvep_exons
   ORDER BY transcript_index, exon_cdna_start"
)
quoted_model <- vapply(
  c("readme", duckvep_model_queries, duckvep_reference),
  function(x) as.character(dbQuoteString(con, x)),
  character(1)
)
dbGetQuery(
  con,
  sprintf(
    "SELECT loaded FROM duckvep_model_load(
       %s, %s, %s, %s, reference_fasta := %s
     )",
    quoted_model[[1]], quoted_model[[2]], quoted_model[[3]],
    quoted_model[[4]], quoted_model[[5]]
  )
)
#>   loaded
#> 1   TRUE

dbExecute(con, "
  CREATE OR REPLACE TABLE readme_duckvep_events AS
  SELECT * FROM (VALUES
    (1::UBIGINT, 1::UINTEGER, 124::UBIGINT, 'A'::VARCHAR, 'G'::VARCHAR,
     NULL::UBIGINT, NULL::VARCHAR, NULL::VARCHAR,
     NULL::UINTEGER, NULL::UBIGINT)
  ) AS e(event_index, seq_region, position, reference, alternate,
         end_position, structural_type, copy_change,
         mate_seq_region, mate_position)
")
#> [1] 1

dbGetQuery(con, "
  SELECT a.event_index, a.transcript_index,
         string_agg(t.consequence, '&' ORDER BY t.severity_rank) AS consequence,
         a.transcript_hgvs, a.protein_hgvs
  FROM duckvep_annotate(
    'readme_duckvep_events', 'readme', hgvs := true,
    upstream_distance := 0, downstream_distance := 0
  ) AS a
  JOIN duckvep_so_terms() AS t
    ON (a.consequence_mask & t.consequence_mask) <> 0
  GROUP BY ALL
  ORDER BY a.event_index, a.transcript_index
")
#>   event_index transcript_index                        consequence
#> 1           1                0 non_coding_transcript_exon_variant
#>   transcript_hgvs protein_hgvs
#> 1         n.25A>G         <NA>

dbGetQuery(con, "SELECT duckvep_model_drop('readme') AS dropped")
#>   dropped
#> 1    TRUE
```

The bundled extension also exposes
`duckvep_transcript_projection(events, annotations, transcripts)`
through DBI. It is the SQL reference for full cDNA/CDS/protein ranges,
transcript-oriented codon and amino-acid strings, exon/intron ordinal
ranges, distance, and CDS quality flags for literal independent alleles.
Supply the same validated prepared model relation used for annotation,
with transcript-ordered nested exons and peptide edits. NULL range
endpoints remain independent; insertions retain an explicit `interbase`
flag. Duplicate annotation rows are preserved. Unavailable coding
sequence leaves codon and amino-acid fields NULL, not `-`.
Canonical/MANE/CCDS/GENCODE/biotype attributes remain relational joins,
not extra hot consequence payload; canonical status is not a VEP quality
flag. The [DuckHTS projection
example](https://github.com/RGenomicsETL/duckhts#canonical-event-relation-and-output-choices)
shows those joins and the literal-allele scope. No additional R wrapper
or reference handle is needed.

### Interval + reference helpers

``` r
bed_path <- system.file("extdata", "targets.bed", package = "Rduckhts")
fai_path <- file.path(tempdir(), "duckhts_readme_00000000000000.fai")
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

`rduckhts_bam_bin_counts()` exposes the native fixed-width read-start
counting kernel. It returns one row per bin across the selected contig
span, including zero-count bins, which makes it suitable as a dense CNV
binning primitive. It can also add one-pass GC and MAPQ summaries on the
same scan.

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

`rduckhts_mosdepth()` writes mosdepth-style outputs to disk and returns
the paths it created. This example writes windowed fragment coverage
from the bundled BAM fixture and previews the generated regions BED.gz.

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

transform(mos_out[, c("summary_path", "regions_path")],
          summary_path = "<tempfile>", regions_path = "<tempfile>")
#>   summary_path regions_path
#> 1   <tempfile>   <tempfile>

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

lift_src_index <- rduckhts_fasta_index(
  con, lift_src, index_path = paste0(lift_src, ".fai")
)
lift_src_index$index_path <- "<tempfile>"
lift_src_index
#>   success index_path
#> 1    TRUE <tempfile>

lift_dst_index <- rduckhts_fasta_index(
  con, lift_dst, index_path = paste0(lift_dst, ".fai")
)
lift_dst_index$index_path <- "<tempfile>"
lift_dst_index
#>   success index_path
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
transform(rduckhts_fasta_index(con, munge_fasta, index_path = paste0(munge_fasta, ".fai")),
          index_path = "<tempfile>")
#>   success index_path
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

`rduckhts_score()` computes per-sample polygenic risk scores (PRS) from
a genotype VCF/BCF and one or more GWAS summary statistics files,
wrapping `bcftools_score`.

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
transform(bgzip_meta[, c("success", "output_path", "bytes_out")],
          output_path = "<tempfile>")
#>   success output_path bytes_out
#> 1    TRUE  <tempfile>       169

bgunzip_meta <- rduckhts_bgunzip(
  con, tmp_bgz,
  output_path = tmp_roundtrip,
  threads = 1,
  keep = TRUE,
  overwrite = TRUE
)
bgunzip_meta$output_path <- "<tempfile>"
bgunzip_meta[, c("success", "output_path", "bytes_out")]
#>   success output_path bytes_out
#> 1    TRUE  <tempfile>       194

bam_index_meta <- rduckhts_bam_index(
  con, bam_src,
  index_path = tmp_bai,
  threads = 1
)
transform(bam_index_meta, index_path = "<tempfile>")
#>   success index_path index_format
#> 1    TRUE <tempfile>          BAI

bcf_index_meta <- rduckhts_bcf_index(
  con, bcf_src,
  index_path = tmp_csi,
  threads = 1
)
transform(bcf_index_meta, index_path = "<tempfile>")
#>   success index_path index_format
#> 1    TRUE <tempfile>          CSI

tabix_meta <- rduckhts_tabix_index(
  con, tmp_bgz,
  preset = "bed",
  index_path = tmp_tbi,
  threads = 1
)
transform(tabix_meta, index_path = "<tempfile>")
#>   success index_path index_format
#> 1    TRUE <tempfile>          TBI

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
fai_path <- file.path(tempdir(), "duckhts_readme_00000000000000.fai")
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

S3 queries require a build with the S3 transport enabled; RTools builds
do not enable it. This remote example is an unevaluated usage snippet.

``` r
# Example S3 URL (1000 Genomes cohort VCF)
s3_base <- "s3://1000genomes-dragen-v3.7.6/data/cohorts/"
s3_path <- "gvcf-genotyper-dragen-3.7.6/hg19/3202-samples-cohort/"
s3_vcf_file <- "3202_samples_cohort_gg_chr22.vcf.gz"
s3_vcf_uri <- paste0(s3_base, s3_path, s3_vcf_file)

rduckhts_bcf(con, "s3_variants", s3_vcf_uri, region = "chr22:16050000-16050500", overwrite = TRUE)
dbGetQuery(con, "SELECT CHROM, COUNT(*) AS n FROM s3_variants GROUP BY CHROM")
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

### FASTQ quality decoding and fused QC

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

For ordinary QC, aggregate the sequence and canonical quality strings
directly. This returns exact global totals plus a small nested per-cycle
relation without generating one SQL row per base.

``` r
fastq_qc <- dbGetQuery(
  con,
  sprintf(
    "WITH q AS (
       SELECT duckhts_fastq_qc(SEQUENCE, QUALITY) AS qc
       FROM read_fastq('%s')
     )
     SELECT qc.reads, qc.bases, qc.q30_bases, qc.max_read_length
     FROM q",
    fastq_r1
  )
)
fastq_qc
#>   reads bases q30_bases max_read_length
#> 1     5   500       475             100
```

Use numeric quality arrays when the query genuinely needs the full
per-position histogram:

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

GFF3 files are read with `rduckhts_gff()` / SQL `read_gff(...)`; GTF
files are read with `rduckhts_gtf()` / SQL `read_gtf(...)`.
`strict = TRUE` enables GFF3 structural validation. Attribute decoding
can be scalar and raw for legacy convenience (`attributes_map`), grouped
and lossless for multi-values (`attributes_list`, a DuckDB
`MAP(VARCHAR, VARCHAR[])`), or exact parser-style pairs
(`attributes_pairs`, a DuckDB `LIST<STRUCT(key, value, idx)>`).

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
depends on CORS policy on both the data object and index object. This
remote example is an unevaluated usage snippet.

``` r
gtex_url <- "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge/Brain_Cerebellar_Hemisphere.tsv.gz" 
rduckhts_tabix(con, "gtex_eqtl", gtex_url, region = "1:11868-14409",
                  header = TRUE, auto_detect = TRUE, overwrite = TRUE)
dbGetQuery(con, "SELECT * FROM gtex_eqtl LIMIT 5")
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
