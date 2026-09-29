
<!-- README.md is generated from README.Rmd using [duckknit](https://github.com/rundel/duckknit). Please edit that file. -->

# DuckHTS

[![CRAN
Status](https://www.r-pkg.org/badges/version/Rduckhts)](https://cran.r-project.org/package=Rduckhts)
[![R-universe
version](https://RGenomicsETL.r-universe.dev/Rduckhts/badges/version)](https://RGenomicsETL.r-universe.dev/Rduckhts)

Read VCF, BCF, BAM, CRAM, FASTA, FASTQ, BigWig, GTF, GFF, GenBank, BED,
and tabix-indexed files directly in [DuckDB](https://duckdb.org/),
locally or over HTTP and S3. DuckHTS uses
[htslib](https://github.com/samtools/htslib) for HTS formats and
provides SQL functions for intervals, coverage, sequence operations,
Somalier-style sample identity, ancestry proportions, compression,
indexing and export. Consequence annotation lives in the sibling
[DuckVEP](https://github.com/RGenomicsETL/DuckVEP) extension and its R
package Rduckvep.

The same extension ships three ways:

- **R**: [`Rduckhts`](https://RGenomicsETL.r-universe.dev/Rduckhts)
  bundles the extension and wraps its readers;
  `install.packages("Rduckhts")` from CRAN, or the development build
  from r-universe:
  `install.packages("Rduckhts", repos = c("https://rgenomicsetl.r-universe.dev", "https://cloud.r-project.org"))`.
- **DuckDB** (CLI, Python, any client):
  `INSTALL duckhts FROM community; LOAD duckhts;`.
- **Browser**: the [`duckhts`](https://www.npmjs.com/package/duckhts)
  npm package for duckdb-wasm.

In R, a connection comes with the extension loaded and every reader is
plain SQL:

``` r
library(DBI)
library(Rduckhts)
con <- rduckhts_connect()
vcf <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
dbGetQuery(con, sprintf(
  "SELECT CHROM, POS, REF, ALT, SAMPLE_ID, FORMAT_GT
   FROM read_bcf(%s, tidy_format := true) LIMIT 4", dbQuoteString(con, vcf)))
#>   CHROM     POS REF ALT SAMPLE_ID FORMAT_GT
#> 1     1 3000150   C   T         A       0/1
#> 2     1 3000150   C   T         B       0/1
#> 3     1 3000151   C   T         A       0/1
#> 4     1 3000151   C   T         B       0/1
dbDisconnect(con, shutdown = TRUE)
```

## Runtime compatibility

DuckHTS requires DuckDB **1.4.0 or newer** for its SQL functions and
macros. The stable **v1.2.0 C API** recorded in `C_STRUCT` extension
metadata is a separate interface version; it does not imply support for
DuckDB 1.2 SQL. For duckdb-wasm, check the embedded DuckDB engine
version rather than comparing npm package version numbers with DuckDB
release numbers.

## Using DuckHTS with a database file

`LOAD` creates the 31 SQL macros only when the default database is
writable and in memory. With a file-backed or read-only default database
it registers native functions and leaves the database file untouched;
the macros are then installed as connection-local `TEMP` macros.
`rduckhts_connect()` does that for you, and `rduckhts_install_macros()`
covers any other DBI connection:

``` r
db <- tempfile(fileext = ".duckdb")
con <- rduckhts_connect(dbdir = db)
dbGetQuery(con, "
  SELECT database_name, count(*) AS duckhts_macros
  FROM duckdb_functions()
  WHERE function_name IN (SELECT name FROM duckhts_macro_definitions())
  GROUP BY ALL")
#>   database_name duckhts_macros
#> 1          temp             33
dbGetQuery(con, "SELECT duckhts_quote_ident('my column') AS quoted")
#>        quoted
#> 1 "my column"
dbDisconnect(con, shutdown = TRUE)
```

With the DuckDB CLI, write the exported statements once and `.read` them
on the connection that needs the macros:

``` bash
cd "$(mktemp -d)"
"$DUCKDB" -unsigned example.duckdb "
  SET allow_extensions_metadata_mismatch = true;
  LOAD '$DUCKHTS_EXTENSION';
  COPY (SELECT sql || ';' FROM duckhts_macro_definitions() ORDER BY install_order)
    TO 'duckhts_macros.sql' (HEADER false, QUOTE '', ESCAPE '');"
printf "SET allow_extensions_metadata_mismatch = true;\nLOAD '%s';\n.read duckhts_macros.sql\nSELECT duckhts_quote_ident('my column') AS quoted;\n" \
  "$DUCKHTS_EXTENSION" | "$DUCKDB" -unsigned example.duckdb
#> ┌─────────────┐
#> │   quoted    │
#> │   varchar   │
#> ├─────────────┤
#> │ "my column" │
#> └─────────────┘
```

`TEMP` macros shadow persistent macros written by older DuckHTS
versions; no persistent macro is deleted. Inspect them with
`duckdb_functions()` where `database_name` equals the file’s catalog
name, and drop them explicitly if wanted. Attaching a database after
`LOAD` does not install macros into it.

## Functions

<details>
<summary>
Show generated function catalog
</summary>

## Extension Function Catalog

This section is generated from `functions.yaml`.

### Utilities

| Function                                                                                               | Kind  | R helper | Description                                                                     |
|--------------------------------------------------------------------------------------------------------|-------|----------|---------------------------------------------------------------------------------|
| [`duckhts_macro_definitions`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_macro_definitions) | table |          | Export the ordered DuckHTS macro definitions for connection-local installation. |

### Diagnostics

| Function                                                                                                                 | Kind         | R helper                              | Description                                                                                                                                                                                                                                                                                                                                                                      |
|--------------------------------------------------------------------------------------------------------------------------|--------------|---------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`duckhts_htslib_version`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_htslib_version)                         | scalar       | `rduckhts_htslib_version`             | Return the runtime version reported by the htslib library loaded with DuckHTS. Rduckhts uses this value to reject a downstream linking receipt whose source/header version does not match the loaded library.                                                                                                                                                                    |
| [`duckhts_htslib_features`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_htslib_features)                       | scalar       |                                       | Return the htslib runtime feature bitfield reported by hts_features(). Use duckhts_htslib_feature_string() for the corresponding build description.                                                                                                                                                                                                                              |
| [`duckhts_htslib_feature_string`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_htslib_feature_string)           | scalar       |                                       | Return htslib’s runtime build-feature description, including configured transports, compression libraries, compiler, and build flags. DuckHTS snapshots it once while loading the extension so parallel SQL calls read immutable text.                                                                                                                                           |
| [`duckhts_simd_backend`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_simd_backend)                             | scalar       | `rduckhts_simd_backend`               | Return the current DuckHTS SIMD dispatch label. For explicit scalar or concrete backend requests this is the requested policy; for auto it is the single selected backend when all logical kernels resolve to the same backend, or mixed when per-kernel auto-dispatch resolves to multiple backends. Use duckhts_simd_kernel_info() for per-kernel details.                     |
| [`duckhts_simd_requested_backend`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_simd_requested_backend)         | scalar       | `rduckhts_simd_requested_backend`     | Return the current explicit SIMD backend request, usually auto unless `SELECT backend FROM duckhts_simd_set_backend('auto'\|'scalar'\|backend)` was called. The selected per-kernel backend may differ under auto-dispatch across x86, ARM, wasm, and scalar-only builds.                                                                                                        |
| [`duckhts_simd_backend_compiled`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_simd_backend_compiled)           | scalar       | `rduckhts_simd_backend_compiled`      | Return whether a concrete DuckHTS SIMD backend was compiled into this build. This is independent of whether the current CPU/runtime supports executing that backend; for example avx512 can be compiled but not CPU-supported on the running host.                                                                                                                               |
| [`duckhts_simd_backend_cpu_supported`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_simd_backend_cpu_supported) | scalar       | `rduckhts_simd_backend_cpu_supported` | Return whether the current CPU/runtime supports a concrete DuckHTS SIMD backend, independent of whether DuckHTS compiled an implementation for it. Availability is the intersection of compiled and CPU-supported.                                                                                                                                                               |
| [`duckhts_simd_backend_available`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_simd_backend_available)         | scalar       | `rduckhts_simd_backend_available`     | Return whether a concrete SIMD backend is usable in the current process. Availability means the backend is compiled into DuckHTS and supported by the current CPU/runtime. auto is a selection request rather than a concrete backend and is not reported as available here.                                                                                                     |
| [`duckhts_simd_info`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_simd_info)                                   | table        | `rduckhts_simd_info`                  | Report compiled, runtime-supported and selected status for each concrete DuckHTS SIMD backend.                                                                                                                                                                                                                                                                                   |
| [`duckhts_simd_kernel_info`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_simd_kernel_info)                     | table        | `rduckhts_simd_kernel_info`           | Return one row per logical DuckHTS SIMD kernel showing the concrete backend selected by the current immutable dispatch table, the selected capability, the requested backend policy, whether scalar was used as a per-kernel fallback, and the dispatch mode. This is the authoritative diagnostic for mixed auto-dispatch when different kernels resolve to different backends. |
| [`duckhts_simd_set_backend`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_simd_set_backend)                     | table        | `rduckhts_simd_set_backend`           | Explicitly select the DuckHTS SIMD dispatch policy for this process using a one-row table-function call and return the current dispatch label in a backend column. Use auto for per-kernel runtime dispatch or scalar for a portable baseline; unavailable platform-specific requests such as avx512 on non-AVX-512 CPUs raise an error instead of silently falling back.        |
| [`duckhts_duckdb_type_supported`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_duckdb_type_supported)           | scalar_macro |                                       | Return whether the currently open DuckDB runtime advertises a logical type with the given name through duckdb_types(). This is a catalog-level runtime probe for feature gating SQL/macros across DuckDB versions.                                                                                                                                                               |
| [`duckhts_duckdb_supports_variant`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_duckdb_supports_variant)       | scalar_macro |                                       | Return whether the currently open DuckDB runtime advertises the VARIANT logical type. Use this to gate optional SQL that depends on DuckDB VARIANT support.                                                                                                                                                                                                                      |
| [`duckhts_duckdb_supports_geometry`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_duckdb_supports_geometry)     | scalar_macro |                                       | Return whether the currently open DuckDB runtime advertises the GEOMETRY logical type. Use this to gate optional SQL that depends on DuckDB GEOMETRY support.                                                                                                                                                                                                                    |

### Readers

| Function                                                                                         | Kind         | R helper                                                                                                                                                               | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
|--------------------------------------------------------------------------------------------------|--------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`read_bcf`](r/Rduckhts/inst/function_catalog/reference.md#read_bcf)                             | table        | `rduckhts_bcf`                                                                                                                                                         | Read VCF/BCF with header-typed INFO/FORMAT, typed CSQ/ANN/BCSQ annotations, sample selection and optional tidy sample rows.                                                                                                                                                                                                                                                                                                                                                                                  |
| [`read_geno`](r/Rduckhts/inst/function_catalog/reference.md#read_geno)                           | table        | `rduckhts_geno`                                                                                                                                                        | Read one row per VCF/BCF record with typed arbitrary-ploidy GT/PS calls, selected FORMAT fields and optional original VCF genotype text.                                                                                                                                                                                                                                                                                                                                                                     |
| [`read_bcf_samples`](r/Rduckhts/inst/function_catalog/reference.md#read_bcf_samples)             | table        | `rduckhts_bcf_samples`                                                                                                                                                 | Read the typed VCF/BCF sample catalog as sample_index UINTEGER and sample_name VARCHAR without reading records. Indices are zero-based positions in the original header, remain stable under selection and join read_geno calls from the same unchanged file. NULL or ‘-’ selects all; an empty string selects none; comma-separated names include samples; ‘^’ excludes them. Names are validated by HTSlib, selected rows retain header order, and unknown names error.                                    |
| [`read_bam`](r/Rduckhts/inst/function_catalog/reference.md#read_bam)                             | table        | `rduckhts_bam`                                                                                                                                                         | Read SAM/BAM/CRAM alignments with optional typed SAM tags, auxiliary maps and packed sequence, quality or CIGAR output.                                                                                                                                                                                                                                                                                                                                                                                      |
| [`read_fasta`](r/Rduckhts/inst/function_catalog/reference.md#read_fasta)                         | table        | `rduckhts_fasta`                                                                                                                                                       | Read full FASTA records or indexed regions with text or packed sequence output.                                                                                                                                                                                                                                                                                                                                                                                                                              |
| [`read_bed`](r/Rduckhts/inst/function_catalog/reference.md#read_bed)                             | table        | `rduckhts_bed`                                                                                                                                                         | Read BED3-BED12 interval files with canonical typed columns and optional tabix-backed region filtering.                                                                                                                                                                                                                                                                                                                                                                                                      |
| [`fasta_nuc`](r/Rduckhts/inst/function_catalog/reference.md#fasta_nuc)                           | table        | `rduckhts_fasta_nuc`                                                                                                                                                   | Compute bedtools nuc-style nucleotide composition for supplied BED intervals or generated fixed-width bins over a FASTA reference. A failed reference fetch fails the query with the file and zero-based half-open interval; requested intervals are not silently omitted. For bgzipped FASTA, gzi_path may point to an explicit .gzi sidecar when it is not colocated with the FASTA.                                                                                                                       |
| [`read_fastq`](r/Rduckhts/inst/function_catalog/reference.md#read_fastq)                         | table        | `rduckhts_fastq`                                                                                                                                                       | Read single-end, paired-end, or interleaved FASTQ files with optional legacy quality decoding. By default, FASTQ qualities are interpreted as modern Phred+33 input. Use sequence_encoding := ‘nt16’ to return SEQUENCE as UTINYINT\[\] and quality_representation := ‘phred’ to return QUALITY as UTINYINT\[\] instead of VARCHAR. input_quality_encoding accepts ‘phred33’, ‘auto’, ‘phred64’, or ‘solexa64’. scan_mode := ‘sequential’ forces raw streaming/counting instead of index-backed count paths. |
| [`read_bigwig`](r/Rduckhts/inst/function_catalog/reference.md#read_bigwig)                       | table        | `rduckhts_bigwig`                                                                                                                                                      | Read stored BigWig signal intervals as CHROM, START0, END0 and VALUE.                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| [`read_gff`](r/Rduckhts/inst/function_catalog/reference.md#read_gff)                             | table        | `rduckhts_gff`                                                                                                                                                         | Read GFF annotations with optional parsed attributes, strict GFF3 validation and indexed region selection.                                                                                                                                                                                                                                                                                                                                                                                                   |
| [`read_gtf`](r/Rduckhts/inst/function_catalog/reference.md#read_gtf)                             | table        | `rduckhts_gtf`                                                                                                                                                         | Read GTF annotations with optional parsed attributes and indexed region selection.                                                                                                                                                                                                                                                                                                                                                                                                                           |
| [`read_genbank`](r/Rduckhts/inst/function_catalog/reference.md#read_genbank)                     | table        | `rduckhts_genbank`                                                                                                                                                     | Read GenBank flat-file features in read_gff’s column shape, with optional parsed qualifier MAP.                                                                                                                                                                                                                                                                                                                                                                                                              |
| [`read_tabix`](r/Rduckhts/inst/function_catalog/reference.md#read_tabix)                         | table        | `rduckhts_tabix`                                                                                                                                                       | Read tabix-indexed text with optional header handling, inferred types and region selection.                                                                                                                                                                                                                                                                                                                                                                                                                  |
| [`fasta_index`](r/Rduckhts/inst/function_catalog/reference.md#fasta_index)                       | table        | `rduckhts_fasta_index`                                                                                                                                                 | Build a FASTA index (.fai) and return a single row with columns success (BOOLEAN) and index_path (VARCHAR).                                                                                                                                                                                                                                                                                                                                                                                                  |
| [`hts_union_query`](r/Rduckhts/inst/function_catalog/reference.md#hts_union_query)               | scalar_macro | `rduckhts_bam_multi, rduckhts_bcf_multi, rduckhts_fastq_multi, rduckhts_fasta_multi, rduckhts_bed_multi, rduckhts_tabix_multi, rduckhts_gff_multi, rduckhts_gtf_multi` | Generate a UNION ALL BY NAME query string that reads every file matching a glob pattern through the named reader function. The result includes a ‘filename’ column identifying the source file for each row. Assign to a variable with SET VARIABLE and execute via query(getvariable(…)). Optional params string is appended to each reader call. In R, use the typed rduckhts\_\*\_multi() helpers instead, which accept file vectors with optional per-file parameters and create DuckDB tables directly. |
| [`hts_region_union_query`](r/Rduckhts/inst/function_catalog/reference.md#hts_region_union_query) | scalar_macro |                                                                                                                                                                        | Generate UNION ALL BY NAME SQL over separate per-region scans of one HTS file.                                                                                                                                                                                                                                                                                                                                                                                                                               |

### Converters

| Function                                                                                                               | Kind         | R helper                         | Description                                                                                                                    |
|------------------------------------------------------------------------------------------------------------------------|--------------|----------------------------------|--------------------------------------------------------------------------------------------------------------------------------|
| [`duckhts_bcf_convert_parquet_sql`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_bcf_convert_parquet_sql)     | scalar_macro | `rduckhts_bcf_convert_parquet`   | Build COPY SQL for read_bcf() output with Parquet metadata, VCF header text and selected columns, filters or partitions.       |
| [`duckhts_bam_convert_parquet_sql`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_bam_convert_parquet_sql)     | scalar_macro | `rduckhts_bam_convert_parquet`   | Build COPY SQL for read_bam() output with Parquet metadata, SAM header text and selected columns, filters or partitions.       |
| [`duckhts_gff_convert_parquet_sql`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_gff_convert_parquet_sql)     | scalar_macro | `rduckhts_gff_convert_parquet`   | Build COPY SQL for read_gff() output with Parquet metadata, GFF/tabix header text and selected columns, filters or partitions. |
| [`duckhts_tabix_convert_parquet_sql`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_tabix_convert_parquet_sql) | scalar_macro | `rduckhts_tabix_convert_parquet` | Build COPY SQL for read_tabix() output with Parquet metadata, header text and selected columns, filters or partitions.         |
| [`genbank_to_fasta`](r/Rduckhts/inst/function_catalog/reference.md#genbank_to_fasta)                                   | table        | `rduckhts_genbank_to_fasta`      | Write the ORIGIN sequence of each GenBank record as FASTA and return success, output_path and records_written.                 |

### Coverage

| Function                                                                                             | Kind  | R helper                    | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
|------------------------------------------------------------------------------------------------------|-------|-----------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`read_pileup`](r/Rduckhts/inst/function_catalog/reference.md#read_pileup)                           | table | `rduckhts_pileup`           | Construct a region-scoped BAM pileup with one row per covered position, emitting chrom, 1-based position, depth, observed bases, and Phred+33 qualities after SAM flag and MAPQ filtering. This is a compact htslib pileup view, not samtools mpileup text parity.                                                                                                                                                                                                                                        |
| [`bam_bin_counts`](r/Rduckhts/inst/function_catalog/reference.md#bam_bin_counts)                     | table | `rduckhts_bam_bin_counts`   | Count BAM or CRAM read starts into fixed-width bins. Returns one row per bin across the selected contig span, including zero-count bins, with total, forward, and reverse counts; `rmdup := 'streaming'` applies the WisecondorX-style larp/larp2 consecutive-position deduplication, `rmdup := 'flag'` drops SAM duplicate-flagged reads, and `stats := 'gc'`, `'mq'`, or `'gc,mq'` adds per-bin pre/post-filter GC and MAPQ sufficient statistics, including reference GC when `reference` is provided. |
| [`duckhts_bam_bed_coverage`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_bam_bed_coverage) | table | `rduckhts_bam_bed_coverage` | Compute samtools coverage-like regional summaries for BAM or CRAM input over a BED target set, returning one row per BED interval with DuckHTS-specific pre/post-filter read counts, covered bases, percentage covered, mean depth, mean baseQ, mean mapQ, and strand-specific post-filter summaries in read mode. Indexed BAM/CRAM input is required in the current implementation. decompression_threads controls htslib worker threads for BAM/CRAM decoding; use 0 to disable them.                   |
| [`duckhts_mosdepth`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_mosdepth)                 | table | `rduckhts_mosdepth`         | Write mosdepth-compatible coverage files from indexed BAM/CRAM.                                                                                                                                                                                                                                                                                                                                                                                                                                           |

### Intervals

| Function                                                                                                             | Kind        | R helper | Description                                                                                                                                                                                                                                                                                                                                                                                                                  |
|----------------------------------------------------------------------------------------------------------------------|-------------|----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`duckhts_cgranges_create`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_cgranges_create)                   | scalar      |          | Create an empty session-scoped cgranges registry entry that can be populated with intervals and finalized for overlap queries.                                                                                                                                                                                                                                                                                               |
| [`duckhts_cgranges_add`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_cgranges_add)                         | scalar      |          | Append an interval to a session-scoped cgranges registry entry before finalization. Labels may be BIGINT-like, DOUBLE, VARCHAR, or BOOLEAN.                                                                                                                                                                                                                                                                                  |
| [`duckhts_cgranges_index`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_cgranges_index)                     | scalar      |          | Finalize a populated cgranges registry entry and build its immutable overlap index for subsequent queries.                                                                                                                                                                                                                                                                                                                   |
| [`duckhts_cgranges_destroy`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_cgranges_destroy)                 | scalar      |          | Destroy a session-scoped cgranges registry entry and release its indexed interval storage when it is not in active use.                                                                                                                                                                                                                                                                                                      |
| [`duckhts_cgranges_from_table`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_cgranges_from_table)           | table_macro |          | Create, populate and finalize a session-scoped cgranges registry entry from the rows of a table or view, on the caller’s connection.                                                                                                                                                                                                                                                                                         |
| [`duckhts_cgranges_has_overlap`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_cgranges_has_overlap)         | scalar      |          | Vectorized scalar predicate for streaming provider rows through a finalized session-scoped cgranges index. Returns TRUE when the query interval overlaps at least one indexed interval, or when mode = ‘contain’ and it fully contains at least one indexed interval; NULL inputs return NULL.                                                                                                                               |
| [`duckhts_cgranges_count_overlaps`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_cgranges_count_overlaps)   | scalar      |          | Vectorized scalar overlap counter for streaming provider rows through a finalized session-scoped cgranges index. Returns the number of indexed intervals that overlap the query interval, or with mode = ‘contain’ the number fully contained by it; NULL inputs return NULL.                                                                                                                                                |
| [`duckhts_cgranges_overlaps_list`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_cgranges_overlaps_list)     | scalar      |          | Vectorized scalar overlap expander for streaming provider rows through a finalized session-scoped cgranges index. Returns a LIST of hit STRUCTs that can be expanded with UNNEST, preserving provider columns while emitting one row per matching indexed interval. Because scalar return types are fixed, labels are returned as text with label_type describing the original cgranges label kind; NULL inputs return NULL. |
| [`duckhts_cgranges_overlaps`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_cgranges_overlaps)               | table       |          | Query a finalized session-scoped cgranges registry entry and return one row per overlapping or containing indexed interval, preserving the original label type and interval coordinates.                                                                                                                                                                                                                                     |
| [`regionkey`](r/Rduckhts/inst/function_catalog/reference.md#regionkey)                                               | scalar      |          | Encode a genomic interval as an official RegionKey-compatible 64-bit unsigned integer. Start and end use 0-based half-open interval semantics, matching BED-style coordinates; strand accepts -1, 0, or 1.                                                                                                                                                                                                                   |
| [`regionkey_hex`](r/Rduckhts/inst/function_catalog/reference.md#regionkey_hex)                                       | scalar      |          | Render a RegionKey as its lowercase 16-character hexadecimal string representation.                                                                                                                                                                                                                                                                                                                                          |
| [`parse_regionkey_hex`](r/Rduckhts/inst/function_catalog/reference.md#parse_regionkey_hex)                           | scalar      |          | Parse a 16-character hexadecimal RegionKey string back into its UBIGINT code. Invalid or non-hex strings return NULL.                                                                                                                                                                                                                                                                                                        |
| [`encode_regionkey`](r/Rduckhts/inst/function_catalog/reference.md#encode_regionkey)                                 | scalar      |          | Encode the raw upstream RegionKey fields directly: chromosome code, 0-based start, 0-based end, and strand code (0 = unknown, 1 = +, 2 = -).                                                                                                                                                                                                                                                                                 |
| [`extract_regionkey_chrom`](r/Rduckhts/inst/function_catalog/reference.md#extract_regionkey_chrom)                   | scalar      |          | Extract the raw upstream RegionKey chromosome code.                                                                                                                                                                                                                                                                                                                                                                          |
| [`extract_regionkey_startpos`](r/Rduckhts/inst/function_catalog/reference.md#extract_regionkey_startpos)             | scalar      |          | Extract the raw upstream RegionKey 0-based start position.                                                                                                                                                                                                                                                                                                                                                                   |
| [`extract_regionkey_endpos`](r/Rduckhts/inst/function_catalog/reference.md#extract_regionkey_endpos)                 | scalar      |          | Extract the raw upstream RegionKey 0-based end position.                                                                                                                                                                                                                                                                                                                                                                     |
| [`extract_regionkey_strand`](r/Rduckhts/inst/function_catalog/reference.md#extract_regionkey_strand)                 | scalar      |          | Extract the raw upstream RegionKey strand code (0 = unknown, 1 = +, 2 = -).                                                                                                                                                                                                                                                                                                                                                  |
| [`decode_regionkey`](r/Rduckhts/inst/function_catalog/reference.md#decode_regionkey)                                 | scalar      |          | Decode a RegionKey into its raw upstream numeric fields: chrom_code, start, end, and strand_code.                                                                                                                                                                                                                                                                                                                            |
| [`reverse_regionkey`](r/Rduckhts/inst/function_catalog/reference.md#reverse_regionkey)                               | scalar      |          | Decode a RegionKey into a STRUCT with chrom, chrom_code, start, end, strand, and strand_code.                                                                                                                                                                                                                                                                                                                                |
| [`extend_regionkey`](r/Rduckhts/inst/function_catalog/reference.md#extend_regionkey)                                 | scalar      |          | Extend a RegionKey interval by a fixed number of bases on both sides, clamping to the official 28-bit RegionKey position range.                                                                                                                                                                                                                                                                                              |
| [`are_overlapping_regions`](r/Rduckhts/inst/function_catalog/reference.md#are_overlapping_regions)                   | scalar      |          | Return TRUE when two explicit 0-based half-open intervals overlap on the same canonical chromosome.                                                                                                                                                                                                                                                                                                                          |
| [`are_overlapping_region_regionkey`](r/Rduckhts/inst/function_catalog/reference.md#are_overlapping_region_regionkey) | scalar      |          | Return TRUE when a 0-based half-open interval overlaps the supplied RegionKey interval.                                                                                                                                                                                                                                                                                                                                      |
| [`are_overlapping_regionkeys`](r/Rduckhts/inst/function_catalog/reference.md#are_overlapping_regionkeys)             | scalar      |          | Return TRUE when two RegionKeys overlap.                                                                                                                                                                                                                                                                                                                                                                                     |

### Quality Control

| Function                                                                             | Kind      | R helper | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
|--------------------------------------------------------------------------------------|-----------|----------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`duckhts_fastq_qc`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_fastq_qc) | aggregate |          | Aggregate canonical sequence and Phred+33 quality strings directly into exact read/base/Q20/Q30/Q40, nucleotide, quality-sum, and per-cycle sufficient statistics. The nested cycles list supports mean-quality, nucleotide-content, GC, and read-length curves without expanding one SQL row per base. Rows with any NULL input are ignored. Per-cycle state defaults to at most 1,048,576 cycles; pass a constant max_cycles per aggregate group to choose a larger explicit limit, up to 16,777,216. |

### Sample Identity

| Function                                                                                                                         | Kind         | R helper                                  | Description                                                                                                 |
|----------------------------------------------------------------------------------------------------------------------------------|--------------|-------------------------------------------|-------------------------------------------------------------------------------------------------------------|
| [`duckhts_somalier_spacing`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_somalier_spacing)                             | scalar       | `rduckhts_somalier_find_sites`            | Greedily select ranked Somalier candidate positions at a minimum genomic distance.                          |
| [`duckhts_somalier_import_sites`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_somalier_import_sites)                   | table_macro  | `rduckhts_somalier_import_sites`          | Import an already selected Somalier sites VCF/BCF as one canonical panel and population-frequency relation. |
| [`duckhts_somalier_vcf_counts`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_somalier_vcf_counts)                       | table_macro  | `rduckhts_somalier_vcf_counts`            | Extract a complete panel-aligned A/B/other count relation from VCF/BCF FORMAT/AD.                           |
| [`duckhts_somalier_bam_counts`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_somalier_bam_counts)                       | table        | `rduckhts_somalier_bam_counts`            | Extract complete panel-aligned A/B/other base counts from one indexed BAM or CRAM source.                   |
| [`duckhts_ancestry_proportions`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_ancestry_proportions)                     | scalar       | `rduckhts_ancestry_proportions`           | Solve a nearest-positive-definite constrained ancestry projection from aggregated PC products.              |
| [`duckhts_somalier_panel_sha256`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_somalier_panel_sha256)                   | scalar_macro |                                           | Derive a stable SHA-256 identity for an ordered biallelic sample-fingerprinting panel.                      |
| [`duckhts_somalier_frequency_sha256`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_somalier_frequency_sha256)           | scalar_macro |                                           | Derive a stable identity for panel-aligned population-B allele frequencies.                                 |
| [`duckhts_somalier_classify`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_somalier_classify)                           | scalar       |                                           | Classify one measured A/B/other count tuple for Somalier-derived autosomal relatedness.                     |
| [`duckhts_somalier_prepare_sketches`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_somalier_prepare_sketches)           | table_macro  | `rduckhts_somalier_sketches`              | Build one panel-verified packed relatedness sketch per sample from typed count evidence.                    |
| [`duckhts_somalier_verify_sketches`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_somalier_verify_sketches)             | scalar_macro |                                           | Verify persisted relatedness sketches against their retained raw count evidence.                            |
| [`duckhts_somalier_relatedness`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_somalier_relatedness)                     | scalar       | `rduckhts_somalier_relatedness`           | Compute fused Somalier-derived relatedness and concordance statistics for two prepared sketches.            |
| [`duckhts_somalier_verify_relatedness`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_somalier_verify_relatedness)       | scalar       |                                           | Verify a typed relatedness result against its two sealed sketches.                                          |
| [`duckhts_somalier_charr`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_somalier_charr)                                 | table_macro  | `rduckhts_somalier_charr`                 | Estimate per-sample contamination with a bounded Somalier-derived CHARR reduction.                          |
| [`duckhts_somalier_matched_contamination`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_somalier_matched_contamination) | table_macro  | `rduckhts_somalier_matched_contamination` | Estimate directional contamination for explicitly selected receiver/anchor sample pairs.                    |

### Metadata

| Function                                                                                               | Kind        | R helper                           | Description                                                                                                                                                                                                                                                                |
|--------------------------------------------------------------------------------------------------------|-------------|------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`detect_quality_encoding`](r/Rduckhts/inst/function_catalog/reference.md#detect_quality_encoding)     | table       | `rduckhts_detect_quality_encoding` | Inspect a FASTQ file’s observed quality ASCII range and report compatible legacy encodings with a heuristic guessed encoding.                                                                                                                                              |
| [`duckhts_samtools_idxstats`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_samtools_idxstats) | table       | `rduckhts_samtools_idxstats`       | Write samtools idxstats-compatible TAB-delimited output for BAM, CRAM, or SAM input. Indexed BAM uses `hts_idx_get_stat(...)` for the fast path; CRAM, SAM, and unindexed BAM fall back to a full scan while preserving samtools-style contig rows plus the final `*` row. |
| [`read_hts_header`](r/Rduckhts/inst/function_catalog/reference.md#read_hts_header)                     | table       | `rduckhts_hts_header`              | Inspect HTS headers in parsed, raw, or combined form across supported formats. Raw VCF/BCF mode includes the final `#CHROM` sample header line so the returned text is suitable for Parquet metadata and future VCF/BCF regeneration.                                      |
| [`read_hts_index`](r/Rduckhts/inst/function_catalog/reference.md#read_hts_index)                       | table       | `rduckhts_hts_index`               | Inspect high-level HTS index metadata such as sequence names and mapped counts.                                                                                                                                                                                            |
| [`read_hts_index_spans`](r/Rduckhts/inst/function_catalog/reference.md#read_hts_index_spans)           | table       | `rduckhts_hts_index_spans`         | Expand index metadata into span and chunk rows suitable for low-level index inspection.                                                                                                                                                                                    |
| [`read_hts_index_raw`](r/Rduckhts/inst/function_catalog/reference.md#read_hts_index_raw)               | table_macro | `rduckhts_hts_index_raw`           | Return the raw on-disk HTS index blob together with basic identifying metadata.                                                                                                                                                                                            |

### Compression

| Function                                                           | Kind  | R helper           | Description                                                                           |
|--------------------------------------------------------------------|-------|--------------------|---------------------------------------------------------------------------------------|
| [`bgzip`](r/Rduckhts/inst/function_catalog/reference.md#bgzip)     | table | `rduckhts_bgzip`   | Compress a plain file to BGZF and return the created output path and byte counts.     |
| [`bgunzip`](r/Rduckhts/inst/function_catalog/reference.md#bgunzip) | table | `rduckhts_bgunzip` | Decompress a BGZF-compressed file and return the created output path and byte counts. |

### Indexing

| Function                                                                   | Kind  | R helper               | Description                                                                                        |
|----------------------------------------------------------------------------|-------|------------------------|----------------------------------------------------------------------------------------------------|
| [`bam_index`](r/Rduckhts/inst/function_catalog/reference.md#bam_index)     | table | `rduckhts_bam_index`   | Build a BAM or CRAM index and report the written index path and format.                            |
| [`bcf_index`](r/Rduckhts/inst/function_catalog/reference.md#bcf_index)     | table | `rduckhts_bcf_index`   | Build a TBI or CSI index for a VCF or BCF file and report the written index path and format.       |
| [`tabix_index`](r/Rduckhts/inst/function_catalog/reference.md#tabix_index) | table | `rduckhts_tabix_index` | Build a tabix index for a BGZF-compressed text file using a preset or explicit coordinate columns. |

### Variants

| Function                                                                                               | Kind        | R helper                 | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
|--------------------------------------------------------------------------------------------------------|-------------|--------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`variantkey`](r/Rduckhts/inst/function_catalog/reference.md#variantkey)                               | scalar      |                          | Encode a normalized biallelic variant as an official VariantKey-compatible 64-bit unsigned integer. This DuckHTS wrapper accepts 1-based VCF/DuckHTS POS to match bcftools `%VKX` / `+add-variantkey`, internally converts to the upstream 0-based field, and preserves the official hashed nonreversible mode for large, ambiguous, and symbolic REF/ALT strings. Only CHROM, POS, REF, and ALT are encoded; END, SVLEN, mate breakend coordinates, and other SV metadata are not.        |
| [`variantkey_hex`](r/Rduckhts/inst/function_catalog/reference.md#variantkey_hex)                       | scalar      |                          | Render a VariantKey as its lowercase 16-character hexadecimal string representation.                                                                                                                                                                                                                                                                                                                                                                                                       |
| [`parse_variantkey_hex`](r/Rduckhts/inst/function_catalog/reference.md#parse_variantkey_hex)           | scalar      |                          | Parse a 16-character hexadecimal VariantKey string back into its UBIGINT code. Invalid or non-hex strings return NULL.                                                                                                                                                                                                                                                                                                                                                                     |
| [`encode_variantkey`](r/Rduckhts/inst/function_catalog/reference.md#encode_variantkey)                 | scalar      |                          | Encode the raw upstream VariantKey fields directly: chromosome code, 0-based position, and 31-bit REF+ALT code.                                                                                                                                                                                                                                                                                                                                                                            |
| [`extract_variantkey_chrom`](r/Rduckhts/inst/function_catalog/reference.md#extract_variantkey_chrom)   | scalar      |                          | Extract the raw upstream VariantKey chromosome code.                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| [`extract_variantkey_pos`](r/Rduckhts/inst/function_catalog/reference.md#extract_variantkey_pos)       | scalar      |                          | Extract the raw upstream VariantKey 0-based position field.                                                                                                                                                                                                                                                                                                                                                                                                                                |
| [`extract_variantkey_refalt`](r/Rduckhts/inst/function_catalog/reference.md#extract_variantkey_refalt) | scalar      |                          | Extract the raw upstream 31-bit VariantKey REF+ALT code.                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| [`decode_variantkey`](r/Rduckhts/inst/function_catalog/reference.md#decode_variantkey)                 | scalar      |                          | Decode a VariantKey into its raw upstream numeric fields: chrom_code, pos0, and refalt_code.                                                                                                                                                                                                                                                                                                                                                                                               |
| [`reverse_variantkey`](r/Rduckhts/inst/function_catalog/reference.md#reverse_variantkey)               | scalar      |                          | Decode a VariantKey into a STRUCT with chrom, chrom_code, 1-based pos, upstream 0-based pos0, ref, alt, refalt_code, and reversible. For hashed nonreversible keys, reversible is FALSE and ref/alt are returned as NULL because DuckHTS v1 does not ship the optional NRVK lookup sidecar.                                                                                                                                                                                                |
| [`variantkey_range`](r/Rduckhts/inst/function_catalog/reference.md#variantkey_range)                   | scalar      |                          | Return the inclusive minimum and maximum VariantKey bounds for a chromosome plus 1-based VCF position range, suitable for numeric range filtering on precomputed VariantKeys.                                                                                                                                                                                                                                                                                                              |
| [`duckhts_contig_key`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_contig_key)               | scalar      |                          | Return a conservative contig join key by removing one non-empty leading chr prefix case-insensitively and normalizing M/MT to MT. X and Y are uppercased; all other suffixes are preserved. This does not map numeric sex chromosomes, accessions, patches, or alternate loci.                                                                                                                                                                                                             |
| [`bcftools_liftover`](r/Rduckhts/inst/function_catalog/reference.md#bcftools_liftover)                 | scalar      | `rduckhts_liftover`      | Row-oriented liftover kernel intended to mirror bcftools +liftover semantics as closely as possible while returning one STRUCT per input row with fields: src_chrom, src_pos, src_ref, src_alt, dest_chrom, dest_pos, dest_end, dest_ref, dest_alt, mapped, reverse_complemented, swap, reject_reason, and note. Set no_left_align := true to skip post-liftover left-alignment of lifted indels (mirrors –no-left-align in bcftools +liftover).                                           |
| [`duckdb_liftover`](r/Rduckhts/inst/function_catalog/reference.md#duckdb_liftover)                     | table_macro | `rduckhts_liftover`      | DuckDB-specific wrapper over bcftools_liftover that takes either a table name or a derived-table expression plus column-name strings for chrom/pos/ref/alt and returns the lifted table. The no_left_align parameter mirrors –no-left-align in bcftools +liftover.                                                                                                                                                                                                                         |
| [`bcftools_norm_row`](r/Rduckhts/inst/function_catalog/reference.md#bcftools_norm_row)                 | scalar      |                          | Normalize one variant against FASTA with bcftools/vt-style left alignment.                                                                                                                                                                                                                                                                                                                                                                                                                 |
| [`duckhts_bcftools_norm`](r/Rduckhts/inst/function_catalog/reference.md#duckhts_bcftools_norm)         | table_macro | `rduckhts_bcftools_norm` | Normalize variants from a table or derived-table expression while preserving input columns.                                                                                                                                                                                                                                                                                                                                                                                                |
| [`bcftools_score`](r/Rduckhts/inst/function_catalog/reference.md#bcftools_score)                       | table       | `rduckhts_score`         | Compute polygenic scores from genotype VCF/BCF and summary statistics using bcftools +score dosage semantics.                                                                                                                                                                                                                                                                                                                                                                              |
| [`bcftools_munge_row`](r/Rduckhts/inst/function_catalog/reference.md#bcftools_munge_row)               | scalar      |                          | Normalize one summary-statistics row into GWAS-VCF-style fields (chrom/pos/ref/alt/effect metrics), resolving REF/ALT orientation against a FASTA reference and applying swap-aware sign/frequency/count transforms. The output flag `alleles_swapped` means REF/ALT orientation was swapped to match the FASTA reference.                                                                                                                                                                 |
| [`duckdb_munge`](r/Rduckhts/inst/function_catalog/reference.md#duckdb_munge)                           | table_macro | `rduckhts_munge`         | DuckDB macro wrapper over bcftools_munge_row that maps source columns (via preset or explicit map) and returns normalized GWAS-VCF-style rows with lean outputs and explicit `alleles_swapped` semantics. Output columns: chrom, pos, id, ref, alt, alleles_swapped, filter, ns, ez, nc, es, se, lp, af, ac, ne (16 columns). For METAL meta-analysis output with SI/I2/CQ/ED columns, use duckdb_munge_metal.                                                                             |
| [`duckdb_munge_metal`](r/Rduckhts/inst/function_catalog/reference.md#duckdb_munge_metal)               | table_macro | `rduckhts_munge`         | Extended munge macro with METAL meta-analysis output columns. Same as duckdb_munge but additionally emits: si (imputation info, from INFO input), i2 (Cochran’s I² heterogeneity, from HET_I2), cq (Cochran’s Q -log10 p, from HET_LP or -log10(HET_P)), and ed (effect direction string, from DIRE; +/- flipped on allele swap). The R wrapper rduckhts_munge() auto-dispatches to this macro when metal keys (INFO, HET_I2, HET_P, HET_LP, DIRE) are present in the resolved column map. |

### Sequence UDFs

| Function                                                                           | Kind   | R helper | Description                                                                                                                                                                                                                                                                                                                                                                                                             |
|------------------------------------------------------------------------------------|--------|----------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`seq_revcomp`](r/Rduckhts/inst/function_catalog/reference.md#seq_revcomp)         | scalar |          | Compute the reverse complement of a DNA sequence using A, C, G, T, and N bases. Overloaded: accepts either a VARCHAR text sequence (returns VARCHAR) or a UTINYINT\[\] of htslib nt16 codes as produced by read_bam(sequence_encoding := ‘nt16’) (returns UTINYINT\[\]); the nt16 overload is bit-identical to the text path after decoding, so BAM pipelines can reverse-complement without leaving the nt16 encoding. |
| [`seq_canonical`](r/Rduckhts/inst/function_catalog/reference.md#seq_canonical)     | scalar |          | Return the lexicographically smaller of a sequence and its reverse complement. Overloaded: accepts either a VARCHAR text sequence (returns VARCHAR) or a UTINYINT\[\] of htslib nt16 codes as produced by read_bam(sequence_encoding := ‘nt16’) (returns UTINYINT\[\]); the nt16 overload compares by decoded base order and is bit-identical to the text path after decoding.                                          |
| [`seq_hash_2bit`](r/Rduckhts/inst/function_catalog/reference.md#seq_hash_2bit)     | scalar |          | Encode a short DNA sequence as a 2-bit unsigned integer hash. Overloaded to also accept a UTINYINT\[\] of htslib nt16 codes (from read_bam(sequence_encoding := ‘nt16’)); non-ACGT codes yield NULL, bit-identical to the text path.                                                                                                                                                                                    |
| [`seq_encode_4bit`](r/Rduckhts/inst/function_catalog/reference.md#seq_encode_4bit) | scalar |          | Encode an IUPAC DNA sequence as a list of 4-bit base codes, preserving ambiguity symbols including N.                                                                                                                                                                                                                                                                                                                   |
| [`seq_decode_4bit`](r/Rduckhts/inst/function_catalog/reference.md#seq_decode_4bit) | scalar |          | Decode a list of 4-bit IUPAC DNA base codes back into a sequence string.                                                                                                                                                                                                                                                                                                                                                |
| [`seq_gc_content`](r/Rduckhts/inst/function_catalog/reference.md#seq_gc_content)   | scalar |          | Compute GC fraction for a DNA sequence as a value between 0 and 1. Overloaded: accepts either a VARCHAR text sequence or a UTINYINT\[\] of htslib nt16 codes as produced by read_bam(sequence_encoding := ‘nt16’); the nt16 overload classifies codes directly and is bit-identical to the text path, so BAM pipelines can compute GC without decoding sequences back to text.                                          |
| [`seq_kmers`](r/Rduckhts/inst/function_catalog/reference.md#seq_kmers)             | table  |          | Expand a sequence into positional k-mers with optional canonicalization.                                                                                                                                                                                                                                                                                                                                                |

### SAM Flag UDFs

| Function                                                                                                                     | Kind   | R helper | Description                                                                                                                                                                        |
|------------------------------------------------------------------------------------------------------------------------------|--------|----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`sam_flag_bits`](r/Rduckhts/inst/function_catalog/reference.md#sam_flag_bits)                                               | scalar |          | Decode a SAM flag into a struct of boolean bit fields using explicit SAM-oriented names such as `is_paired`, `is_proper_pair`, `is_next_segment_unmapped`, and `is_supplementary`. |
| [`sam_flag_has`](r/Rduckhts/inst/function_catalog/reference.md#sam_flag_has)                                                 | scalar |          | Test whether any bits from the provided SAM flag mask are set in a flag value.                                                                                                     |
| [`is_forward_aligned`](r/Rduckhts/inst/function_catalog/reference.md#is_forward_aligned)                                     | scalar |          | Test whether a mapped segment is aligned to the forward strand. Returns `NULL` for unmapped segments because SAM flag `0x10` does not define genomic strand when `0x4` is set.     |
| [`is_paired`](r/Rduckhts/inst/function_catalog/reference.md#is_paired)                                                       | scalar |          | Test whether the SAM flag indicates that the template has multiple segments in sequencing (`0x1`).                                                                                 |
| [`is_proper_pair`](r/Rduckhts/inst/function_catalog/reference.md#is_proper_pair)                                             | scalar |          | Test whether the SAM flag indicates that each segment is properly aligned according to the aligner (`0x2`).                                                                        |
| [`is_unmapped`](r/Rduckhts/inst/function_catalog/reference.md#is_unmapped)                                                   | scalar |          | Test whether the read itself is unmapped according to the SAM flag.                                                                                                                |
| [`is_next_segment_unmapped`](r/Rduckhts/inst/function_catalog/reference.md#is_next_segment_unmapped)                         | scalar |          | Test whether the next segment in the template is flagged as unmapped (`0x8`).                                                                                                      |
| [`is_reverse_complemented`](r/Rduckhts/inst/function_catalog/reference.md#is_reverse_complemented)                           | scalar |          | Test whether `SEQ` is stored reverse complemented (`0x10`); for mapped reads this corresponds to reverse-strand alignment.                                                         |
| [`is_next_segment_reverse_complemented`](r/Rduckhts/inst/function_catalog/reference.md#is_next_segment_reverse_complemented) | scalar |          | Test whether `SEQ` of the next segment in the template is stored reverse complemented (`0x20`).                                                                                    |
| [`is_first_segment`](r/Rduckhts/inst/function_catalog/reference.md#is_first_segment)                                         | scalar |          | Test whether the read is marked as the first segment in the template.                                                                                                              |
| [`is_last_segment`](r/Rduckhts/inst/function_catalog/reference.md#is_last_segment)                                           | scalar |          | Test whether the read is marked as the last segment in the template.                                                                                                               |
| [`is_secondary`](r/Rduckhts/inst/function_catalog/reference.md#is_secondary)                                                 | scalar |          | Test whether the alignment is marked as secondary.                                                                                                                                 |
| [`is_qc_fail`](r/Rduckhts/inst/function_catalog/reference.md#is_qc_fail)                                                     | scalar |          | Test whether the read failed vendor or pipeline quality checks.                                                                                                                    |
| [`is_duplicate`](r/Rduckhts/inst/function_catalog/reference.md#is_duplicate)                                                 | scalar |          | Test whether the alignment is flagged as a duplicate.                                                                                                                              |
| [`is_supplementary`](r/Rduckhts/inst/function_catalog/reference.md#is_supplementary)                                         | scalar |          | Test whether the alignment is marked as supplementary.                                                                                                                             |

### CIGAR Utils

| Function                                                                                                 | Kind   | R helper | Description                                                                                                                                                                                                                                                                                                                             |
|----------------------------------------------------------------------------------------------------------|--------|----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [`cigar_has_soft_clip`](r/Rduckhts/inst/function_catalog/reference.md#cigar_has_soft_clip)               | scalar |          | Test whether a CIGAR string contains any soft-clipped segment (`S`). Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                                                                                          |
| [`cigar_has_hard_clip`](r/Rduckhts/inst/function_catalog/reference.md#cigar_has_hard_clip)               | scalar |          | Test whether a CIGAR string contains any hard-clipped segment (`H`). Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                                                                                          |
| [`cigar_left_soft_clip`](r/Rduckhts/inst/function_catalog/reference.md#cigar_left_soft_clip)             | scalar |          | Return the left-end soft-clipped length from a CIGAR string, or zero if the alignment does not start with `S`. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                                                |
| [`cigar_right_soft_clip`](r/Rduckhts/inst/function_catalog/reference.md#cigar_right_soft_clip)           | scalar |          | Return the right-end soft-clipped length from a CIGAR string, or zero if the alignment does not end with `S`. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                                                 |
| [`cigar_query_length`](r/Rduckhts/inst/function_catalog/reference.md#cigar_query_length)                 | scalar |          | Return the query-consuming length from a CIGAR string, counting `M`, `I`, `S`, `=`, and `X`. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                                                                  |
| [`cigar_aligned_query_length`](r/Rduckhts/inst/function_catalog/reference.md#cigar_aligned_query_length) | scalar |          | Return the aligned query length from a CIGAR string, counting `M`, `=`, and `X` but excluding clips and insertions. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                                           |
| [`cigar_reference_length`](r/Rduckhts/inst/function_catalog/reference.md#cigar_reference_length)         | scalar |          | Return the reference-consuming length from a CIGAR string, counting `M`, `D`, `N`, `=`, and `X`. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                                                              |
| [`cigar_has_op`](r/Rduckhts/inst/function_catalog/reference.md#cigar_has_op)                             | scalar |          | Test whether a CIGAR string contains at least one instance of the requested operator. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path.                                                                         |
| [`cigar_aligned_blocks`](r/Rduckhts/inst/function_catalog/reference.md#cigar_aligned_blocks)             | scalar |          | Return the aligned blocks of a CIGAR as a struct of three parallel BIGINT lists: ref_start, query_start and width, one entry per M, = or X op in CIGAR order. Overloaded to also accept a UINTEGER\[\] binary CIGAR (as produced by read_bam(cigar_representation := ‘binary’)); the binary overload is bit-identical to the text path. |

</details>

`read_fastq` with `mate_path` requires exact QNAME pairing. `read_bam`
supports typed `standard_tags` and `auxiliary_tags` maps. `read_tabix`
supports header-aware parsing (`header`, `header_names`) and optional
type inference (`auto_detect`, `column_types`). Region lists in
comma-separated form are supported by `read_bam`, `read_bcf`,
`read_fasta`, `read_bigwig`, `read_gff`, `read_gtf`, and `read_tabix`.
Indexed `read_bam`, `read_bcf`, `read_bigwig`, `read_gff`, `read_gtf`,
and `read_tabix` multi-region queries emit a matching record once when
requested regions overlap. `read_fasta` retains its separate per-region
sequence-row contract.

## Examples

Every example below runs when this README is rendered, through the
DuckDB CLI with the extension bundled in the installed `Rduckhts`
package. Most read the repository’s test files; the BigWig and
remote-URL examples read public HTTP and S3 data.

### Core readers

``` sql
SELECT CHROM, POS, REF, ALT, SAMPLE_ID
FROM read_bcf('test/data/formatcols.vcf.gz', tidy_format := true)
LIMIT 3;
```

    ┌─────────┬───────┬─────────┬───────────┬───────────┐
    │  CHROM  │  POS  │   REF   │    ALT    │ SAMPLE_ID │
    │ varchar │ int64 │ varchar │ varchar[] │  varchar  │
    ├─────────┼───────┼─────────┼───────────┼───────────┤
    │ 1       │   100 │ A       │ [T]       │ S1        │
    │ 1       │   100 │ A       │ [T]       │ S²        │
    │ 1       │   100 │ A       │ [T]       │ S3        │
    └─────────┴───────┴─────────┴───────────┴───────────┘

``` sql
SELECT count(*) AS n
FROM read_bam('test/data/range.bam', region := 'CHROMOSOME_I:1-1000');
```

    ┌───────┐
    │   n   │
    │ int64 │
    ├───────┤
    │     2 │
    └───────┘

``` sql
SELECT *
FROM fasta_index('test/data/ce.fa');
```

    ┌─────────┬─────────────────────┐
    │ success │     index_path      │
    │ boolean │       varchar       │
    ├─────────┼─────────────────────┤
    │ true    │ test/data/ce.fa.fai │
    └─────────┴─────────────────────┘

``` sql
SELECT NAME, length(SEQUENCE) AS seq_length
FROM read_fasta('test/data/ce.fa', region := 'CHROMOSOME_I:1-25');
```

    ┌──────────────┬────────────┐
    │     NAME     │ seq_length │
    │   varchar    │   int64    │
    ├──────────────┼────────────┤
    │ CHROMOSOME_I │         25 │
    └──────────────┴────────────┘

``` sql
SELECT NAME, MATE, PAIR_ID
FROM read_fastq('test/data/interleaved.fq', interleaved := true)
LIMIT 3;
```

    ┌─────────────────────────────────┬────────┬─────────────────────────────────┐
    │              NAME               │  MATE  │             PAIR_ID             │
    │             varchar             │ uint16 │             varchar             │
    ├─────────────────────────────────┼────────┼─────────────────────────────────┤
    │ HS25_09827:2:1201:1505:59795#49 │      1 │ HS25_09827:2:1201:1505:59795#49 │
    │ HS25_09827:2:1201:1505:59795#49 │      2 │ HS25_09827:2:1201:1505:59795#49 │
    │ HS25_09827:2:1201:1559:70726#49 │      1 │ HS25_09827:2:1201:1559:70726#49 │
    └─────────────────────────────────┴────────┴─────────────────────────────────┘

``` sql
SELECT CHROM, START0, END0, round(VALUE::DOUBLE, 1) AS VALUE
FROM read_bigwig(
  'third_party/libBigWig/test/test.bw',
  region := '1:1-150,10:201-300'
)
ORDER BY CHROM, START0;
```

    ┌─────────┬────────┬────────┬────────┐
    │  CHROM  │ START0 │  END0  │ VALUE  │
    │ varchar │ uint32 │ uint32 │ double │
    ├─────────┼────────┼────────┼────────┤
    │ 1       │      0 │      1 │    0.1 │
    │ 1       │      1 │      2 │    0.2 │
    │ 1       │      2 │      3 │    0.3 │
    │ 1       │    100 │    150 │    1.4 │
    │ 10      │    200 │    300 │    2.0 │
    └─────────┴────────┴────────┴────────┘

### BigWig signal tracks

`read_bigwig()` returns the intervals physically stored in a BigWig as
zero-based, half-open `(CHROM, START0, END0, VALUE)` rows. Its optional
`region` uses the same one-based inclusive, comma-separated syntax as
the indexed HTS readers; overlapping requests are merged and do not
duplicate a stored interval. Local files, native HTTP/S3 paths, and
browser HTTP use the same htslib `hFILE` transport already used by
DuckHTS. A full scan distributes nonempty contigs across DuckDB workers;
a multi-region scan distributes merged ranges. `blocks_per_iteration`
controls indexed block batching inside a worker, not the number of
workers.

This query reads a real 100 kb slice of the UCSC GRCh38 phyloP 100-way
track rather than converting it to an intermediate text file:

``` sql
SELECT count(*) AS stored_intervals,
       round(min(VALUE)::DOUBLE, 3) AS minimum,
       round(max(VALUE)::DOUBLE, 3) AS maximum
FROM read_bigwig(
  'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw',
  region := 'chr22:20000000-20099999'
);
```

    ┌──────────────────┬─────────┬─────────┐
    │ stored_intervals │ minimum │ maximum │
    │      int64       │ double  │ double  │
    ├──────────────────┼─────────┼─────────┤
    │            96783 │ -10.787 │   9.602 │
    └──────────────────┴─────────┴─────────┘

### Variant normalization

`duckhts_bcftools_norm(...)` applies bcftools-style FASTA-backed allele
normalization to a regular table or derived relation while preserving
the original columns. In split mode, multiallelic rows are expanded
first and then normalized one ALT at a time.

``` sql
CREATE OR REPLACE TEMP TABLE readme_norm AS
SELECT *
FROM (VALUES
  ('chrS', 2, 'T', 'TT,TTT'),
  ('chrS', 2, 'T', '*,TT')
) AS t(chrom, pos, ref, alt);
```

``` sql
SELECT chrom, pos, ref, alt, alt_index,
       pos_normed, ref_normed, alt_normed, norm_status
FROM duckhts_bcftools_norm(
  'readme_norm',
  'test/data/liftover_repeat_src.fa',
  split_multiallelic := true
)
ORDER BY alt, alt_index;
```

    ┌─────────┬───────┬─────────┬─────────┬───────────┬────────────┬────────────┬────────────┬──────────────────┐
    │  chrom  │  pos  │   ref   │   alt   │ alt_index │ pos_normed │ ref_normed │ alt_normed │   norm_status    │
    │ varchar │ int32 │ varchar │ varchar │   int64   │   int64    │  varchar   │  varchar   │     varchar      │
    ├─────────┼───────┼─────────┼─────────┼───────────┼────────────┼────────────┼────────────┼──────────────────┤
    │ chrS    │     2 │ T       │ *,TT    │         1 │          2 │ T          │ *          │ SpanningDeletion │
    │ chrS    │     2 │ T       │ *,TT    │         2 │          1 │ G          │ GT         │ Normalized       │
    │ chrS    │     2 │ T       │ TT,TTT  │         1 │          1 │ G          │ GT         │ Normalized       │
    │ chrS    │     2 │ T       │ TT,TTT  │         2 │          1 │ G          │ GTT        │ Normalized       │
    └─────────┴───────┴─────────┴─────────┴───────────┴────────────┴────────────┴────────────┴──────────────────┘

### VariantKey + RegionKey

DuckHTS vendors the official VariantKey / RegionKey C API and exposes
SQL helpers that mirror bcftools `%VKX`-style VariantKey output on VCF
rows. `variantkey(...)` accepts 1-based VCF `POS`, while
`regionkey(...)` uses 0-based half-open interval semantics. Large,
ambiguous, and symbolic alleles still encode through the official hashed
nonreversible VariantKey mode, but those keys do not encode `END`,
`SVLEN`, mate breakend coordinates, or other SV metadata; use RegionKey
explicitly for span-oriented interval work. See Nicola Asuni (2018)
<https://doi.org/10.1101/473744>.

``` sql
SELECT variantkey_hex(variantkey('1', 324684, 'C', 'G')) AS vkx,
       reverse_variantkey(parse_variantkey_hex('08027a2588b00000')) AS reversed;
```

    ┌──────────────────┬─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
    │       vkx        │                                                                  reversed                                                                   │
    │     varchar      │ struct(chrom varchar, chrom_code utinyint, pos bigint, pos0 uinteger, "ref" varchar, alt varchar, refalt_code uinteger, reversible boolean) │
    ├──────────────────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
    │ 08027a2588b00000 │ {'chrom': 1, 'chrom_code': 1, 'pos': 324684, 'pos0': 324683, 'ref': C, 'alt': G, 'refalt_code': 145752064, 'reversible': true}              │
    └──────────────────┴─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘

``` sql
SELECT regionkey_hex(regionkey('X', 1007, 1807, 1)) AS rkx,
       are_overlapping_regionkeys(
         regionkey('X', 1007, 1807, 1),
         parse_regionkey_hex('b80001f78000387a')
       ) AS overlaps;
```

    ┌──────────────────┬──────────┐
    │       rkx        │ overlaps │
    │     varchar      │ boolean  │
    ├──────────────────┼──────────┤
    │ b80001f78000387a │ true     │
    └──────────────────┴──────────┘

### Consequence annotation

Consequence annotation lives in the [DuckVEP
extension](https://github.com/RGenomicsETL/DuckVEP).

### Interval + reference helpers

``` sql
SELECT chrom, start, "end", name, block_count
FROM read_bed('test/data/targets.bed');
```

    ┌────────────────┬───────┬───────┬─────────┬─────────────┐
    │     chrom      │ start │  end  │  name   │ block_count │
    │    varchar     │ int64 │ int64 │ varchar │    int64    │
    ├────────────────┼───────┼───────┼─────────┼─────────────┤
    │ CHROMOSOME_I   │     0 │    10 │ target1 │           2 │
    │ CHROMOSOME_I   │    10 │    20 │ target2 │           1 │
    │ CHROMOSOME_II  │     0 │     8 │ target3 │        NULL │
    │ CHROMOSOME_III │     0 │     6 │ target4 │           1 │
    └────────────────┴───────┴───────┴─────────┴─────────────┘

``` sql
SELECT chrom, start, "end", pct_gc, num_a, num_c, num_g, num_t
FROM fasta_nuc('test/data/ce.fa', bed_path := 'test/data/targets.bed')
ORDER BY chrom, start;
```

    ┌────────────────┬───────┬───────┬────────┬───────┬───────┬───────┬───────┐
    │     chrom      │ start │  end  │ pct_gc │ num_a │ num_c │ num_g │ num_t │
    │    varchar     │ int64 │ int64 │ double │ int64 │ int64 │ int64 │ int64 │
    ├────────────────┼───────┼───────┼────────┼───────┼───────┼───────┼───────┤
    │ CHROMOSOME_I   │     0 │    10 │    0.6 │     2 │     4 │     2 │     2 │
    │ CHROMOSOME_I   │    10 │    20 │    0.5 │     4 │     3 │     2 │     1 │
    │ CHROMOSOME_II  │     0 │     8 │  0.625 │     2 │     4 │     1 │     1 │
    │ CHROMOSOME_III │     0 │     6 │    0.5 │     2 │     2 │     1 │     1 │
    └────────────────┴───────┴───────┴────────┴───────┴───────┴───────┴───────┘

``` sql
SELECT chrom, start, "end", seq_len, pct_gc
FROM fasta_nuc('test/data/ce.fa', bin_width := 10, region := 'CHROMOSOME_I:1-20');
```

    ┌──────────────┬───────┬───────┬─────────┬────────┐
    │    chrom     │ start │  end  │ seq_len │ pct_gc │
    │   varchar    │ int64 │ int64 │  int64  │ double │
    ├──────────────┼───────┼───────┼─────────┼────────┤
    │ CHROMOSOME_I │     0 │    10 │      10 │    0.6 │
    │ CHROMOSOME_I │    10 │    20 │      10 │    0.5 │
    └──────────────┴───────┴───────┴─────────┴────────┘

### cgranges registry entry points

`duckhts_cgranges_*` exposes a session-scoped immutable interval index
for native overlap queries. You can either build it row-wise with
`duckhts_cgranges_create(...)` + `duckhts_cgranges_add(...)`, or
bulk-load it from any table or view with
`duckhts_cgranges_from_table(...)`, which runs on your connection and so
also sees TEMP objects. For row-preserving filters or count annotations
over provider rows, use the vectorized scalar helpers
`duckhts_cgranges_has_overlap(...)` and
`duckhts_cgranges_count_overlaps(...)` directly in queries over
`read_bed(...)`, `read_bam(...)`, `read_bcf(...)`, or regular tables.
For streaming one-row-per-hit expansion while keeping provider columns,
use `duckhts_cgranges_overlaps_list(...)` with `UNNEST(...)` in the
SELECT list, which also covers bulk probing of any relation.

``` sql
SELECT duckhts_cgranges_create('readme_idx');
```

    ┌───────────────────────────────────────┐
    │ duckhts_cgranges_create('readme_idx') │
    │                boolean                │
    ├───────────────────────────────────────┤
    │ true                                  │
    └───────────────────────────────────────┘

``` sql
SELECT duckhts_cgranges_add('readme_idx', 'chr1', 10, 20, 'a');
```

    ┌─────────────────────────────────────────────────────────┐
    │ duckhts_cgranges_add('readme_idx', 'chr1', 10, 20, 'a') │
    │                         boolean                         │
    ├─────────────────────────────────────────────────────────┤
    │ true                                                    │
    └─────────────────────────────────────────────────────────┘

``` sql
SELECT duckhts_cgranges_add('readme_idx', 'chr1', 30, 40, 'b');
```

    ┌─────────────────────────────────────────────────────────┐
    │ duckhts_cgranges_add('readme_idx', 'chr1', 30, 40, 'b') │
    │                         boolean                         │
    ├─────────────────────────────────────────────────────────┤
    │ true                                                    │
    └─────────────────────────────────────────────────────────┘

``` sql
SELECT duckhts_cgranges_index('readme_idx');
```

    ┌──────────────────────────────────────┐
    │ duckhts_cgranges_index('readme_idx') │
    │               boolean                │
    ├──────────────────────────────────────┤
    │ true                                 │
    └──────────────────────────────────────┘

``` sql
SELECT interval_ordinal, label, interval_chrom, interval_start, interval_end
FROM duckhts_cgranges_overlaps('readme_idx', 'chr1', 35, 36, query_row_id := 7);
```

    ┌──────────────────┬─────────┬────────────────┬────────────────┬──────────────┐
    │ interval_ordinal │  label  │ interval_chrom │ interval_start │ interval_end │
    │      int64       │ varchar │    varchar     │     int32      │    int32     │
    ├──────────────────┼─────────┼────────────────┼────────────────┼──────────────┤
    │                1 │ b       │ chr1           │             30 │           40 │
    └──────────────────┴─────────┴────────────────┴────────────────┴──────────────┘

``` sql
CREATE TEMP VIEW readme_targets AS
SELECT * FROM (VALUES ('chr2', 100, 110, 'alpha'), ('chr2', 150, 170, 'beta'))
  AS t(chrom, start, "end", label);
```

``` sql
SELECT * FROM duckhts_cgranges_from_table(
  'readme_qry_idx', 'readme_targets', 'chrom', 'start', 'end', 'label'
);
```

    ┌─────────┐
    │ indexed │
    │ boolean │
    ├─────────┤
    │ true    │
    └─────────┘

``` sql
SELECT interval_ordinal, label, interval_chrom, interval_start, interval_end
FROM duckhts_cgranges_overlaps('readme_qry_idx', 'chr2', 140, 170, mode := 'contain');
```

    ┌──────────────────┬─────────┬────────────────┬────────────────┬──────────────┐
    │ interval_ordinal │  label  │ interval_chrom │ interval_start │ interval_end │
    │      int64       │ varchar │    varchar     │     int32      │    int32     │
    ├──────────────────┼─────────┼────────────────┼────────────────┼──────────────┤
    │                1 │ beta    │ chr2           │            150 │          170 │
    └──────────────────┴─────────┴────────────────┴────────────────┴──────────────┘

``` sql
CREATE TABLE readme_probes AS
SELECT * FROM (VALUES
  (10, 'chr2', 100, 105),
  (20, 'chr2', 160, 161),
  (30, 'chr2', 500, 510)
) AS t(probe_id, chrom, start, "end");
```

``` sql
SELECT
  probe_id,
  hit.interval_ordinal,
  hit.label,
  hit.label_type,
  hit.interval_chrom,
  hit.interval_start,
  hit.interval_end
FROM (
  SELECT
    p.probe_id,
    unnest(duckhts_cgranges_overlaps_list('readme_qry_idx', p.chrom, p.start, p."end")) AS hit
  FROM readme_probes AS p
)
ORDER BY probe_id, hit.interval_ordinal;
```

    ┌──────────┬──────────────────┬─────────┬────────────┬────────────────┬────────────────┬──────────────┐
    │ probe_id │ interval_ordinal │  label  │ label_type │ interval_chrom │ interval_start │ interval_end │
    │  int32   │      int64       │ varchar │  varchar   │    varchar     │     int32      │    int32     │
    ├──────────┼──────────────────┼─────────┼────────────┼────────────────┼────────────────┼──────────────┤
    │       10 │                0 │ alpha   │ VARCHAR    │ chr2           │            100 │          110 │
    │       20 │                1 │ beta    │ VARCHAR    │ chr2           │            150 │          170 │
    └──────────┴──────────────────┴─────────┴────────────┴────────────────┴────────────────┴──────────────┘

``` sql
SELECT duckhts_cgranges_destroy('readme_idx');
```

    ┌────────────────────────────────────────┐
    │ duckhts_cgranges_destroy('readme_idx') │
    │                boolean                 │
    ├────────────────────────────────────────┤
    │ true                                   │
    └────────────────────────────────────────┘

``` sql
SELECT duckhts_cgranges_destroy('readme_qry_idx');
```

    ┌────────────────────────────────────────────┐
    │ duckhts_cgranges_destroy('readme_qry_idx') │
    │                  boolean                   │
    ├────────────────────────────────────────────┤
    │ true                                       │
    └────────────────────────────────────────────┘

``` sql
DROP TABLE readme_probes;
```

``` sql
DROP VIEW readme_targets;
```

### Fixed-bin native counting

`bam_bin_counts()` does fixed-width read-start binning directly in
native code. This is the counting primitive used for WisecondorX-style
workflows: duplicate handling is explicit via `rmdup`, and optional
`stats := 'gc,mq'` adds one-pass GC and MAPQ summaries on the same scan.

``` sql
SELECT
  bin_id,
  count_total,
  count_fwd,
  count_rev,
  count_pre,
  printf('%.2f', gc_perc_pre) AS gc_pre,
  printf('%.2f', gc_perc_post) AS gc_post,
  printf('%.1f', mean_mapq_post) AS mean_mapq_post
FROM bam_bin_counts(
  'test/data/fixture_mixed.cram',
  5000,
  reference := 'test/data/fixture_ref.fa',
  rmdup := 'streaming',
  stats := 'gc,mq'
)
ORDER BY bin_id;
```

    ┌────────┬─────────────┬───────────┬───────────┬───────────┬─────────┬─────────┬────────────────┐
    │ bin_id │ count_total │ count_fwd │ count_rev │ count_pre │ gc_pre  │ gc_post │ mean_mapq_post │
    │ int64  │    int64    │   int64   │   int64   │   int64   │ varchar │ varchar │    varchar     │
    ├────────┼─────────────┼───────────┼───────────┼───────────┼─────────┼─────────┼────────────────┤
    │      0 │           2 │         1 │         1 │         4 │ 0.50    │ 0.00    │ 60.0           │
    │      1 │           2 │         1 │         1 │         2 │ 0.00    │ 0.00    │ 60.0           │
    │      2 │           1 │         1 │         0 │         2 │ 1.00    │ 1.00    │ 60.0           │
    │      3 │           0 │         0 │         0 │         0 │ NULL    │ NULL    │ NULL           │
    │      4 │           0 │         0 │         0 │         0 │ NULL    │ NULL    │ NULL           │
    │      5 │           0 │         0 │         0 │         0 │ NULL    │ NULL    │ NULL           │
    │      6 │           0 │         0 │         0 │         0 │ NULL    │ NULL    │ NULL           │
    │      7 │           0 │         0 │         0 │         0 │ NULL    │ NULL    │ NULL           │
    │      8 │           0 │         0 │         0 │         0 │ NULL    │ NULL    │ NULL           │
    │      9 │           0 │         0 │         0 │         0 │ NULL    │ NULL    │ NULL           │
    └────────┴─────────────┴───────────┴───────────┴───────────┴─────────┴─────────┴────────────────┘
      10 rows                                                                             8 columns

### Mosdepth-compatible coverage outputs

`duckhts_mosdepth()` writes mosdepth-style output files directly from
indexed BAM/CRAM input. The example below writes windowed fragment
coverage and then reads back the generated BED.gz output.

``` sql
SELECT success, summary_path, regions_path
FROM duckhts_mosdepth(
  '/tmp/duckhts_readme_mosdepth',
  'test/data/range.bam',
  chrom := 'CHROMOSOME_II',
  by := '1000',
  no_per_base := TRUE,
  fragment_mode := TRUE,
  use_median := TRUE,
  overwrite := TRUE
);
```

    ┌─────────┬───────────────────────────────────────────────────┬─────────────────────────────────────────────┐
    │ success │                   summary_path                    │                regions_path                 │
    │ boolean │                      varchar                      │                   varchar                   │
    ├─────────┼───────────────────────────────────────────────────┼─────────────────────────────────────────────┤
    │ true    │ /tmp/duckhts_readme_mosdepth.mosdepth.summary.txt │ /tmp/duckhts_readme_mosdepth.regions.bed.gz │
    └─────────┴───────────────────────────────────────────────────┴─────────────────────────────────────────────┘

``` sql
SELECT
  column0 AS chrom,
  CAST(column1 AS BIGINT) AS start,
  CAST(column2 AS BIGINT) AS "end",
  CAST(column3 AS DOUBLE) AS depth
FROM read_csv(
  '/tmp/duckhts_readme_mosdepth.regions.bed.gz',
  delim := '\t',
  header := FALSE,
  compression := 'gzip'
)
LIMIT 3;
```

    ┌───────────────┬───────┬───────┬────────┐
    │     chrom     │ start │  end  │ depth  │
    │    varchar    │ int64 │ int64 │ double │
    ├───────────────┼───────┼───────┼────────┤
    │ CHROMOSOME_II │     0 │  1000 │    0.0 │
    │ CHROMOSOME_II │  1000 │  2000 │    5.0 │
    │ CHROMOSOME_II │  2000 │  3000 │    3.0 │
    └───────────────┴───────┴───────┴────────┘

### Polygenic risk scoring

`bcftools_score` computes per-sample polygenic risk scores (PRS) from a
VCF/BCF and one or more GWAS summary statistics files, mirroring the
`bcftools +score` plugin API.

``` sql
-- Hard-call (GT) PRS — PLINK summary format
-- S1: 0×0.5  + 1×(−0.2) + 2×1.0 = 1.8
-- S2: 1×0.5  + 2×(−0.2) + 0×1.0 = 0.1
SELECT SAMPLE, round(score_summary, 3) AS prs
FROM bcftools_score(
  'test/data/score_input.vcf',
  'test/data/score_summary.tsv',
  use := 'GT',
  columns := 'PLINK'
);
```

    ┌─────────┬────────┐
    │ SAMPLE  │  prs   │
    │ varchar │ double │
    ├─────────┼────────┤
    │ S1      │    1.8 │
    │ S2      │    0.1 │
    └─────────┴────────┘

``` sql
-- Multi-PRS TSV/SSF scoring: multiple summary files in one genotype scan
SELECT SAMPLE,
       round(score_summary, 3) AS prs_a,
       round(score_summary_na, 3) AS prs_b
FROM bcftools_score(
  'test/data/score_input.vcf',
  ['test/data/score_summary.tsv', 'test/data/score_summary_na.tsv'],
  use := 'GT',
  columns := 'PLINK'
);
```

    ┌─────────┬────────┬────────┐
    │ SAMPLE  │ prs_a  │ prs_b  │
    │ varchar │ double │ double │
    ├─────────┼────────┼────────┤
    │ S1      │    1.8 │    2.0 │
    │ S2      │    0.1 │    0.5 │
    └─────────┴────────┴────────┘

``` sql
-- Dosage-based PRS (DS field) — fractional allele dosages from imputed data
-- S1: 0.1×0.5 + 0.8×(−0.2) + 1.8×1.0 = 1.69
-- S2: 1.0×0.5 + 1.9×(−0.2) + 0.2×1.0 = 0.32
SELECT SAMPLE, round(score_summary, 3) AS prs_ds
FROM bcftools_score(
  'test/data/score_dosage.vcf',
  'test/data/score_summary.tsv',
  use := 'DS',
  columns := 'PLINK'
);
```

    ┌─────────┬────────┐
    │ SAMPLE  │ prs_ds │
    │ varchar │ double │
    ├─────────┼────────┤
    │ S1      │   1.69 │
    │ S2      │   0.32 │
    └─────────┴────────┘

``` sql
-- GWAS-VCF multi-PRS: each FORMAT sample column becomes a separate PRS track
SELECT SAMPLE, round(PRS_A, 3) AS prs_a, round(PRS_B, 3) AS prs_b
FROM bcftools_score(
  'test/data/score_input.vcf',
  'test/data/score_gwas_summary.vcf',
  use := 'GT'
);
```

    ┌─────────┬────────┬────────┐
    │ SAMPLE  │ prs_a  │ prs_b  │
    │ varchar │ double │ double │
    ├─────────┼────────┼────────┤
    │ S1      │    1.8 │    1.0 │
    │ S2      │    0.1 │    0.3 │
    └─────────┴────────┴────────┘

### Liftover score-style rows

``` sql
SELECT src_chrom, src_pos, dest_chrom, dest_pos, dest_ref, dest_alt,
       mapped, reverse_complemented, reject_reason, note
FROM duckdb_liftover(
  '(VALUES
     (''chrF'', 2, ''C'', ''T''),
     (''chrR'', 2, ''A'', ''G''),
     (''chrF'', 11, ''A'', ''T'')
   ) AS t(chrom, pos, ref, alt)',
  'chrom',
  'pos',
  ref_col := 'ref',
  alt_col := 'alt',
  chain_path := 'test/data/liftover.chain',
  dst_fasta_ref := 'test/data/liftover_dst.fa',
  src_fasta_ref := 'test/data/liftover_src.fa'
);
```

    ┌───────────┬─────────┬────────────┬──────────┬──────────┬──────────┬─────────┬──────────────────────┬───────────────────┬─────────┐
    │ src_chrom │ src_pos │ dest_chrom │ dest_pos │ dest_ref │ dest_alt │ mapped  │ reverse_complemented │   reject_reason   │  note   │
    │  varchar  │  int64  │  varchar   │  int64   │ varchar  │ varchar  │ boolean │       boolean        │      varchar      │ varchar │
    ├───────────┼─────────┼────────────┼──────────┼──────────┼──────────┼─────────┼──────────────────────┼───────────────────┼─────────┤
    │ chrF      │       2 │ chrLiftF   │        2 │ C        │ T        │ true    │ false                │ NULL              │ NULL    │
    │ chrR      │       2 │ chrLiftR   │        9 │ T        │ C        │ true    │ true                 │ NULL              │ NULL    │
    │ chrF      │      11 │ NULL       │     NULL │ NULL     │ NULL     │ false   │ false                │ SourceRefMismatch │ NULL    │
    └───────────┴─────────┴────────────┴──────────┴──────────┴──────────┴─────────┴──────────────────────┴───────────────────┴─────────┘

### SIMD dispatch flow

DuckHTS uses explicit runtime SIMD dispatch for byte-oriented helper
kernels, starting with `seq_gc_content(...)`. `scalar` is always
available and is the portable baseline. Optional platform backends such
as `avx2` or `avx512` should be checked with
`duckhts_simd_backend_available(...)` before being requested. The `auto`
policy resolves each logical kernel independently from the current
compiled-and-CPU-supported capability mask; use
`duckhts_simd_kernel_info()` for the per-kernel result and
`SELECT backend FROM duckhts_simd_set_backend('auto')` to return to
runtime auto-detection.

``` sql
SELECT backend, selectable, compiled, cpu_supported, available, selected
FROM duckhts_simd_info();
```

    ┌──────────────┬────────────┬──────────┬───────────────┬───────────┬──────────┐
    │   backend    │ selectable │ compiled │ cpu_supported │ available │ selected │
    │   varchar    │  boolean   │ boolean  │    boolean    │  boolean  │ boolean  │
    ├──────────────┼────────────┼──────────┼───────────────┼───────────┼──────────┤
    │ scalar       │ true       │ true     │ true          │ true      │ false    │
    │ sse2         │ false      │ false    │ true          │ false     │ false    │
    │ sse41        │ false      │ false    │ true          │ false     │ false    │
    │ avx2         │ true       │ true     │ true          │ true      │ true     │
    │ avx512       │ true       │ true     │ false         │ false     │ false    │
    │ neon         │ true       │ false    │ false         │ false     │ false    │
    │ wasm_simd128 │ true       │ false    │ false         │ false     │ false    │
    └──────────────┴────────────┴──────────┴───────────────┴───────────┴──────────┘

``` sql
SELECT kernel, selected_backend, scalar_fallback
FROM duckhts_simd_kernel_info();
```

    ┌─────────────────┬──────────────────┬─────────────────┐
    │     kernel      │ selected_backend │ scalar_fallback │
    │     varchar     │     varchar      │     boolean     │
    ├─────────────────┼──────────────────┼─────────────────┤
    │ seq_base_counts │ avx2             │ false           │
    │ bam_nt16_counts │ avx2             │ false           │
    │ nt16_gc_counts  │ avx2             │ false           │
    │ fastq_qc        │ avx2             │ false           │
    └─────────────────┴──────────────────┴─────────────────┘

``` sql
SELECT backend AS selected_backend FROM duckhts_simd_set_backend('scalar');
```

    ┌──────────────────┐
    │ selected_backend │
    │     varchar      │
    ├──────────────────┤
    │ scalar           │
    └──────────────────┘

``` sql
SELECT
  duckhts_simd_requested_backend() AS requested_backend,
  duckhts_simd_backend() AS selected_backend,
  printf('%.3f', seq_gc_content('ACGTNNacgtnn')) AS gc_content;
```

    ┌───────────────────┬──────────────────┬────────────┐
    │ requested_backend │ selected_backend │ gc_content │
    │      varchar      │     varchar      │  varchar   │
    ├───────────────────┼──────────────────┼────────────┤
    │ scalar            │ scalar           │ 0.500      │
    └───────────────────┴──────────────────┴────────────┘

``` sql
SELECT backend IS NOT NULL AS restored_auto FROM duckhts_simd_set_backend('auto');
```

    ┌───────────────┐
    │ restored_auto │
    │    boolean    │
    ├───────────────┤
    │ true          │
    └───────────────┘

### Sequence utilities

``` sql
SELECT
  NAME,
  seq_hash_2bit(substr(SEQUENCE, 1, 12)) AS hash_2bit_prefix,
  seq_encode_4bit(substr(SEQUENCE, 1, 16)) AS codes,
  seq_decode_4bit(seq_encode_4bit(substr(SEQUENCE, 1, 16))) AS roundtrip
FROM read_fasta('test/data/ce.fa')
LIMIT 2;
```

    ┌───────────────┬──────────────────┬──────────────────────────────────────────────────┬──────────────────┐
    │     NAME      │ hash_2bit_prefix │                      codes                       │    roundtrip     │
    │    varchar    │      uint64      │                     uint8[]                      │     varchar      │
    ├───────────────┼──────────────────┼──────────────────────────────────────────────────┼──────────────────┤
    │ CHROMOSOME_I  │          9898352 │ [4, 2, 2, 8, 1, 1, 4, 2, 2, 8, 1, 1, 4, 2, 2, 8] │ GCCTAAGCCTAAGCCT │
    │ CHROMOSOME_II │          6038978 │ [2, 2, 8, 1, 1, 4, 2, 2, 8, 1, 1, 4, 2, 2, 8, 1] │ CCTAAGCCTAAGCCTA │
    └───────────────┴──────────────────┴──────────────────────────────────────────────────┴──────────────────┘

``` sql
SELECT
  NAME,
  MATE,
  seq_encode_4bit(substr(SEQUENCE, 1, 12)) AS codes,
  seq_decode_4bit(seq_encode_4bit(substr(SEQUENCE, 1, 12))) AS roundtrip
FROM read_fastq('test/data/interleaved.fq', interleaved := true)
LIMIT 2;
```

    ┌─────────────────────────────────┬────────┬──────────────────────────────────────┬──────────────┐
    │              NAME               │  MATE  │                codes                 │  roundtrip   │
    │             varchar             │ uint16 │               uint8[]                │   varchar    │
    ├─────────────────────────────────┼────────┼──────────────────────────────────────┼──────────────┤
    │ HS25_09827:2:1201:1505:59795#49 │      1 │ [2, 2, 4, 8, 8, 1, 4, 1, 4, 2, 1, 8] │ CCGTTAGAGCAT │
    │ HS25_09827:2:1201:1505:59795#49 │      2 │ [1, 1, 4, 4, 1, 1, 1, 4, 1, 1, 4, 4] │ AAGGAAAGAAGG │
    └─────────────────────────────────┴────────┴──────────────────────────────────────┴──────────────┘

### FASTQ quality decoding and fused QC

`read_fastq()` separates input interpretation from output
representation:

- `input_quality_encoding` tells DuckHTS how to decode FASTQ ASCII into
  numeric qualities. The default is modern `phred33`. Use `phred64`,
  `solexa64`, or `auto` only for legacy files.
- `quality_representation := 'phred'` returns canonical numeric
  qualities as `UTINYINT[]`.
- `quality_representation := 'string'` returns canonical Phred+33 text.
  For legacy inputs this means decode first, then re-encode as modern
  FASTQ text.

This makes the flow explicit:

1.  FASTQ text input is decoded according to `input_quality_encoding`.
2.  DuckHTS normalizes to numeric Phred values internally.
3.  Output is either raw numeric quality arrays (`phred`) or canonical
    Phred+33 text (`string`).

For BAM/CRAM, qualities are already stored as numeric values, so there
is no FASTQ text-encoding ambiguity on input.

Use `duckhts_fastq_qc(...)` for global and per-cycle quality-control
reductions. It consumes the projected sequence and canonical quality
strings in one bounded aggregate instead of creating one SQL row per
base. Expand only the compact cycle result when plotting or joining
per-cycle statistics.

``` sql
WITH q AS (
  SELECT duckhts_fastq_qc(SEQUENCE, QUALITY) AS qc
  FROM read_fastq('test/data/r1.fq')
)
SELECT
  qc.reads,
  qc.bases,
  qc.q30_bases,
  qc.max_read_length
FROM q;
```

    ┌────────┬────────┬───────────┬─────────────────┐
    │ reads  │ bases  │ q30_bases │ max_read_length │
    │ uint64 │ uint64 │  uint64   │     uint32      │
    ├────────┼────────┼───────────┼─────────────────┤
    │      5 │    500 │       475 │             100 │
    └────────┴────────┴───────────┴─────────────────┘

``` sql
WITH q AS (
  SELECT duckhts_fastq_qc(SEQUENCE, QUALITY) AS qc
  FROM read_fastq('test/data/r1.fq')
)
SELECT cycle.cycle, cycle.bases, cycle.quality_sum
FROM q, UNNEST(qc.cycles) AS u(cycle)
ORDER BY cycle.cycle
LIMIT 5;
```

    ┌────────┬────────┬─────────────┐
    │ cycle  │ bases  │ quality_sum │
    │ uint32 │ uint64 │   uint64    │
    ├────────┼────────┼─────────────┤
    │      1 │      5 │         169 │
    │      2 │      5 │         160 │
    │      3 │      5 │         166 │
    │      4 │      5 │         177 │
    │      5 │      5 │         185 │
    └────────┴────────┴─────────────┘

Use numeric quality arrays when the query genuinely needs the full
quality histogram:

``` sql
SELECT *
FROM detect_quality_encoding('test/data/legacy_phred64.fq');
```

    ┌─────────┬────────────────────┬────────────────────┬─────────────────┬──────────────────────────┬──────────────────┬──────────────┐
    │ format  │ observed_ascii_min │ observed_ascii_max │ records_sampled │   compatible_encodings   │ guessed_encoding │ is_ambiguous │
    │ varchar │       int64        │       int64        │      int64      │         varchar          │     varchar      │   boolean    │
    ├─────────┼────────────────────┼────────────────────┼─────────────────┼──────────────────────────┼──────────────────┼──────────────┤
    │ fastq   │                104 │                104 │               1 │ phred33,phred64,solexa64 │ phred64          │ true         │
    └─────────┴────────────────────┴────────────────────┴─────────────────┴──────────────────────────┴──────────────────┴──────────────┘

``` sql
WITH q AS (
  SELECT NAME, QUALITY
  FROM read_fastq(
    'test/data/r1.fq',
    quality_representation := 'phred'
  )
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
LIMIT 12;
```

    ┌───────┬───────┬─────────┐
    │  pos  │ phred │ n_reads │
    │ int64 │ uint8 │  int64  │
    ├───────┼───────┼─────────┤
    │     1 │    33 │       1 │
    │     1 │    34 │       4 │
    │     2 │    32 │       5 │
    │     3 │    33 │       4 │
    │     3 │    34 │       1 │
    │     4 │    34 │       2 │
    │     4 │    36 │       2 │
    │     4 │    37 │       1 │
    │     5 │    37 │       5 │
    │     6 │    38 │       5 │
    │     7 │    33 │       1 │
    │     7 │    35 │       1 │
    └───────┴───────┴─────────┘
      12 rows       3 columns

### Metadata + export/index helpers

``` sql
SELECT idx, raw
FROM read_hts_header('test/data/formatcols.vcf.gz', mode := 'raw')
LIMIT 3;
```

    ┌───────┬─────────────────────────────────────────────────────┐
    │  idx  │                         raw                         │
    │ int64 │                       varchar                       │
    ├───────┼─────────────────────────────────────────────────────┤
    │     0 │ ##fileformat=VCFv4.3                                │
    │     1 │ ##FILTER=<ID=PASS,Description="All filters passed"> │
    │     2 │ ##contig=<ID=1>                                     │
    └───────┴─────────────────────────────────────────────────────┘

``` sql
SELECT seqname, tid, index_type, chunk_beg_vo, chunk_end_vo
FROM read_hts_index_spans('test/data/formatcols.vcf.gz')
LIMIT 3;
```

    ┌─────────┬───────┬────────────┬──────────────┬──────────────┐
    │ seqname │  tid  │ index_type │ chunk_beg_vo │ chunk_end_vo │
    │ varchar │ int64 │  varchar   │    uint64    │    uint64    │
    ├─────────┼───────┼────────────┼──────────────┼──────────────┤
    │ 1       │     0 │ CSI        │     20381696 │     23789568 │
    └─────────┴───────┴────────────┴──────────────┴──────────────┘

``` sql
SELECT index_type, octet_length(raw) AS raw_bytes
FROM read_hts_index_raw('test/data/formatcols.vcf.gz');
```

    ┌────────────┬───────────┐
    │ index_type │ raw_bytes │
    │  varchar   │   int64   │
    ├────────────┼───────────┤
    │ CSI        │        30 │
    └────────────┴───────────┘

``` sql
COPY (
  SELECT chrom, start, "end", name
  FROM read_bed('test/data/targets.bed')
) TO '/tmp/duckhts_readme_targets.bed' (FORMAT CSV, DELIMITER '\t', HEADER FALSE);
```

``` sql
SELECT success, output_path, bytes_out
FROM bgzip('/tmp/duckhts_readme_targets.bed',
           output_path := '/tmp/duckhts_readme_targets.bed.gz',
           keep := TRUE,
           overwrite := TRUE);
```

    ┌─────────┬────────────────────────────────────┬───────────┐
    │ success │            output_path             │ bytes_out │
    │ boolean │              varchar               │   int64   │
    ├─────────┼────────────────────────────────────┼───────────┤
    │ true    │ /tmp/duckhts_readme_targets.bed.gz │       107 │
    └─────────┴────────────────────────────────────┴───────────┘

``` sql
SELECT success, output_path, bytes_out
FROM bgunzip('/tmp/duckhts_readme_targets.bed.gz',
             output_path := '/tmp/duckhts_readme_targets.roundtrip.bed',
             keep := TRUE,
             overwrite := TRUE);
```

    ┌─────────┬───────────────────────────────────────────┬───────────┐
    │ success │                output_path                │ bytes_out │
    │ boolean │                  varchar                  │   int64   │
    ├─────────┼───────────────────────────────────────────┼───────────┤
    │ true    │ /tmp/duckhts_readme_targets.roundtrip.bed │       106 │
    └─────────┴───────────────────────────────────────────┴───────────┘

``` sql
SELECT success, index_format, index_path
FROM bam_index('test/data/range.bam',
               index_path := '/tmp/duckhts_readme_range.bam.bai',
               threads := 1);
```

    ┌─────────┬──────────────┬───────────────────────────────────┐
    │ success │ index_format │            index_path             │
    │ boolean │   varchar    │              varchar              │
    ├─────────┼──────────────┼───────────────────────────────────┤
    │ true    │ BAI          │ /tmp/duckhts_readme_range.bam.bai │
    └─────────┴──────────────┴───────────────────────────────────┘

``` sql
SELECT success, index_format, index_path
FROM bcf_index('test/data/vcf_file.bcf',
               index_path := '/tmp/duckhts_readme_vcf_file.bcf.csi',
               threads := 1);
```

    ┌─────────┬──────────────┬──────────────────────────────────────┐
    │ success │ index_format │              index_path              │
    │ boolean │   varchar    │               varchar                │
    ├─────────┼──────────────┼──────────────────────────────────────┤
    │ true    │ CSI          │ /tmp/duckhts_readme_vcf_file.bcf.csi │
    └─────────┴──────────────┴──────────────────────────────────────┘

``` sql
SELECT success, index_format, index_path
FROM tabix_index('/tmp/duckhts_readme_targets.bed.gz',
                 preset := 'bed',
                 index_path := '/tmp/duckhts_readme_targets.bed.gz.tbi');
```

    ┌─────────┬──────────────┬────────────────────────────────────────┐
    │ success │ index_format │               index_path               │
    │ boolean │   varchar    │                varchar                 │
    ├─────────┼──────────────┼────────────────────────────────────────┤
    │ true    │ TBI          │ /tmp/duckhts_readme_targets.bed.gz.tbi │
    └─────────┴──────────────┴────────────────────────────────────────┘

### Multi-file queries

`hts_union_query` builds a `UNION ALL BY NAME` across files matching a
glob pattern. Because DuckDB’s `query()` cannot accept subquery
expressions, use the `SET VARIABLE` + `getvariable()` pattern:

``` sql
SET VARIABLE q = hts_union_query('read_fastq', 'test/data/r*.fq');
```

``` sql
SELECT filename, count(*) AS n
FROM query(getvariable('q'))
GROUP BY ALL
ORDER BY filename;
```

    ┌─────────────────┬───────┐
    │    filename     │   n   │
    │     varchar     │ int64 │
    ├─────────────────┼───────┤
    │ test/data/r1.fq │     5 │
    │ test/data/r2.fq │     5 │
    └─────────────────┴───────┘

Per-file parameters can be passed as the third argument (SQL literal):

``` sql
SET VARIABLE q = hts_union_query('read_bam', 'test/data/range.bam',
                                  'region := ''CHROMOSOME_I:1-1000''');
```

``` sql
SELECT count(*) AS n FROM query(getvariable('q'));
```

    ┌───────┐
    │   n   │
    │ int64 │
    ├───────┤
    │     2 │
    └───────┘

## Remote URLs and HTS_PATH

Remote URLs (S3/GCS/HTTP/S) can work in two htslib build modes:

1.  Dynamic plugin mode (`ENABLE_PLUGINS`): remote handlers are loaded
    from `HTS_PATH`.
2.  Static-handler mode (plugins disabled): handlers are compiled into
    `libhts` and `HTS_PATH` is not needed.

Use `HTS_PATH` only when you want dynamic plugin discovery (for example,
to point at an external htslib plugin directory). `Rduckhts` users only
need to call the public helper before the first HTS file is opened:

``` r
library(Rduckhts)
setup_hts_env()
```

Example (works in static-handler mode and plugin mode):

``` sql
SELECT CHROM, COUNT(*) AS n
FROM read_bcf('s3://1000genomes-dragen-v3.7.6/data/cohorts/gvcf-genotyper-dragen-3.7.6/hg19/3202-samples-cohort/3202_samples_cohort_gg_chr22.vcf.gz',
              region := 'chr22:16050000-16050500')
GROUP BY CHROM;
```

    ┌─────────┬───────┐
    │  CHROM  │   n   │
    │ varchar │ int64 │
    ├─────────┼───────┤
    │ chr22   │    11 │
    └─────────┴───────┘

For a direct DuckDB CLI process, set `HTS_PATH` explicitly before its
first HTS read, for example:

``` bash
export HTS_PATH=$(Rscript --quiet -e 'Rduckhts::setup_hts_env(); cat(Sys.getenv("HTS_PATH"),sep="")')
```

If htslib has already opened a file in the process, restart the process
after changing `HTS_PATH`; plugin discovery has already occurred.

If you don’t have htslib plugins installed locally, download the
prebuilt binaries from the r-universe-binaries GitHub release and point
`HTS_PATH` at the extracted htslib/libexec/htslib directory inside the
package bundle.
<https://github.com/RGenomicsETL/duckhts/releases/tag/r-universe-binaries>

### Browser wasm/webR HTTP backend

For browser wasm/webR builds, DuckHTS does **not** use htslib `libcurl`
for remote `http`/`https` access.

- The webR side-module path disables htslib `libcurl`/`S3`/`GCS`
  features because socket-based libcurl calls from a wasm side module
  are not reliable in the current webR runtime model.
- DuckHTS registers a package-owned htslib `hFILE` scheme handler
  implemented in `src/wasm_http_hfile.c` for `http` and `https`.
- This backend uses synchronous `XMLHttpRequest` from the worker for
  range reads, index probes, and seek behavior.

Browser constraints still apply:

- Remote hosts must allow CORS for the main file **and** index sidecars
  (`.tbi`/`.csi`), including range requests.
- Behavior can vary by browser and by server-side CORS policy changes
  over time.
- `ALL_PROXY` / websocket proxy settings do not affect this XHR backend.

Optional header/auth configuration for browser wasm can be provided from
JavaScript before loading/querying:

``` js
Module.duckhtsWasmHttpConfig = {
  headers: {
    Authorization: "Bearer <short-lived-token>",
    "X-Request-Source": "webr-local"
  },
  allowHosts: ["ftp.ebi.ac.uk", ".s3.amazonaws.com"],
  enforceHostAllowlist: true,
  withCredentials: false,
  allowInsecureAuth: false
};
```

Security behavior of this config:

- Headers are only attached when the URL hostname matches `allowHosts`.
- Requests are blocked for non-matching hosts when
  `enforceHostAllowlist: true`.
- `Authorization` is blocked on non-HTTPS URLs unless
  `allowInsecureAuth: true` is set explicitly.
- Credentials/cookies are only sent when `withCredentials: true` is set.

### Browser wasm/duckdb-wasm local setup

DuckHTS is also intended to run as a generic DuckDB community extension
in browser wasm hosts (not only webR).

Use the local duckdb-wasm setup to exercise that path end-to-end:

``` bash
./scripts/start_duckdb_wasm_local_test.sh
```

This setup uses a Docker-only build path via
`scripts/docker/duckdb-wasm-local.Dockerfile`.

The container pre-installs cache-friendly, pinned wasm build
dependencies (`emsdk`, `vcpkg`) so repeated local runs do minimal setup
work.

Builds run in an isolated mirror worktree (`.duckdb_wasm_docker_work`)
and copy back only the wasm extension artifact, so your host native
`build/` and `cmake_build/` trees remain available for normal native
development and testing.

Then open:

``` text
http://127.0.0.1:8001/scripts/duckdb-wasm-local-test.html
```

This setup loads `duckhts.duckdb_extension.wasm` in duckdb-wasm, runs
local HTTP reader checks, and lets you set/clear
`Module.duckhtsWasmHttpConfig` directly in the browser host runtime.

The setup stages `duckdb-browser.mjs`, `duckdb-browser-eh.worker.js`,
and `duckdb-eh.wasm` at the site root for same-origin runtime loading,
while using an import map for `apache-arrow` resolution.

### S3 credentials and configuration

The htslib S3 plugin supports credentials embedded in the URL or
provided via environment variables or standard credentials files. For
AWS-style credentials, the most common variables are:

- `AWS_ACCESS_KEY_ID`
- `AWS_SECRET_ACCESS_KEY`
- `AWS_SESSION_TOKEN` (optional, for temporary credentials)
- `AWS_DEFAULT_REGION`
- `AWS_PROFILE` / `AWS_DEFAULT_PROFILE`
- `AWS_SHARED_CREDENTIALS_FILE` (override credentials file location)

You can also configure htslib-specific settings like
`HTS_S3_ADDRESS_STYLE`, `HTS_S3_HOST`, and `HTS_S3_S3CFG` for
non-default S3 endpoints or path-style access.

See the htslib S3 plugin documentation for full details, URL syntax, and
short‑lived credentials support:
<https://www.htslib.org/doc/htslib-s3-plugin.html>

## Development

This README is rendered with
[`duckknit`](https://github.com/rundel/duckknit) to execute SQL snippets
in a persistent DuckDB session. It discovers `duckdb` on `PATH`;
`DUCKDB_CLI` selects an explicit executable path.

### Clone and environment setup

Clone the pinned extension build tools, then create the local Python
test environment and platform receipt:

``` bash
git clone --recurse-submodules https://github.com/RGenomicsETL/duckhts.git
cd duckhts
make configure
```

`make configure` is an explicit network bootstrap for the Python
SQLLogicTest runner. Normal native builds use the committed DuckDB
headers and vendored C sources; they do not download inputs. Note: MSVC
builds (windows_amd64/windows_arm64) are not supported. Use MinGW/RTools
for Windows.

### Prerequisites

Building the extension requires:

- C compiler (GCC or Clang)
- CMake ≥ 3.5
- Make
- Python 3 + venv
- Git
- [htslib](https://github.com/samtools/htslib) build dependencies: zlib,
  libbz2, liblzma, libdeflate, libcurl, libcrypto (OpenSSL)

Rendering the root documentation additionally requires R with
`rmarkdown` and [`duckknit`](https://github.com/rundel/duckknit).

On Debian/Ubuntu:

``` bash
sudo apt install build-essential cmake python3 python3-venv git \
    zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev libcurl4-openssl-dev libssl-dev
```

On macOS:

``` bash
brew install cmake htslib xz libdeflate
```

The clone already contains the pinned
[htslib](https://github.com/samtools/htslib) source. Vendoring scripts
are dependency-maintainer operations, not development setup.

### Build

``` bash
make configure    # one-time setup (Python venv, platform detection)
make release      # build optimised extension
```

The build runs [htslib](https://github.com/samtools/htslib)’s Makefile
(`make lib-static`) in-tree.

The extension binary is written to
`build/release/duckhts.duckdb_extension`.

### Debug build

``` bash
make debug
```

## Loading

``` sql
-- Unsigned extensions must be loaded with -unsigned flag:
-- duckdb -unsigned

LOAD '/path/to/duckhts.duckdb_extension';
```

## Testing

SQL tests live in `test/sql/` using [DuckDB](https://duckdb.org/)’s
SQLLogicTest format. The small fixtures and required indexes are
committed under `test/data/`; a fresh clone does not prepare or download
test data. The test target removes declared generated outputs after each
file.

``` bash
make test_release
```

### External benchmark and conformance data

Persistent external inputs and derived benchmark relations live outside
the clone under `$DUCKHTS_CACHE_DIR` (default `~/.cache/duckhts`).
Explicit staging commands download only missing inputs from their
recorded public source and write a nearby provenance TSV with the
source, release, transformation, and consumer path. Benchmark rendering
never downloads inputs while measuring.

``` bash
make stage-liftover-references
make stage-giab-v4.2.1
make stage-norm-1000g-dragen-gvcf
```

Use `DUCKHTS_CACHE_DIR=/path/to/cache` to relocate the complete
external-data cache.

## References

- DuckDB: <https://duckdb.org/>
- DuckDB Extension API: <https://duckdb.org/docs/extensions/overview>
- DuckDB extension template (C):
  <https://github.com/duckdb/extension-template-c>
- htslib: <https://github.com/samtools/htslib>
- RBCFTools: <https://github.com/RGenomicsETL/RBCFTools>
- duckknit: <https://github.com/rundel/duckknit>

## License

The DuckHTS DuckDB extension is licensed under the MIT License; see
[LICENSE](LICENSE). The R packages (Rduckhts, duckhtsbench) are licensed
GPL (\>= 2), and the `duckhts` npm package GPL-2.0-or-later. Vendored
and linked third-party code keeps its own licences: see
[r/Rduckhts/inst/COPYRIGHT](r/Rduckhts/inst/COPYRIGHT) and
[js/THIRD_PARTY_NOTICES.md](js/THIRD_PARTY_NOTICES.md).

## Credits

[![Contributors](https://contrib.rocks/image?repo=RGenomicsETL/duckhts)](https://github.com/RGenomicsETL/duckhts/graphs/contributors)

The GenBank reader and FASTA converter were contributed by [Ryan
Ward](https://github.com/ryandward) of [Nurture
Bio](https://github.com/Nurture-Bio).

Thanks to all
[contributors](https://github.com/RGenomicsETL/duckhts/graphs/contributors).
See the [package author credits](r/Rduckhts/DESCRIPTION) and
[third-party notices](r/Rduckhts/inst/COPYRIGHT) for upstream
acknowledgements.
