
<!-- README.md is generated from README.Rmd using [duckknit](https://github.com/rundel/duckknit). Please edit that file. -->

# DuckHTS

A [DuckDB](https://duckdb.org/) extension for reading high-throughput
sequencing (HTS) file formats and running interval-oriented genomics
workflows directly in SQL.

DuckHTS brings sequencing alignments, variants, reference sequences,
annotations, and interval text files into DuckDB through
[htslib](https://github.com/samtools/htslib). Query VCF, BCF, BAM, CRAM,
FASTA, FASTQ, GTF, GFF, BED, and tabix-indexed files with standard SQL,
then compose metadata inspection, compression/indexing, sequence
utilities, and interval/reference helpers in the same session.

It covers the sequencing side of the stack: raw/reference data access,
annotation parsing, BGZF/tabix utilities, and building blocks for
coverage, CNV, QC, and downstream export workflows.

> **Part of a growing genomics stack in DuckDB.** DuckHTS covers
> sequencing and reference-oriented formats.
> [PlinkingDuck](https://github.com/teaguesterling/plinking_duck)
> extends the same idea to PLINK genotype workflows and explicitly
> builds on DuckHTS for the sequencing side of the stack. See also its
> documentation: <https://plinking-duck.readthedocs.io> . Together,
> these tools are pushing DuckDB toward a broader genomics analysis
> environment built around portable file readers, SQL-native analytics,
> and composable local workflows. After all, most of Genomics File
> Formats are Tables/Arrays in disguise.

This extension does not support MSVC builds
(windows_amd64/windows_arm64). Use MinGW/RTools for Windows.

## Functions

<details>
<summary>
Show generated function catalog
</summary>

## Extension Function Catalog

This section is generated from `functions.yaml`.

### Diagnostics

| Function                             | Kind   | Returns                | R helper                              | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|--------------------------------------|--------|------------------------|---------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `duckhts_simd_backend`               | scalar | VARCHAR                | `rduckhts_simd_backend`               | Return the current DuckHTS SIMD dispatch label. For explicit scalar or concrete backend requests this is the requested policy; for auto it is the single selected backend when all logical kernels resolve to the same backend, or mixed when per-kernel auto-dispatch resolves to multiple backends. Use duckhts_simd_kernel_info() for per-kernel details.                                                                                                                                                                                                |
| `duckhts_simd_requested_backend`     | scalar | VARCHAR                | `rduckhts_simd_requested_backend`     | Return the current explicit SIMD backend request, usually auto unless `SELECT backend FROM duckhts_simd_set_backend('auto'\|'scalar'\|backend)` was called. The selected per-kernel backend may differ under auto-dispatch across x86, ARM, wasm, and scalar-only builds.                                                                                                                                                                                                                                                                                   |
| `duckhts_simd_backend_compiled`      | scalar | BOOLEAN                | `rduckhts_simd_backend_compiled`      | Return whether a concrete DuckHTS SIMD backend was compiled into this build. This is independent of whether the current CPU/runtime supports executing that backend; for example avx512 can be compiled but not CPU-supported on the running host.                                                                                                                                                                                                                                                                                                          |
| `duckhts_simd_backend_cpu_supported` | scalar | BOOLEAN                | `rduckhts_simd_backend_cpu_supported` | Return whether the current CPU/runtime supports a concrete DuckHTS SIMD backend, independent of whether DuckHTS compiled an implementation for it. Availability is the intersection of compiled and CPU-supported.                                                                                                                                                                                                                                                                                                                                          |
| `duckhts_simd_backend_available`     | scalar | BOOLEAN                | `rduckhts_simd_backend_available`     | Return whether a concrete SIMD backend is usable in the current process. Availability means the backend is compiled into DuckHTS and supported by the current CPU/runtime. auto is a selection request rather than a concrete backend and is not reported as available here.                                                                                                                                                                                                                                                                                |
| `duckhts_simd_info`                  | table  | table                  | `rduckhts_simd_info`                  | Return one row per known concrete DuckHTS SIMD backend with extension-owned selectable, compiled, CPU-supported, available, selected, requested, and dispatch-mode diagnostics. Availability is the intersection of compiled and CPU/runtime-supported. selectable reports whether the backend has a selectable implementation path; explicit selection still requires available = TRUE. selected is TRUE when the current dispatch table uses that backend for at least one logical kernel. auto is a selection request and is not a concrete backend row. |
| `duckhts_simd_kernel_info`           | table  | table                  | `rduckhts_simd_kernel_info`           | Return one row per logical DuckHTS SIMD kernel showing the concrete backend selected by the current immutable dispatch table, the selected capability, the requested backend policy, whether scalar was used as a per-kernel fallback, and the dispatch mode. This is the authoritative diagnostic for mixed auto-dispatch when different kernels resolve to different backends.                                                                                                                                                                            |
| `duckhts_simd_set_backend`           | table  | table(backend VARCHAR) | `rduckhts_simd_set_backend`           | Explicitly select the DuckHTS SIMD dispatch policy for this process using a one-row table-function call and return the current dispatch label in a backend column. Use auto for per-kernel runtime dispatch or scalar for a portable baseline; unavailable platform-specific requests such as avx512 on non-AVX-512 CPUs raise an error instead of silently falling back.                                                                                                                                                                                   |

### Readers

| Function          | Kind         | Returns | R helper                                                                                                                                                               | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
|-------------------|--------------|---------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `read_bcf`        | table        | table   | `rduckhts_bcf`                                                                                                                                                         | Read VCF and BCF variant data with typed INFO, FORMAT, typed CSQ/ANN/BCSQ subfields, optional tidy sample output, optional bcftools-style CSQ type overrides, and optional htslib decompression worker threads via decompression_threads (default 0 for single-threaded reads).                                                                                                                                                                                                                              |
| `read_bam`        | table        | table   | `rduckhts_bam`                                                                                                                                                         | Read SAM, BAM, and CRAM alignments with optional typed SAMtags and auxiliary tag maps. Use sequence_encoding := ‘nt16’ to return SEQ as UTINYINT\[\], quality_representation := ‘phred’ to return QUAL as UTINYINT\[\], and cigar_representation := ‘binary’ to return packed BAM CIGAR operations as UINTEGER\[\] instead of SAM text. decompression_threads controls per-file htslib worker threads and defaults to 2; use 0 to disable worker threads.                                                    |
| `read_fasta`      | table        | table   | `rduckhts_fasta`                                                                                                                                                       | Read FASTA records or indexed FASTA regions as sequence rows. Use sequence_encoding := ‘nt16’ to return SEQUENCE as UTINYINT\[\] (htslib nt16 4-bit codes) instead of VARCHAR. For bgzipped FASTA, gzi_path may point to an explicit .gzi sidecar when it is not colocated with the FASTA.                                                                                                                                                                                                                   |
| `read_bed`        | table        | table   | `rduckhts_bed`                                                                                                                                                         | Read BED3-BED12 interval files with canonical typed columns and optional tabix-backed region filtering.                                                                                                                                                                                                                                                                                                                                                                                                      |
| `fasta_nuc`       | table        | table   | `rduckhts_fasta_nuc`                                                                                                                                                   | Compute bedtools nuc-style nucleotide composition for supplied BED intervals or generated fixed-width bins over a FASTA reference. For bgzipped FASTA, gzi_path may point to an explicit .gzi sidecar when it is not colocated with the FASTA.                                                                                                                                                                                                                                                               |
| `read_fastq`      | table        | table   | `rduckhts_fastq`                                                                                                                                                       | Read single-end, paired-end, or interleaved FASTQ files with optional legacy quality decoding. By default, FASTQ qualities are interpreted as modern Phred+33 input. Use sequence_encoding := ‘nt16’ to return SEQUENCE as UTINYINT\[\] and quality_representation := ‘phred’ to return QUALITY as UTINYINT\[\] instead of VARCHAR. input_quality_encoding accepts ‘phred33’, ‘auto’, ‘phred64’, or ‘solexa64’.                                                                                              |
| `read_gff`        | table        | table   | `rduckhts_gff`                                                                                                                                                         | Read GFF annotations with optional raw scalar and richer list/pair parsed attribute columns, strict GFF3 structural validation, and indexed region filtering.                                                                                                                                                                                                                                                                                                                                                |
| `read_gtf`        | table        | table   | `rduckhts_gtf`                                                                                                                                                         | Read GTF annotations with optional raw scalar and richer list/pair parsed attribute columns and indexed region filtering.                                                                                                                                                                                                                                                                                                                                                                                    |
| `read_tabix`      | table        | table   | `rduckhts_tabix`                                                                                                                                                       | Read generic tabix-indexed text data with optional header handling and type inference.                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `fasta_index`     | table        | table   | `rduckhts_fasta_index`                                                                                                                                                 | Build a FASTA index (.fai) and return a single row with columns success (BOOLEAN) and index_path (VARCHAR).                                                                                                                                                                                                                                                                                                                                                                                                  |
| `hts_union_query` | scalar_macro | VARCHAR | `rduckhts_bam_multi, rduckhts_bcf_multi, rduckhts_fastq_multi, rduckhts_fasta_multi, rduckhts_bed_multi, rduckhts_tabix_multi, rduckhts_gff_multi, rduckhts_gtf_multi` | Generate a UNION ALL BY NAME query string that reads every file matching a glob pattern through the named reader function. The result includes a ‘filename’ column identifying the source file for each row. Assign to a variable with SET VARIABLE and execute via query(getvariable(…)). Optional params string is appended to each reader call. In R, use the typed rduckhts\_\*\_multi() helpers instead, which accept file vectors with optional per-file parameters and create DuckDB tables directly. |

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
| `read_hts_header`           | table       | table   | `rduckhts_hts_header`              | Inspect HTS headers in parsed, raw, or combined form across supported formats.                                                                                                                                                                                             |
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
| `bcftools_norm_row`         | scalar      | STRUCT   |                          | Normalize one variant row with bcftools/vt-style left-alignment semantics against a FASTA reference. The alt argument may be either a comma-delimited VARCHAR or a VARCHAR\[\] list. The returned STRUCT contains pos_normed, end_pos_normed, ref_normed, alt_normed (always VARCHAR\[\]), normed (TRUE/FALSE/NULL), and norm_status. Symbolic <DEL> rows can use end_pos, and symbolic <DUP> rows can use svlen.                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| `duckhts_bcftools_norm`     | table_macro | table    | `rduckhts_bcftools_norm` | DuckDB table macro wrapper over bcftools_norm_row that normalizes variants from a table or derived-table expression while preserving the original columns. The input ALT column may be either VARCHAR or VARCHAR\[\]. The result appends pos_normed, end_pos_normed, ref_normed, alt_normed, normed, and norm_status; with split_multiallelic := TRUE, multiallelic sites are split before normalization and alt_normed becomes VARCHAR plus alt_index.                                                                                                                                                                                                                                                                                                                                                                                                                      |
| `bcftools_score`            | table       | table    | `rduckhts_score`         | Compute polygenic scores from one genotype BCF/VCF and one or more summary-statistics files with bcftools +score-compatible GT/DS/HDS/AP/GP/AS dosage semantics, sample subsetting, and region/target/FILTER-string controls. The second argument accepts a scalar path or a DuckDB LIST/array of paths; TSV/SSF summaries produce one PRS column per file in a single genotype scan, while GWAS-VCF summaries still produce one PRS column per FORMAT sample. Use summaries_list_file with a NULL second argument to read paths from a file or directory; list-file entries are interpreted as written, matching upstream `bcftools +score --summaries` behavior, while directory inputs scan supported regular summary files in lexicographic order and ignore index sidecars. Use log_path to write per-PRS loaded/matched/allele-mismatch/duplicate-marker audit counts. |
| `bcftools_munge_row`        | scalar      | STRUCT   |                          | Normalize one summary-statistics row into GWAS-VCF-style fields (chrom/pos/ref/alt/effect metrics), resolving REF/ALT orientation against a FASTA reference and applying swap-aware sign/frequency/count transforms. The output flag `alleles_swapped` means REF/ALT orientation was swapped to match the FASTA reference.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `duckdb_munge`              | table_macro | table    | `rduckhts_munge`         | DuckDB macro wrapper over bcftools_munge_row that maps source columns (via preset or explicit map) and returns normalized GWAS-VCF-style rows with lean outputs and explicit `alleles_swapped` semantics. Output columns: chrom, pos, id, ref, alt, alleles_swapped, filter, ns, ez, nc, es, se, lp, af, ac, ne (16 columns). For METAL meta-analysis output with SI/I2/CQ/ED columns, use duckdb_munge_metal.                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `duckdb_munge_metal`        | table_macro | table    | `rduckhts_munge`         | Extended munge macro with METAL meta-analysis output columns. Same as duckdb_munge but additionally emits: si (imputation info, from INFO input), i2 (Cochran’s I² heterogeneity, from HET_I2), cq (Cochran’s Q -log10 p, from HET_LP or -log10(HET_P)), and ed (effect direction string, from DIRE; +/- flipped on allele swap). The R wrapper rduckhts_munge() auto-dispatches to this macro when metal keys (INFO, HET_I2, HET_P, HET_LP, DIRE) are present in the resolved column map.                                                                                                                                                                                                                                                                                                                                                                                   |

### Sequence UDFs

| Function          | Kind   | Returns      | R helper | Description                                                                                           |
|-------------------|--------|--------------|----------|-------------------------------------------------------------------------------------------------------|
| `seq_revcomp`     | scalar | VARCHAR      |          | Compute the reverse complement of a DNA sequence using A, C, G, T, and N bases.                       |
| `seq_canonical`   | scalar | VARCHAR      |          | Return the lexicographically smaller of a sequence and its reverse complement.                        |
| `seq_hash_2bit`   | scalar | UBIGINT      |          | Encode a short DNA sequence as a 2-bit unsigned integer hash.                                         |
| `seq_encode_4bit` | scalar | UTINYINT\[\] |          | Encode an IUPAC DNA sequence as a list of 4-bit base codes, preserving ambiguity symbols including N. |
| `seq_decode_4bit` | scalar | VARCHAR      |          | Decode a list of 4-bit IUPAC DNA base codes back into a sequence string.                              |
| `seq_gc_content`  | scalar | DOUBLE       |          | Compute GC fraction for a DNA sequence as a value between 0 and 1.                                    |
| `seq_kmers`       | table  | table        |          | Expand a sequence into positional k-mers with optional canonicalization.                              |

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

| Function                     | Kind   | Returns | R helper | Description                                                                                                         |
|------------------------------|--------|---------|----------|---------------------------------------------------------------------------------------------------------------------|
| `cigar_has_soft_clip`        | scalar | BOOLEAN |          | Test whether a CIGAR string contains any soft-clipped segment (`S`).                                                |
| `cigar_has_hard_clip`        | scalar | BOOLEAN |          | Test whether a CIGAR string contains any hard-clipped segment (`H`).                                                |
| `cigar_left_soft_clip`       | scalar | BIGINT  |          | Return the left-end soft-clipped length from a CIGAR string, or zero if the alignment does not start with `S`.      |
| `cigar_right_soft_clip`      | scalar | BIGINT  |          | Return the right-end soft-clipped length from a CIGAR string, or zero if the alignment does not end with `S`.       |
| `cigar_query_length`         | scalar | BIGINT  |          | Return the query-consuming length from a CIGAR string, counting `M`, `I`, `S`, `=`, and `X`.                        |
| `cigar_aligned_query_length` | scalar | BIGINT  |          | Return the aligned query length from a CIGAR string, counting `M`, `=`, and `X` but excluding clips and insertions. |
| `cigar_reference_length`     | scalar | BIGINT  |          | Return the reference-consuming length from a CIGAR string, counting `M`, `D`, `N`, `=`, and `X`.                    |
| `cigar_has_op`               | scalar | BOOLEAN |          | Test whether a CIGAR string contains at least one instance of the requested operator.                               |

</details>

`read_fastq` with `mate_path` requires exact QNAME pairing. `read_bam`
supports typed `standard_tags` and `auxiliary_tags` maps. `read_tabix`
supports header-aware parsing (`header`, `header_names`) and optional
type inference (`auto_detect`, `column_types`). Region lists in
comma-separated form are supported by `read_bam`, `read_bcf`,
`read_fasta`, `read_gff`, `read_gtf`, and `read_tabix`. `read_bam`
multi-region queries are deduplicated by htslib, while
`read_bcf`/`read_fasta`/`read_gff`/`read_gtf`/`read_tabix` chain regions
and can return duplicates for overlaps.

## Examples

The examples below run directly against bundled local test files through
the DuckDB CLI and load the built extension from
`build/release/duckhts.duckdb_extension`.

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
bulk-load it from SQL with `duckhts_cgranges_from_query(...)`. For
row-preserving filters or count annotations over provider rows, use the
vectorized scalar helpers `duckhts_cgranges_has_overlap(...)` and
`duckhts_cgranges_count_overlaps(...)` directly in queries over
`read_bed(...)`, `read_bam(...)`, `read_bcf(...)`, or regular tables.
For streaming one-row-per-hit expansion while keeping provider columns,
use `duckhts_cgranges_overlaps_list(...)` with `UNNEST(...)`. The older
`duckhts_cgranges_overlaps_bulk(...)` table function still accepts a
probe query and emits matching indexed intervals in one table-function
call; that bulk query runs on the extension-owned helper connection, so
use a regular table or view rather than a temp table.

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
SELECT duckhts_cgranges_from_query(
  'readme_qry_idx',
  'SELECT * FROM (VALUES (''chr2'', 100, 110, ''alpha''), (''chr2'', 150, 170, ''beta'')) AS t(chrom, start, "end", label)',
  'chrom', 'start', 'end', 'label'
);
```

    ┌────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐
    │ duckhts_cgranges_from_query('readme_qry_idx', 'SELECT * FROM (VALUES (''chr2'', 100, 110, ''alpha''), (''chr2'', 150, 170, ''beta'')) AS t(chrom, start, "end", label)', 'chrom', 'start', 'end', 'label') │
    │                                                                                                  boolean                                                                                                   │
    ├────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
    │ true                                                                                                                                                                                                       │
    └────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘

``` sql
SELECT duckhts_cgranges_index('readme_qry_idx');
```

    ┌──────────────────────────────────────────┐
    │ duckhts_cgranges_index('readme_qry_idx') │
    │                 boolean                  │
    ├──────────────────────────────────────────┤
    │ true                                     │
    └──────────────────────────────────────────┘

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
  p.probe_id,
  hit.interval_ordinal,
  hit.label,
  hit.label_type,
  hit.interval_chrom,
  hit.interval_start,
  hit.interval_end
FROM readme_probes AS p
CROSS JOIN UNNEST(
  duckhts_cgranges_overlaps_list('readme_qry_idx', p.chrom, p.start, p."end")
) AS u(hit)
ORDER BY p.probe_id, hit.interval_ordinal;
```

    ┌──────────┬──────────────────┬─────────┬────────────┬────────────────┬────────────────┬──────────────┐
    │ probe_id │ interval_ordinal │  label  │ label_type │ interval_chrom │ interval_start │ interval_end │
    │  int32   │      int64       │ varchar │  varchar   │    varchar     │     int32      │    int32     │
    ├──────────┼──────────────────┼─────────┼────────────┼────────────────┼────────────────┼──────────────┤
    │       10 │                0 │ alpha   │ VARCHAR    │ chr2           │            100 │          110 │
    │       20 │                1 │ beta    │ VARCHAR    │ chr2           │            150 │          170 │
    └──────────┴──────────────────┴─────────┴────────────┴────────────────┴────────────────┴──────────────┘

``` sql
SELECT query_row_id, interval_ordinal, label, interval_chrom, interval_start, interval_end
FROM duckhts_cgranges_overlaps_bulk(
  'readme_qry_idx',
  'SELECT probe_id, chrom, start, "end" FROM readme_probes',
  'chrom', 'start', 'end',
  query_row_id_col := 'probe_id'
)
ORDER BY query_row_id, interval_ordinal;
```

    ┌──────────────┬──────────────────┬─────────┬────────────────┬────────────────┬──────────────┐
    │ query_row_id │ interval_ordinal │  label  │ interval_chrom │ interval_start │ interval_end │
    │    int64     │      int64       │ varchar │    varchar     │     int32      │    int32     │
    ├──────────────┼──────────────────┼─────────┼────────────────┼────────────────┼──────────────┤
    │           10 │                0 │ alpha   │ chr2           │            100 │          110 │
    │           20 │                1 │ beta    │ chr2           │            150 │          170 │
    └──────────────┴──────────────────┴─────────┴────────────────┴────────────────┴──────────────┘

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

### FASTQ quality decoding and per-position histograms

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
GROUP BY ALL;
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
to point at an external htslib plugin directory).

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

If you need dynamic plugin mode, set `HTS_PATH` before loading the
extension, for example:

``` bash
export HTS_PATH=$(Rscript --quiet -e 'cat(Rduckhts:::duckhts_htslib_plugins_dir(),sep="")')
```

If `HTS_PATH` is changed after loading, restart the session and reload
the extension.

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
local HTTP reader smoke tests, and lets you set/clear
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
in a persistent DuckDB session.

### Environment setup

Run the one-time configure step to create the Python venv and detect
platform settings:

``` bash
make configure
```

Note: MSVC builds (windows_amd64/windows_arm64) are not supported. Use
MinGW/RTools for Windows.

### Prerequisites

- C compiler (GCC or Clang)
- CMake ≥ 3.5
- Make
- Python 3 + venv
- R with the `rmarkdown` and
  [`duckknit`](https://github.com/rundel/duckknit) packages.
- Git
- [htslib](https://github.com/samtools/htslib) build dependencies: zlib,
  libbz2, liblzma, libdeflate, libcurl, libcrypto (OpenSSL)

On Debian/Ubuntu:

``` bash
sudo apt install build-essential cmake python3 python3-venv git \
    zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev libcurl4-openssl-dev libssl-dev
```

On macOS:

``` bash
brew install cmake htslib xz libdeflate
```

### Vendor [htslib](https://github.com/samtools/htslib)

``` bash
./scripts/vendor_htslib.sh
```

This downloads and verifies [htslib](https://github.com/samtools/htslib)
1.23.1 into `third_party/htslib/`.

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
SQLLogicTest format.

Before running tests for the first time, prepare the indexed test data:

``` bash
./test/scripts/prepare_test_data.sh   # requires samtools, bcftools, bgzip, tabix
```

This copies files from the vendored
[htslib](https://github.com/samtools/htslib) test suite into
`test/data/` and builds the required indexes (BAI, CSI, TBI, FAI) so
region queries work.

Then run:

``` bash
make test_release
```

## R demo

The R package lives under `r/Rduckhts` and provides helpers to load the
extension and create [DuckDB](https://duckdb.org/) tables from HTS
files. See its README for R-specific usage:
[r/Rduckhts/README.Rmd](r/Rduckhts/README.Rmd).

``` r
library(DBI)
library(duckdb)

drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
con <- dbConnect(drv, dbdir = ":memory:")
ext_path <- normalizePath("build/release/duckhts.duckdb_extension", mustWork = FALSE)
dbExecute(con, sprintf("LOAD '%s'", ext_path))
#> [1] 0

dbGetQuery(con, "
  SELECT *
  FROM read_bcf('test/data/formatcols.vcf.gz', tidy_format := true)
  LIMIT 5
")
#>   CHROM POS ID REF ALT QUAL FILTER SAMPLE_ID  FORMAT_S
#> 1     1 100  a   A   T   NA   PASS        S1         a
#> 2     1 100  a   A   T   NA   PASS        S²   bbbbbbb
#> 3     1 100  a   A   T   NA   PASS        S3 ccccccccc

parquet_path <- tempfile(fileext = ".parquet")
dbExecute(con, sprintf(
  "COPY (SELECT * FROM read_bcf('test/data/formatcols.vcf.gz', tidy_format := true)) TO '%s' (FORMAT PARQUET)",
  parquet_path
))
#> [1] 3
file.exists(parquet_path)
#> [1] TRUE

dbGetQuery(con, "
  SELECT NAME, SEQUENCE, QUALITY, MATE, PAIR_ID
  FROM read_fastq('test/data/r1.fq', mate_path := 'test/data/r2.fq')
  LIMIT 5
")
#>                              NAME
#> 1 HS25_09827:2:1201:1505:59795#49
#> 2 HS25_09827:2:1201:1505:59795#49
#> 3 HS25_09827:2:1201:1559:70726#49
#> 4 HS25_09827:2:1201:1559:70726#49
#> 5 HS25_09827:2:1201:1564:39627#49
#>                                                                                               SEQUENCE
#> 1 CCGTTAGAGCATTTGTTGAAAATGCTTTCCTTGCTCCATGTGATGACTCTGGTGCCCTTGTCAAAAGCCAGCTGGGCCTATTCGTGTGGGTCTGTTTCTG
#> 2 AAGGAAAGAAGGGAGGGAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAAGTAGGAAGAATTCATCTACCCAATT
#> 3 TTGTTAAAATGACCATACCCAAAGTGATCTACAGACTCAATACAATTTCTATTGAAATACCAATCACACTCTTCACAGAACTAGAAAAACAGTTCTAAAA
#> 4 TTTTCTTTTATTAATTTTATACTTACATTTAAGTCTTTATTCCATTTTGAGTCAATGTTTGTATATGATGAGAGATAGGGGTCTAGTTTCATACTTCTAC
#> 5 ACGCGGCAATCCAATGTGTGAGTTGAGAAGCGGTGAGGAGGGAATCCTAATTTTATGAGCAGGTCAGGACCGTGGGAGATACCTGACACCTGAGATGGTA
#>                                                                                                QUALITY
#> 1 CABCFGDEEFFEFHGHGGFFGDIGIJFIFHHGHEIFGHBCGHDIFBE9GIAICGGICFIBFGGHGDGGGHE?GIGDFGGHEGIEJG>;FG<GGHACEFGH
#> 2 <CBB>DCHFEFBHAGCGACF7CJI8HBIIEFGFEBG?DCGA?ACFGGI=BEDG?EFEHFFFEHFD?HG+DFH>FFHGFBFE4F@I3HF@>A5F?GFH<EG
#> 3 CABEFGFFGFHGGGGJGGFFGKIHHJFIEHHHGIEGGEHJGHDHFGHIGICIJEFIFGIF8GGHKFHGGFEI6GGGFIGHGGIE>EFCFHGGGHEJEAJE
#> 4 ;CBCEFDHDGFGHDGDIGEF@EJIIGEEIECGFHGFHGGGHHHHGGKIFFEHGEGHFIEFFHHGDHHGJEGF?FBHFFGCHHFFII>GCFCFFGGCEBF?
#> 5 BACCFGBFGFHGGJGHGGFEGHIGIJHFEH:HHEHGHHBGGH9IAGHGFHIFJFFAFGIFDIGHKEIG<C>F,CGD66?7EFI5EEG>EGGGGD5=HH6E
#>   MATE                         PAIR_ID
#> 1    1 HS25_09827:2:1201:1505:59795#49
#> 2    2 HS25_09827:2:1201:1505:59795#49
#> 3    1 HS25_09827:2:1201:1559:70726#49
#> 4    2 HS25_09827:2:1201:1559:70726#49
#> 5    1 HS25_09827:2:1201:1564:39627#49

dbGetQuery(con, "
  SELECT QNAME, RNAME, POS, READ_GROUP_ID, SAMPLE_ID
  FROM read_bam('test/data/rg.sam.gz')
  LIMIT 5
")
#>   QNAME RNAME POS READ_GROUP_ID SAMPLE_ID
#> 1    a1    xx   1            x1        x1
#> 2    b1    xx   1            x2        x2
#> 3    c1    xx   1          <NA>      <NA>
#> 4    a2    xx  11            x1        x1
#> 5    b2    xx  11            x2        x2

dbGetQuery(con, "
  SELECT idx, raw
  FROM read_hts_header('test/data/formatcols.vcf.gz', mode := 'raw')
  LIMIT 3
")
#>   idx                                                 raw
#> 1   0                                ##fileformat=VCFv4.3
#> 2   1 ##FILTER=<ID=PASS,Description="All filters passed">
#> 3   2                                     ##contig=<ID=1>

dbGetQuery(con, "
  SELECT seqname, tid, index_type, chunk_beg_vo, chunk_end_vo
  FROM read_hts_index_spans('test/data/formatcols.vcf.gz')
  LIMIT 3
")
#>   seqname tid index_type chunk_beg_vo chunk_end_vo
#> 1       1   0        CSI     20381696     23789568

dbGetQuery(con, "
  SELECT index_type, octet_length(raw) AS raw_bytes
  FROM read_hts_index_raw('test/data/formatcols.vcf.gz')
")
#>   index_type raw_bytes
#> 1        CSI        30
```

### SAMtags + auxiliary tags

Standard SAMtags can be surfaced as typed columns and non-standard tags
captured in a map for ad hoc access:

``` r
dbGetQuery(con, "
  SELECT RG, NM, map_extract(AUXILIARY_TAGS, 'XZ') AS XZ
  FROM read_bam('test/data/aux_tags.sam.gz', standard_tags := true, auxiliary_tags := true)
  LIMIT 1
")
#>   RG NM  XZ
#> 1 x1  2 foo
```

### GFF/GTF annotation attributes

DuckHTS reads GFF3 with `read_gff(...)`, GTF with `read_gtf(...)`, and
generic GFF/GTF-like tabix text with `read_tabix(...)`.
`read_gff(..., strict := true)` performs GFF3 structural validation.
Attribute decoding can be scalar and raw for legacy convenience
(`attributes_map`), grouped and lossless for multi-values
(`attributes_list` as `MAP(VARCHAR, VARCHAR[])`), or exact parser-style
pairs (`attributes_pairs` as `LIST<STRUCT(key, value, idx)>`).

The GFF3 implementation is audited against
[GFFBase](https://github.com/Kuanhao-Chao/gffbase) in
[`benchmarks/benchmark_gffbase_conformance.md`](benchmarks/benchmark_gffbase_conformance.md),
including Rust/Python parser parity checks, local conformance cases,
server specs, and human-scale GENCODE timing.

``` r
dbGetQuery(con, "
  SELECT
    seqname,
    feature,
    start,
    \"end\",
    list_extract(map_extract_value(attributes_list, 'Dbxref'), 1) AS first_dbxref,
    list_count(attributes_pairs) AS n_attr_pairs
  FROM read_gff('test/data/gff_attrs.gff3',
                strict := true,
                attributes_list := true,
                attributes_pairs := true)
")
#>   seqname feature start end first_dbxref n_attr_pairs
#> 1    chr1    gene     1  10     GeneID:1            6

dbGetQuery(con, "
  SELECT p.key, p.value, p.idx
  FROM read_gff('test/data/gff_attrs.gff3', attributes_pairs := true),
       UNNEST(attributes_pairs) AS u(p)
  WHERE p.key IN ('Alias', 'Note')
  ORDER BY p.key, p.idx
")
#>     key       value idx
#> 1 Alias           a   0
#> 2 Alias           b   1
#> 3  Note hello world   0

dbGetQuery(con, "
  SELECT list_extract(map_extract_value(attributes_list, 'note'), 1) AS note
  FROM read_gtf('test/data/gtf_attrs.gtf', attributes_list := true)
")
#>          note
#> 1 weird; semi
```

### Tabix text files

``` r
dbGetQuery(con, "
  SELECT column0, column1
  FROM read_tabix('test/data/meta_tabix.tsv.gz')
  LIMIT 2
")
#>   column0 column1
#> 1    chr1       1
#> 2    chr1       2

dbGetQuery(con, "
  SELECT chrom, pos
  FROM read_tabix('test/data/header_tabix.tsv.gz', header := true)
  LIMIT 2
")
#>   chrom pos
#> 1  chr1   1
#> 2  chr1   2

dbGetQuery(con, "
  SELECT typeof(column1) AS column1_type
  FROM read_tabix('test/data/meta_tabix.tsv.gz', auto_detect := true)
  LIMIT 1
")
#>   column1_type
#> 1       BIGINT

dbGetQuery(con, "
  SELECT pos + 1 AS pos_plus_one
  FROM read_tabix('test/data/header_tabix.tsv.gz', header := true,
                  column_types := ['VARCHAR','BIGINT','VARCHAR'])
  LIMIT 1
")
#>   pos_plus_one
#> 1            2

dbDisconnect(con, shutdown = TRUE)
```

## Project Structure

    src/
      duckhts.c          # Extension entry point
      bcf_reader.c       # VCF/BCF reader (read_bcf)
      bam_reader.c       # SAM/BAM/CRAM reader (read_bam)
      seq_reader.c       # FASTA/FASTQ reader (read_fasta, read_fastq)
      tabix_reader.c     # Tabix/GTF/GFF reader (read_tabix, read_gtf, read_gff)
      vep_parser.c       # VEP/CSQ annotation parser
      ......
      include/
        vcf_types.h
        vep_parser.h
        .....
    third_party/
      htslib/            # Vendored htslib 1.23.1 (built automatically)
    test/
      sql/               # SQL logic tests
    duckdb_capi/
      duckdb.h           # DuckDB C API headers
      duckdb_extension.h
    r/
      Rduckhts/          # R package 

## References

- DuckDB: <https://duckdb.org/>
- DuckDB Extension API: <https://duckdb.org/docs/extensions/overview>
- DuckDB extension template (C):
  <https://github.com/duckdb/extension-template-c>
- htslib: <https://github.com/samtools/htslib>
- RBCFTools: <https://github.com/RGenomicsETL/RBCFTools>
- duckknit: <https://github.com/rundel/duckknit>

## License

MIT
