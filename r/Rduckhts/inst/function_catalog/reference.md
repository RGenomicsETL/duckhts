# Extension Function Reference

Generated from `functions.yaml`.

## duckhts_macro_definitions

Export the ordered DuckHTS macro definitions for connection-local installation.

Signature:

```sql
duckhts_macro_definitions()
```

Returns:

```
TABLE
```

### Installation

LOAD creates macros in a writable in-memory default database. For file-backed or read-only default databases LOAD creates no macros: execute the exported TEMP statements in install_order on each connection that needs them (R: rduckhts_install_macros). The decision concerns the database instance's default database; a connection that USEs another catalog before or after LOAD installs the TEMP statements too. TEMP macros shadow, but do not delete, persistent macros in older files. The definitions_sha256 value identifies the entire ordered definition set.

### Columns

install_order UINTEGER, name VARCHAR, public BOOLEAN (listed in functions.yaml), sql VARCHAR (CREATE OR REPLACE TEMP MACRO statement), definitions_sha256 VARCHAR.

### Examples

```sql
SELECT install_order, name, sql FROM duckhts_macro_definitions() ORDER BY install_order;
```

## duckhts_htslib_version

Return the runtime version reported by the htslib library loaded with DuckHTS. Rduckhts uses this value to reject a downstream linking receipt whose source/header version does not match the loaded library.

Signature:

```sql
duckhts_htslib_version()
```

Returns:

```
VARCHAR
```

### Examples

```sql
SELECT duckhts_htslib_version();
```

## duckhts_htslib_features

Return the htslib runtime feature bitfield reported by hts_features(). Use duckhts_htslib_feature_string() for the corresponding build description.

Signature:

```sql
duckhts_htslib_features()
```

Returns:

```
UINTEGER
```

### Examples

```sql
SELECT duckhts_htslib_features();
```

## duckhts_htslib_feature_string

Return htslib's runtime build-feature description, including configured transports, compression libraries, compiler, and build flags. DuckHTS snapshots it once while loading the extension so parallel SQL calls read immutable text.

Signature:

```sql
duckhts_htslib_feature_string()
```

Returns:

```
VARCHAR
```

### Examples

```sql
SELECT duckhts_htslib_feature_string();
```

## duckhts_simd_backend

Return the current DuckHTS SIMD dispatch label. For explicit scalar or concrete backend requests this is the requested policy; for auto it is the single selected backend when all logical kernels resolve to the same backend, or mixed when per-kernel auto-dispatch resolves to multiple backends. Use duckhts_simd_kernel_info() for per-kernel details.

Signature:

```sql
duckhts_simd_backend()
```

Returns:

```
VARCHAR
```

### Examples

```sql
SELECT duckhts_simd_backend();
```

## duckhts_simd_requested_backend

Return the current explicit SIMD backend request, usually auto unless `SELECT backend FROM duckhts_simd_set_backend('auto'|'scalar'|backend)` was called. The selected per-kernel backend may differ under auto-dispatch across x86, ARM, wasm, and scalar-only builds.

Signature:

```sql
duckhts_simd_requested_backend()
```

Returns:

```
VARCHAR
```

### Examples

```sql
SELECT duckhts_simd_requested_backend();
```

## duckhts_simd_backend_compiled

Return whether a concrete DuckHTS SIMD backend was compiled into this build. This is independent of whether the current CPU/runtime supports executing that backend; for example avx512 can be compiled but not CPU-supported on the running host.

Signature:

```sql
duckhts_simd_backend_compiled(backend)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT duckhts_simd_backend_compiled('scalar');
```

```sql
SELECT duckhts_simd_backend_compiled('avx512');
```

## duckhts_simd_backend_cpu_supported

Return whether the current CPU/runtime supports a concrete DuckHTS SIMD backend, independent of whether DuckHTS compiled an implementation for it. Availability is the intersection of compiled and CPU-supported.

Signature:

```sql
duckhts_simd_backend_cpu_supported(backend)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT duckhts_simd_backend_cpu_supported('avx2');
```

```sql
SELECT duckhts_simd_backend_cpu_supported('avx512');
```

## duckhts_simd_backend_available

Return whether a concrete SIMD backend is usable in the current process. Availability means the backend is compiled into DuckHTS and supported by the current CPU/runtime. auto is a selection request rather than a concrete backend and is not reported as available here.

Signature:

```sql
duckhts_simd_backend_available(backend)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT duckhts_simd_backend_available('scalar');
```

```sql
SELECT duckhts_simd_backend_available('avx2');
```

```sql
SELECT duckhts_simd_backend_available('avx512');
```

## duckhts_simd_info

Report compiled, runtime-supported and selected status for each concrete DuckHTS SIMD backend.

Signature:

```sql
duckhts_simd_info()
```

Returns:

```
table
```

### Diagnostics

Rows include selectable, compiled, CPU-supported, available, selected, requested and dispatch-mode fields. available requires both compiled and CPU/runtime-supported. Explicit selection requires available=TRUE and a selectable implementation path. selected means at least one logical kernel uses that backend. auto is a request, not a concrete backend row.

### Examples

```sql
SELECT * FROM duckhts_simd_info();
```

```sql
SELECT backend FROM duckhts_simd_info() WHERE available;
```

## duckhts_simd_kernel_info

Return one row per logical DuckHTS SIMD kernel showing the concrete backend selected by the current immutable dispatch table, the selected capability, the requested backend policy, whether scalar was used as a per-kernel fallback, and the dispatch mode. This is the authoritative diagnostic for mixed auto-dispatch when different kernels resolve to different backends.

Signature:

```sql
duckhts_simd_kernel_info()
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM duckhts_simd_kernel_info();
```

```sql
SELECT kernel, selected_backend FROM duckhts_simd_kernel_info();
```

## duckhts_simd_set_backend

Explicitly select the DuckHTS SIMD dispatch policy for this process using a one-row table-function call and return the current dispatch label in a backend column. Use auto for per-kernel runtime dispatch or scalar for a portable baseline; unavailable platform-specific requests such as avx512 on non-AVX-512 CPUs raise an error instead of silently falling back.

Signature:

```sql
duckhts_simd_set_backend(backend)
```

Returns:

```
table(backend VARCHAR)
```

### Examples

```sql
SELECT backend FROM duckhts_simd_set_backend('scalar');
```

```sql
SELECT backend FROM duckhts_simd_set_backend('auto');
```

## duckhts_duckdb_type_supported

Return whether the currently open DuckDB runtime advertises a logical type with the given name through duckdb_types(). This is a catalog-level runtime probe for feature gating SQL/macros across DuckDB versions.

Signature:

```sql
duckhts_duckdb_type_supported(candidate_type_name)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT duckhts_duckdb_type_supported('VARIANT');
```

```sql
SELECT duckhts_duckdb_type_supported('GEOMETRY');
```

## duckhts_duckdb_supports_variant

Return whether the currently open DuckDB runtime advertises the VARIANT logical type. Use this to gate optional SQL that depends on DuckDB VARIANT support.

Signature:

```sql
duckhts_duckdb_supports_variant()
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT duckhts_duckdb_supports_variant();
```

## duckhts_duckdb_supports_geometry

Return whether the currently open DuckDB runtime advertises the GEOMETRY logical type. Use this to gate optional SQL that depends on DuckDB GEOMETRY support.

Signature:

```sql
duckhts_duckdb_supports_geometry()
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT duckhts_duckdb_supports_geometry();
```

## read_bcf

Read VCF/BCF with header-typed INFO/FORMAT, typed CSQ/ANN/BCSQ annotations, sample selection and optional tidy sample rows.

Signature:

```sql
read_bcf(path, region := NULL, index_path := NULL, tidy_format := FALSE, additional_csq_column_types := NULL, scan_mode := 'auto', decompression_threads := 0, decode_error_policy := 'null', samples := NULL, regions := NULL)
```

Returns:

```
table
```

### Types

INFO/FORMAT Type and Number determine SQL types and scalar/list shape, including nonstandard declarations for named tags; no schema repair is inferred. additional_csq_column_types overrides typed annotation subfields. Missing numeric/string list positions are NULL; numeric vector-end padding terminates lists. Absent fields are NULL; explicit numeric missing lists are [NULL], without inferred allele/genotype cardinality. Tidy rows repeat complete record annotation across output chunks.

### FILTER

FILTER is VARCHAR[] in header order. PASS is [PASS], an unapplied dot is NULL, and named failures retain their identifiers.

### Errors

decode_error_policy is null, warn or error for header/payload type mismatches and oversized numeric scalars. Multiple scalar elements before vector-end padding count as malformed, including missing elements. FORMAT null/warn withholds that tag for every selected sample on the record.

### Samples

HTSlib selectors: NULL/'-' keeps all, empty string keeps none, comma-separated names include, leading '^' excludes. Unknown names error. Selection retains header order; read_bcf_samples() supplies original-header indices/names.

### Regions

Comma-separated indexed regions use a union iterator, emitting each physical record once across overlaps. NULL/empty region means no filter; empty list items and malformed known-contig intervals error. Unknown contigs follow HTSlib's skip policy.

### Typed regions

Typed regions: regions := STRUCT(chrom VARCHAR, start BIGINT, "end" BIGINT)[] selects 0-based half-open intervals (VCF POS p is [p-1, p)) with one native multi-region iterator built directly from header/index contig ids, not region strings. NULL or omitted regions keeps the ordinary behaviour; a non-NULL empty list selects zero records, including COUNT(*), and never falls back to a scan (list() over no rows is NULL, so wrap it as coalesce(list(...), []::STRUCT(chrom VARCHAR, start BIGINT, "end" BIGINT)[])). region and non-NULL regions are mutually exclusive; regions requires a usable index (.tbi/.csi) and is incompatible with scan_mode := 'sequential'. chrom is a literal contig name, so names such as 'HLA-A*01:01' or 'chr1:2' need no braces; a contig absent from the header/index yields no records for that interval, like an unknown contig in region. NULL fields, start < 0 and end <= start error. Intervals are sorted and overlapping or adjacent ones merged per contig without broadening; a physical record is returned once however many intervals overlap it, and physical duplicates stay separate. Records that merely overlap an interval (a deletion spanning it) are returned: filter POS in SQL for exact positions. Caps: 1,000,000 intervals and 128 MiB of native payload, checked before storage; they do not bound the caller's own list aggregation. Prepare the plan on the caller connection with SET VARIABLE p = (SELECT coalesce(list({'chrom': chrom, 'start': pos - 1, 'end': pos}), []::STRUCT(chrom VARCHAR, start BIGINT, "end" BIGINT)[]) FROM (SELECT DISTINCT chrom, pos FROM requests)) and pass regions := getvariable('p'). An unknown variable name is NULL to DuckDB (a full scan), so rduckhts_bcf(regions_var =) and rduckhts_geno(regions_var =) check that the variable exists.

### Scanning

scan_mode='sequential' streams without loading an index. auto streams if no usable index was available at bind; region queries require one. Index-only counts do not validate data contents. Automatic/region plans retain the parsed bind-time index, including remote/explicit index_path: subsequent index removal, replacement or corruption does not change the plan. Supply an initially matching pair and keep data/header contents unchanged; this is not a snapshot or validation of an unrelated index. Reprepare to use a new index.

### Threads

decompression_threads controls per-file HTSlib workers, default 0; it is separate from DuckDB scan parallelism.

### Examples

```sql
SELECT CHROM, POS, REF, ALT FROM read_bcf('vcf_file.bcf') LIMIT 5;
```

## read_geno

Read one row per VCF/BCF record with typed arbitrary-ploidy GT/PS calls, selected FORMAT fields and optional original VCF genotype text.

Signature:

```sql
read_geno(path, region := NULL, index_path := NULL, samples := NULL, non_reference_only := FALSE, scan_mode := 'auto', decompression_threads := 0, decode_error_policy := 'null', format_fields := NULL, raw_gt := FALSE, include_filter := FALSE, regions := NULL)
```

Returns:

```
table
```

### Output

Return record_index UBIGINT, CHROM, POS, ID, REF, ALT VARCHAR[] and calls STRUCT(sample_index UINTEGER, alleles INTEGER[], phase_before BOOLEAN[], phase_set BIGINT)[]. Sample indices are zero-based original-header positions; join read_bcf_samples() from the same unchanged file for names. Absent GT has NULL allele/phase lists; missing allele slots remain NULL entries. Phase flags follow HTSlib, including its pre-VCF-4.4 leading-slot convention. Absent PS is NULL; PS must declare Number=1,Type=Integer.

### Selection

samples uses read_bcf's selectors. non_reference_only removes calls without a called ALT, never records. Zero selected samples yields empty calls on every row. Calls decode only when projected; no normalization, depth or phase inference is performed.

### Ordering

Full scans stream in physical input order regardless of index availability, with one scan worker and separate HTSlib decompression threads. record_index is zero-based and scan-local, not a persistent file locator. Indexed region unions restart at zero in HTSlib iterator order; overlapping regions visit a physical record once while preserving distinct duplicates. Use ORDER BY record_index when order matters. Sequential mode rejects regions and typed regions.

### FORMAT

format_fields adds header-typed members in calls.format, such as AD/DP/GQ; NULL/[] preserves default schema. Header Type/Number determine scalar/list shape. Missing elements retain ordinals, vector-end padding terminates lists, absent fields are NULL and explicit numeric missing lists are [NULL], without inferred A/R/G padding. Missing GT does not hide selected values unless the call is filtered. Unknown, empty, NULL or case-insensitively duplicate selections error. Lookup uses exact header spelling: GT/PS cannot be reselected; declared lowercase gt/ps are distinct tags.

### FILTER

include_filter := TRUE adds the same physical record's FILTER as VARCHAR[]. PASS is [PASS], an unapplied dot is NULL, and named failing filters retain header order. FALSE preserves the default output schema.

### Errors

Header/type/cardinality mismatches follow decode_error_policy. Scalar cardinality excludes vector-end padding, including after sample selection, but counts missing elements. Multiple scalar elements error or make that FORMAT tag NULL for every selected sample under null/warn; list cardinality is not inferred. Physical read/allocation failures always error.

### Raw GT

raw_gt := TRUE appends calls.raw_gt VARCHAR after any format struct, preserving exact original VCF separators, leading phase markers and allele spelling. Absent GT is NULL; literal '.' stays text. FALSE preserves the default schema. BCF rejects TRUE even when calls is unprojected. Text follows original-header sample selection and scan ordinals, including region unions; it does not emulate a downstream parser.

### Typed regions

regions := STRUCT(chrom VARCHAR, start BIGINT, "end" BIGINT)[] has the same contract as read_bcf: 0-based half-open, NULL vs empty, literal contig names, index required, 1,000,000 interval and 128 MiB caps; see read_bcf. Ordinals restart at zero in HTSlib iterator order.

### Exact alleles

Exact alleles: regions plus SQL restore requested (chrom, pos, REF, ALT) alleles without changing the reader. Keep requests unchanged, plan I/O with SELECT DISTINCT chrom, pos, read with non_reference_only := false and format_fields such as ['GP','DS','HS'] so hom-ref and missing calls are returned, expand ALT with its ordinal, and join on equality: WITH alleles AS (SELECT record_index, CHROM, POS, REF, ALT AS full_alt, calls, unnest(ALT) AS match_alt, generate_subscripts(ALT, 1) AS alt_index FROM read_geno('f.bcf', regions := getvariable('p'), non_reference_only := false)) SELECT s.*, a.record_index, a.full_alt, a.alt_index, a.calls FROM requests s LEFT JOIN alleles a ON s.chrom = a.CHROM AND s.pos = a.POS AND s.ref = a.REF AND s.alt = a.match_alt. Use equality on the expanded ALT: list_contains(a.ALT, s.alt) in the join condition plans a nested-loop join. Multi-ALT records keep every ALT in file order with the matched 1-based alt_index, and a request matching two ALT slots or two physical records returns a row each. GP, DS, HS and phase are full-site values: no allele splitting, GP remapping, dosage selection, trimming, strand flip or liftover, so a shifted indel representation is not matched. An ALT-less record disappears under unnest, so diagnosing unmatched requests needs a second step; rduckhts_geno_sites() composes the complete recipe on the caller connection, labelling each request matched, allele_not_at_site, ref_mismatch or absent (position-level diagnostics only, never a fabricated genotype). record_index is scan-local; order with ORDER BY request_id, record_index, alt_index. The reader does not deduplicate requests.

### Examples

```sql
SELECT record_index, CHROM, POS, unnest(calls) AS call FROM read_geno('geno_calls.bcf', non_reference_only := TRUE) ORDER BY record_index;
```

```sql
SET VARIABLE p = [{'chrom': 'chr1', 'start': 99, 'end': 100}, {'chrom': 'HLA-A*01:01', 'start': 9, 'end': 10}]; SELECT record_index, CHROM, POS, REF, ALT FROM read_geno('geno_sites.bcf', regions := getvariable('p'), non_reference_only := FALSE, format_fields := ['GP', 'DS', 'HS']);
```

## read_bcf_samples

Read the typed VCF/BCF sample catalog as sample_index UINTEGER and sample_name VARCHAR without reading records. Indices are zero-based positions in the original header, remain stable under selection and join read_geno calls from the same unchanged file. NULL or '-' selects all; an empty string selects none; comma-separated names include samples; '^' excludes them. Names are validated by HTSlib, selected rows retain header order, and unknown names error.

Signature:

```sql
read_bcf_samples(path, samples := NULL)
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM read_bcf_samples('geno_calls.bcf', samples := '^S1') ORDER BY sample_index;
```

## read_bam

Read SAM/BAM/CRAM alignments with optional typed SAM tags, auxiliary maps and packed sequence, quality or CIGAR output.

Signature:

```sql
read_bam(path, standard_tags := FALSE, auxiliary_tags := FALSE, region := NULL, index_path := NULL, reference := NULL, sequence_encoding := 'string', quality_representation := 'string', cigar_representation := 'string', scan_mode := 'auto', decompression_threads := 2)
```

Returns:

```
table
```

### Representation

sequence_encoding='nt16' returns SEQ UTINYINT[]; quality_representation='phred' returns QUAL UTINYINT[]; cigar_representation='binary' returns packed BAM CIGAR UINTEGER[] instead of SAM text.

### Offsets

FILE_OFFSET is the BGZF virtual position immediately after each compressed BAM record, not its start. SAM (including compressed SAM), CRAM and non-BGZF BAM return NULL. ORDER BY FILE_OFFSET orders one unchanged BGZF BAM; arrival order is not guaranteed and offsets are not comparable across files.

### Scanning

scan_mode='sequential' streams rather than using indexed count/parallel paths and rejects region. NULL/empty region means no filter; empty comma-separated items and malformed known-contig intervals error. Unknown contigs follow HTSlib's skip policy. An indexed full-file scan is split by reference across DuckDB threads and reports the index's row total to the planner as an estimate. Both need the optional per-reference statistics of a BAI or CSI: when a reference has alignments and no statistics, the scan is one sequential stream and no estimate is reported. A CRAM scan is split by reference and reports no estimate. In a region list, '.' is the whole file and '*' is the reads without coordinates. htslib locates both from the same optional statistics; for a BAM whose index lacks them, read_bam reads those two items sequentially, so the result does not depend on the index statistics.

### Threads

decompression_threads controls per-file HTSlib workers, default 2; use 0 to disable. It does not set DuckDB processing parallelism.

### Examples

```sql
SELECT QNAME, FLAG, RNAME, POS FROM read_bam('range.bam') LIMIT 5;
```

## duckhts_bcf_convert_parquet_sql

Build COPY SQL for read_bcf() output with Parquet metadata, VCF header text and selected columns, filters or partitions.

Signature:

```sql
duckhts_bcf_convert_parquet_sql(path, output, columns := []::VARCHAR[], region := NULL, index_path := NULL, tidy_format := FALSE, additional_csq_column_types := NULL, decompression_threads := 0, where_sql := NULL, compression := 'zstd', row_group_size := 100000, partition_by := []::VARCHAR[], include_metadata := TRUE, header_text := NULL, metadata := map([]::VARCHAR[], []::VARCHAR[]), metadata_json_file := NULL, overwrite := FALSE, write_format_version := '1')
```

Returns:

```
VARCHAR
```

### Metadata

Include DuckHTS key/value metadata, preserve or correct header text, and add user metadata with metadata := map(...). metadata_json_file is caller-managed and requires DuckDB's json extension at builder invocation; use maps for offline/CRAN workflows.

### Execution

The function returns SQL without executing it; the R wrapper executes it through DBI.

### Examples

```sql
SELECT duckhts_bcf_convert_parquet_sql('cohort.vcf.gz', 'cohort.parquet', columns := ['CHROM','POS','REF','ALT'], metadata := map(['project'], ['cohort-a']));
```

## duckhts_bam_convert_parquet_sql

Build COPY SQL for read_bam() output with Parquet metadata, SAM header text and selected columns, filters or partitions.

Signature:

```sql
duckhts_bam_convert_parquet_sql(path, output, columns := []::VARCHAR[], region := NULL, index_path := NULL, reference := NULL, standard_tags := FALSE, auxiliary_tags := FALSE, sequence_encoding := NULL, quality_representation := NULL, cigar_representation := NULL, decompression_threads := 2, where_sql := NULL, compression := 'zstd', row_group_size := 100000, partition_by := []::VARCHAR[], include_metadata := TRUE, header_text := NULL, metadata := map([]::VARCHAR[], []::VARCHAR[]), metadata_json_file := NULL, overwrite := FALSE, write_format_version := '1')
```

Returns:

```
VARCHAR
```

### Metadata

Include DuckHTS key/value metadata, preserve or correct header text, and add user metadata with metadata := map(...). metadata_json_file is caller-managed and requires DuckDB's json extension at builder invocation; use maps for offline/CRAN workflows.

### Execution

The function returns SQL without executing it; the R wrapper executes it through DBI.

### Examples

```sql
SELECT duckhts_bam_convert_parquet_sql('sample.bam', 'sample.parquet', columns := ['QNAME','FLAG','RNAME','POS']);
```

## duckhts_gff_convert_parquet_sql

Build COPY SQL for read_gff() output with Parquet metadata, GFF/tabix header text and selected columns, filters or partitions.

Signature:

```sql
duckhts_gff_convert_parquet_sql(path, output, columns := []::VARCHAR[], region := NULL, index_path := NULL, header := NULL, header_names := []::VARCHAR[], auto_detect := NULL, column_types := []::VARCHAR[], attributes_map := FALSE, attributes_list := FALSE, attributes_pairs := FALSE, strict := FALSE, where_sql := NULL, compression := 'zstd', row_group_size := 100000, partition_by := []::VARCHAR[], include_metadata := TRUE, header_text := NULL, metadata := map([]::VARCHAR[], []::VARCHAR[]), metadata_json_file := NULL, overwrite := FALSE, write_format_version := '1')
```

Returns:

```
VARCHAR
```

### Metadata

Include DuckHTS key/value metadata, preserve or correct header text, and add user metadata with metadata := map(...). metadata_json_file is caller-managed and requires DuckDB's json extension at builder invocation; use maps for offline/CRAN workflows.

### Execution

The function returns SQL without executing it; the R wrapper executes it through DBI.

### Examples

```sql
SELECT duckhts_gff_convert_parquet_sql('annotations.gff3.gz', 'annotations/', columns := ['seqname','feature','start','end'], partition_by := ['feature']);
```

## duckhts_tabix_convert_parquet_sql

Build COPY SQL for read_tabix() output with Parquet metadata, header text and selected columns, filters or partitions.

Signature:

```sql
duckhts_tabix_convert_parquet_sql(path, output, columns := []::VARCHAR[], region := NULL, index_path := NULL, header := NULL, header_names := []::VARCHAR[], auto_detect := NULL, column_types := []::VARCHAR[], where_sql := NULL, compression := 'zstd', row_group_size := 100000, partition_by := []::VARCHAR[], include_metadata := TRUE, header_text := NULL, metadata := map([]::VARCHAR[], []::VARCHAR[]), metadata_json_file := NULL, overwrite := FALSE, write_format_version := '1')
```

Returns:

```
VARCHAR
```

### Metadata

Include DuckHTS key/value metadata, preserve or correct header text, and add user metadata with metadata := map(...). metadata_json_file is caller-managed and requires DuckDB's json extension at builder invocation; use maps for offline/CRAN workflows.

### Execution

The function returns SQL without executing it; the R wrapper executes it through DBI.

### Examples

```sql
SELECT duckhts_tabix_convert_parquet_sql('regions.tsv.gz', 'regions.parquet', header_names := ['chrom','pos','value'], auto_detect := TRUE);
```

## read_pileup

Construct a region-scoped BAM pileup with one row per covered position, emitting chrom, 1-based position, depth, observed bases, and Phred+33 qualities after SAM flag and MAPQ filtering. This is a compact htslib pileup view, not samtools mpileup text parity.

Signature:

```sql
read_pileup(path, region := NULL, index_path := NULL, min_mapq := 0, flag_mask := 1796)
```

Returns:

```
table
```

### Examples

```sql
SELECT chrom, pos, depth FROM read_pileup('range.bam', region := 'CHROMOSOME_I:1-200') LIMIT 5;
```

## read_fasta

Read full FASTA records or indexed regions with text or packed sequence output.

Signature:

```sql
read_fasta(path, region := NULL, index_path := NULL, gzi_path := NULL, sequence_encoding := 'string', scan_mode := 'auto')
```

Returns:

```
table
```

### Output

sequence_encoding='nt16' returns SEQUENCE UTINYINT[] using HTSlib nt16 codes instead of VARCHAR. NAME is the literal indexed header name; HTSlib quoting permits comma/colon names. Repeated requests retain one row per interval.

### Scanning

For bgzipped FASTA, gzi_path selects a non-colocated .gzi sidecar. scan_mode='sequential' streams/counts without indexed count paths and rejects region. NULL/empty region means no filter; empty comma-separated items error.

### Examples

```sql
SELECT NAME, length(SEQUENCE) FROM read_fasta('ce.fa');
```

## read_bed

Read BED3-BED12 interval files with canonical typed columns and optional tabix-backed region filtering.

Signature:

```sql
read_bed(path, region := NULL, index_path := NULL, scan_mode := 'auto', error_policy := 'error')
```

Returns:

```
table
```

### scan_mode

scan_mode := 'sequential' forces full-file streaming/counting instead of index-backed count paths and is incompatible with region.

### validation

read_bed is a lenient reader, not a BEDv1 validator (https://samtools.github.io/hts-specs/BEDv1.pdf): it also reads UCSC track files, skipping track and browser lines, and it only rejects data lines with fewer than three tab-delimited fields. Coordinates are not checked against BEDv1: a non-integer chromStart/chromEnd reads as NULL and chromEnd < chromStart is returned as is, so filter those in SQL when strict BED is required.

### error_policy

error_policy := 'error' aborts on data lines with fewer than three tab-delimited fields; 'skip' drops those lines; 'report' emits them with NULL normal columns and adds error (VARCHAR), line_number (BIGINT, 1-based physical line including headers/comments), and raw_line (VARCHAR, without newline). Good rows have NULL error and raw_line. Report requires a full-file scan, not a region query. Non-integer numeric fields become NULL, reversed intervals remain accepted, and blank/header/track/browser lines are ignored in all modes.

### Examples

```sql
SELECT chrom, start, "end", name FROM read_bed('targets.bed') LIMIT 5;
```

## fasta_nuc

Compute bedtools nuc-style nucleotide composition for supplied BED intervals or generated fixed-width bins over a FASTA reference. A failed reference fetch fails the query with the file and zero-based half-open interval; requested intervals are not silently omitted. For bgzipped FASTA, gzi_path may point to an explicit .gzi sidecar when it is not colocated with the FASTA.

Signature:

```sql
fasta_nuc(path, bed_path := NULL, bin_width := NULL, region := NULL, index_path := NULL, gzi_path := NULL, bed_index_path := NULL, include_seq := FALSE)
```

Returns:

```
table
```

### Examples

```sql
SELECT chrom, start, "end", pct_gc FROM fasta_nuc('ce.fa', bin_width := 1000) LIMIT 5;
```

## duckhts_cgranges_create

Create an empty session-scoped cgranges registry entry that can be populated with intervals and finalized for overlap queries.

Signature:

```sql
duckhts_cgranges_create(name)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT duckhts_cgranges_create('targets_idx');
```

## duckhts_cgranges_add

Append an interval to a session-scoped cgranges registry entry before finalization. Labels may be BIGINT-like, DOUBLE, VARCHAR, or BOOLEAN.

Signature:

```sql
duckhts_cgranges_add(name, chrom, start, end[, label])
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT duckhts_cgranges_add('targets_idx', 'chr1', 10, 20, 'exon1');
```

## duckhts_cgranges_index

Finalize a populated cgranges registry entry and build its immutable overlap index for subsequent queries.

Signature:

```sql
duckhts_cgranges_index(name)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT duckhts_cgranges_index('targets_idx');
```

## duckhts_cgranges_destroy

Destroy a session-scoped cgranges registry entry and release its indexed interval storage when it is not in active use.

Signature:

```sql
duckhts_cgranges_destroy(name)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT duckhts_cgranges_destroy('targets_idx');
```

## duckhts_cgranges_from_table

Create, populate and finalize a session-scoped cgranges registry entry from the rows of a table or view, on the caller's connection.

Signature:

```sql
duckhts_cgranges_from_table(name, table_name, chrom_col, start_col, end_col, label_col)
```

Returns:

```
table(indexed BOOLEAN)
```

### Source

table_name is any relation visible to the calling connection, including TEMP tables, views and uncommitted rows; build it from a query with CREATE TEMP VIEW. chrom_col, start_col, end_col and label_col name its columns; label_col may be omitted (five arguments) for unlabeled intervals. The entry is created, filled and finalized in one statement and returns one TRUE row, so duckhts_cgranges_index(...) is not needed afterwards.

### Errors

The name must be new in the session. A failure while filling (NULL chrom, start or end, an unsupported label type, a coordinate outside int32) leaves the partly filled entry in place; remove it with duckhts_cgranges_destroy(...) before retrying.

### Order

Without label_col the label and interval_ordinal follow insertion order. A parallel scan of a large table does not fix that order, so pass label_col when a stable identity is needed.

### Examples

```sql
SELECT * FROM duckhts_cgranges_from_table('targets_idx', 'targets', 'chrom', 'start', 'end', 'name');
```

## duckhts_cgranges_has_overlap

Vectorized scalar predicate for streaming provider rows through a finalized session-scoped cgranges index. Returns TRUE when the query interval overlaps at least one indexed interval, or when mode = 'contain' and it fully contains at least one indexed interval; NULL inputs return NULL.

Signature:

```sql
duckhts_cgranges_has_overlap(name, chrom, start, end[, mode])
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT * FROM read_bed('queries.bed') WHERE duckhts_cgranges_has_overlap('targets_idx', chrom, start, "end");
```

## duckhts_cgranges_count_overlaps

Vectorized scalar overlap counter for streaming provider rows through a finalized session-scoped cgranges index. Returns the number of indexed intervals that overlap the query interval, or with mode = 'contain' the number fully contained by it; NULL inputs return NULL.

Signature:

```sql
duckhts_cgranges_count_overlaps(name, chrom, start, end[, mode])
```

Returns:

```
BIGINT
```

### Examples

```sql
SELECT chrom, start, "end", duckhts_cgranges_count_overlaps('targets_idx', chrom, start, "end") AS n_targets FROM read_bed('queries.bed');
```

## duckhts_cgranges_overlaps_list

Vectorized scalar overlap expander for streaming provider rows through a finalized session-scoped cgranges index. Returns a LIST of hit STRUCTs that can be expanded with UNNEST, preserving provider columns while emitting one row per matching indexed interval. Because scalar return types are fixed, labels are returned as text with label_type describing the original cgranges label kind; NULL inputs return NULL.

Signature:

```sql
duckhts_cgranges_overlaps_list(name, chrom, start, end[, mode])
```

Returns:

```
STRUCT(interval_ordinal BIGINT, label VARCHAR, label_type VARCHAR, interval_chrom VARCHAR, interval_start INTEGER, interval_end INTEGER)[]
```

### Bulk probing

This is the bulk overlap path for any relation of probes. Expand the list in the SELECT list, as in the example: a select-list UNNEST streams the hits, whereas CROSS JOIN UNNEST over the same list plans a lateral join that was several times slower in measurement. Use row_number() OVER () or an existing key to identify probes.

### Examples

```sql
SELECT q.chrom, q.start, q."end", unnest(duckhts_cgranges_overlaps_list('targets_idx', q.chrom, q.start, q."end")) AS hit FROM read_bed('queries.bed') AS q;
```

## duckhts_cgranges_overlaps

Query a finalized session-scoped cgranges registry entry and return one row per overlapping or containing indexed interval, preserving the original label type and interval coordinates.

Signature:

```sql
duckhts_cgranges_overlaps(name, chrom, start, end, mode := 'overlap', query_row_id := NULL)
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM duckhts_cgranges_overlaps('targets_idx', 'chr1', 100, 150);
```

## read_fastq

Read single-end, paired-end, or interleaved FASTQ files with optional legacy quality decoding. By default, FASTQ qualities are interpreted as modern Phred+33 input. Use sequence_encoding := 'nt16' to return SEQUENCE as UTINYINT[] and quality_representation := 'phred' to return QUALITY as UTINYINT[] instead of VARCHAR. input_quality_encoding accepts 'phred33', 'auto', 'phred64', or 'solexa64'. scan_mode := 'sequential' forces raw streaming/counting instead of index-backed count paths.

Signature:

```sql
read_fastq(path, interleaved := FALSE, mate_path := NULL, sequence_encoding := 'string', quality_representation := 'string', input_quality_encoding := 'phred33', scan_mode := 'auto')
```

Returns:

```
table
```

### Examples

```sql
SELECT NAME, MATE FROM read_fastq('r1.fq', mate_path := 'r2.fq') LIMIT 5;
```

## read_bigwig

Read stored BigWig signal intervals as CHROM, START0, END0 and VALUE.

Signature:

```sql
read_bigwig(path, region := NULL, blocks_per_iteration := 64)
```

Returns:

```
table
```

### Coordinates

Output intervals are zero-based half-open. region uses HTSlib one-based inclusive syntax; comma-separated requests merge per contig and emit each stored interval once.

### Execution

Full scans parallelize over nonempty contigs and region scans over merged ranges, with worker-owned handles/iterators. blocks_per_iteration controls libBigWig batching, not DuckDB worker count. Local, native remote and browser wasm reads use HTSlib hFILE transport.

### Examples

```sql
SELECT * FROM read_bigwig('scores.bw', region := 'chr1:100000-101000,chr2:200000-201000');
```

## duckhts_fastq_qc

Aggregate canonical sequence and Phred+33 quality strings directly into exact read/base/Q20/Q30/Q40, nucleotide, quality-sum, and per-cycle sufficient statistics. The nested cycles list supports mean-quality, nucleotide-content, GC, and read-length curves without expanding one SQL row per base. Rows with any NULL input are ignored. Per-cycle state defaults to at most 1,048,576 cycles; pass a constant max_cycles per aggregate group to choose a larger explicit limit, up to 16,777,216.

Signature:

```sql
duckhts_fastq_qc(sequence, quality [, max_cycles])
```

Returns:

```
STRUCT
```

### Examples

```sql
SELECT duckhts_fastq_qc(SEQUENCE, QUALITY) AS qc FROM read_fastq('reads.fastq.gz');
```

```sql
WITH q AS (SELECT duckhts_fastq_qc(SEQUENCE, QUALITY) AS qc FROM read_fastq('reads.fastq.gz')) SELECT cycle.* FROM q, UNNEST(qc.cycles) AS u(cycle);
```

## duckhts_somalier_spacing

Greedily select ranked Somalier candidate positions at a minimum genomic distance.

Signature:

```sql
duckhts_somalier_spacing(positions UBIGINT[], distance UBIGINT)
```

Returns:

```
BOOLEAN[]
```

### Ordering

positions must be ordered by AF rank within one chromosome. The returned mask has one element per input position; the first position wins each conflict, with equal positions conflicting. Positions are positive one-based integers; distance is positive and candidates are limited to one million per call.

### Recipe

Filter a typed read_bcf-style relation in SQL, group candidates by chromosome with list(CAST(POS AS UBIGINT) ORDER BY score, deterministic_tie_key), call this scalar on the ordered list, expand its Boolean mask, apply class-specific caps, and number the final rows by chromosome and position. rduckhts_somalier_find_sites() composes the complete recipe on the caller connection, including caller TEMP tables. No private connection or init-time macro is involved.

### Examples

```sql
SELECT duckhts_somalier_spacing([100, 110, 109]::UBIGINT[], 10);
```

## duckhts_somalier_import_sites

Import an already selected Somalier sites VCF/BCF as one canonical panel and population-frequency relation, autosomal and X/Y sites together.

Signature:

```sql
duckhts_somalier_import_sites(path, assembly_name, max_sites := 1000000)
```

Returns:

```
table
```

### Orientation

Each biallelic SNV is oriented into lexical allele_a/allele_b order. INFO/AF is the source ALT frequency and is transformed with REF/ALT, so population_b_af always describes allele_b. Dense zero-based site_index numbers every autosomal site first, in Somalier v0.3.4's lexical region then position order, and then the X/Y sites in the same order, so the autosomal ordinals are 0..n-1 whatever the X/Y region names sort as.

### Input

INFO/AF must declare Number=A,Type=Float. Retained records must be distinct uppercase A/C/G/T biallelic SNVs with one finite AF in [0,1]. Records on the exact Somalier v0.3.4 X/Y aliases (X, chrX, NC_000023.10, NC_000023.11, Y, chrY, NC_000024.9, NC_000024.10) are kept as sex-chromosome sites and at least one autosomal record is required to normalise their depth. No pseudo-autosomal filtering happens here or at count time, as in Somalier: choose X sites outside the PAR when selecting the sites file (rduckhts_somalier_find_sites does for human builds). Existing FILTER values are retained as provenance because this imports an already selected sites file; it does not apply selection policy.

### Scope

The returned relation directly satisfies the canonical panel and population-frequency contracts and retains source REF, ALT, ALT frequency and FILTER. This is not Somalier find-sites: population AF/AN and QC filtering, interval exclusions, linkage spacing and target-frequency ranking are a separate panel-selection method. max_sites is an explicit input limit.

### Examples

```sql
CREATE TABLE fingerprint_panel AS SELECT * FROM duckhts_somalier_import_sites('sites.vcf.gz', 'GRCh38');
```

## duckhts_somalier_vcf_counts

Extract a complete panel-aligned A/B/other count relation from VCF/BCF FORMAT/AD.

Signature:

```sql
duckhts_somalier_vcf_counts(path, panel_table, samples := NULL, filter_policy := 'pass_or_unapplied')
```

Returns:

```
table
```

### Panel and samples

panel_table is the canonical typed panel relation with one assembly, dense zero-based site_index, region, one-based position and lexical uppercase allele_a/allele_b. Output has exactly one row per selected header sample and panel ordinal, including absent source records. samples uses HTSlib selection syntax and source_sample_index retains the original zero-based header ordinal.

### AD mapping

FORMAT/AD must declare Number=R,Type=Integer. Exact REF/ALT identity selects A and B; remaining declared-allele AD slots are summed as other. Counts are never inferred from GT or DP. Missing records, alleles or AD produce three NULL counts and a named unavailable status; three zeros are measured evidence. Duplicate source records or allele identities, malformed cardinality and negative counts error.

### Scope and FILTER

Every panel site, autosomal or X/Y, is counted in the one pass. Symbolic alleles are unavailable in this scope. filter_policy='pass_or_unapplied' accepts PASS and dot but makes named failures unavailable; 'include_all' uses their AD; 'error' rejects them. X/Y sites ignore filter_policy and always use their AD, as Somalier extraction does. source_method, count_scope, status, FILTER, record/sample ordinals and matched allele slots preserve extraction provenance. Full scans are sequential and retain physical source-record identity; Parquet is optional downstream storage, not an input requirement.

### Examples

```sql
CREATE TABLE allele_counts AS SELECT * FROM duckhts_somalier_vcf_counts('cohort.bcf', 'fingerprint_panel');
```

## duckhts_somalier_bam_counts

Extract complete panel-aligned A/B/other base counts from one indexed BAM or CRAM source.

Signature:

```sql
duckhts_somalier_bam_counts(source_path, panel_table, sample_id, reference_path, panel_parquet := NULL, index_path := NULL, reference_index_path := NULL, min_mapq := 1, min_baseq := 0, require_flags := 0, exclude_flags := 1796, overlap_policy := 'hileup_v0.1.0', decompression_threads := 0, worker_count := 1, max_depth := 100000, max_overlap_qnames := 100000, max_sites := 1000000, max_region_bytes := 67108864, remote_block_bytes := 1048576, remote_cache_bytes := 67108864, reference_cache_bytes := 67108864)
```

Returns:

```
table
```

### Execution

panel_parquet supplies the canonical typed six-column panel; panel_table is retained for signature compatibility and must be NULL. The panel is read from a local Parquet file, materialized and validated once on a private in-memory DuckDB instance opened and closed for the call, then worker_count partitions it into at most 64 DuckDB scan jobs. Each worker-local job owns one alignment handle, index, multi-region iterator, pileup, reference handle and bounded overlap workspace; no mutable htslib or faidx state is shared. DuckDB's connection thread setting limits concurrent jobs. decompression_threads separately controls htslib workers per alignment handle. Parallel output order is unspecified; use ORDER BY site_index when order matters.

### Evidence

Every panel site, autosomal or X/Y, is counted in the one scan. A and B are exact uppercase panel bases; other counts every other observed base. Valid uncovered sites emit measured 0/0/0. Missing reference contigs or positions, reference mismatches, alignment-header contig absence and positions beyond a declared alignment contig emit NULL counts with a named unavailable status. Deletions and reference skips are not observed bases. MAPQ, base quality, required/excluded flags and overlap suppression are explicit per-call settings.

### Limits and transport

The private reader keeps no connection into the caller's database, so a closed database is released; it also cannot see the caller's tables, so write the panel with COPY ... TO 'panel.parquet' first (rduckhts_somalier_bam_counts does this for a panel_table). panel_parquet must be a local ordinary Parquet file: remote URLs are not read. max_sites bounds the shared panel. max_depth, max_overlap_qnames and max_region_bytes bound each active scan job; concurrent workspace can therefore grow with worker_count. remote_block_bytes, remote_cache_bytes and reference_cache_bytes apply per worker-owned handle, and zero disables that cache override. Explicit non-colocated BAM/CRAM and FASTA indexes are honored; CRAM keeps the alignment file's original @SQ lengths and does not create a default FASTA-index sidecar. overlap_policy='hileup_v0.1.0' pins encounter-order mate suppression only, while 'none' counts every retained observation.

### Examples

```sql
CREATE TABLE allele_counts AS SELECT * FROM duckhts_somalier_bam_counts('sample.bam', NULL, 'sample-1', 'reference.fa', panel_parquet := 'fingerprint_panel.parquet');
```

## duckhts_ancestry_proportions

Solve a nearest-positive-definite constrained ancestry projection from aggregated PC products.

Signature:

```sql
duckhts_ancestry_proportions(x_pc_major, y_pc, group_count, sum_to_one := true)
```

Returns:

```
DOUBLE[]
```

### Input

x_pc_major is a flat DOUBLE[] in PC-major, group-minor order for X = P^T F_ref; y_pc is a DOUBLE[] with one value per PC for (P^T f) times the shrinkage correction. Groups are ordered identically within each PC. 1 to 30 groups and 1 to 64 PCs are supported. All values must be finite. sum_to_one=true enforces nonnegative q summing to one; false allows a sum at most one. The scalar returns full-precision group proportions in input order; the R wrapper rounds displayed proportions to seven decimals after correlation gating.

### Repair

The Gram matrix X^T X follows Matrix::nearPD defaults: Dykstra alternating projections, relative eigen cutoff 1e-6, infinity-norm convergence tolerance 1e-7 and at most 100 iterations, then eigenvalue floor 1e-8 times the largest absolute eigenvalue with diagonal rescaling. A nonconverged repair errors. The native solver is independent of bigsnpr and quadprog code.

### SQL recipe

rduckhts_ancestry_proportions normalises long or keyed wide reference products to one row per locus with PC and group columns, checks long-reference completeness, unique input-matched reference loci and finite values at aligned sites, and matches sample_id, chromosome, positive one-based whole-number position and alleles from the caller's input-frequency relation. It reverses input frequency (1 - f) for reversed alleles, audits strand flips, and drops ambiguous, duplicate, missing and mismatched input sites. A bounded aligned Parquet relation holds the retained physical rows. Per-sample sums of loading * reference_frequency form X in PC-major, group-minor order; sums of loading * aligned_frequency times correction form y in PC order. The native scalar solves for q at full precision; cor(F_ref q, f) and per-group cor_each are evaluated before failed proportions are gated to NULL and returned proportions rounded. The scratch Parquet is removed after the call; no private connection or macro is created.

### Missing data

Input genotypes are diploid dosage / 2; NULL genotypes and frequencies are dropped, not imputed. bigsnpr imputes missing individual genotypes before calling snp_ancestry_summary, so parity comparisons must use the same retained sites. The default min_cor is 0.4 as in bigsnpr 1.12.21; callers may select a stricter gate for cohorts or low-depth count frequencies.

### Examples

```sql
SELECT duckhts_ancestry_proportions([1.0, 0.0, 0.0, 1.0], [0.25, 0.75], 2, true);
```

## duckhts_roh_segments

Decode runs of homozygosity from one sample's sorted per-chromosome site lists, reproducing bcftools roh.

Signature:

```sql
duckhts_roh_segments(positions, af, pl, map_pos, map_cm, rec_rate, hw_to_az, az_to_hw)
```

Returns:

```
STRUCT(start BIGINT, "end" BIGINT, n_markers INTEGER, quality DOUBLE)[]
```

### Overloads

duckhts_roh_segments(positions, af, pl, map_pos, map_cm, rec_rate, hw_to_az, az_to_hw) takes phred-scaled likelihoods as INTEGER[][] (three values per site: RR, RA, AA). duckhts_roh_segments(positions, af, dosage, gt_error, map_pos, map_cm, rec_rate, hw_to_az, az_to_hw) takes INTEGER[] genotype dosages 0, 1 or 2 and a phred gt_error, as bcftools roh -G does. duckhts_roh_segments(positions, af, ref_count, alt_count, seq_error, contamination, map_pos, map_cm, rec_rate, hw_to_az, az_to_hw) takes INTEGER[] read counts of the other allele and of the allele af refers to, a per-read seq_error in (0, 0.5) and a contamination fraction in [0, 1); duckhts_roh_counts describes that read model, which is a DuckHTS extension with no bcftools counterpart. All arguments are positional; the duckhts_roh macro supplies the defaults.

### Input

positions is a BIGINT[] of one-based coordinates, sorted ascending and all in 1..2147483647; af is a DOUBLE[] of ALT allele frequencies and pl or dosage a list of the same length, one entry per site. Build them per sample and chromosome with list(... ORDER BY pos). Lists of unequal length, NULL positions, unsorted positions, an af outside [0, 1], a dosage outside 0..2 and a hw_to_az or az_to_hw outside [0, 1] are errors that name the problem. A NULL list returns NULL; an empty list returns an empty list. One list is decoded at a time, so memory follows the longest list.

### Skipped sites

As in bcftools roh, a site is skipped when its position repeats the previous entry, its af is NULL, NaN or exactly 0, its genotype evidence is NULL or has NULL elements, its PL does not have three values or has a negative value, all three PL values are equal, or the likelihoods sum to zero. PL values above 255 count as 255. Skipped sites are not markers and do not break a run.

### Model

A two-state hidden Markov model, autozygous (AZ) and Hardy-Weinberg (HW), with initial probabilities 1/2, emissions P(D|AZ) = P(D|RR)(1-f) + P(D|AA)f and P(D|HW) = P(D|RR)(1-f)^2 + 2P(D|RA)f(1-f) + P(D|AA)f^2, and transitions hw_to_az (default 6.7e-8) and az_to_hw (default 5e-9) per base pair compounded over the physical distance between sites. Segments are the Viterbi AZ runs; quality is the mean forward-backward phred score of the sites in the run. start and end are the first and last marker of the run, inclusive.

### Genetic map

map_pos (BIGINT[], one-based, strictly ascending) and map_cm (DOUBLE[], cumulative centimorgans) are the nodes of a genetic map, the positions and third column of an IMPUTE2 map file. When both are given, the crossover terms between two sites are scaled by the map distance in Morgans, interpolated as bcftools roh -m does (the slope between the nodes that bracket the two sites, zero outside the map), and rec_rate, when given, scales it further. With no map and a positive rec_rate, the scale is the physical distance times rec_rate, as bcftools roh -M does. With neither, only the per-base-pair compounding applies. NULL map_pos and map_cm mean no map; giving one without the other, lists of unequal length or an empty map are errors.

### Parity

The arithmetic follows the pinned bcftools vcfroh.c and HMM.c operation for operation, so positions and marker counts equal bcftools roh and qualities agree to the one decimal bcftools prints; differences could arise only from the floating-point contraction choices of a different compiler or CPU.

### Examples

```sql
SELECT duckhts_roh_segments([100, 200, 300], [0.3, 0.3, 0.3], [[0, 30, 60], [0, 30, 60], [0, 30, 60]], NULL::BIGINT[], NULL::DOUBLE[], NULL::DOUBLE, 6.7e-8, 5e-9);
```

## duckhts_roh

Find runs of homozygosity in a VCF/BCF with the bcftools roh model, using allele frequencies from an INFO tag.

Signature:

```sql
duckhts_roh(path, af_tag, genetic_map, hw_to_az := 6.7e-8, az_to_hw := 5e-9, gt_error := NULL, rec_rate := NULL, samples := NULL, max_sites := 20000000, max_site_bytes := 4294967296)
```

Returns:

```
table(sample VARCHAR, chrom VARCHAR, start BIGINT, "end" BIGINT, length BIGINT, n_markers INTEGER, quality DOUBLE)
```

### Input

path is a VCF/BCF with FORMAT/PL (default) or FORMAT/GT (gt_error set), read with read_bcf(tidy_format := true) on the caller's connection. af_tag names an INFO tag declared Type=Float, Number=A; its first value is the ALT allele frequency. Sites with an absent or missing tag, or a frequency of exactly 0, are skipped as bcftools roh skips them. Records with more than one ALT, no ALT, or a first ALT that is symbolic (<*>, <NON_REF>) are not analysed. samples uses read_bcf's selectors.

### Output

One row per run: sample, chrom, start and end (one-based, inclusive, the first and last marker), length = end - start + 1, n_markers and quality (mean forward-backward phred score), the RG lines of bcftools roh. Rows are not ordered; ORDER BY sample, chrom, start. Chromosomes are decoded separately and a run never spans two.

### Genotype evidence

gt_error NULL uses PL. A numeric gt_error (phred, bcftools -G, for example 30) uses diploid GT calls 0/0, 0/1 and 1/1 and ignores PL; other or missing genotypes are skipped, and FORMAT/PL need not exist. gt_error is an error of the site: the called genotype has likelihood close to 1 and a genotype one allele away has 10^(-gt_error/10), so one call that contradicts a run counts against it by that factor and no more, whatever caused it. A run can therefore continue through an isolated heterozygous call; a larger gt_error makes that rarer.

### Transitions

hw_to_az and az_to_hw are the per-base-pair transition probabilities (bcftools -a and -H). rec_rate is a constant recombination rate per base pair (bcftools -M). The overload with a third positional argument takes genetic_map, a table name, with columns chrom, pos and cm (cumulative centimorgans at one-based positions, for example an IMPUTE2 map), interpolated as bcftools -m does; chromosomes without map rows are skipped, as bcftools skips a chromosome with no map file.

### Memory

DuckDB manages, and can spill, the scan and the joins. The sites of each sample and chromosome are held in native buffers of 16 bytes per site until that sample and chromosome is decoded, so memory follows samples times sites. max_sites caps one sample and chromosome (default 20,000,000; at most 100,000,000). max_site_bytes (default 4 GiB) bounds the buffers held at once: a decode does not grow a buffer when the bytes held by all ROH decodes in the process would pass its own max_site_bytes. Decodes that run at the same time with different values each apply their own, so the process total can pass the smaller value while the decode with the larger value allocates. Exceeding either limit is an error: decode fewer samples per query with samples, or raise the limit. Decoding one sample and chromosome also needs 38 bytes per site of workspace in each thread. Model parameters and limits are checked even when the input has no record.

### Limits

Duplicate positions on a chromosome keep the first record to arrive, as bcftools roh keeps the first record in the file. That is the file order when one scan thread reads both records; deduplicate upstream when it matters. With a frequency relation (duckhts_roh_af_table, duckhts_roh_ancestry), only the records that match a frequency row take part. bcftools' default AC/AN frequencies, --AF-dflt, --estimate-AF, --include/--exclude, --skip-indels, --ignore-homref, --buffer-size and --viterbi-training are not offered.

### Contig names

The chrom of genetic_map is compared byte for byte with the contig names of the sites. A chromosome that the map spells differently (1 and chr1) is skipped like a chromosome without map rows, without an error. Rename the map's contigs in a view first.

### Examples

```sql
SELECT * FROM duckhts_roh('cohort.vcf.gz', 'AF') ORDER BY sample, chrom, start;
```

```sql
SELECT * FROM duckhts_roh('cohort.vcf.gz', 'AF', 'genetic_map', gt_error := 30) ORDER BY sample, chrom, start;
```

## duckhts_roh_af_table

Find runs of homozygosity in a VCF/BCF with the bcftools roh model, using allele frequencies from a caller relation.

Signature:

```sql
duckhts_roh_af_table(path, af_table, genetic_map, hw_to_az := 6.7e-8, az_to_hw := 5e-9, gt_error := NULL, rec_rate := NULL, samples := NULL, max_sites := 20000000, max_site_bytes := 4294967296)
```

Returns:

```
table(sample VARCHAR, chrom VARCHAR, start BIGINT, "end" BIGINT, length BIGINT, n_markers INTEGER, quality DOUBLE)
```

### Frequencies

af_table is a table or view name with columns chrom, pos (one-based), ref, alt and af, like a bcftools --AF-file. A record takes the frequency of the row matching its chrom, pos, REF and its ALT alleles joined by commas; a record without a matching row is not used and does not claim a repeated position, and a record whose row has a NULL af is skipped. bcftools roh --AF-file differs at a repeated position: it takes the first record in file order and skips the position when that record's alleles differ from the file's. Give one row per key; af must be in [0, 1], and exactly 0 skips the site. This is the entry point for per-sample frequencies, for example ancestry-tuned ones, by running it once per relation.

### Otherwise

Everything else is as duckhts_roh: the same output, genotype evidence (PL, or GT with gt_error), transition parameters, optional genetic_map overload, memory limits and limits.

### Contig names

The chrom of af_table is compared with the VCF's CHROM byte for byte. When the two spell contigs differently (1 and chr1), no record finds a frequency and the function returns no rows, without an error. Rename the contigs of af_table to the VCF's spelling in a view first; duckhts_contig_key() gives a common key.

### Examples

```sql
SELECT * FROM duckhts_roh_af_table('cohort.vcf.gz', 'site_frequencies') ORDER BY sample, chrom, start;
```

## duckhts_roh_ancestry

Find runs of homozygosity using each sample's ancestry-weighted allele frequencies.

Signature:

```sql
duckhts_roh_ancestry(path, reference_table, proportions_table, genetic_map, af_clamp := 1e-3, hw_to_az := 6.7e-8, az_to_hw := 5e-9, gt_error := NULL, rec_rate := NULL, samples := NULL, max_sites := 20000000, max_site_bytes := 4294967296)
```

Returns:

```
table(sample VARCHAR, chrom VARCHAR, start BIGINT, "end" BIGINT, length BIGINT, n_markers INTEGER, quality DOUBLE)
```

### Inputs

reference_table is a long relation with chromosome, position (one-based), allele_a, allele_b, group_id and frequency (frequency of allele_b). proportions_table has sample_id, group_id and proportion. Both relations must use exactly the same group IDs; every called VCF sample needs proportions. Contig names are matched with duckhts_contig_key(), resolved once per distinct name. It removes one leading chr case-insensitively, writes M/MT as MT and uppercases X and Y, so chr1 and 1 match, as do chrX and X. Accessions, patches and numeric sex chromosomes are not mapped. Each site must have one allele orientation across its group rows, and exactly one row per group.

### Frequency model

For each called site and sample, AF is the sum over groups of proportion times group frequency. The reference is aggregated once per site and proportions once per sample, then joined to called sites. REF=allele_a and ALT=allele_b uses the weighted frequency; the reversed orientation uses 1-AF. Reference alleles must be on the forward strand of the VCF's assembly, as in a FASTA-anchored panel; no strand flip is attempted, so palindromic (A/T, C/G) sites are oriented by REF like any other, and a record whose alleles match neither order is not used and does not claim a repeated position. af_clamp in [0, 0.5) clamps AF to [af_clamp, 1-af_clamp]; zero disables clamping. The default avoids interpreting a population frequency of zero as impossible.

### Otherwise

The output, genotype evidence, transitions, optional genetic_map overload, memory limits and limits are as duckhts_roh_af_table.

### Examples

```sql
SELECT * FROM duckhts_roh_ancestry('cohort.vcf.gz', 'reference_long', 'sample_proportions', gt_error := 30) ORDER BY sample, chrom, start;
```

## duckhts_roh_counts

Find runs of homozygosity from allele read counts, with sequencing error and optional contamination.

Signature:

```sql
duckhts_roh_counts(counts_table, genetic_map, seq_error := 1e-3, contamination := 0.0, hw_to_az := 6.7e-8, az_to_hw := 5e-9, rec_rate := NULL, max_sites := 20000000, max_site_bytes := 4294967296)
```

Returns:

```
table(sample VARCHAR, chrom VARCHAR, start BIGINT, "end" BIGINT, length BIGINT, n_markers INTEGER, quality DOUBLE)
```

### Input

counts_table is a table or view with sample_id, chrom, pos (one-based), ref_count, alt_count and af, one row per sample and site. alt_count counts reads showing the allele whose population frequency is af; ref_count counts reads showing the other allele. Any frequency source composes in SQL first: an INFO tag, a frequency table or ancestry-weighted frequencies. Somalier site counts map as ref_count = a, alt_count = b and af = frequency of allele_b.

### Read model

Each read shows the counted allele with probability (1 - contamination) * q + contamination * c, where q is seq_error, 1/2 or 1 - seq_error for zero, one or two copies, and c = af * (1 - seq_error) + (1 - af) * seq_error is the chance that a read from a contaminating individual of the same population shows it. Reads are independent and the binomial coefficient is omitted. The resulting genotype likelihoods enter the bcftools roh model in place of PL. This emission is a DuckHTS extension; bcftools roh has no read-count mode.

### Relation to genotype evidence

seq_error is an error of one read, and the reads of a site multiply. A site with balanced reads excludes both homozygous genotypes whatever seq_error is, and the read model has no term for a site whose reads do not reflect the sample's genotype. The read-count decode therefore corresponds to the GT decode with no tolerated genotype error: it ends a run at a heterozygous site that duckhts_roh with gt_error := 30 would pass through, and reports fewer bases in runs. benchmarks/benchmark_roh_counts_validation.md measures this on three 1000 Genomes samples.

### Skipped sites

Sites with no reads, a NULL count, or an af that is NULL, NaN or 0 are skipped; a repeated position keeps one row, as in the other overloads. The relation need not be sorted. Negative counts, fractional counts, seq_error outside (0, 0.5) and contamination outside [0, 1) are errors.

### Limitations

Biallelic sites only; one error rate for all reads, with no base- or mapping-quality, duplicate or strand modelling. A heterozygote is assumed to show each allele in half its reads, but reference-biased mapping makes the other allele rarer, more so where flanking heterozygosity is high; such heterozygotes can look homozygous and lengthen runs in divergent regions. The contamination fraction is supplied, not estimated, and a contaminant from a different population than af describes is approximated by af.

### Otherwise

Output, transitions, rec_rate, the optional genetic_map overload and the memory limits are as duckhts_roh, with 24 bytes per site.

### Contig names

The chrom of genetic_map is compared byte for byte with the contig names of the sites. A chromosome that the map spells differently (1 and chr1) is skipped like a chromosome without map rows, without an error. Rename the map's contigs in a view first.

### Examples

```sql
SELECT * FROM duckhts_roh_counts('site_counts', contamination := 0.02) ORDER BY sample, chrom, start;
```

## duckhts_somalier_panel_sha256

Derive a stable SHA-256 identity for an ordered biallelic sample-fingerprinting panel.

Signature:

```sql
duckhts_somalier_panel_sha256(panel_table)
```

Returns:

```
VARCHAR
```

### Panel contract

panel_table has assembly, dense zero-based site_index, region, positive one-based position, allele_a and allele_b. Assembly and region are each limited to 1,024 bytes before hashing. One nonempty assembly, unique physical region/position and uppercase single-base A/C/G/T alleles in lexical A < B order are required. The pinned Somalier v0.3.4 X/Y aliases are sex-chromosome sites; every other region is autosomal, since other aliases cannot be classified from a region string. Every autosomal site_index must precede every X/Y site_index. The digest commits to this domain (panel version 3: biallelic SNVs, autosomal then sex-chromosome sites) and every ordered site, independent of physical row order. Panels digested before this version have different identities.

### Examples

```sql
SELECT duckhts_somalier_panel_sha256('fingerprint_panel');
```

## duckhts_somalier_frequency_sha256

Derive a stable identity for panel-aligned population-B allele frequencies.

Signature:

```sql
duckhts_somalier_frequency_sha256(frequency_table, panel_table)
```

Returns:

```
VARCHAR
```

### Frequency contract

frequency_table must cover every autosomal panel site exactly once (X/Y rows are ignored) with the same assembly, ordinal, coordinate and A/B orientation plus a finite population_b_af in [0,1]. The digest commits to the panel identity and every ordered frequency value.

### Examples

```sql
SELECT duckhts_somalier_frequency_sha256('population_frequencies', 'fingerprint_panel');
```

## duckhts_somalier_classify

Classify one measured A/B/other count tuple for Somalier-derived autosomal relatedness.

Signature:

```sql
duckhts_somalier_classify(a, b, other, min_depth, min_het_balance, hom_balance_cutoff)
```

Returns:

```
STRUCT(genotype TINYINT, middling BOOLEAN, unavailable BOOLEAN)
```

### Evidence

a, b and other are all measured or all NULL. Three zeros are measured zero-depth evidence; three NULL values are unavailable. Genotype is -1 unknown, 0 homozygous A, 1 heterozygous or 2 homozygous B.

### Scope

Alleles must already use the panel's ordered A/B orientation. Relatedness classification preserves the pinned Somalier v0.3.4 10% other-read filter; contamination uses separate stricter eligibility.

### Examples

```sql
SELECT duckhts_somalier_classify(20, 20, 0, 7, 0.3, 0.01);
```

## duckhts_somalier_prepare_sketches

Build one panel-verified packed relatedness sketch per sample from typed count evidence.

Signature:

```sql
duckhts_somalier_prepare_sketches(evidence_table, panel_table, min_depth, min_het_balance, hom_balance_cutoff, max_sites := 1000000)
```

Returns:

```
table(sketch STRUCT)
```

### Evidence

evidence_table contains sample_id, the six panel identity columns, and nullable a, b and other counts. Each sample must contain every panel ordinal exactly once. Panel and evidence rows on the X/Y aliases are ignored: this is an autosomal statistic, and results equal those on the autosomal-only panel. The reported panel identity is the whole panel's. Count channels already follow the panel's canonical lexical A/B order; changed coordinates or orientation error. Counts are never inferred from GT or DP. All-NULL tuples are unavailable, all-zero tuples are measured zero depth, and partial NULL tuples error.

### Persistence

The returned struct contains identities, classification settings, counters, a raw-count receipt, content integrity fields and three UBIGINT[] masks. The receipt binds every ordinal, availability state and A/B/other tuple even when changed counts retain the same genotype. It is an ordinary typed value suitable for Parquet, not a serialized native object. max_sites bounds each prepared sketch; sample_id and assembly are each limited to 1,024 bytes.

### Examples

```sql
CREATE TABLE sample_sketches AS SELECT * FROM duckhts_somalier_prepare_sketches('allele_counts', 'fingerprint_panel', 7, 0.3, 0.01);
```

## duckhts_somalier_verify_sketches

Verify persisted relatedness sketches against their retained raw count evidence.

Signature:

```sql
duckhts_somalier_verify_sketches(evidence_table, panel_table, sketches_table, max_sites := 1000000)
```

Returns:

```
BOOLEAN
```

### Integrity

Checks persisted classification settings before using them as rebuild parameters. With valid settings, it rebuilds every sample through the panel and count-validation path and compares the complete typed sketch. It returns false for invalid retained settings, changed A/B/other counts, changed availability, altered masks or receipts, and missing, extra or duplicate sample sketches. Invalid panel/evidence geometry and incomplete or duplicate site ordinals error.

### Scope

evidence_table must contain exactly the samples represented by sketches_table. max_sites is a positive panel-site limit at most 100,000,000. Receipts detect accidental divergence between persisted evidence and sketches; they are not an authenticity mechanism.

### Examples

```sql
SELECT duckhts_somalier_verify_sketches('allele_counts', 'fingerprint_panel', 'sample_sketches');
```

## duckhts_somalier_relatedness

Compute fused Somalier-derived relatedness and concordance statistics for two prepared sketches.

Signature:

```sql
duckhts_somalier_relatedness(sketch_a, sketch_b, max_sites)
```

Returns:

```
STRUCT
```

### Results

The struct retains sample and panel identities, method version, status, jointly-called count, IBS0/IBS2, shared heterozygotes, heterozygote and homozygote denominators, middling/unavailable counters, relatedness and named concordance values. relatedness is 2(shared_hets - 2 IBS0)/max(1,het_ab). inferred_hom_concordance is matching_hom_count/max(1,min(callable_hom_count_a,callable_hom_count_b)). raw_hom_b_concordance is (shared_hom_b - 2 IBS0)/max(1,min(hom_b_count_a,hom_b_count_b)). adjusted_concordance applies pinned v0.3.4's middling and low-homozygote penalties and upper-range transform. Floating results are NULL when no site is jointly callable.

### Execution

Both sketches must have identical panel digests, site counts and classification settings. The contents of both sketches are checked for every call, so comparing every pair of a cohort repeats that check for each pair; use duckhts_somalier_relatedness_all_pairs for all pairs. Mask words are borrowed directly; comparison allocates no pair-sized workspace. SQL chooses the requested pair relation and output ordering.

### Examples

```sql
SELECT unnest(duckhts_somalier_relatedness(a.sketch, b.sketch, 1000000)) FROM selected_pairs p JOIN sample_sketches a ON a.sketch.sample_id = p.sample_a JOIN sample_sketches b ON b.sketch.sample_id = p.sample_b;
```

## duckhts_somalier_relatedness_all_pairs

Compute Somalier-derived relatedness and concordance statistics for every unordered sample pair of one sketch relation.

Signature:

```sql
duckhts_somalier_relatedness_all_pairs(sketches_table, max_sites := 1000000)
```

Returns:

```
table
```

### Input

sketches_table names a relation with one non-NULL sketch column per sample, as duckhts_somalier_prepare_sketches returns. Two sketches of one sample, a NULL sketch, a sketch with more sites than max_sites, or a sketch whose masks do not match its stored content digest is an error. All sketches must share assembly, panel digest, site count and classification settings; a mixed relation is an error. These checks run once over the relation, so they hold for a query that reads no pair column, for example count(*).

### Results

One row per unordered pair with sample_a < sample_b, with the fields of duckhts_somalier_relatedness as columns. The rows equal those of duckhts_somalier_relatedness over the same pairs. No row order is guaranteed.

### Execution

Each sketch is content-checked once and the pairs are compared from the checked sketches, so the content checks grow with the number of samples and the comparisons with pairs times sites. The output has N(N-1)/2 rows for N samples and is not capped. The pair comparison runs in one thread.

### Examples

```sql
SELECT * FROM duckhts_somalier_relatedness_all_pairs('sample_sketches');
```

## duckhts_somalier_verify_relatedness

Verify a typed relatedness result against its two sealed sketches.

Signature:

```sql
duckhts_somalier_verify_relatedness(pair_result, sketch_a, sketch_b, max_sites)
```

Returns:

```
BOOLEAN
```

### Integrity

Returns false when a sketch digest, sample or panel identity, method/status, site denominator, integer statistic, floating statistic or status-dependent NULL value is inconsistent. The check recomputes the pair metrics from borrowed mask words without a per-pair allocation. It detects accidental corruption, not malicious rewriting of both sketches and their integrity fields.

### Limit

max_sites is a positive panel-site limit at most 100,000,000. Invalid or NULL typed inputs return false.

### Examples

```sql
SELECT duckhts_somalier_verify_relatedness(r.pair, a.sketch, b.sketch, 1000000) FROM pair_results r JOIN sample_sketches a ON a.sketch.sample_id = r.pair.sample_a JOIN sample_sketches b ON b.sketch.sample_id = r.pair.sample_b;
```

## duckhts_somalier_charr

Estimate per-sample contamination with a bounded Somalier-derived CHARR reduction.

Signature:

```sql
duckhts_somalier_charr(evidence_table, panel_table, frequency_table, min_depth := 15, max_depth := 1000000, hom_minor_rate := 0.12, hom_tail_alpha := 0.002, max_threshold_work := 16000000, max_sites := 1000000)
```

Returns:

```
table(contamination STRUCT)
```

### Inputs

Count evidence and population_b_af must match the panel's full ordered site identity. Panel and evidence rows on the X/Y aliases are ignored: this is an autosomal statistic, and results equal those on the autosomal-only panel. The reported panel identity is the whole panel's. CHARR uses measured counts, its own homozygous-like binomial eligibility and the pinned 4% other-read filter; relatedness genotype masks are insufficient.

### Results

The struct retains sample, panel and frequency identities, method and numerical status, observed/unavailable/usable and homozygous-site denominators, estimate and every filter/limit. No usable evidence has status no_evidence and a NULL estimate, distinct from measured zero contamination.

### Limits

A+B depth must not exceed 1,000,000. Distinct measured depths are certified once per call; max_threshold_work bounds their cumulative exact recurrence and continued-fraction steps and is at most 100,000,000. Exhaustion errors without publishing a partial result. sample_id, assembly and panel region are each limited to 1,024 bytes. Input order and parallel aggregate reduction order do not change the estimate.

### Compatibility

Eligibility follows Somalier 0.3.4 CHARR except that DuckHTS computes high-depth binomial tails without upstream's numerical underflow. Very deep sites can therefore have different eligibility; a pinned upstream counterexample is retained in the conformance tests.

### Examples

```sql
SELECT unnest(contamination) FROM duckhts_somalier_charr('allele_counts', 'fingerprint_panel', 'population_frequencies');
```

## duckhts_somalier_matched_contamination

Estimate directional contamination for explicitly selected receiver/anchor sample pairs.

Signature:

```sql
duckhts_somalier_matched_contamination(evidence_table, panel_table, frequency_table, pairs_table, min_depth := 15, max_depth := 1000000, hom_minor_rate := 0.05, hom_tail_alpha := 0.001, error_rate := 0.002, min_probability := 1e-10, min_prior_frequency := 1e-6, alpha_min := 0, alpha_max := 1, grid_step := 0.01, refine_tolerance := 1e-10, max_evaluations := 4096, max_threshold_work := 16000000, max_sites := 1000000)
```

Returns:

```
table(contamination STRUCT)
```

### Direction

pairs_table contains distinct receiver_id and anchor_id rows. Panel and evidence rows on the X/Y aliases are ignored: this is an autosomal statistic, and results equal those on the autosomal-only panel. The reported panel identity is the whole panel's. The anchor supplies the receiver's expected uncontaminated homozygous genotype; it is not assumed to identify the contaminating donor. Reversing a pair is a different fit.

### Results

The struct retains ordered sample, panel and frequency identities, method/status, observed/unavailable/usable denominators, alpha, evaluation count and every filter/search limit. No usable evidence returns NULL alpha. relative_log_likelihood omits alpha-independent binomial coefficients and is comparable only across alpha values for the same observations.

### Execution

The query certifies each distinct measured depth once, then prepares one bounded site profile per distinct selected sample and one panel-aligned frequency profile before joining requested ordered pairs. The pair scalar borrows those DuckDB-owned lists and allocates no per-pair workspace. Panel and evidence cardinalities are checked before profile-list construction; max_sites is a per-call panel limit.

### Limits

Panel assembly and region are each limited to 1,024 bytes. Persisted profiles separately limit sample_id and assembly to 1,024 bytes. A+B depth must not exceed max_depth. max_threshold_work bounds cumulative exact binomial certification steps and is at most 100,000,000; exhaustion errors without publishing profiles. max_sites is at most 100,000,000, and max_evaluations must fit the declared search workspace.

### Numerical difference

DuckHTS searches a full 0.01 grid and refines a feasible local optimum rather than reproducing Somalier v0.3.4's fixed coarse/high-refinement sequence. The retained two-site witness fits alpha about 0.39759 versus upstream 0.440983. Both likelihoods and settings are retained in the differential test; results are not guaranteed bitwise identical to the CLI.

### Examples

```sql
SELECT unnest(contamination) FROM duckhts_somalier_matched_contamination('allele_counts', 'fingerprint_panel', 'population_frequencies', 'receiver_anchor_pairs');
```

## duckhts_somalier_sex

Report X/Y dosage evidence and a review-only XX, XY or ambiguous call per sample from panel counts.

Signature:

```sql
duckhts_somalier_sex(counts_table, panel_table, min_depth := 7, min_sex_site_depth := 7, min_het_balance := 0.3, hom_balance_cutoff := 0.01, min_usable_x_sites := 11, xy_max_het_ratio := 0.05, xx_min_het_ratio := 0.4, y_signal_min := 0.4, y_gate := 'sample')
```

Returns:

```
table
```

### Inputs

counts_table is the sample-by-panel relation from duckhts_somalier_vcf_counts() or duckhts_somalier_bam_counts() over a panel that includes X/Y sites; each sample needs one row for every panel site with the panel's geometry and A/B orientation. Coverage is checked per sample from the row count plus the sum and XOR of a 64-bit hash of site_index against the panel's, because a per-sample DISTINCT set would hold samples x sites entries: a missing or repeated ordinal changes the count or the hashes and is rejected, while a different multiset of ordinals with the same count passes only if both hashes collide. Autosomal sites give the depth normaliser, X and Y sites the metrics. Without a usable X site the panel is still accepted and the status says so.

### Metrics

Follows Somalier v0.3.4 relate.nim. autosomal_depth_mean is the mean A+B depth of autosomal sites with other <= 10% of all reads and A+B >= min_depth. An X or Y site is usable when other <= 4%, A+B >= min_sex_site_depth and its B balance is homozygous (below hom_balance_cutoff, above one minus it, or A or B zero) or heterozygous (between min_het_balance and one minus it); x_hom_ref, x_het and x_hom_alt count usable X sites by genotype. x_depth_ratio and y_depth_ratio are 2 * mean A+B depth over usable sites / autosomal_depth_mean, so a diploid X is about 2 and a single X about 1. x_het_hom_alt_ratio is x_het / x_hom_alt (Infinity when there is no hom-alt site and some het site, NULL when both are zero), the statistic Somalier thresholds. Ratios are NULL without usable sites or autosomal depth; the integer counts stay exact zeros. Floating means can differ from Somalier's running mean in the last bits.

### Call

inferred_sex is XY when there are at least min_usable_x_sites usable X sites (Somalier: more than 10) and x_het_hom_alt_ratio < xy_max_het_ratio, XX when it is > xx_min_het_ratio, otherwise ambiguous. y_signal is present when y_depth_ratio > y_signal_min and absent otherwise, NULL without usable Y sites. An XX dosage with a present Y signal is reported as ambiguous with the review flag apparent_y_with_xx_dosage, and an XY dosage with an absent Y signal keeps XY with apparent_y_loss, as Somalier notes. status is ok, no_usable_x_sites, insufficient_x_sites or no_autosomal_depth.

### Y gate

Somalier applies its Y check only when the cohort has Y depth (more than five samples, or more than 10%, with usable Y sites), so a sample's call depends on who else is in the batch. The default y_gate := 'sample' applies the check whenever the sample itself has a Y signal; y_gate := 'cohort' reproduces Somalier's gate, and cohort_has_y reports it. Somalier also leaves the sex unset for samples that fail its autosomal quality gate; DuckHTS reports the X dosage regardless and reports no_autosomal_depth when the normaliser is missing.

### Interpretation

The call is evidence of X/Y dosage for review, not a diagnosis or a legal determination of sex; aneuploidy, mosaicism, poor coverage and sample swaps all change it. Results also carry the whole-panel digest, panel denominators and every threshold used.

### Memory

See benchmarks/benchmark_somalier_sex.md for measured scaling. One pass over the counts relation: state is one aggregate group of a few numbers per sample plus the panel; no product of samples and sites is retained beyond the counts relation itself.

### Examples

```sql
SELECT sample_id, inferred_sex, status, x_depth_ratio, y_depth_ratio FROM duckhts_somalier_sex('allele_counts', 'fingerprint_panel');
```

## detect_quality_encoding

Inspect a FASTQ file's observed quality ASCII range and report compatible legacy encodings with a heuristic guessed encoding.

Signature:

```sql
detect_quality_encoding(path, max_records := 10000)
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM detect_quality_encoding('reads.fq.gz');
```

## read_gff

Read GFF annotations with optional parsed attributes, strict GFF3 validation and indexed region selection.

Signature:

```sql
read_gff(path, header_names := NULL, header := FALSE, column_types := NULL, auto_detect := FALSE, attributes_map := FALSE, attributes_list := FALSE, attributes_pairs := FALSE, attributes := []::VARCHAR[], strict := FALSE, region := NULL, index_path := NULL, scan_mode := 'auto')
```

Returns:

```
table
```

### Attributes

attributes := ['Parent', 'ID'] appends VARCHAR columns named for the requested keys. Each value equals attributes_map[key], including NULL for absent keys; repeated keys use the first value and GFF3 percent encoding is retained. Key lookup is case-sensitive, as GFF3 tags are: attributes := ['id'] on a file that uses ID returns a column of NULLs. Keys must be nonempty and distinct from fixed and optional attribute column names; because column names ignore case, two keys that differ only in case are rejected. Only projected keys are parsed during scanning.

### Scanning

Comma-separated indexed regions emit each row once across overlaps. scan_mode='sequential' streams/counts instead of using indexed count paths and rejects region. NULL/empty region means no filter; empty comma-separated items and malformed known-contig intervals error. Unknown contigs follow HTSlib's skip policy.

### Examples

```sql
SELECT seqname, feature, start, "end" FROM read_gff('gff_file.gff.gz') LIMIT 5;
```

## read_gtf

Read GTF annotations with optional parsed attributes and indexed region selection.

Signature:

```sql
read_gtf(path, header_names := NULL, header := FALSE, column_types := NULL, auto_detect := FALSE, attributes_map := FALSE, attributes_list := FALSE, attributes_pairs := FALSE, attributes := []::VARCHAR[], region := NULL, index_path := NULL, scan_mode := 'auto')
```

Returns:

```
table
```

### Attributes

attributes := ['gene_id', 'transcript_id'] appends VARCHAR columns named for the requested keys. Each value equals attributes_map[key], including NULL for absent keys; repeated keys use the first value. Key lookup is case-sensitive: attributes := ['Gene_ID'] on a file that uses gene_id returns a column of NULLs. Keys must be nonempty and distinct from fixed and optional attribute column names; because column names ignore case, two keys that differ only in case are rejected. Only projected keys are parsed during scanning.

### Scanning

Comma-separated indexed regions emit each row once across overlaps. scan_mode='sequential' streams/counts instead of using indexed count paths and rejects region. NULL/empty region means no filter; empty comma-separated items and malformed known-contig intervals error. Unknown contigs follow HTSlib's skip policy.

### Examples

```sql
SELECT seqname, feature, start, "end" FROM read_gtf('annotations.gtf.gz') LIMIT 5;
```

## read_genbank

Read GenBank flat-file features in read_gff's column shape, with optional parsed qualifier MAP.

Signature:

```sql
read_genbank(path, attributes_map := FALSE, attributes := []::VARCHAR[])
```

Returns:

```
table
```

### Mapping

seqname is VERSION, else ACCESSION, else the LOCUS name; a segment on a remote accession (ACC.1:5..40) is reported under that accession. source is 'GenBank'. join()/order() give one row per segment in biological order, complement(...) sets strand '-', and the GFF3 phase of each CDS segment is carried from /codon_start across segments (absent means 1). A span whose end precedes its start on a circular record wraps the origin into two segments, and a between site n^m is the zero-length site at n. The record-level source feature is dropped, and /translation is omitted as redundant with ORIGIN.

### Attributes

Synthesized GFF3 keys ID, Name and Parent accompany the original qualifiers. ID is gene-<locus_tag> for genes and <key>-<n> otherwise, with n the feature's 0-based position in the file. Parent links a feature to the gene sharing its /locus_tag (or /gene) anywhere in the record. Repeated qualifiers become one key with comma-joined values, valueless qualifiers read 'true', and ; = & , % are percent-encoded. attributes_map := TRUE adds the same pairs as a MAP.

### Named attributes

attributes := ['gene', 'product'] appends VARCHAR columns named for the requested keys, after attributes_map when that is requested. Each value equals attributes_map[key] byte for byte: repeated qualifiers are comma-joined, a valueless qualifier reads 'true', values stay percent-encoded, the synthesized ID, Name and Parent are addressable, and an absent key is NULL. Key lookup is case-sensitive, like qualifier keys themselves: attributes := ['ec_number'] on features that carry /EC_number returns a column of NULLs. Keys must be nonempty and distinct from the fixed columns and attributes_map when it is requested; because column names ignore case, two keys that differ only in case are rejected. Only projected keys are computed during scanning.

### Errors

Records stream one at a time, so memory follows the largest record. A record without a terminating //, a FEATURES table with no sequence section, a malformed location, an unsupported location form (one-of, gap, bond, nested join/order, a.b), an origin-spanning span on a linear record, or a /codon_start outside 1..3 is an error naming the feature and line. Table layout and these rules follow BioPython's GenBank scanner, and test/scripts/genbank_oracle_test.py diffs the reader against it.

### Examples

```sql
SELECT seqname, feature, start, "end" FROM read_genbank('phix174.gb') LIMIT 5;
```

```sql
SELECT feature, locus_tag, product FROM read_genbank('phix174.gb', attributes := ['locus_tag', 'product']) WHERE feature = 'CDS';
```

## genbank_to_fasta

Write the ORIGIN sequence of each GenBank record as FASTA and return success, output_path and records_written.

Signature:

```sql
genbank_to_fasta(path, output_path := NULL, line_width := 70, overwrite := FALSE)
```

Returns:

```
table
```

### Naming

Records are written under the same name read_genbank reports as seqname, so feature coordinates land on the contig of that name; the DEFINITION follows without its trailing period, as in NCBI's FASTA export. A segment on a remote accession is not written.

### Output

output_path defaults to path with .fa appended. The FASTA is written to a temporary file beside output_path and renamed into place only after the input has been read to a clean end, so a failure never leaves a partial output and an existing file is never lost; with overwrite := FALSE an existing output is an error. Records without an ORIGIN block are skipped; zero written records is an error. line_width must be at least 1.

### Examples

```sql
SELECT * FROM genbank_to_fasta('phix174.gb', output_path := 'phix174.fa');
```

## read_tabix

Read tabix-indexed text with optional header handling, inferred types and region selection.

Signature:

```sql
read_tabix(path, header_names := NULL, header := FALSE, column_types := NULL, auto_detect := FALSE, region := NULL, index_path := NULL, scan_mode := 'auto')
```

Returns:

```
table
```

### Scanning

Comma-separated indexed regions emit each row once across overlaps. scan_mode='sequential' streams/counts instead of using indexed count paths and rejects region. NULL/empty region means no filter; empty comma-separated items and malformed known-contig intervals error. Unknown contigs follow HTSlib's skip policy.

### Examples

```sql
SELECT * FROM read_tabix('meta_tabix.tsv.gz') LIMIT 5;
```

## fasta_index

Build a FASTA index (.fai) and return a single row with columns success (BOOLEAN) and index_path (VARCHAR).

Signature:

```sql
fasta_index(path, index_path := NULL)
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM fasta_index('ce.fa');
```

## bgzip

Compress a plain file to BGZF and return the created output path and byte counts.

Signature:

```sql
bgzip(path, output_path := NULL, threads := 4, level := -1, keep := TRUE, overwrite := FALSE)
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM bgzip('regions.bed');
```

## bgunzip

Decompress a BGZF-compressed file and return the created output path and byte counts.

Signature:

```sql
bgunzip(path, output_path := NULL, threads := 4, keep := TRUE, overwrite := FALSE)
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM bgunzip('regions.bed.gz');
```

## bam_index

Build a BAM or CRAM index and report the written index path and format.

Signature:

```sql
bam_index(path, index_path := NULL, min_shift := 0, threads := 4)
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM bam_index('range.bam');
```

## bcf_index

Build a TBI or CSI index for a VCF or BCF file and report the written index path and format.

Signature:

```sql
bcf_index(path, index_path := NULL, min_shift := NULL, threads := 4)
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM bcf_index('formatcols.vcf.gz');
```

## tabix_index

Build a tabix index for a BGZF-compressed text file using a preset or explicit coordinate columns.

Signature:

```sql
tabix_index(path, preset := 'vcf', index_path := NULL, min_shift := 0, threads := 4, seq_col := NULL, start_col := NULL, end_col := NULL, comment_char := NULL, skip_lines := NULL)
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM tabix_index('gff_file.gff.gz', preset := 'gff');
```

## bam_bin_counts

Count BAM or CRAM read starts into fixed-width bins. Returns one row per bin across the selected contig span, including zero-count bins, with total, forward, and reverse counts; `rmdup := 'streaming'` applies the WisecondorX-style larp/larp2 consecutive-position deduplication, `rmdup := 'flag'` drops SAM duplicate-flagged reads, and `stats := 'gc'`, `'mq'`, or `'gc,mq'` adds per-bin pre/post-filter GC and MAPQ sufficient statistics, including reference GC when `reference` is provided.

Signature:

```sql
bam_bin_counts(path, bin_width, chrom := NULL, include_unmapped := FALSE, reference := NULL, index_path := NULL, mapq := 0, require_flags := 0, exclude_flags := 0, rmdup := 'none', stats := NULL)
```

Returns:

```
table
```

### Unmapped records

include_unmapped := TRUE appends one synthetic row with chrom = '*' for no-coordinate records, even when no such records are present; start, end and bin_id are NULL. The default FALSE omits this row. SAM-flag filters, MAPQ filtering and flag-based duplicate removal apply to these records; streaming duplicate removal is not applied to them. Reference-GC fields are NULL for this row.

### Examples

```sql
SELECT * FROM bam_bin_counts('fixture_mixed.cram', 5000, reference := 'fixture_ref.fa', rmdup := 'streaming', stats := 'gc,mq');
```

## duckhts_bam_mismatch_counts

Count aligned read bases against the reference by mate, cycle, base quality and substitution, optionally leaving out masked positions such as known variants.

Signature:

```sql
duckhts_bam_mismatch_counts(path, reference, region := NULL, mask := NULL, index_path := NULL, reference_index_path := NULL, mask_index_path := NULL, min_mapq := 20, require_flags := 0, exclude_flags := 3844, indel_flank := 5)
```

Returns:

```
table
```

### Input

path is a SAM, BAM or CRAM file and reference is its indexed FASTA; fasta_index() builds the index. region is a comma-separated list of contigs or intervals and needs an alignment index; without it the whole file is read. Alignments are kept when MAPQ is at least min_mapq, every bit of require_flags is set and no bit of exclude_flags is set; the default 3844 leaves out unmapped, secondary, failed, duplicate and supplementary alignments.

### Output

One row for each combination that occurs: mate (1 or 2 for a paired read, 0 otherwise), cycle (the one-based position of the base in the read as sequenced, hard clips included), base_quality (NULL when the read stores no qualities), reference_base and read_base (both on the strand of the read, so a reverse alignment is complemented), and bases, the number of aligned bases. Rows with reference_base equal to read_base are the matches, so a rate is a ratio of two sums. Rows come in the order of mate, cycle, base_quality, reference_base, read_base.

### What is counted

Every aligned base (CIGAR M, = or X) of a kept alignment whose read base and reference base are both A, C, G or T. A base within indel_flank read bases of an insertion, a deletion, a reference skip or a soft clip is left out, because bases next to a gap are often misaligned. A region selects alignments by overlap; all aligned bases of a selected alignment are counted, also those outside the region.

### Mask

mask is an indexed BCF or a bgzip-compressed VCF with a tabix index, for example known variants. Every reference position of a mask record is left out. A record with an indel or a symbolic allele also leaves out indel_flank bases on each side; a symbolic allele ends at INFO/END. With a mask of the known variants of the sample or its population, the mismatches that remain are errors of the read, the library or the alignment.

### Contig names

The contig names of the alignment header, the reference and the mask are compared byte for byte. A contig with a kept alignment that the reference does not have is an error. So is a contig that the mask does not know, in its header or in its index: a silent miss would count known variants as errors.

### Limits

One worker. Memory is fixed: 36.5 MB of counters, one reference window of 1 MiB with its mask, and the alignment in hand, plus htslib's state for the file. For CRAM that includes the reference of the slice being decoded; a slice over a sparse region can span many megabases. A cycle above 1000 is counted as 1000 and a base quality above 93 as 93. An alignment that spans more than 64 MiB of reference is an error. Alignments are expected in coordinate order; an unsorted file is counted correctly but fetches a reference window for each alignment.

### Examples

```sql
SELECT base_quality, sum(bases) FILTER (WHERE reference_base != read_base) / sum(bases) AS mismatch_rate FROM duckhts_bam_mismatch_counts('sample.bam', 'reference.fa', region := 'chr20:10000000-12000000', mask := 'known_variants.vcf.gz') GROUP BY base_quality ORDER BY base_quality;
```

## duckhts_bam_bed_coverage

Compute samtools coverage-like regional summaries for BAM or CRAM input over a BED target set, returning one row per BED interval with DuckHTS-specific pre/post-filter read counts, covered bases, percentage covered, mean depth, mean baseQ, mean mapQ, and strand-specific post-filter summaries in read mode. Indexed BAM/CRAM input is required in the current implementation. decompression_threads controls htslib worker threads for BAM/CRAM decoding; use 0 to disable them.

Signature:

```sql
duckhts_bam_bed_coverage(path, bed_path, reference := NULL, index_path := NULL, bed_index_path := NULL, mapq := 0, min_baseq := 0, min_read_len := 0, require_flags := 0, exclude_flags := 1796, min_depth := 1, max_depth := 1000000, decompression_threads := 0, fragment_mode := FALSE, strand_outputs := TRUE, processing_threads := 0)
```

Returns:

```
table
```

### Examples

```sql
SELECT chrom, start, "end", numreads_post, covbases_post, coverage_post FROM duckhts_bam_bed_coverage('fixture_mixed.bam', 'fixture_mixed_regions.bed');
```

## duckhts_mosdepth

Write mosdepth-compatible coverage files from indexed BAM/CRAM.

Signature:

```sql
duckhts_mosdepth(prefix, path, chrom := NULL, by := NULL, fasta := NULL, read_groups := NULL, no_per_base := FALSE, threads := 2, processing_threads := 2, flag := 1796, include_flag := 0, fast_mode := FALSE, fragment_mode := FALSE, use_median := FALSE, mapq := 0, min_frag_len := -1, max_frag_len := -1, precision_digits := 2, quantize := NULL, thresholds := NULL, index_path := NULL, overwrite := FALSE)
```

Returns:

```
table
```

### Output

Produce summary, global distribution, per-base BED.gz/CSI, optional window/BED-region results, quantized BED.gz/CSI and threshold counts for by. precision_digits sets text decimal places.

### Modes

Default fast_mode=FALSE uses CIGAR-aware coverage and mate-overlap correction. fragment_mode counts full insert spans of proper pairs; use_median switches by output from mean to median. read_groups filters comma-separated RG IDs; min_frag_len/max_frag_len filter absolute template length. Supply fasta when CRAM requires reference.

### Threads

processing_threads=0 is sequential; positive values select parallel contig-worker count.

### Examples

```sql
SELECT * FROM duckhts_mosdepth('sample', 'range.cram', fasta := 'ce.fa', by := '1000', fragment_mode := TRUE, read_groups := '1', use_median := TRUE, min_frag_len := 50, max_frag_len := 500, quantize := ':1:4:', thresholds := '1,10,20', precision_digits := 4, overwrite := TRUE);
```

## duckhts_samtools_idxstats

Write samtools idxstats-compatible TAB-delimited output for BAM, CRAM, or SAM input. Indexed BAM uses `hts_idx_get_stat(...)` for the fast path; CRAM, SAM, and unindexed BAM fall back to a full scan while preserving samtools-style contig rows plus the final `*` row.

Signature:

```sql
duckhts_samtools_idxstats(path, output := NULL, index_path := NULL, threads := 0, overwrite := FALSE)
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM duckhts_samtools_idxstats('range.bam', output := 'range.idxstats.txt', overwrite := TRUE);
```

## read_hts_header

Inspect HTS headers in parsed, raw, or combined form across supported formats. Raw VCF/BCF mode includes the final `#CHROM` sample header line so the returned text is suitable for Parquet metadata and future VCF/BCF regeneration.

Signature:

```sql
read_hts_header(path, format := NULL, mode := NULL)
```

Returns:

```
table
```

### Examples

```sql
SELECT record_type, id FROM read_hts_header('formatcols.vcf.gz') LIMIT 10;
```

## read_hts_index

Inspect high-level HTS index metadata such as sequence names and mapped counts.

Signature:

```sql
read_hts_index(path, format := NULL, index_path := NULL)
```

Returns:

```
table
```

### Examples

```sql
SELECT seqname, index_type FROM read_hts_index('vcf_file.bcf');
```

## read_hts_index_spans

Expand index metadata into span and chunk rows suitable for low-level index inspection.

Signature:

```sql
read_hts_index_spans(path, format := NULL, index_path := NULL)
```

Returns:

```
table
```

### Examples

```sql
SELECT seqname, chunk_beg_vo, chunk_end_vo FROM read_hts_index_spans('vcf_file.bcf') LIMIT 5;
```

## read_hts_index_raw

Return the raw on-disk HTS index blob together with basic identifying metadata.

Signature:

```sql
read_hts_index_raw(path, format := NULL, index_path := NULL)
```

Returns:

```
table
```

### Examples

```sql
SELECT length(raw) FROM read_hts_index_raw('formatcols.vcf.gz');
```

## variantkey

Encode a normalized biallelic variant as an official VariantKey-compatible 64-bit unsigned integer. This DuckHTS wrapper accepts 1-based VCF/DuckHTS POS to match bcftools `%VKX` / `+add-variantkey`, internally converts to the upstream 0-based field, and preserves the official hashed nonreversible mode for large, ambiguous, and symbolic REF/ALT strings. Only CHROM, POS, REF, and ALT are encoded; END, SVLEN, mate breakend coordinates, and other SV metadata are not.

Signature:

```sql
variantkey(chrom, pos, ref, alt)
```

Returns:

```
UBIGINT
```

### Chromosome code

The chromosome part of a VariantKey is a code of the upstream format, not a contig key. A leading chr is removed in any letter case. An all-digit name becomes its number modulo 256, with no range check. X, Y, and M or MT, in any letter case, become 23, 24 and 25. Any other name becomes 0. So 25 and MT share a code, and all scaffolds and accessions share code 0. Use duckhts_contig_key() to join on contig names.

### Examples

```sql
SELECT variantkey_hex(variantkey('1', 324684, 'C', 'G'));
```

## variantkey_hex

Render a VariantKey as its lowercase 16-character hexadecimal string representation.

Signature:

```sql
variantkey_hex(vk)
```

Returns:

```
VARCHAR
```

### Examples

```sql
SELECT variantkey_hex(variantkey('1', 324684, 'C', 'G'));
```

## parse_variantkey_hex

Parse a 16-character hexadecimal VariantKey string back into its UBIGINT code. Invalid or non-hex strings return NULL.

Signature:

```sql
parse_variantkey_hex(hex)
```

Returns:

```
UBIGINT
```

### Examples

```sql
SELECT parse_variantkey_hex('08027a2588b00000');
```

## encode_variantkey

Encode the raw upstream VariantKey fields directly: chromosome code, 0-based position, and 31-bit REF+ALT code.

Signature:

```sql
encode_variantkey(chrom_code, pos0, refalt_code)
```

Returns:

```
UBIGINT
```

### Examples

```sql
SELECT variantkey_hex(encode_variantkey(1, 324683, 145752064));
```

## extract_variantkey_chrom

Extract the raw upstream VariantKey chromosome code.

Signature:

```sql
extract_variantkey_chrom(vk)
```

Returns:

```
UTINYINT
```

### Examples

```sql
SELECT extract_variantkey_chrom(parse_variantkey_hex('08027a2588b00000'));
```

## extract_variantkey_pos

Extract the raw upstream VariantKey 0-based position field.

Signature:

```sql
extract_variantkey_pos(vk)
```

Returns:

```
UINTEGER
```

### Examples

```sql
SELECT extract_variantkey_pos(parse_variantkey_hex('08027a2588b00000'));
```

## extract_variantkey_refalt

Extract the raw upstream 31-bit VariantKey REF+ALT code.

Signature:

```sql
extract_variantkey_refalt(vk)
```

Returns:

```
UINTEGER
```

### Examples

```sql
SELECT extract_variantkey_refalt(parse_variantkey_hex('08027a2588b00000'));
```

## decode_variantkey

Decode a VariantKey into its raw upstream numeric fields: chrom_code, pos0, and refalt_code.

Signature:

```sql
decode_variantkey(vk)
```

Returns:

```
STRUCT
```

### Examples

```sql
SELECT (decode_variantkey(parse_variantkey_hex('08027a2588b00000'))).pos0;
```

## reverse_variantkey

Decode a VariantKey into a STRUCT with chrom, chrom_code, 1-based pos, upstream 0-based pos0, ref, alt, refalt_code, and reversible. For hashed nonreversible keys, reversible is FALSE and ref/alt are returned as NULL because DuckHTS v1 does not ship the optional NRVK lookup sidecar.

Signature:

```sql
reverse_variantkey(vk)
```

Returns:

```
STRUCT
```

### Examples

```sql
SELECT (reverse_variantkey(parse_variantkey_hex('08027a2588b00000'))).ref;
```

## variantkey_range

Return the inclusive minimum and maximum VariantKey bounds for a chromosome plus 1-based VCF position range, suitable for numeric range filtering on precomputed VariantKeys.

Signature:

```sql
variantkey_range(chrom, pos_min, pos_max)
```

Returns:

```
STRUCT
```

### Examples

```sql
SELECT variantkey_hex((variantkey_range('1', 100, 100)).min), variantkey_hex((variantkey_range('1', 100, 100)).max);
```

## regionkey

Encode a genomic interval as an official RegionKey-compatible 64-bit unsigned integer. Start and end use 0-based half-open interval semantics, matching BED-style coordinates; strand accepts -1, 0, or 1.

Signature:

```sql
regionkey(chrom, start, end, strand := 0)
```

Returns:

```
UBIGINT
```

### Examples

```sql
SELECT regionkey_hex(regionkey('X', 1007, 1807, 1));
```

## regionkey_hex

Render a RegionKey as its lowercase 16-character hexadecimal string representation.

Signature:

```sql
regionkey_hex(rk)
```

Returns:

```
VARCHAR
```

### Examples

```sql
SELECT regionkey_hex(regionkey('X', 1007, 1807, 1));
```

## parse_regionkey_hex

Parse a 16-character hexadecimal RegionKey string back into its UBIGINT code. Invalid or non-hex strings return NULL.

Signature:

```sql
parse_regionkey_hex(hex)
```

Returns:

```
UBIGINT
```

### Examples

```sql
SELECT parse_regionkey_hex('b80001f78000387a');
```

## encode_regionkey

Encode the raw upstream RegionKey fields directly: chromosome code, 0-based start, 0-based end, and strand code (0 = unknown, 1 = +, 2 = -).

Signature:

```sql
encode_regionkey(chrom_code, start, end, strand_code)
```

Returns:

```
UBIGINT
```

### Examples

```sql
SELECT regionkey_hex(encode_regionkey(23, 1007, 1807, 1));
```

## extract_regionkey_chrom

Extract the raw upstream RegionKey chromosome code.

Signature:

```sql
extract_regionkey_chrom(rk)
```

Returns:

```
UTINYINT
```

### Examples

```sql
SELECT extract_regionkey_chrom(parse_regionkey_hex('b80001f78000387a'));
```

## extract_regionkey_startpos

Extract the raw upstream RegionKey 0-based start position.

Signature:

```sql
extract_regionkey_startpos(rk)
```

Returns:

```
UINTEGER
```

### Examples

```sql
SELECT extract_regionkey_startpos(parse_regionkey_hex('b80001f78000387a'));
```

## extract_regionkey_endpos

Extract the raw upstream RegionKey 0-based end position.

Signature:

```sql
extract_regionkey_endpos(rk)
```

Returns:

```
UINTEGER
```

### Examples

```sql
SELECT extract_regionkey_endpos(parse_regionkey_hex('b80001f78000387a'));
```

## extract_regionkey_strand

Extract the raw upstream RegionKey strand code (0 = unknown, 1 = +, 2 = -).

Signature:

```sql
extract_regionkey_strand(rk)
```

Returns:

```
UTINYINT
```

### Examples

```sql
SELECT extract_regionkey_strand(parse_regionkey_hex('b80001f78000387a'));
```

## decode_regionkey

Decode a RegionKey into its raw upstream numeric fields: chrom_code, start, end, and strand_code.

Signature:

```sql
decode_regionkey(rk)
```

Returns:

```
STRUCT
```

### Examples

```sql
SELECT (decode_regionkey(parse_regionkey_hex('b80001f78000387a'))).strand_code;
```

## reverse_regionkey

Decode a RegionKey into a STRUCT with chrom, chrom_code, start, end, strand, and strand_code.

Signature:

```sql
reverse_regionkey(rk)
```

Returns:

```
STRUCT
```

### Examples

```sql
SELECT (reverse_regionkey(parse_regionkey_hex('b80001f78000387a'))).strand;
```

## extend_regionkey

Extend a RegionKey interval by a fixed number of bases on both sides, clamping to the official 28-bit RegionKey position range.

Signature:

```sql
extend_regionkey(rk, size)
```

Returns:

```
UBIGINT
```

### Examples

```sql
SELECT reverse_regionkey(extend_regionkey(regionkey('X', 10000, 20000, -1), 1000));
```

## duckhts_contig_key

Return a conservative contig join key by removing one non-empty leading chr prefix case-insensitively and normalizing M/MT to MT. X and Y are uppercased; all other suffixes are preserved. This does not map numeric sex chromosomes, accessions, patches, or alternate loci.

Signature:

```sql
duckhts_contig_key(contig)
```

Returns:

```
VARCHAR
```

### Use

Readers and region arguments compare contig names byte for byte. This function is the explicit key for a join of two sources that spell contigs differently: join on duckhts_contig_key(a.chrom) = duckhts_contig_key(b.chrom), or rename one side. Two names of one source can share a key (1 and chr1), so a join on the key repeats rows when one side has both. duckhts_roh_ancestry and the R ancestry wrappers use this key. Functions with another rule state it under Contig names.

### Not mapped

Leading zeros (01 is not 1), numeric sex and mitochondrial codes (23, 24, 25, 26), accessions such as NC_000001.11, patches, alternate loci and unplaced scaffolds. The letter case of other names is kept.

### Examples

```sql
SELECT duckhts_contig_key('chr1'), duckhts_contig_key('chrM');
```

## are_overlapping_regions

Return TRUE when two explicit 0-based half-open intervals overlap on the same canonical chromosome.

Signature:

```sql
are_overlapping_regions(chrom_a, start_a, end_a, chrom_b, start_b, end_b)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT are_overlapping_regions('1', 2, 4, '1', 3, 7);
```

## are_overlapping_region_regionkey

Return TRUE when a 0-based half-open interval overlaps the supplied RegionKey interval.

Signature:

```sql
are_overlapping_region_regionkey(chrom, start, end, rk)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT are_overlapping_region_regionkey('X', 1008, 1800, parse_regionkey_hex('b80001f78000387a'));
```

## are_overlapping_regionkeys

Return TRUE when two RegionKeys overlap.

Signature:

```sql
are_overlapping_regionkeys(rka, rkb)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT are_overlapping_regionkeys(regionkey('X', 1007, 1807, 1), parse_regionkey_hex('b80001f78000387a'));
```

## bcftools_liftover

Row-oriented liftover kernel intended to mirror bcftools +liftover semantics as closely as possible while returning one STRUCT per input row with fields: src_chrom, src_pos, src_ref, src_alt, dest_chrom, dest_pos, dest_end, dest_ref, dest_alt, mapped, reverse_complemented, swap, reject_reason, and note. Set no_left_align := true to skip post-liftover left-alignment of lifted indels (mirrors --no-left-align in bcftools +liftover).

Signature:

```sql
bcftools_liftover(chrom, pos, ref, alt, chain_path, dst_fasta_ref, src_fasta_ref, max_snp_gap, max_indel_inc, lift_mt, end_pos, no_left_align)
```

Returns:

```
STRUCT
```

### Contig names

Input names and chain source names are compared after one rule: a leading chr is removed in any letter case; M, MT and 26 become MT; 23, 25, X, XY, XX, PAR1 and PAR2 become X; 24 and Y become Y. A lifted record carries the first name that the destination FASTA index has among: the chain's destination name; that name with chr added or removed; its form under the rule above, without and with chr; and for the mitochondrion MT, chrM and M. A FASTA that holds two of these names is not an error. A mitochondrial record that passes through without lifting takes the first of chrM, MT and M that the destination FASTA has.

### Examples

```sql
SELECT (bcftools_liftover(chrom, pos, ref, alt, 'hg19ToHg38.over.chain.gz', 'hg38.fa', 'hg19.fa', 1, 250, false, NULL::BIGINT, false)).dest_pos FROM variants;
```

## duckdb_liftover

DuckDB-specific wrapper over bcftools_liftover that takes either a table name or a derived-table expression plus column-name strings for chrom/pos/ref/alt and returns the lifted table. The no_left_align parameter mirrors --no-left-align in bcftools +liftover.

Signature:

```sql
duckdb_liftover(table_name, chrom_col, pos_col, ref_col := NULL, alt_col := NULL, chain_path := NULL, dst_fasta_ref := NULL, src_fasta_ref := NULL, max_snp_gap := 1, max_indel_inc := 250, lift_mt := false, end_pos_col := NULL, no_left_align := false)
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM duckdb_liftover('variants', 'chrom', 'pos', ref_col := 'ref', alt_col := 'alt', chain_path := 'hg19ToHg38.over.chain.gz', dst_fasta_ref := 'hg38.fa');
```

```sql
SELECT * FROM duckdb_liftover('(SELECT chrom, pos, ref, alt FROM variants) AS v', 'chrom', 'pos', ref_col := 'ref', alt_col := 'alt', chain_path := 'hg19ToHg38.over.chain.gz', dst_fasta_ref := 'hg38.fa');
```

## bcftools_norm_row

Normalize one variant against FASTA with bcftools/vt-style left alignment.

Signature:

```sql
bcftools_norm_row(chrom, pos, ref, alt, fasta_ref, end_pos := NULL, svlen := NULL, fasta_index_path := NULL, gzi_path := NULL)
```

Returns:

```
STRUCT
```

### Input

alt accepts comma-delimited VARCHAR or VARCHAR[]. Symbolic `<DEL>` may use end_pos; `<DUP>` may use svlen.

### Output

Return pos_normed, end_pos_normed, ref_normed, alt_normed (always VARCHAR[]), nullable normed and norm_status. gVCF `<NON_REF>`/`<*>` reference blocks pass through with GVCFReferenceBlock. Mixed real/gVCF-symbolic rows normalize real alleles while preserving symbolic alleles and supplied reference-block END.

### Contig names

The reference sequence of a record is looked up in the FASTA index under, in order: the record's name; the name without chr (any letter case), or with chr added when it has none; and for M, MT or chrM also MT, chrM and M. The first name that the index has is used. A FASTA that holds both 1 and chr1 is not an error.

### Examples

```sql
SELECT nr.pos_normed, nr.ref_normed, nr.alt_normed FROM (SELECT bcftools_norm_row('chr1', 100, 'AC', 'A', 'ref.fa', NULL::BIGINT, NULL::BIGINT, NULL, NULL) AS nr);
```

## duckhts_bcftools_norm

Normalize variants from a table or derived-table expression while preserving input columns.

Signature:

```sql
duckhts_bcftools_norm(table_name, fasta_ref, chrom_col := 'chrom', pos_col := 'pos', ref_col := 'ref', alt_col := 'alt', split_multiallelic := FALSE, end_pos_col := NULL, svlen_col := NULL, fasta_index_path := NULL, gzi_path := NULL)
```

Returns:

```
table
```

### Output

ALT accepts VARCHAR or VARCHAR[]. Append pos_normed, end_pos_normed, ref_normed, alt_normed, normed and norm_status. split_multiallelic=TRUE splits sites before normalization; alt_normed becomes VARCHAR and alt_index is added.

### Scope

This wraps bcftools_norm_row, not a full-record VCF/BCF rewrite. GT, PL/GP/DS and PS remain unchanged caller columns unless a separate writer/remapper updates them.

### Examples

```sql
SELECT * FROM duckhts_bcftools_norm('variants', 'ref.fa');
```

```sql
SELECT * FROM duckhts_bcftools_norm('(SELECT CHROM AS chrom, POS AS pos, REF AS ref, ALT AS alt FROM read_bcf(''cohort.vcf.gz'')) AS v', 'ref.fa', split_multiallelic := TRUE);
```

## bcftools_score

Compute polygenic scores from genotype VCF/BCF and summary statistics using bcftools +score dosage semantics.

Signature:

```sql
bcftools_score(bcf_path, summary_path_or_list, use := NULL, columns := 'PLINK', columns_file := NULL, q_score_thr := NULL, summaries_list_file := NULL, log_path := NULL, use_variant_id := FALSE, counts := FALSE, samples := NULL, force_samples := FALSE, regions := NULL, regions_file := NULL, regions_overlap := 1, targets := NULL, targets_file := NULL, targets_overlap := 0, apply_filters := NULL, include := NULL, exclude := NULL)
```

Returns:

```
table
```

### Input

Support GT/DS/HDS/AP/GP/AS dosage, sample subsets and region/target/FILTER-string controls. The second argument accepts one path or a list. TSV/SSF inputs yield one PRS column per file in a single genotype scan; GWAS-VCF yields one per FORMAT sample.

### Summary discovery

With NULL second argument, summaries_list_file reads paths from a file or directory. List entries are interpreted as written; directories scan supported regular files lexicographically and omit index sidecars.

### Audit

log_path writes per-PRS loaded/matched/allele-mismatch/duplicate-marker counts.

### Contig names

Chromosome names of the summary statistics are looked up in the genotype VCF header as the upstream plugin does, in this order: the exact name; the name without a lowercase chr prefix; chr plus a name of at most two characters; 23, 25, XY, XX, PAR1 and PAR2 as X or chrX; 24 as Y or chrY; 26, MT and chrM as MT or chrM. The lookup is case-sensitive, and M alone does not find MT. A marker whose chromosome is not found is skipped without an error.

### Examples

```sql
SELECT * FROM bcftools_score('cohort.bcf', 'gwas.tsv.gz', columns := 'PLINK') LIMIT 5;
```

```sql
SELECT * FROM bcftools_score('cohort.bcf', ['score1.tsv.gz', 'score2.tsv.gz'], columns := 'GWAS-SSF') LIMIT 5;
```

```sql
SELECT * FROM bcftools_score('cohort.bcf', NULL, columns := 'GWAS-SSF', summaries_list_file := 'scores.list', log_path := 'score.log') LIMIT 5;
```

## bcftools_munge_row

Normalize one summary-statistics row into GWAS-VCF-style fields (chrom/pos/ref/alt/effect metrics), resolving REF/ALT orientation against a FASTA reference and applying swap-aware sign/frequency/count transforms. The output flag `alleles_swapped` means REF/ALT orientation was swapped to match the FASTA reference.

Signature:

```sql
bcftools_munge_row(chrom, pos, a1, a2, id, p, z, or, beta, n, n_cas, n_con, info, frq, se, lp, ac, neff, neffdiv2, het_i2, het_p, het_lp, dire, fasta_ref, iffy_tag := 'IFFY', mismatch_tag := 'REF_MISMATCH', ns := NULL, nc := NULL, ne := NULL)
```

Returns:

```
STRUCT
```

### Contig names

The reference sequence of a record is looked up in the FASTA index under, in order: the record's name; the name without chr (any letter case), or with chr added when it has none; and for M, MT or chrM also MT, chrM and M. The first name that the index has is used. A FASTA that holds both 1 and chr1 is not an error.

### Examples

```sql
SELECT (bcftools_munge_row('chr1', 12345, 'A', 'G', 'rs1', 0.01, NULL, NULL, 0.12, 1000, NULL, NULL, NULL, 0.3, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 'ref.fa')).ref;
```

## duckdb_munge

DuckDB macro wrapper over bcftools_munge_row that maps source columns (via preset or explicit map) and returns normalized GWAS-VCF-style rows with lean outputs and explicit `alleles_swapped` semantics. Output columns: chrom, pos, id, ref, alt, alleles_swapped, filter, ns, ez, nc, es, se, lp, af, ac, ne (16 columns). For METAL meta-analysis output with SI/I2/CQ/ED columns, use duckdb_munge_metal.

Signature:

```sql
duckdb_munge(table_name, preset := '', column_map := map([''], ['']), column_map_file := '', fasta_ref := NULL, iffy_tag := 'IFFY', mismatch_tag := 'REF_MISMATCH', ns := NULL, nc := NULL, ne := NULL)
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM duckdb_munge('gwas_table', preset := 'PLINK', fasta_ref := 'ref.fa');
```

```sql
SELECT * FROM duckdb_munge('(SELECT * FROM gwas_table) AS s', column_map := map(['CHR','BP','A1','A2','SNP'], ['chrom','pos','ea','nea','id']), fasta_ref := 'ref.fa');
```

## duckdb_munge_metal

Extended munge macro with METAL meta-analysis output columns. Same as duckdb_munge but additionally emits: si (imputation info, from INFO input), i2 (Cochran's I² heterogeneity, from HET_I2), cq (Cochran's Q -log10 p, from HET_LP or -log10(HET_P)), and ed (effect direction string, from DIRE; +/- flipped on allele swap). The R wrapper rduckhts_munge() auto-dispatches to this macro when metal keys (INFO, HET_I2, HET_P, HET_LP, DIRE) are present in the resolved column map.

Signature:

```sql
duckdb_munge_metal(table_name, preset := '', column_map := map([''], ['']), column_map_file := '', fasta_ref := NULL, iffy_tag := 'IFFY', mismatch_tag := 'REF_MISMATCH', ns := NULL, nc := NULL, ne := NULL)
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM duckdb_munge_metal('metal_results', preset := 'METAL', fasta_ref := 'ref.fa');
```

```sql
SELECT * FROM duckdb_munge_metal('meta', column_map := map(['CHR','BP','A1','A2','SNP','HET_I2','DIRE'], ['chr','pos','ea','nea','snp','HetISq','Direction']), fasta_ref := 'ref.fa');
```

## hts_union_query

Generate a UNION ALL BY NAME query string that reads every file matching a glob pattern through the named reader function. The result includes a 'filename' column identifying the source file for each row. Assign to a variable with SET VARIABLE and execute via query(getvariable(...)). Optional params string is appended to each reader call. In R, use the typed rduckhts_*_multi() helpers instead, which accept file vectors with optional per-file parameters and create DuckDB tables directly.

Signature:

```sql
hts_union_query(reader, pattern, params := '')
```

Returns:

```
VARCHAR
```

### Examples

```sql
SET VARIABLE q = hts_union_query('read_bam', 'samples/*.bam'); SELECT * FROM query(getvariable('q'));
```

```sql
SET VARIABLE q = hts_union_query('read_bcf', 'cohort/*.vcf.gz', 'tidy_format := true'); SELECT * FROM query(getvariable('q'));
```

## hts_region_union_query

Generate UNION ALL BY NAME SQL over separate per-region scans of one HTS file.

Signature:

```sql
hts_region_union_query(reader, path, regions, params := '')
```

Returns:

```
VARCHAR
```

### Input

regions is a list of region strings. params is appended to each reader call and must not include region. Output adds filename, duckhts_region_shard_id and duckhts_region_shard for shard provenance.

### Duplicates

UNION ALL does not deduplicate. Adjacent or overlapping shards can repeat spanning BAM/VCF/BCF records; apply shard-local filters or explicit downstream deduplication when exactly-once output is required.

### Examples

```sql
SET VARIABLE q = hts_region_union_query('read_bam', 'sample.bam', ['chr1:1-1000000','chr1:1000001-2000000']); SELECT * FROM query(getvariable('q'));
```

```sql
SET VARIABLE q = hts_region_union_query('read_bcf', 'cohort.vcf.gz', string_split('22:16000000-16999999,22:17000000-17999999', ','), 'tidy_format := true'); SELECT * FROM query(getvariable('q'));
```

## seq_revcomp

Compute the reverse complement of a DNA sequence using A, C, G, T, and N bases. Overloaded: accepts either a VARCHAR text sequence (returns VARCHAR) or a UTINYINT[] of htslib nt16 codes as produced by read_bam(sequence_encoding := 'nt16') (returns UTINYINT[]); the nt16 overload is bit-identical to the text path after decoding, so BAM pipelines can reverse-complement without leaving the nt16 encoding.

Signature:

```sql
seq_revcomp(sequence)
```

Returns:

```
VARCHAR
```

### Examples

```sql
SELECT seq_revcomp('ACGTN');
```

```sql
SELECT seq_revcomp(SEQ) FROM read_bam('reads.bam', sequence_encoding := 'nt16');
```

## seq_canonical

Return the lexicographically smaller of a sequence and its reverse complement. Overloaded: accepts either a VARCHAR text sequence (returns VARCHAR) or a UTINYINT[] of htslib nt16 codes as produced by read_bam(sequence_encoding := 'nt16') (returns UTINYINT[]); the nt16 overload compares by decoded base order and is bit-identical to the text path after decoding.

Signature:

```sql
seq_canonical(sequence)
```

Returns:

```
VARCHAR
```

### Examples

```sql
SELECT seq_canonical('ACGTN');
```

```sql
SELECT seq_canonical(SEQ) FROM read_bam('reads.bam', sequence_encoding := 'nt16');
```

## seq_hash_2bit

Encode a short DNA sequence as a 2-bit unsigned integer hash. Overloaded to also accept a UTINYINT[] of htslib nt16 codes (from read_bam(sequence_encoding := 'nt16')); non-ACGT codes yield NULL, bit-identical to the text path.

Signature:

```sql
seq_hash_2bit(sequence)
```

Returns:

```
UBIGINT
```

### Examples

```sql
SELECT seq_hash_2bit('ACGT');
```

## seq_encode_4bit

Encode an IUPAC DNA sequence as a list of 4-bit base codes, preserving ambiguity symbols including N.

Signature:

```sql
seq_encode_4bit(sequence)
```

Returns:

```
UTINYINT[]
```

### Examples

```sql
SELECT seq_encode_4bit('ACGTRYSWKMBDHVN');
```

## seq_decode_4bit

Decode a list of 4-bit IUPAC DNA base codes back into a sequence string.

Signature:

```sql
seq_decode_4bit(codes)
```

Returns:

```
VARCHAR
```

### Examples

```sql
SELECT seq_decode_4bit(seq_encode_4bit('ACGTRYSWKMBDHVN'));
```

## seq_gc_content

Compute GC fraction for a DNA sequence as a value between 0 and 1. Overloaded: accepts either a VARCHAR text sequence or a UTINYINT[] of htslib nt16 codes as produced by read_bam(sequence_encoding := 'nt16'); the nt16 overload classifies codes directly and is bit-identical to the text path, so BAM pipelines can compute GC without decoding sequences back to text.

Signature:

```sql
seq_gc_content(sequence)
```

Returns:

```
DOUBLE
```

### Examples

```sql
SELECT seq_gc_content('ACGT');
```

```sql
SELECT seq_gc_content(SEQ) FROM read_bam('reads.bam', sequence_encoding := 'nt16');
```

## seq_kmers

Expand a sequence into positional k-mers with optional canonicalization.

Signature:

```sql
seq_kmers(sequence, k, canonical := FALSE)
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM seq_kmers('ACGT', 2);
```

## sam_flag_bits

Decode a SAM flag into a struct of boolean bit fields using explicit SAM-oriented names such as `is_paired`, `is_proper_pair`, `is_next_segment_unmapped`, and `is_supplementary`.

Signature:

```sql
sam_flag_bits(flag)
```

Returns:

```
STRUCT
```

### Examples

```sql
SELECT (sam_flag_bits(99)).is_proper_pair;
```

## sam_flag_has

Test whether any bits from the provided SAM flag mask are set in a flag value.

Signature:

```sql
sam_flag_has(flag, mask)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT sam_flag_has(99, 2);
```

## is_forward_aligned

Test whether a mapped segment is aligned to the forward strand. Returns `NULL` for unmapped segments because SAM flag `0x10` does not define genomic strand when `0x4` is set.

Signature:

```sql
is_forward_aligned(flag)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT is_forward_aligned(0);
```

## cigar_has_soft_clip

Test whether a CIGAR string contains any soft-clipped segment (`S`). Overloaded to also accept a UINTEGER[] binary CIGAR (as produced by read_bam(cigar_representation := 'binary')); the binary overload is bit-identical to the text path.

Signature:

```sql
cigar_has_soft_clip(cigar[, strict])
```

Returns:

```
BOOLEAN
```

### Input validation

See cigar_query_length for full-input validation and missing-input behavior.

### Examples

```sql
SELECT cigar_has_soft_clip('5S90M5S');
```

## cigar_has_hard_clip

Test whether a CIGAR string contains any hard-clipped segment (`H`). Overloaded to also accept a UINTEGER[] binary CIGAR (as produced by read_bam(cigar_representation := 'binary')); the binary overload is bit-identical to the text path.

Signature:

```sql
cigar_has_hard_clip(cigar[, strict])
```

Returns:

```
BOOLEAN
```

### Input validation

See cigar_query_length for full-input validation and missing-input behavior.

### Examples

```sql
SELECT cigar_has_hard_clip('5H95M');
```

## cigar_left_soft_clip

Return the left-end soft-clipped length from a CIGAR string, or zero if the alignment does not start with `S`. Overloaded to also accept a UINTEGER[] binary CIGAR (as produced by read_bam(cigar_representation := 'binary')); the binary overload is bit-identical to the text path.

Signature:

```sql
cigar_left_soft_clip(cigar[, strict])
```

Returns:

```
BIGINT
```

### Input validation

See cigar_query_length for full-input validation and missing-input behavior. The literal first op determines the left soft clip; a leading H is not skipped.

### Examples

```sql
SELECT cigar_left_soft_clip('5S90M5S');
```

## cigar_right_soft_clip

Return the right-end soft-clipped length from a CIGAR string, or zero if the alignment does not end with `S`. Overloaded to also accept a UINTEGER[] binary CIGAR (as produced by read_bam(cigar_representation := 'binary')); the binary overload is bit-identical to the text path.

Signature:

```sql
cigar_right_soft_clip(cigar[, strict])
```

Returns:

```
BIGINT
```

### Input validation

See cigar_query_length for full-input validation and missing-input behavior. The literal last op determines the right soft clip; a trailing H is not skipped.

### Examples

```sql
SELECT cigar_right_soft_clip('5S90M5S');
```

## cigar_query_length

Return the query-consuming length from a CIGAR string, counting `M`, `I`, `S`, `=`, and `X`. Overloaded to also accept a UINTEGER[] binary CIGAR (as produced by read_bam(cigar_representation := 'binary')); the binary overload is bit-identical to the text path.

Signature:

```sql
cigar_query_length(cigar[, strict])
```

Returns:

```
BIGINT
```

### Input validation

Text and packed CIGARs are validated in full. Supported ops are M, I, D, N, S, H, P, = and X, each with a positive length that fits BIGINT. Consumed query and reference spans must each fit BIGINT. Invalid ops, missing lengths, trailing digits, arithmetic overflow and NULL packed elements are invalid input. This validates operation syntax and numeric ranges, not biological ordering constraints.

### Failure policy

The optional final positional BOOLEAN strict defaults to FALSE: invalid input returns NULL. TRUE uses the same grammar but raises a DuckDB error naming the function and failure. Where an operation can be identified, diagnostics give its 1-based packed-op index or the 1-based byte at which the text operation starts. No read identifier is inferred. A top-level SQL NULL argument, including strict, returns NULL. Empty text, '*' and an empty packed list return NULL in either policy.

### Examples

```sql
SELECT cigar_query_length('5S90M5I');
```

```sql
SELECT cigar_query_length([84, 1440, 81]::UINTEGER[], TRUE);
```

## cigar_aligned_query_length

Return the aligned query length from a CIGAR string, counting `M`, `=`, and `X` but excluding clips and insertions. Overloaded to also accept a UINTEGER[] binary CIGAR (as produced by read_bam(cigar_representation := 'binary')); the binary overload is bit-identical to the text path.

Signature:

```sql
cigar_aligned_query_length(cigar[, strict])
```

Returns:

```
BIGINT
```

### Input validation

See cigar_query_length for full-input validation and missing-input behavior.

### Examples

```sql
SELECT cigar_aligned_query_length('5S90M5I');
```

## cigar_reference_length

Return the reference-consuming length from a CIGAR string, counting `M`, `D`, `N`, `=`, and `X`. Overloaded to also accept a UINTEGER[] binary CIGAR (as produced by read_bam(cigar_representation := 'binary')); the binary overload is bit-identical to the text path.

Signature:

```sql
cigar_reference_length(cigar[, strict])
```

Returns:

```
BIGINT
```

### Input validation

See cigar_query_length for full-input validation and missing-input behavior.

### Examples

```sql
SELECT cigar_reference_length('90M5D');
```

## cigar_has_op

Test whether a CIGAR string contains at least one instance of the requested operator. Overloaded to also accept a UINTEGER[] binary CIGAR (as produced by read_bam(cigar_representation := 'binary')); the binary overload is bit-identical to the text path.

Signature:

```sql
cigar_has_op(cigar, op[, strict])
```

Returns:

```
BOOLEAN
```

### Input validation

Uses the full-input validation of cigar_query_length, including the suffix after any matching op. The requested operator is one supported ASCII character, case-insensitive independently of the process locale.

### Failure policy

The optional final positional BOOLEAN strict defaults to FALSE: an invalid operator or malformed CIGAR returns NULL. TRUE raises an error with the diagnostic conventions of cigar_query_length. A top-level SQL NULL argument, including strict, returns NULL. Empty text, '*' and an empty packed list return false for a valid requested operator in either policy.

### Examples

```sql
SELECT cigar_has_op('5S90M5S', 'S');
```

## cigar_aligned_blocks

Return the aligned blocks of a CIGAR as a struct of three parallel BIGINT lists: ref_start, query_start and width, one entry per M, = or X op in CIGAR order. Overloaded to also accept a UINTEGER[] binary CIGAR (as produced by read_bam(cigar_representation := 'binary')); the binary overload is bit-identical to the text path.

Signature:

```sql
cigar_aligned_blocks(cigar, pos[, strict])
```

Returns:

```
STRUCT
```

### Coordinates

ref_start is pos plus the reference bases consumed before the block, so it carries whatever base pos uses; pass read_bam's 1-based POS for 1-based starts or 0 for offsets from the alignment start. query_start is the 0-based offset into the stored SEQ: soft clips count, hard clips do not. width is the op length. D and N advance the reference and split blocks, I advances the query and splits blocks, H and P consume nothing. Blocks are never merged, as in pysam get_blocks() and GenomicAlignments cigarRangesAlongReferenceSpace over M, = and X.

### Input validation

Uses the full-input CIGAR grammar of cigar_query_length. The half-open reference end, pos + cigar_reference_length(cigar), must also fit BIGINT. A negative pos is permitted. A valid CIGAR with no aligned op returns three empty lists.

### Failure policy

The optional final positional BOOLEAN strict defaults to FALSE: malformed CIGAR or coordinate overflow returns NULL. TRUE raises an error with the diagnostic conventions of cigar_query_length. A top-level SQL NULL argument, including pos or strict, returns NULL. Empty text, '*' and an empty packed list return NULL in either policy.

### Examples

```sql
SELECT (cigar_aligned_blocks('5S90M5S', 100)).ref_start;
```

```sql
SELECT UNNEST((b).ref_start) AS ref_start, UNNEST((b).width) AS width FROM (SELECT cigar_aligned_blocks(CIGAR, POS) AS b FROM read_bam('reads.bam', cigar_representation := 'binary'));
```

## is_paired

Test whether the SAM flag indicates that the template has multiple segments in sequencing (`0x1`).

Signature:

```sql
is_paired(flag)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT is_paired(99);
```

## is_proper_pair

Test whether the SAM flag indicates that each segment is properly aligned according to the aligner (`0x2`).

Signature:

```sql
is_proper_pair(flag)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT is_proper_pair(99);
```

## is_unmapped

Test whether the read itself is unmapped according to the SAM flag.

Signature:

```sql
is_unmapped(flag)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT is_unmapped(4);
```

## is_next_segment_unmapped

Test whether the next segment in the template is flagged as unmapped (`0x8`).

Signature:

```sql
is_next_segment_unmapped(flag)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT is_next_segment_unmapped(9);
```

## is_reverse_complemented

Test whether `SEQ` is stored reverse complemented (`0x10`); for mapped reads this corresponds to reverse-strand alignment.

Signature:

```sql
is_reverse_complemented(flag)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT is_reverse_complemented(16);
```

## is_next_segment_reverse_complemented

Test whether `SEQ` of the next segment in the template is stored reverse complemented (`0x20`).

Signature:

```sql
is_next_segment_reverse_complemented(flag)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT is_next_segment_reverse_complemented(32);
```

## is_first_segment

Test whether the read is marked as the first segment in the template.

Signature:

```sql
is_first_segment(flag)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT is_first_segment(64);
```

## is_last_segment

Test whether the read is marked as the last segment in the template.

Signature:

```sql
is_last_segment(flag)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT is_last_segment(128);
```

## is_secondary

Test whether the alignment is marked as secondary.

Signature:

```sql
is_secondary(flag)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT is_secondary(256);
```

## is_qc_fail

Test whether the read failed vendor or pipeline quality checks.

Signature:

```sql
is_qc_fail(flag)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT is_qc_fail(512);
```

## is_duplicate

Test whether the alignment is flagged as a duplicate.

Signature:

```sql
is_duplicate(flag)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT is_duplicate(1024);
```

## is_supplementary

Test whether the alignment is marked as supplementary.

Signature:

```sql
is_supplementary(flag)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT is_supplementary(2048);
```
