# Extension Function Reference

Generated from `functions.yaml`.

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
duckhts_duckdb_type_supported(type_name)
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

## duckvep_ensembl_regions

Match tiled FASTA sequence to one Ensembl core assembly and assign dense model-local sequence-region ordinals.

Signature:

```sql
duckvep_ensembl_regions(core_schema, reference_chunks_table, assembly, species_id := 1)
```

Returns:

```
table
```

### Input

core_schema supplies seq_region and coord_system. reference_chunks_table supplies chrom, zero-based start, half-open end and seq, normally from fasta_nuc(..., include_seq := TRUE). Chunks must be contiguous from zero with matching sequence lengths; every FASTA contig must match exactly one same-length Ensembl region.

### Model paths

Only supplied FASTA paths enter the model. Alternate haplotypes and patches require their named reference chunks. No assembly_exception projection, X/Y PAR merging or circular-region marking is performed. Ordinary mitochondrial coordinates are accepted; origin-crossing intervals are not.

### Examples

```sql
CREATE TABLE grch38_reference_chunks AS SELECT chrom, start, "end", seq FROM fasta_nuc('GRCh38.primary.fa', bin_width := 1048576, include_seq := TRUE);
```

```sql
CREATE TABLE grch38_regions AS SELECT * FROM duckvep_ensembl_regions('ensembl_core', 'grch38_reference_chunks', 'GRCh38');
```

## duckvep_ensembl_transcripts

Build validated VEP-116 Ensembl core transcript models from core tables and matching tiled FASTA sequence.

Signature:

```sql
duckvep_ensembl_transcripts(core_schema, reference_chunks_table, assembly, species_id := 1)
```

Returns:

```
table
```

### Selection

Before dense ordinals are assigned, retain current transcripts with nonempty stable IDs and exclude artifact biotypes and readthrough_tra. This builds the Ensembl core selection, not VEP --refseq or --merged sets. GENCODE Basic/Primary membership is retained without filtering to those sets.

### Output

Return resident-loader columns, stable/source IDs, biotypes, versioned MANE Select/Plus Clinical RefSeq accessions when present, ranked nested exons, mature-miRNA cDNA attributes projected into genomic exon segments and supported Translation SeqEdits. MANE accessions remain relational; selection flags enter transcript_flags. Sequence contains phase-adjusted CDS and complete transcript-oriented pre-CDS/post-CDS spliced sequence.

### Translation

Codon tables come from seq_region_attrib, defaulting to table 1 as in VEP 116. All BioPerl/VEP-supported NCBI table IDs are accepted; invalid/conflicting attributes reject import. Single-residue initial_met, _selenocysteine, amino_acid_sub and _stop_codon_rt edits form a reference-peptide overlay. Transcript sequence corrections and unsupported Translation SeqEdit shapes withhold sequence with an explicit reason.

### Model paths

Transcripts retain their selected sequence-path identities, including X/Y PAR and alternate haplotypes. No alt_allele collapse, IS_PAR-based result merging or circular-origin lifting is performed.

### Examples

```sql
CREATE TABLE grch38_transcripts AS SELECT * FROM duckvep_ensembl_transcripts('ensembl_core', 'grch38_reference_chunks', 'GRCh38');
```

## duckvep_ensembl_regulation_features

Prepare VEP-116 RegulatoryFeature and MotifFeature intervals for a DuckVEP model.

Signature:

```sql
duckvep_ensembl_regulation_features(funcgen_schema, regions_table)
```

Returns:

```
table
```

### Input

funcgen_schema supplies regulatory_feature, feature_type and motif_feature; regions_table is prepared by duckvep_ensembl_regions(). EMAR RegulatoryFeature rows are excluded as in VEP before dense ordinals are assigned.

### Output

Validate one-based inclusive coordinates and source feature types, map source regions to model ordinals, and retain stable/source IDs, feature metadata, binding-matrix/regulatory-build IDs and motif scores. feature_kind is 1 for regulatory regions or 2 for transcription-factor binding sites.

### Loading

Pass regulation_feature_index, seq_region, feature_start, feature_end and feature_kind to duckvep_model_load(interval_feature_query := ...). Other metadata remains relational for joins.

### Examples

```sql
CREATE TABLE grch38_regulation AS SELECT * FROM duckvep_ensembl_regulation_features('ensembl_funcgen', 'grch38_regions');
```

## duckvep_model_receipt

Create a deterministic provenance receipt and semantic hash for prepared DuckVEP model relations.

Signature:

```sql
duckvep_model_receipt(regions_table, transcripts_table, source_name, source_version, assembly, source_manifest_sha256, reference_sha256, transcript_filter, regulation_features_table := NULL)
```

Returns:

```
table
```

### Output

Return declared source_name, source_version, assembly, source_manifest_sha256, reference_sha256 and transcript_filter; count transcript, sequence, mature-miRNA and peptide-edit content; hash resident-model fields in stable ordinal order. No clock time is included, so identical inputs yield the same model hash.

### Regulation

When regulation_features_table is supplied, validate dense ordinals, region agreement, interval geometry and feature kinds; count regulatory/motif objects and include their five resident columns in the hash.

### Provenance

source_manifest_sha256 should identify every exact core/funcgen input. Record primary versus alternate/patch reference paths and any assembly_exception preprocessing in the source manifest or transcript_filter.

### Examples

```sql
CREATE TABLE grch38_receipt AS SELECT * FROM duckvep_model_receipt('grch38_regions', 'grch38_transcripts', 'Ensembl', '116', 'GRCh38', source_manifest_sha256, reference_sha256, 'VEP 116 core selection on FASTA-covered assembly regions', regulation_features_table := 'grch38_regulation');
```

## duckvep_model_load

Load a validated immutable consequence model under a name in the current DuckDB database; return one TRUE row.

Signature:

```sql
duckvep_model_load(name, sequence_region_query, transcript_query, exon_query, mature_mirna_query := NULL, peptide_edit_query := NULL, interval_feature_query := NULL, reference_fasta := NULL, transcript_coverage_complete := FALSE)
```

Returns:

```
table(loaded BOOLEAN)
```

### Required queries

Three query strings read committed, non-temporary relations: sorted UINTEGER sequence-region ordinals; dense transcripts with genomic span, strand, gene ordinal, flags and optional CDS span/sequence/codon table; and transcript-ordered exons with genomic/cDNA spans and phase. Exon cDNA spans must be contiguous in transcript order and match genomic exon lengths. Models with transcript insertions/deletions relative to the genome are not accepted. Any source catalog satisfying this contract may be used; several named models may coexist.

### Transcript sequence

Use the 11-column CDS-only form or 13 columns ending in complete transcript-oriented pre_cds_sequence/post_cds_sequence BLOBs. CDS-only models lack transcript flanks; short-tail input is unsupported. Complete flanks enable start/stop predicates for length-changing edits crossing a CDS start or end. See the [model input contract](https://github.com/RGenomicsETL/duckhts/blob/main/design/duckvep.md#versioned-model-input-contract) for exact layouts.

### Optional queries

mature_mirna_query supplies transcript_index UINTEGER, mature_mirna_start/end UBIGINT ordered by transcript/start. peptide_edit_query supplies transcript_index UINTEGER, protein_position UINTEGER, alternate_amino_acid VARCHAR ordered uniquely by transcript/position. interval_feature_query supplies regulation_feature_index UINTEGER, seq_region UINTEGER, feature_start/end UINTEGER and feature_kind UTINYINT ordered by region/start/index; kind 1 is RegulatoryFeature and 2 is MotifFeature. Other identifiers and metadata stay relational. Explicit NULL optional queries/reference equal omission.

### Reference

reference_fasta enables reference-validated VEP-style 3-prime HGVS shifting and requires seq_region UINTEGER, sequence_length UBIGINT and seq_region_name VARCHAR in the region query. An existing FASTA index is required; names and exact lengths are verified without modifying it. The model pins validated FASTA/.fai/optional .gzi sources. Workers own independent mutable handles and reject detectable source replacement or in-place mutation around cache-miss fetches. See [reference ownership](https://github.com/RGenomicsETL/duckhts/blob/main/design/duckvep.md#ownership) for platform-specific source retention.

### Coverage

The default partial model returns unresolved when no transcript is loaded, not intergenic. transcript_coverage_complete := TRUE requires sequence_length UBIGINT and permits supported intergenic results after coordinate checks.

### Model paths

Ordinals address exact supplied paths. Apply declared contig-synonym, PAR, patch, haplotype, assembly_exception or circular-origin mappings before loading or constructing events; the loader performs none.

### Examples

```sql
SELECT loaded FROM duckvep_model_load('grch38', 'SELECT seq_region, sequence_length, seq_region_name FROM grch38_regions ORDER BY seq_region', 'SELECT transcript_index, seq_region, transcript_start, transcript_end, strand, gene_index, transcript_flags, cds_start, cds_end, cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence FROM grch38_transcripts ORDER BY seq_region, transcript_start, transcript_index', 'SELECT transcript_index, exon.exon_start, exon.exon_end, exon.exon_cdna_start, exon.exon_cdna_end, exon.phase, exon.end_phase FROM grch38_transcripts, LATERAL unnest(exons) AS u(exon) ORDER BY transcript_index, exon.exon_cdna_start', mature_mirna_query := 'SELECT transcript_index, region.mature_mirna_start, region.mature_mirna_end FROM grch38_transcripts, LATERAL unnest(mature_mirna_regions) AS u(region) ORDER BY transcript_index, region.mature_mirna_start', peptide_edit_query := 'SELECT transcript_index, edit.protein_position, edit.alternate_amino_acid FROM grch38_transcripts, LATERAL unnest(peptide_edits) AS u(edit) ORDER BY transcript_index, edit.protein_position', interval_feature_query := 'SELECT regulation_feature_index, seq_region, feature_start, feature_end, feature_kind FROM grch38_regulation ORDER BY seq_region, feature_start, regulation_feature_index', reference_fasta := '/data/Homo_sapiens.GRCh38.dna.primary_assembly.fa', transcript_coverage_complete := TRUE);
```

## duckvep_model_drop

Remove a named resident DuckVEP consequence model and release its transcript and regulation-feature interval indexes, sequences, and cached worker state. Returns FALSE when the name is absent or the model is in use by an annotation vector.

Signature:

```sql
duckvep_model_drop(name)
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT duckvep_model_drop('grch38');
```

## duckvep_allele_geometry

Separate uploaded, VEP-116 feature and minimized-edit geometry for one literal biallelic allele.

Signature:

```sql
duckvep_allele_geometry(position, reference, alternate)
```

Returns:

```
STRUCT(kind_code UTINYINT, interbase BOOLEAN, anchor_side_code UTINYINT, raw_start0 UBIGINT, raw_end0 UBIGINT, feature_start0 UBIGINT, feature_end0 UBIGINT, edit_start0 UBIGINT, edit_end0 UBIGINT, insertion_boundary0 UBIGINT, reference_difference_offset USMALLINT, reference_difference_length USMALLINT, alternate_difference_offset USMALLINT, alternate_difference_length USMALLINT)
```

### Input

position is positive and one-based. REF/ALT are nonempty A/C/G/T/N strings of at most 65,535 bases and must differ after ASCII case folding.

### Output

kind_code: 0 SNV, 1 insertion, 2 deletion, 3 length-changing replacement, 4 multi-base substitution. anchor_side_code: 0 none, 1 left, 2 right. Coordinates are non-NULL zero-based half-open intervals; insertion_boundary0 is non-NULL only for interbase insertions. Difference offsets/lengths index uploaded REF/ALT. Pure insertions have equal feature/edit start and end.

### Usage

Use feature_start0/feature_end0 for VEP feature overlap and edit_start0/edit_end0 for affected-reference providers, with an explicit insertion-flank policy when interbase is TRUE. Uploaded REF span is not the affected interval. This does not split multiallelic records or left-align: use duckhts_bcftools_norm(..., split_multiallelic := TRUE) when normalized exact keys are required.

### Examples

```sql
SELECT (duckvep_allele_geometry(100, 'A', 'ATG')).*;
```

## duckvep_transcript_projection

Project independent literal alleles and existing DuckVEP annotations into typed, unshifted VEP-116 transcript display fields.

Signature:

```sql
duckvep_transcript_projection(events_table, annotations_table, transcripts_table)
```

Returns:

```
TABLE(event_index UBIGINT, transcript_index UINTEGER, output_allele VARCHAR, interbase BOOLEAN, cdna_start BIGINT, cdna_end BIGINT, cds_start BIGINT, cds_end BIGINT, protein_start BIGINT, protein_end BIGINT, exon_first UINTEGER, exon_last UINTEGER, exon_total UINTEGER, intron_first UINTEGER, intron_last UINTEGER, intron_total UINTEGER, transcript_distance BIGINT, cds_start_nf BOOLEAN, cds_end_nf BOOLEAN, reference_amino_acids VARCHAR, alternate_amino_acids VARCHAR, reference_codons VARCHAR, alternate_codons VARCHAR)
```

### Input

events_table supplies unique event_index UBIGINT, seq_region UINTEGER, position UBIGINT and reference/alternate VARCHAR. annotations_table supplies event_index, nullable transcript_index and consequence_mask from duckvep_annotate(). transcripts_table is the validated model-loading relation: seq_region, transcript spans/strand/flags, CDS spans/sequence, codon_table, complete post_cds_sequence, transcript-ordered nested exons and peptide_edits from duckvep_ensembl_transcripts().

### Validation

Missing events/transcripts, duplicate event/model keys and source/model region disagreement error. Preserve one row per input annotation, including duplicates and transcript-free rows; output order is not guaranteed.

### Coordinates

output_allele uses VEP feature spelling, with '-' for empty ALT. Nullable one-based cDNA/CDS/protein start/end fields retain independently missing endpoints; insertion ranges use ascending flanks with explicit interbase. Exon/intron ordinals follow transcript orientation; totals describe the model even without overlap. Distance is the shortest feature-endpoint to transcript-endpoint distance, present only for upstream/downstream consequences.

### Sequence

Codon strings preserve VEP changed-base casing and variable-length indel forms. Amino-acid fields separately contain reference/alternate strings, including multiple residues. Withheld sequence yields NULL codon/amino-acid fields; '-' means an empty codon. Translation uses the prepared codon table and reference peptide edits.

### Scope

This adds display fields, not consequences, phased/structural projection or HGVS shifting. Join source/model metadata by event_index/transcript_index. CDS quality flags cds_start_nf/cds_end_nf are not canonical flags.

### Examples

```sql
SELECT p.* FROM duckvep_transcript_projection('events', 'annotations', 'grch38_transcripts') p;
```

## duckvep_repeat_alleles

Prepare bounded literal reference and alternate alleles from exact ordered repeat descriptions.

Signature:

```sql
duckvep_repeat_alleles(reference_components, alternate_components, sequence_exact, max_allele_bases := 5000)
```

Returns:

```
STRUCT(reference VARCHAR, alternate VARCHAR, reference_length UBIGINT, alternate_length UBIGINT, length_change BIGINT, length_direction VARCHAR, status VARCHAR)
```

### Input

Each component list contains unit VARCHAR and numeric count fields in sequence order. Required sequence_exact asserts that both complete descriptions, including interruptions, represent the event's literal alleles. CNV:TR summary metadata alone does not establish exactness.

### Output

Exact complete inputs return both literal sequences, their base lengths, signed alternate-minus-reference length_change and length_direction GAIN, LOSS or NEUTRAL. The direction describes expanded base length; it is not a copy-number inference. Empty lists and zero counts produce empty, not NULL, alleles.

### Status

FALSE exactness returns summary_only with all fact fields NULL. With TRUE, missing lists/elements/units/counts return incomplete_input; fractional counts return nonintegral_count without rounding. One unavailable allele withholds the complete event fact.

### Limits

Nonempty IUPAC DNA units preserve case and ambiguity. Invalid units and negative or nonfinite counts error even in summary mode. max_allele_bases is a nonnegative integer no greater than 2147483647 and applies separately to each complete allele; exhaustion errors before either sequence is returned.

### Scope

No RN/RUS/RUC/RB parsing, reference-count inference, FASTA validation, genomic anchoring, normalization or VEP expansion-policy certification is performed. Retain event and ALT ordinals, raw counts, confidence and exactness beside the fact. Use ok literal outputs with annotation or haplotype paths; keep summaries for structural consumers. ok certifies expansion, not the caller's exactness or resolution of ambiguous bases.

### Examples

```sql
SELECT (duckvep_repeat_alleles([{unit: 'CAG', count: 10}], [{unit: 'CAG', count: 5}, {unit: 'CAT', count: 1}, {unit: 'CAG', count: 5}], true, max_allele_bases := 100)).*;
```

## duckvep_breakend_geometry

Parse one raw VCF 4.5 breakend ALT into mate coordinates, orientation and retained replacement sequence.

Signature:

```sql
duckvep_breakend_geometry(alternate)
```

Returns:

```
STRUCT(mate_chrom VARCHAR, mate_position UBIGINT, local_join_after BOOLEAN, mate_extends_right BOOLEAN, replacement_sequence VARCHAR)
```

### Input

Call once per ALT after multiallelic expansion. Paired forms t[p[, t]p], ]p]t and [p[t retain exact mate names, one-based positions and case-preserved replacement t. Split mate names at the last colon; no aliases or reference lookup are applied.

### Orientation

local_join_after is TRUE when brackets follow t; mate_extends_right is TRUE for opening brackets. The mate piece is reverse-complemented when these flags differ. Replacement t includes retained local bases, not necessarily inserted-only sequence. Single forms t. and .t retain replacement/local orientation with NULL mate fields.

### Validation

NULL/non-breakend ALT, including '.', '*' and symbolic alleles, returns NULL. Malformed bracket or leading/trailing-dot forms error. Decimal mate positions from zero through UBIGINT maximum are checked; zero can denote a virtual telomeric mate.

### Scope

Parsing does not make single, telomeric or unknown-contig breakends annotatable: duckvep_annotate requires paired positive model-addressable UINTEGER coordinates after an explicit mate-name join. Retain source ID, MATEID, EVENT, confidence, REF, raw ALT and orientation/payload beside annotations. ALT alone cannot establish reciprocal identity, phase or fusion HGVS. No reciprocal-record merging or rearranged-sequence reconstruction is performed.

### Examples

```sql
SELECT (duckvep_breakend_geometry('C[2:321682[')).*;
```

## duckvep_haplotypes

Replay literal phased CDS/protein paths with carriers, source contributors, coding blocks, aligned differences and optional protein HGVS.

Signature:

```sql
duckvep_haplotypes(calls_query, model_name, phase_policy := 'strict', input_mode := 'alt_events', hgvs := FALSE, max_active_events := 16384, max_active_transcripts := 4096, max_active_carriers := 65536, max_active_prefixes := 262144, max_active_projections := 262144, max_allele_bytes := 8388608, max_leaf_events := 4096, max_leaf_edits := 65536, max_sequence_bases := 1048576, max_ploidy := 64, max_phase_sets := 1024, max_alignment_cells := 16777216, max_leaf_differences := 65536, max_hgvs_operations := 65536, max_hgvs_bytes := 1048576, max_hgvs_reference_bytes := 262144, workspace_limit := 268435456)
```

Returns:

```
TABLE(transcript_index UINTEGER, cds VARCHAR, protein VARCHAR, sequence_flags UINTEGER, evidence_flags UTINYINT, projection_status VARCHAR, sequence_status VARCHAR, edit_count UBIGINT, carrier_count UINTEGER, carriers STRUCT(sample_index UINTEGER, phase_set BIGINT, haplotype_lane USMALLINT, ploidy USMALLINT)[], contributors STRUCT(event_index UBIGINT, seq_region UINTEGER, position UBIGINT, reference VARCHAR, alternate VARCHAR, evidence_flags UTINYINT, projection_status VARCHAR)[], coding_blocks STRUCT(cds_start UINTEGER, reference VARCHAR, alternate VARCHAR, alt_start0 UBIGINT, length_change BIGINT, sequence_flags UINTEGER, coding_status VARCHAR, local_consequence_mask UBIGINT, after_first_stop BOOLEAN, event_indices UBIGINT[])[], cds_differences STRUCT(ref_start0 UBIGINT, alt_start0 UBIGINT, reference VARCHAR, alternate VARCHAR, alignment_start0 UBIGINT)[], protein_differences STRUCT(ref_start0 UBIGINT, alt_start0 UBIGINT, reference VARCHAR, alternate VARCHAR, alignment_start0 UBIGINT)[], stop_in_displaced_frame BOOLEAN, hgvsp VARCHAR, hgvsp_status VARCHAR, nominal_length_diff BIGINT)
```

### Input

calls_query supplies explicit event/transcript/sample candidates for model_name. alt_events uses decoded GT/PS; source_records requires vep116_compat, original gt and complete ordered alternates. Retain reference-only context and equal-position file order in event_index: source-buffer order is planned before genotype filtering. DuckDB sorts input; result/carrier-list order is not guaranteed. See the [phased input contract](https://github.com/RGenomicsETL/duckhts/blob/main/design/duckvep.md#phased-edits) for schemas and coordinates.

### Raw replay

Duplicate identity includes the complete ordered ALT list; shadowed_duplicate contributors execute no replacement. Overlapping records replay ordered full spans, retain overwritten provenance and produce net blocks; local SO/displaced-frame facts are NULL with unsupported_ordered_replacements. Missing-slot replay and validated coding/noncoding-span omissions are conditional. Raw N/U/lowercase ALTs within ACGTUN/acgtun retain source_allele_skipped, conditional evidence and zero edits; unsupported symbols/dashes remain invalid. Other valid sources still replay. Omitted sources retain source_unmapped; model/REF failures withhold sequence. Raw contributors add alt_index: 0 REF, positive ALT, NULL undefined.

### Output

Decoded missing/unphased paths withhold sequence. nominal_length_diff sums projected ALT-minus-REF lengths before current-CDS-end clipping; it is NULL without CDS and may differ from rebuilt length change. Blocks expose local coding status, masks decoded by duckvep_so_terms() and first-stop position facts, not whole-haplotype consequences or protein rescue. Intronic context preserves provenance without changing CDS. Samples without retained exonic genotypes use curated reference peptide; retained exonic calls use mutation translation even with zero edits.

### HGVS

hgvs := TRUE requests a VEP-116-derived predicted protein suffix with VEP local-peptide/unknown-residue presentation. hgvsp_status='ok' means supported computation, not independent HGVS-rule certification; other statuses retain NULL text and sequence/provenance. Phase policy controls assignment, not nomenclature. Conditional/overlapping paths, missing peptide data, unsupported contexts and unrepresentable ends remain explicit. One-original-ALT paths use independent-event VEP HGVS, including shifting and absent results; MNV islands retain that identity. Required unavailable FASTA returns missing_reference. Prepared references keep their own residues/length; reference-only stop-marker loss supplies no extension and insertions need reference flanks. Normalization does not rewrite coding blocks or source IDs.

### Limits

max_hgvs_operations bounds working/final operations; max_hgvs_bytes excludes NUL; max_hgvs_reference_bytes includes NUL and FASTA line-ending scratch. Sequence/edit/allele storage derives from max_sequence_bases, max_leaf_edits and max_allele_bytes. Named capacities error without growth, approximation or dropped rows. workspace_limit includes DuckVEP buffers, not DuckDB input/sort/output or HTSlib handle/transport storage. Disabled HGVS allocates no HGVS workspace/reference handle and returns not_requested.

### Scope

Join identifiers through model ordinals. Whole-haplotype SO/IMPACT/NMD, DNA HGVS, complete protein HGVS, structural composition and splice prediction are unfinished. Preparation reads committed objects on the retained connection, not caller TEMP objects or uncommitted writes; nested/concurrent preparation errors as busy.

### Examples

```sql
SELECT * FROM duckvep_haplotypes('SELECT * FROM prepared_transcript_calls', 'grch38');
```

```sql
SELECT transcript_index, cds, protein, unnest(carriers) AS carrier FROM duckvep_haplotypes('SELECT * FROM prepared_transcript_calls', 'grch38', max_ploidy := 8);
```

```sql
SELECT * FROM duckvep_haplotypes('SELECT * FROM prepared_source_calls', 'grch38', input_mode := 'source_records', phase_policy := 'vep116_compat');
```

## duckvep_phase_call

Assign decoded GT/PS allele slots to haplotype lanes under strict or pinned VEP-116 phase policy.

Signature:

```sql
duckvep_phase_call(alleles, phase_before, phase_set := NULL, phase_policy := 'strict')
```

Returns:

```
STRUCT(input_slot USMALLINT, allele_index INTEGER, haplotype_lane USMALLINT, ploidy USMALLINT, phase_set BIGINT, phase_scope VARCHAR, status VARCHAR)[]
```

### Input

Consume read_geno alleles/phase_before before ALT expansion and retain source record/sample/model identity. Ploidy is 1..65535. NULL GT returns NULL; empty GT or negative called indices error. NULL phase lists/flags mean unavailable phase; otherwise list lengths must match. NULL policy means strict; unknown policy errors.

### Output

One entry preserves each one-based input_slot, source allele_index (0 REF, positive ALT, NULL missing) and ploidy. Scope is local to sample/chromosome. phase_set retains signed BIGINT identity; NULL default differs from zero. Unresolved heterozygous slots have NULL haplotype_lane, unresolved scope and unphased status. Missing alleles stay missing even when their lane is known.

### Strict

VCF 4.4 indicators qualify the following allele: 0|1/2 and /0|1/2 leave slots 1/3 unphased; a later pipe does not phase its preceding allele. Explicitly phased lanes are retained; other slots resolve only when permutation cannot change assignment. Homozygous called/haploid GTs have all_phase_sets scope and apply across sets, not a separate NULL-PS bucket.

### VEP compatibility

vep116_compat ignores separators/PS and assigns allele_slot lanes after removing missing entries; missing entries remain explicit with NULL lane. Decoded ploidy is preserved. This does not emulate raw mixed/prefixed-separator parsing or VEP container ploidy inference.

### Scope

Do not discard missing/unresolved entries. This consumes declared GT/PS, not PSL/PSO/PID or inferred phase; it does not validate VCF PS schema or source ALT bounds, build haplotypes or emit SO/HGVS. Other phase encodings need an explicit adapter.

### Examples

```sql
SELECT unnest(duckvep_phase_call([0, 1], [true, true], phase_set := 10)) AS allele;
```

```sql
SELECT record_index, c.sample_index, unnest(duckvep_phase_call(c.alleles, c.phase_before, phase_set := c.phase_set)) AS allele FROM (SELECT record_index, unnest(calls) AS c FROM read_geno('geno_calls.bcf')) ORDER BY record_index, sample_index, allele.input_slot;
```

## duckvep_annotate

Annotate independent literal alleles, exact typed structural events and paired breakends against a resident VEP-116-compatible model.

Signature:

```sql
duckvep_annotate(events_table, model_name, hgvs := FALSE, upstream_distance := 5000, downstream_distance := 5000, rich := FALSE)
```

Returns:

```
table
```

### Input

events_table supplies event_index UBIGINT, seq_region UINTEGER, positive one-based position UBIGINT, reference/alternate VARCHAR, end_position UBIGINT, structural_type/copy_change VARCHAR and mate_seq_region UINTEGER/mate_position UBIGINT. One row represents one ALT. Keep genotype, phase, confidence, raw ALT, orientation and other provenance in the source relation; join by event_index.

### Geometry

Literal alleles are small variants; explicit structural type, supported symbolic ALT or spans without literal alleles select structural events; a complete mate pair selects BND. Incomplete/contradictory geometry errors, including symbolic/type disagreement. For VCF span SVs use POS+1 through END; insertion position is the interbase site after POS. NULL copy_change implies LOSS for DEL and GAIN for DUP/TDUP.

### Ordering

Supply globally nondecreasing model-region/coordinate order; the macro does not sort. Region ordinals address exact model paths, without contig aliases or PAR/patch/alternate-haplotype/circular-origin projection. Events are independent; GT/PS does not combine them. Zero disables either directional transcript window.

### Output

Compact fields include SO/region masks, IMPACT/status/reason codes, one-based cDNA/CDS/protein coordinates, amino-acid bytes, NMD prediction/escape codes and overlap-object ordinals. protein_position identifies the affected residue/codon. rich := TRUE adds readable consequence/impact/region/amino-acid/NMD/object/audit labels while retaining compact fields; FALSE/NULL leaves labels NULL. Audit status distinguishes unresolved computation from absent consequence and is not an Ensembl CSQ field.

### HGVS

hgvs := TRUE adds independent-event HGVSc/HGVSn/HGVSp for small variants and composes with rich := TRUE. transcript_hgvs and protein_hgvs are accession-free c./n./p. bodies; join versioned transcript/protein identifiers from retained model metadata when a complete serialized VEP-style field is required. Structural/BND HGVS remains NULL.

### gVCF alleles

Expanded `<NON_REF>`, bare `*` and `.` emit no annotation. The distinct `<*>` overlap allele receives coding_sequence_variant or corresponding start/stop-retained terms and no HGVS. Its literal three-character length can cause ablation when long REF contains a feature; filter `<*>` after ALT expansion if unwanted. Record-level END does not turn a literal mixed-gVCF ALT into an SV.

### Metadata

Decode masks with duckvep_so_terms(); join transcript/gene/regulation ordinals to model metadata. This is not the full serialized VEP CSQ vocabulary: exon/intron display, codon strings, distance, CANONICAL/APPRIS/TSL, existing variation and supplementary predictors require their model/annotation relations. See [independent-event execution](https://github.com/RGenomicsETL/duckhts/blob/main/design/duckvep.md#independent-variant-execution) and [structural compatibility](https://github.com/RGenomicsETL/duckhts/blob/main/ERRATA.md#model-inputs-and-structural-events).

### Examples

```sql
SELECT * FROM duckvep_annotate('canonical_events', 'grch38');
```

```sql
SELECT event_index, transcript_index, transcript_hgvs, protein_hgvs FROM duckvep_annotate('canonical_events', 'grch38', hgvs := TRUE, upstream_distance := 5000, downstream_distance := 0);
```

```sql
SELECT event_index, transcript_index, consequence, impact, protein_position, nmd_prediction FROM duckvep_annotate('canonical_events', 'grch38', rich := TRUE);
```

## duckvep_so_terms

Return VEP-116 Sequence Ontology terms, consequence-mask bits, impact, severity rank and evaluator tier.

Signature:

```sql
duckvep_so_terms()
```

Returns:

```
table
```

### Usage

Each row includes bit index and single-bit mask. Filter/aggregate compact annotations before expanding selected terms with (annotation.consequence_mask & term.consequence_mask) <> 0. Metadata comes from the same pinned class model as the native engine.

### Examples

```sql
SELECT a.event_index, t.consequence FROM duckvep_annotate('canonical_events', 'grch38') a JOIN duckvep_so_terms() t ON (a.consequence_mask & t.consequence_mask) <> 0;
```

## read_bcf

Read VCF/BCF with header-typed INFO/FORMAT, typed CSQ/ANN/BCSQ annotations, sample selection and optional tidy sample rows.

Signature:

```sql
read_bcf(path, region := NULL, index_path := NULL, tidy_format := FALSE, additional_csq_column_types := NULL, scan_mode := 'auto', decompression_threads := 0, decode_error_policy := 'null', samples := NULL)
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
read_geno(path, region := NULL, index_path := NULL, samples := NULL, non_reference_only := FALSE, scan_mode := 'auto', decompression_threads := 0, decode_error_policy := 'null', format_fields := NULL, raw_gt := FALSE, include_filter := FALSE)
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

Full scans stream in physical input order regardless of index availability, with one scan worker and separate HTSlib decompression threads. record_index is zero-based and scan-local, not a persistent file locator. Indexed region unions restart at zero in HTSlib iterator order; overlapping regions visit a physical record once while preserving distinct duplicates. Use ORDER BY record_index when order matters. Sequential mode rejects regions.

### FORMAT

format_fields adds header-typed members in calls.format, such as AD/DP/GQ; NULL/[] preserves default schema. Header Type/Number determine scalar/list shape. Missing elements retain ordinals, vector-end padding terminates lists, absent fields are NULL and explicit numeric missing lists are [NULL], without inferred A/R/G padding. Missing GT does not hide selected values unless the call is filtered. Unknown, empty, NULL or case-insensitively duplicate selections error. Lookup uses exact header spelling: GT/PS cannot be reselected; declared lowercase gt/ps are distinct tags.

### FILTER

include_filter := TRUE adds the same physical record's FILTER as VARCHAR[]. PASS is [PASS], an unapplied dot is NULL, and named failing filters retain header order. FALSE preserves the default output schema.

### Errors

Header/type/cardinality mismatches follow decode_error_policy. Scalar cardinality excludes vector-end padding, including after sample selection, but counts missing elements. Multiple scalar elements error or make that FORMAT tag NULL for every selected sample under null/warn; list cardinality is not inferred. Physical read/allocation failures always error.

### Raw GT

raw_gt := TRUE appends calls.raw_gt VARCHAR after any format struct, preserving exact original VCF separators, leading phase markers and allele spelling. Absent GT is NULL; literal '.' stays text. FALSE preserves the default schema. BCF rejects TRUE even when calls is unprojected. Text follows original-header sample selection and scan ordinals, including region unions; it does not emulate a downstream parser.

### Examples

```sql
SELECT record_index, CHROM, POS, unnest(calls) AS call FROM read_geno('geno_calls.bcf', non_reference_only := TRUE) ORDER BY record_index;
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

scan_mode='sequential' streams rather than using indexed count/parallel paths and rejects region. NULL/empty region means no filter; empty comma-separated items and malformed known-contig intervals error. Unknown contigs follow HTSlib's skip policy.

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

Read BED3-BED12 interval files with canonical typed columns and optional tabix-backed region filtering. scan_mode := 'sequential' forces full-file streaming/counting instead of index-backed count paths and is incompatible with region.

Signature:

```sql
read_bed(path, region := NULL, index_path := NULL, scan_mode := 'auto')
```

Returns:

```
table
```

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

## duckhts_cgranges_from_query

Execute a SQL query on an extension-owned DuckDB connection, append its interval rows into a session-scoped cgranges registry entry, and leave the populated index ready for explicit finalization with duckhts_cgranges_index(...).

Signature:

```sql
duckhts_cgranges_from_query(name, query, chrom_col, start_col, end_col[, label_col])
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT duckhts_cgranges_from_query('targets_idx', 'SELECT chrom, start, "end", name FROM targets', 'chrom', 'start', 'end', 'name');
```

## duckhts_cgranges_from_table

Reserved convenience constructor for bulk cgranges population from a table name. The current implementation is intentionally deferred and directs callers to duckhts_cgranges_from_query(...).

Signature:

```sql
duckhts_cgranges_from_table(name, table_name, chrom_col, start_col, end_col[, label_col])
```

Returns:

```
BOOLEAN
```

### Examples

```sql
SELECT duckhts_cgranges_from_table('targets_idx', 'targets', 'chrom', 'start', 'end', 'name');
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

### Examples

```sql
SELECT q.*, hit.interval_chrom, hit.interval_start, hit.interval_end, hit.label FROM read_bed('queries.bed') AS q CROSS JOIN UNNEST(duckhts_cgranges_overlaps_list('targets_idx', q.chrom, q.start, q."end")) AS u(hit);
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

## duckhts_cgranges_overlaps_bulk

Run a SQL query that yields overlap probes, stream those rows through a finalized session-scoped cgranges registry entry, and return one row per matching indexed interval. The probe query runs on the extension-owned helper connection, so it must reference regular tables/views rather than connection-local temp tables. When query_row_id_col is omitted, query_row_id defaults to the 1-based probe row ordinal.

Signature:

```sql
duckhts_cgranges_overlaps_bulk(name, query, chrom_col, start_col, end_col, mode := 'overlap', query_row_id_col := NULL)
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM duckhts_cgranges_overlaps_bulk('targets_idx', 'SELECT probe_id, chrom, start, "end" FROM probes', 'chrom', 'start', 'end', query_row_id_col := 'probe_id');
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

## duckhts_somalier_import_sites

Import an already selected Somalier sites VCF/BCF as one canonical panel and population-frequency relation.

Signature:

```sql
duckhts_somalier_import_sites(path, assembly_name, max_sites := 1000000)
```

Returns:

```
table
```

### Orientation

Each autosomal biallelic SNV is oriented into lexical allele_a/allele_b order. INFO/AF is the source ALT frequency and is transformed with REF/ALT, so population_b_af always describes allele_b. Dense zero-based site_index follows Somalier v0.3.4's lexical region then position order.

### Input

INFO/AF must declare Number=A,Type=Float. Retained records must be distinct uppercase A/C/G/T biallelic SNVs with one finite AF in [0,1]. Exact Somalier v0.3.4 X/Y aliases are excluded. Existing FILTER values are retained as provenance because this imports an already selected sites file; it does not apply selection policy.

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

Symbolic alleles are unavailable in this scope. filter_policy='pass_or_unapplied' accepts PASS and dot but makes named failures unavailable; 'include_all' uses their AD; 'error' rejects them. source_method, count_scope, status, FILTER, record/sample ordinals and matched allele slots preserve extraction provenance. Full scans are sequential and retain physical source-record identity; Parquet is optional downstream storage, not an input requirement.

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

Exactly one of panel_table or panel_parquet supplies the canonical typed six-column panel. The panel is materialized and validated once, then worker_count partitions it into at most 64 DuckDB scan jobs. Each worker-local job owns one alignment handle, index, multi-region iterator, pileup, reference handle and bounded overlap workspace; no mutable htslib or faidx state is shared. DuckDB's connection thread setting limits concurrent jobs. decompression_threads separately controls htslib workers per alignment handle. Parallel output order is unspecified; use ORDER BY site_index when order matters.

### Evidence

A and B are exact uppercase panel bases; other counts every other observed base. Valid uncovered sites emit measured 0/0/0. Missing reference contigs or positions, reference mismatches, alignment-header contig absence and positions beyond a declared alignment contig emit NULL counts with a named unavailable status. Deletions and reference skips are not observed bases. MAPQ, base quality, required/excluded flags and overlap suppression are explicit per-call settings.

### Limits and transport

panel_table must name a committed table or view visible to the database; caller-local TEMP objects and uncommitted changes are not visible to panel preparation. One retained-connection preparation slot is shared by concurrent calls: nested or concurrent preparation returns a busy error and callers may retry. panel_parquet reads an ordinary Parquet file directly. max_sites bounds the shared panel. max_depth, max_overlap_qnames and max_region_bytes bound each active scan job; concurrent workspace can therefore grow with worker_count. remote_block_bytes, remote_cache_bytes and reference_cache_bytes apply per worker-owned handle, and zero disables that cache override. Explicit non-colocated BAM/CRAM and FASTA indexes are honored; CRAM keeps the alignment file's original @SQ lengths and does not create a default FASTA-index sidecar. overlap_policy='hileup_v0.1.0' pins encounter-order mate suppression only, while 'none' counts every retained observation.

### Examples

```sql
CREATE TABLE allele_counts AS SELECT * FROM duckhts_somalier_bam_counts('sample.bam', 'fingerprint_panel', 'sample-1', 'reference.fa');
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

panel_table has assembly, dense zero-based site_index, region, positive one-based position, allele_a and allele_b. Assembly and region are each limited to 1,024 bytes before hashing. One nonempty assembly, unique physical region/position and uppercase single-base A/C/G/T alleles in lexical A < B order are required. The pinned Somalier v0.3.4 X/Y aliases are rejected; other aliases cannot be biologically classified from a region string. The calculation assumes diploid three-state genotypes. The digest commits to this domain and every ordered site, independent of physical row order.

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

frequency_table must cover every panel site exactly once with the same assembly, ordinal, coordinate and A/B orientation plus a finite population_b_af in [0,1]. The digest commits to the panel identity and every ordered frequency value.

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

evidence_table contains sample_id, the six panel identity columns, and nullable a, b and other counts. Each sample must contain every panel ordinal exactly once. Count channels already follow the panel's canonical lexical A/B order; changed coordinates or orientation error. Counts are never inferred from GT or DP. All-NULL tuples are unavailable, all-zero tuples are measured zero depth, and partial NULL tuples error.

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

Both sketches must have identical panel digests, site counts and classification settings. Mask words are borrowed directly; comparison allocates no pair-sized workspace. SQL chooses the requested pair relation and output ordering.

### Examples

```sql
SELECT unnest(duckhts_somalier_relatedness(a.sketch, b.sketch, 1000000)) FROM sample_sketches a JOIN sample_sketches b ON a.sketch.sample_id < b.sketch.sample_id;
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

Count evidence and population_b_af must match the panel's full ordered site identity. CHARR uses measured counts, its own homozygous-like binomial eligibility and the pinned 4% other-read filter; relatedness genotype masks are insufficient.

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

pairs_table contains distinct receiver_id and anchor_id rows. The anchor supplies the receiver's expected uncontaminated homozygous genotype; it is not assumed to identify the contaminating donor. Reversing a pair is a different fit.

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

Read GFF annotations with optional raw scalar and parsed list/pair attributes, strict GFF3 validation and indexed region selection.

Signature:

```sql
read_gff(path, header_names := NULL, header := FALSE, column_types := NULL, auto_detect := FALSE, attributes_map := FALSE, attributes_list := FALSE, attributes_pairs := FALSE, strict := FALSE, region := NULL, index_path := NULL, scan_mode := 'auto')
```

Returns:

```
table
```

### Scanning

Comma-separated indexed regions emit each row once across overlaps. scan_mode='sequential' streams/counts instead of using indexed count paths and rejects region. NULL/empty region means no filter; empty comma-separated items and malformed known-contig intervals error. Unknown contigs follow HTSlib's skip policy.

### Examples

```sql
SELECT seqname, feature, start, "end" FROM read_gff('gff_file.gff.gz') LIMIT 5;
```

## read_gtf

Read GTF annotations with optional raw scalar and parsed list/pair attributes and indexed region selection.

Signature:

```sql
read_gtf(path, header_names := NULL, header := FALSE, column_types := NULL, auto_detect := FALSE, attributes_map := FALSE, attributes_list := FALSE, attributes_pairs := FALSE, region := NULL, index_path := NULL, scan_mode := 'auto')
```

Returns:

```
table
```

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
read_genbank(path, attributes_map := FALSE)
```

Returns:

```
table
```

### Mapping

seqname is VERSION, else ACCESSION, else the LOCUS name; a segment on a remote accession (ACC.1:5..40) is reported under that accession. source is 'GenBank'. join()/order() give one row per segment in biological order, complement(...) sets strand '-', and the GFF3 phase of each CDS segment is carried from /codon_start across segments (absent means 1). A span whose end precedes its start on a circular record wraps the origin into two segments, and a between site n^m is the zero-length site at n. The record-level source feature is dropped, and /translation is omitted as redundant with ORIGIN.

### Attributes

Synthesized GFF3 keys ID, Name and Parent accompany the original qualifiers. ID is gene-<locus_tag> for genes and <key>-<n> otherwise, with n the feature's 0-based position in the file. Parent links a feature to the gene sharing its /locus_tag (or /gene) anywhere in the record. Repeated qualifiers become one key with comma-joined values, valueless qualifiers read 'true', and ; = & , % are percent-encoded. attributes_map := TRUE adds the same pairs as a MAP.

### Errors

Records stream one at a time, so memory follows the largest record. A record without a terminating //, a FEATURES table with no sequence section, a malformed location, an unsupported location form (one-of, gap, bond, nested join/order, a.b), an origin-spanning span on a linear record, or a /codon_start outside 1..3 is an error naming the feature and line. Table layout and these rules follow BioPython's GenBank scanner, and test/scripts/genbank_oracle_test.py diffs the reader against it.

### Examples

```sql
SELECT seqname, feature, start, "end" FROM read_genbank('phix174.gb') LIMIT 5;
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
bam_bin_counts(path, bin_width, chrom := NULL, reference := NULL, index_path := NULL, mapq := 0, require_flags := 0, exclude_flags := 0, rmdup := 'none', stats := NULL)
```

Returns:

```
table
```

### Examples

```sql
SELECT * FROM bam_bin_counts('fixture_mixed.cram', 5000, reference := 'fixture_ref.fa', rmdup := 'streaming', stats := 'gc,mq');
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
