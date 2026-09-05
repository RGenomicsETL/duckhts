
<!-- README.md is generated from README.Rmd using [duckknit](https://github.com/rundel/duckknit). Please edit that file. -->

# DuckHTS

[![CRAN
Status](https://www.r-pkg.org/badges/version/Rduckhts)](https://cran.r-project.org/package=Rduckhts)
[![R-universe
version](https://RGenomicsETL.r-universe.dev/Rduckhts/badges/version)](https://RGenomicsETL.r-universe.dev/Rduckhts)

Read VCF, BCF, BAM, CRAM, FASTA, FASTQ, BigWig, GTF, GFF, BED, and
tabix-indexed files directly in [DuckDB](https://duckdb.org/). DuckHTS
uses [htslib](https://github.com/samtools/htslib) for HTS formats and
provides SQL functions for consequence annotation, intervals, coverage,
sequence operations, compression, indexing, and export.

## Functions

<details>
<summary>
Show generated function catalog
</summary>

## Extension Function Catalog

This section is generated from `functions.yaml`.

### Diagnostics

| Function                             | Kind         | Returns                | R helper                              | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|--------------------------------------|--------------|------------------------|---------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `duckhts_htslib_version`             | scalar       | VARCHAR                | `rduckhts_htslib_version`             | Return the runtime version reported by the htslib library loaded with DuckHTS. Rduckhts uses this value to reject a downstream linking receipt whose source/header version does not match the loaded library.                                                                                                                                                                                                                                                                                                                                               |
| `duckhts_htslib_features`            | scalar       | UINTEGER               |                                       | Return the htslib runtime feature bitfield reported by hts_features(). Use duckhts_htslib_feature_string() for the corresponding build description.                                                                                                                                                                                                                                                                                                                                                                                                         |
| `duckhts_htslib_feature_string`      | scalar       | VARCHAR                |                                       | Return htslib’s runtime build-feature description, including configured transports, compression libraries, compiler, and build flags. DuckHTS snapshots it once while loading the extension so parallel SQL calls read immutable text.                                                                                                                                                                                                                                                                                                                      |
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

| Function                              | Kind        | Returns                                                                                                                                                                                                                                                                                                                                                                                       | R helper | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|---------------------------------------|-------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `duckvep_ensembl_regions`             | table_macro | table                                                                                                                                                                                                                                                                                                                                                                                         |          | Match a tiled reference sequence to one Ensembl core assembly and assign dense DuckVEP sequence-region ordinals. core_schema must contain the Ensembl seq_region and coord_system tables. reference_chunks_table must contain chrom, zero-based start, half-open end, and seq columns, normally materialized from fasta_nuc(…, include_seq := TRUE). Chunks must be contiguous from zero, their sequence lengths must match their intervals, and every FASTA contig must match exactly one same-length Ensembl region. Regions absent from the supplied FASTA are deliberately excluded. The supplied reference therefore defines the exact modeled path set: a primary-assembly FASTA produces a primary-path model, while alternate haplotypes or patches require their named sequence paths in the reference chunks. This macro does not apply Ensembl assembly_exception projections, merge X/Y PAR paths, or mark circular sequence regions; ordinary mitochondrial coordinates can be modeled, but an interval that wraps from the end of a circular sequence to its start is not accepted.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| `duckvep_ensembl_transcripts`         | table_macro | table                                                                                                                                                                                                                                                                                                                                                                                         |          | Build the validated DuckVEP transcript relation from Ensembl core tables and matching tiled FASTA sequence. This supplied builder compiles the VEP-116 Ensembl core selection; it does not build VEP –refseq or –merged transcript sets. Before dense ordinals are assigned, transcripts must be current and have a non-empty stable ID, and artifact-biotype or readthrough_tra transcripts are excluded. The result contains the exact resident-loader columns, stable/source identifiers, biotypes, versioned RefSeq accessions from MANE Select and MANE Plus Clinical attributes when present, a nested ranked-exon projection, mature miRNA cDNA attributes projected into genomic exon segments, and supported Translation SeqEdits. MANE accessions remain cold DuckDB columns while only selection flags enter the resident C model. GENCODE Basic and Primary membership is retained in transcript_flags for late filtering; this macro does not prefilter the model to either set. Prepared sequence contains the phase-adjusted CDS plus the complete transcript-oriented pre-CDS and post-CDS spliced sequence; the legacy post_cds_bases field remains as the first three post-CDS bases. Codon-table IDs come from the Ensembl seq_region_attrib relation and default to table 1 exactly as VEP 116 does; every BioPerl/VEP-supported NCBI table ID is accepted and invalid or conflicting source attributes reject the import. Single-residue initial_met, \_selenocysteine, amino_acid_sub, and \_stop_codon_rt edits are retained as a sparse reference-peptide overlay. Transcript-level sequence corrections, including inserted or deleted transcript bases relative to the genome, and other Translation SeqEdit shapes fail closed by withholding sequence with an explicit reason. Transcript rows remain specific to the exact selected sequence path, including separate X/Y PAR and alternate-haplotype stable IDs. The macro does not collapse alt_allele equivalence groups, use IS_PAR metadata to merge consequence rows, or lift transcripts across a circular sequence origin.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| `duckvep_ensembl_regulation_features` | table_macro | table                                                                                                                                                                                                                                                                                                                                                                                         |          | Build the VEP-116-compatible interval relation for Ensembl RegulatoryFeature and MotifFeature rows on the sequence regions selected for a DuckVEP model. funcgen_schema must contain regulatory_feature, feature_type, and motif_feature; regions_table is the output of duckvep_ensembl_regions(…). Matching VEP’s database and cache source, epigenetically_modified_region (EMAR) RegulatoryFeature rows are excluded before dense ordinals are assigned. The result maps source sequence-region IDs to model ordinals, validates one-based inclusive coordinates and source feature types, and preserves stable/source IDs, feature metadata, binding-matrix IDs, regulatory-build IDs, and motif scores as cold relation columns. feature_kind is the compact resident code: 1 regulatory region or 2 transcription-factor binding site. Pass the five hot columns regulation_feature_index, seq_region, feature_start, feature_end, and feature_kind to duckvep_model_load(interval_feature_query := …); the remaining columns stay relational for late projection.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `duckvep_model_receipt`               | table_macro | table                                                                                                                                                                                                                                                                                                                                                                                         |          | Create one deterministic receipt for prepared DuckVEP region and transcript relations plus an optional prepared regulation-feature relation. The original eight positional arguments remain valid. When regulation_features_table is supplied by name, the receipt also validates its dense ordinals, region agreement, interval geometry, and feature kinds; counts regulatory regions and motif features; and includes their five resident columns in the model hash. It returns the declared provenance as source_name, source_version, assembly, source_manifest_sha256, reference_sha256, and transcript_filter; counts transcript, sequence, mature-miRNA, and peptide-edit content; and hashes every semantic resident-model field in stable ordinal order. source_manifest_sha256 should identify a canonical manifest containing every exact core and funcgen input used. The declared transcript_filter or source manifest should also record whether the reference contains only primary paths or includes alternate/patch paths and whether assembly_exception preprocessing was applied. The receipt contains no clock time, so rebuilding identical inputs yields the same model hash.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| `duckvep_model_load`                  | table       | table(loaded BOOLEAN)                                                                                                                                                                                                                                                                                                                                                                         |          | Load one immutable consequence model into the current DuckDB database instance under a caller-chosen name. The resident kernel is source-neutral: callers may prepare any transcript catalog that satisfies this model contract. Three required query strings read committed, non-temporary DuckDB relations: sorted UINTEGER sequence-region ordinals; dense transcript rows with genomic span, strand, gene ordinal, flags, optional CDS span/sequence/codon table; and transcript-ordered exon rows with genomic and cDNA spans plus phase. Each exon cDNA span must be contiguous in transcript order and the same length as its genomic exon span. Transcript models that insert or delete bases relative to the genome therefore need a future explicit mapping and are not accepted by this interface. The optional mature_mirna_query reads transcript_index UINTEGER, mature_mirna_start UBIGINT, and mature_mirna_end UBIGINT ordered by transcript and start. The optional peptide_edit_query reads transcript_index UINTEGER, protein_position UINTEGER, and alternate_amino_acid VARCHAR ordered uniquely by transcript and protein position. The optional interval_feature_query reads regulation_feature_index UINTEGER, seq_region UINTEGER, feature_start UINTEGER, feature_end UINTEGER, and feature_kind UTINYINT ordered by region, start, and index; kind 1 is RegulatoryFeature and kind 2 is MotifFeature. Loading keeps those hot columns in a compact SoA and builds a separate cgranges seed index, while identifiers and funcgen metadata remain ordinary DuckDB columns. The transcript projection accepts 11 base columns, 12 columns with the legacy post_cds_bases BLOB of up to three bases, or 13 columns ending in complete transcript-oriented pre_cds_sequence and post_cds_sequence BLOBs. Complete flanks enable VEP start/stop predicates for length-changing edits that cross a CDS start or end. reference_fasta enables reference-validated, VEP-style 3-prime HGVS shifting. When supplied, the sequence-region query must return seq_region UINTEGER, sequence_length UBIGINT, and seq_region_name VARCHAR; loading requires an existing FASTA index and verifies every declared name and exact length without creating or modifying the index. The named model pins open read descriptors for the validated FASTA, .fai, and optional .gzi. Linux workers reopen those descriptors through /proc/self/fd, Windows retains the resolved source under deny-write sharing, and other POSIX workers use independent resolved-source handles with identity checks rather than shared /dev/fd seek state. Each annotation worker owns its mutable faidx handle, reuses contained sorted reference requests, and rejects detectable replacement or in-place source mutation around a cache-miss fetch. Sequence-region ordinals are exact model paths: loading performs no contig-synonym, PAR, patch, haplotype, assembly_exception, or circular-origin projection. Callers must apply any declared path mapping before loading or constructing events. Explicit NULL for any optional query or reference parameter is equivalent to omission. The default is a partial transcript model: a variant with no loaded transcript is unresolved, not intergenic. Setting transcript_coverage_complete := TRUE requires sequence_length UBIGINT and permits supported intergenic results after coordinate bounds checks. Loading validates and narrows the model once and returns one TRUE row. Several named models may coexist.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| `duckvep_model_drop`                  | scalar      | BOOLEAN                                                                                                                                                                                                                                                                                                                                                                                       |          | Remove a named resident DuckVEP consequence model and release its transcript and regulation-feature interval indexes, sequences, and cached worker state. Returns FALSE when the name is absent or the model is in use by an annotation vector.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| `duckvep_allele_geometry`             | scalar      | STRUCT(kind_code UTINYINT, interbase BOOLEAN, anchor_side_code UTINYINT, raw_start0 UBIGINT, raw_end0 UBIGINT, feature_start0 UBIGINT, feature_end0 UBIGINT, edit_start0 UBIGINT, edit_end0 UBIGINT, insertion_boundary0 UBIGINT, reference_difference_offset USMALLINT, reference_difference_length USMALLINT, alternate_difference_offset USMALLINT, alternate_difference_length USMALLINT) |          | Interpret one literal biallelic allele with the same lossless geometry authority used by DuckVEP. position is positive and one-based; REF and ALT are non-empty A/C/G/T/N strings of at most 65,535 bases and must differ after ASCII case folding. kind_code is stable: 0 SNV, 1 insertion, 2 deletion, 3 length-changing replacement, and 4 multi-base substitution. anchor_side_code is 0 for no retained anchor, 1 for a left anchor, and 2 for a right anchor. All coordinate fields are non-null zero-based half-open intervals except insertion_boundary0, which is non-null only for an interbase insertion. Difference offsets and lengths are non-null indices into the uploaded REF and ALT. The result separates the uploaded REF span, VEP-116 VariationFeature span, and fully minimized semantic reference edit. A pure insertion has equal feature/edit start and end plus insertion_boundary0. This function does not split multiallelic records or left-align alleles: run duckhts_bcftools_norm(…, split_multiallelic := TRUE) first when canonical exact keys or cross-source joins require normalization. Use feature_start0/feature_end0 for VEP consequence feature overlap, edit_start0/edit_end0 for affected-reference supplementary providers, and an explicit insertion-flank policy when interbase is TRUE; the uploaded raw REF span is not a biological affected interval.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| `duckvep_annotate`                    | table       | table                                                                                                                                                                                                                                                                                                                                                                                         |          | Annotate a narrow canonical event relation against one resident DuckVEP model through one fixed-schema, pipe-friendly surface. events_table names a coordinate-ordered relation with event_index UBIGINT, seq_region UINTEGER, positive one-based position, reference VARCHAR, alternate VARCHAR, end_position UBIGINT, structural_type VARCHAR, copy_change VARCHAR, mate_seq_region UINTEGER, and mate_position UBIGINT. One row represents one ALT allele. Keep genotype, phase-set, confidence-interval, raw-ALT, orientation, and other provenance in the source relation; join selected annotations back by event_index instead of multiplying wide payloads by transcript count. Literal reference/alternate rows are small variants; an explicit structural type, supported symbolic ALT, or span without literal alleles is structural; a complete mate coordinate pair is a breakend. Following VEP 116, expanded gVCF <NON_REF>, bare *, and . ALT rows emit no annotation. The distinct \<*\> catch-all allele remains an overlap allele: coding positions receive coding_sequence_variant, start/stop codon positions receive the corresponding retained consequence, and HGVS fields remain NULL. VEP 116 also compares the literal three-character \<*\> spelling with REF length, so a long-REF catch-all allele that completely contains a feature can receive an ablation term; exclude \<*\> rows after ALT expansion if catch-all consequences are not wanted. A literal ALT from a mixed gVCF record remains a small variant even when the source relation retains record-level INFO/END; keep that END as provenance rather than treating it as an SV type. Contradictory or incomplete geometry fails instead of being guessed, including disagreement between supported symbolic ALT and explicit structural_type values. For VCF span SVs, prepare position as POS + 1 after removing the left anchor; insertion position is the interbase site after POS. DEL implies LOSS and DUP/TDUP imply GAIN when copy_change is NULL. Filtering the globally coordinate-ordered event stream preserves the order required by each native transcript and regulation-feature sweep; the macro does not hide a sort. The event seq_region is an exact ordinal in the loaded model; annotation does not strip contig prefixes or project PAR, patch, alternate-haplotype, or circular-origin coordinates. X/Y PAR and alternate-path transcript rows therefore remain path-specific exactly as prepared. Each event is annotated independently: genotype and phase-set columns are not combined into haplotypes. The stable compact contract includes SO/region masks, IMPACT/status/reason codes, one-based cDNA/CDS/protein positions, amino-acid bytes, NMD prediction/escape codes, and overlap-object ordinals. protein_position is the affected one-based protein residue/codon ordinal used by codon-level evidence such as PS1/PM5. rich := TRUE adds Consequence, IMPACT, region, amino-acid, NMD, overlap-object, duckvep_status, and duckvep_reason text from the same native pass while retaining every compact field; FALSE or explicit NULL leaves those presentation columns NULL. hgvs := TRUE fuses independent-event HGVSc/HGVSn/HGVSp into the small-variant path without repeating candidate discovery; it composes with rich := TRUE in one pass; structural and breakend HGVS remain NULL. DuckHTS audit status/reason fields distinguish an unresolved computation from a genuinely absent consequence and are not Ensembl CSQ fields. Zero disables either directional transcript window. Decode consequence_mask only after filtering by joining duckvep_so_terms(); join transcript_index, gene_index, or regulation_feature_index to the prepared model relations for stable identifiers, biotype, MANE/GENCODE/CCDS flags, and other presentation metadata. This typed rich projection does not claim the complete serialized VEP CSQ vocabulary: exon/intron ordinals, codon strings, distance text, CANONICAL/APPRIS/TSL, existing variation, and supplementary predictors require their declared model or annotation relations. |
| `duckvep_so_terms`                    | table       | table                                                                                                                                                                                                                                                                                                                                                                                         |          | Return the generated Sequence Ontology metadata behind DuckVEP’s stable consequence-mask bits: bit index, single-bit mask, term, impact code and label, VEP severity rank, and evaluator tier. Filter and aggregate compact annotation rows numerically first, then expand only selected mask bits with (annotation.consequence_mask & term.consequence_mask) \<\> 0. The relation is generated from the same pinned VEP-116 class-model metadata as the native consequence engine; it is not a second hand-maintained term inventory.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |

### Readers

| Function                 | Kind         | Returns | R helper                                                                                                                                                               | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
|--------------------------|--------------|---------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `read_bcf`               | table        | table   | `rduckhts_bcf`                                                                                                                                                         | Read VCF and BCF variant data with typed INFO, FORMAT, typed CSQ/ANN/BCSQ subfields, optional tidy sample output, optional bcftools-style CSQ type overrides, explicit scan_mode control (‘auto’ or full-file ‘sequential’), optional htslib decompression worker threads via decompression_threads (default 0 for single-threaded reads), and decode_error_policy (‘null’, ‘warn’, or ‘error’) for corrupt BCF header-vs-payload decode mismatches. The reader retains worker-owned INFO/FORMAT/GT buffers and record-owned annotation across output chunks; tidy sample rows repeat the complete record annotation. Comma-separated indexed regions use htslib’s native multi-region iterator and emit a record once when requested regions overlap.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `read_bcf_appender`      | table        | table   |                                                                                                                                                                        | Experimental side-effecting benchmark helper that reads BCF/VCF with htslib and appends a narrow scalar tidy table through DuckDB’s chunk appender API. The target table is created in the default schema with columns CHROM, POS, REF, ALT, SAMPLE_ID, and FORMAT_GT; set include_file_offset := TRUE to add FILE_OFFSET as a file-order token. Comma-separated regions use htslib’s native deduplicating multi-region iterators. region_threads controls at most 16 worker handles claiming primary-contig jobs iteratively: every contig keeps its complete requested interval set in one native iterator, so repeated, overlapping, or record-spanning intervals produce the same target schema and row multiset for every thread count. Header contig count does not create worker threads or open handles, and a request limited to one contig uses one iterator. Parallel row arrival is unordered; use ORDER BY FILE_OFFSET when file order matters. decompression_threads configures htslib workers per open region-worker handle, so total helper threads can scale as region workers times decompression workers. The operation runs inside an internal transaction and surfaces malformed-record errors as DuckDB errors; decode_error_policy controls corrupt BCF FORMAT/GT header-vs-payload mismatches with ‘null’, ‘warn’, or ‘error’. target_table must be an unqualified SQL identifier. This is intended to compare appender-based materialization against read_bcf(…) query pipelines, not as a full replacement for the typed read_bcf surface. |
| `read_bam`               | table        | table   | `rduckhts_bam`                                                                                                                                                         | Read SAM, BAM, and CRAM alignments with optional typed SAMtags and auxiliary tag maps. Use sequence_encoding := ‘nt16’ to return SEQ as UTINYINT\[\], quality_representation := ‘phred’ to return QUAL as UTINYINT\[\], and cigar_representation := ‘binary’ to return packed BAM CIGAR operations as UINTEGER\[\] instead of SAM text. scan_mode := ‘sequential’ forces a full-file streaming scan instead of index-backed count/parallel paths and is incompatible with region. decompression_threads controls per-file htslib worker threads and defaults to 2; use 0 to disable worker threads.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| `read_fasta`             | table        | table   | `rduckhts_fasta`                                                                                                                                                       | Read FASTA records or indexed FASTA regions as sequence rows. Use sequence_encoding := ‘nt16’ to return SEQUENCE as UTINYINT\[\] (htslib nt16 4-bit codes) instead of VARCHAR. For bgzipped FASTA, gzi_path may point to an explicit .gzi sidecar when it is not colocated with the FASTA. scan_mode := ‘sequential’ forces full-file streaming/counting instead of index-backed count paths and is incompatible with region.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| `read_bed`               | table        | table   | `rduckhts_bed`                                                                                                                                                         | Read BED3-BED12 interval files with canonical typed columns and optional tabix-backed region filtering. scan_mode := ‘sequential’ forces full-file streaming/counting instead of index-backed count paths and is incompatible with region.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| `fasta_nuc`              | table        | table   | `rduckhts_fasta_nuc`                                                                                                                                                   | Compute bedtools nuc-style nucleotide composition for supplied BED intervals or generated fixed-width bins over a FASTA reference. For bgzipped FASTA, gzi_path may point to an explicit .gzi sidecar when it is not colocated with the FASTA.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `read_fastq`             | table        | table   | `rduckhts_fastq`                                                                                                                                                       | Read single-end, paired-end, or interleaved FASTQ files with optional legacy quality decoding. By default, FASTQ qualities are interpreted as modern Phred+33 input. Use sequence_encoding := ‘nt16’ to return SEQUENCE as UTINYINT\[\] and quality_representation := ‘phred’ to return QUALITY as UTINYINT\[\] instead of VARCHAR. input_quality_encoding accepts ‘phred33’, ‘auto’, ‘phred64’, or ‘solexa64’. scan_mode := ‘sequential’ forces raw streaming/counting instead of index-backed count paths.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| `read_bigwig`            | table        | table   | `rduckhts_bigwig`                                                                                                                                                      | Read stored BigWig signal intervals as CHROM, START0, END0, and VALUE. Coordinates are zero-based half-open. region uses htslib’s one-based inclusive syntax; comma-separated requests are merged per contig and emit each stored interval once. Full scans claim nonempty contigs across DuckDB workers, and multi-region scans claim merged ranges. Each worker owns its mutable file handle and iterator. blocks_per_iteration controls libBigWig iterator batching, not the number of DuckDB workers. Local, native remote, and browser wasm reads all use DuckHTS’s htslib hFILE transport rather than a second libcurl stack.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
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

### Quality Control

| Function           | Kind      | Returns | R helper | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
|--------------------|-----------|---------|----------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `duckhts_fastq_qc` | aggregate | STRUCT  |          | Aggregate canonical sequence and Phred+33 quality strings directly into exact read/base/Q20/Q30/Q40, nucleotide, quality-sum, and per-cycle sufficient statistics. The nested cycles list supports mean-quality, nucleotide-content, GC, and read-length curves without expanding one SQL row per base. Rows with any NULL input are ignored. Per-cycle state defaults to at most 1,048,576 cycles; pass a constant max_cycles per aggregate group to choose a larger explicit limit, up to 16,777,216. |

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
| `duckhts_contig_key`        | scalar      | VARCHAR  |                          | Return a conservative contig join key by removing one non-empty leading chr prefix case-insensitively and normalizing M/MT to MT. X and Y are uppercased; all other suffixes are preserved. This does not map numeric sex chromosomes, accessions, patches, or alternate loci.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `bcftools_liftover`         | scalar      | STRUCT   | `rduckhts_liftover`      | Row-oriented liftover kernel intended to mirror bcftools +liftover semantics as closely as possible while returning one STRUCT per input row with fields: src_chrom, src_pos, src_ref, src_alt, dest_chrom, dest_pos, dest_end, dest_ref, dest_alt, mapped, reverse_complemented, swap, reject_reason, and note. Set no_left_align := true to skip post-liftover left-alignment of lifted indels (mirrors –no-left-align in bcftools +liftover).                                                                                                                                                                                                                                                                                                                                                                                                                             |
| `duckdb_liftover`           | table_macro | table    | `rduckhts_liftover`      | DuckDB-specific wrapper over bcftools_liftover that takes either a table name or a derived-table expression plus column-name strings for chrom/pos/ref/alt and returns the lifted table. The no_left_align parameter mirrors –no-left-align in bcftools +liftover.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| `bcftools_norm_row`         | scalar      | STRUCT   |                          | Normalize one variant row with bcftools/vt-style left-alignment semantics against a FASTA reference. The alt argument may be either a comma-delimited VARCHAR or a VARCHAR\[\] list. The returned STRUCT contains pos_normed, end_pos_normed, ref_normed, alt_normed (always VARCHAR\[\]), normed (TRUE/FALSE/NULL), and norm_status. Symbolic <DEL> rows can use end_pos, and symbolic <DUP> rows can use svlen. gVCF <NON_REF>/\<\*\> reference-block alleles pass through with GVCFReferenceBlock, and mixed real-plus-gVCF-symbolic alleles normalize the real alleles while preserving symbolic alleles and caller-supplied reference-block END in site-preserving output.                                                                                                                                                                                              |
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
       min(VALUE)::DOUBLE AS minimum,
       max(VALUE)::DOUBLE AS maximum
FROM read_bigwig(
  'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw',
  region := 'chr22:20000000-20099999'
);
```

    ┌──────────────────┬─────────────────────┬──────────────────┐
    │ stored_intervals │       minimum       │     maximum      │
    │      int64       │       double        │      double      │
    ├──────────────────┼─────────────────────┼──────────────────┤
    │            96783 │ -10.786999702453613 │ 9.60200023651123 │
    └──────────────────┴─────────────────────┴──────────────────┘

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

### DuckVEP: a resident Ensembl VEP consequence engine

#### Design and data ownership

DuckVEP compiles Ensembl release relations into an immutable, shared
transcript and regulation model. Sorted variant alleles then pass
through a native C candidate sweep, transcript projection, consequence
classifier, sequence editor, NMD evaluator, and optional HGVS renderer.
Stable transcript/gene identifiers and supplementary annotations remain
ordinary DuckDB relations and are joined only after the compact
consequence rows have been filtered.

That last sentence is the point of the architecture, not an
implementation detail. DuckVEP is one fast relation producer inside
DuckDB. Anything DuckDB can read–Parquet or GeoParquet, VCF/BCF,
CSV/TSV, JSON, an attached DuckLake, a local DuckDB, HTTP/S3 objects, or
a supported external database–can be an annotation provider without a
new DuckHTS file format or C plugin. Exact allele sources use
collision-safe VariantKey joins, gene sources use stable identifiers,
and interval sources use ordinary inequality joins, bounded key
prefilters, or cgranges. DuckDB can execute two-inequality overlap
predicates with its native
[IEJoin](https://duckdb.org/2022/05/27/iejoin) implementation and can
spill large join state rather than requiring every provider to live in
the resident C model.

The implementation targets executable Ensembl VEP 116 semantics rather
than a new interpretation of the Sequence Ontology rules. The current
differential ledger covers GRCh38, GRCh37, *P. falciparum*, exact
structural variants, paired breakends, regulation/motif features, and
variant-induced NMD. See the rendered [conformance
report](benchmarks/duckvep_conformance.md), the [throughput
report](benchmarks/duckvep_throughput.md), and the implementation
contract in [design/duckvep.md](design/duckvep.md).

#### Building a receipted model

The model compiler reads ordinary DuckDB relations. Stage the exact
Ensembl core and funcgen release tables from their MySQL dump files, or
attach a read-only Ensembl MySQL database through DuckDB’s `mysql`
extension. The compiler is not a downloader: source retrieval and
dump-to-table typing are explicit staging steps, while the same
validation and model receipt apply to either transport.

``` sql
-- Reference chunks are zero-based, half-open and sequence-bearing.
CREATE TABLE grch38_reference_chunks AS
SELECT chrom, start, "end", seq
FROM fasta_nuc(
  'Homo_sapiens.GRCh38.dna.primary_assembly.fa',
  bin_width := 1048576,
  include_seq := TRUE
);

-- ensembl_core and ensembl_funcgen are schemas containing the release-116
-- tables named in design/duckvep.md.
CREATE TABLE grch38_regions AS
SELECT * FROM duckvep_ensembl_regions(
  'ensembl_core', 'grch38_reference_chunks', 'GRCh38'
);

CREATE TABLE grch38_transcripts AS
SELECT * FROM duckvep_ensembl_transcripts(
  'ensembl_core', 'grch38_reference_chunks', 'GRCh38'
);

CREATE TABLE grch38_regulation AS
SELECT * FROM duckvep_ensembl_regulation_features(
  'ensembl_funcgen', 'grch38_regions'
);

CREATE TABLE grch38_receipt AS
SELECT * FROM duckvep_model_receipt(
  'grch38_regions', 'grch38_transcripts',
  'Ensembl', '116', 'GRCh38',
  source_manifest_sha256, reference_sha256,
  'VEP 116 core transcript selection',
  regulation_features_table := 'grch38_regulation'
);
```

`duckvep_ensembl_transcripts(...)` applies the pinned VEP
core-transcript selection, reconstructs ranked exon/cDNA/CDS geometry
and strand-correct sequence, loads the sequence-region genetic code
(including mitochondrial code), prepares complete transcript flanks,
mature-miRNA segments, and supported Ensembl sequence edits, and retains
MANE/GENCODE/CCDS and stable identifiers as cold relational columns.
GRCh37 is built from the separate release-116 GRCh37 core schema and is
intentionally MANE-free. Regulation and motif features are compiled from
the matching funcgen release, not inferred from transcript overlap.

Publication is separate from compilation. Persist the prepared relations
and load their sorted hot projections into a named resident model with
`duckvep_model_load(...)`. `reference_fasta` enables reference
validation and VEP-compatible three-prime HGVS shifting;
`transcript_coverage_complete := TRUE` permits a no-transcript hit to
become a supported `intergenic_variant`.

``` sql
SELECT loaded
FROM duckvep_model_load(
  'human_116_grch38',
  'SELECT seq_region, sequence_length, seq_region_name
     FROM grch38_regions ORDER BY seq_region',
  'SELECT transcript_index, seq_region, transcript_start, transcript_end,
          strand, gene_index, transcript_flags, cds_start, cds_end,
          cds_sequence, codon_table, pre_cds_sequence, post_cds_sequence
     FROM grch38_transcripts
    ORDER BY seq_region, transcript_start, transcript_index',
  'SELECT transcript_index, exon.exon_start, exon.exon_end,
          exon.exon_cdna_start, exon.exon_cdna_end, exon.phase, exon.end_phase
     FROM grch38_transcripts,
          LATERAL unnest(exons) AS u(exon)
    ORDER BY transcript_index, exon.exon_cdna_start',
  mature_mirna_query :=
    'SELECT transcript_index, region.mature_mirna_start,
            region.mature_mirna_end
       FROM grch38_transcripts,
            LATERAL unnest(mature_mirna_regions) AS u(region)
      ORDER BY transcript_index, region.mature_mirna_start',
  peptide_edit_query :=
    'SELECT transcript_index, edit.protein_position,
            edit.alternate_amino_acid
       FROM grch38_transcripts,
            LATERAL unnest(peptide_edits) AS u(edit)
      ORDER BY transcript_index, edit.protein_position',
  interval_feature_query :=
    'SELECT regulation_feature_index, seq_region, feature_start, feature_end,
            feature_kind
       FROM grch38_regulation
      ORDER BY seq_region, feature_start, regulation_feature_index',
  reference_fasta := 'reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa',
  transcript_coverage_complete := TRUE
);
```

#### Why sorted input becomes a sweep

The resident transcript arrays are sorted by
`(seq_region, transcript_start)`. For the first allele in a DuckDB
vector, cgranges supplies a complete seed set. As later alleles advance
in coordinate order, the worker admits newly reachable transcripts from
the front of the sorted array and retires transcripts whose end plus the
requested flank is behind the allele. It classifies only that active
set:

``` text
transcripts ordered by start:  T1------T2----------T3-----T4---->
sorted allele stream:               v1  v2       v3    v4       >

seed once -> admit starts now reachable -> retire expired ends
          -> project/classify active transcripts -> advance
```

This is not a 5,000-record or 5,000-base buffer. The default 5,000 bases
are only the VEP-compatible upstream/downstream candidate distances;
zero, 10,000, 50,000, or another supported `UINTEGER` distance uses the
same algorithm. Long alleles are not clipped to that window. Transcript
and regulation/motif sweeps own different active sets because interval
features do not use transcript flanks. The immutable model is shared
read-only across workers; each worker owns its candidate arrays, exon
cursors, faidx handle, reference window, edit scratch, and output
builder. Consequently memory grows with the shared model plus bounded
per-worker workspace, not with a copy of the model per thread. The
current scalar surface reseeds at each DuckDB vector; a future stateful
table-function surface can preserve the frontier across vectors without
changing the consequence authority.

#### Canonical event relation and output choices

`duckvep_annotate(...)` consumes one coordinate-ordered row per ALT
allele. Its canonical columns cover literal small variants, exact span
SVs, and paired breakends in one relation. Genotypes, phase sets,
imprecision intervals, raw ALT, and sample metadata stay on the source
row and join back by `event_index`.

| `rich`  | `hgvs`  | Result                                                                                                    |
|:--------|:--------|:----------------------------------------------------------------------------------------------------------|
| `FALSE` | `FALSE` | Compact masks, codes, ordinals, cDNA/CDS/protein positions, amino-acid bytes and NMD codes                |
| `TRUE`  | `FALSE` | Compact fields plus Consequence, IMPACT, region, amino-acid, NMD and audit text from the same native pass |
| `FALSE` | `TRUE`  | Compact fields plus independent-event HGVSc/HGVSn/HGVSp                                                   |
| `TRUE`  | `TRUE`  | Rich consequence and HGVS from one fused candidate-discovery pass                                         |

The transcript window defaults to 5,000 bases in each direction,
matching the VEP default. It is an explicit API parameter rather than a
fixed buffer: either direction may be zero or any supported `UINTEGER`
distance. Audit columns named `duckvep_status` and `duckvep_reason`
explain unresolved computation (for example missing sequence or a
reference mismatch); they are not Ensembl CSQ fields.

#### Real human annotation and supplementary providers

The following is the production shape, using the public HG002 40x
PCR-free [DeepVariant GRCh38 WGS
callset](https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/vcf/hiseqx/wgs_pcr_free/40x/HG002.hiseqx.pcr-free.40x.deepvariant-v1.0.grch38.vcf.gz),
an Ensembl 116 GRCh38 model, dated ClinVar/ClinvArbitration,
AlphaMissense v2, gnomAD v2.1.1 gene constraint, and an Ensembl
regulatory Parquet relation. This callset contains chromosomes 1–22, X,
Y, and MT; it is an annotation workload, not a GIAB truth set. The
provider files are ordinary typed Parquet; equivalent HTTP/S3 paths or
DuckLake tables can replace them. Run the statements in one DuckDB CLI
session so `.timer on` reports model load, VCF staging, the consequence
sweep, provider joins, and final Parquet materialization separately.

##### Build the provider relations once

The `providers/*.parquet` paths below are not special DuckHTS files.
They are release-specific projections built with DuckDB from the
declared upstream artifacts:

| Provider relation                       | Upstream artifact                                                                                               |
|:----------------------------------------|:----------------------------------------------------------------------------------------------------------------|
| `clinvar_20260706_grch38_keyed.parquet` | [ClinVar GRCh38 VCF](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/), dated 2026-07-06                    |
| `clinvarbitration_grch38_keyed.parquet` | [ClinvArbitration record 16792026](https://zenodo.org/records/16792026), whose release contract declares GRCh38 |
| `alphamissense_hg38_variantkey.parquet` | [AlphaMissense v2 GRCh38](https://zenodo.org/records/8360242)                                                   |
| UCSC `hg38.phyloP100way.bw`             | [GRCh38 phyloP 100-way BigWig](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/)                   |
| `ensembl116_grch38_regulatory.parquet`  | `grch38_regulation`, compiled above from the Ensembl 116 funcgen tables                                         |
| `gnomad_v211_constraint_gene.parquet`   | [gnomAD v2.1.1 gene constraint](https://gnomad.broadinstitute.org/downloads#v2-constraint)                      |

VariantKey is an exact key over a normalized allele representation.
Normalize arbitrary VCF inputs with `duckhts_bcftools_norm(...)` or an
equivalent pinned normalizer before building either side of the join;
keep the untouched record, genotype arrays and normalization lineage in
their source relations.

ClinVar is already VCF, so `read_bcf()` supplies typed INFO fields while
ALT ordinal expansion produces the allele relation.

``` sql
COPY (
  WITH alleles AS (
    SELECT duckhts_contig_key(CHROM) AS chrom, POS::BIGINT AS pos,
           upper(REF) AS ref, upper(a.alt) AS alt,
           INFO_ALLELEID AS allele_id,
           INFO_CLNSIG AS clinical_significance
    FROM read_bcf('source/clinvar_20260706.vcf.gz', scan_mode := 'sequential')
    CROSS JOIN unnest(ALT) AS a(alt)
  ), keyed AS (
    SELECT *, variantkey(chrom, pos, ref, alt) AS vk FROM alleles
  )
  SELECT chrom, pos, ref, alt, vk,
         mod(extract_variantkey_refalt(vk), 2) = 1 AS is_hash,
         allele_id, clinical_significance
  FROM keyed WHERE vk IS NOT NULL
  ORDER BY vk
) TO 'providers/clinvar_20260706_grch38_keyed.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 122880);
```

ClinvArbitration publishes an aggregate decision table. Its assembly is
a declared GRCh38 release fact; it is not inferred from `chr`-prefixed
strings.

``` sql
COPY (
  WITH source AS (
    SELECT duckhts_contig_key(contig) AS chrom,
           try_cast(position AS BIGINT) AS pos,
           upper(reference) AS ref, upper(alternate) AS alt,
           clinical_significance,
           try_cast(gold_stars AS UTINYINT) AS gold_stars,
           try_cast(allele_id AS BIGINT) AS allele_id
    FROM read_csv_auto(
      'source/clinvarbitration_16792026.tsv',
      delim := '\t', header := TRUE, all_varchar := TRUE
    )
  ), keyed AS (
    SELECT *, variantkey(chrom, pos, ref, alt) AS vk FROM source
  )
  SELECT *, mod(extract_variantkey_refalt(vk), 2) = 1 AS is_hash
  FROM keyed WHERE vk IS NOT NULL
  ORDER BY vk
) TO 'providers/clinvarbitration_grch38_keyed.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 122880);
```

The official AlphaMissense `.tsv.gz` is BGZF-compressed. DuckHTS can
build its coordinate index directly, then use the same file for indexed
regional queries or, as below, stream all 71.7 million rows once to
build the reusable pack.

``` sql
SELECT * FROM tabix_index(
  'source/AlphaMissense_hg38.tsv.gz',
  preset := 'gff', seq_col := 1, start_col := 2, end_col := 2,
  comment_char := '#', skip_lines := 4, threads := 4
);

COPY (
  WITH keyed AS (
    SELECT variantkey(duckhts_contig_key(chrom), pos, ref, alt) AS vk,
           split_part(transcript_id, '.', 1) AS transcript_stable_id,
           try_cast(split_part(transcript_id, '.', 2) AS UINTEGER)
             AS transcript_version,
           transcript_id, protein_variant, uniprot_id,
           am_pathogenicity, am_class
    FROM read_tabix(
      'source/AlphaMissense_hg38.tsv.gz',
      header_names := [
        'chrom', 'pos', 'ref', 'alt', 'genome', 'uniprot_id',
        'transcript_id', 'protein_variant', 'am_pathogenicity', 'am_class'
      ],
      column_types := [
        'VARCHAR', 'BIGINT', 'VARCHAR', 'VARCHAR', 'VARCHAR', 'VARCHAR',
        'VARCHAR', 'VARCHAR', 'DOUBLE', 'VARCHAR'
      ],
      scan_mode := 'sequential'
    )
  )
  SELECT vk, transcript_stable_id, transcript_version,
         transcript_id, protein_variant, uniprot_id,
         am_pathogenicity, am_class
  FROM keyed WHERE vk IS NOT NULL
  ORDER BY vk, transcript_stable_id, transcript_version, protein_variant
) TO 'providers/alphamissense_hg38_variantkey.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 122880);
```

The Ensembl provider is projected from the same release-116 funcgen
relation that supplies the resident core features. Identifiers and SO
metadata remain in Parquet; only compact feature geometry enters the C
model.

``` sql
COPY (
  WITH source AS (
    SELECT f.regulation_feature_index, f.stable_id, f.feature_class,
           coalesce(f.feature_so_term, f.feature_class) AS so_term,
           r.seq_region_name AS chrom,
           f.feature_start::BIGINT - 1 AS start0,
           f.feature_end::BIGINT AS end0
    FROM grch38_regulation f
    JOIN grch38_regions r USING (seq_region)
  ), keyed AS (
    SELECT *, regionkey(chrom, start0, end0) AS rk FROM source
  )
  SELECT *, rk >> 31 AS rk_chrom_start,
         regionkey(chrom, end0, end0) >> 31 AS rk_chrom_end
  FROM keyed WHERE rk IS NOT NULL
  ORDER BY rk
) TO 'providers/ensembl116_grch38_regulatory.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 122880);
```

Gene constraint joins by stable Ensembl gene identifier and therefore
needs no genomic key.

``` sql
COPY (
  SELECT gene_id, gene AS gene_symbol, transcript,
         try_cast(pLI AS DOUBLE) AS pLI,
         try_cast(oe_lof_upper AS DOUBLE) AS oe_lof_upper
  FROM read_csv_auto(
    'source/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz',
    delim := '\t', header := TRUE, nullstr := 'NA'
  )
  WHERE gene_id IS NOT NULL
) TO 'providers/gnomad_v211_constraint_gene.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
```

The resulting files are reusable annotation-pack relations. Sorting
exact providers by `vk` and interval providers by `rk` gives Parquet
zonemaps useful ordering information; DuckDB still performs the final
collision or overlap predicate. The full provider preparation and join
measurements are in the [supplementary annotation
benchmark](benchmarks/benchmark_variantkey_join_overlap.md).

##### Session and resident model

`.timer on` makes the DuckDB CLI report every following statement
separately. Keep the thread count and memory limit fixed for the
complete run. Execute the `duckvep_model_load(...)` statement above in
this session before staging the alleles; the checked run loaded the full
Ensembl 116 model in 2.398 seconds; 644,427 transcripts.

``` sql
.timer on
SET threads = 4;
SET memory_limit = '4GB';
SET temp_directory = 'case/duckdb-tmp';
```

##### Canonical allele relation

The input VCF is a record relation, while the consequence engine accepts
an event relation with one row per ALT. `record_index` and the one-based
`alt_index` preserve the route back to the original record and its
genotype arrays. This small-variant lane admits literal nucleotide
alleles only; exact SVs and paired breakends use the same event schema
through the typed columns that are `NULL` below.

``` sql
-- Preserve the original VCF record and ALT numbering before normalization.
-- The complete sample/FORMAT columns may stay in a wider source relation and
-- join back through (record_index, alt_index).
CREATE TABLE case_records AS
SELECT row_number() OVER ()::UBIGINT AS record_index,
       CHROM AS source_chrom, duckhts_contig_key(CHROM) AS chrom,
       POS::BIGINT AS pos, REF AS ref, ALT AS alt
FROM read_bcf(
  'case/HG002.hiseqx.pcr-free.40x.deepvariant-v1.0.grch38.vcf.gz',
  scan_mode := 'sequential', decompression_threads := 0
);
```

``` sql
-- Split ALT with its original one-based ordinal, validate REF against the
-- release-matched FASTA, and left-align/minimize before any exact key is made.
CREATE TABLE case_normalized AS
SELECT record_index, alt_index, source_chrom, chrom,
       pos_normed AS pos, upper(ref_normed) AS ref,
       upper(alt_normed) AS alt, norm_status
FROM duckhts_bcftools_norm(
  'case_records',
  'reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa',
  split_multiallelic := TRUE
)
WHERE pos_normed IS NOT NULL
  AND regexp_full_match(ref_normed, '[ACGTNacgtn]+')
  AND regexp_full_match(alt_normed, '[ACGTNacgtn]+')
  AND upper(ref_normed) <> upper(alt_normed);
```

``` sql
CREATE TABLE case_alleles AS
WITH interpreted AS (
  -- This calls the same event interpreter as the consequence kernel. It is
  -- not another SQL trim implementation.
  SELECT n.*, duckvep_allele_geometry(pos, ref, alt) AS geometry
  FROM case_normalized n
), keyed AS (
  SELECT i.*, r.seq_region,
         variantkey(chrom, pos, ref, alt) AS vk,
         regionkey(chrom, geometry.feature_start0,
                   geometry.feature_end0) AS feature_rk
  FROM interpreted i
  JOIN grch38_regions r ON r.seq_region_name = i.chrom
)
SELECT row_number() OVER (
         ORDER BY seq_region, pos, record_index, alt_index
       )::UBIGINT AS event_index,
       record_index, alt_index, seq_region, chrom,
       pos::UBIGINT AS position, ref AS reference, alt AS alternate,
       vk, mod(extract_variantkey_refalt(vk), 2) = 1 AS is_hash,
       geometry.kind_code, geometry.interbase,
       geometry.raw_start0, geometry.raw_end0,
       geometry.feature_start0, geometry.feature_end0,
       geometry.edit_start0, geometry.edit_end0,
       geometry.insertion_boundary0,
       feature_rk >> 31 AS feature_rk_chrom_start,
       regionkey(chrom, geometry.feature_end0,
                 geometry.feature_end0) >> 31 AS feature_rk_chrom_end,
       NULL::UBIGINT AS end_position,
       NULL::VARCHAR AS structural_type,
       NULL::VARCHAR AS copy_change,
       NULL::UINTEGER AS mate_seq_region,
       NULL::UBIGINT AS mate_position
FROM keyed
WHERE vk IS NOT NULL AND feature_rk IS NOT NULL
-- Coordinate order is the sweep contract. Normalization can move POS, so sort
-- only after normalization and event interpretation.
ORDER BY seq_region, position, event_index;
```

Measured VCF decoding, ALT expansion, normalization, keying and
ordering: 12.200 seconds; 7,378,240 alleles staged.

##### Consequence and HGVS

The resident model contains only the hot transcript, exon, sequence and
core feature fields. Stable identifiers and supplementary provider
payloads remain in DuckDB and are joined later. `hgvs` and `rich` change
output projection; they do not cause a second candidate-discovery pass.

Do not retain the complete result of `duckvep_annotate(...)` before
writing it: this WGS emits 88 million rows. The final query below
consumes the annotator directly. In isolation, the same four-thread
rich/HGVS stream writes all rows in 24.524 seconds; 88,392,840
annotation rows.

##### Exact allele providers

ClinVar and ClinvArbitration are allele-level providers. VariantKey can
encode short alleles reversibly or store a hash for longer alleles. The
`is_hash` equality prevents mixing those representations; hashed matches
are then checked against the normalized literal contig, position, REF
and ALT so a collision cannot become clinical evidence. Lists of structs
keep each provider record intact instead of independently aggregating
fields that belong together.

``` sql
CREATE TABLE case_clinvar AS
SELECT q.event_index,
       list_distinct(list(struct_pack(
         allele_id := p.allele_id,
         clinical_significance := p.clinical_significance
       ))) AS clinvar_records
FROM case_alleles q
JOIN read_parquet('providers/clinvar_20260706_grch38_keyed.parquet') p
  ON q.vk = p.vk AND q.is_hash = p.is_hash
 AND (NOT q.is_hash OR
      (q.chrom = p.chrom AND q.position = p.pos AND
       q.reference = p.ref AND q.alternate = p.alt))
GROUP BY q.event_index;
```

Measured dated ClinVar join: 0.259 seconds; 50,749 matched alleles.

``` sql
CREATE TABLE case_clinvarbitration AS
SELECT q.event_index,
       list_distinct(list(struct_pack(
         allele_id := p.allele_id,
         classification := p.clinical_significance,
         gold_stars := p.gold_stars
       ))) AS clinvarbitration_records
FROM case_alleles q
JOIN read_parquet('providers/clinvarbitration_grch38_keyed.parquet') p
  ON q.vk = p.vk AND q.is_hash = p.is_hash
 AND (NOT q.is_hash OR
      (q.chrom = p.chrom AND q.position = p.pos AND
       q.reference = p.ref AND q.alternate = p.alt))
GROUP BY q.event_index;
```

Measured ClinvArbitration join: 0.176 seconds; 45,232 matched alleles.

AlphaMissense contains single-nucleotide substitutions, whose
VariantKeys are reversible, so hashed events are excluded rather than
subjected to a literal fallback comparison. Scores are
transcript-specific: the provider pack retains the versioned transcript
accession and the join resolves it to the same resident
`transcript_index` emitted by DuckVEP. Taking the maximum score across
all transcripts would attach one transcript’s prediction to unrelated
transcript consequence rows.

``` sql
CREATE TABLE case_alphamissense AS
SELECT q.event_index, t.transcript_index,
       list_distinct(list(struct_pack(
         transcript_id := p.transcript_id,
         protein_variant := p.protein_variant,
         uniprot_id := p.uniprot_id,
         pathogenicity := p.am_pathogenicity,
         classification := p.am_class
       ))) AS alphamissense_records
FROM case_alleles q
JOIN read_parquet('providers/alphamissense_hg38_variantkey.parquet') p
  ON q.vk = p.vk
JOIN grch38_transcripts t
  ON t.transcript_stable_id = p.transcript_stable_id
 AND t.transcript_version = p.transcript_version
WHERE NOT q.is_hash
GROUP BY q.event_index, t.transcript_index;
```

Measured AlphaMissense join: 0.523 seconds; 14,850 matched alleles.

##### Dense conservation signal

BigWig is already an indexed, compressed genomic signal store, so
conservation does not need a bespoke DuckVEP cache. Read only the case
ranges, preserve missing bases explicitly, and reduce overlaps in SQL.
Conservation over a replacement or deletion uses the minimized affected
reference span, not its VCF padding anchor. A pure insertion affects no
reference base; report its immediate left and right reference flanks
separately instead of fabricating a one-base mean. This example uses the
real UCSC phyloP 100-way track and a chromosome-22 slice; production
orchestration can construct comma-separated batches from the distinct
case coordinates and combine them with `UNION ALL BY NAME`.

``` sql
CREATE TABLE case_phylop_chr22 AS
WITH signal AS (
  SELECT duckhts_contig_key(CHROM) AS chrom,
         START0::BIGINT AS start0, END0::BIGINT AS end0,
         VALUE::DOUBLE AS value
  FROM read_bigwig(
    'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw',
    region := 'chr22:20000000-20099999'
  )
), requests AS (
  SELECT event_index, 'span' AS request_kind,
         edit_start0 AS start0, edit_end0 AS end0
  FROM case_alleles q
  WHERE q.chrom = '22'
    AND NOT interbase
    AND edit_start0 < 20099999 AND edit_end0 > 19999999
  UNION ALL
  SELECT event_index, 'left', insertion_boundary0 - 1, insertion_boundary0
  FROM case_alleles q
  WHERE q.chrom = '22' AND interbase AND insertion_boundary0 > 0
    AND insertion_boundary0 BETWEEN 20000000 AND 20099999
  UNION ALL
  SELECT event_index, 'right', insertion_boundary0, insertion_boundary0 + 1
  FROM case_alleles q
  WHERE q.chrom = '22' AND interbase
    AND insertion_boundary0 BETWEEN 19999999 AND 20099998
), overlaps AS (
  SELECT r.event_index, r.request_kind, r.end0 - r.start0 AS requested_bases,
         greatest(0, least(r.end0, s.end0) - greatest(r.start0, s.start0))
           AS observed_bases, s.value
  FROM requests r
  LEFT JOIN signal s
    ON r.start0 < s.end0 AND r.end0 > s.start0
)
SELECT event_index,
       sum(observed_bases) FILTER (WHERE request_kind = 'span')
         AS observed_bases,
       max(requested_bases) FILTER (WHERE request_kind = 'span')
         AS requested_bases,
       sum(observed_bases * value) FILTER (WHERE request_kind = 'span') /
         nullif(sum(observed_bases) FILTER (WHERE request_kind = 'span'), 0)
         AS phyloP_mean,
       min(value) FILTER (WHERE request_kind = 'span') AS phyloP_min,
       max(value) FILTER (WHERE request_kind = 'span') AS phyloP_max,
       max(value) FILTER (
         WHERE request_kind = 'left' AND observed_bases = 1
       ) AS insertion_left_phyloP,
       max(value) FILTER (
         WHERE request_kind = 'right' AND observed_bases = 1
       ) AS insertion_right_phyloP
FROM overlaps
GROUP BY event_index;
```

The overlap width weights run-length intervals correctly. `NULL` means
the requested span or flank had no observed phyloP value; it is not
silently converted to zero. Repeated cohort workloads may explicitly
materialize release-labelled, coordinate-sorted Parquet or DuckLake
tiles, but BigWig remains a valid provider and precision authority
without such conversion.

##### Interval provider

For RegionKey-supported contigs, shifting a RegionKey right by 31 bits
produces an exact ordered `(chromosome code, start)` value. Constructing
the same value at the half-open end gives an exact
`(chromosome code, end)` value. The two inequalities below therefore
express both chromosome identity and half-open overlap, and DuckDB plans
them as [IEJoin](https://duckdb.org/2022/05/27/iejoin).

This example deliberately uses DuckVEP’s VEP-feature coordinates. For an
insertion they are an empty interval at the interbase site, so the
strict range inequalities require both adjacent reference positions to
lie inside the provider interval, matching VEP’s core-feature test. A
provider whose contract is “affected reference bases” should instead use
`edit_start0/edit_end0` and an explicit insertion policy, as the phyloP
query does above. The uploaded `raw_start0/raw_end0` span is provenance,
not an overlap default.

Do not add `q.chrom = p.chrom` to this plan. DuckDB then chooses a
chromosome hash join and tests the range predicates across enormous
same-chromosome candidate sets. On this workload that physical-plan
change was more than two orders of magnitude slower. Use `EXPLAIN` to
confirm `IE_JOIN`. For accessions, patches, or arbitrary species contigs
outside RegionKey’s encoding domain, use
`duckhts_cgranges_overlaps_bulk(...)` instead; it preserves literal
contig identity and has its own measured memory/speed tradeoff.

``` sql
CREATE TABLE case_interval AS
SELECT q.event_index,
       list_distinct(list(struct_pack(
         regulation_feature_index := p.regulation_feature_index,
         stable_id := p.stable_id,
         so_term := p.so_term
       ))) AS interval_records
FROM case_alleles q
JOIN read_parquet('providers/ensembl116_grch38_regulatory.parquet') p
  ON q.feature_rk_chrom_start < p.rk_chrom_end
 AND q.feature_rk_chrom_end > p.rk_chrom_start
GROUP BY q.event_index;
```

Measured Ensembl regulatory interval join: 0.770 seconds; 414,813
matched alleles. The complete cgranges build, bulk query, aggregation,
and write takes 1.144 seconds; 414,813 matched alleles.

##### Streamed result

Transcript and regulation-feature indices are compact model-local
ordinals. Resolve them to stable identifiers only after annotation, then
attach gene-level and event-level payloads. The interval query above
demonstrates and measures a general relational overlap provider; it is
not needed to rediscover the same Ensembl core feature already selected
by the native sweep. Core feature metadata joins directly by
`regulation_feature_index`, which cannot attach a nearby but different
feature. This keeps provider strings out of the shared C model and lets
DuckDB project or omit them normally.

The query streams annotation into the joins and uses one Parquet writer
per DuckDB worker. `event_index` is the stable identity and ordering
key; forcing a global `ORDER BY` or a single output file serializes the
writer and requires a large sort/materialization. Query the resulting
Parquet dataset with an explicit `ORDER BY` only when a consumer
actually requires ordered rows.

``` sql
COPY (
  SELECT q.chrom, q.position, q.reference, q.alternate,
         a.consequence, a.impact, a.cdna_position, a.cds_position,
         a.protein_position, a.reference_amino_acid,
         a.alternate_amino_acid, a.nmd_prediction,
         a.transcript_hgvs, a.protein_hgvs,
         t.transcript_stable_id, t.gene_stable_id, t.mane_select_refseq,
         rf.stable_id AS regulation_stable_id,
         rf.feature_so_term AS regulation_so_term,
         cv.clinvar_records, ca.clinvarbitration_records,
         am.alphamissense_records,
         gc.pLI, gc.oe_lof_upper,
         a.overlap_object, a.duckvep_status, a.duckvep_reason
  FROM duckvep_annotate(
    'case_alleles', 'human_116_grch38',
    hgvs := TRUE, rich := TRUE,
    upstream_distance := 5000, downstream_distance := 5000
  ) a
  JOIN case_alleles q USING (event_index)
  LEFT JOIN grch38_transcripts t USING (transcript_index)
  LEFT JOIN grch38_regulation rf USING (regulation_feature_index)
  LEFT JOIN read_parquet('providers/gnomad_v211_constraint_gene.parquet') gc
    ON gc.gene_id = t.gene_stable_id
  LEFT JOIN case_clinvar cv USING (event_index)
  LEFT JOIN case_clinvarbitration ca USING (event_index)
  LEFT JOIN case_alphamissense am
    ON am.event_index = a.event_index
   AND am.transcript_index = a.transcript_index
) TO 'case/HG002.deepvariant-v1.0.duckvep'
  (FORMAT PARQUET, COMPRESSION ZSTD, PER_THREAD_OUTPUT TRUE);
```

Measured model load, fused consequence/HGVS, all late joins, and
four-writer ZSTD Parquet output: 28.841 seconds; 88,392,840 annotation
rows written.

##### Whole-genome measurement

The checked real-data integration run retained 7,378,240 alleles from
chromosomes 1–22, X, Y, and MT and emitted 88,392,840
consequence/core-feature rows. With a 4 GB DuckDB memory limit, four
unpinned workers, warm input pages, and four Parquet writers, the
model-load plus annotation/composition stream completed in 28.841
seconds; 88,392,840 annotation rows written. GNU `time -v` recorded 5.31
GiB peak process RSS. That high-water mark includes the transient
full-model load; the immutable C model is allocated outside DuckDB’s
buffer-manager limit. The checked full-row fingerprint matches the
reference output across all 88,392,840 rows.

| measured stage                              | rows       | seconds | peak RSS (GiB) | output bytes |
|:--------------------------------------------|:-----------|--------:|:---------------|:-------------|
| RegionKey IEJoin + interval result          | 414,813    |   0.770 | 1.52           | 4,262,629    |
| cgranges build/query + interval result      | 414,813    |   1.144 | 0.79           | 4,451,147    |
| rich + HGVS + core regulation, four writers | 88,392,840 |  24.524 | 4.98           | 274,046,172  |
| rich + HGVS + all providers, four writers   | 88,392,840 |  28.841 | 5.31           | 285,798,217  |

The exact-source benchmark separately measures the complete
4.0-million-allele GIAB probe against ClinVar, ClinvArbitration,
AlphaMissense, REVEL, gene constraint, and interval providers, including
one-, four-, and twelve-provider plans, Parquet row-group pruning,
materialized output, and peak RSS. See the [supplementary annotation
benchmark](benchmarks/benchmark_variantkey_join_overlap.md).

#### Bundled offline lifecycle model

The deterministic fixture below remains deliberately small so every
README render can execute it offline. It checks the public
load/annotate/drop lifecycle; it is not the human performance or
biological example.

``` sql
CREATE OR REPLACE TABLE readme_duckvep_regions AS
SELECT 1::UINTEGER AS seq_region, 50000::UBIGINT AS sequence_length,
       '11'::VARCHAR AS seq_region_name;
```

``` sql
CREATE OR REPLACE TABLE readme_duckvep_transcripts AS
SELECT 0::UINTEGER AS transcript_index, 1::UINTEGER AS seq_region,
       100::UBIGINT AS transcript_start, 150::UBIGINT AS transcript_end,
       1::TINYINT AS strand, 0::UINTEGER AS gene_index,
       0::UBIGINT AS transcript_flags, NULL::UBIGINT AS cds_start,
       NULL::UBIGINT AS cds_end, NULL::BLOB AS cds_sequence,
       NULL::UTINYINT AS codon_table, NULL::BLOB AS pre_cds_sequence,
       NULL::BLOB AS post_cds_sequence;
```

``` sql
CREATE OR REPLACE TABLE readme_duckvep_exons AS
SELECT 0::UINTEGER AS transcript_index, 100::UBIGINT AS exon_start,
       150::UBIGINT AS exon_end, 1::UBIGINT AS exon_cdna_start,
       51::UBIGINT AS exon_cdna_end, -1::TINYINT AS phase,
       -1::TINYINT AS end_phase;
```

``` sql
SELECT loaded
FROM duckvep_model_load(
  'readme',
  'SELECT * FROM readme_duckvep_regions ORDER BY seq_region',
  'SELECT * FROM readme_duckvep_transcripts ORDER BY seq_region, transcript_start, transcript_index',
  'SELECT * FROM readme_duckvep_exons ORDER BY transcript_index, exon_cdna_start',
  reference_fasta := 'test/data/fixture_ref.fa'
);
```

    ┌─────────┐
    │ loaded  │
    │ boolean │
    ├─────────┤
    │ true    │
    └─────────┘

``` sql
CREATE OR REPLACE TABLE readme_duckvep_events AS
SELECT * FROM (VALUES
  (1::UBIGINT, 1::UINTEGER, 124::UBIGINT, 'A'::VARCHAR, 'G'::VARCHAR,
   NULL::UBIGINT, NULL::VARCHAR, NULL::VARCHAR,
   NULL::UINTEGER, NULL::UBIGINT)
) AS e(event_index, seq_region, position, reference, alternate,
       end_position, structural_type, copy_change,
       mate_seq_region, mate_position);
```

``` sql
SELECT a.event_index, a.transcript_index, a.consequence, a.impact,
       a.cdna_position, a.transcript_hgvs, a.protein_hgvs,
       a.duckvep_status
FROM duckvep_annotate(
  'readme_duckvep_events', 'readme', hgvs := true,
  upstream_distance := 0, downstream_distance := 0, rich := true
) AS a
ORDER BY a.event_index, a.transcript_index;
```

    ┌─────────────┬──────────────────┬────────────────────────────────────┬──────────┬───────────────┬─────────────────┬──────────────┬────────────────┐
    │ event_index │ transcript_index │            consequence             │  impact  │ cdna_position │ transcript_hgvs │ protein_hgvs │ duckvep_status │
    │   uint64    │      uint32      │              varchar               │ varchar  │    uint32     │     varchar     │   varchar    │    varchar     │
    ├─────────────┼──────────────────┼────────────────────────────────────┼──────────┼───────────────┼─────────────────┼──────────────┼────────────────┤
    │           1 │                0 │ non_coding_transcript_exon_variant │ MODIFIER │          NULL │ n.25A>G         │ NULL         │ supported      │
    └─────────────┴──────────────────┴────────────────────────────────────┴──────────┴───────────────┴─────────────────┴──────────────┴────────────────┘

``` sql
SELECT duckvep_model_drop('readme') AS dropped;
```

    ┌─────────┐
    │ dropped │
    │ boolean │
    ├─────────┤
    │ true    │
    └─────────┘

#### Supplementary annotation composition

DuckHTS does not define a closed supplementary-annotation plugin API
because DuckDB already supplies the larger interface: a relation. For
high-cardinality pipelines, filter compact integer masks/codes first and
decode selected SO bits with `duckvep_so_terms()` afterward. Join
`transcript_index`, `gene_index`, and `regulation_feature_index` to the
prepared model relations for stable IDs and cold metadata. Exact or
positional resources such as ClinVar, population frequencies, calibrated
pathogenicity scores, and splice scores use collision-safe allele keys;
protein/domain and regulatory tracks use half-open interval overlap;
constraint, PanelApp, OMIM, and phenotype resources use stable
gene/disease IDs.

Those providers may be local Parquet, partitioned S3 Parquet, DuckLake
snapshots, attached DuckDB databases, tabix-indexed files, or any other
source for which DuckDB has a scanner. The optimizer can push columns
and predicates into scans, prune sorted Parquet row groups, plan several
equality joins together, execute two-inequality overlaps with
[IEJoin](https://duckdb.org/2022/05/27/iejoin), and spill join blocks
when required. cgranges remains useful for a repeatedly queried
immutable interval registry; it is an execution choice, not a provider
format. This makes an annotation-pack release a set of versioned
relations and provenance, while a new case or cohort is only the delta
that must be scanned, annotated, and joined.

#### Testing and conformance infrastructure

The test system is deliberately several independent authorities rather
than one large collection of expected strings:

1.  Pure-C unit and property oracles compare the optimized sweep,
    projection, edit, translation, consequence, HGVS, SV/BND, and
    workspace paths with slower reference implementations. ASan and
    UBSan execute the same properties.
2.  Mass statistical campaigns use a declared seed and explicit
    rare-state strata: strand, exon rank, splice offset,
    CDS/UTR/start/stop geometry, phase, length change and modulo three,
    incomplete codons, NMD thresholds, HGVS shift limits, structural
    containment, BND orientation, and transcript-window distances.
    Coverage counters are part of the result. A passing campaign with a
    required zero-count state fails, and every counterexample is
    retained and minimized into a fixed witness.
3.  The generated finite rule program is exhaustively checked over its
    Boolean consequence-predicate quotient. A transition manifest then
    names the wider geometry/state paths and the fixed, statistical, or
    executable evidence that currently covers each one. This is a
    compact proof obligation, not a claim that sampling proves an
    unobserved state.
4.  Executable differentials run the pinned VEP 116/BioPerl environment
    on the identical allele relation and full-outer-join unique
    allele/object pairs. Consequence, NMD, HGVSc, HGVSn, HGVSp,
    unresolved rows, and missing/extra objects remain separate fatal
    metrics; there is no release option that permits HGVS discordance.
5.  Public corpora search different state distributions: complete/dense
    ClinVar, whole GIAB, selected GRCh37, non-human codon tables,
    regulation/motif, exact SV, paired BND, and HPRC/pangenome-derived
    long alleles. A sparse WGS average cannot replace a dense
    transcript/HGVS or generated structural campaign.

The repository mirrors the relevant Ensembl VEP and Ensembl Variation
self-test cases with exact source lineage. The optional `{targets}`
project orchestrates coarse campaign invalidation, branching, retry, and
retained artifacts; it does not replace the biological comparison logic.
Perl/BioPerl plus micromamba provide the pinned oracle environment, and
[`blit`](https://github.com/WangLabCSU/blit) is the non-interactive
process runner. Cache-backed campaigns create one acquisition receipt
and cheaply revalidate the complete installed-file inventory before
reuse, instead of hashing tens of gigabytes on every invocation. The
[corpus workflow](design/duckvep_corpus_workflow.md) defines evidence
identity, denominators, and release transitions; the [conformance
README](test/duckvep/conformance/README.md) owns commands and
external-campaign details.

#### Measured performance and FastVEP comparison

The resident-engine benchmark uses the complete GIAB HG002 v4.2.1 input
and the complete Ensembl 116 GRCh38 transcript and core-feature model.

| measure                            | value       |
|:-----------------------------------|:------------|
| source ALT alleles                 | 4,096,123   |
| model-addressable literal alleles  | 4,095,611   |
| annotation rows                    | 47,835,851  |
| resident regulatory/motif features | 1,383,580   |
| transcript flank distance          | 5,000 bases |
| source revision                    | ca35fd7b    |

Throughput is input alleles per second:

| projection  | one pinned core | four pinned cores |
|:------------|----------------:|------------------:|
| compact     |         961,411 |         3,172,433 |
| rich        |         430,211 |         1,527,643 |
| HGVS        |         236,385 |           838,749 |
| rich + HGVS |         179,861 |           631,357 |

Input staging and model load are outside this resident-engine
measurement. Each pass validates the complete output relation, and a
full-row fingerprint checks one-core/four-core equality for every
projection. Commands, CPU affinity and the complete evidence are in
[benchmarks/duckvep_throughput.md](benchmarks/duckvep_throughput.md).

The end-to-end comparison includes VCF decoding, DuckVEP’s coordinate
sort, annotation, and uncompressed tabular output for both DuckVEP and
the current native-compiled
[FastVEP](https://github.com/Huang-lab/fastVEP) checkout.

| pinned physical cores | DuckVEP seconds | FastVEP seconds | FastVEP / DuckVEP |
|----------------------:|----------------:|----------------:|------------------:|
|                     1 |           64.54 |          164.38 |             2.55x |
|                     4 |           32.06 |           68.09 |             2.12x |

The tools emit different native tabular schemas. The comparison
therefore keeps the row and byte denominators, credits FastVEP for
native presentation fields left as placeholders in this DuckVEP speed
projection, and separately checks each result against VEP.

The historical supplementary result compares the same ClinVar 2026-07-06
release, 4,095,611 HG002 literal ALT queries, eight requested
clinical/frequency fields, and 44,561 allele hits. DuckDB reads typed
payload columns from collision-safe Parquet joins; FastVEP returns the
complete fastSA JSON payload.

| threads | DuckDB typed Parquet join | FastVEP fastSA lookup | fastSA / DuckDB |
|--------:|--------------------------:|----------------------:|----------------:|
|       1 |                      0.90 |                  3.42 |           3.81x |
|       4 |                      0.26 |                  4.09 |          15.69x |

The illustrated [*DuckVEP: the fastest Ensembl VEP-compatible
consequence predictor in the
West?*](benchmarks/benchmark_duckvep_fastvep.md) report retains the
historical source-comparison evidence and its denominators. It is not an
executable cross-machine benchmark: a renewed campaign must first
declare cache staging for every model and provider input. Supplementary
sources remain ordinary SQL relations rather than provider-specific C or
Rust code.

#### Current scope and explicit gaps

The released claim is independent-event consequence and HGVS parity for
the tested VEP 116 state space. Genotypes and phase sets are preserved
by the HTS relations, but combined phased transcript edit sets,
haplotype consequences, and haplotype HGVS remain separate work. Exact
span SVs and paired BNDs are implemented; confidence-interval
interpretation, STR-specific classification, and arbitrary
producer-specific symbolic records are not silently guessed. GRCh37 uses
its native Ensembl release-116 transcript model and carries no invented
native MANE designation; any MANE projection back to GRCh37 must be a
separate provenance-labelled mapping.

The resident kernel accepts any transcript catalog that satisfies its
validated model contract. The supplied builder currently compiles the
Ensembl core transcript set, not VEP’s `--refseq` or `--merged` sets. It
retains GENCODE Basic and Primary flags for SQL filtering but does not
prefilter the resident model to either set. The builder withholds
sequence for source rows marked with transcript-level sequence
corrections. A model with inserted or deleted transcript bases relative
to its genomic exons also needs a richer coordinate map than the current
resident interface. Mitochondrial genetic codes and ordinary MT
coordinates are supported; events, transcripts, or features that cross a
circular sequence origin are not.

Supplementary clinical/population scores are deliberately SQL relations
rather than hard-wired C fields. The measured joins establish the
storage and execution primitives, not a universal licensed annotation
pack. Broader species/codon-table campaigns, mapped GRCh37 MANE,
combined haplotypes, and the remaining imprecise-SV/STR states stay
visible in the issue tracker and are not included in the 1.5.0 parity
claim.

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
in a persistent DuckDB session.

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
Rebuilding DuckVEP’s executable VEP conformance evidence additionally
requires R with `{targets}`, Perl/BioPerl, micromamba for the pinned
oracle environment, and [`blit`](https://github.com/WangLabCSU/blit) for
non-interactive process execution. Those conformance tools are not
runtime or extension-build dependencies.

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

MIT
