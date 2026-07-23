
# Rduckhts 1.5.0-0.1.0
- make the in-memory downstream htslib consumer test compatible with both the
  original and argument-capable `Rtinycc::tcc_call_symbol()` APIs, and bundle
  bounds-checked BAQ tag renaming so GCC 11 package builds do not emit the
  former zero-size-region warning
- clarify the bundled DuckVEP catalog and README: the resident kernel is
  transcript-source-neutral, the package builder currently compiles the Ensembl core set,
  and transcript-to-genome sequence corrections, circular-origin geometry, implicit path
  projection, and phased composition remain outside the package surface
- make the Fedora CRAN-reproduction flow install the suggested `Rtinycc`
  downstream-linking test dependency before strict package checks
- clarify in the bundled DuckVEP function catalog that the supplied reference defines the
  exact modeled sequence paths and that annotation does not implicitly project X/Y PAR,
  patch, or alternate-haplotype coordinates. Package documentation now distinguishes
  official Ensembl Variation release-VCF product audits from executable-VEP compatibility
- bundle the VEP-116 protein-HGVS endpoint fix: a coding insertion that copies
  the final translated residue no longer becomes a protein duplication when
  VEP's post-variant peptide lookup would suppress the 3-prime shift
- bundle `duckvep_allele_geometry(...)`, which returns the raw VCF span, VEP
  feature span, minimized edit span, and insertion point used by the consequence
  engine. The bundled cumulative consequence/HGVS path now reuses one projected
  transcript edit and one coding fact object per event/transcript pair, with
  release-specific VEP-116 HGVS behavior selected through a named compatibility
  policy rather than formatter-local constants. Add DBI coverage for the geometry
  contract and regenerate the package function catalog
- make the bundled `duckhts_bcftools_norm(...)` derived-query macro accept the
  `VARCHAR[]` ALT column returned by `read_bcf()`. The scalar and list forms of
  `duckhts_alt_to_list(...)` and `bcftools_norm_row(...)` are now registered as
  DuckDB overload sets, matching the package's documented normalization workflow
- bundle peak-capacity preflight for exact-alias multi-edit CDS application. A
  length-neutral edit set whose first applied insertion temporarily grows the CDS now
  returns `BUFFER_TOO_SMALL` before changing the caller buffer when that intermediate
  sequence exceeds its capacity
- add `rduckhts_htslib_config()` as the versioned installed-package contract for
  downstream C/C++ packages linking to the bundled htslib. It resolves exact headers,
  shared/static libraries, static dependencies, loader flags, plugins and compiled
  features from the installed package, and rejects receipt/header/runtime version
  disagreement through the extension-owned htslib diagnostics. Add thin
  `rduckhts_htslib_info()` and `rduckhts_htslib_version()` helpers plus an in-memory
  downstream compile/link test, built through Rtinycc with active-SDK system headers
  on macOS, that opens bundled BAM, CRAM, BCF, and VCF files. When `link` is omitted,
  select the shared or static contract chosen during package configuration. The
  compatibility `htslib_rpath()` helper follows that configured contract, so static-only
  installations return an empty loader path instead of validating a missing shared library
- bundle the projection-aware `read_bigwig(...)` table function and
  `rduckhts_bigwig()` materialization/view helper. Region vectors use htslib's
  one-based inclusive syntax while returned intervals retain stored zero-based,
  half-open coordinates. Add the pinned upstream libBigWig fixture, DBI tests, package
  examples, Devon Ryan's copyright-holder attribution, and the vendored MIT license.
  Oversized chromosome tables, malformed chromosome-tree IDs, truncated zoom headers,
  unreadable R-tree indexes, oversized block allocations, and truncated data blocks return
  reader errors before indexed storage or interval decoding
- reorganize the rendered package DuckVEP example into explicit design/validation,
  supported-scope, current-gap, and executable-example sections. The documentation now
  distinguishes independent-event consequence/HGVS support from combined haplotypes,
  imprecise-SV/STR policy, and supplementary provider relations
- bundle VEP-116 compatibility fixes for complete first-codon in-frame deletions on
  transcripts without a 5-prime UTR and for pure insertions in incomplete terminal codons.
  The latter can now return the exact VEP combination
  `incomplete_terminal_codon_variant&inframe_insertion&stop_gained` together with
  `c.280_281insAGT` and the apparently contradictory `p.Ter94=` protein rendering.
  A bundled DBI regression exercises both rare states through the public
  `duckvep_annotate(..., rich = TRUE, hgvs = TRUE)` relation; its terminal-codon case
  uses an exact remapped GRCh38 transcript/reference fixture with enough downstream
  sequence to pin the one-base VEP 3-prime shift
- make bundled mitochondrial HGVSp reproduce VEP 116's late termination search, which
  uses the standard codon table even when ordinary transcript translation uses
  mitochondrial table 2, and render a combined frame-changing terminal deletion plus
  cached stop loss through VEP's deletion-extension precedence
- bundle rare HGVS compatibility fixes. Shifted in-frame start-loss
  insertions perform VEP's sequence-dependent peptide-level 3-prime rotation;
  position-1 shifted frameshifts retain `Ter?`, while later frameshifts can
  repeat VEP's late original-coordinate, original-allele stop search; and
  terminal complete-feature clamping can render a duplicated
  transcript base even when minimization moved the differing base outside the
  transcript. Bundled SQL behavior is pinned by exact HGVSc/HGVSp witnesses
- bundle `duckvep_annotate(..., rich := TRUE)` so DBI workflows can request
  VEP-facing consequence, IMPACT, region, amino-acid, NMD, overlap-object, and
  explicit DuckHTS audit text without losing the compact numeric fields.
  `rich := TRUE, hgvs := TRUE` returns those fields and independent-event HGVS
  from one fused native pass. Add DBI coverage for the public rich/HGVS schema
  and regenerate the package function catalog from the public SQL authority
- make the bundled `duckvep_annotate(...)` relation reject conflicting
  supported symbolic ALT and explicit structural-type values, and treat
  explicit `hgvs := NULL` as the documented `FALSE` default instead of
  filtering small variants from both annotation lanes
- bundle the unified `duckvep_annotate(events_table, model_name, hgvs := false,
  ...)` relation and generated `duckvep_so_terms()` lookup. DBI callers can keep
  one narrow mixed small/SV/BND event table, request independent-event HGVS
  without changing schemas, decode consequence masks after filtering, and join
  selected rows back by event identity. Add a bundled public-surface tinytest and
  a deterministic model-to-annotation example in the rendered package README
- bundle `duckhts_contig_key(...)` for DBI workflows that need a conservative join key
  across `chr`-prefixed VCF/dbSNP contigs and Ensembl-style region names. The helper
  normalizes mitochondrial `M`/`MT` to `MT` and uppercase `X`/`Y` but deliberately does
  not guess numeric sex chromosomes, accessions, patches, or alternate loci
- reset bundled DuckVEP exon cursors before a worker workspace is reused after a
  non-monotone normalized-event vector, preventing a transcript skipped after a
  forward jump from carrying an ahead exon rank into the next DBI vector
- fix bundled cumulative HGVS output for strings larger than the native
  renderer's initial scratch capacity. The bundled adapter retries with the
  reported exact capacity instead of exposing truncated text, retains a failed
  UTF-8 assignment as NULL, and rejects an invalid internal text slice. A DBI
  regression exercises a transcript HGVS string longer than 1,400 bytes and an
  independently rendered protein HGVS string beyond the initial capacity.
  Cumulative annotation preserves resident regulatory/motif rows from the same
  pass with NULL HGVS fields, so DBI callers do not need a second annotation scan
- speed up the bundled DuckVEP consequence and cumulative HGVS paths by reusing one
  classification pass and compact prepared facts. Bundled annotation now preserves
  deterministic `annotation_index` values across DuckDB vector and disjoint ordered
  partition starts. The resident model remains immutable and shared.
  Tests also cover zero-, 10,000-, and 50,000-base transcript distances so the bundled
  behavior is not specialized to the default distance
- bundle `duckvep_annotate_hgvs(...)` for DBI workflows that need compact independent-event
  consequence rows together with transcript `c.`/`n.` and default-VEP protein `p.` HGVS
  suffixes, the applied 3-prime shift, and explicit supported/unresolved/not-applicable
  states. Bundled `duckvep_model_load(...)` can bind an existing indexed reference FASTA
  through an exact ordinal/name/length relation; it does not create an index and retains
  open read descriptors for the validated FASTA, `.fai`, and optional `.gzi`. Linux workers
  reopen those descriptors, Windows keeps a resolved source under deny-write sharing, and
  other POSIX workers use independent resolved-source handles with identity checks rather
  than sharing `/dev/fd` seek state. Annotation
  workers own separate faidx handles, reuse contained sorted reference windows, and reject
  detectable in-place source mutation around a fetch. Explicit NULL optional model queries
  and reference parameters behave like omission. Bundled tests cover reference-backed substitution and
  insertion rendering through the public DBI surface. Transcript rows admitted only by an
  upstream/downstream distance report HGVS `not_applicable`, matching VEP's absent HGVSc,
  rather than an unresolved projection. Literal exonic SNP HGVSc also retains VEP 116's
  phase-aware CDS-start fast path without shifting intronic SNP, indel, or multi-base
  feature coordinates. Protein HGVS reports `not_applicable` when the VEP feature has a
  leading or trailing genomic-to-peptide mapper Gap. Reference failures use VEP's cached
  complete-feature coding predicate, so a 5-prime-UTR insertion that could shift into CDS
  remains protein-unresolved rather than falsely not applicable. A bundled model without FASTA now
  reports `missing_reference` when retained uploaded REF padding or an anchor cannot be
  checked by the prepared CDS, instead of validating only the minimized differing REF.
  The bundled reference path keeps VEP's exact +/-1000 shift slice separate from complete
  uploaded-REF validation and adjacent duplication-source lookup, so retained padding does
  not change the shift and copied sources longer than 1000 bases still render as `dup`.
  Endpoint-overlapping transcript edits retain VEP's clipped transcript-slice coordinates,
  and protein replay reproduces VEP's one/two-base alternate-CDS trimming assignment bug
  before appending the 3-prime UTR. Bundled HGVS execution now reuses the kernel-prepared
  model, derives each reference first stop once, defers CDS projection until protein HGVS
  needs it, scans a virtual single-edit CDS instead of rebuilding the complete alternate
  sequence, and renders through a reusable worker buffer. The bundled consequence sidecar
  shortcut remains fail-closed for length-changing splice overlaps: without a positive
  frameshift fact, those edits run the complete peptide delta so their VEP-compatible
  frameshift HGVSp is retained
- bundled DuckVEP now gives ordinary long literal deletions the same complete-
  transcript `transcript_ablation` semantics as symbolic deletions. Equal-length
  containing alleles also preserve VEP 116's otherwise-empty endpoint UTR terms
  without leaking an internal unknown-coding fact into
  `coding_sequence_variant`; bundled DBI tests pin both public results
- speed up bundled paired-BND annotation by preserving the already sorted
  transcript stream and linearly merging only the much smaller resident
  regulation-feature stream; transcript-only models no longer sort the full
  expanded result before returning one list per input event
- bundled paired-BND annotation now keeps VEP's fixed 5000-base overlap-allele admission
  independent of a caller-selected zero transcript window: an admitted local allele whose
  directional predicate is disabled contributes default `intergenic_variant` beside
  mate-derived `feature_truncation`. Add bundled DBI regressions at zero and wider caller
  distances
- bundled DuckVEP exact structural annotation now accepts `STR`/`TANDEM_REPEAT`/`CNV:TR` spans and preserves structural tandem-repeat identity while matching VEP 116's tandem-duplication consequence predicates. Bundled paired-BND annotation now returns RegulatoryFeature and MotifFeature rows found at either endpoint, once per feature, in addition to transcript consequences; local exact hits keep the base object term, mate-only exact hits use VEP's generic HIGH-impact `feature_truncation`, and a shifted local point on the same contig but outside and within 5000 bases adds `intergenic_variant` to that mate-discovered object. The local base term wins when both points hit one object exactly, and mixed transcript/regulation batches preserve one result list per input event. A caller-selected transcript distance above 5000 does not widen VEP's separate fixed structural-breakend allele-admission cap. Add DBI tinytests for these surfaces and update generated function documentation
  Structural adapters consume nominal `POS`/`END` geometry, matching VEP 116's registered
  consequence predicates; callers retain `CIPOS`/`CIEND`, inserted sequence, and raw repeat
  metadata beside the result for provenance, later HGVS, and round-trip rendering
- fix bundled projected `read_fastq(...)` scans so FASTQ headers longer than the BAM query-name limit remain usable when `NAME` / `PAIR_ID` are not requested and no paired-file comparison needs them; name-producing scans retain the htslib-compatible limit. An explicit bundled `duckhts_fastq_qc(..., max_cycles)` value below 128 now also caps the initial per-group cycle allocation
- bundle `duckhts_fastq_qc(sequence, quality [, max_cycles])` for DBI workflows that need exact global and per-cycle FASTQ quality statistics without expanding one row per base; the bundled scalar, AVX2, ARM NEON, and wasm SIMD128 implementations share one dispatch contract, and tinytests compare forced scalar with automatic selection over grouped and malformed inputs
- accelerate bundled `read_fastq(...)` scans by parsing FASTQ directly over htslib transport into projected DBI columns instead of round-tripping through temporary BAM records; multiline records, htslib-style query names, paired/interleaved checks, string and packed nt16/Phred output, quality conversion, and count behavior are retained, while truncated quality blocks now surface as R/DBI errors. Add bundled multiline, packed-output, and truncation tinytests
- fix `duckhts_build(...)` to obtain the DuckVEP kernel source list from the package's shared source-list helper; the fallback no longer references a variable local to `duckhts_bootstrap(...)`
- bundle htslib 1.24 and use its native deduplicating BCF/tabix multi-region iterators in the DBI-visible indexed readers. Repeated or overlapping regions passed to bundled `read_bcf(...)`, `read_bcf_v2(...)`, `read_gff(...)`, `read_gtf(...)`, `read_tabix(...)`, and `read_bcf_appender(...)` calls now return each matching record once. The appender partitions `region_threads` work by primary contig while retaining one native iterator for that contig's complete interval set, so thread count no longer changes its target schema or row multiset; remove the thread-only `duckhts_region_idx` column and retain opt-in `FILE_OFFSET` for restoring file order. Bundled `read_bam(...)` already uses the equivalent deduplicating alignment iterator. Keep BCF work claiming iterative and bounded to 16 workers, with bundled regressions over 256 leading header-only contigs and empty requested contigs. Build htslib as GNU C17 so GCC 16 does not apply its GNU C23 default to this pre-C23 dependency. Add bundled tinytests and a Fedora GCC 16 R-devel package check matching CRAN's Fedora compiler environment
- bundle release-matched Ensembl RegulatoryFeature and MotifFeature preparation and integrated consequence output for DBI workflows. `duckvep_model_load(...)` accepts the compact feature projection, while `duckvep_annotate(...)`, its compact form, and exact structural adapters return typed transcript or regulation/motif rows from one call; cold funcgen metadata remains joinable by feature ordinal. The existing eight-argument model-receipt call remains valid and accepts regulation_features_table as an optional named relation. Bundled model receipts validate, count, and hash feature geometry, exclude the same EMAR rows as VEP 116, and the offline fixture includes exact core/funcgen provenance plus three real MotifFeature rows on `KI270395.1`
- bundle dedicated `duckvep_annotate_breakend(...)` and
  `duckvep_annotate_breakend_compact(...)` functions. DBI callers provide both raw VCF
  loci in one row and retain raw ALT, bracket orientation, event identity, and provenance
  as ordinary columns; the bundled engine queries transcripts around both endpoints and
  returns the VEP-116 consequence-set union once per transcript. Mate-only consequences
  have NULL rich region / zero compact region mask
- bundle typed exact structural-span consequence functions `duckvep_annotate_sv(...)` and
  `duckvep_annotate_sv_compact(...)`; DBI callers provide one-based inclusive spans for
  DEL/DUP/TDUP/INV/CNV/UNKNOWN, or `start = end = P` for an INS after reference base P,
  plus explicit copy direction. Contradictory metadata fails as an error. The bundled
  consequence mask retains stable assignments for all 41 terms in VEP 116's registry,
  including regulatory-region and transcription-factor-binding-site terms
- bundled DuckVEP now matches the pinned VEP Plugins release/116 NMD implementation's
  separate feature geometry: consequence and sequence changes use the minimized edit,
  while NMD CDS and exon-end rules use the complete VEP feature. This fixes padded
  equal-length alleles crossing an early-CDS or penultimate-exon threshold and insertions
  at coding exon edges, whose parent VEP object keeps a reversed CDS range—including the
  defined `1,0` range immediately before CDS base 1; the cached first-101-CDS-base result
  is reused only when both spans are identical
- bundled DuckVEP now reuses resolved intron coordinates while evaluating non-SNV mismatch islands and skips predicate evaluation when no intron, splice, or short-intron fact can change; the sorted-SNV classifier remains inlined
- bundled DuckVEP now prepares validated CDS-to-cDNA projection coordinates once per resident transcript, avoiding repeated CDS-endpoint searches while annotating coding variants; topology-only transcripts retain the existing exhaustive projection behavior, and the common sorted-SNV lane keeps a stable instruction-cache-aligned entry
- bundled DuckVEP now rebuilds multi-edit CDS haplotypes in one linear
  reverse-coordinate pass over distinct reference and worker-scratch buffers, replacing
  one CDS-tail move per edit while preserving edit, strand, reference-validation, and
  exact-alias behavior
- bundled `duckvep_annotate(...)` and `duckvep_annotate_compact(...)` now accept VEP-style separate upstream and downstream distances in addition to the existing shared distance; omitted values retain the 5,000-base default and zero disables that direction
- bundled DuckVEP now reproduces VEP 116's reverse-strand insertion-length terminal-stop fallback, reconstructing the original CDS endpoint after the edit before deciding whether the terminal stop is retained
- bundled DuckVEP now reproduces VEP 116's independent unknown-coding and missense predicates for equal-length edits on incomplete-start transcripts
- bundled DuckVEP now matches VEP 116 for an insertion after the final base of a mature-miRNA segment: VEP's minimized reversed insertion interval does not overlap that segment, so annotation returns `non_coding_transcript_exon_variant` instead of a false `mature_miRNA_variant`; the regression covers both transcript strands
- bundled DuckVEP model preparation now reads Ensembl sequence-region codon-table attributes and supports the full VEP 116/BioPerl codon-table set. Valid single-residue initial-methionine, selenocysteine, curated amino-acid-substitution, and stop-readthrough edits can be packed by `duckvep_model_load(...)` and are applied to the reference peptide during annotation; unsupported edit shapes remain fail-closed. Numeric-only attribute columns inferred as integers by dump staging are accepted at the importer boundary. Model receipts hash and count mature-miRNA segments and peptide edits, and the bundled GRCh38 and GRCh37 fixtures exercise real mitochondrial codon-table and peptide-edit behavior
- `duckvep_ensembl_transcripts(...)` now applies VEP's core-transcript selection and
  projects mature-miRNA cDNA attributes into genomic exon segments;
  `duckvep_model_load(...)` accepts those segments through optional `mature_mirna_query`,
  and annotation returns `mature_miRNA_variant` with VEP's precedence over generic
  non-coding exon terms. The bundled engine also fixes insertion placement at
  transcript/exon edges, VEP's empty-UTR endpoint behavior, and coding edits on
  `CDS_END_NF` transcripts with a one- or two-base partial terminal codon; DBI tests cover
  the complete Ensembl-import-to-resident-model path
- bundle `duckvep_annotate_compact(...)`, returning the same DuckVEP transcript consequences as numeric SO/region masks and stable status, reason, amino-acid, and NMD codes so DBI workflows can filter and join before rendering strings. The bundled rich annotation path also caches DuckDB vector validity and string lengths per chunk; SQL and tinytests check numeric-code semantics and parity with `duckvep_annotate(...)`
- bundled DuckVEP now uses one edit/CDS/peptide context for length-changing small
  variants instead of retrying failed contexts through a shape-specific classifier.
  Complete transcript-oriented sequence before and after the CDS supports VEP 116's
  independent start/stop predicates, endpoint-mapper-Gap behavior, and
  transcript-associated default `intergenic_variant`; `duckvep_model_load(...)` accepts
  the complete 13-column projection while retaining its 11/12-column compatibility, and
  older models return `missing_transcript_flank` rather than guessing
- bundled `duckvep_annotate(...)` now matches VEP 116 for equal-length features crossing from 5-prime UTR into the start codon: changing the start returns `start_lost` without a co-occurring stop term, while preserving `ATG` returns `start_retained_variant` even when later coding bases in the same uploaded feature change; DBI regressions pin each state
- bundled `duckvep_annotate(...)` now matches VEP 116 when unchanged bases retained in an equal-length uploaded feature widen the peptide-predicate window beyond the trimmed sequence edit, including representation-dependent `start_lost&start_retained_variant`, `stop_retained_variant`, `missense_variant`, and `stop_gained` results pinned by paired DBI regressions; retained REF mismatches and ambiguous widened windows now return explicit unresolved reasons instead of retrying the smaller edit
- bundled DuckVEP annotation now matches VEP 116 when an equal-length uploaded feature crosses from CDS into the transcript-oriented 3-prime UTR: it preserves VEP's unavailable peptide-mapping state and returns `coding_sequence_variant` rather than a stop or missense consequence derived from a smaller trimmed edit
- make bundled `duckvep_annotate(...)` reproduce VEP 116 when ALT-only mismatch bases extend from an exonic REF-shaped feature into an intron, including VEP's three-base interval-tree cache boundary and co-emission with coding and splice terms
- make bundled `duckvep_annotate(...)` classify length-changing coding edits from VEP 116's codon-local predicate inputs rather than whole-protein or net-length shortcuts. This fixes stop-gained insertions whose preserved peptide flank lies after the new stop, and delins that change length by a complete codon but are `protein_altering_variant` rather than in-frame insertion/deletion
- make the bundled `duckvep_annotate(...)` match VEP 116's terminal-stop insertion
  predicate order, including non-modulo-three insertions that VEP reports as
  `inframe_insertion` together with `coding_sequence_variant`, `stop_lost`, or
  `stop_retained_variant`; unsupported contexts remain explicit
- bundled `duckvep_annotate(...)` now keeps uploaded VCF geometry, VEP parser feature geometry, and the semantic sequence edit distinct, fixing candidate, splice, start-codon, and terminal-stop consequences for padded and multi-base alleles
- bundled `duckvep_annotate(...)` now predicts whether eligible stop-gained, frameshift, splice-donor, and splice-acceptor consequences trigger or escape nonsense-mediated decay under the pinned VEP Plugins release/116 policy, with separate boolean escape reasons and explicit unresolved state when coding placement is unavailable. This remains distinct from `NMD_transcript_variant`, which identifies an already curated nonsense-mediated-decay transcript biotype
- bundled `read_bcf()` and `read_bcf_v2()` now recognize the `Format=...` CSQ
  schema spelling in Ensembl variation release VCFs and expose typed `VEP_*`
  columns through DBI queries
- retain versioned RefSeq accessions from MANE Select and MANE Plus Clinical in bundled prepared transcript relations while keeping the resident model compact; reject empty or conflicting MANE mappings. Fix bundled `duckvep_model_receipt(...)` provenance column names and add roughly 116 KiB of offline Ensembl-116 GRCh38 and GRCh37/GENCODE-19 fixtures to package tests, covering real MANE, ordinary coding, and mitochondrial missing-sequence behavior without network access during CRAN builds or checks
- bundle `duckvep_ensembl_regions(...)`, `duckvep_ensembl_transcripts(...)`, and `duckvep_model_receipt(...)`, allowing DBI workflows to prepare a validated, provenance-hashed resident DuckVEP model directly from Ensembl core relations and matching tiled FASTA sequence. Unsupported Ensembl RNA/peptide edits keep their model flags but return explicit missing-sequence state instead of ordinary coding predictions, including `_rna_edit` records carried by either transcript or translation attributes; bundled SQL and tinytests cover both strands, nested exon projection, receipt generation, and resident loading
- accept Ensembl exon phase `-1` in bundled sequence-backed models when translation begins after 5-prime UTR within that exon, matching VEP's zero-prefix interpretation while retaining exact prepared-CDS validation
- bundled DuckVEP annotation now follows VEP 116 for frame-changing edits that begin in the terminal stop, returning `stop_lost` or `stop_retained_variant` instead of a false frameshift. The resident transcript query accepts an optional twelfth `post_cds_bases` BLOB with up to three transcript-oriented bases; when a terminal deletion needs absent tail sequence, the result is explicitly unresolved as `missing_transcript_tail`
- bundled DuckVEP annotation now matches VEP 116 for insertions at exon, CDS, and transcript boundaries, including coding+splice consequences when the VCF padding base is intronic and UTR/downstream placement at right-hand boundaries
- expose specific DuckVEP unresolved reasons through the bundled extension and return NULL, rather than empty strings, when a protein-positioned frameshift or in-frame indel has no scalar one-letter amino-acid value
- reject malformed bundled DuckVEP resident models at load time when transcript/exon/cDNA/CDS/phase/prepared-sequence coordinates disagree; the SQL loader and standalone C engine now share the same final validator
- make bundled DuckVEP intergenic output fail closed for partial resident models. Ordinary `duckvep_model_load(...)` calls report `no_feature_in_loaded_model` as unresolved; the named `transcript_coverage_complete := TRUE` parameter requires contig lengths before returning supported `intergenic_variant`, and rejects coordinates beyond those lengths
- fix bundled `duckvep_annotate(...)` handling of padded small variants: whole-allele padding no longer turns a one-base substitution into an MNV, and VCF position-1 insertions/deletions with a following padding base are accepted
- speed up the bundled DuckVEP engine for coordinate-sorted SNVs by retaining compact per-transcript exon cursors across adjacent DuckDB chunks and by avoiding full Sequence Ontology scans for the common single-consequence rows; this does not change the existing MNV/indel/SV span classifiers or grouped-haplotype edit core
- bundled DuckVEP annotation now emits VEP-116 `NMD_transcript_variant` for variants inside transcripts imported with the `nonsense_mediated_decay` biotype
- bundle the DuckVEP SQL surface: `duckvep_model_load(...)` loads and validates one of several named resident transcript models from committed DuckDB relations, `duckvep_annotate(...)` returns explicit per-transcript consequence rows for biallelic small variants, and `duckvep_model_drop(...)` releases a model. Add an end-to-end DBI tinytest using the bundled extension and a sequence-backed coding transcript

# Rduckhts 1.4.0-0.1.0 (2026-07-10)
- accelerate the bundled `bam_nt16_counts` and `nt16_gc_counts` logical kernels on AArch64 with NEON, so Apple Silicon and other arm64 package builds no longer fall back to scalar for BAM-bin GC counting or nt16 `seq_gc_content(...)`; add conditional tinytests that require every reported complete backend to own all three sequence-kernel slots and match scalar results. Emscripten package builds now keep `bam_bin_counts(...)` on its synchronous no-index path by disabling unavailable htslib worker threads
- bundled `cigar_*` helpers (`cigar_has_soft_clip`, `cigar_has_hard_clip`, `cigar_left_soft_clip`, `cigar_right_soft_clip`, `cigar_query_length`, `cigar_aligned_query_length`, `cigar_reference_length`, `cigar_has_op`) are overloaded to accept a `UINTEGER[]` binary CIGAR from `read_bam(cigar_representation = 'binary')`, and bundled `seq_hash_2bit(...)` is overloaded to accept a `UTINYINT[]` of htslib nt16 codes from `read_bam(sequence_encoding = 'nt16')`. Both are bit-identical to the text path, so DBI pipelines analyze the projection-pushed binary columns without decoding to text; add bundled tinytest coverage
- bundled `rduckhts_simd_kernel_info()` now reports an additional `nt16_gc_counts` logical kernel, and the bundled nt16 `seq_gc_content(...)` overload is routed through the SIMD dispatch framework; GC results are unchanged (a scalar reference backend is the correctness oracle), and forced `scalar` vs `auto` agree (bundled tinytest coverage added)
- bundled `seq_revcomp(...)` and `seq_canonical(...)` are overloaded to accept BAM nt16 sequences: in addition to text they take the `UTINYINT[]` of htslib nt16 codes from `read_bam(sequence_encoding = 'nt16')` and return `UTINYINT[]`, so DBI pipelines can reverse-complement and canonicalize directly on the nt16 column. Non-ACGTN codes yield NULL as in the text path; results are bit-identical to the text path after decoding; add bundled tinytest coverage
- bundled `seq_gc_content(...)` is overloaded to accept BAM nt16 sequences: in addition to a text sequence it now takes the `UTINYINT[]` of htslib nt16 codes returned by `read_bam(sequence_encoding = 'nt16')`, so DBI pipelines can compute GC directly on the projection-pushed BAM column without decoding to text. The nt16 result is bit-identical to the text path; add bundled tinytest coverage
- bundled `read_bcf(...)`, `read_bcf_v2(...)`, and `read_bcf_appender(...)` gain `decode_error_policy := 'null'|'warn'|'error'` for corrupt BCF FORMAT/INFO header-vs-payload type clashes; the default `null` policy materializes NULLs, `warn` emits a DuckHTS warning and materializes NULLs, and `error` raises a DuckDB/R error. FORMAT and INFO decode now preflight under every policy (including the default `null`), so corrupt inputs neither trigger htslib `exit(1)` termination nor leak raw bytes; an INFO field whose header claims `Type=String` over a numeric payload previously returned truncated garbage under the default policy and now materializes NULL. Add bundled corrupt-BCF fixtures (including a reverse String-header/numeric-payload clash) plus tinytest coverage while keeping valid mixed-ploidy `Number=G` FORMAT records accepted
- bundled `rduckhts_simd_kernel_info()` now reports an additional `bam_nt16_counts` logical kernel, and bundled `bam_bin_counts(...)` GC base counting is routed through the SIMD dispatch framework; GC results are unchanged (a scalar reference backend is the correctness oracle)
- surface malformed-record BCF/VCF scan failures from bundled `read_bcf(...)`, projected `read_bcf_v2(...)` scans, and non-parallel `read_bcf_appender(...)` as DuckDB/R errors instead of treating htslib parse/read failures as EOF; run bundled `read_bcf_appender(...)` writes in an internal transaction so malformed input rolls back target-table side effects; add a malformed-POS fixture plus tinytest coverage
- bundle the restored experimental `read_bcf_v2(...)` table function for DBI users, preserving `read_bcf(...)` schema compatibility while adding sample pushdown, INFO/FORMAT/VEP field filters, projection-aware VCF unpacking, persistent decode caches, and count-only shortcuts for benchmarking
- bundle the initial experimental `read_bcf_appender(...)` benchmark helper for DBI users; at introduction, its single-stream path was transactional while the `region_threads > 1` path was explicitly best-effort and lacked rollback (the current development entry above records its later transactional, thread-invariant contract)
- bundle the experimental `hts_region_union_query(...)` scalar macro for DBI users, which generates a `UNION ALL BY NAME` query string reading one HTS file through separate per-region scans; the generated query does not deduplicate boundary-spanning records
- fix bundled `read_bcf_appender(...)` so a region query matching no records still creates and commits the empty target table instead of rolling back its `DROP`/`CREATE TABLE`, which previously could leave a stale prior table under `overwrite := true` or never create the target
- add bundled tinytest coverage for malformed `additional_csq_column_types` so invalid CSQ/ANN/BCSQ type override rules fail as DuckDB/R errors instead of proceeding with ambiguous parsed-annotation typing
- update the bundled DuckDB C API headers to DuckDB v1.5.3 while keeping stable extension ABI metadata at v1.2.0; bundled SQL sessions loaded through `rduckhts_load()` now expose DuckDB runtime type-support probes for `VARIANT` and `GEOMETRY`
- simplify the bundled `rduckhts_bcftools_norm()` / `duckhts_bcftools_norm(...)` site-preserving table-macro query shape by removing the extra correlated scalar `LATERAL` subquery around `bcftools_norm_row(...)`, eliminating the site-preserving `LEFT_DELIM_JOIN` plan overhead while preserving split-mode ALT row semantics and caller columns whose names collide with DuckHTS helper-column names used internally by earlier macro forms; add tinytest coverage for DuckDB's suffixed behavior when callers already have normalized-output column names
- expose bundled reader `scan_mode = "auto"|"sequential"` controls through the R wrappers and multi-file helpers for `read_bcf`, `read_bam`, `read_fasta`, `read_fastq`, `read_bed`, `read_gff`, `read_gtf`, and `read_tabix`, so callers can force full-file streaming/counting instead of index-backed count or parallel scan paths where applicable; sequential mode is rejected for region queries
- optimize bundled `bcftools_norm_row(...)` / `rduckhts_bcftools_norm()` for already-normalized plain ACGTN allele rows by skipping kstring left-realignment setup when trim predicates prove the row is unchanged; avoid per-row FASTA path duplication after the vector-local cache is established, reuse larger bounded per-thread reference windows, and document/defensively serialize htslib FASTA fetches while keeping normalization reference caches thread-local to avoid the faidx cache race class fixed in https://github.com/RGenomicsETL/duckhts/issues/17 / https://github.com/RGenomicsETL/duckhts/pull/18
- make bundled `rduckhts_bcftools_norm()` / `duckhts_bcftools_norm(...)` gVCF-aware for vt/vcfnorm-style row normalization: `<NON_REF>` and `<*>` reference-block alleles now pass through with `GVCFReferenceBlock`, and mixed real-plus-gVCF-symbolic alleles normalize the real alleles while preserving symbolic alleles and caller-supplied reference-block `END` in site-preserving output; mixed `*` plus real alleles now follows the same ignored-symbolic path, while `*`-only rows remain `SpanningDeletion`; bundled phased GT/PL/GP/DS/PS FORMAT fixtures, including haploid/triploid/tetraploid `Number=G` cardinality cases, and tinytests pin phase-separator preservation through `read_bcf(...)`
- add thin DBI wrappers `rduckhts_bcf_convert_parquet()`, `rduckhts_bam_convert_parquet()`, `rduckhts_gff_convert_parquet()`, and `rduckhts_tabix_convert_parquet()` around the bundled extension SQL builders `duckhts_*_convert_parquet_sql(...)`; these convert DuckHTS scans to Parquet with DuckHTS write-format metadata, preserved raw headers, optional corrected header text, SQL-filter provenance, selected-column/partition metadata, arbitrary user metadata via R named lists/extension `metadata := map(...)`, optional caller-managed JSON-file metadata when DuckDB's `json` extension is available, and partitioned-output support for DuckLake-style registration of premade Parquet files
- include the final VCF `#CHROM`/sample header line in bundled `read_hts_header(..., mode := 'raw')`, so Parquet metadata written from VCF/BCF inputs has the complete header needed for future VCF/BCF regeneration

# Rduckhts 1.3.0-0.1.0 (2026-05-29)
- expose the rebuilt capability-mask SIMD dispatch diagnostics through `rduckhts_simd_kernel_info()`, keeping R wrappers thin while reporting one row per logical kernel and preserving backend-agnostic SQL/R conformance tests for `seq_gc_content(...)`
- harden bundled SIMD backend helpers: retain extension-owned backend-name validation while restoring R-side scalar/non-missing argument shape checks, clarify `selectable` versus `available` diagnostics in generated docs, and preserve ASCII SQL quotes for the rendered `duckhts_simd_set_backend('auto'|'scalar'|backend)` catalog call
- remove htslib autoconf `HAVE_*` macro guards from all bundled SIMD backend translation units; compile-time gate is now `defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))` for x86 backends, available without autoconf; runtime dispatch and scalar fallback behavior are unchanged; add scalar-vs-auto backend R correctness tests for `rduckhts_simd_set_backend()` / `seq_gc_content(...)` covering GC=0/0.5/1.0, embedded-N calling, and soft-masked lowercase bases
- drop internal `.validate_simd_backend()` R helper from SIMD wrapper functions; backend-name normalization and validation now belong to the extension, while the R wrappers only enforce that `backend` is a single non-missing character string

# Rduckhts 1.2.1-0.1.0 (2026-05-07)
- expose bundled SIMD diagnostics and explicit backend selection through SQL functions and R helpers `rduckhts_simd_backend()`, `rduckhts_simd_requested_backend()`, `rduckhts_simd_backend_available()`, and `rduckhts_simd_set_backend()`, route bundled `seq_gc_content(...)` through the new eager scalar/optional-AVX2 runtime dispatch scaffold while preserving scalar fallback behavior on ARM, wasm, and scalar-only builds, add runtime-gated AVX-512, ARM NEON, and wasm SIMD128 backend translation units where compiler-supported, keep the manual `duckhts_build()` rebuild path wired to the SIMD sources, and add README examples for the scalar/auto SIMD flow
- fix bundled `rduckhts_liftover()` / `bcftools_liftover(...)` FASTA contig alias handling during source/destination reference validation and sequence fetches, and align bundled spanning-deletion `*` allele handling with upstream `bcftools +liftover`: inputs such as `23`, `24`, `26`, `X`, `Y`, `MT`, and `chr*` aliases now resolve through the same canonical path, avoiding spurious `SourceRefMismatch` rejects for X/Y/MT indels when the bundled chain names and FASTA names differ only by canonical aliasing; bundled `*`-allele rows now follow upstream swap/ref-add semantics instead of taking the symbolic short-circuit path, full-file GIAB conformance against installed `bcftools +liftover` is now exact, and bundled SQL/tinytest coverage now pins the `23 -> chrX`, `SWAP=2`, and `SWAP=-1` regressions
- bundle the official VariantKey / RegionKey C API (Nicola Asuni, 2018; <https://doi.org/10.1101/473744>) and expose new SQL helpers through `rduckhts_load()` sessions for both the bcftools-style and raw upstream numeric surfaces: `variantkey(...)` now matches bcftools `%VKX` / `+add-variantkey` on 1-based VCF rows, large/ambiguous/symbolic alleles keep the official hashed nonreversible mode, `regionkey(...)` adds 0-based half-open span keys plus overlap helpers, bundled tinytests pin reversible and hashed cases, and the package README now includes concrete DBI examples for VariantKey / RegionKey usage
- fix bundled `rduckhts_bcftools_norm(..., split_multiallelic = TRUE)` row preservation for ref-only and empty-ALT inputs: rows with `ALT='.'`, NULL ALT values, empty ALT lists, or NULL ALT list elements no longer disappear from split-mode DBI results, bundled tinytest coverage now pins the expected `RefOnly` / `NullInput` statuses and `alt_index` behavior, bundled `rduckhts_bcf()` / `rduckhts_bcf_multi()` now expose `decompression_threads = 0` for explicit htslib worker-thread control on bgzipped VCF/BCF reads, and the package README now includes a concrete normalization example
- fix bundled helper-return metadata for omitted output paths: `rduckhts_fasta_index()` now returns the generated `.fai` path when `index_path = NULL` instead of an empty string, bundled regression coverage now also pins default-path returns for BGZF compression/decompression and BAM/BCF/tabix index builders, and the `rduckhts_bgzip()` / `rduckhts_bgunzip()` wrappers now correctly propagate `keep = FALSE` instead of silently falling back to the extension default `keep := TRUE`
- add `rduckhts_bcftools_norm()` and bundle `bcftools_norm_row(...)` / `duckhts_bcftools_norm(...)` for bcftools/vt-style FASTA-backed variant normalization from DBI queries: ALT inputs may be either comma-delimited `VARCHAR` or `VARCHAR[]`, the bundled result appends `pos_normed`, `end_pos_normed`, `ref_normed`, `alt_normed`, `normed`, and `norm_status`, split mode emits one row per ALT with `alt_index`, and bundled SQL/tinytest coverage now exercises sequence, multiallelic, symbolic `<DEL>`/`<DUP>`, and missing-contig rows
- fix bundled `rduckhts_liftover()` / `bcftools_liftover(...)` indel parity in two exact upstream rewrite points: repeat-run source extension now keeps extending across the cached source-reference window boundary when needed, and the bundled clip-pad `Needleman-Wunsch` path now keeps the best shift even when candidate alignment scores are negative instead of leaving padded intervals unshifted; bundled SQL/tinytest coverage now includes dedicated repeat-run and clip-pad regression fixtures, and the real-data conformance workflow reaches exact parity with installed `bcftools +liftover` on GIAB HG001 chr20 plus the full HG006 GRCh37 benchmark VCF
- fix bundled `rduckhts_liftover()` / `bcftools_liftover(...)` row rejection for invalid source-reference indel and difficult-SNP inputs: rows that fail the source-FASTA validation path now stay in the result with `mapped = FALSE` and `reject_reason = 'SourceRefMismatch'` instead of fabricating padded lifted alleles or aborting the query; bundled tests and README examples now reflect the reject-row behavior
- add `rduckhts_pileup()` and bundle native `read_pileup(...)` for region-scoped BAM pileups with per-position `chrom`, `pos`, `depth`, `bases`, and `quals`; expose bundled `read_bam(..., cigar_representation := 'binary')` through `rduckhts_bam(..., cigar_representation = "binary")` and multi-file BAM wrappers, returning packed BAM CIGAR ops as `UINTEGER[]`; and expose explicit `gzi_path` arguments in `rduckhts_fasta()`, `rduckhts_fasta_multi()`, and `rduckhts_fasta_nuc()` so packaged bgzipped FASTA workflows can use relocated `.gzi` sidecars
- speed up bundled `rduckhts_fasta_nuc()` / `fasta_nuc(...)` nucleotide counting on capable x86_64 hosts with an AVX2+popcnt fast path selected via htslib-style runtime dispatch, while preserving the scalar fallback everywhere else
- improve bundled remote HTS performance for long-running scans and `rduckhts_bam_index()`: native remote BAM/BCF/tabix/FASTA/BED reads now apply htslib block/cache tuning by access pattern, while wasm/browser builds use the same policy with smaller budgets appropriate for the XHR-backed worker runtime; the bundled vendored htslib also now exposes a pre-opened `sam_index_build4(...)` entry point so `bam_index(...)` can be tuned before remote index construction begins
- fix bundled `rduckhts_bcf()` / `read_bcf(...)` scanning stability for records where `FILTER` lists were emitted without reserving list-vector capacity, which could crash with allocator corruption (`double free`/`invalid pointer`) during full-table reads; `FILTER` entries now reserve child-list space before writes and scans are stable on files previously triggering crashes
- compile bundled DuckHTS extension sources with `-Wpedantic` during Unix and Windows package builds while leaving vendored `htslib` on its upstream warning flags
- fix the bundled non-Emscripten `wasm_http_hfile.c` translation unit so native package builds do not warn about an empty source file under pedantic C diagnostics
- harden Windows `configure.win` libcurl detection: the package now requires a successful `curl_easy_init` link using the detected `pkg-config` libcurl dependency closure before enabling htslib remote URL support, and otherwise disables libcurl/S3/GCS cleanly

# Rduckhts 1.2.0-0.1.0 (2026-05-07)
- expose richer bundled GFF/GTF parsed attribute outputs through `rduckhts_gff()` / `rduckhts_gtf()` and multi-file wrappers: `attributes_list = TRUE` returns `MAP(VARCHAR, VARCHAR[])` with grouped multi-values and GFF3 percent-decoding, while `attributes_pairs = TRUE` returns `LIST<STRUCT(key VARCHAR, value VARCHAR, idx INTEGER)>` for exact key/value/index records; `attributes_map = TRUE` remains the backward-compatible raw scalar map
- expose bundled `read_gff(..., strict := true)` through `rduckhts_gff(strict = TRUE)` and `rduckhts_gff_multi(strict = TRUE)`, enabling GFF3 structural validation from R/DBI workflows, including wrong field counts and malformed attribute segments, while keeping the default GFF reader permissive for existing ingestion pipelines
- extend bundled `rduckhts_score()` / `bcftools_score(...)` so `summary_path` can be a character vector or callers can use `summaries_list_file`; multiple TSV/SSF summaries are scored in one genotype scan, `log_path` can write per-PRS matching/audit counts for loaded, matched, allele-mismatch, and duplicate markers, `summaries_list_file` directory scans are deterministic and ignore index sidecars, generated score/count column names are validated for uniqueness, and score accumulation now follows upstream `bcftools +score` float32 summation more closely
- collapse the generated package README function catalog behind a disclosure widget so package users can jump to quick-start and workflow examples more easily
- refresh the package README release docs: clarify the bundled `htslib` 1.23.1/system-requirements wording and redact transient temp-file paths in rendered example output so regenerated README diffs stay deterministic
- add bundled `duckhts_cgranges_overlaps_list(...)`, a vectorized scalar overlap expander that returns LIST-of-STRUCT hit records so DBI queries can expand provider rows with `UNNEST(...)` without generated bulk-probe SQL; package tests cover one-row-per-hit expansion over regular tables and bundled BED data, and the existing `duckhts_cgranges_overlaps_bulk(...)` probe path now also handles DuckDB string vector lengths safely
- fix bundled `duckhts_cgranges_from_query(...)` ingestion of DuckDB string vectors by respecting string lengths instead of assuming NUL-terminated buffers; this fixes cgranges construction from providers such as `read_bed(...)` with longer chromosome names and adds package regression coverage
- add bundled vectorized scalar cgranges probe helpers `duckhts_cgranges_has_overlap(...)` and `duckhts_cgranges_count_overlaps(...)`, enabling DBI queries to stream provider rows through an already-finalized session cgranges index for filtering/count annotations without the materializing `overlaps_bulk` query-string path; add package-level coverage for overlap, contain, and NULL probe semantics
- add bundled `duckhts_cgranges_overlaps_bulk(...)` for SQL-first bulk cgranges probing from R/DBI sessions: one table-function call now streams a query of probe intervals through a finalized cgranges index, supports `mode = 'overlap'|'contain'`, accepts an optional `query_row_id_col`, and otherwise emits 1-based probe ordinals as `query_row_id`; add package-level regression coverage for the new bulk path
- document bundled `duckhts_cgranges_*` entry points in the generated function catalog and package README, add bundled DBI smoke coverage for the session-scoped cgranges registry API, and include a packaged overlap-conformance script reference for `bedtk`-style parity checks
- fix bundled `rduckhts_fasta_nuc()` / `fasta_nuc(...)` GC and AT percentages for intervals containing `N`: `pct_gc` and `pct_at` now use only informative `A/C/G/T` bases in the denominator, so ambiguous bases no longer depress reported bin/interval composition percentages; add bundled regression coverage
- add bundled C-built cgranges bulk-ingest support via `duckhts_cgranges_from_query(...)`, which runs the source query on an extension-owned DuckDB connection and builds the cgranges index in C before publishing it to the session registry; `duckhts_cgranges_from_table(...)` remains deferred for now
- bundle htslib 1.23.1 in the package for the upstream CRAM decoder and GZI validation security fixes, including the wasm/browser-exposed parsing path shipped through `Rduckhts`
- add `rduckhts_bam_bed_coverage()`, bundling native `duckhts_bam_bed_coverage(...)` for samtools coverage-like regional summaries over BED targets with DuckHTS-specific pre/post-filter columns and read-mode strand-specific post summaries; bundled SQL/tinytest coverage now checks expected outputs on the packaged mixed BAM fixture, and `fragment_mode` / `processing_threads` are exposed but currently reserved for later phases
- reduce bundled `rduckhts_bam_bed_coverage()` / `duckhts_bam_bed_coverage(...)` peak memory by allocating and freeing per-region working depth buffers during scan processing instead of retaining them for the whole BED, tile large target intervals internally when computing covered-base breadth, keep the tiled implementation single-pass, align `min_depth > 1` mean-depth behavior with `samtools coverage`, and expose `decompression_threads` so package callers can set htslib BAM/CRAM decode worker counts explicitly
- add `rduckhts_samtools_idxstats()`, bundling native `duckhts_samtools_idxstats(...)` for samtools idxstats-compatible BAM/CRAM/SAM summaries with indexed BAM fast-paths and scan fallback; package SQL/tinytest coverage now checks BAM fast-path output, CRAM fallback output, explicit `index_path`, and overwrite errors
- improve package-source hygiene for local development: ignore generated `README.html`, `.Rcheck`, staged `duckhts_extension/htslib` build outputs, wasm/webR harness byproducts, and stray root-level index files under `r/Rduckhts/`; add top-level `make clean_local` to purge the reproducible package-side artifacts
- add `processing_threads` parameter to `rduckhts_mosdepth()` and bundled `duckhts_mosdepth(...)` for parallel contig processing: workers claim contigs atomically and write output in header order; on the NA12878 WGS benchmark with 2 processing threads, fast mode is 1.38x faster, default mode 1.40x faster, and fragment mode 1.61x faster than mosdepth v0.3.13, all byte-identical; new default is `processing_threads = 2`
- change `rduckhts_mosdepth()` defaults to `threads = 2` (decompression) and `processing_threads = 2` (parallel contigs) for better out-of-the-box WGS performance
- ship htslib public headers and static library in the installed package under `duckhts_extension/htslib/{include,lib}/`; add `inst/htslib_config.R` (generated from `htslib_config.R.in` at configure time) providing `htslib_cflags()`, `htslib_libs()`, `htslib_rpath()`, and `htslib_version()` for downstream R packages that link against the bundled htslib
- fix `configure.win` to stage htslib headers into `include/htslib/` alongside `lib/`, matching Unix configure
- change bundled `bam_bin_counts(...)` / `rduckhts_bam_bin_counts()` to return a dense fixed-bin layout across each selected contig span, including zero-count bins up to the contig end instead of only observed bins; this gives downstream CNV/sample serializers stable per-contig bin shapes, and the package docs/tests now describe and validate the dense contract
- add `rduckhts_bam_bin_counts()` and bundle native `bam_bin_counts(...)` fixed-width BAM/CRAM binning in the package. The new wrapper exposes `mapq`, `require_flags`, `exclude_flags`, and `rmdup = "none"|"flag"|"streaming"` duplicate handling, always returns per-bin forward/reverse totals, and can add per-bin GC/MAPQ summaries via `stats = "gc"`, `"mq"`, or `"gc,mq"`; bundled extdata now includes the tiny WisecondorX BAM/CRAM fixtures used by the new SQL/R tests, and the package README now includes a native bin-count example
- add `rduckhts_mosdepth()` examples to the package README, including windowed fragment coverage output and preview of the generated BED.gz regions file, and refresh the generated function-catalog text so the packaged mosdepth description matches the current v0.3.13 parity surface
- Expand `rduckhts_mosdepth()` and bundled `duckhts_mosdepth(...)` to cover the pinned local `mosdepth 0.3.13` option surface for indexed BAM/CRAM input: `fragment_mode = TRUE` now matches upstream `--fragment-mode` full-fragment insert coverage for proper pairs, default mode is supported with CIGAR-aware coverage plus mate-overlap correction, `read_groups = "..."` filters RG tags, `min_frag_len` / `max_frag_len` filter absolute template length, and `use_median = TRUE` switches `by = "<window|bed>"` outputs from mean to median; add bundled SQL/R/conformance coverage for BAM and CRAM fast/fragment/default/median cases.
- Expand `rduckhts_mosdepth()` and bundled `duckhts_mosdepth(...)` fast-mode parity with `quantize = "..."`, writing mosdepth-style `.quantized.bed.gz` + CSI output, and add bundled tests for quantized output plus explicit `by = "<bed>"` validation.
- Expand `rduckhts_mosdepth()` and bundled `duckhts_mosdepth(...)` fast-mode parity with `thresholds = "..."` for `by = "<window|bed>"`, writing mosdepth-style `.thresholds.bed.gz` + CSI outputs; also align window/BED mean accumulation and window-region distribution bucketing with upstream mosdepth's current implementation behavior, and add bundled SQL/R/native-conformance coverage for the new outputs.
- Bundle upstream mosdepth edge-case fixtures (`big`, `empty-tids`, `overlapping-pairs`, `ovl`, `nanopore`, and related BED files) in `inst/extdata/` for stronger mosdepth parity testing, and record Brent Pedersen as the original mosdepth author in the package metadata/copyright bundle.
- Expand `rduckhts_mosdepth()` and bundled `duckhts_mosdepth(...)`: the native mosdepth-compatible fast-mode rewrite now accepts indexed CRAM input via `fasta = ...` when required by htslib, and exposes `precision_digits = 2` as an explicit wrapper argument instead of relying on the `MOSDEPTH_PRECISION` environment variable; add bundled BAM/CRAM tests plus explicit precision validation.
- Expand `README.Rmd` with runnable compression/indexing examples covering `rduckhts_bgzip()`, `rduckhts_bgunzip()`, `rduckhts_bam_index()`, `rduckhts_bcf_index()`, and `rduckhts_tabix_index()`, then regenerate the rendered package README outputs.
- Add `decompression_threads` to `rduckhts_bam()` and `rduckhts_bam_multi()`, matching the bundled `read_bam(..., decompression_threads := 2)` SQL parameter. The previous hardcoded htslib worker-thread count is now the documented default, and `0` disables per-file worker threads.
- Speed up bundled zero-column `COUNT(*)` queries across the HTS readers: `read_bam(...)`, `read_bcf(...)`, `read_tabix(...)`, `read_gff(...)`, `read_gtf(...)`, and indexed `read_bed(...)` now use index metadata for full-file count-only scans when DuckDB projects no output columns; `read_fasta(...)` uses `faidx` sequence counts when an index is available and otherwise counts FASTA headers directly; `read_fastq(...)` continues to count raw FASTQ records directly when no projected columns are needed, while preserving paired/interleaved validation errors.
- Add multi-file reading wrappers: `rduckhts_bam_multi`, `rduckhts_bcf_multi`, `rduckhts_fastq_multi`, `rduckhts_fasta_multi`, `rduckhts_bed_multi`, `rduckhts_tabix_multi`, `rduckhts_gff_multi`, `rduckhts_gtf_multi`. Each follows the standard `(con, table_name, files, ..., overwrite)` convention, creates a DuckDB table with a `filename` column, and accepts an optional `.params` data.frame for per-file parameter overrides (e.g. per-sample regions or index paths). File expansion uses DuckDB's `glob()` so S3 URLs work transparently.
- Add bundled `hts_union_query(reader, pattern, params)` SQL scalar macro for pure-SQL multi-file reading via `SELECT * FROM query(hts_union_query('read_bam', '*.bam'))`.
- Clarify the package README's browser/webR documentation: `README.Rmd` now covers the full `Module.duckhtsWasmHttpConfig` parameter set (`headers`, `allowHosts`, `enforceHostAllowlist`, `withCredentials`, `allowInsecureAuth`), explicitly notes that webR consumers can set that config from R via `webr::eval_js()` without editing the host page, and covers practical wasm/browser behaviors such as same-origin setup, CORS requirements, `.csi` to `.tbi` fallback, and non-fatal `Range` warnings under the local `http.server` harness.
- Use one extension-owned Emscripten compatibility header in the package wasm/webR build: `configure` now includes the shared header from `src/include/` via the bootstrapped `inst/duckhts_extension/include/wasm_socket_compat.h` copy, keeping the bundled browser build aligned with the extension sources without changing native package builds.
- Make the bundled wasm extension self-contained with respect to `htslib`: the Emscripten/webR `configure` path now builds only `libhts.a`, links `duckhts.duckdb_extension` directly against that static archive, and no longer relies on runtime loading of bundled `libhts.so*` files in webR/browser environments.
- Add a browser-native wasm `http` / `https` backend in the bundled extension: `src/wasm_http_hfile.c` now registers a synchronous XHR-backed `htslib` scheme handler from the DuckDB extension entry point, so browser wasm builds can read same-origin and CORS-enabled remote HTS URLs without going through libcurl sockets.
- Keep wasm `libcurl` disabled in `configure`: `r-wasm/webr` ships `/opt/webr/wasm/lib/libcurl.a` and the emcc link test against it passes, but libcurl's `connect()` calls from a SIDE_MODULE still trigger a webR Emscripten message-bus error (`resolved is not a function`) on first network use, so the package-owned XHR backend is the supported wasm HTTP path.
- Harden wasm browser HTTP range behavior in the bundled extension: `wasm_http_hfile.c` now caches object sizes from `Content-Range`/`Content-Length`, clamps range requests when size is known, short-circuits reads at/after EOF, and uses a `GET Range: bytes=0-0` fallback for `SEEK_END` size discovery when `HEAD` metadata is unavailable; this avoids cross-origin 416 failures on `.tbi` index EOF probes (including GTEx tabix in webR/browser).
- Harden non-Range wasm/browser HTTP fallback in the bundled extension: when ranged reads receive `200 OK`, `wasm_http_hfile.c` now caches the full object per open handle and serves later reads from that in-memory cache to avoid repeated full downloads, while still emitting one-time warnings when Range is ignored and when large fallback payloads (>=64 MiB) are used.
- Add optional wasm/browser request-header configuration in the bundled extension via `Module.duckhtsWasmHttpConfig`: supports custom headers (including bearer auth), host allowlisting, optional `withCredentials`, and a default HTTPS-only guard that blocks `Authorization` on non-HTTPS URLs unless `allowInsecureAuth` is explicitly enabled.
- Extend `Module.duckhtsWasmHttpConfig` with `enforceHostAllowlist` in the bundled wasm backend: when enabled, requests to hosts outside `allowHosts` are blocked instead of merely omitting configured headers.
- Fix the bundled wasm side-module final link during `configure`: preserve webR/Emscripten `${LDFLAGS}` on the final `duckhts.duckdb_extension` link so the `SIDE_MODULE` settings reach the extension itself, and export `duckhts_init_c_api` explicitly for DuckDB's loader. This fixes webR/browser `rduckhts_load()` failures where DuckDB could not find a usable init export in `duckhts.duckdb_extension`.
- Set the bundled extension metadata platform to `linux_i686_musl` for the Emscripten/webR path in `configure`, matching the platform value you are using for browser-side loading tests.
- Fix Wasm package builds under `rwasm` / r-universe: the package `configure` script now preserves injected `NAME=VALUE` cache overrides, forwards explicit `--build` / `--host` triplets into the vendored `htslib` `./configure`, forwards webR's Emscripten port flags for `zlib`/`bzip2`, seeds wasm-safe Autoconf cache results for `zlib`/`bzip2`/socket probes, injects a tiny Emscripten-only socket compatibility shim for `recv`/`send`/`closesocket`, and disables the optional `htslib` features that are not available in the stock webR/r-universe wasm toolchain (`libcurl`, `S3`, `GCS`, `lzma`, `plugins`); this fixes the original `ac_cv_func_getrandom=no: command not found` failure and the subsequent nested `htslib` cross-compile probe failures without changing native configure behavior.
- Fix bundled wasm extension artifacts: the package/browser wasm build now includes vendored `htslib` in the linked archive, avoiding unresolved symbols such as `bcf_readrec` at `LOAD`.

# Rduckhts 1.1.6-0.0.2 (2026-04-09)

- Fix `test_bam_file_offset`: cast `COUNT(*)` results to `INTEGER` in SQL so the DuckDB driver returns R `integer` rather than `numeric` (BIGINT maps to double in the duckdb R driver), restoring `expect_identical` assertions.

# Rduckhts 1.1.6-0.0.1 (2026-04-09)

- Fix bundled `read_hts_index_spans(...)` / `rduckhts_hts_index_spans()`: the span view now returns real chunk rows from CSI/TBI/BAI indexes, including populated `bin`, `chunk_beg_vo`, `chunk_end_vo`, `chunk_bytes`, `seq_start`, and `seq_end` values instead of placeholder `NA`s; BCF-backed calls also avoid the previous noisy `tbx` probe warning on `.csi` indexes.
- Add `FILE_OFFSET` column to `rduckhts_bam()` / `read_bam(...)`: exposes the BGZF virtual file offset after each record. Zero runtime overhead (macro over already-open struct fields). Enables `ORDER BY FILE_OFFSET` in SQL `LAG()` / `LAST_VALUE()` window functions to reproduce exact BAM file order for streaming deduplication algorithms. Together with the `//` integer-division operator and `LAST_VALUE(... IGNORE NULLS)`, this permits exact replication of WisecondorX's larp/larp2 state machine in pure SQL, confirmed at 0 mismatches across 25,115 non-zero bins on a real NIPT BAM.

# Rduckhts 1.1.5-0.0.1 (2026-04-08)

- Fix bundled `bcftools_liftover(...)` / `rduckhts_liftover()` cache and realignment hardening: per-thread chain/FASTA contexts are now bounded instead of accumulating for the lifetime of worker threads, and scalar left-alignment no longer reuses stale traceback state after failed/empty alignments.
- Fix bundled `read_bam(...)` / `rduckhts_bam()` and `read_bcf(...)` / `rduckhts_bcf()` indexed parallel full scans when headers contain leading empty contigs: contig claiming now retries iteratively instead of recursively, and the BAM reader no longer returns an empty chunk after successfully handing off to the next contig.
- Fix bundled Windows builds under MinGW and Rtools: vendored `htslib` configuration now distinguishes `windows_amd64_mingw` from `windows_amd64_rtools`, keeping the smaller `configure.win`-style library set on MinGW while restoring the fuller static `libcurl` dependency closure needed on Rtools. `CURL_STATICLIB` remains on built objects rather than `./configure` probes.
- Fix bundled Windows `windows_amd64_rtools` builds: the package build now pins `CC`/`AR`/`RANLIB` from `R CMD config`, avoiding mixed compiler/library selection when vendored `htslib` is configured, and keeps the MinGW static-libcurl configuration aligned with Rtools `libcurl.a`.
- Fix bundled `read_bcf(...)` / `rduckhts_bcf()` mapping of fixed-count INFO/FORMAT arrays: exact-cardinality fields such as `Number=2` and `Number=4` now materialize as DuckDB array/list columns instead of silently dropping all but the first value.
- Fix bundled `read_bcf(...)` / `rduckhts_bcf()` handling of string FORMAT lists such as DRAGEN `FORMAT/LAA`: `Number != 1` string FORMAT fields now materialize as `VARCHAR[]` instead of triggering DuckDB internal assertion failures.
- Fix bundled `duckdb_munge(...)` / `rduckhts_munge()` multithreaded FASTA lookups: FASTA index handles are now thread-local and FASTA fetches are synchronized in `munge`, avoiding intermittent `fai_retrieve` failures and aborts when `fasta_ref` is used with `PRAGMA threads > 1`.
- Add `rduckhts_score()`: polygenic risk score computation backed by the `bcftools +score` plugin, supporting GT/DS/HDS/AP/GP/AS dosage modes, all major GWAS summary presets (PLINK, PLINK2, REGENIE, SAIGE, BOLT, METAL, PGS, SSF/GWAS-SSF), GWAS-VCF multi-PRS scoring, p-value thresholding, sample subsetting, and region/filter controls.
- Add `rduckhts_munge()`: GWAS summary statistics normalization backed by `bcftools +munge`, with FASTA reference allele resolution, swap-aware effect/frequency transforms, and METAL meta-analysis column support.
- Add `rduckhts_liftover()`: variant coordinate liftover backed by `bcftools +liftover` using UCSC chain files, with full indel normalization, INFO/END lifting, and MT passthrough.
- Add `rduckhts_bed()` for BED3–BED12 interval files and `rduckhts_fasta_nuc()` for nucleotide composition over BED intervals or fixed-width bins.
- Add compression and index helpers: `rduckhts_bgzip()`, `rduckhts_bgunzip()`, `rduckhts_bam_index()`, `rduckhts_bcf_index()`, and `rduckhts_tabix_index()`.
- Add HTS metadata readers: `rduckhts_hts_header()`, `rduckhts_hts_index()`, `rduckhts_hts_index_spans()`, and `rduckhts_hts_index_raw()`.
- Add quality encoding controls to `rduckhts_bam()` and `rduckhts_fastq()` (`quality_representation`, `input_quality_encoding`) and `rduckhts_detect_quality_encoding()` for heuristic FASTQ encoding detection.
- Add `sequence_encoding := 'nt16'` parameter to `rduckhts_bam()`, `rduckhts_fasta()`, and `rduckhts_fastq()` for raw htslib nt16 sequence output as `UTINYINT[]`.
- Add SAM flag helpers `sam_flag_bits()` and `sam_flag_has()`, CIGAR utility functions, and `is_forward_aligned()`.
- Bundle duckhts 1.1.5 extension.

# Rduckhts 0.1.3-0.0.2.9000

# Rduckhts 0.1.3-0.0.2

- Conditionaly enable plugins in windows

- Updates the configure script to avoid check faillure on CRAN MacOS 

- Update the extension version to 0.1.3

# Rduckhts  0.1.2-0.1.5

- Fixed inadvertant removal of libexec
- Updated the plugin to add header table functions

# Rduckhts 0.1.2-0.1.4

- CRAN Submission

# Rduckhts 0.1.2-0.0.9000

- Different fixes for CRAN submission
    - Updated DESCRIPTION Title/Description formatting and added HTSlib reference.
    - Removed default write paths in bootstrap/build helpers; now require explicit paths.
    - setup_hts_env now accepts an explicit plugins_dir parameter.
    - duckhts_build now accepts a make argument (GNU make required).

- modified configure to attemp to support wasm
- Update bootstrapped extension code to match `duckhts` 0.1.2.
- Add SAMtags + auxiliary tag support (standard_tags, auxiliary_tags).
- Add tabix header/typing options (header, header_names, auto_detect, column_types).


# Rduckhts 0.1.1-0.0.3

- make the build single threaded

# Rduckhts 0.1.1-0.0.3

- misspeling correction

# Rduckhts 0.1.1-0.0.2

- CRAN resubmission: apply DuckDB C API header patch to avoid strict-prototypes warnings.

# Rduckhts 0.1.1-0.0.1

- CRAN Submission

- Bump bundled duckhts extension version to 0.1.1.

- Initial development release.
- Bundles the DuckHTS DuckDB extension and htslib for HTS file readers.
- Adds table-creation helpers for VCF/BCF, BAM/CRAM, FASTA/FASTQ, GFF/GTF, and tabix.
