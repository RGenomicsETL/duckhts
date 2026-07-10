# DuckVEP Model A v2 and allele/transcript context contract

Status: **current implementation contract (M6a, 2026-07-10).** This note fixes the
persistent Model A v2 relations, their provenance and validation rules, and the internal
`allele_transcript_context` consumed independently by SO evaluation and HGVS rendering.
It is the schema companion to `duckvep_m6a_contract_gate.md`. No DuckVEP source should be
folded from `/root/duckvep-c` until its ABI and adapters are reconciled with this contract.

## Decision and boundary

The schema-first next step is correct, but the salvage boundary is narrower than the old
plan implied:

- salvage the consequence-rule generator, classifier, delta, sweep, and property-test
  algorithms;
- redesign the kernel input/result boundary around Model A v2 and the lossless context;
- adapt or rewrite the old model reader and annotation glue rather than lifting them under
  a new prefix; and
- treat the current fused Model A row as a derived execution view, never the persistent
  source of truth.

The old ABI explicitly describes HGVS as a renderer over selected compact consequence rows
and stores only scalar cDNA/CDS/protein positions
(`/root/duckvep-c/kernel/include/duckvep_kernel.h:207-242`). The old model reader is bound
to one nested transcript row plus an optional `cds_seq`
(`/root/duckvep-c/src/include/duckvep_model_reader.h:5-32`). Both encode the architecture
M6a rejects. The ABI's `uint16_t` REF/ALT lengths and exon counts are also transport limits,
not biological facts. The v2 boundary uses wide lengths/counts and permits narrowing only
behind a manifest-bound, checked algorithm adapter; no old limit becomes an undocumented
rejection.

Model A v2 is a **relation pack**, not a new private cache format. DuckDB tables or Parquet
are the interchange. A fused nested view remains useful at execution time, but it is built
from the normalized relations below and is disposable.

## Exact normative pins

Keep compatibility code, annotation data, and nomenclature versions separate. A single
column named `release` is insufficient.

- **VEP compatibility target:** Ensembl VEP `release/116.0`, commit
  `57ea5c52340acc1f156267f810ad162e26597082`. The upstream tag archive has SHA-256
  `618a4b6d37efbe0968d7ad1115bf6b712f8537c4697659be6c41580708eb5167`, as recorded by
  the installed Bioconda recipe; the installed
  `ensembl-vep-116.0-pl5321h2a3209d_0` package record pins distribution-artifact SHA-256
  `a58cfc1209f55b0a6b63f56b3808b963a36290fb0c5b7f719c0205725fe0cdd3`.
- **Ensembl 116 API/schema pins embedded by that VEP distribution:** core
  `c0cf13daa961d80584bad797b2eb0ff3a7500ef3`, variation
  `2fb834b987ede3824e200197a838ce11e91aeb4b`, funcgen
  `90049ea7ee4d8ae3a6d298dca46d6c6ab20538c4`, compara
  `46f80af9c3864fc557dca90a7eeb1e221e2bad70`, and ensembl-io
  `6afb5dc27a5ae6881d75959153fbe6e9a4a7e788`. The exact core DDL is
  [`sql/table.sql`](https://github.com/Ensembl/ensembl/blob/c0cf13daa961d80584bad797b2eb0ff3a7500ef3/sql/table.sql),
  SHA-256 `ec7bb8dd2fcd6a7012bfd67aff8065b912915d541450fa4f3f05334a92944e8c`.
- **HGVS normative snapshot:** HGVS Nomenclature **21.1.4** (2026-05-11), tag commit
  [`6f85311989e76ead95d3547828f97ebaa3802e35`](https://github.com/HGVSnomenclature/hgvs-nomenclature/tree/6f85311989e76ead95d3547828f97ebaa3802e35),
  with no named extensions in the first claim. VEP 116 remains the compatibility oracle;
  it does not override the normative HGVS snapshot.

The local `.sync/ensembl-vep` mirror and `/root/ensembl-variation` checkout are release 115.
They remain useful for the audited 115-to-116 delta, but are not the 116 source oracle.

An annotation model has its own source identity. For example, the local VEP-116 cache
receipts identify GRCh38 as GENCODE 50 / GRCh38.p14, while the GRCh37 merged cache identifies
GENCODE 19 plus RefSeq 105.20220307 / GRCh37.p13. `vep_compat_release = 116` does not make
those annotation sets the same model.

The existing `/root/duckvep/data/cache/homo_sapiens.116.38/ensembl.duckdb` is useful staging
evidence, not a valid v2 pack: it omits required source joins such as `object_xref`,
`external_db`, and `seq_region_attrib`, and contains none of the five sequence-state
receipts. Do not bless it by adding a manifest around an incomplete schema.

## Persistent relation pack

All source objects use `(model_id, source_namespace, source_internal_id)` as identity.
Stable ID plus version is metadata, not a key: the Ensembl DDL does not make it unique, and
real release-116 data contains repeated translation stable-ID/version pairs. Dense runtime
indices are derived after validation and never escape as biological identifiers.

Every `*_key` below is an opaque pack-local surrogate with a checked bijection to
`(source_namespace, source_internal_id)` for source objects. It exists to make Parquet joins
compact; it is not an accession and must not be reconstructed from a stable ID. Relations
without a source row ID use the duplicate-aware raw-row identity specified below.

### Manifest and artifacts

`model_manifest` has exactly one row:

```text
model_manifest(
  model_id, format_version,
  compatibility_target, compatibility_commit,
  species, taxonomy_id, assembly_name, assembly_accession,
  selection_policy_version, canonicalization_version,
  created_at
)

model_source(
  model_id, source_namespace, provider, source_kind, selection_role,
  annotation_set, annotation_source_release, source_database,
  schema_repository, schema_commit, schema_sha256, source_priority
)

model_selection_audit(
  model_id, source_namespace, object_kind,
  decision, reason_code, row_count, rejected_rows_artifact
)

model_build(
  model_id, model_build_id,
  importer_repository, importer_commit, importer_sha256,
  exporter_repository, exporter_commit, exporter_sha256,
  invocation_sha256, environment_sha256, created_at
)
```

`model_source` is one-to-many: a merged cache, GENCODE overlay, MANE crosswalk, and RefSeq
source are distinct namespaces even when they contribute to one assembly. The namespace
also fixes collision and precedence policy; it is never inferred from accession spelling.

`model_id` is a SHA-256 logical identity over the canonical identity fields in
`model_manifest`, canonical `model_source` rows, required source-input receipt tuples, and
canonical logical digests of the biological relations, excluding repeated `model_id`
columns and timestamps. Importer implementation and output physical layout do not change
`model_id` when the logical model and source receipts are identical.

`model_build_id` is a separate SHA-256 build receipt over `model_id`, importer/exporter code
receipts, normalized invocation, environment receipt, and required input receipts, excluding
timestamps and generated output bytes. This distinguishes two production paths without
pretending they are different biological models. Neither ID is a random UUID or a Parquet-
file checksum. Logical row serialization, duplicate handling, NULL handling, text encoding,
sort keys, included relation/role set, and hash algorithm are versioned by
`canonicalization_version`.

`model_artifact` records every input and generated receipt:

```text
model_artifact(
  model_id, model_build_id, role_class, role,
  source_namespace, logical_name, locator,
  sha256, byte_count, row_count, required
)
```

Required roles include core DDL, each database/dump part, reference FASTA and indices,
cache `info.txt` where applicable, GENCODE/MANE or RefSeq overlays, importer/exporter
SQL/code, oracle export, and the final relation artifacts. Generated artifact checksums are
recorded after both IDs exist and feed back into neither ID, avoiding a digest cycle when a
Parquet file contains those IDs. API modules used by the sequence-state exporter are
software-input artifacts; modules loaded by a conformance oracle belong to that run's
manifest. Both record paths and SHA-256 hashes. `vep --version` is not a receipt.

The manifest also pins model selection. Source kind (`core`, `refseq`, or `merged`), VEP
options, coordinate-system/top-level policy, GENCODE subset, and every inclusion/exclusion
reason are data. Preserve excluded rows in staging/audit counts; do not silently discard
scaffolds, LRGs, artifacts, readthrough transcripts, missing stable IDs, or source-specific
duplicates. `model_selection_audit` records mutually exclusive decision counts, and any
row-level rejection ledger is a checksummed artifact rather than an untracked log.

### Physical type and constraint rules

The pseudotypes above map to DuckDB/Parquet as follows; adapters do not choose narrower
types opportunistically:

- `model_id`, `model_build_id`, and SHA-256 fields are lowercase 64-character hex
  `VARCHAR`; import rejects malformed or non-canonical spellings.
- Opaque pack-local `*_key` values are `UBIGINT`, assigned deterministically after sorting
  the authoritative source identities. `source_internal_id` is lossless canonical `VARCHAR`
  so numeric database IDs and string-backed overlay/cache IDs share one representation.
- Source coordinates and byte counts are `BIGINT` with explicit nonnegative/positive
  checks as appropriate. Ranks, ordinals, priorities, and stable-ID versions are `INTEGER`;
  row counts are `UBIGINT`. The C adapter narrows only after a checked bound.
- Sequence payloads are `BLOB`; accessions, alleles, attribute values, enum/status names,
  and provider-specific accession versions are `VARCHAR`. Persistent enums are named
  allowlists, never private numeric constants.
- Absence is SQL `NULL`, never zero, an empty string, or an all-zero key. Empty sequence is
  retained only when biologically distinct from absence.
- Primary identities are `(model_id, source_namespace, *_key)` for source-object relations
  and the full owner/raw-tuple/duplicate ordinal for attribute rows. All foreign-key,
  cardinality, range, and enum checks run explicitly even when a Parquet reader cannot
  enforce SQL constraints.
- Persistent row order has no meaning. Every derived array names its sort key; no importer
  may inherit dump order except through an explicit receipt-backed ordinal whose semantics
  are documented.

The normalized pack contains no persistent nested `LIST`/`STRUCT` transcript object. Only
the disposable execution view may use that shape.

### Source-faithful geometry and identity

The minimum normalized geometry is:

```text
model_seq_region(
  model_id, source_namespace, seq_region_key,
  source_internal_id, name, length,
  coord_system_name, coord_system_version,
  assembly_accession, reference_accession, reference_version,
  sequence_sha256, is_circular
)

model_gene(
  model_id, source_namespace, gene_key, source_internal_id,
  stable_id, stable_id_version, seq_region_key,
  start1, end1, strand, source, biotype, is_current,
  canonical_transcript_key
)

model_transcript(
  model_id, source_namespace, transcript_key, source_internal_id,
  stable_id, stable_id_version, gene_key,
  seq_region_key, start1, end1, strand, source, biotype, is_current,
  canonical_translation_key, is_canonical
)

model_exon(
  model_id, source_namespace, exon_key, source_internal_id,
  stable_id, stable_id_version, seq_region_key,
  start1, end1, strand, phase, end_phase
)

model_transcript_exon(
  model_id, source_namespace, transcript_key, exon_key, exon_rank,
  raw_cdna_start1, raw_cdna_end1
)

model_translation(
  model_id, source_namespace, translation_key, source_internal_id,
  transcript_key, stable_id, stable_id_version,
  start_exon_key, start_offset1, end_exon_key, end_offset1,
  raw_cdna_coding_start1, raw_cdna_coding_end1,
  edited_cdna_coding_start1, edited_cdna_coding_end1,
  start_phase_padding, codon_table, is_canonical
)
```

Coordinates in persistent source relations are 1-based and inclusive because that is the
Ensembl source contract. `start_offset1` and `end_offset1` are exon-relative in the exon's
transcript orientation (on the minus strand, offset 1 is at the genomic exon end), not
forward-genomic offsets. Exons are N:M and `exon_rank` is transcript 5'-to-3' order. Import
validation requires one exon per rank, contiguous ranks beginning at 1, one consistent
strand/sequence region, and non-overlapping geometry in transcript order.

Translation is 0:N per transcript. Preserve every source translation; the execution view
selects `canonical_translation_key`. Do not collapse alternatives or key them by stable ID.
Gene, transcript, and translation rows retain their source-internal chain; symbols and other
display identifiers remain xrefs rather than being promoted into keys. Canonical transcript
and translation pointers must resolve inside the same source namespace; redundant
`is_canonical` fields are checked against those pointers.

`start_phase_padding` is explicit. The translatable CDS includes that many leading `N`
bases, matching the Ensembl API. The current GFF/FASTA conformance helper omits this padding
while the old kernel projection adds it; a phase-1/2 witness is therefore mandatory before
the helper is reused.

Codon table comes from `seq_region_attrib` code `codon_table`, defaulting to table 1 when
absent. Do not infer table 2 from `M`, `MT`, or `chrM` spelling. The first engine claim
accepts only tables 1 and 2 and inventories/rejects anything else.

### Cardinality ledger

| Edge | Cardinality | Import invariant |
| --- | --- | --- |
| model → source namespace | 1:N | Namespace is unique within `model_id`; priority/collision policy is complete. |
| sequence region → gene/transcript/exon | 1:N | Geometry stays in one declared coordinate system and assembly. |
| gene → transcript | 1:N | Orphan transcripts are admitted only by an explicit source policy and audit count. |
| transcript ↔ exon | N:M | Membership rank is unique/contiguous; exon identity is not duplicated per transcript. |
| transcript → translation | 0:N | All translations survive; canonical pointer is NULL or resolves to exactly one owned row. |
| source object → raw attribute | 0:N multiset | Exact duplicates survive through `duplicate_ordinal`. |
| transcript/translation → typed edit | 0:N ordered | Owner, basis, parse status, and replay order are explicit. |
| object → xref | 0:N | Display choice and source namespace survive; stable ID is not substituted for an xref. |
| transcript → mature-miRNA span | 0:N | Spans remain separate rows in the declared coordinate basis. |
| transcript/translation → sequence state | 0:1 per role | Role-appropriate NULL keys and referenced blob hash/length are validated. |

### Raw attributes, typed edits, and annotations

Owner table is semantic. Preserve raw attributes before deriving flags or edits:

```text
model_attribute_type(model_id, source_namespace, attrib_type_key,
  source_internal_id, code, name, description)

model_seq_region_attribute(model_id, source_namespace, attribute_key,
  seq_region_key, attrib_type_key, duplicate_ordinal, source_row_locator, value)

model_gene_attribute(model_id, source_namespace, attribute_key,
  gene_key, attrib_type_key, duplicate_ordinal, source_row_locator, value)

model_transcript_attribute(model_id, source_namespace, attribute_key,
  transcript_key, attrib_type_key, duplicate_ordinal, source_row_locator, value)

model_translation_attribute(model_id, source_namespace, attribute_key,
  translation_key, attrib_type_key, duplicate_ordinal, source_row_locator, value)
```

The Ensembl owner-attrib tables do not provide a row ID. `attribute_key` therefore derives
from the complete raw tuple plus `duplicate_ordinal`; exact duplicate rows remain a multiset
instead of collapsing. `source_row_locator` is a receipt-backed dump-part/line or byte
locator for diagnosis and is excluded from logical identity. Attribute code, name, and
description come from the pinned attrib-type row, not an unversioned application dictionary.

Typed edit relations reference the raw row:

```text
model_transcript_edit(
  model_id, source_namespace, edit_key, transcript_key, attribute_key,
  code, coordinate_basis, start1, end1,
  basis_ref_seq, preapply_ref_seq, alt_seq, apply_ordinal, status
)

model_translation_edit(
  model_id, source_namespace, edit_key, translation_key, attribute_key,
  code, coordinate_basis, start1, end1,
  basis_ref_seq, preapply_ref_seq, alt_seq, apply_ordinal, status
)
```

The nucleotide allowlist is transcript-owned `_rna_edit`. Its coordinates address the raw
spliced transcript and it is applied before CDS extraction. The peptide allowlist is
translation-owned `initial_met`, `_selenocysteine`, `amino_acid_sub`, and
`_stop_codon_rt`; those coordinates address the unedited peptide and are applied after
translation. A code with the wrong owner is preserved and inventoried, not reinterpreted by
name. In particular, translation-owned `_rna_edit` rows observed in human release 116 are
not automatically applied by the Ensembl API and must remain explicit audit rows.
`status` is a named allowlist (`applied`, `source_ignored`,
`quarantined_ambiguous_order`, or `unsupported_code`); only `applied` mutates a sequence.

Edit coordinates are 1-based inclusive; insertion is `start1 = end1 + 1`; deletion has an
empty alternate. `basis_ref_seq` is derived from the named unedited coordinate basis;
`preapply_ref_seq` records the bytes present immediately before this edit in the ordered
replay. They differ only when earlier-applied edits overlap. Recognized malformed rows fail
the import.

Ensembl applies edits in descending start order, but source attribute retrieval has no
defined row order. Equal-start edits therefore have no source-defined tie order. The
exporter may supply an observed, receipt-backed `apply_ordinal`; otherwise equal-start rows
are quarantined. Distinct-start overlapping edits retain the defined descending-start order
and both reference snapshots; an importer that does not replay that interaction quarantines
it. Never invent a hidden code-priority tie break.

Transcript `_rna_edit` can change coding boundaries. Preserve `_transl_start` and
`_transl_end` raw attributes and materialize the edited cDNA coding bounds above. An importer
that neither models those overrides nor proves them absent is invalid. Multiple conflicting
override values are ambiguous source state and quarantine the transcript; do not reproduce
the API's accidental first-row choice.

MANE, GENCODE, CCDS, partiality, readthrough, and similar facts retain their values and
source rows. MANE values may be RefSeq accessions; a boolean bit loses the crosswalk.
`tx_flags` remains a compact runtime predicate cache only. The execution adapter derives it
from biotype, translation presence, typed edits, and named attributes under a versioned bit
mapping, then checks it against any cached value.

Many-valued transcript features remain relations:

```text
model_mature_mirna(
  model_id, source_namespace, feature_key, transcript_key,
  cdna_start1, cdna_end1, coordinate_basis, source_attribute_key
)
```

Do not reduce mature-miRNA spans to one transcript flag: release-116 human data includes
multiple spans on many transcripts.

### External identifiers and sequence identity

```text
model_object_xref(
  model_id, object_source_namespace, xref_source_namespace,
  xref_key, source_internal_id,
  object_kind, object_key,
  external_db_name, external_db_release,
  accession, accession_version, display_label,
  is_display_xref, linkage_annotation, info_type,
  identity_status, verified_sequence_role, verified_sequence_sha256
)
```

Gene, transcript, translation, and sequence-region aliases are many-valued. The object and
xref namespaces remain separate so an overlay cannot masquerade as the owning annotation
source; the source display-xref choice is retained rather than re-selected lexically.
`identity_status` is `unverified`, `exact`, or `mismatch`, and
`verified_sequence_role` names the raw/edited transcript or peptide basis that was compared.
Only an exact sequence-hash comparison promotes a MANE/RefSeq/other crosswalk to sequence
identity. A MANE tag alone is selection metadata. HGVS generation requires an accession
**and version** whose sequence hash matches the reference bytes used.

### Five sequence states

Sequence bytes are content-addressed and deduplicated:

```text
model_sequence_blob(sequence_sha256, alphabet, byte_count, sequence_bytes)

model_sequence_state(
  model_id, source_namespace, transcript_key, translation_key,
  role, sequence_sha256, byte_count, status
)
```

Required roles are:

1. transcript `raw_spliced` — transcript-oriented genomic exon sequence, before RNA edits;
2. transcript `edited_spliced` — after transcript-owned nucleotide edits;
3. translation `translatable_cds` — extracted after RNA edits and coding-boundary overrides,
   including positive start-phase `N` padding;
4. translation `peptide_unedited` — the Ensembl translation result after codon translation,
   terminal-stop handling, and automatic initial-methionine handling, before translation
   SeqEdits; and
5. translation `peptide_final` — after the peptide-scope edits.

The API has one transcript `edits_enabled` switch that controls both transcript RNA edits
and the later `Translation::modify_translation` call; it does not directly expose state 4.
The pinned exporter therefore obtains state 1 from a fresh edits-disabled object, states 2/3
from a fresh edits-enabled object, translates state 3 through the exact pre-
`modify_translation` codon/terminal-stop/initial-Met logic for state 4, and obtains state 5
from the normal edited `Transcript::translate` path. That exporter logic is a checksummed
software artifact and is tested against transcripts with each peptide edit.

Hash exact uppercase API-export bytes without a terminator or newline; the role declares
DNA versus peptide alphabet. Preserve NULL versus empty: a noncoding transcript has no
translation, not an empty peptide that happens to hash. Missing exon sequence that the API
would replace with `N` plus a warning is a provider error unless the model manifest
explicitly admits it.

Export each state from fresh API objects or with explicit cache invalidation. Reusing an
object while toggling edit state can return a cached translation and bless the wrong hash.

## Derived execution view

The runtime preparation step may construct a shape-aligned view with one selected
translation per transcript, nested exons in transcript order, named sequence-blob IDs,
typed edit slices, and derived `tx_flags`. It may then bulk-copy that bounded model into an
immutable SoA and sequence provider.

This view is not persisted as truth, and the old model reader cannot consume it unchanged.
The v2 reader must:

- validate `model_id` equality across every relation;
- select exactly one `model_build_id` and verify every required physical-artifact receipt;
- derive dense indices and offsets with checked narrowing;
- preserve source keys for diagnostics;
- build raw-to-edited transcript coordinate maps from nucleotide edits;
- keep sequence blobs pooled rather than copying full strings per row;
- retain edit IDs and reference hashes for context provenance; and
- give every worker its own mutable sequence-provider handle/cache.

## Immutable allele/transcript context

The context is an internal ephemeral value for one biallelic ALT and one transcript. SO and
HGVS consume it as siblings. Existing `duckvep_event_t`, `duckvep_coding_context_t`, and
`duckvep_consequence_t` remain useful derived views; none is the canonical intermediate.

### Coordinate and edit representation

Internal C coordinates are **0-based half-open**. A pure insertion is the empty span
`[p,p)`; it is never disguised as a one-base point. Conversion to source 1-based coordinates
is confined to adapters and renderers.

The per-allele base preserves:

- original VCF contig, POS, REF, ALT, ALT ordinal, declared kind, variant ID, and supplied
  end;
- the storage-left-normalized VariantKey identity, if present;
- the prefix/suffix-trimmed semantic replacement envelope;
- a VEP-116 compatibility view of the differing regions (`vep_diff_runs[]`), tagged with
  the exact `VariationFeatureOverlapAllele.pm::_get_differing_regions` algorithm and based
  on VEP's minimized overlap-allele strings, including its unequal-length byte-XOR
  behavior; and
- immutable normalization views, each tagged with algorithm, axis, reference identity/hash,
  rotated alleles, signed displacement, and `exact`/`window_exhausted`/`unavailable` status.

The semantic envelope is the lossless edit. The compatibility runs are a separate derived
view, not an attempt to align arbitrary delins. One envelope alone is insufficient for VEP
parity: VEP 116 evaluates separated differing regions independently, so an MNV with a
retained island can put different runs in different splice compartments. The exact VEP
algorithm, including empty-insertion and unequal-length witnesses, is frozen in tests rather
than silently replaced with a more attractive aligner. Runs are stored 0-based half-open;
the compatibility adapter explicitly converts VEP's 0-based inclusive `{s,e}` pairs,
including its `e = -1` empty-insertion sentinel.

The per-transcript context adds a piecewise projection:

```c
typedef struct duckvep_span0 {
	uint64_t begin;
	uint64_t end;
} duckvep_span0_t;

typedef struct duckvep_projection_piece {
	duckvep_span0_t genomic;
	duckvep_span0_t transcript_raw; /* valid for exonic pieces */
	uint32_t feature_idx;
	uint64_t semantic_ref_offset; /* forward-genomic semantic REF */
	uint8_t kind; /* exon, intron, outside */
} duckvep_projection_piece_t;

typedef struct duckvep_tx_anchor {
	uint64_t anchor_cdna0;
	int64_t intron_offset;
	uint32_t exon_idx;
	uint8_t region;
} duckvep_tx_anchor_t;
```

The complete private shape is component-based rather than one object-like bag:

```text
allele_transcript_context(
  allele_base, transcript_reference_bundle,
  projection_pieces[], endpoint_anchors[2],
  transcript_raw_edit, transcript_reference_edit, translation_cds_edit,
  reference_peptide_diff, alternate_peptide_diff,
  alternate_cds, alternate_peptide,
  applied_model_edit_ids[],
  mapping_status, partiality_flags, provenance_flags
)
```

`transcript_reference_bundle` identifies the model/provider/source/assembly, transcript and
translation source keys, accession/version aliases, all five sequence hashes, raw
attributes, and typed model edits. It points into immutable pooled model storage; it does
not copy cDNA or peptide strings per candidate.

Projection pieces partition the semantic REF span into explicit exon, intron, and
outside-transcript segments. They stay ordered in transcript 5'-to-3' order. A partial
transcript overlap keeps the unclipped genomic span plus clip status; VEP-116-compatible
clamping is a derived rendering view. A wholly outside candidate keeps an explicit outside
piece so upstream/downstream consequences remain evaluable; it simply has no fabricated
transcript coordinate. For a pure insertion there is no REF-consuming piece, so the two
endpoint anchors instead describe the left and right genomic flanks and retain both adjacent
compartments at an exon/intron boundary. Intronic endpoints retain nearest-exon anchor plus
signed offset; the equidistant tie rule must be pinned by a VEP witness before implementation.

Keep these axes distinct:

- forward genomic reference;
- raw spliced transcript;
- RNA-edited transcript reference;
- phase-padded translatable CDS; and
- final post-edit peptide.

HGVS `c.` coordinates are not phase-padded CDS indices. Raw and edited transcript
coordinates may diverge after length-changing `_rna_edit` rows. An allele overlapping a
model edit records an explicit conflict/resolution status; it is never silently projected
through the edit.

Preserve input `GIVEN_REF`, assembly/provider `USED_REF`, transcript-oriented used REF/ALT,
and their validation statuses independently. A valid genomic REF can disagree with a
RefSeq transcript without making the genomic record invalid.

### Normalization is a view

No normalizer mutates the allele base. At minimum retain separate views for storage-left,
semantic/unshifted, VEP compatibility shifting, HGVS `g.` 3'-shift, and HGVS `c./n.`
3'-shift. Genomic and transcript HGVS shifts are different algorithms; negative-strand
`c./n.` moves opposite the genomic direction. Window exhaustion triggers a refetch or an
incomplete status, never a falsely maximal result. HGVS views never feed back into
VariantKey.

### Ownership and lifetime

The context is published as a `const` borrowed view by a private mutable builder.
Projection pieces, reverse-complemented alleles, alternate CDS, and peptide bytes live in a
non-relocating per-worker arena. Immutable reference bytes live in the model pack/blob
provider. The context expires when the tile/workspace resets, and SO/HGVS consume it
synchronously before that reset. Public result rows contain no context pointers. Any later
deferred renderer needs an explicit deep-clone API; accidental tile borrowing is forbidden.

## M6b capability and conformance manifest

Zero discord is meaningful only over a named slice. One machine-readable run manifest must
bind all of the following:

```text
contract_version, model_id, model_build_id,
engine_commit, engine_binary_sha256, rule_program_sha256,
duckdb_version, runtime_library_hashes,
compatibility_target, oracle_container_digest, oracle_command_sha256,
loaded_module_hashes, oracle_lane, assembly, annotation_set, source_namespaces,
input_fixture_hashes,
generator_commit, random_seed, sample_count,
variant_shapes, transcript_states, supported_terms,
max_ref_bases, max_alt_bases, max_affected_span_bases,
normalization_policy_versions, distance_policy_version,
upstream_bases, downstream_bases,
required_strata, explicit_error_strata, intentional_divergence_strata, backlog_strata,
errata_ledger_sha256,
output_equality_rule, expected_discord_count, skip_allowed
```

The initial parity slice is exact, biallelic, non-symbolic SNV/MNV/INS/DEL/delins with one
ALT per row. It excludes spanning `*`, BND/SV/CNV, uncertain coordinates, and haplotypes.
The semantic operation is derived from trimmed REF/ALT lengths; the declared transport kind
is validated, not trusted. Every run supplies finite numeric allele/span bounds; a manifest
with an absent or implied size limit is invalid. Those bounds are capability evidence, not
an excuse for a fixed internal allele buffer.

Required strata include both strands; coding and noncoding transcripts; UTR, exon, intron,
splice, upstream, and downstream topology; phase 0/1/2; complete and partial CDS; repeats;
codon tables 1 and 2; separated MNV differing runs; every supported model-edit code; and at
least one explicit transcript-reference mismatch policy witness. Each stratum is `parity`,
`explicit_error`, `intentional_divergence`, or `backlog` and carries a required minimum
witness count. Every intentional divergence names a checksummed ERRATA witness; it cannot
be introduced by a wildcard. A rejection witness is not an oracle skip, and an empty
required stratum fails the run.

Lane 1 (`--gff`) covers synthetic geometry. Lane 2 (pinned cache/DB) covers curated edits,
source identities, and rare transcript states. A parity row fails on any SO-set discord or
missing oracle. Equality covers the complete emitted allele/transcript-pair set and each
pair's SO-term set; filtering to shared pairs is forbidden. The primary sampling unit for
statistical audit is variant/locus; pair rates remain descriptive. The manifest selects one
authoritative run and never aggregates every dated output in a directory.

The rare-state fuzzer should compare context receipts as well as final SO sets: raw/shifted
coordinates, differing runs, projection pieces, given/used alleles, cDNA/CDS/peptide spans,
model-edit IDs, and normalization offsets. Matching SO sets can hide compensating context
bugs. VEP does not expose every internal component, so the receipt oracle is a structurally
independent slow scalar/base-walk projector (not a call into the production projection
helper) plus source sequence/edit receipts, augmented only by diagnostic fields VEP actually
emits; do not invent unobservable VEP ground truth.

### Rare-state generation discipline

The generator's state vector is explicit: allele shape and lengths, strand, feature topology
and exact boundary distance, phase, partiality, codon table, repeat class, transcript source,
model-edit code/owner/collision, and reference-match status. Deterministic boundary witnesses
cover every required cell; stratified random generation then explores interactions instead
of drawing mostly ordinary intronic SNVs. The manifest records requested and realized counts
for each cell, including rejected generation attempts.

Every random failure is persisted with the generator commit/seed and minimized only by a
shrinker that preserves the failing state labels, oracle availability, and context-receipt
discord. Both original and minimized witnesses are checksummed. Statistical intervals
describe residual sampling uncertainty; they never convert an observed parity discord,
missing oracle, or empty required cell into a pass.

## Required validation before M6b can claim the slice

1. Export all five sequence states from the pinned API and compare logical length/SHA-256
   with the importer, stratified by edit code, phase, partiality, codon table, source, and
   assembly.
2. Validate source cardinalities, contiguous exon ranks, canonical-translation selection,
   attribute parse status, edit order, xref multiplicity, and inclusion/exclusion counts.
3. Brute-force piecewise projection against a base-by-base transcript walk over random exon
   models on both strands.
4. Prove application equivalence between genomic and projected edits when exact; report a
   conflict for overlaps with model edits.
5. Exercise every insertion around exon/intron and transcript boundaries, separated MNV
   runs, unequal-length delins compatibility runs, intronic anchor ties, partial-overlap
   clamps, and wholly outside upstream/downstream candidates.
6. Cover phase 0/1/2 and incomplete CDS, including a test that fails against the current
   unpadded GFF/FASTA helper.
7. Mirror genome/model/alleles by reverse complement and strand flip; transcript-axis
   context and SO must remain equivalent.
8. Check padded-VCF metamorphs: raw provenance may differ, but semantic edit, projection,
   SO, and application result must agree.
9. Check normalization idempotence and application equivalence independently for `g.` and
   `c./n.` repeats; VariantKey must remain unchanged.
10. Include long transcripts above 4096 nt, high exon counts, alleles at every declared
    adapter/capability bound plus one over-limit explicit-error witness, and ASan lifetime
    tests that poison/reset tile memory immediately after allowed consumption.

M6a fixes this contract and the pins. M6b implements the importer, frozen oracle export,
context builder, and SO slice, and earns the hash/conformance evidence above. This removes
the previous circular requirement that a no-code M6a milestone already compare a not-yet-
implemented importer.
