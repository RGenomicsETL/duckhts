# DuckVEP annotation model and context contract

Status: **current implementation contract (2026-07-10).** The executable schema is
`scripts/duckvep_model_schema.sql`; `scripts/duckvep_model.R` initializes and validates
relation packs. Private C components live under `src/duckvep/`.

## Compatibility pins

- Ensembl VEP `release/116.0`, commit
  `57ea5c52340acc1f156267f810ad162e26597082`.
- Ensembl core `c0cf13daa961d80584bad797b2eb0ff3a7500ef3`; pinned `sql/table.sql`
  SHA-256 `ec7bb8dd2fcd6a7012bfd67aff8065b912915d541450fa4f3f05334a92944e8c`.
- Ensembl variation `2fb834b987ede3824e200197a838ce11e91aeb4b`.
- HGVS Nomenclature 21.1.4, commit
  `6f85311989e76ead95d3547828f97ebaa3802e35`.
- Initial public scope: human GRCh37 and GRCh38, codon tables 1 and 2.

Compatibility release, annotation release, assembly, source database, and nomenclature
version are independent fields. A VEP version string is not a source receipt; manifests
record commits and SHA-256 hashes.

## Architecture boundary

- The annotation model is a normalized DuckDB/Parquet relation pack, not a private cache
  format. Nested transcript rows are disposable execution views.
- SO consequence evaluation and HGVS rendering consume the same immutable
  `allele_transcript_context` as siblings. HGVS is never rendered from a compact consequence
  row.
- The old `/root/duckvep-c` model reader, annotation glue, compact kernel ABI, and scalar
  projector are rewrite evidence, not production components. There is one piecewise
  projector and an independent scalar/base-walk test oracle.
- Allele, transcript, exon, and sequence lengths are wide and checked. Limits from the old
  `uint16_t` ABI or the duckhgvs 4096-byte buffers are not biological limits.
- `tx_flags` is a derived execution cache. It is not persistent model truth and never stores
  edit payload.

Before translating or bundling Ensembl/VEP source text, add the pinned Apache-2.0 license
and attribution under `third_party/licenses/` and `r/Rduckhts/inst/COPYRIGHT`. The current
differing-region primitive is an independent behavioral rewrite and copies no upstream text.

## Relation pack

The pack contains exactly these 26 normalized relations:

| Group | Relations |
| --- | --- |
| Manifest and receipts | `model_manifest`, `model_source`, `model_selection_audit`, `model_build`, `model_artifact` |
| Geometry and identity | `model_seq_region`, `model_gene`, `model_transcript`, `model_exon`, `model_transcript_exon`, `model_translation` |
| Raw attributes and typed edits | `model_attribute_type`, `model_seq_region_attribute`, `model_gene_attribute`, `model_transcript_attribute`, `model_translation_attribute`, `model_transcript_edit`, `model_translation_edit`, `model_mature_mirna` |
| External identifiers | `model_external_db`, `model_xref`, `model_object_xref`, `model_seq_region_synonym`, `model_xref_identity` |
| Sequence states | `model_sequence_blob`, `model_sequence_state` |

`format_version = 2` is data in `model_manifest`; filenames and SQL namespaces are not
versioned. Schema evolution changes the manifest version and migration rules rather than
creating parallel implementations.

## Identity and physical rules

- A pack has exactly one `model_manifest` row and at least one `model_source` and
  `model_build` row.
- Source objects use `(model_id, source_namespace, source_internal_id)` as identity. Stable
  ID and version are nullable metadata, never keys.
- Pack-local `*_key` values are positive `UBIGINT`s with a checked bijection to source
  identity. Dense runtime indices never escape as biological identifiers.
- Source coordinates are 1-based inclusive `BIGINT`. Internal C coordinates are 0-based
  half-open. Conversion occurs only at adapters/renderers.
- Absence is `NULL`, never zero, empty text, or an all-zero key. Empty sequence remains
  distinct from absent sequence.
- Text is exact UTF-8 `VARCHAR`; sequence and edit payloads are `BLOB`; hashes are lowercase
  64-character hexadecimal SHA-256.
- `external_db`, `xref`, `object_xref`, and `seq_region_synonym` source IDs remain distinct.
  Sequence-region synonyms may have no external database.
- Display-xref pointers on genes and transcripts resolve through the normalized xref chain
  and agree with `model_object_xref.is_display_xref`.
- `model_xref_identity.verified_sequence_sha256` is the external accession sequence hash.
  `exact` means it equals the named present model sequence state; `mismatch` means it does
  not. Gene xrefs cannot claim sequence identity.
- `model_sequence_blob` has no `model_id`. Its primary key is
  `(sequence_sha256, alphabet)` because the same bytes may be valid DNA and peptide.

`model_id` and `model_build_id` are intended to be canonical SHA-256 identities. Canonical
row serialization and golden vectors are not yet specified, so current tooling must report
identity authentication as unavailable. Structural or model-candidate validation is not
proof that either identifier was recomputed.

## Geometry and edit invariants

- Transcript/exon membership is N:M. Its source-faithful identity includes transcript,
  exon, and rank. Rank is contiguous from 1 for an admitted transcript; rare repeated
  memberships remain representable until selection policy accepts or quarantines them.
- Exons in one transcript share strand and sequence region and are ordered 5′ to 3′.
- A transcript owns zero or more translations. Canonical transcript/translation pointers
  are nullable or resolve within the same source namespace; every source translation is
  preserved.
- Translation offsets are exon-relative in transcript orientation. Start-phase padding is
  explicit and limited to 0, 1, or 2; accepted codon tables are 1 and 2.
- Attribute tables preserve exact duplicates through `duplicate_ordinal`. Attribute type
  description and source display fields retain source nullability.
- Transcript `_rna_edit` rows address raw spliced nucleotide sequence and apply before CDS
  extraction.
- Translation `initial_met`, `_selenocysteine`, `amino_acid_sub`, and `_stop_codon_rt` rows
  address the unedited peptide and apply after translation.
- Applied edits have a non-null replay ordinal. Equal-start edits without receipt-backed
  order are quarantined. `basis_ref_seq` is the unedited coordinate-basis slice;
  `preapply_ref_seq` is replay-derived and may have a different length after overlapping,
  length-changing edits.
- Mature-miRNA spans remain many-valued rows, not transcript flags.

## Sequence states

Each coding model preserves five independently hashed states:

1. transcript `raw_spliced`;
2. transcript `edited_spliced`;
3. translation `translatable_cds`, including start-phase `N` padding;
4. translation `peptide_unedited`, before translation SeqEdits; and
5. translation `peptide_final`.

`model_sequence_state` uses tagged non-null ownership: `(owner_kind, owner_key)`, where
`owner_kind` is `transcript` or `translation`. This avoids mutually exclusive nullable keys
and makes the role/owner/alphabet matrix explicit. Present states reference
`(sequence_sha256, alphabet)` and their byte count must equal both the blob length and the
role receipt.

Export states from fresh pinned-API objects or explicitly invalidate caches. Toggling edit
state on a cached object can return the wrong peptide. A noncoding transcript has no
translation state; it does not have an empty peptide.

## Allele/transcript context

The implemented private primitives are:

- checked 0-based half-open spans;
- immutable borrowed byte and edit views;
- maximal prefix-then-suffix semantic trimming; and
- a separately tagged compatibility differing-region view.

The generic private builder takes a `duckvep_diff_algorithm_t`. The current algorithm tag is
`DUCKVEP_DIFF_VEP_116_BYTE_XOR`; release identity is data, not part of the function name.
Its exact unequal-length Perl byte-XOR behavior and empty insertion sentinel are frozen by
deterministic tests. The semantic edit remains canonical; compatibility runs are derived and
never replace it.

The complete context adds:

- original VCF allele and immutable normalization views;
- piecewise exon/intron/outside projection ordered in transcript orientation;
- explicit validity for transcript coordinates and feature indices;
- endpoint anchor sets that preserve both equidistant intronic candidates and selection
  status;
- raw, edited-transcript, phase-padded CDS, and peptide axes;
- given, used, and transcript-oriented alleles with validation status;
- applied model-edit IDs and conflict status; and
- model, assembly, source, accession/version, and sequence-hash provenance.

For insertions, the semantic reference span is `[p,p)` and projection uses the two genomic
flanks. A partial overlap preserves the unclipped span; VEP-compatible clamping is a derived
view. Missing coordinates use validity flags, never zero sentinels.

Builders own mutable per-worker arenas; published contexts are const borrowed views.
References point into immutable pooled model storage. A context expires when its worker
workspace resets, and public result rows contain no context pointers.

## Validation profiles

`scripts/duckvep_model.R validate` exposes three explicit profiles:

- `structural`: relation shape and local constraints; an empty initialized schema may pass;
- `model_candidate`: nonempty normalized biology, receipts, geometry/cardinality, edit
  replay, five sequence states, and exact sequence hashes; identity authentication remains
  unavailable; and
- `conformance_gate` (default): all model-candidate checks plus canonical identity
  authentication. It fails closed until canonical serialization and golden vectors land.

Diagnostics are deterministic machine-readable JSON. Parquet validation repeats every
constraint explicitly because Parquet views do not enforce DuckDB table constraints.
Missing relations, invalid enum/null coupling, orphaned references, empty required strata,
missing oracle rows, or any declared-slice discord are failures rather than skips.

## Implementation order

1. Land and validate the normalized relation schema and five-state exporter.
2. Complete the lossless context builder and independent scalar projection oracle.
3. Port the VEP-116 consequence rules behind the new context and gate the declared
   biallelic small-variant slice on zero discord.
4. Add DNA HGVS generation with separate `g.` and `c./n.` 3′-shift routines and
   application-equivalence tests.
5. Add protein, compound, mitochondrial, haplotype, and structural HGVS as separately gated
   capabilities.
