# DuckVEP → DuckHTS reconciliation & salvage plan

Status: **future proposal / open design.** Records how to fold the stalled sibling
extension `/root/duckvep-c` (commit `9f922c8`) into this tree as `src/duckvep/`
before that repo is archived, what to salvage, why the principal functional stall was the
*transcript model import*, why the old ABI/glue still require M6a adaptation, and the
conformance/format discipline
that differentiates it from `fastVEP`. Companion to `design/roadmap.md` M6. Nothing
here is implemented in DuckHTS yet.

**Superseded in part (2026-07):** an adversarial gpt-5.6-sol (max) review corrected two
things here — the edit model (edits are split nucleotide-pre / peptide-post translation, not
"all before translation") and the VEP pin (RESOLVED to 116). For those, plus the
SO/HGVS-as-sibling-consumers reframe, the lossless `allele_transcript_context`, the two-lane
oracle, conformance-must-gate, and the `duckhgvs-0.4.0` quarantine,
**`duckvep_m6a_contract_gate.md` plus `duckvep_model_a_v2.md` are authoritative.**

## What duckvep-c actually is (do not under-credit it)

Not a "bcftools csq engine" — a **table-driven port of Ensembl VEP 116's actual
consequence evaluator**, source-grounded in `Bio/EnsEMBL/Variation`
(`%OVERLAP_CONSEQUENCES` in `Constants.pm`, `VariationEffect.pm`,
`BaseTranscriptVariationAllele.pm`) and the haplosaurus path
(`TranscriptHaplotypeContainer::_mutate_sequences`, `CDSHaplotype.pm`,
`ProteinHaplotype.pm`). Shape: `effect_ctx → sequence_delta → compiled rule table`
(tier / rank / include / predicate), the same machinery serving vanilla, haplotype,
and SV cases. See `/root/duckvep-c/design/duckvep_effect_ctx_architecture.md` and
`vep_consequence_state_machine.md`.

## Why it stalled (principal cause plus the ABI gate)

Not incrementalism for its own sake, and **not a weak test story**. The stall is the
**Model A transcript-model import**. The engine reads a `tx_flags` bitset
(`HAS_TRANSLATION, BIOTYPE_PROTEIN_CODING, BIOTYPE_NMD, BIOTYPE_MIRNA,
CDS_START_NF/END_NF, SELENOCYSTEINE, STOP_CODON_READTHROUGH, RNA_EDIT,
AMINO_ACID_SUB, MANE_SELECT, MANE_PLUS_CLINICAL, GENCODE_*, CCDS,
READTHROUGH_TRANSCRIPT, UPSTREAM_ATG`) plus a *prepared, edit-applied* `cds_seq` —
never a biotype string or Ensembl object (`/root/duckvep-c/design/duckvep_annotate_input_contract.md`,
`src/include/duckvep_model_reader.h`). Keeping Ensembl-specific ingestion outside the
kernel is the right cut; the old `cds_seq + tx_flags` data shape is not.

But `model_a_from_gff()` **does not synthesize `cds_seq` yet** (structural-only), and
the evaluator does **not yet emit** `NMD_transcript_variant` / `mature_miRNA_variant`
(with `mature_mirna(...)` a reserved, unconsumed side relation). So the transcript
classes that depend on MySQL-dump-only data — selenoproteins (`UGA`→Sec, not
`stop_gained`), RNA-edited/readthrough transcripts, mature miRNA, MANE-picked
variants — **cannot be represented in the model**, therefore **cannot be traversed
by the fuzzer**, therefore conformance stays a subset statistic no matter how good
the state machine is. Model A v2 is the data pivot. The second M6a gate is the lossless
per-allele/transcript context: the old compact-result ABI cannot be the SO/HGVS boundary.

## The conformance framework is the crown jewel (the anti-slop engine)

Three tiers, each catching what the tier below cannot (all under
`/root/duckvep-c/conformance/`):

1. **Formal spec** — `so_conformance.R` proves the engine's SO vocabulary is a
   subset of VEP's `%OVERLAP_CONSEQUENCES`, extracted by `extract_so_spec.pl` to
   `data/so_consequences.tsv`, dumped from the C engine via `so_dump.c`.
2. **Rare-class witness fuzzer** — `generate_witnesses.R` (D2) *deliberately tiles
   the equivalence classes where consequence bugs live*: each donor / acceptor /
   5th-base / region boundary, start & stop codons, exon edges, same-codon and
   two-codon MNVs. `fuzz_witnesses.R` (D3) runs VEP `--gff` (the geometry lane) vs
   `duckvep_annotate` on the SAME Model A + optional prepared CDS, exact
   SO-term-**set** comparison. Curated edit states additionally require the pinned cache/DB
   lane. **The geometry fuzzer already found a real VEP-state-machine bug**
   (`intron_variant` at essential splice sites → the `PRE_WITHIN_INTRON` fix).
3. **Real-data statistical audit** — `stratified_conformance.R` (D4) samples real
   VCF (GIAB default, ClinVar via `--corpus`), runs VEP vs duckvep, reports SO-set
   agreement per stratum; `statistical_conformance.R` attaches Clopper-Pearson 95%
   bounds. Has a built-in `--fastvep` comparator.

Rule table is **generated**, not hand-written: `generate_effect_rules.pl` compiles
VEP `Constants.pm` → the engine's rule program.

### Why this beats fastVEP's accuracy claim

fastVEP (`/root/fastVEP`, github.com/Huang-lab/fastVEP, `785922e`) claims
**"100% concordance across 23 fields"** — but on **173 real chr22 example variants**
(2,340 shared allele-transcript pairs), **shared pairs only** (emission misses/extras
dropped; `validation/compare_vep.py`, `manuscript/manuscript.md §3.3`). 173 common
variants never hit a selenocysteine `UGA`, a splice-donor 5th base, an incomplete
terminal codon, a start/stop-retained, or a same-codon MNV — the states **rare in
real BCF and decisive in a clinical report**. Tier 2 exists precisely to falsify
that class of claim. duckvep's differentiator is *the fuzzer + honest real-data CI
bounds*, not raw speed.

## Format philosophy — the architectural win (make it non-negotiable)

fastVEP reinvents formats: `fastSA` is an echtvar-inspired bespoke binary (Var32 /
delta / chunked-ZIP; `crates/fastvep-sa/src/lib.rs` "echtvar-inspired"), and the
cache is a serde port of VEP's Sereal structure. **DuckVEP-in-DuckHTS must not
reinvent formats:** Model A is a **DuckDB relation pack / Parquet** with a derived fused
execution view; supplementary
annotations (ClinVar/gnomAD/dbSNP/COSMIC) are **joins** to Parquet / DuckLake; the
Ensembl MySQL dump is read by DuckDB's `mysql` scanner. Storage is solved by the
engine; the consequence logic is the only risk we take on. This aligns with
`design/north_star.md` ("no private cache format").

## Salvage inventory (before archiving `/root/duckvep-c`)

Salvage into `src/duckvep/` only after applying the dispositions below. Priority = unique
algorithms, generators, source-grounded tests, and provenance—not old ABI surface area.

**Verified file-by-file (2026-07, Sonnet 5 audit against the trees): see
`duckvep_salvage_inventory_audit.md` for the corrected, evidence-backed inventory table
(existence, line counts, per-file dedupe diffs, and gaps).** Key corrections folded below.

**Tier A — irreplaceable algorithms, generators, and tracked seed fixtures:**
- From commit `9f922c8`, salvage the 28 tracked `conformance/` programs, binding tables, and
  seed fixtures. Do **not** copy the dirty directory "in full": the 39-file filesystem count
  includes modified, untracked, and ignored result artifacts. Freeze generated oracle data
  only when an authoritative manifest names its checksum.
- From `kernel/`, salvage the pure consequence algorithms and generated-rule machinery:
  `duckvep_classify.c`, `duckvep_effect.c`, `duckvep_sweep.c`, `duckvep_sv.c`,
  `duckvep_delta.c`, `duckvep_codon.c`, `duckvep_coding.c`,
  `duckvep_haplotype.c` (the haplosaurus lane, slice 5 open), `duckvep_projection.c`,
  `duckvep_kernel.c`, and `duckvep_so.c`. The old `duckvep_kernel.h` is retained as rewrite
  evidence, **not lifted as the v2 ABI**: its compact scalar result and HGVS handoff are the
  architecture M6a replaces.
- The design docs: `duckvep_effect_ctx_architecture.md`,
  `vep_consequence_state_machine.md`, `duckvep_annotate_input_contract.md`,
  `duckvep_span_sv_cnv_consolidation_2026-06-23.md`,
  `duckvep_current_design_state_machines.md` (+ its `.qmd` source).
  **Already lifted (do not re-copy — reconcile if drifted):**
  `duckvep_bcftools_csq_port_plan_2026-06-09.md`, `duckvep_layer_keys.md` are already
  present under `/root/duckhts/design/`.
- **Also required for Tier A to compile (audit gap):** the 11 private kernel headers in
  `kernel/src/*.h`, the 2 generated `kernel/src/*.inc` (`duckvep_effect_rules.inc`,
  `duckvep_so_metadata.inc`), and `kernel/CMakeLists.txt`.

**Tier B — engine glue, split by disposition:**
- Reuse/adapt the narrow transport helpers where their contracts survive:
  `src/consequence_udf.c`, `src/duckvep_variant_tile.c`,
  `src/duckvep_relation_cursor.c`, `src/duckvep_chunk_reader.c`, `src/chrom_interner.c` and
  matching headers.
- Rewrite `src/duckvep_model_reader.c`/header for the normalized Model A v2 relation pack,
  five sequence states, typed edit slices, and immutable reference bundles. Adapt or rewrite
  `src/annotate_table.c` around the new context/result lifetime. Neither is a prefix rename.
- **Note (audit):** `src/duckvep.c` (the `dv_`-prefixed extension entry point) is NOT
  lifted verbatim — it is superseded by the new `duckhts_csq` registration (fuse step 5).

**Tier C — DEDUPLICATE against DuckHTS (do not copy; use the DuckHTS copy):**
- `cgranges_api.c`, `variantkey_udf.c`, `vep_parser.c`,
  `bgzip.c`, `tabix_reader.c`, `seq_reader.c`, `quality_encoding.c`,
  `bcftools_norm_udf.c`, `bcftools_shim.c`, `wasm_http_hfile.c`,
  **`hts_index_builder.c`, `hts_meta_reader.c`, `interval_udf.c`** (audit-added; also
  near-dups), `include/seq_encoding.h`, `include/wasm_socket_compat.h` + the 9 paired
  headers the audit enumerates. These are near-duplicates (prefix-only diffs) of files
  already in DuckHTS; unify to one copy.
- **`bcf_reader.c` is NOT a same-generation near-dup — keep DuckHTS's, discard
  duckvep-c's.** DuckHTS's copy (5028 vs 4551 lines) carries the `decode_error_policy`
  (`null`/`warn`/`error`) + GT-preflight feature that duckvep-c's fork predates.
- **SIMD: discard `src/simd/duckvep_simd_*` entirely — do NOT merge.** DuckHTS's
  `src/simd/duckhts_simd_*` is strictly ahead (it added `BAM_NT16_COUNTS` +
  `NT16_GC_COUNTS`); duckvep-c's SIMD tree contributes zero unique content.

## Fuse mechanics

1. Freeze the `9f922c8` tracked-source inventory and selected, checksummed conformance seeds;
   do not move dirty run outputs.
2. Land the **Model A v2 import first**: the normalized relation pack, separate logical-model
   and build identities, artifact/selection manifests, five sequence-state exporter/hashes,
   validation, and derived execution view specified by `duckvep_model_a_v2.md`. The schema
   is species-neutral; first supported packs are human GRCh38 and GRCh37.
3. Create `src/duckvep/` and port the unique classifier/delta/sweep/rule algorithms behind
   the new context and result boundary. Rewrite the old model reader; adapt annotation glue.
4. Deduplicate Tier C against DuckHTS (keep DuckHTS's newer `bcf_reader.c`) and discard
   duckvep-c's `src/simd/` entirely.
5. Wire CMake plus R package Unix/Windows bootstrap/configure paths and both SQL/R test
   levels before exposing a public surface.
6. Add the narrow capability-manifest-bound annotation surface only after the five-state
   hash gate and zero-discord frozen fixtures pass. Reuse DuckHTS readers, cgranges,
   VariantKey, and sequence providers; do not create a private cache format.

## Model A data provenance, schema source, species, and assemblies

The Model A import is a data-engineering problem; get its *inputs* right.

- **Schema from the pinned Ensembl API source, not a live-instance guess.** Model A's core
  DDL is [`sql/table.sql`](https://github.com/Ensembl/ensembl/blob/c0cf13daa961d80584bad797b2eb0ff3a7500ef3/sql/table.sql)
  at commit `c0cf13daa961d80584bad797b2eb0ff3a7500ef3`. Variation DDL is not needed to
  define the transcript pack; pin it separately when importing variation/supplementary
  data. The scanner or server-free dump supplies rows; the repository supplies shape.
- **Species-neutral shape, human-first claim.** Keep species-specific vocabulary in the
  importer and named attributes, not engine branches. Do not claim all-species parity until
  codon-table, attribute, and source-selection inventories pass; M6b supports human
  GRCh37/38 first.
- **GRCh38 AND GRCh37 — GRCh37 is clinically load-bearing.** Clinical genomics still
  runs on GRCh37/hg19. Ensembl serves GRCh37 from a separate archive DB
  (grch37.ensembl.org); VEP has `--grch37`. Model A import must be assembly-parameterized
  and both assemblies are first-class, not GRCh38-only. Reference FASTA for GRCh37 is
  `human_g1k_v37.fasta` (see `/root/GRCh37/`).
- **MANE values from GENCODE, not only Ensembl attribs.** GENCODE ships MANE-tagged
  GTF (`tag "MANE_Select"` / `"MANE_Plus_Clinical"`); `v46lift37` carries MANE onto
  GRCh37. Preserve the accession value and source receipt; then hash-verify ENST/NM sequence
  identity before promoting the crosswalk. The runtime MANE bits are derived selection
  metadata, not identity proof.
  Reference: `/root/bioconnect-sprint-py/docs/data_versions.md` (GENCODE v46lift37,
  `--mane_select --pick`).

## The clinical layer above consequence: a configurable SQL ACMG engine

The endgame is not consequence annotation alone — it is a clinical stack:
**consequence (duckvep) → supplementary-annotation joins (gnomAD / REVEL / SpliceAI /
ClinVar as Parquet/DuckLake) → ACMG/AMP classification.** The classification layer is
*already pure SQL* and composes natively on DuckDB — a direct validation of the
"everything composes in SQL" thesis. Reference implementation:
`/root/bioconnect-acmg/manifests/{00_spec,01_acmg_rules,02_combine}.sql`.

Design principles to adopt (from that engine):
- **Rules as one SQL kernel** over an `annotations_spec` relation (annotations +
  per-gene resolved thresholds) emitting `applied_criteria(variant_key, criterion,
  direction, points, evidence)`; combine = signed-point sum (Tavtigian 2020:
  Supporting/Moderate/Strong/Very-Strong = ±1/2/4/8) → band → gates.
- **Configurability is DATA, not code.** Per-gene / VCEP threshold exceptions
  (`ba1_maf`/`bs1_maf`/`pm2_maf`), known-pathogenic common-variant exception lists,
  and gene-disease validity caps are *rows*, so a panel/VCEP override never touches
  the engine. This is the "configurable SQL ACMG engine."
- **Discipline the engine encodes** (keep these): NULL frequency → PM2 *abstains*
  ("no data" ≠ "rare"); a *single* calibrated predictor per axis (REVEL for missense
  PP3/BP4, SpliceAI for splicing) with ClinGen-calibrated thresholds — never a
  consensus stack (documented bias); CNV/SV are a *separate* kernel (Riggs 2020), not
  this one; gene-disease validity caps the strength.

This becomes roadmap **M8** (clinical SQL layer), sitting on M6 (consequence) +
supplementary-annotation Parquet joins. It is SQL-native, so most of it is manifests,
not C — but the supplementary-annotation ingestion (dbNSFP/gnomAD/ClinVar/SpliceAI →
Parquet, filtering-AF max-pop computation) is real prep work, and calibrated-threshold
provenance (Pejaver 2022, Walker 2023) belongs in the M1 manifest, versioned.

## The ERRATA + manifest discipline (why a flat ERRATA.md is not enough)

A VEP port needs riker's ERRATA instinct but a richer structure, because SO-term
divergences are categorical/clinical, the oracle is versioned and forks
(release, **cache vs DB**), and three oracles (VEP / haplosaurus / bcftools csq)
disagree with each other.

- **(1) Pinned conformance manifest** — *which* model build, engine/rule binary, oracle
  container/modules/command, cache-vs-DB lane, annotation source, assembly, numeric allele
  bounds, and required rare-state witness counts. This separates version skew and empty
  strata from real divergence; emitted-pair equality is checked before SO-set equality.
- **(2) `ERRATA.md`** — three reason-classes, each entry naming *which oracle* it
  diverges from and *backed by a regression witness*: (a) intentional design
  divergence (cite the Ensembl/VEP source mechanism, riker-style); (b) VEP/haplosaurus
  bugs we deliberately do not replicate (the direct analogue of riker's chimeric
  overlap-clip entry — VEP known bugs: multi-allelic HGVS, exon/intron boundaries);
  (c) not-yet-implemented scope = the honest coverage boundary.
- **(3) Every errata entry ⟺ a conformance witness.** No prose-only divergences.

**Trap:** until Model A import exists, selenocysteine/edit/miRNA/MANE mismatches are
*provenance bugs*, not intentional divergence. Keep them as open coverage items;
they must not enter `ERRATA.md` as class-(a). The manifest + import must precede a
trustworthy errata. The accepted-mismatch ledger (roadmap M1) **is** the class-(c)
coverage boundary — shrinking it is the breakthrough metric.

## Version-pin decision — RESOLVED to VEP 116 (2026-07)

Verified against the working trees: `Bio/EnsEMBL/Variation/Utils/Constants.pm` is
**byte-identical** between 115 (`/root/ensembl-variation` @ `release/115`) and 116
(`ensembl-vep-116.0-0`), so regenerating the rule table (`generate_effect_rules.pl`)
would produce the *same* SO vocabulary/rank/tier either way — the pin is **not** a
rule-table decision. The real 115→116 difference is **155 changed lines in
`VariationEffect.pm`**, concentrated on `stop_lost`, `partial_codon`, `stop_retained`,
`_overlaps_stop_codon`→`_overlaps_stop_codon_cil`, `_ins_del_stop_altered`→`…_cil`, and
final-stop-length logic — the stop/inframe/partial-codon predicates at our coding
frontier. **Decision: pin VEP 116** (the salvaged kernel is 116-shaped). If GRCh37/clinical
needs 115, it is a *named* compatibility target with consciously back-ported predicates,
never a 115 wrapper on 116 assumptions. Record loaded-module hashes, not `vep --version`.
Full detail + the sibling-consumer / Model-A-v2 / oracle-lane / conformance-gating decisions
live in the two M6a contract notes.

## duckhgvs-0.4.0 — quarantine (not the head start it looks like)

`/root/duckhgvs-0.4.0` (v0.4.0) is a design spike, NOT a validated engine, despite the
handoff (`/root/duckvep-c/docs/HANDOFF.md:33`) calling it "already exists." It hard-codes
`DUCKHGVS_SEQ_MAX 4096` / `DUCKHGVS_TRANSCRIPT_EXON_MAX 512`, only projects alleles wholly
inside one exon, and its README's "not yet" list covers accession resolution, intronic
`c./n.`, protein normalization, SV, and full grammar. On its bundled GRCh37 RefSeq model,
28,901/108,760 transcripts (26.6%) exceed 4096 nt spliced length; the older 9,742 figure was
a CDS count from a different VEP-116/GRCh38 dataset. Salvage the apply-then-diff test pattern
and seed vectors; **do not fold its ABI.** Its
render-HGVS-from-`duckvep_consequence_t` architecture is exactly the design error M6a
corrects (see `duckvep_m6a_contract_gate.md`).
