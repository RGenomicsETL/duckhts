# DuckVEP → DuckHTS reconciliation & salvage plan

Status: **future proposal / open design.** Records how to fold the stalled sibling
extension `/root/duckvep-c` (commit `9f922c8`) into this tree as `src/duckvep/`
before that repo is archived, what to salvage, why the stall was the *transcript
model import* (not the engine or the tests), and the conformance/format discipline
that differentiates it from `fastVEP`. Companion to `design/roadmap.md` M6. Nothing
here is implemented in DuckHTS yet.

**Superseded in part (2026-07):** an adversarial gpt-5.6-sol (max) review corrected two
things here — the edit model (edits are split nucleotide-pre / peptide-post translation, not
"all before translation") and the VEP pin (RESOLVED to 116). For those, plus the
SO/HGVS-as-sibling-consumers reframe, the lossless `allele_transcript_context`, the two-lane
oracle, conformance-must-gate, and the `duckhgvs-0.4.0` quarantine, **`duckvep_m6a_contract_gate.md`
is authoritative.**

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

## Why it stalled (the correct diagnosis)

Not incrementalism for its own sake, and **not a weak test story**. The stall is the
**Model A transcript-model import**. The engine reads a `tx_flags` bitset
(`HAS_TRANSLATION, BIOTYPE_PROTEIN_CODING, BIOTYPE_NMD, BIOTYPE_MIRNA,
CDS_START_NF/END_NF, SELENOCYSTEINE, STOP_CODON_READTHROUGH, RNA_EDIT,
AMINO_ACID_SUB, MANE_SELECT, MANE_PLUS_CLINICAL, GENCODE_*, CCDS,
READTHROUGH_TRANSCRIPT, UPSTREAM_ATG`) plus a *prepared, edit-applied* `cds_seq` —
never a biotype string or Ensembl object (`/root/duckvep-c/design/duckvep_annotate_input_contract.md`,
`src/include/duckvep_model_reader.h`). That is the *right* cut: Ensembl specifics
live in an import step, the kernel stays uniform.

But `model_a_from_gff()` **does not synthesize `cds_seq` yet** (structural-only), and
the evaluator does **not yet emit** `NMD_transcript_variant` / `mature_miRNA_variant`
(with `mature_mirna(...)` a reserved, unconsumed side relation). So the transcript
classes that depend on MySQL-dump-only data — selenoproteins (`UGA`→Sec, not
`stop_gained`), RNA-edited/readthrough transcripts, mature miRNA, MANE-picked
variants — **cannot be represented in the model**, therefore **cannot be traversed
by the fuzzer**, therefore conformance stays a subset statistic no matter how good
the state machine is. **Model A import is the single pivot.**

## The conformance framework is the crown jewel (the anti-slop engine)

Three tiers, each catching what the tier below cannot (all under
`/root/duckvep-c/conformance/`):

1. **Formal spec** — `so_conformance.R` proves the engine's SO vocabulary is a
   subset of VEP's `%OVERLAP_CONSEQUENCES`, extracted by `extract_so_spec.pl` to
   `data/so_consequences.tsv`, dumped from the C engine via `so_dump.c`.
2. **Rare-class witness fuzzer** — `generate_witnesses.R` (D2) *deliberately tiles
   the equivalence classes where consequence bugs live*: each donor / acceptor /
   5th-base / region boundary, start & stop codons, exon edges, same-codon and
   two-codon MNVs. `fuzz_witnesses.R` (D3) runs VEP `--gff` (sole oracle) vs
   `duckvep_annotate` on the SAME Model A + optional prepared CDS, exact
   SO-term-**set** comparison. **This already found a real VEP-state-machine bug**
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
reinvent formats:** Model A is a **DuckDB relation / Parquet**; supplementary
annotations (ClinVar/gnomAD/dbSNP/COSMIC) are **joins** to Parquet / DuckLake; the
Ensembl MySQL dump is read by DuckDB's `mysql` scanner. Storage is solved by the
engine; the consequence logic is the only risk we take on. This aligns with
`design/north_star.md` ("no private cache format").

## Salvage inventory (before archiving `/root/duckvep-c`)

Lift into `src/duckvep/`. Priority = uniqueness (can't be rebuilt cheaply).

**Tier A — irreplaceable:**
- `conformance/` in full — the three tiers above, `generate_effect_rules.pl`,
  `extract_so_spec.pl`, `generate_so_metadata.pl`, `so_dump.c`,
  `model_a_from_gff.sql`, `model_a_with_cds.R`, `data/`.
- `kernel/` — the pure, DuckDB-free consequence engine:
  `duckvep_classify.c`, `duckvep_effect.c`, `duckvep_sweep.c`, `duckvep_sv.c`,
  `duckvep_delta.c`, `duckvep_codon.c`, `duckvep_coding.c`,
  `duckvep_haplotype.c` (the haplosaurus lane, slice 5 open), `duckvep_projection.c`,
  `duckvep_kernel.c`, `duckvep_so.c` + `include/{duckvep_kernel,duckvep_so}.h`.
- The design docs: `duckvep_effect_ctx_architecture.md`,
  `vep_consequence_state_machine.md`, `duckvep_annotate_input_contract.md`,
  `duckvep_span_sv_cnv_consolidation_2026-06-23.md`,
  `duckvep_current_design_state_machines.md`,
  `duckvep_bcftools_csq_port_plan_2026-06-09.md`, `duckvep_layer_keys.md`.

**Tier B — engine glue to lift + rename:**
- `src/consequence_udf.c`, `src/annotate_table.c`, `src/duckvep_variant_tile.c`,
  `src/duckvep_model_reader.c`, `src/duckvep_relation_cursor.c`,
  `src/duckvep_chunk_reader.c`, `src/chrom_interner.c` + their `include/` headers.

**Tier C — DEDUPLICATE against DuckHTS (do not copy; use the DuckHTS copy):**
- `bcf_reader.c`, `cgranges_api.c`, `variantkey_udf.c`, `vep_parser.c`,
  `bgzip.c`, `tabix_reader.c`, `seq_reader.c`, `quality_encoding.c`,
  `bcftools_norm_udf.c`, `bcftools_shim.c`, `wasm_http_hfile.c`,
  `include/seq_encoding.h`, `include/wasm_socket_compat.h`. These are near-duplicates
  of files already in DuckHTS; unify to one copy.
- **Merge the SIMD dispatchers:** `src/simd/duckvep_simd_*` → the DuckHTS
  `src/simd/duckhts_simd_*` framework. One extension, one dispatch table.

## Fuse mechanics

1. Create `src/duckvep/` (self-contained) + `src/duckvep/kernel/`.
2. Deduplicate Tier C; delete duckvep's copies; repoint includes to DuckHTS.
3. Merge SIMD dispatch (Tier C).
4. Wire the build: `CMakeLists.txt`, `r/Rduckhts/R/bootstrap.R`,
   `r/Rduckhts/configure`, `r/Rduckhts/configure.win` (per `AGENTS.md` — a new
   public function is incomplete until wired on Unix + Windows + R package).
5. Public surface: one `duckhts_csq(...)` / annotate table function reusing DuckHTS
   `read_gtf` / `fasta_nuc` / cgranges / VariantKey.
6. Land the **Model A import** (the pivot) as a first-class DuckHTS surface —
   schema pinned from `github.com/Ensembl/ensembl` (see "Model A data provenance"
   below), assembly-parameterized (GRCh38 + GRCh37), species-agnostic. Ensembl
   MySQL-dump/scanner → DuckDB SQL (`transcript_attrib ⋈ attrib_type`,
   `translation_attrib ⋈ translation`) → flat Model A relation + `tx_flags`
   (MANE flags may come from GENCODE GTF tags instead). **CORRECTION (2026-07, see
   `duckvep_m6a_contract_gate.md`): edits are NOT all "before translation."** Ensembl
   splits them by scope — nucleotide-scope `_rna_edit` changes the *spliced nucleotide*
   sequence pre-translation, but peptide-scope `_selenocysteine` (UGA→Sec, not
   `stop_gained`), `amino_acid_sub`, `initial_met`, `_stop_codon_rt` modify the *peptide
   after translation*. Model A v2 carries **typed edit relations** `(scope, tx/translation
   id, code, start, end, alt_seq)`, not flat flags, and stores/derives both the edited
   spliced transcript sequence and the final edited reference peptide; emit side relations
   (`mature_mirna`, regulatory). Only then does the fuzzer reach the rare clinical classes.

## Model A data provenance, schema source, species, and assemblies

The Model A import is a data-engineering problem; get its *inputs* right.

- **Schema from the Ensembl API repo, not a live MySQL instance.** The authoritative
  DDL is versioned in `github.com/Ensembl/ensembl` (`sql/table.sql`, core: transcript
  / translation / exon / exon_transcript / *_attrib / attrib_type / seq_region) and
  `github.com/Ensembl/ensembl-variation` (`sql/table.sql`, variation; local copy at
  `/root/ensembl-variation/sql/table.sql`). Pin the schema to a release tag and build
  the Model A SQL against *that* file — reverse-engineering a live DB is fragile and
  version-ambiguous. The DuckDB `mysql` scanner (or a server-free dump) supplies the
  *rows*; the repo supplies the *shape*.
- **All Ensembl species, by construction.** The kernel already reads `tx_flags`, never
  a species biotype string (`duckvep_annotate_input_contract.md`), so a species is
  just another Model A relation. Keep it that way: no per-species branches in the
  engine, only per-species import runs.
- **GRCh38 AND GRCh37 — GRCh37 is clinically load-bearing.** Clinical genomics still
  runs on GRCh37/hg19. Ensembl serves GRCh37 from a separate archive DB
  (grch37.ensembl.org); VEP has `--grch37`. Model A import must be assembly-parameterized
  and both assemblies are first-class, not GRCh38-only. Reference FASTA for GRCh37 is
  `human_g1k_v37.fasta` (see `/root/GRCh37/`).
- **MANE tx_flags from GENCODE, not only Ensembl attribs.** GENCODE ships MANE-tagged
  GTF (`tag "MANE_Select"` / `"MANE_Plus_Clinical"`); `v46lift37` carries MANE onto
  GRCh37 (MANE ENST/NM share identical CDS by construction, so it is robust across an
  Ensembl-release gap). This is a simpler, version-stable source for `MANE_SELECT` /
  `MANE_PLUS_CLINICAL` than the core `transcript_attrib` join, and it works on GRCh37.
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

- **(1) Pinned conformance manifest** — *which* oracle: VEP release + cache-vs-DB,
  GENCODE-vs-RefSeq, assembly, bcftools/haplosaurus versions. Separates version
  skew from real divergence.
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
live in `duckvep_m6a_contract_gate.md`.

## duckhgvs-0.4.0 — quarantine (not the head start it looks like)

`/root/duckhgvs-0.4.0` (v0.4.0) is a design spike, NOT a validated engine, despite the
handoff (`/root/duckvep-c/docs/HANDOFF.md:33`) calling it "already exists." It hard-codes
`DUCKHGVS_SEQ_MAX 4096` / `DUCKHGVS_TRANSCRIPT_EXON_MAX 512`, only projects alleles wholly
inside one exon, and its README's "not yet" list covers accession resolution, intronic
`c./n.`, protein normalization, SV, and full grammar; the 4096 buffer already excludes 9,742
real transcripts. Salvage ideas + test vectors; **do not fold its ABI.** Its
render-HGVS-from-`duckvep_consequence_t` architecture is exactly the design error M6a
corrects (see `duckvep_m6a_contract_gate.md`).
