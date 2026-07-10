# DuckVEP M6a: the contract gate before folding duckvep-c

Status: **open design / implementation contract.** M6a is the "contract gate" that must
land *before* any `/root/duckvep-c` code is folded into `src/duckvep/`. It records the
reframe and the verified findings of an adversarial gpt-5.6-sol (max) review (2026-07-10)
of the duckvep + duckhgvs plan, cross-checked against actual source in `/root/duckvep-c`,
`/root/duckhgvs-0.4.0`, and Ensembl VEP 115/116. Companion to `roadmap.md` (M6a–M6d) and
`duckvep_reconciliation.md` (salvage plan). **Where this note disagrees with those two on
the edit model or the VEP pin, THIS note wins** — it corrects them.

## The reframe: SO and HGVS are sibling consumers, not a pipeline

The stalled handoff (`/root/duckvep-c/docs/HANDOFF.md:62-63`) wires HGVS to render *from*
the compact `duckvep_consequence_t` "only on request." That is the core design error: the
consequence row carries scalar cDNA/CDS/protein positions and one ref/alt amino acid (and
several delta paths deliberately set positions to `-1`). HGVS needs *ranges* + full
sequence context: general indels, exon-boundary crossings, intronic coordinates, peptide
delins, frameshift/extension termini, duplication recognition, uncertainty.

**Decision:** introduce an internal immutable **`allele_transcript_context`** that both SO
evaluation and HGVS rendering consume as *siblings*. It carries: original VCF allele;
trimmed semantic edit; raw + projected genomic/transcript ranges; transcript-oriented
alleles; cDNA/CDS/peptide ranges; full sequence-provider identity; edit provenance;
normalization offsets. HGVS is **not** derived from `duckvep_consequence_t`.

## Model A v2: represent the right biological object

`cds_seq + tx_flags` cannot encode Ensembl's edits, and `duckvep_reconciliation.md` was
wrong to say "apply … BEFORE translation." Ensembl splits edits by scope:
- **Nucleotide-scope, pre-translation:** `_rna_edit` changes the spliced nucleotide sequence.
- **Peptide-scope, post-translation:** `_selenocysteine` (UGA→Sec, *not* `stop_gained`),
  `amino_acid_sub`, `initial_met`, `_stop_codon_rt` modify the peptide *after* translation.
A boolean flag like `SELENOCYSTEINE` carries neither position nor replacement, so it cannot
drive the peptide edit; rare-class SO terms and every HGVSp would use the wrong reference
protein.

**Decision:** **typed edit relations** `(scope, transcript_or_translation_id, code, start,
end, alt_seq)`. Keep `tx_flags` for boolean predicates + selection metadata only, never edit
payloads. Store or reproducibly derive BOTH the edited spliced transcript sequence AND the
final edited reference peptide. Model A v2 also carries identity: `model_id`, provider,
release, assembly, transcript stable-ID/version, translation ID/version, RefSeq aliases,
reference accession/version (e.g. `NC_…`), and sequence hashes. GENCODE MANE tags are
*selection annotations*, not proof of sequence identity — crosswalk and hash-verify.

**M6a exit evidence:** for release 116 export each transcript's raw spliced, edited spliced,
translatable CDS, unedited translation, and final edited peptide from the Ensembl API;
compare hashes to the importer, stratified by edit code / phase / partiality / codon table.

## VEP pin: RESOLVED to 116

The roadmap left this open (115 for ClinVar CLNHGVS vs 116 kernel grounding). Verified
2026-07-10:
- `Bio/EnsEMBL/Variation/Utils/Constants.pm` is **byte-identical** between 115
  (`/root/ensembl-variation` @ `release/115`) and 116 (`ensembl-vep-116.0-0`). Regenerating
  the rule table from `Constants.pm` (`generate_effect_rules.pl`) therefore does NOT change
  the SO vocabulary / rank / tier — the pin is *not* a rule-table decision.
- `VariationEffect.pm` differs by **155 lines**, concentrated on `stop_lost`,
  `partial_codon`, `stop_retained`, `_overlaps_stop_codon` → `_overlaps_stop_codon_cil`,
  `_ins_del_stop_altered` → `_ins_del_stop_altered_cil`, and final-stop-length logic — i.e.
  exactly the stop / inframe / partial-codon predicates at our coding frontier.

**Full per-hunk catalog: `vep_115_to_116_predicate_delta.md`** (17 meaningful hunks, 13 on the
coding frontier). The "155 lines" headline is directionally right but *undercounts* — a
116-targeted engine must also implement: `inframe_insertion`/`inframe_deletion` guards that
exclude already-stop reference codons; `stop_gained` now deferring to `stop_lost` (mutual
exclusion by call order); a `frameshift` post-stop guard; an `X` (incomplete/ambiguous peptide)
routing change shared by `stop_lost`/`stop_retained`; plus one real behavioral hunk each in
`Sequence.pm` (the `get_3prime_seq_offset` 3′-shift cap — 115 under-shifts indels before long
repeats) and `TranscriptVariationAllele.pm` (`_var2transcript_slice_coords` now *clamps*
boundary-crossing variants instead of returning `undef`, so 116 emits `c./n.` where 115 emitted
none, for **all** transcripts). Two hunks (the dead-looking `$consider_ins_len` param; whether
the 3′-shift cap change touches `c./n.` or only `g.`) are escalated to gpt-5.6-sol.

**Decision:** **pin VEP 116.** The salvaged kernel, conformance data, and controlled GFF are
116-shaped. If GRCh37/clinical requires 115, make it an explicitly *named* compatibility
target with consciously back-ported predicate behavior — never a 115 wrapper on 116
assumptions. Record loaded-module hashes, not `vep --version`.

## Two oracle lanes

Once Model A carries MySQL-derived translation edits, VEP `--gff` no longer sees the same
transcript model (the GFF source ignores translation-edit payloads). So "Model A import
unblocks rare-class fuzzing" unblocks the *engine*, not its *oracle*.
- **Lane 1 — `--gff`:** synthetic geometry / topology coverage (the witness fuzzer).
- **Lane 2 — pinned `--cache`/DB mode:** curated transcript + translation facts
  (selenocysteine / readthrough / AA-sub). Rare-class parity claims **require Lane 2**, with
  Model A exported from that exact release/source.

## Conformance must GATE, not just report

The current harness reports discord but does not fail on it and skips clean when VEP is
absent (`fuzz_witnesses.R`: `quit(status=0)` skip path; discord surfaced as counts, not a
non-zero exit). `statistical_conformance.R` applies Clopper–Pearson to non-independent
`(variant, transcript)` pairs, and `conformance_progress.R` folds every dated + root-level
stratified file rather than one manifest-bound run.

**Decisions:**
- The declared supported slice **fails its release gate on any discord or skipped oracle.**
- Live VEP lives in a pinned nightly/release job (network-disabled container, forked,
  module + input hashes recorded), **never CRAN**; ordinary CI uses frozen, checksummed
  oracle fixtures.
- Primary uncertainty unit = independently sampled variant/locus (discordant if *any*
  supported transcript pair disagrees); report pair rates descriptively with cluster-aware
  (design-effect / bootstrap) intervals.
- One manifest selects one authoritative run.

## duckhgvs-0.4.0: quarantine as a spike

`/root/duckhgvs-0.4.0` (v0.4.0) is a design spike, **not** a validated engine, despite the
handoff calling it "already exists." Full verified audit: `duckhgvs_0_4_0_audit.md`. Key facts:
`include/duckhgvs.h` hard-codes `DUCKHGVS_SEQ_MAX 4096`; genomic→transcript projection is
verified to accept only alleles wholly inside one exon (`src/duckhgvs.c:1734-1742` — both
intronic and exon-boundary cases → `ENOSUP`); the README "not yet" list = no accession
resolution, no intronic `c./n.`, no protein normalization, no SV, no full grammar.

- **Buffer break (spliced length is the honest metric — apply/`n.`/equivalence all route
  through the 4096 buffer):** ~26.6% of transcripts exceed 4096 nt spliced length on the
  bundled GRCh37 RefSeq GFF (28,901/108,760; max 109,224); 5.8% exceed it on CDS length
  (4,316/74,880; max CDS 107,976 — the figure that corroborates the data source). The "9,742"
  figure came from codex's VEP-116 (GRCh38) CDS scan — a *different dataset*; both agree the
  exclusion is double-digit-percent, not an edge case. `DUCKHGVS_TRANSCRIPT_EXON_MAX 512` does
  NOT break on GRCh37 (max 363, *TTN*) — a latent risk, not a demonstrated one.
- **Its "conformance" is mostly not real:** only the ClinVar `CLNHGVS` lane is an external-
  oracle comparison, and it is restricted to pure SNVs (the trivial `g.` case). The VEP-CSQ
  "oracle" is a parse→format string round-trip on VEP's own output text (never touches
  REF/ALT/transcript); `vcfcheck` has no expected value at all; the TSV corpora are hand-typed
  expectations on 9–14 nt toy sequences. "Tests pass" ≠ HGVS correctness.
- **Salvage (M6c):** lift the **apply-then-diff** equivalence *pattern* + its ~20 `APPLY`/`EQUIV`
  seed vectors (`duckhgvs_apply_to_sequence` / `_equivalent_on_sequence`) — exactly the
  application-equivalence test shape M6c needs. **Discard** `duckhgvs_projected_allele_t` as a
  data model (single scalar tx position, no ranges/provenance, hard-fails on intronic/boundary
  cases): it is a *negative example* for the lossless `allele_transcript_context`. Do NOT fold
  the ABI. The bundled 48 MB RefSeq GFF is orphaned (unreferenced, no NOTICE, dir not even a git
  repo) — discard; re-vendor deliberately with checksum + NOTICE if M6c needs it.

## Scope the first claim: human GRCh37/38 only

The kernel's codon-table enum defines only tables 1 (standard) + 2 (vert-mito)
(`kernel/src/duckvep_codon.h`). The Model A *shape* is species-agnostic, but the *engine* is
not table-complete, so "all species by construction" is premature. First public claim =
human GRCh37/38; promote species only after codon-table / attribute inventories pass.

## HGVS ground-truth: adjudicated panel, not a single oracle

VEP / VariantValidator / Mutalyzer / biocommons-hgvs / ClinVar disagree (algorithmic
choices, nomenclature revisions, reference providers, transcript versions). **Decision:**
pin a dated **HGVS recommendations snapshot** as normative; VEP-116 as the operational
compatibility target; other validators as a differential panel; **ClinVar as corpus, not
oracle.** Every intentional VEP deviation → an ERRATA witness. Build the adjudicated corpus
(repeats, exon boundaries, UTRs, introns, partial transcripts, frameshifts, extensions,
duplications, reference mismatches) before M6c.

## Identity vs HGVS coordinates

Keep the immutable left-normalized VariantKey for storage/joins (`bcftools_norm` left-aligns;
`functions.yaml`). Derive HGVS 3′-shift SEPARATELY per reference sequence and emit the
offset; for `c./n.` the shift follows transcript orientation (negative strand → opposite
genomic direction), `g.` follows genomic direction. **Never** feed HGVS coordinates back into
VariantKey — one genomic variant would acquire transcript-dependent identities and corrupt
joins.

## HGVS↔VCF is not a bijection

Many-to-one (padded / left-aligned VCF forms collapse to one 3′ HGVS) and one-to-many (one
genomic variant → many transcript HGVS; HGVSp discards nucleotide identity). Ship
**generation first**; promise *semantic application-equivalence* only for a declared
`g./c./n.` subset; never promise `p.→VCF` inversion or string identity after a VCF
round-trip. Test `apply(ref, VCF) == apply(ref, parsed HGVS)`, not string round-trips.

## M6a–M6d sequencing

- **M6a — contract gate (this note):** pin VEP 116 + HGVS snapshot; design Model A v2 +
  typed edit relations + provenance/hashes + the lossless `allele_transcript_context`.
  **No code folding.**
- **M6b — Model A + SO MVP:** reproducible GRCh37/38 imports proven by exact
  spliced/CDS/peptide hashes; VEP-116 SO for a declared biallelic small-variant subset, both
  strands, coding + noncoding. No HGVS / haplotype claim.
- **M6c — DNA HGVS:** generation-only `g./c./n.`, independent per-reference 3′ norm,
  application-equivalence tests.
- **M6d — protein + compound:** `p.` after the protein state machine passes independent
  differential tests; then translation edits, partial CDS, haplotypes, `m.`, structural
  HGVS as separately gated expansions.

Smallest defensible MVP matrix (a single plus-strand transcript is a *smoke test*, not an
MVP): both strands, phase, coding/noncoding, complete/partial CDS, repeats, and ≥1
transcript-reference mismatch.

## Memory (open — needs measurement)

HGVS needs cDNA + peptide + edit context; ~458 MB summed CDS / 1.45 GB summed exon bases
(116 GFF, pre-materialization). "gnomAD-flat RSS" is unproven until the relation cursor
genuinely streams. Pool / dedup sequences behind immutable blobs or bounded per-worker
caches; do not copy full transcript strings into every Model A row. Measure model-pack size /
peak RSS / per-worker cache / cold + warm throughput before promising memory bounds.

## Provenance of this note

Adversarial review: `codex exec` (`gpt-5.6-sol`, reasoning effort `max`), session
`019f4b64-b979-7283-a2df-2a4ac972f605`, 2026-07-10. All file-level claims above were
re-verified against the working trees before landing (Constants/VariationEffect diff,
`duckhgvs.h` limits, `duckvep_codon.h` tables, `fuzz_witnesses.R` exit path).
