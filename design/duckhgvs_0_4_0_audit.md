# Audit: `/root/duckhgvs-0.4.0` capability + provenance (verifying the M6a quarantine)

Status: **reference audit (2026-07-10).** Verified capability + provenance audit of the
stalled `/root/duckhgvs-0.4.0` spike, confirming and calibrating the "duckhgvs-0.4.0
quarantine" in `duckvep_m6a_contract_gate.md`. Sonnet 5 subagent; the buffer-break rate and
max-CDS figure were independently reproduced against the bundled GFF. Read when planning M6c
(DNA HGVS) — §5 says what to lift (the apply-then-diff test pattern + seed vectors) vs discard.

Scope: read-only audit of the stalled `/root/duckhgvs-0.4.0` spike against the claims in
`/root/duckhts/design/duckvep_m6a_contract_gate.md` ("duckhgvs-0.4.0: quarantine" + M6c
sections). All line numbers verified by direct read on 2026-07-10 unless marked **[inferred]**.
The target directory was not modified; only prebuilt binaries already present under
`build/` were executed (no `make`/rebuild was run).

---

## 1. Capability matrix

| Capability | Status | Evidence |
|---|---|---|
| Parse `g.` / `c.` / `n.` / `r.` / `p.` / `m.` / `o.` HGVS text | Implemented (grammar only) | `include/duckhgvs.h:44-53` (`duckhgvs_coord_t`), `src/duckhgvs.c:90-101` (`coord_from_char`), `src/duckhgvs.c:416-464` (`duckhgvs_parse`) |
| Format struct back to HGVS text | Implemented | `src/duckhgvs.c:518-563` (`duckhgvs_format`) |
| Edit grammar: sub/del/ins/delins/dup/inv/repeat/no-change | Implemented (text-level) | enum `include/duckhgvs.h:55-66`; parsers `src/duckhgvs.c:245-414` |
| Generate `g.`/`m.` from VCF-like allele | Implemented | `duckhgvs_from_vcf_like`, `src/duckhgvs.c:873-1005`. `o.` (circular) is declared in the enum but **never emitted** — `local_opts.coord` is restricted to `GENOMIC`/`MITO` only (`src/duckhgvs.c:906-909`); no `o.` test vector exists anywhere in `tests/`. |
| Generate `c.`/`n.`/`r.` from a transcript-coordinate allele | Implemented, **sequence-relative only** | `duckhgvs_from_transcript_vcf_like`, `src/duckhgvs.c:1037-1095`. Explicitly refuses intronic offsets: `duckhgvs_normalize` returns `ENOSUP` for `CDS/TX/RNA` coords with offsets (`src/duckhgvs.c:795-800`). |
| Generate `p.` from a transcript-coordinate coding allele | Implemented for the fully-in-CDS case only | `duckhgvs_protein_from_transcript_vcf_like`, `src/duckhgvs.c:1257-1328`. Alleles spanning the CDS boundary are rejected with `ENOSUP` (`src/duckhgvs.c:1290-1293`). Alleles **entirely outside the CDS** (UTR) silently return `p.=` (`src/duckhgvs.c:1287-1289`) — **[inferred]** this conflates "no protein change" (synonymous) with "not in the coding region at all," which is not standard HGVS usage; the built-in oracle rows `txp_utr5`/`txp_utr3` (`tests/oracle/hgvs_semantics.tsv:29-30`) codify this behavior as "expected," so it is a tested design choice, not an oversight, but it is a correctness gap if surfaced as-is. |
| 3′/5′ shuffling (normalization) | Implemented, base-at-a-time via the ref-provider callback | `shuffle_deletion_3p/5p`, `shuffle_insertion_3p/5p`, `src/duckhgvs.c:644-737`. One provider call per flanking base — architecturally fine, but chatty for a real FASTA/htslib backend over long repeats. |
| delins → sub/del/ins minimization | Implemented | `canonicalize_delins_minimal`, `src/duckhgvs.c:739-774` |
| Insertion→dup recognition | Implemented, opportunistic | `maybe_convert_insertion_to_dup`, `src/duckhgvs.c:844-871` |
| Reference-sequence validation | Implemented | `validate_ref_if_requested`, `src/duckhgvs.c:621-635` |
| HGVS → VCF-like conversion | Implemented for `g.`/`m.`, sub/del/ins/delins only | `duckhgvs_to_vcf_like`, `src/duckhgvs.c:1521-1609`; dup/inv/repeat explicitly `ENOSUP` (`src/duckhgvs.c:1607-1608) |
| **Apply edit to an explicit reference string** (`duckhgvs_apply_to_sequence`) | Implemented for sub/del/ins/delins/dup/inv; **not** repeat, **not** protein | `src/duckhgvs.c:1384-1507`; repeat rejected `1399-1402`; protein rejected `1395-1398` |
| **Sequence-level equivalence** (`duckhgvs_equivalent_on_sequence`) | Implemented (apply-then-diff) | `src/duckhgvs.c:1509-1519` |
| Genomic→transcript allele projection | Implemented, **exon-contained alleles only** | see §1a below |
| Intronic `c.`/`n.` genomic resolution | **Not implemented** — matches README's own "not yet" list | README.md "What this library is not yet" |
| Protein normalization / uncertain / predicted notation | **Not implemented** | `duckhgvs_normalize` hard-rejects protein coord (`src/duckhgvs.c:791-794`) |
| Structural variants / symbolic alleles (`<DEL>`, breakends) | Explicitly rejected everywhere | `sequence_allele_is_supported` (`src/duckhgvs.c:1030-1035`), symbolic-allele guards in `duckhgvs_from_vcf_like` (`891-894`), `duckhgvs_from_transcript_vcf_like` (`1051-1054`), `duckhgvs_protein_from_transcript_vcf_like` (`1275-1278`) |
| Codon table coverage | Standard (table 1) + vertebrate-mitochondrial (table 2) exceptions hardcoded in a `switch` | `codon_to_aa`, `src/duckhgvs.c:1097-1133`. No data-driven table registry — same species-limited smell flagged for `duckvep-c`'s `duckvep_codon.h` in the M6a note. |
| External dependencies | None — pure C99 + libc (`ctype/stdarg/stddef/stdint/stdio/stdlib/string/strings`, plus `time.h`/`windows.h` in the bench tool only) | verified via `grep '#include'` across `src/`, `include/`, `tools/`, `tests/` |
| Memory-safety discipline | All fixed buffers are guarded by explicit capacity checks (`safe_copy`, `appendf`, `append_n` all check `n >= cap` and return `DUCKHGVS_EBUFSIZE` rather than truncating/overflowing) | `src/duckhgvs.c:103-115, 466-481, 1354-1363` | No overflow bugs found; the failure mode on oversized input is a clean error return, not memory corruption. This is a genuine positive worth crediting. |

### 1a. The "wholly contained in one exon" restriction — verified

`duckhgvs_project_genomic_vcf_to_transcript` (`src/duckhgvs.c:1714-1757`) maps a genomic
allele's start and end independently through `model_map_exonic_pos` (`src/duckhgvs.c:1672-1681`,
itself walking a linear `model_find_exon` scan, `src/duckhgvs.c:1655-1660`):

```c
if (model_map_exonic_pos(model, genomic_pos1, &tx_s, &sidx) != 0 ||
    model_map_exonic_pos(model, genomic_end1, &tx_e, &eidx) != 0) {
    set_err(err, DUCKHGVS_ENOSUP, "allele is not wholly contained in an exon; ...");
    return -1;
}
if (sidx != eidx) {
    set_err(err, DUCKHGVS_ENOSUP, "allele spans exon boundary; split consequence context is required");
    return -1;
}
```
(`src/duckhgvs.c:1734-1742`)

So: (a) if either allele endpoint falls in an intron, projection fails; (b) if the two
endpoints resolve to *different* exon indices (the allele straddles a splice junction),
projection also fails. Both cases return `DUCKHGVS_ENOSUP`. The claim in the M6a note is
**verified precisely as stated** — projection only succeeds for alleles wholly inside a
single exon.

---

## 2. Hard limits

| Constant | Value | Where used | Where it breaks on real data |
|---|---|---|---|
| `DUCKHGVS_SEQ_MAX` | 4096 | `duckhgvs.h:25`. Backs: `variant_t.ref/alt/repeat_unit` (`duckhgvs.h:98-100`); `rbuf/abuf` allele scratch in `duckhgvs_from_vcf_like` (`duckhgvs.c:916-917`); reference-validation scratch (`duckhgvs.c:628`); dup-detection scratch (`duckhgvs.c:856`); and — the binding one — `ref_cds`/`obs_cds`/`ref_aa`/`alt_aa`/`rbuf`/`abuf`/`protein_edit` in `duckhgvs_protein_from_transcript_vcf_like` (`duckhgvs.c:1295,1298,1304,1309`), which fetches the **entire CDS** into a 4096-byte stack buffer (`duckhgvs.c:1296`) | **Independently reproduced against the bundled `GCF_000001405.25_GRCh37.p13_genomic.gff.gz`.** Grouping `CDS` features by `Parent`: 74,880 transcripts have a CDS; **4,316 (5.8%) have CDS length > 4096 nt**; **max CDS length = 107,976 nt** — this exact figure matches the M6a note's "max CDS 107,976" verbatim, corroborating the note's data source. Any of those 4,316 transcripts hard-fails `duckhgvs_protein_from_transcript_vcf_like` with `EBUFSIZE` from `duckhgvs_memory_ref_fetch` (`duckhgvs.c:587-588`) the moment the whole-CDS fetch is attempted — this is not a precision loss, it is a hard `p.` generation failure for those transcripts. If instead the metric is **spliced transcript length** (the more honest bound, since `n.`/apply/equivalence paths and any full-transcript-window fetch also route through `DUCKHGVS_SEQ_MAX`/`DUCKHGVS_FORMAT_MAX`-sized buffers), grouping `exon` features by `Parent` gives **28,901 of 108,760 transcripts (26.6%) with spliced length > 4096 nt**, max spliced length 109,224 nt. I could not reproduce the note's exact "9,742" figure from either metric alone on this file — it may come from a different GFF build (e.g. GRCh38/116) or a different counting rule — but under either metric I *did* compute, the practical conclusion is the same or worse than stated: the 4096 buffer excludes a double-digit percentage of real human transcripts, not an edge case. |
| `DUCKHGVS_TRANSCRIPT_EXON_MAX` | 512 | `duckhgvs.h:143`. Backs `duckhgvs_transcript_model_t.exons[512]` (`duckhgvs.h:161`); `duckhgvs_transcript_model_prepare` rejects (`EINVAL`) any model with `exon_count > 512` (`duckhgvs.c:1625-1628`) | **Does not break on this specific GRCh37 RefSeq GFF**: max observed exon count per transcript = 363 (consistent with *TTN*, the classic max-exon-count gene). This refines the M6a note, which lists this constant alongside `SEQ_MAX` without a concrete break — on the bundled dataset it is a latent risk (71% utilized at the extreme), not a demonstrated one. It would need a denser annotation (e.g. some GENCODE/Ensembl transcript models, or any caller building one combined "model" across multiple transcripts) to actually overflow. |
| `DUCKHGVS_ACCESSION_MAX` / `DUCKHGVS_SELECTOR_MAX` | 128 | accession/selector/gene_symbol/transcript_accession/protein_accession fields | Not a practical concern — real RefSeq/Ensembl versioned accessions (`NM_000088.3`, `ENST00000357654.9`) are well under 128 bytes. |
| `DUCKHGVS_MSG_MAX` | 512 | error message buffer only | Not a correctness limit. |
| `DUCKHGVS_FORMAT_MAX` | 8192 | `variant_t.raw`, `duckhgvs_format` output, and (critically) the caller-side `observed`/`lhs`/`rhs` buffers used by every built-in `APPLY`/`EQUIV` oracle test (`tools/duckhgvs_oracle.c:88,237,514-515`) | `duckhgvs_apply_to_sequence`'s `observed_cap` is caller-supplied so it is not an intrinsic library ceiling, but every shipped caller (CLI, oracle, tests) sizes it at 8192 — so *as actually used today*, no whole chromosome or even a long transcript can be run through the apply/equivalence oracle without the caller first slicing a local window. |
| Struct size / stack-copy cost | `duckhgvs_variant_t` ≈ 8192 (`raw`) + 3×128 (accession/selector/gene) + 3×4096 (`ref`/`alt`/`repeat_unit`) + interval/edit fields ≈ **~21 KB per variant** | `include/duckhgvs.h:90-104` | **[inferred]** Not a hard limit but a real perf/memory smell: the code repeatedly copies this struct by value (`*out = v;` `duckhgvs.c:1002`, `duckhgvs_variant_t tmp = *var;` `duckhgvs.c:1569,1597`), which is expensive at DuckDB row/chunk scale; nothing here is designed for vectorized bulk processing. |

---

## 3. Test-reality assessment — this is the crux

**Verdict: with one narrow exception, every "conformance"/"oracle" test in this repo is a
parse→format round-trip, a self-consistency check, or a hand-authored expected-value
comparison against toy synthetic sequences. It is not validation against an independent
HGVS-correctness oracle.** Quotes below, with file:line.

### 3a. VEP CSQ "oracle" — round-trip only, not correctness

`tools/duckhgvs_oracle.c:538-558` (`parse_one_hgvs`, shared verbatim with the legacy
`tools/duckhgvs_vepcheck.c:66-82`):

```c
static int parse_one_hgvs(const char *kind, const char *s, int lineno, long long *checked) {
    ...
    if (duckhgvs_parse(s, &v, &err) != 0) { ... return 1; }
    if (duckhgvs_format(&v, out, sizeof(out), &err) != 0) { ... return 1; }
    if (strcmp(out, s) != 0) {
        fprintf(stderr, "vep-csq:%d: %s roundtrip changed expected '%s' got '%s'\n", lineno, kind, s, out);
        return 1;
    }
    (*checked)++;
    return 0;
}
```

This takes VEP's `HGVSg`/`HGVSc`/`HGVSp` **output strings** from CSQ, parses them, formats
them back, and asserts the formatted text equals the original text. It never looks at the
record's own `REF`/`ALT`/`POS`, never computes an HGVS string independently and compares it
to VEP's, and never calls `duckhgvs_apply_to_sequence` against a reference. It proves the
grammar can round-trip real-world VEP output text; it proves **nothing** about whether
duckhgvs's own generation logic would produce the same (or a biologically correct) answer.
The README is honest about this in one place: "This is real-world grammar/round-trip
coverage, not proof of biological c./p. equivalence" (README.md, "VEP CSQ test" section) —
but elsewhere the same README calls 0.4.0 an "explicit C oracle layer" and says the package
"no longer treats 'tests' as parser smoke tests only," which overstates what the VEP lane
actually checks.

### 3b. ClinVar CLNHGVS — the one real external-oracle comparison, but SNV-only and trivial

`tools/duckhgvs_oracle.c:447-486` (`clinvar_clnhgvs`) *does* independently regenerate a `g.`
HGVS string from ClinVar's own `CHROM/POS/REF/ALT` via `duckhgvs_from_vcf_like` and string-
compares it against ClinVar's `INFO/CLNHGVS`:

```c
if (!simple_seq(f[3]) || !simple_seq(f[4]) || strchr(f[4], ',') || strlen(f[3]) != 1 || strlen(f[4]) != 1) { skipped++; continue; }
...
if (duckhgvs_from_vcf_like(acc, pos, f[3], f[4], NULL, &opts, &v, &err) != 0 || duckhgvs_format(&v, out, sizeof(out), &err) != 0) { ... failures++; continue; }
if (strcmp(out, cln) != 0) { fprintf(stderr, "clinvar:%d: CLNHGVS mismatch ...\n", lineno); failures++; }
```

The upstream DuckHTS-side extraction template hard-codes the same restriction:
`tests/duckhts/clinvar_clnhgvs_oracle.sql:23` — `WHERE length(REF) = 1 AND length(alt) = 1`.
So this is a real oracle comparison, but scoped to **pure SNVs only**, and only `g.`
coordinates — i.e. exactly the case where HGVS generation is closest to trivial
(`ACCESSION:g.POSREF>ALT`, no shifting/dup/delins/exon-projection logic exercised at all).
It says nothing about the harder logic (indel normalization, dup recognition, delins
minimization, transcript projection, protein translation) that is the actual value
proposition of an HGVS engine.

### 3c. `duckhgvs_vcfcheck` — crash smoke test, no expected-value comparison at all

`tools/duckhgvs_vcfcheck.c:35-51` (`check_alt`):

```c
if (duckhgvs_from_vcf_like(chrom, (int64_t)pos, ref, alt, NULL, &opts, &v, &err) != 0) {
    fprintf(stderr, "vcfcheck:%d: generation failed ...\n", lineno, ...); return 1;
}
if (duckhgvs_format(&v, out, sizeof(out), &err) != 0) { ... return 1; }
(*checked)++;
return 0;
```

There is no `expected` value anywhere in this file. "Passing" means only "did not return an
error code." This tool cannot detect a wrong-but-well-formed HGVS string.

### 3d. Local TSV "oracle" corpora — hand-authored expected values on toy sequences

`tests/oracle/core.tsv`, `tests/oracle/hgvs_semantics.tsv`, `tests/oracle/semantics.tsv`,
`tests/gold/generation.tsv` (25 + 41 + 30 + 11 = 107 data rows total) use small synthetic
reference strings such as `AACTTTTGA` (9 nt), `CCCATGGCTTAAGG` (14 nt), `ATGGTTGAA` (9 nt).
Example rows:

```
TX  coding_dup_frameshift  NM_TEST.1  NP_TEST.1  ATGGTTGAA  1  9  4  G  GG  none  NM_TEST.1:c.4dupG  NP_TEST.1:p.Val2fsTer2
```
(`tests/oracle/hgvs_semantics.tsv:37`)

```
APPLY  apply_inv  REF:g.3_5inv  AACTTTTGA  AAAAGTTGA
```
(`tests/oracle/semantics.tsv:18`)

The README's own framing is candid about the mechanism: "`APPLY` and `EQUIV` are the local
semantic oracle: they execute the parsed edit against an explicit reference sequence and
compare observed sequence" (README.md, "Semantic oracle implementation"). This is a
legitimately useful *self-consistency* check (it does catch, e.g., an off-by-one in
`duckhgvs_apply_to_sequence` if the author's hand-computed `expected` differs from what the
code produces) — but the "expected" column throughout is authored by whoever wrote the test
row, not computed by an independent HGVS engine (Mutalyzer/VariantValidator/biocommons-hgvs)
or validated against a real annotated dataset. Running the two prebuilt binaries confirms
they currently pass: `./build/test_duckhgvs` → `duckhgvs tests: ok`; `./build/duckhgvs_oracle
conformance tests/oracle/core.tsv` → `duckhgvs oracle conformance: 24 ok`. "All tests pass"
here means "the code agrees with what its own author wrote down it should produce."

### 3e. Unit tests (`tests/test_duckhgvs.c`, `tests/test_gold.c`)

Same pattern: `CHECK(strcmp(out, "REF:g.7dupT") == 0);` (`tests/test_duckhgvs.c:110`) — every
assertion is a hand-written expected literal compared against the library's own output, on
hand-picked short synthetic sequences (`AACTTTTGA`, `CCCATGGCTTAAGG`, `GAT`, `ATGGTTGAA`).
No case draws from an external, independently-computed answer.

### Bottom line for §3

Out of every test surface in the repo, exactly one (`clinvar-clnhgvs`) compares duckhgvs's
own output against an independently produced ground truth, and it is restricted to the
easiest possible case (SNV, genomic coordinate, no normalization). Everything else —
including the file literally named `duckhgvs_oracle` and the README's claim of an "explicit
C oracle layer" — is round-trip/self-consistency testing dressed in oracle terminology. This
matches and sharpens the M6a note's framing; if anything the note is generous in not
spelling out the VEP-lane round-trip mechanism explicitly.

---

## 4. Provenance findings

- **License**: `LICENSE` is plain MIT, "Copyright (c) 2026 duckhgvs contributors" (`LICENSE:1-3`). It covers the *software* only; there is no NOTICE file and no mention anywhere of the bundled data file's separate origin/terms.
- **The bundled GFF is provenance-unattributed and functionally orphaned.** `GCF_000001405.25_GRCh37.p13_genomic.gff.gz` (48,542,794 bytes) is present at repo root. Its own header identifies it precisely: `#!annotation-source NCBI RefSeq GCF_000001405.25-RS_2024_09`, `#!genome-build GRCh37.p13`, `#!genome-build-accession NCBI_Assembly:GCF_000001405.25` (verified via `zcat | head`). **No source file, test, Makefile target, or CMake target anywhere in the repository references this filename or `GCF_000001405` at all** (verified by `grep -rn` across every `.c/.h/.md/.sql/.tsv/Makefile/CMakeLists.txt`, excluding `build/`). It is not decompressed, not parsed, not wired into any oracle or CLI path — it is 48 MB of dead weight in the working tree with zero functional connection to the code that ships alongside it.
- **Licensing assessment**: RefSeq/NCBI annotation data is a US-government work and is not copyrighted in the US (17 U.S.C. §105), so redistributing it is not a copyright *violation* in the strict legal sense. But it is still a provenance gap worth fixing before folding anything from this repo: there is no NOTICE distinguishing "MIT license (our code)" from "NCBI RefSeq data (their data, cite/attribute per NCBI's usage expectations)," and best practice for any tool bundling reference annotation is to say so explicitly. As it stands the LICENSE file is silent about the data file, which reads as an oversight rather than a considered decision.
- **Not a git repository**: `/root/duckhgvs-0.4.0` has no `.git` directory (`git rev-parse --is-inside-work-tree` fails with "not a git repository"). There is no commit history to audit provenance against, no way to see when/why the GFF was added, and no way to distinguish accidental staging from intentional bundling. This alone is consistent with "stalled spike," not a maintained project.
- **No third-party code found.** Grep for `Copyright`, `SPDX`, `Ensembl`, `biocommons`, `Mutalyzer`, `hgvs\b` across all `.c`/`.h` files (excluding `build/`) turns up nothing except the README/header's own *disclaimer* that it does **not** embed VEP/bcftools/Mutalyzer/biocommons code (`include/duckhgvs.h:7-9`). I found no unattributed Ensembl or other third-party source inside `src/`/`tools/`/`tests/` — the code itself appears to be original, dependency-free C99. The provenance problem is entirely the orphaned data file, not embedded foreign code.
- `tests/duckhts/grch38_accessions.tsv` is a hand-maintained chrom→`NC_0000XX.NN` accession map for **GRCh38**, used only by the (SNV-only) ClinVar oracle template — it is unrelated to, and not cross-checked against, the bundled GRCh37 GFF. Nothing in the repo ties the bundled GRCh37 GFF to any test that actually runs.

---

## 5. Salvage vs. discard recommendation, mapped to M6c needs

M6c needs (from the M6a note): independent per-reference 3′ normalization; accession
versions; application-equivalence testing (`apply(ref, VCF) == apply(ref, parsed HGVS)`);
and, separately (M6a scope), the lossless `allele_transcript_context` object.

| M6c/M6a need | What duckhgvs-0.4.0 has | Verdict |
|---|---|---|
| Independent per-reference 3′ normalization | `duckhgvs_ref_provider_t` callback (`duckhgvs.h:107-119`) + `shuffle_deletion_3p/5p`/`shuffle_insertion_3p/5p` (`duckhgvs.c:644-737`) decouple "which reference sequence" from "how to shift," via an opaque `(accession, start1, end1) -> sequence` callback per call. | **Salvage the design idea, not the code.** The interface shape (accession-scoped provider callback, separate 3′/5′ shuffle entry points) is a reasonable, testable pattern for keeping genomic vs. transcript normalization independent, which is exactly M6c's requirement. The implementation is not: it fetches one base at a time in a loop (`simple_fetch(..., e+1, e+1, ...)` per iteration), which is fine for a 9-nt test string and would be a chatty, slow anti-pattern against a real htslib/FASTA-backed provider over a real homopolymer/STR run. Rewrite against windowed fetches; keep the "provider is scoped per accession, shuffle is provider-driven" shape as a reference. |
| Accession versions | Accession is stored as an opaque, unvalidated string field (`DUCKHGVS_ACCESSION_MAX=128`). No version parsing, no comparison, no crosswalk. | **Nothing to salvage.** This is a plain string field; M6c's actual need (MANE/RefSeq crosswalk, sequence-hash-verified accession/version identity per the M6a Model A v2 decision) is entirely absent here and would have to be built from scratch regardless. |
| Application-equivalence testing (`apply(ref,VCF)==apply(ref,parsedHGVS)`) | `duckhgvs_apply_to_sequence` + `duckhgvs_equivalent_on_sequence` (`duckhgvs.c:1384-1519`), exercised by the `APPLY`/`EQUIV` TSV rows in `tests/oracle/semantics.tsv`. | **This is the single most salvageable asset in the repo.** It is precisely the test *shape* the M6a note asks for: execute a parsed edit against an explicit reference string and diff the observed sequence, rather than comparing formatted strings. Recommend lifting the ~20 `APPLY`/`EQUIV` rows (homopolymer-shift equivalence, sub/del/ins/delins/dup/inv application) as a day-one regression seed for M6c's own apply-based test harness, and reusing the "apply-then-diff, not format-then-diff" design principle. Do not reuse the C functions themselves (fixed 4096/8192-byte buffers, `strcmp`-only equivalence, no repeat/protein support) — reimplement against DuckHTS's own memory model. |
| Lossless `allele_transcript_context` | `duckhgvs_projected_allele_t` (`duckhgvs.h:150-156`: `tx_pos1`, `ref`, `alt`, `exon_index`, `is_intronic`) and `duckhgvs_transcript_model_t` (`duckhgvs.h:158-163`). | **Discard as a data-model reference — it is the opposite of what M6a asks for.** It carries one scalar transcript position and no genomic range, no CDS/peptide ranges, no sequence-provider identity, no edit provenance, no normalization offsets — and it *hard-fails* (`ENOSUP`, `duckhgvs.c:1734-1742`) on exactly the cases Model A v2 must represent losslessly: intronic-adjacent and exon-boundary-crossing alleles. If anything this is a useful **negative example**: a concrete illustration of the too-narrow, single-scalar-position design the M6a note is explicitly steering away from. |
| HGVS text grammar (parse/format `g./c./n./r./p.`) | `duckhgvs_parse`/`duckhgvs_format` (`duckhgvs.c:416-563`), dependency-free C99. | **Modest, conditional salvage.** If DuckHTS needs a standalone HGVS *text* parser/formatter independent of correctness-generation (e.g. to ingest/round-trip ClinVar `CLNHGVS` or VEP CSQ `HGVSc`/`HGVSp` strings), this grammar layer is a plausible starting reference or source of test strings (it already round-trips real BRAF/ClinVar-style accessions). Known gaps to fix before any reuse: `p.=` misused for UTR/non-coding alleles (`duckhgvs.c:1287-1289`, tested-for in `tests/oracle/hgvs_semantics.tsv:29-30`); `summarize_protein_change` (`duckhgvs.c:1195-1255`) has dead no-op branches ("keep compiler quiet" comment at `duckhgvs.c:1207`) and no stop-lost/start-lost handling; codon table is a hardcoded `switch`, not data-driven (same species-limited smell flagged for duckvep-c). Per the M6a note's own instruction, **do not fold the ABI** — at most, port ideas/test strings. |
| Provenance / bundled data | 48 MB unattributed, unreferenced RefSeq GFF; no `.git` history. | **Discard outright.** It is not wired into anything, has no NOTICE, and there is no history explaining why it's there. Do not carry it forward into DuckHTS in any form; if a GRCh37 RefSeq GFF is needed for M6c fixtures, re-derive/re-vendor it deliberately with a NOTICE and a pinned URL/checksum per the repo's own vendoring rules (`AGENTS.md` "Pin exact upstream versions/commits and checksum-verify downloads"). |

### Overall

The M6a note's characterization — "a design spike, not a validated engine" — holds up under
direct source reading, including on points the note stated only qualitatively. Where I have
new information beyond the note: (1) the `SEQ_MAX=4096` CDS-length break is real and I
independently reproduced its most extreme data point (max CDS 107,976 nt) exactly, which
corroborates the note's methodology even though I could not reproduce its specific
"9,742 transcripts" count from the bundled GFF using either a CDS-length or spliced-length
metric — under the spliced-length metric (arguably the more relevant one, since apply/`n.`/
equivalence paths also route through fixed-size buffers) the real number is worse: 26.6% of
transcripts, not under 10%; (2) the `TRANSCRIPT_EXON_MAX=512` limit, listed by the note
alongside `SEQ_MAX` without a concrete break, does **not** actually break on this GRCh37
dataset (max 363 exons) — a place where I'd soften the note slightly; (3) the "not a
validated engine" conclusion is, if anything, understated for the test layer specifically:
the VEP-lane code is a byte-for-byte parse→format round-trip with no reference to REF/ALT/
transcript at all, and the repo's own naming ("oracle," "semantic oracle layer") oversells
that. The quarantine recommendation in the M6a note is correctly calibrated, not too harsh.
