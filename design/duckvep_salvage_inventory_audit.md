# DuckVEP salvage inventory audit

Status: **reference audit (2026-07-10).** Verified, evidence-backed companion to
`duckvep_reconciliation.md`'s "Salvage inventory" section — file-by-file existence, line
counts, per-file dedupe diffs against DuckHTS, and gaps. Produced by a Sonnet 5 subagent and
spot-verified against the trees (`bcf_reader.c` divergence, the three gap files, the SIMD
kernel delta). Use when planning the actual duckvep-c fold (roadmap M6b); corrections from
this audit are already folded into `duckvep_reconciliation.md`.

Audited against `/root/duckvep-c` at commit `9f922c88f4feb618165778e77c3158ed591e76d7` (2026-07-07,
matches the `9f922c8` short hash cited in `duckvep_reconciliation.md:4`), working tree has minor
untracked/modified conformance data files only (no source changes). Compared against
`/root/duckhts` (this repo, branch `develop`). All commands and diffs below were run directly
against the filesystem; nothing here is inferred from memory of the design docs.

Methodology for Tier C dedupe verdicts: raw `diff`, then a second pass with
`sed 's/duckvep/duckhts/g; s/DUCKVEP/DUCKHTS/g'` applied to the duckvep-c copy before re-diffing,
since nearly every file differs *only* by the `duckvep_`/`DUCKVEP_` vs `duckhts_`/`DUCKHTS_`
identifier prefix (and the SQL-visible `dv_` function-name prefix, kept for co-load safety per
`functions.yaml:7-8` in duckvep-c). A file is called NEAR-DUP only if the post-rename diff is
empty or trivial; SUBSTANTIALLY DIVERGED if real logic differs after normalizing the prefix.

## (a) Corrected inventory table

### Tier A — irreplaceable (conformance/ + kernel/ + design docs)

| File | Tier | Exists? | Lines | Dedupe verdict | Notes |
|---|---|---|---|---|---|
| `conformance/` (all 39 files, incl. `data/`) | A | yes | n/a | N/A (unique to duckvep-c) | Verified full tree; matches "in full" framing. See explicit sub-list below. |
| `conformance/generate_effect_rules.pl` | A | yes | 97 | N/A | |
| `conformance/extract_so_spec.pl` | A | yes | 28 | N/A | |
| `conformance/generate_so_metadata.pl` | A | yes | 101 | N/A | |
| `conformance/so_dump.c` | A | yes | 34 | N/A | |
| `conformance/model_a_from_gff.sql` | A | yes | 131 | N/A | |
| `conformance/model_a_with_cds.R` | A | yes | 181 | N/A | |
| `kernel/src/duckvep_classify.c` | A | yes | 298 | N/A (no DuckHTS counterpart) | |
| `kernel/src/duckvep_effect.c` | A | yes | 212 | N/A | |
| `kernel/src/duckvep_sweep.c` | A | yes | 197 | N/A | |
| `kernel/src/duckvep_sv.c` | A | yes | 85 | N/A | |
| `kernel/src/duckvep_delta.c` | A | yes | 2045 | N/A | |
| `kernel/src/duckvep_codon.c` | A | yes | 59 | N/A | |
| `kernel/src/duckvep_coding.c` | A | yes | 136 | N/A | |
| `kernel/src/duckvep_haplotype.c` | A | yes | 254 | N/A | |
| `kernel/src/duckvep_projection.c` | A | yes | 193 | N/A | |
| `kernel/src/duckvep_kernel.c` | A | yes | 1054 | N/A | |
| `kernel/src/duckvep_so.c` | A | yes | 92 | N/A | |
| `kernel/include/duckvep_kernel.h` | A | yes | 396 | N/A | |
| `kernel/include/duckvep_so.h` | A | yes | 122 | N/A | |
| `design/duckvep_effect_ctx_architecture.md` | A | yes | 277 | N/A (unique) | |
| `design/vep_consequence_state_machine.md` | A | yes | 122 | N/A | |
| `design/duckvep_annotate_input_contract.md` | A | yes | 283 | N/A | |
| `design/duckvep_span_sv_cnv_consolidation_2026-06-23.md` | A | yes | 303 | N/A | |
| `design/duckvep_current_design_state_machines.md` | A | yes | 323 | N/A | Companion `.qmd` (323 source lines) exists but is unlisted — see gap list. |
| `design/duckvep_bcftools_csq_port_plan_2026-06-09.md` | A | yes | 194 | Already copied into duckhts | `/root/duckhts/design/duckvep_bcftools_csq_port_plan_2026-06-09.md` already exists (per `design/README.md` map) — confirm it's the same content before "salvaging" again. |
| `design/duckvep_layer_keys.md` | A | yes | 163 | Already copied into duckhts | `/root/duckhts/design/duckvep_layer_keys.md` already exists — same note. |

Verdict: **every Tier A file exists exactly as named.** No misnamed or missing files in Tier A.
Two of the seven "design docs" (`duckvep_bcftools_csq_port_plan_2026-06-09.md`,
`duckvep_layer_keys.md`) are **already present in `/root/duckhts/design/`** (confirmed by listing
`/root/duckhts/design/`), so calling them "salvage" targets is stale bookkeeping — they were
already lifted in an earlier pass. (Not diffed byte-for-byte in this audit; flagged for a
follow-up content check before assuming they're current.)

### Tier B — engine glue

| File | Tier | Exists? | Lines | Dedupe verdict | Notes |
|---|---|---|---|---|---|
| `src/consequence_udf.c` | B | yes | 95 | N/A (unique) | No dedicated header (verified via `#include`); "their include/ headers" doesn't apply to this file. |
| `src/annotate_table.c` | B | yes | 511 | N/A | No dedicated header either; includes `duckvep_chunk_reader.h`, `duckvep_model_reader.h`, `duckvep_relation_cursor.h`, `duckvep_variant_tile.h`, `duckvep_kernel.h`, `duckvep_so.h` (all separately listed). |
| `src/duckvep_variant_tile.c` | B | yes | 399 | N/A | |
| `src/duckvep_model_reader.c` | B | yes | 552 | N/A | |
| `src/duckvep_relation_cursor.c` | B | yes | 120 | N/A | |
| `src/duckvep_chunk_reader.c` | B | yes | 52 | N/A | |
| `src/chrom_interner.c` | B | yes | 150 | N/A | |
| `src/include/duckvep_variant_tile.h` | B | yes | 88 | N/A | |
| `src/include/duckvep_model_reader.h` | B | yes | 77 | N/A | |
| `src/include/duckvep_relation_cursor.h` | B | yes | 65 | N/A | |
| `src/include/duckvep_chunk_reader.h` | B | yes | 41 | N/A | |
| `src/include/duckvep_chrom_interner.h` | B | yes | 46 | N/A | Header is named `duckvep_chrom_interner.h`, not `chrom_interner.h` — the reconciliation doc's phrase "+ their include/ headers" doesn't spell out the name, so this isn't a doc error, just worth recording exactly. |

Verdict: **every Tier B file exists.** No misnamed/missing files.

### Tier C — claimed near-duplicates (verified by diff against DuckHTS)

| File | Exists in duckvep-c? | duckvep-c lines | duckhts lines | Verdict | Notes |
|---|---|---|---|---|---|
| `src/bcf_reader.c` | yes | 4551 | 5028 | **SUBSTANTIALLY DIVERGED** | Raw diff 559 changed lines even after normalizing `duckvep_`→`duckhts_`/`DUCKVEP_`→`DUCKHTS_`. DuckHTS's copy has a whole extra feature: `bcf_decode_error_policy_t` (`null`/`warn`/`error`), `bcf_check_decode_ret`, `bcf_encoded_type_matches_*`, `bcf_reader_input_is_bcf`, GT-preflight checks, and a `decode_error_policy` named parameter on `read_bcf`/`read_bcf_appender`. duckvep-c's copy is an **older fork missing this DuckHTS feature**, not a same-generation near-duplicate. **Correction needed: this is not safe to treat as "just drop the duckvep-c copy," but the DIRECTION is still "keep DuckHTS's, discard duckvep-c's" — duckhts is strictly ahead.** |
| `src/cgranges_api.c` | yes | 2664 | 2664 | **IDENTICAL** (after prefix rename) | 0 diff lines post-normalization. True near-dup. |
| `src/variantkey_udf.c` | yes | 1096 | 1096 | **NEAR-DUP** (52 lines differ) | All remaining diffs are `dv_`-prefixed SQL function names (`dv_decode_variantkey` etc.) vs unprefixed DuckHTS names (`decode_variantkey`) — cosmetic/naming only, same logic. |
| `src/vep_parser.c` | yes | 583 | 583 | **IDENTICAL** | 0 diff even *without* any renaming — this file has zero `duckvep_`-prefixed identifiers; it's byte-for-byte the same file in both repos already. |
| `src/bgzip.c` | yes | 380 | 380 | **NEAR-DUP** (8 lines) | Only the `dv_bgzip`/`dv_bgunzip` vs `bgzip`/`bgunzip` SQL names differ. |
| `src/tabix_reader.c` | yes | 1884 | 1884 | **NEAR-DUP** (14 lines) | Only `dv_read_tabix`/`dv_read_gtf`/`dv_read_gff` vs unprefixed names differ. |
| `src/seq_reader.c` | yes | 1334 | 1334 | **NEAR-DUP** (10 lines) | Only `dv_read_fasta`/`dv_fasta_index`/`dv_read_fastq` vs unprefixed names, plus one error-message string differ. |
| `src/quality_encoding.c` | yes | 231 | 231 | **IDENTICAL** (after prefix rename) | 0 diff post-normalization. |
| `src/bcftools_norm_udf.c` | yes | 1360 | 1360 | **NEAR-DUP** (4 lines) | Only `dv_bcftools_norm_row` vs `bcftools_norm_row` differs. |
| `src/bcftools_shim.c` | yes | 94 | 94 | **IDENTICAL** (after prefix rename) | 0 diff post-normalization. |
| `src/wasm_http_hfile.c` | yes | 556 | 556 | **IDENTICAL** (after prefix rename) | 0 diff post-normalization. |
| `src/include/seq_encoding.h` | yes | 87 | 87 | **IDENTICAL** (raw, no rename needed) | 0 diff, no `duckvep_`-prefixed symbols at all. |
| `src/include/wasm_socket_compat.h` | yes | 54 | 54 | **NEAR-DUP** (46 lines, prefix-only) | Every diff line is `duckvep_wasm_*`/`DUCKVEP_*` → `duckhts_wasm_*`/`DUCKHTS_*` renaming. |
| SIMD dispatchers (`src/simd/duckvep_simd_*.c` → `duckhts_simd_*`) | yes | see below | see below | **MIXED — 3 identical, 3 diverged** | See SIMD sub-table. |

**SIMD dispatcher detail** (Tier C says "merge the SIMD dispatchers... one dispatch table"):

| Backend | duckvep-c lines | duckhts lines | Post-rename diff | Verdict |
|---|---|---|---|---|
| `duckvep_simd_avx512.c` | 113 | 113 | 0 | IDENTICAL |
| `duckvep_simd_neon.c` | 97 | 97 | 0 | IDENTICAL |
| `duckvep_simd_wasm_simd128.c` | 105 | 105 | 0 | IDENTICAL |
| `duckvep_simd_dispatch.c` | 915 | 985 | 72 | DIVERGED |
| `duckvep_simd_scalar.c` | 46 | 123 | 79 | DIVERGED |
| `duckvep_simd_avx2.c` | 111 | 238 | 129 | DIVERGED |
| `include/duckvep_simd_kernels.def` | 12 | 14 | 2 kernel rows | DIVERGED |

`duckhts_simd_kernels.def` in DuckHTS defines **3** kernels (`SEQ_BASE_COUNTS`, `BAM_NT16_COUNTS`,
`NT16_GC_COUNTS` — the latter two added per the most recent DuckHTS commit,
`4222c22 docs(simd): add backend x kernel coverage ledger; record the NEON/wasm nt16 gap`).
duckvep-c's `.def` only has `SEQ_BASE_COUNTS`. **Correction:** the doc's "merge the SIMD
dispatchers" undersells that DuckHTS's `dispatch`/`scalar`/`avx2` backends are **strictly ahead**
of duckvep-c's (2 extra kernels DuckHTS already carries that duckvep-c never added), and
duckvep-c's simd/ directory contributes **zero unique kernels**. The correct action is "discard
duckvep-c's simd/ entirely, keep DuckHTS's" — not a symmetric merge.

### Missing-header / unmentioned-header check inside Tier C's own file list

Tier C names only two headers (`seq_encoding.h`, `wasm_socket_compat.h`) even though five of its
eleven listed `.c` files have paired headers that DID have to be checked to make the dedupe claim
meaningful. All exist and are near-dup/identical:

| Header | duckhts counterpart | Verdict |
|---|---|---|
| `src/include/quality_encoding.h` | `/root/duckhts/src/include/quality_encoding.h` | IDENTICAL after prefix rename |
| `src/include/vep_parser.h` | `/root/duckhts/src/include/vep_parser.h` | IDENTICAL (0 diff, raw) |
| `src/include/wasm_http_hfile.h` | `/root/duckhts/src/include/wasm_http_hfile.h` | IDENTICAL (0 diff, raw) |
| `src/include/bcftools_shim.h` | `/root/duckhts/src/include/bcftools_shim.h` | NEAR-DUP, prefix-only |
| `src/include/variantkey_compat.h` | `/root/duckhts/src/include/variantkey_compat.h` | NEAR-DUP, prefix-only |
| `src/include/filter.h` | `/root/duckhts/src/include/filter.h` | IDENTICAL (0 diff, raw) |
| `src/include/hts_io_tuning.h` | `/root/duckhts/src/include/hts_io_tuning.h` | NEAR-DUP, prefix-only |
| `src/include/bcftools_filter_config.h` | `/root/duckhts/src/include/bcftools_filter_config.h` | NEAR-DUP, prefix-only (2-line file, both lines differ) |
| `src/include/vcf_types.h` | `/root/duckhts/src/include/vcf_types.h` | IDENTICAL (0 diff, raw) |

None of these nine headers are missing or misnamed; they're simply **absent from the Tier C
bullet list** even though the `.c` files that need them are present. This is an inventory
completeness gap, not a correctness error (they will obviously be needed for a working dedupe).

## (b) Corrections needed in `duckvep_reconciliation.md`

- **`src/bcf_reader.c` is miscategorized.** It's listed in Tier C as a same-generation
  "near-duplicate," but DuckHTS's copy has a `decode_error_policy` (`null`/`warn`/`error`)
  feature, GT-preflight/type-checking helpers, and a new `read_bcf`/`read_bcf_appender`
  named parameter that duckvep-c's copy entirely lacks (559 diff lines survive prefix
  normalization, vs. 0–14 for every other Tier C `.c` file). State explicitly: "keep DuckHTS's
  `bcf_reader.c`, discard duckvep-c's — duckvep-c's fork predates the decode-error-policy work."
- **"Merge the SIMD dispatchers" overstates the effort/symmetry.** `duckvep_simd_avx512.c`,
  `_neon.c`, `_wasm_simd128.c` are byte-identical to DuckHTS's after prefix rename (nothing to
  merge). `_dispatch.c`, `_scalar.c`, `_avx2.c`, and `duckvep_simd_kernels.def` differ only because
  **DuckHTS added `BAM_NT16_COUNTS` and `NT16_GC_COUNTS`** after duckvep-c forked — duckvep-c
  contributes zero unique kernels. Rewrite as: "discard duckvep-c's `src/simd/`; it has nothing
  DuckHTS's `src/simd/duckhts_simd_*` doesn't already have."
- **Tier C's header list is incomplete.** It names `include/seq_encoding.h` and
  `include/wasm_socket_compat.h` only, but the eleven `.c` files it lists pull in seven more
  paired headers (`quality_encoding.h`, `vep_parser.h`, `wasm_http_hfile.h`, `bcftools_shim.h`,
  `variantkey_compat.h`, `filter.h`, `hts_io_tuning.h`, `bcftools_filter_config.h`, `vcf_types.h`
  — nine total). All verified present and near-dup/identical in DuckHTS, so this is a documentation
  gap, not a blocking discovery, but the fold plan should enumerate them so "dedupe Tier C" isn't
  read as "dedupe just these 13 files."
- **Tier C omits three files that are also verified near-duplicates and belong in it:**
  `src/hts_index_builder.c` (403/403 lines, 10-line raw diff, prefix-only), `src/hts_meta_reader.c`
  (1825/1825 lines, 36-line raw diff, prefix-only), and `src/interval_udf.c` (1078/1078 lines,
  30-line raw diff, prefix-only) all exist under the same name in `/root/duckhts/src/` and are
  clearly the same near-duplicate class as the rest of Tier C. See gap list below — this is the
  single biggest correction (three substantial files with a clean, mechanical dedupe story that
  the plan currently doesn't mention at all).
- **`src/duckvep.c` (402 lines, the extension's `DUCKDB_EXTENSION_EXTERN` entry point /
  registration file) is not mentioned in any tier.** It's the file that wires every Tier A/B/C
  function into the `dv_`/`duckvep_`-prefixed SQL surface. Since fuse-mechanics step 5 says the
  eventual surface is "one `duckhts_csq(...)` / annotate table function," this file is presumably
  meant to be superseded/rewritten rather than lifted verbatim — but the doc should say so
  explicitly instead of leaving it as a silent omission; a reader following the tier lists
  literally would not know this file exists or needs handling.
- **Two of the seven Tier A "design docs" are stale entries** —
  `design/duckvep_bcftools_csq_port_plan_2026-06-09.md` and `design/duckvep_layer_keys.md` already
  exist under `/root/duckhts/design/` (confirmed present by listing; not diffed for content
  drift in this pass). Listing them under "before archiving duckvep-c, lift into src/duckvep/"
  reads as if they still need to be copied; they don't (or if content has drifted, the doc should
  say "reconcile," not "lift").
- **Kernel private headers/includes are absent from Tier A**, even though the seven kernel `.c`
  files Tier A names cannot compile without them: `kernel/src/duckvep_classify.h`,
  `duckvep_coding.h`, `duckvep_codon.h`, `duckvep_delta.h`, `duckvep_effect.h`, `duckvep_haplotype.h`,
  `duckvep_projection.h`, `duckvep_sv.h`, `duckvep_sweep.h`, `duckvep_event.h`,
  `duckvep_workspace_internal.h`, plus two generated `.inc` files (`duckvep_effect_rules.inc`,
  `duckvep_so_metadata.inc`) and `kernel/CMakeLists.txt`. Tier A's header line ("+
  `include/{duckvep_kernel,duckvep_so}.h`") only covers the two *public* headers in
  `kernel/include/`, not the eleven *private* headers colocated in `kernel/src/`. This is worth
  fixing explicitly since "kernel/ — ... `duckvep_classify.c`, ... + include/{...}.h" reads as a
  complete file list and isn't.

## (c) Gap list — files under `/root/duckvep-c/{src,kernel,conformance}` not named in any tier

| File | Guessed tier | Basis |
|---|---|---|
| `src/duckvep.c` (402 lines) | Not salvageable as-is / needs its own decision | Main extension entry point; registers every Tier A/B/C function under `dv_`/`duckvep_` names. Superseded by the new `duckhts_csq` registration per fuse-mechanics step 5, but the doc never says this explicitly. |
| `src/hts_index_builder.c` (403 lines) | **Tier C** (near-dup, prefix-only, 10-line diff) | Same name and near-identical content already in `/root/duckhts/src/hts_index_builder.c`. |
| `src/hts_meta_reader.c` (1825 lines) | **Tier C** (near-dup, prefix-only, 36-line diff) | Same name/near-identical content in `/root/duckhts/src/hts_meta_reader.c`. |
| `src/interval_udf.c` (1078 lines) | **Tier C** (near-dup, prefix-only, 30-line diff) | Same name/near-identical content in `/root/duckhts/src/interval_udf.c`. |
| `src/include/bcftools_filter_config.h`, `bcftools_shim.h`, `filter.h`, `hts_io_tuning.h`, `quality_encoding.h`, `variantkey_compat.h`, `vcf_types.h`, `vep_parser.h`, `wasm_http_hfile.h` (9 headers) | **Tier C** (paired headers for already-listed `.c` files) | See "(b)" — all verified near-dup/identical, just missing from the bullet list. |
| `kernel/src/*.h` (11 private headers) + `kernel/src/*.inc` (2 generated includes) + `kernel/CMakeLists.txt` | **Tier A** (compile-time dependencies of the named kernel `.c` files) | No DuckHTS counterpart exists (unique to duckvep-c); needed to build the kernel at all. |
| `design/duckvep_current_design_state_machines.qmd` (Quarto source, companion to the `.md` already in Tier A) | **Tier A** (source of an already-listed doc) | Rendered `.md` is in Tier A; the `.qmd` it's rendered from is not mentioned. |
| `docs/HANDOFF.md` (239 lines), `docs/VEP_SOURCE_AUDIT.md` (98 lines) | **Historical/provenance, not source** — recommend archiving as a reference, not folding into `src/duckvep/` | `HANDOFF.md` is directly quoted as evidence in `duckvep_m6a_contract_gate.md:13` ("the stalled handoff") and `duckvep_reconciliation.md:37,256`; it and `VEP_SOURCE_AUDIT.md` have real provenance value (what was tried, why it stalled, VEP-source citations) but aren't code to lift into `src/duckvep/`. Worth a one-line disposition decision (keep a copy under `design/` vs. let it die with the archived repo) rather than silence. |
| `design/README.md`, `better_scans.md`, `duckdb_c_api_deprecation_scan_2026-04-21.md`, `multireading.md`, `simd_dispatch_matrix.md`, `simd_future_kernel_proposals.md`, `tigerstyle_audit_2026-06-08.md`, `read_bcf_v2_optimization_comparison_2026-06-09.md` | **Correctly excluded — verified, not a gap** | All eight already exist under `/root/duckhts/design/` with the same names. Five are byte-identical (`better_scans.md`, `duckdb_c_api_deprecation_scan_2026-04-21.md`, `multireading.md`, `simd_future_kernel_proposals.md`, plus `tigerstyle_audit_2026-06-08.md`/`read_bcf_v2_optimization_comparison_2026-06-09.md` present by name, not diffed). `simd_dispatch_matrix.md` differs by 33 lines (DuckHTS's is newer, post-nt16-gap-ledger). None need salvaging; DuckHTS's copies are current or ahead. Recorded here only so the "no-mention = overlooked" reading doesn't get raised as a false gap later. |

Conformance's `data/` subtree (clinvar parquet/csv snapshots, `so_consequences.tsv`,
`effect_rule_bindings.tsv`, etc.) is covered by Tier A's "conformance/ in full" wording — verified
present, not a gap.

## Summary of verification scope

- Tier A: 13 kernel files + 6 explicitly-named conformance files + 7 design docs = all 26
  individually named files confirmed present with line counts recorded above. `conformance/`
  "in full" verified as 39 files, all present.
- Tier B: 7 `.c` + 5 headers (2 of the 7 `.c` files have no header) = all 12 confirmed present.
- Tier C: 11 `.c` + 2 headers = all 13 confirmed present and diffed against DuckHTS; 1 file
  (`bcf_reader.c`) reclassified from NEAR-DUP to SUBSTANTIALLY DIVERGED; SIMD sub-claim corrected
  (3 of 6 backends have zero unique content, the other 3 have zero unique content *for DuckHTS's
  favor* — duckvep-c owns nothing SIMD-side that isn't already superseded).
- Gap check found 3 substantial `.c` files (~3300 lines combined) and 9 headers that belong in
  Tier C but aren't listed, 13 kernel-internal files that belong in Tier A but aren't listed, and
  one un-tiered 402-line registration file (`duckvep.c`) whose disposition the plan never states.
