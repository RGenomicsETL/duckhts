# DuckHTS Roadmap

Status: **future proposal / open design.** Converged roadmap (reconciled between an
in-session analysis and an independent gpt-5.6-sol pass, 2026-07). Milestones are
sequenced but not yet started; each names its enabling primitive, largest risk, and
the conformance oracle that pins correctness. This is direction, not a contract —
update as milestones land or oracles change.

## Thesis

DuckHTS should own the NGS pipeline end-to-end in SQL — `fastq → align → coverage/QC
→ variants → consequence` — where every stage is a queryable relation, memory is
bounded, and cores saturate. Four workstreams get there; three share one enabling
primitive (a typed CIGAR/aligned-pair walker), so it comes first.

Two guardrails, both from `AGENTS.md` discipline:
- **No headline parity-% while backlog remains.** Coverage is a counted, per-class
  number (see the accepted-mismatch ledger, M1), never a marketing figure.
- **Never use a fused collector as the oracle for its own components.**

## Upstream trees and pinned oracles

| Role | Location | Pin |
|---|---|---|
| Picard-compatible QC (Rust) | `.sync/riker` | `09ebb21` (v0.4.0) |
| RSeQC/Qualimap/featureCounts/dupRadar (Rust) | `.sync/RustQC` | `fee1fce` |
| mosdepth | `.sync/mosdepth` | `52813e0` (v0.3.13) |
| bcftools (BCSQ contract) | system / `.sync/pysam/bcftools` | 1.23 |
| Ensembl VEP | `/root/ensembl-vep` | `release/115` (SO + HGVS oracle); duckvep-c kernel grounded in 116 — see reconciliation doc for the pin decision |
| fastVEP (contrast, NOT a source) | `/root/fastVEP`, github.com/Huang-lab/fastVEP | `785922e` (v0.2.0) |
| duckvep-c (consequence engine to fold in) | `/root/duckvep-c` | `9f922c8` |

## Primitive verdicts (reconciled)

- **CIGAR interval explosion → YES, as a typed C walker, not a table function.**
  The primitive is an allocation-free walker yielding op/length, query+ref
  start/end, consume flags, aligned blocks, clips, indels, and `N` splice gaps —
  collectors call it **directly** on the hot path. A `read → ref-interval` table
  function is an *inspection view* only. Ref-only intervals are insufficient for
  error/RNA/clipping metrics. Existing binary CIGAR in `src/bam_reader.c`.
- **Per-cycle SIMD aggregates → YES, scalar oracle first.** base × quality
  histograms stratified by cycle, read number, orientation. SIMD only the
  byte-oriented counting/reduction. Existing kernels: `src/include/duckhts_simd_kernels.def`
  (base/nt16/GC only today).
- **MD/aligned-pair walking → YES, MD optional.** One aligned-pair iterator with
  MD-driven *and* FASTA-driven modes; emit match/mismatch/ins/del, coords, base
  quality, reference base. Define missing/inconsistent-MD behavior.
- **`FILE_OFFSET` deduplication → NO.** `FILE_OFFSET` is an order/identity token,
  not a dup key. Real duplicate signatures need unclipped 5′ coords + orientation
  + mate coords + library/RG policy — a metric-level policy, not a primitive.
  `FILE_OFFSET` is exposed/tested in `src/bam_reader.c` / `test/sql/duckhts.test`;
  the WisecondorX adjacent-suppression path is narrower (`src/bam_bin_counts.c`).

## `FILE_OFFSET` solves the parallel-produce / ordered-consume problem

Do **not** hand-roll a C reorder buffer for QC row order. Workers produce
out-of-order; every row carries `FILE_OFFSET` (the BGZF virtual offset from
`read_bam`); the consumer recovers input order with `ORDER BY FILE_OFFSET`
downstream, and DuckDB's **spilling sort** makes it memory-safe. Cost: storing a
(possibly temporary) offset column. Leverage the engine, don't rebuild it.

Caveat: mosdepth's **file-output** order is a separate exact contract that still
needs per-contig header-order spooling. QC *row* order is now just a sort key.

## Milestones

### M1 — Freeze contracts, oracle pins, parity matrix; minibwa gate
The conformance-first unblock. Publish: metric-by-metric parity matrix; filtering
and dedup policy contracts; the pinned oracle versions above; fixture manifests;
and an **accepted-mismatch ledger** (see `design/duckvep_reconciliation.md` for the
ERRATA + manifest discipline: manifest = *which* oracle, errata = *why* we differ
per-oracle in three reason-classes, each entry backed by a regression witness;
provenance gaps must not hide as "intentional divergence"). Also: time-box a
minibwa native/Windows/WASM portability spike → **promote/reject**. minibwa is
gated/optional, never a dependency.
Risk: "parity" is many incompatible policy contracts, not one algorithm.
Oracle: the pins table; per-metric exact or documented-tolerance.

### M2 — Shared typed CIGAR/aligned-pair walker
The real primitive layer (see verdicts). Per-cycle scalar histograms → SIMD
backends via `duckhts_simd_kernels.def`. Exit: SQL/R tests, forced `scalar` == `auto`.
Feeds QC-RNA and the duckvep annotation joins.
Risk: field-requirement discipline (decode once, only what's needed).
Oracle: bedtools / RSeQC on the same BAM.

### M3 — Ordered streaming infrastructure
mosdepth tiled/event streaming (per-worker coverage scratch ~few MB, not
`int32[chrom_len]` — see `design/coverage_memory_footprint.md`); `FILE_OFFSET`-sort
order recovery (above); streaming `bam_bin_counts`; per-contig header-order spool
for the mosdepth file contract. Uses weighted contig claiming from
`design/better_scans.md`.
Risk: premature tile emission breaking mate-overlap correction / region medians /
exact output order.
Oracle: byte-identical to `duckhts_mosdepth` + mosdepth 0.3.13
(`scripts/mosdepth_conformance.py`, `test/sql/mosdepth.test`); thread-invariance;
≤ ~4 MB coverage scratch/worker.

### M4 — riker parity, then `duckhts_riker` fused
basic/alignment/isize → error/GC/WGS/hybcap, each an independently testable typed
function. Release `duckhts_riker` (one BAM decode, union of collector field
requirements, fan immutable batches to collector-owned states, mergeable collectors
by chunk + one ordered consumer per ordered collector) **only when fused output ==
separate runs**.
Risk: filtering/fragment/dedup/orientation policy divergence.
Oracle: riker `09ebb21` (+ the upstream tool each riker command rewrites);
never self-oracle.

### M5 — RNA + RustQC analytical parity
featureCounts / dupRadar / preseq / RSeQC / Qualimap. RNA is *not* reducible to
CIGAR explosion: transcript-space coverage, strand inference, gene assignment,
junction policy, TIN, transcript-space insert size are metric-level algorithms.
See `.sync/riker/src/commands/rna.rs`, `.sync/RustQC/src/rna/`. Keep plotting/report
rendering out of native kernels.
Oracle: RustQC `fee1fce` + the actual upstream tools (`.sync/RustQC/src/citations.rs`).

### M6 — Fuse duckvep-c consequence/HGVS engine
See `design/duckvep_reconciliation.md` (the detailed salvage + Model A import +
three-tier conformance + format-as-relations plan). Lift into `src/duckvep/`; the
three-tier conformance framework is the deliverable's core (the anti-slop
differentiator vs fastVEP). HGVS g→c→p after SO agreement.
Risk: transcript-model preparation (Model A import) — the pivot.
Oracle: VEP `release/115` (SO/HGVS) + ClinVar `CLNHGVS`; bcftools 1.23 (BCSQ subset);
brute-force genetic-code property tests for the kernel.

### M7 — minibwa `seq_align` (GATED / optional)
Only if M1's portability spike promotes; never a dependency, never delays M2–M6.
Vendor lh3 minibwa; fold ksw2 banded Smith-Waterman as a new SIMD kernel family
(scalar + SSE/NEON, wasm/CRAN-safe scalar fallback). Ship `seq_align(query, target)`
SW UDF first (no index); the FM-index read-mapper (`read_fastq → align(ref)`) is a
separate later arc. Credit lh3.
Oracle: bwa / minimap2 SAM + brute-force SW property tests.

## Dependency graph

```
             ┌─ CIGAR/pair walker (M2) ──┬─→ QC-RNA metrics (M4/M5)
             │                           └─→ duckvep annotation joins (M6)
M1 contracts ┼─ per-cycle SIMD aggs (M2) ─→ QC-DNA metrics (M4)
+ oracles    ├─ mosdepth streaming (M3) ──→ duckhts_riker fused pass (M4)
             ├─ Model A import (M6 pivot) → duckvep fuzzer can reach rare classes
             └─ minibwa ksw2 spike (M1) ──→ [gated] seq_align (M7)
```
