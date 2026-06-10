# DuckHTS north star: architecture spine and roadmap

Status: current north-star architecture and roadmap. This is the apex design note. It defines the layered spine the project builds on, the through-line that connects the per-thread notes, and the ordered roadmap. Individual notes listed in `design/README.md` are pillar documents that hang off one layer of this spine; when they disagree with `functions.yaml` or current source, source wins and the note must be corrected.

## Thesis

DuckHTS exists to prove a specific claim: **C is not dead, htslib is not a dead end, and the right way to make genomics fast is to put domain-aware bytes through a universal relational engine instead of inventing private formats.** DuckDB (the engine DuckDB Labs maintains) gives us a vectorized execution model, a portable build, Parquet/DuckLake as a durable columnar substrate, and SQL as a universal, composable surface. htslib gives us correct, battle-tested transport and codecs. Our job is the thin, high-value layer between them: byte-oriented kernels, packed domain encodings, projection-aware readers, and fused analytics — exposed as ordinary SQL functions and a high-level R API.

We do not invent make-believe formats. We do not absorb other tools' `main()` loops. We rewrite the *operations that fuse with the engine* and we lean on the engine and htslib for everything else. The payoff for users is libraries they can build on and strive with — DuckVEP is the flagship proof that this layer can carry a real, hard genomics workload end to end.

## The spine (bottom-up)

The whole project is one stack. Each layer has a single responsibility and a pillar note.

### L0 — Storage fundamentals and bytes

The foundation is *bytes in a domain space*, and the recurring primitive is **a packed small alphabet with a missing/masked state**:

- nucleotides: `{A,C,G,T}` + `N`/IUPAC ambiguity + soft-mask (case) — 2-bit core, wider with classes.
- BAM nt16: `A=1,C=2,G=4,T=8,…,N=15` — one nibble per base, IUPAC already encoded.
- genotypes: `{hom-ref, het, hom-alt, missing}` — 2 bits per call (PLINK `.bed` layout).
- quality scores: small integer range, delta/bin friendly.

These are the same shape, and they reduce by the same operation: **masked popcount over bitplanes**. That observation is load-bearing for L1 and L2.

Durable storage is ordinary, portable columnar data, never a private cache:

- Parquet owns columnar bytes + file-level key-value metadata (`duckhts_parquet_lake.md`).
- DuckLake owns catalog, snapshots/time-travel, schema evolution, file statistics, registration of premade Parquet.
- htslib owns HTS transport, BGZF, indexes, and codecs; BGZF virtual offsets (`FILE_OFFSET`) are first-class for streaming order and dedup.
- DuckHTS owns parsing, header fidelity, source-format metadata, packed domain encodings, and (future) HTS writers.

Pillars: `duckhts_parquet_lake.md`, `fastq_throughput.md` (transport), `read_bcf_v2_optimization_comparison_2026-06-09.md` (avoid-work-first).

### L1 — Portable SIMD kernel layer

A capability-mask, per-logical-kernel dispatch framework (`src/simd/`, `simd_dispatch_matrix.md`). The manifest (`duckhts_simd_kernels.def`) is the single source of truth; scalar is the correctness oracle; `auto` resolves the best implementation per kernel from the runtime capability mask; callers snapshot one immutable table per chunk.

The kernel idiom this project standardizes on is **LUT-classify + masked popcount**, portable across every backend we ship:

- table lookup: `vpshufb` (AVX2/AVX-512), `vqtbl1q_u8` (NEON), `i8x16.swizzle` (wasm128).
- classify a byte/nibble into a class bitfield (`isGC/isAT/isN/isIUPAC/isMasked/other`), then `movemask` + `popcount` per class.

This idiom is non-destructive (no `& 0xDF` case-fold that erases soft-mask), total over the alphabet (never abort a sequence on one ambiguous base), and the same code shape for text bases, BAM nt16, and genotype reductions. The kernel is a byte oracle that returns per-class counts; **policy lives in the SQL layer above** (does soft-mask count as called? do `S`/`W` contribute to GC? is a stray IUPAC code an error?).

High-value kernels, in priority order:

1. `bam_nt16_counts` — feeds `bam_bin_counts`, markduplicate, QC.
2. quality-score histogram/bin reduction — feeds FASTQ/BAM QC.
3. `base_counts` rewrite — LUT classifier replacing the 5-compare/`0xDF` shortcut; resolves N / soft-mask / IUPAC explicitly.
4. CSQ/INFO delimiter-and-field scan — feeds DuckVEP and the BCF reader hot path.
5. genotype-matrix reductions — AC/AN/missingness as masked popcount over packed rows.

Pillars: `simd_dispatch_matrix.md`, `simd_future_kernel_proposals.md`, `review_simd_dispatch_2026-05-29.md` (historical).

SIMD discipline: it is a hot-loop reduction and parsing/encoding accelerator at the I/O rind. It is **never** the core-logic play — branchy consequence/transcript logic stays scalar.

### L2 — SQL-native primitives and readers

Projection-aware, header-faithful readers and fused analytics, all ordinary table/scalar functions (`functions.yaml` is the catalog of record).

- variant/sites: `read_bcf` (general typed INFO/FORMAT + CSQ/ANN/BCSQ). Population-scale fixes are field/sample pushdown + bind-time projection (`read_bcf_v2_optimization_comparison_2026-06-09.md`). **The dense genotype matrix is a separate reader, not `read_bcf`.**
- genotype matrix (new): variant-major, 2-bit packed `read_bcf_gt`-style surface for matrix ops (AF, missingness, GRM, PCA, LD). Decodes to a contiguous numeric matrix on demand; sample axis is metadata, not columns. On-ramp to L4 SVCR.
- coverage/interval: `read_pileup`, `bam_bin_counts`, `duckhts_bam_bed_coverage`, `duckhts_mosdepth`, `duckhts_samtools_idxstats`, cgranges overlap. One counting model per surface (`coverage_memory_footprint.md`, `duckhts_mosdepth.md`).
- new readers (planned): bigWig / bigBed / bedGraph / PAF. BED-compatible output is the interoperability contract.
- keys: VariantKey (exact allele identity), RegionKey + cgranges (span/overlap) (`duckvep_layer_keys.md`).
- fused bcftools ops: `bcftools_norm_row`, `bcftools_score`, `bcftools_munge_row`, the `bcftools_filter` expression engine. We rewrite operations that fuse; full bcftools CLI fidelity is RBCFTools' job, not ours.
- scan planning: `scan_mode := 'sequential'`, indexed full scans, multi-file via `hts_union_query` (`better_scans.md`, `multireading.md`).
- alignment kernels (future): banded SW / edit-distance / PAF-oriented primitives as SIMD kernels — heavier, gated behind a real consumer.

### L3 — Annotation and consequence engine (DuckVEP) — the flagship proof

A DuckDB-native VEP / `bcftools csq` layer built entirely from L0–L2 (relations, Parquet, VariantKey/RegionKey, cgranges, `read_gff`, `fasta_nuc`, typed CSQ parsing). **No private cache format** (`.fastSA`/Sereal are explicitly not adopted). Storage is solved; the consequence engine is the entire risk. Target `bcftools csq` first (bounded, validatable against an executable), Ensembl VEP `release/115` as semantic comparator, HGVS as a later layer.

Pillars: `duckvep_bcftools_csq_port_plan_2026-06-09.md`, `duckvep_layer_keys.md`.

### L4 — Lakehouse / VariantLake

DuckLake as the catalog/derived/control plane; DuckHTS as the genotype/annotation producer. The differentiated bet is a **local-allele sparse long table (SVCR-style: `LGT/LAD/LPL` + `LA`, reference-block collapsed)** that is N+1 append-native and time-travels via DuckLake snapshots — the only one of the major systems (TileDB-VCF, GenomicsDB, Hail VDS, plain lakehouse) with no good open SQL-native implementation. Dense packed hardcalls (L2) serve common-variant matrix ops; sparse local-allele rows serve cohort-merge / rare variants. Frozen releases export to VCF-Zarr / Parquet / BGEN. **Do not** clone TileDB anchor-cell arrays; **do not** fork DuckLake onto FastLanes — track FastLanes encodings via the converter layer and contribute upstream instead.

Pillar: `duckhts_parquet_lake.md` (+ a future `variant_lake_svcr.md` when the design is written).

### L5 — High-level R API and testing

- A tidy, verb-based R layer (tidyomics-style) above the DBI wrappers in `r/Rduckhts/R/`, CRAN-compatible, composing the L2/L3 surfaces into ergonomic pipelines.
- Testing is two-level and Rmd-driven: SQL conformance in `test/sql/` + tinytest wrappers, plus `.Rmd`-driven benchmark/conformance reports under `benchmarks/` that **run locally first** (we do not pay for GitHub minutes) and degrade gracefully on GitHub CI (skip unavailable backends/data rather than fail). Benchmark in fresh processes; record `duckhts_simd_kernel_info()`.

## Cross-cutting contracts

- **Portability is non-negotiable**: every kernel and reader must survive native x86/ARM, wasm (`-fwasm-exceptions`, init-symbol exports), musl, and scalar-only/RISC-V. Unknown ISAs stay scalar until a tested backend with real probes lands.
- **Compatibility-rewrite discipline** (`rewrites.bio`): pin upstream version/commit, define the supported subset, validate continuously against the executable, document gaps, credit authors. Applies to mosdepth/bcftools/samtools/VEP/WisecondorX ports.
- **TigerStyle**: small C helpers, explicit `cleanup:` ownership, checked size arithmetic, no runtime `abort()`/`exit()` — report through DuckDB errors (`tigerstyle_audit_2026-06-08.md`).
- **Source-of-truth discipline**: `functions.yaml` for the public surface; `src/*.c` authoritative, `r/Rduckhts/inst/duckhts_extension/` generated; dual changelog (`NEWS.md` + `r/Rduckhts/NEWS.md`) mandatory.

## Roadmap (ordered bets)

The sequencing is chosen so each step makes the next easier to land.

1. **markduplicate + `bam_nt16_counts`** — shovel-ready (builds on `FILE_OFFSET` streaming-dedup); the QC consumer that justifies the SIMD framework and exercises the LUT-classify idiom. Lands L1+L2 together.
2. **`read_bcf` pushdown** — sample/field/bind-time projection from the v2 comparison; high leverage for most workloads, prerequisite decode plumbing for the genotype reader.
3. **Packed-state design note** — unify base-class / nt16 / genotype under one kernel family and the `base_counts` N/mask/IUPAC rewrite; review before touching the manifest.
4. **Genotype matrix reader** (`read_bcf_gt`) — variant-major 2-bit packed surface for matrix ops; on-ramp to L4.
5. **DuckVEP phase 2** — narrow `bcftools csq` slice; the flagship proof. Its transcript/exon packs are also L4 annotation substrate, so it de-risks the lake.
6. **VariantLake / SVCR design note** — the big bet; write the design before any code; depends on 4 and 5.
7. **New readers** (bigWig/bigBed/bedGraph/PAF) and **tidyomics R layer** — breadth and ergonomics once the spine is proven.

Explicitly shelved as sprawl traps: absorbing all of bcftools via an in-process hFILE/`main()` shim (that is RBCFTools' product), and forking DuckLake onto FastLanes.