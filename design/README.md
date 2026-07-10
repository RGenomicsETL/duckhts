# DuckHTS design note index

Status: current index for developer orientation. Keep this file in sync when design notes are added, retired, or converted from plans into implementation contracts.

Use this file to decide which design note to read. Individual notes may contain historical context; their status line tells whether they are current guidance or background.

**Start here:** [`north_star.md`](north_star.md) is the apex architecture-and-roadmap note. It defines the layered spine (L0 storage/bytes → L1 portable SIMD → L2 readers/primitives → L3 DuckVEP → L4 VariantLake → L5 R/testing) and the ordered roadmap. Every note below is a pillar hanging off one layer; the **Layer** column says which.

| Design note | Layer | Status | Read when touching |
| --- | --- | --- | --- |
| [`north_star.md`](north_star.md) | apex | Current north-star architecture + roadmap; the through-line for every note below. | Project direction, where a new feature belongs, roadmap sequencing, what is explicitly shelved. |
| [`roadmap.md`](roadmap.md) | apex | Future proposal. Ordered implementation sequence (QC/coverage foundations → riker/RustQC parity → DuckVEP fuse → gated minibwa) with oracle pins, the `FILE_OFFSET` sort-order model, and primitive verdicts. Elaborates `north_star.md` sequencing. | Implementation order, upstream oracle pins, QC-stack primitives, coverage memory/order model, minibwa gate. |
| [`duckvep_model.md`](duckvep_model.md) | L3 | Current implementation contract. Binding normalized annotation relation pack, logical/build receipts, edit and sequence-state rules, immutable `allele_transcript_context`, ownership, and validation gates. | Implementing the annotation importer/context builder or checking whether a model is reproducible. |
| [`duckhts_parquet_lake.md`](duckhts_parquet_lake.md) | L0/L4 | Implemented Parquet converter-helper surface plus future DuckLake/native-writer design. | `duckhts_*_convert_parquet_sql`, R Parquet converter wrappers, metadata keys, DuckLake registration patterns, future HTS writers. |
| [`fastq_throughput.md`](fastq_throughput.md) | L0/L1 | Open performance investigation. Direct htslib-transport FASTQ parser remains future work. | `read_fastq(...)` throughput, FASTQ parser refactors, htslib transport/remote-safe parsing. |
| [`packed_state_kernels.md`](packed_state_kernels.md) | L0/L1 | Open implementation contract; review-before-code gate for the SIMD manifest. Unifies base-class / BAM nt16 / genotype reductions under LUT-classify + masked-popcount; `base_counts` N/mask/IUPAC rewrite. | Adding SIMD kernels, `base_counts`/`seq_gc_content` semantics, `bam_nt16_counts`, genotype AC/AN/missing reductions, soft-mask/IUPAC handling. |
| [`simd_dispatch_matrix.md`](simd_dispatch_matrix.md) | L1 | Current SIMD dispatch contract. | SIMD backend selection, `duckhts_simd_*` diagnostics, logical kernel dispatch behavior. |
| [`simd_future_kernel_proposals.md`](simd_future_kernel_proposals.md) | L1 | Future proposals. | New SIMD kernels, parser/depth/normalization SIMD planning. |
| [`review_simd_dispatch_2026-05-29.md`](review_simd_dispatch_2026-05-29.md) | L1 | Historical code review of the SIMD dispatch landing. | Background on why the dispatch layer is shaped as it is; superseded by `simd_dispatch_matrix.md` for current contract. |
| [`better_scans.md`](better_scans.md) | L2 | Open design/backlog. `scan_mode := 'sequential'` has landed; weighted contig claiming and BCF/VCF offset semantics remain open. | BAM/BCF scan planning, index-backed full scans, explicit sequential scans, future record-offset columns. |
| [`read_bcf_v2_optimization_comparison_2026-06-09.md`](read_bcf_v2_optimization_comparison_2026-06-09.md) | L2 | Current optimization comparison/backlog. Captures rendered `read_bcf_v2` benchmark evidence, current `read_bcf` gaps, field/sample pushdown plan, SIMD guidance, and test plumbing. | `src/bcf_reader.c`, BCF/VCF reader performance, projection/sample pushdown, CSQ/VEP reader materialization costs, SIMD decisions for reader hot paths. |
| [`coverage_memory_footprint.md`](coverage_memory_footprint.md) | L2 | Current backlog. `duckhts_bam_bed_coverage` tiling landed; mosdepth tiling and `bam_bin_counts` streaming remain open; `bcf_reader.c` memory concern was closed as a false alarm. | Coverage-reader memory use, tiling, `bam_bin_counts`, `duckhts_bam_bed_coverage`, `duckhts_mosdepth`. |
| [`duckhts_mosdepth.md`](duckhts_mosdepth.md) | L2 | Compatibility contract plus backlog for the implemented native `duckhts_mosdepth(...)` rewrite. Some implementation-plan sections are historical. | Mosdepth parity, mosdepth outputs/options, coverage-engine changes that affect mosdepth behavior. |
| [`multireading.md`](multireading.md) | L2 | Current guidance for generated SQL multi-file reads plus research notes for native multi-file readers. | `hts_union_query(...)`, multi-file reader ergonomics, `UNION ALL BY NAME`, future native glob/list readers. |
| [`tigerstyle_audit_2026-06-08.md`](tigerstyle_audit_2026-06-08.md) | cross-cut | Open audit note. Structural TigerStyle conformance findings; remediation backlog. | Broad style conformance review and code quality refactor planning across DuckHTS + Rduckhts. |
| [`duckdb_c_api_deprecation_scan_2026-04-21.md`](duckdb_c_api_deprecation_scan_2026-04-21.md) | cross-cut | Historical scan. Re-run/update after DuckDB header/runtime bumps or broad C API refactors. | DuckDB C API upgrades, deprecation cleanup, extension load/bind API changes. |

### Planned pillar notes (not yet written)

Named in `north_star.md`; create when the work starts, then add a row above:

- `genotype_matrix_reader.md` (L2) — variant-major 2-bit packed `read_bcf_gt` surface for matrix ops.
- `variant_lake_svcr.md` (L4) — local-allele sparse long-table cohort store on DuckLake.
- `new_readers_bigwig_bigbed_bedgraph_paf.md` (L2) and `alignment_kernels.md` (L1/L2).
- `r_highlevel_tidyomics.md` (L5) and `testing_rmd_suite.md` (cross-cut).

Do not add `AGENTS.md` references to nonexistent design files. Prefer updating this index and linking here from `AGENTS.md` rather than adding many direct stale-prone pointers.
