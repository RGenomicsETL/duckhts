# DuckHTS design note index

Status: current index for agent/developer orientation. Keep this file in sync when design notes are added, retired, or converted from plans into implementation contracts.

Use this file to decide which design note to read. Individual notes may contain historical context; their status line tells whether they are current guidance or background.

| Design note | Status | Read when touching |
| --- | --- | --- |
| [`better_scans.md`](better_scans.md) | Open design/backlog. `scan_mode := 'sequential'` has landed; weighted contig claiming and BCF/VCF offset semantics remain open. | BAM/BCF scan planning, index-backed full scans, explicit sequential scans, future record-offset columns. |
| [`coverage_memory_footprint.md`](coverage_memory_footprint.md) | Current backlog. `duckhts_bam_bed_coverage` tiling landed; mosdepth tiling and `bam_bin_counts` streaming remain open; `bcf_reader.c` memory concern was closed as a false alarm. | Coverage-reader memory use, tiling, `bam_bin_counts`, `duckhts_bam_bed_coverage`, `duckhts_mosdepth`. |
| [`duckdb_c_api_deprecation_scan_2026-04-21.md`](duckdb_c_api_deprecation_scan_2026-04-21.md) | Historical scan. Re-run/update after DuckDB header/runtime bumps or broad C API refactors. | DuckDB C API upgrades, deprecation cleanup, extension load/bind API changes. |
| [`duckhts_mosdepth.md`](duckhts_mosdepth.md) | Compatibility contract plus backlog for the implemented native `duckhts_mosdepth(...)` rewrite. Some implementation-plan sections are historical. | Mosdepth parity, mosdepth outputs/options, coverage-engine changes that affect mosdepth behavior. |
| [`duckhts_parquet_lake.md`](duckhts_parquet_lake.md) | Implemented Parquet converter-helper surface plus future DuckLake/native-writer design. | `duckhts_*_convert_parquet_sql`, R Parquet converter wrappers, metadata keys, DuckLake registration patterns, future HTS writers. |
| [`duckvep_layer_keys.md`](duckvep_layer_keys.md) | Future DuckVEP/keying design. | VariantKey/RegionKey annotation joins, VEP-style consequence-layer planning, interval/exact annotation tables. |
| [`fastq_throughput.md`](fastq_throughput.md) | Open performance investigation. Direct htslib-transport FASTQ parser remains future work. | `read_fastq(...)` throughput, FASTQ parser refactors, htslib transport/remote-safe parsing. |
| [`multireading.md`](multireading.md) | Current guidance for generated SQL multi-file reads plus research notes for native multi-file readers. | `hts_union_query(...)`, multi-file reader ergonomics, `UNION ALL BY NAME`, future native glob/list readers. |
| [`simd_dispatch_matrix.md`](simd_dispatch_matrix.md) | Current SIMD dispatch contract. | SIMD backend selection, `duckhts_simd_*` diagnostics, logical kernel dispatch behavior. |
| [`simd_future_kernel_proposals.md`](simd_future_kernel_proposals.md) | Future proposals. | New SIMD kernels, parser/depth/normalization SIMD planning. |

Do not add `AGENTS.md` references to nonexistent design files. Prefer updating this index and linking here from `AGENTS.md` rather than adding many direct stale-prone pointers.
