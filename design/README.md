# DuckHTS design note index

Status: current index. Keep notes that describe a live contract or a concrete open
investigation; delete completed planning narratives when code, tests, and history have
replaced them.

| Design note | Status and use |
| --- | --- |
| [`duckvep_model.md`](duckvep_model.md) | Current DuckVEP model, sorted C sweep, supplementary-track, ownership, and validation contract. |
| [`duckhts_parquet_lake.md`](duckhts_parquet_lake.md) | Implemented Parquet converter helpers and open DuckLake/native-writer work. |
| [`fastq_throughput.md`](fastq_throughput.md) | Open FASTQ throughput investigation. |
| [`packed_state_kernels.md`](packed_state_kernels.md) | Open contract for additional packed-state SIMD kernels. |
| [`simd_dispatch_matrix.md`](simd_dispatch_matrix.md) | Current SIMD dispatch contract. |
| [`simd_future_kernel_proposals.md`](simd_future_kernel_proposals.md) | Candidate SIMD kernels requiring measurement before implementation. |
| [`review_simd_dispatch_2026-05-29.md`](review_simd_dispatch_2026-05-29.md) | Historical review; use the dispatch matrix for current behavior. |
| [`better_scans.md`](better_scans.md) | Open scan-planning work. |
| [`read_bcf_v2_optimization_comparison_2026-06-09.md`](read_bcf_v2_optimization_comparison_2026-06-09.md) | Current BCF projection and sample-pushdown comparison. |
| [`coverage_memory_footprint.md`](coverage_memory_footprint.md) | Current coverage-reader memory backlog. |
| [`duckhts_mosdepth.md`](duckhts_mosdepth.md) | Mosdepth compatibility contract and remaining backlog. |
| [`multireading.md`](multireading.md) | Current generated-SQL multi-file guidance and native-reader research. |
| [`tigerstyle_audit_2026-06-08.md`](tigerstyle_audit_2026-06-08.md) | Open cross-project code-quality audit. |
| [`duckdb_c_api_deprecation_scan_2026-04-21.md`](duckdb_c_api_deprecation_scan_2026-04-21.md) | Historical DuckDB C API scan; rerun after dependency bumps. |
