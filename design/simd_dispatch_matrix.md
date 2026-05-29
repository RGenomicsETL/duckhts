# DuckHTS SIMD dispatch matrix

This note tracks the long-term SIMD dispatch infrastructure introduced for byte-oriented DuckHTS kernels.

## Architecture

DuckHTS resolves logical kernels independently from concrete CPU/runtime capabilities:

1. `duckhts_simd_kernels.def` is the source of truth for logical SIMD kernels.
2. `duckhts_kernel_kind_t` is generated from that manifest and may grow to dozens or hundreds of logical kernels.
3. Concrete platform features are represented by a `duckhts_simd_cap_t` bitmask (`scalar`, `avx2`, `avx512`, `neon`, `wasm_simd128`, etc.).
4. Each backend translation unit registers only the kernels it implements with a typed builder helper.
5. Initialization builds immutable dispatch tables for `auto`, `scalar`, and concrete selectable backend policies. Runtime selection swaps an atomic pointer to one immutable table.
6. `auto` uses every currently available selectable capability and chooses the best implementation per logical kernel by priority. Concrete backend requests use that backend where implemented and scalar fallback for missing kernel slots.

The initial manifest contains only `seq_base_counts`, used by `seq_gc_content(...)`, but the shape is intended for future quality scans, delimiter scans, BAM packed-sequence counts, depth reductions, and allele normalization helpers. Unknown architectures, musl-based systems, and RISC-V currently remain scalar-only unless a backend is explicitly added with compiler and runtime probes; this avoids adding libc-specific detection paths before DuckHTS has a tested RVV kernel.

## Conformance test plan

Conformance is backend-agnostic:

- SQL tests assert backend inventory invariants: scalar is always available, `available == compiled && cpu_supported`, `auto` is not a concrete backend row, and known placeholder capabilities such as `sse2` are recognized but not selectable.
- SQL tests assert kernel-dispatch invariants through `duckhts_simd_kernel_info()`: every logical kernel has exactly one selected backend in the current immutable table, the `seq_base_counts` row exists, capability names match backend names for current kernels, and the dispatch mode is `capability-mask-dispatch-table`.
- SQL tests force `scalar`, compute `seq_gc_content(...)` over long deterministic sequences, then force `auto` and require bit-identical decimal results.
- R tinytests repeat the SQL diagnostics through thin DBI wrappers, check R-side scalar/non-missing argument validation, and verify `rduckhts_simd_kernel_info()` exposes the same kernel-level diagnostics.
- Platform-specific assertions must never assume AVX2, AVX512, NEON, wasm SIMD128, or future RVV availability. They may only assert implications from `duckhts_simd_info()` and `duckhts_simd_kernel_info()`.

## Benchmark setup

The benchmark path remains `make bench-simd-seq-gc`:

- `scripts/benchmark_simd_seq_gc.py` launches a fresh process for each backend request.
- The script records both the dispatch label (`duckhts_simd_backend()`) and the concrete `seq_base_counts` kernel backend from `duckhts_simd_kernel_info()`.
- Unavailable concrete backend requests are skipped instead of failing so one report can run across x86, ARM, wasm, and scalar-only builds.
- `benchmarks/benchmark_simd_seq_gc.Rmd` renders the configuration, per-request timings, scalar-relative speedups, and raw timing vectors.

Future benchmarks should follow the same pattern: benchmark logical kernels, record kernel-level dispatch diagnostics, and keep SQL-visible backend selection out of the timed loop.
