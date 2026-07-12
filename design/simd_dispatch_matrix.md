# DuckHTS SIMD dispatch

Status: current implementation contract for byte-oriented SIMD dispatch. Update this file when
backend-selection semantics, diagnostics, or backend coverage change.

## Architecture

1. `src/include/duckhts_simd_kernels.def` is the source of truth for logical kernels.
2. Concrete platform features are represented by a capability mask.
3. Each backend translation unit registers only the kernels it implements through typed
   builder helpers.
4. Initialization builds immutable tables for `auto`, `scalar`, and selectable concrete
   policies. Runtime selection atomically swaps a pointer to one table.
5. `auto` chooses the highest-priority available implementation independently for each logical
   kernel.
6. A concrete backend request uses that backend where implemented and the scalar oracle for a
   missing slot.

Callers snapshot the table once per DuckDB vector, chunk, or worker phase. They do not perform
private ISA checks or load dispatch state once per row.

## Backend coverage

A dash means that an explicit request for that backend uses scalar for the missing slot. Under
`auto`, another available concrete backend may win instead.

| kernel | scalar | avx2 | avx512 | neon | wasm_simd128 |
| --- | :---: | :---: | :---: | :---: | :---: |
| `seq_base_counts` | yes | yes | yes | yes | yes |
| `bam_nt16_counts` | yes | yes | — | yes | yes |
| `nt16_gc_counts` | yes | yes | — | yes | yes |

AVX-512 currently accelerates only `seq_base_counts`; `auto` may select AVX2 for the nt16 slots
on a capable x86 host. The DuckDB-Wasm target opts into SIMD128 and requires all three slots;
webR is a separate toolchain and remains scalar unless its package build enables SIMD128.
Unknown or untested architectures remain scalar until they have compiler gates, runtime
probes, implementations, and conformance coverage.

## Conformance contract

- Scalar is always compiled and is the per-kernel correctness oracle.
- `available` means both compiled and CPU-supported; `auto` is a policy, not a concrete backend.
- Every logical kernel resolves to exactly one selected implementation in an immutable table.
- Portable SQL and R tests compare forced `scalar` with `auto` and validate diagnostic
  relationships without assuming a host ISA.
- Target-specific CI forces the backend promised by that build and fails if compilation,
  runtime support, or required registrations are absent.
- Backend tests cover deterministic boundaries, tails, packing alignment, invalid values, and
  randomized cases against independent slow oracles.
- Browser tests exercise real wasm output; a successful host build or exported symbol name is
  not evidence that SIMD128 executes in the browser.

## Benchmark contract

Use `make bench-simd-kernels` for isolated kernels and the relevant report under `benchmarks/`
for end-to-end materialization. Run process-wide backend selections in fresh processes, record
`duckhts_simd_kernel_info()`, and report scalar-relative results per logical kernel. A
microbenchmark multiplier is not an end-to-end reader or SQL speedup.
