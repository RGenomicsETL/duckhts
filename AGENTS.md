# DuckHTS Agent Guidelines
write C as a BSD kernel programmer rather than a Java programmer that failed upwards
write R as a r-lib programmer rather than a Python programmer that failed upwards
## Project Goal
DuckDB extension that reads HTS file formats (VCF/BCF, BAM/CRAM, FASTA/FASTQ, BED, GTF/GFF, tabix) via htslib, with coverage/interval analytics for CNV and QC workflows. The R package (`Rduckhts`) is on CRAN.

## Source Layout
- Extension sources: `src/` (C source), `src/include/` (headers)
- SQL tests: `test/sql/`
- R package: `r/Rduckhts/`
- Benchmark notebooks and rendered benchmark reports: `benchmarks/`
- Design notes and planning docs: `design/`
- Vendor sources: `third_party/`
- Upstream mirrors (read-only reference): `.sync/`
- `CLAUDE.md` is a symlink to `AGENTS.md`; edit `AGENTS.md` only.

## Agent Working Instructions
1. Read relevant source files before making changes.
2. Check existing tests to understand expected behavior.
3. When mirroring external tool behavior, consult `.sync/` mirrors before secondary sources.

## Changelog — Mandatory
Every user-visible change must update **both**:
- `NEWS.md` (extension, newest first)
- `r/Rduckhts/NEWS.md` (R package, newest first)

Never consider a feature complete until both are updated.

R package changelog scope is strict:
- `r/Rduckhts/NEWS.md` must only describe R package user-facing changes: wrapper behavior, bundled extension behavior as exposed through the package, package docs/examples, package tests, CRAN/build/install behavior.
- Do **not** add repo-only developer tooling, extension-only benchmarking scripts, upstream comparison notes, or other non-package workflow items to `r/Rduckhts/NEWS.md`.
- Put repo-level tooling or extension-only notes in `NEWS.md` instead.

## Documentation
- `functions.yaml` is the source of truth for public function documentation and the community-extension descriptor.
- After adding/removing/renaming functions: update `functions.yaml`, run `python3 scripts/render_function_catalog.py`, bootstrap the R package, verify generated files in `r/Rduckhts/inst/function_catalog/` and `community-extensions/extensions/duckhts/description.yml` are in sync.
- `community-extensions/` is a local sync copy only — do not commit it; copy `description.yml` manually to the community-extensions repo after regenerating.
- `r/Rduckhts/README.Rmd` is wired to the generated function catalog; do not duplicate function lists by hand.
- Keep new benchmark `.Rmd` / `.md` files under `benchmarks/`, not the repo root.

- Scan-planning and deferred record-offset design work is tracked in `design/better_scans.md`. If a task touches weighted contig claiming or a future BCF/VCF offset column, read that file before changing reader code.
- FASTQ reader throughput and remote-safe parser design notes are tracked in `design/fastq_throughput.md`. If a task touches `read_fastq(...)` performance, read that file before changing the reader path.
- Native mosdepth rewrite planning is tracked in `design/duckhts_mosdepth.md`. If a task touches `duckhts_mosdepth(...)`, mosdepth compatibility, or native coverage rewrite work, read that file before implementing reader or coverage code.
- BAM/BED regional coverage, future `fragment_mode := TRUE` work, and irregular-interval coverage design are tracked in `design/duckhts_idxstats_and_interval_coverage.md` and `design/coverage_memory_footprint.md`. If a task touches `duckhts_bam_bed_coverage(...)`, cgranges-backed irregular interval work, or BED regional coverage memory/semantics, read those files first. For `duckhts_bam_bed_coverage(..., fragment_mode := TRUE)`, do not invent new pair semantics from scratch: mirror the existing `duckhts_mosdepth(..., fragment_mode := TRUE)` / upstream `mosdepth --fragment-mode` policy as the primary semantic anchor, and treat `meanbaseq` / `meanmapq` as read-level summaries unless the design doc is explicitly updated.
- Samtools-compatible kernel planning is tracked in `design/duckhts_samtools.md`. If a task touches new `duckhts_samtools_*` functions, scan-fused summary kernels, or future samtools-compatible rewrites beyond `duckhts_samtools_idxstats(...)`, read that file before implementing.
- DuckDB C API deprecation tracking is documented in `design/duckdb_c_api_deprecation_scan_2026-04-21.md`. If a task touches DuckDB C API usage, bundled DuckDB header bumps, or extension/runtime integration code, read that file first and re-run the deprecation scan after updating DuckDB.

## SIMD Kernel Development Workflow
DuckHTS SIMD dispatch is a capability-mask, per-logical-kernel system under `src/simd/`.  Do not add ad-hoc ISA checks inside individual SQL functions when the logic can be represented as a reusable byte-oriented kernel.

Use this flow when adding a SIMD operation:

1. **Define the logical kernel first**
   - Add one entry to `src/include/duckhts_simd_kernels.def` using `DUCKHTS_SIMD_KERNEL(enum_id, dispatch_table_field, signature_tag, sql_name)`.
   - Keep kernel ids separate from capability bits.  Never use enum values as bit positions.
   - If the function signature is new, add the typedef and `DUCKHTS_SIMD_FN_<TAG>` mapping in `src/include/duckhts_simd_internal.h`.
   - Add a typed builder-consider helper in `duckhts_simd_dispatch.c`; avoid untyped function-pointer punning.

2. **Scalar is mandatory**
   - Implement the scalar reference in `src/simd/duckhts_simd_scalar.c` before any ISA backend.
   - Register scalar for the kernel.  Scalar is the correctness oracle and must resolve every logical kernel.
   - `duckhts_simd_require_table_resolved()` must still prove that every dispatch table slot has a fallback.

3. **Register optional backends independently**
   - Put ISA-specific implementations in the existing backend files (`duckhts_simd_avx2.c`, `duckhts_simd_avx512.c`, `duckhts_simd_neon.c`, `duckhts_simd_wasm_simd128.c`) or add a new backend translation unit only when needed.
   - Compile/runtime-gate with compiler/platform-native checks, not htslib autoconf `HAVE_*` macros.
   - x86 kernels that require an ISA should use function-level target attributes where possible so CPU probing remains safe on baseline hosts.
   - Backend registration functions should register only kernels they actually implement.  Public trampolines must tolerate missing non-scalar slots by falling back through the resolved table.
   - Lower priority numbers win under `auto`; choose priorities deliberately and document surprising ordering.

4. **Snapshot dispatch in consumers**
   - Consumers such as UDFs should snapshot once per DuckDB vector/chunk with `duckhts_simd_dispatch_snapshot()` and call `*_with_table(...)` helpers inside the row loop.
   - Table functions should capture the dispatch table in init data when a scan needs a stable view.
   - Do not reload the global dispatch pointer per row unless there is a strong reason.

5. **Wire both extension and R package builds**
   - If adding new source files, update `CMakeLists.txt`, `r/Rduckhts/R/bootstrap.R`, `r/Rduckhts/configure`, and `r/Rduckhts/configure.win`.
   - Run `cd r/Rduckhts && Rscript bootstrap.R /root/duckhts` after source/header changes and verify bundled copies under `r/Rduckhts/inst/duckhts_extension/` match the root sources.

6. **Validation requirements**
   - Add SQL conformance tests that compare forced `scalar` with `auto` for representative inputs and edge cases.
   - Add R tinytests for wrappers or R-visible behavior.
   - Tests must be backend-agnostic: never require AVX2/AVX512/NEON/wasm availability on the runner.
   - Assert diagnostics semantics: `available = compiled && cpu_supported`; `selectable` only means a selectable implementation path exists; `auto` is a request, not a concrete backend row.
   - Cover invalid input, lowercase/uppercase behavior, `N`/missing handling, and any NULL/empty semantics for the kernel.

7. **Benchmarking requirements**
   - Add benchmark drivers under `scripts/` and rendered reports under `benchmarks/`.
   - Benchmarks should run each backend request in a fresh process when measuring process-wide backend selection.
   - Always record `duckhts_simd_kernel_info()` so reports show the concrete kernel backend selected by `auto`.
   - Include both a kernel-isolating benchmark and, when relevant, a real-data end-to-end workload (for example BAM/FASTQ/FASTA input) so I/O and materialization effects are visible.
   - Skip unavailable concrete backends instead of failing portable benchmark runs.

8. **Documentation and changelog**
   - Public SQL/R surfaces require `functions.yaml`, regenerated catalogs, README examples when useful, and both `NEWS.md` files.
   - Internal-only SIMD kernels should still be documented in design notes or benchmark reports when they affect performance claims.
   - Keep `design/simd_dispatch_matrix.md` current when changing dispatch architecture, semantics, or validation strategy.

Required validation after SIMD changes:

```
make release -j2
./configure/venv/bin/python3 -m duckdb_sqllogictest --test-dir test/sql --file-path test/sql/duckhts.test --external-extension build/release/duckhts.duckdb_extension
cd r/Rduckhts && Rscript bootstrap.R /root/duckhts && THREADS=4 make test
git diff --check
```

## Build & CI Rules
- Build scripts must be deterministic and non-interactive.
- No network access in the extension build step; only in explicit vendoring scripts.
- No vcpkg — use CMake and traditional autotools (CRAN requirement).

## Wasm / webR Debugging
- Reproduce `Rduckhts` wasm failures in the webR container, not from host build directories: use `ghcr.io/r-wasm/webr:main` and the `rwasm::build()` path.
- Do not assume `r/Rduckhts/inst/duckhts_extension/build/` in the repo is a wasm build artifact. Host `R CMD build` / `.Rcheck` trees often contain native ELF binaries with the same filename.
- For wasm package debugging, inspect the built `Rduckhts_<Version>.tgz` produced by `rwasm::build()`. The real browser payload is the bundled `Rduckhts/duckhts_extension/build/duckhts.duckdb_extension` inside that `.tgz`.
- For browser-real testing, prefer the local harness in `scripts/start_webr_local_test.sh` plus `scripts/webr-local-test.html`. This path now exercises the actual browser webR worker, installs the locally built `.tgz` from a tiny staged binary repo, runs `rduckhts_load()`, runs a bundled FASTQ smoke test, and can run the installed `tinytest` suite in-browser.
- When checking a wasm extension, verify both:
- the symbol table, e.g. with `emnm`
- the actual wasm export section, because a symbol name may exist in the module but still not be exported for DuckDB's loader
- A missing DuckDB init symbol in webR/browser can come from the final extension link dropping Emscripten `LDFLAGS`. In `r/Rduckhts/configure`, the final `duckhts.duckdb_extension` link must preserve `${LDFLAGS}` so `SIDE_MODULE` and related wasm flags reach the extension itself.
- If DuckDB reports that `duckhts.duckdb_extension` does not contain `duckhts_init_c_api`, do not stop at `strings`/`emnm`. Parse the wasm export section and confirm that `duckhts_init_c_api` is actually exported.
- Keep wasm-specific behavior gated on the real Emscripten target detection (`--host=wasm32-unknown-emscripten`, `emcc`, etc.). Do not let wasm workarounds change non-wasm targets accidentally.
- The DuckDB metadata trailer appended by `append_extension_metadata.R` is separate from `configure`/`htslib` probing. It affects DuckDB's runtime loader selection, not whether the package build succeeds.
- Do not trust the outer host `pkg-config` curl probe on wasm. The webR Docker image exposes host Linux `curl/curl.h` headers via Linux paths, but those are not a usable wasm libcurl backend. The `configure` script skips curl detection entirely on wasm; do not re-enable host pkg-config for this path.
- `r-wasm/webr` ships `/opt/webr/wasm/lib/libcurl.a` (a real wasm-compiled libcurl), and the emcc link test passes against it. However, libcurl's `connect()` / socket calls from a SIDE_MODULE hit a webR Emscripten message-bus limitation (`resolved is not a function` at `prop` in `R.js`) when the first network connection is attempted. Keep libcurl disabled on wasm in the package build.
- Current wasm/browser status: local-file workflows and the installed `tinytest` suite pass end-to-end in a real browser webR worker. HTTP/HTTPS access on wasm now goes through the package-owned `src/wasm_http_hfile.c` backend, which uses synchronous browser XHR from the worker instead of libcurl sockets.
- Browser constraints still apply on wasm: same-origin URLs work, remote URLs require permissive CORS headers, and `ALL_PROXY` / ws-proxy do not affect the XHR backend. Keep upstream `htslib` `libcurl`, `S3`, and `GCS` disabled on wasm unless the runtime model changes.
- Treat webR and duckdb-wasm as different runtimes and different artifacts. The working webR/browser payload is the extension bundled inside the built `Rduckhts_<Version>.tgz`; the duckdb-wasm local harness loads `/duckdb-wasm/duckhts.duckdb_extension.wasm` directly. Do not assume that a symbol/ABI result from one runtime applies to the other.
- Keep the Emscripten socket/i64 compat shim in one canonical header at `src/include/wasm_socket_compat.h`. Bootstrap that header into `r/Rduckhts/inst/duckhts_extension/include/` for the package build; do not keep a second copy under `r/Rduckhts/tools/`. The duckdb-wasm-only `orig$...` remap for native i64 libc imports (`lseek`, `ftruncate`, `strtoll`, `strtoull`, `atoll`, `time`) must stay gated behind `DUCKHTS_WASM_DUCKDB_RUNTIME` so it does not affect the webR/package build.
- In the duckdb-wasm harness, `db.registerFileURL(...)` only maps a DuckDB virtual path to an HTTP URL. If `LOAD` reaches the wasm loader and then traps, the problem is in the wasm side-module ABI/runtime, not in file registration.
- Older duckdb-wasm runtimes can fail with misleading load traps even after the side module is otherwise correct. The local duckdb-wasm harness currently defaults to `@duckdb/duckdb-wasm` `1.31.0`; the older `1.29.0` line (DuckDB `v1.1.1`) produced loader failures while the newer runtime loaded the same extension and passed the HTTP smoke tests. When debugging, log `SELECT version()` in the page before `LOAD`.
- The local `python3 -m http.server` harness ignores HTTP `Range` and returns `200 OK` full-object responses. The `wasm_http_hfile.c` backend is expected to warn and fall back to full-object fetch + local slicing in this setup; those `[W::wasm_http_read] server ignored HTTP Range ...` messages are noisy but not fatal.
- For tabix over HTTP, an initial `.csi` probe and `404` can be expected before htslib falls back to `.tbi`; a `.csi` fetch failure is not by itself a regression if the `.tbi` path succeeds and the query returns rows.
- The browser-side `XrayWrapper`/`content-script.js` error, the WebAssembly `try` deprecation warning, and missing `duckdb-browser.mjs.map` source maps were all non-blocking in this debugging session. Do not mistake them for the root cause if the extension loads and the same-origin HTS smoke tests pass.

## R Package Workflow — Mandatory After Any Extension Source Change
```
cd ~/duckhts/r/Rduckhts/
Rscript bootstrap.R ~/duckhts/
THREADS=4 make test
```

Additional R package rules:
- **Never** run `R CMD INSTALL .` from `r/Rduckhts/` — it mutates `inst/duckhts_extension/htslib` in place. Always build a tarball (`R CMD build .`) and install the tarball.
- When new source files are added under `src/`, update: `CMakeLists.txt`, `r/Rduckhts/R/bootstrap.R`, `r/Rduckhts/configure`, `r/Rduckhts/configure.win`.
- Version scheme: `duckhtsVersion-x` (e.g., `1.1.5-0.0.1`).
- All R changes must maintain CRAN compatibility.
- A new public function is incomplete until its C source is wired through both the extension build and the R package build (Unix and Windows).

## Testing
Two levels required for every public feature:
1. **SQL conformance** in `test/sql/` — schema, semantics, region/index behavior
2. **R package tests** in `r/Rduckhts/inst/tinytest/` — wrapper signatures, argument validation, bundled-extdata access, end-to-end DBI workflows

Guidelines:
- One `.test` file per feature family; one `tinytest` file per wrapper family.
- VCF/BCF fixtures: generate with `test/scripts/vcfpp.R` via `vcfppR`, not by hand. Record in `vcfpp_manifest.tsv`. Copy paired fixtures to `test/data/` and `r/Rduckhts/inst/extdata/`.
- New fixtures needed by SQL tests: document in `test/scripts/prepare_test_data.sh`, then copy into the R package bundle.
- README examples must be deterministic, short, and use bundled `inst/extdata/`.

## Vendoring
- Pin exact versions/commits; checksum-verify all downloads.
- Support offline rebuilds once sources are fetched.
- Do not modify vendored code directly — use explicit patch files.

## CSQ/ANN Typing Strategy
- Builtin CSQ/BCSQ typing rules live in `src/vep_parser.c`; seed from `.sync/` bcftools mirror.
- Default for unresolved fields is `String`; prefer typed allowlist plus string fallback.
- Extending builtin rules requires updating both SQL tests and R wrapper tests.

## Phase 10 — Coverage & Interval Primitives
Layer 0 (SQL-first baseline via `read_bam` + GROUP BY) is available now. Remaining layers:
- **Layer 1**: `read_bed`, `fasta_nuc`, interval UDFs (`src/interval_udf.c`)
- **Layer 2**: native counting kernels — `bam_bin_counts`, `bam_bedcov`, `bam_coverage`, `bam_depth`, `bam_pileup`
- **Layer 3**: indexed BED export / tabix round-trip

Key design rules:
- One function per counting model; BED-compatible output as the interoperability contract.
- `bam_bin_counts` output: `chrom`, `start`, `end`, `bin_id`, `count_total`, `count_fwd`, `count_rev`.
- `count_model := 'read_start'` for WisecondorX/NIPTeR-style; pileup/depth are separate functions.
- Fixed-width bins use arithmetic (`pos // bin_width`); `cgranges` (at `third_party/cgranges/`) is for irregular interval joins only.
- All Phase 10 functions need R wrappers, `functions.yaml` entries, and benchmarks against the SQL-first baseline.
- Parallelism: follow the contig-claiming pattern from `bam_reader.c` (atomic counter, per-thread file handles, `hts_set_threads(fp, 2)`).

## WisecondorX / BAM Dedup Context
`FILE_OFFSET UBIGINT` in `read_bam()` exposes the BGZF virtual offset after each record (`bgzf_tell` — zero-cost macro). Use `ORDER BY FILE_OFFSET` in `LAG()`/`LAST_VALUE()` window functions to reproduce exact BAM file order for streaming dedup. The exact SQL replication of WisecondorX larp/larp2 is in `scripts/wisecondorx_convert_conformance.py`; see `RWisecondorX/AGENTS.md` for the full algorithm description.

## Style
- Keep changes minimal and focused; preserve existing APIs unless the task explicitly requires changes.
- Follow R package conventions for CRAN submission.
