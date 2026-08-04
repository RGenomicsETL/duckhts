# DuckHTS Agent Guidelines


## Prologue

Before landing a change, ask which invariant, interface, test, or shared helper would make
this change and the next one easier to own. Avoid ambiguity, duplicate authorities, code
sprawl, and planning artifacts that outlive their decision.

Read `ARCHITECTURE.md` for the project thesis and `STYLE.md` for the repository-specific
C, R, and SQL discipline.

## Project Goal

DuckHTS is a library-first C extension that makes sequencing data and genomics algorithms
composable in DuckDB. htslib supplies HTS transport and format semantics; reusable native
kernels supply coverage, interval, sequence, normalization, and consequence mechanics;
SQL supplies planning, joins, supplementary annotation, provenance, and explainable
clinical evidence. The same extension is packaged on CRAN as `Rduckhts`.

## Source Layout
- Extension sources: `src/` (C source), `src/include/` (headers), SIMD backends in `src/simd/`.
- SQL tests: `test/sql/`.
- R package: `r/Rduckhts/`.
- Benchmark drivers and rendered reports: `scripts/`, `benchmarks/`.
- Design notes: `design/`; start with `design/README.md` to find current contracts and open investigations.
- Binding architecture and style: `ARCHITECTURE.md` and `STYLE.md`.
- Vendored sources: `third_party/`.
- Upstream mirrors / pinned rewrite references: `.sync/` when present.
- `CLAUDE.md` is a symlink to `AGENTS.md`; edit `AGENTS.md` only.

## Agent Working Instructions
1. Read `ARCHITECTURE.md`, `STYLE.md`, and the relevant source files before making changes.
2. Check existing SQL and R tests to understand expected behavior.
3. When mirroring external tool behavior, consult `.sync/` mirrors before secondary sources.
4. Keep changes focused. Do not create branches, commits, large generated diffs, or workflow sprawl unless explicitly requested.
5. Preserve existing public APIs unless the task explicitly asks for an API change.
6. When referring to public GitHub issues or PRs, use full GitHub URLs.

## Public Surface

`functions.yaml` is the only hand-maintained public SQL catalog. Do not duplicate its
function inventory in agent instructions or design notes.

## Changelog — Mandatory
Every user-visible change must update **both**:
- `NEWS.md` (extension/repo, newest first)
- `r/Rduckhts/NEWS.md` (R package, newest first)

Never consider a feature complete until both are updated.

R package changelog scope is strict:
- `r/Rduckhts/NEWS.md` must only describe R package user-facing changes: wrapper behavior, bundled extension behavior as exposed through the package, package docs/examples, package tests, CRAN/build/install behavior.
- Do **not** add repo-only developer tooling, extension-only benchmarking scripts, upstream comparison notes, or other non-package workflow items to `r/Rduckhts/NEWS.md`.
- Put repo-level tooling or extension-only notes in `NEWS.md` instead.

## Documentation and Catalog Workflow
- `functions.yaml` is the source of truth for public function documentation and the community-extension descriptor.
- After adding/removing/renaming public functions or changing signatures/descriptions: update `functions.yaml`, run `python3 scripts/render_function_catalog.py`, bootstrap the R package, and verify generated files under `r/Rduckhts/inst/function_catalog/` and `community-extensions/extensions/duckhts/description.yml`.
- `community-extensions/` is a local sync copy only — do not commit it here; copy `description.yml` manually to the community-extensions repo after regenerating.
- `r/Rduckhts/README.Rmd` consumes generated catalog output; do not hand-maintain duplicated function lists.
- Render each `README.md` only from its corresponding `README.Rmd`; never hand-edit or cosmetically normalize rendered output.
- Keep benchmark `.Rmd` / `.md` files under `benchmarks/`, not the repo root.
- Design docs must start with a status line that says whether they are current implementation guidance, open design, future proposal, or historical review.

## Design Notes

Read `design/README.md`. Keep a note only while it defines a live contract or a concrete
unresolved choice. Code and tests replace completed implementation plans; git history and
GitHub issues retain the path and backlog.

## Build and CI Rules
- Build scripts must be deterministic and non-interactive.
- No network access in the extension build step; network is allowed only in explicit vendoring/staging scripts.
- No vcpkg for CRAN-oriented package builds; use the repo CMake/autotools paths.
- Pin exact upstream versions/commits and checksum-verify downloads.
- Do not modify vendored code directly; use explicit patch files.

### Pull Request Performance Evidence — Mandatory

Every pull request description must contain a `## Performance` section sourced from a
checked-in rendered benchmark report under `benchmarks/`. Name the report, source
revision, workload, input and output denominators, thread count, and measured result; do
not cite ad hoc terminal output. If the change touches a measured hot path, run and render
the comparable workload at the PR revision and compare it with the nearest identical
recorded workload. Do not claim "no regression" when workloads or environments differ.
If no checked-in benchmark exercises the changed path, state that explicitly in the
section, cite the nearest rendered baseline, and say what measurement is still missing.

### Pull Request Codex Review — Mandatory

Every pull request must receive a Codex review of its current head commit before merge.
After the ready-for-review PR is pushed, comment `@codex review` and wait for the Codex
response. Address every actionable finding, push the fix, and request another review;
repeat until Codex reports no major issue on the current head. Do not merge merely because
GitHub marks the branch mergeable or because CI passes. Record the reviewed commit in the
handoff when the review response names it.

Useful validation commands:

```bash
make release -j2
make test_release
cd r/Rduckhts && Rscript bootstrap.R /root/duckhts && THREADS=4 make test
git diff --check
```

If `duckdb_sqllogictest` is available through the configured venv, focused SQL tests can also be run with:

```bash
./configure/venv/bin/python3 scripts/run_sqllogictest.py \
  --test-dir test/sql \
  --file-path test/sql/duckhts.test \
  --external-extension build/release/duckhts.duckdb_extension
```

## R Package Workflow — Mandatory After Extension Source Changes
```bash
cd /root/duckhts/r/Rduckhts
Rscript bootstrap.R /root/duckhts
THREADS=4 make test
```

Additional R package rules:
- **Never** run `R CMD INSTALL .` from `r/Rduckhts/` — it mutates `inst/duckhts_extension/htslib` in place. Always build a tarball and install the tarball.
- When new source files are added under `src/`, update `CMakeLists.txt`, `r/Rduckhts/R/bootstrap.R`, `r/Rduckhts/configure`, and `r/Rduckhts/configure.win`.
- Version scheme: `duckhtsVersion-x` (for example `1.1.5-0.0.1`).
- All R changes must maintain CRAN compatibility.
- A new public function is incomplete until its C source is wired through both the extension build and the R package build on Unix and Windows.

### Version Bump Workflow
- The authoritative extension version is `version:` in the root `description.yml`.
- The authoritative R package version is `Version:` in `r/Rduckhts/DESCRIPTION`; keep the DuckHTS version before the hyphen and the R packaging revision after it.
- After releasing `X.Y.Z` to both CRAN and the DuckDB community extension repository, start the next development cycle by changing only those two declarations: `X.Y.Z` -> `X.Y.Z.9000` and `X.Y.Z-A` -> `X.Y.Z.9000-A`.
- For a release, perform the inverse transformation while retaining the current R packaging revision `A`.
- Do not hand-edit `configure/extension_version.txt`, rendered README/catalog files, or the local `community-extensions/` descriptor for a version-only bump. Normal configure, bootstrap, and catalog rendering carry the authoritative versions forward.
- A version-only development-cycle bump does not add empty `NEWS.md` sections; changelog entries belong to actual user-visible changes and release finalization.

## Testing and Fixtures
Every public feature needs two levels of testing unless the change is docs-only:
1. SQL conformance in `test/sql/` — schema, semantics, region/index behavior.
2. R package tests in `r/Rduckhts/inst/tinytest/` — wrappers, argument validation, bundled-extdata access, end-to-end DBI workflows.

Guidelines:
- One `.test` file per feature family; one tinytest file per wrapper family.
- README examples must be deterministic, short, and use bundled `inst/extdata/`.
- Prefer fixtures generated with `test/scripts/vcfpp.R` via `vcfppR`; record them in `vcfpp_manifest.tsv` where applicable.
- New fixtures needed by SQL tests should be documented in `test/scripts/prepare_test_data.sh` and copied into the R package bundle when R tests need them.
- Avoid committing large local diagnostic data files.

## Reader and Schema Policies
- Readers should be stable, projection-aware, and header-faithful unless an explicit repair/override API is being added.
- `scan_mode := 'sequential'` is the explicit full-file streaming/counting mode for readers that support it; symlink/no-index tricks are diagnostics, not API semantics.
- Do not share mutable `faidx_t`, reference-window caches, htslib iterators, or other per-reader mutable state across DuckDB worker threads. Use per-thread handles/caches.
- Decompression-thread knobs mean htslib worker threads for a file handle; do not conflate them with DuckDB processing parallelism.
- For multi-file reads, prefer generated `UNION ALL BY NAME` queries through `hts_union_query(...)` unless a task explicitly targets a native multi-file reader.

### VCF/BCF, INFO/FORMAT, and CSQ/ANN/BCSQ
- `read_bcf()` binds INFO/FORMAT columns from declared VCF header metadata. htslib warnings about undeclared record INFO keys do not imply DuckHTS exposes those keys as columns.
- Header `Type` and `Number` drive DuckDB types. For example `Type=String,Number=.` maps to `VARCHAR[]`; DuckHTS should not silently guess that an ExAC-looking string is integer.
- Builtin CSQ/ANN/BCSQ typing rules live in `src/vep_parser.c` and are intentionally allowlist-style with string fallback.
- `additional_csq_column_types` only overrides typed CSQ/ANN/BCSQ subfields after a parseable annotation schema is present in the header.
- Broken or producer-obfuscated headers should be handled by explicit helper macros/repair layers above `read_bcf()` / future reader versions, not by implicit bind-time guessing.
- Full-record VCF/BCF normalization, FORMAT remapping, genotype likelihood remapping, and semantic HTS rewriting are separate writer/remapper problems; row-level helpers must document what they do not rewrite.

## Compatibility Rewrites and Upstream Mirrors
When DuckHTS ports or rewrites an existing tool, follow compatibility-rewrite discipline:

1. identify the upstream tool and exact version/commit;
2. define the supported subset and output contract;
3. consult `.sync/` mirrors first when present;
4. validate continuously against the upstream executable or source semantics;
5. document unsupported features and validation commands;
6. credit original authors and avoid vague compatibility claims.

This applies especially to mosdepth-, bcftools-, samtools-, WisecondorX-, and VEP-inspired behavior.

## Coverage and Interval Work
Current implemented surfaces include `read_pileup`, `bam_bin_counts`, `duckhts_bam_bed_coverage`, `duckhts_mosdepth`, `duckhts_samtools_idxstats`, `fasta_nuc`, and cgranges-backed overlap functions.

Design rules:
- One function per counting model; avoid mixing pileup/depth/bin/region/compatibility semantics in one surface.
- BED-compatible output remains the interoperability contract for interval summaries.
- Fixed-width bins use arithmetic; cgranges is for irregular interval joins.
- `duckhts_mosdepth(...)` is a compatibility rewrite with mosdepth-specific file output semantics; do not fold it into generic coverage APIs.
- Memory backlog: keep `duckhts_bam_bed_coverage` tiled behavior intact, do not regress mosdepth memory parity, and treat `bam_bin_counts` streaming as a separate larger refactor.

## SIMD Kernel Development Workflow
DuckHTS SIMD dispatch is a capability-mask, per-logical-kernel system under `src/simd/`. Do not add ad-hoc ISA checks inside individual SQL functions when the logic can be represented as a reusable byte-oriented kernel.

Use this flow when adding a SIMD operation:
1. Define the logical kernel in `src/include/duckhts_simd_kernels.def`.
2. Implement and register the scalar reference first; scalar is the correctness oracle.
3. Register optional backends independently in the existing backend translation units, compile/runtime-gated by compiler/platform checks, not htslib autoconf macros.
4. Snapshot dispatch once per DuckDB vector/chunk or table-function init where a stable view is needed.
5. Wire new source files through extension and R package builds.
6. Add backend-agnostic SQL and R tests comparing forced `scalar` to `auto`.
7. Benchmark in fresh processes for process-wide backend selection and record `duckhts_simd_kernel_info()`.
8. Keep `design/simd_dispatch_matrix.md` current when changing dispatch architecture or semantics.

## Wasm / webR / Browser Rules
Wasm support must not be treated like native Linux:
- Reproduce `Rduckhts` wasm failures in the webR container/browser workflow, not host build directories.
- Inspect the built `Rduckhts_<Version>.tgz`; host `.Rcheck` trees may contain native ELF artifacts with wasm-looking names.
- Preserve Emscripten linker flags and DuckDB init-symbol exports; verify the wasm export section, not only `strings`/symbol tables.
- Keep wasm behavior gated on real Emscripten target detection.
- Keep htslib `libcurl`, S3, and GCS disabled on wasm unless the runtime model changes; browser HTTP/HTTPS goes through `src/wasm_http_hfile.c`.
- Browser constraints apply: same-origin works, remote URLs require permissive CORS, and proxy environment variables do not affect the XHR backend.
- Treat webR and duckdb-wasm as distinct runtimes and artifacts.
- Keep the Emscripten socket/i64 compatibility shim canonical in `src/include/wasm_socket_compat.h`.

## WisecondorX / BAM Dedup Context
`FILE_OFFSET UBIGINT` in `read_bam()` exposes the BGZF virtual offset after each record. Use `ORDER BY FILE_OFFSET` in window functions to reproduce exact BAM file order for streaming dedup. The SQL replication of WisecondorX larp/larp2 is in `scripts/wisecondorx_convert_conformance.py`; consult the relevant upstream/reference notes before changing these semantics.

## Style

Follow `STYLE.md`. In particular: keep ownership explicit, return errors instead of
aborting DuckDB, bound input-driven allocation, avoid fake object systems and future-only
interfaces, write R as R, preserve SQL composability, and make no benchmark claim beyond a
measured documented workload.

Do not use bare `boundary` as an explanation. Name the concrete axis and location: for
example an exon–intron junction, CDS start or end, transcript endpoint, VCF anchor-removal
case, reference/alternate differing-region edge, allocation limit, DuckDB/C ownership
interface, or DuckDB vector edge. If a genuinely generic boundary is needed, define the
term locally before using it.
