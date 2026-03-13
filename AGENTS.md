# DuckHTS Agent Guidelines

This document provides guidance for AI agents working on the DuckHTS DuckDB extension.

## Project Goal
Build a DuckDB 1.4+ extension that **reads** HTS file formats using htslib, with coverage/interval analytics primitives for CNV and QC workflows. The R package (`Rduckhts`) is published on CRAN. Application code will handle format conversion (e.g., to parquet).

## Current Status
- ✅ All readers implemented: VCF/BCF, SAM/BAM/CRAM, FASTA/FASTQ, GTF/GFF, tabix
- ✅ Sequence UDFs, SAM flag predicates, k-mer UDFs implemented
- ✅ R package published on CRAN (`Rduckhts 0.1.3-0.0.2`, dev `0.1.3-0.0.2.9001`)
- ✅ Comprehensive SQL tests implemented
- 🔄 Phase 10: Coverage & interval primitives for interoperable CNV workflows (in design)

## Target DuckDB Version
- **Minimum**: DuckDB 1.4 (C/C++ extension API)
- Keep all public APIs compatible with DuckDB 1.4+.

## Documentation Conventions
- **Primary documentation**: README.Rmd (R Markdown format)
- **Changelog discipline is mandatory**: every user-visible change must update both top-level [NEWS.md](NEWS.md) for the extension and [r/Rduckhts/NEWS.md](r/Rduckhts/NEWS.md) for the R package in the same change, with **newest entries first**.
- Extension function documentation and community-extension descriptor metadata are sourced from [functions.yaml](functions.yaml); keep it updated whenever functions/macros/UDFs are added, removed, renamed, or re-described.
- Regenerate the package-bundled function catalog and the checked-in community extension descriptor from `functions.yaml` with `python3 scripts/render_function_catalog.py`.
- Keep the generated files in `r/Rduckhts/inst/function_catalog/` and [community-extensions/extensions/duckhts/description.yml](community-extensions/extensions/duckhts/description.yml) in sync with the manifest.
- New specs go in [./.github/SPEC.md](SPEC.md).
- The implementation plan is tracked in [./.github/PLAN.md](PLAN.md).
- Do not hand-edit the generated `community-extensions/extensions/duckhts/description.yml`; update `functions.yaml` and re-render it instead.
- Never consider a feature complete until both changelogs have been updated.

## Agent Working Instructions
**CRITICAL**: Always read existing code and implementation before making changes:
1. Read the relevant source files to understand current implementation
2. Check existing tests to understand expected behavior
3. Review current build system before making changes
4. Preserve existing functionality and API compatibility
5. When mirroring external tool behavior, consult the local upstream mirrors under `.sync/` before relying on secondary summaries

## Source Layout Expectations
- Extension sources live under [src/](src/).
- Public headers live under [src/include/](src/include/).
- Vendor sources live under a dedicated top-level folder (see SPEC.md).
- Test SQL lives under [test/sql/](test/sql/).
- R package structure in [r/duckhts/](r/duckhts/).

## Vendoring Rules
- Vendor scripts must be **reproducible** and pin exact versions/commits.
- All downloads must be checksum-verified.
- Scripts must support **offline rebuilds** once sources are fetched.
- Do not modify vendored code directly; patch with explicit patch files.
- Capture upstream licenses in a dedicated licenses folder and document them.

## Build & CI Rules
- Keep build scripts deterministic and non-interactive.
- Do not add network access to the extension build step.
- Network access is allowed **only** in explicit vendoring scripts.
- **CRAN requirement**: No vcpkg usage - use CMAKE and traditional autotools
- Assume cmake, make, and standard build tools available on R platforms

## R Package Specific Rules
- Package should be self-contained as much as possible
- Use CMAKE and configure/configure.win scripts for cross-platform builds
- Include cleanup and cleanup.win scripts for proper cleanup
- Simplify package bootstrapping to copy necessary extension files
- Keep `r/Rduckhts/README.Rmd` wired to the generated function catalog instead of duplicating extension function lists by hand.
- When upstream extension sources under `src/` change, update the bundled R package copy by running `Rscript bootstrap.R ~/duckhts/` from `r/Rduckhts/`, then reinstall the package before running tests.
- When new extension source files are added under `src/`, update all build/bootstrap paths that enumerate sources:
  - top-level `CMakeLists.txt`
  - `r/Rduckhts/R/bootstrap.R`
  - `r/Rduckhts/configure`
  - `r/Rduckhts/configure.win`
  - any package-side source manifests copied during bootstrap
- Treat a new public function as incomplete until its C source is wired through the extension build and the R package build on Unix and Windows.
- After adding/removing/renaming public functions, update `functions.yaml`, run `python3 scripts/render_function_catalog.py`, bootstrap the R package copy, and verify the generated catalog/descriptor stay in sync.
- Never run `R CMD INSTALL .` from `r/Rduckhts/`: it mutates `inst/duckhts_extension/htslib` in place and removes the vendored htslib source tree needed for subsequent installs. Build a tarball with `R CMD build .` and install the resulting tarball instead.
- Version scheme: duckhtsVersion-x format
- All R package modifications must maintain CRAN compatibility

## Testing Guidelines
- Current tests cover: VCF/BCF, BAM/CRAM, FASTA/FASTQ, GTF/GFF, tabix readers
- Tests use DuckDB `.test` format with type codes: `I` (integer), `T` (text), `R` (real), `B` (boolean)
- Use `__WORKING_DIRECTORY__` placeholder for test data paths
- Plan tests at two levels for every public feature:
  - **Extension-level SQL conformance** in `test/sql/` for schemas, semantics, region/index behavior, and SQL-vs-native parity
  - **R package-level tests** in `r/Rduckhts/inst/tinytest/` for wrapper signatures, argument validation, bundled-extdata access, packaged extension loading, and end-to-end DBI workflows
- Prefer one `.test` file per feature family instead of growing `test/sql/duckhts.test` indefinitely.
- Prefer one `tinytest` file per wrapper family or workflow instead of one monolithic integration script.
- New fixtures required by SQL tests should be generated or documented under `test/scripts/prepare_test_data.sh`, then copied into the R package bundle as needed for package tests.
- Package-level test planning must include:
  - installed-package tests (`tinytest::test_package("Rduckhts")`)
  - tarball build/install path (`R CMD build`, then install tarball)
  - `R CMD check` coverage for any R-facing changes
- Keep README examples deterministic, short, and backed by bundled `inst/extdata` wherever practical; if README rendering is blocked by the harness, note that explicitly.
- Tests should cover:
  - schema correctness
  - record counts  
  - region queries (indexed files)
  - FORMAT/INFO parsing (VCF/BCF)
  - paired FASTQ handling
  - CRAM with reference files
  - attributes_map functionality for GTF/GFF
- Phase 10 coverage tests should additionally cover:
  - BED reader schema inference, passthrough BED4+/extra columns, and indexed `.bed.gz` region reads
  - `fasta_nuc` parity against known interval composition fixtures
  - SQL-vs-native parity (`bam_bin_counts` vs `read_bam` + GROUP BY)
  - strand invariant: `count_total == count_fwd + count_rev`
  - empty regions returning zero rows (not errors)
  - include/exclude flag behavior for duplicates, proper pairs, secondary/supplementary, QC-fail, and unmapped reads
  - overlapping mates, discordant pairs, cross-contig mates
  - coordinate convention correctness (0-based half-open)
  - `bam_bedcov`, `bam_coverage`, and `bam_depth` parity against documented samtools semantics on stable fixtures
  - export round-trip: write `.bed.gz` → `read_tabix` → verify counts
- Test R package loading and function calls

## Best Practices Learned
- Keep Windows builds on MinGW/Rtools; match DuckDB platform strings (e.g., windows_amd64_mingw) in extension metadata.
- Use $ORIGIN rpaths on Linux so libhts.so can be found at runtime.
- Prefer pkg-config (plus Homebrew fallbacks) for OpenSSL/libcurl detection on macOS.
- Avoid plugins on Windows; keep plugin support behind configure checks on Unix.
- Bundle all runtime data under inst/extdata for R and keep README examples runnable with local files.
- Keep README examples deterministic and short; use eval=FALSE only when external network access is required.

## Heng Li Component Scope
- Treat `cgranges` as the interval engine (overlap/contain/index) for BED/GRanges-style operations.
- Treat `seqtk` ideas (`kseq`/buffered stream parsing patterns) as parser implementation guidance for FASTA/FASTQ-like read paths.
- Do not conflate `cgranges` and `seqtk`; they address different layers.
- When adding interval APIs, keep base-level BAM/CRAM pileup/depth logic in htslib-native code paths and use `cgranges` for interval algebra.
- `cgranges` is vendored at `third_party/cgranges/` and needs to be wired into `CMakeLists.txt` for Phase 10 work.

## Phase 10 — Coverage & Interval Primitives

### Architecture Overview
Phase 10 adds coverage/interval analytics as orthogonal building blocks for CNV, fetal-CNV, and QC workflows. The design separates:
- **Immediate**: independent utility functions that should be implemented first (`bgzip`, `bgunzip`, `bam_index`, `bcf_index`, `tabix_index`)
- **Layer 0**: SQL-first baseline using existing `read_bam` + UDFs (already available)
- **Layer 1**: BED/interval primitives (`read_bed`, `fasta_nuc`, `interval_merge`, `interval_overlap`, `interval_nearest`)
- **Layer 2**: Native counting kernels (`bam_bin_counts`, `bam_bedcov`, `bam_coverage`, `bam_depth`, `bam_pileup`)
- **Layer 3**: Indexed export/interoperability over canonical BED-like outputs
- **Later**: candidate/active-site detection and `mpileup`-style richer pileup summaries once the coverage kernels are stable

### Key Design Rules
- **Optimize for a small canonical API**: prefer one function per counting model rather than proliferating wrappers and macros.
- **Use BED-compatible outputs as the main interoperability contract** (`chrom`, `start`, `end`, metrics...).
- **Keep `bam_bin_counts` first-class in the immediate coverage design**: its input contract and stable output schema must cover both WisecondorX-style combined counts and NIPTeR-style strand-split counts without separate public kernels.
- **Do not conflate bin counting with pileup**: read-start binning (WisecondorX/NIPTeR) and per-base depth (Rsamtools) are separate primitives with separate implementations.
- **Do not conflate coverage with candidate-site discovery**: active/candidate-site heuristics are a later analysis layer on top of stable coverage/pileup primitives.
- **Expose filtering policy primarily through explicit `include_flags`, `require_flags`, `exclude_flags`, and quality thresholds**; add higher-level knobs only when an upstream behavior cannot be expressed clearly that way.
- **Prefer symbolic SAM flag names in user-facing docs and R wrappers**; keep raw integer masks as the low-level/compatibility path underneath.
- **Do not hide region merge/combine logic inside BAM scan functions**: use interval preprocessing UDFs.
- **Fixed-width bins use arithmetic** (`pos / bin_width`); `cgranges` is for irregular interval joins only.
- **Avoid SQL macros unless they remove real duplication without obscuring semantics.**
- **Exported `.bed.gz`/`.tsv.gz` files must be tabix-indexable and round-trip through `read_tabix`.**
- **All public Phase 10 functions must have R wrappers and documentation in `functions.yaml`.**
- **All new APIs must be benchmarked against the SQL-first baseline** (`read_bam` + GROUP BY).

### New Source Files
| File | Purpose |
|------|---------|
| `src/hts_index_builder.c` | `bam_index`, `bcf_index`, `tabix_index` — index building (Phase 8) |
| `src/bgzip.c` | `bgzip`, `bgunzip` — BGZF compression/decompression (Phase 8) |
| `src/interval_udf.c` | `read_bed`, `fasta_nuc`, `interval_merge`, `interval_overlap`, `interval_nearest` |
| `src/bam_bin_counts.c` | `bam_bin_counts` table function |
| `src/bam_bedcov.c` | `bam_bedcov` table function |
| `src/bam_coverage.c` | `bam_coverage` table function |
| `src/bam_depth.c` | `bam_depth` table function |
| `src/bam_pileup.c` | `bam_pileup` table function |

### Semantic Parameters
- `count_model := 'read_start'` for fixed-bin CNV counting; pileup/depth are separate functions
- `include_flags`, `require_flags`, and `exclude_flags` are the primary filtering controls for pair, duplicate, secondary, supplementary, QC-fail, and unmapped behavior
- `strand_mode := 'combined' | 'split'`
- overlap/proper-pair behavior must be documented per function; do not assume bin counting, depth, pileup, and BED coverage should share the same mate-handling rules
- `bam_bin_counts` input contract must support regular fixed bins plus optional region restriction, explicit index/reference paths, and optional annotation/mask inputs without changing the canonical output columns
- `bam_bin_counts` output contract must remain stable across WisecondorX-like and NIPTeR-like use:
  - always `chrom`, `start`, `end`, `bin_id`, `count_total`, `count_fwd`, `count_rev`
  - `strand_mode := 'combined'` still emits `count_fwd` / `count_rev` for compatibility, with `count_total = count_fwd + count_rev`
- WisecondorX-like and NIPTeR-like outputs should be expressed as documented parameter sets over one canonical kernel/output contract, not as `profile := 'wisecondorx'` / `profile := 'nipter'` public APIs
- contig selection must be explicit where it matters:
  - support `region`/indexed restriction for targeted scans
  - plan explicit contig include/exclude controls for workflows that want autosomes-only or custom chromosome subsets
- Coordinate contract: 0-based half-open on disk (BED); SQL outputs may add convenience columns.

### DuckDB Parallelism Pattern (CRITICAL)
All new table functions should start from the proven contig-level parallelism pattern from `bam_reader.c`, but the exact threading/iterator strategy must be re-checked per function against htslib 1.23 behavior and workload shape before the API is frozen:

1. **Global init** (`duckdb_table_function_set_init`):
   - Set `duckdb_init_set_max_threads(info, min(n_contigs, 16))` when index is available and no region filter.
   - Store an atomic contig counter in global state for thread coordination.

2. **Local init** (`duckdb_table_function_set_local_init`):
   - Each DuckDB thread opens its **own** `samFile*`, `sam_hdr_t*`, `hts_idx_t*`, `bam1_t*`.
   - Call `hts_set_threads(fp, 2)` per handle for htslib I/O decompression.
   - Multi-region queries use `sam_itr_regarray()` with htslib-internal dedup.

3. **Contig claiming** (in scan or helper):
   - `__sync_fetch_and_add(&global->current_contig, 1)` to claim the next contig.
   - `sam_itr_queryi(idx, tid, 0, HTS_POS_MAX)` for full-contig iteration.
   - No lock contention — each thread processes its own contigs independently.

4. **Key htslib APIs** (vendored 1.23):
   - `sam_open`, `sam_hdr_read`, `sam_hdr_nref`, `sam_hdr_tid2name`, `sam_hdr_tid2len`
   - `sam_index_load3(fp, path, index_path, HTS_IDX_SILENT_FAIL)`
   - `sam_itr_queryi`, `sam_itr_regarray`, `sam_itr_next`
   - `hts_set_opt(fp, CRAM_OPT_REFERENCE, ref)` for CRAM support

**Verify before freezing:**
- Confirm the exact multi-region iterator API and duplicate-suppression semantics in vendored htslib 1.23 before standardizing on `sam_itr_regarray()` in docs or implementation.
- Confirm whether `hts_set_threads(fp, 2)` per DuckDB worker is still the right default once multiple DuckDB threads and multi-region scans are in play; avoid hard-coding a pattern that oversubscribes CPU or regresses small-region workloads.
- Treat contig-level partitioning as the default fast path for full indexed scans, not as a blanket rule for every region-restricted or pileup/depth workload.
- Keep base-level coordinate conventions (`bam_depth`, `bam_pileup`) explicitly under review; do not assume BED-style 0-based half-open semantics without confirming compatibility targets.
- Re-check overlap/proper-pair behavior against upstream mirrors before claiming compatibility:
  - WisecondorX counts read starts, requires proper pairs for paired reads, and does not have pileup-style mate-overlap handling
  - NIPTeR bins read positions from `scanBam()` output and is not explicitly pair-aware in its binning path
- Re-check chromosome include/exclude behavior before freezing APIs for autosome-only or sex-chromosome-aware workflows.

For `bam_bin_counts` specifically, per-contig bin arrays are **independent** — no cross-thread synchronization needed during accumulation. Each thread scans a contig, fills its bin array, emits rows, frees, and claims the next contig.

### Detailed Plan
See [.github/PLAN.md](.github/PLAN.md) Phase 10 for full architecture, code sketches, conformance matrix, and implementation sequence.

## Current Focus Areas
1. **Implement First**: `bgzip`, `bgunzip`, `bam_index`, `bcf_index`, `tabix_index`
2. **Next Interval Layer**: implement `read_bed` and `fasta_nuc` before the BAM coverage/count operations
3. **Coverage Ops**: implement and validate `bam_bin_counts`, `bam_bedcov`, `bam_coverage`, `bam_depth`, and `bam_pileup` against WisecondorX/NIPTeR/samtools/Rsamtools semantics
4. **Phase 10 Design Freeze**: finalize the minimal public API and output contracts for `read_bed`, `fasta_nuc`, `bam_bin_counts`, `bam_bedcov`, `bam_coverage`, `bam_depth`, and `bam_pileup`
4. **cgranges Integration**: wire vendored cgranges into build for irregular BED/interval joins
5. **R Surface**: add and test R wrappers for every public function, and keep `functions.yaml`/generated catalogs synchronized
6. **Later Analysis Layer**: defer candidate/active-site detection and `mpileup`-style richer site summarization until the coverage primitives are stable

## Style
- Keep changes minimal and focused.
- Preserve existing code style and APIs unless the task explicitly requires changes.
- Follow R package conventions for CRAN submission.
- Test changes with R CMD check when modifying R package components.
