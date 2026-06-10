# `read_bcf` vs `read_bcf_v2` optimization comparison (2026-06-09)

Status: current optimization comparison/backlog. `read_bcf_v2(...)` is not present in the current tracked `develop`/`main` source refs, but local ignored benchmark output and Pi command logs show an earlier experimental implementation and benchmark run.

Evidence/provenance clarification:

- `benchmarks/benchmark_read_bcf_v2_current_schema.html` is an ignored local render (`benchmarks/*.html`), not tracked branch evidence.
- Current tracked refs searched with `git grep` did not contain a `read_bcf_v2` implementation.
- Recovery evidence exists in local Pi logs, especially `/tmp/pi-bash-f3081464e17784be.log`, which contains a diff with `register_read_bcf_v2_function(...)`, sample pushdown, field filtering, projected-column descriptors, persistent decode caches, and related SQL/function-catalog additions.

## Scope

This note compares the current top-level `src/bcf_reader.c` implementation with the optimization path captured in the ignored local render `benchmarks/benchmark_read_bcf_v2_current_schema.html` and the recoverable implementation diff in `/tmp/pi-bash-f3081464e17784be.log`.

It is meant to guide the next BCF reader cleanup/optimization pass while preserving:

- current `read_bcf(...)` schema compatibility;
- TigerStyle cleanup goals;
- generated R extension discipline (`r/Rduckhts/inst/duckhts_extension/` is produced by `r/Rduckhts/R/bootstrap.R`);
- no host/DuckDB-session crashes from runtime paths.

## Current `read_bcf(...)` shape

Current source: `src/bcf_reader.c`.

Implemented strengths:

- projection-aware unpack mask via `bcf_projection_unpack_mask()`;
- indexed region scans and index-backed/parallel full scans;
- count-only path using index row counts when valid;
- typed INFO/FORMAT materialization from VCF header metadata;
- typed `CSQ`/`ANN`/`BCSQ` subfield parsing through `src/vep_parser.c`;
- per-record INFO and FORMAT caches inside `bcf_read_function()` to avoid repeated `bcf_get_*()` calls for multiple projected columns in the same record;
- list-vector reserve/set-size fixes for `ALT` and `FILTER` paths;
- htslib decompression worker thread parameter.

Current limits relative to the recovered v2 work:

- no public `read_bcf_v2(...)` function exists in current tracked source;
- no sample pushdown parameter exists (`samples`, `samples_file` in the v2 report);
- no field-filter parameters exist for INFO/FORMAT/VEP (`include_info`, `include_format`, `include_vep`, `info_fields`, `format_fields`, `vep_fields` in the v2 report);
- bind-time schema always includes all detected INFO, FORMAT, and VEP columns;
- `bcf_read_function()` allocates per-chunk cache arrays sized for all bound INFO/FORMAT fields, not just projected fields;
- hot loop still switches over requested columns row-by-row and re-enters generic helpers frequently;
- VEP parsing is per record when VEP projection is needed, but the parser allocates record/transcript/value structures rather than filling projected fields from a streaming split.

## What the rendered `read_bcf_v2` report shows

Report: `benchmarks/benchmark_read_bcf_v2_current_schema.html` (ignored local artifact).

The report describes a non-destructive optimization path whose contract is:

1. retain current wide/tidy schema compatibility when parameters are equivalent;
2. benchmark an optimized current-schema scan path against legacy `read_bcf(...)` and `bcftools view -Ou -o /dev/null`;
3. leave any future array-oriented genotype schema to a separate reader.

Conformance evidence in that report:

- schema/default row conformance against `read_bcf(...)` passed with zero mismatches;
- projected core/FILTER, fixed arrays, tidy arrays, indexed region tidy, VEP projected fields, FILTER list regression, mapping-number families, number-dot INFO/FORMAT, ploidy edge, tidy chunk boundary, and count-only cases all passed;
- v2-only checks passed for sample pushdown, field filtering, VEP field filtering, and no-FORMAT tidy row behavior.

Key benchmark medians from the report:

| Workload | `read_bcf` median | `read_bcf_v2` median | Approx speedup | Notes |
| --- | ---: | ---: | ---: | --- |
| 1000G chr22 full-file core count/checksum | 37.008s | 4.532s | 8.2x | v2 also faster than `bcftools view -Ou` at 30.120s in the report |
| 1000G chr22 region core | 0.556s | 0.178s | 3.1x | indexed region scan |
| 1000G chr22 region tidy GT | 1.771s | 1.497s | 1.2x | tidy expansion dominates, v2 only modestly faster |
| 1000G chr22 region tidy GT selected sample | n/a | 0.213s | n/a | sample pushdown is the main win |
| gnomAD exome chr22 full-file core | 29.129s | 4.879s | 6.0x | sites-only workload |
| gnomAD exome region core | 0.849s | 0.275s | 3.1x | indexed region scan |
| gnomAD exome region INFO+VEP | 1.085s | 0.972s | 1.1x | annotation parsing/materialization dominates |

Interpretation:

- the large gains likely came from avoiding work, not from micro-optimizing individual scalar operations;
- sample pushdown and field filtering are higher-leverage than SIMD for most real BCF workloads;
- once VEP/INFO parsing is requested, speedups shrink because string/list materialization and annotation parsing dominate;
- tidy all-sample output remains expensive because row count explodes by sample count.

## Optimization candidates to port into `read_bcf(...)`

Prefer extending `read_bcf(...)` over adding a public `read_bcf_v2(...)` unless a compatibility break requires a separate name. A temporary internal path is fine, but public API sprawl should be avoided.

### 1. Field filtering at bind time

Add named parameters:

- `include_info := true`
- `include_format := true`
- `include_vep := true`
- `info_fields := NULL`
- `format_fields := NULL`
- `vep_fields := NULL`

Expected benefits:

- smaller DuckDB schema;
- smaller per-chunk cache arrays;
- fewer `bcf_get_info_*()` / `bcf_get_format_*()` calls;
- fewer list child-vector writes;
- much cheaper `DESCRIBE`, Parquet conversion, and projected scans when only core fields are needed.

Implementation notes:

- parse comma-separated field lists at bind time;
- validate requested fields against the header/schema and fail with `duckdb_bind_set_error()`;
- keep default behavior unchanged;
- store only included fields in `bind->info_fields`, `bind->format_fields`, and `bind->vep_schema`/a projected VEP index table;
- update `functions.yaml`, R wrappers, SQL tests, tinytests, and changelogs if public.

### 2. Sample pushdown

Add named parameters:

- `samples := NULL`
- `samples_file := NULL`

Expected benefits:

- very large tidy-mode win when the user needs one/few samples;
- avoids generating wide `FORMAT_*_<sample>` columns for unselected samples;
- reduces GT formatting and FORMAT vector indexing.

Implementation notes:

- bind-time sample selection should build `selected_sample_indices` and `selected_sample_names`;
- for htslib decode, consider `bcf_hdr_set_samples()` if compatible with the current header/record path, otherwise keep local index mapping first;
- reject `samples` + `samples_file` together;
- preserve default all-sample behavior;
- add exact conformance cases against filtering current output by `SAMPLE_ID`/selected columns.

### 3. Projection-indexed cache allocation

Current cache arrays are sized by all bound INFO/FORMAT fields. After field filtering this becomes smaller, but further improvement is possible:

- build a per-init list of projected INFO and FORMAT field indexes;
- allocate caches only for projected fields;
- map `col_id` to compact cache slot in init data;
- keep the code simple enough to satisfy TigerStyle.

This is lower risk than a full column-kernel rewrite and should reduce overhead in projected queries.

### 4. VEP streaming parser for projected fields

Current VEP path builds `vep_record_t` with transcript/value arrays. For annotation-heavy VEP work, the v2 report shows small gains, indicating this remains a hotspot.

Potential optimized path:

- for projected VEP subfields, split the CSQ/ANN/BCSQ string once;
- fill only requested field child vectors;
- avoid allocating `vep_record_t`, `vep_transcript_t`, and every `vep_value_t` when only a few fields are requested;
- retain `vep_parser.c` as the correctness/schema parser.

Do this after field filtering and sample pushdown, not before.

### 5. Hot-loop shape cleanup

The row loop is long and branch-heavy. TigerStyle cleanup and performance align here:

- split core column writing, VEP column writing, INFO writing, and FORMAT writing into small helpers;
- use htslib-style `goto cleanup` for allocation failures inside chunk processing;
- precompute column writer metadata where practical;
- keep one exit cleanup block for chunk-local allocations.

## SIMD assessment

SIMD should be used only where byte-parallel work is significant and reusable. It is not the first lever for BCF reader speed.

Likely poor SIMD targets:

- htslib BCF decoding itself (already library-owned and format-dependent);
- per-row DuckDB string assignment overhead;
- small fixed scalar assignments (`POS`, `QUAL`, flags);
- branchy FORMAT cardinality handling.

Potential SIMD/kernel targets:

1. **Delimited string scans**
   - comma/pipe scanning for `CSQ`/`ANN`/`BCSQ`, `INFO` string lists, and maybe GT text formatting helpers;
   - a reusable byte kernel could find delimiters/invalid bytes faster for long annotations.
2. **Missing/vector-end masks**
   - scanning `int32_t` / `float` buffers for missing and vector-end sentinels before list writes;
   - useful only for large arrays, less so for common short FORMAT fields.
3. **ASCII normalization/validation**
   - uppercasing alleles/bases and lightweight validation in future consequence/HGVS devices.
4. **Sequence context kernels**
   - later DuckVEP/HGVS/codon work may benefit more than `read_bcf()` itself.

If adding SIMD:

- define logical kernels in `src/include/duckhts_simd_kernels.def`;
- scalar reference first;
- backend-specific implementations under `src/simd/`;
- compare forced `scalar` vs `auto` in SQL and R tests;
- update `design/simd_dispatch_matrix.md`.

Do **not** add ad-hoc ISA branches in `bcf_reader.c`.

## Test plumbing to restore before porting v2 optimizations

The report references a benchmark source not currently present in the tree. Recreate source artifacts before implementation:

- `benchmarks/benchmark_read_bcf_v2_current_schema.Rmd` or successor `benchmarks/benchmark_read_bcf_projection_pushdown.Rmd`;
- `make bench-read-bcf-v2` or renamed target;
- SQL conformance cases corresponding to the report;
- R tinytests for wrapper argument validation and generated SQL.

Suggested SQL tests:

- default schema unchanged;
- `include_info := false, include_format := false, include_vep := false` returns only core columns;
- `info_fields := 'SB'` exposes only requested INFO;
- `format_fields := 'GT'` exposes only requested FORMAT;
- `vep_fields := 'Consequence,DISTANCE'` exposes only requested VEP columns;
- malformed/unknown fields fail at bind time without crashing;
- `samples := 'S1'` removes unselected wide columns and emits only selected tidy rows.

Suggested R tests:

- wrapper accepts and quotes the new arguments correctly;
- `samples` and `samples_file` are mutually exclusive;
- bundled fixtures cover one wide and one tidy selected-sample query.

## Proposed order of work

1. Recreate benchmark Rmd/Makefile target from the rendered report so future claims are reproducible.
2. Implement `include_info/include_format/include_vep` only; test core-only schema and count/core scans.
3. Add `info_fields/format_fields/vep_fields` field filtering.
4. Add sample pushdown.
5. Refactor `bcf_read_function()` chunk-local allocations to `goto cleanup` style and compact projected caches.
6. Re-run the cohort benchmark grid and compare to the old rendered report.
7. Only then consider SIMD kernels for CSQ delimiter scanning or missing-mask scans.

## Acceptance gate for merging optimization work

- No default schema or row changes for existing `read_bcf(...)` calls.
- Focused SQL tests pass.
- R package bootstrap regenerates extension copies; no direct generated C edits.
- `git diff --check` passes.
- Benchmark report records dataset paths, row counts, medians, and whether `bcftools` was available.
- NEWS files mention any user-visible parameters or behavior changes.
