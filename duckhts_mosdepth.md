# `duckhts_mosdepth` Rewrite Plan

This file is the implementation plan for a native `duckhts_mosdepth(...)`
rewrite inside duckhts.

The goal is not "a coverage function inspired by mosdepth". The goal is a
rewrite that follows the principles in [`.sync/rewrite.bio.txt`](/root/duckhts/.sync/rewrite.bio.txt):

- credit the original authors
- emulate exactly
- be explicit about validation scope
- build small, validate continuously

This note is intentionally code-oriented. It names the files, structs,
functions, and validation steps that would be needed for an implementation.

## 1. Compatibility target

### Pinned target for v1

Target exact compatibility with:

- upstream tool: `mosdepth`
- validated runtime version: `0.3.9`
- validation date in this repo: `2026-04-14`

Local verification commands run on `2026-04-14`:

```bash
mosdepth --version
python3 scripts/mosdepth_conformance.py \
  HG00106.chrom11.ILLUMINA.bwa.GBR.exome.20130415.bam \
  --extension build/release/duckhts.duckdb_extension
python3 scripts/mosdepth_benchmark.py \
  HG00106.chrom11.ILLUMINA.bwa.GBR.exome.20130415.bam \
  --chrom 11 \
  --modes fast,default \
  --extension build/release/duckhts.duckdb_extension \
  --runs 1 \
  --verify
```

Observed local results on that date:

- the original SQL reconstruction was useful as an early oracle, but it is no
  longer the validation target
- `scripts/mosdepth_conformance.py` and `scripts/mosdepth_benchmark.py` should
  now exercise only the native `duckhts_mosdepth(...)` table function
- any remaining SQL-first coverage logic belongs to generic interval/depth
  primitives, not to the mosdepth rewrite workflow

### Important versioning warning

The `.sync/mosdepth` mirror appears newer than the locally installed
`mosdepth 0.3.9` binary. For example, [`.sync/mosdepth/mosdepth.nim`](/root/duckhts/.sync/mosdepth/mosdepth.nim)
contains a usage string for `mosdepth 0.3.13`.

That means:

- do not claim exact compatibility to "mosdepth" in general
- do not claim exact compatibility to the current `.sync/mosdepth` mirror
- first pin the rewrite to one exact upstream version and record that choice in
  repo docs

Before implementation starts, do one of:

1. replace `.sync/mosdepth` with the exact `0.3.9` source tree, or
2. add a version-pinning note in a validation doc mapping:
   - local binary version used for conformance
   - upstream git tag / commit used as source reference

## 2. What must match exactly

For v1, `duckhts_mosdepth(...)` should match mosdepth's behavior and output
contract for the supported feature subset.

### Output files

Given `prefix`, the rewrite should produce the same file names as mosdepth:

- `{prefix}.mosdepth.global.dist.txt`
- `{prefix}.mosdepth.summary.txt`
- `{prefix}.per-base.bed.gz` unless `no_per_base := TRUE`
- `{prefix}.regions.bed.gz` if `by := ...`
- `{prefix}.mosdepth.region.dist.txt` if `by := ...`
- `{prefix}.quantized.bed.gz` if `quantize := ...`
- `{prefix}.thresholds.bed.gz` if `thresholds := ...`

For the supported modes, "same" means:

- text files: byte-for-byte identical
- gzipped BED outputs: decompressed text must be byte-for-byte identical
- accompanying `.csi` indexes: query-equivalent with tabix; compressed bytes do
  not have to be identical

### Semantics that count

The rewrite must preserve:

- default excluded flags: `1796`
- `threads` meaning: htslib BAM/CRAM decompression threads only
- `fast_mode` semantics
- default-mode CIGAR-aware coverage semantics
- default-mode mate-overlap correction semantics
- `by := <window>` and `by := <bed>` output shape
- summary and distribution formatting
- file naming and header ordering

The rewrite must not silently "improve" outputs.

## 3. Recommended public API

Implement as a DuckDB table function with side effects, similar in spirit to
existing `bgzip(...)` and `bam_index(...)`.

### Proposed SQL signature

```sql
duckhts_mosdepth(
  prefix,
  path,
  chrom := NULL,
  by := NULL,
  no_per_base := FALSE,
  fasta := NULL,
  threads := 0,
  flag := 1796,
  include_flag := 0,
  fast_mode := FALSE,
  fragment_mode := FALSE,
  quantize := NULL,
  mapq := 0,
  min_frag_len := -1,
  max_frag_len := -1,
  thresholds := NULL,
  use_median := FALSE,
  read_groups := NULL,
  index_path := NULL,
  overwrite := FALSE
)
```

### Proposed returned columns

Return one row with:

- `success BOOLEAN`
- `prefix VARCHAR`
- `summary_path VARCHAR`
- `global_dist_path VARCHAR`
- `per_base_path VARCHAR`
- `regions_path VARCHAR`
- `region_dist_path VARCHAR`
- `quantized_path VARCHAR`
- `thresholds_path VARCHAR`

Use `NULL` for outputs that were not requested.

### Why a dedicated function

Do not fold this into `bam_coverage(...)` in v1.

Reason:

- `duckhts_mosdepth(...)` is a compatibility rewrite with file-writing side
  effects and exact mosdepth output contracts
- `bam_coverage(...)` should stay free to be a more SQL-native and possibly
  more parallel coverage primitive later

These are different products and should not share a public surface in v1.

## 4. Source files to add or change

### New extension source files

Add:

- `src/mosdepth_table.c`
- `src/include/mosdepth_table.h`

Optional split if the file gets too large:

- `src/mosdepth_engine.c`
- `src/include/mosdepth_engine.h`
- `src/mosdepth_output.c`
- `src/include/mosdepth_output.h`

### Existing files to update

Build wiring:

- `CMakeLists.txt`
- `r/Rduckhts/bootstrap.R`
- `r/Rduckhts/configure`
- `r/Rduckhts/configure.win`

Docs / wrappers:

- `functions.yaml`
- `r/Rduckhts/R/duckhts.R`
- `r/Rduckhts/inst/function_catalog/functions.yaml`
- `r/Rduckhts/inst/function_catalog/functions.tsv`
- `r/Rduckhts/inst/function_catalog/functions.md`

Tests:

- `test/sql/mosdepth.test`
- `r/Rduckhts/inst/tinytest/test_mosdepth.R`

Release notes:

- `NEWS.md`
- `r/Rduckhts/NEWS.md`

Validation docs:

- `VALIDATION.md` or `benchmark_mosdepth.md`

## 5. Core implementation design

### 5.1 Bind / init / execute model

Implement as a normal DuckDB table function:

- bind validates arguments and constructs immutable options
- init sets `duckdb_init_set_max_threads(info, 1)`
- function body performs the rewrite once and emits a single result row

This function should not use DuckDB scan parallelism.

Reason:

- mosdepth's `threads` parameter is for htslib decompression, not for external
  compute workers
- exact compatibility is easier if one call owns the whole output set and emits
  files in mosdepth-compatible order

### 5.2 Proposed C structs

Add a bind-data struct:

```c
typedef struct {
    char *prefix;
    char *path;
    char *chrom;
    char *by;
    char *fasta;
    char *index_path;
    char *quantize;
    char *thresholds;
    char *read_groups;
    int threads;
    int mapq;
    int min_frag_len;
    int max_frag_len;
    uint16_t exclude_flag;
    uint16_t include_flag;
    bool no_per_base;
    bool fast_mode;
    bool fragment_mode;
    bool use_median;
    bool overwrite;
} duckhts_mosdepth_bind_t;
```

Add a result-data struct:

```c
typedef struct {
    bool done;
    bool success;
    char *prefix;
    char *summary_path;
    char *global_dist_path;
    char *per_base_path;
    char *regions_path;
    char *region_dist_path;
    char *quantized_path;
    char *thresholds_path;
    char *error_message;
} duckhts_mosdepth_result_t;
```

Add a per-contig stats struct mirroring mosdepth summary logic:

```c
typedef struct {
    uint64_t cum_depth;
    uint64_t cum_length;
    uint32_t min_depth;
    uint32_t max_depth;
} duckhts_mosdepth_depth_stat_t;
```

Add a BED region struct:

```c
typedef struct {
    char *chrom;
    uint32_t start;
    uint32_t stop;
    char *name; /* nullable */
} duckhts_mosdepth_region_t;
```

### 5.3 Main engine functions

Suggested function boundaries:

- `mosdepth_bind(...)`
- `mosdepth_init(...)`
- `mosdepth_function(...)`
- `mosdepth_run(const duckhts_mosdepth_bind_t *bind, duckhts_mosdepth_result_t *out)`
- `mosdepth_run_contig(...)`
- `mosdepth_scan_alignment(...)`
- `mosdepth_emit_outputs_for_contig(...)`
- `mosdepth_write_summary_line(...)`
- `mosdepth_write_distribution(...)`
- `mosdepth_write_thresholds(...)`
- `mosdepth_write_quantized(...)`

## 6. Reading and coverage algorithm

### 6.1 htslib APIs to use

Use htslib directly, not `read_bam(...)`:

- `sam_open()`
- `sam_hdr_read()`
- `sam_index_load3()`
- `sam_itr_queryi()`
- `sam_itr_next()`
- `hts_set_threads()`
- `hts_set_opt(fp, CRAM_OPT_REFERENCE, ...)`
- `hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS, ...)`
- `hts_set_opt(fp, CRAM_OPT_DECODE_MD, 0)`

For remote BAM / CRAM support:

- use `sam_open()` for transport
- load indexes with `sam_index_load3(fp, path, index_path, HTS_IDX_SAVE_REMOTE)`

Do not bypass htslib transport with custom readers.

### 6.2 Required-fields optimization

Mirror upstream mosdepth's CRAM field minimization from
[`.sync/mosdepth/mosdepth.nim`](/root/duckhts/.sync/mosdepth/mosdepth.nim:975).

For all modes request at least:

- `SAM_FLAG`
- `SAM_RNAME`
- `SAM_POS`
- `SAM_MAPQ`
- `SAM_CIGAR`
- `SAM_TLEN`

For non-fast mode also request:

- `SAM_QNAME`
- `SAM_RNEXT`
- `SAM_PNEXT`

If `read_groups` filtering is requested, also request:

- `SAM_RGAUX`

This is important for CRAM performance and for fidelity to upstream behavior.

### 6.3 Per-contig processing model

For exact mosdepth-style output, process contigs sequentially in header order.

For each contig:

1. construct a whole-contig iterator with `sam_itr_queryi(idx, tid, 0, HTS_POS_MAX)`
2. allocate a dense `int32_t` delta array of length `target_length + 1`
3. scan records for that contig
4. apply mode-specific delta updates
5. cumulative-sum the array to concrete depth
6. write outputs and update summary/distribution state
7. free or reuse the array for the next contig

Do not use `read_bam(...)` contig-claiming parallelism for this compatibility
rewrite.

Reason:

- mosdepth's public `threads` flag is decompression-only
- a parallel contig-claiming engine changes memory behavior materially
- exact output-order compatibility is simpler with sequential header-order
  processing

### 6.4 Fast-mode update rule

Mirror upstream logic:

```text
arr[rec.start] += 1
arr[rec.stop]  -= 1
```

This uses the aligned record span and ignores internal CIGAR operators.

### 6.5 Default-mode update rule

Mirror upstream logic:

- expand reference-consuming aligned blocks from the CIGAR
- increment block starts
- decrement block ends

Implement a helper:

```c
static void mosdepth_inc_coverage_from_cigar(const bam1_t *rec, int32_t *arr);
```

This should treat only reference-consuming query-consuming operators as depth
contributors, analogous to mosdepth's `gen_start_ends()` / `inc_coverage()`.

### 6.6 Mate-overlap correction

This is the hardest exact-compatibility requirement.

Mirror upstream logic from
[`.sync/mosdepth/mosdepth.nim`](/root/duckhts/.sync/mosdepth/mosdepth.nim:291):

- only in default mode
- only for proper pairs
- not for supplementary reads
- only when both mates are on the same contig
- first overlapping mate is stored in a `seen` map keyed by `QNAME`
- when the partner arrives, subtract the overlap interval(s)

Implementation detail:

- for simple one-CIGAR-block pairs, use the cheap whole-overlap correction
- otherwise merge both mates' block start/end events, track `pair_depth`, and
  subtract intervals where pair depth reaches `2`

Suggested helpers:

- `mosdepth_should_track_overlap(const bam1_t *rec)`
- `mosdepth_record_take_or_store_overlap(...)`
- `mosdepth_apply_overlap_simple(...)`
- `mosdepth_apply_overlap_blocks(...)`

Data structure:

- use a khash table or uthash keyed by QNAME storing a copied `bam1_t *` or a
  compact overlap struct
- free entries immediately when the mate arrives

Do not approximate overlap correction. This is one of the main reasons users
trust mosdepth default mode.

### 6.7 Fragment mode

Implement later, after fast/default parity is complete.

When added, mirror upstream:

- count only read1 from proper pairs
- compute fragment start as `min(start, matepos)`
- increment fragment start
- decrement `fragment_start + abs(isize)`

## 7. Output writers

### 7.1 Summary

Write exactly the same header and columns as mosdepth:

- `chrom`
- `length`
- `bases`
- `mean`
- `min`
- `max`

Formatting:

- use the same decimal precision policy as mosdepth
- support `MOSDEPTH_PRECISION` if we want environment compatibility, or
  document a fixed precision if we deliberately do not mirror that behavior

For v1 compatibility, mirroring `MOSDEPTH_PRECISION` is preferred.

### 7.2 Global and region distributions

Mirror mosdepth's cumulative distribution:

- histogram is per exact depth
- emitted value for depth `d` is proportion of bases covered at least `d`
- write per chromosome and `total`

Use the same sparse-trailing-zero behavior as mosdepth; do not emit a denser
distribution format just because it is easier.

### 7.3 Per-base BED.gz

Generate run-length encoded output from the depth array:

- contiguous runs with same depth become one BED interval
- output order must match mosdepth header-order traversal

Suggested helper:

```c
static int mosdepth_write_rle_depth(BGZF *bgzf, const char *chrom,
                                    const int32_t *depth, int64_t n);
```

### 7.4 `by := <window|bed>` regions output

Implement one region writer that handles both:

- integer window size
- BED file input

For each region, compute:

- mean depth by default
- median if `use_median := TRUE`

For BED input:

- if there is a 4th column, propagate it as the region name
- otherwise emit `unknown`, matching mosdepth threshold behavior

### 7.5 Quantized output

Quantized output is not a new coverage algorithm. It is a post-processing pass
 over the per-base depth array.

Implementation:

- parse quantize breakpoints into sorted integer thresholds
- map each depth to a label/bin
- merge adjacent positions whose quantized label is the same

### 7.6 Thresholds output

Thresholds output is also post-processing over region intervals:

- for each region in `by`
- count bases with depth `>=` each requested threshold
- write one extra output column per threshold

This should be implemented after basic region means are exact.

### 7.7 Compression and indexing

For BED-like outputs:

- write BGZF directly with htslib
- build CSI indexes after close

Use:

- `bgzf_open()` or `hts_open()` in BGZF write mode
- `tbx_index_build3()` or equivalent htslib index builder after the file is
  closed

Do not shell out to `bgzip` or `tabix` from the extension.

## 8. Scope phases

### Phase 1: minimal compatible rewrite

Implement:

- indexed BAM input
- optional CRAM input with explicit reference
- `chrom`
- `threads`
- `flag`
- `include_flag`
- `mapq`
- `fast_mode`
- `no_per_base`
- summary
- global distribution
- `by := <window>`
- `by := <bed>`
- `regions.bed.gz`
- `mosdepth.region.dist.txt`

Validation target:

- exact match to mosdepth `0.3.9` for those features

### Phase 2: default mode exactness

Implement:

- CIGAR-aware default mode
- exact mate-overlap correction

This is the key phase because the benchmark shows SQL fast mode is already
competitive, while SQL default mode is still slower than mosdepth.

### Phase 3: remaining mosdepth features

Implement:

- `fragment_mode`
- `quantize`
- `thresholds`
- `use_median`
- `read_groups`

### Explicitly deferred

Do not implement in v1 unless a real downstream need appears:

- D4 output
- unsupported experimental flags not used in current workflows

## 9. Validation and tests

### 9.1 Keep the current SQL scripts

The current scripts are not throwaway prototypes. They are the executable spec:

- `scripts/mosdepth_conformance.py`
- `scripts/mosdepth_benchmark.py`

Native code should be validated against upstream mosdepth using those scripts
or native equivalents.

### 9.2 New SQL and R tests

Add SQL tests:

- `test/sql/mosdepth.test`

Cover:

- fast mode single chromosome
- default mode single chromosome
- `by := 10000`
- `by := bed`
- output path creation
- bad-argument validation
- CRAM reference requirement

Add R tests:

- `r/Rduckhts/inst/tinytest/test_mosdepth.R`

Cover:

- wrapper signature and defaults
- output path existence
- end-to-end invocation from DBI

### 9.3 Golden-data validation

Keep validated outputs under a repo-controlled location only if they are small
enough and stable.

At minimum, record:

- BAM / CRAM used
- exact mosdepth version
- exact commands
- exact expected output files

This belongs in a repo doc, not just CI logs.

## 10. Documentation and attribution

### README requirements

When `duckhts_mosdepth(...)` is added publicly, README/docs must include:

- that it is a native rewrite of mosdepth behavior
- explicit credit and citation guidance for upstream mosdepth
- an AI assistance disclosure consistent with `rewrites.bio`
- exact validation scope and version

### Suggested attribution text

Document that users of `duckhts_mosdepth(...)` should cite upstream mosdepth
when the output semantics are used in scientific work.

Do not bury this in a footer; make it visible in README and function docs.

## 11. Why this should be done natively

Current status:

- SQL already proves correctness
- SQL fast mode is already very good
- SQL default mode is slower than mosdepth because SQL must:
  - explode CIGAR into rows
  - join mates for overlap correction
  - materialize large intermediate relations

Native code can beat that by:

- keeping per-contig depth in one dense in-memory array
- applying overlap correction during the scan
- computing summary/distribution/region outputs from the same pass
- avoiding relational materialization of internal state

That is the exact "think big, work small" path from `rewrites.bio`:

- big architecture idea: one native pass produces the full mosdepth output set
- small iteration loop: implement fast mode first, validate, then default mode

## 12. Non-goals

For v1, do not:

- silently deviate from mosdepth output formats
- reinterpret `threads` as DuckDB compute threads
- mix this compatibility rewrite into future generic coverage APIs
- claim compatibility against an unpinned upstream version
