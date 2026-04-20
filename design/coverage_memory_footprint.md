# Coverage / Reader Memory Footprint — Audit and Plan

Status: open
Scope: `src/bam_bed_coverage.c`, `src/mosdepth_table.c`, `src/bam_bin_counts.c`, `src/bcf_reader.c`

## Motivation

Benchmark of `duckhts_bam_bed_coverage` vs `samtools coverage` on NA12878 chr11:

| Workload | DuckHTS RSS | samtools RSS |
| --- | --- | --- |
| Capture BED | ~10 MB | ~10 MB |
| WGS full chr11 | **2,100 MB** | 26 MB |
| WGS first 100 × 50k windows | **129 MB** | 26 MB |

Parity on reads/covbases is exact; the gap is purely memory. Walking the
codebase turned up the same structural patterns in three more readers.

## Shared root causes

A small set of patterns recurs across the readers:

1. **Eager, query-lifetime per-region / per-contig scratch.** Peak RSS =
   Σ over all regions, not max. Fix: allocate at the start of processing
   a region and free before moving on.
2. **No tiling for large regions / contigs.** A single 135 Mbp or 249 Mbp
   region still needs one huge buffer. Fix: process in fixed-size tiles
   (~1 Mbp), fold summaries on the fly, reset scratch between tiles.
3. **Oversized element types.** `int64_t` / `uint64_t` / `uint32_t` used
   where `uint32_t` / `uint16_t` / a bitset suffices.
4. **Per-base arrays for quantities that fold on the fly.** Sum-of-depth
   is Σ clipped CIGAR match lengths; does not need an array.
5. **Optional buffers always allocated.** Strand-split arrays, GC/MQ
   caches, FORMAT/INFO caches allocated even when the corresponding
   output is disabled or unprojected.
6. **Dead allocations / stale scratch variables.**

## Punch list

### `src/bam_bed_coverage.c` — step 1 done

All step 1 items are in tree. Summary of what landed:

- [x] Pattern 1: BED-load-time `ensure_region_depth_arrays` call
  removed from `load_bed_regions`; scratch is lazy per tile.
- [x] Pattern 2: tiling added with
  `DUCKHTS_BEDCOV_TILE_SIZE = 1_000_000`.
- [x] Pattern 4: `summed_cov_*` accumulated in
  `flush_completed_tiles`, gated by `>= min_depth` (matches the
  original emit-loop semantics — see test coverage below).
- [x] Pattern 5 (strand gating): `ensure_tile_depth_arrays(tile,
  strand_outputs)` skips `depth_fwd_post` / `depth_rev_post` when
  `strand_outputs := FALSE`.
- [x] Pattern 6: dead `nbytes` / `(void)nbytes;` removed.
- [x] **Single-pass tiled walk.** One `sam_itr_next` loop per region
  at [`:660-705`](../src/bam_bed_coverage.c#L660-L705). Each record
  is decoded once: scalars via `accumulate_record_summary`, depth
  via `accumulate_record_on_tile`.
- [x] **Streaming tile eviction.** Tiles are a per-region linked list
  created lazily in `get_or_create_tile` and flushed in
  `flush_completed_tiles` as soon as `tile->end <= rec->core.pos`.
  Peak scratch is effectively `O(tile_size × live_tiles)` with
  `live_tiles` ≈ 1–2 in the steady state for normal read lengths,
  rather than `O(region_len)`.
- [x] `min_depth > 1` meandepth behaviour pinned to `samtools
  coverage` and covered by new SQL tests in
  [`test/sql/bam_bed_coverage.test`](../test/sql/bam_bed_coverage.test)
  (`min_depth := 2` + `require_flags := 64`; `min_depth := 100`),
  mirrored in
  [`r/Rduckhts/inst/tinytest/test_bam_bed_coverage.R`](../r/Rduckhts/inst/tinytest/test_bam_bed_coverage.R).

Current state: **done except for further type-size tightening and
compaction.**

API note: `duckhts_bam_bed_coverage(...)` now also takes an explicit
`decompression_threads` named parameter (forwarded to
`hts_set_threads` when non-zero; see
[`:1015`](../src/bam_bed_coverage.c#L1015) and
[`:623-624`](../src/bam_bed_coverage.c#L623-L624)). This is separate
from the reserved `processing_threads` stub.

Remaining follow-ups (low priority, not blocking):

- [ ] Pattern 3: bitset / `uint16_t` for `depth_*`. Current peak is
  ~16 MB per live tile-set with `strand_outputs=TRUE`; a bitset
  (sufficient for `min_depth == 1`) would drop that to ~500 KB.
- [ ] Safe pre/post scratch aliasing: when post-filters reduce to the
  pre pass (`mapq == 0 && min_baseq == 0 && min_read_len == 0 &&
  require_flags == 0 && exclude_flags == default`) alias `depth_post`
  onto `depth_pre` instead of a second `calloc`.
- [ ] Optional CIGAR-walk fusion: `accumulate_record_summary` and
  `accumulate_record_on_tile` both iterate the CIGAR; a fused walker
  that takes both `region` and the active tile would halve the parse
  cost per record. Worth doing only if profiling flags it.

Code-polish items (agreed, still to land):

- Make `flush_completed_tiles` return `void` (all callers have dead
  `if (... != 0)` error paths).
- Add a one-line invariant comment on `get_or_create_tile` noting
  that the linear scan is safe because streaming eviction keeps the
  live tile list at ~1–2 entries.

### `src/mosdepth_table.c` — same per-base scratch family, but dominated by per-worker whole-contig coverage

Reference: [`:1485-1498`](../src/mosdepth_table.c#L1485-L1498),
[`:1998-2006`](../src/mosdepth_table.c#L1998-L2006),
[`:96-100`](../src/mosdepth_table.c#L96-L100),
[`:253-260`](../src/mosdepth_table.c#L253-L260),
[`:1246`](../src/mosdepth_table.c#L1246),
[`:1303`](../src/mosdepth_table.c#L1303).

- Worker coverage buffer is `int32_t[chrom_len + 1]`, reused per contig
  but not tiled (Pattern 2). chr1 (~249 Mbp) ≈ 1 GB per worker;
  scales with `processing_threads`. Same allocation shape in the
  sequential path. **Empirical calibration** — whole-WGS `--by 50000`
  default mode, `processing_threads = 0`, via
  `scripts/mosdepth_benchmark.py`: upstream 1960 MB RSS, DuckHTS
  1969 MB RSS (~0.5% higher, 1.33× faster on wall-clock). So in
  sequential mode memory is at upstream parity; the "1 GB per
  worker" cost is a *parallelism tax*, not a structural regression,
  and tiling this buffer is primarily useful when
  `processing_threads > 1`.
- `mosdepth_count_stat_t.counts` is a `uint64_t[65536]` = 512 KB
  scratch histogram, allocated once per `write_window_regions(...)` /
  `write_bed_regions(...)` call when `use_median = TRUE`, then reused
  across windows / BED regions within that call via
  `count_stat_clear(...)`. Not retained one-per-region. `uint32_t`
  suffices → 256 KB (Pattern 3). Worthwhile shrink but not a
  multiplicative footprint like the contig coverage array.

This is *not* the same allocation shape as pre-fix
`bam_bed_coverage.c`: the old bedcov problem was scratch retention
scaling with Σ region_len across the whole BED. For mosdepth the
footprint is instead `max(contig_len) × processing_threads`. That
distinction changes fix priorities — tiling the worker coverage
buffer is the only meaningful lever here.

**Target after fix:** tiled worker scratch ≤ 4 MB/worker; median
histogram 256 KB / call.

#### Aside on `processing_threads` vs upstream `--threads`

Upstream mosdepth's `--threads` is an htslib BGZF/decompression
thread count — `open(bam, …, threads = threads, …)` at
[`.sync/mosdepth/mosdepth.nim:961`](../.sync/mosdepth/mosdepth.nim#L961)
— and `main(...)` iterates contigs sequentially. DuckHTS'
`processing_threads` is different: contig-claiming parallelism that
opens independent BAM handles per worker. The two knobs are
orthogonal; the DuckHTS `duckhts_mosdepth(...)` API exposes both.
Practical scaling of `processing_threads` is bounded by workload
imbalance — one or two workers end up carrying the largest
chromosomes while others go idle — which is the same reason
single-contig workloads (e.g., chr11-only benchmarks) see little
benefit past ~1–2 workers.

### `src/bam_bin_counts.c` — bigger refactor

Reference: [`:316-368`](../src/bam_bin_counts.c#L316-L368),
[`:464`](../src/bam_bin_counts.c#L464).

- `ensure_contig_capacity` grows 3 required + up to 7 optional arrays
  per contig; all contigs stay resident until query end (Pattern 1).
  Peak = Σ across contigs, not max.
- Per-bin counts are `int64_t`; `uint32_t` fits (Pattern 3, 2×).
- GC / MQ arrays allocated whenever `want_gc` / `want_mq` is set, even
  if those columns are not projected (Pattern 5).

This one needs a streaming rewrite: emit a contig's bins as soon as the
iterator crosses to the next contig, then free. Not a drop-in fix.

### `src/bcf_reader.c` — **not a memory footprint issue** (correction)

An earlier draft of this note claimed the FORMAT/INFO caches at
[`:1411-1441`](../src/bcf_reader.c#L1411-L1441) and
[`:1454-1462`](../src/bcf_reader.c#L1454-L1462) should be projection-
gated. That claim was wrong. The reader already has projection
pushdown enabled
([`:2055`](../src/bcf_reader.c#L2055),
[`:1068`](../src/bcf_reader.c#L1068)) and
`bcf_projection_unpack_mask`
([`:338-373`](../src/bcf_reader.c#L338-L373)) drives `bcf_unpack()`
level from the projected column set, so the heavy per-field value
buffers (`fmt_i32[i]` etc.) are only populated when the emit loop
actually emits that column.

What remains at that site is pointer-array scaffolding sized by
*declared* header fields — on a typical multi-sample VCF, roughly
`n_format_fields × 6 × sizeof(ptr)` ≈ 1 KB per scan call, not scaling
with the sample count. Small optimisation: move `vectors[]` and the
pointer scaffolds into `bcf_init` so they are allocated once per query
rather than once per batch. This is on the order of tens of KB saved
in a long-running scan, not a footprint fix. Not worth a dedicated
change; drop to a "while you're in there" cleanup.

`bam_reader.c` has the same projection-pushdown wiring
([`:1259`](../src/bam_reader.c#L1259)) — no parallel concern there.

## Proposed ordering

1. ~~**`bam_bed_coverage.c` step 1**~~ — **done.** Single-pass tiled
   walk, `min_depth` / `summed_cov` semantics pinned with SQL tests,
   strand-array gating, and streaming tile eviction all landed. Only
   follow-ups remaining are type-size tightening / aliasing /
   optional CIGAR-walk fusion, tracked inline in the section above.
2. **`mosdepth_table.c` tiling** — reuse the tiling idiom from step 1
   in the per-worker coverage buffer. Separate PR; shares helpers if
   any fall out of step 1.
3. **`mosdepth_table.c` count_stat type shrink** — mechanical.
4. **`bam_bin_counts.c` streaming rewrite** — largest, defer until the
   above are done. Revisit alongside any work tracked in
   `design/better_scans.md` on contig-claiming, since the streaming
   emit pattern overlaps.

(`bcf_reader.c` dropped from the ordering — see corrected note above.
Projection pushdown already gates the heavy buffers.)

## Non-goals

- No new public columns. No output-schema changes. No changes to
  samtools/mosdepth parity semantics. Floating-point rounding noise in
  the current benchmark is unrelated and stays out of scope.
- No parallelism changes for `bam_bed_coverage.c` in the first pass;
  keep the contig-claim work separate.

## Verification

For each step:

- SQL test in `test/sql/` confirming unchanged output vs. current
  fixtures.
- R package tinytest unchanged.
- Benchmark re-run on the NA12878 chr11 workloads above; record RSS
  in `benchmarks/` alongside the existing numbers.
