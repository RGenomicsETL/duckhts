# Coverage and reader memory

Status: current backlog for bounding coverage scratch. BED coverage tiling is implemented;
mosdepth worker tiling and `bam_bin_counts(...)` streaming remain open.

## Invariants

- Memory should grow with the active tile, contig, or output bound, not the whole input.
- Mutable coverage arrays and htslib handles remain worker-local.
- Refactoring scratch must not change counting, filtering, interval, or output-order semantics.
- htslib decompression workers and DuckHTS processing workers are different budgets.
- A smaller element type is valid only after the supported maximum and overflow behavior are
  explicit.

## `duckhts_bam_bed_coverage(...)`

The reader now walks each BED interval once and keeps lazily allocated one-megabase depth
tiles. Completed tiles are reduced and evicted as coordinate-sorted input advances. Optional
strand arrays are allocated only when strand output is requested.

This replaced query-lifetime scratch proportional to the sum of interval lengths. Preserve the
single-pass tiled walk and streaming eviction.

Remaining low-priority work is limited to measured improvements:

- reduce or alias depth arrays when their numeric and filter contracts permit it; and
- fuse repeated CIGAR walks only if profiling shows material end-to-end cost.

Do not trade the current bound for a faster whole-region allocation.

## `duckhts_mosdepth(...)`

The dominant scratch is a whole-contig `int32_t` coverage array per processing worker. In the
sequential path this is comparable to upstream mosdepth; with multiple processing workers it
becomes a DuckHTS-specific parallelism tax proportional to
`largest_contig_length * processing_threads`.

The open design is a tiled or streaming worker representation that preserves:

- CIGAR-aware updates and mate-overlap correction;
- header-order and byte-compatible output for the supported rewrite contract;
- per-base, quantized, threshold, region, and distribution outputs; and
- deterministic merging of parallel contig results.

The hard part is not allocating a smaller array; it is proving that every output can be
finalized before the relevant coverage state is discarded. Keep the sequential path available
as the compatibility oracle. The median histogram can be narrowed separately if its maximum
count is checked.

The measured WGS calibration and workload details live in
[`benchmark_riker_wgs.md`](../benchmarks/benchmark_riker_wgs.md). Re-run that workload after a
scratch-layout change instead of copying benchmark numbers into this note.

## `bam_bin_counts(...)`

The current implementation retains bin arrays for every contig until the query finishes. GC,
mapping-quality, and other optional statistics add more arrays. The open change is to emit and
release one contig once coordinate-sorted iteration has moved past it.

Before that rewrite, define:

- how unplaced or out-of-order records are handled;
- which projections can avoid optional arrays;
- the checked maximum for narrower count types; and
- how a partially emitted result reports a later malformed-record or I/O error.

This is a streaming/output-lifecycle change, not a small allocation patch, and should remain
separate from reader scheduling work.
