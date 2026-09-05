# Scan planning

Status: open design for indexed full-scan scheduling and a possible VCF/BCF record-location
surface. Explicit `scan_mode := 'sequential'` is implemented for supported readers.

## Current boundary

Indexed full-file scans in `read_bam(...)` and `read_bcf(...)` claim each
contig once in header order. Region scans, native multi-region scans, sequential scans, and
index-backed zero-column counts are separate paths and must keep their existing semantics.

Parallel scans do not promise physical output order. A downstream workflow that needs a
stable order must request one explicitly with `ORDER BY` over a documented key.

Symlinks, hidden indexes, and no-index copies are useful diagnostics, not public scan modes.

## Canonical BCF reader contract

The `read_bcf_v2` experiment is retired, without an alias. Its useful record-lifetime
invariant now belongs to `read_bcf`: projected descriptors and decode buffers are
worker-owned; loaded values and parsed CSQ remain valid until the next input record,
even when one record's tidy samples span several DuckDB output chunks. Reuse is not
sealed allocation: HTSlib decode buffers and annotation storage can still grow.

The experiment does not justify changing these independent policies together:

- Full-file streaming is `scan_mode := 'sequential'`; removing the experiment does
  not change canonical `auto` scheduling. Count-only index statistics are not file
  validation. The experiment's unparsed VCF-line count path is removed.
- DuckDB SQL projection controls materialization. Experimental bind-time field
  lists are removed rather than maintained as a second projection language.
  Earlier VCF-text parsing limits and selected-CSQ decoding need their own
  malformed-input, projection-equivalence, and throughput evidence before adoption.
- Sample selection is a distinct future interface, not equivalent to projecting
  FORMAT columns. HTSlib header selection must be applied consistently at bind and
  worker init; indexed BCF also needs FORMAT subsetting. Unknown/all-excluded
  samples, header order, ploidy, and tidy cardinality must be specified together.
- Batched sample emission must retain the pending record and reacquire list-child
  pointers after a reserve. It is not required for correct resumable tidy emission;
  the experimental sample-run path is removed.
- Explicit region lists use HTSlib 1.24 `bcf_itr_regarray`/`tbx_itr_regarray` for one
  visit per matching stored record. This preserves duplicate records in the file.
  Do not add SQL `DISTINCT` or independent overlapping shard scans as substitutes.
  Input-list validation and efficient scheduling of sparse/dense targets remain
  separate work; native deduplication does not establish those contracts.

The measured full-column, one-sample workload in
[`benchmark_bcf_record_cache.md`](../benchmarks/benchmark_bcf_record_cache.md)
is retained as evidence, not a general endorsement of either scheduling policy.
It does not measure selected samples, dense region lists, remote seeks, or the
appender's narrower schema. A pipeline may preserve batch arrival order without
producing coordinate-sorted output; its consumer must name the required order.

## Weighted contig claiming

The remaining scheduling change is to claim expensive contigs before cheap ones so one large
late contig does not dominate tail latency. This is an internal scheduling change, not a new
reader parameter or result-set change.

The design constraints are:

- retain one-owner-per-contig iteration and the existing atomic claim model;
- derive weights from index statistics when they are available;
- sort by descending weight with original contig ordinal as the deterministic tie-break;
- fall back to header order when statistics cannot be read reliably;
- do not split contigs into overlapping region strings; and
- do not extend this change to explicit or chained region scans.

Before implementation, settle whether a zero-weight contig may be omitted. Keeping it at the
end is the conservative default because an index statistic should not become a silent data
filter. Any shared claim-order helper must have one owner and be used by both readers only if
their fallback behavior is genuinely identical.

Validation should compare row sets and counts with the current paths on uneven and empty
contigs. Tests must not assert incidental parallel output order. A benchmark is useful only if
it records the contig distribution, thread count, and tail-time change.

## VCF/BCF record location

Adding a record-location column to `read_bcf(...)` remains a public-interface decision. The
experimental `read_bcf_appender(...)` can expose a BGZF position, but that does not establish
portable semantics for the general reader.

The public contract must answer all of the following before a column is added:

- whether the value identifies the start of the record or the position after it;
- whether BCF and BGZF-compressed text VCF have the same meaning;
- whether sequential and indexed iterator paths produce comparable values;
- what non-BGZF and remote paths return; and
- which stability guarantees survive parallel scans and htslib upgrades.

The name and nullability follow from those answers. Do not add the column incidentally while
changing scheduling, and do not infer a stronger contract from the appender benchmark helper.
