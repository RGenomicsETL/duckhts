# Scan planning

Status: open design for indexed full-scan scheduling and a possible VCF/BCF record-location
surface. Explicit `scan_mode := 'sequential'` is implemented for supported readers.

## Current boundary

Indexed full-file scans in `read_bam(...)` and the original `read_bcf(...)` path claim each
contig once in header order. Region scans, chained multi-region scans, sequential scans, and
index-backed zero-column counts are separate paths and must keep their existing semantics.

Parallel scans do not promise physical output order. A downstream workflow that needs a
stable order must request one explicitly with `ORDER BY` over a documented key.

Symlinks, hidden indexes, and no-index copies are useful diagnostics, not public scan modes.

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
