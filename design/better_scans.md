# Scan planning

Status: current implementation guidance for byte-balanced BAM scan pieces and explicit
`scan_mode := 'sequential'`; open design for BCF/CRAM scheduling, staged pipelines and a
possible VCF/BCF record-location surface.

## Current scan ownership

Indexed full-file scans claim units of work from a shared atomic counter. A BGZF
BAM claims byte-balanced pieces of each reference (see below); CRAM and BCF claim
each contig once from the header dictionary, and Tabix-backed VCF each contig once
from the index dictionary. VCF headers may omit
indexed contigs or name empty contigs absent from the index, especially with a
non-colocated `index_path`; the work list must use the iterator's dictionary.
Automatic BAM/CRAM scans also claim the no-coordinate tail once.
Region scans, native multi-region scans, sequential scans, and
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

`bcf_scan.h` is the host-neutral authority for file/header/iterator ownership,
physical record advancement and decode diagnostics. The DuckDB adapter owns scan
scheduling, projection, null/warn/error policy and output materialization. A scan
borrows its immutable parsed index; its caller owns the pending `bcf1_t`. Changing
the output chunk does not advance that record. Indexed full-file tasks select
literal contig IDs, not region expressions; explicit user regions retain HTSlib
parsing and brace quoting. Native tests use this same scanner without DuckDB.

An automatic or explicit-region `read_bcf` plan retains its bind-time parsed
CSI/TBI index until the plan is destroyed. This applies to automatic index
discovery and `index_path`, locally and remotely. Removal, in-place corruption,
or replacement of the index after bind cannot change that plan's offsets or
metadata counts; identical rebuilding also leaves the plan valid. SQL prepared
statements can rebind at execution, resolving the current index. DuckDB 2.0
always rebinds C table functions; `PREPARE` does not extend the index snapshot
across executions. Without a usable index at bind, auto mode streams and a
region query errors. Sequential mode does not load an index.

The caller must supply an initially matching data/index pair and keep the data
and header contents unchanged for the plan's lifetime, including remote objects.
This is parsed-index retention, not data-file snapshotting or validation of an
unrelated index supplied before bind. Count statistics do not prove association.
Every worker owns its file handle, mutable header, record, iterator and decoding
storage. Only parsed index state is shared: HTSlib 1.24's index queries are const,
and the Tabix dictionary is fully initialized before publication so query/read
lookups do not lazily allocate it. The shared helper and standalone concurrent
tests own this contract. HTSlib's VCF header reader can still probe a colocated
index independently; retaining the scan index does not promise zero index I/O.

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
  Shared list parsing rejects empty items, and HTSlib validates known-contig
  coordinates before iteration. Efficient scheduling of sparse/dense targets
  remains separate work; native deduplication does not establish that contract.

The measured full-column, one-sample workload in
[`benchmark_bcf_record_cache.md`](../benchmarks/benchmark_bcf_record_cache.md)
is retained as evidence, not a general endorsement of either scheduling policy.
It does not measure selected samples, dense region lists, remote seeks, or the
retired appender's narrower schema.

## Sequential pipelines and materialization

The `read_bcf_appender` experiment is retired, not aliased. It compared narrow
materialization and contig scheduling; it did not establish a second BCF semantic
authority or a general input-order guarantee. `CREATE TABLE AS`, `INSERT SELECT`,
and `COPY` over `read_bcf` remain the materialization interfaces. Historical
appender measurements are retained in `benchmark_bcf_appender_contigs.md` under
`benchmarks/`, and require an extension revision that still contains the experiment.

The useful model for further pipeline work is Heng Li's minibwa
[`kt_pipeline`](https://github.com/lh3/minibwa/blob/d6d9f87d300908622306382cbe17d5ffd2879d2f/kthread.c#L78)
(0.7-r421, pinned by Rminibwa): batches traverse each stage in batch-index order
while different stages overlap. Its mapping stage also uses `kt_for` internally.
Adopting that model here would require these explicit contracts, not another SQL
reader name:

- One sequential producer owns its mutable HTSlib handle. A batch owns or retains
  every record backing its borrowed fields until all consumers finish; a tidy
  record remains live across every output chunk containing its samples.
- Exact record multiplicity, input order, and coordinate sort order are different
  properties. Parallel transforms may finish out of order; an ordered consumer
  must emit by input ordinal if input-order preservation is required. Neither
  input ordinals nor ordered batches sort an unsorted file, and SQL consumers
  requiring order still need an explicit `ORDER BY` contract.
- Bound both the number of in-flight batches and their owned bytes. Cancellation
  and errors reclaim each outstanding batch exactly once; EOF is not a decode,
  capacity, or I/O failure. Worker-local handles/caches remain worker-local.
- Semantic record preservation is not byte-for-byte HTS rewriting. A narrow
  projection cannot preserve unprojected INFO/FORMAT fields; a transport offset
  is not a portable record identifier. A lossless writer needs its own contract.

`scan_mode := 'sequential'` selects full-file streaming; it does not currently
promise a staged parallel decoder. The open question is whether overlapping
decode/transform/materialization beats the ordinary projected DuckDB pipeline
at equal schema, records, output bytes, memory, and thread count.

## Byte-balanced BAM pieces

A scan that claims whole references runs on one thread when one reference holds
most of the data: a single-chromosome BAM, or a prefix of a coordinate-sorted one.
The contract for an indexed full-file BGZF BAM scan is:

- A reference with more indexed compressed data than 16 MiB is cut into at most
  256 pieces with about equal compressed bytes; every other reference is one piece.
  Boundaries come from the BAI/CSI chunk offsets of sampled positions, computed once
  at bind; no record is read to plan the scan.
- A piece is a half-open position range of one reference, read through its own
  iterator. A record belongs to the piece that holds its start position, so a record
  that crosses a boundary is returned exactly once. Records without a reference
  remain one final partition.
- `src/bam_partition.c` is the single owner of this mechanism: the partition
  builder and the ownership rule are host-neutral and tested natively
  (`test/scripts/bam_partition_test.c`, BAI and CSI, spliced reads across
  boundaries, placed-unmapped and unplaced reads). A scanner that accumulates
  fixed-width bins aligns the boundaries to the bin width so no bin is shared.
- CRAM keeps one owner per reference: its index offsets are container addresses,
  not BGZF virtual offsets. BCF keeps one owner per contig until the same evidence
  exists for its iterators. Explicit region lists, sequential scans and
  index-backed zero-column counts keep their own paths.
- Pieces do not change what is returned, only who reads it. Validation compares
  record sets with a sequential scan (count, distinct and summed `FILE_OFFSET`,
  and a record fingerprint) on uneven, empty and single-reference inputs; tests do
  not assert parallel output order.
- DuckDB runs the pieces in parallel only when the plan allows it, for example an
  aggregate, or `CREATE TABLE AS` with `preserve_insertion_order = false`; a
  plain `CREATE TABLE AS` with the default setting runs the scan on one thread.

Byte-balanced pieces replace the earlier plan to claim references in descending
weight: they balance the work directly, so a heavy late reference no longer sets
the tail time. A benchmark of this scheduling must record the reference
distribution, piece count, thread count and the comparison scan.

## VCF/BCF record location

Adding a record-location column to `read_bcf(...)` remains a public-interface decision. The
retired appender experiment exposed a BGZF position, but that did not establish
portable semantics for the general reader.

The public contract must answer all of the following before a column is added:

- whether the value identifies the start of the record or the position after it;
- whether BCF and BGZF-compressed text VCF have the same meaning;
- whether sequential and indexed iterator paths produce comparable values;
- what non-BGZF and remote paths return; and
- which stability guarantees survive parallel scans and htslib upgrades.

The name and nullability follow from those answers. Do not add the column incidentally while
changing scheduling, and do not infer a stronger contract from the appender benchmark helper.
