# DuckHTS architecture

Status: binding project architecture. Public behavior is defined by code, tests, and
`functions.yaml`; this file defines the ideas those artifacts must preserve.

## Thesis

DuckHTS makes sequencing data and genomics algorithms composable inside DuckDB. HTS
records, transcript models, supplementary annotations, and clinical evidence remain
ordinary typed relations or borrowed arrays. DuckDB supplies query planning, columnar
execution, parallelism, memory management, and disk spill. Small C libraries and kernels
own operations whose semantics or cost do not belong in a SQL expression.

The product is not a collection of wrapped command-line tools and it is not a new genomics
file-format ecosystem. Compatibility commands are projections over reusable readers and
kernels. A custom cache or index is justified only after an ordinary DuckDB/Parquet layout
has been measured and found insufficient for a defined workload.

## Durable bets

1. **SQL is the composition layer.** Readers and kernels return typed data that can be
   filtered, joined, grouped, nested, persisted, and served by the rest of DuckDB.
2. **C is the reusable systems layer.** Hot loops and exact biological mechanics live in
   host-neutral C with explicit borrowed views, ownership, and error returns. The same
   core can be exercised outside DuckDB and packaged through R or wasm.
3. **Use mature libraries.** htslib owns HTS transport and format semantics; cgranges,
   SIMDe, and other focused dependencies are reused where they solve the actual problem.
   DuckHTS does not reimplement a library merely to own the code.
4. **Exploit order and representation.** Coordinate sort order, packed alphabets, numeric
   contig keys, projection pushdown, and immutable structure-of-arrays models are part of
   the execution contract. They are not accidental optimizations.
5. **Bound mutable state.** Large records and variants stream. Shared reference models are
   immutable; iterators, caches, workspaces, and output buffers are worker-local. Memory
   grows with a documented active window or result bound, not with the whole input.
6. **One abstraction should serve repeated problems.** The SIMD manifest exists because
   sequence text, BAM nt16, and related reductions share a classify-and-reduce shape. New
   kernels extend that abstraction instead of adding private ISA checks to callers.
7. **Compatibility is evidence.** A rewrite names an upstream version and supported
   contract, then earns claims through differential tests, properties, sanitizers,
   statistical trials, and real data. Similar output is not conformance.
8. **Stable interfaces first.** DuckDB's stable C API is the default boundary. An unstable
   call is acceptable only when a measured requirement cannot be expressed otherwise; it
   must be isolated behind one compatibility module with an explicit removal condition.

## Layers

### Formats and transport

Table functions use htslib to expose VCF/BCF, BAM/CRAM, FASTA/FASTQ, BED, GFF/GTF, and
tabix data faithfully. Header metadata, projection, region selection, and scan mode are
reader concerns. Decompression workers are not DuckDB scan workers.

### Native kernels

Kernels receive typed arrays or compact views and know nothing about SQL strings, R
objects, or file paths. They implement sequence classification, interval operations,
coverage accumulation, normalization, liftover mechanics, and variant consequences.
Scalar reference implementations remain the correctness authority for SIMD backends.

### DuckDB adapters

Adapters bind SQL arguments, narrow checked values, acquire per-query or per-worker state,
call a kernel over vectors or scan batches, and materialize results. Registration code
does not own biological rules or large state machines.

### SQL relations and workflows

DuckDB owns imports, joins, provenance, supplementary annotation, evidence aggregation,
and final presentation. Exact annotations use sorted numeric keys; region annotations use
range planning; wide or serialized output is a final projection rather than the internal
model. Parquet and DuckLake remain ordinary storage/catalog choices, not required DuckHTS
formats.

### Packages and runtimes

`Rduckhts` packages the same extension and provides R-native validation and DBI workflows.
Native DuckDB, CRAN, MinGW, macOS, Linux, and browser runtimes share semantics but may have
different transport and threading capabilities. Platform-specific code stays behind a
small adapter boundary.

## Ownership and concurrency

- Long-lived registries are keyed by DuckDB database instance, never hidden process-wide
  singletons.
- Named transcript/reference models are immutable after publication and may coexist in
  one process.
- Mutable htslib handles, interval cursors, sequence windows, result builders, and scratch
  arenas belong to one scan or worker.
- Locks protect registry structure and lifetime transitions. They do not serialize kernel
  work that can run on isolated state.
- Every allocation has one visible owner and one failure cleanup path. Public input that
  controls allocation has a checked limit.

## Variant annotation and clinical evidence

DuckVEP's deterministic consequence engine is a C kernel over an immutable transcript
model. DuckDB prepares and verifies that model, streams sorted variants, joins population
and clinical annotations, and materializes structured consequence rows. HGVS consumes the
same lossless projected edit facts; it is not string formatting embedded in the hot
consequence loop.

Phased edits and structural events are first-class event shapes, not post-processing of
independent SNV labels. The stateful stream must group phased edits by model, transcript,
sample, phase set, and haplotype, and must retain every contributing variant identifier.

ACMG/AMP reasoning, phenotype evidence, ClinVar, gnomAD, and similar knowledge belong in
typed SQL relations above the deterministic engine. Rules should emit inspectable evidence
and provenance, not only a final opaque classification. The C engine computes genomic
facts; SQL composes the case-specific argument.

## Interface rules

- Public names describe meaning, not implementation generations. Put format versions in
  metadata, not names such as `_v2` or `_new`.
- Do not publish ignored parameters, reserved flags, placeholder result fields, or types
  for code that does not exist.
- Numeric ordinals and bitsets are valid internal contracts when a normal relation maps
  them to stable biological identifiers.
- Sortedness, coordinate conventions, ownership, nullability, and unsupported event shapes
  are explicit at every boundary.
- `functions.yaml` is the public SQL catalog. Internal headers describe only callable
  behavior present in the current tree.

## Change test

Before adding a function, format, cache, or dependency, ask:

1. Is this a reusable primitive or merely another application wrapper?
2. Can DuckDB or an existing focused library already own it?
3. What state owns the memory, and what bounds it?
4. Which sort order, representation, and cache-line behavior does the hot path rely on?
5. What pinned implementation or independent property proves correctness?
6. Does the change strengthen one of the boundaries above, or create a second authority?
