# DuckHTS coding style

Status: binding style and review contract for first-party C, R, and SQL. Vendored code
keeps upstream style and is changed only through explicit patches.

This is DuckHTS's translation of TigerStyle. The objective is bounded behavior, visible
ownership, reviewable state machines, and errors that survive an embedded DuckDB process.
It is not a brace-placement contest and it is not permission to manufacture tiny objects.

## C

### Errors and invariants

- Malformed files, invalid SQL arguments, unsupported runtimes, and resource exhaustion
  return a DuckDB-visible error or an explicit status value. Production code does not
  abort the host process.
- Use `_Static_assert` for compile-time relationships: integer widths, fixed layouts,
  enum/table cardinality, and constant bounds.
- Runtime `assert()` belongs only in opt-in tests, fuzzers, and sanitizer builds. If a
  condition matters in production, check it and return an error.
- Preserve the first useful error. Do not overwrite a precise parser/kernel failure with
  a generic wrapper message.

### Bounds and arithmetic

- Check addition, multiplication, narrowing, and offset-plus-length before allocation or
  byte access. A failed check must not publish a wrapped or partial result.
- External recursion, nesting, queues, active windows, and growable buffers have explicit
  limits. Prefer an iterative cursor when the natural workload is a stream.
- Coordinates name their convention (`pos1`, half-open end, inclusive end) at the type or
  field boundary. Convert once, after validation.

### Ownership and allocation

- Mark views as borrowed in the type or interface comment. Pair every owned object with one
  destroy function and make its allocator family obvious.
- Prefer immutable shared models and worker-owned scratch. Do not share mutable htslib
  iterators, FASTA windows, DuckDB connections, or kernel cursors across workers.
- Allocate or grow at bind/init/tile boundaries and reuse storage. No per-base,
  per-transcript-pair, or per-output-row heap allocation in a hot loop.
- `goto cleanup` is acceptable when it makes one ownership exit auditable. It is better
  than duplicated partial teardown.

### Functions and control flow

- Separate input parsing, validation, state mutation, kernel execution, and output
  materialization when they are distinct phases.
- Function length is a review signal, not a hard numeric law. Split a long function when
  the split names a real invariant or state transition; keep it whole when helpers would
  hide the dataflow.
- Use early returns for failed preconditions and keep the successful path visually clear.
  At a public boundary, separate checks when that produces a more precise error.
- Avoid object-like C layers, generic context bags, callback indirection, and helper
  factories unless multiple real implementations require them.
- Comments explain invariants, ownership, coordinate transforms, or non-obvious bit tricks.
  They do not record implementation phases, agents, old alternatives, or promised code.

### Names and layout

- Prefix externally visible symbols with `duckhts_` or the owning library name such as
  `duckvep_`. Keep file-private helpers `static`; do not export a helper in anticipation of
  another caller.
- Name operations with verbs and stored facts with nouns. Names should expose the unit or
  coordinate convention when confusion is possible, such as `pos1`, `byte_len`, or
  `cds_start`.
- Use `size_t` for in-memory byte counts and indices, fixed-width integers for persisted or
  model contracts, and checked conversion at the boundary between them.
- Follow the surrounding file's indentation and brace style. Keep new lines near 100
  columns when practical, but do not hide a clear expression behind formatting helpers or
  reformat unrelated code.
- Declare ownership-bearing variables near the phase that acquires them. A cleanup block
  should be readable as the reverse of acquisition.

### Data and performance

- Prefer structure-of-arrays for large homogeneous models and streams. Use compact fixed
  widths only after the supported domain has been checked.
- Treat sortedness as an API property. Preserve cursors across batches where the host
  lifecycle permits it; reset explicitly on rewind or model change.
- A branchless mask or lookup table is welcome when it is measured and readable. Add a
  short comment for a non-obvious encoding; do not replace clear scalar code on faith.
- Every SIMD operation has a scalar reference, capability registration, forced-scalar
  tests, and backend-independent semantics. Callers do not perform private ISA dispatch.
- Benchmark the full materialization path and the pure kernel separately. Report input
  rows, expanded work, output rows/bytes, threads, model, and machine.

### Dependencies and portability

- Reuse htslib and focused C libraries before writing a replacement. A new dependency must
  justify its capability, license, vendoring, CRAN/HPC cost, wasm behavior, and maintenance
  surface.
- Keep the core on the DuckDB stable C API. If an unstable entry point is unavoidable,
  isolate it in one `*_compat.c` boundary and document the exact unblock it provides.
- Keep target checks about compiler/platform capability. Do not use unrelated htslib
  configure macros as CPU or DuckDB feature tests.

## R

- Public wrappers should feel native to R: vector inputs, predictable recycling or explicit
  scalar checks, standard conditions, and DBI-compatible return values.
- Validate arguments at the wrapper boundary and preserve the underlying DuckDB error as
  the cause. Do not translate Python control flow, dictionaries, or class hierarchies into
  R syntax.
- Prefer small helpers and base/r-lib conventions over new dependencies. Every dependency
  must be suitable for CRAN and used for more than convenience syntax.
- Keep examples deterministic and runnable from bundled data. Build and install package
  tarballs; never mutate the package source through `R CMD INSTALL .`.

## SQL

- Relations are the public composition mechanism. Keep stable biological identifiers and
  provenance in ordinary tables; use compact numeric keys only inside measured joins and
  kernels.
- Make required ordering explicit at the boundary that consumes it. Do not rely on an
  incidental subquery or file order that the planner is free to change.
- Prefer set operations, projection, and typed nested results over row-by-row macros or
  serialized strings. Render VCF/CSQ/HGVS compatibility text only at the edge.
- Let DuckDB plan filters, projections, joins, row-group pruning, memory, and spill before
  inventing a custom store. Benchmark exact lookup and range workloads before adding an
  index format.
- Keep exact-key, interval, and computed annotations as distinct physical operations. They
  may share prepared variant keys and final output, but not a fake universal lookup path.
- SQL that constructs a resident model or invokes a native kernel must expose its schema,
  assembly/release identity, input hashes, and sort contract.

## Interfaces and documentation

- Public APIs describe behavior implemented today. No reserved flags, ignored arguments,
  version-suffixed experiments, or placeholder fields.
- `functions.yaml`, executable tests, and C headers are contracts. Design notes contain
  durable invariants or one unresolved decision; git history and GitHub issues contain the
  path and backlog.
- Remove obsolete comparisons and completed plans instead of marking them historical
  forever. Never commit chain-of-thought, agent handoffs, or milestone theatre.

## Review checklist

- Is there one authority for each semantic decision?
- Is every allocation bounded and owned?
- Are arithmetic and narrowing checked before use?
- Can malformed input crash DuckDB or read outside a buffer?
- Is mutable state worker-local, and is lock ownership visible?
- Does the hot path exploit its promised order and representation?
- Is the public surface free of code that does not exist?
- Do properties, differential cases, sanitizers, and package tests cover the changed
  boundary?
