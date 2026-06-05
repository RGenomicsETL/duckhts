# FASTQ Throughput Notes

Status: open performance investigation. The direct htslib-transport FASTQ parser and packed pipeline described here remain future work; check `functions.yaml` before assuming any proposed parameter has landed.

This file captures the current understanding of why `.sync/fastp` can achieve
much higher FASTQ throughput than `read_fastq(...)`, how `fastp` uses threads,
and what changes are realistic for duckhts without breaking remote file
support.

## 1. Verified: how `fastp` gets its throughput

Reference mirror: `.sync/fastp`

### Threading model

- `-w/--thread` is the worker-thread count, default `3`.
- Single-end mode uses:
  - `1` reader thread
  - `w` processor threads
  - optional writer thread(s)
- Paired-end mode uses:
  - `2` reader threads, or `1` interleaved reader thread
  - `w` processor threads
  - multiple writer threads

Relevant files:

- `.sync/fastp/src/main.cpp`
- `.sync/fastp/src/options.cpp`
- `.sync/fastp/src/seprocessor.cpp`
- `.sync/fastp/src/peprocessor.cpp`

### Work batching

- `fastp` batches reads into `PACK_SIZE = 1000`.
- Reader threads produce `ReadPack`s.
- Packs are handed to workers through per-worker single-producer /
  single-consumer queues.
- Backpressure limits the number of packs kept in memory.

Relevant files:

- `.sync/fastp/src/common.h`
- `.sync/fastp/src/seprocessor.cpp`
- `.sync/fastp/src/peprocessor.cpp`

### Compression and decompression path

- For ordinary `.gz`, `fastp` uses ISA-L inflate.
- If the input is BGZF, `fastp` switches to a custom `BgzfMtReader`.
- `BgzfMtReader` uses:
  - `1` reader thread
  - a fixed decompression thread pool
  - a ring buffer for backpressure
- `fastp` computes a separate BGZF decompression budget instead of folding it
  into `-w`.

Relevant files:

- `.sync/fastp/src/fastqreader.cpp`
- `.sync/fastp/src/fastqreader.h`
- `.sync/fastp/src/bgzf.h`

### Output path

- For gzipped output and `thread > 1`, `fastp` can use a `pwrite` mode.
- Each worker compresses its own output with `libdeflate`.
- Writes happen to non-overlapping file regions.

Relevant file:

- `.sync/fastp/src/writerthread.cpp`

### Important nuance

`fastp` does not call `read()` concurrently on one reader object. Its
`FastqReader` explicitly states that `read()` on the same object is not
thread-safe. The parallelism comes from the pipeline around the reader, not
from multiple threads sharing one parser object.

## 2. Verified: what duckhts does today

Current file:

- `src/seq_reader.c`

### Current FASTQ scan path

- `read_fastq(...)` opens FASTQ through htslib `sam_open()`.
- It reads one record at a time with `sam_read1()`.
- The htslib FASTQ path parses text records into a `bam1_t`.
- duckhts then decodes that `bam1_t` into DuckDB vectors and strings.

Relevant files:

- `src/seq_reader.c`
- `third_party/htslib/sam.c`

### Current threading behavior

- FASTA/FASTQ scans currently force `duckdb_init_set_max_threads(info, 1)`.
- So the normal FASTQ scan path is single-threaded at the DuckDB scan level.

Relevant file:

- `src/seq_reader.c`

### Current count-only fast path

There is already one important exception:

- zero-column `COUNT(*)` for FASTQ does not use `sam_read1()`
- it uses `hts_open()` plus `hts_getline()` directly
- it skips FASTQ records from the raw text stream

This matters because it proves that duckhts can parse FASTQ directly on top of
htslib streams without going through `bam1_t`.

Relevant file:

- `src/seq_reader.c`

## 3. Verified: what htslib threading does and does not provide

Vendored htslib contains a real `hts_set_threads()` implementation:

- `third_party/htslib/hts.c`

Current behavior:

- if format is `sam`, it dispatches to `sam_set_threads()`
- else if compression is `bgzf`, it dispatches to `bgzf_mt(...)`
- else if format is `cram`, it sets CRAM threading
- otherwise it returns without a FASTQ-specific threaded parse path

Relevant files:

- `third_party/htslib/hts.c`
- `third_party/htslib/sam.c`

Important consequence:

- adding a `decompression_threads` parameter to `read_fastq(...)` would be a
  real and valid change
- but there is no evidence in vendored htslib of a dedicated multithreaded
  FASTQ parser path analogous to `fastp`

Inference:

- an htslib-threading parameter may help BGZF-backed FASTQ inputs
- it is unlikely by itself to close most of the gap to `fastp` for ordinary
  `fastq.gz`

This inference should be treated as a performance expectation, not as a claim
about exact speedup without benchmarking.

## 4. Remote support is a hard design constraint

duckhts is not free to replace htslib transport with a local-file-only reader.

Why:

- current readers rely on `hts_open()` / `hts_getline()` transport
- wasm `http` / `https` support is implemented as an htslib `hFILE` backend
- remote access semantics live below the reader in htslib transport

Relevant files:

- `src/seq_reader.c`
- `src/duckhts.c`
- `src/wasm_http_hfile.c`

Practical consequence:

- directly transplanting `fastp`'s `FILE *`-centric reader architecture into
  duckhts would risk breaking remote file support
- the safe optimization boundary is:
  - keep htslib transport
  - replace only the FASTQ record-parsing layer above that transport

In other words, the preferred architecture is:

- `hts_open()` / `hts_getline()` / hFILE backends
- our FASTQ parser
- DuckDB vector materialization

not:

- local `FILE *`
- custom decompressor / parser that bypasses htslib transport entirely

## 5. Recommended implementation path

### Small change

Add a `decompression_threads` parameter to `read_fastq(...)` and call
`hts_set_threads()` on the handle.

Pros:

- low risk
- simple API extension
- consistent with existing `read_bam(..., decompression_threads := ...)`

Cons:

- likely limited upside for ordinary FASTQ throughput

### Medium change

Keep the current DuckDB table-function shape, but replace the FASTQ
`sam_read1()` path with direct FASTQ parsing over `htsFile`.

Core idea:

- stay on `hts_open()` / `hts_getline()`
- parse FASTQ records ourselves
- decode only projected columns
- write directly into DuckDB vectors
- avoid `fastq_parse1()` and `bam_set1()` overhead

Why this is the best cost / benefit move:

- preserves remote support
- preserves wasm HTTP support
- uses an architecture already proven by the count-only fast path
- avoids the biggest obvious overhead in the current implementation

### Large change

Build a `fastp`-style packed pipeline on top of htslib transport:

- one read / decompress stage
- one or more parse / materialize worker stages
- ordered emit into DuckDB scan chunks

This is much larger because it must preserve:

- DuckDB table function semantics
- projection pushdown behavior
- paired and interleaved validation
- deterministic error handling
- remote transport compatibility

## 6. Size estimate

### Small

- `read_fastq(..., decompression_threads := ...)`
- parameter wiring
- docs and tests

Expected size:

- a few hundred lines across source, docs, and tests

### Medium

- direct FASTQ parser over `htsFile`
- single-threaded first
- projection-aware decode

Expected size:

- one focused feature-sized change in `src/seq_reader.c`
- plus SQL and R tests

### Large

- packed multistage pipeline with worker threads
- still remote-safe via htslib transport

Expected size:

- a subsystem-level change

## 7. Recommended next step

If the goal is to make `read_fastq(...)` stop being much slower than a
specialized FASTQ tool, the right sequence is:

1. add the small htslib-threading parameter if desired
2. implement the medium direct-parser path over `htsFile`
3. benchmark again
4. only then decide whether a large pipeline rewrite is justified

The key conclusion is:

- `fastp` is faster mostly because it is a specialized packed FASTQ pipeline
- the right duckhts response is not to abandon htslib transport
- the right duckhts response is to stop paying the `FASTQ text -> bam1_t ->
  DuckDB` conversion cost on the main scan path
