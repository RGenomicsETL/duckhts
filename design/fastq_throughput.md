# FASTQ throughput

Status: open performance investigation. A direct parser over htslib transport and any packed
worker pipeline remain unimplemented; `functions.yaml` is authoritative for the current
`read_fastq(...)` surface.

## Observed bottleneck

The projected FASTQ scan path in `src/seq_reader.c` reads text through htslib's SAM/sequence
reader into `bam1_t` records and then materializes DuckDB vectors. A zero-column `COUNT(*)`
already bypasses `bam1_t` and skips records with `hts_getline(...)`, showing that direct parsing
above htslib transport is viable.

FASTA and FASTQ scans currently advertise one DuckDB scan worker. htslib can provide threaded
BGZF decompression, but it does not provide a specialized parallel FASTQ parser for ordinary
gzip input. A decompression-thread parameter might help some BGZF workloads; it should not be
presented as the solution without measurements by compression format.

## Comparison caliber

The local fastp mirror at commit `c903553df60fc2e370c565dbc088769fdf4fa5d1` demonstrates the
relevant architecture: bounded packs, one reader per input stream, worker processing, and
separate compression/decompression budgets. It does not make one parser object concurrently
readable.

That comparison is an architectural clue, not a compatibility target. Any throughput claim
must name the input layout, compression format, projected columns, record lengths, thread
budgets, storage, and exact fastp commit.

## Transport boundary

Remote and browser support depends on htslib `hFILE`, including DuckHTS's wasm HTTP/HTTPS
backend. An optimized reader must therefore keep htslib transport below the FASTQ parser. A
local `FILE *` or custom decompressor that bypasses htslib is not an acceptable default path.

The preferred boundary is:

```text
htslib transport -> bounded FASTQ parser -> projected DuckDB vectors
```

## Next experiment

Implement a single-threaded, projection-aware FASTQ parser over `hts_getline(...)` before
designing a worker pipeline. It must preserve:

- multiline sequence and quality handling accepted by the current reader;
- paired-file and interleaved name/count validation;
- configured and auto-detected quality encodings;
- string and packed sequence/quality projections;
- zero-column count behavior;
- remote and wasm transport; and
- precise EOF, truncation, and allocation errors.

Benchmark this parser against the current `sam_read1(...)` path using the same transport and
materialized columns. If parsing remains dominant, the next design may add bounded record packs
with one read/decompress stage and worker parse/materialization stages. It must define ordered
emission, backpressure, error cancellation, and total thread ownership before code is added.

Do not add a public tuning parameter until it controls implemented behavior with a measured
workload and a stable distinction between htslib decompression workers and DuckHTS processing
workers.
