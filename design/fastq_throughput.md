# FASTQ throughput

Status: current implementation guidance and open performance investigation. The direct
single-stream parser is implemented; a parallel parser/QC pipeline remains an evidence-driven
choice, not a public promise. `functions.yaml` defines the current `read_fastq(...)` surface.

## Current data path

`read_fastq(...)` detects FASTQ with htslib, then reads lines through the same htslib `hFILE`
transport used for local, compressed, remote, and browser inputs. It parses each record directly
into projected DuckDB columns:

```text
htslib transport -> multiline FASTQ state machine -> projected DuckDB vectors
```

This removes the previous text -> temporary `bam1_t` -> text/list round trip. The parser retains
only the requested name, sequence, and quality bytes, while still counting every byte needed to
validate record framing. FASTA keeps htslib's SAM-sequence adapter.

For exact read/base/quality and per-cycle QC, the implemented SQL path is:

```text
htslib transport -> projected sequence/quality strings -> duckhts_fastq_qc aggregate
                 -> one compact result plus one small per-cycle relation
```

`duckhts_fastq_qc(sequence, quality [, max_cycles])` consumes canonical Phred+33 quality text
directly from DuckDB string vectors. It updates bounded aggregate state while the scan is live;
it does not create one SQL row or list element per base. Scalar is the semantic oracle, while
AVX2, ARM NEON, and wasm SIMD128 implementations classify blocks and retain the same scalar
per-cycle sufficient statistics. The aggregate returns integer counts and sums, leaving means,
fractions, plotting, and joins to SQL. One DuckDB string materialization remains between the
reader and aggregate; removing that copy would require an explicit fused table-function or
vector interface and must be justified by a measured residual cost.

The direct path preserves the established contract:

- multiline sequence and quality blocks;
- htslib query-name parsing, including terminal `/digit` removal;
- paired-file name/count checks and interleaved count checks;
- string or nt16 sequence output and string or packed-Phred quality output;
- configured or sampled quality encodings;
- indexed count shortcuts and explicit sequential scans; and
- DuckDB-visible truncation, length, allocation, and transport errors.

The parser remains one DuckDB scan worker. This is deliberate: htslib owns transport and any BGZF
decompression workers, while DuckDB owns query parallelism across independent files. Those thread
budgets must not be conflated.

## Measured floor

[`benchmark_fastq_reader.md`](../benchmarks/benchmark_fastq_reader.md) compares the direct parser
with the previous `sam_read1(...)` implementation on the same 2,000,000 uncompressed 150-base
records, one DuckDB thread, and one pinned physical-core sibling. It separately records current
fastp QC-only throughput; fastp performs base statistics and is context, not an equivalent SQL
materialization workload.

The same report includes a deliberately direct SQL formulation that expands every read into one
row per base with `generate_series`/`UNNEST`. On the pinned eight-file HG001 Garvan exome
(168,453,484 reads and 17,013,801,884 bases), that formulation took 1,236.185 seconds while
producing exact fastp global totals. It is retained as a negative control: DuckDB can express the
calculation, but per-base relational materialization is the wrong physical plan for this
reduction. The fused aggregate completed the same corpus in 216.609 seconds (5.71x faster), at
777,685 reads/s and 78.55 million bases/s with 169,588 KiB peak RSS. A final exact-artifact,
equal-priority rerun measured 215.149 seconds (782,962 reads/s; 79.08 million bases/s) against
180.580 seconds for current fastp, leaving DuckHTS 19.1% slower on this one-core gzip workload.
The report records why the paired rerun used elevated scheduler priority on the busy shared host.
A first-file reader-only text scan and the same scan plus fused QC differed by only 0.130 seconds;
compressed transport and parsing, not the compact aggregate, are the next measured target.

Every benchmark must name the input layout, compression, projected columns, input reads/bases,
materialized units, CPU affinity, thread settings, and exact source artifacts. A compressed or
remote result is not interchangeable with this local uncompressed parser measurement.

## Lessons from fastp and paraseq

The pinned mirrors are:

- fastp `d517536b021bca0916cf33cb456f4e4b8aa24456`: one reader feeds bounded 1,000-read packs to
  worker-owned queues; current QC primitives use Highway SIMD; reader, worker, and writer budgets
  are explicit.
- paraseq `c35320c9f2800b7bf680ec79297e4781b96695a9`: reusable bulk byte buffers hold borrowed record
  slices, newline offsets are found in one scan, and only an incomplete trailing record is copied.

Both current parsers assume four physical lines per FASTQ record. DuckHTS accepts multiline FASTQ,
so neither parser can replace its state machine unchanged. The reusable ideas are bounded packs,
bulk buffers with offsets, worker-owned scratch, and one-pass fused QC—not their file interfaces.

## Next measured decision

The next experiment, if profiling still finds parsing or QC dominant, is a bounded ordinary-
four-line fast lane with state-machine fallback for multiline records. It must prove identical
rows, stable input order, bounded memory, prompt cancellation, and no oversubscription before a
worker count becomes public. The implemented QC aggregate already fuses base counts, quality
thresholds, and per-cycle sufficient statistics. Additional length histograms and
adapter/overrepresentation sketches should extend or consume compact batch state rather than
introduce independent full-string rescans.

Keep htslib below the parser. A local `FILE *`, ISA-L-only reader, or custom remote transport would
discard DuckHTS's format, URL, CRAN, and wasm guarantees and is not an acceptable default path.
