Insert-size counts: riker 0.4.1, native and projected SQL
================

## Contract and declared targets

Fulcrum Genomics **riker 0.4.1** (MIT), tag
`40d5e3260c5b1b9d210c28f7c90a75387013b4ce`, is the external oracle.
Picard is context only. Integer counts must match exactly; mean,
population standard deviation, median and MAD must match riker’s
serialized two decimal places exactly. Both directions of each relation
comparison are checked.

Targets are declared before any timed run:

- Every BAM/CRAM, arm, thread count and scale must satisfy exact parity.
- At 1×, each DuckHTS arm’s median elapsed time divided by riker’s must
  be **\<= 1.05 at 1 thread, \<= 1.00 at 2, and \<= 0.90 at 4 and 8**,
  independently for BAM and CRAM. Apply these targets to native and
  projected SQL separately.
- Report the native/projected-SQL elapsed-time ratio without a pass/fail
  gate; this comparison informs the supervisor’s final public-name
  decision.
- For each arm/format/thread count whose 1× median is at least 5
  seconds, the 1×/2×/4× log-time/log-input exponent must be **\<=
  1.25**.
- Peak process RSS must be **\<= 256 MiB overhead + 3× declared live
  state**. Live state includes worst-case overlapping per-worker sparse
  hash tables (16.25 bytes/slot, power-of-two capacity at \<= .77 load),
  8 MiB decoder state per handle, one merged hash table and 24 bytes per
  output count bin. This is a hypothesis to test, not an assertion about
  measured htslib usage.
- No spill is allowed. DuckDB has a 2 GiB buffer limit and zero
  temporary storage allowance; a worker error is a failed run, not a
  discarded trial.

Full workloads use registry artifacts `roh_mixture_hg00403_chr20_cram`
and `insert_size_hg00403_chr20_bam`; the latter is a network-free BAM
transcode of that CRAM using `liftover_grch38_fasta` and RBCFTools’
bundled samtools. Stage with `benchmarks/insert_size_stage.R` after
checking disk space; the function also checks the 8 GiB free-space floor
with a 6 GiB staging reserve. No full HG00188 BAM is staged. The smoke
workload and SQL reference are checksum-pinned committed fixtures staged
through `insert-size-smoke` in the same registry. The SQL arm is the
public `duckhts_bam_insert_size_counts` macro over packed CIGAR,
existing flag predicates and `cigar_reference_length`; the native arm is
`duckhts_bam_insert_size_counts_native`.

Each arm runs in a fresh R process (a fresh DuckDB process for
native/SQL), with thread settings 1/2/4/8. Native settings are histogram
workers, projected SQL settings are DuckDB processing threads, and riker
settings are its own `--threads` semantics; these are not identical
worker models. Native/SQL use zero additional htslib decompression
workers. The full matrix has three replicates and 1×/2×/4× complete file
scans, with output materialization on each scan; scales do not change
bin cardinality. This measures scaling in records, not scaling in
distinct insert sizes. `/usr/bin/time` measures end-to-end worker wall
time and maximum RSS, including R startup, input decoding, reduction,
metrics and output materialization. riker also writes its plot; DuckHTS
emits counts/metrics/histogram relations, not a plot. This output-cost
difference must accompany any speed comparison.

Independent riker raw and default outputs are prepared outside measured
workers. All worker results are compared against them; no failed
comparison is filtered. Record denominators come from registry
identities, retained pair and bin denominators from riker’s untrimmed
histogram. Summary moments trim tails without changing original pair
counts. Smoke runs invoke all 24 format/arm/thread combinations with
**no clock or RSS measurement**.

## Revision and runtime

Source revision: `4e1cc27fd516ebcfaab28ca4ef1102da40c9acd4`.

Oracle: `riker 0.4.1`; SHA-256
`806e0c05ce5993c7fac416cf1e237a13557ca7e7e419832002a31de91ca519e6`.

Runtime: R version 4.6.0 (2026-04-24) DuckDB 1.5.5.

CPU affinity and niceness are supplied by the launcher; this smoke uses
cores 16–19, nice 10.

## Results

**Reduced smoke only; no timed comparison was run.**

24 fresh-process checks passed; exact raw integer counts and rounded
metrics/histograms agree.

| format | arm           | threads | records | read_pairs | bins | metric_rows | histogram_rows | parity |
|:-------|:--------------|--------:|--------:|-----------:|-----:|------------:|---------------:|:-------|
| BAM    | riker         |       1 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | native        |       1 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | projected_sql |       1 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | riker         |       2 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | native        |       2 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | projected_sql |       2 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | riker         |       4 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | native        |       4 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | projected_sql |       4 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | riker         |       8 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | native        |       8 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | projected_sql |       8 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| CRAM   | riker         |       1 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | native        |       1 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | projected_sql |       1 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | riker         |       2 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | native        |       2 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | projected_sql |       2 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | riker         |       4 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | native        |       4 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | projected_sql |       4 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | riker         |       8 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | native        |       8 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | projected_sql |       8 |      26 |         16 |    9 |           3 |              5 | TRUE   |

Speed, RSS and scaling criteria: **not measured**. Spill: disabled for
native/SQL.

## Fixture footprint

The corpus retains all eight alignment cases and all seven riker
argument cases in keyed `metrics.tsv`/`histograms.tsv` files, including
exact upstream decimal text. BAM and CRAM have separate oracle keys. The
R package ships only edges/ empty wrapper fixtures and their oracle
subset, not the full corpus. File payload sizes exclude directory
allocation:

| location   | files | bytes |
|:-----------|------:|------:|
| repository |    42 | 33667 |
| package    |    17 | 15594 |

## Interpretation and next measurement

The native reducer is expected to beat projected SQL because it avoids
emitting one DuckDB row and packed CIGAR per alignment, SQL
flag/orientation expressions, and the per-record hash aggregate. This is
an expectation, not measured evidence. No building-block bottleneck has
been identified with a profile, so no block optimization is justified
here.

Known CRAM cost: `src/bam_reader.c:bam_cram_required_fields` adds
`SAM_SEQ`, `SAM_QUAL`, `SAM_AUX` and `SAM_RGAUX` when CIGAR is
projected, for long-CIGAR CG restoration. Native requests CIGAR but
omits SEQ, QUAL and AUX. The extra decoding’s wall-time and RSS effects
are unmeasured. Profile CRAM `sam_read1`/ decompression separately from
SQL flag predicates, `cigar_reference_length` and packed-CIGAR vector
materialization before attributing an arm difference.

If that decoding dominates, investigate a reader-projection change
separately: htslib 1.24 `sam.c:682-706` (`bam_tag2cigar`) tests the fake
soft-clip length against `l_qseq` before consulting CG;
`cram/cram_decode.c:3231-3233` calls it after constructing the aux
block. Any narrower request must preserve read length and CG
restoration, retaining the long-CIGAR BAM/CRAM fixture. No reader
projection change is made here. The native/SQL ratio informs the
supervisor’s naming decision; the unsuffixed name remains SQL pending
timing. The nearest recorded BAM-reduction baseline is
`benchmarks/benchmark_bam_mismatch_counts.md`; it exercises a different
kernel and cannot establish insert-size throughput or regression status.

Profile scan/decode time, CRAM reference fetches, task-window overlap
decoding, histogram hash operations and merge time separately. At
multiple threads, inspect storage throughput and per-handle
decompression before changing task sizes or adding SIMD. Compare an
explicitly sequential native scan against indexed scans to quantify
index and repeated-overlap costs. Very short inputs will be dominated by
fresh-process startup and riker plotting. The full chr20 matrix,
1×/2×/4× measurements and RSS ceilings remain for the quiet-host run;
none can be inferred from this smoke.
