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
  gate; native owns the unsuffixed public name by the recorded
  comparison.
- For each arm/format/thread count whose 1× median is at least 5
  seconds, each adjacent doubling’s log-time/log-input exponent must be
  **\<= 1.25**.
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
through `insert-size-smoke` in the same registry. The
`insert-size-scale` registry workload names 2×/4× BAM and CRAM inputs.
`stage_insert_size_scale_from_registry()` coordinate-merges two or four
copies of the corresponding 1× input with samtools
`merge -c -p --no-PG`, preserving record flags, read names and RG/PG
identifiers, then indexes the result. It checks alignment counts,
position sums and read-length sums against the registry before
publication. Receipts bind actual input/output hashes and the samtools
version. Scale staging requires 8 GiB free plus 1.2 times the expanded
input bytes. Staging is outside the timed runs. The SQL arm is
`duckhts_bam_insert_size_counts_sql`, a macro over packed CIGAR,
existing flag predicates and `cigar_reference_length`; the native arm is
`duckhts_bam_insert_size_counts`. Both arms remain public.

## Measurement amendment: direct tool processes

This amendment replaces the timing method, not the targets. The first
timed race at commit `676f3a8d` is retained below as **R-worker
harness**, including its original results and verdicts. Those timings
are not pooled with direct-tool timings or corrected by subtraction.

The supervisor measured startup on CPUs 12–15 on 2026-10-14:

| Process                     |       Wall seconds | Peak RSS (MB, as supplied) |
|-----------------------------|-------------------:|---------------------------:|
| `Rscript -e NULL`           |               0.14 |                         63 |
| R + duckdb-R connect        |               0.44 |                        103 |
| R + duckdb-R + LOAD duckhts |               0.51 |                        119 |
| DuckDB CLI + LOAD duckhts   |               0.08 |                         41 |
| riker binary                | approximately 0.00 |                          6 |

Every old arm paid R startup, but only DuckHTS paid duckdb-R package
load and connection setup. The approximately 0.37-second differential
confounds the BAM reducer comparison. Subtracting it estimates
native/riker ratios of 1.00/0.61/0.86/0.70 at 1/2/4/8 threads; these
estimates are diagnostic context, not measurements under this amendment
and not replacement verdicts.

`/usr/bin/time` now wraps tool processes directly, with no R inside the
timed process. DuckHTS runs one `duckdb -unsigned` process and SQL
script per run: LOAD, SET threads, a 2 GiB memory limit, zero spill
allowance, then one temporary counts table, metrics macro and histogram.
COPY writes counts, metrics and histogram TSVs. Every scale uses one
public invocation on one staged file of the corresponding size, in a
fresh tool process for both riker and DuckDB. riker writes metrics,
histogram and its plot. DuckHTS does not write a plot and additionally
writes raw counts; retain these output-cost differences in comparisons.

The CLI version is recorded with the extension’s stable C API target and
header version. A version difference is stated explicitly, including
against the old duckdb-R 1.5.5 runtime; results across runtime versions
do not isolate the harness effect. The available CLI is v1.5.1, the
stable API target is v1.2.0, and build headers are v1.5.3. Loading the
extension is checked before the matrix.

A second, reported-only DuckDB duration uses SQL timestamps after
LOAD/SET and before the first counts query through the last output COPY.
It excludes process startup and extension loading but includes
input/index opening, reduction, summary and serialization. It is not an
acceptance criterion. Parity is checked after timing in R by reading the
output TSVs; no R deserialization or comparison contributes to tool wall
time or RSS. Targets, exact parity rules, scales (1×/2×/4×), thread
settings (1/2/4/8), three replicates and RSS gates are unchanged.
Scaling exponents are checked at both adjacent doublings, as required by
`STYLE.md`. Smoke is untimed. For timed DuckDB runs, JSON profiles of
both temporary-table queries and all three output queries record peak
buffer memory and spill; profiler setup and file-writing costs are
included in elapsed time. Peak RSS remains the process memory gate.
Completed rows and a failed run’s files are retained if execution or
parity fails. Raw parity compares every actual sparse count row against
riker’s positive orientation/size cells by key; duplicate, missing,
extra and changed rows fail. The comparison runs outside measured
processes and does not copy a whole histogram for each individual count
row.

Each arm runs in a fresh tool process, with thread settings 1/2/4/8.
Native settings are histogram workers, projected SQL settings are DuckDB
processing threads, and riker settings are its own `--threads`
semantics; these are not identical worker models. Native/SQL use zero
additional htslib decompression workers. The full matrix has three
replicates and one scan of each 1×/2×/4× file, with output
materialization; scales do not change bin cardinality. This measures
scaling in records, not scaling in distinct insert sizes.
`/usr/bin/time` measures end-to-end tool wall time and maximum RSS,
including tool startup, input decoding, reduction, metrics and output
files, but excluding the R controller.

Independent riker raw and default outputs are prepared outside measured
tool processes. All tool results are compared against them; no failed
comparison is filtered. Record denominators come from registry
identities, retained pair and bin denominators from riker’s untrimmed
histogram. Summary moments trim tails without changing original pair
counts. Smoke runs invoke all 24 format/arm/thread combinations with
**no clock or RSS measurement**.

## Revision and runtime

Source revision: `5abf84a5f0624c92d5b19169b8855a59140ac0d6`; working
tree changed: `TRUE`; driver SHA-256:
`bdb649ff5acdb4c7d85ad886dfaa467b85b436bce3c73f0f5fd474cee9f20f28`.

Oracle: `riker 0.4.1`; SHA-256
`806e0c05ce5993c7fac416cf1e237a13557ca7e7e419832002a31de91ca519e6`.

Runtime: R controller: R version 4.6.0 (2026-04-24) ; timed DuckDB CLI:
v1.5.1 (Variegata) 7dbb2e646f.

CLI executable: `/root/.duckdb/cli/1.5.1/duckdb`. CLI version:
`v1.5.1 (Variegata) 7dbb2e646f`. Repository build defaults: stable C API
target `v1.2.0`, headers `v1.5.3`.

The CLI differs from those versions and from the retained R-worker
runtime (DuckDB 1.5.5). Successful LOAD verifies compatibility, not
identical performance across releases.

Launcher affinity: `pid 3's current affinity list: 0-19`; niceness: `0`.

Host: `Linux; 6.8.0-78-generic; x86_64`; load at start:
`23:20:15 up 373 days,  8:14, 19 users,  load average: 4.04, 3.61, 3.36`.

Cache condition: warm; each input is scanned for independent oracle
outputs before timing. Page cache is not flushed. Concurrent work
declared by launcher: none declared.

## Results

**Reduced smoke only; no timed comparison was run.**

24 fresh-process checks passed; exact raw integer counts and rounded
metrics/histograms agree.

| format | arm                                | threads | records | read_pairs | bins | metric_rows | histogram_rows | parity |
|:-------|:-----------------------------------|--------:|--------:|-----------:|-----:|------------:|---------------:|:-------|
| BAM    | riker                              |       1 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       1 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       1 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | riker                              |       2 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       2 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       2 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | riker                              |       4 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       4 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       4 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | riker                              |       8 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       8 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       8 |      26 |         15 |    9 |           3 |              5 | TRUE   |
| CRAM   | riker                              |       1 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       1 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       1 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | riker                              |       2 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       2 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       2 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | riker                              |       4 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       4 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       4 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | riker                              |       8 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       8 |      26 |         16 |    9 |           3 |              5 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       8 |      26 |         16 |    9 |           3 |              5 | TRUE   |

Speed, RSS and scaling criteria: **not measured**. Spill: disabled for
native/SQL.

## Reported-only in-process context

`work_seconds` excludes CLI startup and LOAD/SET, retaining the query,
profiling and TSV output. The remainder includes startup/setup and
end-of-process work; it is not an isolated startup measurement or a
pass/fail criterion.

Untimed smoke: no process or in-process durations recorded.

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

## Timing verdicts and scale context

The declared timing criterion is the **1× median ratio** against riker,
per format, arm and thread count. The 4× and marginal ratios are
reported context, not alternative acceptance criteria. Marginal cost is
`(4× median - 1× median) / 3` for each arm; its ratio divides by riker’s
corresponding marginal cost.

| format | arm                                | threads | ratio_1x | target | verdict      | ratio_4x | marginal_ratio |
|:-------|:-----------------------------------|--------:|---------:|-------:|:-------------|---------:|---------------:|
| BAM    | duckhts_bam_insert_size_counts     |       1 |       NA |   1.05 | not measured |       NA |             NA |
| BAM    | duckhts_bam_insert_size_counts_sql |       1 |       NA |   1.05 | not measured |       NA |             NA |
| BAM    | duckhts_bam_insert_size_counts     |       2 |       NA |   1.00 | not measured |       NA |             NA |
| BAM    | duckhts_bam_insert_size_counts_sql |       2 |       NA |   1.00 | not measured |       NA |             NA |
| BAM    | duckhts_bam_insert_size_counts     |       4 |       NA |   0.90 | not measured |       NA |             NA |
| BAM    | duckhts_bam_insert_size_counts_sql |       4 |       NA |   0.90 | not measured |       NA |             NA |
| BAM    | duckhts_bam_insert_size_counts     |       8 |       NA |   0.90 | not measured |       NA |             NA |
| BAM    | duckhts_bam_insert_size_counts_sql |       8 |       NA |   0.90 | not measured |       NA |             NA |
| CRAM   | duckhts_bam_insert_size_counts     |       1 |       NA |   1.05 | not measured |       NA |             NA |
| CRAM   | duckhts_bam_insert_size_counts_sql |       1 |       NA |   1.05 | not measured |       NA |             NA |
| CRAM   | duckhts_bam_insert_size_counts     |       2 |       NA |   1.00 | not measured |       NA |             NA |
| CRAM   | duckhts_bam_insert_size_counts_sql |       2 |       NA |   1.00 | not measured |       NA |             NA |
| CRAM   | duckhts_bam_insert_size_counts     |       4 |       NA |   0.90 | not measured |       NA |             NA |
| CRAM   | duckhts_bam_insert_size_counts_sql |       4 |       NA |   0.90 | not measured |       NA |             NA |
| CRAM   | duckhts_bam_insert_size_counts     |       8 |       NA |   0.90 | not measured |       NA |             NA |
| CRAM   | duckhts_bam_insert_size_counts_sql |       8 |       NA |   0.90 | not measured |       NA |             NA |

Native/SQL comparison (reported, not a pass/fail gate):

| format | threads | native_seconds_1x | sql_seconds_1x | native_sql_ratio_1x |
|:-------|--------:|------------------:|---------------:|--------------------:|
| BAM    |       1 |                NA |             NA |                  NA |
| BAM    |       2 |                NA |             NA |                  NA |
| BAM    |       4 |                NA |             NA |                  NA |
| BAM    |       8 |                NA |             NA |                  NA |
| CRAM   |       1 |                NA |             NA |                  NA |
| CRAM   |       2 |                NA |             NA |                  NA |
| CRAM   |       4 |                NA |             NA |                  NA |
| CRAM   |       8 |                NA |             NA |                  NA |

## SQL thread scaling and the next building block

| format | threads | seconds | speedup_vs_1_thread |
|:-------|--------:|--------:|--------------------:|
| BAM    |       1 |      NA |                  NA |
| BAM    |       2 |      NA |                  NA |
| BAM    |       4 |      NA |                  NA |
| BAM    |       8 |      NA |                  NA |
| CRAM   |       1 |      NA |                  NA |
| CRAM   |       2 |      NA |                  NA |
| CRAM   |       4 |      NA |                  NA |
| CRAM   |       8 |      NA |                  NA |

The [`read_bam` catalog,
“Scanning”](../r/Rduckhts/inst/function_catalog/reference.md#read_bam)
specifies: “An indexed full-file scan is split by reference across
DuckDB threads”; “A CRAM scan is split by reference”. The full chr20
workload has one reference, so that split supplies only one scan worker
regardless of the SQL thread setting. Within-reference splitting of
`read_bam` is the next building block to optimise. It is not implemented
here. This scan constraint does not identify the relative costs of flag
predicates, packed-CIGAR materialization or `cigar_reference_length`;
those need a profile.

## CRAM projection context

The comparison table reports native and SQL wall times separately. SQL
projects CIGAR; `src/bam_reader.c:bam_cram_required_fields` consequently
requests `SAM_SEQ`, `SAM_QUAL`, `SAM_AUX` and `SAM_RGAUX` for long-CIGAR
CG restoration. Native requests CIGAR without SEQ/QUAL/AUX. This is a
known difference in work, not a measured attribution of the CRAM timing
gap. Its isolated cost has not been measured. Direct-tool timings also
include raw counts TSV output, unlike riker’s default output, and riker
alone writes a plot.

A narrower reader projection must preserve read length and CG
restoration: htslib 1.24 `sam.c:682-706` (`bam_tag2cigar`) tests the
fake soft-clip length against `l_qseq` before consulting CG;
`cram/cram_decode.c:3231-3233` calls it after constructing the aux
block. Retain the long-CIGAR BAM/CRAM fixture when investigating that
projection separately.

## Startup context

The amendment’s startup measurements identify an asymmetric R-worker
cost. The direct-tool table and reported-only in-process duration keep
tool startup visible without replacing the declared 1× criterion.
Neither 4× nor marginal ratios replace that criterion. Scan/decode,
task-window overlap decoding, histogram updates and merge costs still
require profiling. Smoke results provide no timing or memory evidence.

## Retained R-worker harness

The following is the first timed report, preserved byte-for-byte from
`676f3a8d` in `benchmark_insert_size_counts_r_worker.md`. Its arm names
and verdicts describe that revision and that harness, not the amended
method.

# Insert-size counts: riker 0.4.1, native and projected SQL

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

Source revision: `c010c55e145445e4080247eadc470244567aca6b`.

Oracle: `riker 0.4.1`; SHA-256
`806e0c05ce5993c7fac416cf1e237a13557ca7e7e419832002a31de91ca519e6`.

Runtime: R version 4.6.0 (2026-04-24) DuckDB 1.5.5.

CPU affinity and niceness are supplied by the launcher; this smoke uses
cores 16–19, nice 10.

## Results

| format | arm           | threads | scale | replicate |  records | read_pairs |  bins | metric_rows | histogram_rows | seconds | peak_rss_mib | rss_budget_mib | parity |
|:-------|:--------------|--------:|------:|----------:|---------:|-----------:|------:|------------:|---------------:|--------:|-------------:|---------------:|:-------|
| BAM    | riker         |       1 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    4.08 |     75.92969 |       288.3180 | TRUE   |
| BAM    | riker         |       1 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    4.08 |     75.77734 |       288.3180 | TRUE   |
| BAM    | riker         |       1 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    4.10 |     75.77344 |       288.3180 | TRUE   |
| BAM    | native        |       1 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    4.47 |    143.71094 |       288.3180 | TRUE   |
| BAM    | native        |       1 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    4.46 |    144.03125 |       288.3180 | TRUE   |
| BAM    | native        |       1 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    4.51 |    143.98828 |       288.3180 | TRUE   |
| BAM    | projected_sql |       1 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.35 |    144.33594 |       288.3180 | TRUE   |
| BAM    | projected_sql |       1 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.32 |    144.28516 |       288.3180 | TRUE   |
| BAM    | projected_sql |       1 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.25 |    144.61719 |       288.3180 | TRUE   |
| BAM    | riker         |       1 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    8.02 |     75.77344 |       288.3180 | TRUE   |
| BAM    | riker         |       1 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    8.07 |     75.77344 |       288.3180 | TRUE   |
| BAM    | riker         |       1 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    8.00 |     75.77344 |       288.3180 | TRUE   |
| BAM    | native        |       1 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    8.53 |    145.71875 |       288.3180 | TRUE   |
| BAM    | native        |       1 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    8.55 |    146.32812 |       288.3180 | TRUE   |
| BAM    | native        |       1 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    8.55 |    146.37109 |       288.3180 | TRUE   |
| BAM    | projected_sql |       1 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.42 |    146.21484 |       288.3180 | TRUE   |
| BAM    | projected_sql |       1 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.30 |    146.20312 |       288.3180 | TRUE   |
| BAM    | projected_sql |       1 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.43 |    146.21094 |       288.3180 | TRUE   |
| BAM    | riker         |       1 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   15.83 |     75.92969 |       288.3180 | TRUE   |
| BAM    | riker         |       1 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   15.86 |     75.77344 |       288.3180 | TRUE   |
| BAM    | riker         |       1 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   15.84 |     75.92969 |       288.3180 | TRUE   |
| BAM    | native        |       1 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   16.72 |    150.07812 |       288.3180 | TRUE   |
| BAM    | native        |       1 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   16.61 |    150.16406 |       288.3180 | TRUE   |
| BAM    | native        |       1 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   16.61 |    150.18750 |       288.3180 | TRUE   |
| BAM    | projected_sql |       1 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   28.23 |    150.26562 |       288.3180 | TRUE   |
| BAM    | projected_sql |       1 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   27.88 |    150.60938 |       288.3180 | TRUE   |
| BAM    | projected_sql |       1 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   28.06 |    150.28125 |       288.3180 | TRUE   |
| BAM    | riker         |       2 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    3.50 |     75.77344 |       315.3649 | TRUE   |
| BAM    | riker         |       2 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    3.53 |     75.77344 |       315.3649 | TRUE   |
| BAM    | riker         |       2 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    3.54 |     75.77344 |       315.3649 | TRUE   |
| BAM    | native        |       2 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.52 |    155.29688 |       315.3649 | TRUE   |
| BAM    | native        |       2 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.51 |    155.26562 |       315.3649 | TRUE   |
| BAM    | native        |       2 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.53 |    155.23047 |       315.3649 | TRUE   |
| BAM    | projected_sql |       2 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.28 |    150.68359 |       315.3649 | TRUE   |
| BAM    | projected_sql |       2 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.51 |    150.67969 |       315.3649 | TRUE   |
| BAM    | projected_sql |       2 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.38 |    150.10547 |       315.3649 | TRUE   |
| BAM    | riker         |       2 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    6.90 |     75.92969 |       315.3649 | TRUE   |
| BAM    | riker         |       2 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    6.93 |     75.76953 |       315.3649 | TRUE   |
| BAM    | riker         |       2 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    6.89 |     75.92969 |       315.3649 | TRUE   |
| BAM    | native        |       2 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    4.60 |    157.88281 |       315.3649 | TRUE   |
| BAM    | native        |       2 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    4.64 |    159.26562 |       315.3649 | TRUE   |
| BAM    | native        |       2 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    4.60 |    157.92188 |       315.3649 | TRUE   |
| BAM    | projected_sql |       2 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.79 |    155.18359 |       315.3649 | TRUE   |
| BAM    | projected_sql |       2 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.62 |    153.25000 |       315.3649 | TRUE   |
| BAM    | projected_sql |       2 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.33 |    153.04688 |       315.3649 | TRUE   |
| BAM    | riker         |       2 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   13.72 |     75.77734 |       315.3649 | TRUE   |
| BAM    | riker         |       2 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   13.65 |     75.92969 |       315.3649 | TRUE   |
| BAM    | riker         |       2 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   13.69 |     75.92969 |       315.3649 | TRUE   |
| BAM    | native        |       2 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    8.81 |    163.07031 |       315.3649 | TRUE   |
| BAM    | native        |       2 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    8.79 |    163.45703 |       315.3649 | TRUE   |
| BAM    | native        |       2 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    8.79 |    162.40625 |       315.3649 | TRUE   |
| BAM    | projected_sql |       2 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   28.01 |    157.31250 |       315.3649 | TRUE   |
| BAM    | projected_sql |       2 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   27.91 |    157.64844 |       315.3649 | TRUE   |
| BAM    | projected_sql |       2 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   28.35 |    157.31641 |       315.3649 | TRUE   |
| BAM    | riker         |       4 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.33 |     75.46094 |       369.4586 | TRUE   |
| BAM    | riker         |       4 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.32 |     75.77344 |       369.4586 | TRUE   |
| BAM    | riker         |       4 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.31 |     75.92969 |       369.4586 | TRUE   |
| BAM    | native        |       4 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.51 |    164.27344 |       369.4586 | TRUE   |
| BAM    | native        |       4 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.52 |    163.64844 |       369.4586 | TRUE   |
| BAM    | native        |       4 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.51 |    165.04688 |       369.4586 | TRUE   |
| BAM    | projected_sql |       4 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.23 |    155.27344 |       369.4586 | TRUE   |
| BAM    | projected_sql |       4 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.31 |    156.14062 |       369.4586 | TRUE   |
| BAM    | projected_sql |       4 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.48 |    156.91797 |       369.4586 | TRUE   |
| BAM    | riker         |       4 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.50 |     75.77344 |       369.4586 | TRUE   |
| BAM    | riker         |       4 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.49 |     75.93359 |       369.4586 | TRUE   |
| BAM    | riker         |       4 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.48 |     75.92969 |       369.4586 | TRUE   |
| BAM    | native        |       4 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.62 |    168.19922 |       369.4586 | TRUE   |
| BAM    | native        |       4 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.63 |    168.21875 |       369.4586 | TRUE   |
| BAM    | native        |       4 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.63 |    167.60156 |       369.4586 | TRUE   |
| BAM    | projected_sql |       4 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.22 |    159.59375 |       369.4586 | TRUE   |
| BAM    | projected_sql |       4 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.42 |    159.40234 |       369.4586 | TRUE   |
| BAM    | projected_sql |       4 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.53 |    160.38281 |       369.4586 | TRUE   |
| BAM    | riker         |       4 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    4.83 |     75.77344 |       369.4586 | TRUE   |
| BAM    | riker         |       4 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    4.83 |     75.92969 |       369.4586 | TRUE   |
| BAM    | riker         |       4 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    4.83 |     75.77734 |       369.4586 | TRUE   |
| BAM    | native        |       4 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    4.83 |    173.01562 |       369.4586 | TRUE   |
| BAM    | native        |       4 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    4.86 |    172.46875 |       369.4586 | TRUE   |
| BAM    | native        |       4 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    4.82 |    171.56250 |       369.4586 | TRUE   |
| BAM    | projected_sql |       4 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   28.37 |    166.11328 |       369.4586 | TRUE   |
| BAM    | projected_sql |       4 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   28.18 |    164.99219 |       369.4586 | TRUE   |
| BAM    | projected_sql |       4 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   28.41 |    167.54297 |       369.4586 | TRUE   |
| BAM    | riker         |       8 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.03 |     75.61719 |       477.6461 | TRUE   |
| BAM    | riker         |       8 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.05 |     75.92969 |       477.6461 | TRUE   |
| BAM    | riker         |       8 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.11 |     75.77344 |       477.6461 | TRUE   |
| BAM    | native        |       8 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.11 |    185.80859 |       477.6461 | TRUE   |
| BAM    | native        |       8 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.11 |    182.28516 |       477.6461 | TRUE   |
| BAM    | native        |       8 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.10 |    179.69141 |       477.6461 | TRUE   |
| BAM    | projected_sql |       8 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.66 |    171.36719 |       477.6461 | TRUE   |
| BAM    | projected_sql |       8 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.47 |    170.91797 |       477.6461 | TRUE   |
| BAM    | projected_sql |       8 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.37 |    170.69922 |       477.6461 | TRUE   |
| BAM    | riker         |       8 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.00 |     75.77344 |       477.6461 | TRUE   |
| BAM    | riker         |       8 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.97 |     75.65625 |       477.6461 | TRUE   |
| BAM    | riker         |       8 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.14 |     75.77734 |       477.6461 | TRUE   |
| BAM    | native        |       8 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.87 |    190.05859 |       477.6461 | TRUE   |
| BAM    | native        |       8 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.83 |    187.64062 |       477.6461 | TRUE   |
| BAM    | native        |       8 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.82 |    188.55078 |       477.6461 | TRUE   |
| BAM    | projected_sql |       8 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.52 |    175.85547 |       477.6461 | TRUE   |
| BAM    | projected_sql |       8 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.27 |    175.90625 |       477.6461 | TRUE   |
| BAM    | projected_sql |       8 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.29 |    175.24609 |       477.6461 | TRUE   |
| BAM    | riker         |       8 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    3.78 |     75.93359 |       477.6461 | TRUE   |
| BAM    | riker         |       8 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    3.77 |     75.92969 |       477.6461 | TRUE   |
| BAM    | riker         |       8 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    3.64 |     75.77344 |       477.6461 | TRUE   |
| BAM    | native        |       8 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    3.59 |    196.64453 |       477.6461 | TRUE   |
| BAM    | native        |       8 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    3.26 |    195.65625 |       477.6461 | TRUE   |
| BAM    | native        |       8 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    3.20 |    195.00391 |       477.6461 | TRUE   |
| BAM    | projected_sql |       8 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   29.07 |    182.67578 |       477.6461 | TRUE   |
| BAM    | projected_sql |       8 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   28.22 |    180.26172 |       477.6461 | TRUE   |
| BAM    | projected_sql |       8 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   28.32 |    181.94141 |       477.6461 | TRUE   |
| CRAM   | riker         |       1 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    8.32 |     75.93359 |       288.3180 | TRUE   |
| CRAM   | riker         |       1 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    8.44 |     75.77344 |       288.3180 | TRUE   |
| CRAM   | riker         |       1 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    8.32 |     75.77344 |       288.3180 | TRUE   |
| CRAM   | native        |       1 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.62 |    145.09766 |       288.3180 | TRUE   |
| CRAM   | native        |       1 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.64 |    145.05078 |       288.3180 | TRUE   |
| CRAM   | native        |       1 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.62 |    145.25781 |       288.3180 | TRUE   |
| CRAM   | projected_sql |       1 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |   12.06 |    154.85938 |       288.3180 | TRUE   |
| CRAM   | projected_sql |       1 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |   11.99 |    154.38281 |       288.3180 | TRUE   |
| CRAM   | projected_sql |       1 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |   11.99 |    154.36328 |       288.3180 | TRUE   |
| CRAM   | riker         |       1 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   16.48 |     75.77734 |       288.3180 | TRUE   |
| CRAM   | riker         |       1 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   16.46 |     75.77734 |       288.3180 | TRUE   |
| CRAM   | riker         |       1 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   16.53 |     75.92969 |       288.3180 | TRUE   |
| CRAM   | native        |       1 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    4.89 |    149.26562 |       288.3180 | TRUE   |
| CRAM   | native        |       1 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    4.92 |    148.27344 |       288.3180 | TRUE   |
| CRAM   | native        |       1 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    4.89 |    150.60156 |       288.3180 | TRUE   |
| CRAM   | projected_sql |       1 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.65 |    162.61328 |       288.3180 | TRUE   |
| CRAM   | projected_sql |       1 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.99 |    160.37500 |       288.3180 | TRUE   |
| CRAM   | projected_sql |       1 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.71 |    161.32812 |       288.3180 | TRUE   |
| CRAM   | riker         |       1 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   32.86 |     75.92969 |       288.3180 | TRUE   |
| CRAM   | riker         |       1 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   32.82 |     75.77344 |       288.3180 | TRUE   |
| CRAM   | riker         |       1 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   32.89 |     75.77344 |       288.3180 | TRUE   |
| CRAM   | native        |       1 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    9.27 |    154.02734 |       288.3180 | TRUE   |
| CRAM   | native        |       1 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    9.25 |    152.32812 |       288.3180 | TRUE   |
| CRAM   | native        |       1 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    9.24 |    152.48438 |       288.3180 | TRUE   |
| CRAM   | projected_sql |       1 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   47.57 |    164.27734 |       288.3180 | TRUE   |
| CRAM   | projected_sql |       1 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   46.86 |    164.81641 |       288.3180 | TRUE   |
| CRAM   | projected_sql |       1 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   46.87 |    165.61328 |       288.3180 | TRUE   |
| CRAM   | riker         |       2 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.33 |    113.25000 |       315.3649 | TRUE   |
| CRAM   | riker         |       2 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.30 |    111.29297 |       315.3649 | TRUE   |
| CRAM   | riker         |       2 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.32 |    113.12891 |       315.3649 | TRUE   |
| CRAM   | native        |       2 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.54 |    150.85547 |       315.3649 | TRUE   |
| CRAM   | native        |       2 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.58 |    152.13672 |       315.3649 | TRUE   |
| CRAM   | native        |       2 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.56 |    155.40625 |       315.3649 | TRUE   |
| CRAM   | projected_sql |       2 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |   12.02 |    156.41406 |       315.3649 | TRUE   |
| CRAM   | projected_sql |       2 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |   12.04 |    155.58984 |       315.3649 | TRUE   |
| CRAM   | projected_sql |       2 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |   12.10 |    164.34766 |       315.3649 | TRUE   |
| CRAM   | riker         |       2 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.62 |    112.88281 |       315.3649 | TRUE   |
| CRAM   | riker         |       2 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.55 |    113.13672 |       315.3649 | TRUE   |
| CRAM   | riker         |       2 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.75 |    111.70312 |       315.3649 | TRUE   |
| CRAM   | native        |       2 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.75 |    178.14844 |       315.3649 | TRUE   |
| CRAM   | native        |       2 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.71 |    178.01953 |       315.3649 | TRUE   |
| CRAM   | native        |       2 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.74 |    177.66797 |       315.3649 | TRUE   |
| CRAM   | projected_sql |       2 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.68 |    165.11719 |       315.3649 | TRUE   |
| CRAM   | projected_sql |       2 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.80 |    165.44141 |       315.3649 | TRUE   |
| CRAM   | projected_sql |       2 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.57 |    165.66016 |       315.3649 | TRUE   |
| CRAM   | riker         |       2 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   29.20 |    112.74219 |       315.3649 | TRUE   |
| CRAM   | riker         |       2 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   29.16 |    113.25391 |       315.3649 | TRUE   |
| CRAM   | riker         |       2 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   29.07 |    111.80859 |       315.3649 | TRUE   |
| CRAM   | native        |       2 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    5.10 |    185.21094 |       315.3649 | TRUE   |
| CRAM   | native        |       2 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    5.07 |    185.89844 |       315.3649 | TRUE   |
| CRAM   | native        |       2 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    5.13 |    184.08984 |       315.3649 | TRUE   |
| CRAM   | projected_sql |       2 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   46.81 |    174.10156 |       315.3649 | TRUE   |
| CRAM   | projected_sql |       2 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   47.14 |    171.41797 |       315.3649 | TRUE   |
| CRAM   | projected_sql |       2 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   46.76 |    170.22656 |       315.3649 | TRUE   |
| CRAM   | riker         |       4 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.87 |    135.14453 |       369.4586 | TRUE   |
| CRAM   | riker         |       4 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.80 |    139.55859 |       369.4586 | TRUE   |
| CRAM   | riker         |       4 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.79 |    139.32422 |       369.4586 | TRUE   |
| CRAM   | native        |       4 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.04 |    175.43359 |       369.4586 | TRUE   |
| CRAM   | native        |       4 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.03 |    173.17578 |       369.4586 | TRUE   |
| CRAM   | native        |       4 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.02 |    172.79297 |       369.4586 | TRUE   |
| CRAM   | projected_sql |       4 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |   12.07 |    165.53125 |       369.4586 | TRUE   |
| CRAM   | projected_sql |       4 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |   12.09 |    165.14453 |       369.4586 | TRUE   |
| CRAM   | projected_sql |       4 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |   12.10 |    163.50781 |       369.4586 | TRUE   |
| CRAM   | riker         |       4 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    5.59 |    136.46875 |       369.4586 | TRUE   |
| CRAM   | riker         |       4 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    5.54 |    140.40625 |       369.4586 | TRUE   |
| CRAM   | riker         |       4 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    5.61 |    137.71875 |       369.4586 | TRUE   |
| CRAM   | native        |       4 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.70 |    211.42188 |       369.4586 | TRUE   |
| CRAM   | native        |       4 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.64 |    209.15234 |       369.4586 | TRUE   |
| CRAM   | native        |       4 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.67 |    211.76172 |       369.4586 | TRUE   |
| CRAM   | projected_sql |       4 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.95 |    175.45312 |       369.4586 | TRUE   |
| CRAM   | projected_sql |       4 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   24.40 |    174.71094 |       369.4586 | TRUE   |
| CRAM   | projected_sql |       4 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.78 |    175.85156 |       369.4586 | TRUE   |
| CRAM   | riker         |       4 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   10.87 |    141.26562 |       369.4586 | TRUE   |
| CRAM   | riker         |       4 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   11.04 |    140.94531 |       369.4586 | TRUE   |
| CRAM   | riker         |       4 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   11.27 |    141.64453 |       369.4586 | TRUE   |
| CRAM   | native        |       4 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    3.03 |    219.74219 |       369.4586 | TRUE   |
| CRAM   | native        |       4 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    2.89 |    217.39844 |       369.4586 | TRUE   |
| CRAM   | native        |       4 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    2.93 |    216.43750 |       369.4586 | TRUE   |
| CRAM   | projected_sql |       4 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   47.48 |    182.67578 |       369.4586 | TRUE   |
| CRAM   | projected_sql |       4 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   47.91 |    181.98047 |       369.4586 | TRUE   |
| CRAM   | projected_sql |       4 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   48.29 |    192.50781 |       369.4586 | TRUE   |
| CRAM   | riker         |       8 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.02 |    235.35156 |       477.6461 | TRUE   |
| CRAM   | riker         |       8 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.92 |    243.92188 |       477.6461 | TRUE   |
| CRAM   | riker         |       8 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.08 |    249.27344 |       477.6461 | TRUE   |
| CRAM   | native        |       8 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    0.88 |    215.31250 |       477.6461 | TRUE   |
| CRAM   | native        |       8 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    0.82 |    215.72266 |       477.6461 | TRUE   |
| CRAM   | native        |       8 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    0.81 |    212.91016 |       477.6461 | TRUE   |
| CRAM   | projected_sql |       8 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |   12.32 |    182.80859 |       477.6461 | TRUE   |
| CRAM   | projected_sql |       8 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |   11.83 |    187.90625 |       477.6461 | TRUE   |
| CRAM   | projected_sql |       8 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |   11.95 |    186.83984 |       477.6461 | TRUE   |
| CRAM   | riker         |       8 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    3.98 |    263.53516 |       477.6461 | TRUE   |
| CRAM   | riker         |       8 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    4.07 |    255.49609 |       477.6461 | TRUE   |
| CRAM   | riker         |       8 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    4.23 |    263.97266 |       477.6461 | TRUE   |
| CRAM   | native        |       8 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.36 |    270.20312 |       477.6461 | TRUE   |
| CRAM   | native        |       8 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.24 |    268.25781 |       477.6461 | TRUE   |
| CRAM   | native        |       8 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.22 |    268.37891 |       477.6461 | TRUE   |
| CRAM   | projected_sql |       8 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.64 |    212.58984 |       477.6461 | TRUE   |
| CRAM   | projected_sql |       8 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   24.08 |    213.80078 |       477.6461 | TRUE   |
| CRAM   | projected_sql |       8 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   24.03 |    198.82812 |       477.6461 | TRUE   |
| CRAM   | riker         |       8 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    7.90 |    256.58594 |       477.6461 | TRUE   |
| CRAM   | riker         |       8 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    8.23 |    244.69922 |       477.6461 | TRUE   |
| CRAM   | riker         |       8 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    8.28 |    260.38281 |       477.6461 | TRUE   |
| CRAM   | native        |       8 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    2.43 |    285.20312 |       477.6461 | TRUE   |
| CRAM   | native        |       8 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    2.10 |    285.58594 |       477.6461 | TRUE   |
| CRAM   | native        |       8 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    2.11 |    287.08594 |       477.6461 | TRUE   |
| CRAM   | projected_sql |       8 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   47.73 |    222.26953 |       477.6461 | TRUE   |
| CRAM   | projected_sql |       8 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   47.51 |    223.59375 |       477.6461 | TRUE   |
| CRAM   | projected_sql |       8 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   47.80 |    212.25000 |       477.6461 | TRUE   |

RSS criterion: PASS

BAM 1 native scaling gate not applicable (1x \< 5 s).

BAM 1 projected_sql scaling exponent 0.9692997 PASS

BAM 1 riker scaling gate not applicable (1x \< 5 s).

BAM 2 native scaling gate not applicable (1x \< 5 s).

BAM 2 projected_sql scaling exponent 0.9621246 PASS

BAM 2 riker scaling gate not applicable (1x \< 5 s).

BAM 4 native scaling gate not applicable (1x \< 5 s).

BAM 4 projected_sql scaling exponent 0.9782114 PASS

BAM 4 riker scaling gate not applicable (1x \< 5 s).

BAM 8 native scaling gate not applicable (1x \< 5 s).

BAM 8 projected_sql scaling exponent 0.9613206 PASS

BAM 8 riker scaling gate not applicable (1x \< 5 s).

CRAM 1 native scaling gate not applicable (1x \< 5 s).

CRAM 1 projected_sql scaling exponent 0.9834166 PASS

CRAM 1 riker scaling exponent 0.9908385 PASS

CRAM 2 native scaling gate not applicable (1x \< 5 s).

CRAM 2 projected_sql scaling exponent 0.9794907 PASS

CRAM 2 riker scaling exponent 0.9970376 PASS

CRAM 4 native scaling gate not applicable (1x \< 5 s).

CRAM 4 projected_sql scaling exponent 0.9932563 PASS

CRAM 4 riker scaling gate not applicable (1x \< 5 s).

CRAM 8 native scaling gate not applicable (1x \< 5 s).

CRAM 8 projected_sql scaling exponent 0.9989429 PASS

CRAM 8 riker scaling gate not applicable (1x \< 5 s).

| format | threads | target | native_sql_ratio | native_riker_ratio | sql_riker_ratio | native_pass | sql_pass |
|:-------|--------:|-------:|-----------------:|-------------------:|----------------:|:------------|:---------|
| BAM    |       1 |   1.05 |        0.6106557 |          1.0955882 |        1.794118 | FALSE       | FALSE    |
| BAM    |       2 |   1.00 |        0.3414634 |          0.7138810 |        2.090652 | TRUE        | FALSE    |
| BAM    |       4 |   0.90 |        0.2065663 |          1.1439394 |        5.537879 | FALSE       | FALSE    |
| BAM    |       8 |   0.90 |        0.1485944 |          1.0571429 |        7.114286 | FALSE       | FALSE    |
| CRAM   |       1 |   1.05 |        0.2185154 |          0.3149038 |        1.441106 | TRUE        | FALSE    |
| CRAM   |       2 |   1.00 |        0.1295681 |          0.2131148 |        1.644809 | TRUE        | FALSE    |
| CRAM   |       4 |   0.90 |        0.0851944 |          0.3678571 |        4.317857 | TRUE        | FALSE    |
| CRAM   |       8 |   0.90 |        0.0686192 |          0.4059406 |        5.915842 | TRUE        | FALSE    |

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
