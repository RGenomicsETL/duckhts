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

Source revision: `d29b378fede2741b16359ca41665e29b43c174bd`; working
tree changed: `FALSE`; driver SHA-256:
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

Launcher affinity: `pid 3's current affinity list: 0-7`; niceness: `0`.

Host: `Linux; 6.8.0-78-generic; x86_64`; load at start:
`23:22:49 up 373 days,  8:16, 19 users,  load average: 2.45, 3.10, 3.20`.

Cache condition: warm; each input is scanned for independent oracle
outputs before timing. Page cache is not flushed. Concurrent work
declared by launcher: verified depth-input staging on CPUs 12-19.

## Results

| format | arm                                | threads | scale | replicate |  records | read_pairs |  bins | metric_rows | histogram_rows | seconds | work_seconds | peak_rss_mib | buffer_mib | spill_bytes | rss_budget_mib | parity |
|:-------|:-----------------------------------|--------:|------:|----------:|---------:|-----------:|------:|------------:|---------------:|--------:|-------------:|-------------:|-----------:|------------:|---------------:|:-------|
| BAM    | riker                              |       1 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    3.96 |           NA |     25.93750 |         NA |          NA |       288.3180 | TRUE   |
| BAM    | riker                              |       1 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    4.02 |           NA |     26.01953 |         NA |          NA |       288.3180 | TRUE   |
| BAM    | riker                              |       1 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    3.93 |           NA |     25.92188 |         NA |          NA |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       1 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    4.10 |     4.056327 |     58.61719 |   20.40625 |           0 |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       1 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    4.12 |     4.071332 |     57.81250 |   20.40234 |           0 |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       1 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    4.14 |     4.090580 |     58.44531 |   20.34375 |           0 |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       1 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.04 |     6.993955 |     58.37500 |   20.87109 |           0 |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       1 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.01 |     6.966973 |     57.98828 |   20.88672 |           0 |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       1 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.16 |     7.116711 |     59.01953 |   20.95703 |           0 |       288.3180 | TRUE   |
| BAM    | riker                              |       2 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    3.33 |           NA |     26.56250 |         NA |          NA |       315.3649 | TRUE   |
| BAM    | riker                              |       2 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    3.33 |           NA |     26.50000 |         NA |          NA |       315.3649 | TRUE   |
| BAM    | riker                              |       2 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    3.33 |           NA |     26.93750 |         NA |          NA |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       2 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.12 |     2.071614 |     65.68359 |   29.50781 |           0 |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       2 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.12 |     2.074291 |     65.26562 |   30.49609 |           0 |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       2 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.13 |     2.084040 |     64.53516 |   29.60938 |           0 |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       2 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.10 |     7.054090 |     64.79688 |   29.93750 |           0 |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       2 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.06 |     7.014046 |     64.89062 |   29.50781 |           0 |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       2 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.01 |     6.961507 |     63.98828 |   29.45703 |           0 |       315.3649 | TRUE   |
| BAM    | riker                              |       4 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.26 |           NA |     26.48828 |         NA |          NA |       369.4586 | TRUE   |
| BAM    | riker                              |       4 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.29 |           NA |     26.48047 |         NA |          NA |       369.4586 | TRUE   |
| BAM    | riker                              |       4 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.37 |           NA |     26.67578 |         NA |          NA |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       4 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.15 |     1.103986 |     70.67578 |   36.44141 |           0 |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       4 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.18 |     1.132644 |     67.44141 |   36.45703 |           0 |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       4 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.15 |     1.102791 |     69.68359 |   35.38281 |           0 |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       4 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.01 |     6.964632 |     69.78516 |   36.23047 |           0 |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       4 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.06 |     7.008923 |     69.01953 |   36.56250 |           0 |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       4 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.22 |     7.176481 |     68.49609 |   35.81641 |           0 |       369.4586 | TRUE   |
| BAM    | riker                              |       8 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.16 |           NA |     26.25000 |         NA |          NA |       477.6461 | TRUE   |
| BAM    | riker                              |       8 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.14 |           NA |     26.56250 |         NA |          NA |       477.6461 | TRUE   |
| BAM    | riker                              |       8 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.16 |           NA |     26.52344 |         NA |          NA |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       8 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    0.85 |     0.797218 |     89.05469 |   60.14453 |           0 |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       8 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    0.89 |     0.837646 |     86.71875 |   58.54492 |           0 |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       8 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    0.84 |     0.795681 |     87.43750 |   60.67578 |           0 |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       8 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.32 |     7.266855 |     81.14453 |   58.86328 |           0 |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       8 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.17 |     7.119493 |     82.21484 |   53.89453 |           0 |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       8 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.07 |     7.025916 |     82.49219 |   60.41992 |           0 |       477.6461 | TRUE   |
| BAM    | riker                              |       1 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    5.65 |           NA |     25.91406 |         NA |          NA |       288.3180 | TRUE   |
| BAM    | riker                              |       1 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    5.66 |           NA |     25.78125 |         NA |          NA |       288.3180 | TRUE   |
| BAM    | riker                              |       1 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    5.65 |           NA |     26.02344 |         NA |          NA |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       1 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    6.05 |     5.986829 |     58.56250 |   20.26562 |           0 |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       1 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    6.08 |     6.035486 |     57.78516 |   20.33203 |           0 |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       1 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    6.03 |     5.980312 |     58.44922 |   20.49609 |           0 |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       1 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   11.88 |    11.831133 |     58.69141 |   20.92188 |           0 |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       1 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   12.02 |    11.967886 |     57.95312 |   20.89453 |           0 |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       1 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   11.96 |    11.915165 |     59.01953 |   20.87500 |           0 |       288.3180 | TRUE   |
| BAM    | riker                              |       2 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    4.79 |           NA |     26.38672 |         NA |          NA |       315.3649 | TRUE   |
| BAM    | riker                              |       2 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    4.67 |           NA |     26.47656 |         NA |          NA |       315.3649 | TRUE   |
| BAM    | riker                              |       2 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    4.64 |           NA |     26.66797 |         NA |          NA |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       2 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    3.14 |     3.090817 |     67.12109 |   31.25781 |           0 |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       2 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    3.09 |     3.044762 |     64.77734 |   29.99609 |           0 |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       2 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    3.24 |     3.195256 |     65.00781 |   30.02734 |           0 |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       2 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   12.05 |    11.992365 |     65.60547 |   29.66797 |           0 |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       2 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   11.89 |    11.845322 |     64.23438 |   29.41016 |           0 |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       2 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   12.04 |    11.990423 |     65.06250 |   29.90625 |           0 |       315.3649 | TRUE   |
| BAM    | riker                              |       4 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.11 |           NA |     26.64844 |         NA |          NA |       369.4586 | TRUE   |
| BAM    | riker                              |       4 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.27 |           NA |     26.93750 |         NA |          NA |       369.4586 | TRUE   |
| BAM    | riker                              |       4 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.24 |           NA |     26.48828 |         NA |          NA |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       4 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.65 |     1.605406 |     69.42188 |   38.12500 |           0 |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       4 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.78 |     1.732374 |     69.85547 |   36.98047 |           0 |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       4 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.80 |     1.750127 |     70.82812 |   36.92969 |           0 |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       4 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   11.92 |    11.873872 |     69.37500 |   36.79297 |           0 |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       4 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   12.02 |    11.975398 |     69.38672 |   37.86328 |           0 |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       4 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   11.96 |    11.909132 |     69.53906 |   36.72656 |           0 |       369.4586 | TRUE   |
| BAM    | riker                              |       8 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.97 |           NA |     26.56250 |         NA |          NA |       477.6461 | TRUE   |
| BAM    | riker                              |       8 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.97 |           NA |     26.55859 |         NA |          NA |       477.6461 | TRUE   |
| BAM    | riker                              |       8 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.95 |           NA |     26.31250 |         NA |          NA |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       8 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.22 |     1.167346 |     88.79688 |   62.01367 |           0 |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       8 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.55 |     1.497346 |     85.52734 |   50.24219 |           0 |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       8 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.38 |     1.332457 |     86.00000 |   54.05078 |           0 |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       8 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   11.89 |    11.845427 |     83.14062 |   55.26367 |           0 |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       8 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   11.91 |    11.867784 |     81.90625 |   47.89453 |           0 |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       8 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   11.92 |    11.870086 |     81.59375 |   59.54492 |           0 |       477.6461 | TRUE   |
| BAM    | riker                              |       1 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    8.53 |           NA |     25.78125 |         NA |          NA |       288.3180 | TRUE   |
| BAM    | riker                              |       1 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    8.57 |           NA |     25.91016 |         NA |          NA |       288.3180 | TRUE   |
| BAM    | riker                              |       1 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    8.50 |           NA |     25.62500 |         NA |          NA |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       1 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    9.47 |     9.424880 |     58.54297 |   20.33594 |           0 |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       1 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    9.45 |     9.406099 |     58.97656 |   20.34375 |           0 |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       1 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    9.44 |     9.390906 |     57.62500 |   20.40234 |           0 |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       1 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   21.02 |    20.970418 |     58.16797 |   20.84766 |           0 |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       1 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   21.09 |    21.040650 |     58.06641 |   20.83594 |           0 |       288.3180 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       1 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   21.20 |    21.148479 |     58.13672 |   20.83203 |           0 |       288.3180 | TRUE   |
| BAM    | riker                              |       2 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    6.80 |           NA |     26.46875 |         NA |          NA |       315.3649 | TRUE   |
| BAM    | riker                              |       2 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    6.82 |           NA |     26.40625 |         NA |          NA |       315.3649 | TRUE   |
| BAM    | riker                              |       2 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    6.89 |           NA |     26.22656 |         NA |          NA |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       2 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    4.79 |     4.741569 |     64.94531 |   31.63672 |           0 |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       2 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    4.83 |     4.779691 |     64.90234 |   30.06250 |           0 |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       2 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    4.80 |     4.750222 |     66.87109 |   30.73828 |           0 |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       2 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   21.06 |    21.009905 |     64.87500 |   29.47656 |           0 |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       2 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   21.20 |    21.150367 |     65.67578 |   30.01953 |           0 |       315.3649 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       2 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   21.52 |    21.470591 |     65.58203 |   30.00781 |           0 |       315.3649 | TRUE   |
| BAM    | riker                              |       4 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    3.30 |           NA |     26.38281 |         NA |          NA |       369.4586 | TRUE   |
| BAM    | riker                              |       4 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    3.34 |           NA |     26.06250 |         NA |          NA |       369.4586 | TRUE   |
| BAM    | riker                              |       4 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    2.83 |           NA |     26.33203 |         NA |          NA |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       4 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    3.28 |     3.234968 |     69.61719 |   36.35156 |           0 |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       4 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    3.29 |     3.238943 |     69.98047 |   37.65234 |           0 |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       4 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    2.80 |     2.751978 |     69.94922 |   37.43359 |           0 |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       4 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   21.16 |    21.118285 |     68.31641 |   36.82422 |           0 |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       4 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   21.03 |    20.979728 |     68.66016 |   37.49609 |           0 |       369.4586 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       4 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   21.04 |    20.993071 |     69.90234 |   38.05078 |           0 |       369.4586 | TRUE   |
| BAM    | riker                              |       8 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    3.49 |           NA |     26.56250 |         NA |          NA |       477.6461 | TRUE   |
| BAM    | riker                              |       8 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    3.48 |           NA |     26.06250 |         NA |          NA |       477.6461 | TRUE   |
| BAM    | riker                              |       8 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    3.53 |           NA |     26.15234 |         NA |          NA |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       8 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    2.19 |     2.138877 |     86.53125 |   62.70703 |           0 |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       8 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    2.20 |     2.150745 |     87.56641 |   61.23242 |           0 |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts     |       8 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    2.20 |     2.156754 |     87.12891 |   61.92578 |           0 |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       8 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   21.07 |    21.024245 |     82.93359 |   57.83203 |           0 |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       8 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   21.14 |    21.088479 |     82.41406 |   52.26953 |           0 |       477.6461 | TRUE   |
| BAM    | duckhts_bam_insert_size_counts_sql |       8 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   21.16 |    21.108592 |     83.88281 |   59.70703 |           0 |       477.6461 | TRUE   |
| CRAM   | riker                              |       1 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    8.20 |           NA |     43.88672 |         NA |          NA |       288.3180 | TRUE   |
| CRAM   | riker                              |       1 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    8.23 |           NA |     43.81641 |         NA |          NA |       288.3180 | TRUE   |
| CRAM   | riker                              |       1 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    8.22 |           NA |     43.49609 |         NA |          NA |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       1 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.93 |     2.883073 |     59.72656 |   20.41016 |           0 |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       1 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.98 |     2.934612 |     59.37109 |   20.34375 |           0 |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       1 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.72 |     2.672705 |     59.37109 |   20.34766 |           0 |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       1 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |   12.02 |    11.969781 |     69.54297 |   20.94141 |           0 |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       1 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |   11.84 |    11.789980 |     68.25391 |   20.79297 |           0 |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       1 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |   11.84 |    11.796008 |     66.84766 |   20.91797 |           0 |       288.3180 | TRUE   |
| CRAM   | riker                              |       2 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.16 |           NA |    111.34766 |         NA |          NA |       315.3649 | TRUE   |
| CRAM   | riker                              |       2 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.19 |           NA |    112.24219 |         NA |          NA |       315.3649 | TRUE   |
| CRAM   | riker                              |       2 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    7.27 |           NA |    111.45312 |         NA |          NA |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       2 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.36 |     1.316741 |     65.19922 |   30.10156 |           0 |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       2 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.27 |     1.227175 |     64.51953 |   29.77734 |           0 |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       2 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.27 |     1.225453 |     73.61719 |   31.23438 |           0 |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       2 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |   11.90 |    11.854241 |     71.35938 |   29.48047 |           0 |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       2 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |   11.85 |    11.783469 |     72.02344 |   29.41016 |           0 |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       2 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |   11.94 |    11.874801 |     71.78125 |   29.46875 |           0 |       315.3649 | TRUE   |
| CRAM   | riker                              |       4 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.86 |           NA |    145.23047 |         NA |          NA |       369.4586 | TRUE   |
| CRAM   | riker                              |       4 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    3.03 |           NA |    135.96484 |         NA |          NA |       369.4586 | TRUE   |
| CRAM   | riker                              |       4 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.95 |           NA |    139.48828 |         NA |          NA |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       4 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    0.96 |     0.911215 |     84.00000 |   31.69141 |           0 |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       4 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    0.96 |     0.916812 |     81.76562 |   36.91016 |           0 |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       4 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    1.01 |     0.964989 |     81.69141 |   35.74609 |           0 |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       4 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |   11.86 |    11.816099 |     79.79297 |   34.52734 |           0 |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       4 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |   11.98 |    11.932963 |     78.12500 |   34.20703 |           0 |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       4 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |   11.88 |    11.835177 |     90.17188 |   34.95703 |           0 |       369.4586 | TRUE   |
| CRAM   | riker                              |       8 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.19 |           NA |    244.25781 |         NA |          NA |       477.6461 | TRUE   |
| CRAM   | riker                              |       8 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.59 |           NA |    252.00000 |         NA |          NA |       477.6461 | TRUE   |
| CRAM   | riker                              |       8 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    2.66 |           NA |    251.24609 |         NA |          NA |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       8 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |    0.69 |     0.640418 |    121.82031 |   54.45117 |           0 |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       8 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |    0.69 |     0.635705 |    120.95312 |   67.17578 |           0 |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       8 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |    0.67 |     0.623440 |    120.92578 |   61.98828 |           0 |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       8 |     1 |         1 | 16816430 |    7238229 | 32393 |           1 |           1112 |   12.12 |    12.070252 |     90.55859 |   59.88867 |           0 |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       8 |     1 |         2 | 16816430 |    7238229 | 32393 |           1 |           1112 |   12.06 |    12.011969 |    101.87891 |   59.79492 |           0 |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       8 |     1 |         3 | 16816430 |    7238229 | 32393 |           1 |           1112 |   11.96 |    11.908579 |     91.93750 |   55.11328 |           0 |       477.6461 | TRUE   |
| CRAM   | riker                              |       1 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   16.19 |           NA |     50.46484 |         NA |          NA |       288.3180 | TRUE   |
| CRAM   | riker                              |       1 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   16.48 |           NA |     50.33594 |         NA |          NA |       288.3180 | TRUE   |
| CRAM   | riker                              |       1 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   16.46 |           NA |     50.41406 |         NA |          NA |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       1 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    5.29 |     5.238944 |     59.60547 |   20.33203 |           0 |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       1 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    5.31 |     5.258000 |     59.44922 |   20.42188 |           0 |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       1 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    6.34 |     6.297926 |     58.88281 |   20.34375 |           0 |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       1 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.59 |    23.538686 |     68.68750 |   20.89844 |           0 |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       1 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.52 |    23.474127 |     71.20312 |   20.89453 |           0 |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       1 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.26 |    23.215320 |     67.16406 |   21.00000 |           0 |       288.3180 | TRUE   |
| CRAM   | riker                              |       2 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.25 |           NA |    114.96094 |         NA |          NA |       315.3649 | TRUE   |
| CRAM   | riker                              |       2 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.25 |           NA |    114.46875 |         NA |          NA |       315.3649 | TRUE   |
| CRAM   | riker                              |       2 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   14.33 |           NA |    113.83594 |         NA |          NA |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       2 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.79 |     2.739612 |     72.42188 |   31.26953 |           0 |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       2 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.77 |     2.717091 |     65.01172 |   30.55859 |           0 |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       2 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    2.78 |     2.729367 |     65.32031 |   29.55859 |           0 |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       2 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.35 |    23.293472 |     79.53125 |   29.96875 |           0 |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       2 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.71 |    23.664144 |     74.70703 |   29.38672 |           0 |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       2 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.67 |    23.618956 |     80.13672 |   29.46875 |           0 |       315.3649 | TRUE   |
| CRAM   | riker                              |       4 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    5.20 |           NA |    139.61328 |         NA |          NA |       369.4586 | TRUE   |
| CRAM   | riker                              |       4 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    5.21 |           NA |    141.96094 |         NA |          NA |       369.4586 | TRUE   |
| CRAM   | riker                              |       4 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    5.24 |           NA |    147.14453 |         NA |          NA |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       4 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.51 |     1.463486 |     84.48828 |   37.23828 |           0 |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       4 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.54 |     1.486757 |     88.07422 |   35.96094 |           0 |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       4 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.53 |     1.479767 |     83.42188 |   35.88281 |           0 |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       4 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.70 |    23.652790 |     91.73438 |   37.97266 |           0 |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       4 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.67 |    23.615406 |     88.73828 |   37.10156 |           0 |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       4 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.65 |    23.597279 |     88.09375 |   36.97656 |           0 |       369.4586 | TRUE   |
| CRAM   | riker                              |       8 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    4.39 |           NA |    261.23047 |         NA |          NA |       477.6461 | TRUE   |
| CRAM   | riker                              |       8 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    4.55 |           NA |    243.32812 |         NA |          NA |       477.6461 | TRUE   |
| CRAM   | riker                              |       8 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    4.82 |           NA |    252.28906 |         NA |          NA |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       8 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.42 |     1.366670 |    119.43359 |   61.98828 |           0 |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       8 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.38 |     1.325053 |    118.38672 |   68.07617 |           0 |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       8 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |    1.36 |     1.315307 |    118.71875 |   59.36328 |           0 |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       8 |     2 |         1 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.41 |    23.357089 |    104.08203 |   59.56055 |           0 |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       8 |     2 |         2 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.77 |    23.719608 |     97.51953 |   56.61328 |           0 |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       8 |     2 |         3 | 33632860 |   14476458 | 32393 |           1 |           1112 |   23.67 |    23.615903 |     99.83594 |   61.06055 |           0 |       477.6461 | TRUE   |
| CRAM   | riker                              |       1 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   26.59 |           NA |     48.93750 |         NA |          NA |       288.3180 | TRUE   |
| CRAM   | riker                              |       1 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   26.57 |           NA |     49.48047 |         NA |          NA |       288.3180 | TRUE   |
| CRAM   | riker                              |       1 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   26.67 |           NA |     48.95312 |         NA |          NA |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       1 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    8.93 |     8.881057 |     60.16016 |   20.40234 |           0 |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       1 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    9.18 |     9.127000 |     60.75000 |   20.42188 |           0 |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       1 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    9.22 |     9.174425 |     58.74609 |   20.41406 |           0 |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       1 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   40.81 |    40.754525 |     75.78125 |   20.86719 |           0 |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       1 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   40.84 |    40.788578 |     71.53125 |   20.88672 |           0 |       288.3180 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       1 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   40.80 |    40.753580 |     75.25391 |   20.88281 |           0 |       288.3180 | TRUE   |
| CRAM   | riker                              |       2 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   22.40 |           NA |    113.52734 |         NA |          NA |       315.3649 | TRUE   |
| CRAM   | riker                              |       2 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   22.39 |           NA |    112.83594 |         NA |          NA |       315.3649 | TRUE   |
| CRAM   | riker                              |       2 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   22.21 |           NA |    113.70312 |         NA |          NA |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       2 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    4.72 |     4.676023 |     66.09375 |   30.11719 |           0 |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       2 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    4.40 |     4.355428 |     65.58984 |   30.03125 |           0 |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       2 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    4.57 |     4.521748 |     67.29297 |   30.83594 |           0 |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       2 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   40.43 |    40.383233 |     80.85547 |   29.46875 |           0 |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       2 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   40.35 |    40.303628 |     81.70703 |   29.42969 |           0 |       315.3649 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       2 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   40.32 |    40.267339 |     79.77344 |   28.97656 |           0 |       315.3649 | TRUE   |
| CRAM   | riker                              |       4 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    9.53 |           NA |    150.40625 |         NA |          NA |       369.4586 | TRUE   |
| CRAM   | riker                              |       4 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    9.43 |           NA |    145.03906 |         NA |          NA |       369.4586 | TRUE   |
| CRAM   | riker                              |       4 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    9.05 |           NA |    150.28516 |         NA |          NA |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       4 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    2.61 |     2.560923 |     86.29688 |   36.38281 |           0 |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       4 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    2.56 |     2.511571 |     85.68750 |   35.89453 |           0 |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       4 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    2.85 |     2.802407 |     83.95312 |   37.46875 |           0 |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       4 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   40.57 |    40.525276 |     87.99609 |   38.63672 |           0 |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       4 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   40.62 |    40.569762 |     88.60547 |   37.03516 |           0 |       369.4586 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       4 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   41.07 |    41.016346 |     91.12891 |   39.41016 |           0 |       369.4586 | TRUE   |
| CRAM   | riker                              |       8 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    8.57 |           NA |    270.45703 |         NA |          NA |       477.6461 | TRUE   |
| CRAM   | riker                              |       8 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    8.95 |           NA |    257.88281 |         NA |          NA |       477.6461 | TRUE   |
| CRAM   | riker                              |       8 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    8.88 |           NA |    272.86719 |         NA |          NA |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       8 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |    2.35 |     2.303012 |    127.66797 |   55.40430 |           0 |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       8 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |    2.25 |     2.204624 |    122.90234 |   58.71680 |           0 |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts     |       8 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |    2.29 |     2.238788 |    126.07422 |   61.77930 |           0 |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       8 |     4 |         1 | 67265720 |   28952916 | 32393 |           1 |           1112 |   40.76 |    40.711984 |    104.46094 |   53.92578 |           0 |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       8 |     4 |         2 | 67265720 |   28952916 | 32393 |           1 |           1112 |   40.50 |    40.446442 |    105.87500 |   61.15430 |           0 |       477.6461 | TRUE   |
| CRAM   | duckhts_bam_insert_size_counts_sql |       8 |     4 |         3 | 67265720 |   28952916 | 32393 |           1 |           1112 |   40.81 |    40.759789 |    102.88672 |   59.23438 |           0 |       477.6461 | TRUE   |

RSS criterion: PASS

BAM 1 duckhts_bam_insert_size_counts scaling gate not applicable (1x \<
5 s).

BAM 1 duckhts_bam_insert_size_counts_sql 1x to 2x scaling exponent
0.7645701 PASS

BAM 1 duckhts_bam_insert_size_counts_sql 2x to 4x scaling exponent
0.8183417 PASS

BAM 1 riker scaling gate not applicable (1x \< 5 s).

BAM 2 duckhts_bam_insert_size_counts scaling gate not applicable (1x \<
5 s).

BAM 2 duckhts_bam_insert_size_counts_sql 1x to 2x scaling exponent
0.7700953 PASS

BAM 2 duckhts_bam_insert_size_counts_sql 2x to 4x scaling exponent
0.8162289 PASS

BAM 2 riker scaling gate not applicable (1x \< 5 s).

BAM 4 duckhts_bam_insert_size_counts scaling gate not applicable (1x \<
5 s).

BAM 4 duckhts_bam_insert_size_counts_sql 1x to 2x scaling exponent
0.7604773 PASS

BAM 4 duckhts_bam_insert_size_counts_sql 2x to 4x scaling exponent
0.8149173 PASS

BAM 4 riker scaling gate not applicable (1x \< 5 s).

BAM 8 duckhts_bam_insert_size_counts scaling gate not applicable (1x \<
5 s).

BAM 8 duckhts_bam_insert_size_counts_sql 1x to 2x scaling exponent
0.7321284 PASS

BAM 8 duckhts_bam_insert_size_counts_sql 2x to 4x scaling exponent
0.827802 PASS

BAM 8 riker scaling gate not applicable (1x \< 5 s).

CRAM 1 duckhts_bam_insert_size_counts scaling gate not applicable (1x \<
5 s).

CRAM 1 duckhts_bam_insert_size_counts_sql 1x to 2x scaling exponent
0.990219 PASS

CRAM 1 duckhts_bam_insert_size_counts_sql 2x to 4x scaling exponent
0.7950347 PASS

CRAM 1 riker 1x to 2x scaling exponent 1.001754 PASS

CRAM 1 riker 2x to 4x scaling exponent 0.6919194 PASS

CRAM 2 duckhts_bam_insert_size_counts scaling gate not applicable (1x \<
5 s).

CRAM 2 duckhts_bam_insert_size_counts_sql 1x to 2x scaling exponent
0.9920981 PASS

CRAM 2 duckhts_bam_insert_size_counts_sql 2x to 4x scaling exponent
0.769509 PASS

CRAM 2 riker 1x to 2x scaling exponent 0.9868982 PASS

CRAM 2 riker 2x to 4x scaling exponent 0.6518926 PASS

CRAM 4 duckhts_bam_insert_size_counts scaling gate not applicable (1x \<
5 s).

CRAM 4 duckhts_bam_insert_size_counts_sql 1x to 2x scaling exponent
0.9945249 PASS

CRAM 4 duckhts_bam_insert_size_counts_sql 2x to 4x scaling exponent
0.7791305 PASS

CRAM 4 riker scaling gate not applicable (1x \< 5 s).

CRAM 8 duckhts_bam_insert_size_counts scaling gate not applicable (1x \<
5 s).

CRAM 8 duckhts_bam_insert_size_counts_sql 1x to 2x scaling exponent
0.9728298 PASS

CRAM 8 duckhts_bam_insert_size_counts_sql 2x to 4x scaling exponent
0.7840943 PASS

CRAM 8 riker scaling gate not applicable (1x \< 5 s).

## Reported-only in-process context

`work_seconds` excludes CLI startup and LOAD/SET, retaining the query,
profiling and TSV output. The remainder includes startup/setup and
end-of-process work; it is not an isolated startup measurement or a
pass/fail criterion.

| format | arm                                | threads | scale | seconds | work_seconds | outside_work_seconds |
|:-------|:-----------------------------------|--------:|------:|--------:|-------------:|---------------------:|
| BAM    | duckhts_bam_insert_size_counts     |       1 |     1 |    4.12 |     4.071332 |             0.048668 |
| CRAM   | duckhts_bam_insert_size_counts     |       1 |     1 |    2.93 |     2.883073 |             0.046927 |
| BAM    | duckhts_bam_insert_size_counts_sql |       1 |     1 |    7.04 |     6.993955 |             0.046045 |
| CRAM   | duckhts_bam_insert_size_counts_sql |       1 |     1 |   11.84 |    11.796008 |             0.043992 |
| BAM    | duckhts_bam_insert_size_counts     |       2 |     1 |    2.12 |     2.074291 |             0.045709 |
| CRAM   | duckhts_bam_insert_size_counts     |       2 |     1 |    1.27 |     1.227175 |             0.042825 |
| BAM    | duckhts_bam_insert_size_counts_sql |       2 |     1 |    7.06 |     7.014046 |             0.045954 |
| CRAM   | duckhts_bam_insert_size_counts_sql |       2 |     1 |   11.90 |    11.854241 |             0.045759 |
| BAM    | duckhts_bam_insert_size_counts     |       4 |     1 |    1.15 |     1.103986 |             0.046014 |
| CRAM   | duckhts_bam_insert_size_counts     |       4 |     1 |    0.96 |     0.916812 |             0.043188 |
| BAM    | duckhts_bam_insert_size_counts_sql |       4 |     1 |    7.06 |     7.008923 |             0.051077 |
| CRAM   | duckhts_bam_insert_size_counts_sql |       4 |     1 |   11.88 |    11.835177 |             0.044823 |
| BAM    | duckhts_bam_insert_size_counts     |       8 |     1 |    0.85 |     0.797218 |             0.052782 |
| CRAM   | duckhts_bam_insert_size_counts     |       8 |     1 |    0.69 |     0.635705 |             0.054295 |
| BAM    | duckhts_bam_insert_size_counts_sql |       8 |     1 |    7.17 |     7.119493 |             0.050507 |
| CRAM   | duckhts_bam_insert_size_counts_sql |       8 |     1 |   12.06 |    12.011969 |             0.048031 |
| BAM    | duckhts_bam_insert_size_counts     |       1 |     2 |    6.05 |     5.986829 |             0.063171 |
| CRAM   | duckhts_bam_insert_size_counts     |       1 |     2 |    5.31 |     5.258000 |             0.052000 |
| BAM    | duckhts_bam_insert_size_counts_sql |       1 |     2 |   11.96 |    11.915165 |             0.044835 |
| CRAM   | duckhts_bam_insert_size_counts_sql |       1 |     2 |   23.52 |    23.474127 |             0.045873 |
| BAM    | duckhts_bam_insert_size_counts     |       2 |     2 |    3.14 |     3.090817 |             0.049183 |
| CRAM   | duckhts_bam_insert_size_counts     |       2 |     2 |    2.78 |     2.729367 |             0.050633 |
| BAM    | duckhts_bam_insert_size_counts_sql |       2 |     2 |   12.04 |    11.990423 |             0.049577 |
| CRAM   | duckhts_bam_insert_size_counts_sql |       2 |     2 |   23.67 |    23.618956 |             0.051044 |
| BAM    | duckhts_bam_insert_size_counts     |       4 |     2 |    1.78 |     1.732374 |             0.047626 |
| CRAM   | duckhts_bam_insert_size_counts     |       4 |     2 |    1.53 |     1.479767 |             0.050233 |
| BAM    | duckhts_bam_insert_size_counts_sql |       4 |     2 |   11.96 |    11.909132 |             0.050868 |
| CRAM   | duckhts_bam_insert_size_counts_sql |       4 |     2 |   23.67 |    23.615406 |             0.054594 |
| BAM    | duckhts_bam_insert_size_counts     |       8 |     2 |    1.38 |     1.332457 |             0.047543 |
| CRAM   | duckhts_bam_insert_size_counts     |       8 |     2 |    1.38 |     1.325053 |             0.054947 |
| BAM    | duckhts_bam_insert_size_counts_sql |       8 |     2 |   11.91 |    11.867784 |             0.042216 |
| CRAM   | duckhts_bam_insert_size_counts_sql |       8 |     2 |   23.67 |    23.615903 |             0.054097 |
| BAM    | duckhts_bam_insert_size_counts     |       1 |     4 |    9.45 |     9.406099 |             0.043901 |
| CRAM   | duckhts_bam_insert_size_counts     |       1 |     4 |    9.18 |     9.127000 |             0.053000 |
| BAM    | duckhts_bam_insert_size_counts_sql |       1 |     4 |   21.09 |    21.040650 |             0.049350 |
| CRAM   | duckhts_bam_insert_size_counts_sql |       1 |     4 |   40.81 |    40.754525 |             0.055475 |
| BAM    | duckhts_bam_insert_size_counts     |       2 |     4 |    4.80 |     4.750222 |             0.049778 |
| CRAM   | duckhts_bam_insert_size_counts     |       2 |     4 |    4.57 |     4.521748 |             0.048252 |
| BAM    | duckhts_bam_insert_size_counts_sql |       2 |     4 |   21.20 |    21.150367 |             0.049633 |
| CRAM   | duckhts_bam_insert_size_counts_sql |       2 |     4 |   40.35 |    40.303628 |             0.046372 |
| BAM    | duckhts_bam_insert_size_counts     |       4 |     4 |    3.28 |     3.234968 |             0.045032 |
| CRAM   | duckhts_bam_insert_size_counts     |       4 |     4 |    2.61 |     2.560923 |             0.049077 |
| BAM    | duckhts_bam_insert_size_counts_sql |       4 |     4 |   21.04 |    20.993071 |             0.046929 |
| CRAM   | duckhts_bam_insert_size_counts_sql |       4 |     4 |   40.62 |    40.569762 |             0.050238 |
| BAM    | duckhts_bam_insert_size_counts     |       8 |     4 |    2.20 |     2.150745 |             0.049255 |
| CRAM   | duckhts_bam_insert_size_counts     |       8 |     4 |    2.29 |     2.238788 |             0.051212 |
| BAM    | duckhts_bam_insert_size_counts_sql |       8 |     4 |   21.14 |    21.088479 |             0.051521 |
| CRAM   | duckhts_bam_insert_size_counts_sql |       8 |     4 |   40.76 |    40.711984 |             0.048016 |

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

| format | arm                                | threads | ratio_1x | target | verdict | ratio_4x | marginal_ratio |
|:-------|:-----------------------------------|--------:|---------:|-------:|:--------|---------:|---------------:|
| BAM    | duckhts_bam_insert_size_counts     |       1 |    1.040 |   1.05 | PASS    |    1.108 |          1.166 |
| BAM    | duckhts_bam_insert_size_counts_sql |       1 |    1.778 |   1.05 | FAIL    |    2.472 |          3.074 |
| BAM    | duckhts_bam_insert_size_counts     |       2 |    0.637 |   1.00 | PASS    |    0.704 |          0.768 |
| BAM    | duckhts_bam_insert_size_counts_sql |       2 |    2.120 |   1.00 | FAIL    |    3.109 |          4.052 |
| BAM    | duckhts_bam_insert_size_counts     |       4 |    0.891 |   0.90 | PASS    |    0.994 |          1.060 |
| BAM    | duckhts_bam_insert_size_counts_sql |       4 |    5.473 |   0.90 | FAIL    |    6.376 |          6.955 |
| BAM    | duckhts_bam_insert_size_counts     |       8 |    0.733 |   0.90 | PASS    |    0.630 |          0.579 |
| BAM    | duckhts_bam_insert_size_counts_sql |       8 |    6.181 |   0.90 | FAIL    |    6.057 |          5.996 |
| CRAM   | duckhts_bam_insert_size_counts     |       1 |    0.356 |   1.05 | PASS    |    0.345 |          0.340 |
| CRAM   | duckhts_bam_insert_size_counts_sql |       1 |    1.440 |   1.05 | FAIL    |    1.535 |          1.577 |
| CRAM   | duckhts_bam_insert_size_counts     |       2 |    0.177 |   1.00 | PASS    |    0.204 |          0.217 |
| CRAM   | duckhts_bam_insert_size_counts_sql |       2 |    1.655 |   1.00 | FAIL    |    1.802 |          1.872 |
| CRAM   | duckhts_bam_insert_size_counts     |       4 |    0.325 |   0.90 | PASS    |    0.277 |          0.255 |
| CRAM   | duckhts_bam_insert_size_counts_sql |       4 |    4.027 |   0.90 | FAIL    |    4.308 |          4.435 |
| CRAM   | duckhts_bam_insert_size_counts     |       8 |    0.266 |   0.90 | PASS    |    0.258 |          0.254 |
| CRAM   | duckhts_bam_insert_size_counts_sql |       8 |    4.656 |   0.90 | FAIL    |    4.590 |          4.563 |

Native/SQL comparison (reported, not a pass/fail gate):

| format | threads | native_seconds_1x | sql_seconds_1x | native_sql_ratio_1x |
|:-------|--------:|------------------:|---------------:|--------------------:|
| BAM    |       1 |              4.12 |           7.04 |               0.585 |
| BAM    |       2 |              2.12 |           7.06 |               0.300 |
| BAM    |       4 |              1.15 |           7.06 |               0.163 |
| BAM    |       8 |              0.85 |           7.17 |               0.119 |
| CRAM   |       1 |              2.93 |          11.84 |               0.247 |
| CRAM   |       2 |              1.27 |          11.90 |               0.107 |
| CRAM   |       4 |              0.96 |          11.88 |               0.081 |
| CRAM   |       8 |              0.69 |          12.06 |               0.057 |

## SQL thread scaling and the next building block

BAM SQL 1x medians span 7.040–7.170 s ; largest/smallest = 1.018 .

Thread scaling is flat to within 5% across these measured medians.

CRAM SQL 1x medians span 11.840–12.060 s ; largest/smallest = 1.019 .

Thread scaling is flat to within 5% across these measured medians.

| format | threads | seconds | speedup_vs_1_thread |
|:-------|--------:|--------:|--------------------:|
| BAM    |       1 |    7.04 |               1.000 |
| BAM    |       2 |    7.06 |               0.997 |
| BAM    |       4 |    7.06 |               0.997 |
| BAM    |       8 |    7.17 |               0.982 |
| CRAM   |       1 |   11.84 |               1.000 |
| CRAM   |       2 |   11.90 |               0.995 |
| CRAM   |       4 |   11.88 |               0.997 |
| CRAM   |       8 |   12.06 |               0.982 |

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
