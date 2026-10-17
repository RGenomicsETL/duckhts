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
