# Consequence-extension footprint

## Measurement

Linux amd64, R 4.6.0, DuckDB R 1.5.5; in-tree CMake release build with `make release -j2`. Baseline is `3b90520b443d504286695007bb035ca452fa685b`; candidate native sources are at `42a64cc9`. Both builds used the same host and toolchain. The baseline was built in a fresh detached worktree after `make configure`; the candidate was built after `make -C third_party/htslib -f Makefile distclean` and `make clean` in an already configured worktree. The times are one build each, not a repeated-build speed estimate.

| Measure | Baseline | Candidate | Difference |
|:--|--:|--:|--:|
| Release extension bytes | 2,698,686 | 2,152,798 | −545,888 (−20.2%) |
| Tracked `src/` bytes | 3,792,806 | 2,201,812 | −1,590,994 (−41.9%) |
| `CREATE OR REPLACE MACRO` definitions in C source | 44 | 31 | −13 |
| Clean release build, seconds | 30.46 | 24.43 | −6.03 |

Source bytes are the sum of tracked files under `src/`, including headers. Macro counts use `rg -io 'CREATE OR REPLACE MACRO ' src --glob '*.[ch]' | wc -l`; they count source definitions, not installed catalog functions. The candidate extension SHA-256 is `ca856e7fe520d42cf1a2ea50263840305843ae3d604e403c6ef9c4e9de49e5c6`; baseline is `3cb0c116fcfe2aa163d21f77f3ad9581b82b08fce9c77ed0b8897b0f9a8c66ba`.

## LOAD timing

`Rscript scripts/benchmark_extension_load.R BASELINE CANDIDATE 15 benchmarks/data/consequence_extension_load.csv` alternates builds across 15 fresh R processes per build, loading exactly one unsigned extension into one in-memory DuckDB connection per process. DuckDB is configured with one thread and a private home. The timed interval is only the `LOAD` statement; R startup and connection construction are excluded. The CSV retains all 30 observations.

| Build | Median seconds | Range seconds | Loads |
|:--|--:|:--|--:|
| Baseline | 0.025 | 0.024–0.026 | 15 |
| Candidate | 0.016 | 0.015–0.017 | 15 |

This measures fixed startup cost; there are no input rows or output rows. The nearest checked-in reader baseline is [`benchmark_init_readers.md`](benchmark_init_readers.md), which times `LOAD` and five reader workloads on different revisions and must not be used as a matched reader-performance comparison. No 1×/2×/4× input-scaling claim or memory budget applies to this one-extension load measurement.
