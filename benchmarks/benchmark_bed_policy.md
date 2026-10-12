BED error-policy reader benchmark
================

The public DuckBedQC `duckbedqc_118fc21` corpus supplies
`data/GRCh38_exons.bed` (SHA-256
`75b120e20be18497716c81d81648bb919574b34b90cc15abc474e06f09924073`,
22,802,546 bytes). Stage it with `scripts/stage_duckbedqc_data.sh`. The
baseline is an independent `origin/develop` worktree at a5ab5147
(`BED_POLICY_BASE`, its extension build); the candidate source is
68dfa0332ad9ea34eea556eddea231cd83d842f6 (`BED_POLICY_BRANCH`, its
extension build). Both revisions are available in this repository. The
query forces a full scan with `scan_mode := 'sequential'` and aggregates
`count(*)`, `sum(start)`, and `sum("end")`. Both numeric sums are
checked across policies and builds. Each run reads the 22 MB file and
returns one aggregate row. R/DBI wall time includes query planning and
result materialization. DuckDB R 1.5.5 loads the extension in a fresh R
process for each build. The unindexed file is clean under the reader’s
short-row rule. “default” is the engine’s unmodified thread setting,
recorded alongside the single-thread runs.

The host mixes two core types, so the render records the CPU affinity of
its process, which every child R process inherits, in both observation
files, and compares with the previous render only on the same CPUs.
Choose them with `taskset -c`. `BED_POLICY_SIZE=reduced` runs a driver
check: two scans per cell and two alternating pairs, written to
temporary files instead of `bed_policy_observations.tsv` and
`bed_policy_alternating.tsv`; it is not a measurement.

``` sh
BED_POLICY_BASE=<baseline extension> BED_POLICY_BASE_REV=<baseline revision> \
BED_POLICY_BRANCH=<candidate extension> BED_POLICY_REV=<candidate revision> \
taskset -c 8-11 Rscript -e 'rmarkdown::render("benchmarks/benchmark_bed_policy.Rmd")'
```

``` r
stopifnot(nzchar(Sys.getenv("BED_POLICY_BASE")),
          nzchar(Sys.getenv("BED_POLICY_BRANCH")),
          nzchar(Sys.getenv("BED_POLICY_REV")),
          nzchar(Sys.getenv("BED_POLICY_BASE_REV")))
size <- Sys.getenv("BED_POLICY_SIZE", "full")
if (!size %in% c("full", "reduced")) stop("BED_POLICY_SIZE must be full or reduced")
full_size <- identical(size, "full")
scans_per_cell <- if (full_size) 11L else 2L
alternating_pairs <- if (full_size) 15L else 2L
observations_path <- if (full_size) "benchmarks/bed_policy_observations.tsv" else tempfile(fileext = ".tsv")
alternating_path <- if (full_size) "benchmarks/bed_policy_alternating.tsv" else tempfile(fileext = ".tsv")
# The CPUs this process may run on, as taskset prints them; child processes
# inherit the mask.
affinity_line <- system2("taskset", c("-pc", Sys.getpid()), stdout = TRUE)
affinity <- sub("^.*: *", "", affinity_line[[length(affinity_line)]])
run_benchmark <- function(extension, label, revision) {
  output <- tempfile(fileext = ".tsv")
  result <- system2("Rscript", c("benchmarks/bed_policy_run.R", extension,
                                   label, revision, output, scans_per_cell))
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0L) {
    stop("benchmark failed: ", label)
  }
  read.delim(output, check.names = FALSE)
}
observations <- rbind(
  run_benchmark(Sys.getenv("BED_POLICY_BRANCH"), "candidate",
                Sys.getenv("BED_POLICY_REV")),
  run_benchmark(Sys.getenv("BED_POLICY_BASE"), "baseline",
                Sys.getenv("BED_POLICY_BASE_REV"))
)
stopifnot(length(unique(observations$input_lines)) == 1L,
          length(unique(observations$output_rows)) == 1L,
          length(unique(observations$start_sum)) == 1L,
          length(unique(observations$end_sum)) == 1L,
          all(table(interaction(observations$build, observations$threads,
                                observations$policy, drop = TRUE)) == scans_per_cell))
observations$affinity <- affinity
write.table(observations, observations_path,
            sep = "\t", row.names = FALSE, quote = FALSE)
summary <- aggregate(seconds ~ revision + build + threads + effective_threads +
                       policy + input_lines + output_rows,
                     observations, median)
summary$seconds <- round(summary$seconds, 4)
knitr::kable(summary, caption = sprintf("Median wall seconds over %d independent scans per policy/thread/build",
                                        scans_per_cell))
```

| revision                                 | build     | threads | effective_threads | policy | input_lines | output_rows | seconds |
|:-----------------------------------------|:----------|:--------|------------------:|:-------|------------:|------------:|--------:|
| a5ab5147                                 | baseline  | 1       |                 1 | error  |      389852 |      389852 |   0.038 |
| 68dfa0332ad9ea34eea556eddea231cd83d842f6 | candidate | 1       |                 1 | error  |      389852 |      389852 |   0.038 |
| a5ab5147                                 | baseline  | default |                20 | error  |      389852 |      389852 |   0.039 |
| 68dfa0332ad9ea34eea556eddea231cd83d842f6 | candidate | default |                20 | error  |      389852 |      389852 |   0.038 |
| 68dfa0332ad9ea34eea556eddea231cd83d842f6 | candidate | 1       |                 1 | report |      389852 |      389852 |   0.039 |
| 68dfa0332ad9ea34eea556eddea231cd83d842f6 | candidate | default |                20 | report |      389852 |      389852 |   0.039 |
| 68dfa0332ad9ea34eea556eddea231cd83d842f6 | candidate | 1       |                 1 | skip   |      389852 |      389852 |   0.039 |
| 68dfa0332ad9ea34eea556eddea231cd83d842f6 | candidate | default |                20 | skip   |      389852 |      389852 |   0.039 |

Median wall seconds over 11 independent scans per policy/thread/build

Default policy, 1 threads: 0.038 s to 0.038 s (+0.0%).

Default policy, default threads: 0.039 s to 0.038 s (-2.6%).

The alternating comparison runs the exact default-policy aggregate with
one DuckDB thread in a fresh R process per observation. The baseline and
candidate alternate fifteen times each (twice in a reduced run); process
startup is outside the reported query wall time.

``` r
load_at_start <- system("uptime", intern = TRUE)
run_single <- function(extension, label, run) {
  lines <- system2("Rscript", c("benchmarks/bed_policy_single.R",
                                  extension, label, run), stdout = TRUE)
  if (!is.null(attr(lines, "status")) && attr(lines, "status") != 0L) {
    stop("alternating scan failed: ", label, " run ", run)
  }
  read.delim(text = paste(lines, collapse = "\n"), check.names = FALSE)
}
alternating <- do.call(rbind, lapply(seq_len(alternating_pairs), function(run) {
  rbind(run_single(Sys.getenv("BED_POLICY_BASE"), "baseline", run),
        run_single(Sys.getenv("BED_POLICY_BRANCH"), "candidate", run))
}))
stopifnot(nrow(alternating) == 2L * alternating_pairs,
          length(unique(alternating$result)) == 1L,
          all(table(alternating$build) == alternating_pairs))
alternating$affinity <- affinity
write.table(alternating, alternating_path,
            sep = "\t", row.names = FALSE, quote = FALSE)
alternating_summary <- aggregate(seconds ~ build, alternating,
                                 function(x) c(min = min(x), median = median(x)))
alternating_summary <- data.frame(
  build = alternating_summary$build,
  min = round(alternating_summary$seconds[, "min"], 4),
  median = round(alternating_summary$seconds[, "median"], 4)
)
knitr::kable(alternating_summary,
             caption = sprintf("Alternating fresh-process scans (%d per build, one thread)", alternating_pairs))
```

| build     |   min | median |
|:----------|------:|-------:|
| baseline  | 0.038 |  0.039 |
| candidate | 0.039 |  0.040 |

Alternating fresh-process scans (15 per build, one thread)

Host load before alternating scans: 14:42:11 up 369 days, 23:36, 19
users, load average: 1.18, 1.29, 0.98. CPU affinity: 8-11.

The reader’s only fatal data-line check is fewer than three
tab-delimited fields. Numeric parse failures yield NULL and reversed
intervals are not rejected. This workload does not measure malformed-row
reporting or remote I/O.

## Change from the previous render

The previous render, 8ca12e37, ran the same scans of the same file with
the same driver and host, with candidate
1d87d9e38155596a1cb751e76cfdcbdf203fa2ef and baseline
f741354dc5012f94f59229e922ed5317c07e04fc; “candidate” and “baseline”
name different revisions in the two renders. The previous render
recorded no CPU affinity, so the CPUs it ran on are unknown and the
ratios may include a change of core type. Host load differs between
renders. Median wall seconds per cell:

| build     | threads | policy | effective_threads_before | seconds_before | effective_threads_after | seconds_after | time_ratio |
|:----------|:--------|:-------|:-------------------------|---------------:|------------------------:|--------------:|-----------:|
| baseline  | 1       | error  | 1                        |          0.039 |                       1 |         0.038 |     0.9744 |
| baseline  | default | error  | 20                       |          0.040 |                      20 |         0.039 |     0.9750 |
| candidate | 1       | error  | 1                        |          0.039 |                       1 |         0.038 |     0.9744 |
| candidate | 1       | report | 1                        |          0.039 |                       1 |         0.039 |     1.0000 |
| candidate | 1       | skip   | 1                        |          0.039 |                       1 |         0.039 |     1.0000 |
| candidate | default | error  | 20                       |          0.039 |                      20 |         0.038 |     0.9744 |
| candidate | default | report | 20                       |          0.040 |                      20 |         0.039 |     0.9750 |
| candidate | default | skip   | 20                       |          0.040 |                      20 |         0.039 |     0.9750 |

Alternating fresh-process scans, median per build:

| build     | seconds_before | seconds_after | time_ratio |
|:----------|---------------:|--------------:|-----------:|
| baseline  |          0.041 |         0.039 |     0.9512 |
| candidate |          0.040 |         0.040 |     1.0000 |
