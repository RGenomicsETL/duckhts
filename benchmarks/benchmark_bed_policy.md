BED error-policy reader benchmark
================

The public DuckBedQC `duckbedqc_118fc21` corpus supplies
`data/GRCh38_exons.bed` (SHA-256
`75b120e20be18497716c81d81648bb919574b34b90cc15abc474e06f09924073`,
22,802,546 bytes). Stage it with `scripts/stage_duckbedqc_data.sh`. The
baseline is an independent `git archive` of
`81edf9bff85ab2daabeccd1daaf13e4bfa39ee47` built with
`make configure && make release -j4`; the candidate is the revision
below. The query forces a full scan with `scan_mode := 'sequential'` and
aggregates `count(*)`, `sum(start)`, and `sum("end")`. Both numeric sums
are checked across policies and builds. Each run reads the 22 MB file
and returns one aggregate row. R/DBI wall time includes query planning
and result materialization. DuckDB R 1.5.5 loads the extension in a
fresh R process for each build. The unindexed file is clean under the
reader’s short-row rule. “default” is the engine’s unmodified thread
setting, recorded alongside the single-thread runs.

``` r
stopifnot(nzchar(Sys.getenv("BED_POLICY_BASE")),
          nzchar(Sys.getenv("BED_POLICY_BRANCH")),
          nzchar(Sys.getenv("BED_POLICY_REV")))
run_benchmark <- function(extension, label, revision) {
  output <- tempfile(fileext = ".tsv")
  result <- system2("Rscript", c("benchmarks/bed_policy_run.R", extension,
                                   label, revision, output))
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0L) {
    stop("benchmark failed: ", label)
  }
  read.delim(output, check.names = FALSE)
}
observations <- rbind(
  run_benchmark(Sys.getenv("BED_POLICY_BRANCH"), "candidate",
                Sys.getenv("BED_POLICY_REV")),
  run_benchmark(Sys.getenv("BED_POLICY_BASE"), "baseline",
                "81edf9bff85ab2daabeccd1daaf13e4bfa39ee47")
)
stopifnot(length(unique(observations$input_lines)) == 1L,
          length(unique(observations$output_rows)) == 1L,
          length(unique(observations$start_sum)) == 1L,
          length(unique(observations$end_sum)) == 1L,
          all(table(interaction(observations$build, observations$threads,
                                observations$policy, drop = TRUE)) == 11L))
summary <- aggregate(seconds ~ revision + build + threads + effective_threads +
                       policy + input_lines + output_rows,
                     observations, median)
summary$seconds <- round(summary$seconds, 4)
knitr::kable(summary, caption = "Median wall seconds over eleven independent scans per policy/thread/build")
```

| revision                                 | build     | threads | effective_threads | policy | input_lines | output_rows | seconds |
|:-----------------------------------------|:----------|:--------|------------------:|:-------|------------:|------------:|--------:|
| 81edf9bff85ab2daabeccd1daaf13e4bfa39ee47 | baseline  | 1       |                 1 | error  |      389852 |      389852 |   0.039 |
| 22d4bde23866798a51ea0ab6a4460aed5a326194 | candidate | 1       |                 1 | error  |      389852 |      389852 |   0.042 |
| 81edf9bff85ab2daabeccd1daaf13e4bfa39ee47 | baseline  | default |                20 | error  |      389852 |      389852 |   0.039 |
| 22d4bde23866798a51ea0ab6a4460aed5a326194 | candidate | default |                20 | error  |      389852 |      389852 |   0.041 |
| 22d4bde23866798a51ea0ab6a4460aed5a326194 | candidate | 1       |                 1 | report |      389852 |      389852 |   0.042 |
| 22d4bde23866798a51ea0ab6a4460aed5a326194 | candidate | default |                20 | report |      389852 |      389852 |   0.043 |
| 22d4bde23866798a51ea0ab6a4460aed5a326194 | candidate | 1       |                 1 | skip   |      389852 |      389852 |   0.042 |
| 22d4bde23866798a51ea0ab6a4460aed5a326194 | candidate | default |                20 | skip   |      389852 |      389852 |   0.042 |

Median wall seconds over eleven independent scans per
policy/thread/build

Default policy, 1 threads: 0.039 s to 0.042 s (+7.7%).

Default policy, default threads: 0.039 s to 0.041 s (+5.1%).

The alternating comparison runs the exact default-policy aggregate with
one DuckDB thread in a fresh R process per observation. The baseline and
candidate alternate fifteen times each; process startup is outside the
reported query wall time.

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
alternating <- do.call(rbind, lapply(seq_len(15L), function(run) {
  rbind(run_single(Sys.getenv("BED_POLICY_BASE"), "baseline", run),
        run_single(Sys.getenv("BED_POLICY_BRANCH"), "candidate", run))
}))
stopifnot(nrow(alternating) == 30L,
          length(unique(alternating$result)) == 1L,
          all(table(alternating$build) == 15L))
alternating_summary <- aggregate(seconds ~ build, alternating,
                                 function(x) c(min = min(x), median = median(x)))
alternating_summary <- data.frame(
  build = alternating_summary$build,
  min = round(alternating_summary$seconds[, "min"], 4),
  median = round(alternating_summary$seconds[, "median"], 4)
)
knitr::kable(alternating_summary,
             caption = "Alternating fresh-process scans (15 per build, one thread)")
```

| build     |   min | median |
|:----------|------:|-------:|
| baseline  | 0.039 |  0.041 |
| candidate | 0.039 |  0.042 |

Alternating fresh-process scans (15 per build, one thread)

Host load before alternating scans: 20:23:28 up 353 days, 5:17, 21
users, load average: 0.87, 1.05, 1.14.

The reader’s only fatal data-line check is fewer than three
tab-delimited fields. Numeric parse failures yield NULL and reversed
intervals are not rejected. This workload does not measure malformed-row
reporting or remote I/O.
