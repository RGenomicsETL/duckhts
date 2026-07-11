#!/usr/bin/env Rscript

# Run the pure-C randomized properties and record every theft target reported by
# the executable. A failed process records nothing.

suppressMessages({
  library(blit)
  library(glue)
  library(optparse)
})
options(rlang_backtrace_on_error = "none")

root <- tryCatch(
  system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE),
  error = function(e) "."
)
op <- OptionParser()
op <- add_option(op, "--trials", type = "double", default = 100000)
op <- add_option(op, "--seed", default = "0xd0c0ffee12345678")
op <- add_option(
  op,
  "--history",
  default = file.path(
    root,
    "test",
    "duckvep",
    "conformance",
    "data",
    "property_history.csv"
  )
)
op <- add_option(
  op,
  "--source-revision",
  dest = "source_revision",
  default = ""
)
opt <- parse_args(op)

die <- function(...) stop(glue(..., .envir = parent.frame()), call. = FALSE)
if (
  !is.finite(opt$trials) || opt$trials < 1 || opt$trials != floor(opt$trials)
) {
  die("--trials must be a positive integer")
}
configured_trials <- format(opt$trials, scientific = FALSE, trim = TRUE)
if (!grepl("^(0x[0-9A-Fa-f]+|[0-9]+)$", opt$seed)) {
  die("--seed must be an unsigned decimal or hexadecimal integer")
}
if (utils::packageVersion("blit") < "0.2.0.9000") {
  die("the current WangLabCSU/blit checkout is required")
}

log_path <- tempfile(fileext = ".log")
on.exit(unlink(log_path), add = TRUE)
command <- blit::exec("make", "test-duckvep-kernel") |>
  blit::cmd_wd(root) |>
  blit::cmd_envvar(
    DUCKVEP_PROP_TRIALS = configured_trials,
    DUCKVEP_PROP_SEED = opt$seed
  )
status <- suppressWarnings(blit::cmd_run(
  command,
  stdout = log_path,
  stderr = "2>&1",
  stdin = NULL,
  verbose = FALSE
))
log <- readLines(log_path, warn = FALSE)
if (status != 0L) {
  die(
    "property run failed with exit status {status}:\n",
    "{paste(tail(log, 80L), collapse = '\n')}"
  )
}

capture <- function(lines, pattern, fields) {
  matched <- regexec(pattern, lines, perl = TRUE)
  values <- regmatches(lines, matched)
  values <- values[lengths(values) != 0L]
  if (length(values) == 0L) {
    die("property log did not match: {pattern}")
  }
  out <- as.data.frame(do.call(rbind, lapply(values, `[`, -1L)))
  names(out) <- fields
  out
}
started <- capture(
  log,
  "^== PROP '(.+)': ([0-9]+) trials, seed (0x[0-9a-fA-F]+)$",
  c("target", "trials", "reported_seed")
)
finished <- capture(
  log,
  "^== PASS '(.+)': pass ([0-9]+), fail ([0-9]+), skip ([0-9]+), dup ([0-9]+)$",
  c("target", "passed", "failed", "skipped", "duplicates")
)
properties <- merge(started, finished, by = "target", all = TRUE, sort = FALSE)
if (
  nrow(properties) != nrow(started) ||
    nrow(properties) != nrow(finished) ||
    anyNA(properties)
) {
  die("property start/pass records do not pair exactly")
}
number_columns <- c("trials", "passed", "failed", "skipped", "duplicates")
properties[number_columns] <- lapply(properties[number_columns], as.numeric)
if (
  any(properties$failed != 0) || any(properties$passed != properties$trials)
) {
  die("a property target did not pass every requested trial")
}

total <- capture(
  log,
  "^Total: ([0-9]+) tests .* ([0-9.]+) sec\\), ([0-9]+) assertions$",
  c("suite_tests", "elapsed_seconds", "suite_assertions")
)
suite <- capture(
  log,
  "^Pass: ([0-9]+), fail: ([0-9]+), skip: ([0-9]+)[.]$",
  c("suite_passed", "suite_failed", "suite_skipped")
)
if (
  nrow(total) != 1L || nrow(suite) != 1L || as.integer(suite$suite_failed) != 0L
) {
  die("unexpected greatest summary")
}

if (!nzchar(opt$source_revision)) {
  revision <- suppressWarnings(system2(
    "git",
    c("-C", root, "rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
  ))
  revision_status <- attr(revision, "status")
  if (!is.null(revision_status) && revision_status != 0L) {
    die("cannot determine source revision")
  }
  opt$source_revision <- trimws(revision[[1L]])
}
cc <- Sys.getenv("CC", "cc")
cc_version <- suppressWarnings(system2(
  cc,
  "--version",
  stdout = TRUE,
  stderr = TRUE
))
cc_version <- if (length(cc_version)) cc_version[[1L]] else cc

rows <- data.frame(
  run_date = as.character(Sys.Date()),
  source_revision = opt$source_revision,
  configured_trials = as.numeric(opt$trials),
  seed = properties$reported_seed,
  target = properties$target,
  trials = properties$trials,
  passed = properties$passed,
  failed = properties$failed,
  skipped = properties$skipped,
  duplicates = properties$duplicates,
  suite_tests = as.numeric(total$suite_tests),
  suite_assertions = as.numeric(total$suite_assertions),
  suite_elapsed_seconds = as.numeric(total$elapsed_seconds),
  compiler = cc_version,
  stringsAsFactors = FALSE
)
rows <- rows[order(rows$target), , drop = FALSE]

dir.create(dirname(opt$history), recursive = TRUE, showWarnings = FALSE)
if (file.exists(opt$history)) {
  old <- utils::read.csv(
    opt$history,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = c(source_revision = "character", seed = "character")
  )
  if (!identical(names(old), names(rows))) {
    die("history schema does not match: {opt$history}")
  }
  run_key <- paste(
    rows$source_revision[[1L]],
    rows$seed[[1L]],
    rows$configured_trials[[1L]]
  )
  old_key <- paste(old$source_revision, old$seed, old$configured_trials)
  rows <- rbind(old[old_key != run_key, , drop = FALSE], rows)
}
rows <- rows[
  order(rows$run_date, rows$source_revision, rows$seed, rows$target),
  ,
  drop = FALSE
]
tmp <- tempfile("duckvep-properties-", dirname(opt$history), ".csv")
utils::write.csv(rows, tmp, row.names = FALSE)
if (!file.rename(tmp, opt$history)) {
  unlink(tmp)
  die("cannot replace history: {opt$history}")
}

cat(
  glue(
    "{nrow(properties)} randomized targets; ",
    "{format(sum(properties$trials), big.mark = ',', scientific = FALSE)} trials; ",
    "{total$suite_assertions} assertions; {total$elapsed_seconds} s"
  ),
  "\n",
  sep = ""
)
cat(glue("history -> {opt$history}"), "\n", sep = "")
