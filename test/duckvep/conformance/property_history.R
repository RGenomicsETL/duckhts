#!/usr/bin/env Rscript

# Run the pure-C randomized properties and record every theft target reported by
# the executable. Failed campaigns retain their complete seed-specific log.

suppressMessages({
  library(blit)
  library(glue)
  library(optparse)
})
options(rlang_backtrace_on_error = "none")

main <- function() {
root <- tryCatch(
  system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE),
  error = function(e) "."
)
source(file.path(root, "scripts", "duckvep_evidence.R"), local = TRUE)
campaign_started <- FALSE
campaign_complete <- FALSE
log_path <- NULL
failure_log_path <- NULL
failure_log_reported <- FALSE
preserve_failure_log <- function() {
  if (!isTRUE(campaign_started) || is.null(log_path) || !file.exists(log_path)) {
    return(NULL)
  }
  if (!is.null(failure_log_path)) {
    return(failure_log_path)
  }
  dir.create(opt$failure_log_dir, recursive = TRUE, showWarnings = FALSE)
  safe_seed <- gsub("[^A-Za-z0-9_.-]", "_", opt$seed)
  stamp <- format(Sys.time(), "%Y%m%dT%H%M%S")
  failure_log_path <<- file.path(
    opt$failure_log_dir,
    paste0(
      "property_failure_seed_",
      safe_seed,
      "_trials_",
      configured_trials,
      "_",
      stamp,
      "_pid_",
      Sys.getpid(),
      ".log"
    )
  )
  if (!file.copy(log_path, failure_log_path, overwrite = FALSE)) {
    failure_log_path <<- NULL
  }
  failure_log_path
}
die <- function(...) {
  message <- as.character(glue(..., .envir = parent.frame()))
  preserved <- preserve_failure_log()
  if (!is.null(preserved)) {
    message <- paste0(message, "\nfull property log: ", preserved)
    failure_log_reported <<- TRUE
  }
  stop(message, call. = FALSE)
}

on.exit({
  if (isTRUE(campaign_started) && !isTRUE(campaign_complete)) {
    preserved <- preserve_failure_log()
    if (!is.null(preserved) && !isTRUE(failure_log_reported)) {
      message("full property log: ", preserved)
    }
  }
  if (!is.null(log_path)) {
    unlink(log_path)
  }
}, add = TRUE)

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
  "--coverage-requirements",
  dest = "coverage_requirements",
  default = file.path(
    root,
    "test",
    "duckvep",
    "conformance",
    "data",
    "property_coverage_requirements.tsv"
  )
)
op <- add_option(
  op,
  "--failure-log-dir",
  dest = "failure_log_dir",
  default = file.path(root, "test", "duckvep", "conformance", "results")
)
op <- add_option(
  op,
  "--coverage-history",
  dest = "coverage_history",
  default = file.path(
    root,
    "test",
    "duckvep",
    "conformance",
    "data",
    "property_coverage_history.csv"
  )
)
op <- add_option(
  op,
  "--source-revision",
  dest = "source_revision",
  default = ""
)
opt <- parse_args(op)

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

if (normalizePath(opt$history, mustWork = FALSE) ==
    normalizePath(opt$coverage_history, mustWork = FALSE)) {
  die("property and coverage histories must be different files")
}
dir.create(dirname(opt$history), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$coverage_history), recursive = TRUE, showWarnings = FALSE)
history_dir <- normalizePath(dirname(opt$history), mustWork = TRUE)
coverage_dir <- normalizePath(dirname(opt$coverage_history), mustWork = TRUE)
if (!identical(history_dir, coverage_dir)) {
  die("property and coverage histories must share one publication directory")
}
history_path <- file.path(history_dir, basename(opt$history))
coverage_history_path <- file.path(
  coverage_dir,
  basename(opt$coverage_history)
)
# All publications in one directory share a lock. Pair-specific locks do not
# protect overlapping pairs such as (A, B) and (A, C), both of which replace A.
lock_dir <- file.path(history_dir, ".duckvep-property-history.lock")
recovery_dir <- file.path(
  history_dir,
  ".duckvep-property-history-recovery.lock"
)
owner_path <- file.path(lock_dir, "owner")
journal_path <- file.path(lock_dir, "journal.rds")
property_backup <- file.path(lock_dir, "property_history.backup")
coverage_backup <- file.path(lock_dir, "property_coverage_history.backup")

restore_interrupted_publication <- function(interrupted_lock_dir) {
  interrupted_journal <- file.path(interrupted_lock_dir, "journal.rds")
  interrupted_property_backup <- file.path(
    interrupted_lock_dir,
    "property_history.backup"
  )
  interrupted_coverage_backup <- file.path(
    interrupted_lock_dir,
    "property_coverage_history.backup"
  )
  if (!file.exists(interrupted_journal)) {
    return(invisible(TRUE))
  }
  journal <- tryCatch(
    readRDS(interrupted_journal),
    error = function(e) {
      stop(
        "cannot read interrupted property-history journal in ",
        interrupted_lock_dir,
        call. = FALSE
      )
    }
  )
  expected_journal_names <- c(
    "had_properties",
    "had_coverage",
    "history_path",
    "coverage_history_path"
  )
  if (!identical(names(journal), expected_journal_names) ||
      !is.character(journal$history_path) ||
      length(journal$history_path) != 1L ||
      !is.character(journal$coverage_history_path) ||
      length(journal$coverage_history_path) != 1L) {
    stop(
      "invalid interrupted property-history journal in ",
      interrupted_lock_dir,
      call. = FALSE
    )
  }
  interrupted_history <- normalizePath(journal$history_path, mustWork = FALSE)
  interrupted_coverage <- normalizePath(
    journal$coverage_history_path,
    mustWork = FALSE
  )
  if (!identical(dirname(interrupted_history), history_dir) ||
      !identical(dirname(interrupted_coverage), history_dir) ||
      identical(interrupted_history, interrupted_coverage)) {
    stop(
      "interrupted property-history journal names invalid destinations",
      call. = FALSE
    )
  }
  ok <- TRUE
  if (isTRUE(journal$had_properties)) {
    ok <- file.copy(
      interrupted_property_backup,
      interrupted_history,
      overwrite = TRUE
    ) && ok
  } else if (file.exists(interrupted_history)) {
    ok <- unlink(interrupted_history) == 0L && ok
  }
  if (isTRUE(journal$had_coverage)) {
    ok <- file.copy(
      interrupted_coverage_backup,
      interrupted_coverage,
      overwrite = TRUE
    ) && ok
  } else if (file.exists(interrupted_coverage)) {
    ok <- unlink(interrupted_coverage) == 0L && ok
  }
  if (!ok) {
    stop(
      "cannot recover interrupted property-history publication in ",
      interrupted_lock_dir,
      call. = FALSE
    )
  }
  invisible(TRUE)
}

lock_owner_state <- function(candidate_lock_dir) {
  candidate_owner_path <- file.path(candidate_lock_dir, "owner")
  owner <- suppressWarnings(readLines(
    candidate_owner_path,
    n = 1L,
    warn = FALSE
  ))
  owner_pid <- suppressWarnings(as.integer(owner))
  lock_info <- file.info(candidate_lock_dir)
  lock_age <- suppressWarnings(as.numeric(difftime(
    Sys.time(), lock_info$mtime, units = "secs"
  )))
  owner_alive <- length(owner_pid) == 1L && !is.na(owner_pid) &&
    isTRUE(tryCatch(
      tools::pskill(owner_pid, signal = 0L),
      error = function(e) FALSE,
      warning = function(w) FALSE
    ))
  list(pid = owner_pid, alive = owner_alive, age = lock_age)
}

recover_stale_publication <- function() {
  recovery_owned <- FALSE
  on.exit({
    if (isTRUE(recovery_owned)) {
      unlink(recovery_dir, recursive = TRUE)
    }
  }, add = TRUE)
  if (dir.exists(recovery_dir)) {
    recovery_state <- lock_owner_state(recovery_dir)
    if (recovery_state$alive) {
      die("property-history recovery is active in process {recovery_state$pid}")
    }
    if ((length(recovery_state$pid) != 1L || is.na(recovery_state$pid)) &&
        (length(recovery_state$age) != 1L || is.na(recovery_state$age) ||
         recovery_state$age < 60)) {
      die("property-history recovery lock is still initializing: {recovery_dir}")
    }
    if (unlink(recovery_dir, recursive = TRUE) != 0L) {
      die("cannot remove stale property-history recovery lock: {recovery_dir}")
    }
  }
  if (!dir.exists(lock_dir)) {
    return(invisible(TRUE))
  }
  owner_state <- lock_owner_state(lock_dir)
  if (owner_state$alive) {
    die("property-history publication is locked by process {owner_state$pid}")
  }
  if ((length(owner_state$pid) != 1L || is.na(owner_state$pid)) &&
      (length(owner_state$age) != 1L || is.na(owner_state$age) ||
       owner_state$age < 60)) {
    die("property-history publication lock is still initializing: {lock_dir}")
  }
  if (!dir.create(recovery_dir, showWarnings = FALSE)) {
    die("cannot acquire property-history recovery lock: {recovery_dir}")
  }
  recovery_owned <- TRUE
  writeLines(
    as.character(Sys.getpid()),
    file.path(recovery_dir, "owner"),
    useBytes = TRUE
  )
  owner_state <- lock_owner_state(lock_dir)
  if (owner_state$alive) {
    die("property-history publication is locked by process {owner_state$pid}")
  }
  restore_interrupted_publication(lock_dir)
  if (unlink(lock_dir, recursive = TRUE) != 0L) {
    die("cannot remove stale property-history lock: {lock_dir}")
  }
  invisible(TRUE)
}

# Recovery precedes the cleanliness gate because an interrupted two-file rename
# necessarily leaves at least one tracked ledger dirty until the journal rolls it back.
recover_stale_publication()

git_head <- suppressWarnings(system2(
  "git",
  duckvep_system2_quote(c("-C", root, "rev-parse", "HEAD")),
  stdout = TRUE,
  stderr = FALSE
))
git_head_status <- attr(git_head, "status")
if ((!is.null(git_head_status) && git_head_status != 0L) || length(git_head) != 1L) {
  die("cannot determine source revision")
}
git_head <- trimws(git_head[[1L]])
source_input_paths <- c(
  "Makefile",
  "src/duckvep",
  "test/duckvep/property",
  "test/duckvep/vendor"
)
check_source_tree <- function(expected_head) {
  current_head <- suppressWarnings(system2(
    "git",
    duckvep_system2_quote(c("-C", root, "rev-parse", "HEAD")),
    stdout = TRUE,
    stderr = FALSE
  ))
  current_head_status <- attr(current_head, "status")
  if ((!is.null(current_head_status) && current_head_status != 0L) ||
      length(current_head) != 1L || trimws(current_head[[1L]]) != expected_head) {
    die("source revision changed during the property campaign")
  }
  tracked_changes <- suppressWarnings(system2(
    "git",
    duckvep_system2_quote(c(
      "-C", root, "status", "--porcelain", "--untracked-files=no"
    )),
    stdout = TRUE,
    stderr = FALSE
  ))
  tracked_status <- attr(tracked_changes, "status")
  if (!is.null(tracked_status) && tracked_status != 0L) {
    die("cannot inspect tracked worktree state")
  }
  if (length(tracked_changes)) {
    die(
      "refusing to publish property evidence from a dirty tracked worktree; ",
      "commit the source first and rerun"
    )
  }
  untracked_inputs <- suppressWarnings(system2(
    "git",
    duckvep_system2_quote(c(
      "-C",
      root,
      "ls-files",
      "--others",
      "--",
      source_input_paths
    )),
    stdout = TRUE,
    stderr = FALSE
  ))
  untracked_status <- attr(untracked_inputs, "status")
  if (!is.null(untracked_status) && untracked_status != 0L) {
    die("cannot inspect untracked property-build inputs")
  }
  if (length(untracked_inputs)) {
    die(
      "refusing to publish property evidence with untracked build inputs: ",
      "{paste(untracked_inputs, collapse = ', ')}"
    )
  }
}
check_source_tree(git_head)
if (!nzchar(opt$source_revision)) {
  opt$source_revision <- git_head
} else if (!identical(opt$source_revision, git_head)) {
  die(
    "--source-revision {opt$source_revision} does not match checked-out HEAD ",
    "{git_head}"
  )
}

log_path <- tempfile(fileext = ".log")
command <- do.call(
  blit::exec,
  as.list(duckvep_blit_quote(c(
    Sys.which("make"),
    "-f",
    file.path(root, "Makefile"),
    "DUCKVEP_PROPERTY_ARGS=",
    "DUCKVEP_PROPERTY_CPPFLAGS=",
    "test-duckvep-kernel"
  )))
) |>
  blit::cmd_wd(root) |>
  blit::cmd_envvar(
    MAKEFLAGS = "",
    GNUMAKEFLAGS = "",
    MAKEFILES = "",
    MFLAGS = "",
    DUCKVEP_PROP_TRIALS = configured_trials,
    DUCKVEP_PROP_SEED = opt$seed
  )
campaign_started <- TRUE
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

coverage_lines <- grep("\\[[^]]+ coverage\\]", log, value = TRUE)
coverage <- do.call(
  rbind,
  lapply(coverage_lines, function(line) {
    coverage_name <- sub(
      ".*\\[([^]]+ coverage)\\].*",
      "\\1",
      line,
      perl = TRUE
    )
    payload <- sub(".*\\[[^]]+ coverage\\][[:space:]]*", "", line, perl = TRUE)
    matches <- regmatches(
      payload,
      gregexpr(
        "([A-Za-z_][A-Za-z0-9_]*|[+-][0-9]+)=([0-9]+)",
        payload,
        perl = TRUE
      )
    )[[1L]]
    if (!length(matches) || identical(matches, "")) {
      die("coverage line contains no numeric counters: {line}")
    }
    data.frame(
      coverage = coverage_name,
      metric = sub("=.*", "", matches),
      count = as.numeric(sub("^[^=]+=", "", matches)),
      stringsAsFactors = FALSE
    )
  })
)
if (is.null(coverage) || !nrow(coverage)) {
  die("property log contains no distribution coverage counters")
}
if (anyDuplicated(paste(coverage$coverage, coverage$metric))) {
  die("property log contains duplicate distribution coverage counters")
}

if (!file.exists(opt$coverage_requirements)) {
  die("coverage requirement manifest does not exist: {opt$coverage_requirements}")
}
root_path <- normalizePath(root, mustWork = TRUE)
requirements_path <- normalizePath(opt$coverage_requirements, mustWork = TRUE)
root_prefix <- paste0(root_path, .Platform$file.sep)
if (!startsWith(requirements_path, root_prefix)) {
  die("coverage requirement manifest must be inside the source repository")
}
requirements_relative <- substring(requirements_path, nchar(root_prefix) + 1L)
requirements_tracked <- suppressWarnings(system2(
  "git",
  duckvep_system2_quote(c(
    "-C",
    root,
    "ls-files",
    "--error-unmatch",
    "--",
    requirements_relative
  )),
  stdout = FALSE,
  stderr = FALSE
))
if (!identical(requirements_tracked, 0L)) {
  die(
    "coverage requirement manifest is not tracked at the recorded source revision: ",
    "{requirements_relative}"
  )
}
requirements <- utils::read.delim(
  requirements_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)
required_columns <- c("coverage", "metric", "requirement", "evidence")
if (!identical(names(requirements), required_columns)) {
  die("coverage requirement manifest has an unexpected schema")
}
if (!nrow(requirements) || anyNA(requirements) || any(!nzchar(as.matrix(requirements)))) {
  die("coverage requirement manifest contains an empty field")
}
allowed_requirements <- c("statistical_required", "fixed_witness")
if (any(!requirements$requirement %in% allowed_requirements)) {
  die("coverage requirement manifest contains an unknown requirement class")
}
requirement_key <- paste(requirements$coverage, requirements$metric, sep = "\r")
coverage_key <- paste(coverage$coverage, coverage$metric, sep = "\r")
if (anyDuplicated(requirement_key)) {
  die("coverage requirement manifest contains duplicate coverage/metric keys")
}
requirement_index <- match(requirement_key, coverage_key)
if (anyNA(requirement_index)) {
  missing <- paste(
    paste0(requirements$coverage[is.na(requirement_index)], "/",
           requirements$metric[is.na(requirement_index)]),
    collapse = ", "
  )
  die("property log omitted required coverage counters: {missing}")
}
unmanifested_index <- which(!coverage_key %in% requirement_key)
if (length(unmanifested_index)) {
  unmanifested <- paste(
    paste0(
      coverage$coverage[unmanifested_index],
      "/",
      coverage$metric[unmanifested_index]
    ),
    collapse = ", "
  )
  die("property log contains unmanifested coverage counters: {unmanifested}")
}
required_counts <- coverage$count[requirement_index]
missing_statistical <- requirements$requirement == "statistical_required" &
  required_counts < 1
if (any(missing_statistical)) {
  missing <- paste(
    paste0(requirements$coverage[missing_statistical], "/",
           requirements$metric[missing_statistical]),
    collapse = ", "
  )
  die("statistically required states received zero observations: {missing}")
}

zero_index <- which(coverage$count == 0)
if (length(zero_index)) {
  zero_requirement_index <- match(coverage_key[zero_index], requirement_key)
  allowed_zero <- !is.na(zero_requirement_index) &
    requirements$requirement[zero_requirement_index] == "fixed_witness"
  if (any(!allowed_zero)) {
    missing <- paste(
      paste0(coverage$coverage[zero_index[!allowed_zero]], "/",
             coverage$metric[zero_index[!allowed_zero]]),
      collapse = ", "
    )
    die("coverage counters reached zero without a fixed-witness declaration: {missing}")
  }
}

fixed_requirements <- requirements$requirement == "fixed_witness"
if (any(requirements$requirement == "statistical_required" & requirements$evidence != "-")) {
  die("statistically required counters must use '-' evidence")
}
if (any(fixed_requirements)) {
  evidence <- requirements$evidence[fixed_requirements]
  if (any(!grepl("^[A-Za-z_][A-Za-z0-9_]*$", evidence))) {
    die("fixed-witness evidence names must be C test identifiers")
  }
  property_source_path <- file.path(
    root,
    "test",
    "duckvep",
    "property",
    "duckvep_kernel_prop.c"
  )
  property_source <- paste(readLines(property_source_path, warn = FALSE), collapse = "\n")
  witness_tags <- paste0(
    "COVERAGE_WITNESS: ",
    requirements$coverage[fixed_requirements],
    "/",
    requirements$metric[fixed_requirements]
  )
  for (index in seq_along(evidence)) {
    test_name <- evidence[[index]]
    start <- regexpr(
      paste0(
        "(^|\\n)TEST[[:space:]]+",
        test_name,
        "[[:space:]]*\\([^)]*\\)[[:space:]]*\\{"
      ),
      property_source,
      perl = TRUE
    )
    if (start[[1L]] < 0L) {
      die("coverage requirement manifest names absent fixed witness: {test_name}")
    }
    block_start <- start[[1L]] + attr(start, "match.length")
    remainder <- substring(property_source, block_start)
    next_test <- regexpr(
      "\\nTEST[[:space:]]+[A-Za-z_][A-Za-z0-9_]*[[:space:]]*\\(",
      remainder,
      perl = TRUE
    )
    test_block <- if (next_test[[1L]] < 0L) {
      remainder
    } else {
      substring(remainder, 1L, next_test[[1L]] - 1L)
    }
    if (!grepl(witness_tags[[index]], test_block, fixed = TRUE)) {
      die(
        "fixed witness {test_name} lacks its direct coverage tag: ",
        "{witness_tags[[index]]}"
      )
    }
  }
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

coverage_rows <- data.frame(
  run_date = as.character(Sys.Date()),
  source_revision = opt$source_revision,
  configured_trials = as.numeric(opt$trials),
  seed = properties$reported_seed[[1L]],
  coverage = coverage$coverage,
  metric = coverage$metric,
  count = coverage$count,
  stringsAsFactors = FALSE
)
coverage_rows <- coverage_rows[
  order(coverage_rows$coverage, coverage_rows$metric),
  ,
  drop = FALSE
]

campaign_key <- function(data) {
  paste(
    data$source_revision,
    data$seed,
    format(
      as.numeric(data$configured_trials),
      scientific = FALSE,
      trim = TRUE
    ),
    sep = "\r"
  )
}

read_existing <- function(path, expected_names) {
  if (!file.exists(path)) {
    return(NULL)
  }
  existing <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = c(source_revision = "character", seed = "character")
  )
  if (!identical(names(existing), expected_names)) {
    die("history schema does not match: {path}")
  }
  existing
}

publish_histories <- function() {
  recover_stale_publication()
  if (!dir.create(lock_dir, showWarnings = FALSE)) {
    die("cannot acquire property-history publication lock: {lock_dir}")
  }
  writeLines(as.character(Sys.getpid()), owner_path, useBytes = TRUE)
  publication_complete <- FALSE
  on.exit({
    if (!isTRUE(publication_complete)) {
      restore_interrupted_publication(lock_dir)
    }
    unlink(lock_dir, recursive = TRUE)
  }, add = TRUE)

  old <- read_existing(history_path, names(rows))
  old_coverage <- read_existing(coverage_history_path, names(coverage_rows))
  new_run_key <- unique(campaign_key(rows))
  new_coverage_key <- unique(campaign_key(coverage_rows))
  if (length(new_run_key) != 1L || !identical(new_run_key, new_coverage_key)) {
    die("property and coverage rows do not describe one identical campaign")
  }
  if (!is.null(old) && new_run_key %in% campaign_key(old)) {
    die("property campaign already exists in immutable history: {new_run_key}")
  }
  if (!is.null(old_coverage) && new_run_key %in% campaign_key(old_coverage)) {
    die("property campaign already exists in immutable coverage history: {new_run_key}")
  }
  if (!is.null(old)) {
    rows <- rbind(old, rows)
  }
  rows <- rows[
    order(rows$run_date, rows$source_revision, rows$seed, rows$target),
    ,
    drop = FALSE
  ]
  if (!is.null(old_coverage)) {
    coverage_rows <- rbind(old_coverage, coverage_rows)
  }
  coverage_rows <- coverage_rows[
    order(
      coverage_rows$run_date,
      coverage_rows$source_revision,
      coverage_rows$seed,
      coverage_rows$coverage,
      coverage_rows$metric
    ),
    ,
    drop = FALSE
  ]

  property_tmp <- tempfile("duckvep-properties-", history_dir, ".csv")
  coverage_tmp <- tempfile(
    "duckvep-property-coverage-",
    coverage_dir,
    ".csv"
  )
  on.exit(unlink(c(property_tmp, coverage_tmp)), add = TRUE)

  utils::write.csv(rows, property_tmp, row.names = FALSE)
  utils::write.csv(coverage_rows, coverage_tmp, row.names = FALSE)
  staged_properties <- utils::read.csv(
    property_tmp,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = c(source_revision = "character", seed = "character")
  )
  staged_coverage <- utils::read.csv(
    coverage_tmp,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = c(source_revision = "character", seed = "character")
  )
  if (!identical(names(staged_properties), names(rows)) ||
      nrow(staged_properties) != nrow(rows) ||
      !identical(names(staged_coverage), names(coverage_rows)) ||
      nrow(staged_coverage) != nrow(coverage_rows)) {
    die("staged property histories failed round-trip validation")
  }

  had_properties <- file.exists(history_path)
  had_coverage <- file.exists(coverage_history_path)
  if (had_properties &&
      !file.copy(history_path, property_backup, overwrite = FALSE)) {
    die("cannot back up existing history: {history_path}")
  }
  if (had_coverage &&
      !file.copy(coverage_history_path, coverage_backup, overwrite = FALSE)) {
    die("cannot back up existing coverage history: {coverage_history_path}")
  }
  journal_tmp <- tempfile("journal-", lock_dir, ".rds")
  saveRDS(
    list(
      had_properties = had_properties,
      had_coverage = had_coverage,
      history_path = history_path,
      coverage_history_path = coverage_history_path
    ),
    journal_tmp
  )
  if (!file.rename(journal_tmp, journal_path)) {
    die("cannot publish recovery journal: {journal_path}")
  }

  # The temporary files are in the destination directories. POSIX rename(2)
  # therefore replaces each canonical ledger atomically. The journal and
  # backups recover the pair if the process exits between the two renames;
  # the lock prevents concurrent campaigns from overwriting one another.
  if (!file.rename(property_tmp, history_path)) {
    die("cannot publish history: {history_path}")
  }
  if (!file.rename(coverage_tmp, coverage_history_path)) {
    die("cannot publish coverage history: {coverage_history_path}")
  }
  # Removing the journal commits the pair. If the process dies after this
  # unlink, stale-lock cleanup preserves both published ledgers and discards
  # any remaining backups. Removing backups while the journal still existed
  # would make rollback incomplete after a process crash.
  if (unlink(journal_path) != 0L) {
    die("cannot commit property-history publication: {journal_path}")
  }
  publication_complete <- TRUE
  unlink(c(property_backup, coverage_backup))
}

check_source_tree(git_head)
publish_histories()

campaign_complete <- TRUE

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
cat(glue("coverage history -> {opt$coverage_history}"), "\n", sep = "")
}

main()
