# Comparison of two sets of runs on one chromosome.
#
# A run set is a data frame with `start` and `end`, one-based and inclusive, as
# the ROH macros return them. The runs of one decode do not overlap each other;
# roh_run_overlap_bases() relies on that and checks it.

roh_runs_are_disjoint <- function(runs) {
  if (nrow(runs) < 2L) return(TRUE)
  ordered <- runs[order(runs$start), , drop = FALSE]
  all(ordered$start[-1L] > ordered$end[-nrow(ordered)])
}

roh_run_bases <- function(runs) {
  if (!nrow(runs)) return(0)
  sum(as.numeric(runs$end) - as.numeric(runs$start) + 1)
}

roh_run_overlap_bases <- function(a, b) {
  if (!roh_runs_are_disjoint(a) || !roh_runs_are_disjoint(b)) {
    stop("the runs of one decode must not overlap each other", call. = FALSE)
  }
  if (!nrow(a) || !nrow(b)) return(0)
  first <- outer(as.numeric(a$start), as.numeric(b$start), pmax)
  last <- outer(as.numeric(a$end), as.numeric(b$end), pmin)
  sum(pmax(0, last - first + 1))
}

# One row of metrics for `test` against `reference` over `span_bases`.
# jaccard is the bases in both run sets over the bases in either; it is NA when
# neither set has a run.
roh_run_comparison <- function(test, reference, span_bases) {
  test_bases <- roh_run_bases(test)
  reference_bases <- roh_run_bases(reference)
  shared_bases <- roh_run_overlap_bases(test, reference)
  either_bases <- test_bases + reference_bases - shared_bases
  data.frame(
    reference_runs = nrow(reference), test_runs = nrow(test),
    reference_bases = reference_bases, test_bases = test_bases,
    shared_bases = shared_bases, either_bases = either_bases,
    reference_froh = reference_bases / span_bases, test_froh = test_bases / span_bases,
    froh_difference = (test_bases - reference_bases) / span_bases,
    jaccard = if (either_bases > 0) shared_bases / either_bases else NA_real_)
}

# The two run classes of the report: every run, and runs of at least
# `long_run_bases`. Each decode is filtered by its own run lengths.
roh_run_comparison_by_class <- function(test, reference, span_bases, long_run_bases) {
  long <- function(runs) {
    runs[as.numeric(runs$end) - as.numeric(runs$start) + 1 >= long_run_bases, , drop = FALSE]
  }
  rbind(
    data.frame(run_class = "all", roh_run_comparison(test, reference, span_bases)),
    data.frame(run_class = "long", roh_run_comparison(long(test), long(reference), span_bases)))
}

# The declared verdict. A class in which neither decode has a run agrees.
roh_run_verdict <- function(comparison, declaration) {
  minimum <- ifelse(comparison$run_class == "long", declaration$min_jaccard_long,
                    declaration$min_jaccard_all)
  froh_pass <- abs(comparison$froh_difference) <= declaration$max_froh_difference
  jaccard_pass <- is.na(comparison$jaccard) | comparison$jaccard >= minimum
  data.frame(comparison, froh_pass = froh_pass, jaccard_pass = jaccard_pass,
             pass = froh_pass & jaccard_pass)
}
