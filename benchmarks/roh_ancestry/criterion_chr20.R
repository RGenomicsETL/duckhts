# Base-pair totals print as integers, not in scientific notation.
options(scipen = 999)
cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
metrics <- read.csv(file.path(cache, "chr20.child_metrics.csv"), stringsAsFactors = FALSE)
admixed <- c("ACB", "ASW", "CLM", "MXL", "PEL", "PUR")
controls <- c("YRI", "ESN", "CEU", "CHS")
# A population without truth length has no truth coverage: its coverage check is NA,
# not evaluable, and never counts as a pass. A failed count check still fails.
within_five_percent <- function(value, reference) {
  if (is.na(value) || is.na(reference)) return(NA)
  if (reference == 0) return(value == 0)
  abs(value - reference) <= 0.05 * abs(reference)
}
no_material_coverage_loss <- function(value, reference) {
  if (is.na(value) || is.na(reference)) return(NA)
  if (reference == 0) return(value == 0)
  value >= 0.95 * reference
}
populations <- sort(unique(metrics$Population))
criterion <- do.call(rbind, lapply(populations, function(population) {
  x <- metrics[metrics$Population == population, , drop = FALSE]
  average <- function(arm, field) mean(x[[field]][x$arm == arm])
  coverage <- function(arm) {
    selected <- x[x$arm == arm, , drop = FALSE]
    if (sum(selected$truth_bp) == 0) return(NA_real_)
    sum(selected$truth_covered_bp) / sum(selected$truth_bp)
  }
  truth_bp <- sum(x$truth_bp[x$arm == "ancestry_tuned"])
  ancestry_count <- average("ancestry_tuned", "unsupported_count_1mb")
  ancestry_coverage <- coverage("ancestry_tuned")
  if (population %in% admixed) {
    comparison_arms <- c("pooled", "single_population")
    count_pass <- vapply(comparison_arms, function(arm) {
      ancestry_count < average(arm, "unsupported_count_1mb")
    }, logical(1L))
    coverage_pass <- vapply(comparison_arms, function(arm) {
      no_material_coverage_loss(ancestry_coverage, coverage(arm))
    }, logical(1L))
    data.frame(Population = population, comparison = comparison_arms, truth_bp = truth_bp,
      ancestry_unsupported_count = ancestry_count,
      comparator_unsupported_count = vapply(comparison_arms, function(arm) {
        average(arm, "unsupported_count_1mb")
      }, numeric(1L)), ancestry_truth_coverage = ancestry_coverage,
      comparator_truth_coverage = vapply(comparison_arms, coverage, numeric(1L)),
      count_pass = count_pass, truth_coverage_pass = coverage_pass,
      criterion_pass = count_pass & coverage_pass)
  } else if (population %in% controls) {
    comparison <- "single_population"
    comparison_count <- average(comparison, "unsupported_count_1mb")
    comparison_coverage <- coverage(comparison)
    count_pass <- within_five_percent(ancestry_count, comparison_count)
    coverage_pass <- within_five_percent(ancestry_coverage, comparison_coverage)
    data.frame(Population = population, comparison = comparison, truth_bp = truth_bp,
      ancestry_unsupported_count = ancestry_count,
      comparator_unsupported_count = comparison_count,
      ancestry_truth_coverage = ancestry_coverage,
      comparator_truth_coverage = comparison_coverage,
      count_pass = count_pass, truth_coverage_pass = coverage_pass,
      criterion_pass = count_pass & coverage_pass)
  } else stop("unexpected population", call. = FALSE)
}))
utils::write.csv(criterion, file.path(cache, "chr20.criterion.csv"), row.names = FALSE)
utils::write.csv(criterion, "benchmarks/results/roh-ancestry/chr20/criterion_by_population.csv",
                 row.names = FALSE)
verdict <- if (any(criterion$criterion_pass %in% FALSE)) {
  "not held"
} else if (anyNA(criterion$criterion_pass)) {
  "not evaluable in every population"
} else {
  "held"
}
cat("chr20 criterion:", verdict, "\n")
print(criterion)
