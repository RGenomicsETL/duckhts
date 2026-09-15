#!/usr/bin/env Rscript

# Pinned Somalier v0.3.4 helper outputs are fixture data, not the R oracle.
# This tests helper-level arithmetic only; it does not certify the CLI.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("usage: Rscript somalier_v034_differential.R [native.tsv]")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]])
root <- normalizePath(file.path(dirname(script), "../.."))
cases <- utils::read.delim(file.path(root, "test/data/somalier_v034_helper_cases.tsv"),
  colClasses = "character", check.names = FALSE, na.strings = character())
receipt <- utils::read.delim(file.path(root, "test/data/somalier_v034_helper_receipt.tsv"),
  colClasses = "character", check.names = FALSE, na.strings = character())

case_columns <- c("case_id", "metric", "input", "upstream_value",
  "independent_value", "relation", "abs_tolerance", "expected_native_status",
  "expected_usable_sites", "expected_hom_a", "expected_hom_b")
required_cases <- c("charr_depth_100", "charr_depth_200", "charr_depth_6000",
  "charr_deep_usable", "charr_deep_estimate", "charr_control_estimate",
  "pair_two_sites_usable", "pair_search_alpha", "pair_search_likelihood",
  "pair_likelihood_0", "pair_likelihood_02", "pair_likelihood_04",
  "pair_likelihood_05", "pair_likelihood_1")
if (!identical(names(cases), case_columns) ||
    !identical(cases$case_id, required_cases) ||
    !identical(names(receipt), c("field", "value")) ||
    anyDuplicated(receipt$field)) {
  stop("Somalier helper fixture is incomplete or reordered")
}
fields <- stats::setNames(receipt$value, receipt$field)
if (!identical(fields[["commit"]],
      "ff58fdade8f4f8293d904f10e0a4a13f1fac808d") ||
    !identical(fields[["scope"]],
      "helper-level source-derived witness; not full CLI conformance") ||
    !identical(fields[["source_archive_sha256"]],
      "acd2dc11be6051c80d15628703a9419965cea3fe7a0563b47e14743e7ac6339e") ||
    !identical(fields[["source_archive_bytes"]], "1208978") ||
    !identical(fields[["witness_driver_sha256"]],
      "0f3bbd0c7088e98859221172c73d1399eefb2db33a7afc07052277726d7f0f32") ||
    !identical(fields[["regeneration"]],
      "Rscript test/scripts/somalier_v034_regenerate.R [staged-source-archive]") ||
    !identical(fields[["upstream_result_sha256"]],
      "34e8dff0151c4d556aaaf3c4b3c03bfabdc55e00861f60b058d74303ea3c4aeb") ||
    !identical(fields[["upstream_helper_exit"]], "0") ||
    !identical(fields[["discovery_comparison_exit"]], "1")) {
  stop("Somalier helper provenance or failure status changed")
}
if (!requireNamespace("digest", quietly = TRUE)) {
  stop("digest is required to verify pinned helper output bytes")
}
upstream_bytes <- paste0(paste(c("metric\tinput\tvalue",
  paste(cases$metric, cases$input, cases$upstream_value, sep = "\t")),
  collapse = "\n"), "\n")
if (!identical(digest::digest(charToRaw(upstream_bytes), algo = "sha256",
      serialize = FALSE), fields[["upstream_result_sha256"]])) {
  stop("Pinned helper outputs do not match the retained Nim result SHA-256")
}

parse_value <- function(value) {
  if (identical(value, "NaN")) return(NaN)
  parsed <- suppressWarnings(as.numeric(value))
  if (length(parsed) != 1L || !is.finite(parsed)) {
    stop("Nonfinite or invalid numerical fixture value: ", value)
  }
  parsed
}
near <- function(actual, expected, tolerance) {
  if (is.nan(actual) || is.nan(expected)) return(is.nan(actual) && is.nan(expected))
  is.finite(actual) && is.finite(expected) && abs(actual - expected) <= tolerance
}

binomial_limit <- function(depth, rate = 0.12, tail_alpha = 0.002) {
  candidates <- 0:depth
  tails <- stats::pbinom(candidates - 1L, size = depth, prob = rate,
    lower.tail = FALSE)
  max(candidates[tails >= tail_alpha])
}
charr_one <- function(a, b, other, frequency_b) {
  depth <- a + b
  minor <- min(a, b)
  usable <- depth >= 15 && other / (depth + other) <= 0.04 &&
    minor <= binomial_limit(depth) && a != b
  genotype <- if (usable) if (a > b) "hom_a" else "hom_b" else "unknown"
  infiltrating_frequency <- if (genotype == "hom_a") frequency_b else
    if (genotype == "hom_b") 1 - frequency_b else 0
  usable <- usable && infiltrating_frequency > 1e-5
  list(status = if (usable) "ok" else "no_evidence",
    usable_sites = as.integer(usable),
    hom_a = as.integer(usable && genotype == "hom_a"),
    hom_b = as.integer(usable && genotype == "hom_b"),
    estimate = if (usable) minor / (infiltrating_frequency * depth) else NaN)
}
deep <- charr_one(4000, 2000, 0, 0.5)
control <- charr_one(199, 1, 0, 0.25)

receiver_a <- c(160, 40)
receiver_b <- c(40, 160)
anchor_gt <- c(0, 2)
population_b_af <- c(0.5, 0.5)
pair_log_likelihood <- function(alpha) {
  sum(vapply(seq_along(anchor_gt), function(i) {
    af <- population_b_af[[i]]
    prior <- c((1 - af)^2, 2 * af * (1 - af), af^2)
    probability <- 0.002 + (1 - 2 * 0.002) *
      ((1 - alpha) * anchor_gt[[i]] / 2 + alpha * (0:2) / 2)
    probability <- pmax(1e-10, pmin(1 - 1e-10, probability))
    terms <- log(prior) + receiver_b[[i]] * log(probability) +
      receiver_a[[i]] * log1p(-probability)
    peak <- max(terms)
    peak + log(sum(exp(terms - peak)))
  }, numeric(1)))
}
upstream_alpha <- parse_value(cases$upstream_value[
  cases$case_id == "pair_search_alpha"])
local <- stats::optimize(pair_log_likelihood, c(0.35, 0.45),
  maximum = TRUE, tol = 1e-12)
independent <- c(
  charr_depth_100 = binomial_limit(100),
  charr_depth_200 = binomial_limit(200),
  charr_depth_6000 = binomial_limit(6000),
  charr_deep_usable = deep$usable_sites,
  charr_deep_estimate = deep$estimate,
  charr_control_estimate = control$estimate,
  pair_two_sites_usable = length(anchor_gt),
  pair_search_alpha = local$maximum,
  pair_search_likelihood = pair_log_likelihood(upstream_alpha),
  pair_likelihood_0 = pair_log_likelihood(0),
  pair_likelihood_02 = pair_log_likelihood(0.2),
  pair_likelihood_04 = pair_log_likelihood(0.4),
  pair_likelihood_05 = pair_log_likelihood(0.5),
  pair_likelihood_1 = pair_log_likelihood(1))

for (i in seq_len(nrow(cases))) {
  case <- cases[i, ]
  upstream <- parse_value(case$upstream_value)
  expected <- parse_value(case$independent_value)
  computed <- independent[[case$case_id]]
  tolerance <- parse_value(case$abs_tolerance)
  if (!near(computed, expected, tolerance)) {
    stop("Independent R result changed for ", case$case_id)
  }
  if (case$relation == "equal") {
    if (!near(upstream, computed, tolerance)) {
      stop("Pinned helper and independent R disagree at ", case$case_id)
    }
  } else if (case$relation == "different") {
    if (near(upstream, computed, tolerance)) {
      stop("Retained helper difference disappeared at ", case$case_id)
    }
  } else if (case$relation == "better_feasible") {
    if (abs(upstream - computed) < 1e-3 ||
        local$objective <= pair_log_likelihood(upstream) + 0.5 ||
        pair_log_likelihood(0.4) <= pair_log_likelihood(upstream) + 0.5) {
      stop("Retained pair-search counterexample disappeared")
    }
  } else {
    stop("Unknown helper comparison relation at ", case$case_id)
  }
}
if (deep$status != "no_evidence" || deep$usable_sites != 0L ||
    deep$hom_a != 0L || deep$hom_b != 0L || !is.nan(deep$estimate) ||
    control$status != "ok" || control$usable_sites != 1L ||
    control$hom_a != 1L || control$hom_b != 0L ||
    length(anchor_gt) != 2L) {
  stop("Independent eligibility, named status, or denominators changed")
}
expected_status <- c(rep("ok", 3L), "no_evidence", "no_evidence",
  rep("ok", 9L))
expected_usable <- c(rep("NA", 3L), "0", "0", "1", rep("2", 8L))
expected_hom_a <- c(rep("NA", 3L), "0", "0", "1", rep("NA", 8L))
expected_hom_b <- c(rep("NA", 3L), "0", "0", "0", rep("NA", 8L))
if (!identical(cases$expected_native_status, expected_status) ||
    !identical(cases$expected_usable_sites, expected_usable) ||
    !identical(cases$expected_hom_a, expected_hom_a) ||
    !identical(cases$expected_hom_b, expected_hom_b)) {
  stop("Native status or denominator contract was weakened")
}

if (length(args) == 1L) {
  native <- utils::read.delim(args[[1L]], colClasses = "character",
    check.names = FALSE, na.strings = character())
  native_columns <- c("case_id", "metric", "value", "status", "usable_sites",
    "usable_hom_a", "usable_hom_b", "evaluations")
  if (!identical(names(native), native_columns) ||
      !identical(native$case_id, cases$case_id) ||
      !identical(native$metric, cases$metric)) {
    stop("Native witness rows are missing, extra, or reordered")
  }
  for (i in seq_len(nrow(native))) {
    tolerance <- parse_value(cases$abs_tolerance[[i]])
    if (!near(parse_value(native$value[[i]]), independent[[native$case_id[[i]]]],
        tolerance) ||
        native$status[[i]] != cases$expected_native_status[[i]]) {
      stop("Native value or status disagrees at ", native$case_id[[i]])
    }
    for (column in c("usable_sites", "usable_hom_a", "usable_hom_b")) {
      expected_column <- paste0("expected_", sub("^usable_hom_", "hom_", column))
      expected <- cases[[expected_column]][[i]]
      if (expected != "NA" && native[[column]][[i]] != expected) {
        stop("Native denominator disagrees at ", native$case_id[[i]],
          " (", column, ")")
      }
    }
    evaluations <- native$evaluations[[i]]
    if (evaluations != "NA" &&
        (!grepl("^[0-9]+$", evaluations) ||
         (native$case_id[[i]] == "pair_search_alpha" &&
          as.numeric(evaluations) == 0))) {
      stop("Native evaluation status is invalid at ", native$case_id[[i]])
    }
  }
}

cat("Somalier v0.3.4 helper differential: OK (", nrow(cases),
  " retained cases; pinned helper vs independent R", if (length(args))
  "; current native witness" else "", ")\n", sep = "")
cat("Retained differences: depth-6000 CHARR helper threshold 6000 vs R 793;",
  " high-depth usable sites 1 vs 0; search alpha", format(upstream_alpha),
  "has lower likelihood than feasible alpha 0.4.\n")
