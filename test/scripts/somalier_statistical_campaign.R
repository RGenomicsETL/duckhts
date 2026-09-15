#!/usr/bin/env Rscript

# Numerical conformance of the DuckHTS Somalier kernel against R's binomial
# distribution and an explicit, independently written three-genotype score.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !file.exists(args[[1L]])) {
  stop("usage: Rscript somalier_statistical_campaign.R NATIVE_TEST_BINARY")
}
native_binary <- normalizePath(args[[1L]])
seed <- 20260915L
set.seed(seed)
work_dir <- tempfile("somalier-statistical-", tmpdir = dirname(tempdir()))
dir.create(work_dir)
completed <- FALSE
on.exit(if (completed) unlink(work_dir, recursive = TRUE), add = TRUE)
failures <- data.frame(mode = character(), case_id = character(), reason = character(),
  input = character(), native = character(), expected = character())

row_text <- function(row) paste(names(row), as.character(row), sep = "=", collapse = ";")
record_failure <- function(mode, case_id, reason, input, native, expected) {
  failures <<- rbind(failures, data.frame(mode = mode, case_id = case_id,
    reason = reason, input = row_text(input), native = row_text(native),
    expected = row_text(expected)))
}
real_text <- function(value) sprintf("%.17g", value)
whole_text <- function(value) sprintf("%.0f", value)
numeric_value <- function(value) suppressWarnings(as.numeric(value))
same_number <- function(actual, expected, tolerance) {
  if (is.nan(expected)) return(is.nan(actual))
  is.finite(actual) && is.finite(expected) && abs(actual - expected) <= tolerance
}

run_native <- function(mode, input, expected_columns) {
  input_path <- file.path(work_dir, paste0(mode, "-input.tsv"))
  output_path <- file.path(work_dir, paste0(mode, "-output.tsv"))
  error_path <- file.path(work_dir, paste0(mode, "-stderr.txt"))
  utils::write.table(input, input_path, sep = "\t", row.names = FALSE,
    quote = FALSE, na = "NA")
  status <- system2(native_binary, c("--campaign", mode),
    input = readLines(input_path, warn = FALSE), stdout = output_path,
    stderr = error_path)
  if (status != 0L) {
    stop("Native ", mode, " campaign exited ", status,
      "; retained input/output/stderr at ", work_dir)
  }
  output <- utils::read.delim(output_path, colClasses = "character",
    check.names = FALSE, na.strings = character())
  if (!identical(names(output), expected_columns) ||
      nrow(output) != nrow(input) ||
      anyDuplicated(output$case_id) ||
      !identical(output$case_id, input$case_id)) {
    stop("Native ", mode, " output omitted, reordered, or duplicated cases; ",
      "retained input/output at ", work_dir)
  }
  output
}

# qbinom supplies a candidate; exact pbinom survival checks decide the
# discrete maximum. This does not copy the kernel's beta implementation.
survival <- function(depth, minor_rate, k) {
  if (depth == 1 && k == 1) return(minor_rate)
  # This small rational witness has an exact binary64 tail. Keep it independent
  # of pbinom so adjacent representable cutoffs test the strict comparison.
  if (depth == 3 && minor_rate == 0.25) {
    return(c(1, 0.578125, 0.15625, 0.015625, 0)[pmin(pmax(k, 0), 4) + 1])
  }
  stats::pbinom(k - 1, size = depth, prob = minor_rate, lower.tail = FALSE)
}
binomial_expected <- function(depth, minor_rate, tail_alpha, max_depth) {
  if (depth > max_depth || depth > 1000000) {
    return(list(status = "limit exceeded", max_minor = NA_integer_))
  }
  if (max_depth == 0 || !is.finite(minor_rate) || minor_rate < 0 ||
      minor_rate > 1 || !is.finite(tail_alpha) || tail_alpha <= 0 ||
      tail_alpha >= 1) {
    return(list(status = "invalid argument", max_minor = NA_integer_))
  }
  k <- as.integer(stats::qbinom(1 - tail_alpha, depth, minor_rate))
  while (k > 0L && survival(depth, minor_rate, k) < tail_alpha) k <- k - 1L
  while (k < depth && survival(depth, minor_rate, k + 1L) >= tail_alpha) {
    k <- k + 1L
  }
  list(status = "ok", max_minor = k)
}

depths <- unique(as.integer(round(exp(seq(log(1), log(1000000), length.out = 32L)))))
rates <- c(0, 1e-6, 0.01, 0.05, 0.12, 0.49, 0.75, 0.99, 1)
tails <- c(1e-8, 1e-5, 0.001, 0.002, 0.05, 0.5, 0.95, 1 - 1e-8)
binomial_grid <- expand.grid(depth = depths, minor_rate = rates,
  tail_alpha = tails, KEEP.OUT.ATTRS = FALSE)
binomial_random <- data.frame(
  depth = sample(depths, 160L, replace = TRUE),
  minor_rate = stats::runif(160L, 1e-5, 1 - 1e-5),
  tail_alpha = exp(stats::runif(160L, log(1e-8), log(0.95))))
adjacent <- data.frame(depth = integer(), minor_rate = numeric(),
  tail_alpha = numeric())
adjacent <- rbind(adjacent, data.frame(
  depth = rep(3L, 3L), minor_rate = rep(0.25, 3L),
  tail_alpha = 0.578125 + c(-.Machine$double.eps / 2, 0,
    .Machine$double.eps / 2)))
adjacent <- rbind(adjacent, data.frame(
  depth = rep(1L, 3L), minor_rate = rep(0.05, 3L),
  tail_alpha = 0.05 + c(-.Machine$double.eps / 32, 0,
    .Machine$double.eps / 32)))
for (depth in c(15, 100, 1000, 6000, 100000, 1000000)) {
  for (rate in c(0.05, 0.12, 0.49)) {
    k <- max(1, min(depth, round(depth * rate)))
    target <- survival(depth, rate, k)
    adjacent <- rbind(adjacent, data.frame(depth = rep(depth, 4L),
      minor_rate = rep(rate, 4L),
      tail_alpha = target * c(1 - 1e-8, 1 - 1e-12, 1 + 1e-12, 1 + 1e-8)))
  }
}
binomial_invalid <- data.frame(depth = c(1000001, 101, 100, 100, 100, 100),
  minor_rate = c(0.12, 0.12, -0.01, 1.01, 0.12, 0.12),
  tail_alpha = c(0.002, 0.002, 0.002, 0.002, 0, 1),
  max_depth = c(1000001, 100, 100, 100, 100, 100))
binomial_cases <- rbind(
  transform(binomial_grid, max_depth = 1000000),
  transform(binomial_random, max_depth = 1000000),
  transform(adjacent, max_depth = 1000000), binomial_invalid)
binomial_cases$case_id <- sprintf("binom_%05d", seq_len(nrow(binomial_cases)))
binomial_input <- binomial_cases[c("case_id", "depth", "minor_rate",
  "tail_alpha", "max_depth")]
binomial_input$minor_rate <- real_text(binomial_input$minor_rate)
binomial_input$tail_alpha <- real_text(binomial_input$tail_alpha)
binomial_input$depth <- whole_text(binomial_input$depth)
binomial_input$max_depth <- whole_text(binomial_input$max_depth)
binomial_native <- run_native("binomial", binomial_input,
  c("case_id", "status", "max_minor"))
for (i in seq_len(nrow(binomial_cases))) {
  input <- binomial_cases[i, ]
  observed <- binomial_native[i, ]
  expected <- binomial_expected(input$depth, input$minor_rate,
    input$tail_alpha, input$max_depth)
  reason <- character()
  if (observed$status != expected$status) reason <- c(reason, "status")
  if (expected$status == "ok") {
    reported <- numeric_value(observed$max_minor)
    if (!is.finite(reported) || reported != expected$max_minor) {
      reason <- c(reason, "max_minor")
    }
    if (is.finite(reported)) {
      at_k <- survival(input$depth, input$minor_rate, reported)
      at_next <- survival(input$depth, input$minor_rate, reported + 1)
      if (at_k < input$tail_alpha || at_next >= input$tail_alpha) {
        reason <- c(reason, "k/k+1 survival")
      }
      expected$tail_at_native_k <- at_k
      expected$tail_at_native_k_plus_one <- at_next
    }
  } else if (observed$max_minor != "NA") {
    reason <- c(reason, "non-success result was not NA")
  }
  if (length(reason)) record_failure("binomial", input$case_id,
    paste(reason, collapse = ", "), input, observed, expected)
}

eligibility_expected <- function(input) {
  depth <- input$allele_a + input$allele_b
  empty <- list(classification_status = "ok", genotype = -1L,
    receiver_status = "ok", receiver_usable = 0L,
    charr_status = "no evidence", charr_estimate = NaN,
    charr_usable_sites = 0L, charr_hom_a = 0L, charr_hom_b = 0L)
  if (input$available > 1 || input$min_depth == 0 ||
      input$max_depth < input$min_depth || input$hom_minor_rate <= 0 ||
      input$hom_minor_rate >= 0.5 || input$hom_tail_alpha <= 0 ||
      input$hom_tail_alpha >= 1 || input$population_b_af < 0 ||
      input$population_b_af > 1) {
    empty$classification_status <- "invalid argument"
    empty$receiver_status <- "invalid argument"
    empty$charr_status <- "invalid argument"
    return(empty)
  }
  if (input$available == 0) return(empty)
  if (depth > input$max_depth) {
    empty$classification_status <- "limit exceeded"
    empty$receiver_status <- "limit exceeded"
    empty$charr_status <- "limit exceeded"
    return(empty)
  }
  other_ok <- input$other / (depth + input$other) <= 0.04
  receiver_ok <- depth >= input$min_depth && other_ok
  empty$receiver_usable <- as.integer(receiver_ok)
  if (!receiver_ok) return(empty)
  threshold <- binomial_expected(depth, input$hom_minor_rate,
    input$hom_tail_alpha, input$max_depth)$max_minor
  minor <- min(input$allele_a, input$allele_b)
  if (minor > threshold || input$allele_a == input$allele_b) return(empty)
  empty$genotype <- if (input$allele_a > input$allele_b) 0L else 2L
  infiltrating_frequency <- if (empty$genotype == 0L)
    input$population_b_af else 1 - input$population_b_af
  if (infiltrating_frequency <= 1e-5) return(empty)
  empty$charr_status <- "ok"
  empty$charr_usable_sites <- 1L
  empty$charr_hom_a <- as.integer(empty$genotype == 0L)
  empty$charr_hom_b <- as.integer(empty$genotype == 2L)
  empty$charr_estimate <- minor / (infiltrating_frequency * depth)
  empty
}

eligibility_cases <- data.frame()
for (depth in c(14, 15, 16, 99, 100, 101, 5999, 6000, 6001)) {
  threshold <- binomial_expected(depth, 0.12, 0.002, 1000000)$max_minor
  minors <- unique(pmax(0L, pmin(depth %/% 2L,
    c(threshold - 1L, threshold, threshold + 1L, depth %/% 2L))))
  others <- unique(c(floor(depth / 24), floor(depth / 24) + 1L))
  for (minor in minors) {
    for (other in others) {
      for (frequency in c(0, 1e-5 * (1 - 1e-6), 1e-5 * (1 + 1e-6),
                          0.25, 0.5, 1 - 1e-5 * (1 + 1e-6), 1)) {
        eligibility_cases <- rbind(eligibility_cases, data.frame(
          allele_a = depth - minor, allele_b = minor, other = other,
          available = 1, population_b_af = frequency, min_depth = 15,
          max_depth = 1000000, hom_minor_rate = 0.12, hom_tail_alpha = 0.002))
      }
    }
  }
}
eligibility_controls <- data.frame(
  allele_a = c(20, 0, 20, 20, 14, 15, 20, 20, 20),
  allele_b = c(0, 20, 0, 0, 0, 0, 0, 0, 0),
  other = c(0, 0, 0, 1, 0, 0, 0, 0, 0),
  available = c(0, 1, 1, 1, 1, 1, 1, 2, 1),
  population_b_af = c(0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
                       0.25, 0.25, 0.25),
  min_depth = c(15, 15, 15, 15, 15, 15, 15, 15, 15),
  max_depth = c(1000000, 1000000, 19, 1000000, 1000000,
                1000000, 1000000, 1000000, 1000000),
  hom_minor_rate = c(0.12, 0.12, 0.12, 0.12, 0.12,
                     0.12, 0.12, 0.12, 0.49),
  hom_tail_alpha = c(0.002, 0.002, 0.002, 0.002, 0.002,
                     0.002, 0.002, 0.002, 0.95))
eligibility_exact_tail <- data.frame(
  allele_a = rep(2, 3), allele_b = rep(1, 3), other = rep(0, 3),
  available = rep(1, 3), population_b_af = rep(0.25, 3),
  min_depth = rep(1, 3), max_depth = rep(1000000, 3),
  hom_minor_rate = rep(0.25, 3),
  hom_tail_alpha = 0.578125 + c(-.Machine$double.eps / 2, 0,
    .Machine$double.eps / 2)
)
eligibility_cases <- rbind(eligibility_cases, eligibility_controls,
  eligibility_exact_tail)
eligibility_cases$case_id <- sprintf("elig_%04d", seq_len(nrow(eligibility_cases)))
eligibility_input <- eligibility_cases[c("case_id", "allele_a", "allele_b",
  "other", "available", "population_b_af", "min_depth", "max_depth",
  "hom_minor_rate", "hom_tail_alpha")]
for (column in c("population_b_af", "hom_minor_rate", "hom_tail_alpha")) {
  eligibility_input[[column]] <- real_text(eligibility_input[[column]])
}
for (column in c("allele_a", "allele_b", "other", "available", "min_depth",
                 "max_depth")) eligibility_input[[column]] <- whole_text(eligibility_input[[column]])
eligibility_native <- run_native("eligibility", eligibility_input,
  c("case_id", "classification_status", "genotype", "receiver_status",
    "receiver_usable", "charr_status", "charr_estimate",
    "charr_usable_sites", "charr_hom_a", "charr_hom_b"))
for (i in seq_len(nrow(eligibility_cases))) {
  input <- eligibility_cases[i, ]
  observed <- eligibility_native[i, ]
  expected <- eligibility_expected(input)
  reason <- character()
  for (column in c("classification_status", "receiver_status", "charr_status")) {
    if (observed[[column]] != expected[[column]]) reason <- c(reason, column)
  }
  for (column in c("genotype", "receiver_usable", "charr_usable_sites",
                   "charr_hom_a", "charr_hom_b")) {
    reported <- numeric_value(observed[[column]])
    if (!is.finite(reported) || reported != expected[[column]]) {
      reason <- c(reason, column)
    }
  }
  if (!same_number(numeric_value(observed$charr_estimate),
      expected$charr_estimate, 1e-12)) reason <- c(reason, "charr_estimate")
  if (length(reason)) record_failure("eligibility", input$case_id,
    paste(reason, collapse = ", "), input, observed, expected)
}

matched_score <- function(input, alpha) {
  total <- 0
  for (site in seq_len(input$site_count)) {
    usable <- input[[paste0("receiver_usable_", site)]]
    anchor <- input[[paste0("anchor_gt_", site)]]
    if (usable == 0 || !(anchor %in% c(0, 2))) next
    a <- input[[paste0("receiver_a_", site)]]
    b <- input[[paste0("receiver_b_", site)]]
    af <- input[[paste0("population_b_af_", site)]]
    af <- pmax(input$min_prior_frequency,
      pmin(1 - input$min_prior_frequency, af))
    prior <- c((1 - af)^2, 2 * af * (1 - af), af^2)
    latent <- (1 - alpha) * (anchor / 2) + alpha * (0:2) / 2
    observed_b <- latent * (1 - input$error_rate) +
      (1 - latent) * input$error_rate
    observed_b <- pmax(input$min_probability,
      pmin(1 - input$min_probability, observed_b))
    terms <- log(prior) + b * log(observed_b) + a * log1p(-observed_b)
    peak <- max(terms)
    total <- total + peak + log(sum(exp(terms - peak)))
  }
  total
}
matched_usable <- function(input) {
  sum(vapply(seq_len(input$site_count), function(site) {
    input[[paste0("receiver_usable_", site)]] == 1 &&
      input[[paste0("anchor_gt_", site)]] %in% c(0, 2)
  }, logical(1)))
}
matched_row <- function(case_id, operation, a1, b1, anchor1, af1,
                        a2, b2, anchor2, af2, alpha, usable1 = 1,
                        usable2 = 1, max_depth = 1000000,
                        max_evaluations = 4096) {
  data.frame(case_id, operation, site_count = 2,
    receiver_a_1 = a1, receiver_b_1 = b1, receiver_usable_1 = usable1,
    anchor_gt_1 = anchor1, population_b_af_1 = af1,
    receiver_a_2 = a2, receiver_b_2 = b2, receiver_usable_2 = usable2,
    anchor_gt_2 = anchor2, population_b_af_2 = af2,
    alpha, max_depth, error_rate = 0.002, min_probability = 1e-10,
    min_prior_frequency = 1e-6, alpha_min = 0, alpha_max = 1,
    grid_step = 0.01, refine_tolerance = 1e-10, max_evaluations)
}
matched_cases <- data.frame()
frequency_choices <- c(0, 1e-6 * (1 - 1e-4), 1e-6 * (1 + 1e-4),
  0.01, 0.2, 0.5, 0.8, 1 - 1e-6 * (1 + 1e-4), 1)
alpha_choices <- c(0, 1e-5, 0.2, 0.5, 0.8, 1 - 1e-5, 1)
for (trial in seq_len(240L)) {
  true_alpha <- sample(alpha_choices, 1L)
  frequencies <- sample(frequency_choices, 2L, replace = TRUE)
  anchors <- sample(c(0L, 2L), 2L, replace = TRUE)
  depth <- sample(c(15L, 16L, 100L, 200L, 1000L), 2L, replace = TRUE)
  clean_b <- anchors / 2
  probability_b <- 0.002 + (1 - 2 * 0.002) *
    ((1 - true_alpha) * clean_b + true_alpha * frequencies)
  b <- stats::rbinom(2L, depth, probability_b)
  a <- depth - b
  evaluated_alpha <- sample(alpha_choices, 1L)
  matched_cases <- rbind(matched_cases,
    matched_row(sprintf("matched_%04d_ll", trial), "likelihood",
      a[1], b[1], anchors[1], frequencies[1],
      a[2], b[2], anchors[2], frequencies[2], evaluated_alpha),
    matched_row(sprintf("matched_%04d_fit", trial), "fit",
      a[1], b[1], anchors[1], frequencies[1],
      a[2], b[2], anchors[2], frequencies[2], 0))
}
matched_cases <- rbind(matched_cases,
  # These two rows swap actual receiver and anchor sample roles. The anchor
  # sample's clean counts are 100/0 and 0/100 at the same ordered sites.
  matched_row("direction_receiver_anchor", "likelihood", 90, 10, 0, 0.2,
    10, 90, 2, 0.8, 0.2),
  matched_row("direction_anchor_receiver", "likelihood", 100, 0, 0, 0.2,
    0, 100, 2, 0.8, 0.2),
  matched_row("orientation_original", "likelihood", 90, 10, 0, 0.2,
    10, 90, 2, 0.8, 0.2),
  matched_row("orientation_flipped", "likelihood", 10, 90, 2, 0.8,
    90, 10, 0, 0.2, 0.2),
  matched_row("matched_no_evidence", "likelihood", 90, 10, -1, 0.5,
    10, 90, -1, 0.5, 0.2),
  matched_row("matched_fit_no_evidence", "fit", 90, 10, -1, 0.5,
    10, 90, -1, 0.5, 0),
  matched_row("matched_one_usable", "likelihood", 90, 10, 0, 0.5,
    10, 90, 2, 0.5, 0.2, usable2 = 0),
  matched_row("matched_endpoint_zero", "fit", 100, 0, 0, 0.5,
    0, 100, 2, 0.5, 0),
  matched_row("matched_endpoint_one", "fit", 50, 50, 0, 0.5,
    50, 50, 2, 0.5, 0),
  matched_row("matched_depth_limit", "fit", 90, 10, 0, 0.5,
    10, 90, 2, 0.5, 0, max_depth = 99),
  matched_row("matched_invalid_usable", "fit", 90, 10, 0, 0.5,
    10, 90, 2, 0.5, 0, usable1 = 2),
  matched_row("matched_invalid_zero_depth", "fit", 0, 0, 0, 0.5,
    10, 90, 2, 0.5, 0),
  matched_row("matched_invalid_frequency", "fit", 90, 10, 0, -0.01,
    10, 90, 2, 0.5, 0),
  matched_row("matched_fit_eval_limit", "fit", 90, 10, 0, 0.5,
    10, 90, 2, 0.5, 0, max_evaluations = 10))
matched_input <- matched_cases
for (column in c("population_b_af_1", "population_b_af_2", "alpha",
                 "error_rate", "min_probability", "min_prior_frequency",
                 "alpha_min", "alpha_max", "grid_step", "refine_tolerance")) {
  matched_input[[column]] <- real_text(matched_input[[column]])
}
for (column in c("site_count", "receiver_a_1", "receiver_b_1",
                 "receiver_usable_1", "anchor_gt_1", "receiver_a_2",
                 "receiver_b_2", "receiver_usable_2", "anchor_gt_2",
                 "max_depth", "max_evaluations")) {
  matched_input[[column]] <- whole_text(matched_input[[column]])
}
matched_native <- run_native("matched", matched_input,
  c("case_id", "operation", "status", "alpha", "log_likelihood",
    "usable_sites", "evaluations"))
for (i in seq_len(nrow(matched_cases))) {
  input <- matched_cases[i, ]
  observed <- matched_native[i, ]
  usable <- matched_usable(input)
  expected_status <- if (input$case_id %in% c("matched_depth_limit",
      "matched_fit_eval_limit")) "limit exceeded" else if (
      input$case_id %in% c("matched_invalid_usable",
        "matched_invalid_zero_depth", "matched_invalid_frequency"))
      "invalid argument" else if (usable == 0L) "no evidence" else "ok"
  reason <- character()
  expected <- list(status = expected_status, usable_sites = usable,
    alpha = input$alpha, fixed_score = if (usable) matched_score(input, input$alpha)
      else NaN)
  if (input$case_id %in% c("matched_depth_limit", "matched_invalid_usable",
      "matched_invalid_zero_depth", "matched_invalid_frequency")) {
    # Validation rejects the first malformed/over-depth receiver before a
    # usable-site count is published. The failure result retains zero.
    expected$usable_sites <- 0L
  }
  if (observed$operation != input$operation ||
      observed$status != expected_status) reason <- c(reason, "operation/status")
  reported_usable <- numeric_value(observed$usable_sites)
  if (!is.finite(reported_usable) ||
      reported_usable != expected$usable_sites) {
    reason <- c(reason, "usable_sites")
  }
  if (input$operation == "likelihood") {
    reported_evaluations <- numeric_value(observed$evaluations)
    if (!same_number(numeric_value(observed$alpha), input$alpha, 0) ||
        !same_number(numeric_value(observed$log_likelihood),
          expected$fixed_score, 1e-9) ||
        !is.finite(reported_evaluations) || reported_evaluations != 0) {
      reason <- c(reason, "fixed-alpha score/evaluations")
    }
  } else if (expected_status == "ok") {
    fitted_alpha <- numeric_value(observed$alpha)
    fitted_score <- numeric_value(observed$log_likelihood)
    evaluations <- numeric_value(observed$evaluations)
    grid <- seq(input$alpha_min, input$alpha_max, by = input$grid_step)
    grid_scores <- vapply(grid, function(alpha) matched_score(input, alpha),
      numeric(1))
    expected$grid_max <- max(grid_scores)
    expected$score_at_native_alpha <- if (is.finite(fitted_alpha))
      matched_score(input, fitted_alpha) else NaN
    if (!is.finite(fitted_alpha) || fitted_alpha < input$alpha_min ||
        fitted_alpha > input$alpha_max ||
        !same_number(fitted_score, expected$score_at_native_alpha, 1e-8) ||
        fitted_score < expected$grid_max - 1e-8 ||
        !is.finite(evaluations) || evaluations < length(grid) ||
        evaluations > input$max_evaluations) {
      reason <- c(reason, "fit score/grid/evaluations")
    }
    if (input$case_id == "matched_endpoint_zero" &&
        (!is.finite(fitted_alpha) || abs(fitted_alpha) > 1e-8 ||
         expected$grid_max > matched_score(input, 0) + 1e-8)) {
      reason <- c(reason, "zero endpoint")
    }
    if (input$case_id == "matched_endpoint_one" &&
        (!is.finite(fitted_alpha) || abs(fitted_alpha - 1) > 1e-8 ||
         expected$grid_max > matched_score(input, 1) + 1e-8)) {
      reason <- c(reason, "one endpoint")
    }
  } else {
    failed_likelihood <- numeric_value(observed$log_likelihood)
    if (!is.nan(numeric_value(observed$alpha)) ||
        !is.infinite(failed_likelihood) || failed_likelihood >= 0 ||
        !is.finite(numeric_value(observed$evaluations))) {
      reason <- c(reason, "failed-fit numerical status")
    }
  }
  if (length(reason)) record_failure("matched", input$case_id,
    paste(reason, collapse = ", "), input, observed, expected)
}

direction_forward <- matched_score(matched_cases[
  matched_cases$case_id == "direction_receiver_anchor", ], 0.2)
direction_reverse <- matched_score(matched_cases[
  matched_cases$case_id == "direction_anchor_receiver", ], 0.2)
orientation_original <- matched_score(matched_cases[
  matched_cases$case_id == "orientation_original", ], 0.2)
orientation_flipped <- matched_score(matched_cases[
  matched_cases$case_id == "orientation_flipped", ], 0.2)
if (abs(direction_forward - direction_reverse) <= 1 ||
    abs(orientation_original - orientation_flipped) > 1e-10) {
  record_failure("matched", "direction/orientation_controls",
    "independent directionality or A/B+AF symmetry failed", list(seed = seed),
    list(forward = direction_forward, reverse = direction_reverse),
    list(original = orientation_original, flipped = orientation_flipped))
}

if (nrow(failures)) {
  failure_path <- file.path(work_dir, "unexpected-failures.tsv")
  utils::write.table(failures, failure_path, sep = "\t", row.names = FALSE,
    quote = TRUE)
  stop(nrow(failures), " unexpected observations retained at ", failure_path,
    "; all campaign input/output retained at ", work_dir)
}
completed <- TRUE
cat("Somalier statistical campaign: OK; seed=", seed,
  "; binomial=", nrow(binomial_cases),
  "; eligibility=", nrow(eligibility_cases),
  "; matched=", nrow(matched_cases),
  "; failures=0; k/k+1 and all statuses/denominators checked.\n", sep = "")
