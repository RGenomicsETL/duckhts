suppressPackageStartupMessages(library(tinytest))
setwd(tinytest::get_call_wd())
source("benchmarks/count_mixture_lowpass_genome/collect.R")
setting <- genome_strata()[2L, ]
setting$share <- .1
x <- data.frame(
  status = rep("ok_conditional", 200L), share = .1,
  share_lower = .05, share_upper = .15, relation_diagnostic = "unresolved"
)
summary <- genome_summarize(x, setting)
expect_equal(summary$replicates, 200L)
expect_true(summary$evaluated)
expect_equal(summary$bias, 0, tolerance = 1e-15)
expect_true(summary$bias_pass)
expect_true(summary$accuracy_pass)
expect_equal(summary$coverage, 1)
expect_false(summary$coverage_pass)
expect_false(summary$release_pass)
expect_equal(summary$interval_mean, .1)
x$share_lower[1:10] <- .11
summary <- genome_summarize(x, setting)
expect_equal(summary$coverage, .95)
expect_true(summary$coverage_pass)
expect_true(summary$release_pass)
x$share[1L] <- NA_real_
x$share_lower[2L] <- NA_real_
x$status[3L] <- "weak_pairs"
summary <- genome_summarize(x, setting)
expect_equal(summary$missing_share, 1L)
expect_equal(summary$missing_interval, 1L)
expect_equal(summary$failed_status, 1L)
expect_equal(summary$replicates, 200L)
expect_false(summary$bias_pass)
expect_false(summary$accuracy_pass)
expect_false(summary$release_pass)
summary <- genome_summarize(x[1:199, ], setting)
expect_false(summary$evaluated)
expect_true(is.na(summary$release_pass))
expect_true(is.na(summary$accuracy_pass))
setting$share <- 0
x$share <- 0
x$share_lower <- 0
x$status <- "ok_conditional"
summary <- genome_summarize(x, setting)
expect_equal(summary$detection, 0)
expect_equal(summary$null_upper, 1 - .05^(1 / 200), tolerance = 1e-14)
expect_true(summary$null_pass)
x$share_lower[1:20] <- .01
summary <- genome_summarize(x, setting)
expect_equal(summary$detection, .1)
expect_false(summary$null_pass)
x$share_lower[1L] <- .2
summary <- genome_summarize(x, setting)
expect_equal(summary$missing_interval, 1L)
x$status[2L] <- NA_character_
x$relation_diagnostic[2L] <- NA_character_
summary <- genome_summarize(x, setting)
expect_equal(summary$failed_status, 1L)
expect_false(summary$diagnostic_pass)
