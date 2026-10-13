# Network-free checks of the declared parent/child mixture comparator.
library(tinytest)
root <- Sys.getenv("DUCKHTS_SOURCE_DIR", "../..")
source(file.path(root, "benchmarks/roh_truth_mixtures/declaration_related.R"))
source(file.path(root, "benchmarks/roh_truth_mixtures/mixtures.R"))

alpha <- c(0, 0.01, 0.02, 0.05, 0.10)
fit <- data.frame(alpha = alpha, contamination = alpha / 2,
  contamination_relative = alpha, log_likelihood = rep(-100, length(alpha)),
  log_likelihood_relative = c(-100, -101, -101, -99, -99))
erosion <- data.frame(pass = rep(TRUE, 4))
result <- roh_related_mixture_score(fit, erosion, roh_related_declaration)
expect_equal(result$criteria$criterion, c("M3", "M4", "R1"))
expect_equal(result$criteria$pass, c(TRUE, TRUE, TRUE))
expect_equal(result$fit$relative_tolerance, c(0.005, 0.005, 0.005, 0.0125, 0.025))
expect_equal(result$fit$half_share_tolerance, c(0.003, 0.003, 0.003, 0.00625, 0.0125))
expect_equal(result$fit$relative_loglik_gain, c(0, -1, -1, 1, 1))

# Low-share likelihood preference is descriptive; 5% and 10% must prefer related.
for (index in c(4L, 5L)) {
  wrong <- fit
  wrong$log_likelihood_relative[index] <- wrong$log_likelihood[index]
  expect_equal(roh_related_mixture_score(wrong, erosion, roh_related_declaration)$criteria$pass,
    c(FALSE, TRUE, TRUE))
}
for (index in 2:5) {
  wrong <- fit
  wrong$contamination_relative[index] <- alpha[index] + 0.1
  expect_equal(roh_related_mixture_score(wrong, erosion, roh_related_declaration)$criteria$pass,
    c(FALSE, TRUE, TRUE))
  wrong <- fit
  wrong$contamination[index] <- alpha[index] / 2 + 0.1
  expect_equal(roh_related_mixture_score(wrong, erosion, roh_related_declaration)$criteria$pass,
    c(TRUE, FALSE, TRUE))
}
wrong <- fit
wrong$contamination_relative[2] <- NA_real_
expect_equal(roh_related_mixture_score(wrong, erosion, roh_related_declaration)$criteria$pass,
  c(FALSE, TRUE, TRUE))
wrong <- fit
wrong$log_likelihood_relative[4] <- NA_real_
expect_equal(roh_related_mixture_score(wrong, erosion, roh_related_declaration)$criteria$pass,
  c(FALSE, TRUE, TRUE))
erosion$pass[2] <- FALSE
expect_equal(roh_related_mixture_score(fit, erosion, roh_related_declaration)$criteria$pass,
  c(TRUE, TRUE, FALSE))
erosion$pass[2] <- NA
expect_equal(roh_related_mixture_score(fit, erosion, roh_related_declaration)$criteria$pass,
  c(TRUE, TRUE, FALSE))
