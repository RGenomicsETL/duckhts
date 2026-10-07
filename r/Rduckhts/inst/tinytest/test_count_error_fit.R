library(tinytest)
library(DBI)

# An independent oracle of the count-error model, fitted in R with optim() on
# the same sites. It implements the model of the catalog entry on the sites
# themselves, not on histogram cells, so it does not share the native code's
# binning, tables or optimizer. The native optimum must be at least as good in
# log-likelihood as the oracle's, within the optimizer tolerance, and the two
# estimates must agree within the declared tolerances.
oracle_log_likelihood <- function(par, counts, relation = c("unrelated", "relative")) {
  relation <- match.arg(relation)
  e <- par[["e"]]; c <- par[["c"]]; f_excess <- par[["F"]]; b <- par[["b"]]
  rho <- c(par[["rho_hom"]], par[["rho_het"]], par[["rho_hom"]])
  w <- par[["w"]]
  f <- counts$af
  depth <- counts$ref_count + counts$alt_count
  alt <- counts$alt_count
  q <- c(e, b, 1 - e)
  hwe <- cbind((1 - f)^2, 2 * f * (1 - f), f^2)
  prior <- (1 - f_excess) * hwe + f_excess * cbind(1 - f, 0, f)
  count_probability <- function(p, rho) {
    if (rho < 1e-9) return(stats::dbinom(alt, depth, p))
    s <- (1 - rho) / rho
    exp(lchoose(depth, alt) + lbeta(alt + p * s, depth - alt + (1 - p) * s) - lbeta(p * s, (1 - p) * s))
  }
  mixture <- 0
  for (g in 1:3) {
    given_g <- if (relation == "unrelated") hwe else
      switch(g, cbind(1 - f, f, 0), cbind((1 - f) / 2, 0.5, f / 2), cbind(0, 1 - f, f))
    for (h in 1:3) {
      p <- (1 - c) * q[[g]] + c * q[[h]]
      mixture <- mixture + prior[, g] * given_g[, h] * count_probability(p, rho[[g]])
    }
  }
  sum(log(pmax((1 - w) * mixture + w / (depth + 1), 1e-300)))
}

oracle_fit <- function(counts, relation = "unrelated") {
  scale <- c(e = 0.45, c = 0.45, F = 1, b = 1, rho_hom = 0.5, rho_het = 0.5, w = 0.5)
  natural <- function(x) scale * stats::plogis(x)
  objective <- function(x) {
    value <- -oracle_log_likelihood(as.list(natural(x)), counts, relation)
    if (is.finite(value)) value else 1e300
  }
  best <- NULL
  for (c_start in c(0.005, 0.1)) {
    start <- stats::qlogis(c(e = 1e-3, c = c_start, F = 0.05, b = 0.48, rho_hom = 1e-3, rho_het = 1e-3, w = 1e-3) / scale)
    found <- stats::optim(start, objective, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-12))
    found <- stats::optim(found$par, objective, method = "BFGS", control = list(maxit = 200, reltol = 1e-12))
    if (is.null(best) || found$value < best$value) best <- found
  }
  list(estimate = natural(best$par), log_likelihood = -best$value)
}

test_count_error_fit <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  dbExecute(con, "SET threads = 2")
  # Deterministic synthetic counts, as in test/sql/count_error_fit.test but
  # with 6,000 sites on one chromosome, 2,500 bases apart so that they span
  # 15 Mb, and a share of 0.05. The fit uses 3 Mb blocks so that the sample
  # has five blocks and a block spread.
  dbExecute(con, "CREATE MACRO fit_u(salt, chrom, pos, i) AS (hash(salt, chrom, pos, i) % 1000003) / 1000003.0")
  dbExecute(con, paste(
    "CREATE TABLE counts AS",
    "WITH s AS (SELECT '1' AS chrom, i * 2500 AS pos, 0.05 + 0.9 * fit_u(1, '1', i * 2500, 0) AS af,",
    "25 + (i % 11) AS depth FROM range(1, 6001) t(i)),",
    "g AS (SELECT *, CASE WHEN fit_u(2, chrom, pos, 0) < 0.1 THEN 2 * (fit_u(3, chrom, pos, 0) < af)::INTEGER",
    "ELSE (fit_u(3, chrom, pos, 0) < af)::INTEGER + (fit_u(4, chrom, pos, 0) < af)::INTEGER END AS g,",
    "(fit_u(5, chrom, pos, 0) < af)::INTEGER + (fit_u(6, chrom, pos, 0) < af)::INTEGER AS h FROM s)",
    "SELECT 'S' AS sample_id, chrom, pos, af, depth - alt AS ref_count, alt AS alt_count FROM (",
    "SELECT chrom, pos, af, depth, (SELECT count(*) FROM range(depth) r(i) WHERE fit_u(7, chrom, pos, i) <",
    "0.95 * [1e-3, 0.5, 1 - 1e-3][g + 1] + 0.05 * [1e-3, 0.5, 1 - 1e-3][h + 1]) AS alt FROM g)"))
  counts <- dbGetQuery(con, "SELECT af, ref_count, alt_count FROM counts")

  fit <- rduckhts_count_error_fit(con, "counts", block_bases = 3e6)
  expect_equal(nrow(fit), 1L)
  expect_equal(fit$blocks, 5L)
  expect_equal(fit$sample, "S")
  expect_equal(fit$sites, 6000)
  expect_equal(fit$status, "ok")
  expect_true(abs(fit$contamination - 0.05) < 0.01)
  # The information on the read error comes from the homozygous sites whose
  # frequency is low enough that the second genome adds little, about 300
  # sites and 9,000 reads here, so its standard error is near 7e-4. The
  # precise check of the estimate is the agreement with the oracle below.
  expect_true(abs(fit$seq_error - 1e-3) < 1.5e-3)
  expect_true(abs(fit$allele_balance - 0.5) < 0.02)

  # The oracle, on the sites. The native fit works on cells with the mean
  # frequency of a 1/64 bin, so its log-likelihood is evaluated on the sites
  # here too, at its own estimate, for one comparison.
  oracle <- oracle_fit(counts)
  native_at_sites <- oracle_log_likelihood(list(
    e = fit$seq_error, c = fit$contamination, F = fit$homozygosity_excess, b = fit$allele_balance,
    rho_hom = fit$spread_hom, rho_het = fit$spread_het, w = fit$artefact_weight), counts)
  expect_true(native_at_sites >= oracle$log_likelihood - 0.5)
  expect_true(abs(fit$contamination - oracle$estimate[["c"]]) < 0.003)
  expect_true(abs(fit$seq_error - oracle$estimate[["e"]]) < 3e-4)
  expect_true(abs(fit$allele_balance - oracle$estimate[["b"]]) < 0.01)
  # The relative hypothesis, with the same oracle.
  relative <- oracle_fit(counts, "relative")
  expect_true(abs(fit$contamination_relative - relative$estimate[["c"]]) < 0.005)
  expect_true(fit$log_likelihood_relative < fit$log_likelihood)

  # The wrapper publishes a table and validates its arguments.
  expect_true(rduckhts_count_error_fit(con, "counts", block_bases = 3e6, table_name = "fit_table"))
  expect_equal(dbGetQuery(con, "SELECT count(*) AS n FROM fit_table")$n, 1)
  expect_error(rduckhts_count_error_fit(con, "counts", freq_bins = 0), "freq_bins must be one whole number")
  expect_error(rduckhts_count_error_fit(con, "counts", max_depth = 1.5), "max_depth must be one whole number")
  expect_error(rduckhts_count_error_fit(con, "counts", max_cells = 10), "max_cells")
}

test_count_error_fit()
