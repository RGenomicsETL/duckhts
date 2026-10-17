library(tinytest)
library(DBI)

# Enumerate genotype pairs and the transmitted allele. Neither the native
# coefficient formulas nor its compressed block matrices supply expectations.
lowpass_expectation <- function(f, g, m, e = 0, w = 0, F = 0, ibd1 = FALSE) {
  prior_x <- cbind((1-F)*(1-f)^2 + F*(1-f), (1-F)*2*f*(1-f), (1-F)*f^2 + F*f)
  prior_y <- cbind((1-g)^2, 2*g*(1-g), g^2)
  mean_read <- numeric(length(f))
  discordance <- numeric(length(f))
  for (x in 0:2) {
    if (ibd1) {
      transmission <- x/2
      conditional_y <- cbind((1-transmission)*(1-f),
                            transmission*(1-f)+(1-transmission)*f, transmission*f)
    } else {
      conditional_y <- prior_y
    }
    for (y in 0:2) {
      probability <- prior_x[, x+1L]*conditional_y[, y+1L]
      q <- e + (1-2*e)*((1-m)*x/2 + m*y/2)
      mean_read <- mean_read + probability*q
      discordance <- discordance + probability*2*q*(1-q)
    }
  }
  list(read = (1-w)*mean_read+w/2, pair = (1-w)*discordance+w/3)
}

# The independent objective is a quartic, reconstructed from enumerated site
# expectations at three shares. Evaluate every real stationary point in [0,.5].
lowpass_oracle <- function(counts, e, w, F, ibd1) {
  f <- counts$sample_af
  g <- counts$candidate_af
  n <- counts$ref_count+counts$alt_count
  p <- choose(n, 2)
  u <- cbind(1, f, g)
  x <- cbind(1, f, g, f*f, f*g, g*g)
  block <- (counts$pos-1L) %/% 1024L
  shares <- c(0, .25, .5)
  residuals <- lapply(shares, function(m) {
    expected <- lowpass_expectation(f, g, m, e, w, F, ibd1)
    read <- u*(counts$alt_count-n*expected$read)/sum(n)
    pairs <- x*(counts$ref_count*counts$alt_count-p*expected$pair)/max(1, sum(p))
    rowsum(cbind(read, pairs), block, reorder = FALSE)
  })
  vandermonde <- outer(shares, 0:2, `^`)
  inverse <- solve(vandermonde)
  coefficient <- lapply(seq_len(nrow(residuals[[1L]])), function(b) {
    values <- do.call(cbind, lapply(residuals, function(r) r[b, ]))
    values %*% t(inverse)
  })
  total <- Reduce(`+`, coefficient)
  evaluate <- function(coef, m) drop(coef %*% c(1, m, m*m))
  minimize <- function(weight) {
    polynomial <- numeric(5)
    for (i in 0:2) {
      for (j in 0:2) {
        polynomial[i+j+1L] <- polynomial[i+j+1L] +
          drop(crossprod(total[, i+1L], weight %*% total[, j+1L]))
      }
    }
    roots <- polyroot((1:4)*polynomial[2:5])
    candidates <- sort(c(0, .5, Re(roots[abs(Im(roots)) < 1e-7 & Re(roots) > 0 & Re(roots) < .5])))
    objectives <- vapply(candidates, function(m) {
      r <- evaluate(total, m)
      drop(crossprod(r, weight %*% r))
    }, numeric(1))
    candidates[which.min(objectives)]
  }
  first <- minimize(diag(9))
  block_scores <- t(vapply(coefficient, evaluate, numeric(9), m = first))
  B <- nrow(block_scores)
  centered <- sweep(block_scores, 2, colMeans(block_scores))
  covariance <- crossprod(centered)*B/(B-1)
  eig <- eigen(covariance, symmetric = TRUE)
  keep <- eig$values > max(1e-28, max(eig$values)*1e-10)
  if (!any(keep)) return(first)
  vectors <- eig$vectors[, keep, drop = FALSE]
  weight <- tcrossprod(sweep(vectors, 2, sqrt(eig$values[keep]), `/`))
  minimize(weight)
}

lowpass_simulation <- function(m, depth, model, sites = 32L*1024L) {
  f <- runif(sites, .05, .95)
  g <- if (model == "other_population") .05+.9*(1-f)^2 else f
  gx <- rbinom(sites, 2, f)
  gy <- if (model == "ibd1") rbinom(sites, 1, gx/2)+rbinom(sites, 1, f) else rbinom(sites, 2, g)
  n <- rpois(sites, depth)
  q <- .01+.98*((1-m)*gx/2+m*gy/2)
  alt <- rbinom(sites, n, q)
  artefact <- runif(sites) < .02
  alt[artefact] <- floor(runif(sum(artefact))*(n[artefact]+1))
  data.frame(sample_id = "s", candidate = "c", chrom = "20", pos = seq_len(sites),
             ref_count = n-alt, alt_count = alt, other_count = 0L,
             sample_af = f, candidate_af = g)
}

run_lowpass_tests <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  dbExecute(con, "SET threads=2")
  set.seed(332)
  for (model in c("unrelated", "other_population", "ibd1")) {
    for (depth in c(.1, .4)) {
      for (share in c(0, .05, .2)) {
        counts <- lowpass_simulation(share, depth, model)
        dbWriteTable(con, "counts", counts, overwrite = TRUE)
        fitted <- rduckhts_count_mixture_lowpass(con, "counts", .01, .02, 0, block_bases = 1024)
        relation <- if (model == "ibd1") "ibd1" else "unrelated"
        observed <- fitted[fitted$relationship_model == relation, ]
        expected <- lowpass_oracle(counts, .01, .02, 0, model == "ibd1")
        expect_equal(observed$share, expected, tolerance = 1e-5,
                     info = paste(model, depth, share))
        expect_true(observed$blocks >= 30 && observed$blocks <= 32)
        expect_true(observed$share >= 0 && observed$share <= .5)
        expect_true(is.finite(observed$share_lower) && is.finite(observed$share_upper))
        expect_true(observed$share_lower <= observed$share && observed$share_upper >= observed$share)
        if (observed$share == 0) {
          expect_equal(observed$share_lower, 0)
          expect_true(observed$share_upper > 0)
        }
        if (model != "other_population" && observed$share > 1e-8) {
          expect_equal(observed$relation_diagnostic, "unresolved")
        }
      }
    }
  }
  below_floor <- rduckhts_count_mixture_lowpass(con, "counts", .01, .02, 0,
                                               block_bases = 1024, min_blocks = 33)
  expect_true(all(below_floor$status == "insufficient_blocks"))
  expect_true(all(is.finite(below_floor$share)))
  expect_true(all(is.na(below_floor$share_lower) & is.na(below_floor$share_upper)))
  # The uniform count artefact has discordance 1/3 at every usable depth.
  for (n in c(2, 3, 8)) {
    expect_equal(mean((0:n)*(n-(0:n)))/choose(n, 2), 1/3, tolerance = 1e-15)
  }
  # Exact moment-equivalent count frequencies, obtained by genotype enumeration.
  f <- c(.2, .4, .6, .8)
  ibd <- lowpass_expectation(f, f, .2, ibd1 = TRUE)
  equivalent <- (1-sqrt(.68))/2
  unrelated <- lowpass_expectation(f, f, equivalent)
  expect_true(abs(equivalent-.0876894) <= 5e-8)
  expect_equal(ibd, unrelated, tolerance = 1e-14)
  probabilities <- cbind(1-ibd$read-ibd$pair/2, ibd$pair, ibd$read-ibd$pair/2)
  multiplicity <- round(probabilities*625)
  expect_equal(multiplicity/625, probabilities, tolerance = 1e-14)
  exact <- do.call(rbind, lapply(seq_along(f), function(i) {
    alt <- rep(0:2, multiplicity[i, ])
    data.frame(sample_id = "exact", candidate = "matched", chrom = "20",
               ref_count = 2L-alt, alt_count = alt, other_count = 0L,
               sample_af = f[i], candidate_af = f[i])
  }))
  exact <- exact[rep(seq_len(nrow(exact)), 4L), ]
  exact$pos <- seq_len(nrow(exact))
  dbWriteTable(con, "counts", exact, overwrite = TRUE)
  fitted <- rduckhts_count_mixture_lowpass(con, "counts", 0, 0, 0, block_bases = 2500, min_blocks = 2)
  expect_equal(fitted$share[fitted$relationship_model == "ibd1"], .2, tolerance = 1e-5)
  expect_equal(fitted$share[fitted$relationship_model == "unrelated"], equivalent, tolerance = 1e-5)
  expect_true(all(fitted$relation_diagnostic == "unresolved"))
  expect_true(rduckhts_count_mixture_lowpass(con, "counts", 0, 0, 0, table_name = "result"))
  expect_equal(dbGetQuery(con, "SELECT count(*) n FROM result")$n, 2)
  expect_error(rduckhts_count_mixture_lowpass(con, "counts", 0, 0, 0, table_name = "result"))
  expect_true(rduckhts_count_mixture_lowpass(con, "counts", 0, 0, 0, table_name = "result", overwrite = TRUE))
  expect_error(rduckhts_count_mixture_lowpass(con, "counts"), "seq_error")
  expect_error(rduckhts_count_mixture_lowpass(con, "counts", NaN, 0, 0), "seq_error")
  expect_error(rduckhts_count_mixture_lowpass(con, "counts", 0, 1, 0), "artefact_weight")
  expect_error(rduckhts_count_mixture_lowpass(con, "counts", 0, 0, -1), "homozygosity_excess")
  expect_error(rduckhts_count_mixture_lowpass(con, "counts", 0, 0, 0, min_blocks = 1), "min_blocks")
  expect_error(rduckhts_count_mixture_lowpass(con, "counts", 0, 0, 0, max_candidates = 65), "max_candidates")
  expect_error(rduckhts_count_mixture_lowpass(con, "counts", 0, 0, 0, max_state_bytes = .5), "max_state_bytes")
  expect_error(rduckhts_count_mixture_lowpass(con, "counts", 0, 0, 0, max_depth = 8.5), "max_depth")

  # Parallel reduction preserves fits and statuses at vector-sized input.
  set.seed(332)
  parallel_counts <- data.frame(
    sample_id = "parallel", candidate = "cross", chrom = "20", pos = seq_len(20000),
    ref_count = sample.int(3L, 20000, replace = TRUE) - 1L,
    alt_count = sample.int(3L, 20000, replace = TRUE) - 1L,
    other_count = 0L, sample_af = runif(20000, .01, .99),
    candidate_af = runif(20000, .01, .99)
  )
  dbWriteTable(con, "mixture_parallel", parallel_counts, overwrite = TRUE)
  dbExecute(con, "SET threads=1")
  one <- rduckhts_count_mixture_lowpass(con, "mixture_parallel", .0001, .002, 0,
                                      block_bases = 500, min_blocks = 30)
  dbExecute(con, "SET threads=4")
  four <- rduckhts_count_mixture_lowpass(con, "mixture_parallel", .0001, .002, 0,
                                       block_bases = 500, min_blocks = 30)
  expect_identical(one$status, four$status)
  expect_identical(is.na(one$share), is.na(four$share))
  expect_equal(one$share, four$share, tolerance = 1e-10)
  expect_equal(one$pair_denominator, four$pair_denominator, tolerance = 1e-10)
}
run_lowpass_tests()
