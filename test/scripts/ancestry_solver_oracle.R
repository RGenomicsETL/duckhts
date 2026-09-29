#!/usr/bin/env Rscript
# Independent QP oracle for the native projected-frequency solver.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: ancestry_solver_oracle.R extension private_library")
.libPaths(c(args[[2L]], .libPaths()))
stopifnot(as.character(packageVersion("bigsnpr")) == "1.12.21")
library(DBI)
drv <- duckdb::duckdb(shared_home = FALSE, config = list(allow_unsigned_extensions = "true"))
con <- dbConnect(drv)
on.exit(dbDisconnect(con, shutdown = TRUE))
dbExecute(con, sprintf("LOAD %s", dbQuoteString(con, normalizePath(args[[1L]]))))
sql_list <- function(v) paste0("[", paste(sprintf("%.17g", v), collapse = ","), "]")
solve_native <- function(x, y, groups, equal) {
  query <- sprintf("SELECT duckhts_ancestry_proportions(%s, %s, %d, %s) AS q",
                   sql_list(as.vector(t(x))), sql_list(y), groups,
                   if (equal) "true" else "false")
  dbGetQuery(con, query)$q[[1L]]
}

example_x <- diag(c(1e-6, 1e-6))
example_y <- c(2.5e-7, 7.5e-7)
for (scale in c(1e-12, 1, 1e6)) {
  got <- solve_native(example_x * scale, example_y * scale, 2L, TRUE)
  stopifnot(max(abs(got - c(0.25, 0.75))) < 1e-8)
  cat(sprintf("example scale=%g q=(%.9g, %.9g)\n", scale, got[[1L]], got[[2L]]))
}
set.seed(228)
for (trial in seq_len(32L)) {
  n <- 40L
  pcs <- if (trial > 30L) 16L else if (trial %% 3L == 0L) 3L else 8L
  groups <- if (trial > 30L) 21L else if (trial %% 3L == 0L) 6L else 4L
  f_ref <- matrix(runif(n * groups, 0.1, 0.9), nrow = n)
  projection <- matrix(rnorm(n * pcs), nrow = n)
  freq <- drop(f_ref %*% c(0.6, 0.4, rep(0, groups - 2L))) + rnorm(n, sd = 0.02)
  x <- crossprod(projection, f_ref)
  y <- drop(crossprod(projection, freq))
  gram <- Matrix::nearPD(crossprod(x), base.matrix = TRUE)
  stopifnot(gram$converged)
  equal <- trial %% 2L == 0L
  expected <- quadprog::solve.QP(
    gram$mat, drop(crossprod(y, x)),
    cbind(-1, diag(groups)), c(-1, rep(0, groups)),
    meq = as.integer(equal)
  )$solution
  for (scale in c(1e-12, 1, 1e6)) {
    got <- solve_native(x * scale, y * scale, groups, equal)
    error <- max(abs(got - expected))
    cat(sprintf("trial=%d groups=%d pcs=%d equality=%s scale=%g max_error=%.9g\n",
                trial, groups, pcs, equal, scale, error))
    stopifnot(error < 2e-4)
  }
}

# Group columns whose scales differ by up to 10^6 must fit no worse than the reference:
# compare objectives under the repaired quadratic, since ill-conditioned problems can
# leave the reference slightly short of the optimum.
for (trial in seq_len(200L)) {
  pcs <- c(3L, 8L, 16L)[trial %% 3L + 1L]
  groups <- c(2L, 4L, 6L, 21L)[trial %% 4L + 1L]
  f_ref <- matrix(runif(n * groups, 0.1, 0.9), nrow = n)
  projection <- matrix(rnorm(n * pcs), nrow = n)
  weights <- runif(groups); weights[sample(groups, max(0L, groups - 3L))] <- 0
  weights <- weights / sum(weights)
  equal <- trial %% 2L == 0L
  if (!equal) weights <- weights * runif(1, 0.3, 1)
  x <- crossprod(projection, f_ref) %*% diag(10^runif(groups, -3, 3), groups)
  y <- drop(x %*% weights) + rnorm(pcs, sd = 1e-3 * sqrt(mean((x %*% weights)^2)))
  gram <- Matrix::nearPD(crossprod(x), base.matrix = TRUE)
  if (!gram$converged) next
  linear <- drop(crossprod(y, x))
  expected <- quadprog::solve.QP(gram$mat, linear, cbind(-1, diag(groups)),
                                 c(-1, rep(0, groups)), meq = as.integer(equal))$solution
  got <- solve_native(x, y, groups, equal)
  objective <- function(q) drop(0.5 * crossprod(q, gram$mat %*% q) - sum(linear * q))
  excess <- (objective(got) - objective(expected)) / max(abs(objective(expected)), 1e-300)
  cat(sprintf("column-scaled trial=%d groups=%d pcs=%d equality=%s objective_excess=%.3g\n",
              trial, groups, pcs, equal, excess))
  stopifnot(excess < 1e-9)
}
