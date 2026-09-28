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
