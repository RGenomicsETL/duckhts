library(tinytest)
library(DBI)

test_bcf_filter_list_fetch_regression <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit({
    dbDisconnect(con, shutdown = TRUE)
  }, add = TRUE)

  expect_silent(rduckhts_load(con))

  bcf_path <- system.file("extdata", "bcf_filter_list_regression.vcf", package = "Rduckhts")
  expect_true(file.exists(bcf_path))

  rs <- dbSendQuery(
    con,
    sprintf(
      "SELECT POS, FILTER::VARCHAR AS filter_str FROM read_bcf('%s') ORDER BY POS",
      bcf_path
    )
  )
  on.exit({
    if (DBI::dbIsValid(rs)) {
      DBI::dbClearResult(rs)
    }
  }, add = TRUE)

  n_rows <- 0L
  saw_multi_filter <- FALSE
  repeat {
    chunk <- dbFetch(rs, n = 512)
    if (nrow(chunk) == 0L) {
      break
    }
    n_rows <- n_rows + nrow(chunk)
    if (any(chunk$filter_str == "[q10, q20]")) {
      saw_multi_filter <- TRUE
    }
  }

  expect_equal(n_rows, 5000L)
  expect_true(saw_multi_filter)

  summary <- dbGetQuery(
    con,
    sprintf(
      "SELECT COUNT(*) AS n_rows, SUM(length(FILTER)) AS filter_items FROM read_bcf('%s')",
      bcf_path
    )
  )
  expect_equal(summary$n_rows[1], 5000L)
  expect_equal(summary$filter_items[1], 5384L)
}

test_bcf_filter_list_fetch_regression()
