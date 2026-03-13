library(tinytest)
library(DBI)

test_cigar_utils <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit({
    try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
  }, add = TRUE)

  expect_silent(rduckhts_load(con))

  literal_metrics <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "cigar_left_soft_clip('5S90M5S') AS left_soft,",
      "cigar_right_soft_clip('5S90M5S') AS right_soft,",
      "cigar_query_length('5S90M5I') AS query_len,",
      "cigar_aligned_query_length('5S90M5I') AS aligned_query_len,",
      "cigar_reference_length('90M5D') AS ref_len,",
      "cigar_has_soft_clip('5S90M5S') AS has_soft,",
      "cigar_has_hard_clip('5H95M') AS has_hard,",
      "cigar_has_op('90M5D', 'D') AS has_d"
    )
  )

  expect_equal(literal_metrics$left_soft[[1]], 5)
  expect_equal(literal_metrics$right_soft[[1]], 5)
  expect_equal(literal_metrics$query_len[[1]], 100)
  expect_equal(literal_metrics$aligned_query_len[[1]], 90)
  expect_equal(literal_metrics$ref_len[[1]], 95)
  expect_true(isTRUE(literal_metrics$has_soft[[1]]))
  expect_true(isTRUE(literal_metrics$has_hard[[1]]))
  expect_true(isTRUE(literal_metrics$has_d[[1]]))

  flag_bits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "(sam_flag_bits(99)).is_paired AS is_paired,",
      "(sam_flag_bits(99)).is_proper_pair AS is_proper_pair,",
      "sam_flag_has(99, 2) AS has_proper_pair"
    )
  )

  expect_true(isTRUE(flag_bits$is_paired[[1]]))
  expect_true(isTRUE(flag_bits$is_proper_pair[[1]]))
  expect_true(isTRUE(flag_bits$has_proper_pair[[1]]))

  forward_bits <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT",
      "is_forward_aligned(0) AS forward_zero,",
      "is_forward_aligned(16) AS reverse_sixteen,",
      "is_forward_aligned(4) AS unmapped_four"
    )
  )

  expect_true(isTRUE(forward_bits$forward_zero[[1]]))
  expect_false(isTRUE(forward_bits$reverse_sixteen[[1]]))
  expect_true(is.na(forward_bits$unmapped_four[[1]]))
}

test_cigar_utils()
