# BigWig reader and wrapper tests
library(tinytest)
library(DBI)

test_bigwig <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  expect_silent(rduckhts_load(con))

  bigwig_path <- system.file(
    "extdata", "libbigwig_test.bw", package = "Rduckhts", mustWork = TRUE
  )

  expect_silent(rduckhts_bigwig(
    con,
    "bigwig_intervals",
    bigwig_path,
    region = c("1:1-150", "1:100-250", "10:201-300"),
    blocks_per_iteration = 1L,
    overwrite = TRUE
  ))

  observed <- dbGetQuery(
    con,
    paste(
      "SELECT CHROM, START0, END0, round(VALUE::DOUBLE, 1) AS VALUE",
      "FROM bigwig_intervals ORDER BY CHROM, START0"
    )
  )
  expect_equal(nrow(observed), 6L)
  expect_equal(observed$CHROM, c("1", "1", "1", "1", "1", "10"))
  expect_equal(observed$START0, c(0L, 1L, 2L, 100L, 150L, 200L))
  expect_equal(observed$END0, c(1L, 2L, 3L, 150L, 151L, 300L))
  expect_equal(observed$VALUE, c(0.1, 0.2, 0.3, 1.4, 1.5, 2.0))

  expect_error(
    rduckhts_bigwig(con, "bigwig_intervals", bigwig_path),
    "already exists"
  )
  expect_error(
    rduckhts_bigwig(con, "bad_blocks", bigwig_path, blocks_per_iteration = 0),
    "between 1 and 1048576"
  )
  expect_error(
    rduckhts_bigwig(con, "bad_region", bigwig_path, region = character()),
    "region must contain"
  )
}

test_bigwig()
