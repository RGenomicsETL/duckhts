library(tinytest)
library(DBI)

test_samtools_idxstats <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  mixed_bam <- system.file("extdata", "fixture_mixed.bam", package = "Rduckhts")
  mixed_cram <- system.file("extdata", "fixture_mixed.cram", package = "Rduckhts")
  range_bam <- system.file("extdata", "range.bam", package = "Rduckhts")
  range_bai <- system.file("extdata", "range.bam.bai", package = "Rduckhts")

  bam_out <- tempfile("idxstats_mixed_", fileext = ".txt")
  bam_res <- rduckhts_samtools_idxstats(con, mixed_bam, output = bam_out, overwrite = TRUE)
  expect_true(isTRUE(bam_res$success[[1]]))
  expect_true(isTRUE(bam_res$used_index_fast_path[[1]]))
  bam_lines <- readLines(bam_out, warn = FALSE)
  expect_equal(bam_lines, c("11\t50000\t10\t0", "*\t0\t0\t0"))

  cram_out <- tempfile("idxstats_cram_", fileext = ".txt")
  cram_res <- rduckhts_samtools_idxstats(con, mixed_cram, output = cram_out, overwrite = TRUE)
  expect_true(isTRUE(cram_res$success[[1]]))
  expect_true(!isTRUE(cram_res$used_index_fast_path[[1]]))
  cram_lines <- readLines(cram_out, warn = FALSE)
  expect_equal(cram_lines, c("11\t50000\t10\t0", "*\t0\t0\t0"))

  custom_out <- tempfile("idxstats_range_", fileext = ".txt")
  custom_res <- rduckhts_samtools_idxstats(
    con,
    range_bam,
    output = custom_out,
    index_path = range_bai,
    overwrite = TRUE
  )
  expect_true(isTRUE(custom_res$success[[1]]))
  expect_true(isTRUE(custom_res$used_index_fast_path[[1]]))
  custom_lines <- readLines(custom_out, warn = FALSE)
  expect_equal(
    custom_lines,
    c(
      "CHROMOSOME_I\t1009800\t18\t0",
      "CHROMOSOME_II\t5000\t34\t0",
      "CHROMOSOME_III\t5000\t41\t0",
      "CHROMOSOME_IV\t5000\t19\t0",
      "CHROMOSOME_V\t5000\t0\t0",
      "CHROMOSOME_X\t5000\t0\t0",
      "CHROMOSOME_MtDNA\t5000\t0\t0",
      "*\t0\t0\t0"
    )
  )

  expect_error(
    rduckhts_samtools_idxstats(con, mixed_bam, output = bam_out),
    pattern = "output already exists"
  )
}

test_samtools_idxstats()
