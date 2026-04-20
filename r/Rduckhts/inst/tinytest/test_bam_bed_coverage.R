library(tinytest)
library(DBI)

test_bam_bed_coverage <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  expect_silent(rduckhts_load(con))

  mixed_bam <- system.file("extdata", "fixture_mixed.bam", package = "Rduckhts")
  mixed_bed <- system.file("extdata", "fixture_mixed_regions.bed", package = "Rduckhts")
  min_depth_bed <- system.file("extdata", "fixture_mixed_min_depth_regions.bed", package = "Rduckhts")

  cov <- rduckhts_bam_bed_coverage(con, mixed_bam, mixed_bed)
  expect_equal(cov$chrom, c("11", "11", "11"))
  expect_equal(cov$start, c(1000, 6000, 20000))
  expect_equal(cov$end, c(1200, 6200, 20100))
  expect_equal(cov$name, c("cluster", "dup_only", "zero"))

  expect_equal(cov$numreads_pre, c(4, 2, 0))
  expect_equal(cov$covbases_pre, c(72, 72, 0))
  expect_equal(round(cov$coverage_pre, 2), c(36, 36, 0))
  expect_equal(round(cov$meandepth_pre, 2), c(0.72, 0.36, 0.00))
  expect_equal(round(cov$meanmapq_pre, 1), c(60, 60, 0))
  expect_equal(round(cov$meanbaseq_pre, 1), c(37, 37, 0))

  expect_equal(cov$numreads_post, c(4, 0, 0))
  expect_equal(cov$covbases_post, c(72, 0, 0))
  expect_equal(round(cov$coverage_post, 2), c(36, 0, 0))
  expect_equal(round(cov$meandepth_post, 2), c(0.72, 0.00, 0.00))
  expect_equal(round(cov$meanmapq_post, 1), c(60, 0, 0))
  expect_equal(round(cov$meanbaseq_post, 1), c(37, 0, 0))

  expect_equal(cov$numreads_fwd_post, c(2, 0, 0))
  expect_equal(cov$numreads_rev_post, c(2, 0, 0))
  expect_equal(cov$covbases_fwd_post, c(36, 0, 0))
  expect_equal(cov$covbases_rev_post, c(36, 0, 0))
  expect_equal(round(cov$coverage_fwd_post, 2), c(18, 0, 0))
  expect_equal(round(cov$coverage_rev_post, 2), c(18, 0, 0))
  expect_equal(round(cov$meandepth_fwd_post, 2), c(0.36, 0, 0))
  expect_equal(round(cov$meandepth_rev_post, 2), c(0.36, 0, 0))

  cov_min_depth <- rduckhts_bam_bed_coverage(
    con,
    mixed_bam,
    min_depth_bed,
    require_flags = 64,
    min_depth = 2,
    strand_outputs = FALSE,
    decompression_threads = 1
  )
  expect_equal(cov_min_depth$numreads_post, 1)
  expect_equal(cov_min_depth$covbases_post, 0)
  expect_equal(round(cov_min_depth$coverage_post, 2), 0)
  expect_equal(round(cov_min_depth$meandepth_post, 2), 0)
  expect_equal(cov_min_depth$covbases_fwd_post, 0)
  expect_equal(cov_min_depth$covbases_rev_post, 0)

  cov_high_min_depth <- rduckhts_bam_bed_coverage(
    con,
    mixed_bam,
    mixed_bed,
    min_depth = 100,
    strand_outputs = FALSE
  )
  expect_equal(cov_high_min_depth$name, c("cluster", "dup_only", "zero"))
  expect_equal(cov_high_min_depth$numreads_post, c(4, 0, 0))
  expect_equal(cov_high_min_depth$covbases_post, c(0, 0, 0))
  expect_equal(round(cov_high_min_depth$coverage_post, 2), c(0, 0, 0))
  expect_equal(round(cov_high_min_depth$meandepth_post, 2), c(0, 0, 0))

  expect_error(
    rduckhts_bam_bed_coverage(con, mixed_bam, mixed_bed, decompression_threads = -1),
    pattern = "decompression_threads must be a single whole number"
  )

  expect_error(
    rduckhts_bam_bed_coverage(con, mixed_bam, mixed_bed, fragment_mode = TRUE),
    pattern = "fragment_mode is not implemented yet"
  )
}

test_bam_bed_coverage()
