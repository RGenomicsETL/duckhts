library(tinytest)
library(DBI)

test_bin_counts <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  expect_silent(rduckhts_load(con))

  mixed_bam <- system.file("extdata", "fixture_mixed.bam", package = "Rduckhts")
  mixed_cram <- system.file("extdata", "fixture_mixed.cram", package = "Rduckhts")
  fixture_ref <- system.file("extdata", "fixture_ref.fa", package = "Rduckhts")

  bins_none <- rduckhts_bam_bin_counts(con, mixed_bam, 5000, rmdup = "none")
  bins_streaming <- rduckhts_bam_bin_counts(con, mixed_bam, 5000, rmdup = "streaming")
  bins_flag <- rduckhts_bam_bin_counts(con, mixed_bam, 5000, rmdup = "flag")

  expect_equal(bins_none$count_total, c(4, 2, 2, 2))
  expect_equal(bins_streaming$count_total, c(2, 2, 1))
  expect_equal(bins_flag$count_total, c(4, 2, 2))
  expect_equal(sum(bins_none$count_total), 10)
  expect_equal(sum(bins_streaming$count_total), 5)
  expect_equal(sum(bins_flag$count_total), 8)

  excl_dup <- rduckhts_bam_bin_counts(con, mixed_bam, 5000, exclude_flags = 1024)
  expect_equal(sum(excl_dup$count_total), 8)

  read1_only <- rduckhts_bam_bin_counts(con, mixed_bam, 5000, require_flags = 64)
  expect_equal(sum(read1_only$count_total), 4)

  mq_only <- rduckhts_bam_bin_counts(
    con,
    mixed_bam,
    5000,
    mapq = 61,
    rmdup = "none",
    stats = "mq"
  )
  expect_equal(mq_only$bin_id, c(0, 1, 2, 3))
  expect_equal(mq_only$count_total, c(0, 0, 0, 0))
  expect_equal(mq_only$count_pre, c(4, 2, 2, 2))
  expect_equal(mq_only$mapq_sum_pre, c(240, 120, 120, 120))
  expect_true(all(mq_only$mean_mapq_pre == 60))
  expect_equal(mq_only$mapq_sum_post, c(0, 0, 0, 0))
  expect_true(all(is.na(mq_only$mean_mapq_post)))

  stats_cram <- rduckhts_bam_bin_counts(
    con,
    mixed_cram,
    5000,
    reference = fixture_ref,
    rmdup = "streaming",
    stats = "gc,mq"
  )
  expect_equal(stats_cram$bin_id, c(0, 1, 2))
  expect_equal(stats_cram$count_total, c(2, 2, 1))
  expect_equal(stats_cram$count_pre, c(4, 2, 2))
  expect_equal(stats_cram$gc_bases_pre, c(72, 0, 72))
  expect_equal(stats_cram$bases_pre, c(144, 72, 72))
  expect_equal(stats_cram$gc_perc_pre, c(0.5, 0, 1))
  expect_equal(stats_cram$gc_bases_post, c(0, 0, 36))
  expect_equal(stats_cram$bases_post, c(72, 72, 36))
  expect_equal(stats_cram$gc_perc_post, c(0, 0, 1))
  expect_equal(stats_cram$ref_gc_bases, c(0, 0, 0))
  expect_equal(stats_cram$ref_bases, c(5000, 5000, 5000))
  expect_equal(stats_cram$gc_perc_ref, c(0, 0, 0))
  expect_true(all(stats_cram$mean_mapq_pre == 60))
  expect_true(all(stats_cram$mean_mapq_post == 60))

  expect_error(
    rduckhts_bam_bin_counts(con, mixed_bam, 5000, rmdup = "weird")
  )
  expect_error(
    rduckhts_bam_bin_counts(con, mixed_bam, 5000, stats = "gc,wat")
  )
}

test_bin_counts()
