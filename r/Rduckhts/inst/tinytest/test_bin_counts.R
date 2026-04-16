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
  unmapped_sam <- system.file("extdata", "fixture_unmapped.sam", package = "Rduckhts")

  bins_none <- rduckhts_bam_bin_counts(con, mixed_bam, 5000, rmdup = "none")
  bins_streaming <- rduckhts_bam_bin_counts(con, mixed_bam, 5000, rmdup = "streaming")
  bins_flag <- rduckhts_bam_bin_counts(con, mixed_bam, 5000, rmdup = "flag")

  expect_equal(bins_none$bin_id, 0:9)
  expect_equal(bins_streaming$bin_id, 0:9)
  expect_equal(bins_flag$bin_id, 0:9)
  expect_equal(bins_none$count_total, c(4, 2, 2, 2, 0, 0, 0, 0, 0, 0))
  expect_equal(bins_streaming$count_total, c(2, 2, 1, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(bins_flag$count_total, c(4, 0, 2, 2, 0, 0, 0, 0, 0, 0))
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
  expect_equal(mq_only$bin_id, 0:9)
  expect_equal(mq_only$count_total, rep(0, 10))
  expect_equal(mq_only$count_pre, c(4, 2, 2, 2, 0, 0, 0, 0, 0, 0))
  expect_equal(mq_only$mapq_sum_pre, c(240, 120, 120, 120, 0, 0, 0, 0, 0, 0))
  expect_equal(mq_only$mean_mapq_pre[1:4], rep(60, 4))
  expect_true(all(is.na(mq_only$mean_mapq_pre[5:10])))
  expect_equal(mq_only$mapq_sum_post, rep(0, 10))
  expect_true(all(is.na(mq_only$mean_mapq_post)))

  stats_cram <- rduckhts_bam_bin_counts(
    con,
    mixed_cram,
    5000,
    reference = fixture_ref,
    rmdup = "streaming",
    stats = "gc,mq"
  )
  expect_equal(stats_cram$bin_id, 0:9)
  expect_equal(stats_cram$count_total, c(2, 2, 1, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(stats_cram$count_pre, c(4, 2, 2, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(stats_cram$gc_bases_pre, c(72, 0, 72, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(stats_cram$bases_pre, c(144, 72, 72, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(stats_cram$gc_perc_pre[1:3], c(0.5, 0, 1))
  expect_true(all(is.na(stats_cram$gc_perc_pre[4:10])))
  expect_equal(stats_cram$gc_bases_post, c(0, 0, 36, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(stats_cram$bases_post, c(72, 72, 36, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(stats_cram$gc_perc_post[1:3], c(0, 0, 1))
  expect_true(all(is.na(stats_cram$gc_perc_post[4:10])))
  expect_equal(stats_cram$ref_gc_bases, rep(0, 10))
  expect_equal(stats_cram$ref_bases, rep(5000, 10))
  expect_equal(stats_cram$gc_perc_ref, rep(0, 10))
  expect_equal(stats_cram$mean_mapq_pre[1:3], rep(60, 3))
  expect_true(all(is.na(stats_cram$mean_mapq_pre[4:10])))
  expect_equal(stats_cram$mean_mapq_post[1:3], rep(60, 3))
  expect_true(all(is.na(stats_cram$mean_mapq_post[4:10])))

  samtools <- Sys.which("samtools")
  if (nzchar(samtools)) {
    unmapped_bam <- tempfile("fixture_unmapped_", fileext = ".bam")
    expect_equal(
      system2(samtools, c("view", "-b", "-o", unmapped_bam, unmapped_sam)),
      0L
    )
    expect_equal(system2(samtools, c("index", unmapped_bam)), 0L)

    unmapped_bins <- rduckhts_bam_bin_counts(
      con,
      unmapped_bam,
      5000,
      include_unmapped = TRUE,
      rmdup = "flag",
      stats = "gc,mq"
    )
    unmapped_row <- unmapped_bins[unmapped_bins$chrom == "*", , drop = FALSE]

    expect_equal(nrow(unmapped_row), 1)
    expect_true(is.na(unmapped_row$start[[1]]))
    expect_true(is.na(unmapped_row$end[[1]]))
    expect_true(is.na(unmapped_row$bin_id[[1]]))
    expect_equal(unmapped_row$count_total[[1]], 2)
    expect_equal(unmapped_row$count_fwd[[1]], 1)
    expect_equal(unmapped_row$count_rev[[1]], 1)
    expect_equal(unmapped_row$count_pre[[1]], 3)
    expect_equal(unmapped_row$gc_bases_pre[[1]], 72)
    expect_equal(unmapped_row$bases_pre[[1]], 108)
    expect_equal(unmapped_row$gc_perc_pre[[1]], 2 / 3)
    expect_equal(unmapped_row$gc_bases_post[[1]], 72)
    expect_equal(unmapped_row$bases_post[[1]], 72)
    expect_equal(unmapped_row$gc_perc_post[[1]], 1)
    expect_equal(unmapped_row$mapq_sum_pre[[1]], 180)
    expect_equal(unmapped_row$mean_mapq_pre[[1]], 60)
    expect_equal(unmapped_row$mapq_sum_post[[1]], 120)
    expect_equal(unmapped_row$mean_mapq_post[[1]], 60)
  }

  expect_error(
    rduckhts_bam_bin_counts(con, mixed_bam, 5000, rmdup = "weird")
  )
  expect_error(
    rduckhts_bam_bin_counts(con, mixed_bam, 5000, stats = "gc,wat")
  )
}

test_bin_counts()
