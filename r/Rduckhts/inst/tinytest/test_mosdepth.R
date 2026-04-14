library(tinytest)
library(DBI)

test_mosdepth <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  expect_silent(rduckhts_load(con))

  bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
  bam_index_path <- system.file("extdata", "range.bam.bai", package = "Rduckhts")

  tmp_dir <- tempfile("duckhts_mosdepth_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  prefix_fast <- file.path(tmp_dir, "range_fast")
  out_fast <- rduckhts_mosdepth(
    con,
    prefix_fast,
    bam_path,
    index_path = bam_index_path,
    overwrite = TRUE
  )
  expect_true(isTRUE(out_fast$success[1]))
  expect_true(file.exists(file.path(tmp_dir, "range_fast.mosdepth.summary.txt")))
  expect_true(file.exists(file.path(tmp_dir, "range_fast.mosdepth.global.dist.txt")))
  expect_true(file.exists(file.path(tmp_dir, "range_fast.per-base.bed.gz")))
  expect_true(file.exists(file.path(tmp_dir, "range_fast.per-base.bed.gz.csi")))

  summary_lines <- readLines(file.path(tmp_dir, "range_fast.mosdepth.summary.txt"))
  expect_equal(
    summary_lines[1:3],
    c(
      "chrom\tlength\tbases\tmean\tmin\tmax",
      "CHROMOSOME_I\t1009800\t1730\t0.00\t0\t4",
      "CHROMOSOME_II\t5000\t3400\t0.68\t0\t5"
    )
  )

  per_base <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT chrom, start_pos, end_pos, depth",
        "FROM read_csv(%s,",
        "  columns = {'chrom':'VARCHAR', 'start_pos':'BIGINT', 'end_pos':'BIGINT', 'depth':'INTEGER'},",
        "  delim = '\t',",
        "  header = FALSE",
        ")",
        "WHERE chrom = 'CHROMOSOME_II'",
        "ORDER BY start_pos",
        "LIMIT 3"
      ),
      sprintf("'%s'", file.path(tmp_dir, "range_fast.per-base.bed.gz"))
    )
  )
  expect_equal(
    per_base,
    data.frame(
      chrom = c("CHROMOSOME_II", "CHROMOSOME_II", "CHROMOSOME_II"),
      start_pos = c(0, 1135, 1235),
      end_pos = c(1135, 1235, 1240),
      depth = c(0, 1, 0)
    )
  )

  prefix_win <- file.path(tmp_dir, "range_win")
  out_win <- rduckhts_mosdepth(
    con,
    prefix_win,
    bam_path,
    index_path = bam_index_path,
    by = "1000",
    overwrite = TRUE
  )
  expect_true(isTRUE(out_win$success[1]))
  expect_true(file.exists(file.path(tmp_dir, "range_win.regions.bed.gz")))
  expect_true(file.exists(file.path(tmp_dir, "range_win.regions.bed.gz.csi")))
  expect_true(file.exists(file.path(tmp_dir, "range_win.mosdepth.region.dist.txt")))

  region_dist <- readLines(file.path(tmp_dir, "range_win.mosdepth.region.dist.txt"))
  expect_equal(
    region_dist[1:5],
    c(
      "CHROMOSOME_I\t1\t0.00",
      "CHROMOSOME_I\t0\t1.00",
      "CHROMOSOME_II\t2\t0.20",
      "CHROMOSOME_II\t1\t0.40",
      "CHROMOSOME_II\t0\t1.00"
    )
  )

  prefix_noper <- file.path(tmp_dir, "range_noper")
  out_noper <- rduckhts_mosdepth(
    con,
    prefix_noper,
    bam_path,
    index_path = bam_index_path,
    by = "1000",
    no_per_base = TRUE,
    overwrite = TRUE
  )
  expect_true(is.na(out_noper$per_base_path[1]))
  expect_true(file.exists(file.path(tmp_dir, "range_noper.regions.bed.gz")))
  expect_false(file.exists(file.path(tmp_dir, "range_noper.per-base.bed.gz")))

  expect_error(
    rduckhts_mosdepth(
      con,
      file.path(tmp_dir, "range_default"),
      bam_path,
      index_path = bam_index_path,
      fast_mode = FALSE,
      overwrite = TRUE
    )
  )
}

test_mosdepth()
