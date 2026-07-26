library(tinytest)
library(DBI)

test_mosdepth <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
  bam_index_path <- system.file("extdata", "range.bam.bai", package = "Rduckhts")
  cram_path <- system.file("extdata", "range.cram", package = "Rduckhts")
  cram_index_path <- system.file("extdata", "range.cram.crai", package = "Rduckhts")
  fasta_path <- system.file("extdata", "ce.fa", package = "Rduckhts")

  tmp_dir <- tempfile("duckhts_mosdepth_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  prefix_fast <- file.path(tmp_dir, "range_fast")
  out_fast <- rduckhts_mosdepth(
    con,
    prefix_fast,
    bam_path,
    index_path = bam_index_path,
    fast_mode = TRUE,
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

  prefix_prec <- file.path(tmp_dir, "range_prec")
  out_prec <- rduckhts_mosdepth(
    con,
    prefix_prec,
    bam_path,
    index_path = bam_index_path,
    fast_mode = TRUE,
    precision_digits = 4,
    overwrite = TRUE
  )
  expect_true(isTRUE(out_prec$success[1]))
  summary_prec <- readLines(file.path(tmp_dir, "range_prec.mosdepth.summary.txt"))
  expect_true(any(grepl("^CHROMOSOME_I\\t1009800\\t1730\\t0\\.0017\\t0\\t4$", summary_prec)))
  expect_true(any(grepl("^total\\t1024800\\t11129\\t0\\.0109\\t0\\t6$", summary_prec)))

  prefix_cram <- file.path(tmp_dir, "range_cram")
  out_cram <- rduckhts_mosdepth(
    con,
    prefix_cram,
    cram_path,
    fasta = fasta_path,
    index_path = cram_index_path,
    fast_mode = TRUE,
    overwrite = TRUE
  )
  expect_true(isTRUE(out_cram$success[1]))
  summary_cram <- readLines(file.path(tmp_dir, "range_cram.mosdepth.summary.txt"))
  expect_equal(
    summary_cram[1:3],
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
    fast_mode = TRUE,
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

  prefix_thr <- file.path(tmp_dir, "range_thr")
  out_thr <- rduckhts_mosdepth(
    con,
    prefix_thr,
    bam_path,
    index_path = bam_index_path,
    by = "1000",
    thresholds = "1,2",
    fast_mode = TRUE,
    overwrite = TRUE
  )
  expect_true(isTRUE(out_thr$success[1]))
  expect_true(file.exists(file.path(tmp_dir, "range_thr.thresholds.bed.gz")))
  expect_true(file.exists(file.path(tmp_dir, "range_thr.thresholds.bed.gz.csi")))
  expect_true(!is.na(out_thr$thresholds_path[1]))
  thresholds_rows <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT chrom, region, ge1, ge2",
        "FROM read_csv(%s,",
        "  columns = {'chrom':'VARCHAR', 'start_pos':'BIGINT', 'end_pos':'BIGINT', 'region':'VARCHAR', 'ge1':'BIGINT', 'ge2':'BIGINT'},",
        "  delim = '\t',",
        "  header = FALSE,",
        "  skip = 1",
        ")",
        "WHERE chrom = 'CHROMOSOME_II'",
        "ORDER BY start_pos",
        "LIMIT 5"
      ),
      sprintf("'%s'", file.path(tmp_dir, "range_thr.thresholds.bed.gz"))
    )
  )
  expect_equal(
    thresholds_rows,
    data.frame(
      chrom = rep("CHROMOSOME_II", 5),
      region = rep("unknown", 5),
      ge1 = c(0, 840, 737, 82, 0),
      ge2 = c(0, 555, 368, 75, 0)
    )
  )

  big_bam_path <- system.file("extdata", "big.bam", package = "Rduckhts")
  big_index_path <- system.file("extdata", "big.bam.csi", package = "Rduckhts")
  prefix_big <- file.path(tmp_dir, "big_fast")
  out_big <- rduckhts_mosdepth(
    con,
    prefix_big,
    big_bam_path,
    index_path = big_index_path,
    fast_mode = TRUE,
    overwrite = TRUE
  )
  expect_true(isTRUE(out_big$success[1]))
  summary_big <- readLines(file.path(tmp_dir, "big_fast.mosdepth.summary.txt"))
  expect_equal(
    summary_big[1:3],
    c(
      "chrom\tlength\tbases\tmean\tmin\tmax",
      "ref\t600000045\t16\t0.00\t0\t1",
      "total\t600000045\t16\t0.00\t0\t1"
    )
  )
  expect_true(file.exists(file.path(tmp_dir, "big_fast.per-base.bed.gz.csi")))

  empty_bam_path <- system.file("extdata", "empty-tids.bam", package = "Rduckhts")
  empty_index_path <- system.file("extdata", "empty-tids.bam.bai", package = "Rduckhts")
  empty_bed_path <- system.file("extdata", "empty-tids.bed", package = "Rduckhts")
  prefix_empty <- file.path(tmp_dir, "empty_tids")
  out_empty <- rduckhts_mosdepth(
    con,
    prefix_empty,
    empty_bam_path,
    index_path = empty_index_path,
    by = "1000",
    no_per_base = TRUE,
    fast_mode = TRUE,
    overwrite = TRUE
  )
  expect_true(isTRUE(out_empty$success[1]))
  expect_true(is.na(out_empty$per_base_path[1]))
  summary_empty <- DBI::dbGetQuery(
    con,
    paste0(
      "SELECT chrom, length, bases, round(\"mean\", 2) AS mean_depth, max ",
      "FROM read_csv('", file.path(tmp_dir, "empty_tids.mosdepth.summary.txt"), "', delim = '\t', header = TRUE) ",
      "WHERE chrom = 'HPV18'"
    )
  )
  expect_equal(
    summary_empty,
    data.frame(
      chrom = "HPV18",
      length = 7857,
      bases = 1776620,
      mean_depth = 226.12,
      max = 647
    )
  )
  empty_windows <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT chrom, start_pos, end_pos, round(mean_depth, 2) AS mean_depth",
        "FROM read_csv(%s,",
        "  columns = {'chrom':'VARCHAR', 'start_pos':'BIGINT', 'end_pos':'BIGINT', 'mean_depth':'DOUBLE'},",
        "  delim = '\\t',",
        "  header = FALSE",
        ")",
        "WHERE chrom = 'CMV'",
        "ORDER BY start_pos",
        "LIMIT 3"
      ),
      sprintf("'%s'", file.path(tmp_dir, "empty_tids.regions.bed.gz"))
    )
  )
  expect_equal(
    empty_windows,
    data.frame(
      chrom = rep("CMV", 3),
      start_pos = c(0, 1000, 2000),
      end_pos = c(1000, 2000, 3000),
      mean_depth = c(0, 0, 0)
    )
  )

  prefix_bed <- file.path(tmp_dir, "empty_bed")
  out_bed <- rduckhts_mosdepth(
    con,
    prefix_bed,
    empty_bam_path,
    index_path = empty_index_path,
    by = empty_bed_path,
    no_per_base = TRUE,
    fast_mode = TRUE,
    overwrite = TRUE
  )
  expect_true(isTRUE(out_bed$success[1]))
  expect_true(file.exists(file.path(tmp_dir, "empty_bed.regions.bed.gz")))
  expect_true(file.exists(file.path(tmp_dir, "empty_bed.regions.bed.gz.csi")))
  bed_regions <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT chrom, start_pos, end_pos, round(mean_depth, 2) AS mean_depth",
        "FROM read_csv(%s,",
        "  columns = {'chrom':'VARCHAR', 'start_pos':'BIGINT', 'end_pos':'BIGINT', 'mean_depth':'DOUBLE'},",
        "  delim = '\\t',",
        "  header = FALSE",
        ")",
        "WHERE chrom = 'HPV18'"
      ),
      sprintf("'%s'", file.path(tmp_dir, "empty_bed.regions.bed.gz"))
    )
  )
  expect_equal(
    bed_regions,
    data.frame(
      chrom = "HPV18",
      start_pos = 0,
      end_pos = 7857,
      mean_depth = 226.12
    )
  )

  prefix_quant <- file.path(tmp_dir, "empty_quant")
  out_quant <- rduckhts_mosdepth(
    con,
    prefix_quant,
    empty_bam_path,
    index_path = empty_index_path,
    quantize = ":1:4:",
    no_per_base = TRUE,
    fast_mode = TRUE,
    overwrite = TRUE
  )
  expect_true(isTRUE(out_quant$success[1]))
  expect_true(file.exists(file.path(tmp_dir, "empty_quant.quantized.bed.gz")))
  expect_true(file.exists(file.path(tmp_dir, "empty_quant.quantized.bed.gz.csi")))
  expect_true(!is.na(out_quant$quantized_path[1]))
  quant_rows <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT chrom, start_pos, end_pos, bucket",
        "FROM read_csv(%s,",
        "  columns = {'chrom':'VARCHAR', 'start_pos':'BIGINT', 'end_pos':'BIGINT', 'bucket':'VARCHAR'},",
        "  delim = '\\t',",
        "  header = FALSE",
        ")",
        "WHERE chrom = 'CMV'"
      ),
      sprintf("'%s'", file.path(tmp_dir, "empty_quant.quantized.bed.gz"))
    )
  )
  expect_equal(
    quant_rows,
    data.frame(
      chrom = "CMV",
      start_pos = 0,
      end_pos = 235646,
      bucket = "0:1"
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
    fast_mode = TRUE,
    overwrite = TRUE
  )
  expect_true(is.na(out_noper$per_base_path[1]))
  expect_true(file.exists(file.path(tmp_dir, "range_noper.regions.bed.gz")))
  expect_false(file.exists(file.path(tmp_dir, "range_noper.per-base.bed.gz")))

  ovl_bam_path <- system.file("extdata", "overlapping-pairs.bam", package = "Rduckhts")
  ovl_index_path <- system.file("extdata", "overlapping-pairs.bam.bai", package = "Rduckhts")
  prefix_default <- file.path(tmp_dir, "range_default")
  out_default <- rduckhts_mosdepth(
    con,
    prefix_default,
    ovl_bam_path,
    index_path = ovl_index_path,
    chrom = "1",
    by = "1000",
    fast_mode = FALSE,
    overwrite = TRUE
  )
  expect_true(isTRUE(out_default$success[1]))
  summary_default <- readLines(file.path(tmp_dir, "range_default.mosdepth.summary.txt"))
  expect_equal(
    summary_default[1:4],
    c(
      "chrom\tlength\tbases\tmean\tmin\tmax",
      "1\t249250621\t80\t0.00\t0\t1",
      "1_region\t249250621\t80\t0.00\t0\t1",
      "total\t249250621\t80\t0.00\t0\t1"
    )
  )
  default_per_base <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT chrom, start_pos, end_pos, depth",
        "FROM read_csv(%s,",
        "  columns = {'chrom':'VARCHAR', 'start_pos':'BIGINT', 'end_pos':'BIGINT', 'depth':'INTEGER'},",
        "  delim = '\\t',",
        "  header = FALSE",
        ")",
        "LIMIT 3"
      ),
      sprintf("'%s'", file.path(tmp_dir, "range_default.per-base.bed.gz"))
    )
  )
  expect_equal(
    default_per_base,
    data.frame(
      chrom = rep("1", 3),
      start_pos = c(0, 565173, 565253),
      end_pos = c(565173, 565253, 249250621),
      depth = c(0, 1, 0)
    )
  )

  prefix_rg <- file.path(tmp_dir, "range_rg")
  out_rg <- rduckhts_mosdepth(
    con,
    prefix_rg,
    bam_path,
    index_path = bam_index_path,
    chrom = "CHROMOSOME_I",
    read_groups = "missing",
    no_per_base = TRUE,
    overwrite = TRUE
  )
  expect_true(isTRUE(out_rg$success[1]))
  summary_rg <- readLines(file.path(tmp_dir, "range_rg.mosdepth.summary.txt"))
  expect_equal(
    summary_rg[1:3],
    c(
      "chrom\tlength\tbases\tmean\tmin\tmax",
      "CHROMOSOME_I\t1009800\t0\t0.00\t0\t0",
      "total\t1009800\t0\t0.00\t0\t0"
    )
  )

  prefix_median <- file.path(tmp_dir, "range_median")
  out_median <- rduckhts_mosdepth(
    con,
    prefix_median,
    bam_path,
    index_path = bam_index_path,
    chrom = "CHROMOSOME_II",
    by = "1000",
    use_median = TRUE,
    no_per_base = TRUE,
    overwrite = TRUE
  )
  expect_true(isTRUE(out_median$success[1]))
  median_regions <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT chrom, start_pos, end_pos, value",
        "FROM read_csv(%s,",
        "  columns = {'chrom':'VARCHAR', 'start_pos':'BIGINT', 'end_pos':'BIGINT', 'value':'DOUBLE'},",
        "  delim = '\\t',",
        "  header = FALSE",
        ")",
        "ORDER BY start_pos"
      ),
      sprintf("'%s'", file.path(tmp_dir, "range_median.regions.bed.gz"))
    )
  )
  expect_equal(
    median_regions,
    data.frame(
      chrom = rep("CHROMOSOME_II", 5),
      start_pos = c(0, 1000, 2000, 3000, 4000),
      end_pos = c(1000, 2000, 3000, 4000, 5000),
      value = c(0, 2, 1, 0, 0)
    )
  )
  median_dist <- read.delim(
    file.path(tmp_dir, "range_median.mosdepth.region.dist.txt"),
    sep = "\t",
    header = FALSE,
    col.names = c("chrom", "depth", "fraction")
  )
  median_dist <- subset(median_dist, chrom == "CHROMOSOME_II")
  expect_equal(
    median_dist,
    data.frame(
      chrom = rep("CHROMOSOME_II", 3),
      depth = c(2L, 1L, 0L),
      fraction = c(0.20, 0.40, 1.00)
    )
  )

  prefix_fragment <- file.path(tmp_dir, "range_fragment")
  out_fragment <- rduckhts_mosdepth(
    con,
    prefix_fragment,
    bam_path,
    index_path = bam_index_path,
    chrom = "CHROMOSOME_II",
    by = "1000",
    fragment_mode = TRUE,
    no_per_base = TRUE,
    overwrite = TRUE
  )
  expect_true(isTRUE(out_fragment$success[1]))
  summary_fragment <- readLines(file.path(tmp_dir, "range_fragment.mosdepth.summary.txt"))
  expect_equal(
    summary_fragment[1:5],
    c(
      "chrom\tlength\tbases\tmean\tmin\tmax",
      "CHROMOSOME_II\t5000\t8032\t1.61\t0\t7",
      "CHROMOSOME_II_region\t5000\t8032\t1.61\t0\t7",
      "total\t5000\t8032\t1.61\t0\t7",
      "total_region\t5000\t8032\t1.61\t0\t7"
    )
  )
  fragment_regions <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT chrom, start_pos, end_pos, value",
        "FROM read_csv(%s,",
        "  columns = {'chrom':'VARCHAR', 'start_pos':'BIGINT', 'end_pos':'BIGINT', 'value':'DOUBLE'},",
        "  delim = '\\t',",
        "  header = FALSE",
        ")",
        "ORDER BY start_pos"
      ),
      sprintf("'%s'", file.path(tmp_dir, "range_fragment.regions.bed.gz"))
    )
  )
  expect_equal(
    fragment_regions,
    data.frame(
      chrom = rep("CHROMOSOME_II", 5),
      start_pos = c(0, 1000, 2000, 3000, 4000),
      end_pos = c(1000, 2000, 3000, 4000, 5000),
      value = c(0.00, 3.63, 3.05, 1.35, 0.00)
    )
  )
  fragment_dist <- read.delim(
    file.path(tmp_dir, "range_fragment.mosdepth.region.dist.txt"),
    sep = "\t",
    header = FALSE,
    col.names = c("chrom", "depth", "fraction")
  )
  fragment_dist <- subset(fragment_dist, chrom == "CHROMOSOME_II")
  expect_equal(
    fragment_dist,
    data.frame(
      chrom = rep("CHROMOSOME_II", 5),
      depth = c(4L, 3L, 2L, 1L, 0L),
      fraction = c(0.20, 0.40, 0.40, 0.60, 1.00)
    )
  )

  expect_error(
    rduckhts_mosdepth(
      con,
      file.path(tmp_dir, "range_badprec"),
      bam_path,
      index_path = bam_index_path,
      precision_digits = -1,
      overwrite = TRUE
    ),
    pattern = "precision_digits must be between 0 and 18"
  )

  expect_error(
    rduckhts_mosdepth(
      con,
      file.path(tmp_dir, "range_badthr"),
      bam_path,
      index_path = bam_index_path,
      thresholds = "1,2",
      overwrite = TRUE
    ),
    pattern = "thresholds requires by"
  )

  expect_error(
    rduckhts_mosdepth(
      con,
      file.path(tmp_dir, "range_badthrspec"),
      bam_path,
      index_path = bam_index_path,
      by = "1000",
      thresholds = "1,nope",
      overwrite = TRUE
    ),
    pattern = "invalid thresholds string"
  )

  expect_error(
    rduckhts_mosdepth(
      con,
      file.path(tmp_dir, "range_badquant"),
      bam_path,
      index_path = bam_index_path,
      quantize = "1::4",
      overwrite = TRUE
    ),
    pattern = "invalid quantize string"
  )

  expect_error(
    rduckhts_mosdepth(
      con,
      file.path(tmp_dir, "range_badmode"),
      bam_path,
      index_path = bam_index_path,
      fast_mode = TRUE,
      fragment_mode = TRUE,
      overwrite = TRUE
    ),
    pattern = "only one of fast_mode and fragment_mode can be TRUE"
  )
}

test_mosdepth()
