library(tinytest)
library(DBI)
library(Rduckhts)

test_score <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  expect_silent(rduckhts_load(con))

  extdata <- system.file("extdata", package = "Rduckhts")
  vcf <- file.path(extdata, "score_input.vcf")
  sumf <- file.path(extdata, "score_summary.tsv")
  expect_true(file.exists(vcf))
  expect_true(file.exists(sumf))

  out_gt <- rduckhts_score(
    con,
    bcf_path = vcf,
    summary_path = sumf,
    use = "GT",
    columns = "PLINK"
  )

  expect_true(all(c("SAMPLE", "score_summary") %in% names(out_gt)))
  expect_equal(out_gt$SAMPLE, c("S1", "S2"))
  expect_equal(round(out_gt$score_summary, 3), c(1.8, 0.1))

  out_thr <- rduckhts_score(
    con,
    bcf_path = vcf,
    summary_path = sumf,
    use = "GT",
    columns = "PLINK",
    q_score_thr = "0.01,0.2",
    counts = TRUE
  )

  expect_true(all(c("score_summary_p0.01", "score_summary_p0.01_CNT", "score_summary_p0.2_CNT") %in% names(out_thr)))
  expect_equal(round(out_thr$score_summary_p0.01, 3), c(0.0, 0.5))
  expect_equal(out_thr$score_summary_p0.01_CNT, c(1, 1))
  expect_equal(out_thr$score_summary_p0.2_CNT, c(2, 2))

  out_ds <- rduckhts_score(con, vcf, sumf, use = "DS", columns = "PLINK")
  out_ap <- rduckhts_score(con, vcf, sumf, use = "AP", columns = "PLINK")
  out_gp <- rduckhts_score(con, vcf, sumf, use = "GP", columns = "PLINK")
  out_hds <- rduckhts_score(con, vcf, sumf, use = "HDS", columns = "PLINK")

  expect_equal(out_ds$score_summary, c(0, 0))
  expect_equal(out_ap$score_summary, c(0, 0))
  expect_equal(out_gp$score_summary, c(0, 0))
  expect_equal(out_hds$score_summary, c(0, 0))
}

test_score()
