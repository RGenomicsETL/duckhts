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
  vcf_dosage <- file.path(extdata, "score_dosage.vcf")
  vcf_chrprefix <- file.path(extdata, "score_chrprefix.vcf")
  vcf_missing_gt <- file.path(extdata, "score_missing_gt.vcf")
  sumf <- file.path(extdata, "score_summary.tsv")
  sumf_rsid <- file.path(extdata, "score_summary_rsid.tsv")
  sumf_no_snp <- file.path(extdata, "score_summary_no_snp.tsv")
  sumf_or <- file.path(extdata, "score_summary_or.tsv")
  sumf_plink2 <- file.path(extdata, "score_summary_plink2.tsv")
  sumf_na <- file.path(extdata, "score_summary_na.tsv")
  sumf_mismatch <- file.path(extdata, "score_summary_mismatch.tsv")
  sumf_custom <- file.path(extdata, "score_summary_custom.tsv")
  cols_file <- file.path(extdata, "score_columns.tsv")
  sumf_smallp <- file.path(extdata, "score_summary_smallp.tsv")

  expect_true(file.exists(vcf))
  expect_true(file.exists(vcf_dosage))
  expect_true(file.exists(vcf_chrprefix))
  expect_true(file.exists(vcf_missing_gt))
  expect_true(file.exists(sumf))
  expect_true(file.exists(sumf_rsid))
  expect_true(file.exists(sumf_no_snp))
  expect_true(file.exists(sumf_or))
  expect_true(file.exists(sumf_plink2))
  expect_true(file.exists(sumf_na))
  expect_true(file.exists(sumf_mismatch))
  expect_true(file.exists(sumf_custom))
  expect_true(file.exists(cols_file))
  expect_true(file.exists(sumf_smallp))

  # --- Basic GT scoring ---
  out_gt <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK")
  expect_true(all(c("SAMPLE", "score_summary") %in% names(out_gt)))
  expect_equal(out_gt$SAMPLE, c("S1", "S2"))
  expect_equal(round(out_gt$score_summary, 3), c(1.8, 0.1))

  # --- q_score_thr with counts ---
  # Column naming follows upstream: <prs>_CNT_p<thr> (not <prs>_p<thr>_CNT)
  # Boundary precision: strtof→double promotion means exact-boundary P values
  # are excluded (upstream parity: float LP < double threshold)
  out_thr <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                            q_score_thr = "0.01,0.2", counts = TRUE)
  expect_true(all(c("score_summary_p0.01", "score_summary_CNT_p0.01",
                     "score_summary_CNT_p0.2") %in% names(out_thr)))
  expect_equal(round(out_thr$score_summary_p0.01, 3), c(0.0, 0.0))
  expect_equal(out_thr$score_summary_CNT_p0.01, c(0, 0))
  expect_equal(out_thr$score_summary_CNT_p0.2, c(2, 2))

  # --- rsID matching ---
  out_rsid <- rduckhts_score(con, vcf, sumf_rsid, use = "GT", columns = "PLINK",
                             use_variant_id = TRUE)
  expect_equal(round(out_rsid$score_summary_rsid, 3), c(1.8, 0.1))

  # --- rsID summary without use_variant_id (chr:pos mismatch => zero) ---
  out_no_rsid <- rduckhts_score(con, vcf, sumf_rsid, use = "GT", columns = "PLINK")
  expect_equal(round(out_no_rsid$score_summary_rsid, 3), c(0, 0))

  # --- Error: use_variant_id with no SNP column ---
  expect_error(
    rduckhts_score(con, vcf, sumf_no_snp, use = "GT", columns = "PLINK",
                   use_variant_id = TRUE),
    "marker name column"
  )

  # --- Flexible chromosome name matching (chr prefix) ---
  out_chr <- rduckhts_score(con, vcf_chrprefix, sumf, use = "GT", columns = "PLINK")
  expect_equal(round(out_chr$score_summary, 3), c(1.8, 0.1))

  # --- OR-to-beta conversion ---
  out_or <- rduckhts_score(con, vcf, sumf_or, use = "GT", columns = "PLINK")
  expect_equal(round(out_or$score_summary_or, 3), c(1.8, 0.1))

  # --- PLINK2 preset with LOG10_P ---
  out_plink2 <- rduckhts_score(con, vcf, sumf_plink2, use = "GT", columns = "PLINK2")
  expect_equal(round(out_plink2$score_summary_plink2, 3), c(1.8, 0.1))

  # --- NA/missing value handling (rs2 BETA=NA, skipped) ---
  out_na <- rduckhts_score(con, vcf, sumf_na, use = "GT", columns = "PLINK")
  expect_equal(round(out_na$score_summary_na, 3), c(2.0, 0.5))

  # --- Missing genotype (S1 has ./. at rs1) ---
  out_miss <- rduckhts_score(con, vcf_missing_gt, sumf, use = "GT", columns = "PLINK")
  expect_equal(round(out_miss$score_summary, 3), c(1.8, 0.1))

  # --- Allele mismatch (all markers skipped => zero scores) ---
  out_mm <- rduckhts_score(con, vcf, sumf_mismatch, use = "GT", columns = "PLINK")
  expect_equal(round(out_mm$score_summary_mismatch, 3), c(0.0, 0.0))

  # --- Custom columns_file mapping ---
  out_custom <- rduckhts_score(con, vcf, sumf_custom, columns_file = cols_file)
  expect_equal(round(out_custom$score_summary_custom, 3), c(1.8, 0.1))

  # --- DS dosage mode ---
  out_ds <- rduckhts_score(con, vcf_dosage, sumf, use = "DS", columns = "PLINK")
  expect_equal(round(out_ds$score_summary, 3), c(1.69, 0.32))

  # --- HDS dosage mode ---
  out_hds <- rduckhts_score(con, vcf_dosage, sumf, use = "HDS", columns = "PLINK")
  expect_equal(round(out_hds$score_summary, 3), c(1.69, 0.32))

  # --- AP dosage mode ---
  out_ap <- rduckhts_score(con, vcf_dosage, sumf, use = "AP", columns = "PLINK")
  expect_equal(round(out_ap$score_summary, 3), c(1.69, 0.32))

  # --- GP dosage mode ---
  out_gp <- rduckhts_score(con, vcf_dosage, sumf, use = "GP", columns = "PLINK")
  expect_equal(round(out_gp$score_summary, 3), c(1.65, 0.32))

  # --- Auto-detection on GT-only VCF ---
  out_auto <- rduckhts_score(con, vcf, sumf, columns = "PLINK")
  expect_equal(round(out_auto$score_summary, 3), c(1.8, 0.1))

  # --- Auto-detection on dosage VCF (DS wins priority) ---
  out_auto_ds <- rduckhts_score(con, vcf_dosage, sumf, columns = "PLINK")
  expect_equal(round(out_auto_ds$score_summary, 3), c(1.69, 0.32))

  # --- Small p-value with threshold ---
  out_sp <- rduckhts_score(con, vcf, sumf_smallp, use = "GT", columns = "PLINK",
                           q_score_thr = "1e-200")
  thr_col <- grep("p1e-200$", names(out_sp), value = TRUE)
  expect_true(length(thr_col) > 0)
  expect_equal(round(out_sp[[thr_col[1]]], 3), c(0.0, 0.5))
}

test_score()
