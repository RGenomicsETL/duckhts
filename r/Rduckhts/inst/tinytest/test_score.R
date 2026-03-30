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
  regions_file <- file.path(extdata, "score_regions.txt")
  targets_file <- file.path(extdata, "score_targets.txt")

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
  expect_true(file.exists(regions_file))
  expect_true(file.exists(targets_file))

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

  # --- Sample subsetting: single sample ---
  out_sub <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                            samples = "S1")
  expect_equal(nrow(out_sub), 1L)
  expect_equal(out_sub$SAMPLE, "S1")
  expect_equal(round(out_sub$score_summary, 3), 1.8)

  # --- Sample subsetting: missing sample errors ---
  expect_error(
    rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                   samples = "S1,NOSAMPLE"),
    pattern = "sample"
  )

  # --- Sample subsetting: force_samples ignores missing ---
  out_force <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                              samples = "S1,NOSAMPLE", force_samples = TRUE)
  expect_equal(nrow(out_force), 1L)
  expect_equal(out_force$SAMPLE, "S1")
  expect_equal(round(out_force$score_summary, 3), 1.8)

  # --- Region/target/filter controls ---
  out_regions <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                                regions = "1:150-250")
  expect_equal(round(out_regions$score_summary, 3), c(-0.2, -0.4))

  out_targets <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                                targets = "1:150-250")
  expect_equal(round(out_targets$score_summary, 3), c(-0.2, -0.4))

  out_regions_file <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                                     regions_file = regions_file)
  expect_equal(round(out_regions_file$score_summary, 3), c(-0.2, -0.4))

  out_targets_file <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                                     targets_file = targets_file)
  expect_equal(round(out_targets_file$score_summary, 3), c(-0.2, -0.4))

  out_pass <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                             apply_filters = "PASS")
  expect_equal(round(out_pass$score_summary, 3), c(1.8, 0.1))

  expect_error(
    rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                   include = "POS >>> 100"),
    pattern = "failed to parse include expression"
  )

  out_include <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                                include = "POS > 100")
  expect_equal(round(out_include$score_summary, 3), c(1.8, -0.4))

  out_exclude <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                                exclude = "POS > 100")
  expect_equal(round(out_exclude$score_summary, 3), c(0.0, 0.5))

  out_type <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                             include = "TYPE==\"SNP\" && N_ALT==1")
  expect_equal(round(out_type$score_summary, 3), c(1.8, 0.1))

  out_filter_expr <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                                    include = "FILTER==\"PASS\"")
  expect_equal(round(out_filter_expr$score_summary, 3), c(1.8, 0.1))

  expect_error(
    rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                   include = "INFO/ZZZ==1"),
    pattern = "failed to parse include expression"
  )

  out_format_gt <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                                  include = "FORMAT/GT==\"0/0\"")
  expect_equal(round(out_format_gt$score_summary, 3), c(0.0, 0.5))

  ## ---- S-C4: GWAS-VCF multi-PRS scoring ----

  sumf_gwas <- file.path(extdata, "score_gwas_summary.vcf")
  expect_true(file.exists(sumf_gwas))

  # Basic GWAS-VCF: auto-detected multi-PRS (PRS_A, PRS_B columns in output)
  out_gwas <- rduckhts_score(con, vcf, sumf_gwas, use = "GT")
  expect_true(all(c("SAMPLE", "PRS_A", "PRS_B") %in% names(out_gwas)))
  expect_equal(nrow(out_gwas), 2L)
  gwas_s1 <- out_gwas[out_gwas$SAMPLE == "S1", , drop = FALSE]
  gwas_s2 <- out_gwas[out_gwas$SAMPLE == "S2", , drop = FALSE]
  expect_equal(round(gwas_s1$PRS_A[1], 3), 1.8)
  expect_equal(round(gwas_s1$PRS_B[1], 3), 1.0)
  expect_equal(round(gwas_s2$PRS_A[1], 3), 0.1)
  expect_equal(round(gwas_s2$PRS_B[1], 3), 0.3)

  # GWAS-VCF with counts: adds PRS_A_CNT, PRS_B_CNT columns
  out_gwas_cnt <- rduckhts_score(con, vcf, sumf_gwas, use = "GT", counts = TRUE)
  expect_true(all(c("PRS_A", "PRS_A_CNT", "PRS_B", "PRS_B_CNT") %in% names(out_gwas_cnt)))
  gwas_cnt_s1 <- out_gwas_cnt[out_gwas_cnt$SAMPLE == "S1", , drop = FALSE]
  gwas_cnt_s2 <- out_gwas_cnt[out_gwas_cnt$SAMPLE == "S2", , drop = FALSE]
  expect_equal(gwas_cnt_s1$PRS_A_CNT[1], 3)
  expect_equal(gwas_cnt_s1$PRS_B_CNT[1], 2)
  expect_equal(gwas_cnt_s2$PRS_A_CNT[1], 3)
  expect_equal(gwas_cnt_s2$PRS_B_CNT[1], 2)

  # GWAS-VCF with sample subsetting
  out_gwas_sub <- rduckhts_score(con, vcf, sumf_gwas, use = "GT", samples = "S1")
  expect_equal(nrow(out_gwas_sub), 1L)
  expect_equal(out_gwas_sub$SAMPLE, "S1")
  expect_equal(round(out_gwas_sub$PRS_A[1], 3), 1.8)
  expect_equal(round(out_gwas_sub$PRS_B[1], 3), 1.0)
}

test_score()
