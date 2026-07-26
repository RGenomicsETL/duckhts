library(tinytest)
library(DBI)
library(Rduckhts)

test_score <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))

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
  sumf_cnt <- file.path(extdata, "score_summary_CNT.tsv")
  sumf_mismatch <- file.path(extdata, "score_summary_mismatch.tsv")
  sumf_custom <- file.path(extdata, "score_summary_custom.tsv")
  cols_file <- file.path(extdata, "score_columns.tsv")
  sumf_smallp <- file.path(extdata, "score_summary_smallp.tsv")
  regions_file <- file.path(extdata, "score_regions.txt")
  targets_file <- file.path(extdata, "score_targets.txt")
  summaries_list_file <- file.path(extdata, "score_summaries.list")
  summaries_dir <- file.path(extdata, "score_summary_dir")

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
  expect_true(file.exists(sumf_cnt))
  expect_true(file.exists(sumf_mismatch))
  expect_true(file.exists(sumf_custom))
  expect_true(file.exists(cols_file))
  expect_true(file.exists(sumf_smallp))
  expect_true(file.exists(regions_file))
  expect_true(file.exists(targets_file))
  expect_true(file.exists(summaries_list_file))
  expect_true(dir.exists(summaries_dir))

  # --- Basic GT scoring ---
  out_gt <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK")
  expect_true(all(c("SAMPLE", "score_summary") %in% names(out_gt)))
  expect_equal(out_gt$SAMPLE, c("S1", "S2"))
  expect_equal(round(out_gt$score_summary, 3), c(1.8, 0.1))

  # --- Multi-PRS TSV scoring: character vector and summaries_list_file ---
  out_multi <- rduckhts_score(con, vcf, c(sumf, sumf_na), use = "GT", columns = "PLINK")
  expect_true(all(c("score_summary", "score_summary_na") %in% names(out_multi)))
  expect_equal(round(out_multi$score_summary, 3), c(1.8, 0.1))
  expect_equal(round(out_multi$score_summary_na, 3), c(2.0, 0.5))

  tmp_summaries_list <- tempfile("duckhts_score_summaries_", fileext = ".list")
  writeLines(c(sumf, sumf_na), tmp_summaries_list)
  out_multi_list <- rduckhts_score(con, vcf, summary_path = NULL,
                                   summaries_list_file = tmp_summaries_list,
                                   use = "GT", columns = "PLINK")
  expect_equal(round(out_multi_list$score_summary, 3), c(1.8, 0.1))
  expect_equal(round(out_multi_list$score_summary_na, 3), c(2.0, 0.5))

  old_wd <- setwd(extdata)
  on.exit(setwd(old_wd), add = TRUE)
  out_bundled_list <- rduckhts_score(con, vcf, summary_path = NULL,
                                     summaries_list_file = "score_summaries.list",
                                     use = "GT", columns = "PLINK")
  setwd(old_wd)
  expect_equal(round(out_bundled_list$score_summary, 3), c(1.8, 0.1))
  expect_equal(round(out_bundled_list$score_summary_na, 3), c(2.0, 0.5))

  out_dir_list <- rduckhts_score(con, vcf, summary_path = NULL,
                                 summaries_list_file = summaries_dir,
                                 use = "GT", columns = "PLINK")
  expect_equal(round(out_dir_list$score_summary, 3), c(1.8, 0.1))
  expect_equal(round(out_dir_list$score_summary_na, 3), c(2.0, 0.5))

  expect_error(
    rduckhts_score(con, vcf, c(sumf, sumf_cnt), use = "GT", columns = "PLINK", counts = TRUE),
    "duplicate generated output column name"
  )

  log_path <- tempfile("duckhts_score_", fileext = ".log")
  out_log <- rduckhts_score(con, vcf, c(sumf, sumf_mismatch), use = "GT",
                            columns = "PLINK", log_path = log_path)
  expect_true(file.exists(log_path))
  expect_equal(round(out_log$score_summary_mismatch, 3), c(0, 0))
  log_tbl <- utils::read.delim(log_path, comment.char = "#", check.names = FALSE)
  expect_true(all(c("summary_name", "loaded_markers", "matched_markers",
                    "allele_mismatch_markers") %in% names(log_tbl)))
  mismatch_row <- log_tbl[log_tbl$summary_name == "score_summary_mismatch", , drop = FALSE]
  expect_equal(mismatch_row$loaded_markers[1], 3L)
  expect_equal(mismatch_row$matched_markers[1], 0L)
  expect_equal(mismatch_row$allele_mismatch_markers[1], 3L)

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
    pattern = "failed to evaluate include expression"
  )

  out_format_gt <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                                  include = "FORMAT/GT==\"0/0\"")
  expect_equal(round(out_format_gt$score_summary, 3), c(2.0, 0.5))

  out_shift_like <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                                   include = "POS >> 100")
  expect_equal(round(out_shift_like$score_summary, 3), c(1.8, 0.1))

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

  ## ---- FORMAT/AS integer dosage tag ----
  vcf_as <- file.path(extdata, "score_as.vcf")
  expect_true(file.exists(vcf_as))
  # S1: 2*0.5 + 1*(-0.2) + 2*1.0 = 2.8; S2: 1*0.5 + 2*(-0.2) + 1*1.0 = 1.1
  out_as <- rduckhts_score(con, vcf_as, sumf, use = "AS", columns = "PLINK")
  expect_true("score_summary" %in% names(out_as))
  expect_equal(out_as$SAMPLE, c("S1", "S2"))
  expect_equal(round(out_as$score_summary, 3), c(2.8, 1.1))

  ## ---- Missing DS value (. skipped per variant) ----
  vcf_missing_ds <- file.path(extdata, "score_missing_ds.vcf")
  expect_true(file.exists(vcf_missing_ds))
  # S1: rs2 DS=missing → skip; 0.1*0.5 + 1.8*1.0 = 1.85; S2: 0.32
  out_mds <- rduckhts_score(con, vcf_missing_ds, sumf, use = "DS", columns = "PLINK")
  expect_equal(round(out_mds$score_summary, 3), c(1.85, 0.32))

  ## ---- REGENIE preset ----
  sumf_regenie <- file.path(extdata, "score_summary_regenie.tsv")
  expect_true(file.exists(sumf_regenie))
  out_regenie <- rduckhts_score(con, vcf, sumf_regenie, columns = "REGENIE")
  expect_equal(round(out_regenie$score_summary_regenie, 3), c(1.8, 0.1))

  ## ---- SAIGE preset ----
  sumf_saige <- file.path(extdata, "score_summary_saige.tsv")
  expect_true(file.exists(sumf_saige))
  out_saige <- rduckhts_score(con, vcf, sumf_saige, columns = "SAIGE")
  expect_equal(round(out_saige$score_summary_saige, 3), c(1.8, 0.1))

  ## ---- BOLT preset ----
  sumf_bolt <- file.path(extdata, "score_summary_bolt.tsv")
  expect_true(file.exists(sumf_bolt))
  out_bolt <- rduckhts_score(con, vcf, sumf_bolt, columns = "BOLT")
  expect_equal(round(out_bolt$score_summary_bolt, 3), c(1.8, 0.1))

  ## ---- METAL preset ----
  sumf_metal <- file.path(extdata, "score_summary_metal.tsv")
  expect_true(file.exists(sumf_metal))
  out_metal <- rduckhts_score(con, vcf, sumf_metal, columns = "METAL")
  expect_equal(round(out_metal$score_summary_metal, 3), c(1.8, 0.1))

  ## ---- PGS catalog preset ----
  sumf_pgs <- file.path(extdata, "score_summary_pgs.tsv")
  expect_true(file.exists(sumf_pgs))
  out_pgs <- rduckhts_score(con, vcf, sumf_pgs, columns = "PGS")
  expect_equal(round(out_pgs$score_summary_pgs, 3), c(1.8, 0.1))

  ## ---- SSF / GWAS-SSF preset (same file, two aliases) ----
  sumf_ssf <- file.path(extdata, "score_summary_ssf.tsv")
  expect_true(file.exists(sumf_ssf))
  out_ssf <- rduckhts_score(con, vcf, sumf_ssf, columns = "SSF")
  expect_equal(round(out_ssf$score_summary_ssf, 3), c(1.8, 0.1))
  out_gwas_ssf <- rduckhts_score(con, vcf, sumf_ssf, columns = "GWAS-SSF")
  expect_equal(round(out_gwas_ssf$score_summary_ssf, 3), c(1.8, 0.1))

  ## ---- LOG10_P under PLINK preset ----
  sumf_log10p <- file.path(extdata, "score_summary_log10p.tsv")
  expect_true(file.exists(sumf_log10p))
  out_log10p <- rduckhts_score(con, vcf, sumf_log10p, columns = "PLINK")
  expect_equal(round(out_log10p$score_summary_log10p, 3), c(1.8, 0.1))

  ## ---- TSV counts without threshold ----
  out_cnt <- rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK", counts = TRUE)
  expect_true(all(c("score_summary", "score_summary_CNT") %in% names(out_cnt)))
  expect_equal(round(out_cnt$score_summary, 3), c(1.8, 0.1))
  expect_equal(out_cnt$score_summary_CNT, c(3L, 3L))

  ## ---- GWAS-VCF + q_score_thr: only LP >= threshold variants scored ----
  # q_score_thr=0.1 → LP threshold 1.0; only rs1 qualifies
  out_gwas_thr <- rduckhts_score(con, vcf, sumf_gwas, use = "GT",
                                 q_score_thr = "0.1", counts = TRUE)
  expect_true(all(c("PRS_A_p0.1", "PRS_A_CNT_p0.1",
                    "PRS_B_p0.1", "PRS_B_CNT_p0.1") %in% names(out_gwas_thr)))
  thr_s1 <- out_gwas_thr[out_gwas_thr$SAMPLE == "S1", , drop = FALSE]
  thr_s2 <- out_gwas_thr[out_gwas_thr$SAMPLE == "S2", , drop = FALSE]
  expect_equal(round(thr_s1[["PRS_A_p0.1"]][1], 3), 0.0)
  expect_equal(round(thr_s1[["PRS_B_p0.1"]][1], 3), 0.0)
  expect_equal(thr_s1[["PRS_A_CNT_p0.1"]][1], 1L)
  expect_equal(round(thr_s2[["PRS_A_p0.1"]][1], 3), 0.5)
  expect_equal(round(thr_s2[["PRS_B_p0.1"]][1], 3), 0.3)
  expect_equal(thr_s2[["PRS_B_CNT_p0.1"]][1], 1L)

  ## ---- Error: include + exclude together ----
  expect_error(
    rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                   include = "POS > 100", exclude = "POS > 200"),
    pattern = "only one of include or exclude"
  )

  ## ---- Error: regions + regions_file together ----
  expect_error(
    rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                   regions = "1:100-200", regions_file = regions_file),
    pattern = "only one of regions or regions_file"
  )

  ## ---- Error: targets + targets_file together ----
  expect_error(
    rduckhts_score(con, vcf, sumf, use = "GT", columns = "PLINK",
                   targets = "1:100-200", targets_file = targets_file),
    pattern = "only one of targets or targets_file"
  )
}

test_score()
