library(tinytest)
library(DBI)

test_munge <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))

  tmp_dir <- tempfile("duckhts_munge_")
  dir.create(tmp_dir)

  fasta_ref <- file.path(tmp_dir, "munge.fa")
  writeLines(c(
    ">chrF",
    "ACGTACGTAA"
  ), fasta_ref)

  expect_true(rduckhts_fasta_index(con, fasta_ref, index_path = paste0(fasta_ref, ".fai"))$success[1])

  scalar_sql <- sprintf(
    paste(
      "SELECT",
      "(bcftools_munge_row(",
      "'chrF', 2, 'A', 'C', 'rs1',",
      "NULL, 2.0, NULL, 0.25,",
      "1000, NULL, NULL, NULL, 0.1, NULL, NULL, 200.0, NULL, NULL, NULL, NULL, NULL, NULL,",
      "'%s', 'IFFY', 'REF_MISMATCH', NULL, NULL, NULL",
      ")).chrom AS chrom,",
      "(bcftools_munge_row(",
      "'chrF', 2, 'A', 'C', 'rs1',",
      "NULL, 2.0, NULL, 0.25,",
      "1000, NULL, NULL, NULL, 0.1, NULL, NULL, 200.0, NULL, NULL, NULL, NULL, NULL, NULL,",
      "'%s', 'IFFY', 'REF_MISMATCH', NULL, NULL, NULL",
      ")).ref AS ref,",
      "(bcftools_munge_row(",
      "'chrF', 2, 'A', 'C', 'rs1',",
      "NULL, 2.0, NULL, 0.25,",
      "1000, NULL, NULL, NULL, 0.1, NULL, NULL, 200.0, NULL, NULL, NULL, NULL, NULL, NULL,",
      "'%s', 'IFFY', 'REF_MISMATCH', NULL, NULL, NULL",
      ")).alt AS alt,",
      "(bcftools_munge_row(",
      "'chrF', 2, 'A', 'C', 'rs1',",
      "NULL, 2.0, NULL, 0.25,",
      "1000, NULL, NULL, NULL, 0.1, NULL, NULL, 200.0, NULL, NULL, NULL, NULL, NULL, NULL,",
      "'%s', 'IFFY', 'REF_MISMATCH', NULL, NULL, NULL",
      ")).alleles_swapped AS alleles_swapped,",
      "(bcftools_munge_row(",
      "'chrF', 2, 'A', 'C', 'rs1',",
      "NULL, 2.0, NULL, 0.25,",
      "1000, NULL, NULL, NULL, 0.1, NULL, NULL, 200.0, NULL, NULL, NULL, NULL, NULL, NULL,",
      "'%s', 'IFFY', 'REF_MISMATCH', NULL, NULL, NULL",
      ")).af AS af,",
      "(bcftools_munge_row(",
      "'chrF', 2, 'A', 'C', 'rs1',",
      "NULL, 2.0, NULL, 0.25,",
      "1000, NULL, NULL, NULL, 0.1, NULL, NULL, 200.0, NULL, NULL, NULL, NULL, NULL, NULL,",
      "'%s', 'IFFY', 'REF_MISMATCH', NULL, NULL, NULL",
      ")).ac AS ac,",
      "(bcftools_munge_row(",
      "'chrF', 2, 'A', 'C', 'rs1',",
      "NULL, 2.0, NULL, 0.25,",
      "1000, NULL, NULL, NULL, 0.1, NULL, NULL, 200.0, NULL, NULL, NULL, NULL, NULL, NULL,",
      "'%s', 'IFFY', 'REF_MISMATCH', NULL, NULL, NULL",
      ")).es AS es,",
      "(bcftools_munge_row(",
      "'chrF', 2, 'A', 'C', 'rs1',",
      "NULL, 2.0, NULL, 0.25,",
      "1000, NULL, NULL, NULL, 0.1, NULL, NULL, 200.0, NULL, NULL, NULL, NULL, NULL, NULL,",
      "'%s', 'IFFY', 'REF_MISMATCH', NULL, NULL, NULL",
      ")).ez AS ez"
    ),
    fasta_ref, fasta_ref, fasta_ref, fasta_ref, fasta_ref, fasta_ref, fasta_ref, fasta_ref
  )
  scalar_out <- dbGetQuery(con, scalar_sql)
  expect_equal(scalar_out$chrom[1], "chrF")
  expect_equal(scalar_out$ref[1], "C")
  expect_equal(scalar_out$alt[1], "A")
  expect_false(scalar_out$alleles_swapped[1])
  expect_equal(round(scalar_out$af[1], 3), 0.1)
  expect_equal(round(scalar_out$ac[1], 1), 200)
  expect_equal(round(scalar_out$es[1], 3), 0.25)
  expect_equal(round(scalar_out$ez[1], 3), 2)

  wrapper_map <- rduckhts_munge(
    con,
    query = paste(
      "SELECT * FROM (VALUES",
      "('chrF', 2, 'A', 'C', 'rs2', 0.1, 1000),",
      "('chrF', 2, 'C', 'C', 'rs3', NULL, NULL)",
      ") AS t(CHR, BP, A1, A2, SNP, FRQ, N)"
    ),
    fasta_ref = fasta_ref,
    column_map = c(CHR = "CHR", BP = "BP", A1 = "A1", A2 = "A2", SNP = "SNP", FRQ = "FRQ", N = "N")
  )
  expect_equal(nrow(wrapper_map), 2)
  wrapper_rs2 <- wrapper_map[wrapper_map$id == "rs2", , drop = FALSE]
  wrapper_rs3 <- wrapper_map[wrapper_map$id == "rs3", , drop = FALSE]
  expect_equal(nrow(wrapper_rs2), 1)
  expect_equal(nrow(wrapper_rs3), 1)
  expect_equal(wrapper_rs2$ref[1], "C")
  expect_equal(wrapper_rs2$alt[1], "A")
  expect_false(wrapper_rs2$alleles_swapped[1])
  expect_equal(round(wrapper_rs2$af[1], 3), 0.1)
  expect_equal(round(wrapper_rs2$ns[1], 1), 1000)
  expect_equal(wrapper_rs3$filter[1], "IFFY")

  wrapper_default <- rduckhts_munge(
    con,
    query = "SELECT 'chrF' AS CHR, 2 AS BP, 'A' AS A1, 'C' AS A2, 'rs4' AS SNP",
    fasta_ref = fasta_ref
  )
  expect_equal(nrow(wrapper_default), 1)
  expect_equal(wrapper_default$ref[1], "C")
  expect_equal(wrapper_default$alt[1], "A")

  wrapper_preset <- rduckhts_munge(
    con,
    query = "SELECT 'rs5' AS SNP, 2 AS BP, 'chrF' AS CHR, 'A' AS A1, 'C' AS A2, 0.01 AS P, 0.1 AS FRQ, 1000 AS N",
    fasta_ref = fasta_ref,
    preset = "PLINK"
  )
  expect_equal(nrow(wrapper_preset), 1)
  expect_equal(round(wrapper_preset$lp[1], 3), 2)
  expect_equal(round(wrapper_preset$af[1], 3), 0.1)
  expect_equal(round(wrapper_preset$ns[1], 1), 1000)

  expect_error(
    rduckhts_munge(
      con,
      query = "SELECT 'chrF' AS CHR, 0 AS BP, 'A' AS A1, 'C' AS A2, 'rs6' AS SNP",
      fasta_ref = fasta_ref,
      column_map = c(CHR = "CHR", BP = "BP", A1 = "A1", A2 = "A2", SNP = "SNP")
    ),
    "pos must be >= 1"
  )

  expect_error(
    rduckhts_munge(
      con,
      query = "SELECT NULL AS CHR, 2 AS BP, 'A' AS A1, 'C' AS A2, 'rs7' AS SNP",
      fasta_ref = fasta_ref,
      column_map = c(CHR = "CHR", BP = "BP", A1 = "A1", A2 = "A2", SNP = "SNP")
    ),
    "chrom must be non-null"
  )

  expect_error(
    rduckhts_munge(
      con,
      query = "SELECT 'chrF' AS CHR, 2 AS BP, 'A' AS A1, 'C' AS A2, 'rs8' AS SNP",
      fasta_ref = file.path(tmp_dir, "missing.fa"),
      column_map = c(CHR = "CHR", BP = "BP", A1 = "A1", A2 = "A2", SNP = "SNP")
    ),
    "failed to load FASTA index"
  )

  expect_error(
    rduckhts_munge(
      con,
      query = "SELECT 'chrF' AS CHR, 2 AS BP, 'A' AS A1, 'C' AS A2, 'rs9' AS SNP",
      fasta_ref = fasta_ref,
      preset = "PLINK",
      column_map = c(CHR = "CHR", BP = "BP", A1 = "A1", A2 = "A2", SNP = "SNP")
    ),
    "specify only one of preset, column_map, or column_map_file"
  )

  expect_error(
    rduckhts_munge(
      con,
      query = "SELECT 'chrF' AS CHR, 2 AS BP, 'A' AS A1, 'C' AS A2, 'rs10' AS SNP",
      fasta_ref = fasta_ref,
      column_map = character()
    ),
    "column_map must be non-empty"
  )

  ## ---- METAL preset: SI/I2/CQ/ED output columns ----

  # METAL preset with heterogeneity fields: the R wrapper should detect
  # HET_I2/HET_P/HET_LP/DIRE keys and dispatch to duckdb_munge_metal,
  # returning SI, I2, CQ, ED columns.
  # At pos 2 on chrF (ref=C): Allele1=A (effect allele, A1) != ref,
  # Allele2=C (other allele, A2) = ref => no swap.
  metal_out <- rduckhts_munge(
    con,
    query = paste(
      "SELECT 'rs_m1' AS MarkerName, 2 AS Position, 'chrF' AS Chromosome,",
      "'A' AS Allele1, 'C' AS Allele2, 0.01 AS \"P-value\", 2.5 AS Zscore,",
      "0.3 AS Effect, 1000 AS N, 0.25 AS Freq1, 0.05 AS StdErr,",
      "NULL AS \"log(P)\", 500 AS Weight,",
      "75.2 AS HetISq, 0.001 AS HetPVal, NULL AS logHetP,",
      "'++-?' AS Direction"
    ),
    fasta_ref = fasta_ref,
    preset = "METAL"
  )
  expect_equal(nrow(metal_out), 1)
  # SI/I2/CQ/ED columns should be present (METAL dispatch)
  expect_true("si" %in% names(metal_out))
  expect_true("i2" %in% names(metal_out))
  expect_true("cq" %in% names(metal_out))
  expect_true("ed" %in% names(metal_out))
  # No swap: Allele1=A is the alt, Allele2=C is the ref
  expect_equal(metal_out$ref[1], "C")
  expect_equal(metal_out$alt[1], "A")
  expect_false(metal_out$alleles_swapped[1])
  # SI = NULL (METAL preset doesn't map INFO; no imputation info column)
  expect_true(is.na(metal_out$si[1]))
  # I2 = HetISq passthrough
  expect_equal(round(metal_out$i2[1], 1), 75.2)
  # CQ = -log10(HetPVal) since logHetP is NULL
  expect_equal(round(metal_out$cq[1], 1), 3.0)
  # ED = Direction passthrough (no swap, so unchanged)
  expect_equal(metal_out$ed[1], "++-?")
  # Base fields
  expect_equal(round(metal_out$es[1], 1), 0.3)
  expect_equal(round(metal_out$ez[1], 1), 2.5)
  expect_equal(round(metal_out$ns[1], 0), 1000)
  expect_equal(round(metal_out$ne[1], 0), 500)

  # METAL with allele swap: Allele1=C (effect allele, A1) = ref at pos 2 => swap
  metal_swap <- rduckhts_munge(
    con,
    query = paste(
      "SELECT 'rs_m2' AS MarkerName, 2 AS Position, 'chrF' AS Chromosome,",
      "'C' AS Allele1, 'A' AS Allele2, 0.01 AS \"P-value\", 2.5 AS Zscore,",
      "0.3 AS Effect, 1000 AS N, 0.25 AS Freq1, 0.05 AS StdErr,",
      "NULL AS \"log(P)\", 500 AS Weight,",
      "75.2 AS HetISq, NULL AS HetPVal, 7.5 AS logHetP,",
      "'++-?' AS Direction"
    ),
    fasta_ref = fasta_ref,
    preset = "METAL"
  )
  expect_equal(nrow(metal_swap), 1)
  # Alleles swapped: A1=C matched ref, so ref=C alt=A with swap=true
  expect_equal(metal_swap$ref[1], "C")
  expect_equal(metal_swap$alt[1], "A")
  expect_true(metal_swap$alleles_swapped[1])
  # ED direction flipped on swap: + -> -, - -> +, ? unchanged
  expect_equal(metal_swap$ed[1], "--+?")
  # CQ from logHetP (takes precedence over HetPVal)
  expect_equal(round(metal_swap$cq[1], 1), 7.5)
  # I2 passthrough (unaffected by swap)
  expect_equal(round(metal_swap$i2[1], 1), 75.2)

  # Non-METAL preset should NOT have SI/I2/CQ/ED columns
  plink_out <- rduckhts_munge(
    con,
    query = "SELECT 'rs_p1' AS SNP, 2 AS BP, 'chrF' AS CHR, 'A' AS A1, 'C' AS A2, 0.01 AS P",
    fasta_ref = fasta_ref,
    preset = "PLINK"
  )
  expect_false("si" %in% names(plink_out))
  expect_false("i2" %in% names(plink_out))
  expect_false("cq" %in% names(plink_out))
  expect_false("ed" %in% names(plink_out))

  # Custom column_map with metal keys should dispatch to duckdb_munge_metal
  # A1=A (effect allele), A2=C (other allele) at pos 2 (ref=C): A1 != ref => no swap
  custom_metal <- rduckhts_munge(
    con,
    query = paste(
      "SELECT 'chrF' AS chr_col, 2 AS pos_col, 'A' AS a1_col, 'C' AS a2_col,",
      "'rs_cm1' AS id_col, 0.88 AS info_col, 55.0 AS het_i2_col, '+--' AS dir_col"
    ),
    fasta_ref = fasta_ref,
    column_map = c(
      CHR = "chr_col", BP = "pos_col", A1 = "a1_col", A2 = "a2_col",
      SNP = "id_col", INFO = "info_col", HET_I2 = "het_i2_col", DIRE = "dir_col"
    )
  )
  expect_equal(nrow(custom_metal), 1)
  expect_true("si" %in% names(custom_metal))
  expect_true("i2" %in% names(custom_metal))
  expect_true("ed" %in% names(custom_metal))
  expect_equal(round(custom_metal$si[1], 2), 0.88)
  expect_equal(round(custom_metal$i2[1], 0), 55)
  # No swap, so direction is unchanged
  expect_equal(custom_metal$ed[1], "+--")
  expect_false(custom_metal$alleles_swapped[1])

  ## ---- M-M9: fai-only mode (fasta_ref = NULL) ----

  # Without fasta_ref, alleles pass through as-is (no ref matching, no swap)
  # A1 = effect allele (alt), A2 = other allele (ref): ref=A2, alt=A1
  fai_only_default <- rduckhts_munge(
    con,
    query = "SELECT 'chrF' AS CHR, 2 AS BP, 'A' AS A1, 'C' AS A2, 'rs_fai1' AS SNP, 0.01 AS P",
    column_map = c(CHR = "CHR", BP = "BP", A1 = "A1", A2 = "A2", SNP = "SNP", P = "P")
  )
  expect_equal(nrow(fai_only_default), 1)
  # In fai-only mode: ref=A2=C, alt=A1=A (no reference comparison, no swap)
  expect_equal(fai_only_default$ref[1], "C")
  expect_equal(fai_only_default$alt[1], "A")
  expect_false(fai_only_default$alleles_swapped[1])
  # No FASTA check → filter should be PASS (no IFFY/REF_MISMATCH)
  expect_true(is.na(fai_only_default$filter[1]) || fai_only_default$filter[1] != "IFFY")

  # fai-only with explicit fasta_ref = NULL
  fai_only_explicit <- rduckhts_munge(
    con,
    query = "SELECT 'chrF' AS CHR, 2 AS BP, 'A' AS A1, 'C' AS A2, 'rs_fai2' AS SNP",
    fasta_ref = NULL,
    column_map = c(CHR = "CHR", BP = "BP", A1 = "A1", A2 = "A2", SNP = "SNP")
  )
  expect_equal(nrow(fai_only_explicit), 1)
  expect_equal(fai_only_explicit$ref[1], "C")
  expect_equal(fai_only_explicit$alt[1], "A")
  expect_false(fai_only_explicit$alleles_swapped[1])

  # fai-only with preset: verifies preset column mapping works without FASTA
  # PLINK: A1 = effect allele (alt), A2 = other allele (ref)
  fai_only_preset <- rduckhts_munge(
    con,
    query = "SELECT 'rs_fai3' AS SNP, 2 AS BP, 'chrF' AS CHR, 'A' AS A1, 'C' AS A2, 0.01 AS P, 0.1 AS FRQ, 1000 AS N",
    preset = "PLINK"
  )
  expect_equal(nrow(fai_only_preset), 1)
  expect_equal(fai_only_preset$ref[1], "C")
  expect_equal(fai_only_preset$alt[1], "A")
  expect_false(fai_only_preset$alleles_swapped[1])
  expect_equal(round(fai_only_preset$lp[1], 3), 2)
}

test_munge()
