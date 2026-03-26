library(tinytest)
library(DBI)

test_munge <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  expect_silent(rduckhts_load(con))

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
      ")).swapped AS swapped,",
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
  expect_false(scalar_out$swapped[1])
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
  expect_false(wrapper_rs2$swapped[1])
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
}

test_munge()
