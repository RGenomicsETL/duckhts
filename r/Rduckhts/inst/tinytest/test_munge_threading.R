library(tinytest)
library(DBI)

test_munge_threading <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))

  tmp_dir <- tempfile("duckhts_munge_threading_")
  dir.create(tmp_dir)

  fasta_ref <- file.path(tmp_dir, "munge.fa")
  writeLines(c(
    ">chrF",
    "ACGTACGTAA"
  ), fasta_ref)

  expect_true(rduckhts_fasta_index(con, fasta_ref, index_path = paste0(fasta_ref, ".fai"))$success[1])
  expect_silent(dbExecute(con, "PRAGMA threads=4"))

  out <- rduckhts_munge(
    con,
    query = paste(
      "SELECT 'chrF' AS CHR,",
      "2::BIGINT AS BP,",
      "CASE WHEN i % 2 = 0 THEN 'A' ELSE 'C' END AS A1,",
      "CASE WHEN i % 2 = 0 THEN 'C' ELSE 'A' END AS A2,",
      "'rs' || CAST(i AS VARCHAR) AS SNP,",
      "0.1::DOUBLE AS FRQ,",
      "1000::DOUBLE AS N",
      "FROM range(50000) t(i)"
    ),
    fasta_ref = fasta_ref,
    column_map = c(CHR = "CHR", BP = "BP", A1 = "A1", A2 = "A2", SNP = "SNP", FRQ = "FRQ", N = "N")
  )

  expect_equal(nrow(out), 50000)
  expect_equal(sum(out$alleles_swapped), 25000)
  expect_equal(sum(is.na(out$filter)), 50000)
}

test_munge_threading()
