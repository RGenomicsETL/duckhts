library(tinytest)
library(DBI)

test_liftover <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  expect_silent(rduckhts_load(con))

  tmp_dir <- tempfile("duckhts_liftover_")
  dir.create(tmp_dir)

  src_fa <- file.path(tmp_dir, "liftover_src.fa")
  dst_fa <- file.path(tmp_dir, "liftover_dst.fa")
  chain_path <- file.path(tmp_dir, "liftover.chain")

  writeLines(c(
    ">chrF",
    "ACGTACGTAA",
    ">chrR",
    "AACCGGTTAA"
  ), src_fa)
  writeLines(c(
    ">chrLiftF",
    "ACGTACGTAA",
    ">chrLiftR",
    "TTAACCGGTT"
  ), dst_fa)
  writeLines(c(
    "chain 100 chrF 10 + 0 10 chrLiftF 10 + 0 10 1",
    "10",
    "",
    "chain 100 chrR 10 + 0 10 chrLiftR 10 - 0 10 2",
    "10"
  ), chain_path)

  expect_true(rduckhts_fasta_index(con, src_fa, index_path = paste0(src_fa, ".fai"))$success[1])
  expect_true(rduckhts_fasta_index(con, dst_fa, index_path = paste0(dst_fa, ".fai"))$success[1])

  query <- paste(
    "SELECT * FROM (VALUES",
    "('chrF', 2, 'C', 'T'),",
    "('chrR', 2, 'A', 'G'),",
    "('chrF', 2, NULL, 'T')",
    ") AS t(chrom, pos, ref, alt)"
  )
  out <- rduckhts_liftover(
    con,
    query = query,
    chain_path = chain_path,
    dst_fasta_ref = dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = src_fa
  )

  expect_equal(nrow(out), 3)
  row_forward <- out[out$src_chrom == "chrF" & !is.na(out$src_ref), , drop = FALSE]
  row_reverse <- out[out$src_chrom == "chrR", , drop = FALSE]
  row_iffy <- out[out$src_chrom == "chrF" & is.na(out$src_ref), , drop = FALSE]

  expect_equal(nrow(row_forward), 1)
  expect_equal(row_forward$dest_chrom[1], "chrLiftF")
  expect_equal(row_forward$dest_pos[1], 2)
  expect_equal(row_forward$dest_ref[1], "C")
  expect_equal(row_forward$dest_alt[1], "T")
  expect_false(row_forward$reverse_complemented[1])

  expect_equal(nrow(row_reverse), 1)
  expect_equal(row_reverse$dest_chrom[1], "chrLiftR")
  expect_equal(row_reverse$dest_pos[1], 9)
  expect_equal(row_reverse$dest_ref[1], "T")
  expect_equal(row_reverse$dest_alt[1], "C")
  expect_true(row_reverse$reverse_complemented[1])

  expect_equal(nrow(row_iffy), 1)
  expect_true(row_iffy$mapped[1])
  expect_equal(row_iffy$warning[1], "IFFY")
}

test_liftover()
