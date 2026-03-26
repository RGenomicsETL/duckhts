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

  unmapped <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('chrF', 11, 'A', 'T')) AS t(chrom, pos, ref, alt)",
    chain_path = chain_path,
    dst_fasta_ref = dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = src_fa
  )
  expect_equal(nrow(unmapped), 1)
  expect_false(unmapped$mapped[1])
  expect_equal(unmapped$warning[1], "UNMAPPED")
  expect_true(is.na(unmapped$dest_pos[1]))

  unmapped_indel <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('chrMissing', 2, 'AA', 'A')) AS t(chrom, pos, ref, alt)",
    chain_path = chain_path,
    dst_fasta_ref = dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = src_fa
  )
  expect_equal(nrow(unmapped_indel), 1)
  expect_false(unmapped_indel$mapped[1])
  expect_equal(unmapped_indel$warning[1], "UNMAPPED_ANCHORS")
  expect_true(is.na(unmapped_indel$dest_pos[1]))

  expect_error(
    rduckhts_liftover(
      con,
      query = "SELECT * FROM (VALUES ('chrF', 0, 'C', 'T')) AS t(chrom, pos, ref, alt)",
      chain_path = chain_path,
      dst_fasta_ref = dst_fa,
      ref_col = "ref",
      alt_col = "alt",
      src_fasta_ref = src_fa
    ),
    "pos must be >= 1"
  )

  expect_error(
    rduckhts_liftover(
      con,
      query = "SELECT * FROM (VALUES (NULL, 2, 'C', 'T')) AS t(chrom, pos, ref, alt)",
      chain_path = chain_path,
      dst_fasta_ref = dst_fa,
      ref_col = "ref",
      alt_col = "alt",
      src_fasta_ref = src_fa
    ),
    "chrom must be non-null"
  )

  expect_error(
    rduckhts_liftover(
      con,
      query = "SELECT * FROM (VALUES ('', 2, 'C', 'T')) AS t(chrom, pos, ref, alt)",
      chain_path = chain_path,
      dst_fasta_ref = dst_fa,
      ref_col = "ref",
      alt_col = "alt",
      src_fasta_ref = src_fa
    ),
    "chrom must be non-empty"
  )

  expect_error(
    rduckhts_liftover(
      con,
      query = "SELECT * FROM (VALUES ('chrF', NULL, 'C', 'T')) AS t(chrom, pos, ref, alt)",
      chain_path = chain_path,
      dst_fasta_ref = dst_fa,
      ref_col = "ref",
      alt_col = "alt",
      src_fasta_ref = src_fa
    ),
    "pos must be non-null"
  )

  expect_error(
    rduckhts_liftover(
      con,
      query = "SELECT * FROM (VALUES ('chrF', 2, 'C', 'T')) AS t(chrom, pos, ref, alt)",
      chain_path = chain_path,
      dst_fasta_ref = dst_fa,
      ref_col = "ref",
      alt_col = "alt",
      src_fasta_ref = src_fa,
      max_snp_gap = -1
    ),
    "max_snp_gap must be >= 0"
  )

  expect_error(
    rduckhts_liftover(
      con,
      query = "SELECT * FROM (VALUES ('chrF', 2, 'C', 'T')) AS t(chrom, pos, ref, alt)",
      chain_path = chain_path,
      dst_fasta_ref = dst_fa,
      ref_col = "ref",
      alt_col = "alt",
      src_fasta_ref = src_fa,
      max_indel_inc = -1
    ),
    "max_indel_inc must be >= 0"
  )

  expect_error(
    rduckhts_liftover(
      con,
      query = "SELECT * FROM (VALUES ('chrF', 2, 'C', 'T')) AS t(chrom, pos, ref, alt)",
      chain_path = file.path(tmp_dir, "missing.chain"),
      dst_fasta_ref = dst_fa,
      ref_col = "ref",
      alt_col = "alt",
      src_fasta_ref = src_fa
    ),
    "failed to load chain or FASTA context"
  )

  expect_error(
    rduckhts_liftover(
      con,
      query = "SELECT * FROM (VALUES ('chrF', 2, 'C', 'T')) AS t(chrom, pos, ref, alt)",
      chain_path = chain_path,
      dst_fasta_ref = file.path(tmp_dir, "missing.fa"),
      ref_col = "ref",
      alt_col = "alt",
      src_fasta_ref = src_fa
    ),
    "failed to load chain or FASTA context"
  )
}

test_liftover()
