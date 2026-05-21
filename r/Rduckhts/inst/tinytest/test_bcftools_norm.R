library(tinytest)
library(DBI)

test_bcftools_norm <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  expect_silent(rduckhts_load(con))

  fasta_path <- system.file("extdata", "liftover_repeat_src.fa", package = "Rduckhts", mustWork = TRUE)

  DBI::dbExecute(
    con,
    paste(
      "CREATE OR REPLACE TEMP TABLE norm_seq AS",
      "SELECT * FROM (VALUES",
      "('chrS', 2, 'T', 'TT,TTT'),",
      "('chrS', 2, 'TT', 'T,TTT'),",
      "('chrMissing', 2, 'T', 'TT')",
      ") AS t(chrom, pos, ref, alt)"
    )
  )
  DBI::dbExecute(
    con,
    "CREATE OR REPLACE TEMP TABLE norm_list AS SELECT * FROM (VALUES ('chrS', 2, 'T', ['TT', 'TTT'])) AS t(chrom, pos, ref, alt)"
  )
  DBI::dbExecute(
    con,
    paste(
      "CREATE OR REPLACE TEMP TABLE norm_sv AS",
      "SELECT * FROM (VALUES",
      "('chrS', 2, 'T', '<DEL>', 3, NULL),",
      "('chrS', 3, 'T', '<DUP>', NULL, 2)",
      ") AS t(chrom, pos, ref, alt, end_pos, svlen)"
    )
  )

  out <- rduckhts_bcftools_norm(con, "norm_seq", fasta_path)
  expect_equal(nrow(out), 3)

  row_multi <- out[out$ref == "T" & out$alt == "TT,TTT", , drop = FALSE]
  expect_equal(row_multi$pos_normed[1], 1)
  expect_equal(row_multi$end_pos_normed[1], 1)
  expect_equal(row_multi$ref_normed[1], "G")
  expect_equal(as.character(row_multi$alt_normed[[1]]), c("GT", "GTT"))
  expect_true(row_multi$normed[1])
  expect_equal(row_multi$norm_status[1], "Normalized")

  row_missing <- out[out$chrom == "chrMissing", , drop = FALSE]
  expect_equal(row_missing$pos_normed[1], 2)
  expect_equal(row_missing$end_pos_normed[1], 2)
  expect_equal(row_missing$ref_normed[1], "T")
  expect_equal(as.character(row_missing$alt_normed[[1]]), "TT")
  expect_true(is.na(row_missing$normed[1]))
  expect_equal(row_missing$norm_status[1], "MissingContig")

  out_split <- rduckhts_bcftools_norm(con, "norm_seq", fasta_path, split_multiallelic = TRUE)
  row_split <- out_split[out_split$ref == "T" & out_split$alt == "TT,TTT", , drop = FALSE]
  row_split <- row_split[order(row_split$alt_index), , drop = FALSE]
  expect_equal(nrow(row_split), 2)
  expect_equal(row_split$pos_normed, c(1, 1))
  expect_equal(row_split$ref_normed, c("G", "G"))
  expect_equal(row_split$alt_normed, c("GT", "GTT"))
  expect_equal(row_split$alt_index, c(1, 2))
  expect_true(all(row_split$normed))

  out_list <- rduckhts_bcftools_norm(con, "norm_list", fasta_path)
  expect_equal(nrow(out_list), 1)
  expect_equal(out_list$pos_normed[1], 1)
  expect_equal(out_list$ref_normed[1], "G")
  expect_equal(as.character(out_list$alt_normed[[1]]), c("GT", "GTT"))
  expect_true(out_list$normed[1])

  out_sv <- rduckhts_bcftools_norm(
    con,
    "norm_sv",
    fasta_path,
    end_pos_col = "end_pos",
    svlen_col = "svlen"
  )
  row_del <- out_sv[out_sv$alt == "<DEL>", , drop = FALSE]
  row_dup <- out_sv[out_sv$alt == "<DUP>", , drop = FALSE]
  expect_equal(row_del$pos_normed[1], 1)
  expect_equal(row_del$end_pos_normed[1], 2)
  expect_equal(row_del$ref_normed[1], "G")
  expect_equal(as.character(row_del$alt_normed[[1]]), "<DEL>")
  expect_true(row_del$normed[1])
  expect_equal(row_dup$pos_normed[1], 1)
  expect_equal(row_dup$end_pos_normed[1], 1)
  expect_equal(row_dup$ref_normed[1], "G")
  expect_equal(as.character(row_dup$alt_normed[[1]]), "<DUP>")
  expect_true(row_dup$normed[1])

  out_query <- rduckhts_bcftools_norm(
    con,
    "SELECT chrom, pos, ref, alt FROM norm_list",
    fasta_path
  )
  expect_equal(nrow(out_query), 1)
  expect_equal(out_query$ref_normed[1], "G")
}

test_bcftools_norm()
