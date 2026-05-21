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
  DBI::dbExecute(
    con,
    "CREATE OR REPLACE TEMP TABLE norm_spanning AS SELECT * FROM (VALUES ('chrS', 2, 'T', '*,TT')) AS t(chrom, pos, ref, alt)"
  )
  DBI::dbExecute(
    con,
    paste(
      "CREATE OR REPLACE TEMP TABLE norm_empty_seq AS",
      "SELECT * FROM (VALUES",
      "('dot', 'chrS', 2, 'T', '.'),",
      "('empty_text', 'chrS', 2, 'T', ''),",
      "('null_text', 'chrS', 2, 'T', NULL)",
      ") AS t(case_id, chrom, pos, ref, alt)"
    )
  )
  DBI::dbExecute(
    con,
    paste(
      "CREATE OR REPLACE TEMP TABLE norm_empty_list AS",
      "SELECT * FROM (VALUES",
      "('empty_list', 'chrS', 2, 'T', []::VARCHAR[]),",
      "('null_item', 'chrS', 2, 'T', [NULL]::VARCHAR[]),",
      "('null_list', 'chrS', 2, 'T', NULL::VARCHAR[])",
      ") AS t(case_id, chrom, pos, ref, alt)"
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

  empty_seq_split <- rduckhts_bcftools_norm(con, "norm_empty_seq", fasta_path, split_multiallelic = TRUE)
  empty_seq_split <- empty_seq_split[order(empty_seq_split$case_id), , drop = FALSE]
  expect_equal(empty_seq_split$case_id, c("dot", "empty_text", "null_text"))
  expect_equal(nrow(empty_seq_split), 3)
  expect_true(all(is.na(empty_seq_split$alt_normed[c(1, 3)])))
  expect_equal(nchar(empty_seq_split$alt_normed[2]), 0)
  expect_equal(empty_seq_split$norm_status, c("RefOnly", "NullInput", "RefOnly"))
  expect_equal(as.logical(empty_seq_split$normed), c(FALSE, NA, FALSE))
  expect_true(all(is.na(empty_seq_split$alt_index[c(1, 3)])))
  expect_equal(empty_seq_split$alt_index[2], 1)

  empty_list_split <- rduckhts_bcftools_norm(con, "norm_empty_list", fasta_path, split_multiallelic = TRUE)
  empty_list_split <- empty_list_split[order(empty_list_split$case_id), , drop = FALSE]
  expect_equal(empty_list_split$case_id, c("empty_list", "null_item", "null_list"))
  expect_equal(nrow(empty_list_split), 3)
  expect_true(all(is.na(empty_list_split$alt_normed)))
  expect_equal(empty_list_split$norm_status, c("RefOnly", "NullInput", "RefOnly"))
  expect_equal(as.logical(empty_list_split$normed), c(FALSE, NA, FALSE))
  expect_true(all(is.na(empty_list_split$alt_index[c(1, 3)])))
  expect_equal(empty_list_split$alt_index[2], 1)

  out_spanning <- rduckhts_bcftools_norm(con, "norm_spanning", fasta_path)
  expect_equal(nrow(out_spanning), 1)
  expect_equal(out_spanning$pos_normed[1], 2)
  expect_equal(out_spanning$end_pos_normed[1], 2)
  expect_equal(out_spanning$ref_normed[1], "T")
  expect_equal(as.character(out_spanning$alt_normed[[1]]), c("*", "TT"))
  expect_true(is.na(out_spanning$normed[1]))
  expect_equal(out_spanning$norm_status[1], "SpanningDeletion")

  out_spanning_split <- rduckhts_bcftools_norm(con, "norm_spanning", fasta_path, split_multiallelic = TRUE)
  out_spanning_split <- out_spanning_split[order(out_spanning_split$alt_index), , drop = FALSE]
  expect_equal(nrow(out_spanning_split), 2)
  expect_equal(out_spanning_split$alt_normed, c("*", "GT"))
  expect_equal(out_spanning_split$ref_normed, c("T", "G"))
  expect_equal(out_spanning_split$pos_normed, c(2, 1))
  expect_equal(out_spanning_split$norm_status, c("SpanningDeletion", "Normalized"))
  expect_equal(as.logical(out_spanning_split$normed), c(NA, TRUE))
  expect_equal(out_spanning_split$alt_index, c(1, 2))

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
