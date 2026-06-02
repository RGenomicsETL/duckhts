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
      "CREATE OR REPLACE TEMP TABLE norm_gvcf AS",
      "SELECT * FROM (VALUES",
      "('band_non_ref', 'chrS', 2, 'T', '<NON_REF>', 6::BIGINT),",
      "('mixed_non_ref', 'chrS', 2, 'T', 'TT,<NON_REF>', NULL::BIGINT),",
      "('mixed_symbolic_star', 'chrS', 2, 'T', 'TT,<*>', NULL::BIGINT),",
      "('star_only', 'chrS', 2, 'T', '*', NULL::BIGINT)",
      ") AS t(case_id, chrom, pos, ref, alt, end_pos)"
    )
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
  expect_equal(out_spanning$pos_normed[1], 1)
  expect_equal(out_spanning$end_pos_normed[1], 1)
  expect_equal(out_spanning$ref_normed[1], "G")
  expect_equal(as.character(out_spanning$alt_normed[[1]]), c("*", "GT"))
  expect_true(out_spanning$normed[1])
  expect_equal(out_spanning$norm_status[1], "Normalized")

  out_gvcf <- rduckhts_bcftools_norm(con, "norm_gvcf", fasta_path, end_pos_col = "end_pos")
  out_gvcf <- out_gvcf[order(out_gvcf$case_id), , drop = FALSE]
  expect_equal(out_gvcf$case_id, c("band_non_ref", "mixed_non_ref", "mixed_symbolic_star", "star_only"))
  expect_equal(out_gvcf$pos_normed, c(2, 1, 1, 2))
  expect_equal(out_gvcf$end_pos_normed, c(6, 1, 1, 2))
  expect_equal(out_gvcf$ref_normed, c("T", "G", "G", "T"))
  expect_equal(as.character(out_gvcf$alt_normed[[1]]), "<NON_REF>")
  expect_equal(as.character(out_gvcf$alt_normed[[2]]), c("GT", "<NON_REF>"))
  expect_equal(as.character(out_gvcf$alt_normed[[3]]), c("GT", "<*>"))
  expect_equal(as.character(out_gvcf$alt_normed[[4]]), "*")
  expect_equal(as.logical(out_gvcf$normed), c(FALSE, TRUE, TRUE, NA))
  expect_equal(out_gvcf$norm_status, c("GVCFReferenceBlock", "Normalized", "Normalized", "SpanningDeletion"))

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

  phased_fields_path <- system.file("extdata", "phased_genotype_fields.vcf", package = "Rduckhts", mustWork = TRUE)
  quoted_phased <- DBI::dbQuoteString(con, phased_fields_path)
  phased_gt <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT string_agg(SAMPLE_ID || '=' || FORMAT_GT || ':PS=' || coalesce(CAST(FORMAT_PS AS VARCHAR), 'NA'), ',' ORDER BY POS, SAMPLE_ID) AS gt ",
        "FROM read_bcf(%s, tidy_format := true)"
      ),
      quoted_phased
    )
  )
  expect_equal(
    phased_gt$gt[1],
    paste(
      "S1=0|1:PS=10,S2=1|2:PS=10,S1=1|0:PS=20,S2=0/1:PS=NA",
      "S1=1:PS=30,S2=0|1|1:PS=30,S1=0|1|1|2:PS=40,S2=2|2:PS=40",
      sep = ","
    )
  )
  phased_ploidy_lengths <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT string_agg(SAMPLE_ID || '=' || FORMAT_GT || ':PL=' || length(FORMAT_PL) || ':GP=' || length(FORMAT_GP) || ':DS=' || length(FORMAT_DS), ',' ORDER BY POS, SAMPLE_ID) AS lengths ",
        "FROM read_bcf(%s, tidy_format := true) ",
        "WHERE POS >= 30"
      ),
      quoted_phased
    )
  )
  expect_equal(
    phased_ploidy_lengths$lengths[1],
    "S1=1:PL=2:GP=2:DS=1,S2=0|1|1:PL=4:GP=4:DS=1,S1=0|1|1|2:PL=15:GP=15:DS=2,S2=2|2:PL=6:GP=6:DS=2"
  )
}

test_bcftools_norm()
