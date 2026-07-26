library(tinytest)
library(DBI)

test_bcftools_norm <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))

  fasta_path <- system.file("extdata", "liftover_repeat_src.fa", package = "Rduckhts", mustWork = TRUE)
  quoted_fasta <- DBI::dbQuoteString(con, fasta_path)

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
      "('mixed_non_ref_band', 'chrS', 2, 'T', 'TT,<NON_REF>', 6::BIGINT),",
      "('mixed_symbolic_star', 'chrS', 2, 'T', 'TT,<*>', NULL::BIGINT),",
      "('star_only', 'chrS', 2, 'T', '*', NULL::BIGINT)",
      ") AS t(case_id, chrom, pos, ref, alt, end_pos)"
    )
  )
  DBI::dbExecute(
    con,
    paste(
      "CREATE OR REPLACE TEMP TABLE norm_fast_path AS",
      "SELECT * FROM (VALUES",
      "('plain_mnv_left_trim', 'chrS', 2, 'TT', 'TC', 6::BIGINT),",
      "('plain_mnv_no_common_head', 'chrS', 2, 'TT', 'CC', 6::BIGINT),",
      "('plain_snv_explicit_end', 'chrS', 2, 'T', 'C', 6::BIGINT)",
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
  expect_equal(out_gvcf$case_id, c("band_non_ref", "mixed_non_ref", "mixed_non_ref_band", "mixed_symbolic_star", "star_only"))
  expect_equal(out_gvcf$pos_normed, c(2, 1, 1, 1, 2))
  expect_equal(out_gvcf$end_pos_normed, c(6, 1, 6, 1, 2))
  expect_equal(out_gvcf$ref_normed, c("T", "G", "G", "G", "T"))
  expect_equal(as.character(out_gvcf$alt_normed[[1]]), "<NON_REF>")
  expect_equal(as.character(out_gvcf$alt_normed[[2]]), c("GT", "<NON_REF>"))
  expect_equal(as.character(out_gvcf$alt_normed[[3]]), c("GT", "<NON_REF>"))
  expect_equal(as.character(out_gvcf$alt_normed[[4]]), c("GT", "<*>"))
  expect_equal(as.character(out_gvcf$alt_normed[[5]]), "*")
  expect_equal(as.logical(out_gvcf$normed), c(FALSE, TRUE, TRUE, TRUE, NA))
  expect_equal(out_gvcf$norm_status, c("GVCFReferenceBlock", "Normalized", "Normalized", "Normalized", "SpanningDeletion"))

  out_fast_path <- rduckhts_bcftools_norm(con, "norm_fast_path", fasta_path, end_pos_col = "end_pos")
  out_fast_path <- out_fast_path[order(out_fast_path$case_id), , drop = FALSE]
  expect_equal(out_fast_path$case_id, c("plain_mnv_left_trim", "plain_mnv_no_common_head", "plain_snv_explicit_end"))
  expect_equal(out_fast_path$pos_normed, c(3, 2, 2))
  expect_equal(out_fast_path$end_pos_normed, c(3, 3, 2))
  expect_equal(out_fast_path$ref_normed, c("T", "TT", "T"))
  expect_equal(as.character(out_fast_path$alt_normed[[1]]), "C")
  expect_equal(as.character(out_fast_path$alt_normed[[2]]), "CC")
  expect_equal(as.character(out_fast_path$alt_normed[[3]]), "C")
  expect_equal(as.logical(out_fast_path$normed), c(TRUE, FALSE, FALSE))
  expect_equal(out_fast_path$norm_status, c("Normalized", "Unchanged", "Unchanged"))

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

  DBI::dbExecute(
    con,
    paste(
      "CREATE OR REPLACE TEMP TABLE norm_shadow_text AS",
      "SELECT * FROM (VALUES",
      "('keep_text', 'chrS', 2, 'T', 'TT', 'user_norm', 'user_item', 99)",
      ") AS t(case_id, chrom, pos, ref, alt, __duckhts_norm, __duckhts_alt_item, __duckhts_alt_index)"
    )
  )
  shadow_site <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT __duckhts_norm, pos_normed, ref_normed, alt_normed ",
        "FROM duckhts_bcftools_norm('norm_shadow_text', %s) ORDER BY case_id"
      ),
      quoted_fasta
    )
  )
  expect_equal(shadow_site$`__duckhts_norm`[1], "user_norm")
  expect_equal(shadow_site$pos_normed[1], 1)
  expect_equal(shadow_site$ref_normed[1], "G")
  expect_equal(as.character(shadow_site$alt_normed[[1]]), "GT")

  shadow_split <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT __duckhts_alt_item, __duckhts_alt_index, pos_normed, ref_normed, alt_normed, alt_index ",
        "FROM duckhts_bcftools_norm('norm_shadow_text', %s, split_multiallelic := true) ORDER BY case_id"
      ),
      quoted_fasta
    )
  )
  expect_equal(shadow_split$`__duckhts_alt_item`[1], "user_item")
  expect_equal(shadow_split$`__duckhts_alt_index`[1], 99)
  expect_equal(shadow_split$pos_normed[1], 1)
  expect_equal(shadow_split$ref_normed[1], "G")
  expect_equal(shadow_split$alt_normed[1], "GT")
  expect_equal(shadow_split$alt_index[1], 1)

  DBI::dbExecute(
    con,
    paste(
      "CREATE OR REPLACE TEMP TABLE norm_shadow_helpers AS",
      "SELECT * FROM (VALUES",
      "('chrS', 2, 'T', 'TT', 'bad_contig', 999, 'BADREF', 'BADALT', 123, 456)",
      ") AS t(chrom, pos, ref, alt, __duckhts_chrom, __duckhts_pos, __duckhts_ref, __duckhts_alt, __duckhts_end_pos, __duckhts_svlen)"
    )
  )
  shadow_helpers <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT __duckhts_chrom, __duckhts_pos, __duckhts_ref, __duckhts_alt, __duckhts_end_pos, __duckhts_svlen, ",
        "pos_normed, ref_normed, alt_normed, norm_status ",
        "FROM duckhts_bcftools_norm('norm_shadow_helpers', %s)"
      ),
      quoted_fasta
    )
  )
  expect_equal(shadow_helpers$`__duckhts_chrom`[1], "bad_contig")
  expect_equal(shadow_helpers$`__duckhts_pos`[1], 999)
  expect_equal(shadow_helpers$`__duckhts_ref`[1], "BADREF")
  expect_equal(shadow_helpers$`__duckhts_alt`[1], "BADALT")
  expect_equal(shadow_helpers$`__duckhts_end_pos`[1], 123)
  expect_equal(shadow_helpers$`__duckhts_svlen`[1], 456)
  expect_equal(shadow_helpers$pos_normed[1], 1)
  expect_equal(shadow_helpers$ref_normed[1], "G")
  expect_equal(as.character(shadow_helpers$alt_normed[[1]]), "GT")
  expect_equal(shadow_helpers$norm_status[1], "Normalized")

  shadow_helpers_split <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT __duckhts_chrom, __duckhts_pos, pos_normed, ref_normed, alt_normed, alt_index ",
        "FROM duckhts_bcftools_norm('norm_shadow_helpers', %s, split_multiallelic := true)"
      ),
      quoted_fasta
    )
  )
  expect_equal(shadow_helpers_split$`__duckhts_chrom`[1], "bad_contig")
  expect_equal(shadow_helpers_split$`__duckhts_pos`[1], 999)
  expect_equal(shadow_helpers_split$pos_normed[1], 1)
  expect_equal(shadow_helpers_split$ref_normed[1], "G")
  expect_equal(shadow_helpers_split$alt_normed[1], "GT")
  expect_equal(shadow_helpers_split$alt_index[1], 1)

  DBI::dbExecute(
    con,
    paste(
      "CREATE OR REPLACE TEMP TABLE norm_shadow_struct AS",
      "SELECT * FROM (VALUES",
      "('chrS', 2, 'T', 'TT', struct_pack(pos_normed := 999::BIGINT, end_pos_normed := 999::BIGINT, ref_normed := 'BAD', alt_normed := ['BAD'], normed := false, norm_status := 'BAD'))",
      ") AS t(chrom, pos, ref, alt, __duckhts_norm)"
    )
  )
  shadow_struct <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT __duckhts_norm.pos_normed AS user_pos_normed, pos_normed, ref_normed, alt_normed ",
        "FROM duckhts_bcftools_norm('norm_shadow_struct', %s)"
      ),
      quoted_fasta
    )
  )
  expect_equal(shadow_struct$user_pos_normed[1], 999)
  expect_equal(shadow_struct$pos_normed[1], 1)
  expect_equal(shadow_struct$ref_normed[1], "G")
  expect_equal(as.character(shadow_struct$alt_normed[[1]]), "GT")

  DBI::dbExecute(
    con,
    paste(
      "CREATE OR REPLACE TEMP TABLE norm_shadow_outputs AS",
      "SELECT * FROM (VALUES",
      "('chrS', 2, 'T', 'TT', 900::BIGINT, 901::BIGINT, 'USER_REF', ['USER_ALT']::VARCHAR[], false, 'USER_STATUS', 99::BIGINT)",
      ") AS t(chrom, pos, ref, alt, pos_normed, end_pos_normed, ref_normed, alt_normed, normed, norm_status, alt_index)"
    )
  )
  shadow_outputs <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT pos_normed, pos_normed_1, end_pos_normed, end_pos_normed_1, normed, normed_1, ",
        "ref_normed, ref_normed_1, alt_normed, alt_normed_1, norm_status, norm_status_1 ",
        "FROM duckhts_bcftools_norm('norm_shadow_outputs', %s)"
      ),
      quoted_fasta
    )
  )
  expect_equal(shadow_outputs$pos_normed[1], 900)
  expect_equal(shadow_outputs$pos_normed_1[1], 1)
  expect_equal(shadow_outputs$end_pos_normed[1], 901)
  expect_equal(shadow_outputs$end_pos_normed_1[1], 1)
  expect_false(shadow_outputs$normed[1])
  expect_true(shadow_outputs$normed_1[1])
  expect_equal(shadow_outputs$ref_normed[1], "USER_REF")
  expect_equal(shadow_outputs$ref_normed_1[1], "G")
  expect_equal(as.character(shadow_outputs$alt_normed[[1]]), "USER_ALT")
  expect_equal(as.character(shadow_outputs$alt_normed_1[[1]]), "GT")
  expect_equal(shadow_outputs$norm_status[1], "USER_STATUS")
  expect_equal(shadow_outputs$norm_status_1[1], "Normalized")

  shadow_outputs_split <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT alt_index, alt_index_1, pos_normed, pos_normed_1, end_pos_normed, end_pos_normed_1, ",
        "normed, normed_1, alt_normed, alt_normed_1, norm_status, norm_status_1 ",
        "FROM duckhts_bcftools_norm('norm_shadow_outputs', %s, split_multiallelic := true)"
      ),
      quoted_fasta
    )
  )
  expect_equal(shadow_outputs_split$alt_index[1], 99)
  expect_equal(shadow_outputs_split$alt_index_1[1], 1)
  expect_equal(shadow_outputs_split$pos_normed[1], 900)
  expect_equal(shadow_outputs_split$pos_normed_1[1], 1)
  expect_equal(shadow_outputs_split$end_pos_normed[1], 901)
  expect_equal(shadow_outputs_split$end_pos_normed_1[1], 1)
  expect_false(shadow_outputs_split$normed[1])
  expect_true(shadow_outputs_split$normed_1[1])
  expect_equal(as.character(shadow_outputs_split$alt_normed[[1]]), "USER_ALT")
  expect_equal(shadow_outputs_split$alt_normed_1[1], "GT")
  expect_equal(shadow_outputs_split$norm_status[1], "USER_STATUS")
  expect_equal(shadow_outputs_split$norm_status_1[1], "Normalized")

  DBI::dbExecute(
    con,
    paste(
      "CREATE OR REPLACE TEMP TABLE norm_plan_wide AS",
      "SELECT 'chrS'::VARCHAR AS chrom, 2::BIGINT AS pos, 'T'::VARCHAR AS ref, 'TT'::VARCHAR AS alt,",
      "repeat('x', 1000) AS payload, 'bad_contig' AS __duckhts_chrom"
    )
  )
  plan_site <- DBI::dbGetQuery(
    con,
    sprintf(
      "EXPLAIN SELECT chrom, pos_normed FROM duckhts_bcftools_norm('norm_plan_wide', %s)",
      quoted_fasta
    )
  )
  plan_site_text <- paste(plan_site$explain_value, collapse = "\n")
  expect_false(grepl("LEFT_DELIM_JOIN", plan_site_text, fixed = TRUE))
  expect_false(grepl("payload", plan_site_text, fixed = TRUE))
  expect_false(grepl("__duckhts_chrom", plan_site_text, fixed = TRUE))

  plan_split <- DBI::dbGetQuery(
    con,
    sprintf(
      "EXPLAIN SELECT * FROM duckhts_bcftools_norm('norm_seq', %s, split_multiallelic := true)",
      quoted_fasta
    )
  )
  plan_split_text <- paste(plan_split$explain_value, collapse = "\n")
  fixed_count <- function(pattern, text) {
    match_positions <- gregexpr(pattern, text, fixed = TRUE)[[1]]
    if (length(match_positions) == 1L && match_positions[1] == -1L) 0L else length(match_positions)
  }
  expect_equal(fixed_count("LEFT_DELIM_JOIN", plan_split_text), 1)
  expect_equal(fixed_count("bcftools_norm_row", plan_split_text), 1)

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
