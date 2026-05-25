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

  alias_src_fa <- file.path(tmp_dir, "liftover_chr23_alias_src.fa")
  alias_dst_fa <- file.path(tmp_dir, "liftover_chr23_alias_dst.fa")
  alias_chain_path <- file.path(tmp_dir, "liftover_chr23_alias.chain")
  writeLines(c(
    ">chrX",
    "ACGTACGTAA"
  ), alias_src_fa)
  writeLines(c(
    ">chrLiftX",
    "ACGTACGTAA"
  ), alias_dst_fa)
  writeLines(c(
    "chain 100 X 10 + 0 10 chrLiftX 10 + 0 10 1",
    "10"
  ), alias_chain_path)
  expect_true(rduckhts_fasta_index(con, alias_src_fa, index_path = paste0(alias_src_fa, ".fai"))$success[1])
  expect_true(rduckhts_fasta_index(con, alias_dst_fa, index_path = paste0(alias_dst_fa, ".fai"))$success[1])

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
  row_missing_ref <- out[out$src_chrom == "chrF" & is.na(out$src_ref), , drop = FALSE]

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

  expect_equal(nrow(row_missing_ref), 1)
  expect_true(row_missing_ref$mapped[1])
  expect_true(is.na(row_missing_ref$reject_reason[1]))
  expect_equal(row_missing_ref$note[1], "MissingSourceRef")

  alias_out <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('23', 2, 'CG', 'C')) AS t(chrom, pos, ref, alt)",
    chain_path = alias_chain_path,
    dst_fasta_ref = alias_dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = alias_src_fa
  )
  expect_equal(nrow(alias_out), 1)
  expect_true(alias_out$mapped[1])
  expect_equal(alias_out$dest_chrom[1], "chrLiftX")
  expect_equal(alias_out$dest_pos[1], 2)
  expect_equal(alias_out$dest_ref[1], "CG")
  expect_equal(alias_out$dest_alt[1], "C")
  expect_true(is.na(alias_out$reject_reason[1]))

  star_swap_src_fa <- file.path(tmp_dir, "liftover_star_swap_src.fa")
  star_swap_dst_fa <- file.path(tmp_dir, "liftover_star_swap_dst.fa")
  star_swap_chain <- file.path(tmp_dir, "liftover_star_swap.chain")
  writeLines(c(
    ">chrS",
    "TCCAC"
  ), star_swap_src_fa)
  writeLines(c(
    ">chrD",
    "TCTGC"
  ), star_swap_dst_fa)
  writeLines(c(
    "chain 1 chrS 5 + 0 5 chrD 5 + 0 5 1",
    "5"
  ), star_swap_chain)
  expect_true(rduckhts_fasta_index(con, star_swap_src_fa, index_path = paste0(star_swap_src_fa, ".fai"))$success[1])
  expect_true(rduckhts_fasta_index(con, star_swap_dst_fa, index_path = paste0(star_swap_dst_fa, ".fai"))$success[1])

  star_swap_out <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('chrS', 3, 'C', '*,T')) AS t(chrom, pos, ref, alt)",
    chain_path = star_swap_chain,
    dst_fasta_ref = star_swap_dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = star_swap_src_fa
  )
  expect_equal(nrow(star_swap_out), 1)
  expect_true(star_swap_out$mapped[1])
  expect_equal(star_swap_out$dest_chrom[1], "chrD")
  expect_equal(star_swap_out$dest_pos[1], 3)
  expect_equal(star_swap_out$dest_ref[1], "T")
  expect_equal(star_swap_out$dest_alt[1], "*,C")
  expect_equal(star_swap_out$swap[1], 2)
  expect_true(is.na(star_swap_out$reject_reason[1]))
  expect_true(is.na(star_swap_out$note[1]))

  star_refadd_src_fa <- file.path(tmp_dir, "liftover_star_refadd_src.fa")
  star_refadd_dst_fa <- file.path(tmp_dir, "liftover_star_refadd_dst.fa")
  star_refadd_chain <- file.path(tmp_dir, "liftover_star_refadd.chain")
  writeLines(c(
    ">chrS",
    "GTGCGTGGGTGGGC"
  ), star_refadd_src_fa)
  writeLines(c(
    ">chrD",
    "GTGCGGCCGGGGGGGC"
  ), star_refadd_dst_fa)
  writeLines(c(
    "chain 1 chrS 14 + 0 14 chrD 16 + 0 16 1",
    "5 0 2",
    "9"
  ), star_refadd_chain)
  expect_true(rduckhts_fasta_index(con, star_refadd_src_fa, index_path = paste0(star_refadd_src_fa, ".fai"))$success[1])
  expect_true(rduckhts_fasta_index(con, star_refadd_dst_fa, index_path = paste0(star_refadd_dst_fa, ".fai"))$success[1])

  star_refadd_out <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('chrS', 5, 'G', '*,C')) AS t(chrom, pos, ref, alt)",
    chain_path = star_refadd_chain,
    dst_fasta_ref = star_refadd_dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = star_refadd_src_fa
  )
  expect_equal(nrow(star_refadd_out), 1)
  expect_true(star_refadd_out$mapped[1])
  expect_equal(star_refadd_out$dest_chrom[1], "chrD")
  expect_equal(star_refadd_out$dest_pos[1], 5)
  expect_equal(star_refadd_out$dest_ref[1], "GGC")
  expect_equal(star_refadd_out$dest_alt[1], "G,*,C")
  expect_equal(star_refadd_out$swap[1], -1)
  expect_true(is.na(star_refadd_out$reject_reason[1]))
  expect_true(is.na(star_refadd_out$note[1]))

  # Multi-allelic semantics should preserve all ALT alleles after liftover
  multi_out <- rduckhts_liftover(
    con,
    query = paste(
      "SELECT * FROM (VALUES",
      "('chrF', 2, 'A', 'T,G'),",
      "('chrR', 2, 'A', 'G,T')",
      ") AS t(chrom, pos, ref, alt)"
    ),
    chain_path = chain_path,
    dst_fasta_ref = dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = src_fa
  )
  expect_equal(nrow(multi_out), 2)
  row_multi_f <- multi_out[multi_out$src_chrom == "chrF", , drop = FALSE]
  row_multi_r <- multi_out[multi_out$src_chrom == "chrR", , drop = FALSE]
  expect_equal(row_multi_f$dest_ref[1], "C")
  expect_equal(row_multi_f$dest_alt[1], "A,T,G")
  expect_equal(row_multi_f$swap[1], -1)
  expect_equal(row_multi_r$dest_ref[1], "T")
  expect_equal(row_multi_r$dest_alt[1], "C,A")
  expect_equal(row_multi_r$swap[1], 0)

  # Difficult SNP beyond the source contig should reject without aborting
  rejected_edge_snp <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('chrF', 11, 'A', 'T')) AS t(chrom, pos, ref, alt)",
    chain_path = chain_path,
    dst_fasta_ref = dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = src_fa
  )
  expect_equal(nrow(rejected_edge_snp), 1)
  expect_false(rejected_edge_snp$mapped[1])
  expect_equal(rejected_edge_snp$reject_reason[1], "SourceRefMismatch")
  expect_true(is.na(rejected_edge_snp$note[1]))
  expect_true(is.na(rejected_edge_snp$dest_pos[1]))

  mid_indel <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('chrF', 2, 'CG', 'C')) AS t(chrom, pos, ref, alt)",
    chain_path = chain_path,
    dst_fasta_ref = dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = src_fa
  )
  expect_equal(nrow(mid_indel), 1)
  expect_true(mid_indel$mapped[1])
  expect_equal(mid_indel$dest_chrom[1], "chrLiftF")
  expect_equal(mid_indel$dest_pos[1], 2)
  expect_equal(mid_indel$dest_ref[1], "CG")
  expect_equal(mid_indel$dest_alt[1], "C")
  expect_true(is.na(mid_indel$note[1]))

  rejected_edge_indel <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('chrF', 10, 'AA', 'A')) AS t(chrom, pos, ref, alt)",
    chain_path = chain_path,
    dst_fasta_ref = dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = src_fa
  )
  expect_equal(nrow(rejected_edge_indel), 1)
  expect_false(rejected_edge_indel$mapped[1])
  expect_equal(rejected_edge_indel$reject_reason[1], "SourceRefMismatch")
  expect_true(is.na(rejected_edge_indel$note[1]))
  expect_true(is.na(rejected_edge_indel$dest_pos[1]))

  # Positions beyond the source contig bounds reject before unmapped-anchor handling
  unmapped <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('chrF', 20, 'A', 'T')) AS t(chrom, pos, ref, alt)",
    chain_path = chain_path,
    dst_fasta_ref = dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = src_fa
  )
  expect_equal(nrow(unmapped), 1)
  expect_false(unmapped$mapped[1])
  expect_equal(unmapped$reject_reason[1], "SourceRefMismatch")
  expect_true(is.na(unmapped$note[1]))
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
  expect_equal(unmapped_indel$reject_reason[1], "MissingContig")
  expect_true(is.na(unmapped_indel$note[1]))
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
    DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT (bcftools_liftover(",
          "NULL, 2, 'C', 'T', '%s', '%s', '%s', 1, 250, false, NULL::BIGINT, false",
          ")).src_pos"
        ),
        chain_path, dst_fa, src_fa
      )
    ),
    "chrom must be non-null"
  )

  expect_error(
    DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT (bcftools_liftover(",
          "'chrF', 0, 'C', 'T', '%s', '%s', '%s', 1, 250, false, NULL::BIGINT, false",
          ")).src_pos"
        ),
        chain_path, dst_fa, src_fa
      )
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
    {
      rduckhts_liftover(
        con,
        query = "SELECT * FROM (VALUES ('chrF', 2, 'C', 'T')) AS t(chrom, pos, ref, alt)",
        chain_path = chain_path,
        dst_fasta_ref = dst_fa,
        ref_col = "ref",
        alt_col = "alt",
        src_fasta_ref = src_fa,
        max_snp_gap = -1
      )
    },
    "max_snp_gap must be >= 0"
  )

  expect_error(
    {
      rduckhts_liftover(
        con,
        query = "SELECT * FROM (VALUES ('chrF', 2, 'C', 'T')) AS t(chrom, pos, ref, alt)",
        chain_path = chain_path,
        dst_fasta_ref = dst_fa,
        ref_col = "ref",
        alt_col = "alt",
        src_fasta_ref = src_fa,
        max_indel_inc = -1
      )
    },
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

  ## ---- L-M1: MT passthrough (lift_mt parameter) ----

  mt_src_fa <- file.path(tmp_dir, "liftover_mt_src.fa")
  mt_dst_fa <- file.path(tmp_dir, "liftover_mt_dst.fa")
  mt_chain <- file.path(tmp_dir, "liftover_mt.chain")

  writeLines(c(
    ">chrF",
    "ACGTACGTAA",
    ">MT",
    "AACCGGTTAACCGG"
  ), mt_src_fa)
  writeLines(c(
    ">chrLiftF",
    "ACGTACGTAA",
    ">chrM",
    "AACCGGTTAACCGG"
  ), mt_dst_fa)
  writeLines(c(
    "chain 100 chrF 10 + 0 10 chrLiftF 10 + 0 10 1",
    "10",
    "",
    "chain 50 MT 14 + 0 14 chrM 14 + 0 14 2",
    "14"
  ), mt_chain)

  expect_true(rduckhts_fasta_index(con, mt_src_fa, index_path = paste0(mt_src_fa, ".fai"))$success[1])
  expect_true(rduckhts_fasta_index(con, mt_dst_fa, index_path = paste0(mt_dst_fa, ".fai"))$success[1])

  # lift_mt=FALSE (default): matching MT sizes → contig rename, no chain liftover
  mt_pass <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('MT', 5, 'G', 'A')) AS t(chrom, pos, ref, alt)",
    chain_path = mt_chain,
    dst_fasta_ref = mt_dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = mt_src_fa,
    lift_mt = FALSE
  )
  expect_equal(nrow(mt_pass), 1)
  expect_true(mt_pass$mapped[1])
  expect_equal(mt_pass$dest_chrom[1], "chrM")
  expect_equal(mt_pass$dest_pos[1], 5)
  expect_equal(mt_pass$note[1], "MitochondriaPassthrough")

  # lift_mt=TRUE: force chain liftover (no passthrough)
  mt_chain_lift <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('MT', 5, 'G', 'A')) AS t(chrom, pos, ref, alt)",
    chain_path = mt_chain,
    dst_fasta_ref = mt_dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = mt_src_fa,
    lift_mt = TRUE
  )
  expect_equal(nrow(mt_chain_lift), 1)
  expect_true(mt_chain_lift$mapped[1])
  expect_equal(mt_chain_lift$dest_chrom[1], "chrM")
  expect_equal(mt_chain_lift$dest_pos[1], 5)
  # When lifted through chain, note should be OK/NULL, not MitochondriaPassthrough
  expect_true(is.na(mt_chain_lift$note[1]) || mt_chain_lift$note[1] != "MitochondriaPassthrough")

  # chrM alias also triggers passthrough
  mt_alias <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('chrM', 5, 'G', 'A')) AS t(chrom, pos, ref, alt)",
    chain_path = mt_chain,
    dst_fasta_ref = mt_dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = mt_src_fa,
    lift_mt = FALSE
  )
  expect_equal(nrow(mt_alias), 1)
  expect_true(mt_alias$mapped[1])
  expect_equal(mt_alias$dest_chrom[1], "chrM")
  expect_equal(mt_alias$note[1], "MitochondriaPassthrough")

  ## ---- L-M2: INFO/END liftover (end_pos_col parameter) ----

  # Forward chain with end_pos
  end_fwd <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('chrF', 2, 'C', 'T', 5)) AS t(chrom, pos, ref, alt, end_pos)",
    chain_path = chain_path,
    dst_fasta_ref = dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = src_fa,
    end_pos_col = "end_pos"
  )
  expect_equal(nrow(end_fwd), 1)
  expect_true(end_fwd$mapped[1])
  expect_equal(end_fwd$dest_pos[1], 2)
  expect_equal(end_fwd$dest_end[1], 5)
  expect_true("dest_end" %in% names(end_fwd))

  # Without end_pos_col: dest_end should be NA
  end_none <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('chrF', 2, 'C', 'T')) AS t(chrom, pos, ref, alt)",
    chain_path = chain_path,
    dst_fasta_ref = dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = src_fa
  )
  expect_equal(nrow(end_none), 1)
  expect_true(end_none$mapped[1])
  expect_true("dest_end" %in% names(end_none))
  expect_true(is.na(end_none$dest_end[1]))

  # MT passthrough preserves end_pos
  end_mt <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('MT', 5, 'G', 'A', 8)) AS t(chrom, pos, ref, alt, end_pos)",
    chain_path = mt_chain,
    dst_fasta_ref = mt_dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = mt_src_fa,
    lift_mt = FALSE,
    end_pos_col = "end_pos"
  )
  expect_equal(nrow(end_mt), 1)
  expect_true(end_mt$mapped[1])
  expect_equal(end_mt$dest_end[1], 8)
  expect_equal(end_mt$note[1], "MitochondriaPassthrough")

  ## ---- no_left_align parameter ----
  # For a SNP, no_left_align=TRUE produces identical results to FALSE
  # (left-alignment step only applies to indels)
  out_nla <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('chrF', 2, 'C', 'T')) AS t(chrom, pos, ref, alt)",
    chain_path = chain_path,
    dst_fasta_ref = dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = src_fa,
    no_left_align = TRUE
  )
  expect_equal(nrow(out_nla), 1)
  expect_true(out_nla$mapped[1])
  expect_equal(out_nla$dest_chrom[1], "chrLiftF")
  expect_equal(out_nla$dest_pos[1], 2)
  expect_equal(out_nla$dest_ref[1], "C")
  expect_equal(out_nla$dest_alt[1], "T")

  ## ---- repeat-run indel extension / swap regression ----
  repeat_chain <- system.file("extdata", "liftover_repeat.chain", package = "Rduckhts")
  repeat_src_fa <- system.file("extdata", "liftover_repeat_src.fa", package = "Rduckhts")
  repeat_dst_fa <- system.file("extdata", "liftover_repeat_dst.fa", package = "Rduckhts")
  expect_true(nzchar(repeat_chain))
  expect_true(nzchar(repeat_src_fa))
  expect_true(nzchar(repeat_dst_fa))

  repeat_default <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('chrS', 1, 'G', 'GT')) AS t(chrom, pos, ref, alt)",
    chain_path = repeat_chain,
    dst_fasta_ref = repeat_dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = repeat_src_fa
  )
  expect_equal(nrow(repeat_default), 1)
  expect_true(repeat_default$mapped[1])
  expect_equal(repeat_default$dest_chrom[1], "chrD")
  expect_equal(repeat_default$dest_pos[1], 1)
  expect_equal(repeat_default$dest_ref[1], "GT")
  expect_equal(repeat_default$dest_alt[1], "G")
  expect_equal(repeat_default$swap[1], 1)

  repeat_nla <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('chrS', 1, 'G', 'GT')) AS t(chrom, pos, ref, alt)",
    chain_path = repeat_chain,
    dst_fasta_ref = repeat_dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = repeat_src_fa,
    no_left_align = TRUE
  )
  expect_equal(nrow(repeat_nla), 1)
  expect_true(repeat_nla$mapped[1])
  expect_equal(repeat_nla$dest_chrom[1], "chrD")
  expect_equal(repeat_nla$dest_pos[1], 1)
  expect_equal(repeat_nla$dest_ref[1], "GTTTTTC")
  expect_equal(repeat_nla$dest_alt[1], "GTTTTC")
  expect_equal(repeat_nla$swap[1], 1)

  ## ---- clip-pad negative-score regression ----
  clip_chain <- system.file("extdata", "liftover_clip_pad.chain", package = "Rduckhts")
  clip_src_fa <- system.file("extdata", "liftover_clip_pad_src.fa", package = "Rduckhts")
  clip_dst_fa <- system.file("extdata", "liftover_clip_pad_dst.fa", package = "Rduckhts")
  expect_true(nzchar(clip_chain))
  expect_true(nzchar(clip_src_fa))
  expect_true(nzchar(clip_dst_fa))

  clip_default <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('chrS', 7, 'TC', 'T')) AS t(chrom, pos, ref, alt)",
    chain_path = clip_chain,
    dst_fasta_ref = clip_dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = clip_src_fa
  )
  expect_equal(nrow(clip_default), 1)
  expect_true(clip_default$mapped[1])
  expect_equal(clip_default$dest_chrom[1], "chrD")
  expect_equal(clip_default$dest_pos[1], 7)
  expect_equal(clip_default$dest_ref[1], "T")
  expect_equal(clip_default$dest_alt[1], "TC")
  expect_equal(clip_default$swap[1], 1)
  expect_equal(clip_default$note[1], "Padded")

  clip_nla <- rduckhts_liftover(
    con,
    query = "SELECT * FROM (VALUES ('chrS', 7, 'TC', 'T')) AS t(chrom, pos, ref, alt)",
    chain_path = clip_chain,
    dst_fasta_ref = clip_dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = clip_src_fa,
    no_left_align = TRUE
  )
  expect_equal(nrow(clip_nla), 1)
  expect_true(clip_nla$mapped[1])
  expect_equal(clip_nla$dest_chrom[1], "chrD")
  expect_equal(clip_nla$dest_pos[1], 7)
  expect_equal(clip_nla$dest_ref[1], "TT")
  expect_equal(clip_nla$dest_alt[1], "TCT")
  expect_equal(clip_nla$swap[1], 1)
  expect_equal(clip_nla$note[1], "Padded")

  ## ---- mixed context cache reuse / eviction ----
  chain_copies <- character(12)
  src_copies <- character(12)
  dst_copies <- character(12)
  for (i in seq_len(12)) {
    chain_copies[i] <- file.path(tmp_dir, sprintf("liftover_%02d.chain", i))
    src_copies[i] <- file.path(tmp_dir, sprintf("liftover_src_%02d.fa", i))
    dst_copies[i] <- file.path(tmp_dir, sprintf("liftover_dst_%02d.fa", i))
    file.copy(chain_path, chain_copies[i], overwrite = TRUE)
    file.copy(src_fa, src_copies[i], overwrite = TRUE)
    file.copy(dst_fa, dst_copies[i], overwrite = TRUE)
    expect_true(rduckhts_fasta_index(con, src_copies[i], index_path = paste0(src_copies[i], ".fai"))$success[1])
    expect_true(rduckhts_fasta_index(con, dst_copies[i], index_path = paste0(dst_copies[i], ".fai"))$success[1])
  }

  for (i in c(seq_len(12), 1, 6, 12)) {
    res <- DBI::dbGetQuery(
      con,
      sprintf(
        paste(
          "SELECT lo.dest_chrom, lo.dest_pos, lo.dest_ref, lo.dest_alt, lo.mapped",
          "FROM (SELECT bcftools_liftover(",
          "'chrF', 2, 'C', 'T', '%s', '%s', '%s', 1, 250, false, NULL::BIGINT, false",
          ") AS lo) s"
        ),
        chain_copies[i], dst_copies[i], src_copies[i]
      )
    )
    expect_equal(res$dest_chrom[1], "chrLiftF")
    expect_equal(res$dest_pos[1], 2)
    expect_equal(res$dest_ref[1], "C")
    expect_equal(res$dest_alt[1], "T")
    expect_true(res$mapped[1])
  }
}

test_liftover()
