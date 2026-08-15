library(tinytest)
library(DBI)

test_liftover_stress <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  expect_silent(dbExecute(con, "PRAGMA threads=4"))

  tmp_dir <- tempfile("duckhts_liftover_stress_")
  dir.create(tmp_dir)

  src_fa <- file.path(tmp_dir, "liftover_src.fa")
  dst_fa <- file.path(tmp_dir, "liftover_dst.fa")
  chain_path <- file.path(tmp_dir, "liftover.chain")

  writeLines(
    c(
      ">chrF",
      "ACGTACGTAA",
      ">chrR",
      "AACCGGTTAA"
    ),
    src_fa
  )
  writeLines(
    c(
      ">chrLiftF",
      "ACGTACGTAA",
      ">chrLiftR",
      "TTAACCGGTT"
    ),
    dst_fa
  )
  writeLines(
    c(
      "chain 100 chrF 10 + 0 10 chrLiftF 10 + 0 10 1",
      "10",
      "",
      "chain 100 chrR 10 + 0 10 chrLiftR 10 - 0 10 2",
      "10"
    ),
    chain_path
  )

  expect_true(rduckhts_fasta_index(
    con,
    src_fa,
    index_path = paste0(src_fa, ".fai")
  )$success[1])
  expect_true(rduckhts_fasta_index(
    con,
    dst_fa,
    index_path = paste0(dst_fa, ".fai")
  )$success[1])

  expect_silent(dbExecute(
    con,
    paste(
      "CREATE TEMP TABLE stress_input AS",
      "SELECT CASE WHEN i % 2 = 0 THEN 'chrF' ELSE 'chrR' END AS chrom,",
      "2::BIGINT AS pos,",
      "CASE WHEN i % 2 = 0 THEN 'C' ELSE 'A' END AS ref,",
      "CASE WHEN i % 2 = 0 THEN 'T' ELSE 'G' END AS alt",
      "FROM range(1000000) t(i)"
    )
  ))

  out <- rduckhts_liftover(
    con,
    query = "stress_input",
    chain_path = chain_path,
    dst_fasta_ref = dst_fa,
    ref_col = "ref",
    alt_col = "alt",
    src_fasta_ref = src_fa
  )

  expect_equal(nrow(out), 1000000)
  expect_true(all(out$mapped))
  expect_equal(sum(out$reverse_complemented), 500000)
  expect_equal(sum(out$dest_chrom == "chrLiftF"), 500000)
  expect_equal(sum(out$dest_chrom == "chrLiftR"), 500000)

  nw_chain <- system.file("extdata", "liftover_nw_limit.chain", package = "Rduckhts")
  nw_src <- system.file("extdata", "liftover_nw_limit_src.fa", package = "Rduckhts")
  nw_dst <- system.file("extdata", "liftover_nw_limit_dst.fa", package = "Rduckhts")
  expect_true(nzchar(nw_chain))
  expect_true(nzchar(nw_src))
  expect_true(nzchar(nw_dst))
  expect_error(
    DBI::dbGetQuery(
      con,
      paste0(
        "SELECT bcftools_liftover('chrS', 2, substr(repeat('ACGT', 551), 2, 2200), 'C', ",
        DBI::dbQuoteString(con, nw_chain), ", ",
        DBI::dbQuoteString(con, nw_dst), ", ",
        DBI::dbQuoteString(con, nw_src),
        ", 1, 250, false, NULL::BIGINT, true)"
      )
    ),
    pattern = "clip-pad alignment exceeds the 4194304-cell limit",
    fixed = TRUE
  )
  expect_equal(DBI::dbGetQuery(con, "SELECT 42 AS answer")$answer, 42)
}

test_liftover_stress()
