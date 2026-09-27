library(tinytest)
library(DBI)

test_ancestry_relations <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  dbExecute(con, paste(
    "CREATE TABLE ancestry_ref AS SELECT 'chr1' AS chromosome, (100+i)::INTEGER AS position,",
    "CASE WHEN i % 2 = 0 THEN 'A' ELSE 'C' END AS allele_a,",
    "CASE WHEN i % 2 = 0 THEN 'G' ELSE 'T' END AS allele_b,",
    "g.group_id, CASE WHEN g.group_id = 'A' THEN 0.15 + i*0.06",
    "ELSE 0.8 - i*0.04 END AS frequency",
    "FROM range(8) t(i), (VALUES ('A'), ('B')) g(group_id)"
  ))
  dbExecute(con, paste(
    "CREATE TABLE ancestry_pc AS SELECT DISTINCT chromosome, position, allele_a, allele_b,",
    "pc, CASE WHEN pc = 1 THEN position - 103.0 ELSE sin(position) END AS loading",
    "FROM ancestry_ref, (VALUES (1), (2)) p(pc)"
  ))
  dbExecute(con, "CREATE TABLE ancestry_correction AS SELECT 1 AS pc, 1.0 AS coefficient UNION ALL SELECT 2, 1.0")
  panel <- rduckhts_ancestry_panel(con, "ancestry_ref", "ancestry_pc",
                                   "ancestry_selected", "GRCh38", spacing_bp = 1,
                                   max_sites = 5)
  expect_equal(panel$sites, 5)
  expect_equal(nchar(panel$panel_sha256), 64L)
  selected <- dbGetQuery(con, "SELECT site_index, position FROM ancestry_selected ORDER BY site_index")
  expect_equal(selected$site_index, 0:4)
  expect_equal(selected$position, 100:104)
  dbExecute(con, "CREATE TABLE ancestry_candidates AS SELECT region, position, allele_a, allele_b FROM ancestry_selected")
  matched_panel <- rduckhts_ancestry_panel(con, "ancestry_ref", "ancestry_pc",
                                           "ancestry_matched", "GRCh38",
                                           candidate_table = "ancestry_candidates")
  expect_equal(matched_panel$panel_sha256, panel$panel_sha256)
  dbExecute(con, paste(
    "CREATE TABLE ancestry_input AS SELECT 'S' AS sample_id, chromosome, position, allele_a, allele_b,",
    "sum(CASE WHEN group_id = 'A' THEN 0.7 ELSE 0.3 END * frequency) AS frequency",
    "FROM ancestry_ref GROUP BY chromosome, position, allele_a, allele_b"
  ))
  dbExecute(con, "UPDATE ancestry_input SET allele_a = 'T', allele_b = 'C' WHERE position = 100")
  dbExecute(con, paste("UPDATE ancestry_input SET allele_a = 'T', allele_b = 'C',",
                       "frequency = 1 - frequency WHERE position = 101"))
  dbExecute(con, "UPDATE ancestry_input SET frequency = frequency + 0.04 WHERE position = 102")
  dbExecute(con, "UPDATE ancestry_input SET frequency = NULL WHERE position = 107")
  dbExecute(con, "INSERT INTO ancestry_input SELECT * FROM ancestry_input WHERE position = 106")
  dbExecute(con, "INSERT INTO ancestry_input VALUES ('S', 'chr1', 999, 'A', 'G', 0.5)")
  out <- rduckhts_ancestry_proportions(con, "ancestry_input", "ancestry_ref",
                                       "ancestry_pc", "ancestry_correction", min_cor = 0)
  expect_equal(out$status, rep("ok", 2L))
  expect_equal(out$used_variants, rep(6, 2L))
  expect_equal(out$input_variants, rep(10, 2L))
  expect_equal(out$reversed_variants, rep(1, 2L))
  expect_equal(out$flipped_variants, rep(1, 2L))
  expect_equal(out$duplicate_variants, rep(2, 2L))
  expect_equal(out$unmatched_variants, rep(1, 2L))
  expect_equal(out$missing_variants, rep(1, 2L))
  ref <- dbGetQuery(con, "SELECT position, group_id, frequency FROM ancestry_ref WHERE position BETWEEN 100 AND 105 ORDER BY position, group_id")
  f <- matrix(ref$frequency, ncol = 2L, byrow = TRUE)[1:6, ]
  input <- dbGetQuery(con, "SELECT position, frequency FROM ancestry_input WHERE position BETWEEN 100 AND 105 ORDER BY position")
  input$frequency[input$position == 101L] <- 1 - input$frequency[input$position == 101L]
  p <- cbind((100:105) - 103, sin(100:105))
  x <- crossprod(p, f)
  y <- drop(crossprod(p, input$frequency))
  delta <- x[, 1L] - x[, 2L]
  expected_a <- max(0, min(1, sum(delta * (y - x[, 2L])) / sum(delta^2)))
  expect_true(max(abs(out$proportion - c(expected_a, 1 - expected_a))) < 1e-6)
  gated <- rduckhts_ancestry_proportions(con, "ancestry_input", "ancestry_ref",
                                         "ancestry_pc", "ancestry_correction", min_cor = 1)
  expect_equal(gated$status, rep("low_correlation", 2L))
  expect_true(all(is.na(gated$proportion)))
  dbExecute(con, paste(
    "CREATE TABLE ancestry_scaled AS SELECT 'S' AS sample_id, chromosome, position,",
    "allele_a, allele_b, sum(CASE WHEN group_id = 'A' THEN 0.35 ELSE 0.15 END * frequency) AS frequency",
    "FROM ancestry_ref GROUP BY chromosome, position, allele_a, allele_b"
  ))
  scaled <- rduckhts_ancestry_proportions(con, "ancestry_scaled", "ancestry_ref",
                                          "ancestry_pc", "ancestry_correction",
                                          sum_to_one = FALSE, min_cor = 0)
  expect_true(all(scaled$status == "ok"))
  expect_true(abs(sum(scaled$proportion) - 0.5) < 1e-5)
  dbExecute(con, paste(
    "CREATE TABLE ancestry_geno_calls AS SELECT chromosome AS CHROM, position AS POS,",
    "allele_a AS REF, [allele_b] AS ALT,",
    "[struct_pack(sample_index := 0, alleles := CASE WHEN position = 107",
    "THEN [NULL::INTEGER, 1] ELSE [0, 1] END)] AS calls",
    "FROM ancestry_ref GROUP BY chromosome, position, allele_a, allele_b"
  ))
  dbExecute(con, "CREATE TABLE ancestry_samples AS SELECT 0 AS sample_index, 'S' AS sample_name")
  genotype <- rduckhts_ancestry_geno(con, "ancestry_geno_calls", "ancestry_ref",
                                      "ancestry_pc", "ancestry_correction",
                                      samples_table = "ancestry_samples", min_cor = 0)
  expect_equal(genotype$sample_id, rep("S", 2L))
  expect_equal(genotype$used_variants, rep(7, 2L))
  expect_equal(genotype$missing_variants, rep(1, 2L))
  expect_error(rduckhts_ancestry_proportions(con, "ancestry_input", "ancestry_ref",
                                              "ancestry_pc", "ancestry_correction", min_cor = 2))
}

test_ancestry_relations()

test_ancestry_bam_cram <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  paths <- function(x) system.file("extdata", x, package = "Rduckhts")
  dbExecute(con, paste(
    "CREATE TABLE ancestry_bam_panel AS SELECT * FROM (VALUES",
    "('WBcel235', 0::UBIGINT, 'CHROMOSOME_I', 1::UBIGINT, 'A', 'G'),",
    "('WBcel235', 1::UBIGINT, 'CHROMOSOME_I', 914::UBIGINT, 'A', 'C'),",
    "('WBcel235', 2::UBIGINT, 'CHROMOSOME_I', 2::UBIGINT, 'A', 'G'))",
    "p(assembly, site_index, region, position, allele_a, allele_b)"
  ))
  dbExecute(con, paste(
    "CREATE TABLE ancestry_bam_ref AS SELECT region AS chromosome, position,",
    "allele_a, allele_b, group_id, CASE WHEN group_id='A' THEN 0.2 ELSE 0.8 END AS frequency",
    "FROM ancestry_bam_panel, (VALUES ('A'), ('B')) g(group_id)"
  ))
  dbExecute(con, paste(
    "CREATE TABLE ancestry_bam_pc AS SELECT region AS chromosome, position,",
    "allele_a, allele_b, pc, CASE WHEN pc=1 THEN 1.0 ELSE position/1000.0 END AS loading",
    "FROM ancestry_bam_panel, (VALUES (1), (2)) p(pc)"
  ))
  dbExecute(con, "CREATE TABLE ancestry_bam_correction AS SELECT 1 AS pc, 1.0 AS coefficient UNION ALL SELECT 2, 1.0")
  bam <- rduckhts_ancestry_bam(
    con, paths("range.bam"), "S", paths("ce.fa"), "ancestry_bam_panel",
    "ancestry_bam_ref", "ancestry_bam_pc", "ancestry_bam_correction",
    min_depth = 1, min_cor = 0, index_path = paths("range.bam.bai"),
    reference_index_path = paths("ce.fa.fai")
  )
  cram <- rduckhts_ancestry_bam(
    con, paths("range.cram"), "S", paths("ce.fa"), "ancestry_bam_panel",
    "ancestry_bam_ref", "ancestry_bam_pc", "ancestry_bam_correction",
    frequency_method = "called_genotype", min_depth = 1, min_cor = 0,
    index_path = paths("range.cram.crai"), reference_index_path = paths("ce.fa.fai")
  )
  expect_equal(bam$frequency_method, rep("allele_fraction", 2L))
  expect_equal(cram$frequency_method, rep("called_genotype", 2L))
  expect_equal(bam$used_variants, rep(1, 2L))
  expect_equal(cram$used_variants, rep(1, 2L))
  expect_equal(bam$status, rep("low_correlation", 2L))
  expect_true(all(is.na(cram$proportion)))
}

test_ancestry_bam_cram()
