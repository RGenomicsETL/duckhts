library(tinytest)
library(DBI)

test_somalier_relatedness_wrappers <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  qid <- function(x) as.character(dbQuoteIdentifier(con, x))
  qstr <- function(x) as.character(dbQuoteString(con, x))

  panel_name <- "somalier panel's relation"
  evidence_name <- 'somalier evidence; relation'
  panel <- qid(panel_name)
  evidence <- qid(evidence_name)
  dbExecute(con, sprintf(paste(
    "CREATE TABLE %s AS SELECT 'GRCh38'::VARCHAR AS assembly,",
    "i::UBIGINT AS site_index, 'chr1'::VARCHAR AS region,",
    "(100 + i)::UBIGINT AS position,",
    "CASE WHEN i %% 2 = 0 THEN 'A' ELSE 'C' END::VARCHAR AS allele_a,",
    "CASE WHEN i %% 2 = 0 THEN 'G' ELSE 'T' END::VARCHAR AS allele_b",
    "FROM range(8) sites(i)"), panel))
  dbExecute(con, sprintf(paste(
    "CREATE TABLE %s AS SELECT c.sample_id::VARCHAR AS sample_id,",
    "p.assembly, p.site_index, p.region, p.position, p.allele_a, p.allele_b,",
    "c.a::UBIGINT AS a, c.b::UBIGINT AS b, c.other::UBIGINT AS other",
    "FROM (VALUES",
    "('A', 0, 20, 0, 0), ('A', 1, 20, 0, 0),",
    "('A', 2, 10, 10, 0), ('A', 3, 10, 10, 0),",
    "('A', 4, 0, 20, 0), ('A', 5, NULL, NULL, NULL),",
    "('A', 6, 0, 20, 0), ('A', 7, 0, 0, 0),",
    "('B''s sample', 0, 20, 0, 0), ('B''s sample', 1, 0, 20, 0),",
    "('B''s sample', 2, 10, 10, 0), ('B''s sample', 3, 0, 20, 0),",
    "('B''s sample', 4, 0, 20, 0), ('B''s sample', 5, 20, 0, 0),",
    "('B''s sample', 6, 20, 0, 0), ('B''s sample', 7, NULL, NULL, NULL),",
    "('Z', 0, 0, 0, 0), ('Z', 1, 0, 0, 0),",
    "('Z', 2, 0, 0, 0), ('Z', 3, 0, 0, 0),",
    "('Z', 4, 0, 0, 0), ('Z', 5, 0, 0, 0),",
    "('Z', 6, 0, 0, 0), ('Z', 7, 0, 0, 0),",
    "('U', 0, NULL, NULL, NULL), ('U', 1, NULL, NULL, NULL),",
    "('U', 2, NULL, NULL, NULL), ('U', 3, NULL, NULL, NULL),",
    "('U', 4, NULL, NULL, NULL), ('U', 5, NULL, NULL, NULL),",
    "('U', 6, NULL, NULL, NULL), ('U', 7, NULL, NULL, NULL)",
    ") c(sample_id, site_index, a, b, other)",
    "JOIN %s p USING (site_index)"), evidence, panel))

  selected <- rduckhts_somalier_sketches(
    con, evidence_table = evidence_name, panel_table = panel_name,
    sample_ids = c("A", "B's sample"), max_sites = 8
  )
  expect_equal(nrow(selected), 2L)

  panel_region_limit <- "somalier panel region limit"
  evidence_region_limit <- "somalier evidence region limit"
  panel_region_over_limit <- "somalier panel region over limit"
  evidence_region_over_limit <- "somalier evidence region over limit"
  dbExecute(con, sprintf(
    "CREATE TEMP VIEW %s AS SELECT * REPLACE(repeat('é', 512) AS region) FROM %s",
    qid(panel_region_limit), panel
  ))
  dbExecute(con, sprintf(
    "CREATE TEMP VIEW %s AS SELECT * REPLACE(repeat('é', 512) AS region) FROM %s",
    qid(evidence_region_limit), evidence
  ))
  expect_equal(nrow(rduckhts_somalier_sketches(
    con, evidence_table = evidence_region_limit,
    panel_table = panel_region_limit, max_sites = 8
  )), 4L)
  dbExecute(con, sprintf(paste0(
    "CREATE TEMP VIEW %s AS SELECT * REPLACE(",
    "repeat('é', 512) || 'a' AS region) FROM %s"),
    qid(panel_region_over_limit), panel
  ))
  dbExecute(con, sprintf(paste0(
    "CREATE TEMP VIEW %s AS SELECT * REPLACE(",
    "repeat('é', 512) || 'a' AS region) FROM %s"),
    qid(evidence_region_over_limit), evidence
  ))
  expect_error(
    rduckhts_somalier_sketches(
      con, evidence_table = evidence_region_over_limit,
      panel_table = panel_region_over_limit, max_sites = 8
    ),
    "panel assembly and region must be at most 1024 bytes"
  )
  expect_error(
    rduckhts_somalier_sketches(
      con, evidence_table = evidence_name, panel_table = panel_name,
      sample_ids = c("A", "absent"), max_sites = 8
    ),
    "every sample_ids value"
  )

  sketches_name <- 'somalier sketches "verified"'
  expect_true(rduckhts_somalier_sketches(
    con, evidence_table = evidence_name, panel_table = panel_name,
    table_name = sketches_name, max_sites = 8
  ))
  expect_equal(
    dbGetQuery(con, sprintf(
      "SELECT count(*) AS n, count(DISTINCT sketch.panel_sha256) AS panels FROM %s",
      qid(sketches_name)
    )),
    data.frame(n = 4, panels = 1)
  )

  all_pairs <- rduckhts_somalier_relatedness(
    con, sketches_table = sketches_name, max_sites = 8
  )
  expect_equal(nrow(all_pairs), 6L)
  ab <- all_pairs[all_pairs$sample_a == "A" & all_pairs$sample_b == "B's sample", ]
  expect_equal(
    unname(unlist(ab[c("jointly_called", "ibs0", "ibs2", "shared_hets",
                            "het_ab", "shared_hom_b")])),
    c(6, 2, 3, 1, 3, 1)
  )
  expect_equal(ab$relatedness, -2)

  verified <- dbGetQuery(con, sprintf(paste(
    "WITH comparison AS (SELECT a.sketch AS left_sketch,",
    "b.sketch AS right_sketch,",
    "duckhts_somalier_relatedness(a.sketch, b.sketch, 8) AS pair",
    "FROM %s a, %s b WHERE a.sketch.sample_id = 'A'",
    "AND b.sketch.sample_id = 'B''s sample')",
    "SELECT duckhts_somalier_verify_relatedness(pair, left_sketch,",
    "right_sketch, 8) AS valid,",
    "duckhts_somalier_verify_relatedness(",
    "struct_update(pair, ibs0 := pair.ibs0 + 1::UBIGINT),",
    "left_sketch, right_sketch, 8) AS altered FROM comparison"
  ), qid(sketches_name), qid(sketches_name)))
  expect_equal(verified, data.frame(valid = TRUE, altered = FALSE))
  receipt_sketches <- "somalier_receipt_sketches"
  dbExecute(con, sprintf("CREATE TEMP VIEW %s AS SELECT * FROM %s",
    qid(receipt_sketches), qid(sketches_name)))
  sketch_receipt <- function(evidence_source) {
    dbGetQuery(con, sprintf(
      "SELECT duckhts_somalier_verify_sketches(%s, %s, %s, 8) AS valid",
      qstr(evidence_source), qstr(panel_name), qstr(receipt_sketches)
    ))$valid
  }
  expect_true(sketch_receipt(evidence_name))
  negative_zero_sketches <- "somalier negative zero sketches"
  dbExecute(con, sprintf(paste0(
    "CREATE TEMP TABLE %s AS SELECT * FROM ",
    "duckhts_somalier_prepare_sketches(%s, %s, 7, 0.3, -0.0::DOUBLE, ",
    "max_sites := 8)"), qid(negative_zero_sketches), qstr(evidence_name),
    qstr(panel_name)))
  negative_zero_check <- dbGetQuery(con, sprintf(paste0(
    "SELECT duckhts_somalier_verify_sketches(%s, %s, %s, 8) AS valid, ",
    "count(*) FILTER (WHERE signbit(sketch.hom_balance_cutoff)) AS negative ",
    "FROM %s"), qstr(evidence_name), qstr(panel_name),
    qstr(negative_zero_sketches), qid(negative_zero_sketches)))
  expect_equal(negative_zero_check, data.frame(valid = TRUE, negative = 0))
  expect_equal(nrow(rduckhts_somalier_relatedness(
    con, sketches_table = negative_zero_sketches, max_sites = 8
  )), 6L)
  invalid_settings <- "somalier invalid sketch settings"
  dbExecute(con, sprintf(paste(
    "CREATE TEMP VIEW %s AS SELECT struct_update(",
    "sketch, min_depth := 0::UBIGINT) AS sketch FROM %s"
  ), qid(invalid_settings), qid(receipt_sketches)))
  invalid_receipt <- dbGetQuery(con, sprintf(
    "SELECT duckhts_somalier_verify_sketches(%s, %s, %s, 8) AS valid",
    qstr(evidence_name), qstr(panel_name), qstr(invalid_settings)
  ))$valid
  expect_false(invalid_receipt)
  changed_evidence <- "somalier changed raw count"
  dbExecute(con, sprintf(paste(
    "CREATE TEMP VIEW %s AS SELECT * REPLACE(",
    "CASE WHEN sample_id = 'A' AND site_index = 0 THEN 19::UBIGINT",
    "ELSE a END AS a) FROM %s"), qid(changed_evidence), evidence))
  expect_false(sketch_receipt(changed_evidence))

  pairs_name <- "somalier selected pairs' relation"
  pairs <- qid(pairs_name)
  dbExecute(con, sprintf(
    "CREATE TABLE %s(sample_a VARCHAR, sample_b VARCHAR)", pairs
  ))
  dbExecute(con, sprintf(
    "INSERT INTO %s VALUES ('B''s sample', 'A'), ('Z', 'U')", pairs
  ))
  selected_pairs <- rduckhts_somalier_relatedness(
    con, sketches_table = sketches_name, pairs_table = pairs_name, max_sites = 8
  )
  selected_pairs <- selected_pairs[order(selected_pairs$sample_a,
                                         selected_pairs$sample_b), ]
  rownames(selected_pairs) <- NULL
  expect_equal(nrow(selected_pairs), 2L)
  expect_equal(selected_pairs$sample_a, c("B's sample", "Z"))
  expect_equal(selected_pairs$sample_b, c("A", "U"))
  no_evidence <- selected_pairs[selected_pairs$sample_a == "Z", ]
  expect_equal(no_evidence$status, "no_evidence")
  expect_equal(no_evidence$jointly_called, 0)
  expect_true(is.na(no_evidence$relatedness))

  # Tables and Parquet paths containing SQL punctuation remain data, and packed
  # sketches round-trip through ordinary Parquet without a private serializer.
  guard <- "somalier_injection_guard"
  dbExecute(con, sprintf("CREATE TABLE %s(value INTEGER)", qid(guard)))
  panel_path <- tempfile("somalier panel's ", fileext = ".parquet")
  evidence_path <- tempfile("somalier evidence'); DROP TABLE guard; -- ",
                            fileext = ".parquet")
  sketches_path <- tempfile("somalier sketches' ", fileext = ".parquet")
  on.exit(unlink(c(panel_path, evidence_path, sketches_path)), add = TRUE)
  dbExecute(con, sprintf("COPY %s TO %s (FORMAT PARQUET)", panel, qstr(panel_path)))
  dbExecute(con, sprintf("COPY %s TO %s (FORMAT PARQUET)", evidence, qstr(evidence_path)))
  parquet_sketches <- "somalier parquet sketches; relation"
  expect_true(rduckhts_somalier_sketches(
    con, evidence_parquet = evidence_path, panel_parquet = panel_path,
    table_name = parquet_sketches, max_sites = 8
  ))
  dbExecute(con, sprintf("COPY %s TO %s (FORMAT PARQUET)",
                         qid(parquet_sketches), qstr(sketches_path)))
  parquet_pairs <- rduckhts_somalier_relatedness(
    con, sketches_parquet = sketches_path, pairs_table = pairs_name, max_sites = 8
  )
  parquet_pairs <- parquet_pairs[order(parquet_pairs$sample_a,
                                       parquet_pairs$sample_b), ]
  rownames(parquet_pairs) <- NULL
  expect_equal(
    parquet_pairs[c("sample_a", "sample_b", "status", "jointly_called",
                    "ibs0", "ibs2", "relatedness")],
    selected_pairs[c("sample_a", "sample_b", "status", "jointly_called",
                     "ibs0", "ibs2", "relatedness")]
  )
  expect_true(dbExistsTable(con, guard))

  injected_output <- 'result"; DROP TABLE somalier_injection_guard; --'
  expect_true(rduckhts_somalier_relatedness(
    con, sketches_table = sketches_name, pairs_table = pairs_name,
    table_name = injected_output, max_sites = 8
  ))
  expect_true(dbExistsTable(con, injected_output))
  expect_true(dbExistsTable(con, guard))

  # A failed replacement leaves the prior table intact.
  preserve <- "somalier_atomic_output"
  dbExecute(con, sprintf("CREATE TABLE %s AS SELECT 42 AS sentinel", qid(preserve)))
  bad_evidence <- "somalier_bad_orientation"
  dbExecute(con, sprintf(paste(
    "CREATE TABLE %s AS SELECT * REPLACE(",
    "CASE WHEN site_index = 0 THEN 'T' ELSE allele_a END AS allele_a)",
    "FROM %s WHERE sample_id = 'A'"), qid(bad_evidence), evidence))
  expect_error(
    rduckhts_somalier_sketches(
      con, evidence_table = bad_evidence, panel_table = panel_name,
      table_name = preserve, max_sites = 8, overwrite = TRUE
    ),
    "does not match ordered panel"
  )
  expect_equal(dbGetQuery(con, sprintf("SELECT * FROM %s", qid(preserve))),
               data.frame(sentinel = 42L))

  # Selected-pair and sketch identity checks fail rather than dropping or
  # multiplying rows.
  dbExecute(con, sprintf("INSERT INTO %s VALUES ('A', 'missing')", pairs))
  expect_error(
    rduckhts_somalier_relatedness(
      con, sketches_table = sketches_name, pairs_table = pairs_name, max_sites = 8
    ),
    "different existing sample IDs"
  )
  dbExecute(con, sprintf("DELETE FROM %s WHERE sample_b = 'missing'", pairs))
  dbExecute(con, sprintf("INSERT INTO %s VALUES ('Z', 'U')", pairs))
  expect_error(
    rduckhts_somalier_relatedness(
      con, sketches_table = sketches_name, pairs_table = pairs_name, max_sites = 8
    ),
    "different existing sample IDs"
  )
  duplicate_sketches <- "somalier_duplicate_sketches"
  dbExecute(con, sprintf(
    paste0("CREATE TABLE %s AS SELECT * FROM %s UNION ALL ",
           "SELECT * FROM %s WHERE sketch.sample_id = 'A'"),
    qid(duplicate_sketches), qid(sketches_name), qid(sketches_name)
  ))
  expect_error(
    rduckhts_somalier_relatedness(
      con, sketches_table = duplicate_sketches, max_sites = 8
    ),
    "one non-NULL sketch"
  )

  # The wrapper and native preparation both enforce the explicit capacity.
  expect_error(
    rduckhts_somalier_sketches(
      con, evidence_table = evidence_name, panel_table = panel_name, max_sites = 7
    ),
    "site_count"
  )
  expect_error(
    rduckhts_somalier_relatedness(
      con, sketches_table = sketches_name, max_sites = 7
    ),
    "site_count exceeds max_sites"
  )
  oversized_evidence <- "somalier_oversized_sample_identity"
  dbExecute(con, sprintf(paste(
    "CREATE TEMP VIEW %s AS SELECT * REPLACE(",
    "repeat('s', 1025)::VARCHAR AS sample_id) FROM %s",
    "WHERE sample_id = 'A'"
  ), qid(oversized_evidence), evidence))
  expect_error(
    rduckhts_somalier_sketches(
      con, evidence_table = oversized_evidence, panel_table = panel_name,
      max_sites = 8
    ),
    "sample_id and assembly must be at most 1024 bytes"
  )
  oversized_sketches <- "somalier_oversized_persisted_identity"
  dbExecute(con, sprintf(paste(
    "CREATE TEMP VIEW %s AS SELECT struct_update(sketch,",
    "sample_id := CASE WHEN sketch.sample_id = 'A' THEN repeat('s', 1025)",
    "ELSE sketch.sample_id END) AS sketch FROM %s"
  ), qid(oversized_sketches), qid(sketches_name)))
  expect_error(
    rduckhts_somalier_relatedness(
      con, sketches_table = oversized_sketches, max_sites = 8
    ),
    "persisted sample_id and assembly must be at most 1024 bytes"
  )

  # Word-count transitions exercise the public preparation wrapper, not a
  # hand-constructed mask.
  for (site_count in 63:65) {
    dbExecute(con, paste(
      "CREATE OR REPLACE TEMP TABLE somalier_width_panel AS",
      "SELECT 'GRCh38' AS assembly, i::UBIGINT AS site_index, 'chrW' AS region,",
      "(i + 1)::UBIGINT AS position, 'A' AS allele_a, 'C' AS allele_b",
      sprintf("FROM range(%d) sites(i)", site_count)
    ))
    dbExecute(con, paste(
      "CREATE OR REPLACE TEMP TABLE somalier_width_evidence AS",
      "SELECT 'width' AS sample_id, *, 20::UBIGINT AS a, 0::UBIGINT AS b,",
      "0::UBIGINT AS other FROM somalier_width_panel"
    ))
    rduckhts_somalier_sketches(
      con, evidence_table = "somalier_width_evidence",
      panel_table = "somalier_width_panel", table_name = "somalier_width_sketch",
      max_sites = 65, overwrite = TRUE
    )
    expect_equal(
      dbGetQuery(con, "SELECT len(sketch.hom_a) AS words FROM somalier_width_sketch")$words,
      ceiling(site_count / 64)
    )
  }

  expect_error(
    rduckhts_somalier_sketches(con, evidence_table = evidence_name),
    "panel: supply exactly one"
  )
  expect_error(
    rduckhts_somalier_sketches(
      con, evidence_table = evidence_name, evidence_parquet = evidence_path,
      panel_table = panel_name
    ),
    "evidence: supply exactly one"
  )
  expect_error(
    rduckhts_somalier_sketches(
      con, evidence_table = evidence_name, panel_table = panel_name,
      sample_ids = c("A", "A")
    ),
    "distinct"
  )
  expect_error(
    rduckhts_somalier_sketches(
      con, evidence_table = evidence_name, panel_table = panel_name, max_sites = 0
    ),
    "positive"
  )
  expect_error(
    rduckhts_somalier_sketches(
      con, evidence_table = evidence_name, panel_table = panel_name,
      min_het_balance = 0.2, hom_balance_cutoff = 0.3
    ),
    "hom_balance_cutoff"
  )
  expect_error(
    rduckhts_somalier_relatedness(
      con, sketches_table = sketches_name, max_sites = 100000001
    ),
    "100000000"
  )
}

test_somalier_contamination_wrappers <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  qid <- function(x) as.character(dbQuoteIdentifier(con, x))
  qstr <- function(x) as.character(dbQuoteString(con, x))

  # The aggregate and public macro report the same count constraints.
  expect_error(dbGetQuery(con, paste(
    "SELECT __duckhts_somalier_charr(",
    "'Overflow', 'GRCh38', repeat('a', 64), repeat('b', 64),",
    "0::UBIGINT, 1::UBIGINT, 18446744073709551615::UBIGINT, 1::UBIGINT, 0::UBIGINT,",
    "200::UBIGINT, 11::UBIGINT, 0.25::DOUBLE, 7::UBIGINT,",
    "1000000::UBIGINT, 0.1::DOUBLE, 0.001::DOUBLE,",
    "16000000::UBIGINT, 1::UBIGINT)"
  )), "measured counts must fit UINTEGER")
  expect_error(dbGetQuery(con, paste(
    "SELECT __duckhts_somalier_charr(",
    "'Partial', 'GRCh38', repeat('a', 64), repeat('b', 64),",
    "0::UBIGINT, 1::UBIGINT, 199::UBIGINT, NULL::UBIGINT, 0::UBIGINT,",
    "200::UBIGINT, 11::UBIGINT, 0.25::DOUBLE, 7::UBIGINT,",
    "1000000::UBIGINT, 0.1::DOUBLE, 0.001::DOUBLE,",
    "16000000::UBIGINT, 1::UBIGINT)"
  )), "counts must be all measured or all unavailable")

  panel_name <- "contamination panel's relation"
  evidence_name <- 'contamination evidence; relation'
  frequency_name <- 'population frequency "B"'
  match_frequency_name <- "matched frequency's relation"
  pairs_name <- "ordered contamination pairs' relation"
  dbExecute(con, sprintf(paste(
    "CREATE TABLE %s AS SELECT 'GRCh38' AS assembly, i::UBIGINT AS site_index,",
    "'chr2' AS region, (200 + i)::UBIGINT AS position,",
    "'A' AS allele_a, 'G' AS allele_b FROM range(2) sites(i)"
  ), qid(panel_name)))
  dbExecute(con, sprintf(paste(
    "CREATE TABLE %s AS SELECT p.*, f.population_b_af::DOUBLE AS population_b_af",
    "FROM %s p JOIN (VALUES (0, 0.25), (1, 0.80))",
    "f(site_index, population_b_af) USING (site_index)"
  ), qid(frequency_name), qid(panel_name)))
  dbExecute(con, sprintf(
    "CREATE TABLE %s AS SELECT *, 0.5::DOUBLE AS population_b_af FROM %s",
    qid(match_frequency_name), qid(panel_name)
  ))
  dbExecute(con, sprintf(paste(
    "CREATE TABLE %s AS SELECT c.sample_id::VARCHAR AS sample_id, p.*,",
    "c.a::UBIGINT AS a, c.b::UBIGINT AS b, c.other::UBIGINT AS other",
    "FROM (VALUES",
    "('C', 0, 199, 1, 0), ('C', 1, 2, 198, 0),",
    "('R', 0, 160, 40, 0), ('R', 1, 40, 160, 0),",
    "('K', 0, 199, 1, 0), ('K', 1, 1, 199, 0),",
    "('U', 0, NULL, NULL, NULL), ('U', 1, NULL, NULL, NULL)",
    ") c(sample_id, site_index, a, b, other)",
    "JOIN %s p USING (site_index)"
  ), qid(evidence_name), qid(panel_name)))
  dbExecute(con, sprintf(
    "CREATE TABLE %s(receiver_id VARCHAR, anchor_id VARCHAR)", qid(pairs_name)
  ))
  dbExecute(con, sprintf(
    "INSERT INTO %s VALUES ('R', 'K'), ('K', 'R')", qid(pairs_name)
  ))

  charr <- rduckhts_somalier_charr(
    con, evidence_table = evidence_name, panel_table = panel_name,
    frequency_table = frequency_name, sample_ids = c("C", "U"),
    min_depth = 7, hom_minor_rate = 0.10, hom_tail_alpha = 0.001,
    max_sites = 2
  )
  charr <- charr[order(charr$sample_id), ]
  rownames(charr) <- NULL
  expect_equal(charr$sample_id, c("C", "U"))
  expect_equal(charr$status, c("ok", "no_evidence"))
  expect_equal(charr$usable_sites, c(2, 0))
  expect_equal(charr$usable_hom_a, c(1, 0))
  expect_equal(charr$usable_hom_b, c(1, 0))
  expect_equal(charr$estimate[[1]], 0.035, tolerance = 1e-12)
  expect_true(is.na(charr$estimate[[2]]))
  expect_equal(length(unique(charr$panel_sha256)), 1L)
  expect_equal(length(unique(charr$frequency_sha256)), 1L)

  matched <- rduckhts_somalier_matched_contamination(
    con, evidence_table = evidence_name, panel_table = panel_name,
    frequency_table = match_frequency_name, pairs_table = pairs_name,
    max_sites = 2
  )
  matched <- matched[order(matched$receiver_id, matched$anchor_id), ]
  rownames(matched) <- NULL
  forward <- matched[matched$receiver_id == "R", ]
  reverse <- matched[matched$receiver_id == "K", ]
  expect_equal(forward$status, "ok")
  expect_equal(forward$usable_sites, 2)
  expect_equal(forward$alpha, 0.3975904, tolerance = 2e-6)
  expect_true(forward$evaluations >= 101)
  expect_equal(reverse$status, "no_evidence")
  expect_equal(reverse$usable_sites, 0)
  expect_true(is.na(reverse$alpha))
  expect_true(is.na(reverse$relative_log_likelihood))

  oversized_matched_evidence <- "oversized matched evidence"
  oversized_matched_pairs <- "oversized matched pairs"
  dbExecute(con, sprintf(paste(
    "CREATE TEMP VIEW %s AS SELECT * REPLACE(",
    "CASE WHEN sample_id = 'R' THEN repeat('r', 1025)",
    "ELSE sample_id END AS sample_id) FROM %s"
  ), qid(oversized_matched_evidence), qid(evidence_name)))
  dbExecute(con, sprintf(paste(
    "CREATE TEMP VIEW %s AS SELECT repeat('r', 1025)::VARCHAR AS receiver_id,",
    "'K'::VARCHAR AS anchor_id"
  ), qid(oversized_matched_pairs)))
  expect_error(
    rduckhts_somalier_matched_contamination(
      con, evidence_table = oversized_matched_evidence,
      panel_table = panel_name, frequency_table = match_frequency_name,
      pairs_table = oversized_matched_pairs, max_sites = 2
    ),
    "sample_id and assembly must be at most 1024 bytes"
  )

  # Every contamination relation can be ordinary Parquet. Identity digests and
  # numerical results remain unchanged after the round trip.
  paths <- setNames(
    vapply(c("panel", "evidence", "frequency", "match_frequency", "pairs"), function(kind) {
      tempfile(paste0("somalier ", kind, "'s "), fileext = ".parquet")
    }, character(1)),
    c("panel", "evidence", "frequency", "match_frequency", "pairs")
  )
  on.exit(unlink(paths), add = TRUE)
  sources <- c(panel_name, evidence_name, frequency_name,
               match_frequency_name, pairs_name)
  for (i in seq_along(paths)) {
    dbExecute(con, sprintf("COPY %s TO %s (FORMAT PARQUET)",
                           qid(sources[[i]]), qstr(paths[[i]])))
  }
  charr_parquet <- rduckhts_somalier_charr(
    con, evidence_parquet = paths[["evidence"]],
    panel_parquet = paths[["panel"]],
    frequency_parquet = paths[["frequency"]], sample_ids = c("C", "U"),
    min_depth = 7, hom_minor_rate = 0.10, hom_tail_alpha = 0.001,
    max_sites = 2
  )
  charr_parquet <- charr_parquet[order(charr_parquet$sample_id), ]
  rownames(charr_parquet) <- NULL
  expect_equal(charr_parquet, charr)
  matched_parquet <- rduckhts_somalier_matched_contamination(
    con, evidence_parquet = paths[["evidence"]],
    panel_parquet = paths[["panel"]],
    frequency_parquet = paths[["match_frequency"]],
    pairs_parquet = paths[["pairs"]], max_sites = 2
  )
  matched_parquet <- matched_parquet[
    order(matched_parquet$receiver_id, matched_parquet$anchor_id),
  ]
  rownames(matched_parquet) <- NULL
  expect_equal(matched_parquet, matched)

  charr_table <- "charr output; table"
  matched_table <- 'matched output "table"'
  expect_true(rduckhts_somalier_charr(
    con, evidence_table = evidence_name, panel_table = panel_name,
    frequency_table = frequency_name, table_name = charr_table,
    min_depth = 7, hom_minor_rate = 0.10, hom_tail_alpha = 0.001,
    max_sites = 2
  ))
  expect_true(rduckhts_somalier_matched_contamination(
    con, evidence_table = evidence_name, panel_table = panel_name,
    frequency_table = match_frequency_name, pairs_table = pairs_name,
    table_name = matched_table, max_sites = 2
  ))
  expect_equal(dbGetQuery(con, sprintf("SELECT count(*) AS n FROM %s",
                                       qid(charr_table)))$n, 4)
  expect_equal(dbGetQuery(con, sprintf("SELECT count(*) AS n FROM %s",
                                       qid(matched_table)))$n, 2)

  # A failed overwrite is atomic: malformed frequency orientation does not
  # replace an existing result table.
  preserve <- "contamination_atomic_output"
  bad_frequency <- "frequency_bad_orientation"
  dbExecute(con, sprintf("CREATE TABLE %s AS SELECT 42 AS sentinel", qid(preserve)))
  dbExecute(con, sprintf(paste(
    "CREATE TABLE %s AS SELECT * REPLACE(",
    "CASE WHEN site_index = 0 THEN 'T' ELSE allele_a END AS allele_a)",
    "FROM %s"
  ), qid(bad_frequency), qid(frequency_name)))
  expect_error(
    rduckhts_somalier_charr(
      con, evidence_table = evidence_name, panel_table = panel_name,
      frequency_table = bad_frequency, table_name = preserve,
      min_depth = 7, max_sites = 2, overwrite = TRUE
    ),
    "does not match ordered panel"
  )
  expect_equal(dbGetQuery(con, sprintf("SELECT * FROM %s", qid(preserve))),
               data.frame(sentinel = 42L))

  expect_error(
    rduckhts_somalier_charr(
      con, evidence_table = evidence_name, panel_table = panel_name,
      frequency_table = frequency_name, sample_ids = "missing", max_sites = 2
    ),
    "every sample_ids value"
  )
  expect_error(
    rduckhts_somalier_charr(
      con, evidence_table = evidence_name, panel_table = panel_name,
      frequency_table = frequency_name, min_depth = 16, max_depth = 15,
      max_sites = 2
    ),
    "min_depth"
  )
  expect_error(
    rduckhts_somalier_charr(
      con, evidence_table = evidence_name, panel_table = panel_name,
      frequency_table = frequency_name, hom_minor_rate = 0.5, max_sites = 2
    ),
    "hom_minor_rate"
  )
  expect_error(
    rduckhts_somalier_charr(
      con, evidence_table = evidence_name, panel_table = panel_name,
      frequency_table = frequency_name, max_threshold_work = 0, max_sites = 2
    ),
    "max_threshold_work"
  )
  expect_error(
    rduckhts_somalier_charr(
      con, evidence_table = evidence_name, panel_table = panel_name,
      frequency_table = frequency_name, max_threshold_work = 1, max_sites = 2
    ),
    "binomial threshold certification exceeded max_threshold_work"
  )
  expect_error(
    rduckhts_somalier_matched_contamination(
      con, evidence_table = evidence_name, panel_table = panel_name,
      frequency_table = frequency_name, pairs_table = pairs_name,
      alpha_min = 0.5, alpha_max = 0.5, max_sites = 2
    ),
    "alpha_min"
  )
  expect_error(
    rduckhts_somalier_matched_contamination(
      con, evidence_table = evidence_name, panel_table = panel_name,
      frequency_table = frequency_name, pairs_table = pairs_name,
      error_rate = 0.5, max_sites = 2
    ),
    "error_rate"
  )
  expect_error(
    rduckhts_somalier_matched_contamination(
      con, evidence_table = evidence_name, panel_table = panel_name,
      frequency_table = frequency_name, pairs_table = pairs_name,
      grid_step = 0, max_sites = 2
    ),
    "grid_step"
  )
  expect_error(
    rduckhts_somalier_matched_contamination(
      con, evidence_table = evidence_name, panel_table = panel_name,
      frequency_table = frequency_name, pairs_table = pairs_name,
      max_evaluations = 0, max_sites = 2
    ),
    "max_evaluations"
  )
  expect_error(
    rduckhts_somalier_matched_contamination(
      con, evidence_table = evidence_name, panel_table = panel_name,
      frequency_table = frequency_name, pairs_table = pairs_name,
      max_threshold_work = 100000001, max_sites = 2
    ),
    "max_threshold_work"
  )
  expect_error(
    rduckhts_somalier_matched_contamination(
      con, evidence_table = evidence_name, panel_table = panel_name,
      frequency_table = match_frequency_name, pairs_table = pairs_name,
      max_threshold_work = 1, max_sites = 2
    ),
    "binomial threshold certification exceeded max_threshold_work"
  )
  expect_error(
    rduckhts_somalier_matched_contamination(
      con, evidence_table = evidence_name, panel_table = panel_name,
      frequency_table = frequency_name
    ),
    "pairs: supply exactly one"
  )
}

test_somalier_relatedness_wrappers()
test_somalier_contamination_wrappers()
