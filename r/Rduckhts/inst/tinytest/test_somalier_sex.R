library(tinytest)
library(DBI)

# Ten autosomal, twelve X and three Y sites; profiles are XX, XY and edge cases.
somalier_sex_fixture <- function(con) {
  dbExecute(con, paste(
    "CREATE TEMP TABLE sex_panel AS SELECT 'GRCh38' AS assembly,",
    "i::UBIGINT AS site_index,",
    "CASE WHEN i < 10 THEN 'chr1' WHEN i < 22 THEN 'chrX' ELSE 'chrY' END AS region,",
    "(1000 + i)::UBIGINT AS position, 'A' AS allele_a, 'G' AS allele_b",
    "FROM range(25) sites(i)"
  ))
  dbExecute(con, paste(
    "CREATE TEMP TABLE sex_counts AS SELECT p.sample_id, s.assembly, s.site_index,",
    "s.region, s.position, s.allele_a, s.allele_b,",
    "(CASE WHEN s.region = 'chr1' THEN 15 WHEN s.region = 'chrX' THEN",
    "  CASE WHEN p.kind = 'xy' THEN CASE WHEN s.site_index = 21 THEN 0 ELSE 20 END",
    "       WHEN p.kind = 'none' THEN 1 WHEN (s.site_index - 10) % 2 = 0 THEN 10 ELSE 0 END",
    " ELSE p.y_depth END)::UBIGINT AS a,",
    "(CASE WHEN s.region = 'chr1' THEN 15",
    "  WHEN s.region = 'chrX' THEN CASE WHEN p.kind = 'xy' THEN CASE WHEN s.site_index = 21 THEN 20 ELSE 0 END WHEN p.kind = 'none' THEN 1",
    "       WHEN (s.site_index - 10) % 2 = 0 THEN 10 ELSE 20 END",
    " ELSE 0 END)::UBIGINT AS b, 0::UBIGINT AS other",
    "FROM (VALUES ('xx', 'xx', 0), ('xy', 'xy', 10), ('none', 'none', 0),",
    "  ('xx_y', 'xx', 15)) p(sample_id, kind, y_depth) CROSS JOIN sex_panel s"
  ))
}

test_somalier_sex_wrapper <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  somalier_sex_fixture(con)

  sex <- rduckhts_somalier_sex(con, counts_table = "sex_counts", panel_table = "sex_panel")
  expect_equal(sex$sample_id, c("none", "xx", "xx_y", "xy"))
  expect_equal(sex$inferred_sex, c("ambiguous", "XX", "ambiguous", "XY"))
  expect_equal(sex$status, c("no_usable_x_sites", "ok", "ok", "ok"))
  expect_equal(sex$x_usable_sites, c(0, 12, 12, 12))
  expect_equal(sex$x_het, c(0, 6, 6, 0))
  expect_equal(sex$x_hom_alt, c(0, 6, 6, 1))
  expect_true(is.na(sex$x_depth_ratio[[1L]]))
  # 2 * mean X depth 20 / mean autosomal depth 30
  expect_equal(sex$x_depth_ratio[[2L]], 2 * 20 / 30, tolerance = 1e-12)
  expect_equal(sex$y_signal, c(NA, NA, "present", "present"))
  expect_true(all(grepl("^review_only", sex$interpretation)))
  expect_equal(unique(sex$autosomal_sites), 10)
  expect_equal(
    unique(sex$panel_sha256),
    dbGetQuery(con, "SELECT duckhts_somalier_panel_sha256('sex_panel') AS h")$h
  )

  # The cohort gate reproduces Somalier's cohort-dependent Y check.
  cohort <- rduckhts_somalier_sex(
    con, counts_table = "sex_counts", panel_table = "sex_panel", y_gate = "cohort"
  )
  expect_equal(cohort$cohort_has_y[[1L]], TRUE)

  subset <- rduckhts_somalier_sex(
    con, counts_table = "sex_counts", panel_table = "sex_panel", sample_ids = "xy"
  )
  expect_equal(subset$inferred_sex, "XY")

  expect_true(rduckhts_somalier_sex(
    con, counts_table = "sex_counts", panel_table = "sex_panel",
    table_name = "sex_result"
  ))
  expect_equal(dbGetQuery(con, "SELECT count(*) AS n FROM sex_result")$n, 4)
  expect_error(rduckhts_somalier_sex(
    con, counts_table = "sex_counts", panel_table = "sex_panel",
    table_name = "sex_result"
  ))
}

test_somalier_sex_arguments <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  somalier_sex_fixture(con)
  call <- function(...) {
    rduckhts_somalier_sex(con, counts_table = "sex_counts",
      panel_table = "sex_panel", ...)
  }
  expect_error(call(y_gate = "always"), pattern = "should be one of")
  expect_error(call(min_depth = 0), pattern = "min_depth")
  expect_error(call(min_usable_x_sites = 1.5), pattern = "min_usable_x_sites")
  expect_error(call(min_het_balance = 0.7), pattern = "min_het_balance")
  expect_error(call(hom_balance_cutoff = 0.4), pattern = "must not exceed")
  expect_error(call(xy_max_het_ratio = 0.5, xx_min_het_ratio = 0.1),
    pattern = "must not exceed")
  expect_error(call(y_signal_min = -1), pattern = "non-negative")
  expect_error(call(sample_ids = "missing"), pattern = "must occur")
  expect_error(rduckhts_somalier_sex(con, panel_table = "sex_panel"),
    pattern = "counts: supply exactly one")
}

# The autosomal functions ignore X/Y rows and report the same numbers as on the
# autosomal-only panel.
test_somalier_autosomal_restriction <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  somalier_sex_fixture(con)
  dbExecute(con, "CREATE TEMP TABLE auto_panel AS SELECT * FROM sex_panel WHERE region = 'chr1'")
  dbExecute(con, "CREATE TEMP TABLE auto_counts AS SELECT * FROM sex_counts WHERE region = 'chr1'")
  full <- rduckhts_somalier_sketches(
    con, evidence_table = "sex_counts", panel_table = "sex_panel",
    min_depth = 7, min_het_balance = 0.3, hom_balance_cutoff = 0.01
  )
  auto <- rduckhts_somalier_sketches(
    con, evidence_table = "auto_counts", panel_table = "auto_panel",
    min_depth = 7, min_het_balance = 0.3, hom_balance_cutoff = 0.01
  )
  expect_equal(full$sample_id, auto$sample_id)
  expect_equal(full$site_count, auto$site_count)
  expect_equal(full$het, auto$het)
  expect_equal(full$hom_a, auto$hom_a)
  expect_equal(full$hom_b, auto$hom_b)
}

test_somalier_sex_wrapper()
test_somalier_sex_arguments()
test_somalier_autosomal_restriction()
