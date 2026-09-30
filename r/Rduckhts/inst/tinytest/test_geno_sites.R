library(tinytest)
library(DBI)

test_geno_regions_var_and_sites <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  dbExecute(con, "SET threads=4")
  fixture <- function(name) {
    path <- system.file("extdata", name, package = "Rduckhts")
    stopifnot(nzchar(path), file.exists(path))
    path
  }
  empty_plan <- "[]::STRUCT(chrom VARCHAR, start BIGINT, \"end\" BIGINT)[]"
  variables <- function() dbGetQuery(con, "SELECT name FROM duckdb_variables() WHERE name LIKE 'duckhts_geno_sites%'")$name
  temp_tables <- function() dbGetQuery(con, "SELECT table_name FROM duckdb_tables() WHERE temporary AND table_name LIKE 'duckhts_geno_sites%'")$table_name

  # ---- regions_var on rduckhts_geno() and rduckhts_bcf() -------------------
  for (extension in c("vcf.gz", "bcf")) {
    path <- fixture(paste0("geno_sites.", extension))
    dbExecute(con, paste0("SET VARIABLE rv_plan = [{'chrom': 'chr1', 'start': 99, 'end': 100}, ",
                          "{'chrom': 'chr1', 'start': 500, 'end': 502}, ",
                          "{'chrom': 'HLA-A*01:01', 'start': 9, 'end': 10}, {'chrom': 'chr1:2', 'start': 9, 'end': 10}]"))
    typed <- rduckhts_geno(con, path = path, regions_var = "rv_plan", format_fields = c("GP", "DS", "HS"))
    string <- rduckhts_geno(con, path = path, region = "chr1:100-100,chr1:501-502,{HLA-A*01:01}:10-10,{chr1:2}:10-10",
                            format_fields = c("GP", "DS", "HS"))
    expect_equal(typed$record_index, string$record_index)
    expect_equal(typed$POS, c(100, 100, 501, 10, 10))
    expect_equal(typed$CHROM, c("chr1", "chr1", "chr1", "HLA-A*01:01", "chr1:2"))
    expect_equal(nrow(typed), 5L)
    expect_equal(
      dbGetQuery(con, paste0("SELECT count(*) AS n FROM ((SELECT * FROM read_geno('", path,
        "', regions := getvariable('rv_plan')) EXCEPT ALL SELECT * FROM read_geno('", path,
        "', region := 'chr1:100-100,chr1:501-502,{HLA-A*01:01}:10-10,{chr1:2}:10-10')) UNION ALL ",
        "(SELECT * FROM read_geno('", path, "', region := 'chr1:100-100,chr1:501-502,{HLA-A*01:01}:10-10,{chr1:2}:10-10') ",
        "EXCEPT ALL SELECT * FROM read_geno('", path, "', regions := getvariable('rv_plan'))))"))$n, 0)

    expect_true(rduckhts_bcf(con, "regions_bcf", path, regions_var = "rv_plan", overwrite = TRUE))
    expect_equal(dbGetQuery(con, "SELECT CHROM, POS FROM regions_bcf ORDER BY CHROM, POS"),
                 data.frame(CHROM = c("HLA-A*01:01", "chr1", "chr1", "chr1", "chr1:2"),
                            POS = c(10, 100, 100, 501, 10)))

    # NULL keeps the ordinary scan; an empty list selects nothing, even for COUNT(*).
    dbExecute(con, "SET VARIABLE rv_null = NULL::STRUCT(chrom VARCHAR, start BIGINT, \"end\" BIGINT)[]")
    expect_equal(nrow(rduckhts_geno(con, path = path, regions_var = "rv_null")), 9L)
    dbExecute(con, paste0("SET VARIABLE rv_empty = ", empty_plan))
    expect_equal(nrow(rduckhts_geno(con, path = path, regions_var = "rv_empty")), 0L)
    expect_true(rduckhts_bcf(con, "regions_empty", path, regions_var = "rv_empty", overwrite = TRUE))
    expect_equal(dbGetQuery(con, "SELECT count(*) AS n FROM regions_empty")$n, 0)
    expect_equal(dbGetQuery(con, paste0("SELECT count(*) AS n FROM read_bcf('", path,
                                        "', regions := getvariable('rv_empty'))"))$n, 0)

    # Argument and native errors.
    expect_error(rduckhts_geno(con, path = path, regions_var = "rv_plan", region = "chr1"),
                 pattern = "mutually exclusive")
    expect_error(rduckhts_bcf(con, "x", path, regions_var = "rv_plan", region = "chr1", overwrite = TRUE),
                 pattern = "mutually exclusive")
    expect_error(rduckhts_geno(con, path = path, regions_var = "rv_plan", scan_mode = "sequential"),
                 pattern = "incompatible with regions")
    expect_error(rduckhts_geno(con, path = path, regions_var = NA_character_), pattern = "regions_var")
    expect_error(rduckhts_geno(con, path = path, regions_var = c("a", "b")), pattern = "regions_var")
    expect_error(rduckhts_geno(con, path = path, regions_var = ""), pattern = "regions_var")
    expect_error(rduckhts_geno(con, path = path, regions_var = 1L), pattern = "regions_var")
    expect_error(rduckhts_bcf(con, "x", path, regions_var = NA_character_, overwrite = TRUE), pattern = "regions_var")
    dbExecute(con, "SET VARIABLE rv_bad = [{'chrom': 'chr1', 'start': 5, 'end': 5}]")
    expect_error(rduckhts_geno(con, path = path, regions_var = "rv_bad"), pattern = "end <= start")
    dbExecute(con, "SET VARIABLE rv_wrong_type = 42")
    expect_error(rduckhts_geno(con, path = path, regions_var = "rv_wrong_type"))
  }
  # A variable name is a string literal, never SQL text.
  hostile <- "x'); DROP TABLE regions_bcf; --"
  expect_error(rduckhts_geno(con, path = fixture("geno_sites.bcf"), regions_var = hostile),
               pattern = "no session variable")
  expect_true("regions_bcf" %in% dbListTables(con))
  dbExecute(con, paste0("SET VARIABLE ", dbQuoteIdentifier(con, hostile), " = [{'chrom': 'chr2', 'start': 49, 'end': 50}]"))
  expect_equal(rduckhts_geno(con, path = fixture("geno_sites.bcf"), regions_var = hostile)$POS, 50)
  expect_true("regions_bcf" %in% dbListTables(con))
  # A misspelt variable is an error, never a silent full scan.
  expect_error(rduckhts_geno(con, path = fixture("geno_sites.bcf"), regions_var = "rv_typo"),
               pattern = "no session variable")
  expect_equal(nrow(rduckhts_geno(con, path = fixture("geno_calls.bcf"),
                                  regions_var = "rv_plan", index_path = fixture("geno_calls.bcf.csi"))), 0L)
  # No index: explicit error, not a scan.
  expect_error(rduckhts_geno(con, path = fixture("geno_ps_type.bcf"), regions_var = "rv_plan"),
               pattern = "requires an index")

  # ---- rduckhts_geno_sites() -----------------------------------------------
  requests <- data.frame(
    request_id = 1:13, build = "GRCh38",
    chrom = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "HLA-A*01:01",
              "chrUnknown", "chr1", "chr1:2", "chr2"),
    pos = c(100L, 100L, 100L, 100L, 300L, 502L, 501L, 100L, 10L, 5L, 100L, 10L, 50L),
    ref = c("A", "A", "A", "C", "G", "TA", "AT", "G", "A", "A", "A", "C", "C"),
    alt = c("C", "G", "T", "G", "T", "T", "A", "C", "G", "C", "C", "T", "T"),
    note = paste0("n", 1:13), stringsAsFactors = FALSE)
  dbWriteTable(con, "site_requests", requests, overwrite = TRUE)
  expected <- data.frame(
    request_id = c(1L, 1L, 2L, 2L, 3L, 4L, 5L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 11L, 12L, 13L),
    record_index = c(0, 1, 0, 1, NA, NA, 2, 2, NA, 3, NA, 6, NA, 0, 1, 7, 4),
    full_alt = c("[C, G]", "[C, G]", "[C, G]", "[C, G]", NA, NA, "[T, T]", "[T, T]", NA, "[A]", NA, "[G]",
                 NA, "[C, G]", "[C, G]", "[T]", "[G, T]"),
    alt_index = c(1, 1, 2, 2, NA, NA, 1, 2, NA, 1, NA, 1, NA, 1, 1, 1, 2),
    match_status = c("matched", "matched", "matched", "matched", "allele_not_at_site", "ref_mismatch",
                     "matched", "matched", "absent", "matched", "allele_not_at_site", "matched", "absent",
                     "matched", "matched", "matched", "matched"),
    stringsAsFactors = FALSE)
  summarise <- function(result) {
    result <- result[order(result$request_id, result$record_index, result$alt_index), ]
    data.frame(request_id = result$request_id, record_index = as.numeric(result$record_index),
               full_alt = vapply(result$full_alt, function(x) if (is.null(x)) NA_character_ else
                 paste0("[", paste(x, collapse = ", "), "]"), character(1)),
               alt_index = as.numeric(result$alt_index), match_status = result$match_status,
               stringsAsFactors = FALSE, row.names = NULL)
  }
  for (extension in c("bcf", "vcf.gz")) {
    path <- fixture(paste0("geno_sites.", extension))
    result <- rduckhts_geno_sites(con, path, "site_requests", "GRCh38",
                                  format_fields = c("GP", "DS", "HS"))
    expect_equal(names(result), c(names(requests), "record_index", "full_alt", "alt_index", "calls", "match_status"))
    expect_equal(summarise(result), expected, check.attributes = FALSE)
    expect_equal(sort(result$note[result$match_status == "matched"]), sort(paste0("n", c(1, 1, 2, 2, 5, 5, 7, 9, 11, 11, 12, 13))))
    # Requests keep their identity and columns; unmatched rows have no calls.
    unmatched <- result[result$match_status != "matched", ]
    expect_true(all(vapply(unmatched$calls, is.null, logical(1))))
    expect_true(all(is.na(unmatched$record_index)) && all(is.na(unmatched$alt_index)))
    expect_equal(nrow(result), 17L)
    expect_equal(unique(result$note[result$request_id == 11]), "n11")
    # Full-site calls: GP (Number=G), DS (Number=A), HS and phase are preserved verbatim.
    multi <- dbGetQuery(con, paste0(
      "SELECT c.sample_index AS sample, c.alleles::VARCHAR AS alleles, c.phase_before::VARCHAR AS phase, ",
      "c.format.GP::VARCHAR AS gp, c.format.DS::VARCHAR AS ds, c.format.HS::VARCHAR AS hs ",
      "FROM (SELECT unnest(calls) AS c FROM read_geno('", path, "', regions := [{'chrom': 'chr2', 'start': 49, 'end': 50}], ",
      "format_fields := ['GP', 'DS', 'HS'])) ORDER BY sample"))
    sites_multi <- rduckhts_geno_sites(con, path, DBI::SQL(
      "SELECT 1 AS request_id, 'GRCh38' AS build, 'chr2' AS chrom, 50 AS pos, 'C' AS ref, 'T' AS alt"),
      "GRCh38", format_fields = c("GP", "DS", "HS"))
    expect_equal(nrow(sites_multi), 1L)
    expect_equal(sites_multi$alt_index, 2)
    expect_equal(unlist(sites_multi$full_alt), c("G", "T"))
    expect_equal(sites_multi$calls[[1]]$alleles, list(c(0L, 1L), c(1L, 2L), c(2L, 0L)))
    expect_equal(sites_multi$calls[[1]]$phase_before, list(c(TRUE, TRUE), c(FALSE, FALSE), c(TRUE, TRUE)))
    expect_equal(sites_multi$calls[[1]]$phase_set, c(50, NA, 50))
    expect_equal(sites_multi$calls[[1]]$format$GP[[1]], c(0.1, 0.8, 0.05, 0.02, 0.02, 0.01), tolerance = 1e-6)
    expect_equal(sites_multi$calls[[1]]$format$DS[[3]], c(0.4, 0.7), tolerance = 1e-6)
    expect_equal(sites_multi$calls[[1]]$format$HS[[2]], c(4L, 3L, 2L, 1L))
    expect_equal(multi$hs, c("[1, 2, 3, 4]", "[4, 3, 2, 1]", "[1, 1, 2, 2]"))

    # Missing (./.) and hom-ref calls are calls, not absent records.
    hom <- rduckhts_geno_sites(con, path, DBI::SQL(
      "SELECT 1 AS request_id, 'GRCh38' AS build, 'chr1' AS chrom, 200 AS pos, 'T' AS ref, 'C' AS alt UNION ALL
       SELECT 2, 'GRCh38', 'chr1', 100, 'A', 'C'"), "GRCh38")
    hom <- hom[order(hom$request_id, hom$record_index), ]
    expect_equal(hom$match_status, c("matched", "matched", "matched"))
    expect_equal(hom$calls[[1]]$alleles, list(c(0L, 0L), c(0L, 1L), c(1L, 1L)))
    expect_equal(hom$calls[[2]]$alleles[[3]], c(NA_integer_, NA_integer_))
    expect_equal(hom$calls[[2]]$phase_before[[3]], c(FALSE, FALSE))

    # Sample selection; DBI::Id and relation-name inputs.
    selected <- rduckhts_geno_sites(con, path, DBI::Id(table_name = "site_requests"), "GRCh38", samples = "S2")
    expect_equal(nrow(selected), 17L)
    expect_equal(unique(vapply(selected$calls[selected$match_status == "matched"],
                               function(x) x$sample_index, numeric(1))), 1)
    explicit_index <- rduckhts_geno_sites(con, path, "site_requests", "GRCh38",
      index_path = fixture(paste0("geno_sites.", if (extension == "bcf") "bcf.csi" else "vcf.gz.tbi")))
    expect_equal(summarise(explicit_index), expected, check.attributes = FALSE)
  }

  # table_name / overwrite: a regular table, its name returned invisibly.
  path <- fixture("geno_sites.bcf")
  expect_equal(withVisible(rduckhts_geno_sites(con, path, "site_requests", "GRCh38", table_name = "site_hits"))$visible, FALSE)
  expect_equal(rduckhts_geno_sites(con, path, "site_requests", "GRCh38", table_name = "site_hits", overwrite = TRUE), "site_hits")
  expect_equal(dbGetQuery(con, "SELECT count(*) AS n FROM site_hits")$n, 17)
  expect_true("site_hits" %in% dbListTables(con))
  expect_error(rduckhts_geno_sites(con, path, "site_requests", "GRCh38", table_name = "site_hits"),
               pattern = "already exists")
  expect_equal(dbGetQuery(con, "SELECT count(*) AS n FROM site_hits")$n, 17)
  expect_equal(dbGetQuery(con, "SELECT count(*) AS n FROM site_hits WHERE match_status = 'absent'")$n, 2)
  expect_equal(dbGetQuery(con, "SELECT count(*) AS n FROM duckdb_tables() WHERE table_name = 'site_hits' AND temporary")$n, 0)

  # No requests: an empty plan, an empty result and the requested columns.
  dbWriteTable(con, "no_requests", requests[0, ], overwrite = TRUE)
  none <- rduckhts_geno_sites(con, path, "no_requests", "GRCh38")
  expect_equal(nrow(none), 0L)
  expect_true(all(c("request_id", "record_index", "full_alt", "alt_index", "calls", "match_status") %in% names(none)))

  # ---- validation and cleanup ------------------------------------------------
  bad <- function(sql, pattern, ...) {
    expect_error(rduckhts_geno_sites(con, path, DBI::SQL(sql), "GRCh38", ...), pattern = pattern)
    expect_equal(variables(), character())
    expect_equal(temp_tables(), character())
  }
  head_sql <- "SELECT 1 AS request_id, 'GRCh38' AS build, 'chr1' AS chrom, 100 AS pos, 'A' AS ref, 'C' AS alt"
  bad(paste(head_sql, "UNION ALL SELECT 1, 'GRCh38', 'chr1', 200, 'T', 'C'"), "request_id must be non-NULL and unique")
  bad(paste(head_sql, "UNION ALL SELECT NULL, 'GRCh38', 'chr1', 200, 'T', 'C'"), "request_id must be non-NULL and unique")
  bad(paste(head_sql, "UNION ALL SELECT 2, 'GRCh37', 'chr1', 200, 'T', 'C'"), "build must equal source_build")
  bad(paste(head_sql, "UNION ALL SELECT 2, NULL, 'chr1', 200, 'T', 'C'"), "build must equal source_build")
  bad(paste(head_sql, "UNION ALL SELECT 2, 'GRCh38', NULL, 200, 'T', 'C'"), "chrom must be non-NULL")
  bad(paste(head_sql, "UNION ALL SELECT 2, 'GRCh38', '', 200, 'T', 'C'"), "chrom must be non-NULL")
  bad(paste(head_sql, "UNION ALL SELECT 2, 'GRCh38', 'chr1', 0, 'T', 'C'"), "pos must be non-NULL and at least 1")
  bad(paste(head_sql, "UNION ALL SELECT 2, 'GRCh38', 'chr1', -5, 'T', 'C'"), "pos must be non-NULL and at least 1")
  bad(paste(head_sql, "UNION ALL SELECT 2, 'GRCh38', 'chr1', NULL, 'T', 'C'"), "pos must be non-NULL and at least 1")
  bad(paste(head_sql, "UNION ALL SELECT 2, 'GRCh38', 'chr1', 5, '', 'C'"), "ref must be non-NULL and non-empty")
  bad(paste(head_sql, "UNION ALL SELECT 2, 'GRCh38', 'chr1', 5, NULL, 'C'"), "ref must be non-NULL and non-empty")
  bad(paste(head_sql, "UNION ALL SELECT 2, 'GRCh38', 'chr1', 5, 'T', ''"), "alt must be one non-empty allele")
  bad(paste(head_sql, "UNION ALL SELECT 2, 'GRCh38', 'chr1', 5, 'T', 'C,G'"), "alt must be one non-empty allele")
  bad(paste(head_sql, "UNION ALL SELECT 2, 'GRCh38', 'chr1', 5, 'T', NULL"), "alt must be one non-empty allele")
  bad("SELECT 1 AS request_id, 'GRCh38' AS build, 'chr1' AS chrom, 1.5 AS pos, 'A' AS ref, 'C' AS alt", "integer type")
  bad("SELECT 1 AS request_id, 'GRCh38' AS build, 'chr1' AS chrom, '100' AS pos, 'A' AS ref, 'C' AS alt", "integer type")
  bad("SELECT 1 AS request_id, 'GRCh38' AS build, 'chr1' AS chrom, 100 AS pos, 'A' AS ref", "missing required column")
  bad(paste0(head_sql, ", 1 AS calls"), "reserved output column")
  bad("SELEC nonsense", "Parser Error|syntax")
  # An error raised by the native scan still cleans up the plan variable.
  expect_error(rduckhts_geno_sites(con, fixture("geno_ps_type.bcf"), DBI::SQL(head_sql), "GRCh38"),
               pattern = "requires an index")
  expect_equal(variables(), character())
  expect_equal(temp_tables(), character())
  expect_error(rduckhts_geno_sites(con, path, DBI::SQL(head_sql), "GRCh38", format_fields = "NOPE"),
               pattern = "not declared")
  expect_equal(variables(), character())
  expect_error(rduckhts_geno_sites(con, path, "site_requests", NA_character_), pattern = "source_build")
  expect_error(rduckhts_geno_sites(con, path, "site_requests", c("a", "b")), pattern = "source_build")
  expect_error(rduckhts_geno_sites(con, path, 42, "GRCh38"), pattern = "sites must be")
  expect_error(rduckhts_geno_sites(con, path, "site_requests", "GRCh38", overwrite = NA), pattern = "overwrite")
  expect_error(rduckhts_geno_sites(con, path, "no_such_relation", "GRCh38"))
  expect_equal(variables(), character())
  expect_equal(temp_tables(), character())
  # A cleaned session leaves ordinary variables untouched.
  expect_equal(dbGetQuery(con, "SELECT getvariable('rv_plan') IS NOT NULL AS kept")$kept, TRUE)
}
test_geno_regions_var_and_sites()
