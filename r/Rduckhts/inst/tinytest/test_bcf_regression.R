library(tinytest)
library(DBI)

test_bcf_filter_list_fetch_regression <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit({
    dbDisconnect(con, shutdown = TRUE)
  }, add = TRUE)

  expect_silent(rduckhts_load(con))

  bcf_path <- system.file("extdata", "bcf_filter_list_regression.vcf", package = "Rduckhts")
  expect_true(file.exists(bcf_path))

  quoted_bcf <- DBI::dbQuoteString(con, bcf_path)
  for (reader in c("read_bcf", "read_bcf_v2")) {
    rs <- dbSendQuery(
      con,
      sprintf(
        "SELECT POS, FILTER::VARCHAR AS filter_str FROM %s(%s) ORDER BY POS",
        reader, quoted_bcf
      )
    )

    n_rows <- 0L
    saw_multi_filter <- FALSE
    repeat {
      chunk <- dbFetch(rs, n = 512)
      if (nrow(chunk) == 0L) {
        break
      }
      n_rows <- n_rows + nrow(chunk)
      if (any(chunk$filter_str == "[q10, q20]")) {
        saw_multi_filter <- TRUE
      }
    }

    DBI::dbClearResult(rs)

    expect_equal(n_rows, 5000L)
    expect_true(saw_multi_filter)

    summary <- dbGetQuery(
      con,
      sprintf(
        "SELECT COUNT(*) AS n_rows, SUM(length(FILTER)) AS filter_items FROM %s(%s)",
        reader, quoted_bcf
      )
    )
    expect_equal(summary$n_rows[1], 5000L)
    expect_equal(summary$filter_items[1], 5384L)
  }

  mismatch <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "WITH ",
        "v1 AS (SELECT POS, FILTER FROM read_bcf(%s)), ",
        "v2 AS (SELECT POS, FILTER FROM read_bcf_v2(%s)) ",
        "SELECT count(*) AS n FROM ((SELECT * FROM v1 EXCEPT ALL SELECT * FROM v2) ",
        "UNION ALL (SELECT * FROM v2 EXCEPT ALL SELECT * FROM v1))"
      ),
      quoted_bcf, quoted_bcf
    )
  )
  expect_equal(mismatch$n[1], 0)
}

test_bcf_malformed_record_errors <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit({
    dbDisconnect(con, shutdown = TRUE)
  }, add = TRUE)

  expect_silent(rduckhts_load(con))

  bcf_path <- system.file("extdata", "malformed_bad_pos.vcf", package = "Rduckhts")
  expect_true(file.exists(bcf_path))
  quoted_bcf <- DBI::dbQuoteString(con, bcf_path)

  expect_error(
    dbGetQuery(
      con,
      sprintf("SELECT count(*) FROM read_bcf(%s, tidy_format := true)", quoted_bcf)
    ),
    pattern = "read_bcf: failed to read or parse BCF/VCF record"
  )

  expect_error(
    dbGetQuery(
      con,
      sprintf("SELECT sum(POS) FROM read_bcf_v2(%s, tidy_format := true)", quoted_bcf)
    ),
    pattern = "read_bcf_v2: failed to read or parse BCF/VCF record"
  )

  DBI::dbExecute(con, "DROP TABLE IF EXISTS bcf_appender_malformed_bad_pos")
  expect_error(
    dbGetQuery(
      con,
      sprintf(
        paste0(
          "SELECT rows_written FROM read_bcf_appender(",
          "%s, 'bcf_appender_malformed_bad_pos', tidy_format := true, overwrite := true)"
        ),
        quoted_bcf
      )
    ),
    pattern = "read_bcf_appender: failed to read or parse BCF/VCF record"
  )
  appender_tables <- dbGetQuery(
    con,
    paste0(
      "SELECT count(*) AS n FROM information_schema.tables ",
      "WHERE table_name = 'bcf_appender_malformed_bad_pos'"
    )
  )
  expect_equal(appender_tables$n[1], 0)

  DBI::dbExecute(con, "DROP TABLE IF EXISTS bcf_appender_malformed_bad_pos")
}

test_bcf_filter_list_fetch_regression()
test_bcf_malformed_record_errors()
