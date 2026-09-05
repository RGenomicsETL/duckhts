library(tinytest)
library(DBI)

test_bcf_filter_list_fetch_regression <- function() {
  con <- rduckhts_connect()
  on.exit(
    {
      dbDisconnect(con, shutdown = TRUE)
    },
    add = TRUE
  )

  bcf_path <- system.file(
    "extdata",
    "bcf_filter_list_regression.vcf",
    package = "Rduckhts"
  )
  expect_true(file.exists(bcf_path))

  quoted_bcf <- DBI::dbQuoteString(con, bcf_path)
  rs <- dbSendQuery(
    con,
    sprintf(
      "SELECT POS, FILTER::VARCHAR AS filter_str FROM read_bcf(%s) ORDER BY POS",
      quoted_bcf
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
      "SELECT COUNT(*) AS n_rows, SUM(length(FILTER)) AS filter_items FROM read_bcf(%s)",
      quoted_bcf
    )
  )
  expect_equal(summary$n_rows[1], 5000L)
  expect_equal(summary$filter_items[1], 5384L)

}

test_bcf_malformed_record_errors <- function() {
  con <- rduckhts_connect()
  on.exit(
    {
      dbDisconnect(con, shutdown = TRUE)
    },
    add = TRUE
  )

  bcf_path <- system.file(
    "extdata",
    "malformed_bad_pos.vcf",
    package = "Rduckhts"
  )
  expect_true(file.exists(bcf_path))
  quoted_bcf <- DBI::dbQuoteString(con, bcf_path)

  expect_error(
    dbGetQuery(
      con,
      sprintf(
        "SELECT count(*) FROM read_bcf(%s, tidy_format := true)",
        quoted_bcf
      )
    ),
    pattern = "read_bcf: failed to read or parse BCF/VCF record"
  )

  expect_error(
    dbGetQuery(
      con,
      sprintf(
        "SELECT sum(POS) FROM read_bcf(%s, tidy_format := true)",
        quoted_bcf
      )
    ),
    pattern = "read_bcf: failed to read or parse BCF/VCF record"
  )

  expect_error(dbExecute(con, sprintf(
    "CREATE TABLE bcf_malformed_materialized AS SELECT * FROM read_bcf(%s, tidy_format := true)",
    quoted_bcf)), pattern = "read_bcf: failed to read or parse BCF/VCF record")
  expect_equal(dbGetQuery(con, paste(
    "SELECT count(*) AS n FROM information_schema.tables",
    "WHERE table_name = 'bcf_malformed_materialized'"))$n[[1]], 0)
}

test_bcf_type_clash_errors <- function() {
  con <- rduckhts_connect()
  on.exit(
    {
      dbDisconnect(con, shutdown = TRUE)
    },
    add = TRUE
  )

  info_path <- system.file(
    "extdata",
    "bcf_info_type_clash.bcf",
    package = "Rduckhts"
  )
  format_path <- system.file(
    "extdata",
    "bcf_format_type_clash.bcf",
    package = "Rduckhts"
  )
  str_path <- system.file(
    "extdata",
    "bcf_info_str_clash.bcf",
    package = "Rduckhts"
  )
  expect_true(file.exists(info_path))
  expect_true(file.exists(format_path))
  expect_true(file.exists(str_path))

  quoted_info <- DBI::dbQuoteString(con, info_path)
  quoted_format <- DBI::dbQuoteString(con, format_path)
  quoted_str <- DBI::dbQuoteString(con, str_path)

  null_policy_info <- dbGetQuery(
    con,
    sprintf("SELECT INFO_DP IS NULL AS is_null FROM read_bcf(%s)", quoted_info)
  )
  expect_true(null_policy_info$is_null[1])

  # Reverse clash: String header over numeric payload. htslib's string decode
  # does not type-check, so the reader must NULL this under the default policy
  # rather than leak raw bytes.
  null_str_info <- dbGetQuery(
    con,
    sprintf("SELECT INFO_NN IS NULL AS is_null FROM read_bcf(%s)", quoted_str)
  )
  expect_true(null_str_info$is_null[1])
  null_policy_format <- dbGetQuery(
    con,
    sprintf(
      "SELECT FORMAT_XX_S1 IS NULL AS s1_null, FORMAT_XX_S2 IS NULL AS s2_null FROM read_bcf(%s)",
      quoted_format
    )
  )
  expect_true(null_policy_format$s1_null[1])
  expect_true(null_policy_format$s2_null[1])

  warn_format <- dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT FORMAT_XX_S1 IS NULL AS s1_null, FORMAT_XX_S2 IS NULL AS s2_null ",
        "FROM read_bcf(%s, decode_error_policy := 'warn')"
      ),
      quoted_format
    )
  )
  expect_true(warn_format$s1_null[1])
  expect_true(warn_format$s2_null[1])

  expect_error(
    dbGetQuery(
      con,
      sprintf(
        "SELECT INFO_DP FROM read_bcf(%s, decode_error_policy := 'error')",
        quoted_info
      )
    ),
    pattern = "read_bcf: INFO/DP encoded BCF type CHAR does not match header Type=Integer at chr1:10"
  )

  expect_error(
    dbGetQuery(
      con,
      sprintf(
        "SELECT INFO_NN FROM read_bcf(%s, decode_error_policy := 'error')",
        quoted_str
      )
    ),
    pattern = "read_bcf: INFO/NN encoded BCF type INT8 does not match header Type=String at chr1:10"
  )
  expect_error(
    dbGetQuery(
      con,
      sprintf(
        "SELECT FORMAT_XX_S1 FROM read_bcf(%s, decode_error_policy := 'error')",
        quoted_format
      )
    ),
    pattern = "read_bcf: FORMAT/XX encoded BCF type CHAR does not match header Type=Integer at chr1:10"
  )

}

test_bcf_filter_list_fetch_regression()
test_bcf_malformed_record_errors()
test_bcf_type_clash_errors()
