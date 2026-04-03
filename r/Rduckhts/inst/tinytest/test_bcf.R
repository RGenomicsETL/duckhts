library(tinytest)
library(DBI)

test_bcf_string_format_lists <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit({
    dbDisconnect(con, shutdown = TRUE)
  }, add = TRUE)

  expect_silent(rduckhts_load(con))

  bcf_path <- system.file("extdata", "format_string_list.vcf", package = "Rduckhts")
  fixed_path <- system.file("extdata", "fixed_count_arrays.vcf", package = "Rduckhts")
  expect_true(file.exists(bcf_path))
  expect_true(file.exists(fixed_path))

  expect_silent(rduckhts_bcf(con, "bcf_string_lists", bcf_path, overwrite = TRUE))
  expect_silent(rduckhts_bcf(con, "bcf_fixed_counts", fixed_path, overwrite = TRUE))

  type_row <- DBI::dbGetQuery(
    con,
    "SELECT typeof(FORMAT_LAA_SAMPLE1) AS t FROM bcf_string_lists LIMIT 1"
  )
  expect_equal(type_row$t[1], "VARCHAR[]")

  first_row <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT list_extract(FORMAT_LAA_SAMPLE1, 1) AS laa1,",
      "list_extract(FORMAT_LAA_SAMPLE1, 2) AS laa2",
      "FROM bcf_string_lists WHERE POS = 73"
    )
  )
  expect_equal(first_row$laa1[1], "1")
  expect_equal(first_row$laa2[1], "2")

  missing_row <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT FORMAT_LAA_SAMPLE1 IS NULL AS is_null,",
      "FORMAT_LAA_SAMPLE1::VARCHAR AS laa",
      "FROM bcf_string_lists WHERE POS = 263"
    )
  )
  expect_true(isTRUE(missing_row$is_null[1]))
  expect_true(is.na(missing_row$laa[1]))

  fixed_type <- DBI::dbGetQuery(
    con,
    "SELECT typeof(INFO_SB) AS info_t, typeof(FORMAT_HQ_S1) AS fmt_t FROM bcf_fixed_counts LIMIT 1"
  )
  expect_equal(fixed_type$info_t[1], "INTEGER[]")
  expect_equal(fixed_type$fmt_t[1], "INTEGER[]")

  fixed_row <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT INFO_SB::VARCHAR AS info_sb,",
      "FORMAT_HQ_S1::VARCHAR AS format_hq",
      "FROM bcf_fixed_counts WHERE POS = 100"
    )
  )
  expect_equal(fixed_row$info_sb[1], "[10, 20, 30, 40]")
  expect_equal(fixed_row$format_hq[1], "[5, 9]")
}

test_bcf_string_format_lists()
