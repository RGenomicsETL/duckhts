library(tinytest)
library(DBI)

test_sql_quoting <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  work_dir <- file.path(tempdir(), "Rduckhts issue 170's fixtures")
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  source_bcf <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
  quoted_path <- file.path(work_dir, "quoted sample's variants.bcf")
  expect_true(file.copy(source_bcf, quoted_path, overwrite = TRUE))

  quoted_identifier <- 'select "quoted-table"; punctuation'
  expect_silent(rduckhts_bcf(
    con,
    quoted_identifier,
    quoted_path,
    scan_mode = "sequential"
  ))
  expect_true(dbExistsTable(con, quoted_identifier))
  quoted_table <- as.character(dbQuoteIdentifier(con, quoted_identifier))
  expect_true(dbGetQuery(
    con,
    paste("SELECT count(*) > 0 AS ok FROM", quoted_table)
  )$ok[[1]])

  dbExecute(con, "CREATE TABLE issue170_guard(value INTEGER)")
  injected_identifier <- paste(
    "issue170_result AS SELECT 1;",
    "DROP TABLE issue170_guard; --"
  )
  expect_silent(rduckhts_bcf(
    con,
    injected_identifier,
    quoted_path,
    scan_mode = "sequential"
  ))
  expect_true(dbExistsTable(con, injected_identifier))
  expect_true(dbExistsTable(con, "issue170_guard"))

  plain_path <- file.path(work_dir, "plain.bcf")
  injected_path <- paste0(
    plain_path,
    "') ; DROP TABLE issue170_guard; --.bcf"
  )
  expect_true(file.copy(source_bcf, plain_path, overwrite = TRUE))
  expect_true(file.copy(source_bcf, injected_path, overwrite = TRUE))
  expect_silent(rduckhts_bcf(
    con,
    "injected_path_is_data",
    injected_path,
    scan_mode = "sequential"
  ))
  expect_true(dbExistsTable(con, "issue170_guard"))
}

test_sql_quoting()
