library(tinytest)
library(DBI)

macro_catalog <- function(con, catalog) {
  dbGetQuery(
    con,
    paste0(
      "SELECT function_name FROM duckdb_functions() ",
      "WHERE function_type IN ('macro', 'table_macro') ",
      "AND database_name = ", dbQuoteString(con, catalog),
      " AND function_name IN (SELECT name FROM duckhts_macro_definitions())"
    )
  )$function_name
}

# Release the database file before reopening it: Windows refuses to reopen a file
# that a live instance still holds, and a disconnected connection object can keep
# its driver's instance alive.
close_database <- function(con) {
  driver <- methods::slot(con, "driver")
  dbDisconnect(con)
  duckdb::duckdb_shutdown(driver)
  invisible(gc())
}

file <- tempfile("rduckhts_macros_", fileext = ".duckdb")
con <- rduckhts_connect(dbdir = file)
catalog <- dbGetQuery(con, "SELECT current_database()")[[1L]]
expect_equal(length(macro_catalog(con, "temp")), 31L)
expect_equal(length(macro_catalog(con, catalog)), 0L)
expect_equal(dbGetQuery(con, "SELECT duckhts_quote_ident('a') AS value")$value, '"a"')
second <- dbConnect(methods::slot(con, "driver"))
expect_false("duckhts_quote_ident" %in% macro_catalog(second, "temp"))
expect_true(rduckhts_load(second))
expect_equal(length(rduckhts_install_macros(second)), 1L)
expect_equal(dbGetQuery(second, "SELECT duckhts_quote_ident('a') AS value")$value, '"a"')
dbExecute(second, "CREATE TEMP TABLE inputs AS SELECT 'b' AS label")
expect_equal(
  dbGetQuery(second, "WITH cte AS (SELECT label FROM inputs) SELECT duckhts_quote_ident(label) AS value FROM cte")$value,
  '"b"'
)
dbDisconnect(second)
rm(second)
expect_equal(length(rduckhts_install_macros(con)), 1L)
expect_equal(length(macro_catalog(con, "temp")), 31L)
close_database(con)

readonly_before <- tools::md5sum(file)
con <- rduckhts_connect(dbdir = file, read_only = TRUE)
expect_equal(length(macro_catalog(con, "temp")), 31L)
expect_equal(length(macro_catalog(con, catalog)), 0L)
dbBegin(con)
rduckhts_install_macros(con)
dbRollback(con)
expect_equal(length(macro_catalog(con, "temp")), 31L)
close_database(con)
expect_identical(unname(tools::md5sum(file)), unname(readonly_before))

# Installation inside a caller transaction is rolled back with that transaction.
con <- rduckhts_connect(dbdir = file)
dbExecute(con, "DROP MACRO IF EXISTS temp.duckhts_quote_ident")
expect_false("duckhts_quote_ident" %in% macro_catalog(con, "temp"))
dbBegin(con)
rduckhts_install_macros(con)
expect_true("duckhts_quote_ident" %in% macro_catalog(con, "temp"))
dbRollback(con)
expect_false("duckhts_quote_ident" %in% macro_catalog(con, "temp"))
rduckhts_install_macros(con)
expect_true("duckhts_quote_ident" %in% macro_catalog(con, "temp"))
close_database(con)
unlink(file)
