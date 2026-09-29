library(tinytest)
library(DBI)

# Each phase opens the database file in its own R process so that connection-local
# TEMP macros and installed state cannot leak between phases.
run_phase <- function(file, body) {
  script <- tempfile("rduckhts_macro_phase_", fileext = ".R")
  result <- tempfile("rduckhts_macro_result_", fileext = ".rds")
  on.exit(unlink(c(script, result)), add = TRUE)
  writeLines(c(
    sprintf(".libPaths(%s)", paste(deparse(.libPaths()), collapse = "")),
    "suppressPackageStartupMessages({ library(DBI); library(Rduckhts) })",
    "macro_catalog <- function(con, catalog) unique(dbGetQuery(con, paste0(",
    "  \"SELECT function_name FROM duckdb_functions() \",",
    "  \"WHERE function_type IN ('macro', 'table_macro') \",",
    "  \"AND database_name = \", dbQuoteString(con, catalog),",
    "  \" AND function_name IN (SELECT name FROM duckhts_macro_definitions())\"))$function_name)",
    sprintf("file <- %s", deparse(file)),
    body,
    sprintf("saveRDS(result, %s)", deparse(result))
  ), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script),
                    stdout = TRUE, stderr = TRUE)
  if (!file.exists(result)) stop(paste(status, collapse = "\n"))
  readRDS(result)
}

file <- tempfile("rduckhts_macros_", fileext = ".duckdb")

# A writable file: LOAD writes nothing to the file; R installs TEMP macros per
# connection, and they see the caller's TEMP tables and CTEs.
writable <- run_phase(file, c(
  "con <- rduckhts_connect(dbdir = file)",
  "catalog <- dbGetQuery(con, 'SELECT current_database()')[[1L]]",
  "result <- list(catalog = catalog,",
  "  temp = length(macro_catalog(con, 'temp')),",
  "  persistent = length(macro_catalog(con, catalog)),",
  "  quoted = dbGetQuery(con, \"SELECT duckhts_quote_ident('a') AS value\")$value)",
  "second <- dbConnect(methods::slot(con, 'driver'))",
  "result$second_before <- 'duckhts_quote_ident' %in% macro_catalog(second, 'temp')",
  "result$second_loaded <- rduckhts_load(second)",
  "result$second_installed <- length(rduckhts_install_macros(second))",
  "result$second_quoted <- dbGetQuery(second, \"SELECT duckhts_quote_ident('a') AS value\")$value",
  "dbExecute(second, \"CREATE TEMP TABLE inputs AS SELECT 'b' AS label\")",
  "result$cte <- dbGetQuery(second, paste('WITH cte AS (SELECT label FROM inputs)',",
  "  'SELECT duckhts_quote_ident(label) AS value FROM cte'))$value",
  "dbDisconnect(second)",
  "result$reinstalled <- length(rduckhts_install_macros(con))",
  "result$temp_after <- length(macro_catalog(con, 'temp'))",
  "dbDisconnect(con, shutdown = TRUE)"
))
expect_equal(writable$temp, 32L)
expect_equal(writable$persistent, 0L)
expect_equal(writable$quoted, '"a"')
expect_false(writable$second_before)
expect_true(writable$second_loaded)
expect_equal(writable$second_installed, 1L)
expect_equal(writable$second_quoted, '"a"')
expect_equal(writable$cte, '"b"')
expect_equal(writable$reinstalled, 1L)
expect_equal(writable$temp_after, 32L)

# A read-only file: LOAD succeeds, installation works inside a transaction that is
# rolled back, and the file is left byte-identical.
readonly_before <- unname(tools::md5sum(file))
readonly <- run_phase(file, c(
  sprintf("catalog <- %s", deparse(writable$catalog)),
  "con <- rduckhts_connect(dbdir = file, read_only = TRUE)",
  "result <- list(temp = length(macro_catalog(con, 'temp')),",
  "  persistent = length(macro_catalog(con, catalog)))",
  "dbBegin(con)",
  "rduckhts_install_macros(con)",
  "dbRollback(con)",
  "result$temp_after <- length(macro_catalog(con, 'temp'))",
  "dbDisconnect(con, shutdown = TRUE)"
))
expect_equal(readonly$temp, 32L)
expect_equal(readonly$persistent, 0L)
expect_equal(readonly$temp_after, 32L)
expect_identical(unname(tools::md5sum(file)), readonly_before)

# Installation inside a caller transaction is rolled back with that transaction.
rollback <- run_phase(file, c(
  "con <- rduckhts_connect(dbdir = file)",
  "dbExecute(con, 'DROP MACRO IF EXISTS temp.duckhts_quote_ident')",
  "result <- list(dropped = 'duckhts_quote_ident' %in% macro_catalog(con, 'temp'))",
  "dbBegin(con)",
  "rduckhts_install_macros(con)",
  "result$inside <- 'duckhts_quote_ident' %in% macro_catalog(con, 'temp')",
  "dbRollback(con)",
  "result$rolled_back <- 'duckhts_quote_ident' %in% macro_catalog(con, 'temp')",
  "rduckhts_install_macros(con)",
  "result$reinstalled <- 'duckhts_quote_ident' %in% macro_catalog(con, 'temp')",
  "dbDisconnect(con, shutdown = TRUE)"
))
expect_false(rollback$dropped)
expect_true(rollback$inside)
expect_false(rollback$rolled_back)
expect_true(rollback$reinstalled)
unlink(file)
