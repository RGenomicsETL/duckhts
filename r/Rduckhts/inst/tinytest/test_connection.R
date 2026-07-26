library(tinytest)
library(DBI)

expect_error(rduckhts_connect(config = "threads=1"), "config must be a named list")
expect_error(rduckhts_connect(config = list("1")), "config must be a named list")

.connection_test_driver <- function(dbdir, config) {
  driver_args <- list(dbdir = dbdir, config = config)
  duckdb_formals <- names(formals(duckdb::duckdb))
  if ("allow_extensions" %in% duckdb_formals) {
    driver_args$allow_extensions <- FALSE
  }
  if ("shared_home" %in% duckdb_formals) {
    driver_args$shared_home <- FALSE
  }
  do.call(duckdb::duckdb, driver_args)
}

.connection_test_cleanup <- function(con, drv, dbdir) {
  if (!is.null(con) && DBI::dbIsValid(con)) {
    try(DBI::dbDisconnect(con), silent = TRUE)
  }
  if (!is.null(drv) && DBI::dbIsValid(drv)) {
    try(duckdb::duckdb_shutdown(drv), silent = TRUE)
  }
  unlink(dbdir)
}

test_package_owned_connection <- function() {
  duckdb_formals <- names(formals(duckdb::duckdb))
  if ("allow_extensions" %in% duckdb_formals) {
    old_allow_extensions <- getOption("duckdb.allow_extensions")
    options(duckdb.allow_extensions = FALSE)
    on.exit(
      options(duckdb.allow_extensions = old_allow_extensions),
      add = TRUE
    )
  }

  con <- rduckhts_connect(config = list(
    threads = "1",
    autoinstall_known_extensions = "true",
    autoload_known_extensions = "true"
  ))
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  settings <- dbGetQuery(
    con,
    paste(
      "SELECT name, lower(value) AS value FROM duckdb_settings()",
      "WHERE name IN ('allow_unsigned_extensions',",
      "'autoinstall_known_extensions', 'autoload_known_extensions')"
    )
  )
  values <- setNames(settings$value, settings$name)
  expect_identical(values[["allow_unsigned_extensions"]], "true")
  expect_identical(values[["autoinstall_known_extensions"]], "false")
  expect_identical(values[["autoload_known_extensions"]], "false")
  threads <- dbGetQuery(
    con,
    "SELECT value FROM duckdb_settings() WHERE name = 'threads'"
  )$value[[1L]]
  expect_identical(threads, "1")

  loaded <- dbGetQuery(
    con,
    paste(
      "SELECT count(*) AS n FROM duckdb_functions()",
      "WHERE function_name = 'duckhts_htslib_version'"
    )
  )$n
  expect_equal(loaded, 1)

  json_state <- dbGetQuery(
    con,
    paste(
      "SELECT loaded FROM duckdb_extensions()",
      "WHERE extension_name = 'json'"
    )
  )
  if (nrow(json_state) == 1L) {
    expect_false(json_state$loaded[[1L]])
  }
}

test_reused_file_driver_rejected <- function() {
  dbdir <- tempfile("rduckhts_reused_", fileext = ".duckdb")
  drv <- .connection_test_driver(
    dbdir,
    config = list(allow_unsigned_extensions = "false")
  )
  con <- DBI::dbConnect(drv)
  on.exit(.connection_test_cleanup(con, drv, dbdir), add = TRUE)

  expect_error(
    rduckhts_connect(dbdir = dbdir),
    "already has a live instance"
  )
  expect_true(DBI::dbIsValid(con))
  expect_equal(DBI::dbGetQuery(con, "SELECT 42 AS answer")$answer[[1L]], 42)
}

test_failed_file_connection_releases_driver <- function() {
  dbdir <- tempfile("rduckhts_failed_", fileext = ".duckdb")
  missing_extension <- tempfile("rduckhts_missing_", fileext = ".duckdb_extension")
  expect_error(
    rduckhts_connect(dbdir = dbdir, extension_path = missing_extension),
    "extension not found"
  )

  drv <- .connection_test_driver(
    dbdir,
    config = list(allow_unsigned_extensions = "false")
  )
  con <- DBI::dbConnect(drv)
  on.exit(.connection_test_cleanup(con, drv, dbdir), add = TRUE)
  allow_unsigned <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT lower(value) AS value FROM duckdb_settings()",
      "WHERE name = 'allow_unsigned_extensions'"
    )
  )$value[[1L]]
  expect_identical(allow_unsigned, "false")
}

test_compatibility_connection <- function() {
  con <- duckhts_load()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  expect_equal(
    DBI::dbGetQuery(
      con,
      paste(
        "SELECT count(*) AS n FROM duckdb_functions()",
        "WHERE function_name = 'duckhts_htslib_version'"
      )
    )$n[[1L]],
    1
  )
}

test_package_owned_connection()
test_reused_file_driver_rejected()
test_failed_file_connection_releases_driver()
test_compatibility_connection()
