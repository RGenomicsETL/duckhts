#' Create a DuckDB connection with bundled DuckHTS loaded
#'
#' Creates a DuckDB connection that explicitly permits loading the package-built
#' DuckHTS extension, disables automatic installation and loading of unrelated
#' known DuckDB extensions, and loads the bundled DuckHTS extension.
#'
#' @details
#' Current versions of the \pkg{duckdb} R package can disable extension loading on
#' Linux builds that use a C++ standard library other than \code{libstdc++}.
#' DuckHTS is compiled locally with the same package toolchain, so this helper
#' explicitly enables extension loading for this package-owned connection.
#' Package-owned connections use per-session DuckDB extension/secret storage
#' by default when the installed \pkg{duckdb} version supports that setting.
#'
#' DuckDB reuses a live database instance when the same file-backed \code{dbdir}
#' is opened again, ignoring new driver and configuration settings. This helper
#' rejects such reuse rather than silently weakening its connection policy. Use
#' \code{rduckhts_load()} on the existing connection if it already permits
#' unsigned, driver-level extension loading, or close the existing database
#' instance before calling this helper.
#'
#' The settings \code{allow_unsigned_extensions},
#' \code{autoinstall_known_extensions}, and
#' \code{autoload_known_extensions} are controlled by this helper and override
#' entries with those names in \code{config}. Other named DuckDB settings are
#' passed through unchanged.
#'
#' @param dbdir Path to a DuckDB database, or \code{":memory:"} for an
#'   in-memory database.
#' @param read_only Logical; open a file-backed database read-only.
#' @param bigint How DuckDB 64-bit integers are returned; passed to
#'   \code{duckdb::duckdb()}.
#' @param config Named list of additional DuckDB configuration settings.
#' @param extension_path Optional path to a DuckHTS \code{.duckdb_extension}
#'   file. If \code{NULL}, uses the extension bundled with Rduckhts.
#'
#' @return A DuckDB connection with DuckHTS loaded.
#'
#' @examples
#' con <- rduckhts_connect()
#' DBI::dbGetQuery(con, "SELECT duckhts_htslib_version() AS version")
#' DBI::dbDisconnect(con, shutdown = TRUE)
#'
#' @export
rduckhts_connect <- function(
  dbdir = ":memory:",
  read_only = FALSE,
  bigint = "numeric",
  config = list(),
  extension_path = NULL
) {
  if (!is.list(config)) {
    stop("config must be a named list", call. = FALSE)
  }
  if (length(config) > 0L &&
      (is.null(names(config)) || anyNA(names(config)) || any(!nzchar(names(config))))) {
    stop("config must be a named list", call. = FALSE)
  }

  controlled_config <- list(
    allow_unsigned_extensions = "true",
    autoinstall_known_extensions = "false",
    autoload_known_extensions = "false"
  )
  config <- utils::modifyList(config, controlled_config)

  # duckdb() caches live file-backed drivers and ignores settings on reuse. A
  # private R attribute survives on a newly created driver's config slot but is
  # not a DuckDB setting, letting us reject a cached driver without reaching
  # into duckdb's private driver registry.
  driver_token <- new.env(parent = emptyenv())
  attr(config, "rduckhts_driver_token") <- driver_token

  driver_args <- list(
    dbdir = dbdir,
    read_only = read_only,
    bigint = bigint,
    config = config
  )
  duckdb_args <- names(formals(duckdb::duckdb))
  if ("allow_extensions" %in% duckdb_args) {
    driver_args$allow_extensions <- TRUE
  }
  if ("shared_home" %in% duckdb_args) {
    driver_args$shared_home <- FALSE
  }

  drv <- do.call(duckdb::duckdb, driver_args)
  driver_config <- methods::slot(drv, "config")
  owns_driver <- identical(
    attr(driver_config, "rduckhts_driver_token"),
    driver_token
  )
  if (!owns_driver) {
    stop(
      "DuckDB already has a live instance for dbdir; rduckhts_connect() ",
      "cannot apply its connection policy to a reused driver. Use ",
      "rduckhts_load() on an appropriately configured existing connection, ",
      "or close the existing database instance before retrying.",
      call. = FALSE
    )
  }

  con <- NULL
  loaded <- FALSE
  on.exit({
    if (!loaded) {
      if (!is.null(con)) {
        try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
      }
      try(duckdb::duckdb_shutdown(drv), silent = TRUE)
    }
  }, add = TRUE)
  con <- DBI::dbConnect(drv)
  rduckhts_load(con, extension_path = extension_path)
  loaded <- TRUE
  con
}

#' Load DuckHTS Extension
#'
#' Loads the DuckHTS extension into an existing DuckDB connection. This must be
#' called before using HTS reader functions on a connection not created by
#' \code{rduckhts_connect()}.
#'
#' @details
#' The connection must permit unsigned extension loading. With current versions
#' of the \pkg{duckdb} R package, its driver must also permit extension loading.
#' Prefer \code{rduckhts_connect()} when Rduckhts owns the connection.
#'
#' @param con A DuckDB connection object.
#' @param extension_path Optional path to the DuckHTS extension file. If
#'   \code{NULL}, uses the extension bundled with Rduckhts.
#'
#' @return \code{TRUE} if the extension was loaded successfully.
#'
#' @examples
#' con <- rduckhts_connect()
#' DBI::dbDisconnect(con, shutdown = TRUE)
#'
#' @export
rduckhts_load <- function(con, extension_path = NULL) {
  if (is.null(extension_path)) {
    # Try to use the bundled extension build directory
    ext_dir <- duckhts_extension_dir()
    if (is.null(ext_dir)) {
      stop(
        "duckhts_extension directory not found in installed package.",
        call. = FALSE
      )
    }
    extension_path <- file.path(ext_dir, "build", "duckhts.duckdb_extension")
  }

  if (!file.exists(extension_path)) {
    stop(
      "DuckHTS extension not found at: ",
      extension_path,
      "\nThis suggests the package was not built correctly during installation.",
      "\nTry recompiling the package with R CMD INSTALL --preclean Rduckhts"
    )
  }

  # Ensure unsigned extensions are allowed (must be set at connection creation)
  setting <- DBI::dbGetQuery(
    con,
    "SELECT value FROM duckdb_settings() WHERE name = 'allow_unsigned_extensions'"
  )
  if (nrow(setting) == 1 && tolower(setting$value[1]) != "true") {
    stop(
      "DuckDB connection must allow unsigned extensions. Recreate it with ",
      "rduckhts_connect(), or explicitly enable unsigned and driver-level ",
      "extension loading when constructing an external DuckDB connection.",
      call. = FALSE
    )
  }
  DBI::dbExecute(con, "SET enable_progress_bar = false")

  # Load the extension. dbQuoteString handles package library paths containing
  # quotes or other SQL-sensitive characters.
  quoted_extension <- as.character(DBI::dbQuoteString(con, extension_path))
  result <- DBI::dbExecute(con, paste("LOAD", quoted_extension))
  return(result == 0)
}
