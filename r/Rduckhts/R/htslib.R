#' Setup HTSlib Environment
#'
#' Sets the `HTS_PATH` environment variable to point to the bundled htslib
#' plugins directory. This enables remote file access via libcurl plugins
#' (e.g., s3://, gs://, http://) when plugins are available.
#'
#' @param plugins_dir Optional path to the htslib plugins directory. When NULL,
#'   uses the bundled plugins directory if available.
#'
#' @return Invisibly returns the previous value of `HTS_PATH` (or `NA` if unset).
#'
#' @details
#' Call this before the process opens its first HTS file. htslib discovers
#' dynamic plugins on first file access and does not rescan `HTS_PATH` later.
#'
#' @examples
#' \dontrun{
#' setup_hts_env()
#'
#' plugins_path <- tempfile("hts_plugins_")
#' dir.create(plugins_path)
#' setup_hts_env(plugins_dir = plugins_path)
#' unlink(plugins_path, recursive = TRUE)
#' }
#' @export
setup_hts_env <- function(plugins_dir = NULL) {
  old_value <- Sys.getenv("HTS_PATH", unset = NA)
  if (is.null(plugins_dir)) {
    plugins_dir <- duckhts_htslib_plugins_dir()
  }
  if (nzchar(plugins_dir) && dir.exists(plugins_dir)) {
    Sys.setenv(HTS_PATH = plugins_dir)
  }
  invisible(old_value)
}

#' @keywords internal
duckhts_htslib_plugins_dir <- function() {
  ext_dir <- duckhts_extension_dir()
  candidates <- c(
    file.path(ext_dir, "htslib", "libexec", "htslib"),
    file.path(ext_dir, "htslib", "plugins"),
    file.path(ext_dir, "htslib")
  )
  for (p in candidates) {
    if (dir.exists(p)) {
      return(normalizePath(p))
    }
  }
  ""
}

#' Inspect the Loaded htslib Build
#'
#' Returns the version and build features reported by the htslib library that
#' is actually loaded with DuckHTS.
#'
#' @param con A DuckDB connection with DuckHTS loaded.
#'
#' @return A one-row data frame with `version`, `feature_bits`, and
#'   `feature_string` columns.
#'
#' @examples
#' \dontrun{
#' con <- rduckhts_connect()
#' rduckhts_htslib_info(con)
#' DBI::dbDisconnect(con, shutdown = TRUE)
#' }
#' @export
rduckhts_htslib_info <- function(con) {
  DBI::dbGetQuery(
    con,
    paste(
      "SELECT duckhts_htslib_version() AS version,",
      "duckhts_htslib_features() AS feature_bits,",
      "duckhts_htslib_feature_string() AS feature_string"
    )
  )
}

#' Return the Loaded htslib Version
#'
#' @param con A DuckDB connection with DuckHTS loaded.
#'
#' @return The runtime htslib semantic version as a character scalar.
#'
#' @export
rduckhts_htslib_version <- function(con) {
  as.character(rduckhts_htslib_info(con)$version[[1L]])
}

duckhts_runtime_htslib_info <- function() {
  con <- rduckhts_connect()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  rduckhts_htslib_info(con)
}

#' Get the Installed htslib Linking Contract
#'
#' Resolves headers, an exact shared or static library, linker flags, enabled
#' features, and build identity from the installed Rduckhts package. With
#' validation enabled, the receipt and public headers are also compared with
#' the htslib version reported by the loaded DuckHTS extension.
#'
#' @param link Either `"shared"` or `"static"`. When omitted, use the link
#'   mode selected when this Rduckhts package was configured.
#' @param validate Whether to validate installed files, header identity, and
#'   the loaded htslib runtime version.
#'
#' @return An object of class `rduckhts_htslib_config`. Its `cppflags` and
#'   `ldflags` elements can be consumed by a downstream package configure
#'   script; `duckdb_platform` records the extension footer platform and the
#'   remaining fields form the versioned build receipt.
#'
#' @details
#' The shared contract is currently available on native Unix builds. MinGW and
#' browser-wasm builds expose the static contract unless a shared htslib was
#' explicitly built. Static consumers must review `static_license_note`.
#'
#' @examples
#' \dontrun{
#' config <- rduckhts_htslib_config()
#' config$cppflags
#' config$ldflags
#' config$features
#' }
#' @export
rduckhts_htslib_config <- function(
  link = NULL,
  validate = TRUE
) {
  config_path <- system.file("htslib_config.R", package = "Rduckhts")
  if (!nzchar(config_path) || !file.exists(config_path)) {
    stop("Rduckhts installed htslib contract was not found", call. = FALSE)
  }

  contract <- new.env(parent = baseenv())
  sys.source(config_path, envir = contract)
  if (is.null(link)) {
    link <- contract$htslib_default_link()
  }
  link <- match.arg(link, c("shared", "static"))
  config <- contract$htslib_config(link = link, validate = validate)

  if (isTRUE(validate)) {
    runtime <- duckhts_runtime_htslib_info()
    runtime_version <- as.character(runtime$version[[1L]])
    if (!identical(runtime_version, config$htslib_version)) {
      stop(
        "Rduckhts htslib receipt/runtime version mismatch: receipt=",
        config$htslib_version,
        ", runtime=",
        runtime_version,
        call. = FALSE
      )
    }
    config$runtime_version <- runtime_version
    config$runtime_feature_bits <- runtime$feature_bits[[1L]]
    config$runtime_feature_string <- as.character(runtime$feature_string[[1L]])
  }

  config
}
