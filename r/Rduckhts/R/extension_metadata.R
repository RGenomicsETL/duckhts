# Canonical DuckDB extension metadata writer used by package configure and
# the retained command-line adapter in tools/append_extension_metadata.R.

duckhts_append_extension_metadata <- function(
    library_file,
    out_file,
    extension_name = "duckhts",
    duckdb_platform,
    duckdb_version,
    extension_version,
    abi_type = "C_STRUCT"
) {
  required <- list(
    library_file = library_file,
    out_file = out_file,
    extension_name = extension_name,
    duckdb_platform = duckdb_platform,
    duckdb_version = duckdb_version,
    extension_version = extension_version,
    abi_type = abi_type
  )
  missing <- names(required)[vapply(required, function(x) {
    length(x) != 1L || is.na(x) || !nzchar(as.character(x))
  }, logical(1))]
  if (length(missing)) {
    stop(
      "Missing extension metadata value(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!file.exists(library_file)) {
    stop("Extension library not found: ", library_file, call. = FALSE)
  }

  pad32 <- function(value, name) {
    bytes <- charToRaw(as.character(value))
    if (length(bytes) > 32L) {
      stop(name, " exceeds the 32-byte extension metadata field", call. = FALSE)
    }
    c(bytes, raw(32L - length(bytes)))
  }

  out_tmp <- paste0(out_file, ".tmp")
  unlink(out_tmp)
  on.exit(unlink(out_tmp), add = TRUE)
  if (!file.copy(library_file, out_tmp, overwrite = TRUE)) {
    stop("Could not stage extension library: ", library_file, call. = FALSE)
  }

  con <- file(out_tmp, open = "ab")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  writeBin(as.raw(c(0x00, 0x93, 0x04, 0x10)), con, useBytes = TRUE)
  writeBin(charToRaw("duckdb_signature"), con, useBytes = TRUE)
  writeBin(as.raw(c(0x80, 0x04)), con, useBytes = TRUE)
  writeBin(pad32("", "field 8"), con, useBytes = TRUE)
  writeBin(pad32("", "field 7"), con, useBytes = TRUE)
  writeBin(pad32("", "field 6"), con, useBytes = TRUE)
  writeBin(pad32(abi_type, "ABI type"), con, useBytes = TRUE)
  writeBin(pad32(extension_version, "extension version"), con, useBytes = TRUE)
  writeBin(pad32(duckdb_version, "DuckDB API version"), con, useBytes = TRUE)
  writeBin(pad32(duckdb_platform, "DuckDB platform"), con, useBytes = TRUE)
  writeBin(pad32("4", "metadata magic"), con, useBytes = TRUE)
  writeBin(raw(256L), con, useBytes = TRUE)
  close(con)

  if (file.exists(out_file)) {
    unlink(out_file)
    if (file.exists(out_file)) {
      stop("Could not replace extension library: ", out_file, call. = FALSE)
    }
  }
  if (!file.rename(out_tmp, out_file)) {
    stop("Could not publish extension library: ", out_file, call. = FALSE)
  }
  invisible(out_file)
}
