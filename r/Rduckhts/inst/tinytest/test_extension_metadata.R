library(tinytest)

# DuckDB refuses an extension whose footer names a C API version newer than its
# own. The footer of the bundled extension must therefore name the version the
# bundled sources were compiled for, which the bundled header defines.

extension_dir <- Rduckhts:::duckhts_extension_dir()
header <- file.path(extension_dir, "duckdb_capi", "duckdb_extension.h")
extension_file <- file.path(extension_dir, "build", "duckhts.duckdb_extension")
expect_true(file.exists(header))
expect_true(file.exists(extension_file))

header_version <- Rduckhts:::duckhts_capi_version(header)
# The declared baseline. Raising it also raises the oldest DuckDB that loads the package.
expect_identical(header_version, "v1.2.0")

# Read the footer by DuckDB's layout, not through the package writer: eight
# 32-byte fields, then a 256-byte signature. The API version is the third
# field from the end.
footer_field <- function(path, fields_from_end) {
  size <- file.size(path)
  con <- file(path, open = "rb")
  on.exit(close(con))
  seek(con, where = size - 256L - 32L * fields_from_end)
  bytes <- readBin(con, what = "raw", n = 32L)
  rawToChar(bytes[bytes != as.raw(0L)])
}
expect_identical(footer_field(extension_file, 1L), "4")
expect_identical(footer_field(extension_file, 3L), header_version)
expect_identical(footer_field(extension_file, 5L), "C_STRUCT")

# The header guards its defaults with `#if !defined(...)` lines that name the
# same macros. Only the `#define` lines carry the version.
write_header <- function(lines) {
  path <- tempfile(fileext = ".h")
  writeLines(lines, path)
  path
}
guarded <- write_header(c(
  "#if !defined(DUCKDB_EXTENSION_API_VERSION_MAJOR) && !defined(DUCKDB_EXTENSION_API_VERSION_MINOR) && \\",
  "    !defined(DUCKDB_EXTENSION_API_VERSION_PATCH)",
  "#define DUCKDB_EXTENSION_API_VERSION_MAJOR 1",
  "#define DUCKDB_EXTENSION_API_VERSION_MINOR 5",
  "#define DUCKDB_EXTENSION_API_VERSION_PATCH 6",
  "#endif"
))
expect_identical(Rduckhts:::duckhts_capi_version(guarded), "v1.5.6")

missing_patch <- write_header(c(
  "#define DUCKDB_EXTENSION_API_VERSION_MAJOR 1",
  "#define DUCKDB_EXTENSION_API_VERSION_MINOR 5"
))
expect_error(
  Rduckhts:::duckhts_capi_version(missing_patch),
  "Expected one DUCKDB_EXTENSION_API_VERSION_PATCH definition"
)
expect_error(
  Rduckhts:::duckhts_capi_version(file.path(tempdir(), "no_such_header.h")),
  "DuckDB C API header not found"
)
unlink(c(guarded, missing_patch))
