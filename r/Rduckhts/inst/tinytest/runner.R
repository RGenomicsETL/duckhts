if (requireNamespace("tinytest", quietly = TRUE)) {
  # HTSlib writes directly to process stderr and does not frame concurrent logs.
  diagnostics_file <- tempfile("rduckhts-tinytest-stderr-")
  expression <- "tinytest::test_package('Rduckhts')"
  status <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", "-e", shQuote(expression)),
    stdout = "",
    stderr = diagnostics_file
  ))
  diagnostics <- readLines(diagnostics_file, warn = FALSE)
  unlink(diagnostics_file)

  if (!identical(status, 0L)) {
    if (length(diagnostics) > 0L) writeLines(diagnostics, stderr())
    stop("tinytest subprocess failed with status ", status, call. = FALSE)
  }
} else {
  message("tinytest not available; skipping tests")
}
