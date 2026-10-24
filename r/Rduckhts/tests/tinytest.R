if (requireNamespace("tinytest", quietly = TRUE)) {
  # NOT_CRAN=true (CI, devtools, rcmdcheck) also runs the long tests that exit early
  # on CRAN's machines, where the whole check must stay under 10 minutes.
  tinytest::test_package(
    "Rduckhts",
    at_home = identical(Sys.getenv("NOT_CRAN"), "true") ||
      identical(Sys.getenv("TT_AT_HOME"), "TRUE")
  )
} else {
  message("tinytest not available; skipping tests")
}
