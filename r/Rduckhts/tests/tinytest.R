if (requireNamespace("tinytest", quietly = TRUE)) {
    tinytest::test_package("Rduckhts")
} else {
    message("tinytest not available; skipping tests")
}
