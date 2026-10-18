# Run the two BAM/CRAM test functions with a flushed marker before each statement.
# This diagnostic supplements the full package suite; it does not replace it.
library(Rduckhts)
library(DBI)
library(tinytest)

print(sessionInfo())
print(Sys.getenv(c(
  "_R_CHECK_ELAPSED_TIMEOUT_", "_R_CHECK_TESTS_ELAPSED_TIMEOUT_",
  "_R_CHECK_ONE_TEST_ELAPSED_TIMEOUT_", "_R_WIN_CHECK_INVALID_PARAMETERS_"
)))

targets <- c(
  test_ancestry.R = "test_ancestry_bam_cram",
  test_somalier_extraction.R = "test_somalier_bam_count_extraction"
)

for (file in names(targets)) {
  path <- system.file("tinytest", file, package = "Rduckhts", mustWork = TRUE)
  definitions <- Filter(function(expr) {
    is.call(expr) && identical(expr[[1L]], as.name("<-"))
  }, as.list(parse(path)))
  name <- targets[[file]]
  index <- which(vapply(definitions, function(expr) {
    identical(expr[[2L]], as.name(name))
  }, logical(1)))
  stopifnot(length(index) == 1L)
  test_function <- eval(definitions[[index]][[3L]])
  statements <- as.list(body(test_function))[-1L]
  traced <- list(as.name("{"))
  for (i in seq_along(statements)) {
    statement <- statements[[i]]
    if (is.call(statement) && identical(statement[[1L]], as.name("on.exit"))) {
      cleanup <- statement[[2L]]
      if (is.call(cleanup) && identical(cleanup[[1L]], as.name("dbDisconnect"))) {
        statement[[2L]] <- substitute({
          cat(NAME, "starting dbDisconnect\n")
          flush(stdout())
          DISCONNECT
          cat(NAME, "finished dbDisconnect\n")
          flush(stdout())
        }, list(NAME = name, DISCONNECT = cleanup))
      }
    }
    label <- paste(name, i, deparse(statement)[1L])
    marker <- substitute({
      cat(format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"), LABEL, "\n")
      flush(stdout())
    }, list(LABEL = label))
    traced <- c(traced, list(marker, statement))
  }
  body(test_function) <- as.call(traced)
  definitions[[index]] <- substitute(NAME <- FUN, list(
    NAME = as.name(name), FUN = test_function
  ))
  code <- c(definitions, list(call(name)))
  probe <- tempfile(pattern = "bam-extraction-", fileext = ".R")
  writeLines(unlist(lapply(code, deparse)), probe)
  results <- tinytest::run_test_file(probe, verbose = 0)
  unlink(probe)
  print(results)
  stopifnot(length(results) > 0L, all(vapply(results, isTRUE, logical(1))))
  cat(name, "COMPLETE (including on.exit cleanup)\n")
  flush(stdout())
}
