# Build and verify an in-tree extension for a reproducible benchmark.

benchmark_command <- function(command, args, context) {
  output <- suppressWarnings(system2(command, shQuote(args), stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(context, ":\n", paste(output, collapse = "\n"), call. = FALSE)
  }
  output
}

benchmark_evidence_revision <- function(root) {
  revision <- trimws(benchmark_command(
    "git", c("-C", root, "rev-parse", "HEAD"), "cannot identify source revision"
  )[[1L]])
  if (!grepl("^[0-9a-f]{40}$", revision)) {
    stop("source revision is not a full Git object name", call. = FALSE)
  }
  revision
}

benchmark_extension_path <- function(root, extension) {
  expected <- normalizePath(file.path(root, "build", "release", "duckhts.duckdb_extension"),
                            mustWork = FALSE)
  if (!identical(normalizePath(extension, mustWork = FALSE), expected)) {
    stop("benchmark must use the in-tree release extension", call. = FALSE)
  }
  expected
}

benchmark_evidence_sha256 <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  before <- file.info(path)[c("size", "mtime", "ctime")]
  output <- benchmark_command("sha256sum", path, "cannot hash extension")
  hash <- strsplit(trimws(output[[1L]]), "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", hash) ||
      !identical(before, file.info(path)[c("size", "mtime", "ctime")])) {
    stop("extension changed while hashing", call. = FALSE)
  }
  hash
}

benchmark_assert_checkout <- function(root, revision, allowed_outputs = character()) {
  if (!identical(benchmark_evidence_revision(root), revision)) {
    stop("benchmark source revision changed", call. = FALSE)
  }
  changed <- benchmark_command("git", c("-C", root, "status", "--porcelain=v1",
                                      "--untracked-files=no"), "cannot inspect checkout")
  root <- normalizePath(root, mustWork = TRUE)
  allowed <- normalizePath(allowed_outputs, mustWork = FALSE)
  allowed <- substring(allowed[startsWith(allowed, paste0(root, "/"))],
                       nchar(root) + 2L)
  if (length(setdiff(sub("^.. ", "", changed), allowed)) > 0L) {
    stop("benchmark needs a clean tracked checkout", call. = FALSE)
  }
  inputs <- benchmark_command(
    "git", c("-C", root, "ls-files", "--others", "--", "Makefile", "CMakeLists.txt",
             "description.yml", "src", "third_party/htslib"),
    "cannot inspect untracked build inputs"
  )
  inputs <- inputs[grepl("\\.(c|cc|cpp|h|hpp|inc|def|o|cmake)$|(^|/)(Makefile|CMakeLists[.]txt)$",
                         inputs)]
  inputs <- setdiff(inputs, c("third_party/htslib/config.h",
                              "third_party/htslib/config_vars.h",
                              "third_party/htslib/version.h"))
  if (length(inputs) > 0L) {
    stop("untracked benchmark build inputs: ", paste(inputs, collapse = ", "),
         call. = FALSE)
  }
}

benchmark_evidence_build_extension <- function(root, extension, revision,
                                               allowed_outputs = character()) {
  root <- normalizePath(root, mustWork = TRUE)
  extension <- benchmark_extension_path(root, extension)
  benchmark_assert_checkout(root, revision, allowed_outputs)
  htslib <- file.path(root, "third_party", "htslib")
  if (!file.exists(file.path(htslib, "Makefile")) ||
      !file.exists(file.path(root, "configure", "venv", "bin", "python3"))) {
    stop("benchmark requires a configured DuckHTS checkout", call. = FALSE)
  }
  if (length(benchmark_command("git", c("-C", file.path(root, "extension-ci-tools"),
                                      "status", "--porcelain=v1", "--untracked-files=all"),
                               "cannot inspect build tools")) > 0L) {
    stop("build tools checkout is not clean", call. = FALSE)
  }
  benchmark_command("make", c("-C", htslib, "-f", "Makefile", "distclean"),
                    "cannot distclean vendored htslib")
  benchmark_command("make", c("-C", root, "-f", "Makefile", "clean"),
                    "cannot clean extension")
  benchmark_command("make", c("-C", root, "-f", "Makefile", "-j2", "release"),
                    "cannot build extension")
  benchmark_assert_checkout(root, revision, allowed_outputs)
  if (!file.exists(extension)) stop("release extension was not built", call. = FALSE)
  list(path = extension, binding = "htslib_distclean_make_release",
       sha256 = benchmark_evidence_sha256(extension))
}

benchmark_evidence_write_extension_receipt <- function(path, revision, extension) {
  if (!grepl("^[0-9a-f]{40}$", revision) ||
      !is.list(extension) ||
      !identical(extension$binding, "htslib_distclean_make_release") ||
      !identical(extension$sha256, benchmark_evidence_sha256(extension$path))) {
    stop("extension receipt requires a verified release build", call. = FALSE)
  }
  value <- data.frame(
    field = c("source_revision", "path", "binding", "sha256"),
    value = c(revision, normalizePath(extension$path, mustWork = TRUE),
              extension$binding, extension$sha256)
  )
  path <- normalizePath(path.expand(path), mustWork = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("extension-receipt-", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  utils::write.table(value, temporary, sep = "\t", row.names = FALSE,
                     quote = FALSE, na = "")
  if (!file.rename(temporary, path)) {
    stop("cannot publish extension receipt: ", path, call. = FALSE)
  }
  path
}

benchmark_evidence_read_extension_receipt <- function(path, root, extension, revision) {
  value <- utils::read.delim(path, colClasses = "character", quote = "",
                             comment.char = "", check.names = FALSE)
  fields <- c("source_revision", "path", "binding", "sha256")
  if (!identical(names(value), c("field", "value")) ||
      !identical(value$field, fields)) {
    stop("extension receipt has an invalid schema", call. = FALSE)
  }
  receipt <- stats::setNames(value$value, value$field)
  expected <- benchmark_extension_path(root, extension)
  if (!identical(receipt[["source_revision"]], revision) ||
      !identical(normalizePath(receipt[["path"]], mustWork = FALSE), expected) ||
      !identical(receipt[["binding"]], "htslib_distclean_make_release") ||
      !identical(receipt[["sha256"]], benchmark_evidence_sha256(expected))) {
    stop("extension receipt does not match the source and binary", call. = FALSE)
  }
  list(path = expected, binding = receipt[["binding"]], sha256 = receipt[["sha256"]])
}
