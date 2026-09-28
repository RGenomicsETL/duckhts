# Record and verify the build identity of each builder measured by
# benchmark_ancestry_panel_run.R.
#
# Rscript benchmarks/benchmark_ancestry_panel_identity.R <out.tsv> <name>=<lib>=<commit>=<source_dir> ...
#
# <source_dir> is the Rduckhts source tree the library was built from. Two checks tie
# the installed library to <commit>:
# - every tracked code file of r/Rduckhts at <commit> is byte-identical in <source_dir>,
#   and <source_dir> holds no other files except build outputs (object files, the
#   shared library and the configure-generated htslib_config.R); documentation
#   (NEWS, README, man/) is not compared;
# - every function in the installed Rduckhts namespace equals the function its R
#   sources at <commit> define, compared without source references.
# Builds are not bit-reproducible across build directories, so hashes of the installed
# package code and extension are recorded to tie measured runs to these libraries.
args <- commandArgs(TRUE)

if (identical(args[1], "--functions")) {
  # Child process: compare the installed namespace with the recorded R sources.
  .libPaths(c(args[2], .libPaths()))
  ns <- asNamespace("Rduckhts")
  collate <- read.dcf(file.path(args[3], "DESCRIPTION"), fields = "Collate")[1, 1]
  files <- if (is.na(collate)) sort(list.files(file.path(args[3], "R"), "\\.[Rr]$")) else
    scan(text = gsub("\n", " ", collate), what = "", quiet = TRUE)
  env <- new.env(parent = ns)
  for (f in files) sys.source(file.path(args[3], "R", f), envir = env, keep.source = FALSE)
  names <- Filter(function(n) is.function(env[[n]]), ls(env, all.names = TRUE))
  differ <- Filter(function(n) {
    !exists(n, envir = ns, inherits = FALSE) ||
      !identical(utils::removeSource(env[[n]]), utils::removeSource(get(n, envir = ns)),
                 ignore.environment = TRUE)
  }, names)
  if (length(differ)) stop("installed functions differ from the recorded source: ",
                           paste(head(differ, 10), collapse = ", "))
  cat(length(names), "\n")
  quit(save = "no")
}

stopifnot(length(args) >= 2L)
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
git <- function(...) system2("git", c(...), stdout = TRUE)
build_outputs <- "^(src/.*\\.(o|so|dll)|inst/htslib_config\\.R)$"
documentation <- "^(NEWS\\.md|README\\.(Rmd|md)|man/.*)$"
rows <- lapply(strsplit(args[-1], "=", fixed = TRUE), function(p) {
  stopifnot(length(p) == 4L)
  commit <- git("rev-parse", "--verify", paste0(p[3], "^{commit}"))
  recorded <- tempfile("rduckhts-source-")
  dir.create(recorded)
  on.exit(unlink(recorded, recursive = TRUE), add = TRUE)
  system2("sh", c("-c", shQuote(sprintf("git archive %s r/Rduckhts | tar -x -C %s",
                                        commit, shQuote(recorded)))))
  recorded <- file.path(recorded, "r", "Rduckhts")
  tracked <- list.files(recorded, recursive = TRUE, all.files = TRUE)
  present <- list.files(p[4], recursive = TRUE, all.files = TRUE)
  code <- tracked[!grepl(documentation, tracked)]
  missing <- setdiff(code, present)
  extra <- setdiff(present, tracked)
  extra <- extra[!grepl(build_outputs, extra)]
  changed <- code[tools::md5sum(file.path(recorded, code)) != tools::md5sum(file.path(p[4], code))]
  if (length(c(missing, extra, changed))) {
    stop(p[1], ": build source differs from ", commit, ": ",
         paste(head(c(missing, extra, changed), 10), collapse = ", "))
  }
  functions <- system2("Rscript", shQuote(c(script, "--functions", p[2], recorded)),
                       stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(functions, "status"))) stop(p[1], ": ", paste(functions, collapse = "\n"))
  pkg <- file.path(p[2], "Rduckhts")
  extension <- list.files(file.path(pkg, "duckhts_extension"), pattern = "\\.duckdb_extension$",
                          recursive = TRUE, full.names = TRUE)
  stopifnot(length(extension) == 1L)
  description <- read.dcf(file.path(pkg, "DESCRIPTION"), fields = c("Version", "Packaged"))
  data.frame(implementation = p[1], commit = commit,
             builder_blob = git("rev-parse", paste0(commit, ":r/Rduckhts/R/ancestry_panel.R")),
             verified_code_files = length(code),
             verified_functions = as.integer(trimws(tail(functions, 1))),
             package_version = description[1, "Version"],
             package_code_sha256 = digest::digest(file = file.path(pkg, "R", "Rduckhts.rdb"),
                                                  algo = "sha256"),
             extension_sha256 = digest::digest(file = extension, algo = "sha256"),
             packaged = description[1, "Packaged"])
})
utils::write.table(do.call(rbind, rows), args[1], sep = "\t", quote = FALSE, row.names = FALSE)
