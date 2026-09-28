# Record the build identity of each builder measured by benchmark_ancestry_panel_run.R.
#
# Rscript benchmarks/benchmark_ancestry_panel_identity.R <out.tsv> <name>=<lib>=<commit>=<source_dir> ...
#
# <source_dir> is the Rduckhts source tree the library was built from. The builder
# source in it must equal r/Rduckhts/R/ancestry_panel.R at <commit> (compared as git
# blobs); the installed package code and extension binary are hashed so a rerun can be
# tied to the same builds.
args <- commandArgs(TRUE)
stopifnot(length(args) >= 2L)
git <- function(...) system2("git", c(...), stdout = TRUE)
rows <- lapply(strsplit(args[-1], "=", fixed = TRUE), function(p) {
  stopifnot(length(p) == 4L)
  commit <- git("rev-parse", "--verify", paste0(p[3], "^{commit}"))
  blob <- git("rev-parse", paste0(commit, ":r/Rduckhts/R/ancestry_panel.R"))
  built <- git("hash-object", file.path(p[4], "R", "ancestry_panel.R"))
  if (!identical(blob, built)) stop(p[1], ": builder source differs from ", commit)
  pkg <- file.path(p[2], "Rduckhts")
  extension <- list.files(file.path(pkg, "duckhts_extension"), pattern = "\\.duckdb_extension$",
                          recursive = TRUE, full.names = TRUE)
  stopifnot(length(extension) == 1L)
  description <- read.dcf(file.path(pkg, "DESCRIPTION"), fields = c("Version", "Packaged"))
  data.frame(implementation = p[1], commit = commit, builder_blob = blob,
             package_version = description[1, "Version"],
             package_code_sha256 = digest::digest(file = file.path(pkg, "R", "Rduckhts.rdb"), algo = "sha256"),
             extension_sha256 = digest::digest(file = extension, algo = "sha256"),
             packaged = description[1, "Packaged"])
})
utils::write.table(do.call(rbind, rows), args[1], sep = "\t", quote = FALSE, row.names = FALSE)
