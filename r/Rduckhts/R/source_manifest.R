# The bundled package manifest is the installed copy of src/duckhts_sources.tsv.
# Keep source inventory consumers on this one receipt instead of maintaining
# another R-only list.

duckhts_source_manifest <- function() {
  manifest <- file.path(duckhts_extension_dir(), "duckhts_sources.tsv")
  if (!file.exists(manifest)) {
    stop("DuckHTS source manifest not found: ", manifest, call. = FALSE)
  }
  utils::read.delim(
    manifest,
    sep = "|",
    header = FALSE,
    col.names = c("repo_path", "package_path"),
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
}

duckhts_duckvep_kernel_source_files <- function() {
  paths <- duckhts_source_manifest()$package_path
  prefix <- "duckvep/kernel/src/"
  paths <- paths[startsWith(paths, prefix) & grepl("[.]c$", paths)]
  sub(prefix, "", paths, fixed = TRUE)
}
