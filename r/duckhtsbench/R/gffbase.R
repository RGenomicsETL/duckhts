#' Stage the pinned GFFBase Python package for conformance benchmarks.
#'
#' @param site_dir Python target directory, defaulting to its registry cache path.
#' @param python Path to Python.
#' @return The staged Python site directory, invisibly.
#' @export
duckhts_bench_stage_gffbase <- function(
    site_dir = duckhts_bench_artifact_path("gffbase_010"), python = Sys.which("python3")) {
  if (!nzchar(python)) stop("python3 is required to stage GFFBase", call. = FALSE)
  row <- duckhts_bench_registry()
  row <- row[row$id == "gffbase_010", , drop = FALSE]
  if (nrow(row) != 1L || row$locator != "artifact:gffbase_010_linux_x86_64_wheel") {
    stop("GFFBase derived registry entry is missing or has an unregistered wheel", call. = FALSE)
  }
  wheel <- duckhts_bench_fetch("gffbase_010_linux_x86_64_wheel")
  status <- system2(python, c("-c", shQuote("import duckdb, pyarrow")))
  if (status != 0L) stop("GFFBase requires preinstalled duckdb and pyarrow", call. = FALSE)
  dir.create(site_dir, recursive = TRUE, showWarnings = FALSE)
  status <- system2(
    python,
    c("-m", "pip", "install", "--no-deps", "--no-index", "--upgrade", "--force-reinstall", "--target", shQuote(site_dir), shQuote(wheel))
  )
  if (status != 0L) stop("could not install verified GFFBase wheel", call. = FALSE)
  writeLines(c("field\tvalue", "workload\tgffbase-conformance", "package\tgffbase", "version\t0.1.0", paste0("source\t", row$locator), paste0("site_directory\t", site_dir)), file.path(site_dir, "provenance.tsv"))
  invisible(site_dir)
}
