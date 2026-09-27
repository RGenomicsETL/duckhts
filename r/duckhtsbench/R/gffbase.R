#' Stage the pinned GFFBase Python package for conformance benchmarks.
#'
#' @param site_dir Python target directory, defaulting to its registry cache path.
#' @param python Path to Python.
#' @param artifact_id Registry ID for the pinned Python package.
#' @return The staged Python site directory, invisibly.
#' @export
duckhts_bench_stage_gffbase <- function(
    site_dir = duckhts_bench_artifact_path(artifact_id), python = Sys.which("python3"),
    artifact_id = "gffbase_010") {
  if (!nzchar(python)) stop("python3 is required to stage GFFBase", call. = FALSE)
  if (!artifact_id %in% c("gffbase_010", "gffbase_021")) {
    stop("unsupported GFFBase registry ID", call. = FALSE)
  }
  row <- duckhts_bench_registry()
  row <- row[row$id == artifact_id, , drop = FALSE]
  wheel_id <- paste0(artifact_id, "_linux_x86_64_wheel")
  if (nrow(row) != 1L || row$locator != paste0("artifact:", wheel_id)) {
    stop("GFFBase derived registry entry is missing or has an unregistered wheel", call. = FALSE)
  }
  wheel <- duckhts_bench_fetch(wheel_id)
  status <- system2(python, c("-c", shQuote("import duckdb, pyarrow")))
  if (status != 0L) stop("GFFBase requires preinstalled duckdb and pyarrow", call. = FALSE)
  dir.create(site_dir, recursive = TRUE, showWarnings = FALSE)
  status <- system2(
    python,
    c("-m", "pip", "install", "--no-deps", "--no-index", "--upgrade", "--force-reinstall", "--target", shQuote(site_dir), shQuote(wheel))
  )
  if (status != 0L) stop("could not install verified GFFBase wheel", call. = FALSE)
  version <- if (artifact_id == "gffbase_010") "0.1.0" else "0.2.1"
  writeLines(c("field\tvalue", paste0("workload\t", row$workload), "package\tgffbase",
    paste0("version\t", version), paste0("source\t", row$locator),
    paste0("site_directory\t", site_dir)), file.path(site_dir, "provenance.tsv"))
  invisible(site_dir)
}
