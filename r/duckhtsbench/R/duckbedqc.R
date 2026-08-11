#' Stage DuckBedQC BED providers for the real cgranges benchmark.
#'
#' @param destination Checkout directory, defaulting to the registered cache path.
#' @param git Path to `git`.
#' @return The detached checkout directory, invisibly.
#' @export
duckhts_bench_stage_duckbedqc <- function(
    destination = duckhts_bench_artifact_path("duckbedqc_118fc21"), git = Sys.which("git")) {
  if (!nzchar(git)) stop("git is required to stage DuckBedQC", call. = FALSE)
  row <- duckhts_bench_registry()
  row <- row[row$id == "duckbedqc_118fc21", , drop = FALSE]
  if (nrow(row) != 1L) stop("DuckBedQC registry entry is missing", call. = FALSE)
  commit <- sub("^git_commit=", "", row$supplier_identity)
  run <- function(args, error) {
    status <- system2(git, args)
    if (status != 0L) stop(error, call. = FALSE)
  }
  if (!dir.exists(file.path(destination, ".git"))) {
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    run(c("clone", row$locator, destination), "could not clone registered DuckBedQC repository")
  }
  run(c("-C", destination, "fetch", "--quiet", "origin", commit), "could not fetch registered DuckBedQC commit")
  run(c("-C", destination, "checkout", "--quiet", "--force", "--detach", commit), "could not check out registered DuckBedQC commit")
  actual <- trimws(system2(git, c("-C", destination, "rev-parse", "HEAD"), stdout = TRUE))
  if (!identical(actual, commit)) stop("DuckBedQC checkout does not match registered commit", call. = FALSE)
  required <- file.path(destination, "data", c("GRCh38_exons.bed", "GRCh38_illumina_clinical_regions_v100.39.0.bed"))
  if (!all(file.exists(required))) stop("registered DuckBedQC commit lacks required cgranges providers", call. = FALSE)
  writeLines(c("field\tvalue", "workload\tcgranges-real-benchmark", paste0("repository\t", row$locator), paste0("commit\t", commit), paste0("checkout\t", destination)), file.path(destination, "provenance.tsv"))
  invisible(destination)
}
