args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: stage_remote_chromosome.R CHROMOSOME", call. = FALSE)
chromosome <- as.integer(args[[1L]])
source("benchmarks/roh_ancestry_stage.R")
registry <- duckhtsbench::duckhts_bench_registry()
source_id <- sprintf("roh_ancestry_chr%d_source", chromosome)
url <- registry$locator[registry$id == source_id]
if (length(url) != 1L || !startsWith(url, "https://")) {
  stop("registry has no remote source for ", source_id, call. = FALSE)
}
headers <- system2("curl", shQuote(c("-sIL", "--max-time", "90", url)),
                   stdout = TRUE, stderr = TRUE)
if (!is.null(attr(headers, "status"))) stop("curl failed to retrieve source metadata")
header_value <- function(name) {
  matches <- grep(paste0("^", name, ":"), headers, value = TRUE, ignore.case = TRUE)
  if (!length(matches)) stop("source response lacks ", name)
  trimws(sub("^[^:]+:", "", tail(matches, 1L)))
}
last_modified <- header_value("Last-Modified")
etag <- header_value("ETag")
output <- duckhtsbench::duckhts_bench_artifact_path(sprintf("roh_ancestry_chr%d_children_bcf", chromosome))
# The registry-aware entry point checks the derived BCF's registered record and
# sample counts before publishing it.
result <- stage_roh_children_from_registry(chromosome, registry,
  duckhtsbench::duckhts_bench_artifact_path("roh_ancestry_pedigree"), output,
  bcftools = Sys.getenv("BCFTOOLS", unname(Sys.which("bcftools"))))
receipt <- data.frame(
  artifact = c("source_vcf", "children_bcf"),
  locator = c(url, sprintf("derived:chr%d_source+pedigree", chromosome)),
  remote_last_modified = c(last_modified, ""), remote_etag = c(etag, ""),
  samples = c(result$source_samples, result$samples),
  records = c(NA_integer_, result$records),
  sha256 = c("", result$output_sha256),
  status = c("remote metadata recorded; streamed by bcftools",
             sprintf("bcftools 1.23.1-70-g6dbd8fef; biallelic SNVs; INFO/AF fields retained")),
  stringsAsFactors = FALSE)
result_dir <- file.path("benchmarks/results/roh-ancestry", sprintf("chr%d", chromosome))
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.table(receipt, file.path(result_dir, "input_receipt.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE, na = "")
cat("source", url, "\n")
cat("source_etag", etag, "last_modified", last_modified, "\n")
cat("chromosome", result$chromosome, "samples", result$samples,
    "records", result$records, "output_sha256", result$output_sha256, "\n")
