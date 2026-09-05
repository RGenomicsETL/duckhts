#' Stage the deterministic indexed multi-region reader workload, without network.
#'
#' @return Named paths for the registered VCF, BGZF/TBI, and BCF/CSI artifacts.
#' @export
duckhts_bench_stage_multi_region <- function() {
  plan <- duckhts_bench_stage_plan("multi-region-readers")
  ids <- paste0("multi_region_", c("vcf", "vcf_gz", "tbi", "bcf", "csi"))
  stopifnot(identical(plan$id, ids))
  paths <- stats::setNames(vapply(ids, duckhts_bench_artifact_path, character(1)),
                          c("vcf", "vcf_gz", "tbi", "bcf", "csi"))
  stopifnot(paths[["tbi"]] == paste0(paths[["vcf_gz"]], ".tbi"),
            paths[["csi"]] == paste0(paths[["bcf"]], ".csi"))
  rows <- as.integer(duckhts_bench_identity_fields(plan$supplier_identity[[1]])[["records"]])
  stopifnot(length(rows) == 1L, !is.na(rows), rows > 0L)
  for (tool in c("bgzip", "tabix", "bcftools")) {
    if (!nzchar(Sys.which(tool))) stop(tool, " is required for multi-region staging", call. = FALSE)
  }
  for (path in paths) dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  out <- file(paths[["vcf"]], open = "wt")
  tryCatch({
    writeLines(c("##fileformat=VCFv4.3", sprintf("##contig=<ID=1,length=%d>", rows),
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), out)
    for (first in seq.int(1L, rows, by = 100000L)) {
      pos <- seq.int(first, min(rows, first + 99999L))
      writeLines(sprintf("1\t%d\t.\tA\tC\t.\tPASS\t.", pos), out)
    }
  }, finally = close(out))
  stopifnot(system2("bgzip", c("-@", "1", "-c", shQuote(paths[["vcf"]])),
                    stdout = paths[["vcf_gz"]]) == 0L)
  stopifnot(system2("tabix", c("-f", "-p", "vcf", shQuote(paths[["vcf_gz"]]))) == 0L)
  stopifnot(system2("bcftools", c("view", "--no-version", "-Ob", "-o",
                    shQuote(paths[["bcf"]]), shQuote(paths[["vcf_gz"]]))) == 0L)
  stopifnot(system2("bcftools", c("index", "-f", shQuote(paths[["bcf"]]))) == 0L)
  for (i in seq_along(ids)) duckhts_bench_write_provenance(ids[[i]], paths[[i]])
  invisible(paths)
}
