# Encode the same regional records in both formats, retaining all samples and
# all fields. Also used by the network-free staging test with a local source.
duckhts_bench_genotype_pair <- function(source, region, outputs, bcftools) {
  stopifnot(identical(names(outputs), c("vcf", "bcf")), nzchar(region), nzchar(bcftools))
  run <- function(args) {
    status <- system2(bcftools, shQuote(args))
    if (status != 0L) stop("could not stage genotype cohort", call. = FALSE)
  }
  for (output in outputs) dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary <- setNames(paste0(outputs, ".partial-", Sys.getpid()), names(outputs))
  on.exit(unlink(temporary), add = TRUE)
  run(c("view", "--no-version", "-r", region, "-Oz", "-o", temporary[["vcf"]], source))
  run(c("view", "--no-version", "-Ob", "-o", temporary[["bcf"]], temporary[["vcf"]]))
  query <- function(args) {
    result <- system2(bcftools, shQuote(args), stdout = TRUE)
    if (!is.null(attr(result, "status"))) stop("could not validate staged genotypes", call. = FALSE)
    result
  }
  samples <- query(c("query", "-l", source))
  records <- query(c("view", "-H", temporary[["vcf"]]))
  stopifnot(length(samples) > 0L, length(records) > 0L,
            identical(query(c("query", "-l", temporary[["vcf"]])), samples),
            identical(query(c("query", "-l", temporary[["bcf"]])), samples),
            identical(query(c("view", "-H", temporary[["bcf"]])), records))
  duckhts_bench_duckvep_publish_pair(temporary[["vcf"]], temporary[["bcf"]],
                                    outputs[["vcf"]], outputs[["bcf"]])
  invisible(list(records = length(records), samples = length(samples)))
}

#' Stage a Real Cohort for Genotype Reader Benchmarks
#'
#' Reads one registry-declared region of the immutable HPRC source through its
#' registered index. Unlike consequence-only corpora, this keeps every sample,
#' allele and GT/PS field. Network access occurs only in this explicit staging
#' step, never while rendering a benchmark or building the extension.
#' @return Named local VCF.gz and BCF paths.
#' @export
duckhts_bench_stage_genotypes <- function() {
  plan <- duckhts_bench_stage_plan("genotype-reader")
  stopifnot(identical(plan$id, c("geno_hprc_vcfgz", "geno_hprc_bcf")))
  rows <- duckhts_bench_registry()
  definition <- duckhts_bench_duckvep_corpus_definitions()[["hprc-african4-chr22"]]
  source <- duckhts_bench_duckvep_corpus_row(rows, definition$source)
  stopifnot(identical(plan$locator, c(
    paste0("artifact:", definition$source, ";artifact:", definition$source_index),
    "artifact:geno_hprc_vcfgz")), all(plan$release == source$release))
  identity <- duckhts_bench_identity_fields(plan$supplier_identity[[1]])
  stopifnot(identical(unname(identity[c("region", "all_samples", "genotypes_removed")]),
                      c("chr22:20000000-21000000", "true", "false")))
  index <- duckhts_bench_fetch(definition$source_index)
  duckhts_bench_duckvep_validate_source(rows, definition, index, NULL, Sys.which("curl"))
  outputs <- setNames(vapply(plan$id, duckhts_bench_artifact_path, character(1)), c("vcf", "bcf"))
  counts <- duckhts_bench_genotype_pair(paste0(source$locator, "##idx##", index),
                                       identity[["region"]], outputs, Sys.which("bcftools"))
  for (i in seq_along(outputs)) {
    duckhts_bench_validate_identity(plan$id[[i]], outputs[[i]])
    receipt <- duckhts_bench_write_provenance(plan$id[[i]], outputs[[i]])
    fields <- utils::read.delim(receipt)
    fields <- rbind(fields, data.frame(
      field = c("source_supplier_identity", "source_index_artifact", "bcftools_version",
                "observed_md5", "observed_bytes", "records", "samples"),
      value = c(source$supplier_identity, definition$source_index,
                system2("bcftools", "--version", stdout = TRUE)[[1]],
                unname(tools::md5sum(outputs[[i]])), file.info(outputs[[i]])$size,
                counts$records, counts$samples)))
    duckhts_bench_duckvep_atomic_table(fields, receipt)
  }
  outputs
}
