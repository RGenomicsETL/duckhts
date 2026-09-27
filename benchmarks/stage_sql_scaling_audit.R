# Run from the repository root; source acquisition is explicit with --download.
source("r/duckhtsbench/R/registry.R")
source("r/duckhtsbench/R/stage.R")
Sys.setenv(DUCKHTSBENCH_REGISTRY = normalizePath(
  "r/duckhtsbench/inst/benchmark_registry.tsv", mustWork = TRUE))

download <- identical(commandArgs(trailingOnly = TRUE), "--download")
stopifnot(download || length(commandArgs(trailingOnly = TRUE)) == 0L)
inputs <- c("geno_giab_phased_source", "geno_giab_phased_source_tbi",
            "sql_scaling_1000g_phase3_chr22", "sql_scaling_epilepsy")
for (id in inputs) {
  path <- duckhts_bench_artifact_path(id)
  if (!file.exists(path) && !download) {
    stop("staged input missing (pass --download to fetch): ", id)
  }
  if (download) duckhts_bench_fetch(id)
  duckhts_bench_validate_identity(id)
}

references <- c("liftover_grch37_fasta", "liftover_grch38_fasta",
                "liftover_grch37_grch38_chain")
missing_references <- !vapply(references, function(id) {
  file.exists(duckhts_bench_artifact_path(id))
}, logical(1))
if (any(missing_references)) {
  if (!download) stop("stage the registered liftover reference bundle")
  source("r/duckhtsbench/R/liftover.R")
  duckhts_bench_stage_liftover()
}
for (id in references) duckhts_bench_validate_identity(id)

bcftools <- Sys.which("bcftools")
if (!nzchar(bcftools)) stop("bcftools is required to stage indexed VCF subsets")
phase3 <- duckhts_bench_artifact_path("sql_scaling_1000g_phase3_chr22")
if (!file.exists(paste0(phase3, ".tbi"))) {
  status <- system2(bcftools, c("index", "-t", shQuote(phase3)))
  if (status != 0L) stop("could not index phase-3 source")
}

giab <- duckhts_bench_artifact_path("geno_giab_phased_source")
for (multiplier in c(1L, 2L, 4L)) {
  for (kind in c("giab", "phase3")) {
    id <- sprintf("sql_scaling_%s_%dx", kind, multiplier)
    output <- duckhts_bench_artifact_path(id)
    dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
    candidate <- if (file.exists(output)) output else {
      paste0(output, ".partial-", Sys.getpid(), ".vcf.gz")
    }
    tryCatch({
      if (candidate != output) {
        if (kind == "giab") {
          regions <- paste0("chr", seq_len(multiplier), collapse = ",")
          status <- system2(bcftools, c("view", "--threads", "2", "-r", regions,
                                        "-Oz", "-o", shQuote(candidate), shQuote(giab)))
        } else {
          end <- 16000000L + multiplier * 5000000L
          command <- paste("bcftools view -r", shQuote(sprintf("22:16000000-%d", end)),
                           "-m2 -M2 -v snps -Ou", shQuote(phase3),
                           "| bcftools norm -d all -Oz -o", shQuote(candidate))
          status <- system2("bash", c("-o", "pipefail", "-c", shQuote(command)))
        }
        if (status != 0L) stop("could not derive ", id)
      }
      index <- paste0(candidate, ".tbi")
      if (!file.exists(index)) {
        status <- system2(bcftools, c("index", "-t", shQuote(candidate)))
        if (status != 0L) stop("could not index ", id)
      }
      record_count <- as.integer(system2(bcftools, c("index", "-n", shQuote(candidate)),
                                         stdout = TRUE))
      row <- duckhts_bench_registry()
      row <- row[row$id == id, , drop = FALSE]
      expected <- as.integer(duckhts_bench_identity_fields(row$supplier_identity)[["records"]])
      if (!identical(record_count, expected)) stop("record count mismatch: ", id)
      if (candidate != output) {
        if (!file.rename(candidate, output) || !file.rename(index, paste0(output, ".tbi"))) {
          stop("could not publish ", id)
        }
      }
      duckhts_bench_write_provenance(id, output)
      message(id, ": ", record_count, " records")
    }, finally = {
      if (candidate != output) unlink(c(candidate, paste0(candidate, ".tbi")))
    })
  }
}
