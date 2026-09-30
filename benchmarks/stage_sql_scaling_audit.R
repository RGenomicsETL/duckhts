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

duplicate_id <- "sql_scaling_vcf_counts_duplicate_million"
duplicate_source_id <- "sql_scaling_giab_1x"
duplicate_source <- duckhts_bench_artifact_path(duplicate_source_id)
duplicate_output <- duckhts_bench_artifact_path(duplicate_id)
if (!file.exists(paste0(duplicate_source, ".tbi"))) {
  stop("indexed 1x GIAB input is required to stage the duplicate-error case")
}
header <- system2(bcftools, c("view", "-h", shQuote(duplicate_source)), stdout = TRUE)
if (!any(grepl("^##FORMAT=<ID=AD,", header))) {
  stop("registered 1x GIAB input does not declare FORMAT/AD")
}
seed_rows <- system2(
  bcftools,
  c("view", "-H", "-r", "chr1:1-2000000", "-m2", "-M2", "-v", "snps",
    shQuote(duplicate_source)),
  stdout = TRUE)
seed_fields <- strsplit(seed_rows, "\t", fixed = TRUE)
valid_seed <- vapply(seed_fields, function(fields) {
  if (length(fields) < 10L || fields[[1L]] != "chr1" ||
      !fields[[4L]] %in% c("A", "C", "G", "T") ||
      !fields[[5L]] %in% c("A", "C", "G", "T") ||
      fields[[4L]] == fields[[5L]]) {
    return(FALSE)
  }
  format_fields <- strsplit(fields[[9L]], ":", fixed = TRUE)[[1L]]
  sample_fields <- strsplit(fields[[10L]], ":", fixed = TRUE)[[1L]]
  format_indices <- match(c("GT", "AD"), format_fields)
  !anyNA(format_indices) && length(sample_fields) >= max(format_indices) &&
    !any(sample_fields[format_indices] %in% c("", ".")) &&
    grepl("^[0-9]+,[0-9]+$", sample_fields[[format_indices[[2L]]]])
}, logical(1))
if (!any(valid_seed)) stop("could not find a biallelic chr1 SNP with FORMAT/AD")
seed_fields <- seed_fields[[which(valid_seed)[[1L]]]]
format_fields <- strsplit(seed_fields[[9L]], ":", fixed = TRUE)[[1L]]
sample_fields <- strsplit(seed_fields[[10L]], ":", fixed = TRUE)[[1L]]
format_indices <- match(c("GT", "AD"), format_fields)
if (anyNA(format_indices) || length(sample_fields) < max(format_indices) ||
    any(sample_fields[format_indices] %in% c("", ".")) ||
    !grepl("^[0-9]+,[0-9]+$", sample_fields[[format_indices[[2L]]]])) {
  stop("selected chr1 SNP lacks a usable GT:AD sample call")
}
seed_record <- paste(c(seed_fields[[1L]], seed_fields[[2L]], ".",
                       seed_fields[[4L]], seed_fields[[5L]], ".", "PASS", ".",
                       "GT:AD", paste(sample_fields[format_indices], collapse = ":")),
                     collapse = "\t")
registry <- duckhts_bench_registry()
duplicate_row <- registry[registry$id == duplicate_id, , drop = FALSE]
expected_records <- as.integer(duckhts_bench_identity_fields(
  duplicate_row$supplier_identity[[1L]])[["records"]])
dir.create(dirname(duplicate_output), recursive = TRUE, showWarnings = FALSE)
temporary_output <- paste0(duplicate_output, ".partial-", Sys.getpid())
tryCatch({
  duckhts_bench_write_repeated_vcf(header, seed_record, temporary_output,
                                   expected_records)
  record_count <- as.integer(system2(
    "grep", c("-vc", shQuote("^#"), shQuote(temporary_output)), stdout = TRUE))
  if (!identical(record_count, expected_records)) {
    stop("duplicate-error VCF record count mismatch")
  }
  if (file.exists(duplicate_output)) unlink(duplicate_output)
  if (!file.rename(temporary_output, duplicate_output)) {
    stop("could not publish duplicate-error VCF")
  }
  duckhts_bench_write_provenance(duplicate_id, duplicate_output)
  message(duplicate_id, ": ", record_count, " duplicate records")
}, finally = unlink(temporary_output))
