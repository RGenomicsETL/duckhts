# Identity checks of the staged inputs of the read-count validation.
#
# The cache can hold a stale, damaged or replaced file under the right name.
# Before any decode, each input is compared with the identity its registry row
# records, so the report cannot certify results from other inputs. The
# identities are those staging computes: roh_counts_identity() for the read
# counts and roh_validation_genotype_identity() for the genotype BCF.
#
# Run from the repository root.

source("benchmarks/roh_counts_stage.R")
source("benchmarks/roh_counts_validation/stage_genotypes.R")

roh_validation_registered_identity <- function(id) {
  registry <- duckhtsbench::duckhts_bench_registry()
  row <- registry[registry$id == id, , drop = FALSE]
  if (nrow(row) != 1L) stop("unknown or non-unique benchmark artifact: ", id, call. = FALSE)
  duckhtsbench:::duckhts_bench_identity_fields(row$supplier_identity[[1L]])
}

roh_validation_require_fields <- function(id, registered, required) {
  absent <- setdiff(required, names(registered))
  if (length(absent)) {
    stop("registered identity of ", id, " lacks: ", paste(absent, collapse = ", "), call. = FALSE)
  }
}

# The read counts of one sample: its registered sample name, row count and
# canonical counts digest.
roh_validation_check_counts <- function(con, id, sample) {
  registered <- roh_validation_registered_identity(id)
  roh_validation_require_fields(id, registered, c("sample", "rows", "counts_sha256"))
  observed <- roh_counts_identity(con, duckhtsbench::duckhts_bench_artifact_path(id))
  if (!identical(registered[["sample"]], sample) ||
      observed$rows != as.numeric(registered[["rows"]]) ||
      !identical(observed$counts_sha256, registered[["counts_sha256"]])) {
    stop("read counts of ", sample, " differ from the registered identity of ", id, ": ",
         observed$rows, " rows, counts SHA-256 ", observed$counts_sha256, call. = FALSE)
  }
  invisible(observed)
}

# The genotype BCF: its registered samples in order, record count and record
# digest. `samples` are the samples the run decodes.
roh_validation_check_genotypes <- function(id, samples, bcftools = "/usr/local/bin/bcftools") {
  registered <- roh_validation_registered_identity(id)
  roh_validation_require_fields(id, registered,
                                c("sample_list", "records", "samples", "records_sha256"))
  observed <- roh_validation_genotype_identity(duckhtsbench::duckhts_bench_artifact_path(id),
                                               bcftools)
  registered_samples <- strsplit(registered[["sample_list"]], ",", fixed = TRUE)[[1L]]
  if (!identical(registered_samples, samples) ||
      !identical(observed$sample_list, registered_samples) ||
      observed$records != as.integer(registered[["records"]]) ||
      observed$samples != as.integer(registered[["samples"]]) ||
      !identical(observed$records_sha256, registered[["records_sha256"]])) {
    stop("genotype BCF differs from the registered identity of ", id, ": ",
         observed$records, " records, samples ", paste(observed$sample_list, collapse = ","),
         ", record SHA-256 ", observed$records_sha256, call. = FALSE)
  }
  invisible(observed)
}

# Every input of the declared validation.
roh_validation_check_inputs <- function(con, declaration,
                                        bcftools = "/usr/local/bin/bcftools") {
  for (sample in declaration$samples) {
    roh_validation_check_counts(con, declaration$count_artifacts[[sample]], sample)
  }
  roh_validation_check_genotypes(declaration$genotype_artifact, declaration$samples, bcftools)
  invisible(TRUE)
}
