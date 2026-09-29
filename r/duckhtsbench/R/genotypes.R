duckhts_bench_publish_files <- function(temporaries, targets) {
  if (!length(temporaries) || length(temporaries) != length(targets) ||
      anyDuplicated(targets) || !all(file.exists(temporaries))) {
    stop("publication requires distinct targets and one existing temporary per target",
      call. = FALSE)
  }
  backups <- paste0(targets, ".backup-", Sys.getpid())
  unlink(backups, force = TRUE)
  backed_up <- logical(length(targets))
  for (index in seq_along(targets)) {
    if (file.exists(targets[[index]])) {
      if (!file.rename(targets[[index]], backups[[index]])) {
        for (restore in which(backed_up)) file.rename(backups[[restore]], targets[[restore]])
        stop("could not preserve existing artifact: ", targets[[index]], call. = FALSE)
      }
      backed_up[[index]] <- TRUE
    }
  }
  published <- logical(length(targets))
  for (index in seq_along(targets)) {
    if (!file.rename(temporaries[[index]], targets[[index]])) {
      for (remove in which(published)) unlink(targets[[remove]], force = TRUE)
      for (restore in which(backed_up)) file.rename(backups[[restore]], targets[[restore]])
      stop("could not publish complete artifact set: ", targets[[index]], call. = FALSE)
    }
    published[[index]] <- TRUE
  }
  unlink(backups, force = TRUE)
  invisible(targets)
}

duckhts_bench_write_genotype_receipt <- function(fields, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".partial-", Sys.getpid())
  backup <- paste0(path, ".backup-", Sys.getpid())
  unlink(c(temporary, backup), force = TRUE)
  utils::write.table(fields, temporary, sep = "\t", row.names = FALSE, quote = FALSE)
  if (file.exists(path) && !file.rename(path, backup)) {
    unlink(temporary, force = TRUE)
    stop("could not preserve existing receipt before publication: ", path, call. = FALSE)
  }
  if (!file.rename(temporary, path)) {
    if (file.exists(backup)) file.rename(backup, path)
    unlink(temporary, force = TRUE)
    stop("could not publish receipt: ", path, call. = FALSE)
  }
  unlink(backup, force = TRUE)
  invisible(path)
}

# Encode the same regional records in both formats, retaining all samples and
# all fields. Also used by the network-free staging test with a local source.
duckhts_bench_genotype_pair <- function(source, region, outputs, bcftools) {
  stopifnot(identical(names(outputs), c("vcf", "bcf")), nzchar(region), nzchar(bcftools))
  run <- function(args) {
    status <- system2(bcftools, shQuote(args))
    if (status != 0L) stop("could not stage genotype cohort", call. = FALSE)
  }
  for (output in outputs) dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary <- stats::setNames(paste0(outputs, ".partial-", Sys.getpid()), names(outputs))
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
  duckhts_bench_publish_files(temporary, outputs)
  invisible(list(records = length(records), samples = length(samples)))
}

duckhts_bench_genotype_staging_paths <- function(outputs, label) {
  directory <- tempfile(paste0(label, "-"), tmpdir = dirname(outputs[[1L]]))
  if (!dir.create(directory, recursive = TRUE)) {
    stop("could not create genotype staging directory", call. = FALSE)
  }
  stats::setNames(file.path(directory, basename(outputs)), names(outputs))
}

duckhts_bench_publish_genotype_bundle <- function(staged, outputs) {
  temporary <- c(staged, paste0(staged, ".provenance.tsv"))
  targets <- c(outputs, paste0(outputs, ".provenance.tsv"))
  duckhts_bench_publish_files(temporary, targets)
}

duckhts_bench_genotype_phase_observation <- function(path, bcftools) {
  query <- function(args) {
    result <- system2(bcftools, shQuote(args), stdout = TRUE)
    if (!is.null(attr(result, "status"))) {
      stop("could not inspect staged phase-set genotypes", call. = FALSE)
    }
    result
  }
  header <- query(c("view", "-h", path))
  phase_header <- header[startsWith(header, "##FORMAT=<ID=PS,")]
  if (length(phase_header) != 1L ||
      !grepl("Number=1", phase_header, fixed = TRUE) ||
      !grepl("Type=Integer", phase_header, fixed = TRUE)) {
    stop("phase-set benchmark requires FORMAT/PS Number=1,Type=Integer", call. = FALSE)
  }
  samples <- query(c("query", "-l", path))
  records <- query(c("query", "-f", "%CHROM\\n", path))
  calls <- query(c("query", "-f", "[%GT\\t%PS\\n]", path))
  regions <- unique(records)
  if (!length(samples) || length(regions) != 1L ||
      length(calls) != length(records) * length(samples)) {
    stop("phase-set benchmark record/sample cardinality is inconsistent", call. = FALSE)
  }
  fields <- strsplit(calls, "\t", fixed = TRUE)
  if (any(lengths(fields) != 2L)) {
    stop("phase-set benchmark GT/PS projection is malformed", call. = FALSE)
  }
  gt <- vapply(fields, `[[`, character(1L), 1L)
  ps <- vapply(fields, `[[`, character(1L), 2L)
  list(
    region = regions, ps_type = "Integer", samples = samples, calls = calls,
    counts = c(records = length(records), samples = length(samples), calls = length(calls),
      allele_slots = sum(nchar(gsub("[^|/]", "", gt)) + 1L),
      nonnull_ps = sum(ps != "."))
  )
}

duckhts_bench_genotype_phase_counts <- function(path, bcftools) {
  duckhts_bench_genotype_phase_observation(path, bcftools)$counts
}

duckhts_bench_genotype_phase_set_sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256")
}

duckhts_bench_read_genotype_phase_set_receipt <- function(path) {
  receipt <- utils::read.delim(path, colClasses = "character", check.names = FALSE,
                               stringsAsFactors = FALSE)
  if (!identical(names(receipt), c("field", "value")) || !nrow(receipt) ||
      anyNA(receipt$field) || anyNA(receipt$value) ||
      any(!nzchar(receipt$field)) || any(!nzchar(receipt$value)) ||
      anyDuplicated(receipt$field)) {
    stop("phase-set provenance receipt is malformed", call. = FALSE)
  }
  stats::setNames(receipt$value, receipt$field)
}

duckhts_bench_validate_genotype_phase_set_evidence <- function(
    registry, ids, source_ids, paths, receipts, observations) {
  formats <- c("VCF", "BCF")
  semantic_fields <- c("region", "all_samples", "genotypes_removed", "ps_type",
                       "records", "samples", "calls", "allele_slots", "nonnull_ps")
  count_fields <- c("records", "samples", "calls", "allele_slots", "nonnull_ps")
  registry_fields <- c("workload", "release", "locator", "access", "transform",
                       "supplier_identity", "consumer")
  if (!identical(names(ids), formats) || !identical(names(paths), formats) ||
      !identical(names(receipts), formats) || !identical(names(observations), formats) ||
      !identical(names(source_ids), c("source", "index"))) {
    stop("phase-set evidence inputs must use the declared names", call. = FALSE)
  }
  required_ids <- c(unname(source_ids), unname(ids))
  if (anyDuplicated(required_ids) ||
      any(vapply(required_ids, function(id) sum(registry$id == id) != 1L, logical(1L)))) {
    stop("phase-set registry closure is incomplete or non-unique", call. = FALSE)
  }
  rows <- lapply(required_ids, function(id) registry[registry$id == id, , drop = FALSE])
  names(rows) <- required_ids
  if (!all(registry_fields %in% names(registry))) {
    stop("phase-set registry lacks provenance fields", call. = FALSE)
  }
  derived <- lapply(ids, function(id) rows[[id]])
  source <- rows[[source_ids[["source"]]]]
  source_index <- rows[[source_ids[["index"]]]]
  closure <- c(list(source, source_index), derived)
  if (!all(vapply(closure, function(row)
      identical(as.character(row$workload), "genotype-phase-set"), logical(1L))) ||
      !identical(as.character(derived[["VCF"]]$locator),
                 paste0("artifact:", source_ids[["source"]], ";artifact:",
                        source_ids[["index"]])) ||
      !identical(as.character(derived[["BCF"]]$locator),
                 paste0("artifact:", ids[["VCF"]])) ||
      length(unique(vapply(closure, function(row)
        as.character(row$release), character(1L)))) != 1L) {
    stop("phase-set registry dependency or release contract is inconsistent", call. = FALSE)
  }

  identities <- lapply(derived, function(row)
    duckhts_bench_identity_fields(row$supplier_identity))
  if (!all(vapply(identities, function(identity)
      setequal(names(identity), semantic_fields), logical(1L)))) {
    stop("phase-set registry identity has an incomplete semantic contract", call. = FALSE)
  }
  identities <- lapply(identities, function(identity) identity[semantic_fields])
  if (!identical(identities[["VCF"]], identities[["BCF"]]) ||
      !identical(unname(identities[["VCF"]][c("all_samples", "genotypes_removed")]),
                 c("true", "false"))) {
    stop("phase-set registry identities disagree", call. = FALSE)
  }
  expected_identity <- identities[["VCF"]]
  expected_counts <- suppressWarnings(as.numeric(expected_identity[count_fields]))
  if (anyNA(expected_counts) || any(!is.finite(expected_counts)) ||
      any(expected_counts < 0) || any(expected_counts != floor(expected_counts))) {
    stop("phase-set registry denominators are invalid", call. = FALSE)
  }
  names(expected_counts) <- count_fields

  for (format in formats) {
    row <- derived[[format]]
    receipt <- receipts[[format]]
    observation <- observations[[format]]
    required_receipt <- c("artifact_id", "workload", "release", "source_locator",
                          "access", "transform", "supplier_identity", "cached_output",
                          "consumer", "source_artifact", "source_supplier_identity",
                          "source_index_artifact", "source_index_supplier_identity",
                          "bcftools_version", "observed_sha256", count_fields)
    if (!is.character(receipt) || is.null(names(receipt)) ||
        !all(required_receipt %in% names(receipt)) ||
        any(vapply(receipt[required_receipt], length, integer(1L)) != 1L)) {
      stop("phase-set provenance receipt lacks required fields", call. = FALSE)
    }
    expected_receipt <- c(
      artifact_id = ids[[format]], workload = as.character(row$workload),
      release = as.character(row$release), source_locator = as.character(row$locator),
      access = as.character(row$access), transform = as.character(row$transform),
      supplier_identity = as.character(row$supplier_identity), cached_output = paths[[format]],
      consumer = as.character(row$consumer), source_artifact = source_ids[["source"]],
      source_supplier_identity = as.character(source$supplier_identity),
      source_index_artifact = source_ids[["index"]],
      source_index_supplier_identity = as.character(source_index$supplier_identity)
    )
    if (!identical(unname(receipt[names(expected_receipt)]), unname(expected_receipt)) ||
        !identical(receipt[["observed_sha256"]],
                   duckhts_bench_genotype_phase_set_sha256(paths[[format]])) ||
        length(receipt[["bcftools_version"]]) != 1L ||
        !nzchar(receipt[["bcftools_version"]])) {
      stop("phase-set provenance does not match the registry or artifact", call. = FALSE)
    }
    if (!is.list(observation) ||
        !all(c("region", "ps_type", "counts") %in% names(observation)) ||
        length(observation$region) != 1L || length(observation$ps_type) != 1L ||
        !identical(names(observation$counts), count_fields)) {
      stop("phase-set observation is incomplete", call. = FALSE)
    }
    if (!identical(as.character(observation$region), expected_identity[["region"]]) ||
        !identical(as.character(observation$ps_type), expected_identity[["ps_type"]]) ||
        !identical(as.numeric(observation$counts), unname(expected_counts)) ||
        !identical(unname(receipt[count_fields]), unname(expected_identity[count_fields]))) {
      stop("phase-set observation or receipt disagrees with the registry", call. = FALSE)
    }
  }
  if (!identical(receipts[["VCF"]][["bcftools_version"]],
                 receipts[["BCF"]][["bcftools_version"]])) {
    stop("phase-set encodings were not produced by the same bcftools version", call. = FALSE)
  }
  identity <- do.call(rbind, lapply(formats, function(format) data.frame(
    format = format, artifact = ids[[format]], bytes = file.info(paths[[format]])$size,
    sha256 = receipts[[format]][["observed_sha256"]],
    source = receipts[[format]][["source_artifact"]], stringsAsFactors = FALSE)))
  rownames(identity) <- NULL
  identity
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
  source_id <- "geno_hprc_source"
  index_id <- "geno_hprc_source_tbi"
  source <- rows[rows$id == source_id, , drop = FALSE]
  source_index <- rows[rows$id == index_id, , drop = FALSE]
  stopifnot(nrow(source) == 1L, nrow(source_index) == 1L,
            identical(plan$locator, c(
              paste0("artifact:", source_id, ";artifact:", index_id),
              "artifact:geno_hprc_vcfgz")), all(plan$release == source$release))
  identity <- duckhts_bench_identity_fields(plan$supplier_identity[[1]])
  stopifnot(identical(unname(identity[c("region", "all_samples", "genotypes_removed")]),
                      c("chr22:20000000-21000000", "true", "false")))
  index <- duckhts_bench_fetch(index_id)
  source_identity <- duckhts_bench_identity_fields(source$supplier_identity)
  index_identity <- duckhts_bench_identity_fields(source_index$supplier_identity)
  if (!file.exists(index) || !file.info(index)$size ||
      !all(c("version_id", "bytes") %in% names(source_identity)) ||
      !nzchar(source_identity[["version_id"]]) ||
      !all(c("version_id", "sha256", "bytes") %in% names(index_identity)) ||
      !grepl(paste0("versionId=", source_identity[["version_id"]]),
             source$locator, fixed = TRUE) ||
      !grepl(paste0("versionId=", index_identity[["version_id"]]),
             source_index$locator, fixed = TRUE)) {
    stop("HPRC source and index must retain their registered S3 version IDs",
         call. = FALSE)
  }
  outputs <- stats::setNames(
    vapply(plan$id, duckhts_bench_artifact_path, character(1)), c("vcf", "bcf")
  )
  staged_outputs <- duckhts_bench_genotype_staging_paths(outputs, "genotype-reader")
  on.exit(unlink(dirname(staged_outputs[[1L]]), recursive = TRUE), add = TRUE)
  counts <- duckhts_bench_genotype_pair(paste0(source$locator, "##idx##", index),
                                       identity[["region"]], staged_outputs, Sys.which("bcftools"))
  for (i in seq_along(outputs)) {
    duckhts_bench_validate_identity(plan$id[[i]], staged_outputs[[i]])
    receipt <- duckhts_bench_write_provenance(plan$id[[i]], staged_outputs[[i]])
    fields <- utils::read.delim(receipt)
    fields$value[fields$field == "cached_output"] <- outputs[[i]]
    fields <- rbind(fields, data.frame(
      field = c("source_supplier_identity", "source_index_artifact", "bcftools_version",
                "observed_md5", "observed_bytes", "records", "samples"),
      value = c(source$supplier_identity, index_id,
                system2("bcftools", "--version", stdout = TRUE)[[1]],
                unname(tools::md5sum(staged_outputs[[i]])), file.info(staged_outputs[[i]])$size,
                counts$records, counts$samples)))
    duckhts_bench_write_genotype_receipt(fields, receipt)
  }
  duckhts_bench_publish_genotype_bundle(staged_outputs, outputs)
  outputs
}

#' Stage a Real Integer Phase-Set Benchmark
#'
#' Downloads the pinned GIAB HG002 v4.2.1 phased benchmark and its index, then
#' derives complete chr1 VCF.gz and BCF inputs without rewriting GT or PS.
#' Network access occurs only in this explicit staging step.
#' @param bcftools Path to the bcftools executable.
#' @return Named local VCF.gz and BCF paths.
#' @export
duckhts_bench_stage_genotype_phase_set <- function(bcftools = Sys.which("bcftools")) {
  if (length(bcftools) != 1L || is.na(bcftools) || !nzchar(bcftools)) {
    stop("bcftools is required to stage the phase-set benchmark", call. = FALSE)
  }
  plan <- duckhts_bench_stage_plan("genotype-phase-set")
  expected_ids <- c("geno_giab_phased_source", "geno_giab_phased_source_tbi",
                    "geno_giab_phased_chr1_vcfgz", "geno_giab_phased_chr1_bcf")
  if (!identical(plan$id, expected_ids)) {
    stop("genotype phase-set registry plan is incomplete", call. = FALSE)
  }
  source <- duckhts_bench_fetch(expected_ids[[1L]])
  source_index <- duckhts_bench_fetch(expected_ids[[2L]])
  if (!identical(source_index, paste0(source, ".tbi"))) {
    stop("phase-set source index must be adjacent to its VCF.gz", call. = FALSE)
  }
  outputs <- stats::setNames(
    vapply(expected_ids[3:4], duckhts_bench_artifact_path, character(1L)),
    c("vcf", "bcf")
  )
  staged_outputs <- duckhts_bench_genotype_staging_paths(outputs, "genotype-phase-set")
  on.exit(unlink(dirname(staged_outputs[[1L]]), recursive = TRUE), add = TRUE)
  identities <- lapply(plan$supplier_identity[3:4], duckhts_bench_identity_fields)
  required <- c("region", "all_samples", "genotypes_removed", "ps_type",
                "records", "samples", "calls", "allele_slots", "nonnull_ps")
  if (!all(vapply(identities, function(identity) all(required %in% names(identity)),
                  logical(1L)))) {
    stop("phase-set derived artifact lacks required workload identity fields", call. = FALSE)
  }
  if (!identical(identities[[1L]][required], identities[[2L]][required])) {
    stop("phase-set VCF.gz and BCF registry identities differ", call. = FALSE)
  }
  identity <- identities[[1L]]
  staged <- duckhts_bench_genotype_pair(source, identity[["region"]], staged_outputs, bcftools)
  counts <- duckhts_bench_genotype_phase_counts(staged_outputs[["vcf"]], bcftools)
  expected <- as.numeric(identity[c("records", "samples", "calls", "allele_slots", "nonnull_ps")])
  names(expected) <- names(counts)
  if (!identical(names(counts), names(expected)) || any(counts != expected) ||
      staged$records != counts[["records"]] || staged$samples != counts[["samples"]]) {
    stop("staged phase-set denominators do not match the registry", call. = FALSE)
  }
  bcf_counts <- duckhts_bench_genotype_phase_counts(staged_outputs[["bcf"]], bcftools)
  if (!identical(counts, bcf_counts)) {
    stop("staged phase-set VCF.gz and BCF denominators differ", call. = FALSE)
  }
  source_identity <- plan$supplier_identity[[1L]]
  source_index_identity <- plan$supplier_identity[[2L]]
  version <- system2(bcftools, "--version", stdout = TRUE)[[1L]]
  for (i in seq_along(outputs)) {
    duckhts_bench_validate_identity(expected_ids[[i + 2L]], staged_outputs[[i]])
    receipt <- duckhts_bench_write_provenance(expected_ids[[i + 2L]], staged_outputs[[i]])
    fields <- utils::read.delim(receipt, colClasses = "character", check.names = FALSE)
    fields$value[fields$field == "cached_output"] <- outputs[[i]]
    fields <- rbind(fields, data.frame(
      field = c("source_artifact", "source_supplier_identity", "source_index_artifact",
                "source_index_supplier_identity", "bcftools_version", "observed_sha256",
                names(counts)),
      value = c(expected_ids[[1L]], source_identity, expected_ids[[2L]],
                source_index_identity, version,
                digest::digest(file = staged_outputs[[i]], algo = "sha256"), counts),
      stringsAsFactors = FALSE
    ))
    duckhts_bench_write_genotype_receipt(fields, receipt)
  }
  duckhts_bench_publish_genotype_bundle(staged_outputs, outputs)
  outputs
}
