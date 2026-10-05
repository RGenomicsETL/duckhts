# Stage the genotypes of a few samples at the sites of an AF-site Parquet.
#
# The output is an indexed BCF with FORMAT/GT only, the named samples in the
# given order, and exactly one record per AF site: the source record with the
# site's chrom, pos, REF and ALT. Record IDs are rewritten to
# chrom_pos_ref_alt, which is how the records are matched to the sites.
#
# The identity is the record count, the sample count and the SHA-256 of the
# record lines (`bcftools view -H`). The file's SHA-256 is not compared:
# bcftools writes its command line and date into the header. The record digest
# is of bcftools text output, so it holds for the bcftools version the registry
# records. When an expected identity is given, a mismatch publishes nothing.

roh_validation_site_ids <- function(af_sites) {
  con <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  DBI::dbGetQuery(con, sprintf(
    "SELECT chrom, pos::BIGINT AS pos, ref, alt FROM read_parquet(%s) ORDER BY chrom, pos",
    DBI::dbQuoteString(con, af_sites)))
}

# The bcftools executable of staging and of the validation run: the BCFTOOLS
# environment variable when it is set, otherwise /usr/local/bin/bcftools.
roh_validation_bcftools <- function() {
  Sys.getenv("BCFTOOLS", "/usr/local/bin/bcftools")
}

# The identity of a genotype BCF as it is on disk: its record count, its
# sample list in header order, and the SHA-256 of its records as
# `bcftools view -H` writes them. Staging registers this identity and the
# validation run checks it before decoding.
roh_validation_genotype_identity <- function(path, bcftools = roh_validation_bcftools()) {
  if (length(path) != 1L || !file.exists(path) ||
      length(bcftools) != 1L || !file.exists(bcftools)) {
    stop("a genotype BCF and bcftools must be supplied", call. = FALSE)
  }
  run <- function(command) {
    value <- system2("/bin/bash", c("-o", "pipefail", "-c", shQuote(command)), stdout = TRUE)
    if (!is.null(attr(value, "status"))) {
      stop("bcftools could not read the genotype BCF: ", path, call. = FALSE)
    }
    value
  }
  records <- paste(shQuote(bcftools), "view -H", shQuote(path))
  digest <- run(paste(records, "| sha256sum"))
  count <- run(paste(records, "| wc -l"))
  if (length(digest) != 1L || length(count) != 1L) {
    stop("could not digest the genotype BCF records", call. = FALSE)
  }
  sample_list <- run(paste(shQuote(bcftools), "query -l", shQuote(path)))
  list(records = as.integer(count), samples = length(sample_list),
       sample_list = sample_list, records_sha256 = sub(" .*$", "", digest))
}

stage_roh_validation_genotypes <- function(source_vcf, af_sites, samples, output,
                                           bcftools = roh_validation_bcftools(),
                                           expected_source_sha256 = NULL,
                                           expected = NULL) {
  if (length(source_vcf) != 1L || !nzchar(source_vcf) ||
      length(af_sites) != 1L || !file.exists(af_sites) ||
      length(output) != 1L || !nzchar(output) ||
      length(bcftools) != 1L || !file.exists(bcftools)) {
    stop("source, AF sites, output, and bcftools must be supplied", call. = FALSE)
  }
  if (!is.character(samples) || !length(samples) || anyNA(samples) ||
      anyDuplicated(samples) || any(!grepl("^[A-Za-z0-9_.-]+$", samples))) {
    stop("`samples` must be distinct sample names", call. = FALSE)
  }
  if (!is.null(expected_source_sha256)) {
    source_hash <- if (file.exists(source_vcf)) {
      unname(digest::digest(file = source_vcf, algo = "sha256"))
    } else NA_character_
    if (!identical(source_hash, expected_source_sha256)) {
      stop("local VCF checksum differs from the expected source identity", call. = FALSE)
    }
  }
  sites <- roh_validation_site_ids(af_sites)
  if (!nrow(sites)) stop("the AF-site Parquet has no site", call. = FALSE)
  site_ids <- paste(sites$chrom, sites$pos, sites$ref, sites$alt, sep = "_")
  if (anyDuplicated(site_ids)) stop("the AF-site Parquet repeats a site", call. = FALSE)

  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  targets <- paste0(output, ".targets-", Sys.getpid(), ".tsv")
  ids <- paste0(output, ".ids-", Sys.getpid(), ".txt")
  temporary <- paste0(output, ".partial-", Sys.getpid(), ".bcf")
  temporary_index <- paste0(temporary, ".csi")
  on.exit(unlink(c(targets, ids, temporary, temporary_index), force = TRUE), add = TRUE)
  utils::write.table(sites[c("chrom", "pos")], targets, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  writeLines(site_ids, ids)

  query <- function(arguments) {
    value <- system2(bcftools, shQuote(arguments), stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(value, "status"))) {
      stop("bcftools command failed: ", paste(arguments, collapse = " "), call. = FALSE)
    }
    value
  }
  absent <- setdiff(samples, query(c("query", "-l", source_vcf)))
  if (length(absent)) {
    stop("samples absent from the VCF header: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  # The targets file streams the source, so no source index is needed. Records
  # at a site position with other alleles are removed by the ID filter.
  command <- paste(
    shQuote(bcftools), "view -s", shQuote(paste(samples, collapse = ",")),
    "-T", shQuote(targets), "-m2 -M2 -v snps -Ou", shQuote(source_vcf), "|",
    shQuote(bcftools), "annotate -x INFO,^FORMAT/GT",
    "--set-id", shQuote("%CHROM\\_%POS\\_%REF\\_%FIRST_ALT"), "-Ou |",
    shQuote(bcftools), "view -i", shQuote(paste0("ID=@", ids)), "-Ob -o", shQuote(temporary))
  status <- system2("/bin/bash", c("-o", "pipefail", "-c", shQuote(command)))
  if (status != 0L) stop("bcftools failed to build the genotype BCF", call. = FALSE)
  status <- system2(bcftools, shQuote(c("index", "-f", temporary)))
  if (status != 0L || !file.exists(temporary_index)) {
    stop("bcftools failed to index the genotype BCF", call. = FALSE)
  }
  if (!identical(query(c("query", "-l", temporary)), samples)) {
    stop("genotype BCF sample order differs from the requested samples", call. = FALSE)
  }
  record_ids <- query(c("query", "-f", "%ID\\n", temporary))
  if (!identical(sort(record_ids), sort(site_ids))) {
    stop("genotype BCF has ", length(record_ids), " records for ", length(site_ids),
         " AF sites; each site needs exactly one record", call. = FALSE)
  }
  identity <- roh_validation_genotype_identity(temporary, bcftools)[
    c("records", "samples", "records_sha256")]
  if (identity$records != length(record_ids) || identity$samples != length(samples)) {
    stop("genotype BCF identity does not match its records and samples", call. = FALSE)
  }
  if (!is.null(expected) &&
      (!identical(identity$records, as.integer(expected$records)) ||
       !identical(identity$samples, as.integer(expected$samples)) ||
       !identical(identity$records_sha256, expected$records_sha256))) {
    stop("genotype BCF differs from its registered identity: ", identity$records,
         " records, ", identity$samples, " samples, record SHA-256 ",
         identity$records_sha256, call. = FALSE)
  }
  if (!file.rename(temporary, output)) stop("could not publish the genotype BCF", call. = FALSE)
  if (!file.rename(temporary_index, paste0(output, ".csi"))) {
    unlink(output, force = TRUE)
    stop("could not publish the genotype BCF index", call. = FALSE)
  }
  c(identity, list(output = output))
}

# A genotype artifact's locator names two artifacts in order: the source VCF
# and the AF-site Parquet. Its supplier identity holds the comma-separated
# sample list and the expected records, samples and records_sha256.
stage_roh_validation_genotypes_from_registry <- function(
  id, bcftools = roh_validation_bcftools()
) {
  registry <- duckhtsbench::duckhts_bench_registry()
  row <- registry[registry$id == id, , drop = FALSE]
  if (nrow(row) != 1L) stop("unknown or non-unique benchmark artifact: ", id, call. = FALSE)
  inputs <- sub("^artifact:", "", strsplit(row$locator[[1L]], ";", fixed = TRUE)[[1L]])
  if (length(inputs) != 2L || !all(inputs %in% registry$id)) {
    stop("genotype artifact must name a source VCF and an AF-site Parquet", call. = FALSE)
  }
  fields <- duckhtsbench:::duckhts_bench_identity_fields(row$supplier_identity[[1L]])
  required <- c("sample_list", "records", "samples", "records_sha256")
  if (!all(required %in% names(fields))) {
    stop("genotype artifact identity lacks: ",
         paste(setdiff(required, names(fields)), collapse = ", "), call. = FALSE)
  }
  source_row <- registry[registry$id == inputs[[1L]], , drop = FALSE]
  source_fields <- duckhtsbench:::duckhts_bench_identity_fields(
    source_row$supplier_identity[[1L]])
  # The cached copy of the source is checked against its registered SHA-256.
  # Without a cached copy bcftools streams the registered locator, and only the
  # output identity is checked.
  source_vcf <- duckhtsbench::duckhts_bench_artifact_path(inputs[[1L]])
  expected_source_sha256 <- NULL
  if (file.exists(source_vcf)) {
    if ("sha256" %in% names(source_fields)) expected_source_sha256 <- source_fields[["sha256"]]
  } else {
    source_vcf <- source_row$locator[[1L]]
  }
  stage_roh_validation_genotypes(
    source_vcf = source_vcf,
    af_sites = duckhtsbench::duckhts_bench_artifact_path(inputs[[2L]]),
    samples = strsplit(fields[["sample_list"]], ",", fixed = TRUE)[[1L]],
    output = duckhtsbench::duckhts_bench_artifact_path(id), bcftools = bcftools,
    expected_source_sha256 = expected_source_sha256,
    expected = list(records = fields[["records"]], samples = fields[["samples"]],
                    records_sha256 = fields[["records_sha256"]]))
}
