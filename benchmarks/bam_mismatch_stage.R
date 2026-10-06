# Stage the inputs of benchmark_bam_mismatch_counts.Rmd: a CRAM slice of one
# registered public CRAM, and a sites-only BCF mask of one registered panel
# VCF. samtools and bcftools are the binaries bundled by the RBCFTools package;
# tools on PATH are not used. A source is read from its cached copy when there
# is one, and otherwise where its registered locator says, so a slice of a
# large public file needs no local copy of that file.
#
# The identity of a slice is its number of alignments, the sum of their
# positions and the sum of their stored read lengths. The identity of a mask is
# its number of records and the sum of their positions. Neither depends on the
# encoding. When an expected identity is given, the staged file is compared
# with it before anything is published.

bam_mismatch_tool <- function(tool) {
  if (!requireNamespace("RBCFTools", quietly = TRUE) ||
      utils::packageVersion("RBCFTools") < "1.24.1.1.0") {
    stop("staging needs RBCFTools 1.24-1.1.0 or later for its bundled samtools and bcftools",
         call. = FALSE)
  }
  switch(tool, samtools = RBCFTools::samtools_path(), bcftools = RBCFTools::bcftools_path(),
         stop("unknown tool: ", tool, call. = FALSE))
}

# The bundled tools read a remote file through the package's htslib plugins.
# The plugins link the package's shared htslib and carry no run path to it, so
# the tool process gets both directories.
bam_mismatch_tool_environment <- function() {
  loader_path <- c(RBCFTools::htslib_lib_dir(), Sys.getenv("LD_LIBRARY_PATH", unset = NA))
  c(paste0("HTS_PATH=", shQuote(RBCFTools::htslib_plugins_dir())),
    paste0("LD_LIBRARY_PATH=", shQuote(paste(loader_path[!is.na(loader_path)], collapse = ":"))))
}

bam_mismatch_run <- function(tool, arguments) {
  status <- system2(bam_mismatch_tool(tool), shQuote(arguments), env = bam_mismatch_tool_environment())
  if (status != 0L) stop(tool, " failed: ", paste(arguments, collapse = " "), call. = FALSE)
}

bam_mismatch_connect <- function(extension) {
  driver <- duckdb::duckdb(dbdir = ":memory:", shared_home = FALSE,
    allow_extensions = TRUE, config = list(threads = "4",
      allow_unsigned_extensions = "true", autoinstall_known_extensions = "false",
      autoload_known_extensions = "false"))
  con <- DBI::dbConnect(driver)
  DBI::dbExecute(con, sprintf("LOAD %s", DBI::dbQuoteString(con, normalizePath(extension))))
  con
}

bam_mismatch_same_identity <- function(observed, expected) {
  all(names(observed) %in% names(expected)) &&
    all(vapply(names(observed), function(field) {
      identical(as.numeric(observed[[field]]), as.numeric(expected[[field]]))
    }, logical(1L)))
}

bam_mismatch_slice_identity <- function(path, reference, extension) {
  con <- bam_mismatch_connect(extension)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  identity <- DBI::dbGetQuery(con, sprintf(
    "SELECT count(*)::DOUBLE AS alignments, coalesce(sum(POS), 0)::DOUBLE AS position_sum,
            coalesce(sum(len(SEQ)), 0)::DOUBLE AS base_sum
     FROM read_bam(%s, reference := %s)",
    DBI::dbQuoteString(con, path), DBI::dbQuoteString(con, reference)))
  as.list(identity)
}

bam_mismatch_mask_identity <- function(path, extension) {
  con <- bam_mismatch_connect(extension)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  identity <- DBI::dbGetQuery(con, sprintf(
    "SELECT count(*)::DOUBLE AS records, coalesce(sum(POS), 0)::DOUBLE AS position_sum FROM read_bcf(%s)",
    DBI::dbQuoteString(con, path)))
  as.list(identity)
}

# Publish `temporary` and its index as `output`, after the identity check.
bam_mismatch_publish <- function(temporary, output, index_suffix, identity, expected, what) {
  if (!is.null(expected) && !bam_mismatch_same_identity(identity, expected)) {
    stop(what, " differs from its registered identity: ",
         paste(names(identity), format(unlist(identity), scientific = FALSE, trim = TRUE),
               sep = "=", collapse = ";"), call. = FALSE)
  }
  if (!file.rename(paste0(temporary, index_suffix), paste0(output, index_suffix)) ||
      !file.rename(temporary, output)) {
    stop("could not publish ", output, call. = FALSE)
  }
  c(identity, list(output = output))
}

bam_mismatch_is_url <- function(path) grepl("^[A-Za-z][A-Za-z0-9+.-]*://", path)

bam_mismatch_require_source <- function(source) {
  if (!bam_mismatch_is_url(source) && !file.exists(source)) {
    stop("missing staging input: ", source, call. = FALSE)
  }
}

# One region of `source` (a local file, or a URL read by range requests) as an
# indexed CRAM against `reference`.
stage_bam_mismatch_slice <- function(source, region, reference, output, extension, expected = NULL) {
  bam_mismatch_require_source(source)
  if (!file.exists(reference)) stop("missing staging input: ", reference, call. = FALSE)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  # The working directory changes below, so no path stays relative.
  output <- file.path(normalizePath(dirname(output)), basename(output))
  reference <- normalizePath(reference)
  extension <- normalizePath(extension, mustWork = TRUE)
  if (file.exists(source)) source <- normalizePath(source)
  temporary <- paste0(output, ".partial-", Sys.getpid(), ".cram")
  on.exit(unlink(c(temporary, paste0(temporary, ".crai")), force = TRUE), add = TRUE)
  # htslib saves the index of a remote file in the working directory.
  scratch <- tempfile("bam-mismatch-stage-")
  dir.create(scratch)
  previous <- setwd(scratch)
  on.exit({
    setwd(previous)
    unlink(scratch, recursive = TRUE)
  }, add = TRUE)
  bam_mismatch_run("samtools", c("view", "--no-PG", "-C", "-T", reference, "-o", temporary, source, region))
  bam_mismatch_run("samtools", c("index", temporary))
  identity <- bam_mismatch_slice_identity(temporary, reference, extension)
  bam_mismatch_publish(temporary, output, ".crai", identity, expected, "alignment slice")
}

# The records of one contig of `source` (a local file, or a URL read as a
# stream), without genotypes, as an indexed BCF. A `modulo` above 1 keeps the
# records whose position is a multiple of it: a thinner mask of the same file,
# for scaling the mask records alone.
stage_bam_mismatch_mask <- function(source, contig, output, extension, expected = NULL, modulo = 1L) {
  bam_mismatch_require_source(source)
  if (!is.numeric(modulo) || length(modulo) != 1L || modulo < 1 || modulo != floor(modulo)) {
    stop("modulo must be a whole number of at least 1", call. = FALSE)
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(output, ".partial-", Sys.getpid(), ".bcf")
  on.exit(unlink(c(temporary, paste0(temporary, ".csi")), force = TRUE), add = TRUE)
  thinning <- if (modulo > 1) c("-i", sprintf("POS %% %d == 0", as.integer(modulo))) else character()
  bam_mismatch_run("bcftools", c("view", "--no-version", "-G", "-t", contig, thinning, "-Ob", "-o", temporary, source))
  bam_mismatch_run("bcftools", c("index", "-f", temporary))
  identity <- bam_mismatch_mask_identity(temporary, extension)
  bam_mismatch_publish(temporary, output, ".csi", identity, expected, "mask")
}

bam_mismatch_registry_row <- function(id) {
  registry <- duckhtsbench::duckhts_bench_registry()
  row <- registry[registry$id == id, , drop = FALSE]
  if (nrow(row) != 1L) stop("unknown or non-unique benchmark artifact: ", id, call. = FALSE)
  inputs <- sub("^artifact:", "", strsplit(row$locator[[1L]], ";", fixed = TRUE)[[1L]])
  if (!all(inputs %in% registry$id)) stop("artifact ", id, " names an unknown input", call. = FALSE)
  list(registry = registry, inputs = inputs,
       fields = duckhtsbench:::duckhts_bench_identity_fields(row$supplier_identity[[1L]]))
}

# Where a staging input is read from: its cached copy when there is one,
# otherwise its registered locator when that is a URL.
bam_mismatch_registry_source <- function(registry, id) {
  cached <- duckhtsbench::duckhts_bench_artifact_path(id)
  if (file.exists(cached)) return(cached)
  locator <- registry$locator[registry$id == id]
  if (!bam_mismatch_is_url(locator)) {
    stop("staging input ", id, " is not cached and its locator is not a URL", call. = FALSE)
  }
  locator
}

bam_mismatch_require_fields <- function(id, fields, required) {
  if (!all(required %in% names(fields))) {
    stop("artifact ", id, " identity lacks: ", paste(setdiff(required, names(fields)), collapse = ", "),
         call. = FALSE)
  }
}

# A slice artifact's locator names the alignment source and the reference. Its
# supplier identity holds the region and the expected identity.
stage_bam_mismatch_slice_from_registry <- function(id, extension) {
  row <- bam_mismatch_registry_row(id)
  identity_fields <- c("alignments", "position_sum", "base_sum")
  bam_mismatch_require_fields(id, row$fields, c("region", identity_fields))
  if (length(row$inputs) != 2L) stop("slice artifact must name a source and a reference", call. = FALSE)
  stage_bam_mismatch_slice(
    source = bam_mismatch_registry_source(row$registry, row$inputs[[1L]]), region = row$fields[["region"]],
    reference = duckhtsbench::duckhts_bench_artifact_path(row$inputs[[2L]]),
    output = duckhtsbench::duckhts_bench_artifact_path(id), extension = extension,
    expected = as.list(row$fields[identity_fields]))
}

# A mask artifact's locator names the variant source. Its supplier identity
# holds the contig, the thinning modulo when there is one, and the expected
# identity.
stage_bam_mismatch_mask_from_registry <- function(id, extension) {
  row <- bam_mismatch_registry_row(id)
  identity_fields <- c("records", "position_sum")
  bam_mismatch_require_fields(id, row$fields, c("contig", identity_fields))
  if (length(row$inputs) != 1L) stop("mask artifact must name one variant source", call. = FALSE)
  modulo <- if ("modulo" %in% names(row$fields)) as.integer(row$fields[["modulo"]]) else 1L
  stage_bam_mismatch_mask(
    source = bam_mismatch_registry_source(row$registry, row$inputs[[1L]]),
    contig = row$fields[["contig"]],
    output = duckhtsbench::duckhts_bench_artifact_path(id), extension = extension,
    expected = as.list(row$fields[identity_fields]), modulo = modulo)
}

# The staged file of a registered artifact has its registered identity.
bam_mismatch_check_staged <- function(id, extension) {
  row <- bam_mismatch_registry_row(id)
  path <- duckhtsbench::duckhts_bench_artifact_path(id)
  if (!file.exists(path)) stop("artifact ", id, " is not staged: ", path, call. = FALSE)
  observed <- if ("alignments" %in% names(row$fields)) {
    bam_mismatch_slice_identity(path, duckhtsbench::duckhts_bench_artifact_path(row$inputs[[2L]]), extension)
  } else {
    bam_mismatch_mask_identity(path, extension)
  }
  if (!bam_mismatch_same_identity(observed, as.list(row$fields[names(observed)]))) {
    stop("staged artifact ", id, " differs from its registered identity", call. = FALSE)
  }
  invisible(observed)
}
