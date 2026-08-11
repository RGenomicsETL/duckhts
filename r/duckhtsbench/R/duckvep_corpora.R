duckhts_bench_duckvep_corpus_definitions <- function() {
  list(
    `hprc-african4-chr22` = list(
      source = "duckvep_hprc_v2_grch38_source",
      source_index = "duckvep_hprc_v2_grch38_source_tbi",
      manifest = NULL,
      output = "duckvep_hprc_v2_grch38_chr22_african4",
      output_index = "duckvep_hprc_v2_grch38_chr22_african4_tbi",
      transform = "bcftools_chr22_african4_carried_split_sites",
      kind = "hprc"
    ),
    `sniffles2-chr22` = list(
      source = "duckvep_sniffles2_1kgp_source",
      source_index = "duckvep_sniffles2_1kgp_source_csi",
      manifest = NULL,
      output = "duckvep_sniffles2_1kgp_chr22",
      output_index = "duckvep_sniffles2_1kgp_chr22_tbi",
      transform = "bcftools_chr22_sites",
      kind = "sniffles2"
    ),
    `dbvar-chr22` = list(
      source = "duckvep_dbvar_grch38_20260127_source",
      source_index = "duckvep_dbvar_grch38_20260127_source_tbi",
      manifest = "duckvep_dbvar_grch38_20260127_manifest",
      output = "duckvep_dbvar_grch38_20260127_chr22",
      output_index = "duckvep_dbvar_grch38_20260127_chr22_tbi",
      transform = "bcftools_chr22_region",
      kind = "dbvar"
    )
  )
}

duckhts_bench_duckvep_corpus_rows <- function(definition) {
  registry <- duckhts_bench_registry()
  ids <- unlist(definition[c("source", "source_index", "manifest", "output", "output_index")],
    use.names = FALSE
  )
  ids <- ids[!is.na(ids) & nzchar(ids)]
  rows <- registry[match(ids, registry$id), , drop = FALSE]
  if (nrow(rows) != length(ids) || any(is.na(rows$id)) || anyDuplicated(rows$id)) {
    stop("DuckVEP corpus registry closure is incomplete or non-unique", call. = FALSE)
  }
  if (any(rows$workload != "duckvep-conformance-corpora")) {
    stop("DuckVEP corpus artifacts must belong to duckvep-conformance-corpora", call. = FALSE)
  }
  if (length(unique(rows$release)) != 1L) {
    stop("DuckVEP corpus source and derived artifacts must share one release identity", call. = FALSE)
  }
  names <- stats::setNames(seq_len(nrow(rows)), rows$id)
  source <- rows[names[[definition$source]], , drop = FALSE]
  source_index <- rows[names[[definition$source_index]], , drop = FALSE]
  output <- rows[names[[definition$output]], , drop = FALSE]
  output_index <- rows[names[[definition$output_index]], , drop = FALSE]
  expected_dependencies <- paste0(
    "artifact:", definition$source, ";artifact:", definition$source_index,
    if (is.null(definition$manifest)) "" else paste0(";artifact:", definition$manifest)
  )
  if (source$role != "raw_vcf_source_receipt" ||
      source$transform != "verify_remote_indexed_vcf_source_to_receipt" ||
      source_index$role != "raw_vcf_index" ||
      source_index$transform != "direct_download" ||
      output$role != "derived_chr22_vcf" ||
      output$locator != expected_dependencies ||
      output$transform != definition$transform ||
      output_index$role != "derived_chr22_vcf_index" ||
      output_index$locator != paste0("artifact:", definition$output) ||
      output_index$transform != "tabix_vcf_index" ||
      output_index$cache_relpath != paste0(output$cache_relpath, ".tbi")) {
    stop("DuckVEP corpus registry dependency or transform contract is invalid", call. = FALSE)
  }
  if (!is.null(definition$manifest)) {
    manifest <- rows[names[[definition$manifest]], , drop = FALSE]
    if (manifest$role != "source_manifest" || manifest$transform != "direct_download") {
      stop("DuckVEP corpus source manifest must be a direct registered download", call. = FALSE)
    }
  }
  source_index_identity <- duckhts_bench_identity_fields(source_index$supplier_identity)
  output_identity <- duckhts_bench_identity_fields(output$supplier_identity)
  index_identity <- duckhts_bench_identity_fields(output_index$supplier_identity)
  if (!all(c("sha256", "bytes") %in% names(source_index_identity)) ||
      !all(c("sha256", "bytes", "region") %in% names(output_identity)) ||
      !all(c("sha256", "bytes") %in% names(index_identity))) {
    stop("DuckVEP source index and derived identities must declare checksums and sizes", call. = FALSE)
  }
  rows
}

duckhts_bench_duckvep_corpus_row <- function(rows, id) {
  rows[match(id, rows$id), , drop = FALSE]
}

duckhts_bench_duckvep_atomic_table <- function(fields, path) {
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

duckhts_bench_duckvep_source_receipt <- function(rows, definition, source_index, manifest,
                                                   verification) {
  source <- duckhts_bench_duckvep_corpus_row(rows, definition$source)
  path <- duckhts_bench_artifact_path(definition$source)
  fields <- data.frame(
    field = c(
      "artifact_id", "workload", "role", "release", "source_locator", "access",
      "transform", "supplier_identity", "source_index_artifact", "source_index_path",
      if (!is.null(definition$manifest)) c("source_manifest_artifact", "source_manifest_path") else character(),
      "verification"
    ),
    value = c(
      source$id, source$workload, source$role, source$release, source$locator, source$access,
      source$transform, source$supplier_identity, definition$source_index, source_index,
      if (!is.null(definition$manifest)) c(definition$manifest, manifest) else character(),
      verification
    ),
    stringsAsFactors = FALSE
  )
  duckhts_bench_duckvep_atomic_table(fields, path)
  duckhts_bench_write_provenance(definition$source, path)
  invisible(path)
}

duckhts_bench_duckvep_validate_source <- function(rows, definition, source_index, manifest,
                                                    curl) {
  source <- duckhts_bench_duckvep_corpus_row(rows, definition$source)
  identity <- duckhts_bench_identity_fields(source$supplier_identity)
  if (!file.exists(source_index) || !file.info(source_index)$size) {
    stop("registered DuckVEP corpus source index is missing", call. = FALSE)
  }
  if (definition$kind == "hprc") {
    source_index_row <- duckhts_bench_duckvep_corpus_row(rows, definition$source_index)
    source_index_identity <- duckhts_bench_identity_fields(source_index_row$supplier_identity)
    if (!all(c("version_id", "bytes") %in% names(identity)) ||
        !nzchar(identity[["version_id"]]) ||
        !grepl(paste0("versionId=", identity[["version_id"]]), source$locator, fixed = TRUE) ||
        !all(c("version_id", "sha256", "bytes") %in% names(source_index_identity)) ||
        !grepl(paste0("versionId=", source_index_identity[["version_id"]]),
          source_index_row$locator, fixed = TRUE
        )) {
      stop("HPRC source and index must retain their registered S3 version IDs", call. = FALSE)
    }
    return("registered_s3_version_ids")
  }
  if (definition$kind == "sniffles2") {
    if (!all(c("etag", "bytes", "source_version") %in% names(identity)) ||
        !nzchar(identity[["etag"]])) {
      stop("Sniffles2 source registry row lacks its ETag", call. = FALSE)
    }
    if (!nzchar(curl)) stop("curl is required to validate the Sniffles2 source ETag", call. = FALSE)
    headers <- system2(curl, c("-fsSIL", shQuote(source$locator)), stdout = TRUE, stderr = TRUE)
    status <- attr(headers, "status")
    if (!is.null(status) && status != 0L) {
      stop("could not read Sniffles2 source headers", call. = FALSE)
    }
    etag_headers <- headers[grepl("^[[:space:]]*[Ee][Tt][Aa][Gg]:", headers)]
    etags <- trimws(sub("^[^:]+:", "", etag_headers))
    etags <- sub('^"', "", sub('"[[:space:]\\r]*$', "", etags))
    if (!length(etags) || utils::tail(etags, 1L) != identity[["etag"]]) {
      stop("Sniffles2 source ETag does not match the registry", call. = FALSE)
    }
    return(paste0("http_etag=", identity[["etag"]]))
  }
  if (definition$kind == "dbvar") {
    manifest_row <- duckhts_bench_duckvep_corpus_row(rows, definition$manifest)
    manifest_identity <- duckhts_bench_identity_fields(manifest_row$supplier_identity)
    if (is.null(manifest) || !file.exists(manifest) ||
        !all(c("md5", "file_date", "bytes") %in% names(identity)) ||
        !"required_entry_md5" %in% names(manifest_identity) ||
        manifest_identity[["required_entry_md5"]] != identity[["md5"]]) {
      stop("dbVar source registry row lacks its publisher MD5 manifest", call. = FALSE)
    }
    source_name <- basename(sub("[?].*$", "", source$locator))
    expected <- paste0(identity[["md5"]], "  ", source_name)
    if (!any(grepl(expected, readLines(manifest, warn = FALSE), fixed = TRUE))) {
      stop("dbVar source MD5 is absent from the registered manifest", call. = FALSE)
    }
    return(paste0("publisher_manifest_md5=", identity[["md5"]]))
  }
  stop("unknown DuckVEP corpus derivation", call. = FALSE)
}

duckhts_bench_duckvep_output_path <- function(row, output_dir) {
  if (is.null(output_dir)) return(duckhts_bench_artifact_path(row$id))
  prefix <- "corpora/duckvep/"
  if (!startsWith(row$cache_relpath, prefix)) {
    stop("DuckVEP corpus output path is outside its registered workload directory", call. = FALSE)
  }
  file.path(output_dir, substring(row$cache_relpath, nchar(prefix) + 1L))
}

duckhts_bench_duckvep_run_bcftools <- function(bcftools, args, error) {
  status <- system2(bcftools, args)
  if (status != 0L) stop(error, call. = FALSE)
}

duckhts_bench_duckvep_run_pipeline <- function(bcftools, commands, error) {
  command_lines <- vapply(commands, function(arguments) {
    paste(vapply(c(bcftools, arguments), shQuote, character(1L)), collapse = " ")
  }, character(1L))
  bash <- Sys.which("bash")
  if (!nzchar(bash)) stop("bash is required for fail-closed bcftools pipelines", call. = FALSE)
  pipeline <- paste(command_lines, collapse = " | ")
  status <- system2(bash, c("-o", "pipefail", "-c", shQuote(pipeline)))
  if (status != 0L) stop(error, call. = FALSE)
}

duckhts_bench_duckvep_publish_pair <- function(temporary, temporary_index, output, output_index) {
  targets <- c(output, output_index)
  temporaries <- c(temporary, temporary_index)
  backups <- paste0(targets, ".backup-", Sys.getpid())
  unlink(backups, force = TRUE)
  backed_up <- logical(2L)
  for (index in seq_along(targets)) {
    if (file.exists(targets[[index]])) {
      if (!file.rename(targets[[index]], backups[[index]])) {
        for (restore in which(backed_up)) file.rename(backups[[restore]], targets[[restore]])
        stop("could not preserve existing DuckVEP corpus artifact: ", targets[[index]], call. = FALSE)
      }
      backed_up[[index]] <- TRUE
    }
  }
  published <- logical(2L)
  for (index in seq_along(targets)) {
    if (!file.rename(temporaries[[index]], targets[[index]])) {
      for (remove in which(published)) unlink(targets[[remove]], force = TRUE)
      for (restore in which(backed_up)) file.rename(backups[[restore]], targets[[restore]])
      stop("could not atomically publish DuckVEP corpus artifact: ", targets[[index]], call. = FALSE)
    }
    published[[index]] <- TRUE
  }
  unlink(backups, force = TRUE)
  invisible(targets)
}

duckhts_bench_duckvep_write_derived_receipts <- function(rows, definition, output, output_index,
                                                           bcftools) {
  output_row <- duckhts_bench_duckvep_corpus_row(rows, definition$output)
  index_row <- duckhts_bench_duckvep_corpus_row(rows, definition$output_index)
  dependencies <- c(definition$source, definition$source_index, definition$manifest)
  dependencies <- dependencies[!vapply(dependencies, is.null, logical(1L))]
  dependency_rows <- rows[match(dependencies, rows$id), , drop = FALSE]
  versions <- system2(bcftools, "--version", stdout = TRUE, stderr = TRUE)
  status <- attr(versions, "status")
  if (!is.null(status) && status != 0L) versions <- character()
  fields <- data.frame(
    field = c(
      "artifact_id", "workload", "role", "release", "source_url", "source_id",
      "source_artifacts", "source_releases", "source_locators", "source_identities", "access", "transform",
      "supplier_identity", "output_path", "output_bytes", "output_sha256", "index_artifact",
      "index_path", "index_bytes", "index_sha256", "bcftools", "htslib", "consumer"
    ),
    value = c(
      output_row$id, output_row$workload, output_row$role, output_row$release,
      dependency_rows$locator[[1L]],
      paste(paste0(dependency_rows$id, "[", dependency_rows$supplier_identity, "]"), collapse = " | "),
      paste(dependency_rows$id, collapse = ";"),
      paste(dependency_rows$release, collapse = ";"),
      paste(dependency_rows$locator, collapse = " | "),
      paste(paste0(dependency_rows$id, "=", dependency_rows$supplier_identity), collapse = " | "),
      output_row$access, output_row$transform, output_row$supplier_identity, output,
      file.info(output)$size,
      duckhts_bench_identity_fields(output_row$supplier_identity)[["sha256"]], index_row$id,
      output_index, file.info(output_index)$size,
      duckhts_bench_identity_fields(index_row$supplier_identity)[["sha256"]],
      if (length(versions)) versions[[1L]] else "unknown",
      if (length(versions) > 1L) versions[[2L]] else "unknown", output_row$consumer
    ),
    stringsAsFactors = FALSE
  )
  receipt <- paste0(output, ".provenance.tsv")
  duckhts_bench_duckvep_atomic_table(fields, receipt)
  duckhts_bench_duckvep_atomic_table(fields, paste0(output, ".receipt.tsv"))
  index_fields <- data.frame(
    field = c(
      "artifact_id", "workload", "role", "release", "source_artifact", "source_path",
      "access", "transform", "supplier_identity", "output_path", "output_bytes",
      "output_sha256", "consumer"
    ),
    value = c(
      index_row$id, index_row$workload, index_row$role, index_row$release,
      output_row$id, output, index_row$access, index_row$transform,
      index_row$supplier_identity, output_index, file.info(output_index)$size,
      duckhts_bench_identity_fields(index_row$supplier_identity)[["sha256"]], index_row$consumer
    ),
    stringsAsFactors = FALSE
  )
  duckhts_bench_duckvep_atomic_table(index_fields, paste0(output_index, ".provenance.tsv"))
  invisible(receipt)
}

duckhts_bench_duckvep_derive_corpus <- function(rows, definition, source_index, output,
                                                 bcftools, tabix) {
  source <- duckhts_bench_duckvep_corpus_row(rows, definition$source)
  output_row <- duckhts_bench_duckvep_corpus_row(rows, definition$output)
  output_index_row <- duckhts_bench_duckvep_corpus_row(rows, definition$output_index)
  identity <- duckhts_bench_identity_fields(output_row$supplier_identity)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    pattern = paste0(".", basename(output), ".partial-"), tmpdir = dirname(output),
    fileext = ".vcf.gz"
  )
  temporary_index <- paste0(temporary, ".tbi")
  work <- tempfile("duckhtsbench-duckvep-corpus-")
  dir.create(work)
  on.exit(unlink(c(temporary, temporary_index, work), recursive = TRUE, force = TRUE), add = TRUE)
  source_spec <- paste0(source$locator, "##idx##", source_index)
  run <- function(args, error) duckhts_bench_duckvep_run_bcftools(bcftools, args, error)

  if (definition$kind == "hprc") {
    required <- c(
      region = "chr22", samples = "HG02055,HG02145,HG02723,HG03098",
      carried_alt_only = "true", split_multiallelic = "true",
      genotypes_removed = "true", chrom_alias = "chr22:22"
    )
    if (!identical(unname(identity[names(required)]), unname(required))) {
      stop("HPRC derivation metadata does not match the implemented transform", call. = FALSE)
    }
    aliases <- strsplit(identity[["chrom_alias"]], ":", fixed = TRUE)[[1L]]
    rename <- file.path(work, "rename.tsv")
    writeLines(paste(aliases, collapse = "\t"), rename)
    duckhts_bench_duckvep_run_pipeline(bcftools, list(
      c(
        "view", "--no-version", "-r", identity[["region"]], "-s", identity[["samples"]],
        "-Ou", source_spec
      ),
      c("view", "--no-version", "-i", 'GT="alt"', "-Ou"),
      c("norm", "--no-version", "-m", "-any", "-Ou"),
      c("view", "--no-version", "-i", 'GT="alt"', "-Ou"),
      c("view", "--no-version", "-G", "-Ou"),
      c("annotate", "--no-version", "--rename-chrs", rename, "-Oz", "-o", temporary)
    ), "could not derive the registered HPRC chr22 corpus")
  } else if (definition$kind == "sniffles2") {
    required <- c(region = "chr22", genotypes_removed = "true", chrom_alias = "chr22:22")
    if (!identical(unname(identity[names(required)]), unname(required))) {
      stop("Sniffles2 derivation metadata does not match the implemented transform", call. = FALSE)
    }
    aliases <- strsplit(identity[["chrom_alias"]], ":", fixed = TRUE)[[1L]]
    rename <- file.path(work, "rename.tsv")
    writeLines(paste(aliases, collapse = "\t"), rename)
    duckhts_bench_duckvep_run_pipeline(bcftools, list(
      c("view", "--no-version", "-r", identity[["region"]], "-G", "-Ou", source_spec),
      c("annotate", "--no-version", "--rename-chrs", rename, "-Oz", "-o", temporary)
    ), "could not derive the registered Sniffles2 chr22 corpus")
  } else if (definition$kind == "dbvar") {
    required <- c(region = "22", file_date = "20260127")
    if (!identical(unname(identity[names(required)]), unname(required))) {
      stop("dbVar derivation metadata does not match the implemented transform", call. = FALSE)
    }
    run(c(
      "view", "--no-version", "-r", identity[["region"]], "-Oz", "-o",
      shQuote(temporary), shQuote(source_spec)
    ), "could not select the registered dbVar chr22 corpus")
  }
  if (!file.exists(temporary) || !file.info(temporary)$size) {
    stop("DuckVEP corpus derivation produced no VCF", call. = FALSE)
  }
  status <- system2(tabix, c("-f", "-p", "vcf", shQuote(temporary)))
  if (status != 0L || !file.exists(temporary_index) || !file.info(temporary_index)$size) {
    stop("could not index the derived DuckVEP corpus", call. = FALSE)
  }
  duckhts_bench_validate_identity(output_row$id, temporary)
  duckhts_bench_validate_identity(output_index_row$id, temporary_index)
  duckhts_bench_duckvep_publish_pair(
    temporary, temporary_index, output, paste0(output, ".tbi")
  )
  invisible(output)
}

#' Stage deterministic chr22 corpora for DuckVEP conformance campaigns.
#'
#' Complete HPRC, Sniffles2, and dbVar source VCFs are read through their
#' registered remote indexes; only source indexes, source identity receipts,
#' and deterministic chr22 derivatives are stored in the DuckHTS cache.
#'
#' @param corpus One of `all`, `hprc-african4-chr22`, `sniffles2-chr22`, or
#'   `dbvar-chr22`.
#' @param output_dir Optional compatibility destination. By default every output
#'   uses its registry cache path.
#' @param bcftools Path to `bcftools`.
#' @param tabix Path to `tabix`.
#' @param curl Path to `curl`, used for the registered Sniffles2 ETag check.
#' @return A named character vector of derived VCF paths, invisibly.
#' @export
duckhts_bench_stage_duckvep_corpora <- function(
    corpus = "all", output_dir = NULL, bcftools = Sys.which("bcftools"),
    tabix = Sys.which("tabix"), curl = Sys.which("curl")) {
  definitions <- duckhts_bench_duckvep_corpus_definitions()
  if (length(corpus) != 1L || is.na(corpus) || !corpus %in% c("all", names(definitions))) {
    stop("unknown DuckVEP corpus: ", paste(corpus, collapse = ", "), call. = FALSE)
  }
  if (!is.null(output_dir) &&
      (length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir))) {
    stop("output_dir must be NULL or one non-empty path", call. = FALSE)
  }
  if (!nzchar(bcftools) || !nzchar(tabix)) {
    stop("bcftools and tabix are required to stage DuckVEP corpora", call. = FALSE)
  }
  selected <- if (corpus == "all") definitions else definitions[corpus]
  outputs <- character(length(selected))
  names(outputs) <- names(selected)
  for (name in names(selected)) {
    definition <- selected[[name]]
    rows <- duckhts_bench_duckvep_corpus_rows(definition)
    source_index <- duckhts_bench_fetch(definition$source_index)
    manifest <- NULL
    if (!is.null(definition$manifest)) {
      manifest <- duckhts_bench_fetch(definition$manifest, overwrite = TRUE)
    }
    verification <- duckhts_bench_duckvep_validate_source(
      rows, definition, source_index, manifest, curl
    )
    duckhts_bench_duckvep_source_receipt(
      rows, definition, source_index, manifest, verification
    )
    output_row <- duckhts_bench_duckvep_corpus_row(rows, definition$output)
    output <- duckhts_bench_duckvep_output_path(output_row, output_dir)
    output_index <- paste0(output, ".tbi")
    ready <- file.exists(output) && file.exists(output_index) &&
      isTRUE(tryCatch({
        duckhts_bench_validate_identity(definition$output, output)
        duckhts_bench_validate_identity(definition$output_index, output_index)
        TRUE
      }, error = function(error) FALSE))
    if (!ready) {
      duckhts_bench_duckvep_derive_corpus(
        rows, definition, source_index, output, bcftools, tabix
      )
    }
    duckhts_bench_duckvep_write_derived_receipts(
      rows, definition, output, output_index, bcftools
    )
    outputs[[name]] <- output
  }
  invisible(outputs)
}
