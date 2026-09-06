#' Stage complete region shards from a pinned VEP cache archive.
#'
#' The archive is streamed through GNU tar; only registered region subtrees and
#' every top-level metadata file are retained. No variant or transcript selection
#' occurs here. Publication requires a complete transfer and successful extraction.
#'
#' @param repo DuckHTS checkout supplying the canonical cache receipt functions.
#' @param id Derived cache artifact in the benchmark registry.
#' @param archive Optional local archive with a registry-pinned content checksum.
#' @param curl Curl executable for streaming HTTP acquisition.
#' @return Named cache root, info.txt, and acquisition receipt paths.
#' @export
duckhts_bench_stage_vep_cache <- function(
  repo, id = "vep116_grch38_cache_chr21", archive = NULL, curl = Sys.which("curl")
) {
  registry <- duckhts_bench_registry()
  cache <- registry[registry$id == id, , drop = FALSE]
  if (nrow(cache) != 1L || cache$transform != "extract_vep_cache_regions" ||
      !grepl("^artifact:[A-Za-z0-9_]+$", cache$locator)) {
    stop("expected one registered regional VEP cache", call. = FALSE)
  }
  source_id <- sub("^artifact:", "", cache$locator)
  source_row <- registry[registry$id == source_id, , drop = FALSE]
  if (nrow(source_row) != 1L || source_row$role != "vep_cache_archive") {
    stop("regional VEP cache must name its source archive", call. = FALSE)
  }
  scope <- duckhts_bench_identity_fields(cache$supplier_identity)
  required <- c("species", "cache_version", "assembly", "regions")
  if (!all(required %in% names(scope)) ||
      any(!grepl("^[A-Za-z0-9_-]+$", scope[required]))) {
    stop("invalid VEP cache species/version/assembly/region scope", call. = FALSE)
  }
  region <- scope[["regions"]]
  prefix <- paste(scope[["species"]], paste(scope[["cache_version"]],
    scope[["assembly"]], sep = "_"), sep = "/")
  identity <- duckhts_bench_identity_fields(source_row$supplier_identity)
  source_identity <- if ("sha256" %in% names(identity)) {
    paste0("sha256:", identity[["sha256"]])
  } else if ("md5" %in% names(identity)) {
    paste0("md5:", identity[["md5"]])
  } else if (all(c("etag", "bytes") %in% names(identity))) {
    paste0("http-etag:", identity[["etag"]], ":", identity[["bytes"]])
  } else stop("VEP archive needs a checksum or ETag/byte identity", call. = FALSE)
  evidence <- new.env(parent = baseenv())
  sys.source(file.path(normalizePath(repo), "scripts", "duckvep_evidence.R"), evidence)
  output <- duckhts_bench_artifact_path(id)
  paths <- function(directory) c(
    cache = directory, info = file.path(directory, prefix, "info.txt"),
    receipt = file.path(directory, "acquisition.tsv")
  )
  validate <- function(directory) {
    result <- paths(directory)
    fields <- evidence$duckvep_evidence_read_cache_receipt(result[["receipt"]])
    if (!identical(unname(fields["source_url"]), source_row$locator) ||
        !identical(unname(fields["source_identity"]), source_identity) ||
        !identical(unname(fields["cache_regions"]), region)) {
      stop("VEP cache receipt differs from the registered source or scope", call. = FALSE)
    }
    evidence$duckvep_evidence_validate_cache_receipt(result[["receipt"]], directory,
      scope[["species"]], scope[["cache_version"]], scope[["assembly"]])
    result
  }
  if (file.exists(output)) return(validate(output))
  tar <- Sys.which("tar")
  if (!nzchar(tar) || !any(grepl("GNU tar", system2(tar, "--version", stdout = TRUE)))) {
    stop("regional VEP cache staging requires GNU tar", call. = FALSE)
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  staging <- tempfile(".vep-cache-", tmpdir = dirname(output))
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  headers_path <- file.path(staging, "source.headers")
  log_path <- file.path(staging, "acquisition.log")
  members_path <- file.path(staging, "archive-members.txt")
  tar_args <- c("--extract", "--gzip", "--verbose", "--ignore-zeros",
    "--no-same-owner", "--no-same-permissions", "--keep-old-files",
    "--directory", staging, "--file", if (is.null(archive)) "-" else normalizePath(archive),
    "--no-recursion", "--wildcards", "--no-wildcards-match-slash", paste0(prefix, "/*"),
    "--wildcards-match-slash", paste0(prefix, "/", region, "/*"))
  if (!is.null(archive)) {
    if (!any(c("sha256", "md5") %in% names(identity))) {
      stop("a local VEP archive requires a registry-pinned content checksum", call. = FALSE)
    }
    duckhts_bench_validate_identity(source_id, archive)
    status <- system2(tar, shQuote(tar_args), stdout = members_path, stderr = log_path)
  } else {
    if (any(c("sha256", "md5") %in% names(identity))) {
      stop("supply a checksum-pinned archive locally; streaming verifies ETag/bytes", call. = FALSE)
    }
    if (!all(c("etag", "bytes") %in% names(identity)) ||
        !grepl("^[0-9a-fA-F-]+$", identity[["etag"]]) ||
        !grepl("^[1-9][0-9]*$", identity[["bytes"]]) || !nzchar(curl) ||
        !nzchar(Sys.which("bash")) || !startsWith(source_row$locator, "https://")) {
      stop("streaming VEP acquisition requires HTTPS, pinned ETag/bytes, curl and bash", call. = FALSE)
    }
    curl_args <- c("--fail", "--silent", "--show-error", "--location",
      "--proto", "=https", "--proto-redir", "=https", "--connect-timeout", "30",
      "--max-time", "7200", "--header", paste0('If-Match: "', identity[["etag"]], '"'),
      "--dump-header", headers_path, "--write-out",
      "%{stderr}\nduckhts_cache_bytes=%{size_download}\n", source_row$locator)
    command <- paste(paste(shQuote(c(curl, curl_args)), collapse = " "), "|",
      paste(shQuote(c(tar, tar_args)), collapse = " "))
    status <- system2(Sys.which("bash"), c("-o", "pipefail", "-c", shQuote(command)),
      stdout = members_path, stderr = log_path)
    if (status == 0L) {
      headers <- readLines(headers_path, warn = FALSE)
      etags <- trimws(sub("^[^:]+:", "", headers[grepl("^ETag:", headers, ignore.case = TRUE)]))
      etags <- sub('^"(.*)"$', "\\1", etags)
      transferred <- sub("^duckhts_cache_bytes=", "",
        grep("^duckhts_cache_bytes=", readLines(log_path, warn = FALSE), value = TRUE))
      if (!identical(utils::tail(etags, 1L), identity[["etag"]]) ||
          !identical(transferred, identity[["bytes"]])) {
        stop("VEP archive transfer does not match its registered ETag/byte identity", call. = FALSE)
      }
    }
  }
  if (status != 0L) {
    stop("VEP archive acquisition/extraction failed:\n",
      paste(utils::tail(readLines(log_path, warn = FALSE), 8L), collapse = "\n"), call. = FALSE)
  }
  leaf <- file.path(staging, prefix)
  files <- list.files(leaf, recursive = TRUE, all.files = TRUE, full.names = TRUE)
  relative <- substring(files, nchar(leaf) + 2L)
  # The selector includes all root metadata, but no files in other region trees.
  if (!length(files) || !file.exists(file.path(leaf, "info.txt")) ||
      !any(startsWith(relative, paste0(region, "/"))) ||
      any(grepl("/", relative, fixed = TRUE) & !startsWith(relative, paste0(region, "/")))) {
    stop("VEP archive did not yield exactly the registered region and root metadata", call. = FALSE)
  }
  evidence$duckvep_evidence_write_cache_receipt(paths(staging)[["receipt"]], staging,
    scope[["species"]], scope[["cache_version"]], scope[["assembly"]],
    source_row$locator, source_identity)
  write(paste("cache_regions", region, sep = "\t"), paths(staging)[["receipt"]], append = TRUE)
  validate(staging)
  if (file.exists(output) || !file.rename(staging, output)) {
    stop("could not publish verified VEP cache; existing output was preserved", call. = FALSE)
  }
  duckhts_bench_write_provenance(id, output)
  validate(output)
}
