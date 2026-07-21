# Shared receipt checks for DuckVEP performance and executable-oracle evidence.
# Checked evidence is allowed to create declared output files, but the source
# revision and every other tracked path must remain unchanged while it runs.

duckvep_evidence_command <- function(command, args, context) {
  value <- suppressWarnings(system2(
    command,
    args,
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(value, "status")
  if (!is.null(status) && status != 0L) {
    detail <- paste(value, collapse = "\n")
    if (nzchar(detail)) detail <- paste0(":\n", detail)
    stop(paste0(context, detail), call. = FALSE)
  }
  value
}

duckvep_evidence_revision <- function(root) {
  value <- duckvep_evidence_command(
    "git",
    c("-C", root, "rev-parse", "HEAD"),
    "cannot determine DuckHTS source revision"
  )
  revision <- trimws(value[[1L]])
  if (!grepl("^[0-9a-f]{40}$", revision)) {
    stop("DuckHTS source revision is not a full Git object name", call. = FALSE)
  }
  revision
}

duckvep_evidence_tracked_changes <- function(root) {
  value <- duckvep_evidence_command(
    "git",
    c("-C", root, "status", "--porcelain=v1", "--untracked-files=no"),
    "cannot inspect DuckHTS tracked files"
  )
  if (!length(value)) return(character())
  sub("^.. ", "", value)
}

duckvep_evidence_repo_path <- function(root, path) {
  if (!nzchar(path)) return("")
  root <- normalizePath(root, mustWork = TRUE)
  path <- normalizePath(path, mustWork = FALSE)
  prefix <- paste0(root, .Platform$file.sep)
  if (!startsWith(path, prefix)) return("")
  substring(path, nchar(prefix) + 1L)
}

duckvep_evidence_assert_checkout <- function(
    root,
    revision = NULL,
    allowed_outputs = character(),
    context = "checked evidence") {
  current <- duckvep_evidence_revision(root)
  if (!is.null(revision) && !identical(current, revision)) {
    stop(
      paste0(context, " changed source revision from ", revision, " to ", current),
      call. = FALSE
    )
  }
  allowed <- unique(vapply(
    allowed_outputs,
    function(path) duckvep_evidence_repo_path(root, path),
    character(1L)
  ))
  allowed <- allowed[nzchar(allowed)]
  changed <- duckvep_evidence_tracked_changes(root)
  unexpected <- setdiff(changed, allowed)
  if (length(unexpected)) {
    stop(
      paste0(
        context,
        " requires a clean tracked checkout; changed path(s):\n",
        paste(unexpected, collapse = "\n")
      ),
      call. = FALSE
    )
  }
  current
}

duckvep_evidence_sha256 <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  value <- duckvep_evidence_command(
    "sha256sum",
    path,
    paste0("sha256sum failed for ", path)
  )
  digest <- strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", digest)) {
    stop(paste0("invalid SHA-256 for ", path), call. = FALSE)
  }
  digest
}

duckvep_evidence_build_extension <- function(
    root,
    extension,
    revision,
    allowed_outputs = character()) {
  expected <- normalizePath(
    file.path(root, "build", "release", "duckhts.duckdb_extension"),
    mustWork = FALSE
  )
  requested <- normalizePath(extension, mustWork = FALSE)
  if (!identical(requested, expected)) {
    stop(
      paste0(
        "checked evidence must use the in-tree release extension: ", expected
      ),
      call. = FALSE
    )
  }
  duckvep_evidence_assert_checkout(
    root,
    revision,
    allowed_outputs,
    context = "DuckHTS release build"
  )
  htslib <- file.path(root, "third_party", "htslib")
  if (!file.exists(file.path(htslib, "Makefile"))) {
    stop("vendored htslib Makefile is missing", call. = FALSE)
  }
  duckvep_evidence_command(
    "make",
    c("-C", htslib, "distclean"),
    "failed to distclean vendored htslib"
  )
  htslib_files <- list.files(
    htslib,
    recursive = TRUE,
    all.files = TRUE,
    include.dirs = FALSE
  )
  stale_htslib <- htslib_files[grepl(
    paste0(
      "(^|/)(config[.](h|mk|status|cache)|config_vars[.]h|",
      "hts-object-files|libhts([.][^/]+)?)$|[.](o|pico)$"
    ),
    htslib_files
  )]
  if (length(stale_htslib)) {
    stop(
      paste0(
        "vendored htslib distclean left build artifact(s):\n",
        paste(stale_htslib, collapse = "\n")
      ),
      call. = FALSE
    )
  }
  duckvep_evidence_command(
    "make",
    c("-C", root, "clean"),
    "failed to clean the DuckHTS release build"
  )
  duckvep_evidence_command(
    "make",
    c("-C", root, "-j2", "release"),
    "failed to rebuild the DuckHTS release extension"
  )
  duckvep_evidence_assert_checkout(
    root,
    revision,
    allowed_outputs,
    context = "DuckHTS release build"
  )
  if (!file.exists(expected)) {
    stop("DuckHTS release build did not produce the extension", call. = FALSE)
  }
  list(
    path = expected,
    binding = "htslib_distclean_make_release",
    sha256 = duckvep_evidence_sha256(expected)
  )
}
