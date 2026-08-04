# Shared receipt checks for DuckVEP performance and executable-oracle evidence.
# Checked evidence is allowed to create declared output files, but the source
# revision and every other tracked path must remain unchanged while it runs.

duckvep_system2_quote <- function(value) {
  type <- if (identical(.Platform$OS.type, "windows")) "cmd" else "sh"
  unname(vapply(as.character(value), shQuote, character(1L), type = type))
}

duckvep_evidence_command <- function(
  command,
  args,
  context,
  env = character()
) {
  value <- suppressWarnings(system2(
    command,
    duckvep_system2_quote(args),
    env = env,
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(value, "status")
  if (!is.null(status) && status != 0L) {
    detail <- paste(value, collapse = "\n")
    if (nzchar(detail)) {
      detail <- paste0(":\n", detail)
    }
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
  if (!length(value)) {
    return(character())
  }
  sub("^.. ", "", value)
}

duckvep_evidence_repo_path <- function(root, path) {
  if (!nzchar(path)) {
    return("")
  }
  root <- normalizePath(root, mustWork = TRUE)
  path <- normalizePath(path, mustWork = FALSE)
  prefix <- paste0(root, .Platform$file.sep)
  if (!startsWith(path, prefix)) {
    return("")
  }
  substring(path, nchar(prefix) + 1L)
}

duckvep_evidence_assert_checkout <- function(
  root,
  revision = NULL,
  allowed_outputs = character(),
  context = "checked evidence"
) {
  current <- duckvep_evidence_revision(root)
  if (!is.null(revision) && !identical(current, revision)) {
    stop(
      paste0(
        context,
        " changed source revision from ",
        revision,
        " to ",
        current
      ),
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

duckvep_evidence_file_state <- function(path) {
  info <- file.info(path)
  if (
    nrow(info) != 1L ||
      is.na(info$size) ||
      is.na(info$mtime) ||
      is.na(info$ctime) ||
      isTRUE(info$isdir)
  ) {
    stop(paste0("cannot identify regular artifact ", path), call. = FALSE)
  }
  c(
    size = format(info$size, scientific = FALSE, trim = TRUE),
    mtime = sprintf("%.9f", as.numeric(info$mtime)),
    ctime = sprintf("%.9f", as.numeric(info$ctime))
  )
}

duckvep_evidence_sha256 <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  state <- duckvep_evidence_file_state(path)
  value <- duckvep_evidence_command(
    "sha256sum",
    path,
    paste0("sha256sum failed for ", path)
  )
  digest <- strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", digest)) {
    stop(paste0("invalid SHA-256 for ", path), call. = FALSE)
  }
  if (!identical(state, duckvep_evidence_file_state(path))) {
    stop(paste0("artifact changed while hashing: ", path), call. = FALSE)
  }
  digest
}

duckvep_evidence_cache_info_path <- function(
  cache_dir,
  species,
  cache_version,
  assembly
) {
  normalizePath(
    file.path(
      normalizePath(cache_dir, mustWork = TRUE),
      species,
      paste0(cache_version, "_", assembly),
      "info.txt"
    ),
    mustWork = TRUE
  )
}

duckvep_evidence_cache_leaf <- function(
  cache_dir,
  species,
  cache_version,
  assembly
) {
  components <- c(
    species = species,
    cache_version = cache_version,
    assembly = assembly
  )
  invalid <- names(components)[
    !grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", components) |
      components %in% c(".", "..")
  ]
  if (length(invalid)) {
    stop(
      "invalid VEP cache path component(s): ",
      paste(invalid, collapse = ", "),
      call. = FALSE
    )
  }
  normalizePath(
    file.path(
      normalizePath(cache_dir, mustWork = TRUE),
      species,
      paste0(cache_version, "_", assembly)
    ),
    mustWork = TRUE
  )
}

duckvep_evidence_sha256_lines <- function(lines) {
  path <- tempfile("duckvep-sha256-lines-")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  writeLines(enc2utf8(lines), path, useBytes = TRUE)
  duckvep_evidence_sha256(path)
}

duckvep_evidence_cache_inventory <- function(
  cache_dir,
  species,
  cache_version,
  assembly
) {
  leaf <- duckvep_evidence_cache_leaf(
    cache_dir,
    species,
    cache_version,
    assembly
  )
  entries <- list.files(
    leaf,
    all.files = TRUE,
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = TRUE,
    no.. = TRUE
  )
  if (!length(entries)) {
    stop("VEP cache leaf contains no files: ", leaf, call. = FALSE)
  }
  if (any(nzchar(Sys.readlink(entries)))) {
    stop("VEP cache receipts do not admit symbolic links", call. = FALSE)
  }
  entry_info <- file.info(entries)
  if (any(is.na(entry_info$isdir))) {
    stop("cannot stat every VEP cache entry", call. = FALSE)
  }
  paths <- entries[!entry_info$isdir]
  if (!length(paths)) {
    stop("VEP cache leaf contains no files: ", leaf, call. = FALSE)
  }
  paths <- sort(normalizePath(paths, winslash = "/", mustWork = TRUE))
  prefix <- paste0(chartr("\\", "/", leaf), "/")
  if (any(!startsWith(paths, prefix))) {
    stop("VEP cache inventory escaped its cache leaf", call. = FALSE)
  }
  relative <- substring(paths, nchar(prefix) + 1L)
  if (any(grepl("[\t\r\n]", relative))) {
    stop("VEP cache paths may not contain tabs or newlines", call. = FALSE)
  }
  info <- file.info(paths)
  if (
    any(is.na(info$size)) ||
      any(is.na(info$mtime)) ||
      any(is.na(info$ctime)) ||
      any(info$isdir)
  ) {
    stop("cannot stat every VEP cache file", call. = FALSE)
  }
  size <- format(info$size, scientific = FALSE, trim = TRUE)
  mtime <- sprintf("%.6f", as.numeric(info$mtime))
  ctime <- sprintf("%.6f", as.numeric(info$ctime))
  lines <- paste(relative, size, mtime, ctime, sep = "\t")
  list(
    leaf = leaf,
    entries = length(paths),
    bytes = format(sum(info$size), scientific = FALSE, trim = TRUE),
    sha256 = duckvep_evidence_sha256_lines(lines)
  )
}

duckvep_evidence_read_cache_receipt <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  values <- utils::read.delim(
    path,
    header = FALSE,
    col.names = c("field", "value"),
    colClasses = "character",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )
  if (
    !nrow(values) || any(!nzchar(values$field)) || anyDuplicated(values$field)
  ) {
    stop("invalid VEP cache receipt fields", call. = FALSE)
  }
  required <- c(
    "schema_version",
    "species",
    "cache_version",
    "assembly",
    "source_url",
    "source_identity",
    "cache_info_sha256",
    "inventory_kind",
    "inventory_entries",
    "inventory_bytes",
    "inventory_sha256"
  )
  missing <- setdiff(required, values$field)
  if (length(missing)) {
    stop(
      "VEP cache receipt is missing field(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  stats::setNames(values$value, values$field)
}

duckvep_evidence_write_cache_receipt <- function(
  path,
  cache_dir,
  species,
  cache_version,
  assembly,
  source_url,
  source_identity,
  overwrite = FALSE
) {
  if (!grepl("^[A-Za-z][A-Za-z0-9+.-]*://[^\t\r\n]+$", source_url)) {
    stop("cache source_url must be an absolute URL", call. = FALSE)
  }
  if (
    !grepl(
      paste0(
        "^(?:sha256:[0-9a-fA-F]{64}|md5:[0-9a-fA-F]{32}|",
        "bsd-sum:[0-9]+:[0-9]+|http-etag:[0-9a-fA-F-]+:[0-9]+)$"
      ),
      source_identity,
      perl = TRUE
    )
  ) {
    stop("invalid cache source_identity", call. = FALSE)
  }
  path <- normalizePath(path, mustWork = FALSE)
  leaf <- duckvep_evidence_cache_leaf(
    cache_dir,
    species,
    cache_version,
    assembly
  )
  path_key <- chartr("\\", "/", path)
  leaf_key <- paste0(chartr("\\", "/", leaf), "/")
  if (startsWith(path_key, leaf_key)) {
    stop(
      "cache receipt must live outside the inventoried cache leaf",
      call. = FALSE
    )
  }
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("cache receipt already exists: ", path, call. = FALSE)
  }
  inventory <- duckvep_evidence_cache_inventory(
    cache_dir,
    species,
    cache_version,
    assembly
  )
  info_path <- duckvep_evidence_cache_info_path(
    cache_dir,
    species,
    cache_version,
    assembly
  )
  fields <- c(
    schema_version = "1",
    species = species,
    cache_version = cache_version,
    assembly = assembly,
    source_url = source_url,
    source_identity = tolower(source_identity),
    cache_info_sha256 = duckvep_evidence_sha256(info_path),
    inventory_kind = "relative_path_size_mtime_ctime",
    inventory_entries = as.character(inventory$entries),
    inventory_bytes = inventory$bytes,
    inventory_sha256 = inventory$sha256
  )
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(".cache-receipt-", tmpdir = dirname(path))
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  writeLines(paste(names(fields), fields, sep = "\t"), temporary)
  if (!file.rename(temporary, path)) {
    stop("cannot publish VEP cache receipt: ", path, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

duckvep_evidence_validate_cache_receipt <- function(
  path,
  cache_dir,
  species,
  cache_version,
  assembly
) {
  fields <- duckvep_evidence_read_cache_receipt(path)
  if (
    !grepl(
      "^[A-Za-z][A-Za-z0-9+.-]*://[^\t\r\n]+$",
      fields[["source_url"]]
    ) ||
      !grepl(
        paste0(
          "^(?:sha256:[0-9a-f]{64}|md5:[0-9a-f]{32}|",
          "bsd-sum:[0-9]+:[0-9]+|http-etag:[0-9a-f-]+:[0-9]+)$"
        ),
        fields[["source_identity"]],
        perl = TRUE
      )
  ) {
    stop("invalid VEP cache receipt source identity", call. = FALSE)
  }
  expected <- c(
    schema_version = "1",
    species = species,
    cache_version = cache_version,
    assembly = assembly,
    inventory_kind = "relative_path_size_mtime_ctime"
  )
  mismatch <- names(expected)[fields[names(expected)] != expected]
  if (length(mismatch)) {
    stop(
      "VEP cache receipt mismatch for field(s): ",
      paste(mismatch, collapse = ", "),
      call. = FALSE
    )
  }
  info_path <- duckvep_evidence_cache_info_path(
    cache_dir,
    species,
    cache_version,
    assembly
  )
  info_sha256 <- duckvep_evidence_sha256(info_path)
  if (!identical(unname(fields[["cache_info_sha256"]]), info_sha256)) {
    stop("VEP cache info.txt differs from its cache receipt", call. = FALSE)
  }
  inventory <- duckvep_evidence_cache_inventory(
    cache_dir,
    species,
    cache_version,
    assembly
  )
  observed <- c(
    inventory_entries = as.character(inventory$entries),
    inventory_bytes = inventory$bytes,
    inventory_sha256 = inventory$sha256
  )
  mismatch <- names(observed)[fields[names(observed)] != observed]
  if (length(mismatch)) {
    stop(
      "VEP cache tree differs from its acquisition receipt: ",
      paste(mismatch, collapse = ", "),
      call. = FALSE
    )
  }
  list(
    receipt_sha256 = duckvep_evidence_sha256(path),
    inventory_sha256 = inventory$sha256,
    entries = inventory$entries,
    bytes = inventory$bytes,
    source_url = unname(fields[["source_url"]]),
    source_identity = unname(fields[["source_identity"]])
  )
}

duckvep_evidence_explicit_packages <- function(lines) {
  records <- trimws(lines)
  records <- records[
    nzchar(records) &
      records != "@EXPLICIT" &
      !startsWith(records, "#") &
      !startsWith(records, "List of packages in environment:")
  ]
  if (
    length(records) &&
      any(
        !grepl(
          "^[A-Za-z][A-Za-z0-9+.-]*://",
          records
        )
      )
  ) {
    stop("invalid non-URL record in explicit Conda package set", call. = FALSE)
  }
  sort(unique(records))
}

duckvep_blit_quote <- function(value, os_type = .Platform$OS.type) {
  if (!identical(os_type, "windows")) {
    return(unname(vapply(
      as.character(value),
      shQuote,
      character(1L),
      type = "sh"
    )))
  }
  # blit writes a cmd.exe batch file. Quote first for the target application's
  # argv parser, then escape that quoted token for cmd.exe itself.
  unname(vapply(
    as.character(value),
    function(item) shQuote(shQuote(item, type = "cmd"), type = "cmd2"),
    character(1L)
  ))
}

duckvep_evidence_untracked_build_inputs <- function(root) {
  value <- duckvep_evidence_command(
    "git",
    c(
      "-C",
      root,
      "ls-files",
      "--others",
      "--",
      "Makefile",
      "GNUmakefile",
      "makefile",
      "CMakeLists.txt",
      "extension_config.cmake",
      "description.yml",
      "functions.yaml",
      "duckdb_capi",
      "scripts/duckvep_evidence.R",
      "src",
      "third_party/cgranges",
      "third_party/variantkey/include",
      "third_party/htslib"
    ),
    "cannot inspect untracked DuckHTS build inputs"
  )
  value[grepl(
    paste0(
      "(?:\\.(?:c|cc|cpp|h|hpp|inc|def|gch|pch|patch|cmake)|",
      "(?:^|/)(?:Makefile|GNUmakefile|makefile|configure|CMakeLists[.]txt))$"
    ),
    value,
    perl = TRUE
  )]
}

duckvep_evidence_build_extension <- function(
  root,
  extension,
  revision,
  allowed_outputs = character()
) {
  expected <- normalizePath(
    file.path(root, "build", "release", "duckhts.duckdb_extension"),
    mustWork = FALSE
  )
  requested <- normalizePath(extension, mustWork = FALSE)
  if (!identical(requested, expected)) {
    stop(
      paste0(
        "checked evidence must use the in-tree release extension: ",
        expected
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
  build_environment_names <- names(Sys.getenv())
  build_environment_names <- build_environment_names[grepl(
    paste0(
      "^(?:MAKEFLAGS|GNUMAKEFLAGS|MAKEFILES|MFLAGS|",
      "CFLAGS|CXXFLAGS|CPPFLAGS|LDFLAGS|LIBS|",
      "CC|CXX|AR|RANLIB|LD|NM|STRIP|",
      "CPATH|C_INCLUDE_PATH|CPLUS_INCLUDE_PATH|LIBRARY_PATH|",
      "LD_LIBRARY_PATH|PKG_CONFIG_PATH|PYTHONHOME|PYTHONPATH|",
      "VIRTUAL_ENV|RUSTFLAGS|RUSTUP_TOOLCHAIN|SDKROOT|",
      "MACOSX_DEPLOYMENT_TARGET|GEN|",
      "CMAKE(?:_|$)|VCPKG(?:_|$)|DUCKDB(?:_|$)|",
      "EXTENSION(?:_|$)|EXTRA_(?:CMAKE|EXTENSION)|",
      "BUILD_EXTENSION_|DEFAULT_TEST_|FULL_TEST_|SUBSET_|",
      "USE_MERGED_|LINUX_CI_|OSX_BUILD_)"
    ),
    build_environment_names,
    perl = TRUE
  )]
  make_env <- paste0(unique(build_environment_names), "=")
  extension_ci_status <- duckvep_evidence_command(
    "git",
    c(
      "-C",
      file.path(root, "extension-ci-tools"),
      "status",
      "--porcelain=v1",
      "--untracked-files=all"
    ),
    "cannot inspect extension-ci-tools worktree"
  )
  if (length(extension_ci_status)) {
    stop("extension-ci-tools worktree is not clean", call. = FALSE)
  }
  duckvep_evidence_command(
    "make",
    c("-C", htslib, "-f", "Makefile", "distclean"),
    "failed to distclean vendored htslib",
    env = make_env
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
  untracked_build_inputs <- duckvep_evidence_untracked_build_inputs(root)
  if (length(untracked_build_inputs)) {
    stop(
      paste0(
        "checked evidence found untracked compiler input(s):\n",
        paste(untracked_build_inputs, collapse = "\n")
      ),
      call. = FALSE
    )
  }
  platform_file <- file.path(root, "configure", "platform.txt")
  configure_python <- file.path(root, "configure", "venv", "bin", "python3")
  if (!file.exists(configure_python)) {
    stop(
      paste0(
        "checked evidence needs the configured in-tree Python environment; ",
        "run make -f Makefile configure before the network-free evidence build"
      ),
      call. = FALSE
    )
  }
  duckvep_evidence_command(
    "make",
    c("-C", root, "-f", "Makefile", "platform", "extension_version"),
    "failed to refresh DuckHTS platform/version metadata",
    env = make_env
  )
  platform <- trimws(readLines(platform_file, warn = FALSE))
  if (length(platform) != 1L || !grepl("^[A-Za-z0-9_.-]+$", platform)) {
    stop("DuckHTS configure/platform.txt is invalid", call. = FALSE)
  }
  duckvep_evidence_command(
    "make",
    c("-C", root, "-f", "Makefile", "clean"),
    "failed to clean the DuckHTS release build",
    env = make_env
  )
  duckvep_evidence_command(
    "make",
    c("-C", root, "-f", "Makefile", "-j2", "release"),
    "failed to rebuild the DuckHTS release extension",
    env = make_env
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
  post_build_inputs <- setdiff(
    duckvep_evidence_untracked_build_inputs(root),
    c(
      "third_party/htslib/config.h",
      "third_party/htslib/config_vars.h",
      "third_party/htslib/version.h"
    )
  )
  if (length(post_build_inputs)) {
    stop(
      paste0(
        "checked evidence build created untracked compiler input(s):\n",
        paste(post_build_inputs, collapse = "\n")
      ),
      call. = FALSE
    )
  }
  extension_ci_status <- duckvep_evidence_command(
    "git",
    c(
      "-C",
      file.path(root, "extension-ci-tools"),
      "status",
      "--porcelain=v1",
      "--untracked-files=all"
    ),
    "cannot recheck extension-ci-tools worktree"
  )
  if (length(extension_ci_status)) {
    stop(
      "extension-ci-tools worktree changed during release build",
      call. = FALSE
    )
  }
  list(
    path = expected,
    binding = "htslib_distclean_make_release",
    sha256 = duckvep_evidence_sha256(expected)
  )
}

duckvep_evidence_write_extension_receipt <- function(
  path,
  source_revision,
  extension
) {
  if (!grepl("^[0-9a-f]{40}$", source_revision)) {
    stop("extension receipt needs a full Git object name", call. = FALSE)
  }
  required <- c("path", "binding", "sha256")
  if (!is.list(extension) || !all(required %in% names(extension))) {
    stop("extension receipt is missing build fields", call. = FALSE)
  }
  extension_path <- normalizePath(extension$path, mustWork = TRUE)
  if (!identical(extension$binding, "htslib_distclean_make_release")) {
    stop("extension receipt has an unsupported build binding", call. = FALSE)
  }
  if (!grepl("^[0-9a-f]{64}$", extension$sha256)) {
    stop("extension receipt has an invalid SHA-256", call. = FALSE)
  }
  value <- data.frame(
    field = c("source_revision", "path", "binding", "sha256"),
    value = c(
      source_revision,
      extension_path,
      extension$binding,
      extension$sha256
    ),
    stringsAsFactors = FALSE
  )
  path <- normalizePath(path.expand(path), mustWork = FALSE)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile("duckvep_extension_", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  utils::write.table(
    value,
    temporary,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = ""
  )
  if (!file.rename(temporary, path)) {
    stop(paste0("cannot publish extension receipt ", path), call. = FALSE)
  }
  path
}

duckvep_evidence_read_extension_receipt <- function(
  path,
  root,
  extension,
  source_revision
) {
  path <- normalizePath(path, mustWork = TRUE)
  value <- utils::read.delim(
    path,
    colClasses = "character",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  fields <- c("source_revision", "path", "binding", "sha256")
  if (
    !identical(names(value), c("field", "value")) ||
      nrow(value) != length(fields) ||
      !identical(value$field, fields)
  ) {
    stop("extension receipt has an invalid schema", call. = FALSE)
  }
  receipt <- stats::setNames(value$value, value$field)
  if (!identical(receipt[["source_revision"]], source_revision)) {
    stop("extension receipt belongs to another source revision", call. = FALSE)
  }
  expected <- normalizePath(
    file.path(root, "build", "release", "duckhts.duckdb_extension"),
    mustWork = FALSE
  )
  requested <- normalizePath(extension, mustWork = FALSE)
  recorded <- normalizePath(receipt[["path"]], mustWork = FALSE)
  if (!identical(requested, expected) || !identical(recorded, expected)) {
    stop(
      "extension receipt does not bind the in-tree release extension",
      call. = FALSE
    )
  }
  if (
    !identical(
      receipt[["binding"]],
      "htslib_distclean_make_release"
    )
  ) {
    stop("extension receipt has an unsupported build binding", call. = FALSE)
  }
  observed_sha256 <- duckvep_evidence_sha256(expected)
  if (!identical(receipt[["sha256"]], observed_sha256)) {
    stop("extension receipt does not match the built extension", call. = FALSE)
  }
  list(
    path = expected,
    binding = receipt[["binding"]],
    sha256 = observed_sha256
  )
}
