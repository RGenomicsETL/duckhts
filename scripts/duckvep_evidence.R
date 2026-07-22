# Shared receipt checks for DuckVEP performance and executable-oracle evidence.
# Checked evidence is allowed to create declared output files, but the source
# revision and every other tracked path must remain unchanged while it runs.

duckvep_evidence_command <- function(
    command,
    args,
    context,
    env = character()) {
  value <- suppressWarnings(system2(
    command,
    args,
    env = env,
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

duckvep_evidence_sha256_cache_path <- function() {
  explicit <- Sys.getenv("DUCKVEP_EVIDENCE_DIGEST_CACHE", "")
  if (nzchar(explicit)) return(path.expand(explicit))
  cache_home <- Sys.getenv("XDG_CACHE_HOME", "")
  if (!nzchar(cache_home)) cache_home <- file.path(path.expand("~"), ".cache")
  file.path(cache_home, "duckhts", "evidence_sha256.tsv")
}

duckvep_evidence_sha256_file_state <- function(path) {
  info <- file.info(path)
  if (nrow(info) != 1L || is.na(info$size) || is.na(info$mtime) ||
      is.na(info$ctime) || isTRUE(info$isdir)) {
    stop(paste0("cannot identify regular artifact ", path), call. = FALSE)
  }
  c(
    size = format(info$size, scientific = FALSE, trim = TRUE),
    mtime = sprintf("%.9f", as.numeric(info$mtime)),
    ctime = sprintf("%.9f", as.numeric(info$ctime))
  )
}

duckvep_evidence_sha256_cache_read <- function(path) {
  if (!file.exists(path)) return(NULL)
  value <- tryCatch(
    utils::read.delim(
      path,
      colClasses = "character",
      quote = "",
      comment.char = "",
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    error = function(e) NULL
  )
  expected <- c("path", "size", "mtime", "ctime", "sha256")
  if (is.null(value) || !identical(names(value), expected)) return(NULL)
  value
}

duckvep_evidence_sha256_cache_write <- function(cache_path, row) {
  tryCatch({
    old <- duckvep_evidence_sha256_cache_read(cache_path)
    if (is.null(old)) {
      old <- data.frame(
        path = character(),
        size = character(),
        mtime = character(),
        ctime = character(),
        sha256 = character(),
        stringsAsFactors = FALSE
      )
    }
    old <- old[old$path != row$path, , drop = FALSE]
    value <- rbind(old, row)
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile("evidence_sha256_", tmpdir = dirname(cache_path))
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
    Sys.chmod(temporary, mode = "0600")
    if (file.rename(temporary, cache_path)) return(invisible(NULL))
    warning(
      paste0("could not update DuckVEP digest cache ", cache_path),
      call. = FALSE
    )
  }, error = function(error) {
    warning(
      paste0(
        "could not update DuckVEP digest cache ",
        cache_path,
        ": ",
        conditionMessage(error)
      ),
      call. = FALSE
    )
  })
  invisible(NULL)
}

duckvep_evidence_sha256 <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  mode <- tolower(Sys.getenv("DUCKVEP_ARTIFACT_VERIFY", "reuse"))
  if (!(mode %in% c("reuse", "full"))) {
    stop(
      "DUCKVEP_ARTIFACT_VERIFY must be reuse or full",
      call. = FALSE
    )
  }
  state <- duckvep_evidence_sha256_file_state(path)
  cache_path <- duckvep_evidence_sha256_cache_path()
  if (identical(mode, "reuse")) {
    cache <- duckvep_evidence_sha256_cache_read(cache_path)
    if (!is.null(cache)) {
      hit <- cache$path == path &
        cache$size == unname(state[["size"]]) &
        cache$mtime == unname(state[["mtime"]]) &
        cache$ctime == unname(state[["ctime"]]) &
        grepl("^[0-9a-f]{64}$", cache$sha256)
      hit[is.na(hit)] <- FALSE
      if (sum(hit) == 1L) return(cache$sha256[hit][[1L]])
    }
  }
  value <- duckvep_evidence_command(
    "sha256sum",
    path,
    paste0("sha256sum failed for ", path)
  )
  digest <- strsplit(value[[1L]], "[[:space:]]+")[[1L]][[1L]]
  if (!grepl("^[0-9a-f]{64}$", digest)) {
    stop(paste0("invalid SHA-256 for ", path), call. = FALSE)
  }
  state_after <- duckvep_evidence_sha256_file_state(path)
  if (!identical(state, state_after)) {
    stop(paste0("artifact changed while hashing: ", path), call. = FALSE)
  }
  duckvep_evidence_sha256_cache_write(
    cache_path,
    data.frame(
      path = path,
      size = unname(state[["size"]]),
      mtime = unname(state[["mtime"]]),
      ctime = unname(state[["ctime"]]),
      sha256 = digest,
      stringsAsFactors = FALSE
    )
  )
  digest
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
    stop("extension-ci-tools worktree changed during release build", call. = FALSE)
  }
  list(
    path = expected,
    binding = "htslib_distclean_make_release",
    sha256 = duckvep_evidence_sha256(expected)
  )
}
