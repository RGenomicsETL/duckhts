# ---------------------------------------------------------------------------
# Bootstrap: copy extension sources into inst/duckhts_extension/
# ---------------------------------------------------------------------------

#' Bootstrap the duckhts extension sources into the R package
#'
#' Copies extension source files from the parent duckhts repository into
#' \code{inst/duckhts_extension/} so the R package becomes self-contained.
#' Run this before \code{R CMD build} to prepare the source tarball.
#'
#' @param repo_root Path to the duckhts repository root. Required.
#' @return Invisibly returns the destination directory.
#' @export
duckhts_bootstrap <- function(repo_root = NULL) {
  if (is.null(repo_root)) {
    stop(
      "repo_root is required. Pass the duckhts repository root explicitly.",
      call. = FALSE
    )
  }
  repo_root <- normalizePath(repo_root, mustWork = TRUE)
  message("Repo root: ", repo_root)

  if (!file.exists(file.path(repo_root, "src", "duckhts.c"))) {
    stop("Not a duckhts repo: missing src/duckhts.c", call. = FALSE)
  }

  pkg_src_dir <- file.path(repo_root, "r", "Rduckhts")
  dest <- file.path(pkg_src_dir, "inst", "duckhts_extension")

  if (dir.exists(dest)) {
    unlink(dest, recursive = TRUE)
  }
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  message("Destination: ", dest)

  # Translation units and their package-relative destinations are maintained
  # in one manifest shared with CMake and both package configure scripts.
  src_dir <- file.path(repo_root, "src")
  source_manifest <- file.path(src_dir, "duckhts_sources.tsv")
  if (!file.exists(source_manifest)) {
    stop(
      "DuckHTS source manifest not found at: ",
      source_manifest,
      call. = FALSE
    )
  }
  source_rows <- duckhts_read_source_manifest(source_manifest)
  source_rows <- source_rows[nzchar(source_rows$repo_path), , drop = FALSE]
  for (i in seq_len(nrow(source_rows))) {
    source_path <- file.path(repo_root, source_rows$repo_path[[i]])
    destination_path <- file.path(dest, source_rows$package_path[[i]])
    dir.create(
      dirname(destination_path),
      recursive = TRUE,
      showWarnings = FALSE
    )
    if (!file.copy(source_path, destination_path, overwrite = TRUE)) {
      stop("Failed to copy DuckHTS source: ", source_path, call. = FALSE)
    }
  }
  file.copy(
    source_manifest,
    file.path(dest, "duckhts_sources.tsv"),
    overwrite = TRUE
  )

  duckvep_dest <- file.path(dest, "duckvep")
  dir.create(duckvep_dest, recursive = TRUE, showWarnings = FALSE)
  duckvep_headers <- c("duckvep_model.h", "duckvep_variant_tile.h")
  file.copy(file.path(src_dir, "duckvep", duckvep_headers), duckvep_dest)
  duckvep_kernel_headers <- c("duckvep_kernel.h", "duckvep_so.h")
  duckvep_kernel_private_headers <- list.files(
    file.path(src_dir, "duckvep", "kernel", "src"),
    pattern = "[.](h|inc)$"
  )
  duckvep_kernel_dest <- file.path(duckvep_dest, "kernel")
  dir.create(
    file.path(duckvep_kernel_dest, "include"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  dir.create(
    file.path(duckvep_kernel_dest, "src"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  file.copy(
    file.path(src_dir, "duckvep", "kernel", "include", duckvep_kernel_headers),
    file.path(duckvep_kernel_dest, "include")
  )
  file.copy(
    file.path(
      src_dir,
      "duckvep",
      "kernel",
      "src",
      duckvep_kernel_private_headers
    ),
    file.path(duckvep_kernel_dest, "src")
  )
  message(
    "  Copied ",
    nrow(source_rows),
    " C source files from the shared manifest"
  )

  # Headers
  inc_dest <- file.path(dest, "include")
  dir.create(inc_dest, showWarnings = FALSE)
  inc_files <- list.files(file.path(src_dir, "include"), full.names = FALSE)
  file.copy(file.path(src_dir, "include", inc_files), inc_dest)
  cgranges_dir <- file.path(repo_root, "third_party", "cgranges")
  file.copy(file.path(cgranges_dir, c("cgranges.h", "khash.h")), inc_dest)
  variantkey_inc_dest <- file.path(inc_dest, "variantkey")
  dir.create(variantkey_inc_dest, showWarnings = FALSE)
  variantkey_dir <- file.path(
    repo_root,
    "third_party",
    "variantkey",
    "include",
    "variantkey"
  )
  variantkey_files <- c("hex.h", "variantkey.h", "regionkey.h")
  file.copy(file.path(variantkey_dir, variantkey_files), variantkey_inc_dest)
  message(
    "  Copied ",
    length(inc_files) + 2L + length(variantkey_files),
    " header files"
  )

  # Vendored read-only libBigWig sources and the upstream correctness fixture.
  libbigwig_src <- file.path(repo_root, "third_party", "libBigWig")
  libbigwig_dest <- file.path(dest, "libBigWig")
  dir.create(
    file.path(libbigwig_dest, "test"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  libbigwig_files <- c(
    "LICENSE",
    "README.md",
    "SOURCE_URL",
    "COMMIT",
    "bigWig.h",
    "bigWigIO.h",
    "bwCommon.h",
    "bwValues.h"
  )
  file.copy(file.path(libbigwig_src, libbigwig_files), libbigwig_dest)
  file.copy(
    file.path(libbigwig_src, "test", "test.bw"),
    file.path(libbigwig_dest, "test", "test.bw")
  )
  extdata_dest <- file.path(pkg_src_dir, "inst", "extdata")
  dir.create(extdata_dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(
    file.path(libbigwig_src, "test", "test.bw"),
    file.path(extdata_dest, "libbigwig_test.bw"),
    overwrite = TRUE
  )
  message("  Copied vendored read-only libBigWig sources and test fixture")

  # DuckDB C API headers
  capi_dest <- file.path(dest, "duckdb_capi")
  dir.create(capi_dest, showWarnings = FALSE)
  file.copy(
    file.path(repo_root, "duckdb_capi", c("duckdb.h", "duckdb_extension.h")),
    capi_dest
  )
  message("  Copied DuckDB C API headers")

  # Apply local patch(es) to C API headers for R package only
  patch_file <- file.path(
    pkg_src_dir,
    "inst",
    "patches",
    "duckdb_capi_strict_prototypes.patch"
  )
  if (file.exists(patch_file)) {
    message("  Applying DuckDB C API patch: ", patch_file)
    patch_status <- system2("patch", c("-p1", "-d", dest, "-i", patch_file))
    if (!identical(patch_status, 0L)) {
      stop("Failed to apply DuckDB C API patch: ", patch_file, call. = FALSE)
    }
  } else {
    message("  DuckDB C API patch not found (skipping): ", patch_file)
  }

  # htslib source tree (full copy, then clean)
  htslib_src <- file.path(repo_root, "third_party", "htslib")
  if (!dir.exists(htslib_src)) {
    stop("htslib source not found at: ", htslib_src, call. = FALSE)
  }
  htslib_dest <- file.path(dest, "htslib")
  system2("cp", c("-a", htslib_src, htslib_dest))
  system2(
    "make",
    c("-C", htslib_dest, "distclean"),
    stdout = FALSE,
    stderr = FALSE
  )
  message("  Copied htslib source tree")

  render_script <- file.path(repo_root, "scripts", "render_function_catalog.py")
  if (!file.exists(render_script)) {
    stop(
      "Function catalog renderer not found at: ",
      render_script,
      call. = FALSE
    )
  }
  python <- Sys.which("python3")
  if (!nzchar(python)) {
    python <- Sys.which("python")
  }
  if (!nzchar(python)) {
    stop(
      "python3/python is required to render the function catalog",
      call. = FALSE
    )
  }
  render_status <- system2(python, c(render_script, repo_root))
  if (!identical(render_status, 0L)) {
    stop(
      "Failed to render generated documentation assets from functions.yaml",
      call. = FALSE
    )
  }
  message("  Rendered generated documentation assets from functions.yaml")

  message("Bootstrap complete. Run 'R CMD build .' to create the tarball.")
  invisible(dest)
}
