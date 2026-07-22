# ---------------------------------------------------------------------------
# Bootstrap: copy extension sources into inst/duckhts_extension/
# ---------------------------------------------------------------------------

duckhts_duckvep_kernel_source_files <- function() {
  c(
    "duckvep_kernel.c",
    "duckvep_so.c",
    "duckvep_sweep.c",
    "duckvep_classify.c",
    "duckvep_effect.c",
    "duckvep_sv.c",
    "duckvep_delta.c",
    "duckvep_projection.c",
    "duckvep_transcript_edit.c",
    "duckvep_hgvs.c",
    "duckvep_codon.c",
    "duckvep_coding.c",
    "duckvep_haplotype.c"
  )
}

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

  # C sources
  src_dir <- file.path(repo_root, "src")
  c_files <- c(
    "duckhts.c",
    "bcf_reader.c",
    "bam_reader.c",
    "bam_pileup.c",
    "bgzip.c",
    "hts_index_builder.c",
    "liftover_udf.c",
    "bcftools_norm_udf.c",
    "munge_udf.c",
    "kmer_udf.c",
    "interval_udf.c",
    "seq_reader.c",
    "fastq_qc.c",
    "tabix_reader.c",
    "bigwig_reader.c",
    "libbigwig_hfile_io.c",
    "hts_meta_reader.c",
    "quality_encoding.c",
    "quality_encoding_reader.c",
    "mosdepth_table.c",
    "bam_bin_counts.c",
    "samtools_idxstats_table.c",
    "bam_bed_coverage.c",
    "cgranges_api.c",
    "variantkey_udf.c",
    "bcftools_filter.c",
    "bcftools_shim.c",
    "score_udf.c",
    "vep_parser.c",
    "wasm_http_hfile.c"
  )
  file.copy(file.path(src_dir, c_files), dest)
  duckvep_files <- c(
    "duckvep_model.c",
    "duckvep_model.h",
    "duckvep_variant_tile.c",
    "duckvep_variant_tile.h",
    "duckvep_annotate.c",
    "duckvep_ensembl.c",
    "duckvep_sql.c"
  )
  duckvep_dest <- file.path(dest, "duckvep")
  dir.create(duckvep_dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(src_dir, "duckvep", duckvep_files), duckvep_dest)
  duckvep_kernel_headers <- c("duckvep_kernel.h", "duckvep_so.h")
  duckvep_kernel_sources <- duckhts_duckvep_kernel_source_files()
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
      c(duckvep_kernel_sources, duckvep_kernel_private_headers)
    ),
    file.path(duckvep_kernel_dest, "src")
  )
  simd_files <- c(
    "duckhts_simd_dispatch.c",
    "duckhts_simd_scalar.c",
    "duckhts_simd_avx2.c",
    "duckhts_simd_avx512.c",
    "duckhts_simd_neon.c",
    "duckhts_simd_wasm_simd128.c"
  )
  simd_dest <- file.path(dest, "simd")
  dir.create(simd_dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(src_dir, "simd", simd_files), simd_dest)
  message(
    "  Copied ",
    length(c_files) + length(simd_files) +
      sum(endsWith(duckvep_files, ".c")) + length(duckvep_kernel_sources),
    " C source files and ",
    sum(endsWith(duckvep_files, ".h")) + length(duckvep_kernel_headers) +
      length(duckvep_kernel_private_headers),
    " private DuckVEP headers"
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
  variantkey_dir <- file.path(repo_root, "third_party", "variantkey", "include", "variantkey")
  variantkey_files <- c("hex.h", "variantkey.h", "regionkey.h")
  file.copy(file.path(variantkey_dir, variantkey_files), variantkey_inc_dest)
  message("  Copied ", length(inc_files) + 2L + length(variantkey_files), " header files")

  # Vendored cgranges C source
  file.copy(file.path(repo_root, "third_party", "cgranges", "cgranges.c"), dest)
  message("  Copied vendored cgranges source")

  # Vendored read-only libBigWig sources and the upstream correctness fixture.
  libbigwig_src <- file.path(repo_root, "third_party", "libBigWig")
  libbigwig_dest <- file.path(dest, "libBigWig")
  dir.create(file.path(libbigwig_dest, "test"), recursive = TRUE, showWarnings = FALSE)
  libbigwig_files <- c(
    "LICENSE", "README.md", "SOURCE_URL", "COMMIT",
    "bigWig.h", "bigWigIO.h", "bwCommon.h", "bwValues.h",
    "bwRead.c", "bwValues.c"
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
  patch_file <- file.path(pkg_src_dir, "inst", "patches", "duckdb_capi_strict_prototypes.patch")
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
    stop("Function catalog renderer not found at: ", render_script, call. = FALSE)
  }
  python <- Sys.which("python3")
  if (!nzchar(python)) {
    python <- Sys.which("python")
  }
  if (!nzchar(python)) {
    stop("python3/python is required to render the function catalog", call. = FALSE)
  }
  render_status <- system2(python, c(render_script, repo_root))
  if (!identical(render_status, 0L)) {
    stop("Failed to render generated documentation assets from functions.yaml", call. = FALSE)
  }
  message("  Rendered generated documentation assets from functions.yaml")

  message("Bootstrap complete. Run 'R CMD build .' to create the tarball.")
  invisible(dest)
}

# ---------------------------------------------------------------------------
# Build: compile htslib + extension from the bundled sources
# ---------------------------------------------------------------------------

#' Build the duckhts DuckDB extension
#'
#' Compiles htslib and the duckhts extension from the sources bundled in the
#' installed R package. The built \code{.duckdb_extension} file is placed in
#' the extension directory.
#'
#' @param build_dir Where to build. Required. Use a writable location such as
#'   \code{tempdir()} when the installed package directory is read-only.
#' @param make Optional GNU make command to use (e.g., "gmake" or "make").
#'   When NULL, auto-detects gmake or make. If a non-GNU make is used,
#'   htslib's configure step will fail.
#' @param force Rebuild even if the extension file already exists.
#' @param verbose Print build output.
#' @return Path to the built \code{duckhts.duckdb_extension} file.
#' @export
duckhts_build <- function(build_dir = NULL, make = NULL, force = FALSE, verbose = TRUE) {
  ext_dir <- duckhts_extension_dir()

  if (is.null(build_dir)) {
    stop(
      "build_dir is required. Use a writable location such as tempdir().",
      call. = FALSE
    )
  }
  dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)

  ext_file <- file.path(build_dir, "duckhts.duckdb_extension")
  if (!force && file.exists(ext_file)) {
    if (verbose) {
      message("Extension already built: ", ext_file)
    }
    return(ext_file)
  }

  htslib_dir <- file.path(ext_dir, "htslib")
  hts_a <- file.path(htslib_dir, "libhts.a")

  # --- Build htslib if not already built (configure does this at install) ---
  if (!file.exists(hts_a)) {
    if (verbose) {
      message("Building htslib...")
    }
    if (is.null(make)) {
      make <- Sys.which("gmake")
      if (!nzchar(make)) {
        make <- Sys.which("make")
      }
    }
    if (!nzchar(make)) {
      stop("GNU make not found", call. = FALSE)
    }

    cfg_status <- system2(
      file.path(htslib_dir, "configure"),
      c("CFLAGS=-fPIC -O2 -std=gnu17", "--disable-plugins"),
      stdout = if (verbose) "" else FALSE,
      stderr = if (verbose) "" else FALSE
    )
    if (cfg_status != 0) {
      stop("htslib configure failed", call. = FALSE)
    }

    bld_status <- system2(
      make,
      c("-C", htslib_dir, "-j1", "lib-static"),
      stdout = if (verbose) "" else FALSE,
      stderr = if (verbose) "" else FALSE
    )
    if (bld_status != 0) stop("htslib build failed", call. = FALSE)
  } else {
    if (verbose) message("Using pre-built htslib: ", hts_a)
  }

  # --- Compile extension C files ---
  if (verbose) {
    message("Building duckhts extension...")
  }

  cc <- Sys.getenv("CC", unset = "gcc")
  uname <- system2("uname", stdout = TRUE)
  is_mac <- identical(uname, "Darwin")
  shared_ext <- if (is_mac) ".dylib" else ".so"
  shared_flags <- if (is_mac) {
    "-shared -fPIC -undefined dynamic_lookup"
  } else {
    "-shared -fPIC"
  }

  c_files <- file.path(
    ext_dir,
    c(
      "duckhts.c",
      file.path("simd", "duckhts_simd_dispatch.c"),
      file.path("simd", "duckhts_simd_scalar.c"),
      file.path("simd", "duckhts_simd_avx2.c"),
      file.path("simd", "duckhts_simd_avx512.c"),
      file.path("simd", "duckhts_simd_neon.c"),
      file.path("simd", "duckhts_simd_wasm_simd128.c"),
      "bcf_reader.c",
      "bam_reader.c",
      "bam_pileup.c",
      "bgzip.c",
      "hts_index_builder.c",
      "liftover_udf.c",
      "bcftools_norm_udf.c",
      "munge_udf.c",
      "kmer_udf.c",
      "interval_udf.c",
      "seq_reader.c",
      "fastq_qc.c",
      "tabix_reader.c",
      "bigwig_reader.c",
      "libbigwig_hfile_io.c",
      "hts_meta_reader.c",
      "quality_encoding.c",
      "quality_encoding_reader.c",
      "mosdepth_table.c",
      "bam_bin_counts.c",
      "samtools_idxstats_table.c",
      "bam_bed_coverage.c",
      "cgranges_api.c",
      "variantkey_udf.c",
      "cgranges.c",
      file.path("libBigWig", "bwRead.c"),
      file.path("libBigWig", "bwValues.c"),
      "bcftools_filter.c",
      "bcftools_shim.c",
      "score_udf.c",
      "vep_parser.c",
      file.path(
        "duckvep",
        "kernel",
        "src",
        duckhts_duckvep_kernel_source_files()
      ),
      file.path("duckvep", "duckvep_variant_tile.c"),
      file.path("duckvep", "duckvep_model.c"),
      file.path("duckvep", "duckvep_annotate.c"),
      file.path("duckvep", "duckvep_ensembl.c"),
      file.path("duckvep", "duckvep_sql.c"),
      "wasm_http_hfile.c"
    )
  )
  o_files <- file.path(build_dir, sub("\\.c$", ".o", basename(c_files)))
  includes <- paste(
    paste0("-I", file.path(ext_dir, "include")),
    paste0("-I", file.path(ext_dir, "duckdb_capi")),
    paste0("-I", htslib_dir),
    paste0("-I", file.path(ext_dir, "libBigWig")),
    paste0("-I", file.path(ext_dir, "duckvep", "kernel", "include")),
    paste0("-I", file.path(ext_dir, "duckvep", "kernel", "src"))
  )

  for (i in seq_along(c_files)) {
    cmd <- paste(
      cc, "-O2 -fPIC -Wpedantic -DNOCURL=1", includes,
      "-c", c_files[i], "-o", o_files[i]
    )
    if (verbose) {
      message("  ", basename(c_files[i]))
    }
    if (system(cmd) != 0) stop("Compile failed: ", c_files[i], call. = FALSE)
  }

  # --- Link ---
  shared_lib <- file.path(build_dir, paste0("libduckhts", shared_ext))
  link_libs <- "-lz -lbz2 -llzma -lcurl -lpthread"
  if (!is_mac) {
    link_libs <- paste(link_libs, "-lm")
  }

  link_cmd <- paste(
    cc,
    shared_flags,
    "-o",
    shared_lib,
    paste(o_files, collapse = " "),
    hts_a,
    link_libs
  )
  if (verbose) {
    message("  Linking...")
  }
  if (system(link_cmd) != 0) {
    stop("Link failed", call. = FALSE)
  }

  # Copy the shared library, then append DuckDB extension metadata
  file.copy(shared_lib, ext_file, overwrite = TRUE)
  duckhts_append_metadata(ext_file, verbose = verbose)
  if (verbose) {
    message("Extension built: ", ext_file)
  }

  ext_file
}

# ---------------------------------------------------------------------------
# Load: connect to DuckDB and load the extension
# ---------------------------------------------------------------------------

#' Load the duckhts extension into a DuckDB connection
#'
#' @param con An existing DuckDB connection, or \code{NULL} to create one.
#' @param extension_path Explicit path to the \code{.duckdb_extension} file.
#'   If \code{NULL}, uses the default location in the installed package.
#' @return The DuckDB connection (invisibly).
#' @export
duckhts_load <- function(con = NULL, extension_path = NULL) {
  if (!requireNamespace("duckdb", quietly = TRUE)) {
    stop("duckdb R package is required", call. = FALSE)
  }
  if (is.null(con)) {
    drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
    con <- DBI::dbConnect(drv)
  }

  if (is.null(extension_path)) {
    extension_path <- file.path(
      duckhts_extension_dir(),
      "build",
      "duckhts.duckdb_extension"
    )
  }
  if (!file.exists(extension_path)) {
    stop(
      "Extension not found: ",
      extension_path,
      "\nRun duckhts_build() first.",
      call. = FALSE
    )
  }

  DBI::dbExecute(con, sprintf("LOAD '%s'", extension_path))

  invisible(con)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Append DuckDB extension metadata trailer to a shared library
#' @keywords internal
duckhts_append_metadata <- function(ext_file, verbose = FALSE) {
  platform <- duckhts_detect_platform()
  api_version <- "v1.2.0"
  ext_version <- as.character(utils::packageVersion("Rduckhts"))

  if (verbose) {
    message(
      "  Appending metadata: platform=",
      platform,
      " api=",
      api_version,
      " ext=",
      ext_version
    )
  }

  # Helper: pad a string to exactly 32 bytes with null bytes
  pad32 <- function(s) {
    raw_s <- charToRaw(s)
    n <- length(raw_s)
    if (n > 32) {
      stop("Metadata field too long: ", s, call. = FALSE)
    }
    c(raw_s, raw(32L - n))
  }

  con <- file(ext_file, open = "ab")
  on.exit(close(con))

  # 1. Start signature (23 bytes): Wasm custom section header
  writeBin(
    as.raw(c(
      0x00, # Wasm custom section marker
      0x93,
      0x04, # LEB128(531)
      0x10 # name length = 16
    )),
    con
  )
  writeBin(charToRaw("duckdb_signature"), con)
  writeBin(as.raw(c(0x80, 0x04)), con) # LEB128(512)

  # 2. Eight 32-byte fields in reverse order (FIELD8..FIELD1)
  writeBin(pad32(""), con) # FIELD8 unused
  writeBin(pad32(""), con) # FIELD7 unused
  writeBin(pad32(""), con) # FIELD6 unused
  writeBin(pad32("C_STRUCT"), con) # FIELD5 abi_type
  writeBin(pad32(ext_version), con) # FIELD4 extension version
  writeBin(pad32(api_version), con) # FIELD3 DuckDB API version
  writeBin(pad32(platform), con) # FIELD2 platform
  writeBin(pad32("4"), con) # FIELD1 magic

  # 3. Signature space (256 null bytes)
  writeBin(raw(256L), con)
}

#' Detect DuckDB platform string
#' @keywords internal
duckhts_detect_platform <- function() {
  sysname <- tolower(Sys.info()[["sysname"]])
  machine <- Sys.info()[["machine"]]

  os <- switch(sysname,
    linux = "linux",
    darwin = "osx",
    windows = "windows",
    sysname
  )
  arch <- switch(machine,
    x86_64 = "amd64",
    aarch64 = "aarch64",
    arm64 = "aarch64",
    machine
  )
  paste0(os, "_", arch)
}

#' @keywords internal
duckhts_extension_dir <- function() {
  d <- system.file("duckhts_extension", package = "Rduckhts", mustWork = FALSE)
  if (!nzchar(d) || !dir.exists(d)) {
    # Fallback for development: look relative to working directory
    wd <- getwd()
    candidate <- file.path(wd, "inst", "duckhts_extension")
    if (dir.exists(candidate)) {
      return(normalizePath(candidate))
    }
    stop(
      "duckhts_extension directory not found in installed package.\n",
      "Run duckhts_bootstrap() first, then reinstall.",
      call. = FALSE
    )
  }
  normalizePath(d)
}

#' @keywords internal
duckhts_find_repo <- function() {
  # Try from working directory
  wd <- getwd()
  # If in r/Rduckhts/
  if (basename(wd) == "Rduckhts" && basename(dirname(wd)) == "r") {
    candidate <- normalizePath(file.path(wd, "..", ".."), mustWork = FALSE)
    if (file.exists(file.path(candidate, "CMakeLists.txt"))) {
      return(candidate)
    }
  }
  # If in repo root
  if (
    file.exists(file.path(wd, "CMakeLists.txt")) &&
      file.exists(file.path(wd, "src", "duckhts.c"))
  ) {
    return(wd)
  }
  stop(
    "Cannot auto-detect repo root. Pass repo_root= explicitly or ",
    "run from the duckhts repo directory.",
    call. = FALSE
  )
}
