# ---------------------------------------------------------------------------
# Bootstrap: copy extension sources into inst/duckhts_extension/
# ---------------------------------------------------------------------------

#' Bootstrap the duckhts extension sources into the R package
#'
#' Copies extension source files from the parent duckhts repository into
#' \code{inst/duckhts_extension/} so the R package becomes self-contained.
#' Run this before \code{R CMD build} to prepare the source tarball.
#'
#' @param repo_root Path to the duckhts repository root. Default auto-detects
#'   from the package source directory (assumes \code{<repo>/r/duckhts}).
#' @return Invisibly returns the destination directory.
#' @export
duckhts_bootstrap <- function(repo_root = NULL) {
    if (is.null(repo_root)) {
        repo_root <- duckhts_find_repo()
    }
    repo_root <- normalizePath(repo_root, mustWork = TRUE)
    message("Repo root: ", repo_root)

    if (!file.exists(file.path(repo_root, "src", "duckhts.c"))) {
        stop("Not a duckhts repo: missing src/duckhts.c", call. = FALSE)
    }

    pkg_src_dir <- file.path(repo_root, "r", "duckhts")
    dest <- file.path(pkg_src_dir, "inst", "duckhts_extension")

    if (dir.exists(dest)) unlink(dest, recursive = TRUE)
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    message("Destination: ", dest)

    # C sources
    src_dir <- file.path(repo_root, "src")
    c_files <- c(
        "duckhts.c", "bcf_reader.c", "bam_reader.c",
        "seq_reader.c", "tabix_reader.c", "vep_parser.c"
    )
    file.copy(file.path(src_dir, c_files), dest)
    message("  Copied ", length(c_files), " C source files")

    # Headers
    inc_dest <- file.path(dest, "include")
    dir.create(inc_dest, showWarnings = FALSE)
    inc_files <- list.files(file.path(src_dir, "include"), full.names = FALSE)
    file.copy(file.path(src_dir, "include", inc_files), inc_dest)
    message("  Copied ", length(inc_files), " header files")

    # DuckDB C API headers
    capi_dest <- file.path(dest, "duckdb_capi")
    dir.create(capi_dest, showWarnings = FALSE)
    file.copy(
        file.path(repo_root, "duckdb_capi", c("duckdb.h", "duckdb_extension.h")),
        capi_dest
    )
    message("  Copied DuckDB C API headers")

    # htslib source tree (full copy, then clean)
    htslib_src <- file.path(repo_root, "third_party", "htslib")
    if (!dir.exists(htslib_src)) {
        stop("htslib source not found at: ", htslib_src, call. = FALSE)
    }
    htslib_dest <- file.path(dest, "htslib")
    system2("cp", c("-a", htslib_src, htslib_dest))
    system2("make", c("-C", htslib_dest, "distclean"),
        stdout = FALSE, stderr = FALSE
    )
    message("  Copied htslib source tree")

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
#' @param build_dir Where to build. Defaults to a \code{build/} sub-directory
#'   inside the bundled extension sources. When the installed package directory
#'   is read-only, set this to a writable location such as \code{tempdir()}.
#' @param force Rebuild even if the extension file already exists.
#' @param verbose Print build output.
#' @return Path to the built \code{duckhts.duckdb_extension} file.
#' @export
duckhts_build <- function(build_dir = NULL, force = FALSE, verbose = TRUE) {
    ext_dir <- duckhts_extension_dir()

    if (is.null(build_dir)) {
        build_dir <- file.path(ext_dir, "build")
    }
    dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)

    ext_file <- file.path(build_dir, "duckhts.duckdb_extension")
    if (!force && file.exists(ext_file)) {
        if (verbose) message("Extension already built: ", ext_file)
        return(ext_file)
    }

    htslib_dir <- file.path(ext_dir, "htslib")
    hts_a <- file.path(htslib_dir, "libhts.a")

    # --- Build htslib if not already built (configure does this at install) ---
    if (!file.exists(hts_a)) {
        if (verbose) message("Building htslib...")
        make <- Sys.which("gmake")
        if (!nzchar(make)) make <- Sys.which("make")
        if (!nzchar(make)) stop("GNU make not found", call. = FALSE)

        cfg_status <- system2(
            file.path(htslib_dir, "configure"),
            c("CFLAGS=-fPIC -O2", "--disable-plugins"),
            stdout = if (verbose) "" else FALSE,
            stderr = if (verbose) "" else FALSE
        )
        if (cfg_status != 0) stop("htslib configure failed", call. = FALSE)

        bld_status <- system2(
            make, c("-C", htslib_dir, "-j", "lib-static"),
            stdout = if (verbose) "" else FALSE,
            stderr = if (verbose) "" else FALSE
        )
        if (bld_status != 0) stop("htslib build failed", call. = FALSE)
    } else {
        if (verbose) message("Using pre-built htslib: ", hts_a)
    }

    # --- Compile extension C files ---
    if (verbose) message("Building duckhts extension...")

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
            "duckhts.c", "bcf_reader.c", "bam_reader.c",
            "seq_reader.c", "tabix_reader.c", "vep_parser.c"
        )
    )
    o_files <- file.path(build_dir, sub("\\.c$", ".o", basename(c_files)))
    includes <- paste(
        paste0("-I", file.path(ext_dir, "include")),
        paste0("-I", file.path(ext_dir, "duckdb_capi")),
        paste0("-I", htslib_dir)
    )

    for (i in seq_along(c_files)) {
        cmd <- paste(cc, "-O2 -fPIC", includes, "-c", c_files[i], "-o", o_files[i])
        if (verbose) message("  ", basename(c_files[i]))
        if (system(cmd) != 0) stop("Compile failed: ", c_files[i], call. = FALSE)
    }

    # --- Link ---
    shared_lib <- file.path(build_dir, paste0("libduckhts", shared_ext))
    link_libs <- "-lz -lbz2 -llzma -lcurl -lpthread"
    if (!is_mac) link_libs <- paste(link_libs, "-lm")

    link_cmd <- paste(
        cc, shared_flags, "-o", shared_lib,
        paste(o_files, collapse = " "),
        hts_a, link_libs
    )
    if (verbose) message("  Linking...")
    if (system(link_cmd) != 0) stop("Link failed", call. = FALSE)

    # The .so IS the extension (unsigned, no metadata append)
    file.copy(shared_lib, ext_file, overwrite = TRUE)
    if (verbose) message("Extension built: ", ext_file)

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
        con <- DBI::dbConnect(duckdb::duckdb())
    }

    if (is.null(extension_path)) {
        extension_path <- file.path(
            duckhts_extension_dir(),
            "build", "duckhts.duckdb_extension"
        )
    }
    if (!file.exists(extension_path)) {
        stop("Extension not found: ", extension_path,
            "\nRun duckhts_build() first.",
            call. = FALSE
        )
    }

    DBI::dbExecute(con, "SET allow_unsigned_extensions = true")
    DBI::dbExecute(con, sprintf("LOAD '%s'", extension_path))

    invisible(con)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' @keywords internal
duckhts_extension_dir <- function() {
    d <- system.file("duckhts_extension", package = "duckhts", mustWork = FALSE)
    if (!nzchar(d) || !dir.exists(d)) {
        # Fallback for development: look relative to working directory
        wd <- getwd()
        candidate <- file.path(wd, "inst", "duckhts_extension")
        if (dir.exists(candidate)) {
            return(normalizePath(candidate))
        }
        stop("duckhts_extension directory not found in installed package.\n",
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
    # If in r/duckhts/
    if (basename(wd) == "duckhts" && basename(dirname(wd)) == "r") {
        candidate <- normalizePath(file.path(wd, "..", ".."), mustWork = FALSE)
        if (file.exists(file.path(candidate, "CMakeLists.txt"))) {
            return(candidate)
        }
    }
    # If in repo root
    if (file.exists(file.path(wd, "CMakeLists.txt")) &&
        file.exists(file.path(wd, "src", "duckhts.c"))) {
        return(wd)
    }
    stop("Cannot auto-detect repo root. Pass repo_root= explicitly or ",
        "run from the duckhts repo directory.",
        call. = FALSE
    )
}
