duckhts_bootstrap <- function(vendor_conformance = FALSE) {
    root_dir <- duckhts_root_dir()
    script <- file.path(root_dir, "scripts", "vendor_all.sh")
    if (!file.exists(script)) {
        stop("vendor_all.sh not found: ", script)
    }

    env <- character()
    if (isTRUE(vendor_conformance)) {
        env <- c(env, "VENDOR_CONFORMANCE=1")
    }

    status <- system2("bash", c(script), env = env)
    if (status != 0) {
        stop("Vendoring failed with status ", status)
    }
    invisible(TRUE)
}

duckhts_build <- function(target = "release") {
    root_dir <- duckhts_root_dir()
    makefile <- file.path(root_dir, "Makefile")
    if (!file.exists(makefile)) {
        stop("Makefile not found: ", makefile)
    }

    status <- system2("make", c(target), env = c(sprintf("PWD=%s", root_dir)))
    if (status != 0) {
        stop("Build failed with status ", status)
    }
    invisible(TRUE)
}

duckhts_load <- function(con = NULL) {
    if (!requireNamespace("duckdb", quietly = TRUE)) {
        stop("duckdb package not available")
    }
    if (is.null(con)) {
        con <- duckdb::dbConnect(duckdb::duckdb())
    }

    root_dir <- duckhts_root_dir()
    ext_path <- file.path(root_dir, "build", "release", "duckhts.duckdb_extension")
    if (!file.exists(ext_path)) {
        stop("Extension not found: ", ext_path)
    }

    duckdb::duckdb_load_extension(con, ext_path)
    invisible(con)
}

# internal
# Find repo root from the package location or working directory
# This expects the package to live at <root>/r/duckhts

duckhts_root_dir <- function() {
    pkg_dir <- normalizePath(system.file(package = "duckhts"), mustWork = FALSE)
    if (nzchar(pkg_dir)) {
        return(normalizePath(file.path(pkg_dir, "..", "..")))
    }

    wd <- normalizePath(getwd())
    # If run from within the package directory
    if (basename(wd) == "duckhts" && basename(dirname(wd)) == "r") {
        return(normalizePath(file.path(wd, "..", "..")))
    }

    # Fallback to current working directory
    return(wd)
}
