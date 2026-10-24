#!/usr/bin/env bash
# Shared steps of the CRAN-flavour check scripts (check_cran_fedora.sh,
# check_cran_nold.sh). Source this file; it defines functions only.
#
# Each lane checks one source tarball, records what it ran on, and keeps its
# evidence under $DUCKHTS_CHECK_OUT so a failure can be read after the container
# is gone.

# Install or upgrade the check dependencies to the newest version CRAN offers.
# A cached library must not keep an older DuckDB: CRAN checks against the
# current one, and a stale cache hid exactly that difference before.
cran_flavour_refresh_deps() {
    local log="$1"
    shift
    local ncpus="${DUCKHTS_NCPUS:-2}"
    if ! DUCKHTS_DEPS="$*" DUCKHTS_NCPUS="$ncpus" Rscript --vanilla -e '
        packages <- strsplit(Sys.getenv("DUCKHTS_DEPS"), " ", fixed = TRUE)[[1L]]
        repos <- "https://cloud.r-project.org"
        library_path <- Sys.getenv("R_LIBS_USER")
        installed <- installed.packages(lib.loc = library_path)
        available <- available.packages(repos = repos)
        stale <- vapply(packages, function(package) {
            if (!package %in% rownames(installed)) return(TRUE)
            if (!package %in% rownames(available)) return(FALSE)
            package_version(installed[package, "Version"]) <
                package_version(available[package, "Version"])
        }, logical(1L))
        if (any(stale)) {
            install.packages(packages[stale], lib = library_path, repos = repos,
                             Ncpus = as.integer(Sys.getenv("DUCKHTS_NCPUS")))
        }
        installed <- installed.packages(lib.loc = library_path)
        missing <- setdiff(packages, rownames(installed))
        if (length(missing)) stop("could not install: ", paste(missing, collapse = ", "))
    ' >"$log" 2>&1; then
        tail -n 200 "$log" >&2
        return 1
    fi
    tail -n 20 "$log"
}

# Use the tarball built once for the whole workflow when DUCKHTS_TARBALL names
# one, so every lane checks the same bytes; otherwise build it here. Prints the
# path of the tarball inside $1.
cran_flavour_stage_tarball() {
    local check_root="$1"
    local repo_root="$2"
    local tarball
    if [[ -z "${DUCKHTS_TARBALL:-}" && -n "${DUCKHTS_TARBALL_DIR:-}" ]]; then
        DUCKHTS_TARBALL="$(find "$DUCKHTS_TARBALL_DIR" -maxdepth 1 -name 'Rduckhts_*.tar.gz' -print -quit)"
    fi
    if [[ -n "${DUCKHTS_TARBALL:-}" ]]; then
        test -f "$DUCKHTS_TARBALL"
        if [[ -n "${DUCKHTS_TARBALL_SHA256:-}" ]]; then
            echo "${DUCKHTS_TARBALL_SHA256}  ${DUCKHTS_TARBALL}" | sha256sum -c - >&2
        fi
        cp "$DUCKHTS_TARBALL" "$check_root/"
    else
        (cd "$check_root" && R CMD build --no-manual "$repo_root/r/Rduckhts") >&2
    fi
    tarball="$(find "$check_root" -maxdepth 1 -name 'Rduckhts_*.tar.gz' -print -quit)"
    test -n "$tarball"
    printf '%s\n' "$tarball"
}

# Write the facts that decide whether two runs are comparable.
cran_flavour_receipt() {
    local out="$1"
    local tarball="$2"
    local repo_root="$3"
    local label="$4"
    mkdir -p "$out"
    {
        echo "lane: $label"
        echo "image: ${DUCKHTS_IMAGE_REF:-unknown}"
        echo "commit: $(git -C "$repo_root" rev-parse HEAD 2>/dev/null || echo unknown)"
        echo "tarball: $(basename "$tarball")"
        echo "tarball_sha256: $(sha256sum "$tarball" | cut -d' ' -f1)"
        echo "r_version: $(R --version | head -n 1)"
        echo "r_svn_revision: $(Rscript --vanilla -e 'cat(R.version[["svn rev"]])' 2>/dev/null)"
        echo "long_double: $(Rscript --vanilla -e 'cat(capabilities()[["long.double"]])' 2>/dev/null)"
        echo "cc: $(R CMD config CC)"
        echo "cc_version: $($(R CMD config CC | awk '{print $1}') --version 2>/dev/null | head -n 1)"
        echo "cflags: $(R CMD config CFLAGS)"
        echo "cxx: $(R CMD config CXX)"
        echo "cxxflags: $(R CMD config CXXFLAGS)"
        echo "packages:"
        Rscript --vanilla -e '
            ip <- installed.packages(lib.loc = Sys.getenv("R_LIBS_USER"))
            for (p in intersect(c("DBI", "duckdb", "Rtinycc", "tinytest"), rownames(ip)))
                cat("  ", p, " ", ip[p, "Version"], "\n", sep = "")
        '
    } >"$out/receipt.txt"
    cat "$out/receipt.txt"
}

# Keep the check output, the install log and every HTSlib configure log. The
# package cleanup script removes config.log when installation succeeds, and a
# failed configure is otherwise opaque.
cran_flavour_keep_evidence() {
    local out="$1"
    local check_root="$2"
    mkdir -p "$out/configure"
    if [[ -d "$check_root/Rduckhts.Rcheck" ]]; then
        cp -a "$check_root/Rduckhts.Rcheck" "$out/Rduckhts.Rcheck" 2>/dev/null || true
        grep -E '^Status:' "$check_root/Rduckhts.Rcheck/00check.log" >"$out/status.txt" 2>/dev/null || true
    fi
    find "$check_root" /tmp -maxdepth 8 \( -name config.log -o -name config.status \) \
        -path '*duckhts*' -print0 2>/dev/null |
        while IFS= read -r -d '' file; do
            cp --parents "$file" "$out/configure/" 2>/dev/null || true
        done
}
