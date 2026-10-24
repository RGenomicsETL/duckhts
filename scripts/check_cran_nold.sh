#!/usr/bin/env bash
# Reproduce CRAN's R-devel check with long-double support disabled.
set -euo pipefail

IMAGE="ghcr.io/r-hub/containers/nold@sha256:089518d59395f2f36c357e0fd2325d35c9066438fbaf0854d2204e8603415476"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ "${DUCKHTS_NOLD_CONTAINER:-0}" != "1" ]]; then
    RLIBS_CACHE="${DUCKHTS_NOLD_RLIBS:-${HOME}/.cache/duckhts-nold-r-libs}"
    mkdir -p "$RLIBS_CACHE"
    exec docker run --rm \
        -v "$REPO_ROOT:/repo" \
        -v "$RLIBS_CACHE:/root/R-libs" \
        -w /repo \
        -e DUCKHTS_NOLD_CONTAINER=1 \
        "$IMAGE" \
        bash scripts/check_cran_nold.sh
fi

export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    libbz2-dev \
    libcurl4-openssl-dev \
    libdeflate-dev \
    liblzma-dev \
    libssl-dev \
    pandoc \
    pkg-config \
    qpdf \
    zlib1g-dev

Rscript --vanilla -e '
    stopifnot(!capabilities()[["long.double"]])
    print(R.version.string)
    print(capabilities()["long.double"])
    cat("decimal predecessor: ", sprintf("%.17g", 0.49999999999999994), "\n", sep = "")
'

mkdir -p /root/R-libs
export R_LIBS_USER=/root/R-libs
export MAKEFLAGS=-j2

DEPS_LOG=/tmp/rduckhts-nold-dependencies.log
if ! Rscript --vanilla -e '
    packages <- c("DBI", "duckdb", "Rtinycc", "tinytest")
    installed <- rownames(installed.packages(lib.loc = Sys.getenv("R_LIBS_USER")))
    needed <- setdiff(packages, installed)
    if (length(needed)) {
        install.packages(
            needed,
            lib = Sys.getenv("R_LIBS_USER"),
            repos = "https://cloud.r-project.org",
            Ncpus = 2L
        )
    }
' >"$DEPS_LOG" 2>&1; then
    tail -n 200 "$DEPS_LOG" >&2
    exit 1
fi
tail -n 20 "$DEPS_LOG"

CHECK_ROOT="$(mktemp -d /tmp/rduckhts-nold-check.XXXXXX)"
cd "$CHECK_ROOT"
R CMD build --no-manual "$REPO_ROOT/r/Rduckhts"
PACKAGE_TARBALL="$(find "$CHECK_ROOT" -maxdepth 1 -name 'Rduckhts_*.tar.gz' -print -quit)"
test -n "$PACKAGE_TARBALL"
set +e
R CMD check --as-cran --no-manual "$PACKAGE_TARBALL"
check_status=$?
set -e

CHECK_LOG="$CHECK_ROOT/Rduckhts.Rcheck/00check.log"
test -f "$CHECK_LOG"
if [[ "$check_status" -ne 0 ]]; then
    INSTALL_LOG="$CHECK_ROOT/Rduckhts.Rcheck/00install.out"
    if [[ -f "$INSTALL_LOG" ]]; then
        echo "--- $INSTALL_LOG (last 240 lines) ---" >&2
        tail -n 240 "$INSTALL_LOG" >&2
    fi
fi
if grep -Eq '^Status:.*(ERROR|WARNING)' "$CHECK_LOG"; then
    echo "R-devel no-long-double check produced an ERROR or WARNING" >&2
    exit 1
fi
exit "$check_status"
