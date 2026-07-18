#!/usr/bin/env bash
# Reproduce CRAN's r-devel-linux-x86_64-fedora-gcc package check.
#
# Local use launches the same R-hub Fedora/GCC container used by CI:
#   scripts/check_cran_fedora.sh
#
# The container path is also the implementation called by GitHub Actions;
# keeping both entry points here prevents the local and CI recipes drifting.
set -euo pipefail

IMAGE="${DUCKHTS_FEDORA_IMAGE:-ghcr.io/r-hub/containers/gcc16:latest}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ "${DUCKHTS_FEDORA_CONTAINER:-0}" != "1" ]]; then
    RLIBS_CACHE="${DUCKHTS_FEDORA_RLIBS:-${HOME}/.cache/duckhts-fedora-r-libs}"
    mkdir -p "$RLIBS_CACHE"
    exec docker run --rm \
        -v "$REPO_ROOT:/repo" \
        -v "$RLIBS_CACHE:/root/R-libs" \
        -w /repo \
        -e DUCKHTS_FEDORA_CONTAINER=1 \
        "$IMAGE" \
        bash scripts/check_cran_fedora.sh
fi

dnf install -y --setopt=install_weak_deps=False \
    bzip2-devel \
    cmake \
    libcurl-devel \
    libdeflate-devel \
    openssl-devel \
    pandoc \
    patch \
    pkgconf-pkg-config \
    xz-devel \
    zlib-devel

mkdir -p /root/R-libs
export R_LIBS_USER=/root/R-libs

DEPS_LOG=/tmp/rduckhts-fedora-dependencies.log
if ! R -q -e '
    packages <- c("DBI", "duckdb", "tinytest")
    installed <- rownames(installed.packages(lib.loc = Sys.getenv("R_LIBS_USER")))
    needed <- setdiff(packages, installed)
    if (length(needed)) {
        install.packages(
            needed,
            lib = Sys.getenv("R_LIBS_USER"),
            repos = "https://cloud.r-project.org"
        )
    }
' >"$DEPS_LOG" 2>&1; then
    tail -n 200 "$DEPS_LOG" >&2
    exit 1
fi
tail -n 20 "$DEPS_LOG"

CHECK_ROOT="$(mktemp -d /tmp/rduckhts-fedora-check.XXXXXX)"
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
if grep -Eq '^Status:.*(ERROR|WARNING)' "$CHECK_LOG"; then
    echo "Fedora GCC 16 check produced an ERROR or WARNING" >&2
    exit 1
fi
exit "$check_status"
