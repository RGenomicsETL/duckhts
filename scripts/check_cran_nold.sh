#!/usr/bin/env bash
# Reproduce CRAN's R-devel check with long-double support disabled.
set -euo pipefail

IMAGE="ghcr.io/r-hub/containers/nold@sha256:089518d59395f2f36c357e0fd2325d35c9066438fbaf0854d2204e8603415476"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$REPO_ROOT/scripts/cran_flavour_common.sh"
DUCKHTS_CHECK_OUT="${DUCKHTS_CHECK_OUT:-$REPO_ROOT/check-artifacts/nold}"
export DUCKHTS_IMAGE_REF="$IMAGE"

if [[ "${DUCKHTS_NOLD_CONTAINER:-0}" != "1" ]]; then
    RLIBS_CACHE="${DUCKHTS_NOLD_RLIBS:-${HOME}/.cache/duckhts-nold-r-libs}"
    mkdir -p "$RLIBS_CACHE"
    # Pass a pre-built tarball and an output directory through to the container.
    extra=()
    if [[ -n "${DUCKHTS_TARBALL:-}" ]]; then
        extra+=(-v "$(cd "$(dirname "$DUCKHTS_TARBALL")" && pwd):/in:ro"
                -e "DUCKHTS_TARBALL=/in/$(basename "$DUCKHTS_TARBALL")"
                -e DUCKHTS_TARBALL_SHA256)
    fi
    if [[ -n "${DUCKHTS_CHECK_OUT:-}" ]]; then
        mkdir -p "$DUCKHTS_CHECK_OUT"
        extra+=(-v "$(cd "$DUCKHTS_CHECK_OUT" && pwd):/out" -e DUCKHTS_CHECK_OUT=/out)
    fi
    exec docker run --rm "${extra[@]}" \
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
DUCKHTS_NCPUS=2 cran_flavour_refresh_deps "$DEPS_LOG" DBI duckdb Rtinycc tinytest

CHECK_ROOT="$(mktemp -d /tmp/rduckhts-nold-check.XXXXXX)"
cd "$CHECK_ROOT"
PACKAGE_TARBALL="$(cran_flavour_stage_tarball "$CHECK_ROOT" "$REPO_ROOT")"
cran_flavour_receipt "$DUCKHTS_CHECK_OUT" "$PACKAGE_TARBALL" "$REPO_ROOT" "nold"
cp "$DEPS_LOG" "$DUCKHTS_CHECK_OUT/dependencies.log"
set +e
R CMD check --as-cran --no-manual "$PACKAGE_TARBALL"
check_status=$?
set -e
cran_flavour_keep_evidence "$DUCKHTS_CHECK_OUT" "$CHECK_ROOT"

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
