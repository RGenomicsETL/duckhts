#!/usr/bin/env bash
# Check one source tarball inside a local CRAN-flavour image. The runner,
# cran_flavours_local.sh, mounts the repository at /repo and the tarball at /in.
set -euo pipefail

REPO_ROOT=/repo
git config --global --add safe.directory /repo 2>/dev/null || true
source "$REPO_ROOT/scripts/cran_flavour_common.sh"
LANE="${DUCKHTS_LANE:?lane name}"
DUCKHTS_CHECK_OUT="${DUCKHTS_CHECK_OUT:-/out}"
export DUCKHTS_TARBALL="${DUCKHTS_TARBALL:?tarball}"
mkdir -p /root/R-libs
export R_LIBS_USER=/root/R-libs MAKEFLAGS="-j${DUCKHTS_JOBS:-4}"

# The check controls CRAN's Fedora flavours use that matter to an extension.
export _R_CHECK_COMPILATION_FLAGS_=true _R_CHECK_INSTALL_DEPENDS_=true \
    _R_CHECK_SUGGESTS_ONLY_=true _R_CHECK_NO_RECOMMENDED_=true \
    _R_CHECK_DEPRECATED_DEFUNCT_=true

DUCKHTS_NCPUS="${DUCKHTS_JOBS:-4}" cran_flavour_refresh_deps /tmp/dependencies.log DBI duckdb Rtinycc tinytest
CHECK_ROOT="$(mktemp -d /tmp/check.XXXXXX)"
cd "$CHECK_ROOT"
PACKAGE_TARBALL="$(cran_flavour_stage_tarball "$CHECK_ROOT" "$REPO_ROOT")"
cran_flavour_receipt "$DUCKHTS_CHECK_OUT" "$PACKAGE_TARBALL" "$REPO_ROOT" "$LANE"
cp /tmp/dependencies.log "$DUCKHTS_CHECK_OUT/dependencies.log"
set +e
R CMD check --as-cran --no-manual "$PACKAGE_TARBALL" >"$DUCKHTS_CHECK_OUT/check.out" 2>&1
check_status=$?
set -e
cran_flavour_keep_evidence "$DUCKHTS_CHECK_OUT" "$CHECK_ROOT"
cat "$DUCKHTS_CHECK_OUT/status.txt" 2>/dev/null || true
exit "$check_status"
