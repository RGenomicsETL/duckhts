#!/usr/bin/env bash
set -eu

ROOT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
# shellcheck source=duckhts_cache.sh
source "${ROOT_DIR}/scripts/duckhts_cache.sh"
PORT=${PORT:-8001}
HOST=${HOST:-127.0.0.1}
LOCAL_WASM_IMAGE=${LOCAL_WASM_IMAGE:-duckhts/duckdb-wasm-local:latest}
LOCAL_WASM_DOCKERFILE=${LOCAL_WASM_DOCKERFILE:-scripts/docker/duckdb-wasm-local.Dockerfile}
# duckdb-wasm 1.29.0 is based on DuckDB v1.1.1 and traps while loading our
# current extension build; use the first 1.4.x-aligned runtime by default.
DUCKDB_WASM_NPM_VERSION=${DUCKDB_WASM_NPM_VERSION:-1.31.0}
# SERVE=1 (default) builds the site then serves it with python http.server.
# SERVE=0 builds the site and exits, printing SITE_ROOT/PORT so an external
# driver (e.g. Playwright's webServer) can own the HTTP server lifecycle.
SERVE=${SERVE:-1}
ARTIFACT_ROOT=${ARTIFACT_ROOT:-$(duckhts_cache_subdir wasm/local-artifacts)}
SITE_ROOT=${SITE_ROOT:-${ARTIFACT_ROOT}/site}
WASM_BUILD_DIR=${WASM_BUILD_DIR:-${ARTIFACT_ROOT}/build/wasm_eh/extension/duckhts}
DOCKER_WORK_ROOT=${DOCKER_WORK_ROOT:-$(duckhts_cache_subdir wasm/docker-work)}
DOCKER_REBUILD_IMAGE=${DOCKER_REBUILD_IMAGE:-0}

echo "Building DuckHTS wasm extension (DuckDB community extension path)..."

if [ "${DOCKER_REBUILD_IMAGE}" = "1" ] || ! docker image inspect "${LOCAL_WASM_IMAGE}" >/dev/null 2>&1; then
  docker build -f "${ROOT_DIR}/${LOCAL_WASM_DOCKERFILE}" -t "${LOCAL_WASM_IMAGE}" "${ROOT_DIR}"
fi

mkdir -p "${DOCKER_WORK_ROOT}"
rsync -a --delete \
  --exclude '.git/' \
  --exclude 'build/' \
  --exclude 'cmake_build/' \
  --exclude 'configure/venv/' \
  --exclude 'duckdb-wasm-local-site/' \
  --exclude '.duckdb-wasm-local-artifacts/' \
  --exclude '.webr-local-artifacts/' \
  --exclude '.duckdb_wasm_docker_work/' \
  "${ROOT_DIR}/" "${DOCKER_WORK_ROOT}/"
docker run --rm -v "${DOCKER_WORK_ROOT}:/work" -w /work \
  "${LOCAL_WASM_IMAGE}" bash -lc '
set -eu
source /opt/emsdk/emsdk_env.sh
export VCPKG_TOOLCHAIN_PATH=/opt/vcpkg/scripts/buildsystems/vcpkg.cmake
export VCPKG_TARGET_TRIPLET=wasm32-emscripten
export VCPKG_HOST_TRIPLET=x64-linux
export VCPKG_OVERLAY_PORTS=/work/extension-ci-tools/vcpkg_ports
export VCPKG_OVERLAY_TRIPLETS=/work/extension-ci-tools/toolchains
# Avoid mixed-arch htslib artifacts from host copies contaminating wasm build.
# cmake_build/ is preserved across runs so vcpkg packages are not re-downloaded.
find third_party/htslib -type f \( \
  -name "*.o" -o -name "*.a" -o -name "*.pico" -o -name "*.lo" -o -name "*.la" -o \
  -name "*.so" -o -name "*.dylib" -o -name "*.dll" -o \
  -name "config.h" -o -name "config.mk" -o -name "config_vars.h" -o \
  -name "config.status" -o -name "config.log" -o -name "config.cache" -o \
  -name "htslib.pc.tmp" -o -name "a.wasm" \
\) -delete || true
# Always build the venv fresh inside the container: a host configure/venv (its
# python symlinks point at host paths absent here) is excluded from the rsync,
# but rsync --delete cannot remove an excluded path, so a stale copy from a
# prior run could shadow recreation.  Removing it lets `make configure` build a
# working venv with the container python.
rm -rf configure/venv
make DUCKDB_PLATFORM=wasm_eh configure
# cmake_build is retained for fast local rebuilds.  Invalidate the htslib
# ExternalProject completion stamps explicitly: otherwise a host-side
# config.status copied by rsync can regenerate a curl-enabled wasm htslib after
# the configure step is skipped.  duckdb-wasm exports no libcurl functions, so
# that stale mix produces a side module which compiles but fails during LOAD.
rm -f cmake_build/release/htslib_build-prefix/src/htslib_build-stamp/htslib_build-configure
rm -f cmake_build/release/htslib_build-prefix/src/htslib_build-stamp/htslib_build-build
make DUCKDB_PLATFORM=wasm_eh release move_wasm_extension
'

mkdir -p "${WASM_BUILD_DIR}"
cp -f "${DOCKER_WORK_ROOT}/build/wasm_eh/extension/duckhts/duckhts.duckdb_extension.wasm" \
      "${WASM_BUILD_DIR}/duckhts.duckdb_extension.wasm"

WASM_EXT="${WASM_BUILD_DIR}/duckhts.duckdb_extension.wasm"
if [ ! -f "${WASM_EXT}" ]; then
  echo "Expected wasm extension not found: ${WASM_EXT}" >&2
  exit 1
fi

rm -rf "${SITE_ROOT}"
mkdir -p "${SITE_ROOT}/scripts" "${SITE_ROOT}/duckdb-wasm" "${SITE_ROOT}/extdata"

RUNTIME_BASE="https://cdn.jsdelivr.net/npm/@duckdb/duckdb-wasm@${DUCKDB_WASM_NPM_VERSION}/dist"
RUNTIME_CACHE="$(duckhts_cache_subdir "runtime/duckdb-wasm/${DUCKDB_WASM_NPM_VERSION}")"
if [[ "$DUCKDB_WASM_NPM_VERSION" != "1.31.0" ]]; then
  echo "DuckDB wasm runtime version $DUCKDB_WASM_NPM_VERSION has no pinned asset manifest" >&2
  exit 2
fi
mkdir -p "$RUNTIME_CACHE"
runtime_sha256() { # filename
  case "$1" in
    duckdb-browser.mjs) printf '%s\n' '660ee2979e878cfd1d70122ab6d511d4454c5e687f59977e77ff70097730a8a0' ;;
    duckdb-browser-eh.worker.js) printf '%s\n' 'fb692cd56e87c71849ff545e14fef54c91ed4cdef295f16172ab27be8de76b5d' ;;
    duckdb-eh.wasm) printf '%s\n' '07993a5cda534ebb303d476cbdf3d1f7271841c1298709a4b4a5713d8c78b156' ;;
    duckdb-browser-eh.worker.js.map) printf '%s\n' 'e3750fb2ea26e4e1e1baf49f3a2b2f3d722d3a0bd840830f2abd431f8546dc90' ;;
    *) return 2 ;;
  esac
}
verify_runtime() { # filename path
  [[ "$(sha256sum "$2" | awk '{print $1}')" == "$(runtime_sha256 "$1")" ]]
}
fetch_runtime() { # filename [optional]
  local name="$1"
  local optional="${2:-0}"
  local cached="$RUNTIME_CACHE/$name"
  local temporary="${cached}.partial.$$"
  if [[ ! -s "$cached" ]] || ! verify_runtime "$name" "$cached"; then
    rm -f "$cached" "$temporary"
    if ! curl -fsSL "${RUNTIME_BASE}/$name" -o "$temporary" || ! verify_runtime "$name" "$temporary"; then
      rm -f "$temporary"
      if [[ "$optional" = "1" ]]; then
        return 0
      fi
      echo "could not stage verified duckdb-wasm runtime asset: $name" >&2
      return 1
    fi
    mv -f "$temporary" "$cached"
  fi
  cp -f "$cached" "$SITE_ROOT/$name"
}
fetch_runtime duckdb-browser.mjs
fetch_runtime duckdb-browser-eh.worker.js
fetch_runtime duckdb-eh.wasm
fetch_runtime duckdb-browser-eh.worker.js.map 1
printf '%s\n' '<!doctype html><title>duckhts</title>' > "${SITE_ROOT}/favicon.ico"

ln -sf "${ROOT_DIR}/scripts/duckdb-wasm-local-test.html" "${SITE_ROOT}/scripts/duckdb-wasm-local-test.html"
cp -f "${WASM_EXT}" "${SITE_ROOT}/duckdb-wasm/duckhts.duckdb_extension.wasm"

# Expose local test data over HTTP for real browser http(s) reads
cp -f "${ROOT_DIR}/r/Rduckhts/inst/extdata/r1.fq" "${SITE_ROOT}/extdata/"
cp -f "${ROOT_DIR}/r/Rduckhts/inst/extdata/header_tabix.tsv.gz" "${SITE_ROOT}/extdata/"
cp -f "${ROOT_DIR}/r/Rduckhts/inst/extdata/header_tabix.tsv.gz.tbi" "${SITE_ROOT}/extdata/"
cp -f "${ROOT_DIR}/r/Rduckhts/inst/extdata/fixture_mixed.bam" "${SITE_ROOT}/extdata/"

cat <<EOF
Built:
  ${WASM_EXT}

Open in a browser:
  http://${HOST}:${PORT}/scripts/duckdb-wasm-local-test.html

Artifact roots:
  ${ARTIFACT_ROOT}
  ${DOCKER_WORK_ROOT}

Served files:
  /duckdb-browser.mjs
  /duckdb-browser-eh.worker.js
  /duckdb-eh.wasm
  /duckdb-wasm/duckhts.duckdb_extension.wasm
  /extdata/r1.fq
  /extdata/header_tabix.tsv.gz
  /extdata/header_tabix.tsv.gz.tbi
  /extdata/fixture_mixed.bam
EOF

if [ "${SERVE}" = "0" ]; then
  echo "SERVE=0: site built, not starting HTTP server."
  echo "DUCKHTS_WASM_SITE_ROOT=${SITE_ROOT}"
  echo "DUCKHTS_WASM_PORT=${PORT}"
  echo "DUCKHTS_WASM_HOST=${HOST}"
  exit 0
fi

cd "${SITE_ROOT}"
exec python3 -m http.server "${PORT}" --bind "${HOST}"
