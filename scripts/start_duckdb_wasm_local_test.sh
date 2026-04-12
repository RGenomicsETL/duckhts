#!/usr/bin/env bash
set -eu

ROOT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
PORT=${PORT:-8001}
HOST=${HOST:-127.0.0.1}
EMSDK_IMAGE=${EMSDK_IMAGE:-emscripten/emsdk:3.1.74}
LOCAL_WASM_IMAGE=${LOCAL_WASM_IMAGE:-duckhts/duckdb-wasm-local:latest}
LOCAL_WASM_DOCKERFILE=${LOCAL_WASM_DOCKERFILE:-scripts/docker/duckdb-wasm-local.Dockerfile}
DUCKDB_WASM_NPM_VERSION=${DUCKDB_WASM_NPM_VERSION:-1.29.0}
SITE_ROOT=${SITE_ROOT:-${ROOT_DIR}/duckdb-wasm-local-site}
DOCKER_WORK_ROOT=${DOCKER_WORK_ROOT:-${ROOT_DIR}/.duckdb_wasm_docker_work}
VCPKG_TOOLCHAIN_PATH=${VCPKG_TOOLCHAIN_PATH:-${ROOT_DIR}/local_vcpkg_installation/scripts/buildsystems/vcpkg.cmake}
VCPKG_TARGET_TRIPLET=${VCPKG_TARGET_TRIPLET:-wasm32-emscripten}
VCPKG_HOST_TRIPLET=${VCPKG_HOST_TRIPLET:-x64-linux}
VCPKG_OVERLAY_PORTS=${VCPKG_OVERLAY_PORTS:-${ROOT_DIR}/extension-ci-tools/vcpkg_ports}
VCPKG_OVERLAY_TRIPLETS=${VCPKG_OVERLAY_TRIPLETS:-${ROOT_DIR}/extension-ci-tools/toolchains}

echo "Building DuckHTS wasm extension (DuckDB community extension path)..."

if command -v emcmake >/dev/null 2>&1; then
  if [ ! -f "${VCPKG_TOOLCHAIN_PATH}" ]; then
    echo "Bootstrapping local vcpkg (missing ${VCPKG_TOOLCHAIN_PATH})..."
    rm -rf "${ROOT_DIR}/local_vcpkg_installation"
    mkdir -p "${ROOT_DIR}/local_vcpkg_installation"
    git -C "${ROOT_DIR}/local_vcpkg_installation" init
    git -C "${ROOT_DIR}/local_vcpkg_installation" remote add origin https://github.com/microsoft/vcpkg.git
    git -C "${ROOT_DIR}/local_vcpkg_installation" fetch --depth=1 origin 84bab45d415d22042bd0b9081aea57f362da3f35
    git -C "${ROOT_DIR}/local_vcpkg_installation" checkout FETCH_HEAD
    "${ROOT_DIR}/local_vcpkg_installation/bootstrap-vcpkg.sh"
  fi

  export VCPKG_TOOLCHAIN_PATH
  export VCPKG_TARGET_TRIPLET
  export VCPKG_HOST_TRIPLET
  export VCPKG_OVERLAY_PORTS
  export VCPKG_OVERLAY_TRIPLETS
  make -C "${ROOT_DIR}" DUCKDB_PLATFORM=wasm_eh configure release move_wasm_extension
else
  echo "emcmake not found on host; using local DuckDB-wasm build image (${LOCAL_WASM_IMAGE})"
  docker build -f "${ROOT_DIR}/${LOCAL_WASM_DOCKERFILE}" -t "${LOCAL_WASM_IMAGE}" "${ROOT_DIR}"
  mkdir -p "${DOCKER_WORK_ROOT}"
  rsync -a --delete \
    --exclude '.git/' \
    --exclude 'build/' \
    --exclude 'cmake_build/' \
    --exclude 'duckdb-wasm-local-site/' \
    --exclude '.duckdb_wasm_docker_work/' \
    "${ROOT_DIR}/" "${DOCKER_WORK_ROOT}/"
  docker run --rm -v "${DOCKER_WORK_ROOT}:/work" -w /work "${LOCAL_WASM_IMAGE}" bash -lc '
set -eu
source /opt/emsdk/emsdk_env.sh
export VCPKG_TOOLCHAIN_PATH=/opt/vcpkg/scripts/buildsystems/vcpkg.cmake
export VCPKG_TARGET_TRIPLET=wasm32-emscripten
export VCPKG_HOST_TRIPLET=x64-linux
export VCPKG_OVERLAY_PORTS=/work/extension-ci-tools/vcpkg_ports
export VCPKG_OVERLAY_TRIPLETS=/work/extension-ci-tools/toolchains
rm -rf cmake_build
make DUCKDB_PLATFORM=wasm_eh configure release move_wasm_extension
'
  mkdir -p "${ROOT_DIR}/build/wasm_eh/extension/duckhts"
  cp -f "${DOCKER_WORK_ROOT}/build/wasm_eh/extension/duckhts/duckhts.duckdb_extension.wasm" \
        "${ROOT_DIR}/build/wasm_eh/extension/duckhts/duckhts.duckdb_extension.wasm"
fi

WASM_EXT="${ROOT_DIR}/build/wasm_eh/extension/duckhts/duckhts.duckdb_extension.wasm"
if [ ! -f "${WASM_EXT}" ]; then
  echo "Expected wasm extension not found: ${WASM_EXT}" >&2
  exit 1
fi

rm -rf "${SITE_ROOT}"
mkdir -p "${SITE_ROOT}/scripts" "${SITE_ROOT}/duckdb-wasm" "${SITE_ROOT}/extdata"

# Host duckdb-wasm runtime locally to avoid cross-origin worker loading issues
mkdir -p "${SITE_ROOT}/duckdb-wasm-runtime"
RUNTIME_BASE="https://cdn.jsdelivr.net/npm/@duckdb/duckdb-wasm@${DUCKDB_WASM_NPM_VERSION}/dist"
curl -fsSL "${RUNTIME_BASE}/duckdb-browser.mjs" -o "${SITE_ROOT}/duckdb-wasm-runtime/duckdb-browser.mjs"
curl -fsSL "${RUNTIME_BASE}/duckdb-browser-eh.worker.js" -o "${SITE_ROOT}/duckdb-wasm-runtime/duckdb-browser-eh.worker.js"
curl -fsSL "${RUNTIME_BASE}/duckdb-eh.wasm" -o "${SITE_ROOT}/duckdb-wasm-runtime/duckdb-eh.wasm"

cp -f "${ROOT_DIR}/scripts/duckdb-wasm-local-test.html" "${SITE_ROOT}/scripts/duckdb-wasm-local-test.html"
cp -f "${WASM_EXT}" "${SITE_ROOT}/duckdb-wasm/duckhts.duckdb_extension.wasm"

# Expose local test data over HTTP for real browser http(s) reads
cp -f "${ROOT_DIR}/r/Rduckhts/inst/extdata/r1.fq" "${SITE_ROOT}/extdata/"
cp -f "${ROOT_DIR}/r/Rduckhts/inst/extdata/header_tabix.tsv.gz" "${SITE_ROOT}/extdata/"
cp -f "${ROOT_DIR}/r/Rduckhts/inst/extdata/header_tabix.tsv.gz.tbi" "${SITE_ROOT}/extdata/"

cat <<EOF
Built:
  ${WASM_EXT}

Open in a browser:
  http://${HOST}:${PORT}/scripts/duckdb-wasm-local-test.html

Served files:
  /duckdb-wasm-runtime/duckdb-browser.mjs
  /duckdb-wasm-runtime/duckdb-browser-eh.worker.js
  /duckdb-wasm-runtime/duckdb-eh.wasm
  /duckdb-wasm/duckhts.duckdb_extension.wasm
  /extdata/r1.fq
  /extdata/header_tabix.tsv.gz
  /extdata/header_tabix.tsv.gz.tbi
EOF

cd "${SITE_ROOT}"
exec python3 -m http.server "${PORT}" --bind "${HOST}"
