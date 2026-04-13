#!/usr/bin/env bash
set -eu

ROOT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
PORT=${PORT:-8001}
HOST=${HOST:-127.0.0.1}
LOCAL_WASM_IMAGE=${LOCAL_WASM_IMAGE:-duckhts/duckdb-wasm-local:latest}
LOCAL_WASM_DOCKERFILE=${LOCAL_WASM_DOCKERFILE:-scripts/docker/duckdb-wasm-local.Dockerfile}
DUCKDB_WASM_NPM_VERSION=${DUCKDB_WASM_NPM_VERSION:-1.29.0}
ARTIFACT_ROOT=${ARTIFACT_ROOT:-${ROOT_DIR}/.duckdb-wasm-local-artifacts}
SITE_ROOT=${SITE_ROOT:-${ARTIFACT_ROOT}/site}
WASM_BUILD_DIR=${WASM_BUILD_DIR:-${ARTIFACT_ROOT}/build/wasm_eh/extension/duckhts}
DOCKER_WORK_ROOT=${DOCKER_WORK_ROOT:-${ROOT_DIR}/.duckdb_wasm_docker_work}
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
  -name "htslib.pc.tmp" -o -name "a.wasm" \
\) -delete || true
make DUCKDB_PLATFORM=wasm_eh configure release move_wasm_extension
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
curl -fsSL "${RUNTIME_BASE}/duckdb-browser.mjs" -o "${SITE_ROOT}/duckdb-browser.mjs"
curl -fsSL "${RUNTIME_BASE}/duckdb-browser-eh.worker.js" -o "${SITE_ROOT}/duckdb-browser-eh.worker.js"
curl -fsSL "${RUNTIME_BASE}/duckdb-eh.wasm" -o "${SITE_ROOT}/duckdb-eh.wasm"
curl -fsSL "${RUNTIME_BASE}/duckdb-browser-eh.worker.js.map" -o "${SITE_ROOT}/duckdb-browser-eh.worker.js.map" || true
printf '%s\n' '<!doctype html><title>duckhts</title>' > "${SITE_ROOT}/favicon.ico"

ln -sf "${ROOT_DIR}/scripts/duckdb-wasm-local-test.html" "${SITE_ROOT}/scripts/duckdb-wasm-local-test.html"
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
EOF

cd "${SITE_ROOT}"
exec python3 -m http.server "${PORT}" --bind "${HOST}"
