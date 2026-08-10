#!/usr/bin/env bash
set -eu

ROOT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
# shellcheck source=duckhts_cache.sh
source "${ROOT_DIR}/scripts/duckhts_cache.sh"
PORT=${PORT:-8000}
HOST=${HOST:-127.0.0.1}
WEBR_IMAGE=${WEBR_IMAGE:-ghcr.io/r-wasm/webr:main}
ARTIFACT_ROOT=${ARTIFACT_ROOT:-$(duckhts_cache_subdir webr/local-artifacts)}

echo "Building local webR artifact with ${WEBR_IMAGE}..."

docker run --rm -v "${ROOT_DIR}:/work" -w /work/r/Rduckhts "${WEBR_IMAGE}" bash -lc '
set -eu

export THREADS=5

Rscript bootstrap.R /work

mkdir -p /tmp/Rduckhts-src
cp -a . /tmp/Rduckhts-src
rm -rf /tmp/Rduckhts-src/.github \
       /tmp/Rduckhts-src/Rduckhts.Rcheck \
       /tmp/Rduckhts-src/..Rcheck \
       /tmp/Rduckhts-src/README.html

R CMD build /tmp/Rduckhts-src

Rscript -e "if (!requireNamespace(\"pak\", quietly = TRUE)) install.packages(\"pak\", repos = \"https://repo.r-wasm.org/\"); pak::pak(\"r-wasm/rwasm\", ask = FALSE)"
Rscript -e "ver <- read.dcf(\"DESCRIPTION\")[1, \"Version\"]; rwasm::build(sprintf(\"./Rduckhts_%s.tar.gz\", ver))"
'

VERSION=$(Rscript -e 'cat(read.dcf("r/Rduckhts/DESCRIPTION")[1, "Version"])')
BUILD_TGZ_PATH="${ROOT_DIR}/r/Rduckhts/Rduckhts_${VERSION}.tgz"
TGZ_PATH="${ARTIFACT_ROOT}/Rduckhts_${VERSION}.tgz"
REPO_ROOT="${ARTIFACT_ROOT}/repo"
LOCAL_REPO_DIR="${REPO_ROOT}/src/contrib"
SITE_ROOT="${ARTIFACT_ROOT}/site"
WEBR_R_SERIES=$(
  docker run --rm "${WEBR_IMAGE}" bash -lc \
    "Rscript -e 'cat(paste(R.version\$major, sub(\"\\\\..*$\", \"\", R.version\$minor), sep = \".\"))'"
)
BINARY_REPO_DIR="${REPO_ROOT}/bin/emscripten/contrib/${WEBR_R_SERIES}"

if [ ! -f "${BUILD_TGZ_PATH}" ]; then
  echo "Expected wasm package artifact not found: ${BUILD_TGZ_PATH}" >&2
  exit 1
fi

mkdir -p "${ARTIFACT_ROOT}"
cp -f "${BUILD_TGZ_PATH}" "${TGZ_PATH}"

rm -rf "${REPO_ROOT}"
mkdir -p "${LOCAL_REPO_DIR}"
mkdir -p "${BINARY_REPO_DIR}"
cp -f "${TGZ_PATH}" "${LOCAL_REPO_DIR}/"
cp -f "${TGZ_PATH}" "${BINARY_REPO_DIR}/"
Rscript -e "
desc <- read.dcf('r/Rduckhts/DESCRIPTION')
fields <- as.list(desc[1, , drop = TRUE])
fields\$File <- sprintf('Rduckhts_%s.tgz', desc[1, 'Version'])
db <- as.data.frame(fields, stringsAsFactors = FALSE, check.names = FALSE)
mat <- as.matrix(db)
rownames(mat) <- mat[, 'Package']
repo_src <- '${LOCAL_REPO_DIR}'
repo_bin <- '${BINARY_REPO_DIR}'
write_index <- function(repo) {
  write.dcf(db, file = file.path(repo, 'PACKAGES'))
  con <- gzfile(file.path(repo, 'PACKAGES.gz'), 'wb')
  writeLines(readLines(file.path(repo, 'PACKAGES')), con)
  close(con)
  saveRDS(mat, file.path(repo, 'PACKAGES.rds'))
}
write_index(repo_src)
write_index(repo_bin)
"

rm -rf "${SITE_ROOT}"
mkdir -p "${SITE_ROOT}/scripts" "${SITE_ROOT}/r/Rduckhts" "${SITE_ROOT}/extdata"
cp -f "${ROOT_DIR}/scripts/webr-local-test.html" "${SITE_ROOT}/scripts/webr-local-test.html"
cp -f "${ROOT_DIR}/r/Rduckhts/DESCRIPTION" "${SITE_ROOT}/r/Rduckhts/DESCRIPTION"
cp -a "${REPO_ROOT}" "${SITE_ROOT}/r/Rduckhts/webr-repo"
# Expose tabix-indexed test files under /extdata so the browser REPL can open
# them via a real http:// URL without needing a remote server or ws-proxy.
cp -f "${ROOT_DIR}/r/Rduckhts/inst/extdata/header_tabix.tsv.gz" "${SITE_ROOT}/extdata/"
cp -f "${ROOT_DIR}/r/Rduckhts/inst/extdata/header_tabix.tsv.gz.tbi" "${SITE_ROOT}/extdata/"
cp -f "${ROOT_DIR}/r/Rduckhts/inst/extdata/gff_file.gff.gz" "${SITE_ROOT}/extdata/"
cp -f "${ROOT_DIR}/r/Rduckhts/inst/extdata/gff_file.gff.gz.tbi" "${SITE_ROOT}/extdata/"

cat <<EOF
Built:
  ${TGZ_PATH}

Local webR repo:
  ${REPO_ROOT}

Open in a browser:
  http://${HOST}:${PORT}/scripts/webr-local-test.html

The page will load:
  /r/Rduckhts/DESCRIPTION
  /r/Rduckhts/webr-repo/src/contrib/PACKAGES
  /r/Rduckhts/webr-repo/bin/emscripten/contrib/${WEBR_R_SERIES}/PACKAGES
  /r/Rduckhts/webr-repo/bin/emscripten/contrib/${WEBR_R_SERIES}/Rduckhts_${VERSION}.tgz

Artifact root:
  ${ARTIFACT_ROOT}
EOF

cd "${SITE_ROOT}"
exec python3 -m http.server "${PORT}" --bind "${HOST}"
