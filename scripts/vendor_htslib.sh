#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=vendor_common.sh
source "$SCRIPT_DIR/vendor_common.sh"

VERSION="1.23.1"
TARBALL="htslib-${VERSION}.tar.bz2"
URL="https://github.com/samtools/htslib/releases/download/${VERSION}/${TARBALL}"
SHA256="f8a3f36effeec38f043c53ab1f2d9ed45064f14205c5ef8e3c815763b90803c4"

archive_path="$DIST_DIR/$TARBALL"

download_if_missing "$URL" "$archive_path" "$SHA256"

rm -rf "$TP_DIR/htslib" "$TP_DIR/htslib-${VERSION}"
extract_tar_bz2 "$archive_path" "$TP_DIR"

mv "$TP_DIR/htslib-${VERSION}" "$TP_DIR/htslib"
echo "$VERSION" > "$TP_DIR/htslib/VERSION"

capture_licenses "$TP_DIR/htslib" "htslib"

echo "Vendored HTSlib ${VERSION} into $TP_DIR/htslib"
