#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=vendor_common.sh
source "$SCRIPT_DIR/vendor_common.sh"

VERSION="1.23"
TARBALL="htslib-${VERSION}.tar.bz2"
URL="https://github.com/samtools/htslib/releases/download/${VERSION}/${TARBALL}"
SHA256="63927199ef9cea03096345b95d96cb600ae10385248b2ef670b0496c2ab7e4cd"

archive_path="$DIST_DIR/$TARBALL"

download_if_missing "$URL" "$archive_path" "$SHA256"

rm -rf "$TP_DIR/htslib" "$TP_DIR/htslib-${VERSION}"
extract_tar_bz2 "$archive_path" "$TP_DIR"

mv "$TP_DIR/htslib-${VERSION}" "$TP_DIR/htslib"
echo "$VERSION" > "$TP_DIR/htslib/VERSION"

capture_licenses "$TP_DIR/htslib" "htslib"

echo "Vendored HTSlib ${VERSION} into $TP_DIR/htslib"
