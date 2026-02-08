#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=vendor_common.sh
source "$SCRIPT_DIR/vendor_common.sh"

VERSION="1.23"
TARBALL="bcftools-${VERSION}.tar.bz2"
URL="https://github.com/samtools/bcftools/releases/download/${VERSION}/${TARBALL}"
SHA256="5acde0ac38f7981da1b89d8851a1a425d1c275e1eb76581925c04ca4252c0778"

archive_path="$DIST_DIR/$TARBALL"

download_if_missing "$URL" "$archive_path" "$SHA256"

rm -rf "$TP_DIR/bcftools" "$TP_DIR/bcftools-${VERSION}"
extract_tar_bz2 "$archive_path" "$TP_DIR"

mv "$TP_DIR/bcftools-${VERSION}" "$TP_DIR/bcftools"
echo "$VERSION" > "$TP_DIR/bcftools/VERSION"

# bcftools release tarballs bundle an htslib copy; remove to avoid duplication
if [[ -d "$TP_DIR/bcftools/htslib" ]]; then
  rm -rf "$TP_DIR/bcftools/htslib"
fi

capture_licenses "$TP_DIR/bcftools" "bcftools"

echo "Vendored bcftools ${VERSION} into $TP_DIR/bcftools"
