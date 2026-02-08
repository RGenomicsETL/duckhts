#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=vendor_common.sh
source "$SCRIPT_DIR/vendor_common.sh"

VERSION="1.23"
TARBALL="samtools-${VERSION}.tar.bz2"
URL="https://github.com/samtools/samtools/releases/download/${VERSION}/${TARBALL}"
SHA256="f228db57d25b724ea26fe55c1c91529f084ef564888865fb190dd87bd04ee74c"

archive_path="$DIST_DIR/$TARBALL"

download_if_missing "$URL" "$archive_path" "$SHA256"

rm -rf "$TP_DIR/samtools" "$TP_DIR/samtools-${VERSION}"
extract_tar_bz2 "$archive_path" "$TP_DIR"

mv "$TP_DIR/samtools-${VERSION}" "$TP_DIR/samtools"
echo "$VERSION" > "$TP_DIR/samtools/VERSION"

# samtools release tarballs bundle an htslib copy; remove to avoid duplication
if [[ -d "$TP_DIR/samtools/htslib" ]]; then
  rm -rf "$TP_DIR/samtools/htslib"
fi

capture_licenses "$TP_DIR/samtools" "samtools"

echo "Vendored samtools ${VERSION} into $TP_DIR/samtools"
