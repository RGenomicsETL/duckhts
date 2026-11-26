#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=vendor_common.sh
source "$SCRIPT_DIR/vendor_common.sh"

VERSION="1.25"
TARBALL="htslib-${VERSION}.tar.bz2"
URL="https://github.com/samtools/htslib/releases/download/${VERSION}/${TARBALL}"
SHA256="158c271e4e2986bb2dd11283b640edc10c56695b293df77a109ce99480303730"

archive_path="$DIST_DIR/$TARBALL"

download_if_missing "$URL" "$archive_path" "$SHA256"

rm -rf "$TP_DIR/htslib" "$TP_DIR/htslib-${VERSION}"
extract_tar_bz2 "$archive_path" "$TP_DIR"

mv "$TP_DIR/htslib-${VERSION}" "$TP_DIR/htslib"
echo "$VERSION" > "$TP_DIR/htslib/VERSION"

# The release contributor guide is not a build, runtime, attribution, or API
# input, and its upstream trailing whitespace conflicts with this repository's
# mandatory diff check. Keep the vendored distribution focused on owned inputs.
rm -f "$TP_DIR/htslib/CONTRIBUTING.md"

apply_patches "$TP_DIR/htslib" "htslib"
capture_licenses "$TP_DIR/htslib" "htslib"

echo "Vendored HTSlib ${VERSION} into $TP_DIR/htslib"
