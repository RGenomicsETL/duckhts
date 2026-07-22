#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=vendor_common.sh
source "$SCRIPT_DIR/vendor_common.sh"

LIBBIGWIG_COMMIT="43c294ef1721a73b760803ca5e9410d581b98f17"
LIBBIGWIG_ARCHIVE="libBigWig-${LIBBIGWIG_COMMIT}.tar.gz"
LIBBIGWIG_URL="https://github.com/dpryan79/libBigWig/archive/${LIBBIGWIG_COMMIT}.tar.gz"
LIBBIGWIG_SHA256="4d1311f4ca8eba81c9ded3042667bd02020548c1b538fb165ecea94a9f2a7ed3"
ARCHIVE="$DIST_DIR/$LIBBIGWIG_ARCHIVE"
TARGET_DIR="$TP_DIR/libBigWig"

download_if_missing "$LIBBIGWIG_URL" "$ARCHIVE" "$LIBBIGWIG_SHA256"

tmp="$(mktemp -d "${TMPDIR:-/tmp}/duckhts-libbigwig.XXXXXX")"
trap 'rm -rf "$tmp"' EXIT HUP INT TERM
tar -xzf "$ARCHIVE" -C "$tmp"
source_dir="$tmp/libBigWig-${LIBBIGWIG_COMMIT}"

reset_dir "$TARGET_DIR"
for file in LICENSE README.md bigWig.h bigWigIO.h bwCommon.h bwValues.h \
    bwRead.c bwStats.c bwValues.c bwWrite.c io.c; do
  cp "$source_dir/$file" "$TARGET_DIR/$file"
done
mkdir -p "$TARGET_DIR/test"
cp "$source_dir/test/test.bw" "$TARGET_DIR/test/test.bw"
apply_patches "$TARGET_DIR" libBigWig
rm -f "$TARGET_DIR/bwRead.c.orig"

# Keep imported text compatible with the repository-wide whitespace check.
# This changes formatting only and is replayed on every vendoring run.
while IFS= read -r -d '' file; do
  perl -0pi -e 's/[ \t]+$//mg; s/(?:\r?\n)+\z/\n/' "$file"
done < <(find "$TARGET_DIR" -maxdepth 1 -type f -print0)

echo "$LIBBIGWIG_COMMIT" > "$TARGET_DIR/COMMIT"
echo "https://github.com/dpryan79/libBigWig/tree/$LIBBIGWIG_COMMIT" \
  > "$TARGET_DIR/SOURCE_URL"

mkdir -p "$TP_DIR/licenses/libBigWig"
cp "$TARGET_DIR/LICENSE" "$TP_DIR/licenses/libBigWig/LICENSE"

echo "Vendored libBigWig ${LIBBIGWIG_COMMIT} into $TARGET_DIR"
