#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=vendor_common.sh
source "$SCRIPT_DIR/vendor_common.sh"

CGRANGES_COMMIT="b3d5e2c5b9a0a379f54592ab85f6cff5d58c387e"
CGRANGES_BASE_URL="https://raw.githubusercontent.com/lh3/cgranges/${CGRANGES_COMMIT}"

LOCAL_SOURCE_DEFAULT="$ROOT_DIR/../RCGRanges/src"
LOCAL_SOURCE="${CGRANGES_LOCAL_SOURCE:-$LOCAL_SOURCE_DEFAULT}"

CACHE_DIR="$DIST_DIR/cgranges-${CGRANGES_COMMIT}"
TARGET_DIR="$TP_DIR/cgranges"

expected_sha256() {
  case "$1" in
    cgranges.h) echo "dfb6223abc07e98bd054ea31dea297e5e493ba95d77bca0804dad739f5e19997" ;;
    cgranges.c) echo "bb969a661f31f9d9c5bdd68be3d5d9bae265c679d48f1aea031a0c0e0b02fc99" ;;
    khash.h)    echo "4582c7fd4b98f61339ce677bb1daa7a45a53075b112a494a335525c4bda72074" ;;
    README.md)  echo "33c33093b86253e11233ec17d7d90e20fa2f9301b9a8d369140f1441de67dcbd" ;;
    *)
      echo "Unknown CGRanges file: $1" >&2
      exit 1
      ;;
  esac
}

stage_file() {
  local file="$1"
  local sha
  local cached
  sha="$(expected_sha256 "$file")"
  cached="$CACHE_DIR/$file"

  mkdir -p "$CACHE_DIR"
  if [[ -f "$cached" ]]; then
    echo "${sha}  ${cached}" | sha256sum -c -
    return
  fi

  if [[ -f "$LOCAL_SOURCE/$file" ]]; then
    cp "$LOCAL_SOURCE/$file" "$cached"
    echo "${sha}  ${cached}" | sha256sum -c -
    echo "Cached $file from local source: $LOCAL_SOURCE"
    return
  fi

  download_if_missing "$CGRANGES_BASE_URL/$file" "$cached" "$sha"
}

FILES=("cgranges.h" "cgranges.c" "khash.h" "README.md")

for f in "${FILES[@]}"; do
  stage_file "$f"
done

reset_dir "$TARGET_DIR"
for f in "${FILES[@]}"; do
  cp "$CACHE_DIR/$f" "$TARGET_DIR/$f"
done

echo "$CGRANGES_COMMIT" > "$TARGET_DIR/COMMIT"
echo "https://github.com/lh3/cgranges/tree/$CGRANGES_COMMIT" > "$TARGET_DIR/SOURCE_URL"

mkdir -p "$TP_DIR/licenses/cgranges"
cp "$TARGET_DIR/cgranges.h" "$TP_DIR/licenses/cgranges/MIT_LICENSE_HEADER.txt"
cp "$TARGET_DIR/README.md" "$TP_DIR/licenses/cgranges/README.md"

echo "Vendored cgranges ${CGRANGES_COMMIT} into $TARGET_DIR"
