#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=vendor_common.sh
source "$SCRIPT_DIR/vendor_common.sh"

VARIANTKEY_COMMIT="4ee491c11cd1990ca628eab207054da028f8910a"
VARIANTKEY_BASE_URL="https://raw.githubusercontent.com/genomics-dev/variantkey/${VARIANTKEY_COMMIT}"

LOCAL_SOURCE_DEFAULT="$ROOT_DIR/.sync/variantkey"
LOCAL_SOURCE="${VARIANTKEY_LOCAL_SOURCE:-$LOCAL_SOURCE_DEFAULT}"

CACHE_DIR="$DIST_DIR/variantkey-${VARIANTKEY_COMMIT}"
TARGET_DIR="$TP_DIR/variantkey"

expected_sha256() {
  case "$1" in
    c/src/variantkey/hex.h) echo "a7f71ccb312865c670ae9f1d269b49ee6eb0074fa704f5d08d0b947d8b7a7bc9" ;;
    c/src/variantkey/variantkey.h) echo "6bef9b79e5ef63e467092968b700e9004a052720bff67ca3619330ba80bb06da" ;;
    c/src/variantkey/regionkey.h) echo "32f9432615b57abc2c77daca835ae810e17982eec69c7da7311997f5c7aae71a" ;;
    LICENSE) echo "ae6d86d0d22adbe6ab061c6ec76748eb78bb5b4250c4fdc8643002825fd5ecab" ;;
    README.md) echo "86d60ce8286ab002e5bffa939eb5cd80f1ce5ceefd4c546b5744d1e9d0d6067b" ;;
    *)
      echo "Unknown VariantKey file: $1" >&2
      exit 1
      ;;
  esac
}

stage_file() {
  local rel="$1"
  local sha
  local cached

  sha="$(expected_sha256 "$rel")"
  cached="$CACHE_DIR/$rel"

  mkdir -p "$(dirname "$cached")"
  if [[ -f "$cached" ]]; then
    echo "${sha}  ${cached}" | sha256sum -c -
    return
  fi

  if [[ -f "$LOCAL_SOURCE/$rel" ]]; then
    cp "$LOCAL_SOURCE/$rel" "$cached"
    echo "${sha}  ${cached}" | sha256sum -c -
    echo "Cached $rel from local source: $LOCAL_SOURCE"
    return
  fi

  download_if_missing "$VARIANTKEY_BASE_URL/$rel" "$cached" "$sha"
}

FILES=(
  "c/src/variantkey/hex.h"
  "c/src/variantkey/variantkey.h"
  "c/src/variantkey/regionkey.h"
  "LICENSE"
  "README.md"
)

for f in "${FILES[@]}"; do
  stage_file "$f"
done

reset_dir "$TARGET_DIR"
mkdir -p "$TARGET_DIR/include/variantkey"
cp "$CACHE_DIR/c/src/variantkey/hex.h" "$TARGET_DIR/include/variantkey/hex.h"
cp "$CACHE_DIR/c/src/variantkey/variantkey.h" "$TARGET_DIR/include/variantkey/variantkey.h"
cp "$CACHE_DIR/c/src/variantkey/regionkey.h" "$TARGET_DIR/include/variantkey/regionkey.h"
cp "$CACHE_DIR/LICENSE" "$TARGET_DIR/LICENSE"
cp "$CACHE_DIR/README.md" "$TARGET_DIR/README.md"

echo "$VARIANTKEY_COMMIT" > "$TARGET_DIR/COMMIT"
echo "https://github.com/genomics-dev/variantkey/tree/$VARIANTKEY_COMMIT" > "$TARGET_DIR/SOURCE_URL"

apply_patches "$TARGET_DIR" "variantkey"

mkdir -p "$TP_DIR/licenses/variantkey"
cp "$TARGET_DIR/LICENSE" "$TP_DIR/licenses/variantkey/LICENSE"
cp "$TARGET_DIR/README.md" "$TP_DIR/licenses/variantkey/README.md"

echo "Vendored VariantKey $VARIANTKEY_COMMIT into $TARGET_DIR"
