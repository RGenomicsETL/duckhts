#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TP_DIR="$ROOT_DIR/third_party"
DIST_DIR="$TP_DIR/distfiles"

require_cmd() {
  command -v "$1" >/dev/null 2>&1 || {
    echo "Missing required command: $1" >&2
    exit 1
  }
}

require_cmd curl
require_cmd sha256sum
require_cmd tar

reset_dir() {
  local dir="$1"
  rm -rf "$dir"
  mkdir -p "$dir"
}

download_if_missing() {
  local url="$1"
  local dest="$2"
  local sha="$3"

  mkdir -p "$(dirname "$dest")"
  if [[ -f "$dest" ]]; then
    echo "Using cached $dest"
  else
    echo "Downloading $url"
    curl -L --fail -o "$dest" "$url"
  fi

  echo "${sha}  ${dest}" | sha256sum -c -
}

extract_tar_bz2() {
  local archive="$1"
  local target="$2"
  tar -xjf "$archive" -C "$target"
}

capture_licenses() {
  local src="$1"
  local name="$2"
  local out="$TP_DIR/licenses/$name"
  mkdir -p "$out"

  for f in LICENSE LICENSE.txt COPYING COPYING.LESSER COPYING.LIB COPYING.LGPL COPYING.GPL NOTICE; do
    if [[ -f "$src/$f" ]]; then
      cp "$src/$f" "$out/"
    fi
  done
}
