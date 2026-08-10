#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CACHE_HELPER="$ROOT_DIR/scripts/duckhts_cache.sh"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/duckhts-cache-test.XXXXXX")"
trap 'rm -rf "$TMP_DIR"' EXIT

HOME="$TMP_DIR/home" XDG_CACHE_HOME="$TMP_DIR/xdg" bash -c '
  set -euo pipefail
  source "$1"
  test "$DUCKHTS_CACHE_DIR" = "$XDG_CACHE_HOME/duckhts"
  test "$(duckhts_cache_subdir benchmarks/norm)" = "$XDG_CACHE_HOME/duckhts/benchmarks/norm"
  ! duckhts_cache_subdir ../outside >/dev/null 2>&1
' bash "$CACHE_HELPER"

XDG_CACHE_HOME="$TMP_DIR/python-xdg" python3 -c '
import sys
from pathlib import Path
sys.path.insert(0, sys.argv[1])
from duckhts_cache import duckhts_cache_dir, duckhts_cache_subdir
assert duckhts_cache_dir() == Path(sys.argv[2]) / "duckhts"
assert duckhts_cache_subdir("benchmarks", "cgranges") == Path(sys.argv[2]) / "duckhts" / "benchmarks" / "cgranges"
' "$ROOT_DIR/scripts" "$TMP_DIR/python-xdg"

DUCKHTS_CACHE_DIR="$TMP_DIR/custom-cache" bash -c '
  set -euo pipefail
  source "$1"
  test "$DUCKHTS_CACHE_DIR" = "$2"
  duckhts_write_provenance "$DUCKHTS_CACHE_DIR/example.provenance.tsv" \
    "workload=cache_test" "source_url=https://example.test/data" "output=example"
  test -s "$DUCKHTS_CACHE_DIR/example.provenance.tsv"
  grep -Fq $'"'"'source_url\thttps://example.test/data'"'"' "$DUCKHTS_CACHE_DIR/example.provenance.tsv"
' bash "$CACHE_HELPER" "$TMP_DIR/custom-cache"

echo "DuckHTS cache helpers: OK"
