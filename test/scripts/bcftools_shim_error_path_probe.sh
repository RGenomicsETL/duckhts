#!/usr/bin/env bash
# Compile and run a focused probe for duckhts_bcftools_error* recovery behavior.
#
# This intentionally exercises two paths:
#   1. guarded calls inside duckhts_filter_try_begin()/duckhts_filter_try_end(),
#      which must longjmp back to the recovery scope;
#   2. unguarded calls, which currently log to stderr and return to the caller.
#
# The unguarded-return result documents the invalid-state risk: caller code after
# duckhts_bcftools_error* continues executing unless it explicitly checks an
# errors-as-values status of its own.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CC="${CC:-cc}"
BUILD_DIR="${TMPDIR:-/tmp}/duckhts_bcftools_shim_probe.$$"

cleanup() {
  rm -rf "$BUILD_DIR"
}
trap cleanup EXIT

mkdir -p "$BUILD_DIR"

"$CC" \
  -std=c99 \
  -Wall \
  -Wextra \
  -I"$REPO_ROOT/src/include" \
  -I"$REPO_ROOT/third_party/htslib" \
  "$REPO_ROOT/test/scripts/bcftools_shim_error_path_probe.c" \
  "$REPO_ROOT/src/bcftools_shim.c" \
  -o "$BUILD_DIR/bcftools_shim_error_path_probe"

"$BUILD_DIR/bcftools_shim_error_path_probe" \
  >"$BUILD_DIR/stdout.txt" \
  2>"$BUILD_DIR/stderr.txt"

cat "$BUILD_DIR/stdout.txt"
cat "$BUILD_DIR/stderr.txt" >&2

grep -q '^GUARDED_LONGJMP=1$' "$BUILD_DIR/stdout.txt"
grep -q '^GUARDED_ERRNO_LONGJMP=1$' "$BUILD_DIR/stdout.txt"
grep -q '^UNGUARDED_RETURNED=1$' "$BUILD_DIR/stdout.txt"
grep -q '^UNGUARDED_CALLER_STATE=2$' "$BUILD_DIR/stdout.txt"
grep -q '^UNGUARDED_ERRNO_RETURNED=1$' "$BUILD_DIR/stdout.txt"
grep -q '^UNGUARDED_ERRNO_CALLER_STATE=2$' "$BUILD_DIR/stdout.txt"
grep -q '^BCFTOOLS_SHIM_ERROR_PATH_PROBE_OK=1$' "$BUILD_DIR/stdout.txt"

grep -q '^unguarded boom$' "$BUILD_DIR/stderr.txt"
grep -q '^unguarded errno: No such file' "$BUILD_DIR/stderr.txt"
