#!/usr/bin/env bash
# Loads the built extension into a process that already carries a host library
# with its own hopen() and seq_nt16_str, and checks that the bundled htslib
# uses neither.
set -eu
cd "$(dirname "$0")/../.."
extension=$(realpath "${1:-build/release/duckhts.duckdb_extension}")
python=${DUCKHTS_TEST_PYTHON:-python3}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
"${CC:-cc}" -std=gnu11 -O1 -Wall -Wextra -Werror -fPIC -shared \
    -Ithird_party/htslib \
    test/scripts/htslib_binding_host.c \
    -o "$work/htslib_binding_host.so"
LD_PRELOAD="$work/htslib_binding_host.so" timeout 60 "$python" \
    test/scripts/htslib_binding_test.py \
    --extension "$extension" --host "$work/htslib_binding_host.so"
