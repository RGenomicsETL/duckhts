#!/bin/sh
# Run native fault-injection tests against the production translation units.
set -eu
cd "$(dirname "$0")/../.."
mkdir -p build
tmp=$(mktemp -d build/error-tests.XXXXXX)
trap 'rm -rf "$tmp"' EXIT HUP INT TERM
if [ "$#" -eq 0 ]; then
    set -- registration_errors_test bigwig_errors_test norm_errors_test
fi
for name in "$@"; do
    "${CC:-cc}" -std=gnu11 -Wall -Wextra -Werror \
        -ffunction-sections -fdata-sections \
        -Iduckdb_capi -Isrc/include -Ithird_party/htslib -Ithird_party/libBigWig \
        "test/scripts/$name.c" -Wl,--gc-sections -pthread -o "$tmp/$name"
    timeout 30 "$tmp/$name"
    printf '%s: PASS\n' "$name"
done
