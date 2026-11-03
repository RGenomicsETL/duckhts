#!/usr/bin/env bash
set -eu
cd "$(dirname "$0")/../.."
work=$(mktemp -d "${PWD}/test/scripts/vector-size.XXXXXX")
trap 'rm -rf "$work"' EXIT
for i in $(seq 1 64); do
    printf 'chr1\t0\t1\n'
done > "$work/input.bed"
failed=0
for mode in tabix interval; do
    flags=()
    if [ "$mode" = tabix ]; then
        flags=(-DTEST_TABIX)
    fi
    "${CC:-cc}" -std=gnu11 -O1 -ffunction-sections -fdata-sections \
        -Iduckdb_capi -Isrc/include -Ithird_party/htslib -Ithird_party/cgranges \
        "${flags[@]}" test/scripts/reader_vector_size_test.c \
        -Wl,--gc-sections -Wl,-rpath,"${PWD}/build/release" \
        -Lbuild/release -l:duckhts.duckdb_extension -lm -o "$work/$mode"
    "$work/$mode" "$work/input.bed" || failed=1
done
exit "$failed"
