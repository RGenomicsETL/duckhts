#!/usr/bin/env bash
set -eu
cd "$(dirname "$0")/../.."
work=$(mktemp -d "${PWD}/test/scripts/bam-open.XXXXXX")
trap 'rm -rf "$work"' EXIT
"${CC:-cc}" -std=gnu11 -O1 -ffunction-sections -fdata-sections \
    -Iduckdb_capi -Isrc/include -Ithird_party/htslib \
    test/scripts/bam_open_error_test.c \
    -Wl,--gc-sections -Wl,-rpath,"${PWD}/build/release" \
    -Lbuild/release -l:duckhts.duckdb_extension -lm -o "$work/bam_open_error_test"
timeout 30 "$work/bam_open_error_test" "$work/missing.bam"
