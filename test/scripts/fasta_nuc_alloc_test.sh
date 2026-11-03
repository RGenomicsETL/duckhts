#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
work=$(mktemp -d test/fasta-nuc-alloc-XXXXXX)
trap 'rm -rf "$work"' EXIT
"${CC:-cc}" -std=c11 -D_POSIX_C_SOURCE=200809L -O1 \
    -ffunction-sections -fdata-sections \
    -Iduckdb_capi -Isrc/include -Ithird_party/htslib \
    test/scripts/fasta_nuc_alloc_test.c src/bed_line.c -Wl,--gc-sections \
    third_party/htslib/libhts.a \
    -lz -lbz2 -llzma -lcurl -ldeflate -lcrypto -lpthread -lm \
    -o "$work/test"
"$work/test"
