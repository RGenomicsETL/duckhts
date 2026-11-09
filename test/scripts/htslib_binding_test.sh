#!/usr/bin/env bash
# Links a host program that defines its own hopen() against the built extension
# and checks that the bundled htslib does not call it.
set -eu
cd "$(dirname "$0")/../.."
library_dir="${1:-build/release}"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
"${CC:-cc}" -std=gnu11 -O1 -Wall -Wextra -Werror \
    -Ithird_party/htslib \
    test/scripts/htslib_binding_test.c \
    -Wl,-rpath,"$(cd "$library_dir" && pwd)" \
    -L"$library_dir" -l:duckhts.duckdb_extension -o "$work/htslib_binding_test"
timeout 30 "$work/htslib_binding_test" "$work/missing.bam"
