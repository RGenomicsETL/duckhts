#!/usr/bin/env bash
# Run from the repository root after make release.
set -euo pipefail

work=$(mktemp -d "$PWD/test/bgzip-identity.XXXXXX")
trap 'rm -rf "$work"' EXIT
extension="$PWD/build/release/duckhts.duckdb_extension"
printf 'hello\n' > "$work/input"
cp "$work/input" "$work/expected"
ln "$work/input" "$work/hardlink"
ln -s input "$work/symlink"

for alias in hardlink symlink; do
    if timeout 30 duckdb -unsigned -c "LOAD '$extension'; SELECT * FROM bgzip('$work/input', output_path := '$work/$alias', overwrite := TRUE, keep := FALSE);" > "$work/result" 2>&1; then
        printf 'bgzip accepted %s as output\n' "$alias" >&2
        exit 1
    fi
    grep -q 'same file' "$work/result"
    cmp "$work/input" "$work/expected"
    cmp "$work/$alias" "$work/expected"
done

timeout 30 duckdb -unsigned -c "LOAD '$extension'; SELECT * FROM bgzip('$work/input', output_path := '$work/input.gz');" > /dev/null
ln "$work/input.gz" "$work/compressed-link"
cp "$work/input.gz" "$work/expected.gz"
if timeout 30 duckdb -unsigned -c "LOAD '$extension'; SELECT * FROM bgunzip('$work/input.gz', output_path := '$work/compressed-link', overwrite := TRUE, keep := FALSE);" > "$work/result" 2>&1; then
    printf 'bgunzip accepted a hard link as output\n' >&2
    exit 1
fi
grep -q 'same file' "$work/result"
cmp "$work/input.gz" "$work/expected.gz"
cmp "$work/compressed-link" "$work/expected.gz"
