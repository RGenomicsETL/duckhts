#!/usr/bin/env bash
# Synchronized Linux regression for output created after the initial path check.
set -euo pipefail
if [[ $(uname -s) != Linux ]]; then
    echo 'writer race test: Linux only'
    exit 0
fi
extension=$(realpath "${1:-build/release/duckhts.duckdb_extension}")
cli=${DUCKDB_CLI:-duckdb}
root=$(mktemp -d)
reader=
cleanup() {
    if [[ -n $reader ]]; then
        kill "$reader" 2>/dev/null || true
        wait "$reader" 2>/dev/null || true
    fi
    rm -rf "$root"
}
trap cleanup EXIT

run_case() {
    local name=$1
    local dir="$root/$name"
    local sql status=0 blocked=false
    mkdir "$dir"
    mkfifo "$dir/input"
    if [[ $name == bgzip ]]; then
        sql="SELECT * FROM bgzip('$dir/input', output_path := '$dir/output', threads := 1, overwrite := false);"
    else
        sql="SELECT * FROM duckhts_samtools_idxstats('$dir/input', output := '$dir/output', overwrite := false);"
    fi
    "$cli" -unsigned -c "LOAD '$extension'; $sql" >"$dir/query.log" 2>&1 &
    reader=$!
    for ((i = 0; i < 400; i++)); do
        if grep -q wait_for_partner "/proc/$reader/wchan" 2>/dev/null; then
            blocked=true
            break
        fi
        if ! kill -0 "$reader" 2>/dev/null; then
            break
        fi
        sleep 0.05
    done
    if [[ $blocked != true ]]; then
        echo "$name: input-open synchronization failed" >&2
        return 1
    fi
    printf 'independently published data\n' >"$dir/expected"
    cp "$dir/expected" "$dir/output"
    if [[ $name == bgzip ]]; then
        printf 'payload\n' >"$dir/input"
    else
        printf '@HD\tVN:1.6\n@SQ\tSN:chr1\n' >"$dir/input"
    fi
    wait "$reader" || status=$?
    reader=
    if [[ $status -eq 0 ]] || ! cmp -s "$dir/expected" "$dir/output"; then
        echo "$name: query must fail and preserve the concurrent output" >&2
        return 1
    fi
    if [[ $name == idxstats ]] && ! grep -q 'failed to read alignment header' "$dir/query.log"; then
        echo 'idxstats: expected the malformed SAM header error' >&2
        return 1
    fi
    echo "$name: failed without clobbering concurrent output"
}
run_case bgzip
run_case idxstats
