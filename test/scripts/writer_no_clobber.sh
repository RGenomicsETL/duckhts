#!/usr/bin/env bash
# Linux regressions for concurrent creators and symlink output ownership.
set -euo pipefail
if [[ $(uname -s) != Linux ]]; then
    echo 'writer race and symlink tests: Linux only (skipped)'
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

run_symlink_case() {
    local name=$1
    local outcome=$2
    local dir="$root/${name}_${outcome}"
    local input sql status=0
    mkdir "$dir"
    printf 'referent must survive\n' >"$dir/referent"
    ln -s referent "$dir/output"
    case "$name" in
        bgzip)
            if [[ $outcome == success ]]; then
                printf 'payload\n' >"$dir/input"
            else
                mkdir "$dir/input"
            fi
            input="$dir/input"
            sql="SELECT * FROM bgzip('$input', output_path := '$dir/output', threads := 1, overwrite := true);"
            ;;
        bgunzip)
            printf 'payload\n' >"$dir/plain"
            "$cli" -unsigned -c "LOAD '$extension'; SELECT * FROM bgzip('$dir/plain', output_path := '$dir/input.gz', threads := 1);" >"$dir/prepare.log" 2>&1
            if [[ $outcome == failure ]]; then
                printf '\377' | dd of="$dir/input.gz" bs=1 seek=20 conv=notrunc status=none
            fi
            input="$dir/input.gz"
            sql="SELECT * FROM bgunzip('$input', output_path := '$dir/output', threads := 1, overwrite := true);"
            ;;
        idxstats)
            if [[ $outcome == success ]]; then
                input="$(pwd)/test/data/range.bam"
            else
                printf '@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\nbad\t0\tchr1\tinvalid\t60\t1M\t*\t0\t0\tA\tI\n' >"$dir/input.sam"
                input="$dir/input.sam"
            fi
            sql="SELECT * FROM duckhts_samtools_idxstats('$input', output := '$dir/output', overwrite := true);"
            ;;
    esac
    "$cli" -unsigned -c "LOAD '$extension'; $sql" >"$dir/query.log" 2>&1 || status=$?
    if [[ $outcome == success ]]; then
        if [[ $status -ne 0 || ! -f $dir/output || -L $dir/output ]]; then
            echo "$name: overwrite must replace symlink with a regular output" >&2
            return 1
        fi
        if [[ $name == bgzip ]]; then
            gzip -dc "$dir/output" >"$dir/result"
            cmp "$dir/input" "$dir/result"
        elif [[ $name == bgunzip ]]; then
            cmp "$dir/plain" "$dir/output"
        else
            grep -q $'CHROMOSOME_I\t' "$dir/output"
        fi
    else
        if [[ $status -eq 0 || -e $dir/output || -L $dir/output ]]; then
            echo "$name: failed overwrite must remove only its own output" >&2
            return 1
        fi
        if [[ $name == idxstats ]]; then
            grep -q 'failed while scanning input' "$dir/query.log"
        else
            grep -q 'read error' "$dir/query.log"
        fi
    fi
    if [[ -L $dir/output ]] || ! cmp -s "$dir/referent" <(printf 'referent must survive\n'); then
        echo "$name: overwrite touched the symlink referent" >&2
        return 1
    fi
    echo "$name: symlink overwrite $outcome preserved the referent"
}
for writer in bgzip bgunzip idxstats; do
    run_symlink_case "$writer" success
    run_symlink_case "$writer" failure
done
