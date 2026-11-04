#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
# Built by CMake with the extension's htslib link libraries (Makefile test target).
work=$(mktemp -d test/score-summary-error-test-XXXXXX)
trap 'rm -rf "$work"' EXIT
cp cmake_build/release/duckhts_score_summary_error_test "$work/test"
if [ "$#" -gt 0 ]; then
    "$work/test" "$@"
    exit
fi
for mode in tsv rsid gwas; do
    case "$mode" in
        tsv) path=test/data/score_summary_skips.tsv ;;
        rsid) path=test/data/score_summary.tsv ;;
        gwas) path=test/data/score_gwas_summary.vcf ;;
    esac
    output=$("$work/test" "$mode" "$path")
    printf '%s\n' "$output"
    count=$(printf '%s\n' "$output" | grep -o 'allocations=[0-9]*' | cut -d= -f2)
    for ((i=1; i<=count; i++)); do
        "$work/test" "$mode" "$path" "$i"
    done
done
for failure in record ES LP; do
    "$work/test" gwas test/data/score_gwas_summary.vcf 0 "$failure"
done
"$work/test" tsv test/data/score_truncated.tsv.gz 0 read
"$work/test" tsv test/data/score_truncated.header.tsv.gz 0 read
"$work/test" gwas test/data/score_truncated.vcf.gz 0 read
