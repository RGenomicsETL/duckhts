#!/usr/bin/env bash
# Retain a complete first BGZF block and truncate the following block.
set -euo pipefail
cd "$(dirname "$0")/../.."
bgzip=${BGZIP:-bgzip}
work=$(mktemp -d test/score-truncated-fixtures-XXXXXX)
trap 'rm -rf "$work"' EXIT
for format in tsv vcf header.tsv; do
    case "$format" in
        tsv) source=test/data/score_summary.tsv; prefix_lines=2 ;;
        vcf) source=test/data/score_gwas_summary.vcf; prefix_lines=6 ;;
        header.tsv) source=test/data/score_summary_comments.tsv; prefix_lines=1 ;;
    esac
    head -n "$prefix_lines" "$source" | "$bgzip" -c > "$work/prefix.gz"
    tail -n +"$((prefix_lines + 1))" "$source" | "$bgzip" -c > "$work/suffix.gz"
    # Drop the 28-byte EOF marker from the prefix. Drop the suffix EOF marker
    # and 12 bytes of its data block, leaving a physical read failure.
    head -c -28 "$work/prefix.gz" > "test/data/score_truncated.$format.gz"
    head -c -40 "$work/suffix.gz" >> "test/data/score_truncated.$format.gz"
done
