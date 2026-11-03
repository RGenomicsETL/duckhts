#!/usr/bin/env bash
# Reconstruct tiny, truncated BGZF fixtures from their checked-in text sources.
set -euo pipefail
cd "$(dirname "$0")/../.."
bgzip=${BGZIP:-/usr/local/lib/R/site-library/RBCFTools/htslib/bin/bgzip}
tabix=${TABIX:-/usr/local/lib/R/site-library/RBCFTools/htslib/bin/tabix}
work=$(mktemp -d test/data/reader-io-XXXXXX)
trap 'rm -rf "$work"' EXIT

# Put the final record in a separate BGZF block. Index the complete BED, then
# remove the EOF block and the last eight bytes of the data block's trailer.
for suffix in bed sam; do
    source="test/data/reader_io_error.$suffix"
    head -n -1 "$source" | "$bgzip" -c > "$work/first.gz"
    tail -n 1 "$source" | "$bgzip" -c > "$work/last.gz"
    head -c -28 "$work/first.gz" > "$work/complete.$suffix.gz"
    cat "$work/last.gz" >> "$work/complete.$suffix.gz"
    if [[ $suffix == bed ]]; then
        "$tabix" -p bed "$work/complete.bed.gz"
        cp "$work/complete.bed.gz.tbi" test/data/reader_io_error.bed.gz.tbi
    fi
    head -c -36 "$work/complete.$suffix.gz" > "test/data/reader_io_error.$suffix.gz"
done

# No complete line is available to schema discovery in this file.
head -c 24 "$work/first.gz" > test/data/reader_io_error_first.gz
