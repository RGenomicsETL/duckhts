#!/usr/bin/env bash

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
work_dir=$(mktemp -d)
trap 'rm -rf "$work_dir"' EXIT

cd "$repo_dir"

"${CC:-cc}" -std=c11 -O2 -g -Wall -Wextra -Werror -pedantic \
    -Isrc/include src/somalier.c test/scripts/somalier_native_test.c \
    -lm -o "$work_dir/somalier_native_test"

"$work_dir/somalier_native_test" --campaign pinned > "$work_dir/pinned.tsv"
Rscript test/scripts/somalier_v034_differential.R "$work_dir/pinned.tsv"
Rscript test/scripts/somalier_statistical_campaign.R "$work_dir/somalier_native_test"
Rscript test/scripts/test_somalier_upstream_staging.R
