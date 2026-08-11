#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/duckhts-cgranges-r.XXXXXX")"
trap 'rm -rf "$TMP_DIR"' EXIT

Rscript "$ROOT_DIR/r/duckhtsbench/scripts/cgranges_benchmark.R" \
  --extension "$ROOT_DIR/build/release/duckhts.duckdb_extension" \
  --bedtk duckhts-test-no-bedtk --bedtools duckhts-test-no-bedtools \
  --subjects 20 --queries 10 --passes 1 --out-dir "$TMP_DIR/out" >/dev/null

Rscript - "$TMP_DIR/out" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
summary <- utils::read.csv(file.path(args[[1L]], "summary.csv"), check.names = FALSE)
metadata <- jsonlite::fromJSON(file.path(args[[1L]], "metadata.json"))
stopifnot(identical(summary$variant, c("scalar_filter", "scalar_count", "scalar_expand")))
stopifnot(identical(as.integer(summary$subject_intervals), rep.int(20L, 3L)))
stopifnot(identical(as.integer(summary$query_intervals), rep.int(10L, 3L)))
stopifnot(nzchar(metadata$source_revision), nzchar(metadata$extension_md5))
RS

echo "R cgranges benchmark: OK"
