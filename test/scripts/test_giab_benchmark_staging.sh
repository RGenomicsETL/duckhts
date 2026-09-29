#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/duckhts-giab-benchmark-stage.XXXXXX")"
trap 'rm -rf "$TMP_DIR"' EXIT
mkdir -p "$TMP_DIR/cache"

Rscript - "$ROOT_DIR" "$TMP_DIR" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
root <- args[[1L]]
work <- args[[2L]]
registry <- utils::read.delim(
  file.path(root, "r", "duckhtsbench", "inst", "benchmark_registry.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
row <- registry[registry$workload == "giab-benchmark-vcf", , drop = FALSE]
stopifnot(nrow(row) == 1L, identical(row$transform, "direct_download"))
source <- file.path(work, "fixture.vcf.gz")
writeLines("##fileformat=VCFv4.2", source)
row$locator <- paste0("file://", source)
row$supplier_identity <- paste0(
  "bytes=", file.info(source)$size, ";md5=", unname(tools::md5sum(source))
)
utils::write.table(row, file.path(work, "registry.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE)
RS

export DUCKHTS_CACHE_DIR="$TMP_DIR/cache"
export DUCKHTSBENCH_REGISTRY="$TMP_DIR/registry.tsv"
Rscript "$ROOT_DIR/r/duckhtsbench/scripts/stage_giab_benchmark_vcf.R" >/dev/null
Rscript - "$TMP_DIR" <<'RS'
work <- commandArgs(trailingOnly = TRUE)[[1L]]
registry <- utils::read.delim(file.path(work, "registry.tsv"),
                              stringsAsFactors = FALSE, check.names = FALSE)
output <- file.path(work, "cache", registry$cache_relpath)
stopifnot(identical(readLines(output), "##fileformat=VCFv4.2"),
          file.exists(paste0(output, ".provenance.tsv")))
writeLines("stale cache", output)
RS
Rscript "$ROOT_DIR/r/duckhtsbench/scripts/stage_giab_benchmark_vcf.R" >/dev/null
Rscript - "$TMP_DIR" <<'RS'
work <- commandArgs(trailingOnly = TRUE)[[1L]]
registry <- utils::read.delim(file.path(work, "registry.tsv"),
                              stringsAsFactors = FALSE, check.names = FALSE)
output <- file.path(work, "cache", registry$cache_relpath)
stopifnot(identical(readLines(output), "##fileformat=VCFv4.2"),
          file.exists(paste0(output, ".provenance.tsv")))
RS

echo "GIAB benchmark staging: OK"
