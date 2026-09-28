#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/duckhts-variantkey-stage-test.XXXXXX")"
trap 'rm -rf "$TMP_DIR"' EXIT
mkdir -p "$TMP_DIR/source" "$TMP_DIR/cache"

Rscript - "$ROOT_DIR" "$TMP_DIR" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
root <- args[[1L]]
work <- args[[2L]]
source_dir <- file.path(work, "source")
registry <- utils::read.delim(
  file.path(root, "r", "duckhtsbench", "inst", "benchmark_registry.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
writeLines("fixture", file.path(source_dir, "plain"))
connection <- gzfile(file.path(source_dir, "fasta.gz"), "wb")
writeLines(c(">chr1", "ACGT"), connection)
close(connection)
writeLines(c("chr,hg19_pos,ref,alt,REVEL", "1,1,A,C,0.5", "X,2,A,C,0.3"),
           file.path(source_dir, "revel_with_transcript_ids"))
old <- getwd()
setwd(source_dir)
utils::zip("revel.zip", "revel_with_transcript_ids")
writeLines(c("contig\tposition\treference\talternate", "1\t1\tA\tC"),
           "clinvar_decisions.tsv")
utils::tar("clinvarbitration.tar.gz", files = "clinvar_decisions.tsv",
           compression = "gzip", tar = "internal")
setwd(old)
source_map <- c(
  ensembl116_grch38_fasta = "fasta.gz",
  revel_v13_source_zip = "revel.zip",
  clinvarbitration_202508_source_archive = "clinvarbitration.tar.gz"
)
for (id in registry$id[registry$workload == "variantkey-providers" &
                       registry$transform == "direct_download"]) {
  source <- source_map[id]
  if (is.na(source)) source <- "plain"
  registry$locator[registry$id == id] <- paste0("file://", file.path(source_dir, source))
  registry$supplier_identity[registry$id == id] <- ""
}
utils::write.table(registry, file.path(work, "registry.tsv"), sep = "\t",
                   row.names = FALSE, quote = FALSE)
RS

export DUCKHTS_CACHE_DIR="$TMP_DIR/cache"
export DUCKHTSBENCH_REGISTRY="$TMP_DIR/registry.tsv"
Rscript "$ROOT_DIR/r/duckhtsbench/scripts/stage_variantkey_providers.R" >/dev/null
printf 'stale parquet' > "$TMP_DIR/cache/benchmarks/variantkey-providers/raw/revel_grch37.parquet"
Rscript "$ROOT_DIR/r/duckhtsbench/scripts/stage_variantkey_providers.R" >/dev/null

Rscript - "$TMP_DIR/registry.tsv" "$TMP_DIR/cache" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
registry <- utils::read.delim(args[[1L]], stringsAsFactors = FALSE, check.names = FALSE)
rows <- registry[registry$workload == "variantkey-providers", , drop = FALSE]
paths <- file.path(args[[2L]], rows$cache_relpath)
stopifnot(length(paths) > 0L, all(file.exists(paths)),
          all(file.info(paths)$size > 0L),
          all(file.exists(paste0(paths, ".provenance.tsv"))))
revel <- paths[rows$id == "revel_v13_grch37"]
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(duckdb))
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
result <- DBI::dbGetQuery(con, sprintf(
  "SELECT chrom, pos FROM read_parquet('%s') ORDER BY pos",
  gsub("'", "''", revel, fixed = TRUE)
))
stopifnot(identical(result$chrom, c("1", "X")), identical(result$pos, c(1, 2)))
RS

echo "VariantKey provider staging: OK"
