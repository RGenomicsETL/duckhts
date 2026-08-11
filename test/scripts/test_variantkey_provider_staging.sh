#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/duckhts-variantkey-stage-test.XXXXXX")"
SOURCE_DIR="$TMP_DIR/source"
FAKE_BIN="$TMP_DIR/bin"
CACHE_DIR="$TMP_DIR/cache"
REGISTRY="$TMP_DIR/registry.tsv"
trap 'rm -rf "$TMP_DIR"' EXIT
mkdir -p "$SOURCE_DIR" "$FAKE_BIN"

cat >"$FAKE_BIN/samtools" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
[[ "$1" == "faidx" ]]
printf '1\t4\t3\t4\t5\n' >"${2}.fai"
EOF
chmod +x "$FAKE_BIN/samtools"

Rscript - "$ROOT_DIR" "$SOURCE_DIR" "$CACHE_DIR" "$REGISTRY" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
root <- args[[1L]]
source_dir <- args[[2L]]
cache_dir <- args[[3L]]
registry_path <- args[[4L]]
registry <- utils::read.delim(
  file.path(root, "r", "duckhtsbench", "inst", "benchmark_registry.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

writeLines("fixture", file.path(source_dir, "plain"))
connection <- gzfile(file.path(source_dir, "fasta.gz"), "wb")
writeLines(c(">1", "ACGT"), connection)
close(connection)
work <- file.path(source_dir, "revel")
dir.create(work)
writeLines(c("chr,hg19_pos,ref,alt,REVEL", "1,1,A,C,0.5"), file.path(work, "revel_with_transcript_ids"))
old <- getwd()
setwd(work)
utils::zip(file.path(source_dir, "revel.zip"), "revel_with_transcript_ids")
setwd(old)
work <- file.path(source_dir, "clinvarbitration")
dir.create(work)
writeLines(c("contig\tposition\treference\talternate", "1\t1\tA\tC"), file.path(work, "clinvar_decisions.tsv"))
old <- getwd()
setwd(work)
utils::tar(file.path(source_dir, "clinvarbitration.tar.gz"), files = "clinvar_decisions.tsv", compression = "gzip", tar = "internal")
setwd(old)

source_map <- stats::setNames(rep(file.path(source_dir, "plain"), nrow(registry)), registry$id)
source_map[["ensembl116_grch38_fasta"]] <- file.path(source_dir, "fasta.gz")
source_map[["revel_v13_source_zip"]] <- file.path(source_dir, "revel.zip")
source_map[["clinvarbitration_202508_source_archive"]] <- file.path(source_dir, "clinvarbitration.tar.gz")
for (id in registry$id[registry$transform == "direct_download"]) {
  registry$locator[registry$id == id] <- paste0("file://", source_map[[id]])
  registry$supplier_identity[registry$id == id] <- ""
}
utils::write.table(registry, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)

model <- file.path(cache_dir, registry$cache_relpath[registry$id == "duckvep_ensembl116_model"])
dir.create(dirname(model), recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(duckdb))
con <- DBI::dbConnect(duckdb::duckdb(), model)
DBI::dbExecute(con, "CREATE TABLE model_regions (seq_region UINTEGER, seq_region_name VARCHAR)")
DBI::dbExecute(con, "INSERT INTO model_regions VALUES (1, '1')")
DBI::dbExecute(con, paste(
  "CREATE TABLE duckvep_regulation_features (seq_region UINTEGER, stable_id VARCHAR,",
  "feature_start UINTEGER, feature_end UINTEGER, feature_class VARCHAR,",
  "feature_so_term VARCHAR, feature_so_accession VARCHAR)"
))
DBI::dbExecute(con, "INSERT INTO duckvep_regulation_features VALUES (1, 'ENSR1', 10, 20, 'regulatory_region', 'enhancer', 'SO:0000165')")
DBI::dbExecute(con, "CREATE TABLE model_transcripts (gene_stable_id VARCHAR)")
DBI::dbExecute(con, "INSERT INTO model_transcripts VALUES ('ENSG000001')")
DBI::dbExecute(con, paste(
  "CREATE TABLE model_receipt (source_name VARCHAR, source_version VARCHAR, assembly VARCHAR,",
  "source_manifest_sha256 VARCHAR, reference_sha256 VARCHAR, model_sha256 VARCHAR)"
))
DBI::dbExecute(con, paste(
  "INSERT INTO model_receipt VALUES ('Ensembl', '116', 'GRCh38',",
  "'34eaee64f47916a2d9d3d864377f911424741c7ead12883e47dcf998f9677703',",
  "'d8c3af0094a7bba6125763bad779ec18a81483c739c6ed122094bdf86c187b92',",
  "'296bc90633562da49e900a723675c7a0b4aba53b8cbdec027f107d46972aec68')"
))
DBI::dbDisconnect(con, shutdown = TRUE)
registry <- utils::read.delim(registry_path, stringsAsFactors = FALSE, check.names = FALSE)
model_row <- registry$id == "duckvep_ensembl116_model"
identity <- registry$supplier_identity[model_row]
identity <- sub("sha256=[^;]+", paste0("sha256=", unname(tools::sha256sum(model))), identity)
identity <- sub("bytes=[0-9]+", paste0("bytes=", file.info(model)$size), identity)
registry$supplier_identity[model_row] <- identity
utils::write.table(registry, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)
RS

PATH="$FAKE_BIN:$PATH" DUCKHTS_CACHE_DIR="$CACHE_DIR" DUCKHTSBENCH_REGISTRY="$REGISTRY" \
  Rscript "$ROOT_DIR/r/duckhtsbench/scripts/stage_variantkey_providers.R" >/dev/null

printf 'stale parquet' >"$CACHE_DIR/benchmarks/variantkey-providers/raw/revel_grch37.parquet"
PATH="$FAKE_BIN:$PATH" DUCKHTS_CACHE_DIR="$CACHE_DIR" DUCKHTSBENCH_REGISTRY="$REGISTRY" \
  Rscript "$ROOT_DIR/r/duckhtsbench/scripts/stage_variantkey_providers.R" >/dev/null

Rscript - "$CACHE_DIR/models/duckvep/ensembl-116-grch38.duckdb" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(duckdb))
con <- DBI::dbConnect(duckdb::duckdb(), args[[1L]])
DBI::dbExecute(con, "INSERT INTO model_transcripts VALUES ('ENSG_MUTATED')")
DBI::dbDisconnect(con, shutdown = TRUE)
RS
if PATH="$FAKE_BIN:$PATH" DUCKHTS_CACHE_DIR="$CACHE_DIR" DUCKHTSBENCH_REGISTRY="$REGISTRY" \
  Rscript "$ROOT_DIR/r/duckhtsbench/scripts/stage_variantkey_providers.R" >/dev/null 2>&1; then
  echo "mutated DuckVEP model unexpectedly passed identity validation" >&2
  exit 1
fi

Rscript - "$REGISTRY" "$CACHE_DIR" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
registry <- utils::read.delim(args[[1]], stringsAsFactors = FALSE, check.names = FALSE)
paths <- file.path(args[[2]], registry$cache_relpath)
expected <- registry$id[registry$workload == "variantkey-providers"]
for (id in expected) {
  output <- paths[registry$id == id]
  stopifnot(length(output) == 1L, file.exists(output), file.info(output)$size > 0L)
  stopifnot(file.exists(paste0(output, ".provenance.tsv")))
}
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(duckdb))
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
for (id in c("duckvep_ensembl116_regulatory_parquet", "duckvep_ensembl116_transcripts_parquet")) {
  output <- gsub("'", "''", paths[registry$id == id], fixed = TRUE)
  stopifnot(DBI::dbGetQuery(con, sprintf("SELECT count(*) AS n FROM read_parquet('%s')", output))$n == 1L)
}
DBI::dbDisconnect(con, shutdown = TRUE)
RS

echo "VariantKey provider staging: OK"
