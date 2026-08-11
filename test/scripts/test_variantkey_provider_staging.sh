#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/duckhts-variantkey-stage-test.XXXXXX")"
SOURCE_DIR="$TMP_DIR/source"
FAKE_BIN="$TMP_DIR/bin"
CACHE_DIR="$TMP_DIR/cache root"
REGISTRY="$TMP_DIR/registry.tsv"
trap 'rm -rf "$TMP_DIR"' EXIT
mkdir -p "$SOURCE_DIR" "$FAKE_BIN"

command -v samtools >/dev/null

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
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(duckdb))
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
fixture_dir <- file.path(root, "test", "data", "duckvep", "ensembl_core", "grch38")

reference <- DBI::dbGetQuery(
  con,
  sprintf("SELECT * FROM read_parquet('%s') ORDER BY chrom, start", file.path(fixture_dir, "reference_chunks.parquet"))
)
connection <- gzfile(file.path(source_dir, "fasta.gz"), "wb")
for (chrom in unique(reference$chrom)) {
  writeLines(c(paste0(">", chrom), paste0(reference$seq[reference$chrom == chrom], collapse = "")), connection)
}
close(connection)

work <- file.path(source_dir, "revel")
dir.create(work)
writeLines(c("chr,hg19_pos,ref,alt,REVEL", "1,1,A,C,0.5", "X,2,A,C,0.3"), file.path(work, "revel_with_transcript_ids"))
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
registry$supplier_identity[registry$id == "duckvep_ensembl116_model"] <- ""

source_rows <- registry[grepl("^source_(manifest|schema|table):(core|funcgen)(:|$)", registry$role), , drop = FALSE]
source_rows$group <- vapply(strsplit(source_rows$role, ":", fixed = TRUE), `[[`, character(1L), 2L)
source_rows$kind <- vapply(strsplit(source_rows$role, ":", fixed = TRUE), `[[`, character(1L), 1L)
source_rows$table <- vapply(
  strsplit(source_rows$role, ":", fixed = TRUE),
  function(x) if (length(x) == 3L) x[[3L]] else NA_character_,
  character(1L)
)
source_rows$database <- ifelse(
  source_rows$group == "core", "homo_sapiens_core_116_38", "homo_sapiens_funcgen_116_38"
)

mysql_type <- function(x) {
  if (grepl("INT|BOOL", x)) return("bigint")
  if (grepl("DOUBLE|FLOAT|DECIMAL|NUMERIC", x)) return("double")
  "varchar(255)"
}
write_dump <- function(data, path) {
  for (name in names(data)) {
    if (inherits(data[[name]], "POSIXt")) data[[name]] <- format(data[[name]], "%Y-%m-%d %H:%M:%S")
    if (is.logical(data[[name]])) data[[name]] <- as.integer(data[[name]])
    if (is.character(data[[name]])) {
      data[[name]] <- gsub("\t", "\\\\t", data[[name]], fixed = TRUE)
      data[[name]] <- gsub("\n", "\\\\n", data[[name]], fixed = TRUE)
    }
  }
  connection <- gzfile(path, "wb")
  utils::write.table(data, connection, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, na = "\\N")
  close(connection)
}

for (group in c("core", "funcgen")) {
  selected <- source_rows[source_rows$group == group, , drop = FALSE]
  database <- unique(selected$database)
  group_dir <- file.path(source_dir, group)
  dir.create(group_dir)
  table_rows <- selected[selected$kind == "source_table", , drop = FALSE]
  definitions <- character()
  for (index in seq_len(nrow(table_rows))) {
    table <- table_rows$table[[index]]
    if (table == "meta") {
      data <- if (group == "core") {
        data.frame(
          meta_id = 1:4, species_id = c(NA, 1, 1, 1),
          meta_key = c("schema_version", "assembly.default", "species.production_name", "species.scientific_name"),
          meta_value = c("116", "GRCh38", "homo_sapiens", "Homo sapiens")
        )
      } else {
        data.frame(
          meta_id = 1:2, species_id = c(NA, 1),
          meta_key = c("schema_version", "species.production_name"),
          meta_value = c("116", "homo_sapiens")
        )
      }
      types <- c("BIGINT", "BIGINT", "VARCHAR", "VARCHAR")
    } else {
      parquet_name <- if (group == "funcgen") paste0("funcgen_", table, ".parquet") else paste0(table, ".parquet")
      parquet <- file.path(fixture_dir, parquet_name)
      data <- DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s')", parquet))
      types <- DBI::dbGetQuery(con, sprintf("DESCRIBE SELECT * FROM read_parquet('%s')", parquet))$column_type
    }
    path <- file.path(group_dir, basename(table_rows$cache_relpath[[index]]))
    write_dump(data, path)
    source_map[[table_rows$id[[index]]]] <- path
    columns <- paste0("  `", names(data), "` ", vapply(types, mysql_type, character(1L)))
    definitions <- c(
      definitions,
      paste0("CREATE TABLE `", table, "` (\n", paste(columns, collapse = ",\n"), "\n) ENGINE=MyISAM;")
    )
  }
  schema_row <- selected[selected$kind == "source_schema", , drop = FALSE]
  schema_path <- file.path(group_dir, basename(schema_row$cache_relpath[[1L]]))
  connection <- gzfile(schema_path, "wb")
  writeLines(c(paste0("-- Host: fixture    Database: ", database), definitions), connection)
  close(connection)
  source_map[[schema_row$id[[1L]]]] <- schema_path

  manifest_inputs <- selected[selected$kind != "source_manifest", , drop = FALSE]
  manifest_lines <- character(nrow(manifest_inputs))
  for (index in seq_len(nrow(manifest_inputs))) {
    id <- manifest_inputs$id[[index]]
    fields <- strsplit(trimws(system2("sum", c("-r", shQuote(source_map[[id]])), stdout = TRUE)), "[[:space:]]+")[[1L]]
    manifest_lines[[index]] <- paste(fields[[1L]], fields[[2L]], basename(source_map[[id]]))
    registry$supplier_identity[registry$id == id] <- paste0("sum=", fields[[1L]], ";blocks=", fields[[2L]])
  }
  manifest_row <- selected[selected$kind == "source_manifest", , drop = FALSE]
  manifest_path <- file.path(group_dir, "CHECKSUMS")
  writeLines(manifest_lines, manifest_path)
  source_map[[manifest_row$id[[1L]]]] <- manifest_path
}
DBI::dbDisconnect(con, shutdown = TRUE)

for (id in registry$id[registry$transform == "direct_download"]) {
  registry$locator[registry$id == id] <- paste0("file://", source_map[[id]])
}
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
invisible(DBI::dbExecute(con, "UPDATE model_transcripts SET transcript_flags = transcript_flags + 1 WHERE transcript_index = 0"))
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
model <- paths[registry$id == "duckvep_ensembl116_model"]
model_con <- DBI::dbConnect(duckdb::duckdb(), model, read_only = TRUE)
metadata <- DBI::dbGetQuery(model_con, "SELECT source_release, vep_release, species, species_id, assembly FROM model_metadata")
stopifnot(
  nrow(metadata) == 1L,
  identical(unlist(metadata[1, ], use.names = FALSE), c("116", "116", "homo_sapiens", "1", "GRCh38"))
)
required_tables <- DBI::dbGetQuery(
  model_con,
  paste(
    "SELECT table_schema, table_name FROM information_schema.tables",
    "WHERE table_schema IN ('ensembl_core', 'ensembl_funcgen')"
  )
)
stopifnot(nrow(required_tables) == 16L)
receipt <- DBI::dbGetQuery(model_con, "SELECT model_sha256, reference_sha256 FROM model_receipt")
stopifnot(
  grepl("^[0-9a-f]{64}$", receipt$model_sha256),
  receipt$reference_sha256 == "4fb2dcf3cc1844c54e8105695826d099be075c8cd5b18f74ccd80fb0d38d2cdf"
)
DBI::dbDisconnect(model_con, shutdown = TRUE)

con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
expected_rows <- c(
  duckvep_ensembl116_regulatory_parquet = 3L,
  duckvep_ensembl116_transcripts_parquet = 39L
)
for (id in names(expected_rows)) {
  output <- gsub("'", "''", paths[registry$id == id], fixed = TRUE)
  stopifnot(DBI::dbGetQuery(con, sprintf("SELECT count(*) AS n FROM read_parquet('%s')", output))$n == expected_rows[[id]])
}
revel <- gsub("'", "''", paths[registry$id == "revel_v13_grch37"], fixed = TRUE)
revel_rows <- DBI::dbGetQuery(con, sprintf("SELECT chrom, pos FROM read_parquet('%s') ORDER BY pos", revel))
stopifnot(identical(revel_rows$chrom, c("1", "X")), identical(revel_rows$pos, c(1, 2)))
DBI::dbDisconnect(con, shutdown = TRUE)
RS

echo "VariantKey provider staging: OK"
