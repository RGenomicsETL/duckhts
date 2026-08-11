#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/duckhts-corpus-stage-test.XXXXXX")"
SOURCE_DIR="$TMP_DIR/sources"
FAKE_BIN="$TMP_DIR/bin"
CACHE_DIR="$TMP_DIR/cache root"
REGISTRY="$TMP_DIR/registry.tsv"
TOOL_LOG="$TMP_DIR/tools.log"
trap 'rm -rf "$TMP_DIR"' EXIT
mkdir -p "$SOURCE_DIR" "$FAKE_BIN"

cat >"$FAKE_BIN/bcftools" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" >>"${FAKE_TOOL_LOG:?}"
if [[ "${1:-}" == "--version" ]]; then
  printf 'bcftools fixture 1.0\nUsing htslib fixture 1.0\n'
  exit 0
fi
command="$1"
out=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o) out="$2"; shift 2 ;;
    *) shift ;;
  esac
done
if [[ -z "$out" ]]; then
  cat >/dev/null
  printf 'BCF\n'
  exit 0
fi
if [[ "$command" == "annotate" ]]; then cat >/dev/null; fi
mkdir -p "$(dirname "$out")"
case "$out" in
  *.bcf) printf 'BCF\n' >"$out" ;;
  *hprc_v2_grch38*) printf 'HPRC\n' >"$out" ;;
  *sniffles2_joint*) printf 'SNIFFLES2\n' >"$out" ;;
  *GRCh38.variant_region*) printf 'DBVAR\n' >"$out" ;;
  *) printf 'unexpected %s output\n' "$command" >&2; exit 2 ;;
esac
EOF

cat >"$FAKE_BIN/tabix" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" >>"${FAKE_TOOL_LOG:?}"
for argument in "$@"; do target="$argument"; done
printf 'INDEX\n' >"$target.tbi"
EOF

cat >"$FAKE_BIN/curl" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" >>"${FAKE_TOOL_LOG:?}"
printf 'HTTP/1.1 200 OK\r\nETag: "%s"\r\n' "${FAKE_SNIFFLES_ETAG:-fixture-sniffles-etag}"
EOF
chmod +x "$FAKE_BIN"/*

mkdir -p "$SOURCE_DIR/hprc"
printf 'HPRC SOURCE INDEX\n' >"$SOURCE_DIR/hprc/versionId=fixture-hprc-index"
printf 'SNIFFLES SOURCE INDEX\n' >"$SOURCE_DIR/sniffles.csi"
printf 'DBVAR SOURCE INDEX\n' >"$SOURCE_DIR/dbvar.tbi"
printf 'fixture-dbvar-md5  GRCh38.variant_region.all.vcf.gz\n' >"$SOURCE_DIR/dbvar.md5"

Rscript - "$ROOT_DIR" "$SOURCE_DIR" "$REGISTRY" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
root <- args[[1L]]
source_dir <- args[[2L]]
registry_path <- args[[3L]]
registry <- utils::read.delim(
  file.path(root, "r", "duckhtsbench", "inst", "benchmark_registry.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
registry <- registry[registry$workload == "duckvep-conformance-corpora", , drop = FALSE]
sha256 <- function(path) {
  command <- Sys.which("sha256sum")
  arguments <- shQuote(path)
  if (!nzchar(command)) {
    command <- Sys.which("shasum")
    arguments <- c("-a", "256", shQuote(path))
  }
  strsplit(trimws(system2(command, arguments, stdout = TRUE)), "[[:space:]]+")[[1L]][[1L]]
}
identity <- function(path, extra = character()) {
  paste(c(paste0("sha256=", sha256(path)), paste0("bytes=", file.info(path)$size), extra), collapse = ";")
}
file_url <- function(path) paste0("file://", normalizePath(path))

hprc_source <- registry$id == "duckvep_hprc_v2_grch38_source"
hprc_index <- registry$id == "duckvep_hprc_v2_grch38_source_tbi"
registry$locator[hprc_source] <- paste0(file_url(source_dir), "/hprc/versionId=fixture-hprc-vcf")
registry$supplier_identity[hprc_source] <- "version_id=fixture-hprc-vcf;bytes=999"
registry$locator[hprc_index] <- file_url(file.path(source_dir, "hprc", "versionId=fixture-hprc-index"))
registry$supplier_identity[hprc_index] <- identity(
  file.path(source_dir, "hprc", "versionId=fixture-hprc-index"),
  "version_id=fixture-hprc-index"
)

sniffles_source <- registry$id == "duckvep_sniffles2_1kgp_source"
sniffles_index <- registry$id == "duckvep_sniffles2_1kgp_source_csi"
registry$locator[sniffles_source] <- "https://fixtures.invalid/sniffles2_joint_sv_calls.vcf.gz"
registry$supplier_identity[sniffles_source] <- "etag=fixture-sniffles-etag;bytes=999;source_version=Sniffles2_2.0.7"
registry$locator[sniffles_index] <- file_url(file.path(source_dir, "sniffles.csi"))
registry$supplier_identity[sniffles_index] <- identity(file.path(source_dir, "sniffles.csi"))

dbvar_source <- registry$id == "duckvep_dbvar_grch38_20260127_source"
dbvar_index <- registry$id == "duckvep_dbvar_grch38_20260127_source_tbi"
dbvar_manifest <- registry$id == "duckvep_dbvar_grch38_20260127_manifest"
registry$locator[dbvar_source] <- "https://fixtures.invalid/GRCh38.variant_region.all.vcf.gz"
registry$supplier_identity[dbvar_source] <- "md5=fixture-dbvar-md5;file_date=20260127;bytes=999"
registry$locator[dbvar_index] <- file_url(file.path(source_dir, "dbvar.tbi"))
registry$supplier_identity[dbvar_index] <- identity(file.path(source_dir, "dbvar.tbi"))
registry$locator[dbvar_manifest] <- file_url(file.path(source_dir, "dbvar.md5"))
registry$supplier_identity[dbvar_manifest] <- "required_entry_md5=fixture-dbvar-md5"

outputs <- c(
  duckvep_hprc_v2_grch38_chr22_african4 = "HPRC\n",
  duckvep_sniffles2_1kgp_chr22 = "SNIFFLES2\n",
  duckvep_dbvar_grch38_20260127_chr22 = "DBVAR\n"
)
for (id in names(outputs)) {
  path <- file.path(source_dir, paste0(id, ".expected"))
  writeChar(outputs[[id]], path, eos = NULL, useBytes = TRUE)
  fields <- strsplit(registry$supplier_identity[registry$id == id], ";", fixed = TRUE)[[1L]]
  metadata <- fields[!grepl("^(sha256|bytes)=", fields)]
  registry$supplier_identity[registry$id == id] <- identity(path, metadata)
}
index_expected <- file.path(source_dir, "index.expected")
writeChar("INDEX\n", index_expected, eos = NULL, useBytes = TRUE)
for (id in registry$id[registry$role == "derived_chr22_vcf_index"]) {
  registry$supplier_identity[registry$id == id] <- identity(index_expected)
}
utils::write.table(registry, registry_path, sep = "\t", row.names = FALSE, quote = FALSE)

wrong_version <- registry
wrong_version$supplier_identity[wrong_version$id == "duckvep_hprc_v2_grch38_source"] <-
  "version_id=wrong-version;bytes=999"
utils::write.table(wrong_version, paste0(registry_path, ".wrong-version"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
wrong_output <- registry
wrong_output$supplier_identity[wrong_output$id == "duckvep_hprc_v2_grch38_chr22_african4"] <- sub(
  "sha256=[0-9a-f]+", paste0("sha256=", paste(rep("0", 64L), collapse = "")),
  wrong_output$supplier_identity[wrong_output$id == "duckvep_hprc_v2_grch38_chr22_african4"]
)
utils::write.table(wrong_output, paste0(registry_path, ".wrong-output"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
RS

STAGE="$ROOT_DIR/r/duckhtsbench/scripts/stage_duckvep_conformance_corpora.R"
FAKE_TOOL_LOG="$TOOL_LOG" DUCKHTS_CACHE_DIR="$CACHE_DIR" DUCKHTSBENCH_REGISTRY="$REGISTRY" \
  Rscript "$STAGE" --corpus all --bcftools "$FAKE_BIN/bcftools" \
  --tabix "$FAKE_BIN/tabix" --curl "$FAKE_BIN/curl" >/dev/null

Rscript - "$REGISTRY" "$CACHE_DIR" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
registry <- utils::read.delim(args[[1L]], stringsAsFactors = FALSE, check.names = FALSE)
paths <- file.path(args[[2L]], registry$cache_relpath)
stopifnot(all(file.exists(paths)), all(file.info(paths)$size > 0L))
stopifnot(all(file.exists(paste0(paths, ".provenance.tsv"))))
derived <- paths[registry$role == "derived_chr22_vcf"]
stopifnot(all(file.exists(paste0(derived, ".receipt.tsv"))))
for (path in derived) {
  receipt <- utils::read.delim(paste0(path, ".provenance.tsv"), stringsAsFactors = FALSE)
  stopifnot(all(c(
    "source_url", "source_id", "source_artifacts", "source_releases",
    "source_locators", "source_identities"
  ) %in% receipt$field))
}
RS

HPRC_INDEX="$CACHE_DIR/corpora/duckvep/source-indexes/hprc-v2.0-mc-grch38.vcf.gz.tbi"
grep -F "versionId=fixture-hprc-vcf##idx##$HPRC_INDEX" "$TOOL_LOG" >/dev/null
grep -F -- "-s HG02055,HG02145,HG02723,HG03098" "$TOOL_LOG" >/dev/null
grep -F "https://fixtures.invalid/sniffles2_joint_sv_calls.vcf.gz##idx##$CACHE_DIR/corpora/duckvep/source-indexes/sniffles2_joint_sv_calls.vcf.gz.csi" "$TOOL_LOG" >/dev/null
grep -F "https://fixtures.invalid/GRCh38.variant_region.all.vcf.gz##idx##$CACHE_DIR/corpora/duckvep/source-indexes/GRCh38.variant_region.all.vcf.gz.tbi" "$TOOL_LOG" >/dev/null

printf 'stale source index\n' >"$HPRC_INDEX"
FAKE_TOOL_LOG="$TOOL_LOG" DUCKHTS_CACHE_DIR="$CACHE_DIR" DUCKHTSBENCH_REGISTRY="$REGISTRY" \
  Rscript "$STAGE" --corpus hprc-african4-chr22 --bcftools "$FAKE_BIN/bcftools" \
  --tabix "$FAKE_BIN/tabix" --curl "$FAKE_BIN/curl" >/dev/null
cmp "$HPRC_INDEX" "$SOURCE_DIR/hprc/versionId=fixture-hprc-index"

if FAKE_TOOL_LOG="$TOOL_LOG" DUCKHTS_CACHE_DIR="$CACHE_DIR" \
  DUCKHTSBENCH_REGISTRY="$REGISTRY.wrong-version" Rscript "$STAGE" \
  --corpus hprc-african4-chr22 --bcftools "$FAKE_BIN/bcftools" \
  --tabix "$FAKE_BIN/tabix" --curl "$FAKE_BIN/curl" >/dev/null 2>&1; then
  echo "wrong HPRC source version unexpectedly passed" >&2
  exit 1
fi
if FAKE_SNIFFLES_ETAG=wrong-etag FAKE_TOOL_LOG="$TOOL_LOG" DUCKHTS_CACHE_DIR="$CACHE_DIR" \
  DUCKHTSBENCH_REGISTRY="$REGISTRY" Rscript "$STAGE" --corpus sniffles2-chr22 \
  --bcftools "$FAKE_BIN/bcftools" --tabix "$FAKE_BIN/tabix" \
  --curl "$FAKE_BIN/curl" >/dev/null 2>&1; then
  echo "wrong Sniffles2 ETag unexpectedly passed" >&2
  exit 1
fi
printf 'wrong  GRCh38.variant_region.all.vcf.gz\n' >"$SOURCE_DIR/dbvar.md5"
if FAKE_TOOL_LOG="$TOOL_LOG" DUCKHTS_CACHE_DIR="$CACHE_DIR" DUCKHTSBENCH_REGISTRY="$REGISTRY" \
  Rscript "$STAGE" --corpus dbvar-chr22 --bcftools "$FAKE_BIN/bcftools" \
  --tabix "$FAKE_BIN/tabix" --curl "$FAKE_BIN/curl" >/dev/null 2>&1; then
  echo "wrong dbVar manifest unexpectedly passed" >&2
  exit 1
fi
printf 'fixture-dbvar-md5  GRCh38.variant_region.all.vcf.gz\n' >"$SOURCE_DIR/dbvar.md5"

HPRC_OUTPUT="$CACHE_DIR/corpora/duckvep/hprc_v2_grch38/hprc_v2_grch38_chr22_african4_carried.vcf.gz"
BEFORE="$(sha256sum "$HPRC_OUTPUT" | awk '{print $1}')"
if FAKE_TOOL_LOG="$TOOL_LOG" DUCKHTS_CACHE_DIR="$CACHE_DIR" \
  DUCKHTSBENCH_REGISTRY="$REGISTRY.wrong-output" Rscript "$STAGE" \
  --corpus hprc-african4-chr22 --bcftools "$FAKE_BIN/bcftools" \
  --tabix "$FAKE_BIN/tabix" --curl "$FAKE_BIN/curl" >/dev/null 2>&1; then
  echo "wrong derived identity unexpectedly passed" >&2
  exit 1
fi
AFTER="$(sha256sum "$HPRC_OUTPUT" | awk '{print $1}')"
[[ "$BEFORE" == "$AFTER" ]]

LEGACY_DIR="$TMP_DIR/legacy output"
FAKE_TOOL_LOG="$TOOL_LOG" DUCKHTS_CACHE_DIR="$CACHE_DIR" DUCKHTSBENCH_REGISTRY="$REGISTRY" \
  BCFTOOLS_BIN="$FAKE_BIN/bcftools" TABIX_BIN="$FAKE_BIN/tabix" CURL_BIN="$FAKE_BIN/curl" \
  bash "$ROOT_DIR/scripts/stage_duckvep_conformance_corpora.sh" \
  "$LEGACY_DIR" dbvar-chr22 >/dev/null
[[ -s "$LEGACY_DIR/dbvar/GRCh38_20260127/GRCh38.variant_region.all.chr22.vcf.gz" ]]
[[ -s "$LEGACY_DIR/dbvar/GRCh38_20260127/GRCh38.variant_region.all.chr22.vcf.gz.tbi" ]]

echo "DuckVEP corpus registry staging: OK"
