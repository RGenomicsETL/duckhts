#!/usr/bin/env bash
set -euo pipefail

# Run scripts/liftover_conformance.sh across a curated case table.
# Default case table currently contains GIAB benchmark VCF slices but the TSV
# schema is generic enough for future known-sites/site-only VCF sources.
#
# Usage:
#   bash scripts/liftover_conformance_batch.sh [CASE_ID ...]
#
# Examples:
#   bash scripts/liftover_conformance_batch.sh
#   bash scripts/liftover_conformance_batch.sh giab_hg001_grch37_chr20 giab_hg006_grch37_chr20
#   LIFTOVER_REGION_OVERRIDE=1:1-5000000 bash scripts/liftover_conformance_batch.sh giab_hg002_grch37_chr20
#
# Environment variables:
#   LIFTOVER_CASES_TSV     case table (default: scripts/conformance_case_table.tsv)
#   LIFTOVER_OUT_DIR       output directory (default: $DUCKHTS_CACHE_DIR/conformance/liftover)
#   LIFTOVER_CHAIN_PATH    chain file (default: staged liftover cache path)
#   LIFTOVER_SRC_FASTA     source FASTA (default: staged liftover cache path)
#   LIFTOVER_DST_FASTA     destination FASTA (default: staged liftover cache path)
#                           Run scripts/stage_liftover_references.sh first.
#   LIFTOVER_REGION_OVERRIDE optional region string applied to every selected case
#   LIFTOVER_STAGE_ONLY=1 stage selected registered inputs without running comparison
#   DUCKHTS_EXT, BCFTOOLS_BIN, BCFTOOLS_PLUGIN_DIR, KEEP_SLICE are forwarded to
#   scripts/liftover_conformance.sh unchanged.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=duckhts_cache.sh
source "$SCRIPT_DIR/duckhts_cache.sh"
LIFTOVER_REFERENCE_DIR=${LIFTOVER_REFERENCE_DIR:-$(duckhts_cache_subdir "references/liftover/grch37-to-grch38")}
CASES_TSV="${LIFTOVER_CASES_TSV:-$SCRIPT_DIR/conformance_case_table.tsv}"
REGISTRY_TSV="${DUCKHTSBENCH_REGISTRY:-$REPO_ROOT/r/duckhtsbench/inst/benchmark_registry.tsv}"
OUT_DIR="${LIFTOVER_OUT_DIR:-$(duckhts_cache_subdir conformance/liftover)}"
CHAIN_PATH="${LIFTOVER_CHAIN_PATH:-$LIFTOVER_REFERENCE_DIR/GRCh37_to_GRCh38.chain.gz}"
SRC_FASTA="${LIFTOVER_SRC_FASTA:-$LIFTOVER_REFERENCE_DIR/human_g1k_v37.fasta}"
DST_FASTA="${LIFTOVER_DST_FASTA:-$LIFTOVER_REFERENCE_DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna}"
REGION_OVERRIDE="${LIFTOVER_REGION_OVERRIDE:-}"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
SUMMARY_TSV="$OUT_DIR/summary.tsv"
run_count=0

if [[ ! -f "$CASES_TSV" ]]; then
  echo "Case table not found: $CASES_TSV" >&2
  exit 1
fi
if [[ ! -f "$CHAIN_PATH" ]]; then
  echo "Chain file not found: $CHAIN_PATH" >&2
  exit 1
fi
if [[ ! -f "$SRC_FASTA" ]]; then
  echo "Source FASTA not found: $SRC_FASTA" >&2
  exit 1
fi
if [[ ! -f "$DST_FASTA" ]]; then
  echo "Destination FASTA not found: $DST_FASTA" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
printf 'case_id\tdataset\tsample\tregion\tmapped_statuses\treject_statuses\tmapped_mismatches\treject_mismatches\tcompare_tsv\treject_compare_tsv\n' > "$SUMMARY_TSV"

have_selection=0
if [[ $# -gt 0 ]]; then
  have_selection=1
fi

case_requested() {
  local wanted="$1"
  if [[ "$have_selection" -eq 0 ]]; then
    return 0
  fi
  shift || true
  for arg in "$@"; do
    if [[ "$arg" == "$wanted" ]]; then
      return 0
    fi
  done
  return 1
}

stage_case_input() { # case_id dataset sample artifact_id
  local case_id="$1"
  local dataset="$2"
  local sample="$3"
  local artifact_id="$4"
  local source cache_relpath transform cached index provenance row

  [[ -f "$REGISTRY_TSV" ]] || { echo "Benchmark registry not found: $REGISTRY_TSV" >&2; return 1; }
  row="$(awk -F '\t' -v id="$artifact_id" '$1 == id {print $5 "\t" $7 "\t" $8; exit}' "$REGISTRY_TSV")"
  [[ -n "$row" ]] || { echo "unknown liftover artifact: $artifact_id" >&2; return 1; }
  IFS=$'\t' read -r source cache_relpath transform <<<"$row"
  [[ "$transform" == "direct_download" ]] || { echo "liftover artifact is not a direct VCF download: $artifact_id" >&2; return 1; }
  cached="$(duckhts_cache_subdir "$cache_relpath")"
  index="${cached}.tbi"
  provenance="${cached}.provenance.tsv"
  mkdir -p "$(dirname "$cached")"

  if [[ ! -s "$cached" ]]; then
    command -v curl >/dev/null 2>&1 || { echo "curl is required to stage $case_id" >&2; return 1; }
    echo "Downloading $case_id from registered artifact $artifact_id" >&2
    curl --fail --location --retry 5 --retry-delay 5 --continue-at - --output "$cached" "$source"
  else
    echo "Using cached $case_id input: $cached" >&2
  fi

  if [[ ! -s "$index" ]]; then "$BCFTOOLS_BIN" index -f -t "$cached"; fi
  if [[ ! -f "$provenance" ]]; then
    duckhts_write_provenance "$provenance" \
      "workload=liftover_conformance" "case_id=$case_id" "dataset=$dataset" \
      "sample=$sample" "artifact_id=$artifact_id" "source=$source" \
      "cached_vcf=$cached" "cached_index=$index" \
      "staging_command=curl; bcftools index -t" \
      "bcftools=$($BCFTOOLS_BIN --version | head -n 1)"
  fi
  printf '%s\n' "$cached"
}

while IFS=$'\t' read -r case_id dataset sample description default_region artifact_id; do
  [[ -n "$case_id" ]] || continue
  [[ "$case_id" == "case_id" ]] && continue
  if ! case_requested "$case_id" "$@"; then
    continue
  fi

  region="$default_region"
  if [[ -n "$REGION_OVERRIDE" ]]; then
    region="$REGION_OVERRIDE"
  fi

  cached_vcf="$(stage_case_input "$case_id" "$dataset" "$sample" "$artifact_id")"
  run_count=$((run_count + 1))
  if [[ "${LIFTOVER_STAGE_ONLY:-0}" == "1" ]]; then
    echo "==> [$case_id] staged $cached_vcf"
    continue
  fi
  out_prefix="$OUT_DIR/$case_id"
  echo "==> [$case_id] $description"
  if [[ -n "$region" ]]; then
    LIFTOVER_REGION="$region" bash "$SCRIPT_DIR/liftover_conformance.sh" \
      "$cached_vcf" "$CHAIN_PATH" "$SRC_FASTA" "$DST_FASTA" "$out_prefix"
  else
    bash "$SCRIPT_DIR/liftover_conformance.sh" \
      "$cached_vcf" "$CHAIN_PATH" "$SRC_FASTA" "$DST_FASTA" "$out_prefix"
  fi

  compare_tsv="$out_prefix.compare.tsv"
  reject_compare_tsv="$out_prefix.reject_compare.tsv"
  mapped_statuses="$(duckdb -unsigned -csv -noheader -c "SELECT string_agg(status || ':' || CAST(n AS VARCHAR), ', ' ORDER BY status) FROM (SELECT status, COUNT(*) AS n FROM read_csv_auto('$compare_tsv', delim := '\t', header := true) GROUP BY status);")"
  reject_statuses="$(duckdb -unsigned -csv -noheader -c "SELECT string_agg(status || ':' || CAST(n AS VARCHAR), ', ' ORDER BY status) FROM (SELECT status, COUNT(*) AS n FROM read_csv_auto('$reject_compare_tsv', delim := '\t', header := true) GROUP BY status);")"
  mapped_mismatches="$(duckdb -unsigned -csv -noheader -c "SELECT COUNT(*) FROM read_csv_auto('$compare_tsv', delim := '\t', header := true) WHERE status <> 'match';")"
  reject_mismatches="$(duckdb -unsigned -csv -noheader -c "SELECT COUNT(*) FROM read_csv_auto('$reject_compare_tsv', delim := '\t', header := true) WHERE status <> 'match';")"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$case_id" "$dataset" "$sample" "${region:-full_input}" \
    "$mapped_statuses" "$reject_statuses" \
    "$mapped_mismatches" "$reject_mismatches" \
    "$compare_tsv" "$reject_compare_tsv" >> "$SUMMARY_TSV"
done < "$CASES_TSV"

if [[ "$run_count" -eq 0 ]]; then
  echo "No cases selected from: $CASES_TSV" >&2
  exit 1
fi

echo
echo "Batch summary:"
if command -v column >/dev/null 2>&1; then
  column -t -s $'\t' "$SUMMARY_TSV"
else
  cat "$SUMMARY_TSV"
fi
