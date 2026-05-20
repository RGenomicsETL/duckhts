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
#   LIFTOVER_OUT_DIR       output directory (default: /tmp/duckhts_liftover_cases)
#   LIFTOVER_CHAIN_PATH    chain file (default: /root/GRCh37/GRCh37_to_GRCh38.chain.gz)
#   LIFTOVER_SRC_FASTA     source FASTA (default: /root/GRCh37/human_g1k_v37.fasta)
#   LIFTOVER_DST_FASTA     destination FASTA (default: /root/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna)
#   LIFTOVER_REGION_OVERRIDE optional region string applied to every selected case
#   DUCKHTS_EXT, BCFTOOLS_BIN, BCFTOOLS_PLUGIN_DIR, KEEP_SLICE are forwarded to
#   scripts/liftover_conformance.sh unchanged.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CASES_TSV="${LIFTOVER_CASES_TSV:-$SCRIPT_DIR/conformance_case_table.tsv}"
OUT_DIR="${LIFTOVER_OUT_DIR:-/tmp/duckhts_liftover_cases}"
CHAIN_PATH="${LIFTOVER_CHAIN_PATH:-/root/GRCh37/GRCh37_to_GRCh38.chain.gz}"
SRC_FASTA="${LIFTOVER_SRC_FASTA:-/root/GRCh37/human_g1k_v37.fasta}"
DST_FASTA="${LIFTOVER_DST_FASTA:-/root/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna}"
REGION_OVERRIDE="${LIFTOVER_REGION_OVERRIDE:-}"
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

while IFS=$'\t' read -r case_id dataset sample description default_region input_vcf; do
  [[ -n "$case_id" ]] || continue
  [[ "$case_id" == "case_id" ]] && continue
  if ! case_requested "$case_id" "$@"; then
    continue
  fi

  region="$default_region"
  if [[ -n "$REGION_OVERRIDE" ]]; then
    region="$REGION_OVERRIDE"
  fi

  out_prefix="$OUT_DIR/$case_id"
  run_count=$((run_count + 1))
  echo "==> [$case_id] $description"
  if [[ -n "$region" ]]; then
    LIFTOVER_REGION="$region" bash "$SCRIPT_DIR/liftover_conformance.sh" \
      "$input_vcf" "$CHAIN_PATH" "$SRC_FASTA" "$DST_FASTA" "$out_prefix"
  else
    bash "$SCRIPT_DIR/liftover_conformance.sh" \
      "$input_vcf" "$CHAIN_PATH" "$SRC_FASTA" "$DST_FASTA" "$out_prefix"
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
