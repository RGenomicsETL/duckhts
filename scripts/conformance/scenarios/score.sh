#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:?output directory required}"
ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
source "${ROOT_DIR}/scripts/conformance/common.sh"

conformance_prepare_dirs "$OUT_DIR"
SUMMARY_TSV="${OUT_DIR}/report/summary.tsv"
conformance_write_summary_header "$SUMMARY_TSV"

run_score_case() {
  local case_name="$1"
  local vcf_rel="$2"
  local summary_rel="$3"
  local use_tag="$4"
  local columns_preset="$5"
  local q_score_thr="$6"
  local counts="$7"
  local use_variant_id="$8"

  local prefix="${OUT_DIR}/raw/${case_name}"
  local compare_tsv="${prefix}.compare.tsv"
  local detail_tsv="${OUT_DIR}/report/${case_name}.details.tsv"
  local log_path="${OUT_DIR}/raw/${case_name}.log"
  local mismatch_count
  local reason

  (
    cd "$ROOT_DIR"
    env \
      DUCKHTS_EXT="${DUCKHTS_EXT}" \
      BCFTOOLS_BIN="${BCFTOOLS_BIN}" \
      BCFTOOLS_PLUGINS="${BCFTOOLS_PLUGINS}" \
      USE_TAG="${use_tag}" \
      COLUMNS_PRESET="${columns_preset}" \
      Q_SCORE_THR="${q_score_thr}" \
      COUNTS="${counts}" \
      USE_VARIANT_ID="${use_variant_id}" \
      bash scripts/score_conformance.sh "$vcf_rel" "$summary_rel" "$prefix"
  ) >"$log_path" 2>&1

  duckdb -unsigned <<SQL
COPY (
  SELECT *,
         CASE
           WHEN status = 'diff_values'
             AND abs(duck_value - bcf_value) <= 1e-6 THEN 'float_precision'
           ELSE 'unclassified'
         END AS reason
  FROM read_csv_auto('${compare_tsv}', delim := '\t', header := true)
  WHERE status NOT IN ('match', 'both_null')
  ORDER BY SAMPLE, metric
)
TO '${detail_tsv}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
SQL

  mismatch_count="$(duckdb -csv -noheader -unsigned -c "SELECT COUNT(*) FROM read_csv_auto('${detail_tsv}', delim := '\t', header := true)")"
  if [[ "$mismatch_count" == "0" ]]; then
    conformance_append_summary_row "$SUMMARY_TSV" "$case_name" "match" "match" "0" "${detail_tsv}"
    return
  fi

  reason="$(duckdb -csv -noheader -unsigned -c "SELECT CASE WHEN COUNT(DISTINCT reason) = 1 THEN MIN(reason) ELSE 'mixed' END FROM read_csv_auto('${detail_tsv}', delim := '\t', header := true)")"
  conformance_append_summary_row "$SUMMARY_TSV" "$case_name" "diff" "$reason" "$mismatch_count" "${detail_tsv}"
}

run_score_case "basic_gt" "test/data/score_input.vcf" "test/data/score_summary.tsv" "GT" "PLINK" "" "0" "0"
run_score_case "qthr_counts" "test/data/score_input.vcf" "test/data/score_summary.tsv" "GT" "PLINK" "0.01,0.2" "1" "0"
run_score_case "gp_dosage" "test/data/score_dosage.vcf" "test/data/score_summary.tsv" "GP" "PLINK" "" "0" "0"
