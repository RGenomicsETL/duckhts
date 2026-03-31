#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${ROOT_DIR}/scripts/conformance/common.sh"

SCENARIO="${1:-all}"
OUT_DIR="${2:-${ROOT_DIR}/.tmp/conformance/${SCENARIO}}"
SCENARIO_DIR="${ROOT_DIR}/scripts/conformance/scenarios"

conformance_require_cmd duckdb
conformance_require_cmd bcftools

DUCKHTS_EXT="${DUCKHTS_EXT:-${ROOT_DIR}/build/release/duckhts.duckdb_extension}"
conformance_require_file "$DUCKHTS_EXT"

PLUGIN_DIR="${BCFTOOLS_PLUGINS:-$(conformance_resolve_plugin_dir "$ROOT_DIR")}"

mkdir -p "$OUT_DIR"

run_one() {
  local scenario_name="$1"
  local script_path="${SCENARIO_DIR}/${scenario_name}.sh"
  local scenario_out="${OUT_DIR}/${scenario_name}"

  conformance_require_file "$script_path"
  mkdir -p "$scenario_out"

  DUCKHTS_EXT="$DUCKHTS_EXT" \
  BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}" \
  BCFTOOLS_PLUGINS="$PLUGIN_DIR" \
  ROOT_DIR="$ROOT_DIR" \
  "$script_path" "$scenario_out"
}

if [[ "$SCENARIO" == "all" ]]; then
  while IFS= read -r scenario_file; do
    run_one "$(basename "$scenario_file" .sh)"
  done < <(find "$SCENARIO_DIR" -maxdepth 1 -type f -name '*.sh' | sort)
else
  run_one "$SCENARIO"
fi

SUMMARY_OUT="${OUT_DIR}/summary.tsv"
printf 'scenario\tcase\tstatus\treason\tmismatch_count\tdetail_path\n' > "$SUMMARY_OUT"
while IFS= read -r summary_file; do
  scenario_name="$(basename "$(dirname "$(dirname "$summary_file")")")"
  awk -F'\t' -v OFS='\t' -v scenario="$scenario_name" 'NR > 1 { print scenario, $0 }' "$summary_file" >> "$SUMMARY_OUT"
done < <(find "$OUT_DIR" -mindepth 3 -maxdepth 3 -path '*/report/summary.tsv' | sort)

echo "Conformance summary: ${SUMMARY_OUT}"
column -t -s $'\t' "$SUMMARY_OUT" || cat "$SUMMARY_OUT"
