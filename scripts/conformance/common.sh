#!/usr/bin/env bash
set -euo pipefail

conformance_repo_root() {
  local script_dir
  script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  (cd "${script_dir}/../.." && pwd)
}

conformance_require_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Missing required command: $cmd" >&2
    exit 1
  fi
}

conformance_require_file() {
  local path="$1"
  if [[ ! -f "$path" ]]; then
    echo "Missing required file: $path" >&2
    exit 1
  fi
}

conformance_resolve_plugin_dir() {
  local root="$1"
  local zip_path="${root}/score_1.22-20250819.zip"
  local out_dir="${root}/.tmp/score-plugin"

  mkdir -p "$out_dir"
  if [[ ! -f "${out_dir}/score.so" || ! -f "${out_dir}/munge.so" || ! -f "${out_dir}/liftover.so" ]]; then
    conformance_require_file "$zip_path"
    unzip -oq "$zip_path" -d "$out_dir"
  fi
  printf '%s\n' "$out_dir"
}

conformance_prepare_dirs() {
  local base_dir="$1"
  mkdir -p "${base_dir}/raw" "${base_dir}/normalized" "${base_dir}/report" "${base_dir}/tmp"
}

conformance_write_summary_header() {
  local summary_tsv="$1"
  printf 'case\tstatus\treason\tmismatch_count\tdetail_path\n' > "$summary_tsv"
}

conformance_append_summary_row() {
  local summary_tsv="$1"
  local case_name="$2"
  local status="$3"
  local reason="$4"
  local mismatch_count="$5"
  local detail_path="$6"
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$case_name" "$status" "$reason" "$mismatch_count" "$detail_path" >> "$summary_tsv"
}

