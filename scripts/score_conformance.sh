#!/usr/bin/env bash
set -euo pipefail

# Score conformance test: compare DuckHTS bcftools_score() against
# upstream bcftools +score on the same VCF + summary input.
#
# Usage:
#   bash scripts/score_conformance.sh [VCF_PATH] [SUMMARY_PATH]
#
# The script:
#  1. Prepares an indexed VCF (bgzip + tabix) for bcftools
#  2. Runs DuckHTS bcftools_score() -> TSV
#  3. Runs bcftools +score -> TSV
#  4. Compares per-sample scores with numeric tolerance via FULL OUTER JOIN
#
# Environment variables:
#   DUCKHTS_EXT          path to duckhts.duckdb_extension
#   BCFTOOLS_BIN         path to bcftools binary
#   SCORE_PLUGIN_DIR     directory containing score.so
#   USE_TAG              GT / DS / HDS / AP / GP  (default: GT)
#   COLUMNS_PRESET       PLINK / PLINK2 / METAL / SSF / … (default: PLINK)
#   Q_SCORE_THR          comma-separated p-value thresholds (default: empty)
#   USE_VARIANT_ID       1 to match on rsID (default: 0)
#   COUNTS               1 to include CNT columns (default: 0)
#   NUM_TOL              absolute numeric tolerance (default: 1e-4)
#   REL_TOL              relative numeric tolerance (default: 1e-5)
#   OUT_PREFIX           output file prefix (default: score_conformance)

VCF_PATH="${1:-test/data/score_input.vcf}"
SUMMARY_PATH="${2:-test/data/score_summary.tsv}"
OUT_PREFIX="${3:-score_conformance}"
DUCKHTS_EXT="${DUCKHTS_EXT:-build/release/duckhts.duckdb_extension}"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
SCORE_PLUGIN_DIR="${SCORE_PLUGIN_DIR:-}"
USE_TAG="${USE_TAG:-GT}"
COLUMNS_PRESET="${COLUMNS_PRESET:-PLINK}"
Q_SCORE_THR="${Q_SCORE_THR:-}"
USE_VARIANT_ID="${USE_VARIANT_ID:-0}"
COUNTS="${COUNTS:-0}"
NUM_TOL="${NUM_TOL:-1e-4}"
REL_TOL="${REL_TOL:-1e-5}"

DUCK_TSV="${OUT_PREFIX}.duckhts.tsv"
BCF_TSV="${OUT_PREFIX}.bcftools.tsv"
COMPARE_TSV="${OUT_PREFIX}.compare.tsv"

# --- Validation ---

if [[ ! -f "$VCF_PATH" ]]; then
  echo "Input VCF not found: $VCF_PATH" >&2; exit 1
fi
if [[ ! -f "$SUMMARY_PATH" ]]; then
  echo "Summary file not found: $SUMMARY_PATH" >&2; exit 1
fi
if [[ ! -f "$DUCKHTS_EXT" ]]; then
  echo "DuckHTS extension not found: $DUCKHTS_EXT" >&2; exit 1
fi
if ! command -v "$BCFTOOLS_BIN" >/dev/null 2>&1; then
  echo "bcftools not found: $BCFTOOLS_BIN" >&2; exit 1
fi

# --- Resolve score plugin ---

if [[ -z "$SCORE_PLUGIN_DIR" ]]; then
  if [[ -f "score_1.22-20250819.zip" ]]; then
    mkdir -p .tmp/score-plugin
    if [[ ! -f .tmp/score-plugin/score.so ]]; then
      unzip -oq score_1.22-20250819.zip -d .tmp/score-plugin
    fi
    SCORE_PLUGIN_DIR="$(realpath .tmp/score-plugin)"
  else
    echo "Set SCORE_PLUGIN_DIR or provide score_1.22-20250819.zip" >&2
    exit 1
  fi
fi

if [[ ! -f "$SCORE_PLUGIN_DIR/score.so" ]]; then
  echo "score.so not found under SCORE_PLUGIN_DIR: $SCORE_PLUGIN_DIR" >&2
  exit 1
fi

# --- Prepare indexed VCF for bcftools ---
# bcftools +score requires indexed VCF/BCF (BCF_SR_REQUIRE_IDX).

VCF_GZ=""
CLEANUP_GZ=0

if [[ "$VCF_PATH" == *.vcf.gz ]]; then
  VCF_GZ="$VCF_PATH"
  if [[ ! -f "${VCF_GZ}.tbi" && ! -f "${VCF_GZ}.csi" ]]; then
    "$BCFTOOLS_BIN" index -t "$VCF_GZ"
  fi
elif [[ "$VCF_PATH" == *.bcf ]]; then
  VCF_GZ="$VCF_PATH"
  if [[ ! -f "${VCF_GZ}.csi" ]]; then
    "$BCFTOOLS_BIN" index "$VCF_GZ"
  fi
else
  VCF_GZ="${OUT_PREFIX}.input.vcf.gz"
  "$BCFTOOLS_BIN" view -Oz -o "$VCF_GZ" "$VCF_PATH"
  "$BCFTOOLS_BIN" index -t "$VCF_GZ"
  CLEANUP_GZ=1
fi

cleanup() {
  if [[ "$CLEANUP_GZ" == "1" ]]; then
    rm -f "$VCF_GZ" "${VCF_GZ}.tbi" "${VCF_GZ}.csi"
  fi
}
trap cleanup EXIT

# --- Build DuckHTS named parameters ---

DUCK_OPTS="use := '${USE_TAG}', columns := '${COLUMNS_PRESET}'"
if [[ "$USE_VARIANT_ID" == "1" ]]; then
  DUCK_OPTS="${DUCK_OPTS}, use_variant_id := true"
fi
if [[ -n "$Q_SCORE_THR" ]]; then
  DUCK_OPTS="${DUCK_OPTS}, q_score_thr := '${Q_SCORE_THR}'"
fi
if [[ "$COUNTS" == "1" ]]; then
  DUCK_OPTS="${DUCK_OPTS}, counts := true"
fi

# --- Build bcftools CLI args ---

BCF_ARGS=(--use "$USE_TAG" -c "$COLUMNS_PRESET")
if [[ "$USE_VARIANT_ID" == "1" ]]; then
  BCF_ARGS+=(--use-variant-id)
fi
if [[ -n "$Q_SCORE_THR" ]]; then
  BCF_ARGS+=(--q-score-thr "$Q_SCORE_THR")
fi
if [[ "$COUNTS" == "1" ]]; then
  BCF_ARGS+=(--counts)
fi

# Resolve to absolute paths for SQL
VCF_ABS="$(realpath "$VCF_PATH")"
SUMMARY_ABS="$(realpath "$SUMMARY_PATH")"

# --- [1/4] Run DuckHTS bcftools_score ---

echo "[1/4] DuckHTS bcftools_score -> $DUCK_TSV"
duckdb -unsigned <<SQL
LOAD '${DUCKHTS_EXT}';

COPY (
  SELECT *
  FROM bcftools_score(
    '${VCF_ABS}',
    '${SUMMARY_ABS}',
    ${DUCK_OPTS}
  )
  ORDER BY SAMPLE
)
TO '${DUCK_TSV}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
SQL

# --- [2/4] Run bcftools +score ---

echo "[2/4] bcftools +score -> $BCF_TSV"
BCFTOOLS_PLUGINS="$SCORE_PLUGIN_DIR" \
  "$BCFTOOLS_BIN" +score "${BCF_ARGS[@]}" -o "$BCF_TSV" "$VCF_GZ" "$SUMMARY_ABS"

# --- [3/4] Normalize and compare ---
# Both produce TSV: rows = samples, columns = SAMPLE + score metrics.
# Both derive PRS name from summary filename (basename minus extensions),
# so column names should already match. UNPIVOT both to long format and
# FULL OUTER JOIN on SAMPLE + metric with numeric tolerance.

echo "[3/4] Compare -> $COMPARE_TSV"
duckdb -unsigned <<SQL
LOAD '${DUCKHTS_EXT}';

CREATE TEMP TABLE duck_raw AS
  SELECT * FROM read_csv_auto('${DUCK_TSV}', delim := '\t', header := true);

CREATE TEMP TABLE bcf_raw AS
  SELECT * FROM read_csv_auto('${BCF_TSV}', delim := '\t', header := true);

-- UNPIVOT to long format: (SAMPLE, metric, value)
CREATE TEMP TABLE duck_long AS
  SELECT SAMPLE, metric, CAST(value AS DOUBLE) AS value
  FROM (UNPIVOT duck_raw ON COLUMNS(* EXCLUDE(SAMPLE)) INTO NAME metric VALUE value);

CREATE TEMP TABLE bcf_long AS
  SELECT SAMPLE, metric, CAST(value AS DOUBLE) AS value
  FROM (UNPIVOT bcf_raw ON COLUMNS(* EXCLUDE(SAMPLE)) INTO NAME metric VALUE value);

-- Column names match between engines (both derive PRS name from filename).
-- No normalization needed.

-- FULL OUTER JOIN with numeric tolerance
COPY (
  SELECT
    COALESCE(d.SAMPLE, b.SAMPLE) AS SAMPLE,
    COALESCE(d.metric, b.metric) AS metric,
    d.value AS duck_value,
    b.value AS bcf_value,
    CASE
      WHEN d.SAMPLE IS NULL THEN 'only_bcftools'
      WHEN b.SAMPLE IS NULL THEN 'only_duckhts'
      WHEN d.value IS NULL AND b.value IS NULL THEN 'both_null'
      WHEN d.value IS NULL OR b.value IS NULL THEN 'null_mismatch'
      WHEN abs(d.value - b.value) <= ${NUM_TOL}
        OR abs(d.value - b.value) <= ${REL_TOL} * greatest(abs(d.value), abs(b.value), 1.0)
        THEN 'match'
      ELSE 'diff_values'
    END AS status
  FROM duck_long d
  FULL OUTER JOIN bcf_long b
    ON d.SAMPLE = b.SAMPLE AND d.metric = b.metric
  ORDER BY SAMPLE, metric
)
TO '${COMPARE_TSV}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);

-- Summary
.mode duckbox
SELECT status, COUNT(*) AS n
FROM read_csv_auto('${COMPARE_TSV}', delim := '\t', header := true)
GROUP BY status
ORDER BY status;

-- Show mismatches (if any)
SELECT *
FROM read_csv_auto('${COMPARE_TSV}', delim := '\t', header := true)
WHERE status NOT IN ('match', 'both_null')
LIMIT 20;
SQL

echo ""
echo "Settings:"
echo "  use_tag:          ${USE_TAG}"
echo "  columns_preset:   ${COLUMNS_PRESET}"
echo "  q_score_thr:      ${Q_SCORE_THR:-none}"
echo "  use_variant_id:   ${USE_VARIANT_ID}"
echo "  counts:           ${COUNTS}"
echo "  numeric_tolerance: ${NUM_TOL}"
echo "  relative_tolerance: ${REL_TOL}"
echo ""
echo "Output files:"
echo "  $DUCK_TSV"
echo "  $BCF_TSV"
echo "  $COMPARE_TSV"
