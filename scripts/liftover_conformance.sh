#!/usr/bin/env bash
set -euo pipefail

# Liftover conformance test: compare DuckHTS duckdb_liftover() against
# upstream bcftools +liftover on the same input VCF/BCF.
#
# Usage:
#   bash scripts/liftover_conformance.sh INPUT_VCF CHAIN_PATH SRC_FASTA DST_FASTA [OUT_PREFIX]
#
# Stage references first with scripts/stage_liftover_references.sh. The
# batch driver materializes each public GIAB slice below the DuckHTS cache.
#
# Environment variables:
#   DUCKHTS_EXT         path to duckhts.duckdb_extension
#   BCFTOOLS_BIN        path to bcftools binary (defaults to RBCFTools::bcftools_path())
#   BCFTOOLS_PLUGIN_DIR path to bcftools plugin directory
#                       (defaults to RBCFTools::bcftools_plugins_dir();
#                        SCORE_PLUGIN_DIR is also honored as a fallback)
#   LIFTOVER_REGION     optional bcftools region string to slice before comparison
#   KEEP_SLICE          1 to keep the materialized local VCF slice

if [[ $# -lt 4 || $# -gt 5 ]]; then
  echo "Usage: bash scripts/liftover_conformance.sh INPUT_VCF CHAIN_PATH SRC_FASTA DST_FASTA [OUT_PREFIX]" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=duckhts_cache.sh
source "$SCRIPT_DIR/duckhts_cache.sh"

INPUT_VCF="$1"
CHAIN_PATH="$2"
SRC_FASTA="$3"
DST_FASTA="$4"
OUT_PREFIX="${5:-$(duckhts_cache_subdir conformance/liftover/single/liftover_conformance)}"
DUCKHTS_EXT="${DUCKHTS_EXT:-build/release/duckhts.duckdb_extension}"
LIFTOVER_REGION="${LIFTOVER_REGION:-}"
KEEP_SLICE="${KEEP_SLICE:-0}"

mkdir -p "$(dirname "$OUT_PREFIX")"

if [[ -z "${BCFTOOLS_BIN:-}" ]]; then
  BCFTOOLS_BIN="$(Rscript -e "if (!requireNamespace('RBCFTools', quietly=TRUE)) quit(status=1); cat(RBCFTools::bcftools_path())")" || {
    echo "RBCFTools is required to resolve the bundled bcftools binary; set BCFTOOLS_BIN to override" >&2
    exit 1
  }
fi

PLUGIN_ENV_FALLBACK="${SCORE_PLUGIN_DIR:-}"
if [[ -z "${BCFTOOLS_PLUGIN_DIR:-}" ]]; then
  if [[ -n "$PLUGIN_ENV_FALLBACK" ]]; then
    BCFTOOLS_PLUGIN_DIR="$PLUGIN_ENV_FALLBACK"
  else
    BCFTOOLS_PLUGIN_DIR="$(Rscript -e "if (!requireNamespace('RBCFTools', quietly=TRUE)) quit(status=1); cat(RBCFTools::bcftools_plugins_dir())")" || {
      echo "RBCFTools is required to resolve the bundled bcftools plugin directory; set BCFTOOLS_PLUGIN_DIR to override" >&2
      exit 1
    }
  fi
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
if [[ ! -f "$DUCKHTS_EXT" ]]; then
  echo "DuckHTS extension not found: $DUCKHTS_EXT" >&2
  exit 1
fi
if [[ ! -f "$BCFTOOLS_PLUGIN_DIR/liftover.so" ]]; then
  echo "liftover.so not found under BCFTOOLS_PLUGIN_DIR: $BCFTOOLS_PLUGIN_DIR" >&2
  exit 1
fi
if ! command -v "$BCFTOOLS_BIN" >/dev/null 2>&1; then
  echo "bcftools not found: $BCFTOOLS_BIN" >&2
  exit 1
fi

INPUT_SLICE_GZ="${OUT_PREFIX}.input.vcf.gz"
DUCK_TSV="${OUT_PREFIX}.duckhts.tsv"
BCF_TSV="${OUT_PREFIX}.bcftools.tsv"
COMPARE_TSV="${OUT_PREFIX}.compare.tsv"
DUCK_REJECT_TSV="${OUT_PREFIX}.duckhts.rejects.tsv"
BCF_REJECT_TSV="${OUT_PREFIX}.bcftools.rejects.tsv"
REJECT_COMPARE_TSV="${OUT_PREFIX}.reject_compare.tsv"
BCF_OUT="${OUT_PREFIX}.bcftools.lifted.bcf"
BCF_REJECT="${OUT_PREFIX}.bcftools.reject.vcf"

cleanup() {
  if [[ "$KEEP_SLICE" != "1" ]]; then
    rm -f "$INPUT_SLICE_GZ" "${INPUT_SLICE_GZ}.tbi" "${INPUT_SLICE_GZ}.csi"
  fi
}
trap cleanup EXIT

VIEW_ARGS=()
if [[ -n "$LIFTOVER_REGION" ]]; then
  VIEW_ARGS+=( -r "$LIFTOVER_REGION" )
fi

echo "[1/5] Materialize comparison slice -> $INPUT_SLICE_GZ"
"$BCFTOOLS_BIN" view "${VIEW_ARGS[@]}" -Oz -o "$INPUT_SLICE_GZ" "$INPUT_VCF"
"$BCFTOOLS_BIN" index -f -t "$INPUT_SLICE_GZ"

echo "[2/5] DuckHTS duckdb_liftover -> $DUCK_TSV"
duckdb -unsigned <<SQL
LOAD '${DUCKHTS_EXT}';

CREATE TEMP TABLE src AS
SELECT CHROM AS chrom,
       POS AS pos,
       REF AS ref,
       array_to_string(ALT, ',') AS alt
FROM read_bcf('${INPUT_SLICE_GZ}');

CREATE TEMP TABLE lifted AS
SELECT *
FROM duckdb_liftover(
  'src', 'chrom', 'pos',
  ref_col := 'ref',
  alt_col := 'alt',
  chain_path := '${CHAIN_PATH}',
  dst_fasta_ref := '${DST_FASTA}',
  src_fasta_ref := '${SRC_FASTA}'
);

COPY (
  SELECT dest_chrom, dest_pos, dest_ref, dest_alt, COUNT(*) AS duck_n
  FROM lifted
  WHERE mapped
  GROUP BY ALL
  ORDER BY dest_chrom, dest_pos, dest_ref, dest_alt
) TO '${DUCK_TSV}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);

COPY (
  SELECT reject_reason, COUNT(*) AS duck_n
  FROM lifted
  WHERE NOT mapped
  GROUP BY ALL
  ORDER BY reject_reason
) TO '${DUCK_REJECT_TSV}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
SQL

echo "[3/5] bcftools +liftover -> $BCF_OUT"
BCFTOOLS_PLUGINS="$BCFTOOLS_PLUGIN_DIR" \
  "$BCFTOOLS_BIN" +liftover --no-version -Ob -o "$BCF_OUT" "$INPUT_SLICE_GZ" -- \
  -s "$SRC_FASTA" \
  -f "$DST_FASTA" \
  -c "$CHAIN_PATH" \
  --reject "$BCF_REJECT" \
  --reject-type v \
  --write-src \
  --write-reject

{
  printf 'dest_chrom\tdest_pos\tdest_ref\tdest_alt\tbcf_n\n'
  "$BCFTOOLS_BIN" query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$BCF_OUT" \
    | sort | uniq -c \
    | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$5,$1}'
} > "$BCF_TSV"

{
  printf 'reject_reason\tbcf_n\n'
  awk '!/^#/ {print $7}' "$BCF_REJECT" \
    | sort | uniq -c \
    | awk 'BEGIN{OFS="\t"} {print $2,$1}'
} > "$BCF_REJECT_TSV"

echo "[4/5] Compare mapped outputs -> $COMPARE_TSV"
duckdb -unsigned <<SQL
CREATE TEMP TABLE duck AS
  SELECT * FROM read_csv('${DUCK_TSV}', delim := '\t', header := true,
    columns = {'dest_chrom':'VARCHAR','dest_pos':'BIGINT','dest_ref':'VARCHAR','dest_alt':'VARCHAR','duck_n':'BIGINT'});
CREATE TEMP TABLE bcf AS
  SELECT * FROM read_csv('${BCF_TSV}', delim := '\t', header := true,
    columns = {'dest_chrom':'VARCHAR','dest_pos':'BIGINT','dest_ref':'VARCHAR','dest_alt':'VARCHAR','bcf_n':'BIGINT'});

COPY (
  SELECT
    COALESCE(d.dest_chrom, b.dest_chrom) AS dest_chrom,
    COALESCE(d.dest_pos, b.dest_pos) AS dest_pos,
    COALESCE(d.dest_ref, b.dest_ref) AS dest_ref,
    COALESCE(d.dest_alt, b.dest_alt) AS dest_alt,
    d.duck_n,
    b.bcf_n,
    CASE
      WHEN d.dest_chrom IS NULL THEN 'only_bcftools'
      WHEN b.dest_chrom IS NULL THEN 'only_duckhts'
      WHEN d.duck_n <> b.bcf_n THEN 'diff_counts'
      ELSE 'match'
    END AS status
  FROM duck d
  FULL OUTER JOIN bcf b USING (dest_chrom, dest_pos, dest_ref, dest_alt)
  ORDER BY dest_chrom, dest_pos, dest_ref, dest_alt
) TO '${COMPARE_TSV}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
SQL

echo "[5/5] Compare reject outputs -> $REJECT_COMPARE_TSV"
duckdb -unsigned <<SQL
CREATE TEMP TABLE duck AS
  SELECT * FROM read_csv('${DUCK_REJECT_TSV}', delim := '\t', header := true,
    columns = {'reject_reason':'VARCHAR','duck_n':'BIGINT'});
CREATE TEMP TABLE bcf AS
  SELECT * FROM read_csv('${BCF_REJECT_TSV}', delim := '\t', header := true,
    columns = {'reject_reason':'VARCHAR','bcf_n':'BIGINT'});

COPY (
  SELECT
    COALESCE(d.reject_reason, b.reject_reason) AS reject_reason,
    d.duck_n,
    b.bcf_n,
    CASE
      WHEN d.reject_reason IS NULL THEN 'only_bcftools'
      WHEN b.reject_reason IS NULL THEN 'only_duckhts'
      WHEN d.duck_n <> b.bcf_n THEN 'diff_counts'
      ELSE 'match'
    END AS status
  FROM duck d
  FULL OUTER JOIN bcf b USING (reject_reason)
  ORDER BY reject_reason
) TO '${REJECT_COMPARE_TSV}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
SQL

echo

echo "Mapped output status summary:"
duckdb -unsigned -c "SELECT status, COUNT(*) AS n FROM read_csv_auto('${COMPARE_TSV}', delim := '\t', header := true) GROUP BY status ORDER BY status;"

echo
echo "Reject output status summary:"
duckdb -unsigned -c "SELECT status, COUNT(*) AS n FROM read_csv_auto('${REJECT_COMPARE_TSV}', delim := '\t', header := true) GROUP BY status ORDER BY status;"

echo
echo "Mapped mismatch examples:"
duckdb -unsigned -c "SELECT * FROM read_csv_auto('${COMPARE_TSV}', delim := '\t', header := true) WHERE status <> 'match' LIMIT 20;"

echo
echo "Reject mismatch examples:"
duckdb -unsigned -c "SELECT * FROM read_csv_auto('${REJECT_COMPARE_TSV}', delim := '\t', header := true) WHERE status <> 'match' LIMIT 20;"

echo
echo "Settings:"
echo "  input_vcf:          ${INPUT_VCF}"
echo "  region:             ${LIFTOVER_REGION:-full_input}"
echo "  chain_path:         ${CHAIN_PATH}"
echo "  src_fasta:          ${SRC_FASTA}"
echo "  dst_fasta:          ${DST_FASTA}"
echo "  bcftools_bin:       ${BCFTOOLS_BIN}"
echo "  bcftools_plugin_dir:${BCFTOOLS_PLUGIN_DIR}"
echo ""
echo "Output files:"
echo "  ${DUCK_TSV}"
echo "  ${BCF_TSV}"
echo "  ${COMPARE_TSV}"
echo "  ${DUCK_REJECT_TSV}"
echo "  ${BCF_REJECT_TSV}"
echo "  ${REJECT_COMPARE_TSV}"
