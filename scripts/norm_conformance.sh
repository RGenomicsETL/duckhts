#!/usr/bin/env bash
set -euo pipefail

# Compare DuckHTS duckhts_bcftools_norm(...) against installed bcftools norm.
#
# Usage:
#   bash scripts/norm_conformance.sh INPUT_VCF FASTA [OUT_PREFIX]
#
# Environment variables:
#   DUCKHTS_EXT   path to duckhts.duckdb_extension
#   BCFTOOLS_BIN  path to bcftools binary (defaults to RBCFTools::bcftools_path())
#   NORM_REGION   optional region string to slice before comparison
#   NORM_SPLIT    1 => compare split multiallelic normalization (`bcftools norm -m -any`)
#                 0 => compare site-preserving normalization (default)
#   KEEP_SLICE    1 to keep the materialized local VCF slice

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: bash scripts/norm_conformance.sh INPUT_VCF FASTA [OUT_PREFIX]" >&2
  exit 1
fi

INPUT_VCF="$1"
FASTA="$2"
OUT_PREFIX="${3:-norm_conformance}"
DUCKHTS_EXT="${DUCKHTS_EXT:-build/release/duckhts.duckdb_extension}"
NORM_REGION="${NORM_REGION:-}"
NORM_SPLIT="${NORM_SPLIT:-0}"
KEEP_SLICE="${KEEP_SLICE:-0}"

if [[ -z "${BCFTOOLS_BIN:-}" ]]; then
  BCFTOOLS_BIN="$(Rscript -e "if (!requireNamespace('RBCFTools', quietly=TRUE)) quit(status=1); cat(RBCFTools::bcftools_path())")" || {
    echo "RBCFTools is required to resolve the bundled bcftools binary; set BCFTOOLS_BIN to override" >&2
    exit 1
  }
fi

if [[ ! -f "$FASTA" ]]; then
  echo "FASTA not found: $FASTA" >&2
  exit 1
fi
if [[ ! -f "$DUCKHTS_EXT" ]]; then
  echo "DuckHTS extension not found: $DUCKHTS_EXT" >&2
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
DUCK_STATUS_TSV="${OUT_PREFIX}.duckhts.status.tsv"
BCF_OUT="${OUT_PREFIX}.bcftools.norm.vcf.gz"

cleanup() {
  if [[ "$KEEP_SLICE" != "1" ]]; then
    rm -f "$INPUT_SLICE_GZ" "${INPUT_SLICE_GZ}.tbi" "${INPUT_SLICE_GZ}.csi"
  fi
}
trap cleanup EXIT

VIEW_ARGS=()
if [[ -n "$NORM_REGION" ]]; then
  VIEW_ARGS+=( -r "$NORM_REGION" )
fi

NORM_ARGS=(-f "$FASTA")
DUCK_SPLIT_SQL="false"
if [[ "$NORM_SPLIT" == "1" ]]; then
  NORM_ARGS=(-m -any -f "$FASTA")
  DUCK_SPLIT_SQL="true"
fi

echo "[1/4] Materialize comparison slice -> $INPUT_SLICE_GZ"
"$BCFTOOLS_BIN" view "${VIEW_ARGS[@]}" -Oz -o "$INPUT_SLICE_GZ" "$INPUT_VCF"
"$BCFTOOLS_BIN" index -f -t "$INPUT_SLICE_GZ"

echo "[2/4] DuckHTS duckhts_bcftools_norm -> $DUCK_TSV"
duckdb -unsigned <<SQL
LOAD '${DUCKHTS_EXT}';

CREATE TEMP TABLE src AS
SELECT CHROM AS chrom,
       POS AS pos,
       REF AS ref,
       ALT
FROM read_bcf('${INPUT_SLICE_GZ}');

CREATE TEMP TABLE normed AS
SELECT *
FROM duckhts_bcftools_norm(
  'src',
  '${FASTA}',
  split_multiallelic := ${DUCK_SPLIT_SQL}
);

COPY (
  SELECT chrom,
         pos_normed,
         end_pos_normed,
         ref_normed,
         CASE
           WHEN typeof(alt_normed) = 'VARCHAR[]'
             THEN replace(replace(replace(alt_normed::VARCHAR, '[', ''), ']', ''), ', ', ',')
           ELSE alt_normed::VARCHAR
         END AS alt_normed,
         COUNT(*) AS duck_n
  FROM normed
  WHERE coalesce(norm_status, '') IN ('Normalized', 'Unchanged', 'RefOnly', 'SpanningDeletion', 'Breakend')
  GROUP BY ALL
  ORDER BY chrom, pos_normed, ref_normed, alt_normed
) TO '${DUCK_TSV}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);

COPY (
  SELECT norm_status, COUNT(*) AS duck_n
  FROM normed
  GROUP BY ALL
  ORDER BY norm_status
) TO '${DUCK_STATUS_TSV}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
SQL

echo "[3/4] bcftools norm -> $BCF_OUT"
"$BCFTOOLS_BIN" norm "${NORM_ARGS[@]}" -Oz -o "$BCF_OUT" "$INPUT_SLICE_GZ"
"$BCFTOOLS_BIN" index -f -t "$BCF_OUT"

HAS_END_TAG=0
if "$BCFTOOLS_BIN" view -h "$BCF_OUT" | grep -q '^##INFO=<ID=END,'; then
  HAS_END_TAG=1
fi

if [[ "$HAS_END_TAG" == "1" ]]; then
  {
    printf 'chrom\tpos_normed\tend_pos_normed\tref_normed\talt_normed\tbcf_n\n'
    "$BCFTOOLS_BIN" query -f '%CHROM\t%POS\t%INFO/END\t%REF\t%ALT\n' "$BCF_OUT" \
      | awk 'BEGIN{OFS="\t"} { if ($3=="." || $3=="") $3=$2+length($4)-1; print $1,$2,$3,$4,$5 }' \
      | sort | uniq -c \
      | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$5,$6,$1}'
  } > "$BCF_TSV"
else
  {
    printf 'chrom\tpos_normed\tend_pos_normed\tref_normed\talt_normed\tbcf_n\n'
    "$BCFTOOLS_BIN" query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$BCF_OUT" \
      | awk 'BEGIN{OFS="\t"} {print $1,$2,$2+length($3)-1,$3,$4}' \
      | sort | uniq -c \
      | awk 'BEGIN{OFS="\t"} {print $2,$3,$4,$5,$6,$1}'
  } > "$BCF_TSV"
fi

echo "[4/4] Compare outputs -> $COMPARE_TSV"
duckdb -unsigned <<SQL
CREATE TEMP TABLE duck AS
  SELECT * FROM read_csv('${DUCK_TSV}', delim := '\t', header := true,
    columns = {'chrom':'VARCHAR','pos_normed':'BIGINT','end_pos_normed':'BIGINT','ref_normed':'VARCHAR','alt_normed':'VARCHAR','duck_n':'BIGINT'});
CREATE TEMP TABLE bcf AS
  SELECT * FROM read_csv('${BCF_TSV}', delim := '\t', header := true,
    columns = {'chrom':'VARCHAR','pos_normed':'BIGINT','end_pos_normed':'BIGINT','ref_normed':'VARCHAR','alt_normed':'VARCHAR','bcf_n':'BIGINT'});

COPY (
  SELECT
    COALESCE(d.chrom, b.chrom) AS chrom,
    COALESCE(d.pos_normed, b.pos_normed) AS pos_normed,
    COALESCE(d.end_pos_normed, b.end_pos_normed) AS end_pos_normed,
    COALESCE(d.ref_normed, b.ref_normed) AS ref_normed,
    COALESCE(d.alt_normed, b.alt_normed) AS alt_normed,
    d.duck_n,
    b.bcf_n,
    CASE
      WHEN d.chrom IS NULL THEN 'only_bcftools'
      WHEN b.chrom IS NULL THEN 'only_duckhts'
      WHEN d.duck_n <> b.bcf_n THEN 'diff_counts'
      ELSE 'match'
    END AS status
  FROM duck d
  FULL OUTER JOIN bcf b USING (chrom, pos_normed, end_pos_normed, ref_normed, alt_normed)
  ORDER BY chrom, pos_normed, end_pos_normed, ref_normed, alt_normed
) TO '${COMPARE_TSV}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
SQL

echo
echo "DuckHTS norm_status summary:"
duckdb -unsigned -c "SELECT * FROM read_csv_auto('${DUCK_STATUS_TSV}', delim := '\t', header := true);"

echo
echo "Comparison status summary:"
duckdb -unsigned -c "SELECT status, COUNT(*) AS n FROM read_csv_auto('${COMPARE_TSV}', delim := '\t', header := true) GROUP BY status ORDER BY status;"

echo
echo "Mismatch examples:"
duckdb -unsigned -c "SELECT * FROM read_csv_auto('${COMPARE_TSV}', delim := '\t', header := true) WHERE status <> 'match' LIMIT 20;"

echo
echo "Settings:"
echo "  input_vcf:    ${INPUT_VCF}"
echo "  fasta:        ${FASTA}"
echo "  region:       ${NORM_REGION:-full_input}"
echo "  split:        ${NORM_SPLIT}"
echo "  bcftools_bin: ${BCFTOOLS_BIN}"
echo ""
echo "Output files:"
echo "  ${DUCK_TSV}"
echo "  ${BCF_TSV}"
echo "  ${COMPARE_TSV}"
echo "  ${DUCK_STATUS_TSV}"
