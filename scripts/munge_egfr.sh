#!/usr/bin/env bash
set -euo pipefail

IN_GZ="${1:-hg38_ma_egfr_7_files_maf0.01_rsid.txt.gz}"
HG38_FA="${2:-$HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna}"
OUT_PREFIX="${3:-munge_egfr}"
DUCKHTS_EXT="${4:-build/release/duckhts.duckdb_extension}"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
SCORE_PLUGIN_DIR="${SCORE_PLUGIN_DIR:-}"
COMPARE_FIELDS="${COMPARE_FIELDS:-core}"
COLUMNS_FILE="${COLUMNS_FILE:-scripts/munge_egfr.columns.tsv}"
NUM_TOL="${NUM_TOL:-1e-6}"
REL_TOL="${REL_TOL:-1e-7}"

HG38_FA="$(realpath "$HG38_FA")"

DUCK_TSV="${OUT_PREFIX}.duckhts.tsv"
BCF_PATH="${OUT_PREFIX}.bcftools.bcf"
BCF_NORM_TSV="${OUT_PREFIX}.bcftools.norm.tsv"
DUCK_NORM_TSV="${OUT_PREFIX}.duckhts.norm.tsv"
COMPARE_TSV="${OUT_PREFIX}.compare.tsv"

if [[ ! -f "$IN_GZ" ]]; then
  echo "Input file not found: $IN_GZ" >&2
  exit 1
fi
if [[ ! -f "$HG38_FA" ]]; then
  echo "Reference FASTA not found: $HG38_FA" >&2
  exit 1
fi
if [[ ! -f "$DUCKHTS_EXT" ]]; then
  echo "DuckHTS extension not found: $DUCKHTS_EXT" >&2
  exit 1
fi
if [[ ! -f "$COLUMNS_FILE" ]]; then
  echo "columns file not found: $COLUMNS_FILE" >&2
  exit 1
fi
if ! command -v "$BCFTOOLS_BIN" >/dev/null 2>&1; then
  echo "bcftools not found: $BCFTOOLS_BIN" >&2
  exit 1
fi

if [[ -z "$SCORE_PLUGIN_DIR" ]]; then
  if [[ -f "score_1.22-20250819.zip" ]]; then
    mkdir -p .tmp/score-plugin
    if [[ ! -f .tmp/score-plugin/munge.so ]]; then
      unzip -oq score_1.22-20250819.zip -d .tmp/score-plugin
    fi
    SCORE_PLUGIN_DIR="$(realpath .tmp/score-plugin)"
  else
    echo "Set SCORE_PLUGIN_DIR or provide score_1.22-20250819.zip" >&2
    exit 1
  fi
fi

if [[ ! -f "$SCORE_PLUGIN_DIR/munge.so" ]]; then
  echo "munge.so not found under SCORE_PLUGIN_DIR: $SCORE_PLUGIN_DIR" >&2
  exit 1
fi

echo "[1/4] Run DuckHTS munge -> $DUCK_TSV"
duckdb -unsigned <<SQL
LOAD '${DUCKHTS_EXT}';

SELECT * FROM fasta_index('${HG38_FA}', index_path := '${HG38_FA}.fai');

CREATE OR REPLACE TEMP TABLE munge_input AS
SELECT
  CAST(snpid AS VARCHAR) AS snpid,
  CASE
    WHEN left(CAST(CHR AS VARCHAR), 3) = 'chr' THEN CAST(CHR AS VARCHAR)
    ELSE 'chr' || CAST(CHR AS VARCHAR)
  END AS CHR,
  CAST(BP AS BIGINT) AS BP,
  CAST(A1 AS VARCHAR) AS A1,
  CAST(A2 AS VARCHAR) AS A2,
  CAST(P AS DOUBLE) AS P,
  CAST(b AS DOUBLE) AS BETA,
  CAST(Freq1 AS DOUBLE) AS Freq1,
  CAST(se AS DOUBLE) AS SE,
  CAST(N AS DOUBLE) AS N
FROM read_csv_auto('${IN_GZ}', delim := '\t', header := true);

COPY (
  SELECT
    chrom,
    pos,
    id,
    ref,
    alt,
    alleles_swapped,
    COALESCE(filter, 'PASS') AS filter,
    ns,
    es,
    se,
    lp,
    af,
    ac
  FROM duckdb_munge(
    'munge_input',
    column_map := map(
      ['SNP','CHR','BP','A1','A2','P','BETA','FRQ','SE','N'],
      ['snpid','CHR','BP','A1','A2','P','BETA','Freq1','SE','N']
    ),
    fasta_ref := '${HG38_FA}',
    iffy_tag := 'IFFY',
    mismatch_tag := 'REF_MISMATCH'
  )
)
TO '${DUCK_TSV}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
SQL

echo "[2/4] Run bcftools +munge -> $BCF_PATH"
BCFTOOLS_PLUGINS="$SCORE_PLUGIN_DIR" \
  "$BCFTOOLS_BIN" +munge --no-version -C "$COLUMNS_FILE" -f "$HG38_FA" -Ob -o "$BCF_PATH" "$IN_GZ"

echo "[3/4] Export normalized debug TSVs"
duckdb -unsigned <<SQL
LOAD '${DUCKHTS_EXT}';

COPY (
  SELECT
    CHROM AS chrom,
    POS AS pos,
    COALESCE(ID, CHROM || ':' || CAST(POS AS VARCHAR)) AS id,
    REF AS ref,
    ALT[1] AS alt,
    CASE
      WHEN FILTER IS NULL OR list_count(FILTER) = 0 THEN 'PASS'
      WHEN list_count(FILTER) = 1 THEN FILTER[1]
      ELSE array_to_string(FILTER, ';')
    END AS filter,
    CAST(NULL AS DOUBLE) AS ns,
    CAST(NULL AS DOUBLE) AS es,
    FORMAT_SE_SAMPLE[1]::DOUBLE AS se,
    FORMAT_LP_SAMPLE[1]::DOUBLE AS lp,
    CAST(NULL AS DOUBLE) AS af,
    CAST(NULL AS DOUBLE) AS ac
  FROM read_bcf('${BCF_PATH}')
)
TO '${BCF_NORM_TSV}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);

COPY (
  SELECT chrom, pos, id, ref, alt, filter, ns, es, se, lp, af, ac
  FROM read_csv_auto('${DUCK_TSV}', delim := '\t', header := true)
)
TO '${DUCK_NORM_TSV}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
SQL

echo "[4/4] Compare normalized outputs -> $COMPARE_TSV (fields: $COMPARE_FIELDS)"
if [[ "${COMPARE_FIELDS,,}" == "all" ]]; then
  COMPARE_PREDICATE="
      coalesce(d.filter, '') = coalesce(b.filter, '')
      AND (d.ns IS NULL AND b.ns IS NULL OR d.ns IS NOT NULL AND b.ns IS NOT NULL AND (abs(d.ns - b.ns) <= ${NUM_TOL} OR abs(d.ns - b.ns) <= ${REL_TOL} * greatest(abs(d.ns), abs(b.ns), 1.0)))
      AND (d.es IS NULL AND b.es IS NULL OR d.es IS NOT NULL AND b.es IS NOT NULL AND (abs(d.es - b.es) <= ${NUM_TOL} OR abs(d.es - b.es) <= ${REL_TOL} * greatest(abs(d.es), abs(b.es), 1.0)))
      AND (d.se IS NULL AND b.se IS NULL OR d.se IS NOT NULL AND b.se IS NOT NULL AND (abs(d.se - b.se) <= ${NUM_TOL} OR abs(d.se - b.se) <= ${REL_TOL} * greatest(abs(d.se), abs(b.se), 1.0)))
      AND (d.lp IS NULL AND b.lp IS NULL OR d.lp IS NOT NULL AND b.lp IS NOT NULL AND (abs(d.lp - b.lp) <= ${NUM_TOL} OR abs(d.lp - b.lp) <= ${REL_TOL} * greatest(abs(d.lp), abs(b.lp), 1.0)))
      AND (d.af IS NULL AND b.af IS NULL OR d.af IS NOT NULL AND b.af IS NOT NULL AND (abs(d.af - b.af) <= ${NUM_TOL} OR abs(d.af - b.af) <= ${REL_TOL} * greatest(abs(d.af), abs(b.af), 1.0)))
      AND (d.ac IS NULL AND b.ac IS NULL OR d.ac IS NOT NULL AND b.ac IS NOT NULL AND (abs(d.ac - b.ac) <= ${NUM_TOL} OR abs(d.ac - b.ac) <= ${REL_TOL} * greatest(abs(d.ac), abs(b.ac), 1.0)))"
  COMPARE_FIELD_LIST="filter,ns,es,se,lp,af,ac"
else
  COMPARE_PREDICATE="
      coalesce(d.filter, '') = coalesce(b.filter, '')
      AND (d.se IS NULL AND b.se IS NULL OR d.se IS NOT NULL AND b.se IS NOT NULL AND (abs(d.se - b.se) <= ${NUM_TOL} OR abs(d.se - b.se) <= ${REL_TOL} * greatest(abs(d.se), abs(b.se), 1.0)))
      AND (d.lp IS NULL AND b.lp IS NULL OR d.lp IS NOT NULL AND b.lp IS NOT NULL AND (abs(d.lp - b.lp) <= ${NUM_TOL} OR abs(d.lp - b.lp) <= ${REL_TOL} * greatest(abs(d.lp), abs(b.lp), 1.0)))"
  COMPARE_FIELD_LIST="filter,se,lp"
fi

duckdb -unsigned <<SQL
LOAD '${DUCKHTS_EXT}';

COPY (
  WITH duck_norm AS (
    SELECT chrom, pos, id, ref, alt, filter, ns, es, se, lp, af, ac
    FROM read_csv_auto('${DUCK_TSV}', delim := '\t', header := true)
  ),
  bcf_norm AS (
    SELECT
      CHROM AS chrom,
      POS AS pos,
      COALESCE(ID, CHROM || ':' || CAST(POS AS VARCHAR)) AS id,
      REF AS ref,
      ALT[1] AS alt,
      CASE
        WHEN FILTER IS NULL OR list_count(FILTER) = 0 THEN 'PASS'
        WHEN list_count(FILTER) = 1 THEN FILTER[1]
        ELSE array_to_string(FILTER, ';')
      END AS filter,
      CAST(NULL AS DOUBLE) AS ns,
      CAST(NULL AS DOUBLE) AS es,
      FORMAT_SE_SAMPLE[1]::DOUBLE AS se,
      FORMAT_LP_SAMPLE[1]::DOUBLE AS lp,
      CAST(NULL AS DOUBLE) AS af,
      CAST(NULL AS DOUBLE) AS ac
    FROM read_bcf('${BCF_PATH}')
  ),
  j AS (
    SELECT
      coalesce(d.chrom, b.chrom) AS chrom,
      coalesce(d.pos, b.pos) AS pos,
      d.id AS duck_id,
      b.id AS bcf_id,
      coalesce(d.ref, b.ref) AS ref,
      coalesce(d.alt, b.alt) AS alt,
      d.filter AS duck_filter,
      b.filter AS bcf_filter,
      d.ns AS duck_ns,
      b.ns AS bcf_ns,
      d.es AS duck_es,
      b.es AS bcf_es,
      d.se AS duck_se,
      b.se AS bcf_se,
      d.lp AS duck_lp,
      b.lp AS bcf_lp,
      d.af AS duck_af,
      b.af AS bcf_af,
      d.ac AS duck_ac,
      b.ac AS bcf_ac,
      CASE
        WHEN d.chrom IS NULL THEN 'only_bcftools'
        WHEN b.chrom IS NULL THEN 'only_duckhts'
        WHEN ${COMPARE_PREDICATE} THEN 'match'
        ELSE 'diff_values'
      END AS status
    FROM duck_norm d
    FULL OUTER JOIN bcf_norm b
      ON d.chrom = b.chrom
     AND d.pos = b.pos
     AND d.ref = b.ref
     AND d.alt = b.alt
  )
  SELECT *
  FROM j
  ORDER BY chrom, pos, ref, alt
)
TO '${COMPARE_TSV}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);

SELECT status, COUNT(*) AS n
FROM read_csv_auto('${COMPARE_TSV}', delim := '\t', header := true)
GROUP BY status
ORDER BY status;
SQL

echo "compare_fields: ${COMPARE_FIELD_LIST}"
echo "numeric_tolerance: ${NUM_TOL}"
echo "relative_tolerance: ${REL_TOL}"
echo "compare file: ${COMPARE_TSV}"

echo "Done. Files:"
echo "  $DUCK_TSV"
echo "  $BCF_PATH"
echo "  $BCF_NORM_TSV"
echo "  $DUCK_NORM_TSV"
echo "  $COMPARE_TSV"
