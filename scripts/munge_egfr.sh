#!/usr/bin/env bash
set -euo pipefail

IN_GZ="${1:-hg38_ma_egfr_7_files_maf0.01_rsid.txt.gz}"
HG38_FA="${2:-$HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna}"
OUT_PREFIX="${3:-munge_egfr}"
DUCKHTS_EXT="${4:-build/release/duckhts.duckdb_extension}"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
SCORE_PLUGIN_DIR="${SCORE_PLUGIN_DIR:-}"

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
if ! command -v "$BCFTOOLS_BIN" >/dev/null 2>&1; then
  echo "bcftools not found: $BCFTOOLS_BIN" >&2
  exit 1
fi

if [[ -z "$SCORE_PLUGIN_DIR" ]]; then
  if [[ -f "score_1.22-20250819.zip" ]]; then
    mkdir -p .tmp/score-plugin
    if [[ ! -f .tmp/score-plugin/plugins/munge.so ]]; then
      unzip -oq score_1.22-20250819.zip -d .tmp/score-plugin
    fi
    SCORE_PLUGIN_DIR="$(realpath .tmp/score-plugin/plugins)"
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
  "$BCFTOOLS_BIN" +munge --no-version -c PLINK -f "$HG38_FA" -Ob -o "$BCF_PATH" "$IN_GZ"

echo "[3/4] Read bcftools BCF with DuckHTS and normalize"
duckdb -unsigned <<SQL
LOAD '${DUCKHTS_EXT}';

COPY (
  SELECT
    CHROM AS chrom,
    POS AS pos,
    ID AS id,
    REF AS ref,
    ALT[1] AS alt,
    COALESCE(NULLIF(FILTER, ''), 'PASS') AS filter,
    INFO_NS AS ns,
    INFO_ES AS es,
    INFO_SE AS se,
    INFO_LP AS lp,
    INFO_AF AS af,
    INFO_AC AS ac
  FROM read_bcf('${BCF_PATH}')
)
TO '${BCF_NORM_TSV}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);

COPY (
  SELECT chrom, pos, id, ref, alt, filter, ns, es, se, lp, af, ac
  FROM read_csv_auto('${DUCK_TSV}', delim := '\t', header := true)
)
TO '${DUCK_NORM_TSV}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
SQL

echo "[4/4] Compare normalized outputs -> $COMPARE_TSV"
DUCK_NORM_TSV="$DUCK_NORM_TSV" BCF_NORM_TSV="$BCF_NORM_TSV" COMPARE_TSV="$COMPARE_TSV" python3 - <<'PY'
import csv
import os
from pathlib import Path

duck = Path(os.environ['DUCK_NORM_TSV'])
bcf = Path(os.environ['BCF_NORM_TSV'])
out = Path(os.environ['COMPARE_TSV'])

def load(path):
    rows = {}
    with path.open() as fh:
        r = csv.DictReader(fh, delimiter='\t')
        for row in r:
            k = (row['chrom'], row['pos'], row['id'], row['ref'], row['alt'])
            rows[k] = row
    return rows

d = load(duck)
b = load(bcf)
keys = sorted(set(d) | set(b))

fields = ['chrom','pos','id','ref','alt','status','duck_filter','bcf_filter','duck_ns','bcf_ns','duck_es','bcf_es','duck_se','bcf_se','duck_lp','bcf_lp','duck_af','bcf_af','duck_ac','bcf_ac']
with out.open('w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter='\t')
    w.writeheader()
    for k in keys:
        dr = d.get(k)
        br = b.get(k)
        status = 'match'
        if dr is None:
            status = 'only_bcftools'
        elif br is None:
            status = 'only_duckhts'
        else:
            for f in ('filter','ns','es','se','lp','af','ac'):
                if (dr.get(f,'') or '') != (br.get(f,'') or ''):
                    status = 'diff_values'
                    break
        w.writerow({
            'chrom': k[0], 'pos': k[1], 'id': k[2], 'ref': k[3], 'alt': k[4], 'status': status,
            'duck_filter': '' if dr is None else dr.get('filter',''),
            'bcf_filter': '' if br is None else br.get('filter',''),
            'duck_ns': '' if dr is None else dr.get('ns',''),
            'bcf_ns': '' if br is None else br.get('ns',''),
            'duck_es': '' if dr is None else dr.get('es',''),
            'bcf_es': '' if br is None else br.get('es',''),
            'duck_se': '' if dr is None else dr.get('se',''),
            'bcf_se': '' if br is None else br.get('se',''),
            'duck_lp': '' if dr is None else dr.get('lp',''),
            'bcf_lp': '' if br is None else br.get('lp',''),
            'duck_af': '' if dr is None else dr.get('af',''),
            'bcf_af': '' if br is None else br.get('af',''),
            'duck_ac': '' if dr is None else dr.get('ac',''),
            'bcf_ac': '' if br is None else br.get('ac',''),
        })

from collections import Counter
c = Counter()
with out.open() as fh:
    r = csv.DictReader(fh, delimiter='\t')
    for row in r:
        c[row['status']] += 1

print('comparison summary:', dict(c))
print('compare file:', out)
PY

echo "Done. Files:"
echo "  $DUCK_TSV"
echo "  $BCF_PATH"
echo "  $BCF_NORM_TSV"
echo "  $DUCK_NORM_TSV"
echo "  $COMPARE_TSV"
