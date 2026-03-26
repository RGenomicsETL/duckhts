#!/usr/bin/env bash
set -euo pipefail

IN_GZ="${1:-hg38_ma_egfr_7_files_maf0.01_rsid.txt.gz}"
HG38_FA="${2:-"$HOME/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"}"
HG38_FA="$(realpath "$HG38_FA")"
OUT_TSV="${3:-munge_out.tsv}"
DUCKHTS_EXT="${4:-build/release/duckhts.duckdb_extension}"

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

time duckdb -unsigned <<SQL
LOAD '${DUCKHTS_EXT}';

SELECT *
FROM fasta_index('${HG38_FA}', index_path := '${HG38_FA}.fai');

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
  SELECT *
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
TO '${OUT_TSV}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
SQL

echo "Wrote: ${OUT_TSV}"
