#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VCF_PATH="${VARIANTKEY_VCF:-}"
REGION="${VARIANTKEY_REGION:-}"
BCFTOOLS_BIN="${VARIANTKEY_BCFTOOLS_BIN:-bcftools}"
EXT_PATH="${DUCKHTS_EXTENSION:-$ROOT_DIR/build/release/duckhts.duckdb_extension}"
OUT_DIR="${VARIANTKEY_OUT_DIR:-$(mktemp -d /tmp/duckhts_variantkey_conf.XXXXXX)}"
LABEL="${VARIANTKEY_LABEL:-$(basename "${VCF_PATH:-variantkey}")}" 

if [[ -z "$VCF_PATH" ]]; then
  echo "ERROR: set VARIANTKEY_VCF to a local or remote indexed VCF/BCF path" >&2
  exit 1
fi

if [[ ! -f "$EXT_PATH" ]]; then
  echo "ERROR: DuckHTS extension not found at $EXT_PATH" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"
BCFTOOLS_TSV="$OUT_DIR/${LABEL}.bcftools.tsv"
DUCKHTS_TSV="$OUT_DIR/${LABEL}.duckhts.tsv"
COMPARE_TSV="$OUT_DIR/${LABEL}.compare.tsv"
SUMMARY_TSV="$OUT_DIR/${LABEL}.summary.tsv"

BCFTOOLS_ARGS=(query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%VKX\n')
if [[ -n "$REGION" ]]; then
  BCFTOOLS_ARGS+=(-r "$REGION")
fi
BCFTOOLS_ARGS+=("$VCF_PATH")

"$BCFTOOLS_BIN" "${BCFTOOLS_ARGS[@]}" > "$BCFTOOLS_TSV"

python3 - <<'PY' "$VCF_PATH" "$REGION" "$EXT_PATH" "$DUCKHTS_TSV" "$BCFTOOLS_TSV" "$COMPARE_TSV" "$SUMMARY_TSV"
import csv
import duckdb
import os
import sys

vcf_path, region, ext_path, duckhts_tsv, bcftools_tsv, compare_tsv, summary_tsv = sys.argv[1:8]

con = duckdb.connect(config={"allow_unsigned_extensions": "true"})
con.execute(f"LOAD '{ext_path}'")

def sql_quote(text: str) -> str:
    return "'" + text.replace("'", "''") + "'"

region_sql = ""
if region:
    region_sql = f", region := {sql_quote(region)}"

con.execute(
    f"""
    COPY (
      SELECT CHROM,
             POS,
             REF,
             coalesce(ALT[1], '.') AS ALT,
             variantkey_hex(variantkey(CHROM, POS, REF, coalesce(ALT[1], '.'))) AS VKX
      FROM read_bcf({sql_quote(vcf_path)}{region_sql})
      ORDER BY CHROM, POS, REF, ALT
    ) TO {sql_quote(duckhts_tsv)} (FORMAT CSV, DELIMITER '\t', HEADER TRUE)
    """
)

for path in (bcftools_tsv, duckhts_tsv):
    if os.path.getsize(path) == 0:
        with open(path, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter='\t')
            writer.writerow(["CHROM", "POS", "REF", "ALT", "VKX"])
    else:
        with open(path, "r", newline="") as fh:
            first = fh.readline()
        if not first.startswith("CHROM\tPOS\tREF\tALT\tVKX"):
            original = open(path, "r", newline="").read()
            with open(path, "w", newline="") as fh:
                fh.write("CHROM\tPOS\tREF\tALT\tVKX\n")
                fh.write(original)

compare_sql = f"""
COPY (
  WITH bcftools_rows AS (
    SELECT CHROM, POS, REF, ALT, VKX, COUNT(*) AS n
    FROM read_csv({sql_quote(bcftools_tsv)}, header := TRUE, delim := '\t',
                  columns := {{'CHROM':'VARCHAR','POS':'BIGINT','REF':'VARCHAR','ALT':'VARCHAR','VKX':'VARCHAR'}})
    GROUP BY ALL
  ),
  duckhts_rows AS (
    SELECT CHROM, POS, REF, ALT, VKX, COUNT(*) AS n
    FROM read_csv({sql_quote(duckhts_tsv)}, header := TRUE, delim := '\t',
                  columns := {{'CHROM':'VARCHAR','POS':'BIGINT','REF':'VARCHAR','ALT':'VARCHAR','VKX':'VARCHAR'}})
    GROUP BY ALL
  )
  SELECT coalesce(b.CHROM, d.CHROM) AS CHROM,
         coalesce(b.POS, d.POS) AS POS,
         coalesce(b.REF, d.REF) AS REF,
         coalesce(b.ALT, d.ALT) AS ALT,
         coalesce(b.VKX, d.VKX) AS VKX,
         coalesce(b.n, 0) AS bcftools_n,
         coalesce(d.n, 0) AS duckhts_n
  FROM bcftools_rows b
  FULL OUTER JOIN duckhts_rows d USING (CHROM, POS, REF, ALT, VKX)
  WHERE coalesce(b.n, 0) <> coalesce(d.n, 0)
  ORDER BY 1, 2, 3, 4, 5
) TO {sql_quote(compare_tsv)} (FORMAT CSV, DELIMITER '\t', HEADER TRUE)
"""
con.execute(compare_sql)

summary = con.execute(
    f"""
    WITH bcftools_rows AS (
      SELECT COUNT(*) AS n
      FROM read_csv({sql_quote(bcftools_tsv)}, header := TRUE, delim := '\t',
                    columns := {{'CHROM':'VARCHAR','POS':'BIGINT','REF':'VARCHAR','ALT':'VARCHAR','VKX':'VARCHAR'}})
    ),
    duckhts_rows AS (
      SELECT COUNT(*) AS n
      FROM read_csv({sql_quote(duckhts_tsv)}, header := TRUE, delim := '\t',
                    columns := {{'CHROM':'VARCHAR','POS':'BIGINT','REF':'VARCHAR','ALT':'VARCHAR','VKX':'VARCHAR'}})
    ),
    mismatch_rows AS (
      SELECT COUNT(*) AS n
      FROM read_csv({sql_quote(compare_tsv)}, header := TRUE, delim := '\t',
                    columns := {{'CHROM':'VARCHAR','POS':'BIGINT','REF':'VARCHAR','ALT':'VARCHAR','VKX':'VARCHAR','bcftools_n':'BIGINT','duckhts_n':'BIGINT'}})
    )
    SELECT (SELECT n FROM bcftools_rows) AS bcftools_rows,
           (SELECT n FROM duckhts_rows) AS duckhts_rows,
           (SELECT n FROM mismatch_rows) AS mismatch_groups
    """
).fetchone()

with open(summary_tsv, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter='\t')
    writer.writerow(["bcftools_rows", "duckhts_rows", "mismatch_groups"])
    writer.writerow(summary)

print(f"bcftools_rows\tduckhts_rows\tmismatch_groups")
print(f"{summary[0]}\t{summary[1]}\t{summary[2]}")

if summary[2] != 0:
    sys.exit(1)
PY

echo "Wrote:"
echo "  $BCFTOOLS_TSV"
echo "  $DUCKHTS_TSV"
echo "  $COMPARE_TSV"
echo "  $SUMMARY_TSV"
