#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:?output directory required}"
ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
source "${ROOT_DIR}/scripts/conformance/common.sh"

conformance_prepare_dirs "$OUT_DIR"
SUMMARY_TSV="${OUT_DIR}/report/summary.tsv"
conformance_write_summary_header "$SUMMARY_TSV"

REF_FASTA="${ROOT_DIR}/test/data/liftover_src.fa"
TMP_DIR="${OUT_DIR}/tmp"

core_input="${TMP_DIR}/core.tsv"
core_columns="${TMP_DIR}/core.columns.tsv"
cat > "$core_input" <<'EOF'
SNP	CHR	BP	A1	A2	P	BETA	FRQ	SE	N
rs_swap	chrF	2	C	A	0.01	0.25	0.1	0.05	1000
rs_pass	chrF	2	A	C	0.01	0.25	0.1	0.05	1000
rs_iffy	chrF	2	C	C	0.01	0.25	0.1	0.05	1000
rs_mismatch	chrF	2	G	T	0.01	0.25	0.1	0.05	1000
EOF
cat > "$core_columns" <<'EOF'
SNP	SNP
CHR	CHR
BP	BP
A1	A1
A2	A2
P	P
BETA	BETA
FRQ	FRQ
SE	SE
N	N
EOF

duckdb -unsigned >"${OUT_DIR}/raw/core.duck.log" 2>&1 <<SQL
LOAD '${DUCKHTS_EXT}';
CREATE OR REPLACE TEMP TABLE munge_input AS
SELECT * FROM read_csv_auto('${core_input}', delim := '\t', header := true);
COPY (
  SELECT chrom, pos, id, ref, alt, coalesce(filter, 'PASS') AS filter,
         ns, es, se, lp, af
  FROM duckdb_munge(
    'munge_input',
    column_map := map(
      ['SNP','CHR','BP','A1','A2','P','BETA','FRQ','SE','N'],
      ['SNP','CHR','BP','A1','A2','P','BETA','FRQ','SE','N']
    ),
    fasta_ref := '${REF_FASTA}'
  )
  ORDER BY id
)
TO '${OUT_DIR}/normalized/core.duck.tsv' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
SQL

"${BCFTOOLS_BIN}" +munge --no-version -C "$core_columns" -f "$REF_FASTA" -Ov \
  -o "${OUT_DIR}/raw/core.upstream.vcf" "$core_input" \
  >"${OUT_DIR}/raw/core.upstream.log" 2>&1

duckdb -unsigned <<SQL
LOAD '${DUCKHTS_EXT}';
COPY (
  SELECT CHROM AS chrom,
         POS AS pos,
         ID AS id,
         REF AS ref,
         ALT[1] AS alt,
         CASE
           WHEN FILTER IS NULL OR list_count(FILTER) = 0 THEN 'PASS'
           WHEN list_count(FILTER) = 1 THEN FILTER[1]
           ELSE array_to_string(FILTER, ';')
         END AS filter,
         FORMAT_NS_SAMPLE[1]::DOUBLE AS ns,
         FORMAT_ES_SAMPLE[1]::DOUBLE AS es,
         FORMAT_SE_SAMPLE[1]::DOUBLE AS se,
         FORMAT_LP_SAMPLE[1]::DOUBLE AS lp,
         FORMAT_AF_SAMPLE[1]::DOUBLE AS af
  FROM read_bcf('${OUT_DIR}/raw/core.upstream.vcf')
  ORDER BY id
)
TO '${OUT_DIR}/normalized/core.upstream.tsv' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);

COPY (
  WITH duck_norm AS (
    SELECT * FROM read_csv_auto('${OUT_DIR}/normalized/core.duck.tsv', delim := '\t', header := true)
  ),
  upstream_norm AS (
    SELECT * FROM read_csv_auto('${OUT_DIR}/normalized/core.upstream.tsv', delim := '\t', header := true)
  )
  SELECT
    coalesce(d.id, u.id) AS id,
    d.chrom AS duck_chrom,
    u.chrom AS upstream_chrom,
    d.ref AS duck_ref,
    u.ref AS upstream_ref,
    d.alt AS duck_alt,
    u.alt AS upstream_alt,
    d.filter AS duck_filter,
    u.filter AS upstream_filter,
    d.ns AS duck_ns,
    u.ns AS upstream_ns,
    d.es AS duck_es,
    u.es AS upstream_es,
    d.se AS duck_se,
    u.se AS upstream_se,
    d.lp AS duck_lp,
    u.lp AS upstream_lp,
    d.af AS duck_af,
    u.af AS upstream_af,
    CASE
      WHEN d.id IS NULL THEN 'only_bcftools'
      WHEN u.id IS NULL THEN 'only_duckhts'
      WHEN coalesce(d.chrom, '') = coalesce(u.chrom, '')
       AND coalesce(d.ref, '') = coalesce(u.ref, '')
       AND coalesce(d.alt, '') = coalesce(u.alt, '')
       AND coalesce(d.filter, '') = coalesce(u.filter, '')
       AND abs(coalesce(d.ns, 0.0) - coalesce(u.ns, 0.0)) <= 1e-6
       AND abs(coalesce(d.es, 0.0) - coalesce(u.es, 0.0)) <= 1e-6
       AND abs(coalesce(d.se, 0.0) - coalesce(u.se, 0.0)) <= 1e-6
       AND abs(coalesce(d.lp, 0.0) - coalesce(u.lp, 0.0)) <= 1e-6
       AND abs(coalesce(d.af, 0.0) - coalesce(u.af, 0.0)) <= 1e-6
        THEN 'match'
      ELSE 'diff_values'
    END AS status,
    CASE
      WHEN d.id IS NULL OR u.id IS NULL THEN 'surface_difference'
      WHEN abs(coalesce(d.ns, 0.0) - coalesce(u.ns, 0.0)) <= 1e-6
       AND abs(coalesce(d.es, 0.0) - coalesce(u.es, 0.0)) <= 1e-6
       AND abs(coalesce(d.se, 0.0) - coalesce(u.se, 0.0)) <= 1e-6
       AND abs(coalesce(d.lp, 0.0) - coalesce(u.lp, 0.0)) <= 1e-6
       AND abs(coalesce(d.af, 0.0) - coalesce(u.af, 0.0)) <= 1e-6
        THEN 'match'
      ELSE 'unclassified'
    END AS reason
  FROM duck_norm d
  FULL OUTER JOIN upstream_norm u USING (id)
  WHERE NOT (
    d.id IS NOT NULL AND u.id IS NOT NULL
    AND coalesce(d.chrom, '') = coalesce(u.chrom, '')
    AND coalesce(d.ref, '') = coalesce(u.ref, '')
    AND coalesce(d.alt, '') = coalesce(u.alt, '')
    AND coalesce(d.filter, '') = coalesce(u.filter, '')
    AND abs(coalesce(d.ns, 0.0) - coalesce(u.ns, 0.0)) <= 1e-6
    AND abs(coalesce(d.es, 0.0) - coalesce(u.es, 0.0)) <= 1e-6
    AND abs(coalesce(d.se, 0.0) - coalesce(u.se, 0.0)) <= 1e-6
    AND abs(coalesce(d.lp, 0.0) - coalesce(u.lp, 0.0)) <= 1e-6
    AND abs(coalesce(d.af, 0.0) - coalesce(u.af, 0.0)) <= 1e-6
  )
  ORDER BY id
)
TO '${OUT_DIR}/report/core.details.tsv' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
SQL

core_mismatch_count="$(duckdb -csv -noheader -unsigned -c "SELECT COUNT(*) FROM read_csv_auto('${OUT_DIR}/report/core.details.tsv', delim := '\t', header := true)")"
if [[ "$core_mismatch_count" == "0" ]]; then
  conformance_append_summary_row "$SUMMARY_TSV" "core_columns_file" "match" "match" "0" "${OUT_DIR}/report/core.details.tsv"
else
  conformance_append_summary_row "$SUMMARY_TSV" "core_columns_file" "diff" "unclassified" "$core_mismatch_count" "${OUT_DIR}/report/core.details.tsv"
fi

bolt_input="${TMP_DIR}/bolt.tsv"
cat > "$bolt_input" <<'EOF'
SNP	BP	CHR	ALLELE1	ALLELE0	P_LINREG	BETA	A1FREQ	SE
rs_bolt	2	chrF	A	C	0.01	0.25	0.1	0.05
EOF

set +e
duckdb -unsigned >"${OUT_DIR}/raw/bolt.duck.stderr" 2>&1 <<SQL
LOAD '${DUCKHTS_EXT}';
CREATE OR REPLACE TEMP TABLE bolt_input AS
SELECT * FROM read_csv_auto('${bolt_input}', delim := '\t', header := true);
SELECT * FROM duckdb_munge('bolt_input', preset := 'BOLT', fasta_ref := '${REF_FASTA}');
SQL
duck_status=$?
set -e

"${BCFTOOLS_BIN}" +munge --no-version -c BOLT -f "$REF_FASTA" -Ov \
  -o "${OUT_DIR}/raw/bolt.upstream.vcf" "$bolt_input" \
  >"${OUT_DIR}/raw/bolt.upstream.log" 2>&1

printf 'engine\tstatus\treason\n' > "${OUT_DIR}/report/bolt.details.tsv"
printf 'duckhts\t%s\tpreset_alias_gap\n' "$duck_status" >> "${OUT_DIR}/report/bolt.details.tsv"
printf 'bcftools\t0\tpreset_alias_gap\n' >> "${OUT_DIR}/report/bolt.details.tsv"
if [[ "$duck_status" == "0" ]]; then
  conformance_append_summary_row "$SUMMARY_TSV" "bolt_preset_alias" "match" "match" "0" "${OUT_DIR}/report/bolt.details.tsv"
else
  conformance_append_summary_row "$SUMMARY_TSV" "bolt_preset_alias" "diff" "preset_alias_gap" "1" "${OUT_DIR}/report/bolt.details.tsv"
fi

missing_input="${TMP_DIR}/missing.tsv"
missing_columns="${TMP_DIR}/missing.columns.tsv"
cat > "$missing_input" <<'EOF'
SNP	CHR	BP	A1	A2	P	BETA	FRQ	SE	N
rs_missing	chrNONE	100	A	C	0.01	0.25	0.1	0.05	1000
EOF
cp "$core_columns" "$missing_columns"

duckdb -unsigned >"${OUT_DIR}/raw/missing.duck.log" 2>&1 <<SQL
LOAD '${DUCKHTS_EXT}';
CREATE OR REPLACE TEMP TABLE missing_input AS
SELECT * FROM read_csv_auto('${missing_input}', delim := '\t', header := true);
COPY (
  SELECT chrom, pos, id, ref, alt, filter
  FROM duckdb_munge(
    'missing_input',
    column_map := map(
      ['SNP','CHR','BP','A1','A2','P','BETA','FRQ','SE','N'],
      ['SNP','CHR','BP','A1','A2','P','BETA','FRQ','SE','N']
    ),
    fasta_ref := '${REF_FASTA}'
  )
)
TO '${OUT_DIR}/normalized/missing.duck.tsv' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
SQL

set +e
"${BCFTOOLS_BIN}" +munge --no-version -C "$missing_columns" -f "$REF_FASTA" -Ov \
  -o "${OUT_DIR}/raw/missing.upstream.vcf" "$missing_input" \
  >"${OUT_DIR}/raw/missing.upstream.log" 2>&1
missing_status=$?
set -e

printf 'engine\tstatus\treason\tnote\n' > "${OUT_DIR}/report/missing.details.tsv"
printf 'duckhts\t0\terror_handling_divergence\temits structured row with filter=MissingContig\n' >> "${OUT_DIR}/report/missing.details.tsv"
printf 'bcftools\t%s\terror_handling_divergence\tskips record and reports warning only\n' "$missing_status" >> "${OUT_DIR}/report/missing.details.tsv"
conformance_append_summary_row "$SUMMARY_TSV" "missingcontig_behavior" "diff" "error_handling_divergence" "1" "${OUT_DIR}/report/missing.details.tsv"
