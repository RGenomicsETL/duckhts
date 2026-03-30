#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${1:?output directory required}"
ROOT_DIR="${ROOT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
source "${ROOT_DIR}/scripts/conformance/common.sh"

conformance_prepare_dirs "$OUT_DIR"
SUMMARY_TSV="${OUT_DIR}/report/summary.tsv"
conformance_write_summary_header "$SUMMARY_TSV"

CORE_SRC_FA="${ROOT_DIR}/test/data/liftover_src.fa"
CORE_DST_FA="${ROOT_DIR}/test/data/liftover_dst.fa"
CORE_CHAIN="${ROOT_DIR}/test/data/liftover.chain"
MT_SRC_FA="${ROOT_DIR}/test/data/liftover_mt_src.fa"
MT_DST_FA="${ROOT_DIR}/test/data/liftover_mt_dst.fa"
MT_CHAIN="${ROOT_DIR}/test/data/liftover_mt.chain"

run_liftover_compare_case() {
  local case_name="$1"
  local input_vcf="$2"
  local src_fa="$3"
  local dst_fa="$4"
  local chain_path="$5"
  local plugin_opts="$6"
  local duck_sql="$7"
  local upstream_norm_sql="$8"
  local detail_tsv="${OUT_DIR}/report/${case_name}.details.tsv"
  local mismatch_count

  "${BCFTOOLS_BIN}" +liftover -Ov -o "${OUT_DIR}/raw/${case_name}.upstream.vcf" "$input_vcf" -- \
    -s "$src_fa" -f "$dst_fa" -c "$chain_path" $plugin_opts \
    >"${OUT_DIR}/raw/${case_name}.upstream.log" 2>&1

  duckdb -unsigned >"${OUT_DIR}/raw/${case_name}.duck.log" 2>&1 <<SQL
LOAD '${DUCKHTS_EXT}';
COPY (
${duck_sql}
)
TO '${OUT_DIR}/normalized/${case_name}.duck.tsv' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
${upstream_norm_sql}
COPY (
  WITH duck_norm AS (
    SELECT src_chrom,
           CAST(src_pos AS BIGINT) AS src_pos,
           dest_chrom,
           CAST(dest_pos AS BIGINT) AS dest_pos,
           CAST(dest_end AS BIGINT) AS dest_end,
           dest_ref,
           dest_alt,
           CAST(reverse_complemented AS BOOLEAN) AS reverse_complemented,
           CAST(swap AS INTEGER) AS swap,
           note
    FROM read_csv_auto('${OUT_DIR}/normalized/${case_name}.duck.tsv', delim := '\t', header := true)
  ),
  upstream_norm AS (
    SELECT src_chrom,
           CAST(src_pos AS BIGINT) AS src_pos,
           dest_chrom,
           CAST(dest_pos AS BIGINT) AS dest_pos,
           CAST(dest_end AS BIGINT) AS dest_end,
           dest_ref,
           dest_alt,
           CAST(reverse_complemented AS BOOLEAN) AS reverse_complemented,
           CAST(swap AS INTEGER) AS swap,
           note
    FROM read_csv_auto('${OUT_DIR}/normalized/${case_name}.upstream.tsv', delim := '\t', header := true)
  )
  SELECT
    coalesce(d.src_chrom, u.src_chrom) AS src_chrom,
    coalesce(d.src_pos, u.src_pos) AS src_pos,
    d.dest_chrom AS duck_dest_chrom,
    u.dest_chrom AS upstream_dest_chrom,
    d.dest_pos AS duck_dest_pos,
    u.dest_pos AS upstream_dest_pos,
    d.dest_end AS duck_dest_end,
    u.dest_end AS upstream_dest_end,
    d.dest_ref AS duck_dest_ref,
    u.dest_ref AS upstream_dest_ref,
    d.dest_alt AS duck_dest_alt,
    u.dest_alt AS upstream_dest_alt,
    d.reverse_complemented AS duck_reverse_complemented,
    u.reverse_complemented AS upstream_reverse_complemented,
    d.swap AS duck_swap,
    u.swap AS upstream_swap,
    d.note AS duck_note,
    u.note AS upstream_note,
    CASE
      WHEN d.src_chrom IS NULL THEN 'only_bcftools'
      WHEN u.src_chrom IS NULL THEN 'only_duckhts'
      WHEN coalesce(d.dest_chrom, '') = coalesce(u.dest_chrom, '')
       AND coalesce(d.dest_pos, -1) = coalesce(u.dest_pos, -1)
       AND coalesce(d.dest_end, -1) = coalesce(u.dest_end, -1)
       AND coalesce(d.dest_ref, '') = coalesce(u.dest_ref, '')
       AND coalesce(d.dest_alt, '') = coalesce(u.dest_alt, '')
       AND coalesce(d.reverse_complemented, false) = coalesce(u.reverse_complemented, false)
       AND coalesce(d.swap, 0) = coalesce(u.swap, 0)
        THEN 'match'
      ELSE 'diff_values'
    END AS status,
    CASE
      WHEN d.src_chrom IS NULL OR u.src_chrom IS NULL THEN 'surface_difference'
      ELSE 'unclassified'
    END AS reason
  FROM duck_norm d
  FULL OUTER JOIN upstream_norm u USING (src_chrom, src_pos)
  WHERE NOT (
    d.src_chrom IS NOT NULL AND u.src_chrom IS NOT NULL
    AND coalesce(d.dest_chrom, '') = coalesce(u.dest_chrom, '')
    AND coalesce(d.dest_pos, -1) = coalesce(u.dest_pos, -1)
    AND coalesce(d.dest_end, -1) = coalesce(u.dest_end, -1)
    AND coalesce(d.dest_ref, '') = coalesce(u.dest_ref, '')
    AND coalesce(d.dest_alt, '') = coalesce(u.dest_alt, '')
    AND coalesce(d.reverse_complemented, false) = coalesce(u.reverse_complemented, false)
    AND coalesce(d.swap, 0) = coalesce(u.swap, 0)
  )
  ORDER BY src_chrom, src_pos
)
TO '${detail_tsv}' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);
SQL

  mismatch_count="$(duckdb -csv -noheader -unsigned -c "SELECT COUNT(*) FROM read_csv_auto('${detail_tsv}', delim := '\t', header := true)")"
  if [[ "$mismatch_count" == "0" ]]; then
    conformance_append_summary_row "$SUMMARY_TSV" "$case_name" "match" "match" "0" "${detail_tsv}"
  else
    conformance_append_summary_row "$SUMMARY_TSV" "$case_name" "diff" "unclassified" "$mismatch_count" "${detail_tsv}"
  fi
}

core_len="$(awk 'NR == 1 { print $2 }' "${CORE_SRC_FA}.fai")"
cat > "${OUT_DIR}/tmp/core.vcf" <<EOF
##fileformat=VCFv4.2
##contig=<ID=chrF,length=${core_len}>
##contig=<ID=chrR,length=${core_len}>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chrF	2	rs1	C	T	.	PASS	.
chrR	2	rs2	A	G	.	PASS	.
EOF

run_liftover_compare_case \
  "core" \
  "${OUT_DIR}/tmp/core.vcf" \
  "${CORE_SRC_FA}" \
  "${CORE_DST_FA}" \
  "${CORE_CHAIN}" \
  "" \
  "SELECT lo.src_chrom, lo.src_pos, lo.dest_chrom, lo.dest_pos,
          CAST(NULL AS BIGINT) AS dest_end,
          lo.dest_ref, lo.dest_alt, lo.reverse_complemented, lo.swap, lo.note
   FROM (
     SELECT bcftools_liftover(chrom, pos, ref, alt,
       '${CORE_CHAIN}', '${CORE_DST_FA}', '${CORE_SRC_FA}', 1, 250, false, NULL::BIGINT) AS lo
     FROM (VALUES ('chrF', 2, 'C', 'T'), ('chrR', 2, 'A', 'G')) AS t(chrom, pos, ref, alt)
   ) s
   ORDER BY lo.src_chrom, lo.src_pos" \
  "COPY (
     SELECT
       CASE WHEN ID = 'rs1' THEN 'chrF' ELSE 'chrR' END AS src_chrom,
       2::BIGINT AS src_pos,
       CHROM AS dest_chrom,
       POS AS dest_pos,
       CAST(NULL AS BIGINT) AS dest_end,
       REF AS dest_ref,
       ALT[1] AS dest_alt,
       coalesce(INFO_FLIP, false) AS reverse_complemented,
       coalesce(INFO_SWAP, 0) AS swap,
       CAST(NULL AS VARCHAR) AS note
     FROM read_bcf('${OUT_DIR}/raw/core.upstream.vcf')
     ORDER BY src_chrom, src_pos
   )
   TO '${OUT_DIR}/normalized/core.upstream.tsv' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);"

mt_len="$(awk 'NR == 1 { print $2 }' "${MT_SRC_FA}.fai")"
cat > "${OUT_DIR}/tmp/mt.vcf" <<EOF
##fileformat=VCFv4.2
##contig=<ID=MT,length=${mt_len}>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
MT	5	rs_mt	G	A	.	PASS	.
EOF

run_liftover_compare_case \
  "mt_passthrough" \
  "${OUT_DIR}/tmp/mt.vcf" \
  "${MT_SRC_FA}" \
  "${MT_DST_FA}" \
  "${MT_CHAIN}" \
  "" \
  "SELECT lo.src_chrom, lo.src_pos, lo.dest_chrom, lo.dest_pos,
          CAST(NULL AS BIGINT) AS dest_end,
          lo.dest_ref, lo.dest_alt, lo.reverse_complemented, lo.swap, lo.note
   FROM (
     SELECT bcftools_liftover('MT', 5, 'G', 'A',
       '${MT_CHAIN}', '${MT_DST_FA}', '${MT_SRC_FA}', 1, 250, false, NULL::BIGINT) AS lo
   ) s" \
  "COPY (
     SELECT
       'MT' AS src_chrom,
       5::BIGINT AS src_pos,
       CHROM AS dest_chrom,
       POS AS dest_pos,
       CAST(NULL AS BIGINT) AS dest_end,
       REF AS dest_ref,
       ALT[1] AS dest_alt,
       coalesce(INFO_FLIP, false) AS reverse_complemented,
       coalesce(INFO_SWAP, 0) AS swap,
       'MitochondriaPassthrough' AS note
     FROM read_bcf('${OUT_DIR}/raw/mt_passthrough.upstream.vcf')
   )
   TO '${OUT_DIR}/normalized/mt_passthrough.upstream.tsv' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);"

cat > "${OUT_DIR}/tmp/end.vcf" <<EOF
##fileformat=VCFv4.2
##contig=<ID=chrF,length=${core_len}>
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chrF	2	rs_end	C	T	.	PASS	END=5
EOF

run_liftover_compare_case \
  "lift_end" \
  "${OUT_DIR}/tmp/end.vcf" \
  "${CORE_SRC_FA}" \
  "${CORE_DST_FA}" \
  "${CORE_CHAIN}" \
  "--lift-end" \
  "SELECT lo.src_chrom, lo.src_pos, lo.dest_chrom, lo.dest_pos,
          lo.dest_end,
          lo.dest_ref, lo.dest_alt, lo.reverse_complemented, lo.swap, lo.note
   FROM (
     SELECT bcftools_liftover('chrF', 2, 'C', 'T',
       '${CORE_CHAIN}', '${CORE_DST_FA}', '${CORE_SRC_FA}', 1, 250, false, 5::BIGINT) AS lo
   ) s" \
  "COPY (
     SELECT
       'chrF' AS src_chrom,
       2::BIGINT AS src_pos,
       CHROM AS dest_chrom,
       POS AS dest_pos,
       INFO_END::BIGINT AS dest_end,
       REF AS dest_ref,
       ALT[1] AS dest_alt,
       coalesce(INFO_FLIP, false) AS reverse_complemented,
       coalesce(INFO_SWAP, 0) AS swap,
       CAST(NULL AS VARCHAR) AS note
     FROM read_bcf('${OUT_DIR}/raw/lift_end.upstream.vcf')
   )
   TO '${OUT_DIR}/normalized/lift_end.upstream.tsv' (FORMAT CSV, DELIMITER '\t', HEADER TRUE);"
