#!/usr/bin/env bash
# DuckDB-CLI cgranges benchmark over BED providers.
#
# This intentionally does NOT generate one duckhts_cgranges_overlaps(...) query
# per probe. DuckHTS variants stream query BED rows through read_bed(...) and the
# vectorized scalar helpers:
#   * duckhts_cgranges_has_overlap(...)
#   * duckhts_cgranges_count_overlaps(...)
#   * duckhts_cgranges_overlaps_list(...) + UNNEST for one-row-per-hit expansion
#
# Baselines:
#   * bedtk flt
#   * bedtools intersect -u / -c, when bedtools is available

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

EXT="${EXT:-$REPO_ROOT/build/release/duckhts.duckdb_extension}"
BEDTK="${BEDTK:-$REPO_ROOT/.sync/bedtk/bedtk}"
BEDTOOLS="${BEDTOOLS:-bedtools}"
SUBJECT="${SUBJECT:-$REPO_ROOT/DuckBedQC/data/GRCh38_exons.bed}"
QUERY="${QUERY:-$REPO_ROOT/DuckBedQC/data/GRCh38_illumina_clinical_regions_v100.39.0.bed}"
OUT_DIR="${OUT_DIR:-$REPO_ROOT/.tmp/cgranges_benchmark_cli}"

mkdir -p "$OUT_DIR"
SUMMARY="$OUT_DIR/summary.tsv"
printf 'tool\tvariant\telapsed_sec\tpeak_rss_mb\tmatched_query_intervals\ttotal_hits\n' > "$SUMMARY"

sql_quote() {
  local s="$1"
  printf "%s" "${s//\'/\'\'}"
}

elapsed_seconds() {
  awk -F': ' '/Elapsed \(wall clock\)/ {print $NF}' "$1" | awk -F: '{
    n=NF; s=0;
    for (i=1;i<=n;i++) s = s*60 + $i;
    printf "%.3f", s
  }'
}

rss_mb() {
  awk -F': ' '/Maximum resident set size/ {printf "%.1f", $NF/1024.0}' "$1"
}

append_summary() {
  local tool="$1" variant="$2" log="$3" matched="$4" hits="$5"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$tool" "$variant" "$(elapsed_seconds "$log")" "$(rss_mb "$log")" "$matched" "$hits" >> "$SUMMARY"
}

run_duckdb_variant() {
  local variant="$1"
  local sql="$OUT_DIR/duckhts_${variant}.sql"
  local stdout_log="$OUT_DIR/duckhts_${variant}.stdout"
  local stderr_log="$OUT_DIR/duckhts_${variant}.stderr"
  local time_log="$OUT_DIR/duckhts_${variant}.time"
  local subject_q query_q ext_q
  subject_q="$(sql_quote "$SUBJECT")"
  query_q="$(sql_quote "$QUERY")"
  ext_q="$(sql_quote "$EXT")"

  cat > "$sql" <<SQL
LOAD '$ext_q';
CREATE TABLE cgr_build AS
  SELECT duckhts_cgranges_from_query(
    'bench',
    'SELECT chrom, start, "end" FROM read_bed(''$subject_q'')',
    'chrom', 'start', 'end'
  ) AS ok;
CREATE TABLE cgr_index AS SELECT duckhts_cgranges_index('bench') AS ok;
SQL

  if [[ "$variant" == "scalar_filter" ]]; then
    cat >> "$sql" <<SQL
COPY (
  SELECT count(*)::BIGINT AS matched_query_intervals, NULL::BIGINT AS total_hits
  FROM read_bed('$query_q') q
  WHERE duckhts_cgranges_has_overlap('bench', q.chrom, q.start, q."end")
) TO STDOUT (HEADER false, DELIMITER '\t');
SQL
  elif [[ "$variant" == "scalar_count" ]]; then
    cat >> "$sql" <<SQL
COPY (
  WITH counts AS (
    SELECT duckhts_cgranges_count_overlaps('bench', q.chrom, q.start, q."end") AS n
    FROM read_bed('$query_q') q
  )
  SELECT count(*) FILTER (WHERE n > 0)::BIGINT AS matched_query_intervals,
         COALESCE(sum(n), 0)::BIGINT AS total_hits
  FROM counts
) TO STDOUT (HEADER false, DELIMITER '\t');
SQL
  elif [[ "$variant" == "scalar_expand" ]]; then
    cat >> "$sql" <<SQL
COPY (
  WITH q AS (
    SELECT row_number() OVER () AS qid, chrom, start, "end"
    FROM read_bed('$query_q')
  ), hits AS (
    SELECT qid
    FROM q
    CROSS JOIN UNNEST(duckhts_cgranges_overlaps_list('bench', q.chrom, q.start, q."end")) AS u(hit)
  )
  SELECT count(DISTINCT qid)::BIGINT AS matched_query_intervals,
         count(*)::BIGINT AS total_hits
  FROM hits
) TO STDOUT (HEADER false, DELIMITER '\t');
SQL
  else
    echo "unknown DuckHTS variant: $variant" >&2
    return 1
  fi

  /usr/bin/time -v -o "$time_log" duckdb -unsigned -batch < "$sql" > "$stdout_log" 2> "$stderr_log"
  local result matched hits
  result="$(tail -n 1 "$stdout_log")"
  matched="$(printf '%s' "$result" | awk -F'\t' '{print $1}')"
  hits="$(printf '%s' "$result" | awk -F'\t' '{print $2}')"
  append_summary "duckhts" "$variant" "$time_log" "$matched" "${hits:-NA}"
  echo "duckhts $variant: matched=$matched hits=${hits:-NA}"
}

run_bedtk() {
  if [[ ! -x "$BEDTK" ]]; then
    echo "Skipping bedtk: executable not found at $BEDTK" >&2
    return 0
  fi
  local time_log="$OUT_DIR/bedtk_flt.time"
  local out="$OUT_DIR/bedtk_flt.count"
  /usr/bin/time -v -o "$time_log" bash -c "'$BEDTK' flt '$SUBJECT' '$QUERY' | awk 'END { print NR }' > '$out'"
  append_summary "bedtk" "flt" "$time_log" "$(cat "$out")" "NA"
  echo "bedtk flt: matched=$(cat "$out")"
}

run_bedtools() {
  local bedtools_bin
  bedtools_bin="$(command -v "$BEDTOOLS" || true)"
  if [[ -z "$bedtools_bin" ]]; then
    echo "Skipping bedtools: executable not found ($BEDTOOLS)" >&2
    return 0
  fi

  local time_log out matched hits
  time_log="$OUT_DIR/bedtools_intersect_u.time"
  out="$OUT_DIR/bedtools_intersect_u.count"
  /usr/bin/time -v -o "$time_log" bash -c "'$bedtools_bin' intersect -a '$QUERY' -b '$SUBJECT' -u | awk 'END { print NR }' > '$out'"
  matched="$(cat "$out")"
  append_summary "bedtools" "intersect_u" "$time_log" "$matched" "NA"
  echo "bedtools intersect -u: matched=$matched"

  time_log="$OUT_DIR/bedtools_intersect_c.time"
  out="$OUT_DIR/bedtools_intersect_c.counts"
  /usr/bin/time -v -o "$time_log" bash -c "'$bedtools_bin' intersect -a '$QUERY' -b '$SUBJECT' -c | awk -F '\t' '{ if (\$NF > 0) matched++; hits += \$NF } END { print matched+0 \"\t\" hits+0 }' > '$out'"
  matched="$(awk -F'\t' '{print $1}' "$out")"
  hits="$(awk -F'\t' '{print $2}' "$out")"
  append_summary "bedtools" "intersect_c" "$time_log" "$matched" "$hits"
  echo "bedtools intersect -c: matched=$matched hits=$hits"
}

echo "Subject: $SUBJECT ($(grep -cv '^#' "$SUBJECT") intervals)"
echo "Query  : $QUERY ($(grep -cv '^#' "$QUERY") intervals)"
echo "Output : $OUT_DIR"
echo

run_duckdb_variant scalar_filter
run_duckdb_variant scalar_count
run_duckdb_variant scalar_expand
run_bedtk
run_bedtools

echo
echo "=== Summary ($SUMMARY) ==="
if command -v column >/dev/null 2>&1; then
  column -s $'\t' -t "$SUMMARY"
else
  cat "$SUMMARY"
fi
