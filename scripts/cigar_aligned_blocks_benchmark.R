# Independent quadratic packed-CIGAR oracle. Missing CIGARs produce NULL;
# nonempty CIGARs without M, = or X produce empty lists.
sql_blocks <- paste0(
  "CASE WHEN len(CIGAR) > 0 AND POS IS NOT NULL THEN ",
  "list_transform(list_filter(list_transform(CIGAR, lambda c, i: {'kind': c & 15, ",
  "'v': POS + COALESCE(list_sum(list_transform(list_slice(CIGAR, 1, i - 1), ",
  "lambda p: CASE WHEN (p & 15) IN (0, 2, 3, 7, 8) THEN (p >> 4)::BIGINT ELSE 0 END)), 0)}), ",
  "lambda x: x.kind IN (0, 7, 8)), lambda x: x.v) END AS ref_start, ",
  "CASE WHEN len(CIGAR) > 0 AND POS IS NOT NULL THEN ",
  "list_transform(list_filter(list_transform(CIGAR, lambda c, i: {'kind': c & 15, ",
  "'v': COALESCE(list_sum(list_transform(list_slice(CIGAR, 1, i - 1), ",
  "lambda p: CASE WHEN (p & 15) IN (0, 1, 4, 7, 8) THEN (p >> 4)::BIGINT ELSE 0 END)), 0)}), ",
  "lambda x: x.kind IN (0, 7, 8)), lambda x: x.v) END AS query_start, ",
  "CASE WHEN len(CIGAR) > 0 AND POS IS NOT NULL THEN ",
  "list_transform(list_filter(CIGAR, lambda c: (c & 15) IN (0, 7, 8)), ",
  "lambda c: (c >> 4)::BIGINT) END AS width"
)
scalar_blocks <- "(b).ref_start AS ref_start, (b).query_start AS query_start, (b).width AS width"

# Compare every physical input row before timing; logical read IDs need not be unique.
validate_blocks <- function(con, source) {
  result <- DBI::dbGetQuery(con, sprintf(paste0(
    "SELECT count(*)::VARCHAR AS records, count(*) FILTER (WHERE ",
    "(b IS NULL) IS DISTINCT FROM missing OR ",
    "ref_start IS DISTINCT FROM (b).ref_start OR ",
    "query_start IS DISTINCT FROM (b).query_start OR ",
    "width IS DISTINCT FROM (b).width)::VARCHAR AS mismatches ",
    "FROM (SELECT %s, CIGAR IS NULL OR len(CIGAR) = 0 OR POS IS NULL AS missing, ",
    "cigar_aligned_blocks(CIGAR, POS) AS b FROM %s)"
  ), sql_blocks, source))
  if (result$mismatches != "0") {
    stop(result$mismatches, " physical-record geometry mismatches out of ", result$records)
  }
  invisible(result)
}

# Timed proofs include the physical ordinal, explicit NULL markers, XOR and
# an independent sum reduction. They supplement the exact untimed comparison.
workload_sql <- function(source, workload) {
  records <- paste0("(SELECT row_number() OVER () AS record_id, QNAME, FLAG, POS, CIGAR FROM ",
                    source, ")")
  if (workload == "read_bam only") {
    projection <- paste0(
      "SELECT len(CIGAR) AS ops, 0::BIGINT AS blocks, ",
      "hash(record_id, QNAME, FLAG, POS, CIGAR IS NULL, CIGAR) AS fingerprint FROM ", records
    )
  } else {
    geometry <- switch(workload,
      "blocks in SQL" = paste0(
        "SELECT record_id, QNAME, FLAG, POS, len(CIGAR) AS ops, ",
        "CIGAR IS NULL OR len(CIGAR) = 0 OR POS IS NULL AS missing, ",
        sql_blocks, " FROM ", records
      ),
      "cigar_aligned_blocks" = paste0(
        "SELECT record_id, QNAME, FLAG, POS, len(CIGAR) AS ops, b IS NULL AS missing, ", scalar_blocks,
        " FROM (SELECT *, cigar_aligned_blocks(CIGAR, POS) AS b FROM ", records, ")"
      ),
      stop("unknown workload")
    )
    projection <- paste0(
      "SELECT ops, coalesce(len(width), 0) AS blocks, ",
      "hash(record_id, QNAME, FLAG, POS, missing, ref_start IS NULL, ref_start, ",
      "query_start IS NULL, query_start, width IS NULL, width) AS fingerprint FROM (", geometry, ")"
    )
  }
  paste0(
    "SELECT count(*)::VARCHAR AS reads, coalesce(sum(ops), 0)::VARCHAR AS ops, ",
    "sum(blocks)::VARCHAR AS blocks, hex(bit_xor(fingerprint)) AS fingerprint, ",
    "sum(fingerprint::HUGEINT)::VARCHAR AS fingerprint_sum FROM (", projection, ")"
  )
}
