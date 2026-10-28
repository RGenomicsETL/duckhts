.bcf_snapshot_statements <- function(con, reader, path, index_paths, replacement,
                                     min_shift, count = FALSE) {
  quote <- function(x) as.character(DBI::dbQuoteString(con, x))
  rewrite_index <- function(data, prefix) {
    sprintf("bcf_index(%s,index_path:=%s,min_shift:=%d,threads:=1) AS %s%d",
            quote(data), quote(index_paths), min_shift, prefix,
            seq_along(index_paths))
  }
  rewrite <- function(sources) {
    paste("SELECT count(*) AS n FROM", paste(sources, collapse = " CROSS JOIN "))
  }
  replace <- rewrite_index(replacement, "replaced")
  replaced <- paste0("replaced", seq_along(index_paths), ".success", collapse = " AND ")
  projection <- if (count) "count(*) AS n" else "records.POS, records.* EXCLUDE (POS)"
  # bcf_index builds during execution. DuckDB finishes the build side of a cross
  # product before the probe side scans, and the one-row bcf_index results are
  # the build sides, so the replacement runs after the reader binds and before
  # it scans. `window` checks that order: read_blob reads the index bytes when
  # it scans, so it must see the replacement.
  list(
    restore = rewrite(rewrite_index(path, "restored")),
    replace = rewrite(replace),
    after_bind = paste("SELECT", projection, "FROM",
                       paste(c(paste(reader, "AS records"), replace), collapse = " CROSS JOIN "),
                       "WHERE", replaced, "ORDER BY ALL"),
    window = paste("SELECT md5(blob.content) AS digest FROM",
                   paste(c(sprintf("read_blob(%s) AS blob", quote(index_paths[[1]])), replace),
                         collapse = " CROSS JOIN "),
                   "WHERE", replaced)
  )
}
