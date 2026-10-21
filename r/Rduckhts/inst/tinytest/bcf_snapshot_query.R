.bcf_snapshot_query <- function(con, reader, path, index_paths, replacement,
                                min_shift, count = FALSE,
                                replace_before_bind = FALSE) {
  quote <- function(x) as.character(DBI::dbQuoteString(con, x))
  rewrite_index <- function(data, prefix) {
    sprintf("bcf_index(%s,index_path:=%s,min_shift:=%d,threads:=1) AS %s%d",
            quote(data), quote(index_paths), min_shift, prefix,
            seq_along(index_paths))
  }
  # bcf_index writes in its bind callback. FROM items bind left to right, so
  # each bind restores the pair, binds the reader, then replaces its sidecars.
  restore <- rewrite_index(path, "restored")
  replace <- rewrite_index(replacement, "replaced")
  read <- paste(reader, "AS records")
  sources <- if (replace_before_bind) c(restore, replace, read) else
    c(restore, read, replace)
  aliases <- c(paste0("restored", seq_along(index_paths)),
               paste0("replaced", seq_along(index_paths)))
  projection <- if (count) "count(*) AS n" else "records.POS, records.* EXCLUDE (POS)"
  paste("SELECT", projection, "FROM", paste(sources, collapse = " CROSS JOIN "),
        "WHERE", paste0(aliases, ".success", collapse = " AND "), "ORDER BY ALL")
}
