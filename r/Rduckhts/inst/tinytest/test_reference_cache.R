library(tinytest)
library(DBI)

test_reference_cache <- function() {
  directory <- tempfile("duckhts_reference_cache_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE, after = FALSE)
  motifs <- c("ACGT", "CGTA", "GTAC", "TACG")
  paths <- character(12)
  for (i in seq_along(paths)) {
    path <- file.path(directory, paste0("ref-", i, ".fa"))
    writeLines(c(">chrCache", strrep(motifs[(i - 1L) %% 4L + 1L], 65536)), path)
    if (i %% 2L == 0L) {
      expect_true(rduckhts_bgzip(con, path, output_path = paste0(path, ".gz"),
        keep = TRUE, overwrite = TRUE)$success[1])
      path <- paste0(path, ".gz")
    }
    expect_true(rduckhts_fasta_index(con, path, index_path = paste0(path, ".fai"))$success[1])
    paths[i] <- path
  }
  dbWriteTable(con, "cache_references", data.frame(
    file_id = 0:11, path = paths, motif = rep(motifs, 3)))
  dbExecute(con, paste(
    "CREATE TABLE cache_input AS SELECT i, path, 1 + i % 200000 AS pos,",
    "substr(motif, 1 + i % 4, 1) AS ref FROM range(1000000) t(i)",
    "JOIN cache_references ON file_id = (i // 8192) % 12 ORDER BY i"))
  dbExecute(con, paste(
    "CREATE MACRO cache_run(input_table) AS TABLE",
    "SELECT i, pos, ref, CASE WHEN ref = 'A' THEN 'C' ELSE 'A' END AS alt,",
    "bcftools_norm_row('chrCache',pos,ref,[alt],path,NULL,NULL,NULL,NULL) AS norm,",
    "bcftools_munge_row('chrCache',pos,alt,ref,NULL,",
    "NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,",
    "path,'IFFY','REF_MISMATCH',NULL,NULL,NULL) AS munge FROM query_table(input_table)"))
  check_sql <- paste(
    "SELECT count(*) AS n, sum(i)::DOUBLE AS input_ids, count(*) FILTER (WHERE",
    "norm.pos_normed = pos AND norm.ref_normed = ref AND norm.alt_normed = [alt]",
    "AND norm.norm_status = 'Unchanged' AND NOT norm.normed",
    "AND munge.pos = pos AND munge.ref = ref AND munge.alt = alt",
    "AND munge.filter IS NULL AND NOT munge.alleles_swapped) AS correct",
    "FROM cache_run('%s')")
  for (threads in c(1L, 4L, 8L)) {
    dbExecute(con, paste("SET threads =", threads))
    out <- dbGetQuery(con, sprintf(check_sql, "cache_input"))
    expect_equal(as.numeric(out$n), 1000000)
    expect_equal(out$input_ids, 999999 * 1000000 / 2)
    expect_equal(as.numeric(out$correct), 1000000)
  }
  dbExecute(con, "CREATE TABLE cache_reverse AS SELECT * FROM cache_input ORDER BY i DESC")
  reverse <- dbGetQuery(con, sprintf(check_sql, "cache_reverse"))
  expect_equal(reverse, out)
}

test_reference_cache()
