library(DBI)

# A shut-down file database must be released: the same file reopens in this
# process (it could not on Windows while LOAD-time connections pinned it) and,
# where /proc exists, no descriptor for it remains.
test_database_lifetime <- function() {
  database_file <- tempfile("rduckhts_lifetime_", fileext = ".duckdb")
  on.exit(unlink(database_file), add = TRUE)
  extdata <- function(name) system.file("extdata", name, package = "Rduckhts")
  paths <- vapply(
    c("range.bam", "range.bam.bai", "ce.fa", "ce.fa.fai"),
    extdata, character(1)
  )
  expect_true(all(nzchar(paths)) && all(file.exists(paths)))

  con <- rduckhts_connect(dbdir = database_file)
  dbExecute(con, paste(
    "CREATE TEMP TABLE lifetime_targets AS SELECT * FROM (VALUES",
    "('chr1', 10, 20, 'a'), ('chr1', 15, 30, 'b'))",
    "t(chrom, start, \"end\", label)"
  ))
  expect_equal(
    dbGetQuery(con, paste(
      "SELECT * FROM duckhts_cgranges_from_table('lifetime_idx',",
      "'lifetime_targets', 'chrom', 'start', 'end', 'label') AS t(ok)"
    ))$ok[1],
    TRUE
  )
  expect_equal(
    dbGetQuery(con, paste(
      "SELECT count(*) AS n FROM (SELECT unnest(",
      "duckhts_cgranges_overlaps_list('lifetime_idx', 'chr1', 16, 18)))"
    ))$n[1],
    2
  )
  dbExecute(con, "SELECT duckhts_cgranges_destroy('lifetime_idx')")
  dbExecute(con, paste(
    "CREATE TABLE lifetime_panel AS SELECT * FROM (VALUES",
    "('WBcel235', 0::UBIGINT, 'CHROMOSOME_I', 914::UBIGINT, 'A', 'C'))",
    "p(assembly, site_index, region, position, allele_a, allele_b)"
  ))
  counts <- rduckhts_somalier_bam_counts(
    con, paths[["range.bam"]], "sample-1", paths[["ce.fa"]],
    panel_table = "lifetime_panel", index_path = paths[["range.bam.bai"]],
    reference_index_path = paths[["ce.fa.fai"]]
  )
  expect_equal(counts$status, "measured")
  dbDisconnect(con, shutdown = TRUE)
  rm(con)
  invisible(gc())

  if (dir.exists("/proc/self/fd")) {
    target <- normalizePath(database_file, mustWork = TRUE)
    links <- Sys.readlink(list.files("/proc/self/fd", full.names = TRUE))
    expect_false(target %in% normalizePath(links[nzchar(links)], mustWork = FALSE))
  }

  reopened <- rduckhts_connect(dbdir = database_file)
  expect_equal(
    dbGetQuery(reopened, "SELECT count(*) AS n FROM lifetime_panel")$n[1],
    1
  )
  dbDisconnect(reopened, shutdown = TRUE)
}

test_database_lifetime()
