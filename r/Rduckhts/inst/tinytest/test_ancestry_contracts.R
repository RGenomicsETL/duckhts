library(tinytest)
library(DBI)

test_ancestry_checked_sites <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  dbExecute(con, "SET threads=1")
  pcs <- paste0("sin(position * ", seq_len(16L), " / 10.0) AS PC", seq_len(16L))
  dbExecute(con, paste0(
    "CREATE TABLE ref AS SELECT 1 AS chromosome, position, 'A' AS allele_a, ",
    "'G' AS allele_b, ", paste(pcs, collapse = ", "),
    ", 0.1 + position * 0.03 AS A, 0.9 - position * 0.02 AS B, ",
    "0.5 + cos(position)/5 AS C FROM range(1, 13) t(position)"))
  dbExecute(con, paste0(
    "CREATE TABLE input AS SELECT 'S' AS sample_id, 'chr1' AS chromosome, ",
    "position, allele_a, allele_b, ",
    "0.712345645*A+0.2*B+0.087654355*C+sin(position*0.4)*0.015 ",
    "AS frequency FROM ref"))
  dbExecute(con, "CREATE TABLE correction AS SELECT pc, 1.0 AS coefficient FROM range(1, 17) t(pc)")
  run <- function(input_table = "input", reference_table = "ref", min_cor = -1) {
    rduckhts_ancestry_proportions(con, input_table, reference_table,
                                  "correction", min_cor = min_cor)
  }
  baseline <- run()
  expect_equal(unique(baseline$input_variants), 12)
  expect_equal(unique(baseline$used_variants), 12)
  extra_pcs <- paste0("sin(position * ", 17:64, " / 10.0) AS PC", 17:64)
  dbExecute(con, paste0("CREATE VIEW ref64 AS SELECT *, ",
                        paste(extra_pcs, collapse = ", "), " FROM ref"))
  dbExecute(con, "CREATE VIEW correction64 AS SELECT pc, 1.0 AS coefficient FROM range(1, 65) t(pc)")
  full_pc <- rduckhts_ancestry_proportions(con, "input", "ref64", "correction64",
                                            min_cor = -1)
  expect_equal(nrow(full_pc), 3L)
  expect_equal(unique(full_pc$used_variants), 12)
  dbExecute(con, paste0("CREATE VIEW alias_input AS SELECT * FROM input UNION ALL ",
                        "SELECT * REPLACE ('1' AS chromosome) FROM input WHERE position=1"))
  alias <- run("alias_input")
  expect_equal(unique(alias$input_variants), 13)
  expect_equal(unique(alias$used_variants), 11)
  expect_equal(unique(alias$duplicate_variants), 2)
  dbExecute(con, "CREATE VIEW duplicated_ref AS SELECT * FROM ref UNION ALL SELECT * FROM ref WHERE position=1")
  expect_equal(dbGetQuery(con, paste0("SELECT count(*) AS n FROM input i JOIN ",
    "duplicated_ref r ON regexp_replace(i.chromosome::VARCHAR, '^chr', '') = ",
    "r.chromosome::VARCHAR AND i.position = r.position"))$n, 13)
  expect_error(run(reference_table = "duplicated_ref"), "duplicate reference locus")
  dbExecute(con, paste0("CREATE VIEW discarded_missing AS SELECT * REPLACE (",
    "CASE WHEN position=1 THEN NULL ELSE frequency END AS frequency) FROM input"))
  dbExecute(con, paste0("CREATE VIEW discarded_ambiguous AS SELECT * REPLACE (",
    "CASE WHEN position=1 THEN 'T' ELSE allele_b END AS allele_b) FROM input"))
  dbExecute(con, paste0("CREATE VIEW discarded_unmatched AS SELECT * REPLACE (",
    "CASE WHEN position=1 THEN 'C' ELSE allele_b END AS allele_b) FROM input"))
  for (view in c("discarded_missing", "discarded_ambiguous", "discarded_unmatched")) {
    expect_equal(dbGetQuery(con, paste0("SELECT count(*) AS n FROM ", view))$n, 12)
    expect_error(run(view, "duplicated_ref"), "duplicate reference locus")
  }
  dbExecute(con, paste0("CREATE VIEW duplicated_input AS SELECT * FROM input UNION ALL ",
                        "SELECT * FROM input WHERE position=1"))
  expect_equal(dbGetQuery(con, "SELECT count(*) AS n FROM duplicated_input")$n, 13)
  expect_equal(dbGetQuery(con, paste0("SELECT count(*) AS n FROM duplicated_input i ",
    "JOIN duplicated_ref r ON regexp_replace(i.chromosome::VARCHAR, '^chr', '') = ",
    "r.chromosome::VARCHAR AND i.position = r.position"))$n, 15)
  expect_error(run("duplicated_input", "duplicated_ref"), "duplicate reference locus")
  dbExecute(con, paste0("CREATE VIEW aliased_ref AS SELECT * FROM ref UNION ALL ",
                        "SELECT * REPLACE ('chr1' AS chromosome) FROM ref WHERE position=1"))
  expect_error(run(reference_table = "aliased_ref"), "duplicate reference locus")
  dbExecute(con, paste0("CREATE VIEW alternate_ref AS SELECT * FROM ref UNION ALL ",
                        "SELECT * REPLACE ('T' AS allele_b) FROM ref WHERE position=1"))
  expect_error(run(reference_table = "alternate_ref"), "duplicate reference locus")
  dbExecute(con, paste0("CREATE VIEW unused_invalid_ref AS SELECT * FROM ref UNION ALL ",
                        "SELECT * REPLACE (13 AS position, NULL::DOUBLE AS A) ",
                        "FROM ref WHERE position=1"))
  expect_equal(run(reference_table = "unused_invalid_ref")$proportion,
               baseline$proportion)
  dbExecute(con, paste0("CREATE VIEW unused_duplicated_ref AS SELECT * FROM ref ",
                        "UNION ALL SELECT * REPLACE (13 AS position) FROM ref ",
                        "WHERE position=1 UNION ALL SELECT * REPLACE (13 AS position) ",
                        "FROM ref WHERE position=1"))
  expect_equal(run(reference_table = "unused_duplicated_ref")$proportion,
               baseline$proportion)
  dbExecute(con, "CREATE VIEW integer_ref AS SELECT * REPLACE (1::INTEGER AS chromosome) FROM ref")
  expect_equal(run(reference_table = "integer_ref")$proportion,
               baseline$proportion)
  dbExecute(con, "CREATE VIEW missing_ref AS SELECT * REPLACE (CASE WHEN position=1 THEN NULL ELSE A END AS A) FROM ref")
  expect_error(run(reference_table = "missing_ref"), "finite frequencies and loadings")
  dbExecute(con, "CREATE VIEW missing_pc AS SELECT * REPLACE (CASE WHEN position=1 THEN NULL ELSE PC1 END AS PC1) FROM ref")
  expect_error(run(reference_table = "missing_pc"), "finite frequencies and loadings")
  dbExecute(con, "CREATE VIEW infinite_ref AS SELECT * REPLACE ('Inf'::DOUBLE AS A) FROM ref")
  expect_error(run(reference_table = "infinite_ref"), "finite frequencies and loadings")
  dbExecute(con, "CREATE VIEW null_sample AS SELECT * REPLACE (NULL::VARCHAR AS sample_id) FROM input")
  expect_error(run("null_sample"), "sample_id must not be NULL")
  dbExecute(con, paste0(
    "CREATE VIEW long_ref AS SELECT chromosome, position, allele_a, allele_b, ",
    "g.group_id, CASE g.group_id WHEN 'A' THEN A WHEN 'B' THEN B ELSE C END ",
    "AS frequency FROM ref, (VALUES ('A'), ('B'), ('C')) g(group_id)"))
  dbExecute(con, paste0(
    "CREATE VIEW long_pc AS SELECT chromosome, position, allele_a, allele_b, ",
    "pc, CASE pc WHEN 1 THEN PC1 ELSE PC2 END AS loading ",
    "FROM ref, (VALUES (1), (2)) p(pc)"))
  dbExecute(con, "CREATE VIEW correction2 AS SELECT * FROM correction WHERE pc <= 2")
  dbExecute(con, paste0("CREATE VIEW incomplete_long_ref AS SELECT * FROM ",
                        "long_ref WHERE NOT (position=1 AND group_id='A')"))
  expect_error(rduckhts_ancestry_proportions(
    con, "input", "incomplete_long_ref", "long_pc", "correction2"),
    "one finite frequency per group")
  dbExecute(con, paste0("CREATE VIEW incomplete_long_pc AS SELECT * FROM ",
                        "long_pc WHERE NOT (position=1 AND pc=1)"))
  expect_error(rduckhts_ancestry_proportions(
    con, "input", "long_ref", "incomplete_long_pc", "correction2"),
    "one finite loading per PC")
  dbExecute(con, paste0("CREATE VIEW zero_pcs AS SELECT * REPLACE (",
                        paste0("0.0::DOUBLE AS PC", seq_len(16L), collapse = ", "),
                        ") FROM ref"))
  dbBegin(con)
  error <- tryCatch(run(reference_table = "zero_pcs"), error = conditionMessage)
  expect_true(grepl("ancestry|solver|converge|singular", error, ignore.case = TRUE))
  expect_false(grepl("transaction is aborted", error, ignore.case = TRUE))
  dbRollback(con)
  expect_equal(length(list.files(tempdir(), "^ancestry_aligned_")), 0L)
  expect_equal(dbGetQuery(con, paste0("SELECT count(*) AS n FROM duckdb_tables() ",
    "WHERE temporary AND table_name LIKE 'ancestry_%'"))$n, 0)
  gate <- 0.99503379597967712
  gated <- run(min_cor = gate)
  expect_equal(unique(gated$status), "low_correlation")
  expect_true(all(is.na(gated$proportion)))
}

test_ancestry_checked_sites()
