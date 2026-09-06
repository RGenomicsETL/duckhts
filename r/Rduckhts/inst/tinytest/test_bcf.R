library(tinytest)
library(DBI)

test_bcf_string_format_lists <- function() {
  con <- rduckhts_connect()
  on.exit(
    {
      dbDisconnect(con, shutdown = TRUE)
    },
    add = TRUE
  )

  bcf_path <- system.file(
    "extdata",
    "format_string_list.vcf",
    package = "Rduckhts"
  )
  fixed_path <- system.file(
    "extdata",
    "fixed_count_arrays.vcf",
    package = "Rduckhts"
  )
  mapping_path <- system.file(
    "extdata",
    "mapping_number_families.vcf",
    package = "Rduckhts"
  )
  vep_path <- system.file("extdata", "test_vep.vcf", package = "Rduckhts")
  expect_true(file.exists(bcf_path))
  expect_true(file.exists(fixed_path))
  expect_true(file.exists(mapping_path))
  expect_true(file.exists(vep_path))

  expect_silent(rduckhts_bcf(
    con,
    "bcf_string_lists",
    bcf_path,
    decompression_threads = 0,
    overwrite = TRUE
  ))
  expect_silent(rduckhts_bcf(
    con,
    "bcf_fixed_counts",
    fixed_path,
    overwrite = TRUE
  ))
  expect_error(
    rduckhts_bcf(
      con,
      "bcf_bad_threads",
      bcf_path,
      decompression_threads = -1,
      overwrite = TRUE
    ),
    pattern = "decompression_threads must be a single whole number"
  )
  expect_error(
    rduckhts_bcf(
      con,
      "bcf_bad_csq_types",
      vep_path,
      additional_csq_column_types = "DISTANCE BogusType",
      overwrite = TRUE
    ),
    pattern = "additional_csq_column_types must use bcftools-style"
  )

  type_row <- DBI::dbGetQuery(
    con,
    "SELECT typeof(FORMAT_LAA_SAMPLE1) AS t FROM bcf_string_lists LIMIT 1"
  )
  expect_equal(type_row$t[1], "VARCHAR[]")

  first_row <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT list_extract(FORMAT_LAA_SAMPLE1, 1) AS laa1,",
      "list_extract(FORMAT_LAA_SAMPLE1, 2) AS laa2",
      "FROM bcf_string_lists WHERE POS = 73"
    )
  )
  expect_equal(first_row$laa1[1], "1")
  expect_equal(first_row$laa2[1], "2")

  missing_row <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT FORMAT_LAA_SAMPLE1 IS NULL AS is_null,",
      "FORMAT_LAA_SAMPLE1::VARCHAR AS laa",
      "FROM bcf_string_lists WHERE POS = 263"
    )
  )
  expect_true(isTRUE(missing_row$is_null[1]))
  expect_true(is.na(missing_row$laa[1]))

  fixed_type <- DBI::dbGetQuery(
    con,
    "SELECT typeof(INFO_SB) AS info_t, typeof(FORMAT_HQ_S1) AS fmt_t FROM bcf_fixed_counts LIMIT 1"
  )
  expect_equal(fixed_type$info_t[1], "INTEGER[]")
  expect_equal(fixed_type$fmt_t[1], "INTEGER[]")

  fixed_row <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT INFO_SB::VARCHAR AS info_sb,",
      "FORMAT_HQ_S1::VARCHAR AS format_hq",
      "FROM bcf_fixed_counts WHERE POS = 100"
    )
  )
  expect_equal(fixed_row$info_sb[1], "[10, 20, 30, 40]")
  expect_equal(fixed_row$format_hq[1], "[5, 9]")

  quoted_mapping <- DBI::dbQuoteString(con, mapping_path)

  number_dot_row <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT typeof(INFO_VSTR) AS info_vstr_type, INFO_VSTR::VARCHAR AS info_vstr,",
        "typeof(INFO_VINT) AS info_vint_type, INFO_VINT::VARCHAR AS info_vint,",
        "typeof(FORMAT_LAA_S1) AS fmt_laa_type, FORMAT_LAA_S1::VARCHAR AS fmt_laa",
        "FROM read_bcf(%s) WHERE POS = 200"
      ),
      quoted_mapping
    )
  )
  expect_equal(number_dot_row$info_vstr_type[1], "VARCHAR[]")
  expect_equal(number_dot_row$info_vstr[1], "[red, blue, green]")
  expect_equal(number_dot_row$info_vint_type[1], "INTEGER[]")
  expect_equal(number_dot_row$info_vint[1], "[1, 2, 3]")
  expect_equal(number_dot_row$fmt_laa_type[1], "VARCHAR[]")
  expect_equal(number_dot_row$fmt_laa[1], "[1, 2]")

}

test_bcf_projection_sql <- function() {
  con <- rduckhts_connect()
  on.exit(
    {
      dbDisconnect(con, shutdown = TRUE)
    },
    add = TRUE
  )
  bcf_path <- system.file("extdata", "formatcols.vcf.gz", package = "Rduckhts")
  expect_true(file.exists(bcf_path))
  quoted_path <- DBI::dbQuoteString(con, bcf_path)

  out <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT count(*) AS n, min(POS) AS min_pos, max(POS) AS max_pos, ",
        "min(FORMAT_S) AS fmt_s ",
        "FROM read_bcf(%s, tidy_format := true)"
      ),
      quoted_path
    )
  )
  expect_equal(out$n[1], 3)
  expect_equal(out$min_pos[1], 100)
  expect_equal(out$max_pos[1], 100)
  expect_equal(out$fmt_s[1], "a")

  expect_equal(dbGetQuery(con, paste(
    "SELECT count(*) AS n FROM duckdb_functions()",
    "WHERE function_name = 'read_bcf_v2'"
  ))$n[[1]], 0)

  vep_tidy_path <- system.file(
    "extdata",
    "test_vep_tidy.vcf",
    package = "Rduckhts"
  )
  expect_true(file.exists(vep_tidy_path))
  quoted_vep_tidy <- DBI::dbQuoteString(con, vep_tidy_path)
  vep_tidy <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT SAMPLE_ID, VEP_SYMBOL::VARCHAR AS symbol,",
        "VEP_DISTANCE::VARCHAR AS distance, FORMAT_GT",
        "FROM read_bcf(%s, tidy_format := true)",
        "ORDER BY SAMPLE_ID"
      ),
      quoted_vep_tidy
    )
  )
  expect_equal(vep_tidy$SAMPLE_ID, c("S1", "S2"))
  expect_equal(vep_tidy$symbol, c("[GENE1]", "[GENE1]"))
  expect_equal(vep_tidy$distance, c("[12]", "[12]"))
  expect_equal(vep_tidy$FORMAT_GT, c("0/1", "1/1"))

  indexed_bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
  quoted_indexed_bcf <- DBI::dbQuoteString(con, indexed_bcf_path)

  gt_edge_cases <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT string_agg(FORMAT_GT, ',' ORDER BY POS, SAMPLE_ID) AS gt ",
        "FROM read_bcf(%s, tidy_format := true) ",
        "WHERE FORMAT_GT IN ('2', '0/300', '240/260')"
      ),
      quoted_indexed_bcf
    )
  )
  expect_equal(gt_edge_cases$gt[1], "2,0/300,240/260")

  ploidy_edge_path <- system.file(
    "extdata",
    "genotype_ploidy_edge_cases.vcf",
    package = "Rduckhts"
  )
  expect_true(file.exists(ploidy_edge_path))
  ploidy_edges <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT string_agg(SAMPLE_ID || '=' || FORMAT_GT, ',' ORDER BY SAMPLE_ID) AS gt ",
        "FROM read_bcf(%s, tidy_format := true)"
      ),
      DBI::dbQuoteString(con, ploidy_edge_path)
    )
  )
  expect_equal(
    ploidy_edges$gt[1],
    "BIG=0/10,DIP=0|1,HAP=1,TET=0/1/2/3,TRI=0/1/2"
  )

  phased_fields_path <- system.file(
    "extdata",
    "phased_genotype_fields.vcf",
    package = "Rduckhts"
  )
  expect_true(file.exists(phased_fields_path))
  quoted_phased <- DBI::dbQuoteString(con, phased_fields_path)
  phased_gt <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT string_agg(SAMPLE_ID || '=' || FORMAT_GT || ':PS=' || coalesce(CAST(FORMAT_PS AS VARCHAR), 'NA'), ',' ORDER BY POS, SAMPLE_ID) AS gt ",
        "FROM read_bcf(%s, tidy_format := true)"
      ),
      quoted_phased
    )
  )
  expect_equal(
    phased_gt$gt[1],
    paste(
      "S1=0|1:PS=10,S2=1|2:PS=10,S1=1|0:PS=20,S2=0/1:PS=NA",
      "S1=1:PS=30,S2=0|1|1:PS=30,S1=0|1|1|2:PS=40,S2=2|2:PS=40",
      sep = ","
    )
  )
  phased_ploidy_lengths <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT string_agg(SAMPLE_ID || '=' || FORMAT_GT || ':PL=' || length(FORMAT_PL) || ':GP=' || length(FORMAT_GP) || ':DS=' || length(FORMAT_DS), ',' ORDER BY POS, SAMPLE_ID) AS lengths ",
        "FROM read_bcf(%s, tidy_format := true) ",
        "WHERE POS >= 30"
      ),
      quoted_phased
    )
  )
  expect_equal(
    phased_ploidy_lengths$lengths[1],
    "S1=1:PL=2:GP=2:DS=1,S2=0|1|1:PL=4:GP=4:DS=1,S1=0|1|1|2:PL=15:GP=15:DS=2,S2=2|2:PL=6:GP=6:DS=2"
  )

  tidy_chunk_path <- system.file(
    "extdata",
    "tidy_chunk_boundary.vcf",
    package = "Rduckhts"
  )
  expect_true(file.exists(tidy_chunk_path))
  tidy_chunk <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "WITH t AS (",
        "SELECT SAMPLE_ID, FORMAT_GT ",
        "FROM read_bcf(%s, tidy_format := true)",
        ") ",
        "SELECT ",
        "(SELECT count(*) FROM t) AS n, ",
        "(SELECT count(DISTINCT SAMPLE_ID) FROM t) AS samples, ",
        "(SELECT string_agg(SAMPLE_ID || '=' || FORMAT_GT, ',' ORDER BY SAMPLE_ID) ",
        " FROM t WHERE SAMPLE_ID IN ('S0001', 'S2048', 'S2049', 'S2053')) AS edge_gt"
      ),
      DBI::dbQuoteString(con, tidy_chunk_path)
    )
  )
  expect_equal(tidy_chunk$n[1], 2053)
  expect_equal(tidy_chunk$samples[1], 2053)
  expect_equal(tidy_chunk$edge_gt[1], "S0001=0/1,S2048=0/1,S2049=0/1,S2053=0/1")

  ensembl_path <- system.file(
    "extdata",
    "ensembl_release_consequences.vcf",
    package = "Rduckhts"
  )
  expect_true(file.exists(ensembl_path))
  ensembl_counts <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT count(*) AS records, ",
        "sum(list_count(VEP_Consequence)) AS csq_entries, ",
        "sum(list_count(INFO_VE)) AS ve_entries ",
        "FROM read_bcf(%s)"
      ),
      DBI::dbQuoteString(con, ensembl_path)
    )
  )
  expect_equal(ensembl_counts$records[1], 2)
  expect_equal(ensembl_counts$csq_entries[1], 6)
  expect_equal(ensembl_counts$ve_entries[1], 9)

}

test_bcf_materialization <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
  quoted <- dbQuoteString(con, path)
  expect_equal(dbGetQuery(con, paste(
    "SELECT count(*) AS n FROM duckdb_functions()",
    "WHERE function_name IN ('read_bcf_v2','read_bcf_appender')"))$n[[1]], 0)
  dbExecute(con, sprintf(paste(
    "CREATE TABLE bcf_materialized AS SELECT * FROM read_bcf(%s,",
    "region := '1:3062915-3062915,1:3062918-3062918,1:3062915-3062915,2:3199812-3199812',",
    "tidy_format := true)"), quoted))
  actual <- dbGetQuery(con, "SELECT * FROM bcf_materialized ORDER BY CHROM,POS,SAMPLE_ID")
  expected <- dbGetQuery(con, sprintf(paste(
    "SELECT * FROM read_bcf(%s, scan_mode := 'sequential', tidy_format := true)",
    "WHERE (CHROM='1' AND POS IN (3062915,3062918)) OR (CHROM='2' AND POS=3199812)",
    "ORDER BY CHROM,POS,SAMPLE_ID"), quoted))
  expect_equal(nrow(actual), 6L)
  expect_equal(actual, expected)
  dbExecute(con, sprintf(paste(
    "CREATE OR REPLACE TABLE bcf_materialized AS SELECT * FROM read_bcf(%s,",
    "region := 'absent_contig:1-1000', tidy_format := true)"), quoted))
  expect_equal(dbGetQuery(con, "SELECT count(*) AS n FROM bcf_materialized")$n[[1]], 0)
}

test_bcf_record_cache <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  path <- system.file("extdata", "bcf_cache_lifecycle.vcf", package = "Rduckhts")
  expect_true(file.exists(path))
  out <- dbGetQuery(con, sprintf(paste(
    "SELECT POS, SAMPLE_ID, VEP_SYMBOL::VARCHAR AS symbol, INFO_DP,",
    "FORMAT_GT, FORMAT_AD::VARCHAR AS ad, FORMAT_ST::VARCHAR AS st",
    "FROM read_bcf(%s, tidy_format := true) ORDER BY POS, SAMPLE_ID"
  ), dbQuoteString(con, path)))
  expect_equal(as.integer(out$POS), rep(c(10L, 20L, 30L), each = 2L))
  expect_equal(out$SAMPLE_ID, rep(c("S1", "S2"), 3L))
  expect_equal(out$symbol, rep(c("[GENE1]", NA, "[GENE2, GENE3]"), each = 2L))
  expect_equal(out$INFO_DP, rep(c(7L, NA_integer_, 1234L), each = 2L))
  expect_equal(out$FORMAT_GT, c("0/1", "1/1", "./.", "0", "0|1|2", "2/2"))
  expect_equal(out$ad, c("[6, 1]", "[0, 7]", NA, NA, "[4, 5, 6]", "[0, 0, 9]"))
  expect_equal(out$st, c("[x, y]", "[z]", NA, NA, "[longer, value, three]", "[last]"))

  path <- system.file("extdata", "tidy_chunk_boundary.vcf", package = "Rduckhts")
  out <- dbGetQuery(con, sprintf(paste(
    "SELECT count(*) AS n, count(*) FILTER (WHERE VEP_SYMBOL = ['GENE1']",
    "AND VEP_DISTANCE = [12] AND FORMAT_GT = '0/1') AS correct",
    "FROM read_bcf(%s, tidy_format := true)"
  ), dbQuoteString(con, path)))
  expect_equal(out$n[[1]], 2053)
  expect_equal(out$correct[[1]], 2053)
}

test_bcf_index_reload <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  index <- tempfile("duckhts-bcf-index-")
  hidden <- paste0(index, ".hidden")
  on.exit(unlink(c(index, hidden)), add = TRUE)
  for (file in c("parallel_empty_contigs.bcf", "parallel_empty_contigs.vcf.gz")) {
    path <- system.file("extdata", file, package = "Rduckhts")
    expect_true(nzchar(path))
    build <- sprintf("SELECT * FROM bcf_index(%s, index_path := %s, threads := 1)",
                     dbQuoteString(con, path), dbQuoteString(con, index))
    for (workers in c(1L, 4L)) for (tidy in c(FALSE, TRUE)) for (damage in c("missing", "corrupt")) {
      dbExecute(con, sprintf("SET threads=%d", workers))
      dbGetQuery(con, build)
      source <- sprintf("read_bcf(%s, index_path := %s, tidy_format := %s, decompression_threads := 0",
                        dbQuoteString(con, path), dbQuoteString(con, index), tolower(tidy))
      # R DuckDB's EXECUTE bookkeeping coerces the first column to numeric.
      # Keep POS first while retaining every column for exact row comparison.
      projection <- "SELECT POS, * EXCLUDE (POS) FROM "
      sql <- paste0(projection, source, ") ORDER BY CHROM, POS")
      expected <- dbGetQuery(con, sql)
      expect_equal(nrow(expected), 2L)
      dbExecute(con, paste("PREPARE index_scan AS", sql))
      if (damage == "missing") {
        expect_true(file.rename(index, hidden))
      } else writeLines("broken index", index)
      expect_error(dbGetQuery(con, "EXECUTE index_scan"),
                   pattern = "read_bcf: failed to reload index for parallel scan")
      expect_equal(dbGetQuery(con, "SELECT 4242 AS n")$n, 4242L)
      expect_equal(dbGetQuery(con, sql), expected) # New bind: no-index fallback.
      sequential <- paste0(projection, source, ", scan_mode := 'sequential') ORDER BY CHROM, POS")
      expect_equal(dbGetQuery(con, sequential), expected)
      dbGetQuery(con, build)
      expect_equal(dbGetQuery(con, "EXECUTE index_scan"), expected)
      dbExecute(con, "DEALLOCATE index_scan")
      if (file.exists(hidden)) unlink(hidden)
    }
  }
}

test_bcf_index_reload()
test_bcf_record_cache()
test_bcf_string_format_lists()
test_bcf_projection_sql()
test_bcf_materialization()
