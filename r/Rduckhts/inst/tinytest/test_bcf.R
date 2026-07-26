library(tinytest)
library(DBI)

test_bcf_string_format_lists <- function() {
  con <- rduckhts_connect()
  on.exit({
    dbDisconnect(con, shutdown = TRUE)
  }, add = TRUE)

  bcf_path <- system.file("extdata", "format_string_list.vcf", package = "Rduckhts")
  fixed_path <- system.file("extdata", "fixed_count_arrays.vcf", package = "Rduckhts")
  mapping_path <- system.file("extdata", "mapping_number_families.vcf", package = "Rduckhts")
  vep_path <- system.file("extdata", "test_vep.vcf", package = "Rduckhts")
  expect_true(file.exists(bcf_path))
  expect_true(file.exists(fixed_path))
  expect_true(file.exists(mapping_path))
  expect_true(file.exists(vep_path))

  expect_silent(rduckhts_bcf(con, "bcf_string_lists", bcf_path, decompression_threads = 0, overwrite = TRUE))
  expect_silent(rduckhts_bcf(con, "bcf_fixed_counts", fixed_path, overwrite = TRUE))
  expect_error(
    rduckhts_bcf(con, "bcf_bad_threads", bcf_path, decompression_threads = -1, overwrite = TRUE),
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

  quoted_bcf <- DBI::dbQuoteString(con, bcf_path)
  quoted_fixed <- DBI::dbQuoteString(con, fixed_path)
  quoted_mapping <- DBI::dbQuoteString(con, mapping_path)

  v2_type_row <- DBI::dbGetQuery(
    con,
    sprintf("SELECT typeof(FORMAT_LAA_SAMPLE1) AS t FROM read_bcf_v2(%s) LIMIT 1", quoted_bcf)
  )
  expect_equal(v2_type_row$t[1], "VARCHAR[]")

  v2_missing_row <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT FORMAT_LAA_SAMPLE1 IS NULL AS is_null,",
        "FORMAT_LAA_SAMPLE1::VARCHAR AS laa",
        "FROM read_bcf_v2(%s) WHERE POS = 263"
      ),
      quoted_bcf
    )
  )
  expect_true(isTRUE(v2_missing_row$is_null[1]))
  expect_true(is.na(v2_missing_row$laa[1]))

  v2_fixed_row <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT INFO_SB::VARCHAR AS info_sb,",
        "FORMAT_HQ_S1::VARCHAR AS format_hq",
        "FROM read_bcf_v2(%s) WHERE POS = 100"
      ),
      quoted_fixed
    )
  )
  expect_equal(v2_fixed_row$info_sb[1], "[10, 20, 30, 40]")
  expect_equal(v2_fixed_row$format_hq[1], "[5, 9]")

  v2_number_dot_row <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT typeof(INFO_VSTR) AS info_vstr_type, INFO_VSTR::VARCHAR AS info_vstr,",
        "typeof(INFO_VINT) AS info_vint_type, INFO_VINT::VARCHAR AS info_vint,",
        "typeof(FORMAT_LAA_S1) AS fmt_laa_type, FORMAT_LAA_S1::VARCHAR AS fmt_laa",
        "FROM read_bcf_v2(%s) WHERE POS = 200"
      ),
      quoted_mapping
    )
  )
  expect_equal(v2_number_dot_row$info_vstr_type[1], "VARCHAR[]")
  expect_equal(v2_number_dot_row$info_vstr[1], "[red, blue, green]")
  expect_equal(v2_number_dot_row$info_vint_type[1], "INTEGER[]")
  expect_equal(v2_number_dot_row$info_vint[1], "[1, 2, 3]")
  expect_equal(v2_number_dot_row$fmt_laa_type[1], "VARCHAR[]")
  expect_equal(v2_number_dot_row$fmt_laa[1], "[1, 2]")

  mapping_mismatch <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "WITH ",
        "v1 AS (SELECT * FROM read_bcf(%s)), ",
        "v2 AS (SELECT * FROM read_bcf_v2(%s)) ",
        "SELECT count(*) AS n FROM ((SELECT * FROM v1 EXCEPT ALL SELECT * FROM v2) ",
        "UNION ALL (SELECT * FROM v2 EXCEPT ALL SELECT * FROM v1))"
      ),
      quoted_mapping, quoted_mapping
    )
  )
  expect_equal(mapping_mismatch$n[1], 0)
}

test_bcf_v2_sql <- function() {
  con <- rduckhts_connect()
  on.exit({
    dbDisconnect(con, shutdown = TRUE)
  }, add = TRUE)
  bcf_path <- system.file("extdata", "formatcols.vcf.gz", package = "Rduckhts")
  expect_true(file.exists(bcf_path))
  quoted_path <- DBI::dbQuoteString(con, bcf_path)

  out <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT count(*) AS n, min(POS) AS min_pos, max(POS) AS max_pos, ",
        "min(FORMAT_S) AS fmt_s ",
        "FROM read_bcf_v2(%s, tidy_format := true)"
      ),
      quoted_path
    )
  )
  expect_equal(out$n[1], 3)
  expect_equal(out$min_pos[1], 100)
  expect_equal(out$max_pos[1], 100)
  expect_equal(out$fmt_s[1], "a")

  schema_mismatch <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "WITH ",
        "v1 AS (SELECT row_number() OVER () AS ord, column_name, column_type ",
        "       FROM (DESCRIBE SELECT * FROM read_bcf(%s))), ",
        "v2 AS (SELECT row_number() OVER () AS ord, column_name, column_type ",
        "       FROM (DESCRIBE SELECT * FROM read_bcf_v2(%s))) ",
        "SELECT count(*) AS n FROM ((SELECT * FROM v1 EXCEPT ALL SELECT * FROM v2) ",
        "UNION ALL (SELECT * FROM v2 EXCEPT ALL SELECT * FROM v1))"
      ),
      quoted_path, quoted_path
    )
  )
  expect_equal(schema_mismatch$n[1], 0)

  fixed_path <- system.file("extdata", "fixed_count_arrays.vcf", package = "Rduckhts")
  expect_true(file.exists(fixed_path))
  quoted_fixed <- DBI::dbQuoteString(con, fixed_path)

  default_mismatch <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "WITH ",
        "v1 AS (SELECT * FROM read_bcf(%s)), ",
        "v2 AS (SELECT * FROM read_bcf_v2(%s)) ",
        "SELECT count(*) AS n FROM ((SELECT * FROM v1 EXCEPT ALL SELECT * FROM v2) ",
        "UNION ALL (SELECT * FROM v2 EXCEPT ALL SELECT * FROM v1))"
      ),
      quoted_fixed, quoted_fixed
    )
  )
  expect_equal(default_mismatch$n[1], 0)

  tidy_mismatch <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "WITH ",
        "v1 AS (SELECT * FROM read_bcf(%s, tidy_format := true)), ",
        "v2 AS (SELECT * FROM read_bcf_v2(%s, tidy_format := true)) ",
        "SELECT count(*) AS n FROM ((SELECT * FROM v1 EXCEPT ALL SELECT * FROM v2) ",
        "UNION ALL (SELECT * FROM v2 EXCEPT ALL SELECT * FROM v1))"
      ),
      quoted_fixed, quoted_fixed
    )
  )
  expect_equal(tidy_mismatch$n[1], 0)

  vep_tidy_path <- system.file("extdata", "test_vep_tidy.vcf", package = "Rduckhts")
  expect_true(file.exists(vep_tidy_path))
  quoted_vep_tidy <- DBI::dbQuoteString(con, vep_tidy_path)
  v2_vep_tidy <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT SAMPLE_ID, VEP_SYMBOL::VARCHAR AS symbol,",
        "VEP_DISTANCE::VARCHAR AS distance, FORMAT_GT",
        "FROM read_bcf_v2(%s, tidy_format := true)",
        "ORDER BY SAMPLE_ID"
      ),
      quoted_vep_tidy
    )
  )
  expect_equal(v2_vep_tidy$SAMPLE_ID, c("S1", "S2"))
  expect_equal(v2_vep_tidy$symbol, c("[GENE1]", "[GENE1]"))
  expect_equal(v2_vep_tidy$distance, c("[12]", "[12]"))
  expect_equal(v2_vep_tidy$FORMAT_GT, c("0/1", "1/1"))

  samples_file <- system.file("extdata", "samples_s1.txt", package = "Rduckhts")
  expect_true(file.exists(samples_file))
  quoted_samples_file <- DBI::dbQuoteString(con, samples_file)

  sample_schema <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT count(*) AS n FROM ",
        "(DESCRIBE SELECT * FROM read_bcf_v2(%s, samples := 'S1')) ",
        "WHERE column_name = 'FORMAT_GT_S2'"
      ),
      quoted_fixed
    )
  )
  expect_equal(sample_schema$n[1], 0)

  sample_pushdown <- DBI::dbGetQuery(
    con,
    sprintf(
      paste(
        "SELECT count(*) AS n, count(DISTINCT SAMPLE_ID) AS samples,",
        "string_agg(FORMAT_GT, ',' ORDER BY POS) AS gt",
        "FROM read_bcf_v2(%s, samples_file := %s, tidy_format := true)"
      ),
      quoted_fixed, quoted_samples_file
    )
  )
  expect_equal(sample_pushdown$n[1], 2)
  expect_equal(sample_pushdown$samples[1], 1)
  expect_equal(sample_pushdown$gt[1], "0/1,0/0")

  indexed_bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
  expect_true(file.exists(indexed_bcf_path))
  quoted_indexed_bcf <- DBI::dbQuoteString(con, indexed_bcf_path)
  zero_sample_tidy <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT ",
        "(SELECT count(*) FROM read_bcf_v2(%s, tidy_format := true, exclude_samples := 'A,B')) AS count_only, ",
        "(SELECT count(CHROM) FROM read_bcf_v2(%s, tidy_format := true, exclude_samples := 'A,B')) AS materialized, ",
        "(SELECT count(*) FROM (DESCRIBE SELECT * FROM read_bcf_v2(%s, tidy_format := true, exclude_samples := 'A,B')) ",
        " WHERE column_name = 'SAMPLE_ID') AS sample_id_cols"
      ),
      quoted_indexed_bcf, quoted_indexed_bcf, quoted_indexed_bcf
    )
  )
  expect_equal(zero_sample_tidy$count_only[1], zero_sample_tidy$materialized[1])
  expect_equal(zero_sample_tidy$sample_id_cols[1], 0)

  minimal_schema <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT count(*) AS n, string_agg(column_name, ',' ORDER BY column_name) AS cols ",
        "FROM (DESCRIBE SELECT * FROM read_bcf_v2(%s, ",
        "include_info := false, include_format := false, include_vep := false))"
      ),
      quoted_fixed
    )
  )
  expect_equal(minimal_schema$n[1], 7)
  expect_equal(minimal_schema$cols[1], "ALT,CHROM,FILTER,ID,POS,QUAL,REF")

  filtered_schema <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT string_agg(column_name, ',' ORDER BY column_name) AS cols ",
        "FROM (DESCRIBE SELECT * FROM read_bcf_v2(%s, info_fields := 'SB', format_fields := 'GT'))"
      ),
      quoted_fixed
    )
  )
  expect_equal(
    filtered_schema$cols[1],
    "ALT,CHROM,FILTER,FORMAT_GT_S1,FORMAT_GT_S2,ID,INFO_SB,POS,QUAL,REF"
  )

  filtered_rows <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT count(*) AS n, sum(list_count(INFO_SB)) AS sb, ",
        "sum(length(coalesce(FORMAT_GT_S1, '')) + length(coalesce(FORMAT_GT_S2, ''))) AS gt_len ",
        "FROM read_bcf_v2(%s, info_fields := 'SB', format_fields := 'GT')"
      ),
      quoted_fixed
    )
  )
  expect_equal(filtered_rows$n[1], 2)
  expect_equal(filtered_rows$sb[1], 4)
  expect_equal(filtered_rows$gt_len[1], 12)

  gt_projection_mismatch <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "WITH ",
        "v1 AS (SELECT POS, SAMPLE_ID, FORMAT_GT FROM read_bcf(%s, tidy_format := true)), ",
        "v2 AS (SELECT POS, SAMPLE_ID, FORMAT_GT FROM read_bcf_v2(%s, tidy_format := true, format_fields := 'GT')) ",
        "SELECT count(*) AS n FROM ((SELECT * FROM v1 EXCEPT ALL SELECT * FROM v2) ",
        "UNION ALL (SELECT * FROM v2 EXCEPT ALL SELECT * FROM v1))"
      ),
      quoted_indexed_bcf, quoted_indexed_bcf
    )
  )
  expect_equal(gt_projection_mismatch$n[1], 0)

  gt_edge_cases <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT string_agg(FORMAT_GT, ',' ORDER BY POS, SAMPLE_ID) AS gt ",
        "FROM read_bcf_v2(%s, tidy_format := true, format_fields := 'GT') ",
        "WHERE FORMAT_GT IN ('2', '0/300', '240/260')"
      ),
      quoted_indexed_bcf
    )
  )
  expect_equal(gt_edge_cases$gt[1], "2,0/300,240/260")

  ploidy_edge_path <- system.file("extdata", "genotype_ploidy_edge_cases.vcf", package = "Rduckhts")
  expect_true(file.exists(ploidy_edge_path))
  ploidy_edges <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT string_agg(SAMPLE_ID || '=' || FORMAT_GT, ',' ORDER BY SAMPLE_ID) AS gt ",
        "FROM read_bcf_v2(%s, tidy_format := true, format_fields := 'GT')"
      ),
      DBI::dbQuoteString(con, ploidy_edge_path)
    )
  )
  expect_equal(ploidy_edges$gt[1], "BIG=0/10,DIP=0|1,HAP=1,TET=0/1/2/3,TRI=0/1/2")

  phased_fields_path <- system.file("extdata", "phased_genotype_fields.vcf", package = "Rduckhts")
  expect_true(file.exists(phased_fields_path))
  quoted_phased <- DBI::dbQuoteString(con, phased_fields_path)
  phased_mismatch <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "WITH ",
        "v1 AS (SELECT POS, SAMPLE_ID, FORMAT_GT, FORMAT_PL, FORMAT_GP, FORMAT_DS, FORMAT_PS ",
        "       FROM read_bcf(%s, tidy_format := true)), ",
        "v2 AS (SELECT POS, SAMPLE_ID, FORMAT_GT, FORMAT_PL, FORMAT_GP, FORMAT_DS, FORMAT_PS ",
        "       FROM read_bcf_v2(%s, tidy_format := true, format_fields := 'GT,PL,GP,DS,PS')) ",
        "SELECT count(*) AS n FROM ((SELECT * FROM v1 EXCEPT ALL SELECT * FROM v2) ",
        "UNION ALL (SELECT * FROM v2 EXCEPT ALL SELECT * FROM v1))"
      ),
      quoted_phased, quoted_phased
    )
  )
  expect_equal(phased_mismatch$n[1], 0)
  phased_gt <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT string_agg(SAMPLE_ID || '=' || FORMAT_GT || ':PS=' || coalesce(CAST(FORMAT_PS AS VARCHAR), 'NA'), ',' ORDER BY POS, SAMPLE_ID) AS gt ",
        "FROM read_bcf_v2(%s, tidy_format := true, format_fields := 'GT,PL,GP,DS,PS')"
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
        "FROM read_bcf_v2(%s, tidy_format := true, format_fields := 'GT,PL,GP,DS,PS') ",
        "WHERE POS >= 30"
      ),
      quoted_phased
    )
  )
  expect_equal(
    phased_ploidy_lengths$lengths[1],
    "S1=1:PL=2:GP=2:DS=1,S2=0|1|1:PL=4:GP=4:DS=1,S1=0|1|1|2:PL=15:GP=15:DS=2,S2=2|2:PL=6:GP=6:DS=2"
  )

  tidy_chunk_path <- system.file("extdata", "tidy_chunk_boundary.vcf", package = "Rduckhts")
  expect_true(file.exists(tidy_chunk_path))
  tidy_chunk <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "WITH t AS (",
        "SELECT SAMPLE_ID, FORMAT_GT ",
        "FROM read_bcf_v2(%s, tidy_format := true, format_fields := 'GT')",
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

  no_format_tidy <- DBI::dbGetQuery(
    con,
    sprintf("SELECT count(*) AS n FROM read_bcf_v2(%s, tidy_format := true, include_format := false)", quoted_fixed)
  )
  expect_equal(no_format_tidy$n[1], 2)

  vep_path <- system.file("extdata", "test_vep.vcf", package = "Rduckhts")
  expect_true(file.exists(vep_path))
  quoted_vep <- DBI::dbQuoteString(con, vep_path)
  vep_schema <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT string_agg(column_name, ',' ORDER BY column_name) AS cols ",
        "FROM (DESCRIBE SELECT * FROM read_bcf_v2(%s, ",
        "include_info := false, include_format := false, vep_fields := 'Consequence,DISTANCE'))"
      ),
      quoted_vep
    )
  )
  expect_equal(vep_schema$cols[1], "ALT,CHROM,FILTER,ID,POS,QUAL,REF,VEP_Consequence,VEP_DISTANCE")

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
        "FROM read_bcf_v2(%s, info_fields := 'VE', include_format := false, ",
        "vep_fields := 'Allele,Consequence,Feature_type,Feature,Amino_acids,SIFT')"
      ),
      DBI::dbQuoteString(con, ensembl_path)
    )
  )
  expect_equal(ensembl_counts$records[1], 2)
  expect_equal(ensembl_counts$csq_entries[1], 6)
  expect_equal(ensembl_counts$ve_entries[1], 9)

  expect_error(
    DBI::dbGetQuery(
      con,
      sprintf("SELECT count(*) FROM read_bcf_v2(%s, info_fields := 'NOPE')", quoted_fixed)
    ),
    pattern = "unknown field"
  )

  score_gwas_path <- system.file("extdata", "score_gwas_summary.vcf", package = "Rduckhts")
  expect_true(file.exists(score_gwas_path))
  expect_error(
    DBI::dbGetQuery(
      con,
      sprintf(
        "SELECT count(*) FROM read_bcf_v2(%s, format_fields := 'GT')",
        DBI::dbQuoteString(con, score_gwas_path)
      )
    ),
    pattern = "unknown field"
  )

  expect_error(
    DBI::dbGetQuery(
      con,
      sprintf("SELECT count(*) FROM read_bcf_v2(%s, format_fields := 'GT')", quoted_vep)
    ),
    pattern = "no samples"
  )

  expect_error(
    DBI::dbGetQuery(
      con,
      sprintf(
        "SELECT count(*) FROM read_bcf_v2(%s, exclude_samples := 'S1,S2', format_fields := 'GT')",
        quoted_fixed
      )
    ),
    pattern = "no samples"
  )

  expect_error(
    DBI::dbGetQuery(
      con,
      sprintf("SELECT count(*) FROM read_bcf_v2(%s, samples := 'NOPE')", quoted_fixed)
    ),
    pattern = "unknown sample"
  )
}

test_bcf_appender <- function() {
  con <- rduckhts_connect()
  on.exit({
    dbDisconnect(con, shutdown = TRUE)
  }, add = TRUE)
  bcf_path <- system.file("extdata", "vcf_file.bcf", package = "Rduckhts")
  expect_true(file.exists(bcf_path))

  quoted_path <- DBI::dbQuoteString(con, bcf_path)
  out <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT rows_written FROM read_bcf_appender(",
        "%s, 'bcf_appender_test', ",
        "region := '1:3000150-3000151', tidy_format := true, overwrite := true)"
      ),
      quoted_path
    )
  )
  expect_equal(out$rows_written[1], 4)

  schema <- DBI::dbGetQuery(con, "PRAGMA table_info('bcf_appender_test')")
  expect_false("FILE_OFFSET" %in% schema$name)

  counts <- DBI::dbGetQuery(
    con,
    "SELECT count(*) AS n, count(DISTINCT SAMPLE_ID) AS samples FROM bcf_appender_test"
  )
  expect_equal(counts$n[1], 4)
  expect_equal(counts$samples[1], 2)

  out_offset <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT rows_written FROM read_bcf_appender(",
        "%s, 'bcf_appender_offset_test', ",
        "region := '1:3000150-3000151', tidy_format := true, ",
        "overwrite := true, include_file_offset := true)"
      ),
      quoted_path
    )
  )
  expect_equal(out_offset$rows_written[1], 4)

  offset_counts <- DBI::dbGetQuery(
    con,
    paste0(
      "SELECT count(*) AS n, count(DISTINCT SAMPLE_ID) AS samples, ",
      "count(DISTINCT FILE_OFFSET) AS offsets FROM bcf_appender_offset_test"
    )
  )
  expect_equal(offset_counts$n[1], 4)
  expect_equal(offset_counts$samples[1], 2)
  expect_equal(offset_counts$offsets[1], 2)

  region_spec <- paste0(
    "1:3062915-3062915,1:3062918-3062918,",
    "1:3062915-3062915,2:3199812-3199812"
  )
  out_serial_regions <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT rows_written FROM read_bcf_appender(",
        "%s, 'bcf_appender_serial_regions_test', ",
        "region := %s, tidy_format := true, overwrite := true, ",
        "include_file_offset := true, region_threads := 1)"
      ),
      quoted_path,
      DBI::dbQuoteString(con, region_spec)
    )
  )
  expect_equal(out_serial_regions$rows_written[1], 6)

  out_parallel <- DBI::dbGetQuery(
    con,
    sprintf(
      paste0(
        "SELECT rows_written FROM read_bcf_appender(",
        "%s, 'bcf_appender_parallel_regions_test', ",
        "region := %s, tidy_format := true, overwrite := true, ",
        "include_file_offset := true, region_threads := 2)"
      ),
      quoted_path,
      DBI::dbQuoteString(con, region_spec)
    )
  )
  expect_equal(out_parallel$rows_written[1], 6)

  parallel_schema <- DBI::dbGetQuery(
    con,
    "PRAGMA table_info('bcf_appender_parallel_regions_test')"
  )
  expect_false("duckhts_region_idx" %in% parallel_schema$name)

  delta <- DBI::dbGetQuery(
    con,
    paste0(
      "SELECT count(*) AS n FROM (",
      "(SELECT * FROM bcf_appender_serial_regions_test EXCEPT ALL ",
      " SELECT * FROM bcf_appender_parallel_regions_test) UNION ALL ",
      "(SELECT * FROM bcf_appender_parallel_regions_test EXCEPT ALL ",
      " SELECT * FROM bcf_appender_serial_regions_test))"
    )
  )
  expect_equal(delta$n[1], 0)
}

test_bcf_string_format_lists()
test_bcf_v2_sql()
test_bcf_appender()
