library(tinytest)
library(DBI)

test_somalier_find_sites_selection <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  dbExecute(con, paste(
    "CREATE TEMP TABLE population AS SELECT * FROM (VALUES",
    "('chr1',100,'A',['G'],['PASS'],[0.48::FLOAT],120000),",
    "('chr1',110,'G',['A'],['PASS'],[0.52::FLOAT],120000),",
    "('chr1',200,'G',['A'],['PASS'],[0.25::FLOAT],120000),",
    "('chr1',230,'A',['G'],['PASS'],[0.15::FLOAT],120000),",
    "('chr1',300,'G',['A'],['PASS'],[0.85::FLOAT],120000),",
    "('chr1',400,'A',['G'],['PASS'],[NULL::FLOAT],120000),",
    "('chr1',500,'A',['G'],['PASS'],[0.48::FLOAT],99),",
    "('chr1',600,'A',['G'],['q10'],[0.48::FLOAT],120000),",
    "('chr1',700,'C',['T'],['PASS'],[0.48::FLOAT],120000),",
    "('chrX',2781480,'A',['G'],['PASS'],[0.48::FLOAT],99),",
    "('chrX',2781479,'A',['G'],['PASS'],[0.48::FLOAT],99),",
    "('chrY',900,'A',['G'],['PASS'],[0.05::FLOAT],99)",
    ") v(CHROM,POS,REF,ALT,FILTER,INFO_AF,INFO_AN)"
  ))
  dbExecute(con, "CREATE TEMP TABLE region_exclude(chrom VARCHAR, start BIGINT, stop BIGINT)")
  dbExecute(con, "INSERT INTO region_exclude VALUES ('chr1', 299, 300)")
  dbExecute(con, "CREATE TEMP TABLE region_include(chrom VARCHAR, start BIGINT, stop BIGINT)")
  dbExecute(con, paste(
    "INSERT INTO region_include VALUES ('chr1', 99, 400),",
    "('chrX', 2781479, 2781480), ('chrY', 899, 900)"
  ))
  selected <- rduckhts_somalier_find_sites(
    con, source_table = "population", assembly = "GRCh38", snp_dist = 20,
    exclude_table = "region_exclude", include_table = "region_include"
  )
  expect_equal(selected$region, c("chr1", "chr1", "chr1", "chrX", "chrY"))
  expect_equal(selected$position, c(100, 200, 230, 2781480, 900))
  expect_equal(selected$site_index, 0:4)
  expect_equal(selected$allele_a[2], "A")
  expect_equal(selected$allele_b[2], "G")
  expect_equal(round(selected$population_b_af[2], 2), 0.75)
  expect_equal(selected$source_an[4], 99)
  expect_equal(selected$source_alt_af[1], 0.48, tolerance = 1e-6)
  limited <- rduckhts_somalier_find_sites(
    con, source_table = "population", assembly = "GRCh38", snp_dist = 20,
    exclude_table = "region_exclude", include_table = "region_include",
    max_autosomal = 1, max_x = 1, max_y = 1
  )
  expect_equal(limited$position, c(100, 2781480, 900))
  dbExecute(con, paste(
    "CREATE TEMP VIEW custom_fields AS SELECT *, INFO_AF AS INFO_POP_AF,",
    "INFO_AN AS INFO_POP_AN FROM population"
  ))
  custom <- rduckhts_somalier_find_sites(
    con, source_table = "custom_fields", assembly = "GRCh38", snp_dist = 20,
    af_field = "POP_AF", an_field = "POP_AN",
    exclude_table = "region_exclude", include_table = "region_include"
  )
  expect_equal(custom[, names(custom) != "source_path"],
               selected[, names(selected) != "source_path"])

  dbExecute(con, "CREATE TEMP TABLE reversed AS SELECT * FROM population ORDER BY POS DESC")
  dbExecute(con, "SET threads=1")
  first <- rduckhts_somalier_find_sites(
    con, source_table = "reversed", assembly = "GRCh38", snp_dist = 20,
    exclude_table = "region_exclude", include_table = "region_include",
    tie_order = "lexical"
  )
  dbExecute(con, "SET threads=4")
  second <- rduckhts_somalier_find_sites(
    con, source_table = "population", assembly = "GRCh38", snp_dist = 20,
    exclude_table = "region_exclude", include_table = "region_include",
    tie_order = "lexical"
  )
  expect_equal(first[, names(first) != "source_path"],
               second[, names(second) != "source_path"])
  expect_error(rduckhts_somalier_find_sites(
    con, source_table = "population", assembly = "GRCh38", min_af = -1
  ), pattern = "min_af")
}

test_somalier_find_sites_ties_and_sex_spacing <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  dbExecute(con, "SET threads=1")
  dbExecute(con, paste(
    "CREATE TEMP TABLE population AS SELECT * FROM (VALUES",
    "('chr1',200,'A',['G'],['PASS'],[0.48::FLOAT],120000),",
    "('chr1',100,'A',['G'],['PASS'],[0.48::FLOAT],120000),",
    "('chrX',2781500,'A',['G'],['PASS'],[0.48::FLOAT],120000),",
    "('chrX',2781600,'A',['G'],['PASS'],[0.48::FLOAT],120000),",
    "('chrY',1000,'A',['G'],['PASS'],[0.48::FLOAT],120000),",
    "('chrY',1100,'A',['G'],['PASS'],[0.48::FLOAT],120000)",
    ") v(CHROM,POS,REF,ALT,FILTER,INFO_AF,INFO_AN)"
  ))
  input <- rduckhts_somalier_find_sites(
    con, source_table = "population", assembly = "GRCh38", snp_dist = 1000
  )
  expect_equal(input$position, c(200, 2781500, 2781600, 1000, 1100))
  lexical <- rduckhts_somalier_find_sites(
    con, source_table = "population", assembly = "GRCh38", snp_dist = 1000,
    tie_order = "lexical", sex_spacing = "enforced"
  )
  expect_equal(lexical$position, c(100, 2781500, 1000))
}

test_somalier_find_sites_nearby <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  dbExecute(con, paste(
    "CREATE TEMP TABLE population AS SELECT * FROM (VALUES",
    "('chr1',100,'A',['G'],['PASS'],[0.48::FLOAT],120000),",
    "('chr1',102,'G',['A'],['PASS'],[0.48::FLOAT],120000),",
    "('chr1',200,'A',['G'],['PASS'],[0.48::FLOAT],120000),",
    "('chr1',208,'A',['AT'],['q10'],[0.03::FLOAT],120000),",
    "('chr1',400,'A',['G'],['PASS'],[0.48::FLOAT],120000),",
    "('chr1',430,'G',['A'],['PASS'],[0.48::FLOAT],120000)",
    ") v(CHROM,POS,REF,ALT,FILTER,INFO_AF,INFO_AN)"
  ))
  selected <- rduckhts_somalier_find_sites(
    con, source_table = "population", assembly = "GRCh38", snp_dist = 40
  )
  expect_equal(selected$position, 400)
  dbExecute(con, "CREATE TEMP TABLE gno(chrom VARCHAR, pos BIGINT, ref VARCHAR, alt VARCHAR)")
  dbExecute(con, "INSERT INTO gno VALUES ('chr1', 400, 'A', 'G')")
  excluded <- rduckhts_somalier_find_sites(
    con, source_table = "population", assembly = "GRCh38", snp_dist = 40,
    gnotate_exclude_table = "gno"
  )
  expect_equal(nrow(excluded), 1L)
  expect_equal(excluded$position, 430)
}

test_somalier_find_sites_sorted_neighbors <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  dbExecute(con, paste(
    "CREATE TEMP TABLE population AS SELECT * FROM (VALUES",
    "('chr1',100,'A',['G'],['PASS'],[0.48::FLOAT],120000),",
    "('chr1',100,'A',['G'],['PASS'],[0.48::FLOAT],120000),",
    "('chr1',160,repeat('A',50),['A'],['q10'],[0.03::FLOAT],120000),",
    "('chr1',190,'AT',['A'],['q10'],[0.03::FLOAT],120000),",
    "('chr1',200,'A',['G'],['PASS'],[0.48::FLOAT],120000),",
    "('chr1',300,'A',['G'],['PASS'],[0.48::FLOAT],120000),",
    "('chr1',302,'G',['A'],['PASS'],[0.03::FLOAT],120000),",
    "('chr1',400,'A',['G'],['PASS'],[0.48::FLOAT],120000)",
    ") v(CHROM,POS,REF,ALT,FILTER,INFO_AF,INFO_AN)"
  ))
  selected <- rduckhts_somalier_find_sites(
    con, source_table = "population", assembly = "GRCh38"
  )
  expect_equal(selected$position, 400)
  counts <- dbGetQuery(con, Rduckhts:::.somalier_find_sites_query(
    con, "population", "GRCh38", 0.15, 115000, "AF", "AN", 10000, 0.48,
    list(include = NULL, exclude = NULL, gnotate = NULL),
    "somalier_v0.3.4", 65535, 10001, 5001, diagnostics = TRUE
  ))
  expect_equal(counts$records[counts$gate == "indels"], 2)
  expect_equal(counts$records[counts$gate == "snps"], 6)
  expect_equal(counts$records[counts$gate == "indel_clear"], 4)
  expect_equal(counts$records[counts$gate == "neighbors"], 1)
}

test_somalier_find_sites_sparse_events <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  dbExecute(con, paste(
    "CREATE TEMP TABLE population AS SELECT * FROM (VALUES",
    "('chr1',100,'C',['T'],['PASS'],[0.48::FLOAT],120000),",
    "('chr1',200,'A',['G'],['q10'],[0.48::FLOAT],120000),",
    "('chr1',300,'A',['G'],['PASS'],[0.0::FLOAT],120000)",
    ") v(CHROM,POS,REF,ALT,FILTER,INFO_AF,INFO_AN)"
  ))
  selected <- rduckhts_somalier_find_sites(
    con, source_table = "population", assembly = "GRCh38", min_af = 0
  )
  expect_equal(selected$position, 300)
  counts <- dbGetQuery(con, Rduckhts:::.somalier_find_sites_query(
    con, "population", "GRCh38", 0, 115000, "AF", "AN", 10000, 0.48,
    list(include = NULL, exclude = NULL, gnotate = NULL),
    "somalier_v0.3.4", 65535, 10001, 5001, diagnostics = TRUE
  ))
  expect_equal(counts$records[counts$gate == "src"], 3)
  expect_equal(counts$records[counts$gate == "eligible"], 2)
  expect_equal(counts$records[counts$gate == "snps"], 0)
  expect_equal(counts$records[counts$gate == "selected"], 1)
  dbExecute(con, "DELETE FROM population WHERE REF != 'C'")
  empty <- rduckhts_somalier_find_sites(
    con, source_table = "population", assembly = "GRCh38"
  )
  expect_equal(nrow(empty), 0L)
  counts <- dbGetQuery(con, Rduckhts:::.somalier_find_sites_query(
    con, "population", "GRCh38", 0.15, 115000, "AF", "AN", 10000, 0.48,
    list(include = NULL, exclude = NULL, gnotate = NULL),
    "somalier_v0.3.4", 65535, 10001, 5001, diagnostics = TRUE
  ))
  expect_equal(counts$records[counts$gate == "src"], 1)
  expect_equal(counts$records[counts$gate == "eligible"], 0)
}

test_somalier_find_sites_flags_and_endpoints <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  dbExecute(con, paste(
    "CREATE TEMP TABLE annotated AS SELECT * FROM (VALUES",
    "('chr2',100,'A',['G'],['PASS'],[0.16::FLOAT],100,'PASS',false,false,false,false,0.0,20.0,60.0),",
    "('chr2',200,'G',['A'],['PASS'],[0.84::FLOAT],99,'PASS',false,false,false,false,0.0,20.0,60.0),",
    "('chr2',300,'G',['A'],['PASS'],[NULL::FLOAT],NULL::BIGINT,'PASS',false,false,false,false,0.0,20.0,60.0),",
    "('chr2',400,'A',['G'],['PASS'],[0.48::FLOAT],100,'FAIL',false,false,false,false,0.0,20.0,60.0),",
    "('chr2',500,'A',['G'],['PASS'],[0.48::FLOAT],100,'PASS',true,false,false,false,0.0,20.0,60.0),",
    "('chr2',600,'A',['G'],['PASS'],[0.48::FLOAT],100,'PASS',false,true,false,false,0.0,20.0,60.0),",
    "('chr2',700,'A',['G'],['PASS'],[0.48::FLOAT],100,'PASS',false,false,true,false,0.0,20.0,60.0),",
    "('chr2',800,'A',['G'],['PASS'],[0.48::FLOAT],100,'PASS',false,false,false,false,2.5,20.0,60.0),",
    "('chr2',900,'A',['G'],['PASS'],[0.48::FLOAT],100,'PASS',false,false,false,false,0.0,11.0,60.0),",
    "('chr2',1000,'A',['G'],['PASS'],[0.48::FLOAT],100,'PASS',false,false,false,false,0.0,20.0,49.0),",
    "('chrX',2781479,'A',['G'],['PASS'],[0.48::FLOAT],100,'PASS',false,false,false,false,0.0,20.0,60.0),",
    "('chrX',2781480,'A',['G'],['PASS'],[0.48::FLOAT],100,'PASS',false,false,false,false,0.0,20.0,60.0),",
    "('chrX',154931045,'A',['G'],['PASS'],[0.48::FLOAT],100,'PASS',false,false,false,false,0.0,20.0,60.0),",
    "('chrX',154931046,'A',['G'],['PASS'],[0.48::FLOAT],100,'PASS',false,false,false,false,0.0,20.0,60.0),",
    "('chrY',1200,'A',['G'],['PASS'],[0.48::FLOAT],100,'PASS',false,false,true,false,0.0,20.0,60.0)",
    ") v(CHROM,POS,REF,ALT,FILTER,INFO_AF,INFO_AN,INFO_AS_FilterStatus,",
    "INFO_OLD_MULTIALLELIC,INFO_OLD_VARIANT,INFO_segdup,INFO_lcr,",
    "INFO_BaseQRankSum,INFO_QD,INFO_MQ)"
  ))
  dbExecute(con, "ALTER TABLE annotated ADD COLUMN INFO_FS DOUBLE")
  dbExecute(con, paste(
    "INSERT INTO annotated SELECT CHROM, 1100, REF, ALT, FILTER, INFO_AF,",
    "INFO_AN, INFO_AS_FilterStatus, INFO_OLD_MULTIALLELIC, INFO_OLD_VARIANT,",
    "INFO_segdup, INFO_lcr, INFO_BaseQRankSum, INFO_QD, INFO_MQ, 2.5",
    "FROM annotated WHERE CHROM = 'chr2' AND POS = 100"
  ))
  selected <- rduckhts_somalier_find_sites(
    con, source_table = "annotated", assembly = "GRCh38", min_an = 100,
    snp_dist = 10
  )
  expect_equal(selected$position, c(100, 2781480, 154931045, 1200))
  expect_equal(selected$source_filter[[1]], "PASS")
  dbExecute(con, "CREATE TEMP TABLE excluded(chrom VARCHAR, start BIGINT, stop BIGINT)")
  dbExecute(con, "INSERT INTO excluded VALUES ('chr2', 105, 106)")
  not_overlapping <- rduckhts_somalier_find_sites(
    con, source_table = "annotated", assembly = "GRCh38", min_an = 100,
    snp_dist = 10, exclude_table = "excluded"
  )
  expect_equal(not_overlapping$position, selected$position)
  dbExecute(con, "UPDATE excluded SET start = 104, stop = 105")
  overlapping <- rduckhts_somalier_find_sites(
    con, source_table = "annotated", assembly = "GRCh38", min_an = 100,
    snp_dist = 10, exclude_table = "excluded"
  )
  expect_equal(overlapping$position, selected$position[-1])
  missing_af <- rduckhts_somalier_find_sites(
    con, source_table = "annotated", assembly = "GRCh38", min_an = 100,
    min_af = 0, snp_dist = 10
  )
  expect_true(300 %in% missing_af$position)
  expect_true(is.na(missing_af$source_alt_af[missing_af$position == 300]))
  dbExecute(con, paste(
    "CREATE TEMP VIEW aliases AS SELECT CASE WHEN CHROM = 'chrX' THEN",
    "'NC_000023.11' ELSE 'NC_000024.10' END AS CHROM,",
    "POS, REF, ALT, FILTER, INFO_AF, INFO_AN, INFO_AS_FilterStatus,",
    "INFO_OLD_MULTIALLELIC, INFO_OLD_VARIANT, INFO_segdup, INFO_lcr,",
    "INFO_BaseQRankSum, INFO_QD, INFO_MQ, INFO_FS",
    "FROM annotated WHERE CHROM IN ('chrX', 'chrY')"
  ))
  alias_sites <- rduckhts_somalier_find_sites(
    con, source_table = "aliases", assembly = "GRCh38", min_an = 100,
    snp_dist = 10
  )
  expect_equal(alias_sites$position, c(2781480, 154931045, 1200))
}

test_somalier_find_sites_af_edges <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  dbExecute(con, paste(
    "CREATE TEMP TABLE frequency_edges AS SELECT * FROM (VALUES",
    "('1',100,'A',['G'],['PASS'],[0.15::FLOAT],115000),",
    "('1',200,'A',['G'],['PASS'],[0.85::FLOAT],115000),",
    "('1',300,'G',['A'],['PASS'],[0.16::FLOAT],115000),",
    "('1',400,'G',['A'],['PASS'],[0.84::FLOAT],115000),",
    "('1',500,'A',['G'],['PASS'],[0.48::FLOAT],114999),",
    "('Y',700,'A',['G'],['PASS'],[0.04::FLOAT],1),",
    "('Y',1100,'A',['G'],['PASS'],[0.05::FLOAT],1)",
    ") v(CHROM,POS,REF,ALT,FILTER,INFO_AF,INFO_AN)"
  ))
  selected <- rduckhts_somalier_find_sites(
    con, source_table = "frequency_edges", assembly = "GRCh38",
    snp_dist = 10
  )
  expect_equal(selected$position, c(100, 300, 400, 1100))
  expect_equal(round(selected$population_b_af[selected$position == 300], 2),
               0.84)
}

test_somalier_find_sites_vcf <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  path <- system.file("extdata", "somalier_population_1000g.vcf.gz",
                      package = "Rduckhts", mustWork = TRUE)
  selected <- rduckhts_somalier_find_sites(
    con, source_vcf = path, assembly = "GRCh37", min_an = 6, snp_dist = 100
  )
  expect_equal(nrow(selected), 3L)
  expect_equal(selected$position, c(10583, 11508, 16378))
  dbExecute(con, sprintf(
    "CREATE TEMP VIEW population_vcf AS SELECT * FROM read_bcf(%s, samples := '', scan_mode := 'sequential')",
    as.character(dbQuoteString(con, path))
  ))
  from_view <- rduckhts_somalier_find_sites(
    con, source_table = "population_vcf", assembly = "GRCh37", min_an = 6,
    snp_dist = 100
  )
  expect_equal(from_view[, names(from_view) != "source_path"],
               selected[, names(selected) != "source_path"])
  parquet <- tempfile(fileext = ".parquet")
  on.exit(unlink(parquet), add = TRUE)
  dbExecute(con, sprintf("COPY population_vcf TO %s (FORMAT PARQUET)",
                         as.character(dbQuoteString(con, parquet))))
  from_parquet <- rduckhts_somalier_find_sites(
    con, source_parquet = parquet, assembly = "GRCh37", min_an = 6,
    snp_dist = 100
  )
  expect_equal(from_parquet[, names(from_parquet) != "source_path"],
               selected[, names(selected) != "source_path"])
}

test_somalier_find_sites_selection()
test_somalier_find_sites_ties_and_sex_spacing()
test_somalier_find_sites_nearby()
test_somalier_find_sites_sorted_neighbors()
test_somalier_find_sites_sparse_events()
test_somalier_find_sites_flags_and_endpoints()
test_somalier_find_sites_af_edges()
test_somalier_find_sites_vcf()
