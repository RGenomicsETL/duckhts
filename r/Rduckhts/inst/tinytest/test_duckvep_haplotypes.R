library(tinytest)
library(DBI)

local({
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  tx <- paste("SELECT 0::UINTEGER transcript_index, 0::UINTEGER seq_region,",
    "100::UBIGINT transcript_start,111::UBIGINT transcript_end,1::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,",
    "111::UBIGINT cds_end,'AAAAAAAAAAAA'::BLOB cds_sequence,1::UTINYINT codon_table")
  exons <- paste("SELECT 0::UINTEGER transcript_index,100::UBIGINT exon_start,111::UBIGINT exon_end,",
    "1::UBIGINT exon_cdna_start,12::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('haps',",
    dbQuoteString(con, "SELECT 0::UINTEGER seq_region"), ",",
    dbQuoteString(con, tx), ",", dbQuoteString(con, exons), ")"))$loaded)
  empty <- dbGetQuery(con, paste("SELECT a.consequence,a.status,a.reason FROM",
    "unnest(_duckvep_annotate_small_rich('haps',0::UINTEGER,1::UBIGINT,'A','C',0::UBIGINT)) u(a)"))
  expect_identical(empty, data.frame(consequence = "sequence_variant", status = "unresolved",
    reason = "no_feature_in_loaded_model"))
  calls <- paste("SELECT event_index,0 seq_region,position,'A' reference,alternate,1 alt_index,",
    "0 transcript_index,0 sample_index,alleles,[true,true] phase_before,phase_set FROM",
    "(VALUES (1,100,'C',[1,1],NULL),(2,101,'G',[1,0],10),(3,102,'C',[0,1],20))",
    "v(event_index,position,alternate,alleles,phase_set)")
  actual <- rduckhts_haplotypes(con, calls, "haps")
  expect_equal(sort(actual$cds), c("CAAAAAAAAAAA", "CACAAAAAAAAA", "CGAAAAAAAAAA"))
  expect_equal(sort(actual$carrier_count), c(1,1,2))
  expect_equal(sort(actual$protein), c("HKKK", "QKKK", "RKKK"))
  expect_true(all(actual$sequence_status == "ok"))
  expect_true(all(actual$projection_status == "ok"))
  expect_equal(sum(lengths(lapply(actual$contributors, function(x) x$event_index))), 5L)
  expect_equal(nrow(rduckhts_haplotypes(con, calls, "haps", "vep116_compat")), 2L)
  noncoding <- paste("SELECT * REPLACE(0::UBIGINT AS transcript_flags,",
    "NULL::UBIGINT AS cds_start,NULL::UBIGINT AS cds_end,NULL::BLOB AS cds_sequence,",
    "NULL::UTINYINT AS codon_table) FROM (", tx, ")")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('noncoding',",
    dbQuoteString(con, "SELECT 0::UINTEGER seq_region"), ",",
    dbQuoteString(con, noncoding), ",", dbQuoteString(con, exons), ")"))$loaded)
  unavailable <- rduckhts_haplotypes(con, calls, "noncoding")
  expect_true(nrow(unavailable) == 3L && sum(unavailable$carrier_count) == 4L &&
    all(unavailable$projection_status == "outside_cds") &&
    all(is.na(unavailable$cds)) && all(is.na(unavailable$protein)))
  expect_error(rduckhts_haplotypes(con, calls, "haps", max_ploidy = 1), pattern = "max_ploidy")
  expect_error(rduckhts_haplotypes(con, calls, "haps", workspace_limit = 1), pattern = "workspace")
  expect_error(rduckhts_haplotypes(con, calls, "missing"), pattern = "loaded model")
  expect_error(rduckhts_haplotypes(con, calls, "haps", max_ploidy = -1), pattern = "positive")
  expect_error(rduckhts_haplotypes(con, calls, "haps", max_ploidy = 1.5), pattern = "integer")
  expect_error(rduckhts_haplotypes(con, NA_character_, "haps"), pattern = "nonempty")
  expect_error(rduckhts_haplotypes(con, calls, "haps", "strict", 12), pattern = "parameter names")
  rduckhts_haplotypes(con, calls, "haps", table_name = "hap_output")
  expect_equal(dbGetQuery(con, "SELECT count(*) n FROM hap_output")$n, 3)
  expect_error(rduckhts_haplotypes(con, calls, "haps", table_name = "hap_output"))
  rduckhts_haplotypes(con, calls, "haps", "vep116_compat", table_name = "hap_output", overwrite = TRUE)
  expect_equal(dbGetQuery(con, "SELECT count(*) n FROM hap_output")$n, 2)
  expect_equal(dbGetQuery(con, "SELECT 42 n")$n, 42L)
})
