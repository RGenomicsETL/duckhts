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
  blocks <- do.call(rbind, actual$coding_blocks)
  blocks <- blocks[order(blocks$alternate), ]
  expect_equal(blocks$cds_start, rep(1, 3))
  expect_equal(blocks$reference, c("A", "AAA", "AA"))
  expect_equal(blocks$alternate, c("C", "CAC", "CG"))
  expect_equal(blocks$alt_start0, rep(0, 3))
  expect_equal(blocks$edit_count, c(1, 2, 2))
  expect_true(all(blocks$length_change == 0 & blocks$sequence_flags == 0))
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
  null_blocks <- dbGetQuery(con, paste0("SELECT bool_and(coding_blocks IS NULL) ok FROM ",
    "duckvep_haplotypes(", dbQuoteString(con, calls), ",'noncoding')"))
  expect_true(null_blocks$ok)
  # One restored-frame block followed by an independent substitution. Both
  # input records in the block stay visible in contributor provenance.
  indels <- paste("SELECT event_index,0 seq_region,position,reference,alternate,1 alt_index,",
    "0 transcript_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set FROM",
    "(VALUES (1,100,'A','AT'),(2,106,'AA','A'),(3,110,'A','G'))",
    "v(event_index,position,reference,alternate)")
  restored <- rduckhts_haplotypes(con, indels, "haps", max_leaf_edits = 3)
  blocks <- restored$coding_blocks[[1]]
  expect_equal(nrow(restored), 1L)
  expect_equal(blocks$cds_start, c(2, 11))
  expect_equal(blocks$reference, c("AAAAAAA", "A"))
  expect_equal(blocks$alternate, c("TAAAAAA", "G"))
  expect_equal(blocks$alt_start0, c(1, 10))
  expect_equal(blocks$sequence_flags, c(5, 0))
  expect_equal(blocks$edit_count, c(2, 1))
  expect_equal(nrow(restored$contributors[[1]]), 3L)
  stop_tx <- paste("SELECT * REPLACE('ATGAAATAACCC'::BLOB AS cds_sequence) FROM (", tx, ")")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('stops',",
    dbQuoteString(con, "SELECT 0::UINTEGER seq_region"), ",",
    dbQuoteString(con, stop_tx), ",", dbQuoteString(con, exons), ")"))$loaded)
  stop_calls <- paste("SELECT event_index,0 seq_region,position,reference,alternate,1 alt_index,",
    "0 transcript_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set",
    "FROM (VALUES (1,103,'A','G'),(2,109,'C','A')) v(event_index,position,reference,alternate)")
  stopped <- rduckhts_haplotypes(con, stop_calls, "stops")
  expect_equal(stopped$cds, "ATGGAATAAACC")
  expect_equal(stopped$protein, "ME*")
  expect_equal(stopped$sequence_flags, 8)
  expect_equal(nrow(stopped$coding_blocks[[1]]), 2L)
  expect_equal(stopped$contributors[[1]]$event_index, 1:2)
  expect_error(rduckhts_haplotypes(con, calls, "haps", max_ploidy = 1), pattern = "max_ploidy")
  expect_error(rduckhts_haplotypes(con, calls, "haps", workspace_limit = 1), pattern = "workspace")
  expect_error(rduckhts_haplotypes(con, calls, "missing"), pattern = "loaded model")
  expect_error(rduckhts_haplotypes(con, calls, "haps", max_ploidy = -1), pattern = "positive")
  expect_error(rduckhts_haplotypes(con, calls, "haps", max_ploidy = 1.5), pattern = "integer")
  expect_error(rduckhts_haplotypes(con, NA_character_, "haps"), pattern = "nonempty")
  expect_error(rduckhts_haplotypes(con, calls, "haps", "strict", 12), pattern = "parameter names")
  rduckhts_haplotypes(con, calls, "haps", table_name = "hap_output")
  expect_equal(dbGetQuery(con, "SELECT count(*) n FROM hap_output")$n, 3)
  expect_equal(dbGetQuery(con, "SELECT sum(len(coding_blocks)) n FROM hap_output")$n, 3)
  expect_error(rduckhts_haplotypes(con, calls, "haps", table_name = "hap_output"))
  rduckhts_haplotypes(con, calls, "haps", "vep116_compat", table_name = "hap_output", overwrite = TRUE)
  expect_equal(dbGetQuery(con, "SELECT count(*) n FROM hap_output")$n, 2)
  expect_equal(dbGetQuery(con, "SELECT 42 n")$n, 42L)
})
