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
  overlapping <- paste("SELECT event_index,0 seq_region,position,reference,alternates,",
    "0 transcript_index,0 sample_index,gt FROM",
    "(VALUES (1,100,'AAA',['CAA'],'0|1'),(2,101,'A',['G'],'1|1'))",
    "v(event_index,position,reference,alternates,gt)")
  composed <- rduckhts_haplotypes(con, overlapping, "haps", "vep116_compat",
    input_mode = "source_records")
  composed <- composed[order(composed$cds), ]
  expect_equal(composed$cds, c("AAAAAAAAAAAA", "CAAAAAAAAAAA"))
  expect_equal(composed$protein, c("KKKK", "QKKK"))
  expect_equal(composed$edit_count, c(2, 2))
  expect_true(all(is.na(composed$stop_in_displaced_frame)))
  for (i in 1:2) {
    expect_equal(composed$coding_blocks[[i]]$coding_status, "unsupported_ordered_replacements")
    expect_true(is.na(composed$coding_blocks[[i]]$local_consequence_mask))
    expect_equal(composed$coding_blocks[[i]]$event_indices[[1]], c(1, 2))
    expect_equal(composed$contributors[[i]]$alt_index, c(i - 1L, 1L))
    expect_equal(composed$contributors[[i]]$evidence_flags, c(1L, 1L))
  }
  uncertain <- rduckhts_haplotypes(con, sub("0|1", ".|.", overlapping, fixed = TRUE),
    "haps", "vep116_compat", input_mode = "source_records")
  expect_equal(uncertain$cds, "AGAAAAAAAAAA")
  expect_equal(uncertain$sequence_status, "conditional")
  expect_equal(uncertain$edit_count, 1)
  expect_equal(uncertain$carrier_count, 2L)
  expect_error(rduckhts_haplotypes(con, overlapping, "haps", "vep116_compat",
    input_mode = "source_records", max_leaf_edits = 1), pattern = "max_leaf_edits")
  tied <- paste("SELECT event_index,0 seq_region,position,'A' AS reference,alternates,",
    "0 transcript_index,0 sample_index,gt FROM (VALUES",
    "(1,100,['C'],'1|1'),(2,100,['G'],'1|1'),(3,109,['T'],'0|0'))",
    "r(event_index,position,alternates,gt)")
  for (keep_context in c(FALSE, TRUE)) {
    query <- paste(tied, if (keep_context) "" else "WHERE event_index<3", "ORDER BY event_index DESC")
    result <- rduckhts_haplotypes(con, query, "haps", "vep116_compat", input_mode = "source_records")
    expect_equal(result$cds, if (keep_context) "CAAAAAAAAAAA" else "GAAAAAAAAAAA")
    expect_equal(result$edit_count, 2)
    expect_equal(result$carrier_count, 2L)
  }
  duplicates <- paste("SELECT event_index,0 seq_region,100 AS position,'A' AS reference,",
    "['C'] alternates,0 transcript_index,0 sample_index,gt",
    "FROM (VALUES (1,'1|0'),(2,'0|1')) r(event_index,gt)")
  result <- rduckhts_haplotypes(con, duplicates, "haps", "vep116_compat", input_mode = "source_records")
  result <- result[order(result$cds), ]
  expect_equal(result$cds, c("AAAAAAAAAAAA", "CAAAAAAAAAAA"))
  expect_equal(result$contributors[[1]]$projection_status, "shadowed_duplicate")
  expect_equal(result$contributors[[1]]$event_index, 1)
  expect_equal(result$edit_count, c(0, 1))
  expect_equal(result$coding_blocks[[2]]$event_indices[[1]], 2)
  raw_calls <- paste("SELECT i event_index,0 seq_region,99+i AS position,'A' AS reference,",
    "CASE i WHEN 1 THEN ['C','T'] ELSE ['G'] END alternates,0 transcript_index,0 sample_index,",
    "CASE i WHEN 1 THEN '.|1' ELSE '1|1' END gt FROM range(1,3) r(i)")
  raw <- rduckhts_haplotypes(con, raw_calls, "haps", "vep116_compat", input_mode = "source_records")
  raw <- raw[order(raw$cds), ]
  expect_equal(raw$cds, c("CGAAAAAAAAAA", "GAAAAAAAAAA"))
  expect_equal(raw$protein, c("RKKK", "EKK"))
  expect_equal(raw$sequence_status, rep("conditional", 2L))
  expect_equal(raw$evidence_flags, c(11L, 11L))
  expect_equal(raw$edit_count, c(2, 2))
  expect_equal(raw$carrier_count, c(1L, 1L))
  expect_equal(raw$contributors[[1]]$alt_index, c(1L, 1L))
  expect_equal(raw$contributors[[2]]$alt_index, c(NA_integer_, 1L))
  expect_equal(raw$contributors[[2]]$evidence_flags, c(10L, 1L))
  expect_equal(raw$contributors[[2]]$alternate, c("", "G"))
  expect_equal(raw$contributors[[2]]$event_index, c(1, 2))
  expect_equal(do.call(rbind, raw$carriers)$haplotype_lane, 1:2)
  expect_true(all(is.na(do.call(rbind, raw$carriers)$phase_set)))
  for (spelling in c("0|1", "|0|1", ".", "1")) {
    result <- rduckhts_haplotypes(con, sub(".|1", spelling, raw_calls, fixed = TRUE), "haps",
      "vep116_compat", input_mode = "source_records")
    expect_equal(sum(result$carrier_count), 2L)
    expect_equal(nrow(result), if (spelling %in% c("0|1", "1")) 2L else 1L)
  }
  for (replacement in c("NULL AS alternates", "['C',NULL] AS alternates", "[''] AS alternates",
                        "NULL AS gt", "'3|1' AS gt", "'1||1' AS gt")) {
    malformed <- paste("SELECT * REPLACE(", replacement, ") FROM (", raw_calls, ")")
    expect_error(rduckhts_haplotypes(con, malformed, "haps", "vep116_compat",
      input_mode = "source_records"))
  }
  for (name in c("event_index", "seq_region", "position", "reference", "transcript_index", "sample_index")) {
    null_input <- paste("SELECT * REPLACE(NULL AS", name, ") FROM (", raw_calls, ")")
    expect_error(rduckhts_haplotypes(con, null_input, "haps", "vep116_compat",
      input_mode = "source_records"), pattern = "required input")
  }
  for (replacement in c("'T' AS reference", "['T','C'] AS alternates", "101 AS position")) {
    inconsistent <- paste("SELECT * FROM (", raw_calls, ") UNION ALL SELECT * REPLACE(",
      "1 AS sample_index,", replacement, ") FROM (", raw_calls, ")")
    expect_error(rduckhts_haplotypes(con, inconsistent, "haps", "vep116_compat",
      input_mode = "source_records"), pattern = "inconsistent source record identity")
  }
  different_gt <- paste("SELECT * FROM (", raw_calls, ") UNION ALL SELECT * REPLACE(",
    "1 AS transcript_index,'1|1' AS gt) FROM (", raw_calls, ")")
  expect_error(rduckhts_haplotypes(con, different_gt, "haps", "vep116_compat",
    input_mode = "source_records"), pattern = "source GT")
  expect_error(rduckhts_haplotypes(con, paste(raw_calls, "UNION ALL", raw_calls), "haps",
    "vep116_compat", input_mode = "source_records"), pattern = "duplicate call")
  expect_error(rduckhts_haplotypes(con, raw_calls, "haps", input_mode = "source_records"),
    pattern = "requires phase_policy")
  expect_error(rduckhts_haplotypes(con, raw_calls, "haps", "vep116_compat", input_mode = "unknown"))
  expect_error(rduckhts_haplotypes(con, raw_calls, "haps", "vep116_compat",
    input_mode = "source_records", max_ploidy = 1), pattern = "max_ploidy")
  rduckhts_haplotypes(con, raw_calls, "haps", "vep116_compat", input_mode = "source_records",
    table_name = "raw_haplotypes")
  expect_equal(dbGetQuery(con, "SELECT sum(carrier_count) n FROM raw_haplotypes")$n, 2)
  for (policy in c("strict", "vep116_compat")) {
    # Changing any source ALT identity across samples must fail even without
    # duplicate (event, transcript, sample) keys. Input order cannot hide it.
    for (replacement in c("1 AS seq_region", "position+1 AS position", "'T' AS reference",
                          "'T' AS alternate", "2 AS alt_index")) {
      inconsistent <- paste("SELECT * FROM (", calls, ") UNION ALL SELECT * REPLACE(",
        "1 AS sample_index,", replacement, ") FROM (", calls, ") ORDER BY sample_index DESC")
      expect_error(rduckhts_haplotypes(con, inconsistent, "haps", phase_policy = policy),
        pattern = "inconsistent event identity")
    }
    for (column in c("event_index", "transcript_index", "sample_index")) {
      null_key <- paste("SELECT * REPLACE(NULL AS", column, ") FROM (", calls, ")")
      expect_error(rduckhts_haplotypes(con, null_key, "haps", phase_policy = policy),
        pattern = "required input")
    }
    ref_only <- paste("SELECT * REPLACE([0,0] AS alleles) FROM (", calls, ")")
    expect_error(rduckhts_haplotypes(con, paste(ref_only, "UNION ALL", ref_only),
      "haps", phase_policy = policy), pattern = "duplicate call")
    changed_ploidy <- paste("SELECT * REPLACE(CASE WHEN event_index=3 THEN [1] ELSE alleles END",
      "AS alleles,NULL::BOOLEAN[] AS phase_before) FROM (", calls, ")")
    expect_error(rduckhts_haplotypes(con, changed_ploidy, "haps", phase_policy = policy),
      pattern = "ploidy")
    cohort <- paste("SELECT * FROM (", calls, ") UNION ALL SELECT * REPLACE(",
      "1 AS sample_index,[1] AS alleles,[true] AS phase_before) FROM (", calls, ")")
    mixed <- rduckhts_haplotypes(con, cohort, "haps", phase_policy = policy)
    carriers <- do.call(rbind, mixed$carriers)
    expect_equal(sum(mixed$carrier_count), if (policy == "strict") 5 else 3)
    expect_true(all(carriers$ploidy == ifelse(carriers$sample_index == 0, 2, 1)))
  }
  n_tx <- sub("AAAAAAAAAAAA", "ATGGCNTGNGCC", tx, fixed = TRUE)
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('n_codons',",
    dbQuoteString(con, "SELECT 0::UINTEGER seq_region"), ",", dbQuoteString(con, n_tx), ",",
    dbQuoteString(con, exons), ")"))$loaded)
  n_call <- paste("SELECT 1 event_index,0 seq_region,102 AS position,'G' AS reference,'A' alternate,",
    "1 alt_index,0 transcript_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set")
  for (policy in c("strict", "vep116_compat")) {
    n_result <- rduckhts_haplotypes(con, n_call, "n_codons", phase_policy = policy)
    expect_equal(n_result$cds, "ATAGCNTGNGCC")
    expect_equal(n_result$protein, "IAXA")
    expect_identical(n_result$stop_in_displaced_frame, FALSE)
    expect_equal(n_result$sequence_status, "ok")
    expect_equal(nrow(n_result$contributors[[1L]]), 1L)
    n_changed <- rduckhts_haplotypes(con, sub("102 AS position", "103 AS position", n_call),
      "n_codons", phase_policy = policy)
    expect_equal(n_changed$cds, "ATGACNTGNGCC")
    expect_equal(n_changed$protein, "MTXA")
    expect_equal(n_changed$sequence_status, "ok")
    expect_equal(n_changed$coding_blocks[[1]]$coding_status, "unsupported")
    expect_true(is.na(n_changed$coding_blocks[[1]]$local_consequence_mask))
  }
  expect_equal(sort(actual$cds), c("CAAAAAAAAAAA", "CACAAAAAAAAA", "CGAAAAAAAAAA"))
  expect_equal(sort(actual$carrier_count), c(1,1,2))
  expect_equal(sort(actual$protein), c("HKKK", "QKKK", "RKKK"))
  expect_true(all(actual$sequence_status == "ok"))
  expect_identical(actual$stop_in_displaced_frame, rep(FALSE, nrow(actual)))
  expect_true(all(actual$projection_status == "ok"))
  expect_equal(sum(lengths(lapply(actual$contributors, function(x) x$event_index))), 5L)
  blocks <- do.call(rbind, actual$coding_blocks)
  blocks <- blocks[order(blocks$alternate), ]
  missense_mask <- dbGetQuery(con,
    "SELECT consequence_mask FROM duckvep_so_terms() WHERE consequence='missense_variant'")$consequence_mask
  synonymous_mask <- dbGetQuery(con,
    "SELECT consequence_mask FROM duckvep_so_terms() WHERE consequence='synonymous_variant'")$consequence_mask
  start_lost_mask <- dbGetQuery(con,
    "SELECT consequence_mask FROM duckvep_so_terms() WHERE consequence='start_lost'")$consequence_mask
  expect_true(all(blocks$coding_status == "ok"))
  expect_equal(blocks$local_consequence_mask, rep(start_lost_mask, 3L))
  expect_true(all(!blocks$after_first_stop))
  expect_equal(blocks$cds_start, rep(1, 3))
  expect_equal(blocks$reference, c("A", "AAA", "AA"))
  expect_equal(blocks$alternate, c("C", "CAC", "CG"))
  expect_equal(blocks$alt_start0, rep(0, 3))
  expect_equal(lengths(blocks$event_indices), c(1, 2, 2))
  expect_equal(blocks$event_indices, list(1, c(1, 3), c(1, 2)))
  cis_tx <- sub("AAAAAAAAAAAA", "ATGTCTGCCTAA", tx, fixed = TRUE)
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('cis',",
    dbQuoteString(con, "SELECT 0::UINTEGER seq_region"), ",", dbQuoteString(con, cis_tx), ",",
    dbQuoteString(con, exons), ")"))$loaded)
  cis_calls <- paste("SELECT event_index,0 seq_region,position,reference,alternate,1 alt_index,",
    "0 transcript_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set",
    "FROM (VALUES (1,103,'T','A'),(2,104,'C','G')) v(event_index,position,reference,alternate)")
  for (policy in c("strict", "vep116_compat")) {
    cis <- rduckhts_haplotypes(con, cis_calls, "cis", phase_policy = policy)
    expect_equal(cis$cds, "ATGAGTGCCTAA")
    expect_equal(cis$protein, "MSA*")
    expect_equal(cis$coding_blocks[[1]]$coding_status, "ok")
    expect_equal(cis$coding_blocks[[1]]$local_consequence_mask, synonymous_mask)
    expect_equal(cis$coding_blocks[[1]]$event_indices, list(c(1, 2)))
    expect_equal(nrow(cis$protein_differences[[1]]), 0L)
  }
  expect_true(all(blocks$length_change == 0 & blocks$sequence_flags == 0))
  expect_equal(nrow(rduckhts_haplotypes(con, calls, "haps", "vep116_compat")), 2L)
  noncoding <- paste("SELECT * REPLACE(0::UBIGINT AS transcript_flags,",
    "NULL::UBIGINT AS cds_start,NULL::UBIGINT AS cds_end,NULL::BLOB AS cds_sequence,",
    "NULL::UTINYINT AS codon_table) FROM (", tx, ")")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('noncoding',",
    dbQuoteString(con, "SELECT 0::UINTEGER seq_region"), ",",
    dbQuoteString(con, noncoding), ",", dbQuoteString(con, exons), ")"))$loaded)
  unavailable <- rduckhts_haplotypes(con, calls, "noncoding")
  expect_true(all(is.na(unavailable$stop_in_displaced_frame)))
  expect_true(nrow(unavailable) == 3L && sum(unavailable$carrier_count) == 4L &&
    all(unavailable$projection_status == "outside_cds") &&
    all(is.na(unavailable$cds)) && all(is.na(unavailable$protein)))
  null_blocks <- dbGetQuery(con, paste0("SELECT bool_and(coding_blocks IS NULL) ok FROM ",
    "duckvep_haplotypes(", dbQuoteString(con, calls), ",'noncoding')"))
  expect_true(null_blocks$ok)
  mixed_tx <- paste("SELECT i::UINTEGER transcript_index,0::UINTEGER seq_region,",
    "100::UBIGINT transcript_start,151::UBIGINT transcript_end,",
    "(CASE i WHEN 0 THEN 1 ELSE -1 END)::TINYINT strand,0::UINTEGER gene_index,",
    "3::UBIGINT transcript_flags,103::UBIGINT cds_start,148::UBIGINT cds_end,",
    "(CASE i WHEN 0 THEN 'AAAAAAAAAAAA' ELSE 'TTTTTTTTTTTT' END)::BLOB cds_sequence,",
    "1::UTINYINT codon_table FROM range(2) t(i)")
  dbExecute(con, paste("CREATE TABLE mixed_tx AS", mixed_tx))
  mixed_exons <- paste("SELECT transcript_index,s::UBIGINT exon_start,e::UBIGINT exon_end,",
    "(CASE WHEN (strand=1 AND s=100) OR (strand=-1 AND s=143) THEN 1 ELSE 10 END)::UBIGINT exon_cdna_start,",
    "(exon_cdna_start+8)::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase",
    "FROM mixed_tx CROSS JOIN (VALUES (100,108),(143,151)) v(s,e)",
    "ORDER BY transcript_index,exon_cdna_start")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('mixed',",
    dbQuoteString(con, "SELECT 0::UINTEGER seq_region"), ",",
    dbQuoteString(con, "SELECT * FROM mixed_tx"), ",", dbQuoteString(con, mixed_exons), ")"))$loaded)
  mixed_calls <- paste("SELECT event_index,0 seq_region,position,'A' reference,alternate,1 alt_index,",
    "transcript_index,s.i sample_index,[CASE WHEN s.i=1 AND event_index=2 THEN 0 ELSE 1 END,0] alleles,",
    "[true,true] phase_before,10 phase_set FROM (VALUES (1,100,'C'),(2,105,'G'),",
    "(3,115,'C'),(4,149,'C')) v(event_index,position,alternate)",
    "CROSS JOIN mixed_tx CROSS JOIN range(2) s(i)")
  for (policy in c("strict", "vep116_compat")) {
    mixed <- rduckhts_haplotypes(con, mixed_calls, "mixed", phase_policy = policy)
    mixed <- mixed[order(mixed$transcript_index, mixed$cds), ]
    expect_equal(mixed$cds, c("AAAAAAAAAAAA", "AAGAAAAAAAAA", "TTTTTTTTTCTT", "TTTTTTTTTTTT"))
    expect_true(all(mixed$projection_status == "ok" & mixed$sequence_status == "ok"))
    expect_equal(mixed$edit_count, c(0, 1, 1, 0))
    expect_equal(vapply(mixed$contributors, nrow, 1L), c(3L, 4L, 4L, 3L))
    contributors <- do.call(rbind, mixed$contributors)
    expect_equal(as.integer(table(contributors$projection_status)[c("ok", "outside_cds")]), c(2L, 12L))
    expect_equal(vapply(mixed$coding_blocks, nrow, 1L), c(0L, 1L, 1L, 0L))
    for (site in c(115, 142)) {
      insertion <- paste("SELECT * REPLACE(CASE WHEN event_index=3 THEN", site,
        "ELSE position END AS position,CASE WHEN event_index=3 THEN 'AC' ELSE alternate END AS alternate)",
        "FROM (", mixed_calls, ")")
      insertion <- rduckhts_haplotypes(con, insertion, "mixed", phase_policy = policy)
      expect_equal(nrow(insertion), 4L)
      expect_equal(nchar(insertion$cds), rep(if (site == 115) 12L else 13L, 4L))
      expect_equal(sum(insertion$edit_count), if (site == 115) 2 else 6)
    }
    uncertain <- paste("SELECT * REPLACE(CASE WHEN event_index=3 THEN [NULL,0] ELSE alleles END AS alleles)",
      "FROM (", mixed_calls, ")")
    uncertain <- rduckhts_haplotypes(con, uncertain, "mixed", phase_policy = policy)
    # Compatibility compacts the called REF into lane 1, preserving the missing
    # slot as a separate lane-2 prefix shared by the two samples per transcript.
    expect_equal(nrow(uncertain), if (policy == "strict") 4L else 6L)
    expect_equal(sum(uncertain$carrier_count), if (policy == "strict") 4 else 8)
    expect_true(all(is.na(uncertain$cds)) &&
      all(uncertain$sequence_status == "incomplete_input"))
    crossing <- paste("SELECT * REPLACE(CASE WHEN event_index=1 THEN 102 ELSE position END AS position,",
      "CASE WHEN event_index=1 THEN 'AA' ELSE reference END AS reference,",
      "CASE WHEN event_index=1 THEN 'CC' ELSE alternate END AS alternate) FROM (", mixed_calls, ")")
    crossing <- rduckhts_haplotypes(con, crossing, "mixed", phase_policy = policy)
    expect_true(nrow(crossing) == 4L && all(is.na(crossing$cds)) &&
      all(crossing$projection_status == "outside_cds"))
  }
  unmapped_calls <- paste("SELECT event_index,0 seq_region,position,reference,alternates,",
    "transcript_index,s.i sample_index,CASE WHEN event_index=2 THEN '1|1'",
    "ELSE ['1|1','0|1','.','.|1','0|0'][s.i+1] END gt FROM",
    "(VALUES (1,107,repeat('A',38),['C']),(2,147,'A',['G']))",
    "v(event_index,position,reference,alternates) CROSS JOIN mixed_tx CROSS JOIN range(5) s(i)")
  unmapped <- rduckhts_haplotypes(con, unmapped_calls, "mixed", "vep116_compat",
    input_mode = "source_records", max_leaf_edits = 1)
  expect_equal(nrow(unmapped), 12L)
  expect_equal(sum(unmapped$carrier_count), 20L)
  expect_equal(sum(unmapped$carrier_count[unmapped$sequence_status == "conditional"]), 16L)
  expect_true(all(unmapped$projection_status == "ok" & unmapped$edit_count == 1))
  expect_true(all(unmapped$cds[unmapped$transcript_index == 0] == "AAAAAAAAAAGA"))
  expect_true(all(unmapped$cds[unmapped$transcript_index == 1] == "TCTTTTTTTTTT"))
  for (i in seq_len(nrow(unmapped))) {
    expect_equal(unmapped$coding_blocks[[i]]$event_indices[[1]], 2)
    source <- subset(unmapped$contributors[[i]], event_index == 1)
    if (nrow(source)) {
      expect_equal(source$projection_status, "source_unmapped")
      expect_equal(bitwAnd(source$evidence_flags, 8L), 8L)
      expect_equal(source$reference, strrep("A", 38L))
      expect_equal(source$position, 107)
    }
  }
  mismatch <- paste("SELECT * REPLACE(CASE WHEN event_index=1 THEN 'C'||repeat('A',37)",
    "ELSE reference END AS reference) FROM (", unmapped_calls, ")")
  mismatch <- rduckhts_haplotypes(con, mismatch, "mixed", "vep116_compat", input_mode = "source_records")
  expect_equal(sum(mismatch$carrier_count[is.na(mismatch$cds)]), 16L)
  expect_equal(sum(mismatch$carrier_count[mismatch$projection_status == "reference_mismatch"]), 16L)
  decoded <- paste("SELECT *,alternates[1] alternate,1 alt_index,[1,1] alleles,",
    "[true,true] phase_before,NULL::BIGINT phase_set FROM (", unmapped_calls, ") WHERE sample_index=0")
  decoded <- rduckhts_haplotypes(con, decoded, "mixed")
  expect_true(all(is.na(decoded$cds)) && all(decoded$projection_status == "outside_cds"))
  expect_equal(sum(decoded$carrier_count), 4L)
  expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('mixed') dropped")$dropped)
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
  expect_equal(lengths(blocks$event_indices), c(2, 1))
  expect_equal(blocks$event_indices, list(c(1, 2), 3))
  expect_equal(nrow(restored$contributors[[1]]), 3L)
  identical_indels <- paste("SELECT * REPLACE('AA' AS alternate) FROM (", indels,
    ") WHERE event_index=1 UNION ALL SELECT * FROM (", indels, ") WHERE event_index=2")
  for (policy in c("strict", "vep116_compat")) {
    identity <- rduckhts_haplotypes(con, identical_indels, "haps", phase_policy = policy)
    expect_equal(identity$cds, "AAAAAAAAAAAA")
    expect_equal(identity$edit_count, 2)
    expect_equal(identity$sequence_flags, 5)
    expect_equal(identity$contributors[[1]]$event_index, c(1, 2))
    expect_equal(nrow(identity$cds_differences[[1]]), 0L)
    expect_equal(nrow(identity$protein_differences[[1]]), 0L)
    blocks <- identity$coding_blocks[[1]]
    expect_equal(nrow(blocks), 1L)
    expect_equal(blocks$reference, blocks$alternate)
    expect_equal(blocks$sequence_flags, 5)
    expect_equal(blocks$event_indices[[1]], c(1, 2))
    expect_equal(blocks$coding_status, "ok")
    expect_equal(blocks$local_consequence_mask, synonymous_mask)
    expect_identical(blocks$after_first_stop, FALSE)
  }
  reverse_tx <- paste("SELECT * REPLACE(1::UINTEGER AS transcript_index,-1::TINYINT AS strand,",
    "'TTTTTTTTTTTT'::BLOB AS cds_sequence) FROM (", tx, ")")
  reverse_exons <- paste("SELECT * REPLACE(1::UINTEGER AS transcript_index) FROM (", exons, ")")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('islands',",
    dbQuoteString(con, "SELECT 0::UINTEGER seq_region"), ",",
    dbQuoteString(con, paste(tx, "UNION ALL", reverse_tx)), ",",
    dbQuoteString(con, paste(exons, "UNION ALL", reverse_exons)), ")"))$loaded)
  island_calls <- paste("SELECT event_index,0 seq_region,position,reference,alternate,1 alt_index,",
    "t.transcript_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set",
    "FROM (VALUES (17,100,'AAAAAAAA','CACAAAAC'),(19,101,'A','G'),(18,104,'A','G'))",
    "e(event_index,position,reference,alternate) CROSS JOIN (VALUES (0),(1)) t(transcript_index)")
  for (policy in c("strict", "vep116_compat")) {
    islands <- rduckhts_haplotypes(con, island_calls, "islands", phase_policy = policy, max_leaf_edits = 5)
    islands <- islands[order(islands$transcript_index), ]
    expect_equal(islands$edit_count, c(5, 5))
    expect_equal(lapply(islands$contributors, function(x) x$event_index), rep(list(c(17, 19, 18)), 2))
    expect_equal(islands$coding_blocks[[1]]$event_indices, list(c(17, 19, 17), 18, 17))
    expect_equal(islands$coding_blocks[[2]]$event_indices, list(17, 18, c(17, 19, 17)))
    expect_equal(islands$coding_blocks[[1]]$cds_start, c(1, 5, 8))
    expect_equal(islands$coding_blocks[[2]]$cds_start, c(5, 8, 10))
    expect_error(rduckhts_haplotypes(con, island_calls, "islands", phase_policy = policy,
      max_leaf_edits = 4), pattern = "max_leaf_edits")
  }
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
  expect_identical(stopped$stop_in_displaced_frame, FALSE)
  expect_equal(nrow(stopped$coding_blocks[[1]]), 2L)
  expect_equal(stopped$coding_blocks[[1]]$coding_status, c("ok", "ok"))
  expect_equal(stopped$coding_blocks[[1]]$local_consequence_mask, rep(missense_mask, 2L))
  expect_identical(stopped$coding_blocks[[1]]$after_first_stop, c(FALSE, TRUE))
  expect_equal(stopped$contributors[[1]]$event_index, 1:2)
  expect_equal(stopped$cds_differences[[1]]$ref_start0, c(3, 9))
  expect_equal(stopped$cds_differences[[1]]$reference, c("A", "C"))
  expect_equal(stopped$cds_differences[[1]]$alternate, c("G", "A"))
  after_frame_calls <- paste("SELECT event_index,0 seq_region,position,reference,alternate,1 alt_index,",
    "0 transcript_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set",
    "FROM (VALUES (1,103,'A','AC'),(2,104,'AA','A')) v(event_index,position,reference,alternate)")
  for (policy in c("strict", "vep116_compat")) {
    after_frame <- rduckhts_haplotypes(con, after_frame_calls, "stops", phase_policy = policy)
    expect_equal(after_frame$cds, "ATGACATAACCC")
    expect_equal(after_frame$protein, "MT*")
    expect_equal(after_frame$sequence_flags, 13)
    expect_identical(after_frame$stop_in_displaced_frame, FALSE)
  }
  missing_calls <- sub("alleles,[true,true]", "[NULL,NULL]::INTEGER[] alleles,[true,true]", calls, fixed = TRUE)
  missing_frame <- rduckhts_haplotypes(con, missing_calls, "haps")
  expect_true(all(is.na(missing_frame$stop_in_displaced_frame)))
  expect_true(all(missing_frame$sequence_status == "incomplete_input"))
  # The DNA frame restores after translation has already reached a stop. These
  # are the unchanged seed-173 DHT000002 edits, also expressed on the forward strand.
  early_cds <- paste0("ATGGGCTTGCCTAGCTTAAAACACTTGGTAAACTTTCTCGTGGTTACAACCTCGTTAGGGCTGCAT",
    "CTTCTCTACATGGAGAAGTGCGAAGTCCGGGTAAGGACCGGGTGTTATAACCCAGCAGACAAGCT",
    "GAAGGGGCGCGTGTACTTAGCTGCCCGCTCGAGGCATTGGTCCCCGTAA")
  early_tx <- paste("SELECT i::UINTEGER transcript_index,0::UINTEGER seq_region,",
    "11::UBIGINT transcript_start,190::UBIGINT transcript_end,s::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,11::UBIGINT cds_start,",
    "190::UBIGINT cds_end,", dbQuoteString(con, early_cds), "::BLOB cds_sequence,",
    "1::UTINYINT codon_table FROM (VALUES (0,1),(1,-1)) t(i,s)")
  early_exons <- paste("SELECT transcript_index,transcript_start exon_start,transcript_end exon_end,",
    "1::UBIGINT exon_cdna_start,180::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase",
    "FROM (", early_tx, ")")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('early_stop',",
    dbQuoteString(con, "SELECT 0::UINTEGER seq_region"), ",", dbQuoteString(con, early_tx), ",",
    dbQuoteString(con, early_exons), ")"))$loaded)
  early_calls <- paste("SELECT event_index,0 seq_region,position,reference,alternate,1 alt_index,",
    "transcript_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set",
    "FROM (VALUES (1,0,19,'G','GT'),(2,0,49,'CG','C'),(3,1,150,'AC','A'),(4,1,181,'G','GA'))",
    "v(event_index,transcript_index,position,reference,alternate)")
  for (policy in c("strict", "vep116_compat")) {
    early <- rduckhts_haplotypes(con, early_calls, "early_stop", phase_policy = policy)
    early <- early[order(early$transcript_index), ]
    expect_equal(early$protein, c("MGLS*", "MGLS*"))
    expect_equal(nchar(early$cds), c(180L, 180L))
    expect_equal(early$sequence_flags, c(13, 13))
    expect_identical(early$stop_in_displaced_frame, c(TRUE, TRUE))
    expect_equal(lapply(early$contributors, function(x) x$event_index), list(c(1, 2), c(3, 4)))
    expect_equal(lapply(early$coding_blocks, function(x) x$event_indices), list(list(c(1, 2)), list(c(4, 3))))
    expect_equal(vapply(early$coding_blocks, function(x) x$length_change, 0), c(0, 0))
  }
  spelling_tx <- sub("ATGAAATAACCC", "atgaaataaccc", stop_tx, fixed = TRUE)
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('spelling',",
    dbQuoteString(con, "SELECT 0::UINTEGER seq_region"), ",", dbQuoteString(con, spelling_tx), ",",
    dbQuoteString(con, exons), ")"))$loaded)
  spelled <- rduckhts_haplotypes(con, stop_calls, "spelling")
  expect_equal(spelled$cds_differences, stopped$cds_differences)
  insertion <- paste("SELECT 1 event_index,0 seq_region,105 AS position,'A' AS reference,'AA' alternate,",
    "1 alt_index,0 transcript_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set")
  inserted <- rduckhts_haplotypes(con, insertion, "haps", max_alignment_cells = 39)
  expect_equal(inserted$cds_differences[[1]]$ref_start0, 12)
  expect_equal(inserted$cds_differences[[1]]$alt_start0, 12)
  expect_equal(inserted$cds_differences[[1]]$reference, "")
  expect_equal(inserted$cds_differences[[1]]$alternate, "A")
  expect_error(rduckhts_haplotypes(con, insertion, "haps", max_alignment_cells = 38),
    pattern = "max_alignment_cells=38, required=39")
  cancelling <- paste("SELECT event_index,0 seq_region,position,reference,alternate,",
    "1 alt_index,0 transcript_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set",
    "FROM (VALUES (1,101,'A','AA'),(2,108,'AA','A')) v(event_index,position,reference,alternate)")
  unchanged <- rduckhts_haplotypes(con, cancelling, "haps", max_alignment_cells = 1, max_leaf_differences = 1)
  expect_equal(nrow(unchanged$cds_differences[[1]]), 0L)
  expect_equal(nrow(unchanged$protein_differences[[1]]), 0L)
  expect_equal(nrow(unchanged$contributors[[1]]), 2L)
  ref_edits <- paste("SELECT 0::UINTEGER transcript_index,(i+1)::UINTEGER protein_position,",
    "'G'::VARCHAR alternate_amino_acid FROM range(4) r(i) ORDER BY protein_position")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('reference_limit',",
    dbQuoteString(con, "SELECT 0::UINTEGER seq_region"), ",", dbQuoteString(con, tx), ",",
    dbQuoteString(con, exons), ",peptide_edit_query := ", dbQuoteString(con, ref_edits), ")"))$loaded)
  reference_aligned <- rduckhts_haplotypes(con, cancelling, "reference_limit", max_alignment_cells = 25)
  expect_equal(nrow(reference_aligned$cds_differences[[1L]]), 0L)
  expect_equal(reference_aligned$protein_differences[[1L]]$reference, "GGGG")
  expect_equal(reference_aligned$protein_differences[[1L]]$alternate, "KKKK")
  expect_error(rduckhts_haplotypes(con, cancelling, "reference_limit", max_alignment_cells = 24),
    pattern = "protein difference status 3, max_alignment_cells=24, required=25")

  dbWriteTable(con, "reference_protein_source", data.frame(i = 0:5,
    cds = c("CTGGCCTAA", "ATGTGAGCCTAA", "ATGGCCTAA", "ATGGCCTGA", "ctggcctaa", "AT"),
    code = c(1L, 1L, 1L, 2L, 1L, 1L)))
  p_tx <- paste("SELECT i::UINTEGER transcript_index,0::UINTEGER seq_region,",
    "(100+20*i)::UBIGINT transcript_start,(99+20*i+len(cds))::UBIGINT transcript_end,",
    "1::TINYINT strand,0::UINTEGER gene_index,3::UBIGINT transcript_flags,",
    "transcript_start cds_start,transcript_end cds_end,cds::BLOB cds_sequence,",
    "code::UTINYINT codon_table FROM reference_protein_source")
  p_exons <- paste("SELECT transcript_index,transcript_start exon_start,transcript_end exon_end,",
    "1::UBIGINT exon_cdna_start,(transcript_end-transcript_start+1)::UBIGINT exon_cdna_end,",
    "0::TINYINT phase,0::TINYINT end_phase FROM (", p_tx, ")")
  p_edits <- paste("SELECT * FROM (VALUES (1::UINTEGER,2::UINTEGER,'U'::VARCHAR),",
    "(2::UINTEGER,3::UINTEGER,'W'::VARCHAR)) e(transcript_index,protein_position,alternate_amino_acid)")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('reference_proteins',",
    dbQuoteString(con, "SELECT 0::UINTEGER seq_region"), ",", dbQuoteString(con, p_tx), ",",
    dbQuoteString(con, p_exons), ",peptide_edit_query := ", dbQuoteString(con, p_edits), ")"))$loaded)
  p_calls <- paste("SELECT transcript_index+1 event_index,0 seq_region,transcript_start+",
    "CASE transcript_index WHEN 1 THEN 8 WHEN 5 THEN 1 ELSE 5 END AS position,",
    "CASE transcript_index WHEN 5 THEN 'T' ELSE 'C' END AS reference,",
    "CASE transcript_index WHEN 5 THEN 'C' ELSE 'T' END alternate,1 alt_index,transcript_index,",
    "0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set FROM (", p_tx, ")")
  for (policy in c("strict", "vep116_compat")) {
    p <- rduckhts_haplotypes(con, p_calls, "reference_proteins", phase_policy = policy)
    p <- p[order(p$transcript_index), ]
    expect_equal(p$protein, c("LA*", "M*", "MA*", "MAW", "LA*", ""))
    expect_true(is.null(p$protein_differences[[6L]]))
    expect_equal(vapply(p$protein_differences[1:5], nrow, 1L), c(1L, 2L, 2L, 1L, 2L))
    differences <- do.call(rbind, p$protein_differences[1:5])
    expect_equal(differences$reference, c("M", "U", "A*", "W", "*", "*", "M", ""))
    expect_equal(differences$alternate, c("L", "*", "", "*", "", "", "L", "*"))
    expect_equal(differences$ref_start0, c(0, 1, 2, 2, 3, 3, 0, 2))
    expect_equal(differences$alignment_start0, differences$ref_start0)
    expect_true(all(vapply(p$contributors, nrow, 1L) == 1L))
  }
  expect_error(rduckhts_haplotypes(con, paste("SELECT * FROM (", p_calls,
    ") WHERE transcript_index=1"), "reference_proteins", max_leaf_differences = 1),
    pattern = "protein difference status 4, max_leaf_differences=1, required=2")
  p_raw <- paste("SELECT event_index,seq_region,position,reference,[alternate] alternates,",
    "transcript_index,sample_index,'.' gt FROM (", p_calls, ")")
  p_missing <- rduckhts_haplotypes(con, p_raw, "reference_proteins", "vep116_compat",
    input_mode = "source_records")
  p_missing <- p_missing[order(p_missing$transcript_index), ]
  expect_equal(p_missing$protein, c("MA*", "MUA*", "MAW*", "MAW*", "MA", ""))
  expect_true(all(p_missing$edit_count == 0 & p_missing$carrier_count == 2 &
    p_missing$sequence_flags == 0 & p_missing$sequence_status == "conditional"))
  expect_true(all(vapply(p_missing$protein_differences[1:5], nrow, 1L) == 0L))
  expect_true(is.null(p_missing$protein_differences[[6L]]))
  p_retained <- rduckhts_haplotypes(con, sub("'.' gt", "'0|1' gt", p_raw, fixed = TRUE),
    "reference_proteins", "vep116_compat", input_mode = "source_records")
  p_retained <- p_retained[vapply(p_retained$carriers, function(x) 1L %in% x$haplotype_lane, TRUE), ]
  p_retained <- p_retained[order(p_retained$transcript_index), ]
  expect_equal(p_retained$protein, c("LA*", "M*", "MA*", "MAW", "LA*", ""))
  expect_true(all(p_retained$edit_count == 0))
  route_tx <- paste("SELECT i::UINTEGER transcript_index,i::UINTEGER seq_region,",
    "11::UBIGINT transcript_start,(22+gap)::UBIGINT transcript_end,1::TINYINT strand,",
    "i::UINTEGER gene_index,3::UBIGINT transcript_flags,transcript_start cds_start,",
    "transcript_end cds_end,cds::BLOB cds_sequence,1::UTINYINT codon_table,gap FROM",
    "(VALUES (0,12,'ATGAAACCCTAA'),(1,12,'ATGAAATAACCC'),(2,12,'CTGAAACCCTAA'),",
    "(3,24,'ATGAAACCCTAA'),(4,24,'ATGAAATAACCC'),(5,24,'CTGAAACCCTAA')) v(i,gap,cds)")
  dbExecute(con, paste("CREATE TABLE route_tx AS", route_tx))
  route_queries <- c("SELECT seq_region FROM route_tx", "SELECT * EXCLUDE(gap) FROM route_tx",
    paste("SELECT transcript_index,(11+(6+gap)*i)::UBIGINT exon_start,(16+(6+gap)*i)::UBIGINT exon_end,",
      "(1+6*i)::UBIGINT exon_cdna_start,(6+6*i)::UBIGINT exon_cdna_end,",
      "0::TINYINT phase,0::TINYINT end_phase FROM route_tx,range(2) e(i) ORDER BY transcript_index,i"))
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('intronic_routes',",
    paste(dbQuoteString(con, route_queries), collapse = ","), ")"))$loaded)
  route_calls <- paste("SELECT transcript_index*2+e.i event_index,seq_region,",
    "CASE e.i WHEN 1 THEN 14 ELSE 22 END AS position,CASE e.i WHEN 1 THEN 'A' ELSE 'C' END AS reference,",
    "CASE e.i WHEN 1 THEN ['C'] ELSE ['T'] END alternates,transcript_index,s.i sample_index,",
    "CASE WHEN s.i=0 THEN CASE e.i WHEN 1 THEN '.' ELSE '0|1' END WHEN s.i=1 THEN '0|0'",
    "ELSE CASE e.i WHEN 1 THEN '1|1' ELSE '0|0' END END gt FROM route_tx,range(1,3) e(i),range(3) s(i)")
  routes <- rduckhts_haplotypes(con, route_calls, "intronic_routes", "vep116_compat",
    input_mode = "source_records")
  routes <- routes[vapply(routes$carriers, function(x) 0 %in% x$sample_index, TRUE), ]
  expect_equal(nrow(routes), 12L)
  expect_equal(routes$protein, c("MKP*", "MK*P", "MKP*")[routes$transcript_index %% 3 + 1L])
  expect_true(all(routes$edit_count == 0 & routes$sequence_flags == 0 &
    routes$sequence_status == "conditional" & routes$carrier_count == 1))
  expect_true(all(vapply(routes$protein_differences, nrow, 1L) == 0L))
  intronic_sources <- do.call(rbind, routes$contributors)
  intronic_sources <- intronic_sources[intronic_sources$position == 22, ]
  expect_equal(nrow(intronic_sources), 6L)
  expect_true(all(intronic_sources$projection_status == "outside_cds" &
    intronic_sources$reference == "C" & intronic_sources$alternate == "T" & intronic_sources$evidence_flags == 1L))
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
