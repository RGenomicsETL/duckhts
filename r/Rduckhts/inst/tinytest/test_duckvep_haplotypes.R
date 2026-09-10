library(tinytest)
library(DBI)

local({
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  tx <- paste("SELECT i::UINTEGER transcript_index,i::UINTEGER seq_region,",
    "100::UBIGINT transcript_start,111::UBIGINT transcript_end,(1-2*(i%2))::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,transcript_start cds_start,",
    "transcript_end cds_end,'AAAAAAAAAAAA'::BLOB cds_sequence,1::UTINYINT codon_table",
    "FROM range(8) t(i)")
  exons <- paste("SELECT i::UINTEGER transcript_index,100::UBIGINT exon_start,",
    "111::UBIGINT exon_end,1::UBIGINT exon_cdna_start,12::UBIGINT exon_cdna_end,",
    "0::TINYINT phase,0::TINYINT end_phase FROM range(8) t(i)")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('nominal',",
    "'SELECT i::UINTEGER seq_region FROM range(8) t(i)',", dbQuoteString(con, tx), ",",
    dbQuoteString(con, exons), ")"))$loaded)
  calls <- paste("WITH edits(scenario,ordinal,cds_start,ref,alt,gt) AS (VALUES",
    "(0,0,12,'A','AGC','.|1'),(0,1,12,'A','C','1|1'),",
    "(1,0,7,'AAAAAA','AAA','1|1'),(1,1,9,'AAAA','A','1|1'),",
    "(2,0,9,'AAAA','A','1|1'),(3,0,4,'A','AAA','1|1'),(3,1,9,'AA','A','1|1')),",
    "oriented AS (SELECT i transcript_index,i seq_region,e.ordinal,",
    "CASE WHEN i%2=0 THEN 99+cds_start ELSE 113-cds_start-len(ref) END AS position,",
    "CASE WHEN i%2=0 THEN ref ELSE reverse(translate(ref,'ACGT','TGCA')) END AS reference,",
    "[CASE WHEN i%2=0 THEN alt ELSE reverse(translate(alt,'ACGT','TGCA')) END] alternates,gt",
    "FROM range(8) t(i) JOIN edits e ON i//2=e.scenario)",
    "SELECT 2*transcript_index+row_number() OVER(PARTITION BY transcript_index",
    "ORDER BY position,ordinal) event_index,seq_region,position,reference,alternates,",
    "transcript_index,0 sample_index,gt FROM oriented ORDER BY event_index DESC")
  expected <- data.frame(transcript_index = rep(0:7, each = 2L), lane = rep(1:2, 8L),
    nominal = c(2, -1, 2, -1, rep(-6, 4L), rep(-3, 4L), rep(1, 4L)),
    net = c(2, 0, 2, 0, rep(-3, 8L), rep(1, 4L)),
    flags = c(rep(3L, 4L), rep(1L, 8L), rep(3L, 4L)))
  for (threads in c(1L, 4L)) {
    dbExecute(con, paste("SET threads =", threads))
    result <- rduckhts_haplotypes(con, calls, "nominal", phase_policy = "vep116_compat",
      input_mode = "source_records")
    expect_identical(tail(names(result), 1L), "nominal_length_diff")
    expect_equal(nrow(result), 10L)
    expect_equal(sum(result$carrier_count), 16L)
    actual <- do.call(rbind, lapply(seq_len(nrow(result)), function(i) {
      expect_equal(result$sequence_status[i], if (result$transcript_index[i] < 2L) "conditional" else "ok")
      expect_equal(nrow(result$contributors[[i]]), if (result$transcript_index[i] %in% 4:5) 1L else 2L)
      data.frame(transcript_index = result$transcript_index[i],
        lane = result$carriers[[i]]$haplotype_lane, nominal = result$nominal_length_diff[i],
        net = nchar(result$cds[i]) - 12, flags = result$sequence_flags[i])
    }))
    actual <- actual[order(actual$transcript_index, actual$lane), ]
    rownames(actual) <- NULL
    expect_equal(actual, expected)
    # Clipping changes actual component lengths, not the nominal source sum.
    block_net <- vapply(result$coding_blocks, function(x) sum(x$length_change), 0)
    expect_equal(sum(result$nominal_length_diff != block_net), 4L)
  }
})

local({
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  tx <- paste("SELECT i::UINTEGER transcript_index,i::UINTEGER seq_region,100::UBIGINT transcript_start,",
    "120::UBIGINT transcript_end,CASE WHEN i<3 THEN 1 ELSE -1 END::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,120::UBIGINT cds_end,",
    "(repeat('N',i%3)||'ATGGGTCCTGCTGAACAATAA')::BLOB cds_sequence,1::UTINYINT codon_table,",
    "''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence FROM range(6) t(i)")
  ex <- paste("SELECT i::UINTEGER transcript_index,100::UBIGINT exon_start,120::UBIGINT exon_end,",
    "1::UBIGINT exon_cdna_start,21::UBIGINT exon_cdna_end,(i%3)::TINYINT phase,",
    "0::TINYINT end_phase FROM range(6) t(i)")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('start_padding',",
    "'SELECT i::UINTEGER seq_region FROM range(6) t(i)',", dbQuoteString(con, tx), ",",
    dbQuoteString(con, ex), ")"))$loaded)
  calls <- paste("SELECT 2*i+j+1 AS event_index,i AS seq_region,",
    "CASE WHEN i<3 THEN CASE WHEN j=0 THEN 102 ELSE 109 END",
    "ELSE CASE WHEN j=0 THEN 118 ELSE 111 END END AS position,",
    "CASE WHEN i<3 THEN 'G' ELSE 'C' END AS reference,",
    "CASE WHEN i<3 THEN CASE WHEN j=0 THEN 'AAAA' ELSE 'T' END",
    "ELSE CASE WHEN j=0 THEN 'TTTT' ELSE 'A' END END AS alternate,",
    "1 alt_index,i transcript_index,0 sample_index,[1,1] alleles,[true,true] phase_before,",
    "NULL::BIGINT phase_set FROM range(6) t(i) CROSS JOIN range(2) e(j)")
  for (threads in c(1L, 4L)) {
    dbExecute(con, paste("SET threads =", threads))
    result <- rduckhts_haplotypes(con, calls, "start_padding", hgvs = TRUE)
    result <- result[order(result$transcript_index), ]
    expect_equal(result$transcript_index, 0:5)
    expect_equal(result$hgvsp, rep(c("p.[(Met1?;Ala4Ser)]", "p.[(Gly2?;Cys4Phe)]", "p.(Trp2?)"), 2))
    expect_true(all(result$hgvsp_status == "ok" & result$carrier_count == 2))
    expect_equal(result$cds, paste0(strrep("N", (0:5) %% 3L), "ATAAAAGGTCCTTCTGAACAATAA"))
    expect_equal(result$nominal_length_diff, rep(3, 6))
  }
})

local({
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  tx <- paste("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,",
    "100::UBIGINT transcript_start,111::UBIGINT transcript_end,1::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,100::UBIGINT cds_start,",
    "111::UBIGINT cds_end,'ATGGCTGCTTAA'::BLOB cds_sequence,1::UTINYINT codon_table,",
    "''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence")
  ex <- paste("SELECT 0::UINTEGER transcript_index,100::UBIGINT exon_start,111::UBIGINT exon_end,",
    "1::UBIGINT exon_cdna_start,12::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('partial_phase',",
    "'SELECT 0::UINTEGER seq_region',", dbQuoteString(con, tx), ",", dbQuoteString(con, ex), ")"))$loaded)
  path <- system.file("extdata", "geno_phase_partial.vcf", package = "Rduckhts")
  stopifnot(nzchar(path))
  calls <- paste0("SELECT i AS event_index,0 AS seq_region,POS AS position,REF AS reference,",
    "ALT[i] AS alternate,i AS alt_index,0 AS transcript_index,record_index AS sample_index,",
    "calls[1].alleles alleles,calls[1].phase_before phase_before,calls[1].phase_set phase_set ",
    "FROM read_geno(", dbQuoteString(con, path), ") CROSS JOIN range(1,3) alts(i) ",
    "WHERE record_index IN (3,5,6,8)")
  expected <- data.frame(sample = c(3L,3L,5L,5L,6L,6L,6L,8L,8L,8L),
    lane = c(1L,3L,1L,3L,1L,2L,3L,1L,2L,3L),
    cds = c(NA,NA,"ATGGATGCTTAA","ATGGGTGCTTAA",NA,"ATGGATGCTTAA",NA,
      NA,"ATGGATGCTTAA","ATGGGTGCTTAA"),
    evidence = c(4L,4L,1L,1L,6L,1L,6L,2L,1L,1L),
    contributors = c(2L,2L,1L,1L,2L,1L,2L,2L,1L,1L))
  for (threads in c(1L, 4L)) {
    dbExecute(con, paste("SET threads =", threads))
    result <- rduckhts_haplotypes(con, calls, "partial_phase")
    expect_equal(result$nominal_length_diff, ifelse(is.na(result$cds), NA_real_, 0))
    actual <- do.call(rbind, lapply(seq_len(nrow(result)), function(i) {
      carriers <- result$carriers[[i]]
      expect_true(all(carriers$phase_set == 10 & carriers$ploidy == 3))
      if (is.na(result$cds[i])) expect_equal(result$sequence_status[i], "incomplete_input")
      data.frame(sample = carriers$sample_index, lane = carriers$haplotype_lane,
        cds = result$cds[i], evidence = result$evidence_flags[i],
        contributors = nrow(result$contributors[[i]]))
    }))
    actual <- actual[order(actual$sample, actual$lane), ]
    rownames(actual) <- NULL
    expect_equal(actual, expected)
    expect_equal(sum(result$carrier_count), 10)
  }
})

local({
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  tx <- paste("SELECT i::UINTEGER transcript_index,0::UINTEGER seq_region,",
    "(100+20*i)::UBIGINT transcript_start,(111+20*i)::UBIGINT transcript_end,",
    "1::TINYINT strand,0::UINTEGER gene_index,3::UBIGINT transcript_flags,",
    "transcript_start AS cds_start,transcript_end AS cds_end,",
    "'AAAAAAAAAAAA'::BLOB cds_sequence,1::UTINYINT codon_table,",
    "''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence FROM range(2051) t(i)")
  ex <- paste("SELECT i::UINTEGER transcript_index,(100+20*i)::UBIGINT exon_start,",
    "(111+20*i)::UBIGINT exon_end,1::UBIGINT exon_cdna_start,12::UBIGINT exon_cdna_end,",
    "0::TINYINT phase,0::TINYINT end_phase FROM range(2051) t(i)")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('borrowed',",
    "'SELECT 0::UINTEGER seq_region',", dbQuoteString(con, tx), ",",
    dbQuoteString(con, ex), ")"))$loaded)
  calls <- paste("SELECT i+1 AS event_index,0 AS seq_region,100+20*i AS position,",
    "CASE WHEN i%2=0 THEN 'A' ELSE 'AAA' END AS reference,",
    "CASE WHEN i%2=0 THEN 'C' ELSE 'CAA' END AS alternate,1 AS alt_index,",
    "i AS transcript_index,s AS sample_index,[1,1] AS alleles,[true,true] AS phase_before,",
    "NULL::BIGINT AS phase_set FROM range(2051) t(i) CROSS JOIN range(3) samples(s)")
  for (threads in c(1L, 4L)) for (mode in c("alt_events", "source_records")) {
    dbExecute(con, paste("SET threads =", threads))
    query <- if (mode == "alt_events") calls else paste("SELECT event_index,seq_region,",
      "position,reference,[alternate] alternates,transcript_index,sample_index,'1|1' gt FROM (",
      calls, ")")
    # Raw input retains REF and undefined-slot descriptors as well as ALT.
    raw <- mode == "source_records"
    result <- rduckhts_haplotypes(con, query, "borrowed", hgvs = TRUE, input_mode = mode,
      phase_policy = if (raw) "vep116_compat" else "strict",
      max_active_events = if (raw) 3 else 2, max_active_transcripts = 1,
      max_active_projections = if (raw) 3 else 2, max_allele_bytes = if (raw) 24 else 8)
    expect_equal(nrow(result), 2051L)
    expect_equal(sort(result$transcript_index), 0:2050)
    expect_true(all(result$carrier_count == 6 & result$cds == "CAAAAAAAAAAA" &
      result$protein == "QKKK"))
    # Independent annotation uses uncertain-start HGVS and needs FASTA for the
    # padded allele. No missing value may be replaced by a guessed protein label.
    expect_identical(result$hgvsp, ifelse(result$transcript_index %% 2 == 0, "p.(Lys1?)", NA_character_))
    expect_identical(result$hgvsp_status, ifelse(result$transcript_index %% 2 == 0, "ok", "missing_reference"))
    expect_true(all(vapply(result$contributors, nrow, 0L) == 1L))
    provenance <- do.call(rbind, result$contributors)
    expect_equal(provenance$event_index, result$transcript_index + 1)
    expect_equal(provenance$position, 100 + 20 * result$transcript_index)
    expect_identical(provenance$reference, ifelse(result$transcript_index %% 2 == 0, "A", "AAA"))
    expect_identical(provenance$alternate, ifelse(result$transcript_index %% 2 == 0, "C", "CAA"))
    expect_true(all(vapply(result$carriers, function(x)
      identical(sort(as.integer(x$sample_index)), rep(0:2, each = 2L)), FALSE)))
  }
})

local({
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  tx <- paste("SELECT i::UINTEGER transcript_index,i::UINTEGER seq_region,",
    "100::UBIGINT transcript_start,125::UBIGINT transcript_end,",
    "CASE WHEN i<3 THEN 1 ELSE -1 END::TINYINT strand,0::UINTEGER gene_index,",
    "3::UBIGINT transcript_flags,CASE WHEN i<3 THEN 110 ELSE 104 END::UBIGINT cds_start,",
    "CASE WHEN i<3 THEN 121 ELSE 115 END::UBIGINT cds_end,",
    "(repeat('N',i%3)||'ATGGCTGCTTAA')::BLOB cds_sequence,1::UTINYINT codon_table,",
    "'AAAAAA'::BLOB pre_cds_sequence,'AAAA'::BLOB post_cds_sequence FROM range(6) t(i)")
  ex <- paste("SELECT i::UINTEGER transcript_index,",
    "CASE WHEN i<3 THEN CASE WHEN e=0 THEN 100 ELSE 110 END ELSE CASE WHEN e=0 THEN 120 ELSE 100 END END::UBIGINT exon_start,",
    "CASE WHEN i<3 THEN CASE WHEN e=0 THEN 105 ELSE 125 END ELSE CASE WHEN e=0 THEN 125 ELSE 115 END END::UBIGINT exon_end,",
    "CASE WHEN e=0 THEN 1 ELSE 7 END::UBIGINT exon_cdna_start,",
    "CASE WHEN e=0 THEN 6 ELSE 22 END::UBIGINT exon_cdna_end,",
    "CASE WHEN e=0 THEN -1 ELSE i%3 END::TINYINT phase,-1::TINYINT end_phase",
    "FROM range(6) t(i) CROSS JOIN range(2) ex(e) ORDER BY i,e")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('phase',",
    "'SELECT i::UINTEGER seq_region FROM range(6) t(i)',", dbQuoteString(con, tx), ",",
    dbQuoteString(con, ex), ")"))$loaded)
  dbExecute(con, paste("CREATE TABLE phase_events AS SELECT i+1 event_index,i seq_region,",
    "CASE WHEN i<3 THEN 114 ELSE 111 END AS position,CASE WHEN i<3 THEN 'C' ELSE 'G' END AS reference,",
    "CASE WHEN i<3 THEN 'A' ELSE 'T' END AS alternate,NULL::UBIGINT end_position,",
    "NULL::VARCHAR structural_type,NULL::VARCHAR copy_change,NULL::UINTEGER mate_seq_region,",
    "NULL::UBIGINT mate_position FROM range(6) t(i)"))
  expected <- dbGetQuery(con, paste("SELECT event_index,protein_hgvs FROM duckvep_annotate(",
    "'phase_events','phase',hgvs:=true,upstream_distance:=0,downstream_distance:=0) ORDER BY event_index"))
  expect_equal(nrow(expected), 6L)
  expect_true(any(!is.na(expected$protein_hgvs)))
  for (threads in c(1L, 4L)) for (mode in c("alt_events", "source_records")) {
    dbExecute(con, paste("SET threads =", threads))
    calls <- if (mode == "alt_events") paste("SELECT *,event_index-1 transcript_index,",
      "0 sample_index,1 alt_index,[1,1] alleles,[true,true] phase_before,",
      "NULL::BIGINT phase_set FROM phase_events") else paste("SELECT *,[alternate] alternates,",
      "event_index-1 transcript_index,0 sample_index,'1|1' gt FROM phase_events")
    result <- rduckhts_haplotypes(con, calls, "phase", hgvs = TRUE, input_mode = mode,
      phase_policy = if (mode == "source_records") "vep116_compat" else "strict")
    result <- result[order(result$transcript_index), ]
    expect_equal(nrow(result), 6L)
    expect_equal(result$transcript_index, 0:5)
    expect_identical(gsub("[()]", "", result$hgvsp), expected$protein_hgvs)
    expect_true(all(result$cds == paste0(strrep("N", (0:5) %% 3), "ATGGATGCTTAA")))
    expect_true(all(result$carrier_count == 2 & vapply(result$contributors, nrow, 0L) == 1L))
    expect_equal(do.call(rbind, result$contributors)$event_index, expected$event_index)
  }
})

local({
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  tx <- paste("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,",
    "11::UBIGINT transcript_start,45::UBIGINT transcript_end,1::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,11::UBIGINT cds_start,",
    "22::UBIGINT cds_end,'ATGGGTCCTTAA'::BLOB cds_sequence,1::UTINYINT codon_table,",
    "''::BLOB pre_cds_sequence,'AAAGAACAATAATAACTAGCTGA'::BLOB post_cds_sequence")
  ex <- paste("SELECT 0::UINTEGER transcript_index,11::UBIGINT exon_start,45::UBIGINT exon_end,",
    "1::UBIGINT exon_cdna_start,35::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase")
  regions <- "SELECT 0::UINTEGER seq_region,55::UBIGINT sequence_length,'chrA1'::VARCHAR seq_region_name"
  reference <- system.file("extdata", "duckvep_hgvs_anchor.fa", package = "Rduckhts")
  expect_true(nzchar(reference))
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('anchor',",
    dbQuoteString(con, regions), ",", dbQuoteString(con, tx), ",",
    dbQuoteString(con, ex), ",reference_fasta:=", dbQuoteString(con, reference), ")"))$loaded)
  calls <- paste("SELECT event_index,0 seq_region,position,reference,alternate,1 alt_index,",
    "0 transcript_index,event_index-1 sample_index,[1,1] alleles,",
    "[true,true] phase_before,NULL::BIGINT phase_set FROM",
    "(VALUES (1,19,'T','TT'),(2,20,'T','TT'),(3,20,'T','TT'),(4,21,'A','TA'))",
    "v(event_index,position,reference,alternate)")
  inputs <- dbGetQuery(con, calls)
  # These are per-record executable VEP-116 suffixes, not HGVS-rule corrections.
  expected <- c("p.(Ter4LeufsTer9)", rep("p.(Ter4delinsLeuTer)", 3L))
  baseline <- list()
  for (threads in c(1L, 4L)) for (route in c("strict", "vep116_compat", "source_records")) {
    dbExecute(con, paste("SET threads =", threads))
    raw <- route == "source_records"
    query <- if (raw) paste("SELECT event_index,seq_region,position,reference,",
      "[alternate] alternates,transcript_index,sample_index,'1|1' gt FROM (", calls, ")") else calls
    # Independent record comparisons must not trigger raw-file duplicate retention.
    queries <- if (raw) paste("SELECT * FROM (", query, ") WHERE event_index=", 1:4) else query
    result <- do.call(rbind, lapply(queries, function(sql) rduckhts_haplotypes(con, sql, "anchor",
      phase_policy = if (raw) "vep116_compat" else route,
      input_mode = if (raw) "source_records" else "alt_events", hgvs = TRUE)))
    expect_equal(nrow(result), 4L)
    ids <- vapply(result$contributors, function(x) x$event_index, 0)
    expect_equal(sort(ids), 1:4)
    result <- result[order(ids), ]
    rownames(result) <- NULL
    expect_identical(result$hgvsp, expected)
    expect_true(all(result$hgvsp_status == "ok"))
    expect_true(all(result$cds == "ATGGGTCCTTTAA" & result$protein == "MGPL"))
    expect_true(all(result$carrier_count == 2L))
    for (i in 1:4) {
      expect_equal(nrow(result$contributors[[i]]), 1L)
      expect_equal(result$contributors[[i]]$position, inputs$position[i])
      expect_identical(result$contributors[[i]]$reference, inputs$reference[i])
      expect_identical(result$contributors[[i]]$alternate, inputs$alternate[i])
      expect_equal(result$carriers[[i]]$sample_index, rep(i-1L, 2L))
    }
    if (is.null(baseline[[route]])) baseline[[route]] <- result else
      expect_identical(result, baseline[[route]])
  }
  expect_error(rduckhts_haplotypes(con, calls, "anchor", hgvs = TRUE,
    max_hgvs_reference_bytes = 56), pattern = "reference workspace bytes=56, required=57")
  bounded <- rduckhts_haplotypes(con, calls, "anchor", hgvs = TRUE, max_hgvs_reference_bytes = 57)
  expect_equal(sum(!is.na(bounded$hgvsp)), 4L)
  disabled <- rduckhts_haplotypes(con, calls, "anchor", hgvs = FALSE, max_hgvs_reference_bytes = 1)
  expect_true(all(is.na(disabled$hgvsp) & disabled$hgvsp_status == "not_requested"))
})

local({
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  tx <- paste("SELECT 0::UINTEGER transcript_index,0::UINTEGER seq_region,",
    "101::UBIGINT transcript_start,136::UBIGINT transcript_end,1::TINYINT strand,",
    "0::UINTEGER gene_index,3::UBIGINT transcript_flags,101::UBIGINT cds_start,",
    "136::UBIGINT cds_end,('ATG'||repeat('CAG',10)||'TAA')::BLOB cds_sequence,",
    "1::UTINYINT codon_table,''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence")
  ex <- paste("SELECT 0::UINTEGER transcript_index,101::UBIGINT exon_start,",
    "136::UBIGINT exon_end,1::UBIGINT exon_cdna_start,36::UBIGINT exon_cdna_end,",
    "0::TINYINT phase,0::TINYINT end_phase")
  reference <- system.file("extdata", "duckvep_repeat.fa", package = "Rduckhts")
  expect_true(nzchar(reference))
  regions <- "SELECT 0::UINTEGER seq_region,236::UBIGINT sequence_length,'chr1'::VARCHAR seq_region_name"
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('repeat',",
    dbQuoteString(con, regions), ",", dbQuoteString(con, tx), ",",
    dbQuoteString(con, ex), ",reference_fasta:=", dbQuoteString(con, reference), ")"))$loaded)
  # A VCF repeat summary may omit interruptions. Literal alleles determine replay.
  calls <- paste("SELECT event_index,0 seq_region,103 AS position,",
    "'G'||repeat('CAG',10) AS reference,'G'||sequence AS alternate,1 alt_index,",
    "0 transcript_index,event_index-1 sample_index,[1,1] alleles,[true,true] phase_before,",
    "NULL::BIGINT phase_set,'CAG' repeat_unit,11 repeat_count FROM",
    "(VALUES (1,repeat('CAG',11)),(2,repeat('CAG',5)||'CAT'||repeat('CAG',5)))",
    "r(event_index,sequence)")
  exact <- dbGetQuery(con, calls)
  prepared <- paste("SELECT * REPLACE('G'||(duckvep_repeat_sequence(",
    "[{unit:'CAG',count:10}],true)).sequence AS reference,",
    "'G'||(duckvep_repeat_sequence(CASE event_index WHEN 1 THEN [{unit:'CAG',count:11}]",
    "ELSE [{unit:'CAG',count:5},{unit:'CAT',count:1},{unit:'CAG',count:5}] END,true,",
    "max_sequence_bases:=33)).sequence AS alternate) FROM (", calls, ")")
  expect_identical(dbGetQuery(con, prepared), exact)
  annotations <- lapply(list(literal = calls, prepared = prepared), function(query) {
    dbExecute(con, paste("CREATE OR REPLACE TABLE repeat_events AS SELECT event_index,",
      "seq_region,position,reference,alternate,NULL::UBIGINT end_position,",
      "NULL::VARCHAR structural_type,NULL::VARCHAR copy_change,NULL::UINTEGER mate_seq_region,",
      "NULL::UBIGINT mate_position FROM (", query, ")"))
    dbGetQuery(con, paste("SELECT * FROM duckvep_annotate('repeat_events','repeat',",
      "hgvs:=true,upstream_distance:=0,downstream_distance:=0) ORDER BY event_index"))
  })
  expect_equal(nrow(annotations$prepared), 2L)
  expect_identical(annotations$prepared, annotations$literal)
  expect_true(all(!is.na(annotations$prepared$protein_hgvs)))
  for (threads in c(1L, 4L)) for (mode in c("alt_events", "source_records")) {
    dbExecute(con, paste("SET threads =", threads))
    query <- if (mode == "alt_events") calls else paste(
      "SELECT event_index,seq_region,position,reference,[alternate] alternates,",
      "transcript_index,sample_index,'1|1' gt FROM (", calls, ")")
    result <- rduckhts_haplotypes(con, query, "repeat",
      phase_policy = if (mode == "alt_events") "strict" else "vep116_compat",
      input_mode = mode, hgvs = TRUE)
    expect_equal(nrow(result), 2L)
    result <- result[order(vapply(result$contributors, function(x) x$event_index, 0)), ]
    prepared_query <- if (mode == "alt_events") prepared else paste(
      "SELECT event_index,seq_region,position,reference,[alternate] alternates,",
      "transcript_index,sample_index,'1|1' gt FROM (", prepared, ")")
    expanded <- rduckhts_haplotypes(con, prepared_query, "repeat",
      phase_policy = if (mode == "alt_events") "strict" else "vep116_compat",
      input_mode = mode, hgvs = TRUE)
    expanded <- expanded[order(vapply(expanded$contributors, function(x) x$event_index, 0)), ]
    rownames(expanded) <- rownames(result) <- NULL
    expect_identical(expanded, result)
    expect_identical(result$protein, c("MQQQQQQQQQQQ*", "MQQQQQHQQQQQ*"))
    expect_identical(result$hgvsp, c("p.(Gln11dup)", "p.(Gln6_Gln7insHis)"))
    expect_identical(result$hgvsp_status, c("ok", "ok"))
    expect_equal(result$carrier_count, c(2L, 2L))
    expect_identical(result$cds, paste0("AT", exact$alternate, "TAA"))
    for (i in 1:2) {
      expect_equal(result$contributors[[i]]$event_index, i)
      expect_identical(result$contributors[[i]]$reference, exact$reference[i])
      expect_identical(result$contributors[[i]]$alternate, exact$alternate[i])
      expect_equal(result$carriers[[i]]$sample_index, rep(i - 1L, 2L))
    }
  }
})

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
  expect_true(all(is.na(actual$hgvsp) & actual$hgvsp_status == "not_requested"))
  for (invalid in list(NA, 1, "true", logical()))
    expect_error(rduckhts_haplotypes(con, calls, "haps", hgvs = invalid), pattern = "hgvs must")
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
  expect_equal(composed$nominal_length_diff, c(0, 0))
  expect_true(all(is.na(composed$stop_in_displaced_frame)))
  for (i in 1:2) {
    expect_equal(composed$coding_blocks[[i]]$coding_status, "unsupported_ordered_replacements")
    expect_true(is.na(composed$coding_blocks[[i]]$local_consequence_mask))
    expect_equal(composed$coding_blocks[[i]]$event_indices[[1]], c(1, 2))
    expect_equal(composed$contributors[[i]]$alt_index, c(i - 1L, 1L))
    expect_equal(composed$contributors[[i]]$evidence_flags, c(1L, 1L))
  }
  uncertain <- rduckhts_haplotypes(con, sub("0|1", ".|.", overlapping, fixed = TRUE),
    "haps", "vep116_compat", input_mode = "source_records", hgvs = TRUE)
  expect_equal(uncertain$cds, "AGAAAAAAAAAA")
  expect_equal(uncertain$sequence_status, "conditional")
  expect_equal(uncertain$edit_count, 1)
  expect_equal(uncertain$carrier_count, 2L)
  expect_equal(uncertain$nominal_length_diff, 0)
  expect_true(all(is.na(uncertain$hgvsp) & uncertain$hgvsp_status == "incomplete_input"))
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
  missense_mask <- dbGetQuery(con,
    "SELECT consequence_mask FROM duckvep_so_terms() WHERE consequence='missense_variant'")$consequence_mask
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
    expect_equal(n_changed$coding_blocks[[1]]$coding_status, "ok")
    expect_equal(n_changed$coding_blocks[[1]]$local_consequence_mask, missense_mask)
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
  restore_tx <- paste("SELECT i::UINTEGER transcript_index,0::UINTEGER seq_region,",
    "(100+100*i)::UBIGINT transcript_start,(135+100*i)::UBIGINT transcript_end,",
    "(1-2*i)::TINYINT strand,0::UINTEGER gene_index,3::UBIGINT transcript_flags,",
    "transcript_start cds_start,transcript_end cds_end,",
    "'ATGGGTGGTGCTGATGATGCTGATGCTGATGGTTAA'::BLOB cds_sequence,1::UTINYINT codon_table",
    "FROM range(2) r(i)")
  restore_exons <- paste("SELECT transcript_index,transcript_start exon_start,",
    "transcript_end exon_end,1::UBIGINT exon_cdna_start,36::UBIGINT exon_cdna_end,",
    "0::TINYINT phase,0::TINYINT end_phase FROM (", restore_tx, ")")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('restore',",
    dbQuoteString(con, "SELECT 0::UINTEGER seq_region"), ",",
    dbQuoteString(con, restore_tx), ",", dbQuoteString(con, restore_exons), ")"))$loaded)
  restore_calls <- paste("SELECT event_index,0 seq_region,position,reference,alternate,",
    "1 alt_index,transcript_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set",
    "FROM (VALUES (1,0,102,'GGGT','G'),(2,0,109,'GCT','GGT'),(3,0,111,'T','TGCT'),",
    "(11,1,229,'CACC','C'),(12,1,224,'AGC','ACC'),(13,1,223,'C','CAGC'))",
    "v(event_index,transcript_index,position,reference,alternate) ORDER BY position DESC")
  for (policy in c("strict", "vep116_compat")) {
    result <- rduckhts_haplotypes(con, restore_calls, "restore", policy, max_leaf_edits = 3)
    expect_equal(nrow(result), 2L)
    expect_equal(result$cds, rep("ATGGGTGGTGCTGATGATGCTGATGCTGATGGTTAA", 2L))
    expect_equal(result$protein, rep("MGGADDADADG*", 2L))
    expect_equal(result$edit_count, c(3, 3))
    expect_equal(result$sequence_flags, c(1L, 1L))
    expect_true(all(!result$stop_in_displaced_frame))
    for (i in seq_len(nrow(result))) {
      blocks <- result$coding_blocks[[i]]
      expect_equal(blocks$cds_start, c(4, 11, 13))
      expect_equal(blocks$coding_status, rep("ok", 3L))
      expect_equal(blocks$local_consequence_mask[2], missense_mask)
      expected_events <- 10 * result$transcript_index[i] + 1:3
      expect_equal(unlist(blocks$event_indices), expected_events)
      expect_equal(sort(result$contributors[[i]]$event_index), expected_events)
      expect_equal(nrow(result$cds_differences[[i]]), 0L)
      expect_equal(nrow(result$protein_differences[[i]]), 0L)
    }
  }
  expect_error(rduckhts_haplotypes(con, restore_calls, "restore", max_leaf_edits = 2),
    pattern = "max_leaf_edits")
  # Finalize R-owned failed statements before checking model release.
  gc()
  expect_true(dbGetQuery(con, "SELECT duckvep_model_drop('restore') dropped")$dropped)
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

  views_tx <- paste("SELECT i::UINTEGER transcript_index,0::UINTEGER seq_region,",
    "(100+20*i)::UBIGINT transcript_start,(108+20*i)::UBIGINT transcript_end,",
    "1::TINYINT strand,0::UINTEGER gene_index,3::UBIGINT transcript_flags,",
    "transcript_start cds_start,transcript_end cds_end,cds::BLOB cds_sequence,1::UTINYINT codon_table",
    "FROM (VALUES (0,'ATGGCNTAA'),(1,'ATGGCCTAA')) v(i,cds)")
  views_exons <- paste("SELECT transcript_index,transcript_start exon_start,transcript_end exon_end,",
    "1::UBIGINT exon_cdna_start,9::UBIGINT exon_cdna_end,",
    "0::TINYINT phase,0::TINYINT end_phase FROM (", views_tx, ")")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('reference_views',",
    dbQuoteString(con, "SELECT 0::UINTEGER seq_region"), ",", dbQuoteString(con, views_tx), ",",
    dbQuoteString(con, views_exons), ")"))$loaded)
  views_calls <- paste("SELECT transcript_index+1 event_index,0 seq_region,",
    "transcript_start+3 AS position,'G' reference,'A' alternate,1 alt_index,transcript_index,",
    "0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set FROM (", views_tx, ")")
  for (workers in c(1L, 4L)) {
    dbExecute(con, paste("SET threads=", workers))
    for (policy in c("strict", "vep116_compat")) {
      views <- rduckhts_haplotypes(con, views_calls, "reference_views", phase_policy = policy)
      views <- views[order(views$transcript_index), ]
      expect_equal(views$cds, c("ATGACNTAA", "ATGACCTAA"))
      expect_equal(views$protein, rep("MT*", 2L))
      expect_equal(views$carrier_count, rep(1, 2L))
      expect_equal(vapply(views$coding_blocks, function(x) x$coding_status, ""),
        c("ok", "ok"))
      for (difference in views$protein_differences) {
        expect_equal(difference$reference, "A")
        expect_equal(difference$alternate, "T")
        expect_equal(difference$ref_start0, 1)
        expect_equal(difference$alt_start0, 1)
      }
    }
  }

  # Original VEP-116 witnesses: surrounding N, first/terminal codons and uploaded REF N.
  # Unknown-residue HGVSp uses Ter; these expectations do not infer a stop consequence.
  n_codon_source <- data.frame(i = 0:10,
    cds = c("ATGGCNTAA", "ATGAANTAA", "ATGNCNTAA", "GCNGCCTAA", "AANGCCTAA",
      "CTNGCCTAA", "ATGGCN", "ATGTAN", "ATGTCN", "GCNGCCTAA", "ATNGCCTAA"),
    code = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 29L, 22L, 1L, 1L),
    position = c(4L, 5L, 4L, 2L, 2L, 1L, 4L, 5L, 5L, 1L, 1L),
    reference = c("G", "A", "N", "C", "A", "C", "G", "A", "C", "G", "A"),
    alternate = c("A", "C", "A", "A", "C", "A", "A", "C", "A", "A", "C"))
  dbWriteTable(con, "n_codon_source", n_codon_source)
  n_codon_tx <- paste("SELECT i::UINTEGER transcript_index,i::UINTEGER seq_region,",
    "1::UBIGINT transcript_start,length(cds)::UBIGINT transcript_end,1::TINYINT strand,",
    "i::UINTEGER gene_index,3::UBIGINT transcript_flags,1::UBIGINT cds_start,",
    "transcript_end cds_end,cds::BLOB cds_sequence,code::UTINYINT codon_table",
    "FROM n_codon_source ORDER BY transcript_index")
  n_codon_exons <- paste("SELECT i::UINTEGER transcript_index,1::UBIGINT exon_start,",
    "length(cds)::UBIGINT exon_end,1::UBIGINT exon_cdna_start,",
    "length(cds)::UBIGINT exon_cdna_end,0::TINYINT phase,0::TINYINT end_phase",
    "FROM n_codon_source ORDER BY transcript_index")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('n_codon_witnesses',",
    "'SELECT i::UINTEGER seq_region FROM n_codon_source ORDER BY i',",
    dbQuoteString(con, n_codon_tx), ",", dbQuoteString(con, n_codon_exons), ")"))$loaded)
  n_codon_calls <- paste("SELECT i::UBIGINT event_index,i::UINTEGER seq_region,",
    "position::UBIGINT AS position,reference,alternate,1::UINTEGER alt_index,",
    "i::UINTEGER transcript_index,0::UINTEGER sample_index,[1]::INTEGER[] alleles,",
    "[true]::BOOLEAN[] phase_before,NULL::BIGINT phase_set FROM n_codon_source")
  dbExecute(con, paste("CREATE TABLE n_codon_events AS SELECT event_index,seq_region,",
    "position,reference,alternate,NULL::UBIGINT end_position,NULL::VARCHAR structural_type,",
    "NULL::VARCHAR copy_change,NULL::UINTEGER mate_seq_region,NULL::UBIGINT mate_position",
    "FROM (", n_codon_calls, ")"))
  n_codon_expected_so <- c("missense_variant", "coding_sequence_variant&missense_variant",
    "coding_sequence_variant", "coding_sequence_variant&missense_variant", "start_lost",
    "coding_sequence_variant&missense_variant", "missense_variant", "missense_variant",
    "coding_sequence_variant", "start_lost", "start_lost")
  n_codon_expected_independent <- c("p.Ala2Thr", "p.Ter2Thr", NA_character_, "p.Ala1Ter",
    "p.Ter1?", "p.Leu1Ter", "p.Ala2Thr", "p.Tyr2Ser", "p.Ter2=", "p.Ala1?", "p.Ter1?")
  n_codon_expected_singleton <- c("p.(Ala2Thr)", "p.(Ter2Thr)", NA_character_, "p.(Ala1Ter)",
    "p.(Ter1?)", "p.(Leu1Ter)", "p.(Ala2Thr)", "p.(Tyr2Ser)", "p.(Ter2=)",
    "p.(Ala1?)", "p.(Ter1?)")
  n_codon_independent <- dbGetQuery(con, paste("SELECT event_index,",
    "(SELECT string_agg(t.consequence,'&' ORDER BY t.consequence) FROM duckvep_so_terms() t",
    "WHERE (a.consequence_mask & t.consequence_mask)<>0) consequences,protein_hgvs",
    "FROM duckvep_annotate('n_codon_events','n_codon_witnesses',hgvs:=true,",
    "upstream_distance:=0,downstream_distance:=0) a ORDER BY event_index"))
  expect_equal(n_codon_independent$event_index, n_codon_source$i)
  expect_identical(n_codon_independent$consequences, n_codon_expected_so)
  expect_identical(n_codon_independent$protein_hgvs, n_codon_expected_independent)
  for (policy in c("strict", "vep116_compat")) {
    n_codon_singletons <- rduckhts_haplotypes(con, n_codon_calls, "n_codon_witnesses",
      phase_policy = policy, hgvs = TRUE)
    n_codon_singletons <- n_codon_singletons[order(n_codon_singletons$transcript_index), ]
    expect_equal(n_codon_singletons$transcript_index, n_codon_source$i)
    expect_identical(n_codon_singletons$hgvsp, n_codon_expected_singleton)
    expect_equal(n_codon_singletons$carrier_count, rep(1, nrow(n_codon_source)))
    expect_equal(vapply(n_codon_singletons$contributors, nrow, 0L), rep(1L, nrow(n_codon_source)))
    n_codon_contributors <- do.call(rbind, n_codon_singletons$contributors)
    expect_equal(n_codon_contributors$event_index, n_codon_source$i)
    expect_equal(n_codon_contributors$position, n_codon_source$position)
    expect_identical(n_codon_contributors$reference, n_codon_source$reference)
    expect_identical(n_codon_contributors$alternate, n_codon_source$alternate)
  }

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
    p_hgvs <- rduckhts_haplotypes(con, p_calls, "reference_proteins",
      phase_policy = policy, hgvs = TRUE)
    p_hgvs <- p_hgvs[order(p_hgvs$transcript_index), ]
    # Pinned TranscriptVariationAllele observes each synonymous ALT separately
    # from Haplosaurus's contrast against the curated reference protein.
    expect_equal(p_hgvs$hgvsp, c("p.(Ala2=)", "p.(Ala3=)", rep("p.(Ala2=)", 3L), NA_character_))
    expect_equal(p_hgvs$hgvsp_status, c(rep("ok", 5L), "missing_reference_protein"))
    fields <- setdiff(names(p), c("hgvsp", "hgvsp_status"))
    expect_equal(p_hgvs[fields], p[fields])
  }
  p_equal <- rduckhts_haplotypes(con, paste("SELECT * REPLACE ('CTGG' AS alternate) FROM (",
    p_calls, ") WHERE transcript_index=2"), "reference_proteins", hgvs = TRUE)
  expect_equal(p_equal$protein, "MAW*")
  expect_true(is.na(p_equal$hgvsp))
  expect_equal(p_equal$hgvsp_status, "missing_reference")
  expect_equal(p_equal$edit_count, 1)
  expect_equal(nrow(p_equal$coding_blocks[[1L]]), 1L)
  expect_equal(nrow(p_equal$contributors[[1L]]), 1L)
  p_terminal <- rduckhts_haplotypes(con, paste("SELECT * REPLACE ('CTGGGCC' AS alternate) FROM (",
    p_calls, ") WHERE transcript_index=2"), "reference_proteins", hgvs = TRUE)
  expect_equal(p_terminal$protein, "MAWA*")
  expect_true(is.na(p_terminal$hgvsp))
  expect_equal(p_terminal$hgvsp_status, "missing_reference")
  expect_equal(p_terminal$edit_count, 1)
  expect_equal(nrow(p_terminal$coding_blocks[[1L]]), 1L)
  expect_equal(nrow(p_terminal$contributors[[1L]]), 1L)
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
    "reference_proteins", "vep116_compat", input_mode = "source_records", hgvs = TRUE)
  p_retained <- p_retained[vapply(p_retained$carriers, function(x) 1L %in% x$haplotype_lane, TRUE), ]
  p_retained <- p_retained[order(p_retained$transcript_index), ]
  expect_equal(p_retained$protein, c("LA*", "M*", "MA*", "MAW", "LA*", ""))
  expect_true(all(p_retained$edit_count == 0))
  expect_equal(p_retained$hgvsp[1:3], c("p.(Met1Leu)", "p.(Sec2Ter)", "p.(Trp3Ter)"))
  expect_equal(p_retained$hgvsp_status[1:3], rep("ok", 3L))
  expect_equal(nrow(p_retained$coding_blocks[[1L]]), 0L)
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

local({
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
  dbExecute(con, paste(
    "CREATE TABLE protein_tx AS SELECT i::UINTEGER transcript_index,0::UINTEGER seq_region,",
    "(100+100*i)::UBIGINT transcript_start,(99+100*i+len(cds))::UBIGINT transcript_end,",
    "(CASE WHEN i<3 THEN 1 ELSE -1 END)::TINYINT strand,0::UINTEGER gene_index,",
    "3::UBIGINT transcript_flags,transcript_start cds_start,transcript_end cds_end,",
    "cds::BLOB cds_sequence,1::UTINYINT codon_table,''::BLOB pre_cds_sequence,''::BLOB post_cds_sequence",
    "FROM range(6) r(i) JOIN (VALUES (0,'ATGGCTGCTGCTGCTGCTGAATAA'),",
    "(1,'ATGCGGCATTTCTATGAATAA'),(2,'ATGGGTCCTGCTGAACAATAA')) c(k,cds) ON i%3=k"))
  exons <- paste("SELECT transcript_index,transcript_start exon_start,transcript_end exon_end,",
    "1::UBIGINT exon_cdna_start,(transcript_end-transcript_start+1)::UBIGINT exon_cdna_end,",
    "0::TINYINT phase,0::TINYINT end_phase FROM protein_tx ORDER BY transcript_index")
  expect_true(dbGetQuery(con, paste0("SELECT loaded FROM duckvep_model_load('protein',",
    dbQuoteString(con, "SELECT 0::UINTEGER seq_region"), ",",
    dbQuoteString(con, "SELECT * FROM protein_tx ORDER BY transcript_index"), ",",
    dbQuoteString(con, exons), ")"))$loaded)
  dbExecute(con, paste(
    "CREATE TABLE protein_calls AS SELECT 10*transcript_index+j AS event_index,0 seq_region,",
    "CASE WHEN strand=1 THEN transcript_start+p-1 ELSE transcript_end-p-len(ref)+2 END AS position,",
    "CASE WHEN strand=1 THEN ref ELSE seq_revcomp(ref) END AS reference,",
    "CASE WHEN strand=1 THEN alt ELSE seq_revcomp(alt) END AS alternate,1 alt_index,",
    "transcript_index,0 sample_index,[1] alleles,[true] phase_before,NULL::BIGINT phase_set",
    "FROM protein_tx JOIN (VALUES (0,1,3,'G','GGCT'),(0,2,11,'C','A'),",
    "(1,1,3,'G','GC'),(1,2,12,'CT','C'),(2,1,3,'GGGT','G'),(2,2,9,'T','TT'))",
    "v(k,j,p,ref,alt) ON transcript_index%3=k"))
  expected <- rep(c("p.(Ala4_Ala5insAsp)", "p.[(Arg2_His3delinsProAla;Tyr5His)]",
    "p.[(Gly2del;Ala4CysfsTer2)]"), 2L)
  for (policy in c("strict", "vep116_compat")) {
    result <- rduckhts_haplotypes(con, "SELECT * FROM protein_calls ORDER BY position DESC",
      "protein", phase_policy = policy, hgvs = TRUE)
    result <- result[order(result$transcript_index), ]
    expect_identical(result$hgvsp, expected)
    expect_identical(result$hgvsp_status, rep("ok", 6L))
    expect_identical(result$protein, rep(c("MAAADAAE*", "MPAFHE*", "MPC*"), 2L))
    expect_equal(result$edit_count, rep(2, 6L))
    expect_equal(vapply(result$contributors, nrow, 1L), rep(2L, 6L))
    # Codon-aligned deletion plus insertion has different physical grouping
    # on the two source strands, without changing the protein operation set.
    expect_equal(vapply(result$coding_blocks, nrow, 1L)[c(3L, 6L)], c(2L, 1L))
  }
  one <- "SELECT * FROM protein_calls WHERE transcript_index=0"
  source_calls <- paste("SELECT event_index,seq_region,position,reference,[alternate] alternates,",
    "transcript_index,sample_index,'1|1' gt FROM protein_calls")
  raw <- rduckhts_haplotypes(con, source_calls, "protein", "vep116_compat",
    input_mode = "source_records", hgvs = TRUE)
  raw <- raw[order(raw$transcript_index), ]
  expect_identical(raw$hgvsp, expected)
  expect_identical(raw$hgvsp_status, rep("ok", 6L))
  expect_equal(raw$carrier_count, rep(2L, 6L))
  expect_error(rduckhts_haplotypes(con, one, "protein", hgvs = TRUE, max_hgvs_bytes = 18),
    pattern = "max_hgvs_bytes=18, required=19")
  exact <- rduckhts_haplotypes(con, one, "protein", hgvs = TRUE,
    max_hgvs_bytes = 19, max_hgvs_operations = 2)
  expect_identical(exact$hgvsp, expected[1L])
  expect_error(rduckhts_haplotypes(con, "SELECT * FROM protein_calls WHERE transcript_index=1",
    "protein", hgvs = TRUE, max_hgvs_operations = 1), pattern = "max_hgvs_operations=1")
  expect_error(rduckhts_haplotypes(con, one, "protein", hgvs = TRUE, workspace_limit = 1),
    pattern = "workspace")
  rduckhts_haplotypes(con, one, "protein", hgvs = TRUE, table_name = "protein_output")
  expect_identical(dbGetQuery(con, "SELECT hgvsp FROM protein_output")$hgvsp, expected[1L])
  expect_equal(dbGetQuery(con, "SELECT 42 n")$n, 42L)
})
