# Sequence encoding and UDF tests for Rduckhts
library(tinytest)
library(DBI)

test_seq_ops <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  expect_silent(rduckhts_load(con))

  bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
  fasta_path <- system.file("extdata", "ce.fa", package = "Rduckhts")
  fastq_r1 <- system.file("extdata", "r1.fq", package = "Rduckhts")
  fasta_index_path <- tempfile("duckhts_seq_ops_", fileext = ".fai")
  on.exit(unlink(fasta_index_path), add = TRUE)
  expect_silent(rduckhts_fasta_index(con, fasta_path, index_path = fasta_index_path))

  # =========================================================================
  # seq_encode_4bit / seq_decode_4bit UDF tests
  # =========================================================================

  # roundtrip: all IUPAC bases
  rt <- DBI::dbGetQuery(con,
    "SELECT seq_decode_4bit(seq_encode_4bit('ACGTRYSWKMBDHVN')) AS seq")
  expect_equal(rt$seq[1], "ACGTRYSWKMBDHVN")

  # U (RNA) normalizes to T (code 8)
  u_enc <- DBI::dbGetQuery(con,
    "SELECT seq_encode_4bit('ACGU')::VARCHAR AS codes")
  expect_equal(u_enc$codes[1], "[1, 2, 4, 8]")

  # code 0 ('=') is valid nt16
  eq_dec <- DBI::dbGetQuery(con,
    "SELECT seq_decode_4bit([1::UTINYINT, 0::UTINYINT]) AS decoded")
  expect_equal(eq_dec$decoded[1], "A=")

  # unknown char '!' maps to N (15) — permissive like htslib
  inv <- DBI::dbGetQuery(con,
    "SELECT seq_encode_4bit('ACG!')::VARCHAR AS codes")
  expect_equal(inv$codes[1], "[1, 2, 4, 15]")

  # out-of-range code 16 yields NULL
  inv2 <- DBI::dbGetQuery(con,
    "SELECT seq_decode_4bit([1::UTINYINT, 16::UTINYINT]) IS NULL AS is_null")
  expect_true(inv2$is_null[1])

  # empty string roundtrips
  empty <- DBI::dbGetQuery(con,
    "SELECT length(seq_decode_4bit(seq_encode_4bit(''))) AS len")
  expect_equal(empty$len[1], 0L)

  # seq_gc_content
  gc <- DBI::dbGetQuery(con,
    "SELECT seq_gc_content('ACGTNN') AS gc, seq_gc_content('NNNN') IS NULL AS all_n")
  expect_true(abs(gc$gc[1] - 0.5) < 1e-6)
  expect_true(gc$all_n[1])

  # seq_revcomp
  rc <- DBI::dbGetQuery(con, "SELECT seq_revcomp('ACGT') AS rc")
  expect_equal(rc$rc[1], "ACGT")  # ACGT is its own reverse complement

  rc2 <- DBI::dbGetQuery(con, "SELECT seq_revcomp('AAACCC') AS rc")
  expect_equal(rc2$rc[1], "GGGTTT")

  # seq_kmers (k is positional, not named)
  km <- DBI::dbGetQuery(con,
    "SELECT kmer FROM seq_kmers('ACGTAC', 3) ORDER BY pos")
  expect_equal(km$kmer, c("ACG", "CGT", "GTA", "TAC"))

  # =========================================================================
  # read_bam with sequence_encoding := 'nt16'
  # =========================================================================

  # wrapper: default encoding (string) returns VARCHAR SEQ
  rduckhts_bam(con, "bam_str", bam_path, overwrite = TRUE)
  bam_str_type <- DBI::dbGetQuery(con, "SELECT typeof(SEQ) AS t FROM bam_str LIMIT 1")
  expect_equal(bam_str_type$t[1], "VARCHAR")

  # wrapper: configurable htslib decompression workers
  rduckhts_bam(con, "bam_threads_off", bam_path,
    decompression_threads = 0, overwrite = TRUE)
  bam_threads_off_n <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM bam_threads_off")
  expect_equal(bam_threads_off_n$n[1], 112)

  # wrapper: nt16 encoding returns UTINYINT[]
  rduckhts_bam(con, "bam_nt16", bam_path,
    sequence_encoding = "nt16", overwrite = TRUE)
  bam_nt16_type <- DBI::dbGetQuery(con, "SELECT typeof(SEQ) AS t FROM bam_nt16 LIMIT 1")
  expect_equal(bam_nt16_type$t[1], "UTINYINT[]")

  # first 5 codes match expected AGCTA -> [1, 4, 2, 8, 1]
  bam_codes <- DBI::dbGetQuery(con,
    "SELECT SEQ[1:5]::VARCHAR AS codes FROM bam_nt16 LIMIT 1")
  expect_equal(bam_codes$codes[1], "[1, 4, 2, 8, 1]")

  # roundtrip: nt16 decodes back to original string
  bam_rt <- DBI::dbGetQuery(con, paste(
    "SELECT seq_decode_4bit(n.SEQ) = s.SEQ AS match",
    "FROM (SELECT SEQ FROM bam_nt16 LIMIT 1) n,",
    "     (SELECT SEQ FROM bam_str LIMIT 1) s"))
  expect_true(bam_rt$match[1])

  # row counts match between encodings
  bam_str_n <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM bam_str")
  bam_nt16_n <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM bam_nt16")
  expect_equal(bam_str_n$n[1], bam_nt16_n$n[1])

  # direct SQL COUNT(*) should stay on the zero-column BAM fast path
  bam_count_only <- DBI::dbGetQuery(con, sprintf(paste(
    "SELECT COUNT(*) AS n FROM read_bam(",
    "'%s', index_path := '%s', sequence_encoding := 'nt16', quality_representation := 'phred', decompression_threads := 1)"
  ), bam_path, paste0(bam_path, ".bai")))
  expect_equal(bam_count_only$n[1], bam_str_n$n[1])

  # direct SQL: invalid encoding errors
  expect_error(DBI::dbGetQuery(con, sprintf(
    "SELECT * FROM read_bam('%s', sequence_encoding := 'invalid') LIMIT 1",
    bam_path)))

  # direct SQL: invalid decompression thread count errors
  expect_error(DBI::dbGetQuery(con, sprintf(
    "SELECT * FROM read_bam('%s', decompression_threads := -1) LIMIT 1",
    bam_path)))

  # =========================================================================
  # read_fasta with sequence_encoding := 'nt16'
  # =========================================================================

  # wrapper: default encoding returns VARCHAR
  rduckhts_fasta(con, "fa_str", fasta_path,
    region = "CHROMOSOME_I:1-100",
    index_path = fasta_index_path,
    overwrite = TRUE)
  fa_str_type <- DBI::dbGetQuery(con,
    "SELECT typeof(SEQUENCE) AS t FROM fa_str LIMIT 1")
  expect_equal(fa_str_type$t[1], "VARCHAR")

  # wrapper: nt16 encoding returns UTINYINT[]
  rduckhts_fasta(con, "fa_nt16", fasta_path,
    region = "CHROMOSOME_I:1-100",
    index_path = fasta_index_path,
    sequence_encoding = "nt16",
    overwrite = TRUE)
  fa_nt16_type <- DBI::dbGetQuery(con,
    "SELECT typeof(SEQUENCE) AS t FROM fa_nt16 LIMIT 1")
  expect_equal(fa_nt16_type$t[1], "UTINYINT[]")

  # first 5 codes: GCCTA -> [4, 2, 2, 8, 1]
  fa_codes <- DBI::dbGetQuery(con,
    "SELECT SEQUENCE[1:5]::VARCHAR AS codes FROM fa_nt16 LIMIT 1")
  expect_equal(fa_codes$codes[1], "[4, 2, 2, 8, 1]")

  # roundtrip: nt16 decodes back to original string
  fa_rt <- DBI::dbGetQuery(con, paste(
    "SELECT seq_decode_4bit(n.SEQUENCE) = s.SEQUENCE AS match",
    "FROM fa_nt16 n, fa_str s"))
  expect_true(fa_rt$match[1])

  # direct SQL COUNT(*) should stay on the zero-column FASTA fast path
  fa_count_only <- DBI::dbGetQuery(con, sprintf(paste(
    "SELECT COUNT(*) AS n FROM read_fasta(",
    "'%s', index_path := '%s', sequence_encoding := 'nt16')"
  ), fasta_path, fasta_index_path))
  expect_equal(fa_count_only$n[1], 7)

  # direct SQL: invalid encoding errors
  expect_error(DBI::dbGetQuery(con, sprintf(
    "SELECT * FROM read_fasta('%s', sequence_encoding := 'bad')", fasta_path)))

  # =========================================================================
  # read_fastq with sequence_encoding := 'nt16'
  # =========================================================================

  # wrapper: default encoding returns VARCHAR
  rduckhts_fastq(con, "fq_str", fastq_r1, overwrite = TRUE)
  fq_str_type <- DBI::dbGetQuery(con,
    "SELECT typeof(SEQUENCE) AS t FROM fq_str LIMIT 1")
  expect_equal(fq_str_type$t[1], "VARCHAR")

  # wrapper: nt16 encoding returns UTINYINT[]
  rduckhts_fastq(con, "fq_nt16", fastq_r1,
    sequence_encoding = "nt16", overwrite = TRUE)
  fq_nt16_type <- DBI::dbGetQuery(con,
    "SELECT typeof(SEQUENCE) AS t FROM fq_nt16 LIMIT 1")
  expect_equal(fq_nt16_type$t[1], "UTINYINT[]")

  # first 5 codes: CCGTT -> [2, 2, 4, 8, 8]
  fq_codes <- DBI::dbGetQuery(con,
    "SELECT SEQUENCE[1:5]::VARCHAR AS codes FROM fq_nt16 LIMIT 1")
  expect_equal(fq_codes$codes[1], "[2, 2, 4, 8, 8]")

  # roundtrip: nt16 decodes back to original string
  fq_rt <- DBI::dbGetQuery(con, paste(
    "SELECT seq_decode_4bit(n.SEQUENCE) = s.SEQUENCE AS match",
    "FROM (SELECT SEQUENCE FROM fq_nt16 LIMIT 1) n,",
    "     (SELECT SEQUENCE FROM fq_str LIMIT 1) s"))
  expect_true(fq_rt$match[1])

  # row counts match between encodings
  fq_str_n <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM fq_str")
  fq_nt16_n <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM fq_nt16")
  expect_equal(fq_str_n$n[1], fq_nt16_n$n[1])

  # direct SQL COUNT(*) should stay on the zero-column FASTQ fast path
  fq_count_only <- DBI::dbGetQuery(con, sprintf(paste(
    "SELECT COUNT(*) AS n FROM read_fastq(",
    "'%s', sequence_encoding := 'nt16', quality_representation := 'phred')"
  ), fastq_r1))
  expect_equal(fq_count_only$n[1], fq_str_n$n[1])

  # direct SQL: invalid encoding errors
  expect_error(DBI::dbGetQuery(con, sprintf(
    "SELECT * FROM read_fastq('%s', sequence_encoding := 'xyz')", fastq_r1)))

  dbDisconnect(con, shutdown = TRUE)
  message("Sequence operation tests passed!")
}

test_seq_ops()

message("All sequence operation tests completed!")
