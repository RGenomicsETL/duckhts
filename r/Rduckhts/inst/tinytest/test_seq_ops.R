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

  # seq_gc_content and SIMD backend diagnostics
  gc <- DBI::dbGetQuery(con,
    "SELECT seq_gc_content('ACGTNN') AS gc, seq_gc_content('NNNN') IS NULL AS all_n")
  expect_true(abs(gc$gc[1] - 0.5) < 1e-6)
  expect_true(gc$all_n[1])
  gc_lower <- DBI::dbGetQuery(con, "SELECT seq_gc_content('acgtnn') AS gc")
  expect_true(abs(gc_lower$gc[1] - 0.5) < 1e-6)
  simd <- DBI::dbGetQuery(con,
    "SELECT duckhts_simd_backend() AS selected, duckhts_simd_requested_backend() AS requested, duckhts_simd_backend_available('scalar') AS has_scalar")
  expect_true(nzchar(simd$selected[1]))
  expect_true(nzchar(simd$requested[1]))
  expect_true(simd$has_scalar[1])
  forced <- DBI::dbGetQuery(con, "SELECT backend AS selected FROM duckhts_simd_set_backend('scalar')")
  expect_identical(forced$selected[1], "scalar")
  auto <- DBI::dbGetQuery(con, "SELECT backend AS selected FROM duckhts_simd_set_backend('auto')")
  expect_true(nzchar(auto$selected[1]))
  auto_req <- DBI::dbGetQuery(con, "SELECT duckhts_simd_requested_backend() AS requested")
  expect_identical(auto_req$requested[1], "auto")
  expect_true(nzchar(rduckhts_simd_backend(con)))
  expect_identical(rduckhts_simd_requested_backend(con), "auto")
  expect_true(rduckhts_simd_backend_available(con, "scalar"))
  expect_identical(rduckhts_simd_set_backend(con, "scalar"), "scalar")
  expect_identical(rduckhts_simd_backend(con), "scalar")
  expect_true(nzchar(rduckhts_simd_set_backend(con, "auto")))
  expect_identical(rduckhts_simd_requested_backend(con), "auto")
  expect_false(rduckhts_simd_backend_available(con, "not-a-backend"))
  expect_true(rduckhts_simd_backend_compiled(con, "scalar"))
  expect_false(rduckhts_simd_backend_compiled(con, "sse2"))
  expect_true(rduckhts_simd_backend_cpu_supported(con, "scalar"))
  expect_false(rduckhts_simd_backend_cpu_supported(con, "not-a-backend"))
  expect_true(is.logical(rduckhts_simd_backend_available(con, "avx512")))
  simd_info <- rduckhts_simd_info(con)
  expect_true(is.data.frame(simd_info))
  expect_true(all(c("backend", "selectable", "compiled", "cpu_supported", "available", "selected") %in% names(simd_info)))
  expect_true("scalar" %in% simd_info$backend)
  expect_true("avx512" %in% simd_info$backend)
  expect_false("auto" %in% simd_info$backend)
  expect_true(all(simd_info$available == (simd_info$compiled & simd_info$cpu_supported)))
  kernel_info <- rduckhts_simd_kernel_info(con)
  expect_true(is.data.frame(kernel_info))
  expect_true(all(c("kernel", "selected_backend", "selected_capability", "requested_backend", "scalar_fallback") %in% names(kernel_info)))
  expect_true("seq_base_counts" %in% kernel_info$kernel)
  expect_true(all(nzchar(kernel_info$selected_backend)))
  expect_error(rduckhts_simd_backend_available(con, character()), "single non-missing")
  expect_error(rduckhts_simd_set_backend(con, c("scalar", "auto")), "single non-missing")
  expect_error(rduckhts_simd_set_backend(con, NA_character_), "single non-missing")
  expect_error(rduckhts_simd_set_backend(con, ""), "backend")
  expect_true(nzchar(rduckhts_simd_set_backend(con, "auto")))
  expect_error(
    DBI::dbGetQuery(con, paste(
      "SELECT * FROM duckhts_simd_set_backend('scalar')",
      "UNION ALL SELECT * FROM duckhts_simd_set_backend('not-a-backend')"
    )),
    "unknown SIMD backend"
  )
  expect_identical(rduckhts_simd_requested_backend(con), "auto")
  param_set <- DBI::dbGetQuery(con, "SELECT backend FROM duckhts_simd_set_backend(?)", params = list("scalar"))
  expect_identical(param_set$backend[1], "scalar")
  expect_true(nzchar(rduckhts_simd_set_backend(con, "auto")))
  expect_error(rduckhts_simd_set_backend(con, "sse2"), "no selectable implementation")
  expect_error(rduckhts_simd_set_backend(con, "not-a-backend"), "unknown SIMD backend")

  # scalar vs auto-selected backend correctness on sequences long enough to
  # exercise the SIMD vector loop (>32 bytes for AVX2, >16 for NEON/wasm).
  # The test is ISA-agnostic: it forces scalar, records results, forces auto,
  # records results again, and asserts they are identical.
  seqs <- list(
    gc50_64   = strrep("GCATGCAT", 8L),   # 64 chars, GC = 0.5
    gc100_64  = strrep("GC",       32L),  # 64 chars, GC = 1.0
    gc0_64    = strrep("AT",       32L),  # 64 chars, GC = 0.0
    gc50_96   = strrep("GCGCATAT", 12L),  # 96 chars, GC = 0.5
    gc_n      = strrep("GCNN",     24L),  # 96 chars, GC/called = 1.0
    gc_lower  = strrep("gcgcatat", 12L)   # soft-masked lowercase, GC = 0.5
  )
  sql_vals <- function(xs) paste(sprintf("seq_gc_content('%s')", xs), collapse = ", ")

  expect_identical(rduckhts_simd_set_backend(con, "scalar"), "scalar")
  scalar_kernel_info <- rduckhts_simd_kernel_info(con)
  expect_identical(scalar_kernel_info$selected_backend[scalar_kernel_info$kernel == "seq_base_counts"], "scalar")
  expect_false(scalar_kernel_info$scalar_fallback[scalar_kernel_info$kernel == "seq_base_counts"])
  scalar_gc <- DBI::dbGetQuery(con, sprintf("SELECT %s", sql_vals(seqs)))

  expect_true(nzchar(rduckhts_simd_set_backend(con, "auto")))
  auto_gc <- DBI::dbGetQuery(con, sprintf("SELECT %s", sql_vals(seqs)))

  for (col in names(scalar_gc)) {
    expect_true(
      abs(scalar_gc[[col]][1] - auto_gc[[col]][1]) < 1e-9,
      info = sprintf("scalar vs auto GC mismatch for sequence '%s'", col)
    )
  }
  # Leave backend in auto state
  expect_true(nzchar(rduckhts_simd_set_backend(con, "auto")))

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

  # wrapper: binary CIGAR representation returns packed BAM ops
  rduckhts_bam(con, "bam_cigar_bin", bam_path,
    cigar_representation = "binary", overwrite = TRUE)
  bam_cigar_bin_type <- DBI::dbGetQuery(con,
    "SELECT typeof(CIGAR) AS t FROM bam_cigar_bin LIMIT 1")
  expect_equal(bam_cigar_bin_type$t[1], "UINTEGER[]")
  bam_cigar_bin_txt <- DBI::dbGetQuery(con,
    "SELECT CIGAR::VARCHAR AS cigar FROM bam_cigar_bin ORDER BY POS LIMIT 1")
  expect_equal(bam_cigar_bin_txt$cigar[1], "[1248, 18, 352]")

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
  expect_error(DBI::dbGetQuery(con, sprintf(
    "SELECT * FROM read_bam('%s', cigar_representation := 'invalid') LIMIT 1",
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

  # bgzipped FASTA can use an explicit relocated .gzi sidecar
  gz_fasta_path <- tempfile("duckhts_ce_bgz_", fileext = ".fa.gz")
  gz_fasta_index_path <- tempfile("duckhts_ce_bgz_", fileext = ".fa.gz.fai")
  custom_gzi_path <- tempfile("duckhts_ce_bgz_", fileext = ".gzi")
  on.exit(unlink(c(gz_fasta_path, gz_fasta_index_path, custom_gzi_path), force = TRUE), add = TRUE)
  expect_true(rduckhts_bgzip(
    con,
    fasta_path,
    output_path = gz_fasta_path,
    keep = TRUE,
    overwrite = TRUE
  )$success[1])
  expect_true(rduckhts_fasta_index(con, gz_fasta_path, index_path = gz_fasta_index_path)$success[1])
  default_gzi_path <- paste0(gz_fasta_path, ".gzi")
  expect_true(file.rename(default_gzi_path, custom_gzi_path))
  expect_error(DBI::dbGetQuery(con, sprintf(paste(
    "SELECT * FROM read_fasta(",
    "'%s', region := 'CHROMOSOME_I:1-10', index_path := '%s')"
  ), gz_fasta_path, gz_fasta_index_path)))
  rduckhts_fasta(con, "fa_bgz_gzi", gz_fasta_path,
    region = "CHROMOSOME_I:1-10",
    index_path = gz_fasta_index_path,
    gzi_path = custom_gzi_path,
    overwrite = TRUE)
  fa_bgz_gzi <- DBI::dbGetQuery(con,
    "SELECT NAME, length(SEQUENCE) AS len FROM fa_bgz_gzi")
  expect_equal(fa_bgz_gzi$NAME[1], "CHROMOSOME_I")
  expect_equal(fa_bgz_gzi$len[1], 10)
  fa_bgz_nuc <- rduckhts_fasta_nuc(
    con,
    gz_fasta_path,
    bin_width = 10,
    region = "CHROMOSOME_I:1-20",
    index_path = gz_fasta_index_path,
    gzi_path = custom_gzi_path
  )
  expect_equal(sum(fa_bgz_nuc$seq_len), 20)

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
