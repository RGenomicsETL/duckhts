library(tinytest)
library(DBI)

test_interval_readers <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  expect_silent(rduckhts_load(con))

  fasta_path <- system.file("extdata", "ce.fa", package = "Rduckhts")
  bed_path <- system.file("extdata", "targets.bed", package = "Rduckhts")
  fasta_index_path <- tempfile("duckhts_ce_", fileext = ".fai")

  expect_true(file.exists(fasta_path))
  expect_true(file.exists(bed_path))
  expect_true(rduckhts_fasta_index(con, fasta_path, index_path = fasta_index_path)$success[1])

  expect_silent(rduckhts_bed(con, "targets", bed_path, overwrite = TRUE))
  expect_equal(DBI::dbGetQuery(con, "SELECT count(*) AS n FROM targets")$n[1], 4)

  typed_row <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT chrom, start, \"end\", name, score, strand, thick_start, block_count",
      "FROM targets ORDER BY chrom, start LIMIT 1"
    )
  )
  expect_equal(typed_row$chrom[1], "CHROMOSOME_I")
  expect_equal(typed_row$start[1], 0)
  expect_equal(typed_row$end[1], 10)
  expect_equal(typed_row$name[1], "target1")
  expect_equal(typed_row$score[1], "100")
  expect_equal(typed_row$strand[1], "+")
  expect_equal(typed_row$thick_start[1], 0)
  expect_equal(typed_row$block_count[1], 2)

  nuc <- rduckhts_fasta_nuc(con, fasta_path, bed_path = bed_path, index_path = fasta_index_path)
  first_nuc <- nuc[nuc$chrom == "CHROMOSOME_I" & nuc$start == 0, , drop = FALSE]
  expect_equal(nrow(first_nuc), 1)
  expect_equal(first_nuc$pct_at[1], 0.4)
  expect_equal(first_nuc$pct_gc[1], 0.6)
  expect_equal(first_nuc$num_a[1], 2)
  expect_equal(first_nuc$num_c[1], 4)
  expect_equal(first_nuc$num_g[1], 2)
  expect_equal(first_nuc$num_t[1], 2)
  expect_equal(first_nuc$seq_len[1], 10)

  bins <- rduckhts_fasta_nuc(
    con,
    fasta_path,
    bin_width = 10,
    region = "CHROMOSOME_I:1-20",
    index_path = fasta_index_path
  )
  expect_equal(nrow(bins), 2)
  expect_equal(sum(bins$seq_len), 20)

  seq_nuc <- rduckhts_fasta_nuc(
    con,
    fasta_path,
    bed_path = bed_path,
    index_path = fasta_index_path,
    include_seq = TRUE
  )
  first_seq <- seq_nuc[seq_nuc$chrom == "CHROMOSOME_I" & seq_nuc$start == 0, "seq", drop = TRUE]
  expect_equal(first_seq[1], "GCCTAAGCCT")

  gz_path <- tempfile("targets", fileext = ".bed.gz")
  tbi_path <- paste0(gz_path, ".tbi")
  expect_true(rduckhts_bgzip(con, bed_path, output_path = gz_path, keep = TRUE, overwrite = TRUE)$success[1])
  expect_true(rduckhts_tabix_index(con, gz_path, preset = "bed", index_path = tbi_path, threads = 1)$success[1])
  expect_silent(rduckhts_bed(con, "targets_idx", gz_path, region = "CHROMOSOME_I:1-20", index_path = tbi_path, overwrite = TRUE))
  expect_equal(DBI::dbGetQuery(con, "SELECT count(*) AS n FROM targets_idx")$n[1], 2)
}

test_interval_readers()
