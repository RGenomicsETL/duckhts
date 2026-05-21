library(tinytest)
library(DBI)

test_indexes_bgzip <- function() {
  drv <- duckdb::duckdb(config = list(allow_unsigned_extensions = "true"))
  con <- dbConnect(drv)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  expect_silent(rduckhts_load(con))

  bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
  vcf_path <- system.file("extdata", "formatcols.vcf.gz", package = "Rduckhts")
  fasta_path <- system.file("extdata", "ce.fa", package = "Rduckhts")

  tmp_dir <- tempfile("duckhts_index_build_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  bam_index_path <- file.path(tmp_dir, "range.bam.bai")
  vcf_index_path <- file.path(tmp_dir, "formatcols.vcf.gz.tbi")
  bed_path <- file.path(tmp_dir, "targets.bed")
   auto_fasta_path <- file.path(tmp_dir, "ce.fa")
  auto_fasta_index_path <- paste0(auto_fasta_path, ".fai")
  auto_bam_path <- file.path(tmp_dir, "range.auto.bam")
  auto_bam_index_path <- paste0(auto_bam_path, ".bai")
  auto_vcf_path <- file.path(tmp_dir, "formatcols.auto.vcf.gz")
  auto_vcf_index_path <- paste0(auto_vcf_path, ".tbi")
  auto_bed_path <- file.path(tmp_dir, "targets.auto.bed")
  bed_gz_path <- paste0(bed_path, ".gz")
  bed_tbi_path <- paste0(bed_gz_path, ".tbi")
  bed_roundtrip_path <- file.path(tmp_dir, "targets.roundtrip.bed")

  writeLines(
    c(
      "chr1\t0\t10\ta",
      "chr1\t10\t20\tb",
      "chr2\t0\t5\tc"
    ),
    bed_path
  )
  writeLines(
    c(
      "chr1\t0\t10\ta",
      "chr1\t10\t20\tb",
      "chr2\t0\t5\tc"
    ),
    auto_bed_path
  )
  expect_true(file.copy(fasta_path, auto_fasta_path, overwrite = TRUE))
  expect_true(file.copy(bam_path, auto_bam_path, overwrite = TRUE))
  expect_true(file.copy(vcf_path, auto_vcf_path, overwrite = TRUE))

  bam_idx <- rduckhts_bam_index(con, bam_path, index_path = bam_index_path, threads = 1)
  expect_true(isTRUE(bam_idx$success[1]))
  expect_equal(bam_idx$index_format[1], "BAI")
  expect_true(file.exists(bam_index_path))

  vcf_idx <- rduckhts_bcf_index(con, vcf_path, index_path = vcf_index_path, threads = 1)
  expect_true(isTRUE(vcf_idx$success[1]))
  expect_equal(vcf_idx$index_format[1], "TBI")
  expect_true(file.exists(vcf_index_path))

  bgz <- rduckhts_bgzip(
    con,
    bed_path,
    output_path = bed_gz_path,
    threads = 1,
    keep = TRUE,
    overwrite = TRUE
  )
  expect_true(isTRUE(bgz$success[1]))
  expect_true(file.exists(bed_gz_path))

  tbx <- rduckhts_tabix_index(
    con,
    bed_gz_path,
    preset = "bed",
    index_path = bed_tbi_path,
    threads = 1
  )
  expect_true(isTRUE(tbx$success[1]))
  expect_equal(tbx$index_format[1], "TBI")
  expect_true(file.exists(bed_tbi_path))

  gunz <- rduckhts_bgunzip(
    con,
    bed_gz_path,
    output_path = bed_roundtrip_path,
    threads = 1,
    keep = TRUE,
    overwrite = TRUE
  )
  expect_true(isTRUE(gunz$success[1]))
  expect_true(file.exists(bed_roundtrip_path))
  expect_equal(readLines(bed_roundtrip_path), readLines(bed_path))

  fai_auto <- rduckhts_fasta_index(con, auto_fasta_path)
  expect_true(isTRUE(fai_auto$success[1]))
  expect_equal(fai_auto$index_path[1], auto_fasta_index_path)
  expect_true(file.exists(auto_fasta_index_path))

  bam_auto <- rduckhts_bam_index(con, auto_bam_path, threads = 1)
  expect_true(isTRUE(bam_auto$success[1]))
  expect_equal(bam_auto$index_path[1], auto_bam_index_path)
  expect_true(file.exists(auto_bam_index_path))

  vcf_auto <- rduckhts_bcf_index(con, auto_vcf_path, threads = 1)
  expect_true(isTRUE(vcf_auto$success[1]))
  expect_equal(vcf_auto$index_path[1], auto_vcf_index_path)
  expect_true(file.exists(auto_vcf_index_path))

  bgz_auto <- rduckhts_bgzip(
    con,
    auto_bed_path,
    threads = 1,
    keep = FALSE,
    overwrite = TRUE
  )
  expect_true(isTRUE(bgz_auto$success[1]))
  expect_equal(bgz_auto$output_path[1], paste0(auto_bed_path, ".gz"))
  expect_true(file.exists(bgz_auto$output_path[1]))
  expect_false(file.exists(auto_bed_path))

  tbx_auto <- rduckhts_tabix_index(
    con,
    bgz_auto$output_path[1],
    preset = "bed",
    threads = 1
  )
  expect_true(isTRUE(tbx_auto$success[1]))
  expect_equal(tbx_auto$index_path[1], paste0(bgz_auto$output_path[1], ".tbi"))
  expect_true(file.exists(tbx_auto$index_path[1]))

  gunz_auto <- rduckhts_bgunzip(
    con,
    bgz_auto$output_path[1],
    threads = 1,
    keep = FALSE,
    overwrite = TRUE
  )
  expect_true(isTRUE(gunz_auto$success[1]))
  expect_equal(gunz_auto$output_path[1], auto_bed_path)
  expect_true(file.exists(auto_bed_path))
  expect_equal(readLines(auto_bed_path), readLines(bed_path))

  expect_silent(rduckhts_bam(
    con,
    "reads_idx_tmp",
    bam_path,
    region = "CHROMOSOME_I:1-1000",
    index_path = bam_index_path,
    overwrite = TRUE
  ))
  expect_equal(DBI::dbGetQuery(con, "SELECT count(*) AS n FROM reads_idx_tmp")$n[1], 2)

  expect_silent(rduckhts_bcf(
    con,
    "variants_idx_tmp",
    vcf_path,
    region = "1:100-100",
    index_path = vcf_index_path,
    overwrite = TRUE
  ))
  expect_equal(DBI::dbGetQuery(con, "SELECT count(*) AS n FROM variants_idx_tmp")$n[1], 1)

  expect_silent(rduckhts_tabix(
    con,
    "tabix_idx_tmp",
    bed_gz_path,
    header_names = c("chrom", "start", "end", "name"),
    column_types = c("VARCHAR", "BIGINT", "BIGINT", "VARCHAR"),
    region = "chr1",
    index_path = bed_tbi_path,
    overwrite = TRUE
  ))
  expect_equal(DBI::dbGetQuery(con, "SELECT count(*) AS n FROM tabix_idx_tmp")$n[1], 2)
}

test_indexes_bgzip()
