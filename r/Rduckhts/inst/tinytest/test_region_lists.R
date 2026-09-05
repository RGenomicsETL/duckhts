library(tinytest)
library(DBI)

test_region_lists <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  file <- function(name) dbQuoteString(con, system.file("extdata", name, package = "Rduckhts"))
  explicit_index <- sprintf("read_bcf(%s, index_path := %s, region := %%s)",
                           file("region_union_no_contig.vcf.gz"), file("region_union_no_contig.index.tbi"))
  expect_error(dbGetQuery(con, paste("SELECT count(*) FROM",
    sprintf(explicit_index, "'chr1:10-20,chr1:20-10'"))), pattern = "region list: invalid item")
  expect_equal(dbGetQuery(con, paste("SELECT count(*) AS n FROM",
    sprintf(explicit_index, "'chr1:11-11,chr1:19-19'")))$n[[1]], 2)
  readers <- c(read_bcf = "vcf_file.bcf", read_bam = "range.bam",
               read_tabix = "gff_file.gff.gz", read_fasta = "ce.fa")
  contigs <- c("1", "CHROMOSOME_I", "X", "CHROMOSOME_I")
  for (i in seq_along(readers)) {
    reader <- names(readers)[[i]]
    path <- file(readers[[i]])
    regions <- c(",", " , , ", paste0(contigs[[i]], ","),
                 paste0(",", contigs[[i]]), paste0(contigs[[i]], ", ,", contigs[[i]]))
    for (region in regions) {
      expect_error(dbGetQuery(con, sprintf("SELECT count(*) FROM %s(%s, region := %s)",
        reader, path, dbQuoteString(con, region))), pattern = "region list: empty item")
      expect_equal(dbGetQuery(con, "SELECT 4242 AS n")$n, 4242L)
    }
    if (reader != "read_fasta") {
      region <- paste0(contigs[[i]], ":10-20,", contigs[[i]], ":20-10")
      expect_error(dbGetQuery(con, sprintf("SELECT count(*) FROM %s(%s, region := %s)",
        reader, path, dbQuoteString(con, region))), pattern = "region list: invalid item")
    }
  }
  for (format in c("bcf", "vcf.gz")) for (tidy in c(FALSE, TRUE)) {
    path <- file(paste0("region_union.", format))
    mode <- if (tidy) "true" else "false"
    suffix <- if (tidy) "SAMPLE_ID,FORMAT_GT" else "FORMAT_GT_S1,FORMAT_GT_S2"
    columns <- paste("ID,POS,REF,ALT,INFO_END,", suffix)
    for (region in c("chr1:19-19,chr1:11-11", "chr1:11-11,chr1:11-11")) {
      actual <- dbGetQuery(con, sprintf(
        "SELECT %s FROM read_bcf(%s, region := %s, tidy_format := %s) ORDER BY ALL",
        columns, path, dbQuoteString(con, region), mode))
      expected <- dbGetQuery(con, sprintf(paste(
        "SELECT %s FROM read_bcf(%s, scan_mode := 'sequential', tidy_format := %s)",
        "WHERE ID IN ('long_ref','long_end') ORDER BY ALL"), columns, path, mode))
      expect_equal(nrow(actual), if (tidy) 4L else 2L)
      expect_equal(actual, expected)
    }
    duplicate <- dbGetQuery(con, sprintf(paste(
      "SELECT count(*) AS n FROM read_bcf(%s,",
      "region := 'chr1:18-18,chr1:18-18', tidy_format := %s) WHERE ID='duplicate'"), path, mode))
    expect_equal(duplicate$n[[1]], if (tidy) 4 else 2)
  }
  fasta <- dbGetQuery(con, sprintf(paste(
    "SELECT count(*) AS n, min(length(SEQUENCE)) AS bases FROM read_fasta(%s,",
    "region := 'CHROMOSOME_I:1-10,CHROMOSOME_I:1-10')"), file("ce.fa")))
  expect_equal(fasta$n[[1]], 2)
  expect_equal(fasta$bases[[1]], 10)
}

test_region_lists()
