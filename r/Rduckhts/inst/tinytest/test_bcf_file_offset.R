# Tests for the opt-in FILE_OFFSET column of read_bcf() through rduckhts_bcf().
#
# FILE_OFFSET is the BGZF virtual position just after each record. The expected
# positions in bcf_record_ends.tsv are read from the BGZF blocks and record
# lengths of each file, not from HTSlib, and are listed in file order
# (test/scripts/prepare_bcf_offset_fixtures.R in the DuckHTS repository).
library(tinytest)
library(DBI)

con <- rduckhts_connect()
dbExecute(con, "SET threads=4")

extdata <- function(name) {
  path <- system.file("extdata", name, package = "Rduckhts")
  expect_true(nzchar(path), info = name)
  path
}
oracle <- read.delim(extdata("bcf_record_ends.tsv"), stringsAsFactors = FALSE)

# The default schema is unchanged; on request FILE_OFFSET follows FILTER.
rduckhts_bcf(con, "bcf_default", extdata("bcf_offset.bcf"), overwrite = TRUE)
expect_false("FILE_OFFSET" %in% dbListFields(con, "bcf_default"))
rduckhts_bcf(con, "bcf_offsets", extdata("bcf_offset.bcf"), include_file_offset = TRUE,
             overwrite = TRUE)
fields <- dbListFields(con, "bcf_offsets")
expect_identical(which(fields == "FILE_OFFSET"), which(fields == "FILTER") + 1L)

# Every scan mode gives each record its expected position, and ordering by it
# gives the records in file order, including the two at chr2:20000.
for (name in c("bcf_offset.bcf", "bcf_offset.vcf.gz")) {
  expected <- oracle[oracle$file == name, ]
  for (mode in c("auto", "sequential")) {
    rduckhts_bcf(con, "bcf_offsets", extdata(name), include_file_offset = TRUE,
                 scan_mode = mode, overwrite = TRUE)
    rows <- dbGetQuery(con, "SELECT ID, FILE_OFFSET FROM bcf_offsets ORDER BY FILE_OFFSET")
    expect_identical(rows$ID, expected$ID, info = paste(name, mode))
    expect_equal(as.numeric(rows$FILE_OFFSET), as.numeric(expected$FILE_OFFSET),
                 info = paste(name, mode))
  }
}

# Tidy rows of one record share its position.
rduckhts_bcf(con, "bcf_tidy", extdata("bcf_offset.vcf.gz"), include_file_offset = TRUE,
             tidy_format = TRUE, overwrite = TRUE)
tidy <- dbGetQuery(con, paste(
  "SELECT ID, count(*) AS n, count(DISTINCT FILE_OFFSET) AS positions",
  "FROM bcf_tidy GROUP BY ID"))
expect_identical(nrow(tidy), 120L)
expect_true(all(tidy$n == 2 & tidy$positions == 1))

# An uncompressed VCF has no BGZF position.
rduckhts_bcf(con, "bcf_plain", extdata("bcf_offset_plain.vcf"), include_file_offset = TRUE,
             overwrite = TRUE)
plain <- dbGetQuery(con, "SELECT FILE_OFFSET FROM bcf_plain")
expect_identical(nrow(plain), 8L)
expect_true(all(is.na(plain$FILE_OFFSET)))

for (bad in list(NA, "yes", c(TRUE, FALSE), 1)) {
  expect_error(
    rduckhts_bcf(con, "bcf_bad", extdata("bcf_offset.bcf"), include_file_offset = bad,
                 overwrite = TRUE),
    "include_file_offset must be TRUE or FALSE"
  )
}
expect_false(dbExistsTable(con, "bcf_bad"))

dbDisconnect(con, shutdown = TRUE)
