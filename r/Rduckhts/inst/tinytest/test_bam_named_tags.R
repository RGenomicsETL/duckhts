library(tinytest)
library(DBI)

# auxiliary_tags accepts tag names and returns one VARCHAR column per tag,
# formatted like the AUXILIARY_TAGS map value and NULL when a record lacks it.
test_bam_named_tags <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  path <- system.file("extdata", "bam_materialize.sam", package = "Rduckhts")
  expect_true(nzchar(path))

  rduckhts_bam(con, "named", path, auxiliary_tags = c("XZ", "X7"), overwrite = TRUE)
  named <- dbGetQuery(con, "SELECT QNAME, XZ, X7 FROM named ORDER BY QNAME")
  expect_equal(named$XZ, c("previous", "", NA_character_))
  expect_equal(named$X7, c("f,0,1.25,1e+30", "f,-1.5,2.25", NA_character_))
  expect_false("AUXILIARY_TAGS" %in% dbListFields(con, "named"))

  rduckhts_bam(con, "mapped", path, auxiliary_tags = TRUE, overwrite = TRUE)
  mapped <- dbGetQuery(con, paste(
    "SELECT QNAME, map_extract_value(AUXILIARY_TAGS, 'XZ') AS XZ,",
    "map_extract_value(AUXILIARY_TAGS, 'X7') AS X7 FROM mapped ORDER BY QNAME"))
  expect_equal(named, mapped)

  # One name is still a list, not a string.
  rduckhts_bam(con, "one", path, auxiliary_tags = "XA", overwrite = TRUE)
  expect_equal(dbGetQuery(con, "SELECT XA FROM one ORDER BY QNAME")$XA,
               c("q", NA_character_, NA_character_))

  # The multi-file wrapper renders a single name as a list too.
  rduckhts_bam_multi(con, "multi", path, auxiliary_tags = "XZ", overwrite = TRUE)
  multi <- dbGetQuery(con, "SELECT XZ FROM multi ORDER BY QNAME")
  expect_equal(multi$XZ, c("previous", "", NA_character_))

  for (bad in list(NA, NA_character_, 1L, character(0), rep("XZ", 65L), c(TRUE, FALSE))) {
    expect_error(rduckhts_bam(con, "bad", path, auxiliary_tags = bad, overwrite = TRUE),
                 "auxiliary_tags must be TRUE, FALSE or a character vector")
  }
  expect_error(rduckhts_bam(con, "bad", path, auxiliary_tags = c("xz", "XZ"), overwrite = TRUE),
               "the same column")
  expect_error(rduckhts_bam(con, "bad", path, auxiliary_tags = "X", overwrite = TRUE),
               "two characters")
}

test_bam_named_tags()
