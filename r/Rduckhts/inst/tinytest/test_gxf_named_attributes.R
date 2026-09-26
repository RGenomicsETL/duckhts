library(tinytest)
library(DBI)

(function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE))
  cases <- list(
    gff = list(wrapper = rduckhts_gff, file = "gff_named_attributes.gff3.gz",
               keys = c("ID", "Parent")),
    gtf = list(wrapper = rduckhts_gtf, file = "gtf_named_attributes.gtf.gz",
               keys = c("gene_id", "transcript_id"))
  )
  for (kind in names(cases)) {
    case <- cases[[kind]]
    path <- system.file("extdata", case$file, package = "Rduckhts", mustWork = TRUE)
    keys <- c(case$keys, "encoded", "spaced", "empty", "missing", "odd-key", "quote'key")
    expect_silent(case$wrapper(con, "named", path, attributes = keys,
                              attributes_map = TRUE, overwrite = TRUE))
    actual <- dbGetQuery(con, "SELECT * FROM named ORDER BY seqname, start")
    expect_identical(tail(names(actual), length(keys)), keys)
    for (key in keys) {
      comparison <- dbGetQuery(con, sprintf(
        "SELECT bool_and(%s IS NOT DISTINCT FROM attributes_map[%s]) AS same FROM named",
        dbQuoteIdentifier(con, key), dbQuoteString(con, key)
      ))
      expect_true(comparison$same, info = paste(kind, key))
    }
    expect_identical(actual[[case$keys[1]]], c("first", "third", NA, "fourth", "last"))
    expect_identical(actual$encoded, c("a%3Bb%20c", "%25", NA, "a=b", NA))
    expect_identical(actual$empty, c("", "", NA, NA, NA))
    expect_true(all(is.na(actual$missing)))

    case$wrapper(con, "named_seq", path, attributes = keys, scan_mode = "sequential",
                 overwrite = TRUE)
    case$wrapper(con, "named_region", path, attributes = keys,
                 region = "chr1:5-15,chr1:14-20", overwrite = TRUE)
    expected <- dbGetQuery(con, paste(
      "SELECT * FROM named_seq WHERE seqname = 'chr1' AND start <= 20 AND \"end\" >= 5",
      "ORDER BY seqname, start"
    ))
    indexed <- dbGetQuery(con, "SELECT * FROM named_region ORDER BY seqname, start")
    expect_equal(nrow(indexed), 2L)
    expect_identical(indexed, expected)

    expect_silent(case$wrapper(con, "named_empty", path, attributes = character(),
                              overwrite = TRUE))
    expect_equal(length(dbListFields(con, "named_empty")), 9L)
    expect_error(case$wrapper(con, "bad", path, attributes = 1), "character vector")
    expect_error(case$wrapper(con, "bad", path, attributes = NA_character_), "NA")
    expect_error(case$wrapper(con, "bad", path, attributes = ""), "empty")
    expect_error(case$wrapper(con, "bad", path, attributes = "START"), "collides")
    expect_error(case$wrapper(con, "bad", path, attributes = c("ID", "id")), "distinct")
    expect_error(case$wrapper(con, "bad", path, attributes = paste0("k", 1:257)), "256")

    for (invalid in list(NA_character_, "START")) {
      dbExecute(con, "CREATE OR REPLACE TABLE preserved AS SELECT 1 AS marker UNION ALL SELECT 2")
      expected <- dbGetQuery(con, "SELECT marker FROM preserved ORDER BY marker")
      expect_error(
        case$wrapper(con, "preserved", path, attributes = invalid, overwrite = TRUE),
        if (is.na(invalid)) "NA" else "collides",
        info = kind
      )
      expect_true(dbExistsTable(con, "preserved"), info = paste(kind, invalid))
      if (dbExistsTable(con, "preserved")) {
        expect_identical(
          dbGetQuery(con, "SELECT marker FROM preserved ORDER BY marker"),
          expected,
          info = paste(kind, invalid)
        )
      }
    }
  }
})()
