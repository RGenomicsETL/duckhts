library(tinytest)
library(DBI)

test_pileup <- function() {
  con <- rduckhts_connect()
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
  bam_index_path <- system.file("extdata", "range.bam.bai", package = "Rduckhts")

  rduckhts_pileup(con, "pileup_rows", bam_path,
    region = "CHROMOSOME_I:900-1100",
    index_path = bam_index_path,
    overwrite = TRUE)

  types <- DBI::dbGetQuery(con, paste(
    "SELECT typeof(chrom) AS chrom_t, typeof(pos) AS pos_t,",
    "typeof(depth) AS depth_t, typeof(bases) AS bases_t, typeof(quals) AS quals_t",
    "FROM pileup_rows LIMIT 1"
  ))
  expect_equal(types$chrom_t[1], "VARCHAR")
  expect_equal(types$pos_t[1], "BIGINT")
  expect_equal(types$depth_t[1], "INTEGER")
  expect_equal(types$bases_t[1], "VARCHAR")
  expect_equal(types$quals_t[1], "VARCHAR")

  stats <- DBI::dbGetQuery(con, paste(
    "SELECT COUNT(*) AS n, MIN(pos) AS min_pos, MAX(pos) AS max_pos,",
    "COUNT(*) FILTER (WHERE depth = length(bases) AND depth = length(quals)) AS ok_rows,",
    "MAX(depth) AS max_depth",
    "FROM pileup_rows"
  ))
  expect_equal(stats$n[1], 133)
  expect_equal(stats$min_pos[1], 914)
  expect_equal(stats$max_pos[1], 1047)
  expect_equal(stats$ok_rows[1], stats$n[1])
  expect_equal(stats$max_depth[1], 2)

  head_rows <- DBI::dbGetQuery(con,
    "SELECT chrom, pos, depth, bases, quals FROM pileup_rows ORDER BY pos LIMIT 3")
  expect_identical(head_rows$chrom, rep("CHROMOSOME_I", 3))
  expect_identical(head_rows$pos, c(914, 915, 916))
  expect_identical(head_rows$depth, c(1L, 1L, 1L))
  expect_identical(head_rows$bases, c("a", "g", "c"))
  expect_identical(head_rows$quals, c("G", "G", "G"))

  strict_mapq <- DBI::dbGetQuery(con, sprintf(paste(
    "SELECT COUNT(*) AS n FROM read_pileup(",
    "'%s', region := 'CHROMOSOME_I:900-1100', index_path := '%s', min_mapq := 255)"
  ), bam_path, bam_index_path))
  expect_equal(strict_mapq$n[1], 0)

  clipped <- DBI::dbGetQuery(con, sprintf(paste(
    "SELECT COUNT(*) AS n, MIN(pos) AS min_pos, MAX(pos) AS max_pos",
    "FROM read_pileup(",
    "'%s', region := 'CHROMOSOME_I:930-980', index_path := '%s', min_mapq := 20)"
  ), bam_path, bam_index_path))
  expect_equal(clipped$n[1], 51)
  expect_equal(clipped$min_pos[1], 930)
  expect_equal(clipped$max_pos[1], 980)

  expect_error(DBI::dbGetQuery(con, sprintf(
    "SELECT * FROM read_pileup('%s', region := 'bad region', index_path := '%s')",
    bam_path, bam_index_path
  )))

}

test_pileup()
