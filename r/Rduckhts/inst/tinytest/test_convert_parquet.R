library(tinytest)
library(DBI)

qstr <- function(x) sprintf("'%s'", gsub("'", "''", x, fixed = TRUE))
parquet_metadata <- function(con, path) {
  DBI::dbGetQuery(
    con,
    sprintf(
      "SELECT key::VARCHAR AS key, value::VARCHAR AS value FROM parquet_kv_metadata(%s)",
      qstr(path)
    )
  )
}
meta_value <- function(meta, key) {
  vals <- meta$value[meta$key == key]
  if (length(vals) == 0L) return(NA_character_)
  vals[[1]]
}

with_duckhts_con <- function(code) {
  code <- substitute(code)
  con <- rduckhts_connect()
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  env <- new.env(parent = parent.frame())
  env$con <- con
  eval(code, envir = env)
}

expect_error(
  rduckhts_bcf_convert_parquet(NULL, "input.vcf", tempfile(fileext = ".parquet"), tidy_format = NA),
  "tidy_format must be TRUE or FALSE"
)
expect_error(
  rduckhts_bam_convert_parquet(NULL, "input.bam", tempfile(fileext = ".parquet"), standard_tags = c(TRUE, FALSE)),
  "standard_tags must be TRUE or FALSE"
)
expect_error(
  rduckhts_gff_convert_parquet(NULL, "input.gff3", tempfile(fileext = ".parquet"), header = "yes"),
  "header must be TRUE or FALSE"
)
expect_error(
  rduckhts_tabix_convert_parquet(NULL, "input.tsv.gz", tempfile(fileext = ".parquet"), auto_detect = NA),
  "auto_detect must be TRUE or FALSE"
)

with_duckhts_con({
  bcf_path <- system.file("extdata", "formatcols.vcf.gz", package = "Rduckhts")
  expect_true(file.exists(bcf_path))
  out <- tempfile(fileext = ".parquet")
  sql_preview <- DBI::dbGetQuery(con, sprintf(
    "SELECT duckhts_bcf_convert_parquet_sql(%s, %s, metadata := map(['project'], ['cohort-a'])) AS sql",
    qstr(bcf_path), qstr(out)
  ))$sql
  expect_true(grepl("KV_METADATA", sql_preview, fixed = TRUE))
  expect_equal(
    rduckhts_bcf_convert_parquet(
      con,
      bcf_path,
      out,
      columns = c("CHROM", "POS", "REF", "ALT"),
      metadata = list(
        duckhts_test_marker = "bcf",
        project = "cohort-a",
        `sample.id` = "dot-ok",
        nested = '{"batch":1}'
      ),
      overwrite = TRUE
    ),
    out
  )
  expect_true(file.exists(out))
  meta <- parquet_metadata(con, out)
  expect_equal(meta_value(meta, "duckhts_write_format_version"), "1")
  expect_equal(meta_value(meta, "duckhts_reader"), "read_bcf")
  expect_equal(meta_value(meta, "duckhts_tidy_format"), "false")
  expect_equal(meta_value(meta, "duckhts_test_marker"), "bcf")
  expect_equal(meta_value(meta, "project"), "cohort-a")
  expect_equal(meta_value(meta, "sample.id"), "dot-ok")
  expect_equal(meta_value(meta, "nested"), '{\\x22batch\\x22:1}')
  expect_true(grepl("##fileformat=VCF", meta_value(meta, "vcf_header"), fixed = TRUE))
  expect_true(grepl("#CHROM", meta_value(meta, "vcf_header"), fixed = TRUE))
  expect_error(
    rduckhts_bcf_convert_parquet(con, bcf_path, out, columns = c("CHROM", "POS")),
    "Output already exists"
  )
  rows <- DBI::dbGetQuery(con, sprintf("SELECT count(*) AS n FROM read_parquet(%s)", qstr(out)))$n
  expected <- DBI::dbGetQuery(
    con,
    sprintf("SELECT count(*) AS n FROM read_bcf(%s)", qstr(bcf_path))
  )$n
  expect_equal(rows, expected)
  json_state <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT loaded FROM duckdb_extensions()",
      "WHERE extension_name = 'json'"
    )
  )
  if (nrow(json_state) == 1L) {
    expect_false(json_state$loaded[[1L]])
  }
})

with_duckhts_con({
  bam_path <- system.file("extdata", "range.bam", package = "Rduckhts")
  expect_true(file.exists(bam_path))
  out <- tempfile(fileext = ".parquet")
  expect_silent(rduckhts_bam_convert_parquet(
    con,
    bam_path,
    out,
    columns = c("QNAME", "FLAG", "RNAME", "POS"),
    where = "FLAG >= 0",
    overwrite = TRUE
  ))
  meta <- parquet_metadata(con, out)
  expect_equal(meta_value(meta, "duckhts_reader"), "read_bam")
  expect_true(grepl("@HD", meta_value(meta, "sam_header"), fixed = TRUE))
  expect_equal(meta_value(meta, "duckhts_filter"), "FLAG >= 0")
  expect_true(DBI::dbGetQuery(con, sprintf("SELECT count(*) > 0 AS ok FROM read_parquet(%s)", qstr(out)))$ok)
})

with_duckhts_con({
  gff_path <- system.file("extdata", "gff_file.gff.gz", package = "Rduckhts")
  expect_true(file.exists(gff_path))
  out_dir <- tempfile("duckhts_gff_parquet_")
  expect_silent(rduckhts_gff_convert_parquet(
    con,
    gff_path,
    out_dir,
    columns = c("seqname", "source", "feature", "start", "end"),
    partition_by = "feature",
    overwrite = TRUE
  ))
  files <- list.files(out_dir, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)
  expect_true(length(files) > 0L)
  meta <- parquet_metadata(con, files[[1]])
  expect_equal(meta_value(meta, "duckhts_reader"), "read_gff")
  expect_true(grepl("##gff-version 3", meta_value(meta, "gff_header"), fixed = TRUE))
  expect_equal(meta_value(meta, "duckhts_partition_by"), "feature")
  expect_true(DBI::dbGetQuery(
    con,
    sprintf("SELECT count(*) > 0 AS ok FROM read_parquet(%s, hive_partitioning = true)", qstr(file.path(out_dir, "**", "*.parquet")))
  )$ok)
})

with_duckhts_con({
  tabix_path <- system.file("extdata", "meta_tabix.tsv.gz", package = "Rduckhts")
  expect_true(file.exists(tabix_path))
  out <- tempfile(fileext = ".parquet")
  expect_silent(rduckhts_tabix_convert_parquet(
    con,
    tabix_path,
    out,
    header_names = c("chrom", "pos", "value"),
    auto_detect = TRUE,
    overwrite = TRUE
  ))
  meta <- parquet_metadata(con, out)
  expect_equal(meta_value(meta, "duckhts_reader"), "read_tabix")
  expect_true(grepl("#", meta_value(meta, "tabix_header"), fixed = TRUE))
  expect_equal(meta_value(meta, "duckhts_header_names"), "chrom,pos,value")
  expect_identical(
    DBI::dbGetQuery(con, sprintf("SELECT typeof(pos) AS type FROM read_parquet(%s) LIMIT 1", qstr(out)))$type,
    "BIGINT"
  )
})
