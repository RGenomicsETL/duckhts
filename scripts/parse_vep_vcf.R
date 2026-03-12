#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
})

vep_vcf <- Sys.getenv("VEP_VCF", unset = "")
if (!nzchar(vep_vcf)) {
  stop("Set VEP_VCF to a local VCF/BCF path.", call. = FALSE)
}
if (!file.exists(vep_vcf)) {
  stop(sprintf("VEP_VCF does not exist: %s", vep_vcf), call. = FALSE)
}

ext_path <- Sys.getenv(
  "DUCKHTS_EXTENSION",
  unset = "build/release/duckhts.duckdb_extension"
)
ext_path <- normalizePath(ext_path, mustWork = TRUE)

sql_quote <- function(x) {
  sprintf("'%s'", gsub("'", "''", x, fixed = TRUE))
}
drv <- duckdb::duckdb(
  dbdir = ":memory:",
  config = list(allow_unsigned_extensions = "true")
)
con <- dbConnect(drv)
on.exit(
  {
    try(dbDisconnect(con, shutdown = TRUE), silent = TRUE)
  },
  add = TRUE
)

invisible(dbExecute(con, sprintf("LOAD '%s';", gsub("\\\\", "/", ext_path))))
invisible(dbExecute(con, "PRAGMA threads=4;"))

src <- sprintf(
  "read_bcf(%s)",
  sql_quote(normalizePath(vep_vcf, mustWork = TRUE))
)

cat("VEP_VCF:", vep_vcf, "\n")
cat("Extension:", ext_path, "\n\n")

header_sql <- sprintf(
  paste(
    "SELECT record_type, id, value_type",
    "FROM read_hts_header(%s)",
    "WHERE record_type = 'INFO' AND id IN ('CSQ', 'ANN', 'BCSQ', 'VEP', 'vep')"
  ),
  sql_quote(normalizePath(vep_vcf, mustWork = TRUE))
)
print(dbGetQuery(con, header_sql))
cat("\n")

count_sql <- sprintf("SELECT COUNT(*) AS n FROM %s", src)
print(dbGetQuery(con, count_sql))
cat("\n")

describe_sql <- sprintf("DESCRIBE SELECT * FROM %s", src)
cols <- dbGetQuery(con, describe_sql)
print(head(cols, 25))
cat("\n")

vep_cols <- grep("^VEP_", cols$column_name, value = TRUE)
if (length(vep_cols) > 0) {
  picked <- paste(c("CHROM", "POS", head(vep_cols, 8)), collapse = ", ")
  sample_sql <- sprintf("SELECT %s FROM %s LIMIT 1000000000", picked, src)
  print(head(dbGetQuery(con, sample_sql)))
} else {
  cat("No VEP-style consequence columns were materialized.\n")
  sample_sql <- sprintf(
    "SELECT CHROM, POS, REF, ALT FROM %s LIMIT 1000000000",
    src
  )
  print(head(dbGetQuery(con, sample_sql)))
}
