# Per-child chr20 truth windows: for every 100 kb window, the number of called
# genotypes and of heterozygous genotypes. A window without called genotypes keeps
# sites = 0, so truth_intervals.R can tell missing evidence apart from homozygous
# evidence.
cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
path <- file.path(cache, "chr20.children.bcf")
output <- duckhtsbench::duckhts_bench_artifact_path("roh_ancestry_chr20_truth_windows")
# The extension built from this tree, loaded as the arm driver loads it.
driver <- duckdb::duckdb(dbdir = ":memory:", shared_home = FALSE,
  allow_extensions = TRUE, config = list(memory_limit = "12GB", threads = "4",
    temp_directory = file.path(cache, "duckdb-tmp"),
    allow_unsigned_extensions = "true", autoinstall_known_extensions = "false",
    autoload_known_extensions = "false"))
con <- DBI::dbConnect(driver)
q <- function(x) as.character(DBI::dbQuoteString(con, x))
DBI::dbExecute(con, sprintf("LOAD %s", q(normalizePath("build/release/duckhts.duckdb_extension"))))

children <- readLines(file.path(cache, "chr20.children.present.txt"))
length_bp <- 64444167L
windows <- data.frame(window_id = seq.int(0L, (length_bp - 1L) %/% 100000L))
DBI::dbWriteTable(con, "children", data.frame(sample_id = children), temporary = TRUE)
DBI::dbWriteTable(con, "windows", windows, temporary = TRUE)

# sites counts the sample's fully called diploid genotypes in the window: a missing
# or partly missing genotype is no evidence either way. The input is biallelic.
called <- "('0/0', '0/1', '1/0', '1/1', '0|0', '0|1', '1|0', '1|1')"
DBI::dbExecute(con, sprintf(paste0(
  "CREATE TEMP TABLE observed AS ",
  "SELECT SAMPLE_ID AS sample_id, floor((POS - 1) / 100000)::INTEGER AS window_id, ",
  "count(*) FILTER (WHERE FORMAT_GT IN ", called, ")::INTEGER AS sites, ",
  "count(*) FILTER (WHERE FORMAT_GT IN ('0/1', '1/0', '0|1', '1|0'))::INTEGER AS heterozygotes ",
  "FROM read_bcf(%s, tidy_format := true) GROUP BY SAMPLE_ID, window_id"), q(path)))
DBI::dbExecute(con, sprintf(paste0(
  "COPY (SELECT c.sample_id, w.window_id, ",
  "coalesce(o.sites, 0)::INTEGER AS sites, ",
  "coalesce(o.heterozygotes, 0)::INTEGER AS heterozygotes ",
  "FROM children c CROSS JOIN windows w ",
  "LEFT JOIN observed o USING (sample_id, window_id) ",
  "ORDER BY c.sample_id, w.window_id) TO %s (FORMAT PARQUET, COMPRESSION ZSTD)"), q(output)))

written <- DBI::dbGetQuery(con, sprintf(paste0(
  "SELECT count(*) AS window_rows, count(DISTINCT sample_id) AS children, ",
  "count(*) FILTER (WHERE sites = 0) AS empty_child_windows FROM read_parquet(%s)"), q(output)))
if (written$window_rows != length(children) * nrow(windows) ||
    written$children != length(children)) {
  stop("chr20 truth windows do not cover every child and window", call. = FALSE)
}
print(written)
# The registry records this digest as the artifact's supplier identity.
cat(system2("sha256sum", shQuote(output), stdout = TRUE), "\n")
DBI::dbDisconnect(con, shutdown = TRUE)
