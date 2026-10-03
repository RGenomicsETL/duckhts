# The installed Rduckhts; load_all() would compile r/Rduckhts in place.
library(Rduckhts)
cache <- file.path(duckhtsbench::duckhts_bench_cache_dir(), "benchmarks", "roh-ancestry")
reference_path <- duckhtsbench::duckhts_bench_artifact_path("ancestry_reference_grch38_parquet")
correction_path <- normalizePath("test/data/ancestry_correction_bigsnpr_1.12.21.tsv")
bcf_path <- file.path(cache, "chr20.children.bcf")
con <- rduckhts_connect(extension_path = normalizePath("build/release/duckhts.duckdb_extension"),
                        config = list(memory_limit = "12GB", threads = "4",
                                      temp_directory = file.path(cache, "duckdb-tmp")))
q <- function(x) as.character(DBI::dbQuoteString(con, x))
DBI::dbExecute(con, sprintf(
  "CREATE VIEW ref20 AS SELECT * FROM read_parquet(%s) WHERE try_cast(regexp_replace(CAST(chromosome AS VARCHAR), '^chr', '') AS INTEGER) = 20",
  q(reference_path)))
DBI::dbExecute(con, sprintf(
  "CREATE VIEW ancestry_correction AS SELECT * FROM read_csv(%s, delim='\\t', header=true)",
  q(correction_path)))
DBI::dbExecute(con, sprintf(paste0(
  "CREATE TABLE child_dosage AS SELECT b.SAMPLE_ID AS sample_id, b.CHROM AS chromosome, b.POS AS position, b.REF AS allele_a, b.ALT[1] AS allele_b, ",
  "CASE WHEN regexp_matches(b.FORMAT_GT, '^[01][/|][01]$') THEN ",
  "(length(b.FORMAT_GT) - length(replace(b.FORMAT_GT, '1', '')))::DOUBLE END AS dosage ",
  "FROM read_bcf(%s, tidy_format := true) b JOIN ref20 r ON b.POS = r.position AND ",
  "((b.REF = r.allele_a AND b.ALT[1] = r.allele_b) OR (b.REF = r.allele_b AND b.ALT[1] = r.allele_a)) ",
  "WHERE len(b.ALT) = 1 AND NOT ((r.allele_a IN ('A','T') AND r.allele_b IN ('A','T')) OR ",
  "(r.allele_a IN ('C','G') AND r.allele_b IN ('C','G')))"), q(bcf_path)))
cat("Dosage input rows:", DBI::dbGetQuery(con, "SELECT count(*) AS n FROM child_dosage")$n, "\n")
proportions <- rduckhts_ancestry_proportions(con, "child_dosage", "ref20",
  "ancestry_correction", input_kind = "dosage", min_cor = 0.4)
utils::write.csv(proportions, file.path(cache, "chr20.q.csv"), row.names = FALSE)
cat("Per-child q status and top groups:\n")
for (sample in unique(proportions$sample_id)) {
  x <- proportions[proportions$sample_id == sample, ]
  cat(sample, unique(x$status), unique(x$used_variants), unique(x$cor_pred),
      paste(utils::head(x$group_id[order(-x$proportion)], 3), collapse = ","),
      paste(format(utils::head(sort(x$proportion, decreasing = TRUE), 3), digits = 4), collapse = ","), "\n")
}
DBI::dbDisconnect(con, shutdown = TRUE)
